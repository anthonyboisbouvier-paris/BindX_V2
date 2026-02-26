"""
DockIt pipeline -- Scaffold analysis and R-group decomposition.

Identifies the Murcko scaffold core, BRICS retrosynthetic bonds, and
secondary modifiable positions (halogens, O/S, aromatic C-H).  Returns
annotated SVG and per-position metadata for the frontend ScaffoldAnalyzer.

Uses RDKit where available; falls back to a minimal stub otherwise.
"""

from __future__ import annotations

import logging
import re
from typing import Optional

logger = logging.getLogger(__name__)


def analyze_scaffold(smiles: str) -> dict:
    """Analyze a molecule's scaffold, R-group positions and return annotated SVG.

    Parameters
    ----------
    smiles : str
        Input SMILES string.

    Returns
    -------
    dict
        Keys: smiles, scaffold_smiles, positions (list), annotated_svg,
        core_atom_indices, brics_bond_count, stats.
    """
    try:
        from rdkit import Chem
        return _analyze_rdkit(smiles)
    except ImportError:
        logger.warning("RDKit not available; returning empty scaffold analysis")
        return _fallback(smiles)


def _fallback(smiles: str) -> dict:
    """Return a minimal response when RDKit is unavailable."""
    return {
        "smiles": smiles,
        "scaffold_smiles": None,
        "positions": [],
        "annotated_svg": '<svg xmlns="http://www.w3.org/2000/svg" width="300" height="200">'
                         '<text x="50%" y="50%" text-anchor="middle" fill="#888" font-size="14">'
                         'RDKit required</text></svg>',
        "core_atom_indices": [],
        "brics_bond_count": 0,
        "stats": {},
    }


# ---------------------------------------------------------------------------
# RDKit implementation
# ---------------------------------------------------------------------------

def _analyze_rdkit(smiles: str) -> dict:
    """Full scaffold analysis using RDKit."""
    from rdkit import Chem
    from rdkit.Chem import Descriptors, RingInfo
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.Chem import BRICS

    if not smiles or not smiles.strip():
        return _fallback(smiles)

    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        return _fallback(smiles)

    # --- Murcko scaffold ---
    try:
        core_mol = MurckoScaffold.GetScaffoldForMol(mol)
        scaffold_smiles = Chem.MolToSmiles(core_mol)
        # Map core atoms: match scaffold back onto mol
        core_match = mol.GetSubstructMatch(core_mol)
        core_set = set(core_match) if core_match else set()
    except Exception:
        scaffold_smiles = None
        core_set = set()

    # --- BRICS bonds ---
    brics_bonds = []
    brics_atom_set = set()
    try:
        raw_brics = BRICS.FindBRICSBonds(mol)
        for (i, j), _ in raw_brics:
            brics_bonds.append((i, j))
            brics_atom_set.add(i)
            brics_atom_set.add(j)
    except Exception:
        pass

    # --- Build position list ---
    positions = []
    label_counter = 1

    # 1) BRICS bond sites
    for (begin_idx, end_idx) in brics_bonds:
        # The "R-group side" is the atom NOT in the core
        if begin_idx in core_set and end_idx not in core_set:
            rgroup_idx, anchor_idx = end_idx, begin_idx
        elif end_idx in core_set and begin_idx not in core_set:
            rgroup_idx, anchor_idx = begin_idx, end_idx
        else:
            # Both in or both out — pick the one with fewer neighbors
            a1 = mol.GetAtomWithIdx(begin_idx)
            a2 = mol.GetAtomWithIdx(end_idx)
            if a1.GetDegree() <= a2.GetDegree():
                rgroup_idx, anchor_idx = begin_idx, end_idx
            else:
                rgroup_idx, anchor_idx = end_idx, begin_idx

        atom = mol.GetAtomWithIdx(rgroup_idx)
        fragment_smi = _get_fragment_at_bond(mol, anchor_idx, rgroup_idx, core_set)
        strategies = _strategies_for_atom(atom)
        suggestions = _suggested_groups(atom, strategies)

        positions.append({
            "label": f"R{label_counter}",
            "position_idx": rgroup_idx,
            "atom_symbol": atom.GetSymbol(),
            "current_group": fragment_smi,
            "applicable_strategies": strategies,
            "suggested_replacements": suggestions,
            "is_brics_site": True,
        })
        label_counter += 1

    # 2) Secondary modifiable positions (non-BRICS)
    seen_indices = {p["position_idx"] for p in positions}
    secondary = _classify_modifiable_atoms(mol, core_set, brics_atom_set)
    for info in secondary:
        if info["idx"] in seen_indices:
            continue
        seen_indices.add(info["idx"])
        atom = mol.GetAtomWithIdx(info["idx"])
        strategies = info["strategies"]
        suggestions = _suggested_groups(atom, strategies)

        positions.append({
            "label": f"R{label_counter}",
            "position_idx": info["idx"],
            "atom_symbol": atom.GetSymbol(),
            "current_group": atom.GetSymbol(),
            "applicable_strategies": strategies,
            "suggested_replacements": suggestions,
            "is_brics_site": False,
        })
        label_counter += 1

    # --- Annotated SVG ---
    annotated_svg = _generate_annotated_svg(mol, core_set, positions, brics_bonds)

    # --- Stats ---
    stats = {
        "n_positions": len(positions),
        "n_brics_bonds": len(brics_bonds),
        "scaffold_mw": round(Descriptors.ExactMolWt(core_mol), 1) if scaffold_smiles else None,
        "n_rings": mol.GetRingInfo().NumRings(),
    }

    return {
        "smiles": smiles,
        "scaffold_smiles": scaffold_smiles,
        "positions": positions,
        "annotated_svg": annotated_svg,
        "core_atom_indices": sorted(core_set),
        "brics_bond_count": len(brics_bonds),
        "stats": stats,
    }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _get_fragment_at_bond(mol, anchor_idx: int, rgroup_idx: int, core_set: set) -> str:
    """Cut bond between anchor and rgroup_idx, return SMILES of the R-group side."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        bond = mol.GetBondBetweenAtoms(anchor_idx, rgroup_idx)
        if bond is None:
            return mol.GetAtomWithIdx(rgroup_idx).GetSymbol()

        frags = Chem.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=True)
        frag_smiles = Chem.MolToSmiles(frags)
        parts = frag_smiles.split(".")

        # Return the fragment that does NOT contain the majority of core atoms
        if len(parts) >= 2:
            # Pick the smaller fragment (likely the R-group)
            parts.sort(key=len)
            # Remove dummy atom markers [*] or [1*] etc
            frag = re.sub(r'\[\d*\*\]', '[H]', parts[0])
            return frag if frag else parts[0]
        return frag_smiles
    except Exception:
        return mol.GetAtomWithIdx(rgroup_idx).GetSymbol()


def _classify_modifiable_atoms(mol, core_set: set, brics_atoms: set) -> list[dict]:
    """Find secondary modifiable positions: terminal halogens, O/S, aromatic C-H."""
    from rdkit import Chem

    results = []
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in brics_atoms:
            continue

        anum = atom.GetAtomicNum()

        # Terminal halogens (F, Cl, Br)
        if anum in (9, 17, 35) and atom.GetDegree() == 1:
            results.append({"idx": idx, "strategies": ["swap_halogen"]})

        # O or S not in ring (potential swap targets)
        elif anum in (8, 16) and not atom.IsInRing() and atom.GetDegree() <= 2:
            results.append({"idx": idx, "strategies": ["swap_atom"]})

        # Aromatic C-H (could add functional group)
        elif (anum == 6 and atom.GetIsAromatic()
              and atom.GetTotalNumHs() > 0 and idx not in core_set):
            results.append({"idx": idx, "strategies": ["add_fg"]})

    return results


def _strategies_for_atom(atom) -> list[str]:
    """Determine applicable strategies for a given atom."""
    anum = atom.GetAtomicNum()
    strategies = []

    if anum in (9, 17, 35):
        strategies.append("swap_halogen")
    if anum in (8, 16):
        strategies.append("swap_atom")
    if anum == 6 and atom.GetIsAromatic():
        strategies.append("add_fg")
    if anum == 6 and not atom.GetIsAromatic():
        strategies.append("modify_chain")

    # All atoms can be targets of functional group addition in theory
    if "add_fg" not in strategies:
        strategies.append("add_fg")

    return strategies


def _suggested_groups(atom, strategies: list[str]) -> list[str]:
    """Suggest replacement groups based on atom type and strategies."""
    anum = atom.GetAtomicNum()
    groups = []

    if "swap_halogen" in strategies:
        groups.extend(["F", "Cl", "Br"])
    if "add_fg" in strategies:
        groups.extend(["C", "OC", "N", "C(F)(F)F", "C#N", "C(=O)N", "O"])
    if "swap_atom" in strategies:
        if anum == 8:
            groups.extend(["S", "N"])
        elif anum == 16:
            groups.extend(["O", "N"])
    if "modify_chain" in strategies:
        groups.extend(["C", "CC", "C(C)C"])

    # Deduplicate while preserving order
    seen = set()
    unique = []
    for g in groups:
        if g not in seen:
            seen.add(g)
            unique.append(g)
    return unique


def _generate_annotated_svg(
    mol,
    core_atoms: set,
    positions: list[dict],
    brics_bonds: list[tuple],
    width: int = 400,
    height: int = 300,
) -> Optional[str]:
    """Generate an SVG with core atoms highlighted orange, R-groups green, BRICS bonds red."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Chem.Draw import rdMolDraw2D

        # Compute 2D coords
        AllChem.Compute2DCoords(mol)

        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        opts = drawer.drawOptions()
        opts.addAtomIndices = False
        opts.addStereoAnnotation = False

        # Build highlight maps
        highlight_atoms = {}
        highlight_bonds = {}
        highlight_radii = {}

        # Core atoms in orange
        for idx in core_atoms:
            highlight_atoms[idx] = (0.85, 0.6, 0.2, 0.3)
            highlight_radii[idx] = 0.3

        # R-group positions in green
        position_indices = set()
        for pos in positions:
            pidx = pos["position_idx"]
            position_indices.add(pidx)
            highlight_atoms[pidx] = (0.13, 0.77, 0.37, 0.5)
            highlight_radii[pidx] = 0.4

        # BRICS bonds in red
        for (i, j) in brics_bonds:
            bond = mol.GetBondBetweenAtoms(i, j)
            if bond is not None:
                highlight_bonds[bond.GetIdx()] = (0.9, 0.2, 0.2, 0.8)

        # Draw with highlights
        atom_colors = {idx: color for idx, color in highlight_atoms.items()}
        bond_colors = {idx: color for idx, color in highlight_bonds.items()}

        drawer.DrawMolecule(
            mol,
            highlightAtoms=list(atom_colors.keys()),
            highlightAtomColors=atom_colors,
            highlightBonds=list(bond_colors.keys()),
            highlightBondColors=bond_colors,
            highlightAtomRadii=highlight_radii,
        )
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()

        # Post-process: inject R-group labels using drawer's 2D coordinates
        svg = _inject_labels(svg, mol, positions, drawer)

        return svg

    except Exception as exc:
        logger.warning("Annotated SVG generation failed: %s", exc)
        return None


def _inject_labels(svg: str, mol, positions: list[dict], drawer=None) -> str:
    """Inject R-group label text into the SVG XML using drawer's 2D coordinates."""
    try:
        if drawer is None:
            return svg

        labels_xml = []
        for pos in positions:
            idx = pos["position_idx"]
            label = pos["label"]

            # Use the drawer's coordinate mapping (molecule-space -> SVG-space)
            try:
                point = drawer.GetDrawCoords(idx)
                x = point.x
                y = point.y
            except Exception:
                continue

            # Offset the label slightly above and to the right of the atom
            lx = x + 12
            ly = y - 12

            labels_xml.append(
                f'<g>'
                f'<rect x="{lx - 2}" y="{ly - 10}" width="{len(label) * 7 + 4}" height="14" rx="3" '
                f'fill="white" stroke="#1e3a5f" stroke-width="0.8" opacity="0.9"/>'
                f'<text x="{lx}" y="{ly}" font-size="10" fill="#1e3a5f" font-weight="bold"'
                f' font-family="Arial,sans-serif" dominant-baseline="auto"'
                f' data-atom-idx="{idx}" class="rgroup-label">{label}</text>'
                f'</g>'
            )

        if labels_xml:
            labels_group = '<g class="rgroup-labels">' + ''.join(labels_xml) + '</g>'
            svg = svg.replace('</svg>', labels_group + '</svg>')

        return svg
    except Exception:
        return svg
