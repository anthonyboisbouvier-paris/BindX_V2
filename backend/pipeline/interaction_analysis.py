"""
DockIt pipeline -- V5bis Protein-ligand interaction analysis using ProLIF.

Analyzes specific interactions (H-bonds, pi-stacking, hydrophobic, etc.)
between docked ligands and the target protein. Crosses with functional
residues from UniProt to assess binding quality.

Strategies (in order):
  1. ProLIF for interaction fingerprinting (requires prolif + MDAnalysis).
  2. RDKit-based distance analysis from PDB coordinates.
  3. Deterministic mock based on known target data.
"""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Known functional residues for common drug targets
# ---------------------------------------------------------------------------

KNOWN_FUNCTIONAL_RESIDUES: dict[str, dict] = {
    "P00533": {  # EGFR
        "residues": [719, 721, 743, 745, 790, 791, 793, 797, 854, 855],
        "key_hbond_residues": [793, 790],  # MET793 (hinge), THR790 (gatekeeper)
        "description": "ATP-binding site key residues",
    },
    "P04626": {  # HER2/ERBB2
        "residues": [726, 753, 798, 862, 863],
        "key_hbond_residues": [798],
        "description": "ATP-binding site key residues",
    },
    "P07550": {  # ADRB2
        "residues": [113, 203, 207, 286, 289, 290, 312, 316],
        "key_hbond_residues": [113, 207],
        "description": "Orthosteric binding site residues",
    },
    "P00519": {  # ABL1
        "residues": [253, 255, 271, 286, 315, 318, 381, 382],
        "key_hbond_residues": [318, 315],
        "description": "ATP-binding site key residues",
    },
}


# ---------------------------------------------------------------------------
# UniProt -> PDB numbering offset helpers
# ---------------------------------------------------------------------------

def _get_uniprot_pdb_offset(pdb_path: Path, uniprot_id: str) -> int:
    """Parse PDB DBREF records to find numbering offset between UniProt and PDB.

    The DBREF record maps PDB residue numbers to database (UniProt) residue
    numbers.  The offset is ``pdb_start - db_start``, so that
    ``uniprot_residue + offset == pdb_residue``.

    Parameters
    ----------
    pdb_path : Path
        Path to the protein PDB file.
    uniprot_id : str
        UniProt accession (e.g. ``"P00533"``).

    Returns
    -------
    int
        Offset to add to UniProt residue numbers to obtain PDB residue
        numbers.  Returns 0 if no DBREF record is found.
    """
    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("DBREF") and uniprot_id.upper() in line.upper():
                    # DBREF format (PDB standard):
                    #   cols 15-18 = PDB sequence begin
                    #   cols 55-60 = database sequence begin
                    try:
                        pdb_start = int(line[14:18].strip())
                        db_start = int(line[55:60].strip())
                        offset = pdb_start - db_start
                        logger.info(
                            "DBREF offset for %s: pdb_start=%d, db_start=%d, offset=%d",
                            uniprot_id, pdb_start, db_start, offset,
                        )
                        return offset
                    except (ValueError, IndexError):
                        continue
    except Exception as exc:
        logger.debug("Could not read DBREF from %s: %s", pdb_path, exc)
    return 0  # No offset found, assume same numbering


def _map_functional_residues(known_residues: list[int], offset: int) -> set[int]:
    """Map UniProt residue numbers to PDB numbering.

    Parameters
    ----------
    known_residues : list[int]
        Residue numbers in UniProt numbering.
    offset : int
        Offset to add (``pdb_start - db_start``).

    Returns
    -------
    set[int]
        Residue numbers in PDB numbering.
    """
    return {r + offset for r in known_residues}


# ---------------------------------------------------------------------------
# VDW redundancy filter (PLIP / Discovery Studio convention)
# ---------------------------------------------------------------------------

_VDW_TYPES = {"VDW", "VdWContact"}


def _filter_redundant_vdw(interactions: list[dict]) -> list[dict]:
    """Remove VDW entries when the same residue has a more specific interaction.

    Follows PLIP / Discovery Studio convention: VDW contacts are only
    reported for residues where no specific interaction (H-bond, hydrophobic,
    pi-stacking, ionic, etc.) is detected.  This avoids noise in the table
    and keeps the residue count consistent with the interaction diagram.
    """
    has_specific: set[int] = set()
    for i in interactions:
        if i.get("type") not in _VDW_TYPES:
            has_specific.add(i.get("residue_number"))
    return [
        i for i in interactions
        if i.get("type") not in _VDW_TYPES or i.get("residue_number") not in has_specific
    ]


def _compute_quality_metrics(
    interactions: list[dict],
    uniprot_id: str,
    pdb_path: Optional[str] = None,
    functional_residues: Optional[dict] = None,
) -> dict:
    """Compute interaction quality metrics with offset mapping and fuzzy matching.

    This is the shared quality computation used by all three analysis
    strategies.  It handles:

    1. UniProt -> PDB numbering offset via DBREF parsing.
    2. Exact residue matching against mapped functional residues.
    3. Fuzzy matching (+-2 residues) if exact matching yields 0.
    4. Safety net: if there are interactions but no functional contacts,
       assign a floor quality proportional to the number of interactions.
    5. Re-tagging ``is_functional`` on each interaction dict.

    Parameters
    ----------
    interactions : list[dict]
        Interaction dicts, each containing at least ``residue_number``.
    uniprot_id : str
        UniProt accession for functional residue lookup.
    pdb_path : str | None
        Path to the protein PDB file (for DBREF offset).  May be ``None``
        for the mock strategy.
    functional_residues : dict | None
        External functional residues from project target_preview.
        Prioritized over the hardcoded KNOWN_FUNCTIONAL_RESIDUES dict.

    Returns
    -------
    dict
        Keys: ``functional_contacts``, ``total_functional``,
        ``interaction_quality``, ``key_hbonds``.
    """
    # Roles that are biologically meaningful for quality metrics.
    # Exclude "pocket" (P2Rank geometry) and "site" (too generic UniProt).
    _QUALITY_ROLES = {"active_site", "binding", "metal_binding", "custom", "key"}

    # Priority: external functional_residues > hardcoded KNOWN_FUNCTIONAL_RESIDUES
    if functional_residues and functional_residues.get("residues"):
        known_res_uniprot = [
            r["number"] if isinstance(r, dict) else r
            for r in functional_residues["residues"]
            if (
                (isinstance(r, dict) and r.get("enabled", True)
                 and r.get("role", "key") in _QUALITY_ROLES)
                or not isinstance(r, dict)
            )
        ]
        key_hbond_uniprot = functional_residues.get("key_hbond_residues", [])
    else:
        functional = KNOWN_FUNCTIONAL_RESIDUES.get(uniprot_id, {})
        known_res_uniprot = functional.get("residues", [])
        key_hbond_uniprot = functional.get("key_hbond_residues", [])

    # --- Compute offset ---
    offset = 0
    if pdb_path:
        offset = _get_uniprot_pdb_offset(Path(pdb_path), uniprot_id)

    # --- Map residues to PDB numbering ---
    mapped_residues = _map_functional_residues(known_res_uniprot, offset)
    mapped_key_hbond = _map_functional_residues(key_hbond_uniprot, offset)

    # --- Build PDB-numbered → info mapping for origin tagging ---
    mapped_to_info: dict[int, dict] = {}
    if functional_residues and functional_residues.get("residues"):
        for r in functional_residues["residues"]:
            if isinstance(r, dict) and r.get("enabled", True) and r.get("role", "key") in _QUALITY_ROLES:
                pdb_num = r["number"] + offset
                mapped_to_info[pdb_num] = {
                    "role": r.get("role", "key"),
                    "description": r.get("description", ""),
                    "aa": r.get("aa", ""),
                    "number": pdb_num,
                    "uniprot_number": r["number"],
                }

    # --- Build fuzzy set (+-2 around each mapped residue) ---
    fuzzy_residues: set[int] = set()
    for r in mapped_residues:
        for delta in range(-2, 3):
            fuzzy_residues.add(r + delta)

    # --- Detected residue numbers ---
    detected_residue_numbers = set(i["residue_number"] for i in interactions)

    # --- Exact matching ---
    func_contacted = len(detected_residue_numbers & mapped_residues)

    # --- Fuzzy fallback ---
    if func_contacted == 0 and mapped_residues:
        func_contacted = len(detected_residue_numbers & fuzzy_residues)
        if func_contacted > 0:
            logger.info(
                "Fuzzy matching found %d functional contacts (exact was 0)",
                func_contacted,
            )

    total_func = len(mapped_residues) if mapped_residues else 1
    quality = func_contacted / total_func if total_func > 0 else 0.5

    # --- Safety net ---
    n_interactions = len(interactions)
    if quality == 0.0 and n_interactions > 0:
        quality = min(n_interactions / 10.0, 0.5)
        logger.info(
            "Safety net: %d interactions but 0 functional contacts -> quality=%.2f",
            n_interactions, quality,
        )

    # --- Key H-bonds ---
    key_hb = sum(
        1 for i in interactions
        if i["residue_number"] in mapped_key_hbond and "HB" in i.get("type", "")
    )

    # Also count via fuzzy key H-bond set
    if key_hb == 0 and mapped_key_hbond:
        fuzzy_key_hbond: set[int] = set()
        for r in mapped_key_hbond:
            for delta in range(-2, 3):
                fuzzy_key_hbond.add(r + delta)
        key_hb = sum(
            1 for i in interactions
            if i["residue_number"] in fuzzy_key_hbond and "HB" in i.get("type", "")
        )

    # --- Re-tag is_functional + origin on interactions ---
    # Use the same set that was used for functional_contacts counting
    # so that checkmark count == KPI number
    func_set = mapped_residues if func_contacted > 0 else fuzzy_residues
    for i in interactions:
        i["is_functional"] = i["residue_number"] in func_set
        info = mapped_to_info.get(i["residue_number"])
        if info:
            i["origin_role"] = info["role"]
            i["origin_label"] = info.get("description") or info["role"]

    return {
        "functional_contacts": func_contacted,
        "total_functional": len(known_res_uniprot),
        "interaction_quality": round(quality, 3),
        "key_hbonds": key_hb,
        "mapped_functional": list(mapped_to_info.values()),
    }


# =====================================================================
# PUBLIC API
# =====================================================================

def analyze_interactions(
    protein_path: str,
    ligand_path: str,
    uniprot_id: str = "",
    smiles: str = "",
    functional_residues: Optional[dict] = None,
) -> dict:
    """Analyze protein-ligand interactions.

    Tries ProLIF first, falls back to RDKit distance-based analysis,
    and finally to a deterministic mock.

    Parameters
    ----------
    protein_path : str
        Path to the protein PDB file.
    ligand_path : str
        Path to the docked ligand file (PDB, SDF, or MOL2).
    uniprot_id : str
        UniProt accession (e.g. "P00533") for functional residue lookup.
    smiles : str
        SMILES of the ligand (used for deterministic mock fallback).
    functional_residues : dict | None
        External functional residues from project target_preview.

    Returns
    -------
    dict
        Interaction analysis result with keys:
        - ``interactions``: list of interaction dicts.
        - ``functional_contacts``: int, number of functional residues contacted.
        - ``total_functional``: int, total known functional residues.
        - ``interaction_quality``: float in [0, 1].
        - ``key_hbonds``: int, number of key H-bond contacts.
        - ``method``: "prolif", "rdkit_distance", or "mock".
        - ``summary``: human-readable summary string.
    """
    # Try ProLIF
    try:
        return _analyze_prolif(protein_path, ligand_path, uniprot_id, functional_residues=functional_residues)
    except ImportError:
        logger.info("ProLIF not available, using RDKit distance analysis")
    except Exception as exc:
        logger.warning("ProLIF failed: %s, falling back", exc)

    # Try RDKit distance-based
    try:
        return _analyze_rdkit_distance(protein_path, ligand_path, uniprot_id, functional_residues=functional_residues)
    except ImportError:
        raise RuntimeError(
            "Interaction analysis requires RDKit. Install rdkit-pypi."
        )
    except Exception as exc:
        raise RuntimeError(
            f"Interaction analysis failed for {uniprot_id}: {exc}"
        ) from exc


def compute_pose_quality(
    protein_path: str,
    ligand_path: str,
    uniprot_id: str = "",
    smiles: str = "",
    functional_residues: Optional[dict] = None,
) -> dict:
    """Compute structured pose quality metrics from interaction analysis.

    Wraps ``analyze_interactions()`` and extracts a structured dict of
    quality metrics suitable for the ``PoseQuality`` Pydantic model.

    Parameters
    ----------
    protein_path : str
        Path to the protein PDB file.
    ligand_path : str
        Path to the docked ligand file (PDB, SDF, PDBQT).
    uniprot_id : str
        UniProt accession for functional residue lookup.
    smiles : str
        SMILES of the ligand (for mock fallback).
    functional_residues : dict | None
        External functional residues from project target_preview.

    Returns
    -------
    dict
        Keys: ``n_contacts_4A``, ``n_hbonds``, ``has_clashes``, ``n_clashes``,
        ``key_residue_distances``, ``interaction_quality``.
    """
    raw = analyze_interactions(protein_path, ligand_path, uniprot_id, smiles, functional_residues=functional_residues)
    interactions = raw.get("interactions", [])

    contacts_4A = 0
    hbonds = 0
    clashes = 0

    for itype in interactions:
        dist = itype.get("distance", 3.5)  # default for methods without distance
        if dist < 4.0:
            contacts_4A += 1
        if dist < 1.5:
            clashes += 1
        if itype.get("type", "") in ("HBDonor", "HBAcceptor", "Hbond"):
            hbonds += 1

    # Key residue distances from known functional residues
    key_residue_distances: dict[str, float] = {}
    functional = KNOWN_FUNCTIONAL_RESIDUES.get(uniprot_id, {})
    known_res = functional.get("residues", [])
    if known_res:
        # Build a map from residue_number -> min distance
        res_min_dist: dict[int, float] = {}
        for i in interactions:
            rn = i.get("residue_number", 0)
            d = i.get("distance", 3.5)
            if rn not in res_min_dist or d < res_min_dist[rn]:
                res_min_dist[rn] = d

        # Common residue name lookup (simplified)
        aa3 = {
            "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
            "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
            "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
            "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        }
        # Build residue label from interaction data
        res_labels: dict[int, str] = {}
        for i in interactions:
            rn = i.get("residue_number", 0)
            res_name = i.get("residue", "")
            if rn and res_name:
                res_labels[rn] = res_name

        for rn in known_res:
            label = res_labels.get(rn, f"RES{rn}")
            if rn in res_min_dist:
                key_residue_distances[label] = round(res_min_dist[rn], 1)

    return {
        "n_contacts_4A": contacts_4A,
        "n_hbonds": hbonds,
        "has_clashes": clashes > 0,
        "n_clashes": clashes,
        "key_residue_distances": key_residue_distances if key_residue_distances else None,
        "interaction_quality": raw.get("interaction_quality", 0.0),
    }


def score_interaction_quality(interactions_result: dict) -> float:
    """Extract interaction_quality score (0-1) for composite scoring.

    Parameters
    ----------
    interactions_result : dict
        Output from ``analyze_interactions()``.

    Returns
    -------
    float
        Interaction quality score in [0, 1].
    """
    return interactions_result.get("interaction_quality", 0.5)


# =====================================================================
# STRATEGY 1: ProLIF
# =====================================================================

def _analyze_prolif(
    protein_path: str,
    ligand_path: str,
    uniprot_id: str,
    functional_residues: Optional[dict] = None,
) -> dict:
    """Use ProLIF for interaction fingerprinting.

    Requires ``prolif`` and ``MDAnalysis`` to be installed.

    Parameters
    ----------
    protein_path : str
        Path to protein PDB.
    ligand_path : str
        Path to docked ligand file.
    uniprot_id : str
        UniProt accession for functional residue annotation.
    functional_residues : dict | None
        External functional residues from project target_preview.

    Returns
    -------
    dict
        Interaction analysis result.
    """
    import prolif  # type: ignore[import-untyped]
    import MDAnalysis as mda  # type: ignore[import-untyped]

    # GNINA output files may have .pdbqt extension but contain SDF/MOL or
    # real PDBQT data.  MDAnalysis cannot parse either format directly for
    # ProLIF, so convert to standard PDB via Open Babel first.
    import subprocess as _sp

    actual_ligand_path = ligand_path
    needs_conversion = Path(ligand_path).suffix.lower() in (".pdbqt", ".sdf", ".mol2")

    if needs_conversion or Path(ligand_path).suffix.lower() == ".pdb":
        # Detect actual format by inspecting first line
        try:
            with open(ligand_path, "r") as _fh:
                first_lines = _fh.read(500)
        except Exception:
            first_lines = ""

        if "V2000" in first_lines or "V3000" in first_lines or "$$$$" in first_lines:
            in_fmt = "sdf"
        elif "BRANCH" in first_lines or "ENDROOT" in first_lines or "TORSDOF" in first_lines:
            in_fmt = "pdbqt"
        elif "ATOM" in first_lines or "HETATM" in first_lines:
            in_fmt = "pdb"
            needs_conversion = False  # already PDB
        else:
            in_fmt = "pdbqt"  # default guess

    if needs_conversion:
        pdb_tmp = Path(ligand_path).with_suffix(".prolif.pdb")
        try:
            _sp.run(
                ["obabel", f"-i{in_fmt}", ligand_path, "-opdb", "-O", str(pdb_tmp), "-h"],
                capture_output=True, text=True, timeout=30, check=True,
            )
            if pdb_tmp.stat().st_size > 0:
                actual_ligand_path = str(pdb_tmp)
                logger.info("Converted %s -> PDB for ProLIF: %s", in_fmt.upper(), pdb_tmp)
            else:
                logger.warning("Conversion produced empty file, trying generic detection")
                _sp.run(
                    ["obabel", ligand_path, "-opdb", "-O", str(pdb_tmp), "-h"],
                    capture_output=True, text=True, timeout=30, check=True,
                )
                if pdb_tmp.stat().st_size > 0:
                    actual_ligand_path = str(pdb_tmp)
                    logger.info("Generic obabel conversion succeeded: %s", pdb_tmp)
                else:
                    raise RuntimeError("obabel produced empty PDB")
        except Exception as conv_exc:
            logger.warning("Ligand->PDB conversion failed: %s", conv_exc)
            raise

    # Ensure protein PDB has hydrogens (ProLIF requires explicit H)
    protein_h_path = protein_path
    prot_pdb = Path(protein_path)
    if prot_pdb.suffix.lower() == ".pdb":
        prot_h_tmp = prot_pdb.with_suffix(".prolif_h.pdb")
        if not prot_h_tmp.exists():
            try:
                _sp.run(
                    ["obabel", "-ipdb", protein_path, "-opdb", "-O", str(prot_h_tmp), "-h"],
                    capture_output=True, text=True, timeout=60, check=True,
                )
                if prot_h_tmp.stat().st_size > 0:
                    protein_h_path = str(prot_h_tmp)
                    logger.info("Added H to protein for ProLIF: %s", prot_h_tmp)
            except Exception:
                logger.warning("Could not add H to protein, using original")
        elif prot_h_tmp.stat().st_size > 0:
            protein_h_path = str(prot_h_tmp)

    # Load protein and ligand
    prot = mda.Universe(protein_h_path)
    lig = mda.Universe(actual_ligand_path)

    prot_mol = prolif.Molecule.from_mda(prot, force=True)
    lig_mol = prolif.Molecule.from_mda(lig, force=True)

    # Compute interaction fingerprint
    fp = prolif.Fingerprint()
    fp.run_from_iterable([lig_mol], prot_mol)

    # Extract interactions (is_functional will be re-tagged by _compute_quality_metrics)
    interactions: list[dict] = []
    ifp = fp.ifp[0]  # First (only) frame

    # Compute per-residue min heavy-atom distance (ligand ↔ protein)
    # ProLIF adds hydrogens, but we report heavy-atom distances (convention)
    _residue_distances: dict[int, float] = {}
    try:
        # Filter to heavy atoms only (exclude H) for distance calculation
        lig_heavy = lig.select_atoms("not element H")
        lig_positions = lig_heavy.positions if len(lig_heavy) > 0 else lig.atoms.positions
        for (lig_res, prot_res), interaction_dict in ifp.items():
            if not any(interaction_dict.values()):
                continue
            rn = prot_res.number
            if rn in _residue_distances:
                continue
            try:
                res_heavy = prot.select_atoms(f"resid {rn} and not element H")
                res_atoms = res_heavy if len(res_heavy) > 0 else prot.select_atoms(f"resid {rn}")
                if len(res_atoms) > 0:
                    prot_positions = res_atoms.positions  # (M, 3)
                    # Pairwise distances: (N, M) — find min
                    diff = lig_positions[:, None, :] - prot_positions[None, :, :]
                    dists = (diff ** 2).sum(axis=2) ** 0.5
                    _residue_distances[rn] = round(float(dists.min()), 2)
            except Exception:
                pass
    except Exception:
        logger.debug("Could not compute ProLIF residue distances")

    for (lig_res, prot_res), interaction_dict in ifp.items():
        res_num = prot_res.number
        res_name = f"{prot_res.name}{res_num}"

        for int_type, is_present in interaction_dict.items():
            if is_present:
                entry: dict = {
                    "residue": res_name,
                    "residue_number": res_num,
                    "type": int_type,
                    "is_functional": False,  # re-tagged below
                }
                if res_num in _residue_distances:
                    entry["distance"] = _residue_distances[res_num]
                interactions.append(entry)

    # Filter redundant VDW (PLIP convention: only keep VDW when no specific interaction)
    interactions = _filter_redundant_vdw(interactions)

    # Compute quality with offset mapping + fuzzy matching + safety net
    metrics = _compute_quality_metrics(interactions, uniprot_id, pdb_path=protein_path, functional_residues=functional_residues)

    return {
        "interactions": interactions,
        **metrics,
        "method": "prolif",
        "summary": (
            f"Functional contacts: {metrics['functional_contacts']}/{metrics['total_functional']} "
            f"key residues ({metrics['interaction_quality']:.0%}) - "
            f"{metrics['key_hbonds']} key H-bonds"
        ),
    }


# =====================================================================
# STRATEGY 2: RDKit distance-based analysis
# =====================================================================

def _analyze_rdkit_distance(
    protein_path: str,
    ligand_path: str,
    uniprot_id: str,
    functional_residues: Optional[dict] = None,
) -> dict:
    """Simplified interaction analysis using RDKit and coordinate distances.

    Parses PDB files, finds contacts within 4 Angstroms, and classifies
    them by atom types.

    Parameters
    ----------
    protein_path : str
        Path to protein PDB.
    ligand_path : str
        Path to docked ligand PDB/SDF.
    uniprot_id : str
        UniProt accession for functional residue annotation.
    functional_residues : dict | None
        External functional residues from project target_preview.

    Returns
    -------
    dict
        Interaction analysis result.
    """
    from rdkit import Chem  # type: ignore[import-untyped]

    # Parse protein PDB to extract atom coordinates (heavy atoms only for distances)
    protein_atoms = [a for a in _parse_pdb_atoms(protein_path) if a["element"].upper() not in ("H", "D")]
    if not protein_atoms:
        raise ValueError(f"Could not parse protein atoms from {protein_path}")

    # Parse ligand atoms (heavy atoms only)
    ligand_atoms = [a for a in _parse_ligand_atoms(ligand_path) if a["element"].upper() not in ("H", "D")]
    if not ligand_atoms:
        raise ValueError(f"Could not parse ligand atoms from {ligand_path}")

    # Find contacts within 4.0 Angstroms
    contact_distance = 4.0
    interactions: list[dict] = []
    seen_contacts: set[tuple[int, str]] = set()

    for lig_atom in ligand_atoms:
        for prot_atom in protein_atoms:
            dist = _euclidean_distance(lig_atom["coords"], prot_atom["coords"])
            if dist <= contact_distance:
                res_num = prot_atom["res_num"]
                res_name = prot_atom["res_name"]
                contact_key = (res_num, prot_atom["atom_name"])
                if contact_key in seen_contacts:
                    continue
                seen_contacts.add(contact_key)

                # Classify interaction type
                int_type = _classify_interaction(
                    lig_atom, prot_atom, dist,
                )

                interactions.append({
                    "residue": f"{res_name}{res_num}",
                    "residue_number": res_num,
                    "type": int_type,
                    "distance": round(dist, 2),
                    "is_functional": False,  # re-tagged below
                })

    # Deduplicate by (residue, type) pair — keep closest contact per pair
    # This preserves multiple interaction types per residue (like ProLIF/PLIP)
    pair_best: dict[tuple[int, str], dict] = {}
    for interaction in interactions:
        key = (interaction["residue_number"], interaction["type"])
        if key not in pair_best or interaction.get("distance", 99) < pair_best[key].get("distance", 99):
            pair_best[key] = interaction
    interactions = list(pair_best.values())

    # Filter redundant VDW (PLIP convention: only keep VDW when no specific interaction)
    interactions = _filter_redundant_vdw(interactions)

    # Compute quality with offset mapping + fuzzy matching + safety net
    metrics = _compute_quality_metrics(interactions, uniprot_id, pdb_path=protein_path, functional_residues=functional_residues)

    return {
        "interactions": interactions,
        **metrics,
        "method": "rdkit_distance",
        "summary": (
            f"Functional contacts: {metrics['functional_contacts']}/{metrics['total_functional']} "
            f"key residues ({metrics['interaction_quality']:.0%}) - "
            f"{metrics['key_hbonds']} key H-bonds"
        ),
    }


def _parse_pdb_atoms(pdb_path: str) -> list[dict]:
    """Parse ATOM records from a PDB file.

    Parameters
    ----------
    pdb_path : str
        Path to the PDB file.

    Returns
    -------
    list[dict]
        List of atom dicts with keys: atom_name, element, res_name,
        res_num, coords.
    """
    atoms: list[dict] = []
    try:
        with open(pdb_path, "r") as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        atom_name = line[12:16].strip()
                        res_name = line[17:20].strip()
                        res_num = int(line[22:26].strip())
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        element = line[76:78].strip() if len(line) > 76 else atom_name[0]
                        atoms.append({
                            "atom_name": atom_name,
                            "element": element,
                            "res_name": res_name,
                            "res_num": res_num,
                            "coords": (x, y, z),
                        })
                    except (ValueError, IndexError):
                        continue
    except Exception as exc:
        logger.warning("Failed to parse PDB %s: %s", pdb_path, exc)
    return atoms


def _parse_ligand_atoms(ligand_path: str) -> list[dict]:
    """Parse atoms from a ligand file (PDB or SDF).

    Parameters
    ----------
    ligand_path : str
        Path to the ligand file.

    Returns
    -------
    list[dict]
        List of atom dicts with keys: atom_name, element, coords.
    """
    path = Path(ligand_path)
    if path.suffix.lower() in (".pdb", ".pdbqt"):
        return _parse_pdb_atoms(ligand_path)
    elif path.suffix.lower() in (".sdf", ".mol"):
        try:
            from rdkit import Chem
            supplier = Chem.SDMolSupplier(ligand_path)
            atoms: list[dict] = []
            for mol in supplier:
                if mol is None:
                    continue
                conf = mol.GetConformer()
                for atom in mol.GetAtoms():
                    pos = conf.GetAtomPosition(atom.GetIdx())
                    atoms.append({
                        "atom_name": atom.GetSymbol(),
                        "element": atom.GetSymbol(),
                        "coords": (pos.x, pos.y, pos.z),
                    })
                break  # Only first molecule
            return atoms
        except Exception as exc:
            logger.warning("Failed to parse SDF %s: %s", ligand_path, exc)
            return []
    else:
        # Try as PDB
        return _parse_pdb_atoms(ligand_path)


def _euclidean_distance(
    coords1: tuple[float, float, float],
    coords2: tuple[float, float, float],
) -> float:
    """Compute Euclidean distance between two 3D points."""
    return (
        (coords1[0] - coords2[0]) ** 2
        + (coords1[1] - coords2[1]) ** 2
        + (coords1[2] - coords2[2]) ** 2
    ) ** 0.5


def _classify_interaction(
    lig_atom: dict,
    prot_atom: dict,
    distance: float,
) -> str:
    """Classify an interaction type based on atom types and distance.

    Parameters
    ----------
    lig_atom : dict
        Ligand atom dict with ``element`` key.
    prot_atom : dict
        Protein atom dict with ``element`` key.
    distance : float
        Distance in Angstroms.

    Returns
    -------
    str
        Interaction type string.
    """
    lig_elem = lig_atom.get("element", "C").upper()
    prot_elem = prot_atom.get("element", "C").upper()

    polar = {"N", "O", "F", "S"}
    hydrophobic = {"C"}

    # H-bond: polar...polar within 3.5 A
    # Convention: HBDonor = ligand provides H (N > O as donor)
    #             HBAcceptor = ligand accepts H from protein
    if distance <= 3.5 and lig_elem in polar and prot_elem in polar:
        # N is primarily a donor (NH, NH2), O/F primarily acceptors (C=O, COO-)
        if lig_elem == "N" and prot_elem in {"O", "F", "S"}:
            return "HBDonor"
        if lig_elem in {"O", "F", "S"} and prot_elem == "N":
            return "HBAcceptor"
        # Same element or ambiguous → use context: N-N=donor, O-O=acceptor
        if lig_elem == "N":
            return "HBDonor"
        return "HBAcceptor"

    # Hydrophobic: C...C within 4.0 A
    if lig_elem in hydrophobic and prot_elem in hydrophobic and distance <= 4.0:
        return "Hydrophobic"

    # Salt bridge: charged N...O
    if distance <= 4.0:
        if (lig_elem == "N" and prot_elem == "O") or (lig_elem == "O" and prot_elem == "N"):
            return "Ionic"

    return "VDW"


# =====================================================================
# STRATEGY 3: Mock interactions
# =====================================================================

def _mock_interactions(uniprot_id: str, smiles: str = "") -> dict:
    """Mock interaction analysis — removed. Requires ProLIF or RDKit."""
    raise RuntimeError(
        f"Mock interaction analysis is not available for {uniprot_id}. "
        "Install ProLIF or RDKit for real interaction analysis."
    )
    """Original docstring:
    Generates deterministic, reproducible interaction data using a hash
    of the SMILES/UniProt ID. Useful for development and testing when
    neither ProLIF nor RDKit is available.

    Parameters
    ----------
    uniprot_id : str
        UniProt accession.
    smiles : str
        Ligand SMILES for deterministic hash.

    Returns
    -------
    dict
        Mock interaction analysis result.
    """
    functional = KNOWN_FUNCTIONAL_RESIDUES.get(uniprot_id, {})
    known_res = functional.get("residues", [718, 745, 793, 855])

    # Generate deterministic mock interactions based on smiles hash
    h = int(hashlib.md5((smiles or uniprot_id).encode()).hexdigest()[:8], 16)

    interaction_types = ["HBDonor", "HBAcceptor", "Hydrophobic", "PiStacking", "CationPi", "VDW"]
    residue_names = ["MET", "THR", "LYS", "LEU", "ASP", "GLU", "ALA", "VAL", "PHE", "TYR"]

    interactions: list[dict] = []
    n_contacts = 3 + (h % 5)  # 3-7 contacts

    for i in range(n_contacts):
        res_idx = (h + i * 7) % len(known_res)
        res_num = known_res[res_idx] if res_idx < len(known_res) else 700 + (h + i) % 200
        res_name = residue_names[(h + i) % len(residue_names)]
        int_type = interaction_types[(h + i * 3) % len(interaction_types)]

        interactions.append({
            "residue": f"{res_name}{res_num}",
            "residue_number": res_num,
            "type": int_type,
            "is_functional": False,  # re-tagged below
        })

    # Mock has no PDB file, so no offset mapping -- use _compute_quality_metrics
    # with pdb_path=None (offset will be 0, matching UniProt numbering directly)
    metrics = _compute_quality_metrics(interactions, uniprot_id, pdb_path=None)

    return {
        "interactions": interactions,
        **metrics,
        "method": "mock",
        "summary": (
            f"Functional contacts: {metrics['functional_contacts']}/{metrics['total_functional']} "
            f"key residues ({metrics['interaction_quality']:.0%}) - "
            f"{metrics['key_hbonds']} key H-bonds"
        ),
    }


# =====================================================================
# CLI / SELF-TEST
# =====================================================================

if __name__ == "__main__":
    import json
    import sys

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    # Test with mock for EGFR + Erlotinib
    erlotinib_smiles = "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OC)cc23)c1"

    print("=" * 72)
    print("DockIt Interaction Analysis Self-Test")
    print("=" * 72)

    result = _mock_interactions("P00533", smiles=erlotinib_smiles)
    print(f"\nMethod: {result['method']}")
    print(f"Summary: {result['summary']}")
    print(f"Interaction quality: {result['interaction_quality']}")
    print(f"Functional contacts: {result['functional_contacts']}/{result['total_functional']}")
    print(f"Key H-bonds: {result['key_hbonds']}")
    print(f"\nInteractions:")
    for inter in result["interactions"]:
        func_tag = " [FUNCTIONAL]" if inter["is_functional"] else ""
        print(f"  {inter['residue']:>10s} -- {inter['type']:<15s}{func_tag}")

    print(f"\nFull JSON:")
    print(json.dumps(result, indent=2))
    print("\n--- All tests passed ---")
    sys.exit(0)
