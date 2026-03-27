"""
BindX — SAR (Structure-Activity Relationship) analysis engine.

Scientifically rigorous implementation based on established methods:
- Murcko scaffold extraction (Bemis-Murcko, J. Med. Chem. 1996)
- R-group decomposition (RDKit RGroupDecompose)
- Activity cliff detection via SALI (Guha & Van Drie, JCIM 2008)
- Matched Molecular Pairs (Hussain-Rea algorithm via rdMMPA)
- R-group contribution analysis (mean ± std per substituent)

No Free-Wilson regression (requires scikit-learn). Instead, we report
per-substituent mean activity — simpler, more interpretable, no overfitting risk.
"""

from __future__ import annotations

import logging
import re
from collections import defaultdict

import numpy as np

logger = logging.getLogger(__name__)

# In-memory cache: (phase_id, property_key, scaffold_key) -> (n_molecules, result)
_sar_cache: dict[tuple[str, str, str], tuple[int, dict]] = {}

# Properties eligible for SAR analysis (priority order)
ANALYZABLE_PROPERTIES = [
    "docking_score", "composite_score", "weighted_score",
    "QED", "logP", "TPSA", "MW", "solubility",
]

# Internal pipeline keys that should NOT appear in Activity Explorer chips.
# These come from the generation pipeline and get flattened by _deep_flatten()
# alongside legitimate dashboard properties. Blacklist approach = any new
# legitimate property passes automatically without needing a whitelist update.
_EXCLUDED_PROPERTIES = {
    # generation pipeline internals
    "estimated_affinity", "novelty_score", "lipinski_violations",
    "diversity_score", "drug_likeness_score", "generation_score",
    "parent_similarity", "validity_score",
    # ranking / internal metadata
    "_rank", "_score", "rank",
    # string/boolean fields that sometimes appear as numeric
    "safety_color_code", "lipinski_pass", "ro3_pass",
    "inchikey", "reagents_available", "is_representative",
    # synth cost is a string ($, $$, $$$)
    "synth_cost_estimate",
    # pharmacophore feature_counts sub-keys (leaked by _deep_flatten)
    # these are internal feature counts, not independent molecular properties
    "acceptor", "aromatic", "donor", "hydrophobic", "negative", "positive",
}


# ═══════════════════════════════════════════════════════════════════════════
# Public API
# ═══════════════════════════════════════════════════════════════════════════

def analyze_sar(
    molecules: list[dict],
    property_key: str = "auto",
    scaffold_filter: str | None = None,
    matrix_r1: str | None = None,
    matrix_r2: str | None = None,
    r_filters: dict[str, str | None] | None = None,
) -> dict:
    """Run SAR analysis on a list of flat molecule dicts.

    Parameters
    ----------
    molecules : list of dicts with at least {id, smiles} and numeric properties.
    property_key : which property to analyze. "auto" picks the best available.
    scaffold_filter : if set, scope analysis to molecules matching this Murcko scaffold SMILES.
                      If None, auto-selects the largest series.
    matrix_r1 : R-group key for matrix rows (overrides auto-selection).
    matrix_r2 : R-group key for matrix columns (overrides auto-selection).
    r_filters : dict of {R-group: value} to filter non-displayed R-groups.

    Returns a JSON-serializable dict with all SAR analysis results.
    """
    if not molecules:
        return _error_result("No molecules provided")

    valid = [m for m in molecules if m.get("smiles")]
    if len(valid) < 2:
        return _error_result("Need at least 2 molecules with SMILES")

    try:
        from rdkit import Chem
    except ImportError:
        return _error_result("RDKit not available")

    # Parse molecules — ensure each has an 'id' for internal tracking
    parsed = []
    for i, m in enumerate(valid):
        rdmol = Chem.MolFromSmiles(m["smiles"])
        if rdmol is not None:
            if "id" not in m:
                m = {**m, "id": m.get("name", f"mol_{i}")}
            parsed.append({"d": m, "mol": rdmol})

    if len(parsed) < 2:
        return _error_result("Less than 2 valid molecules after SMILES parsing")

    n_molecules_total = len(parsed)

    # ── Cluster by scaffold (always, O(n)) ──
    series_list = _cluster_by_scaffold(parsed)

    # Determine which scaffold to use
    selected_scaffold = scaffold_filter
    if selected_scaffold is None and series_list:
        # Auto-select largest series, but only if it has ≥2 molecules
        if series_list[0]["n_molecules"] >= 2:
            selected_scaffold = series_list[0]["scaffold_smiles"]
        else:
            # All scaffolds are singletons → analyze all molecules (no filtering)
            selected_scaffold = None
            logger.info("No scaffold with ≥2 molecules, analyzing all %d molecules", len(parsed))

    # Filter parsed to selected series
    if selected_scaffold is not None:
        series_entry = next(
            (s for s in series_list if s["scaffold_smiles"] == selected_scaffold), None
        )
        if series_entry:
            idx_set = set(series_entry["molecule_indices"])
            parsed = [parsed[i] for i in range(len(parsed)) if i in idx_set]
        else:
            logger.warning("Scaffold filter %s not found, using all molecules", selected_scaffold)
            selected_scaffold = None

    if len(parsed) < 2:
        # Return error but still include series metadata so frontend can offer alternatives
        err = _error_result("Less than 2 molecules in selected series")
        err["series"] = [
            {"scaffold_smiles": s["scaffold_smiles"], "scaffold_svg": s["scaffold_svg"], "n_molecules": s["n_molecules"]}
            for s in series_list
        ]
        err["selected_scaffold"] = selected_scaffold
        err["n_molecules_total"] = n_molecules_total
        err["n_molecules"] = len(parsed)
        return err

    # Auto-detect property key
    if property_key == "auto":
        property_key = _detect_best_property(parsed)

    # List all available numeric properties
    available_props = _list_available_properties(parsed)

    # Extract activities for the selected property
    activities = {}  # mol_id -> float
    for md in parsed:
        val = md["d"].get(property_key)
        if val is not None and isinstance(val, (int, float)):
            activities[md["d"]["id"]] = float(val)

    has_activity = len(activities) >= 2

    # ── Step 1: Murcko Scaffold ──
    scaffold_smiles, scaffold_svg = _extract_murcko_scaffold(parsed)

    # ── Step 2: R-group decomposition ──
    rgroup_result = _rgroup_decomposition(
        parsed, scaffold_smiles,
        matrix_r1=matrix_r1, matrix_r2=matrix_r2, r_filters=r_filters,
    )

    # ── Step 3: Activity cliffs (SALI) ──
    cliff_result = _compute_activity_cliffs(parsed, activities, property_key) if has_activity else None

    # ── Step 4: MMP analysis ──
    try:
        mmp_result = _compute_mmps(parsed, activities, property_key)
    except Exception as e:
        logger.warning("MMP analysis failed: %s", e)
        mmp_result = {"rules": [], "n_pairs_total": 0, "n_rules": 0, "property_key": property_key}

    # ── Step 5: R-group contributions ──
    contribution_result = None
    if rgroup_result and has_activity:
        contribution_result = _compute_rgroup_contributions(
            rgroup_result["molecules"], activities, property_key,
        )

    # ── Series metadata (strip molecule_indices for JSON response) ──
    series_meta = [
        {"scaffold_smiles": s["scaffold_smiles"], "scaffold_svg": s["scaffold_svg"], "n_molecules": s["n_molecules"]}
        for s in series_list
    ]

    # ── Build molecule_properties lookup (all numeric props per mol) ──
    molecule_properties = {}
    for md in parsed:
        mid = md["d"].get("id", "")
        if mid:
            molecule_properties[mid] = {
                k: v for k, v in md["d"].items()
                if isinstance(v, (int, float)) and k not in ("id",)
            }

    # ── Assemble result ──
    return {
        "property_key": property_key,
        "available_properties": available_props,
        "n_molecules": len(parsed),
        "n_molecules_total": n_molecules_total,
        "n_with_activity": len(activities),

        # Series
        "series": series_meta,
        "selected_scaffold": selected_scaffold,

        # Scaffold
        "scaffold_smiles": scaffold_smiles,
        "scaffold_svg": scaffold_svg,

        # R-groups
        "rgroup": rgroup_result,

        # Activity cliffs
        "cliffs": cliff_result,

        # MMPs
        "mmps": mmp_result,

        # R-group contributions
        "contributions": contribution_result,

        # All numeric properties per molecule (for frontend multi-prop display)
        "molecule_properties": molecule_properties,
    }


# ═══════════════════════════════════════════════════════════════════════════
# Scaffold Clustering (series detection)
# ═══════════════════════════════════════════════════════════════════════════

def _cluster_by_scaffold(parsed: list[dict], min_series_size: int = 3) -> list[dict]:
    """Group molecules by their *generic* Murcko scaffold → one series per scaffold.

    Uses MakeScaffoldGeneric (ring topology only, no side-chains) to reduce
    fragmentation — typically yields ~10-20 series instead of ~95 for 100 mols.
    Series with fewer than `min_series_size` molecules are merged into "Other".

    Returns a list sorted by n_molecules desc:
    [{"scaffold_smiles": str|None, "scaffold_svg": str|None, "n_molecules": int, "molecule_indices": [int]}]
    """
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    scaffold_groups: dict[str | None, list[int]] = defaultdict(list)
    scaffold_mols: dict[str, object] = {}

    for i, md in enumerate(parsed):
        try:
            core = MurckoScaffold.GetScaffoldForMol(md["mol"])
            if core and core.GetNumAtoms() >= 3:
                # Make generic: strip side chains, keep only ring topology
                generic = MurckoScaffold.MakeScaffoldGeneric(core)
                smi = Chem.MolToSmiles(generic)
                scaffold_groups[smi].append(i)
                if smi not in scaffold_mols:
                    scaffold_mols[smi] = generic
            else:
                scaffold_groups[None].append(i)
        except Exception:
            scaffold_groups[None].append(i)

    # Merge small series (< min_series_size) into "Other"
    other_indices = list(scaffold_groups.pop(None, []))
    small_keys = [k for k, v in scaffold_groups.items() if len(v) < min_series_size]
    for k in small_keys:
        other_indices.extend(scaffold_groups.pop(k))

    result = []
    for smi, indices in scaffold_groups.items():
        svg = _mol_to_svg(scaffold_mols[smi], 200, 140) if smi in scaffold_mols else None
        result.append({
            "scaffold_smiles": smi,
            "scaffold_svg": svg,
            "n_molecules": len(indices),
            "molecule_indices": indices,
        })

    if other_indices:
        result.append({
            "scaffold_smiles": None,
            "scaffold_svg": None,
            "n_molecules": len(other_indices),
            "molecule_indices": other_indices,
        })

    # Sort by count desc, "Other" (None) last
    result.sort(key=lambda s: (s["scaffold_smiles"] is not None, s["n_molecules"]), reverse=True)

    logger.info("Scaffold clustering: %d series from %d molecules (generic scaffolds, min_size=%d)",
                len(result), len(parsed), min_series_size)
    return result


# ═══════════════════════════════════════════════════════════════════════════
# Step 1 — Murcko Scaffold
# ═══════════════════════════════════════════════════════════════════════════

def _extract_murcko_scaffold(parsed: list[dict]) -> tuple[str | None, str | None]:
    """Extract the most common Murcko framework from the dataset."""
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    scaffold_counts: dict[str, int] = defaultdict(int)
    scaffold_mols: dict[str, object] = {}

    for md in parsed:
        try:
            core = MurckoScaffold.GetScaffoldForMol(md["mol"])
            if core and core.GetNumAtoms() >= 3:
                smi = Chem.MolToSmiles(core)
                scaffold_counts[smi] += 1
                scaffold_mols[smi] = core
        except Exception:
            continue

    if not scaffold_counts:
        return None, None

    # Most common scaffold
    best_smi = max(scaffold_counts, key=scaffold_counts.get)
    best_mol = scaffold_mols[best_smi]

    svg = _mol_to_svg(best_mol, 300, 200)
    logger.info("Murcko scaffold: %s (covers %d/%d molecules)",
                best_smi, scaffold_counts[best_smi], len(parsed))

    return best_smi, svg


def _mol_to_svg(mol, width=300, height=200) -> str | None:
    try:
        from rdkit.Chem import AllChem
        from rdkit.Chem.Draw import rdMolDraw2D
        AllChem.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.drawOptions().addAtomIndices = False
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception as e:
        logger.warning("SVG generation failed: %s", e)
        return None


# ═══════════════════════════════════════════════════════════════════════════
# Step 2 — R-Group Decomposition
# ═══════════════════════════════════════════════════════════════════════════

def _rgroup_decomposition(
    parsed: list[dict],
    scaffold_smiles: str | None,
    matrix_r1: str | None = None,
    matrix_r2: str | None = None,
    r_filters: dict[str, str | None] | None = None,
) -> dict | None:
    """Decompose molecules against the Murcko scaffold."""
    if not scaffold_smiles:
        return None

    try:
        from rdkit import Chem
        from rdkit.Chem import rdRGroupDecomposition

        core_mol = Chem.MolFromSmiles(scaffold_smiles)
        if core_mol is None:
            return None

        rdmols = [md["mol"] for md in parsed]
        res, unmatched = rdRGroupDecomposition.RGroupDecompose(
            [core_mol], rdmols, asSmiles=True,
        )

        if not res:
            return None

        n_matched = len(res)
        match_rate = n_matched / len(parsed)

        # Collect R-group keys
        rg_keys = sorted(
            k for k in set().union(*(entry.keys() for entry in res))
            if k.startswith("R")
        )

        # Build R-group summary
        rg_summary = {}
        for rk in rg_keys:
            counts = defaultdict(int)
            for entry in res:
                val = entry.get(rk, "")
                if val:
                    counts[val] += 1

            sorted_vals = sorted(counts.items(), key=lambda x: -x[1])
            n_unique = len(sorted_vals)

            # Skip R-groups with only 1 unique value (invariant)
            rg_summary[rk] = {
                "n_unique": n_unique,
                "is_variable": n_unique > 1,
                "values": [
                    {"smiles": s, "display": _clean_rg(s), "count": c}
                    for s, c in sorted_vals[:30]
                ],
            }

        # Build molecule-level R-group assignments
        mol_rgroups = []
        for i, entry in enumerate(res):
            if i >= len(parsed):
                break
            md = parsed[i]["d"]
            rgs = {}
            rgs_display = {}
            for rk in rg_keys:
                raw = entry.get(rk, "")
                rgs[rk] = raw
                rgs_display[rk] = _clean_rg(raw)
            # Collect all numeric properties for this molecule
            mol_props = {
                k: v for k, v in md.items()
                if isinstance(v, (int, float)) and k not in ("id",)
            }
            mol_rgroups.append({
                "id": md.get("id", ""),
                "name": md.get("name", ""),
                "smiles": md.get("smiles", ""),
                "r_groups": rgs,
                "r_groups_display": rgs_display,
                "properties": mol_props,
            })

        variable_rgs = [k for k, v in rg_summary.items() if v["is_variable"]]

        # Build SAR matrix if ≥2 variable R-groups
        sar_matrix = None
        if len(variable_rgs) >= 2:
            # Use custom axes if provided, otherwise default to first two variable R-groups
            r1 = matrix_r1 if matrix_r1 and matrix_r1 in rg_summary else variable_rgs[0]
            r2 = matrix_r2 if matrix_r2 and matrix_r2 in rg_summary else variable_rgs[1]
            # Ensure r1 != r2
            if r1 == r2 and len(variable_rgs) >= 2:
                r2 = next((v for v in variable_rgs if v != r1), variable_rgs[1])
            sar_matrix = _build_2d_matrix(r1, r2, mol_rgroups, r_filters=r_filters)

        return {
            "n_matched": n_matched,
            "n_total": len(parsed),
            "match_rate": round(match_rate, 3),
            "r_groups": rg_summary,
            "variable_r_groups": variable_rgs,
            "molecules": mol_rgroups,
            "sar_matrix": sar_matrix,
        }

    except Exception as e:
        logger.warning("R-group decomposition failed: %s", e)
        return None


def _build_2d_matrix(
    rg1: str,
    rg2: str,
    mol_rgroups: list[dict],
    r_filters: dict[str, str | None] | None = None,
) -> dict:
    """Build 2D SAR matrix for any 2 selected R-groups.

    Parameters
    ----------
    rg1 : R-group key for rows (e.g. "R1")
    rg2 : R-group key for columns (e.g. "R2")
    mol_rgroups : list of molecule dicts with r_groups_display and properties
    r_filters : optional dict of {R-group_key: value} to filter molecules
                by non-displayed R-groups (e.g. {"R3": "F"}). None or "All" = no filter.

    Returns dict with rows, cols, cells (each cell can hold multiple molecules).
    """
    # Filter molecules by r_filters (for R-groups not on axes)
    filtered = mol_rgroups
    if r_filters:
        for rk, rv in r_filters.items():
            if rv is None or rv == "" or rv.lower() == "all":
                continue
            filtered = [
                m for m in filtered
                if m.get("r_groups_display", {}).get(rk, "") == rv
            ]

    r1_vals = sorted(set(
        m["r_groups_display"].get(rg1, "") for m in filtered
        if m["r_groups_display"].get(rg1)
    ))
    r2_vals = sorted(set(
        m["r_groups_display"].get(rg2, "") for m in filtered
        if m["r_groups_display"].get(rg2)
    ))

    # Group molecules by (r1, r2) — multiple molecules possible per cell
    lookup: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for m in filtered:
        r1 = m["r_groups_display"].get(rg1, "")
        r2 = m["r_groups_display"].get(rg2, "")
        if r1 and r2:
            lookup[(r1, r2)].append(m)

    cells = []
    for r1 in r1_vals:
        row = []
        for r2 in r2_vals:
            mols = lookup.get((r1, r2), [])
            if not mols:
                row.append(None)
            elif len(mols) == 1:
                m = mols[0]
                row.append({
                    "mol_id": m["id"],
                    "smiles": m["smiles"],
                    "name": m.get("name", ""),
                    "properties": m.get("properties", {}),
                    "n_molecules": 1,
                })
            else:
                # Multiple molecules in same cell — return all
                row.append({
                    "mol_id": mols[0]["id"],  # representative
                    "smiles": mols[0]["smiles"],
                    "name": mols[0].get("name", ""),
                    "properties": mols[0].get("properties", {}),
                    "n_molecules": len(mols),
                    "molecules": [
                        {
                            "mol_id": m["id"],
                            "smiles": m["smiles"],
                            "name": m.get("name", ""),
                            "properties": m.get("properties", {}),
                        }
                        for m in mols
                    ],
                })
        cells.append(row)

    return {
        "r1_key": rg1,
        "r2_key": rg2,
        "rows": r1_vals,
        "cols": r2_vals,
        "cells": cells,
    }


def _clean_rg(smi: str) -> str:
    """Strip RDKit attachment-point notation for display."""
    if not smi:
        return ""
    cleaned = re.sub(r'\[\*:\d+\]', '', smi)
    cleaned = cleaned.strip('.')
    if not cleaned or cleaned == '[H]':
        return "H"
    return cleaned


# ═══════════════════════════════════════════════════════════════════════════
# Step 3 — Activity Cliffs (SALI)
# ═══════════════════════════════════════════════════════════════════════════

def _compute_activity_cliffs(
    parsed: list[dict],
    activities: dict[str, float],
    property_key: str,
) -> dict:
    """Compute activity cliffs using SALI (Structure-Activity Landscape Index).

    SALI(i,j) = |A_i - A_j| / (1 - sim(i,j))

    Reference: Guha & Van Drie, JCIM 2008.
    Fingerprint: ECFP4 (Morgan radius 2, 2048 bits) — standard for SALI.
    """
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem

    # Filter to molecules with activity
    active_mols = []
    for md in parsed:
        mid = md["d"].get("id", "")
        if mid in activities:
            active_mols.append(md)

    n = len(active_mols)
    if n < 2:
        return {"pairs": [], "n_cliffs": 0, "sali_threshold": 0}

    # Cap at 200 molecules for O(n²) computation
    if n > 200:
        # Keep molecules with most extreme activities
        active_mols.sort(key=lambda m: activities.get(m["d"]["id"], 0))
        # Take top 100 and bottom 100
        active_mols = active_mols[:100] + active_mols[-100:]
        n = len(active_mols)

    # Compute ECFP4 fingerprints
    fps = []
    mids = []
    for md in active_mols:
        try:
            fp = AllChem.GetMorganFingerprintAsBitVect(md["mol"], 2, nBits=2048)
            fps.append(fp)
            mids.append(md["d"]["id"])
        except Exception:
            continue

    n = len(fps)
    if n < 2:
        return {"pairs": [], "n_cliffs": 0, "sali_threshold": 0}

    # Compute all pairwise SALI
    sali_pairs = []
    all_sali = []

    for i in range(n):
        for j in range(i + 1, n):
            sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            act_i = activities[mids[i]]
            act_j = activities[mids[j]]
            delta_act = abs(act_i - act_j)

            if sim >= 0.999:
                # Nearly identical molecules — not a cliff
                continue

            sali = delta_act / (1.0 - sim)
            all_sali.append(sali)

            sali_pairs.append({
                "mol_a_id": mids[i],
                "mol_b_id": mids[j],
                "similarity": round(sim, 3),
                "delta_activity": round(delta_act, 4),
                "sali": round(sali, 2),
                # For scatter plot
                "activity_a": round(act_i, 4),
                "activity_b": round(act_j, 4),
            })

    if not all_sali:
        return {"pairs": [], "n_cliffs": 0, "sali_threshold": 0}

    # Threshold: 95th percentile (standard)
    threshold = float(np.percentile(all_sali, 95))

    # Sort by SALI desc, keep top 50 cliff pairs
    sali_pairs.sort(key=lambda p: p["sali"], reverse=True)
    cliff_pairs = [p for p in sali_pairs if p["sali"] >= threshold]

    # Enrich with molecule names
    name_map = {md["d"]["id"]: md["d"].get("name", "") for md in parsed}
    smiles_map = {md["d"]["id"]: md["d"].get("smiles", "") for md in parsed}
    for p in cliff_pairs[:50]:
        p["mol_a_name"] = name_map.get(p["mol_a_id"], "")
        p["mol_b_name"] = name_map.get(p["mol_b_id"], "")
        p["mol_a_smiles"] = smiles_map.get(p["mol_a_id"], "")
        p["mol_b_smiles"] = smiles_map.get(p["mol_b_id"], "")

    # Scatter data (subsample for performance)
    scatter = sali_pairs[:2000]
    for p in scatter:
        p["is_cliff"] = p["sali"] >= threshold

    # Global landscape stats
    arr = np.array(all_sali)

    return {
        "pairs": cliff_pairs[:50],
        "n_cliffs": len(cliff_pairs),
        "n_pairs_analyzed": len(all_sali),
        "sali_threshold": round(threshold, 2),
        "sali_stats": {
            "mean": round(float(arr.mean()), 2),
            "median": round(float(np.median(arr)), 2),
            "std": round(float(arr.std()), 2),
            "max": round(float(arr.max()), 2),
        },
        "scatter": scatter[:500],
        "property_key": property_key,
    }


# ═══════════════════════════════════════════════════════════════════════════
# Step 4 — Matched Molecular Pairs
# ═══════════════════════════════════════════════════════════════════════════

def _compute_mmps(
    parsed: list[dict],
    activities: dict[str, float],
    property_key: str,
    max_rules: int = 50,
) -> dict:
    """Compute MMPs using rdMMPA (Hussain-Rea algorithm).

    Fragments molecules, groups by shared context (core), computes
    per-transformation statistics: count, mean Δ, std Δ.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMMPA
    except ImportError:
        return {"rules": [], "n_pairs": 0}

    # Fragment each molecule
    context_index: dict[str, list[tuple[int, str]]] = defaultdict(list)

    for i, md in enumerate(parsed):
        try:
            frags = rdMMPA.FragmentMol(md["mol"], maxCuts=1, resultsAsMols=False)
        except Exception:
            continue

        if not frags:
            continue

        for _core, chains in frags:
            if not chains:
                continue
            try:
                chain_mol = Chem.MolFromSmiles(chains)
                if chain_mol is None:
                    continue
                parts = Chem.MolToSmiles(chain_mol).split('.')
                if len(parts) < 2:
                    continue
                parts.sort(key=len, reverse=True)
                context = parts[0]
                variable = parts[1]
            except Exception:
                continue

            context_index[context].append((i, variable))

    # Build transformation rules
    transform_data: dict[str, list[dict]] = defaultdict(list)
    seen_pairs = set()

    for context, mol_vars in context_index.items():
        if len(mol_vars) < 2:
            continue
        for ii in range(len(mol_vars)):
            for jj in range(ii + 1, len(mol_vars)):
                a_idx, var_a = mol_vars[ii]
                b_idx, var_b = mol_vars[jj]
                if var_a == var_b:
                    continue

                pair_key = (min(a_idx, b_idx), max(a_idx, b_idx))
                if pair_key in seen_pairs:
                    continue
                seen_pairs.add(pair_key)

                ma = parsed[a_idx]["d"]
                mb = parsed[b_idx]["d"]
                act_a = activities.get(ma.get("id"))
                act_b = activities.get(mb.get("id"))

                # Canonical direction: alphabetical order of variable SMILES
                if var_a > var_b:
                    var_a, var_b = var_b, var_a
                    act_a, act_b = act_b, act_a

                transform_key = f"{_clean_rg(var_a)} → {_clean_rg(var_b)}"
                delta = (act_b - act_a) if (act_a is not None and act_b is not None) else None

                transform_data[transform_key].append({
                    "delta": delta,
                    "mol_a": ma.get("name", ma.get("id", "")),
                    "mol_b": mb.get("name", mb.get("id", "")),
                    "from_smiles": var_a,  # raw fragment with [*:1]
                    "to_smiles": var_b,    # raw fragment with [*:1]
                })

    # Aggregate into rules with statistics
    rules = []
    total_pairs = 0

    for transform, pairs in transform_data.items():
        n = len(pairs)
        total_pairs += n
        deltas = [p["delta"] for p in pairs if p["delta"] is not None]

        # Extract raw fragment SMILES from first pair (same for all pairs of this rule)
        first_pair = pairs[0]
        rule = {
            "transform": transform,
            "from_smiles": first_pair.get("from_smiles"),
            "to_smiles": first_pair.get("to_smiles"),
            "n_pairs": n,
            "n_with_activity": len(deltas),
        }

        if deltas:
            arr = np.array(deltas)
            rule["mean_delta"] = round(float(arr.mean()), 4)
            rule["std_delta"] = round(float(arr.std()), 4) if len(deltas) > 1 else None
            rule["median_delta"] = round(float(np.median(arr)), 4)
            rule["min_delta"] = round(float(arr.min()), 4)
            rule["max_delta"] = round(float(arr.max()), 4)

            # Confidence: "high" if ≥5 pairs and consistent sign
            if n >= 5 and len(deltas) >= 5:
                sign_consistency = abs(arr.mean()) / (arr.std() + 1e-9)
                rule["confidence"] = "high" if sign_consistency > 1.0 else "medium"
            elif n >= 3:
                rule["confidence"] = "medium"
            else:
                rule["confidence"] = "low"
        else:
            rule["mean_delta"] = None
            rule["std_delta"] = None
            rule["confidence"] = "low"

        rules.append(rule)

    # Sort by n_pairs desc, then by |mean_delta|
    rules.sort(key=lambda r: (
        r["n_pairs"],
        abs(r["mean_delta"]) if r["mean_delta"] is not None else 0,
    ), reverse=True)

    return {
        "rules": rules[:max_rules],
        "n_pairs_total": total_pairs,
        "n_rules": len(rules),
        "property_key": property_key,
    }


# ═══════════════════════════════════════════════════════════════════════════
# Step 5 — R-Group Contributions
# ═══════════════════════════════════════════════════════════════════════════

def _compute_rgroup_contributions(
    mol_rgroups: list[dict],
    activities: dict[str, float],
    property_key: str,
) -> dict:
    """Compute mean activity per substituent at each R-group position.

    Returns per-position, per-substituent: mean, std, n, vs. global mean.
    """
    global_acts = list(activities.values())
    if not global_acts:
        return None

    global_mean = float(np.mean(global_acts))
    global_std = float(np.std(global_acts)) if len(global_acts) > 1 else 0.0

    positions = {}

    # Collect all R-group keys across molecules
    all_rg_keys = set()
    for m in mol_rgroups:
        all_rg_keys.update(m.get("r_groups_display", {}).keys())

    for rk in sorted(all_rg_keys):
        # Group activities by substituent
        sub_activities: dict[str, list[float]] = defaultdict(list)
        for m in mol_rgroups:
            mid = m.get("id", "")
            rg_val = m.get("r_groups_display", {}).get(rk, "")
            if rg_val and mid in activities:
                sub_activities[rg_val].append(activities[mid])

        if not sub_activities:
            continue

        substituents = []
        for sub_smi, acts in sub_activities.items():
            arr = np.array(acts)
            n = len(acts)
            mean = float(arr.mean())
            std = float(arr.std()) if n > 1 else None

            # Effect vs global mean
            effect = mean - global_mean

            substituents.append({
                "smiles": sub_smi,
                "n": n,
                "mean": round(mean, 4),
                "std": round(std, 4) if std is not None else None,
                "effect": round(effect, 4),
            })

        # Sort by effect (best first — depends on property direction)
        substituents.sort(key=lambda s: s["effect"], reverse=True)

        positions[rk] = {
            "substituents": substituents,
            "n_unique": len(substituents),
        }

    return {
        "positions": positions,
        "global_mean": round(global_mean, 4),
        "global_std": round(global_std, 4),
        "property_key": property_key,
        "n_molecules": len(global_acts),
    }


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

def _detect_best_property(parsed: list[dict]) -> str:
    """Auto-detect the best property for SAR analysis."""
    for key in ANALYZABLE_PROPERTIES:
        count = sum(
            1 for md in parsed
            if md["d"].get(key) is not None and isinstance(md["d"].get(key), (int, float))
        )
        if count >= len(parsed) * 0.3:
            return key
    return "docking_score"  # default


def _list_available_properties(parsed: list[dict]) -> list[dict]:
    """List ALL numeric properties with their coverage.

    Returns every key that has at least 1 numeric value.
    `auto_eligible` = True when coverage ≥ 30% (recommended for analysis).
    """
    # Collect all numeric keys across molecules
    key_counts: dict[str, int] = defaultdict(int)
    total = len(parsed)
    for md in parsed:
        for k, v in md["d"].items():
            if k in ("id", "name", "smiles", "bookmarked"):
                continue
            if isinstance(v, (int, float)) and k not in _EXCLUDED_PROPERTIES:
                key_counts[k] += 1

    # Build result, prioritizing ANALYZABLE_PROPERTIES order first
    seen = set()
    result = []
    for key in ANALYZABLE_PROPERTIES:
        if key in key_counts:
            count = key_counts[key]
            coverage = round(count / total, 2)
            result.append({
                "key": key,
                "label": _prop_label(key),
                "n_values": count,
                "coverage": coverage,
                "auto_eligible": coverage >= 0.3,
            })
            seen.add(key)

    # Then add remaining numeric keys alphabetically
    for key in sorted(key_counts.keys()):
        if key in seen:
            continue
        count = key_counts[key]
        coverage = round(count / total, 2)
        result.append({
            "key": key,
            "label": _prop_label(key),
            "n_values": count,
            "coverage": coverage,
            "auto_eligible": coverage >= 0.3,
        })

    return result


def _prop_label(key: str) -> str:
    labels = {
        "docking_score": "Docking Score",
        "composite_score": "Composite Score",
        "weighted_score": "Weighted Score",
        "QED": "QED",
        "logP": "LogP",
        "TPSA": "TPSA",
        "MW": "Mol. Weight",
        "HBD": "H-Bond Donors",
        "HBA": "H-Bond Acceptors",
        "rotatable_bonds": "Rotatable Bonds",
        "heavy_atom_count": "Heavy Atoms",
        "sa_score": "SA Score",
        "ligand_efficiency": "Ligand Efficiency",
        "cnn_score": "CNN Score",
        "cnn_affinity": "CNN Affinity",
        "cnn_vs": "CNN VS",
        "consensus_ecr": "Consensus ECR",
        "pocket_distance": "Pocket Distance",
        "solubility": "Solubility",
        "oral_bioavailability": "Oral Bioavail.",
        "intestinal_permeability": "Intestinal Perm.",
        "BBB": "BBB Permeability",
        "vd": "Volume Dist.",
        "plasma_protein_binding": "Plasma Binding",
        "half_life": "Half-Life",
        "hERG": "hERG",
        "hepatotoxicity": "Hepatotoxicity",
        "carcinogenicity": "Carcinogenicity",
        "ames_mutagenicity": "Ames Mutagenicity",
        "cns_mpo": "CNS MPO",
        "confidence_score": "Confidence",
        "tanimoto_to_centroid": "Tanimoto Centroid",
        "selectivity_score": "Selectivity",
        "off_target_hits": "Off-Target Hits",
        "selectivity_ratio": "Selectivity Ratio",
        "n_interactions": "Interactions",
        "functional_contacts": "Contacts",
        "interaction_quality": "Interaction Quality",
        "key_hbonds": "Key H-Bonds",
        "n_synth_steps": "Synth. Steps",
        "synth_confidence": "Synth. Confidence",
        "synth_cost_estimate": "Synth. Cost",
        "pharmacophore_fit_1": "Ref. 1 Fit",
        "pharmacophore_fit_2": "Ref. 2 Fit",
        "pharmacophore_fit_3": "Ref. 3 Fit",
        "pharmacophore_fit_4": "Ref. 4 Fit",
        "pharmacophore_fit_5": "Ref. 5 Fit",
        "pharmacophore_fit_best": "Best Ref. Fit",
        "pharmacophore_fit_avg": "Avg Ref. Fit",
        "cyp1a2_inhibitor": "CYP1A2",
        "cyp2c9_inhibitor": "CYP2C9",
        "cyp2c19_inhibitor": "CYP2C19",
        "cyp2d6_inhibitor": "CYP2D6",
        "cyp3a4_inhibitor": "CYP3A4",
        "pgp_substrate": "P-gp Substrate",
        "skin_sensitization": "Skin Sensitization",
        "cluster_id": "Cluster ID",
        "cluster_size": "Cluster Size",
        "pchembl_value": "pChEMBL",
        "activity_value_nM": "Activity (nM)",
    }
    return labels.get(key, key)


def _error_result(msg: str) -> dict:
    return {
        "error": msg,
        "fallback": True,
        "property_key": None,
        "available_properties": [],
        "n_molecules": 0,
        "n_molecules_total": 0,
        "n_with_activity": 0,
        "series": [],
        "selected_scaffold": None,
        "scaffold_smiles": None,
        "scaffold_svg": None,
        "rgroup": None,
        "cliffs": None,
        "mmps": None,
        "contributions": None,
    }


# ═══════════════════════════════════════════════════════════════════════════
# MMP-Guided Suggestion Generation
# ═══════════════════════════════════════════════════════════════════════════

def apply_mmp_transform(
    molecules: list[dict],
    from_smiles: str,
    to_smiles: str,
    max_results: int = 50,
) -> list[dict]:
    """Apply an MMP transformation to a list of molecules.

    Builds a reaction SMARTS from the from/to fragments, applies it to each
    molecule, validates products, and returns up to max_results candidates.

    Parameters
    ----------
    molecules : list of dicts with at least {id, smiles, name}.
    from_smiles : raw fragment SMILES with attachment point [*:1].
    to_smiles : raw fragment SMILES with attachment point [*:1].
    max_results : cap on returned suggestions.

    Returns list of suggestion dicts.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
    except ImportError:
        logger.error("RDKit not available for MMP transform")
        return []

    if not from_smiles or not to_smiles:
        return []

    # Normalize attachment point notation for reaction SMARTS
    def _to_rxn_smarts(frag: str) -> str:
        """Convert rdMMPA fragment to reaction SMARTS atom."""
        # Ensure [*:1] notation
        if "[*:1]" not in frag:
            frag = frag.replace("[*]", "[*:1]")
        return frag

    from_rxn = _to_rxn_smarts(from_smiles)
    to_rxn = _to_rxn_smarts(to_smiles)

    rxn_smarts = f"[*:1]{from_rxn.replace('[*:1]', '')}>>[*:1]{to_rxn.replace('[*:1]', '')}"

    try:
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
    except Exception:
        # Fallback: try substructure replacement approach
        logger.warning("Reaction SMARTS failed for %s, trying substructure replacement", rxn_smarts)
        return _apply_substructure_replacement(molecules, from_smiles, to_smiles, max_results)

    if rxn is None:
        return _apply_substructure_replacement(molecules, from_smiles, to_smiles, max_results)

    seen_smiles = set()
    suggestions = []

    for mol_dict in molecules:
        if len(suggestions) >= max_results:
            break

        smiles = mol_dict.get("smiles", "")
        rdmol = Chem.MolFromSmiles(smiles)
        if rdmol is None:
            continue

        try:
            products = rxn.RunReactants((rdmol,))
        except Exception:
            continue

        for product_set in products:
            if len(suggestions) >= max_results:
                break
            for product in product_set:
                try:
                    Chem.SanitizeMol(product)
                    product_smiles = Chem.MolToSmiles(product)
                except Exception:
                    continue

                if not product_smiles or product_smiles in seen_smiles:
                    continue
                # Skip if product is same as input
                if product_smiles == Chem.MolToSmiles(rdmol):
                    continue

                seen_smiles.add(product_smiles)

                # Compute basic properties
                try:
                    props = _compute_basic_props(product)
                except Exception:
                    props = {}

                suggestions.append({
                    "smiles": product_smiles,
                    "source_mol_id": str(mol_dict.get("id", "")),
                    "source_mol_name": mol_dict.get("name", ""),
                    "qed": props.get("qed"),
                    "mw": props.get("mw"),
                    "logp": props.get("logp"),
                    "hbd": props.get("hbd"),
                    "hba": props.get("hba"),
                    "valid": props.get("valid", True),
                })

    logger.info("MMP transform %s: %d suggestions from %d molecules", rxn_smarts, len(suggestions), len(molecules))
    return suggestions


def _apply_substructure_replacement(
    molecules: list[dict],
    from_smiles: str,
    to_smiles: str,
    max_results: int = 50,
) -> list[dict]:
    """Fallback: simple substructure search+replace when reaction SMARTS fails."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    # Clean fragments (remove attachment points for substructure match)
    clean_from = _clean_rg(from_smiles)
    clean_to = _clean_rg(to_smiles)

    from_mol = Chem.MolFromSmiles(clean_from)
    if from_mol is None:
        return []

    seen_smiles = set()
    suggestions = []

    for mol_dict in molecules:
        if len(suggestions) >= max_results:
            break

        smiles = mol_dict.get("smiles", "")
        rdmol = Chem.MolFromSmiles(smiles)
        if rdmol is None:
            continue

        if not rdmol.HasSubstructMatch(from_mol):
            continue

        try:
            to_mol = Chem.MolFromSmiles(clean_to)
            if to_mol is None:
                continue
            product = AllChem.ReplaceSubstructs(rdmol, from_mol, to_mol, replaceAll=False)[0]
            Chem.SanitizeMol(product)
            product_smiles = Chem.MolToSmiles(product)
        except Exception:
            continue

        if not product_smiles or product_smiles in seen_smiles:
            continue
        if product_smiles == Chem.MolToSmiles(rdmol):
            continue

        seen_smiles.add(product_smiles)

        try:
            props = _compute_basic_props(product)
        except Exception:
            props = {}

        suggestions.append({
            "smiles": product_smiles,
            "source_mol_id": str(mol_dict.get("id", "")),
            "source_mol_name": mol_dict.get("name", ""),
            "qed": props.get("qed"),
            "mw": props.get("mw"),
            "logp": props.get("logp"),
            "hbd": props.get("hbd"),
            "hba": props.get("hba"),
            "valid": props.get("valid", True),
        })

    return suggestions


def _compute_basic_props(rdmol) -> dict:
    """Compute QED, MW, logP, HBD, HBA for a molecule."""
    from rdkit.Chem import Descriptors, rdMolDescriptors, QED as QEDModule

    mw = Descriptors.MolWt(rdmol)
    logp = Descriptors.MolLogP(rdmol)
    hbd = rdMolDescriptors.CalcNumHBD(rdmol)
    hba = rdMolDescriptors.CalcNumHBA(rdmol)
    try:
        qed = QEDModule.qed(rdmol)
    except Exception:
        qed = None

    return {
        "qed": round(qed, 3) if qed is not None else None,
        "mw": round(mw, 1),
        "logp": round(logp, 2),
        "hbd": hbd,
        "hba": hba,
        "valid": True,
    }


def get_cached_sar(phase_id: str, prop_key: str, n_molecules: int, scaffold_key: str = "__auto__") -> dict | None:
    cached = _sar_cache.get((phase_id, prop_key, scaffold_key))
    if cached and cached[0] == n_molecules:
        return cached[1]
    return None


def set_cached_sar(phase_id: str, prop_key: str, n_molecules: int, result: dict, scaffold_key: str = "__auto__"):
    _sar_cache[(phase_id, prop_key, scaffold_key)] = (n_molecules, result)
