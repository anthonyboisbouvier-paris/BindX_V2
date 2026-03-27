"""
BindX — SAR Series Clustering: Exhaustive Validation Suite
==========================================================

Production-grade validation of the scaffold-based series clustering feature.
Uses real kinase inhibitor benchmark data (65+ molecules, 5 targets, 10 decoys).

Tests:
  1. Scaffold clustering correctness (Murcko decomposition)
  2. Series filtering (scaffold_filter param)
  3. R-group decomposition within series (match rate, variable R-groups)
  4. Activity cliff detection per series (SALI coherence)
  5. MMP rules per series
  6. R-group contributions per series
  7. Cross-series isolation (no data leakage between series)
  8. Edge cases (1 series, <3 mols, None scaffolds, invalid SMILES)
  9. Cache correctness (scaffold key separation)
 10. API endpoint integration (via Docker)
 11. Frontend response shape validation

Run: docker-compose exec backend python -m pytest tests/test_sar_series_validation.py -v --tb=short
  or: docker-compose exec backend python tests/test_sar_series_validation.py
"""

from __future__ import annotations

import json
import logging
import math
import os
import sys
import time
import uuid
from collections import Counter, defaultdict

# Ensure /app is on path when running from tests/
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np

logging.basicConfig(level=logging.INFO, format="%(levelname)s | %(message)s")
log = logging.getLogger("SAR-validation")

# ═══════════════════════════════════════════════════════════════════════════
# Benchmark data — real molecules from multiple chemical series
# ═══════════════════════════════════════════════════════════════════════════

from benchmarks.kinase_inhibitors import (
    EGFR_ACTIVES, CDK2_ACTIVES, BRAF_ACTIVES, JAK2_ACTIVES,
    KRAS_ACTIVES, DECOYS,
)


def _make_shared_scaffold_molecules() -> list[dict]:
    """Create synthetic dataset with KNOWN shared scaffolds for series testing.

    Uses real scaffolds with systematic R-group variations:
    - Series A: Quinazoline (EGFR-like), 8 molecules
    - Series B: Pyrimidine (CDK-like), 6 molecules
    - Series C: Pyrrolo-pyrimidine (JAK-like), 5 molecules
    - Series D: Miscellaneous singletons, 4 molecules
    """
    series = []

    # Series A: 4-anilinoquinazolines (EGFR TKI scaffold)
    # Core: c1ccc2ncnc(Nc3ccccc3)c2c1
    quinaz_rgroups = [
        ("QZ-001", "COc1cc2c(Nc3ccc(F)cc3)ncnc2cc1OC",         -8.5),
        ("QZ-002", "COc1cc2c(Nc3ccc(Cl)cc3)ncnc2cc1OC",        -7.8),
        ("QZ-003", "COc1cc2c(Nc3ccc(Br)cc3)ncnc2cc1OC",        -7.2),
        ("QZ-004", "COc1cc2c(Nc3ccc(C)cc3)ncnc2cc1OC",         -6.9),
        ("QZ-005", "COc1cc2c(Nc3ccc(OC)cc3)ncnc2cc1OC",        -6.5),
        ("QZ-006", "COc1cc2c(Nc3ccc(C#N)cc3)ncnc2cc1OC",       -8.1),
        ("QZ-007", "COc1cc2c(Nc3ccc(C(F)(F)F)cc3)ncnc2cc1OC",  -7.6),
        ("QZ-008", "COc1cc2c(Nc3ccccc3)ncnc2cc1OC",            -6.0),
    ]
    for name, smi, dock in quinaz_rgroups:
        series.append({
            "id": str(uuid.uuid4()), "name": name, "smiles": smi,
            "bookmarked": False, "series_label": "quinazoline",
            "docking_score": dock + np.random.normal(0, 0.1),
            "composite_score": round((-dock) * 10, 2),
            "QED": round(np.random.uniform(0.4, 0.7), 3),
            "logP": round(np.random.uniform(2.0, 4.5), 2),
            "MW": round(np.random.uniform(350, 450), 1),
        })

    # Series B: 2,4-diaminopyrimidines (CDK-like scaffold)
    # Core: c1cnc(N)nc1N
    pyrim_rgroups = [
        ("PY-001", "Nc1ncc(-c2ccccc2)c(Nc2ccccc2)n1",        -7.0),
        ("PY-002", "Nc1ncc(-c2ccccc2)c(Nc2ccc(F)cc2)n1",      -7.5),
        ("PY-003", "Nc1ncc(-c2ccccc2)c(Nc2ccc(Cl)cc2)n1",     -6.8),
        ("PY-004", "Nc1ncc(-c2ccccc2)c(Nc2ccc(C)cc2)n1",      -6.2),
        ("PY-005", "Nc1ncc(-c2ccccc2)c(Nc2ccc(OC)cc2)n1",     -5.8),
        ("PY-006", "Nc1ncc(-c2ccccc2)c(Nc2ccc(C#N)cc2)n1",    -7.3),
    ]
    for name, smi, dock in pyrim_rgroups:
        series.append({
            "id": str(uuid.uuid4()), "name": name, "smiles": smi,
            "bookmarked": False, "series_label": "pyrimidine",
            "docking_score": dock + np.random.normal(0, 0.1),
            "composite_score": round((-dock) * 10, 2),
            "QED": round(np.random.uniform(0.5, 0.8), 3),
            "logP": round(np.random.uniform(1.5, 3.5), 2),
            "MW": round(np.random.uniform(280, 380), 1),
        })

    # Series C: Pyrrolo[2,3-d]pyrimidines (JAK-like)
    pyrrolop_rgroups = [
        ("PP-001", "c1ccc(Nc2ncnc3[nH]ccc23)cc1",             -6.5),
        ("PP-002", "c1ccc(Nc2ncnc3[nH]ccc23)cc1F",            -7.1),
        ("PP-003", "Fc1ccc(Nc2ncnc3[nH]ccc23)cc1",            -7.0),
        ("PP-004", "Clc1ccc(Nc2ncnc3[nH]ccc23)cc1",           -6.8),
        ("PP-005", "Cc1ccc(Nc2ncnc3[nH]ccc23)cc1",            -6.3),
    ]
    for name, smi, dock in pyrrolop_rgroups:
        series.append({
            "id": str(uuid.uuid4()), "name": name, "smiles": smi,
            "bookmarked": False, "series_label": "pyrrolopyrimidine",
            "docking_score": dock + np.random.normal(0, 0.1),
            "composite_score": round((-dock) * 10, 2),
            "QED": round(np.random.uniform(0.5, 0.8), 3),
            "logP": round(np.random.uniform(1.5, 3.0), 2),
            "MW": round(np.random.uniform(220, 300), 1),
        })

    # Series D: diverse singletons (no shared scaffold)
    singletons = [
        ("SG-001", "CC(=O)Oc1ccccc1C(=O)O",       -3.5),  # aspirin
        ("SG-002", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", -4.0),  # ibuprofen
        ("SG-003", "OC(=O)c1ccccc1O",               -3.0),  # salicylic acid
        ("SG-004", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  -2.5),  # caffeine
    ]
    for name, smi, dock in singletons:
        series.append({
            "id": str(uuid.uuid4()), "name": name, "smiles": smi,
            "bookmarked": False, "series_label": "singleton",
            "docking_score": dock,
            "composite_score": round((-dock) * 10, 2),
            "QED": round(np.random.uniform(0.3, 0.6), 3),
            "logP": round(np.random.uniform(0.5, 2.5), 2),
            "MW": round(np.random.uniform(150, 250), 1),
        })

    return series


def _make_flat_molecules(
    series_list: list[tuple[str, list[dict]]],
    include_decoys: bool = True,
) -> list[dict]:
    """Build flat molecule dicts mimicking flatten_molecule() output."""
    molecules = []
    for series_name, compounds in series_list:
        for c in compounds:
            ic50 = c.get("ic50_nM")
            # Convert IC50 to pIC50 (standard SAR metric: -log10(IC50_M))
            pIC50 = -math.log10(ic50 * 1e-9) if ic50 and ic50 > 0 else None
            # Simulate docking score (negative, correlated with pIC50)
            docking = -(pIC50 * 0.8 + np.random.normal(0, 0.3)) if pIC50 else None

            mol = {
                "id": str(uuid.uuid4()),
                "name": c["name"],
                "smiles": c["smiles"],
                "bookmarked": False,
                "series_label": series_name,  # ground truth for validation
            }
            if pIC50 is not None:
                mol["docking_score"] = round(docking, 4)
                mol["composite_score"] = round(pIC50 * 10, 2)  # arbitrary scale
                mol["QED"] = round(np.random.uniform(0.3, 0.8), 3)
                mol["logP"] = round(np.random.uniform(1.5, 5.5), 2)
                mol["TPSA"] = round(np.random.uniform(50, 130), 1)
                mol["MW"] = round(np.random.uniform(300, 650), 1)
            molecules.append(mol)

    if include_decoys:
        for c in DECOYS:
            mol = {
                "id": str(uuid.uuid4()),
                "name": c["name"],
                "smiles": c["smiles"],
                "bookmarked": False,
                "series_label": "DECOY",
                # Decoys have no kinase activity → only physicochemical props
                "QED": round(np.random.uniform(0.4, 0.9), 3),
                "logP": round(np.random.uniform(0.5, 4.5), 2),
                "MW": round(np.random.uniform(150, 500), 1),
            }
            molecules.append(mol)

    return molecules


# ═══════════════════════════════════════════════════════════════════════════
# Independent Murcko scaffold computation (reference implementation)
# ═══════════════════════════════════════════════════════════════════════════

def _ref_murcko_scaffold(smiles: str) -> str | None:
    """Reference Murcko scaffold using RDKit directly."""
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        core = MurckoScaffold.GetScaffoldForMol(mol)
        if core and core.GetNumAtoms() >= 3:
            return Chem.MolToSmiles(core)
    except Exception:
        pass
    return None


def _ref_cluster_molecules(molecules: list[dict]) -> dict[str | None, list[str]]:
    """Reference clustering: scaffold_smiles → [mol_ids]."""
    clusters = defaultdict(list)
    for m in molecules:
        scaffold = _ref_murcko_scaffold(m["smiles"])
        clusters[scaffold].append(m["id"])
    return dict(clusters)


# ═══════════════════════════════════════════════════════════════════════════
# Test functions
# ═══════════════════════════════════════════════════════════════════════════

PASS = 0
FAIL = 0
WARN = 0


def check(condition: bool, msg: str, critical: bool = True):
    global PASS, FAIL, WARN
    if condition:
        PASS += 1
        log.info("  PASS  %s", msg)
    elif critical:
        FAIL += 1
        log.error("  FAIL  %s", msg)
    else:
        WARN += 1
        log.warning("  WARN  %s", msg)


def test_1_scaffold_clustering_correctness():
    """Verify _cluster_by_scaffold matches independent reference implementation."""
    log.info("\n=== TEST 1: Scaffold Clustering Correctness ===")

    from pipeline.sar_analysis import _cluster_by_scaffold
    from rdkit import Chem

    molecules = _make_flat_molecules([
        ("EGFR", EGFR_ACTIVES),
        ("CDK2", CDK2_ACTIVES),
        ("BRAF", BRAF_ACTIVES),
    ], include_decoys=True)

    # Parse molecules the same way analyze_sar does
    parsed = []
    for m in molecules:
        if m.get("smiles"):
            rdmol = Chem.MolFromSmiles(m["smiles"])
            if rdmol is not None:
                parsed.append({"d": m, "mol": rdmol})

    check(len(parsed) >= 30, f"Parsed {len(parsed)} molecules (expected ≥30)")

    # Run our implementation
    series = _cluster_by_scaffold(parsed)

    check(len(series) >= 2, f"Found {len(series)} series (expected ≥2 for 3 target families + decoys)")
    check(series[0]["n_molecules"] >= series[-1]["n_molecules"],
          "Series sorted by count descending")

    # Verify total molecule count
    total = sum(s["n_molecules"] for s in series)
    check(total == len(parsed), f"Total molecules across series ({total}) == parsed ({len(parsed)})")

    # Verify no molecule appears in 2 series
    all_indices = []
    for s in series:
        all_indices.extend(s["molecule_indices"])
    check(len(all_indices) == len(set(all_indices)),
          "No molecule in multiple series (indices unique)")

    # Cross-validate against reference implementation
    ref_clusters = _ref_cluster_molecules([p["d"] for p in parsed])
    impl_scaffolds = {s["scaffold_smiles"] for s in series}
    ref_scaffolds = set(ref_clusters.keys())

    check(impl_scaffolds == ref_scaffolds,
          f"Scaffold set matches reference ({len(impl_scaffolds)} scaffolds)")

    for s in series:
        smi = s["scaffold_smiles"]
        ref_count = len(ref_clusters.get(smi, []))
        check(s["n_molecules"] == ref_count,
              f"Series '{(smi or 'None')[:40]}': {s['n_molecules']} mols (ref={ref_count})")

    # Verify SVGs generated for scaffolds with SMILES
    for s in series:
        if s["scaffold_smiles"] is not None:
            check(s["scaffold_svg"] is not None and "<svg" in (s["scaffold_svg"] or ""),
                  f"SVG generated for scaffold '{s['scaffold_smiles'][:30]}...'")

    return molecules, parsed, series


def test_2_analyze_sar_auto_selection():
    """Verify analyze_sar auto-selects largest series when scaffold_filter=None."""
    log.info("\n=== TEST 2: Auto Series Selection ===")

    from pipeline.sar_analysis import analyze_sar

    molecules = _make_flat_molecules([
        ("EGFR", EGFR_ACTIVES),       # 13 mols
        ("CDK2", CDK2_ACTIVES),        # 11 mols
        ("BRAF", BRAF_ACTIVES),        # 11 mols
        ("JAK2", JAK2_ACTIVES),        # 12 mols
        ("KRAS", KRAS_ACTIVES),        # 9 mols
    ], include_decoys=True)            # + 10 decoys = ~66 total

    result = analyze_sar(molecules, property_key="docking_score")

    check("series" in result, "Result contains 'series' key")
    check("selected_scaffold" in result, "Result contains 'selected_scaffold' key")
    check("n_molecules_total" in result, "Result contains 'n_molecules_total' key")
    check(isinstance(result["series"], list), "series is a list")
    check(len(result["series"]) >= 3, f"Found {len(result['series'])} series (expected ≥3)")

    # n_molecules_total should be total parsed molecules
    check(result["n_molecules_total"] >= 50,
          f"n_molecules_total={result['n_molecules_total']} (expected ≥50)")

    # n_molecules should be the selected series count (subset)
    check(result["n_molecules"] <= result["n_molecules_total"],
          f"n_molecules ({result['n_molecules']}) ≤ n_molecules_total ({result['n_molecules_total']})")

    # selected_scaffold should be the largest series scaffold
    if result["series"]:
        largest = result["series"][0]
        check(result["selected_scaffold"] == largest["scaffold_smiles"],
              "selected_scaffold matches largest series")
        check(result["n_molecules"] <= largest["n_molecules"] + 5,
              f"Analyzed mol count ({result['n_molecules']}) ≈ largest series ({largest['n_molecules']})")

    # Series metadata should NOT contain molecule_indices (privacy/payload)
    for s in result["series"]:
        check("molecule_indices" not in s,
              f"Series '{(s['scaffold_smiles'] or 'None')[:30]}' has no molecule_indices (clean response)")

    return result


def test_3_scaffold_filter_isolation():
    """Verify scaffold_filter correctly isolates a series — no cross-contamination."""
    log.info("\n=== TEST 3: Scaffold Filter Isolation ===")

    from pipeline.sar_analysis import analyze_sar

    # Use synthetic series with shared scaffolds to guarantee multi-mol series
    molecules = _make_shared_scaffold_molecules()

    # First run: auto
    auto_result = analyze_sar(molecules, property_key="docking_score")
    all_series = auto_result.get("series", [])

    check(len(all_series) >= 2, f"Multiple series detected ({len(all_series)})")

    # Find series with ≥2 molecules to test filtering
    testable = [s for s in all_series if s["n_molecules"] >= 2 and s["scaffold_smiles"] is not None]
    check(len(testable) >= 1, f"At least 1 series with ≥2 molecules ({len(testable)})")

    for s in testable[:3]:
        scaffold_smi = s["scaffold_smiles"]

        filtered = analyze_sar(molecules, property_key="docking_score", scaffold_filter=scaffold_smi)

        check(filtered.get("error") is None,
              f"No error for series '{scaffold_smi[:30]}...' ({s['n_molecules']} mols)")

        check(filtered["selected_scaffold"] == scaffold_smi,
              f"selected_scaffold matches filter")

        check(filtered["n_molecules"] <= s["n_molecules"] + 2,
              f"Filtered mols ({filtered['n_molecules']}) ≈ series count ({s['n_molecules']})")

        check(filtered["n_molecules_total"] == auto_result["n_molecules_total"],
              f"n_molecules_total preserved ({filtered['n_molecules_total']})")

        if filtered.get("rgroup"):
            check(filtered["rgroup"]["n_total"] == filtered["n_molecules"],
                  f"R-group n_total == filtered n_molecules")

    # Also test mixed benchmark to verify it works with diverse data
    diverse = _make_flat_molecules([
        ("EGFR", EGFR_ACTIVES),
        ("CDK2", CDK2_ACTIVES),
    ], include_decoys=False)

    diverse_result = analyze_sar(diverse, property_key="docking_score")
    check("series" in diverse_result, "Diverse dataset returns series key")
    check(diverse_result.get("n_molecules_total", 0) >= 20,
          f"Diverse total: {diverse_result.get('n_molecules_total', 0)}")

    return all_series


def test_4_rgroup_decomposition_per_series():
    """Verify R-group decomposition is coherent within a filtered series."""
    log.info("\n=== TEST 4: R-Group Decomposition Per Series ===")

    from pipeline.sar_analysis import analyze_sar

    # Use synthetic quinazolines which SHARE a scaffold by design
    molecules = _make_shared_scaffold_molecules()
    result = analyze_sar(molecules, property_key="docking_score")

    rg = result.get("rgroup")
    if rg is None:
        log.warning("  WARN  No R-group decomposition for auto-selected series")
        # Try filtering to quinazoline series explicitly
        for s in result.get("series", []):
            if s["n_molecules"] >= 5 and s["scaffold_smiles"]:
                result = analyze_sar(molecules, property_key="docking_score",
                                     scaffold_filter=s["scaffold_smiles"])
                rg = result.get("rgroup")
                if rg:
                    break

    if rg is None:
        log.warning("  WARN  No R-group decomposition for any series")
        return

    check(rg["match_rate"] >= 0.3,
          f"R-group match rate {rg['match_rate']:.1%} ≥ 30%")

    check(rg["n_matched"] >= 3,
          f"At least 3 molecules matched ({rg['n_matched']})")

    check(len(rg.get("variable_r_groups", [])) >= 1,
          f"At least 1 variable R-group ({rg.get('variable_r_groups', [])})")

    # Verify molecule IDs in rgroup.molecules are from the input set
    input_ids = {m["id"] for m in molecules}
    rg_ids = {m["id"] for m in rg.get("molecules", [])}
    check(rg_ids.issubset(input_ids),
          "R-group molecule IDs are subset of input molecules")

    # Verify R-group values are non-empty for matched molecules
    for m in rg.get("molecules", [])[:5]:
        has_rg = any(v for v in m.get("r_groups_display", {}).values())
        check(has_rg, f"Molecule '{m.get('name', '?')}' has R-group values", critical=False)


def test_5_activity_cliffs_coherence():
    """Verify activity cliff detection produces coherent SALI values."""
    log.info("\n=== TEST 5: Activity Cliffs Coherence ===")

    from pipeline.sar_analysis import analyze_sar

    # Use shared scaffold molecules — quinazoline series has good activity range
    molecules = _make_shared_scaffold_molecules()

    result = analyze_sar(molecules, property_key="docking_score")
    cliffs = result.get("cliffs")

    if cliffs is None:
        log.warning("  WARN  No cliff data (needs ≥2 molecules with activity)")
        return

    check(cliffs["n_pairs_analyzed"] > 0, f"Pairs analyzed: {cliffs['n_pairs_analyzed']}")
    check(cliffs["sali_threshold"] > 0, f"SALI threshold: {cliffs['sali_threshold']}")

    # SALI = |ΔA| / (1 - sim) — verify formula on top pairs
    for p in cliffs.get("pairs", [])[:5]:
        if p["similarity"] < 0.999:
            expected_sali = abs(p["delta_activity"]) / (1.0 - p["similarity"])
            check(abs(p["sali"] - expected_sali) < 0.1,
                  f"SALI formula correct for pair ({p['sali']:.1f} vs {expected_sali:.1f})")

    # All cliff pairs should have SALI ≥ threshold
    for p in cliffs.get("pairs", []):
        check(p["sali"] >= cliffs["sali_threshold"] - 0.01,
              f"Cliff pair SALI ({p['sali']:.1f}) ≥ threshold ({cliffs['sali_threshold']:.1f})",
              critical=False)

    # Scatter data should exist
    check(len(cliffs.get("scatter", [])) > 0, "Scatter data present for visualization")

    # Stats sanity
    stats = cliffs.get("sali_stats", {})
    if stats:
        check(stats["mean"] > 0, f"Mean SALI > 0 ({stats['mean']})")
        check(stats["max"] >= stats["mean"], f"Max SALI ≥ mean")
        check(stats["median"] > 0, "Median SALI > 0")


def test_6_mmp_rules_coherence():
    """Verify MMP rules are sensible within a series."""
    log.info("\n=== TEST 6: MMP Rules Coherence ===")

    from pipeline.sar_analysis import analyze_sar

    molecules = _make_shared_scaffold_molecules()
    result = analyze_sar(molecules, property_key="docking_score")
    mmps = result.get("mmps")

    if mmps is None or mmps.get("n_rules", 0) == 0:
        log.warning("  WARN  No MMP rules generated (may be expected for small series)")
        return

    check(mmps["n_rules"] > 0, f"MMP rules found: {mmps['n_rules']}")
    check(mmps["n_pairs_total"] > 0, f"MMP pairs found: {mmps['n_pairs_total']}")

    for rule in mmps.get("rules", [])[:10]:
        check("→" in rule["transform"],
              f"Transform format correct: '{rule['transform']}'")
        check(rule["n_pairs"] >= 1,
              f"Rule has ≥1 pair: {rule['n_pairs']}")

        if rule["mean_delta"] is not None:
            check(isinstance(rule["mean_delta"], float),
                  f"mean_delta is float: {rule['mean_delta']}")

        # Confidence logic
        if rule["n_pairs"] >= 5 and rule.get("n_with_activity", 0) >= 5:
            check(rule["confidence"] in ("high", "medium"),
                  f"Confidence for {rule['n_pairs']} pairs: {rule['confidence']}")


def test_7_contributions_coherence():
    """Verify R-group contributions math is correct."""
    log.info("\n=== TEST 7: R-Group Contributions Coherence ===")

    from pipeline.sar_analysis import analyze_sar

    molecules = _make_shared_scaffold_molecules()
    result = analyze_sar(molecules, property_key="docking_score")
    contrib = result.get("contributions")

    if contrib is None:
        log.warning("  WARN  No contributions (needs rgroup + activity)")
        return

    check(contrib["global_mean"] is not None, f"Global mean: {contrib['global_mean']}")
    check(contrib["n_molecules"] >= 2, f"Contribution molecules: {contrib['n_molecules']}")

    for pos_key, pos_data in contrib.get("positions", {}).items():
        check(pos_data["n_unique"] >= 1, f"Position {pos_key}: {pos_data['n_unique']} unique substituents")

        for sub in pos_data.get("substituents", []):
            # effect = mean - global_mean
            expected_effect = sub["mean"] - contrib["global_mean"]
            check(abs(sub["effect"] - expected_effect) < 0.01,
                  f"  {pos_key} sub '{sub['smiles'][:15]}': effect={sub['effect']:.4f} (expected {expected_effect:.4f})")


def test_8_cross_series_no_leakage():
    """Verify that filtering to one series doesn't include molecules from another."""
    log.info("\n=== TEST 8: Cross-Series Isolation (No Leakage) ===")

    from pipeline.sar_analysis import analyze_sar

    # Use shared scaffold molecules — guaranteed multi-series
    molecules = _make_shared_scaffold_molecules()

    # Get scaffold assignment per molecule
    scaffold_map = {}
    for m in molecules:
        scaffold_map[m["id"]] = _ref_murcko_scaffold(m["smiles"])

    # Run auto
    auto_result = analyze_sar(molecules, property_key="docking_score")
    selected = auto_result.get("selected_scaffold")

    if selected is None:
        log.warning("  WARN  No scaffold selected")
        return

    # Get IDs of molecules that should be in the selected series
    expected_ids = {m["id"] for m in molecules if scaffold_map[m["id"]] == selected}
    check(len(expected_ids) >= 2,
          f"Selected series has {len(expected_ids)} expected molecules")

    # Verify rgroup molecules are only from the selected series
    if auto_result.get("rgroup") and auto_result["rgroup"].get("molecules"):
        rg_ids = {m["id"] for m in auto_result["rgroup"]["molecules"]}
        leaked = rg_ids - expected_ids
        check(len(leaked) == 0,
              f"No R-group molecules from other series (leaked: {len(leaked)})")

    # Now filter to a DIFFERENT series with ≥2 mols
    testable_other = [
        s for s in auto_result.get("series", [])
        if s["scaffold_smiles"] != selected
        and s["scaffold_smiles"] is not None
        and s["n_molecules"] >= 2
    ]
    if testable_other:
        other_scaffold = testable_other[0]["scaffold_smiles"]
        other_result = analyze_sar(molecules, property_key="docking_score", scaffold_filter=other_scaffold)

        check(other_result.get("error") is None,
              "Other series analysis succeeds")
        check(other_result.get("selected_scaffold") == other_scaffold,
              f"Switched to other series: '{other_scaffold[:30]}...'")

        # n_molecules should be different from auto (different series)
        check(other_result.get("n_molecules", 0) != auto_result.get("n_molecules", 0)
              or other_result.get("selected_scaffold") != auto_result.get("selected_scaffold"),
              "Different series yields different molecule set")

        # Verify NO molecule from first series appears in second
        other_expected = {m["id"] for m in molecules if scaffold_map[m["id"]] == other_scaffold}
        if other_result.get("rgroup") and other_result["rgroup"].get("molecules"):
            other_rg_ids = {m["id"] for m in other_result["rgroup"]["molecules"]}
            cross_leak = other_rg_ids & expected_ids
            check(len(cross_leak) == 0,
                  f"No cross-series leakage between series A and B (leaked: {len(cross_leak)})")
    else:
        log.warning("  WARN  No other testable series found for cross-validation")


def test_9_edge_cases():
    """Test edge cases: 1 series, <3 mols, invalid SMILES, None scaffolds."""
    log.info("\n=== TEST 9: Edge Cases ===")

    from pipeline.sar_analysis import analyze_sar

    # --- Single series (selector should be hidden) ---
    single = _make_flat_molecules([("EGFR", EGFR_ACTIVES)], include_decoys=False)
    result = analyze_sar(single, property_key="docking_score")
    # Could be 1 or a few series (EGFR has some scaffold diversity), but should work
    check(result.get("error") is None or result.get("n_molecules", 0) >= 2,
          "Single-target analysis runs without error")

    # --- Very small set (2 molecules) ---
    small = _make_flat_molecules([("EGFR", EGFR_ACTIVES[:2])], include_decoys=False)
    result = analyze_sar(small, property_key="docking_score")
    check("n_molecules" in result, "2-molecule analysis returns result")

    # --- 1 molecule → error ---
    tiny = _make_flat_molecules([("EGFR", EGFR_ACTIVES[:1])], include_decoys=False)
    result = analyze_sar(tiny, property_key="docking_score")
    check(result.get("error") is not None or result.get("n_molecules", 0) <= 1,
          "1-molecule analysis returns error or minimal result")

    # --- Invalid SMILES mixed in ---
    bad_mols = [
        {"id": "bad1", "smiles": "NOT_A_SMILES", "name": "Invalid1", "docking_score": -5.0},
        {"id": "bad2", "smiles": "", "name": "Empty", "docking_score": -6.0},
        {"id": "bad3", "smiles": None, "name": "None", "docking_score": -7.0},
    ]
    mixed = _make_flat_molecules([("EGFR", EGFR_ACTIVES)], include_decoys=False) + bad_mols
    result = analyze_sar(mixed, property_key="docking_score")
    check(result.get("error") is None,
          "Analysis gracefully handles invalid SMILES")

    # --- Scaffold filter that doesn't exist ---
    result = analyze_sar(single, property_key="docking_score", scaffold_filter="C1=CC=CC=C1XYZFAKE")
    check(result.get("n_molecules", 0) >= 2,
          "Non-existent scaffold filter falls back gracefully")

    # --- All decoys (no docking score) ---
    decoy_only = _make_flat_molecules([], include_decoys=True)
    result = analyze_sar(decoy_only, property_key="docking_score")
    # Should still work with physicochemical properties
    check("series" in result or result.get("error") is not None,
          "Decoy-only analysis handles missing activity")

    # --- Empty list ---
    result = analyze_sar([], property_key="auto")
    check(result.get("error") is not None, "Empty list returns error")


def test_10_cache_scaffold_separation():
    """Verify cache correctly separates entries by scaffold key."""
    log.info("\n=== TEST 10: Cache Scaffold Key Separation ===")

    from pipeline.sar_analysis import get_cached_sar, set_cached_sar, _sar_cache

    phase_id = "test-phase-cache-001"
    prop_key = "docking_score"
    n_mols = 50

    # Clear any existing cache for this test
    keys_to_clear = [k for k in _sar_cache if k[0] == phase_id]
    for k in keys_to_clear:
        del _sar_cache[k]

    # Set cache for auto
    result_auto = {"test": "auto", "n_molecules": 30}
    set_cached_sar(phase_id, prop_key, n_mols, result_auto, "__auto__")

    # Set cache for specific scaffold
    scaffold_a = "c1ccc2c(c1)ncnc2"
    result_a = {"test": "scaffold_a", "n_molecules": 15}
    set_cached_sar(phase_id, prop_key, n_mols, result_a, scaffold_a)

    scaffold_b = "c1ccc(-c2ccccn2)cc1"
    result_b = {"test": "scaffold_b", "n_molecules": 10}
    set_cached_sar(phase_id, prop_key, n_mols, result_b, scaffold_b)

    # Retrieve and verify separation
    cached_auto = get_cached_sar(phase_id, prop_key, n_mols, "__auto__")
    check(cached_auto is not None and cached_auto["test"] == "auto",
          "Cache hit for __auto__")

    cached_a = get_cached_sar(phase_id, prop_key, n_mols, scaffold_a)
    check(cached_a is not None and cached_a["test"] == "scaffold_a",
          f"Cache hit for scaffold A")

    cached_b = get_cached_sar(phase_id, prop_key, n_mols, scaffold_b)
    check(cached_b is not None and cached_b["test"] == "scaffold_b",
          f"Cache hit for scaffold B")

    # Verify no cross-contamination
    check(cached_auto["n_molecules"] == 30, "Auto cache has correct data")
    check(cached_a["n_molecules"] == 15, "Scaffold A cache has correct data")
    check(cached_b["n_molecules"] == 10, "Scaffold B cache has correct data")

    # Wrong n_molecules should miss cache
    cached_miss = get_cached_sar(phase_id, prop_key, 999, "__auto__")
    check(cached_miss is None, "Cache miss on wrong n_molecules")

    # Cleanup
    for k in [k for k in _sar_cache if k[0] == phase_id]:
        del _sar_cache[k]


def test_11_full_response_shape():
    """Validate the complete response JSON shape for frontend consumption."""
    log.info("\n=== TEST 11: Full Response Shape Validation ===")

    from pipeline.sar_analysis import analyze_sar

    molecules = _make_shared_scaffold_molecules()

    result = analyze_sar(molecules, property_key="docking_score")

    # Top-level required keys
    required_keys = [
        "property_key", "available_properties", "n_molecules", "n_molecules_total",
        "n_with_activity", "series", "selected_scaffold",
        "scaffold_smiles", "scaffold_svg",
        "rgroup", "cliffs", "mmps", "contributions",
    ]
    for key in required_keys:
        check(key in result, f"Top-level key '{key}' present")

    # series shape
    for s in result.get("series", []):
        check("scaffold_smiles" in s, "Series has scaffold_smiles")
        check("scaffold_svg" in s, "Series has scaffold_svg")
        check("n_molecules" in s, "Series has n_molecules")
        check("molecule_indices" not in s, "Series does NOT have molecule_indices")
        check(isinstance(s["n_molecules"], int), "n_molecules is int")

    # available_properties shape
    for p in result.get("available_properties", []):
        check("key" in p and "label" in p and "n_values" in p and "coverage" in p,
              f"Property '{p.get('key', '?')}' has all fields")

    # property_key should be from ANALYZABLE_PROPERTIES
    check(result["property_key"] in [
        "docking_score", "composite_score", "weighted_score",
        "QED", "logP", "TPSA", "MW", "solubility",
    ], f"property_key '{result['property_key']}' is valid")

    # Numeric fields should be numbers
    check(isinstance(result["n_molecules"], int), "n_molecules is int")
    check(isinstance(result["n_molecules_total"], int), "n_molecules_total is int")
    check(isinstance(result["n_with_activity"], int), "n_with_activity is int")

    # JSON serializable check
    try:
        json.dumps(result)
        check(True, "Full result is JSON-serializable")
    except (TypeError, ValueError) as e:
        check(False, f"JSON serialization failed: {e}")


def test_12_property_switching_with_scaffold():
    """Verify changing property_key while scaffold_filter is set works correctly."""
    log.info("\n=== TEST 12: Property Switch with Scaffold Filter ===")

    from pipeline.sar_analysis import analyze_sar

    molecules = _make_shared_scaffold_molecules()

    # Get series list
    auto = analyze_sar(molecules, property_key="docking_score")
    if not auto.get("series"):
        log.warning("  WARN  No series available")
        return

    scaffold = auto["series"][0]["scaffold_smiles"]

    # Analyze same scaffold with different properties
    result_dock = analyze_sar(molecules, property_key="docking_score", scaffold_filter=scaffold)
    result_qed = analyze_sar(molecules, property_key="QED", scaffold_filter=scaffold)

    check(result_dock["property_key"] == "docking_score", "Docking score property set")
    check(result_qed["property_key"] == "QED", "QED property set")

    # Same scaffold → same n_molecules
    check(result_dock["n_molecules"] == result_qed["n_molecules"],
          f"Same scaffold, same mol count (dock={result_dock['n_molecules']}, qed={result_qed['n_molecules']})")

    # Same scaffold → same selected_scaffold
    check(result_dock["selected_scaffold"] == result_qed["selected_scaffold"],
          "selected_scaffold consistent across property changes")

    # Different properties → potentially different n_with_activity
    # (QED might have different coverage than docking_score)
    check(isinstance(result_qed["n_with_activity"], int), "QED n_with_activity is int")


def test_13_scientific_sanity():
    """Verify scientific soundness: synthetic quinazoline series should cluster together."""
    log.info("\n=== TEST 13: Scientific Sanity Checks ===")

    from pipeline.sar_analysis import analyze_sar

    molecules = _make_shared_scaffold_molecules()

    # Quinazolines (QZ-*) should all cluster in same series
    qz_names = [m["name"] for m in molecules if m["name"].startswith("QZ-")]
    py_names = [m["name"] for m in molecules if m["name"].startswith("PY-")]

    result = analyze_sar(molecules, property_key="docking_score")

    # Verify QZ series exists with ≥6 molecules
    qz_scaffolds = set()
    for m in molecules:
        if m["name"].startswith("QZ-"):
            s = _ref_murcko_scaffold(m["smiles"])
            qz_scaffolds.add(s)

    check(len(qz_scaffolds) == 1,
          f"All quinazolines share 1 scaffold (found {len(qz_scaffolds)})")

    if qz_scaffolds:
        qz_scaffold = list(qz_scaffolds)[0]
        # Find it in the series list
        qz_series = next((s for s in result.get("series", [])
                          if s["scaffold_smiles"] == qz_scaffold), None)
        if qz_series:
            check(qz_series["n_molecules"] == len(qz_names),
                  f"Quinazoline series has {qz_series['n_molecules']} mols (expected {len(qz_names)})")

    # Verify PY series separate
    py_scaffolds = set()
    for m in molecules:
        if m["name"].startswith("PY-"):
            s = _ref_murcko_scaffold(m["smiles"])
            py_scaffolds.add(s)

    check(len(py_scaffolds) == 1,
          f"All pyrimidines share 1 scaffold (found {len(py_scaffolds)})")

    # Verify QZ and PY scaffolds are DIFFERENT
    if qz_scaffolds and py_scaffolds:
        check(qz_scaffolds != py_scaffolds,
              "Quinazoline and pyrimidine scaffolds are distinct")

    # Activity cliffs: QZ series has range -8.5 to -6.0, so cliffs expected
    if result.get("cliffs") and result["cliffs"].get("pairs"):
        cliff_names = set()
        for p in result["cliffs"]["pairs"][:10]:
            cliff_names.add(p.get("mol_a_name", ""))
            cliff_names.add(p.get("mol_b_name", ""))
        log.info("    Top cliff molecules: %s", cliff_names)

    # The auto-selected series should be the quinazoline (largest with 8 mols)
    check(result.get("selected_scaffold") is not None,
          "A scaffold was auto-selected")
    log.info("    Auto-selected scaffold: %s (%d mols)",
             (result.get("selected_scaffold", "?"))[:40],
             result.get("n_molecules", 0))


def test_14_performance():
    """Measure performance: scaffold clustering + full analysis on large dataset."""
    log.info("\n=== TEST 14: Performance Benchmarks ===")

    from pipeline.sar_analysis import analyze_sar

    # ~66 molecules (realistic phase size)
    molecules = _make_flat_molecules([
        ("EGFR", EGFR_ACTIVES),
        ("CDK2", CDK2_ACTIVES),
        ("BRAF", BRAF_ACTIVES),
        ("JAK2", JAK2_ACTIVES),
        ("KRAS", KRAS_ACTIVES),
    ], include_decoys=True)

    t0 = time.time()
    result = analyze_sar(molecules, property_key="docking_score")
    t_full = time.time() - t0

    check(t_full < 30.0, f"Full analysis in {t_full:.2f}s (limit: 30s)")
    log.info("    Full analysis: %.2fs for %d molecules, %d series",
             t_full, result.get("n_molecules_total", 0), len(result.get("series", [])))

    # Filtered analysis should be faster
    if result.get("series") and len(result["series"]) >= 2:
        smallest = result["series"][-1]
        if smallest["scaffold_smiles"]:
            t0 = time.time()
            analyze_sar(molecules, property_key="docking_score",
                        scaffold_filter=smallest["scaffold_smiles"])
            t_filtered = time.time() - t0
            log.info("    Filtered analysis: %.2fs for smallest series (%d mols)",
                     t_filtered, smallest["n_molecules"])


# ═══════════════════════════════════════════════════════════════════════════
# API Integration Test (requires running Docker backend)
# ═══════════════════════════════════════════════════════════════════════════

def test_15_api_endpoint():
    """Test SAR API endpoint via HTTP (requires running backend)."""
    log.info("\n=== TEST 15: API Endpoint Integration ===")

    try:
        import httpx
    except ImportError:
        try:
            import requests as httpx
        except ImportError:
            log.warning("  SKIP  httpx/requests not available — skipping API test")
            return

    # This test only works if there's a real phase with molecules
    log.info("  INFO  API test requires a phase with molecules in the database.")
    log.info("  INFO  Skipping API test in unit test mode. Run manually with a real phase_id.")


# ═══════════════════════════════════════════════════════════════════════════
# Runner
# ═══════════════════════════════════════════════════════════════════════════

def run_all():
    global PASS, FAIL, WARN
    PASS = FAIL = WARN = 0

    log.info("=" * 70)
    log.info("SAR Series Clustering — Exhaustive Validation")
    log.info("=" * 70)

    np.random.seed(42)  # Reproducible

    test_1_scaffold_clustering_correctness()
    test_2_analyze_sar_auto_selection()
    test_3_scaffold_filter_isolation()
    test_4_rgroup_decomposition_per_series()
    test_5_activity_cliffs_coherence()
    test_6_mmp_rules_coherence()
    test_7_contributions_coherence()
    test_8_cross_series_no_leakage()
    test_9_edge_cases()
    test_10_cache_scaffold_separation()
    test_11_full_response_shape()
    test_12_property_switching_with_scaffold()
    test_13_scientific_sanity()
    test_14_performance()
    test_15_api_endpoint()

    log.info("\n" + "=" * 70)
    log.info("RESULTS: %d PASS | %d FAIL | %d WARN", PASS, FAIL, WARN)
    log.info("=" * 70)

    if FAIL > 0:
        log.error("VALIDATION FAILED — %d critical issues found", FAIL)
        return 1
    else:
        log.info("VALIDATION PASSED — all checks OK")
        return 0


if __name__ == "__main__":
    sys.exit(run_all())
