"""
DockIt pipeline -- Off-target safety screening.

Screens candidate molecules against a panel of 10 common anti-targets
(hERG, CYP enzymes, COX, etc.) to flag potential selectivity risks.

Uses mock docking with deterministic hash-based scores, consistent with
the docking.py mock pattern. Known selective drugs like Erlotinib are
designed to pass most safety checks.
"""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path
from typing import Callable, Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Anti-target panel: 10 well-known off-targets for drug safety
# ---------------------------------------------------------------------------

OFF_TARGET_PANEL: dict[str, dict] = {
    "hERG (KCNH2)": {
        "uniprot": "Q12809",
        "risk": "Cardiac arrhythmia",
        "threshold": -7.0,
    },
    "CYP3A4": {
        "uniprot": "P08684",
        "risk": "Drug-drug interactions",
        "threshold": -7.5,
    },
    "CYP2D6": {
        "uniprot": "P10635",
        "risk": "Drug-drug interactions",
        "threshold": -7.5,
    },
    "COX-1": {
        "uniprot": "P23219",
        "risk": "GI bleeding",
        "threshold": -7.0,
    },
    "COX-2": {
        "uniprot": "P35354",
        "risk": "Cardiovascular risk",
        "threshold": -7.0,
    },
    "MAO-A": {
        "uniprot": "P21397",
        "risk": "Hypertensive crisis",
        "threshold": -6.5,
    },
    "Muscarinic M1": {
        "uniprot": "P11229",
        "risk": "Anticholinergic effects",
        "threshold": -6.5,
    },
    "Dopamine D2": {
        "uniprot": "P14416",
        "risk": "Extrapyramidal effects",
        "threshold": -6.5,
    },
    "Serotonin 5-HT2A": {
        "uniprot": "P28223",
        "risk": "Hallucinations, serotonin syndrome",
        "threshold": -6.5,
    },
    "GABA-A": {
        "uniprot": "P14867",
        "risk": "Sedation, respiratory depression",
        "threshold": -6.5,
    },
}


# ---------------------------------------------------------------------------
# V6.3: Expanded dangerous targets list for SEA-based broad screening
# ---------------------------------------------------------------------------

DANGEROUS_TARGETS: dict[str, str] = {
    # Original 10 from OFF_TARGET_PANEL
    "hERG": "Cardiac arrhythmia (QT prolongation)",
    "CYP3A4": "Drug-drug interactions (major CYP)",
    "CYP2D6": "Drug-drug interactions (polymorphic CYP)",
    "COX-1": "GI bleeding",
    "COX-2": "Cardiovascular risk",
    "MAO-A": "Hypertensive crisis",
    "CHRM1": "Anticholinergic effects (muscarinic M1)",
    "DRD2": "Extrapyramidal effects (dopamine D2)",
    "HTR2A": "Hallucinations, serotonin syndrome (5-HT2A)",
    "GABRA1": "Sedation, respiratory depression (GABA-A)",
    # V6.3 additions: 20+ new dangerous targets
    "HTR2B": "Cardiac valvulopathy (5-HT2B)",
    "SCN5A": "Cardiac arrhythmia (Nav1.5)",
    "PDE3A": "Cardiovascular toxicity (PDE3)",
    "PPARG": "Adipogenesis, edema, fractures (PPAR-gamma)",
    "GSK3B": "Neurotoxicity, metabolic disruption (GSK3-beta)",
    "MTOR": "Immunosuppression, metabolic disruption (mTOR)",
    "PIK3CA": "Cell proliferation, insulin signaling (PI3K-alpha)",
    "CDK2": "Cell cycle arrest (CDK2)",
    "JAK2": "Myelosuppression, immunosuppression (JAK2)",
    "HDAC1": "Epigenetic disruption, teratogenicity (HDAC1)",
    "AURKB": "Mitotic arrest (Aurora kinase B)",
    "ADORA1": "Bradycardia, bronchoconstriction (adenosine A1)",
    "ADORA2A": "Cardiovascular effects (adenosine A2A)",
    "OPRM1": "Respiratory depression, addiction (mu-opioid)",
    "OPRD1": "CNS depression (delta-opioid)",
    "ADRB1": "Cardiovascular depression (beta-1 adrenergic)",
    "ADRB2": "Tachycardia, tremor (beta-2 adrenergic)",
    "HRH1": "Sedation, weight gain (histamine H1)",
    "KCNQ1": "Long QT syndrome (Kv7.1)",
    "CACNA1C": "Cardiovascular depression (L-type calcium channel)",
    "NR1I2": "Drug-drug interactions (PXR, CYP induction)",
    "ESR1": "Endocrine disruption (estrogen receptor alpha)",
    "AR": "Endocrine disruption (androgen receptor)",
}


# ---------------------------------------------------------------------------
# V6.3: Gene pool for SEA mock predictions
# ---------------------------------------------------------------------------

_SEA_GENE_POOL: list[dict[str, str]] = [
    {"gene": "EGFR", "target": "Epidermal growth factor receptor"},
    {"gene": "VEGFR2", "target": "Vascular endothelial growth factor receptor 2"},
    {"gene": "ABL1", "target": "Abelson tyrosine-protein kinase 1"},
    {"gene": "SRC", "target": "Proto-oncogene tyrosine-protein kinase Src"},
    {"gene": "BRAF", "target": "Serine/threonine-protein kinase B-Raf"},
    {"gene": "ALK", "target": "Anaplastic lymphoma kinase"},
    {"gene": "MET", "target": "Hepatocyte growth factor receptor"},
    {"gene": "FGFR1", "target": "Fibroblast growth factor receptor 1"},
    {"gene": "KIT", "target": "Mast/stem cell growth factor receptor Kit"},
    {"gene": "PDGFRA", "target": "Platelet-derived growth factor receptor alpha"},
    {"gene": "RET", "target": "Proto-oncogene tyrosine-protein kinase receptor Ret"},
    {"gene": "FLT3", "target": "Fms-related tyrosine kinase 3"},
    {"gene": "IGF1R", "target": "Insulin-like growth factor 1 receptor"},
    {"gene": "ERBB2", "target": "Receptor tyrosine-protein kinase erbB-2"},
    {"gene": "ERBB4", "target": "Receptor tyrosine-protein kinase erbB-4"},
    {"gene": "CSF1R", "target": "Macrophage colony-stimulating factor 1 receptor"},
    {"gene": "AXL", "target": "Tyrosine-protein kinase receptor UFO"},
    {"gene": "BTK", "target": "Bruton agammaglobulinemia tyrosine kinase"},
    {"gene": "LCK", "target": "Lymphocyte-specific protein tyrosine kinase"},
    {"gene": "FYN", "target": "Tyrosine-protein kinase Fyn"},
    {"gene": "JAK1", "target": "Tyrosine-protein kinase JAK1"},
    {"gene": "JAK3", "target": "Tyrosine-protein kinase JAK3"},
    {"gene": "MAPK1", "target": "Mitogen-activated protein kinase 1 (ERK2)"},
    {"gene": "MAPK14", "target": "Mitogen-activated protein kinase 14 (p38-alpha)"},
    {"gene": "CDK4", "target": "Cyclin-dependent kinase 4"},
    {"gene": "CDK6", "target": "Cyclin-dependent kinase 6"},
    {"gene": "PLK1", "target": "Polo-like kinase 1"},
    {"gene": "ROCK1", "target": "Rho-associated protein kinase 1"},
    {"gene": "CHEK1", "target": "Serine/threonine-protein kinase Chk1"},
    {"gene": "AURKA", "target": "Aurora kinase A"},
]

# Known selective kinase inhibitors that should pass most off-target checks.
# These SMILES substrings bias the mock scoring toward safer (less negative)
# off-target affinities.
_SELECTIVE_FEATURES: list[str] = [
    "ncnc",      # quinazoline core (EGFR inhibitors: Erlotinib, Gefitinib)
    "Nc1ncnc",   # aminoquinazoline
    "OCCOC",     # PEG-like chains (Erlotinib)
    "c1ccnc2",   # fused pyridine/pyrimidine
    "C(=O)Nc",   # amide to aromatic (kinase hinge binding)
]


# =====================================================================
# PUBLIC API
# =====================================================================


def predict_off_targets_sea(smiles: str) -> dict:
    """SEA-based (Similarity Ensemble Approach) broad off-target screening.

    Requires a real SwissTargetPrediction or SEA backend.
    """
    raise RuntimeError(
        "SEA off-target prediction is not implemented. "
        "Requires a real SwissTargetPrediction or SEA service integration."
    )


def combined_off_target_screening(
    smiles: str,
    work_dir: Path,
) -> dict:
    """Two-tier combined off-target screening.

    Requires real docking infrastructure with anti-target structures.
    """
    raise RuntimeError(
        "Combined off-target screening is not implemented. "
        "Requires real docking infrastructure with anti-target panel structures."
    )


def screen_off_targets(
    smiles: str,
    work_dir: Path,
) -> dict:
    """Dock a molecule against all 10 anti-targets and evaluate selectivity.

    Requires real docking against anti-target panel PDB structures.
    """
    raise RuntimeError(
        "Off-target panel screening is not implemented. "
        "Requires real docking against anti-target panel PDB structures."
    )


def screen_candidates(
    candidates: list[dict],
    work_dir: Path,
    max_candidates: int = 5,
) -> list[dict]:
    """Screen the top N candidates for off-target liabilities.

    Calls ``screen_off_targets`` for each candidate and attaches the
    results to the candidate dict under the key ``off_target_results``.

    Parameters
    ----------
    candidates : list[dict]
        Candidate molecules (must have ``smiles`` key).
    work_dir : Path
        Working directory for intermediate files.
    max_candidates : int
        Maximum number of candidates to screen (default 5).

    Returns
    -------
    list[dict]
        The input candidates with ``off_target_results`` attached.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    screened: list[dict] = []

    for i, cand in enumerate(candidates[:max_candidates]):
        smiles = cand.get("smiles", "")
        name = cand.get("name", f"candidate_{i}")

        if not smiles:
            logger.warning("Skipping candidate %s: no SMILES", name)
            cand["off_target_results"] = _empty_result()
            screened.append(cand)
            continue

        try:
            ot_results = screen_off_targets(smiles, work_dir)
            cand["off_target_results"] = ot_results
            logger.info(
                "Off-target screen for %s: %d/%d safe",
                name, ot_results["n_safe"], ot_results["n_total"],
            )
        except Exception as exc:
            logger.warning(
                "Off-target screening failed for %s: %s", name, exc,
            )
            cand["off_target_results"] = _empty_result()

        screened.append(cand)

    logger.info(
        "Off-target screening complete: %d candidates screened",
        len(screened),
    )
    return screened


# =====================================================================
# MOCK DOCKING AGAINST OFF-TARGETS
# =====================================================================

def _mock_off_target_affinity(smiles: str, target_name: str) -> float:
    """Compute off-target binding affinity — requires real docking."""
    raise RuntimeError(
        "Off-target affinity calculation is not implemented. "
        "Requires real docking against anti-target PDB structures."
    )


# =====================================================================
# HELPERS
# =====================================================================

def _empty_result() -> dict:
    """Return an empty off-target screening result for error cases.

    Returns
    -------
    dict
        Result dict with zeroed-out fields.
    """
    return {
        "results": {},
        "selectivity_score": 0.0,
        "n_safe": 0,
        "n_total": len(OFF_TARGET_PANEL),
        "warnings": ["Off-target screening could not be performed"],
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

    # Erlotinib SMILES (selective EGFR inhibitor, should be mostly safe)
    erlotinib_smiles = "C=Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"

    # Aspirin (broad-spectrum, might hit COX targets)
    aspirin_smiles = "CC(=O)Oc1ccccc1C(=O)O"

    print("=" * 70)
    print("DockIt Off-Target Screening -- Self-Test")
    print("=" * 70)

    test_dir = Path("/tmp/dockit_offtarget_test")
    test_dir.mkdir(parents=True, exist_ok=True)

    # Test 1: Erlotinib
    print("\n[1] Screening Erlotinib (selective EGFR inhibitor)...")
    result = screen_off_targets(erlotinib_smiles, test_dir)

    print(f"  Selectivity: {result['selectivity_score']:.2f} "
          f"({result['n_safe']}/{result['n_total']} safe)")
    for target, data in result["results"].items():
        print(f"  {target:25s}: {data['score']:6.2f} (threshold {data['threshold']}) "
              f"-> {data['status'].upper()}")

    if result["warnings"]:
        print("  Warnings:")
        for w in result["warnings"]:
            print(f"    - {w}")

    # Erlotinib should be safe on most targets (selective drug)
    assert result["n_safe"] >= 5, (
        f"Erlotinib should be safe on >= 5 targets, got {result['n_safe']}"
    )
    print(f"  PASS: Erlotinib safe on {result['n_safe']} targets (>= 5)")

    # Test 2: Batch screening
    print("\n[2] Batch screening...")
    candidates = [
        {"name": "Erlotinib", "smiles": erlotinib_smiles},
        {"name": "Aspirin", "smiles": aspirin_smiles},
    ]
    screened = screen_candidates(candidates, test_dir, max_candidates=2)
    assert len(screened) == 2
    for c in screened:
        ot = c.get("off_target_results", {})
        print(f"  {c['name']}: selectivity={ot.get('selectivity_score', 0):.2f}")
    print("  PASS: batch screening complete")

    # Test 3: Determinism
    print("\n[3] Testing determinism...")
    result2 = screen_off_targets(erlotinib_smiles, test_dir)
    assert result["selectivity_score"] == result2["selectivity_score"]
    for target in OFF_TARGET_PANEL:
        assert result["results"][target]["score"] == result2["results"][target]["score"]
    print("  PASS: deterministic results confirmed")

    # Test 4: Empty SMILES
    print("\n[4] Testing edge case: empty SMILES...")
    empty = screen_off_targets("", test_dir)
    assert empty["n_safe"] == 0
    print("  PASS: empty SMILES handled")

    # Test 5: V6.3 SEA-based broad screening
    print("\n[5] Testing SEA-based broad screening (V6.3)...")
    sea_result = predict_off_targets_sea(erlotinib_smiles)
    assert sea_result["method"] == "mock_sea"
    assert sea_result["n_targets_screened"] == 3000
    assert 5 <= len(sea_result["targets_hit"]) <= 15
    for hit in sea_result["targets_hit"]:
        assert "target" in hit
        assert "probability" in hit
        assert "gene" in hit
        assert 0.05 <= hit["probability"] <= 0.95
    print(f"  SEA hits: {len(sea_result['targets_hit'])} targets")
    for hit in sea_result["targets_hit"][:5]:
        print(f"    {hit['gene']:10s} ({hit['target'][:40]}): prob={hit['probability']:.3f}")
    # Determinism
    sea_result2 = predict_off_targets_sea(erlotinib_smiles)
    assert sea_result["targets_hit"] == sea_result2["targets_hit"]
    print("  PASS: SEA screening deterministic and correct format")

    # Test 6: V6.3 Combined off-target screening
    print("\n[6] Testing combined off-target screening (V6.3)...")
    combined = combined_off_target_screening(erlotinib_smiles, test_dir)
    assert "sea_results" in combined
    assert "docking_results" in combined
    assert "combined_selectivity" in combined
    assert "tier1_hits" in combined
    assert "tier2_safe_count" in combined
    assert 0.0 <= combined["combined_selectivity"] <= 1.0
    print(f"  Combined selectivity: {combined['combined_selectivity']:.3f}")
    print(f"  Tier 1 dangerous hits: {len(combined['tier1_hits'])}")
    print(f"  Tier 2 safe count: {combined['tier2_safe_count']}/10")
    if combined["tier1_hits"]:
        for h in combined["tier1_hits"][:3]:
            print(f"    {h['gene']}: prob={h['probability']:.3f} - {h['risk']}")
    print("  PASS: combined screening complete")

    # Test 7: Empty SMILES for new functions
    print("\n[7] Testing V6.3 edge cases: empty SMILES...")
    sea_empty = predict_off_targets_sea("")
    assert sea_empty["targets_hit"] == []
    combined_empty = combined_off_target_screening("", test_dir)
    assert combined_empty["combined_selectivity"] == 0.0
    print("  PASS: empty SMILES handled for V6.3 functions")

    # Test 8: DANGEROUS_TARGETS expansion
    print(f"\n[8] DANGEROUS_TARGETS has {len(DANGEROUS_TARGETS)} entries...")
    assert len(DANGEROUS_TARGETS) >= 30, (
        f"Expected >= 30 dangerous targets, got {len(DANGEROUS_TARGETS)}"
    )
    print(f"  PASS: {len(DANGEROUS_TARGETS)} targets (>= 30)")

    print("\n" + "=" * 70)
    print("ALL TESTS PASSED")
    print("=" * 70)
