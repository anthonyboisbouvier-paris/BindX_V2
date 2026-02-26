"""
DockIt pipeline -- Confidence scoring for docking results.

Computes a multi-component confidence score for each molecule result,
reflecting the reliability of the entire pipeline: structure quality,
pocket detection, docking method, ADMET availability, and synthesis
feasibility.

Each component is scored on a [0, 1] scale. The overall confidence
is the weighted mean of all components.

V6: Updated structure/pocket confidence hierarchy.
    Added disorder penalty (fraction_disordered > 0.3 => -0.10).
"""

from __future__ import annotations

import logging
from typing import Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Component weights (must sum to 1.0)
# ---------------------------------------------------------------------------

_WEIGHTS: dict[str, float] = {
    "structure": 0.25,
    "pocket": 0.15,
    "docking": 0.25,
    "admet": 0.20,
    "synthesis": 0.15,
}


# =====================================================================
# PUBLIC API
# =====================================================================

def calculate_confidence(
    molecule: dict,
    structure_source: str = "alphafold",
    disorder_info: Optional[dict] = None,
) -> dict:
    """Compute a multi-component confidence score for a single molecule.

    The confidence reflects how much trust we can place in the pipeline
    results for this molecule, based on the quality and availability of
    each pipeline component.

    Parameters
    ----------
    molecule : dict
        Molecule result dict. Expected keys (all optional):
        - ``source``: str ("chembl", "zinc", "mock_generation", etc.)
        - ``admet``: dict or None
        - ``synthesis_route``: dict or None
        - ``docking_method``: str or None
    structure_source : str
        Source of the protein structure: "alphafold", "esmfold",
        "pdb_experimental", "mock", or "unknown".
    disorder_info : dict, optional
        Disorder prediction result from ``predict_disorder()``.
        If provided and ``fraction_disordered > 0.3``, a penalty
        of 0.10 is applied to the overall confidence.

    Returns
    -------
    dict
        Keys:
        - ``overall``: float in [0.0, 1.0]
        - ``components``: dict mapping component names to
          ``{"score": float, "note": str}``
        - ``disorder_penalty``: float (0.0 or 0.10)
    """
    components: dict[str, dict] = {}

    # --- Structure confidence ---
    struct_score, struct_note = _score_structure(structure_source)
    components["structure"] = {"score": struct_score, "note": struct_note}

    # --- Pocket confidence ---
    pocket_score, pocket_note = _score_pocket(molecule)
    components["pocket"] = {"score": pocket_score, "note": pocket_note}

    # --- Docking confidence ---
    docking_score, docking_note = _score_docking(molecule)
    components["docking"] = {"score": docking_score, "note": docking_note}

    # --- ADMET confidence ---
    admet_score, admet_note = _score_admet(molecule)
    components["admet"] = {"score": admet_score, "note": admet_note}

    # --- Synthesis confidence ---
    synth_score, synth_note = _score_synthesis(molecule)
    components["synthesis"] = {"score": synth_score, "note": synth_note}

    # --- Weighted overall ---
    overall = 0.0
    for comp_name, weight in _WEIGHTS.items():
        comp_data = components.get(comp_name, {"score": 0.0})
        overall += weight * comp_data["score"]

    # --- V6: Disorder penalty ---
    disorder_penalty = 0.0
    if disorder_info and isinstance(disorder_info, dict):
        fraction_disordered = disorder_info.get("fraction_disordered", 0.0)
        if fraction_disordered > 0.3:
            disorder_penalty = 0.10
            logger.info(
                "Disorder penalty applied: fraction_disordered=%.3f (> 0.3), -0.10",
                fraction_disordered,
            )

    overall -= disorder_penalty
    overall = round(max(0.0, min(1.0, overall)), 3)

    return {
        "overall": overall,
        "components": components,
        "disorder_penalty": disorder_penalty,
    }


def calculate_confidence_batch(
    molecules: list[dict],
    structure_source: str = "alphafold",
    disorder_info: Optional[dict] = None,
) -> list[dict]:
    """Compute confidence scores for a list of molecules.

    Parameters
    ----------
    molecules : list[dict]
        List of molecule result dicts.
    structure_source : str
        Source of the protein structure.
    disorder_info : dict, optional
        Disorder prediction result from ``predict_disorder()``.

    Returns
    -------
    list[dict]
        The input molecules with a ``confidence`` key added to each.
    """
    for mol in molecules:
        try:
            mol["confidence"] = calculate_confidence(mol, structure_source, disorder_info=disorder_info)
        except Exception as exc:
            logger.warning(
                "Confidence calculation failed for %s: %s",
                mol.get("name", "unknown"), exc,
            )
            mol["confidence"] = {
                "overall": 0.0,
                "components": {},
                "disorder_penalty": 0.0,
            }
    return molecules


# =====================================================================
# COMPONENT SCORING FUNCTIONS
# =====================================================================

def _score_structure(structure_source: str) -> tuple[float, str]:
    """Score confidence based on protein structure source.

    Parameters
    ----------
    structure_source : str
        One of: "alphafold", "pdb_experimental", "esmfold", "mock", "unknown".

    Returns
    -------
    tuple[float, str]
        (score, explanatory note)
    """
    source_lower = structure_source.lower().strip()

    # V6: Updated hierarchy -- PDB holo/apo distinction handled via pdb_info
    # in the caller; here we score by broad source category.
    if source_lower in ("pdb_holo",):
        return 0.98, "Experimental structure with co-crystallized ligand (holo) -- very high confidence"
    elif source_lower in ("pdb_experimental", "experimental", "pdb_apo"):
        return 0.90, "Experimental structure (X-ray/cryo-EM, apo) -- high confidence"
    elif source_lower == "alphafold":
        return 0.85, "AlphaFold predicted structure -- high confidence"
    elif source_lower == "esmfold":
        return 0.60, "ESMFold predicted structure -- moderate confidence"
    elif source_lower == "mock":
        return 0.20, "Mock structure -- low confidence (development only)"
    else:
        return 0.40, f"Unknown structure source '{structure_source}' -- limited confidence"


def _score_pocket(molecule: dict) -> tuple[float, str]:
    """Score confidence based on pocket detection quality.

    V5bis: supports co-crystallized ligand pockets, P2Rank ML, and fpocket.

    Parameters
    ----------
    molecule : dict
        Molecule result dict. Checked for ``pocket_method`` key.

    Returns
    -------
    tuple[float, str]
        (score, explanatory note)
    """
    pocket_method = (molecule.get("pocket_method") or "").lower()

    # V6: Updated pocket confidence hierarchy
    if "co-crystal" in pocket_method or "ligand" in pocket_method:
        return 0.98, "Pocket from co-crystallized ligand -- very high confidence"
    elif "p2rank" in pocket_method:
        return 0.85, "Pocket detected by P2Rank (ML-based) -- high confidence"
    elif "fpocket" in pocket_method:
        return 0.75, "Pocket detected by fpocket (geometric) -- moderate confidence"
    elif "mock" in pocket_method or "geometric" in pocket_method:
        return 0.50, "Mock/geometric pocket detection -- low confidence (development only)"
    else:
        return 0.75, "Pocket detected (method not specified)"


def _score_docking(molecule: dict) -> tuple[float, str]:
    """Score confidence based on docking method and molecule source.

    ChEMBL compounds (known bioactive molecules) are docked with higher
    confidence than de novo generated molecules, because their chemical
    space is better explored. The docking engine also affects confidence.

    Parameters
    ----------
    molecule : dict
        Must have ``source`` and optionally ``docking_method`` keys.

    Returns
    -------
    tuple[float, str]
        (score, explanatory note)
    """
    source = (molecule.get("source") or "").lower()
    method = (molecule.get("docking_method") or "vina").lower()

    # Source-based scoring
    if source in ("chembl", "user"):
        source_score = 0.70
        source_note = "Known bioactive compound"
    elif source in ("zinc", "zinc20"):
        source_score = 0.60
        source_note = "ZINC drug-like compound"
    elif source in ("mock_generation", "reinvent4", "generated"):
        source_score = 0.45
        source_note = "AI-generated compound (novel, less validated)"
    else:
        source_score = 0.50
        source_note = f"Source: {source or 'unknown'}"

    # Method-based adjustment (V5bis: GNINA > Vina)
    engine = (molecule.get("docking_engine") or method).lower()
    if engine == "gnina":
        method_adj = 0.10
        method_note = "GNINA (CNN-scored docking)"
    elif engine == "diffdock":
        method_adj = 0.05
        method_note = "DiffDock (AI docking)"
    elif engine == "vina":
        method_adj = 0.0
        method_note = "AutoDock Vina"
    elif engine == "mock":
        method_adj = -0.10
        method_note = "Mock docking (development only)"
    else:
        method_adj = -0.05
        method_note = f"Unknown method: {engine}"

    # V5bis: CNN score bonus (pose quality indicator)
    cnn_score = molecule.get("cnn_score", 0)
    if cnn_score and cnn_score > 0.7:
        method_adj += 0.05  # High CNN confidence bonus

    score = max(0.0, min(1.0, source_score + method_adj))
    note = f"{source_note} docked with {method_note}"

    return round(score, 3), note


def _score_admet(molecule: dict) -> tuple[float, str]:
    """Score confidence based on ADMET data availability and quality.

    Parameters
    ----------
    molecule : dict
        Must have optional ``admet`` key (dict or None).

    Returns
    -------
    tuple[float, str]
        (score, explanatory note)
    """
    admet = molecule.get("admet")

    if admet is None or not isinstance(admet, dict):
        return 0.30, "No ADMET data available -- limited safety assessment"

    # Check how many ADMET categories are populated
    categories = ["absorption", "distribution", "metabolism", "excretion", "toxicity"]
    populated = sum(
        1 for cat in categories
        if isinstance(admet.get(cat), dict) and len(admet[cat]) > 0
    )

    if populated >= 4:
        return 0.70, f"ADMET predictions available ({populated}/5 categories)"
    elif populated >= 2:
        return 0.55, f"Partial ADMET data ({populated}/5 categories)"
    else:
        # Has some data but incomplete
        composite = admet.get("composite_score")
        if composite is not None:
            return 0.50, "ADMET composite score available (limited detail)"
        return 0.35, "Minimal ADMET data available"


def _score_synthesis(molecule: dict) -> tuple[float, str]:
    """Score confidence based on retrosynthesis route availability.

    Parameters
    ----------
    molecule : dict
        Must have optional ``synthesis_route`` key (dict or None).

    Returns
    -------
    tuple[float, str]
        (score, explanatory note)
    """
    route = molecule.get("synthesis_route")

    if route is None or not isinstance(route, dict):
        return 0.20, "No synthesis route planned -- feasibility unknown"

    n_steps = route.get("n_steps", 0)
    confidence = route.get("confidence", 0.0)
    all_available = route.get("all_reagents_available", False)

    if n_steps == 0:
        return 0.20, "Synthesis route has 0 steps -- likely failed"

    if all_available and confidence >= 0.5:
        return 0.80, (
            f"Synthesis route: {n_steps} steps, confidence {confidence:.2f}, "
            f"all reagents available"
        )
    elif n_steps > 0 and confidence >= 0.3:
        return 0.50, (
            f"Synthesis route: {n_steps} steps, confidence {confidence:.2f}, "
            f"some reagents may not be available"
        )
    else:
        return 0.30, (
            f"Synthesis route: {n_steps} steps, low confidence ({confidence:.2f})"
        )


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

    print("=" * 70)
    print("DockIt Confidence Scoring -- Self-Test")
    print("=" * 70)

    # Test 1: High-confidence molecule (ChEMBL, AlphaFold, full data)
    print("\n[1] High-confidence molecule (ChEMBL + AlphaFold + ADMET + synthesis)...")
    mol_high = {
        "name": "Erlotinib",
        "smiles": "C=Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1",
        "source": "chembl",
        "docking_method": "vina",
        "admet": {
            "absorption": {"oral_bioavailability": 0.85},
            "distribution": {"bbb_permeability": 0.55},
            "metabolism": {"cyp3a4_inhibitor": 0.3},
            "excretion": {"clearance": 0.4},
            "toxicity": {"herg_inhibition": 0.1},
            "composite_score": 0.65,
        },
        "synthesis_route": {
            "n_steps": 3,
            "confidence": 0.75,
            "all_reagents_available": True,
        },
    }
    result_high = calculate_confidence(mol_high, structure_source="alphafold")
    print(f"  Overall confidence: {result_high['overall']:.3f}")
    for comp_name, comp_data in result_high["components"].items():
        print(f"  {comp_name:12s}: {comp_data['score']:.3f} -- {comp_data['note']}")
    assert result_high["overall"] >= 0.6, f"Expected >= 0.6, got {result_high['overall']}"
    print("  PASS: high-confidence molecule scores >= 0.6")

    # Test 2: Low-confidence molecule (generated, mock structure, no ADMET)
    print("\n[2] Low-confidence molecule (generated, mock, no ADMET, no synthesis)...")
    mol_low = {
        "name": "GEN_001",
        "smiles": "Cc1ccnc(Nc2ccc(CN3CCNCC3)cc2)c1C(=O)Nc1ccccc1F",
        "source": "mock_generation",
        "docking_method": "vina",
    }
    result_low = calculate_confidence(mol_low, structure_source="mock")
    print(f"  Overall confidence: {result_low['overall']:.3f}")
    for comp_name, comp_data in result_low["components"].items():
        print(f"  {comp_name:12s}: {comp_data['score']:.3f} -- {comp_data['note']}")
    assert result_low["overall"] < result_high["overall"], (
        "Low-confidence should be lower than high-confidence"
    )
    print("  PASS: low-confidence molecule scores lower than high-confidence")

    # Test 3: Batch calculation
    print("\n[3] Batch confidence calculation...")
    batch = [mol_high.copy(), mol_low.copy()]
    batch_result = calculate_confidence_batch(batch, structure_source="alphafold")
    assert len(batch_result) == 2
    assert "confidence" in batch_result[0]
    assert "confidence" in batch_result[1]
    print(f"  Molecule 1 confidence: {batch_result[0]['confidence']['overall']:.3f}")
    print(f"  Molecule 2 confidence: {batch_result[1]['confidence']['overall']:.3f}")
    print("  PASS: batch calculation works")

    # Test 4: ESMFold vs AlphaFold
    print("\n[4] ESMFold vs AlphaFold structure confidence...")
    result_esmfold = calculate_confidence(mol_high, structure_source="esmfold")
    result_alphafold = calculate_confidence(mol_high, structure_source="alphafold")
    print(f"  AlphaFold confidence: {result_alphafold['overall']:.3f}")
    print(f"  ESMFold confidence:   {result_esmfold['overall']:.3f}")
    assert result_alphafold["overall"] > result_esmfold["overall"], (
        "AlphaFold should give higher confidence than ESMFold"
    )
    print("  PASS: AlphaFold > ESMFold confidence")

    # Test 5: V6 disorder penalty
    print("\n[5] V6 disorder penalty...")
    disorder_high = {"fraction_disordered": 0.45, "idr_regions": [(10, 100)], "method": "mock"}
    disorder_low = {"fraction_disordered": 0.10, "idr_regions": [], "method": "mock"}
    result_disordered = calculate_confidence(mol_high, structure_source="alphafold", disorder_info=disorder_high)
    result_ordered = calculate_confidence(mol_high, structure_source="alphafold", disorder_info=disorder_low)
    print(f"  High disorder confidence: {result_disordered['overall']:.3f} (penalty={result_disordered['disorder_penalty']})")
    print(f"  Low disorder confidence:  {result_ordered['overall']:.3f} (penalty={result_ordered['disorder_penalty']})")
    assert result_disordered["disorder_penalty"] == 0.10, "Expected 0.10 penalty for high disorder"
    assert result_ordered["disorder_penalty"] == 0.0, "Expected no penalty for low disorder"
    assert result_disordered["overall"] < result_ordered["overall"], (
        "High-disorder protein should have lower confidence"
    )
    print("  PASS: disorder penalty applied correctly")

    # Test 6: V6 PDB holo score
    print("\n[6] V6 PDB holo vs apo structure confidence...")
    result_holo = calculate_confidence(mol_high, structure_source="pdb_holo")
    result_apo = calculate_confidence(mol_high, structure_source="pdb_experimental")
    print(f"  PDB holo confidence:  {result_holo['overall']:.3f}")
    print(f"  PDB apo confidence:   {result_apo['overall']:.3f}")
    assert result_holo["overall"] >= result_apo["overall"], (
        "PDB holo should give >= confidence than PDB apo"
    )
    print("  PASS: PDB holo >= PDB apo confidence")

    print("\n" + "=" * 70)
    print("ALL TESTS PASSED")
    print("=" * 70)
