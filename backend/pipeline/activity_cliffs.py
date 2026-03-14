"""
BindX V9 — Activity Cliff Detection (SALI).

Identifies structure-activity landscape index (SALI) pairs:
molecules that are structurally similar but have large activity differences.

SALI = |activity_i - activity_j| / (1 - Tc(i,j))

Reference: Guha & Van Drie, J. Chem. Inf. Model. 2008, 48, 646-658.
Thresholds from bindx_reference_visuel.html §P17:
  - Tc >= 0.40 on ECFP4 for structural similarity
  - delta_activity >= 2 log units (100x difference)
  - SALI percentile 90-95 = cliff threshold
"""

from __future__ import annotations

import logging
from typing import Optional

logger = logging.getLogger(__name__)


def detect_activity_cliffs(
    molecules: list[dict],
    activity_key: str = "docking_score",
    tc_threshold: float = 0.40,
    activity_threshold: float = 2.0,
    top_n: int = 20,
) -> list[dict]:
    """Detect activity cliffs in a set of molecules.

    Parameters
    ----------
    molecules : list[dict]
        Each must have 'smiles' and the specified activity_key.
    activity_key : str
        Property key to use as activity (default: docking_score).
    tc_threshold : float
        Minimum Tanimoto similarity to be considered similar (default: 0.40).
    activity_threshold : float
        Minimum activity difference to qualify as cliff (default: 2.0).
    top_n : int
        Return at most top_n cliffs per molecule.

    Returns
    -------
    list[dict]
        Per-molecule cliff annotations:
        {
            is_cliff: bool,
            sali_max: float,
            cliff_partner: str (SMILES of most extreme partner),
            cliff_delta: float (activity difference with partner),
            n_cliffs: int (number of cliff pairs),
        }
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, DataStructs
    except ImportError:
        logger.warning("RDKit not available for activity cliff detection")
        return [_nd_result() for _ in molecules]

    # Filter molecules with valid SMILES and activity
    valid = []
    has_activity = set()  # Track which molecules have activity data
    for i, mol in enumerate(molecules):
        smi = mol.get("smiles")
        act = mol.get(activity_key)
        if smi and act is not None:
            has_activity.add(i)
            try:
                act_val = float(act)
                rdmol = Chem.MolFromSmiles(smi)
                if rdmol is not None:
                    fp = AllChem.GetMorganFingerprintAsBitVect(rdmol, 2, nBits=2048)
                    valid.append({"idx": i, "smiles": smi, "activity": act_val, "fp": fp})
            except (ValueError, TypeError):
                continue

    if len(valid) < 2:
        # Distinguish: molecules WITH activity but < 2 valid → "not a cliff"
        #              molecules WITHOUT activity → N/D
        return [
            _no_cliff_result() if i in has_activity else _nd_result()
            for i in range(len(molecules))
        ]

    # Compute pairwise SALI
    n = len(valid)
    cliff_data = [[] for _ in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            tc = DataStructs.TanimotoSimilarity(valid[i]["fp"], valid[j]["fp"])
            if tc < tc_threshold:
                continue
            if tc >= 1.0:
                continue  # identical molecules

            delta = abs(valid[i]["activity"] - valid[j]["activity"])
            if delta < activity_threshold:
                continue

            sali = delta / (1.0 - tc)
            cliff_data[i].append({"partner_idx": j, "sali": sali, "delta": delta, "tc": tc})
            cliff_data[j].append({"partner_idx": i, "sali": sali, "delta": delta, "tc": tc})

    # Build per-molecule results
    # Molecules with activity data → evaluated (is_cliff: True/False)
    # Molecules without activity data → N/D (is_cliff: None)
    results = [
        _no_cliff_result() if i in has_activity else _nd_result()
        for i in range(len(molecules))
    ]

    for vi, cliffs in enumerate(cliff_data):
        orig_idx = valid[vi]["idx"]
        if not cliffs:
            continue

        cliffs.sort(key=lambda c: c["sali"], reverse=True)
        best = cliffs[0]
        partner = valid[best["partner_idx"]]

        results[orig_idx] = {
            "is_cliff": True,
            "sali_max": round(best["sali"], 2),
            "cliff_partner": partner["smiles"],
            "cliff_delta": round(best["delta"], 2),
            "n_cliffs": min(len(cliffs), top_n),
        }

    return results


def _no_cliff_result() -> dict:
    """Molecule was evaluated but is not an activity cliff."""
    return {
        "is_cliff": False,
        "sali_max": None,
        "cliff_partner": None,
        "cliff_delta": None,
        "n_cliffs": 0,
    }


def _nd_result() -> dict:
    """Molecule could not be evaluated (missing activity data / docking not run)."""
    return {
        "is_cliff": None,
        "sali_max": None,
        "cliff_partner": None,
        "cliff_delta": None,
        "n_cliffs": None,
    }
