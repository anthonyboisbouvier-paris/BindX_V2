"""Tests for NaN safety in scoring.py helpers and consensus detail."""

import math
import sys
import os

# Ensure the backend package is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pipeline.scoring import (
    _safe_float,
    _safe_zscore,
    enrich_consensus_detail,
    _extract_pareto_objectives,
)


# ---------------------------------------------------------------------------
# _safe_float
# ---------------------------------------------------------------------------

class TestSafeFloat:
    def test_normal_value(self):
        assert _safe_float(3.14) == 3.14

    def test_none_returns_zero(self):
        assert _safe_float(None) == 0.0

    def test_nan_returns_zero(self):
        assert _safe_float(float("nan")) == 0.0

    def test_string_returns_zero(self):
        assert _safe_float("not_a_number") == 0.0

    def test_zero(self):
        assert _safe_float(0.0) == 0.0

    def test_negative(self):
        assert _safe_float(-7.5) == -7.5

    def test_int(self):
        assert _safe_float(42) == 42.0

    def test_inf(self):
        # inf is a valid float, not NaN, so it should pass through
        assert _safe_float(float("inf")) == float("inf")


# ---------------------------------------------------------------------------
# _safe_zscore
# ---------------------------------------------------------------------------

class TestSafeZscore:
    def test_normal_computation(self):
        z = _safe_zscore(10.0, 5.0, 2.5)
        assert z == 2.0

    def test_zero_sd_returns_zero(self):
        assert _safe_zscore(10.0, 5.0, 0.0) == 0.0

    def test_negative_sd_returns_zero(self):
        assert _safe_zscore(10.0, 5.0, -1.0) == 0.0

    def test_nan_val_returns_zero(self):
        assert _safe_zscore(float("nan"), 0.0, 1.0) == 0.0

    def test_nan_mu_returns_zero(self):
        assert _safe_zscore(5.0, float("nan"), 1.0) == 0.0

    def test_result_is_rounded(self):
        z = _safe_zscore(1.0, 0.0, 3.0)
        assert z == 0.333


# ---------------------------------------------------------------------------
# enrich_consensus_detail with NaN inputs
# ---------------------------------------------------------------------------

class TestEnrichConsensusDetailNaN:
    def test_nan_values_produce_zero_zscores(self):
        """Molecules with NaN scores should get 0.0 z-scores, not NaN."""
        molecules = [
            {"vina_score": float("nan"), "cnn_score": 0.5, "cnn_affinity": -6.0},
            {"vina_score": -7.0, "cnn_score": float("nan"), "cnn_affinity": -8.0},
            {"vina_score": -9.0, "cnn_score": 0.8, "cnn_affinity": float("nan")},
        ]
        result = enrich_consensus_detail(molecules)
        for mol in result:
            cd = mol["consensus_detail"]
            assert not math.isnan(cd["z_vina"]), f"z_vina is NaN for {mol}"
            assert not math.isnan(cd["z_cnn_score"]), f"z_cnn_score is NaN for {mol}"
            assert not math.isnan(cd["z_cnn_affinity"]), f"z_cnn_affinity is NaN for {mol}"

    def test_all_nan_values(self):
        """All NaN should produce valid (zero) z-scores."""
        molecules = [
            {"vina_score": float("nan"), "cnn_score": float("nan"), "cnn_affinity": float("nan")},
            {"vina_score": float("nan"), "cnn_score": float("nan"), "cnn_affinity": float("nan")},
        ]
        result = enrich_consensus_detail(molecules)
        for mol in result:
            cd = mol["consensus_detail"]
            assert cd["z_vina"] == 0.0
            assert cd["z_cnn_score"] == 0.0
            assert cd["z_cnn_affinity"] == 0.0
            # Agreement should also be valid
            assert cd["agreement"] in ("0/3", "1/3", "2/3", "3/3")

    def test_none_values(self):
        """Missing keys (None) should also be handled safely."""
        molecules = [
            {},  # all keys missing
            {"vina_score": None, "cnn_score": None, "cnn_affinity": None},
        ]
        result = enrich_consensus_detail(molecules)
        for mol in result:
            cd = mol["consensus_detail"]
            assert not math.isnan(cd["z_vina"])
            assert not math.isnan(cd["z_cnn_score"])
            assert not math.isnan(cd["z_cnn_affinity"])

    def test_normal_values_still_work(self):
        """Normal numeric values should produce meaningful z-scores."""
        molecules = [
            {"vina_score": -10.0, "cnn_score": 0.9, "cnn_affinity": -9.0},
            {"vina_score": -5.0, "cnn_score": 0.3, "cnn_affinity": -4.0},
        ]
        result = enrich_consensus_detail(molecules)
        # The better molecule (lower vina, higher cnn) should have positive z_vina
        cd0 = result[0]["consensus_detail"]
        cd1 = result[1]["consensus_detail"]
        assert cd0["z_vina"] > cd1["z_vina"]  # -10 is better
        assert cd0["z_cnn_score"] > cd1["z_cnn_score"]  # 0.9 > 0.3

    def test_empty_list(self):
        """Empty list should be a no-op."""
        assert enrich_consensus_detail([]) == []


# ---------------------------------------------------------------------------
# _extract_pareto_objectives with NaN inputs
# ---------------------------------------------------------------------------

class TestExtractParetoObjectivesNaN:
    def test_nan_vina_score(self):
        mol = {"vina_score": float("nan")}
        obj = _extract_pareto_objectives(mol, 10)
        assert not math.isnan(obj["affinity"])

    def test_nan_affinity_fallback(self):
        mol = {"affinity": float("nan")}
        obj = _extract_pareto_objectives(mol, 10)
        assert not math.isnan(obj["affinity"])

    def test_normal_vina_score(self):
        mol = {"vina_score": -8.0}
        obj = _extract_pareto_objectives(mol, 10)
        expected = round(8.0 / 12.0, 4)
        assert obj["affinity"] == expected
