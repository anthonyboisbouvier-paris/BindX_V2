"""
DockIt pipeline -- Pharmacological filter (Pass 1 of massive screening).

Filters a large SMILES pool (up to ~2.4M molecules) down to ~200K candidates
using drug-likeness criteria. This is the first and broadest filter, designed
to eliminate molecules that have no chance of being drug-like before more
expensive 3D and docking calculations.

Filter criteria (RDKit full mode):
  - Lipinski Rule of 5 with max 1 violation:
      MW <= 500, LogP <= 5, HBD <= 5, HBA <= 10
  - QED (Quantitative Estimate of Drug-likeness) > 0.3
  - PAINS substructure rejection (10+ SMARTS patterns)
  - MW compatible with pocket volume (heuristic: MW < pocket_volume * 0.8)

Strategies (3-tier fallback):
  1. RDKit full: Lipinski + QED + PAINS + pocket MW constraint.
  2. RDKit basic: Lipinski only (when QED/PAINS imports fail).
  3. Hash-based mock: deterministic ~10% pass rate (no cheminformatics).
"""

from __future__ import annotations

import hashlib
import logging
import time
from typing import Callable, Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Availability flags (resolved lazily, once)
# ---------------------------------------------------------------------------
_RDKIT_AVAILABLE: Optional[bool] = None
_RDKIT_QED_AVAILABLE: Optional[bool] = None


def _check_rdkit() -> bool:
    """Check whether RDKit core can be imported."""
    global _RDKIT_AVAILABLE
    if _RDKIT_AVAILABLE is None:
        try:
            from rdkit import Chem  # noqa: F401
            from rdkit.Chem import Descriptors  # noqa: F401
            _RDKIT_AVAILABLE = True
        except ImportError:
            _RDKIT_AVAILABLE = False
            logger.warning(
                "RDKit not available; pharmacological filter will use "
                "hash-based mock (~10%% pass rate)"
            )
    return _RDKIT_AVAILABLE


def _check_rdkit_qed() -> bool:
    """Check whether RDKit QED module is available."""
    global _RDKIT_QED_AVAILABLE
    if _RDKIT_QED_AVAILABLE is None:
        try:
            from rdkit.Chem import QED  # noqa: F401
            _RDKIT_QED_AVAILABLE = True
        except ImportError:
            _RDKIT_QED_AVAILABLE = False
            logger.info("RDKit QED not available; will skip QED filter")
    return _RDKIT_QED_AVAILABLE


# ---------------------------------------------------------------------------
# PAINS SMARTS patterns (curated subset of the most common alerts)
# ---------------------------------------------------------------------------

# At least 10 SMARTS patterns covering the most frequent PAINS families.
# These are consistent with the patterns in generation.py plus extras.
_PAINS_SMARTS: list[str] = [
    # 1. Azo compounds (diazo dyes)
    "[#6]1:[#6]:[#6](:[#6]:[#6]:[#6]:1)-[#7]=[#7]-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2",
    # 2. Disulfide-like
    "[#6]-[#16]-[#16]-[#6]",
    # 3. Rhodanine-like (thiazolidine-2,4-dione)
    "[#6]=[#6](-[#7])-[#16]",
    # 4. Rhodanine core
    "[#8]=[#6]-1-[#6]=[#6]-[#16]-[#7]-1",
    # 5. Quinone (oxidized catechol)
    "[#8]=[#6]1[#6]=[#6][#6](=[#8])[#6]=[#6]1",
    # 6. Maleic anhydride
    "[#6]1=[#6]-[#6](=[#8])-[#8]-[#6](=[#8])-1",
    # 7. Michael acceptor: alpha-beta unsaturated carbonyl
    "[#6]=[#6]-[#6](=[#8])",
    # 8. Acyl hydrazide
    "[#6](=[#8])-[#7]-[#7]",
    # 9. Alkylidene barbiturate / cyanoacrylate
    "[#6]=[#6](-[#6]#[#7])-[#6](=[#8])",
    # 10. Catechol (1,2-dihydroxybenzene) -- prone to oxidation
    "[OH]c1ccccc1[OH]",
    # 11. Nitro aromatic
    "[c][N+](=O)[O-]",
    # 12. Acyl halide (reactive electrophile)
    "[#6](=[#8])-[F,Cl,Br,I]",
    # 13. Epoxide (strained reactive ring)
    "C1OC1",
    # 14. Isothiocyanate / isocyanate
    "[#7]=[#6]=[#16]",
    # 15. Aldehyde (often reactive, metabolic liability)
    "[CX3H1](=O)[#6]",
]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def filter_pharmacological(
    smiles_list: list[str],
    pocket_volume: float = 800.0,
    progress_callback: Optional[Callable[[dict], None]] = None,
) -> list[str]:
    """Filter a SMILES list by pharmacological drug-likeness criteria.

    This is Pass 1 of the massive screening pipeline. It eliminates
    molecules that are clearly not drug-like based on physicochemical
    properties, quality scores, and reactive substructure alerts.

    Parameters
    ----------
    smiles_list : list[str]
        Input SMILES strings (may be millions).
    pocket_volume : float
        Estimated binding pocket volume in Angstrom^3. Used to set a
        soft upper bound on molecular weight via the heuristic
        MW < pocket_volume * 0.8.
    progress_callback : callable, optional
        Called with a dict ``{"pass_name": str, "processed": int,
        "total": int, "kept": int}`` after each batch.

    Returns
    -------
    list[str]
        Filtered SMILES that pass all drug-likeness criteria.

    Notes
    -----
    Processing is done in batches of 10000 to limit memory usage when
    dealing with millions of molecules.
    """
    total = len(smiles_list)
    if total == 0:
        logger.info("filter_pharmacological: empty input, nothing to filter")
        return []

    logger.info(
        "filter_pharmacological: starting with %d molecules, "
        "pocket_volume=%.1f",
        total, pocket_volume,
    )
    start_time = time.monotonic()

    # Determine the MW ceiling from pocket volume
    mw_ceiling = pocket_volume * 0.8
    # Clamp to a sensible range: at least 200, at most 800
    mw_ceiling = max(200.0, min(800.0, mw_ceiling))
    logger.info("MW ceiling from pocket volume: %.1f Da", mw_ceiling)

    # --- Strategy 1: RDKit full (Lipinski + QED + PAINS + pocket MW) ---
    if _check_rdkit():
        try:
            result = _filter_rdkit_full(
                smiles_list, mw_ceiling, progress_callback,
            )
            elapsed = time.monotonic() - start_time
            logger.info(
                "filter_pharmacological (RDKit full): %d -> %d molecules "
                "(%.1f%% pass rate) in %.1fs",
                total, len(result),
                len(result) / max(total, 1) * 100,
                elapsed,
            )
            return result
        except Exception as exc:
            logger.warning(
                "RDKit full filter failed: %s; trying basic fallback", exc
            )

    # --- Strategy 2: RDKit basic (Lipinski only) ---
    if _check_rdkit():
        try:
            result = _filter_rdkit_basic(
                smiles_list, mw_ceiling, progress_callback,
            )
            elapsed = time.monotonic() - start_time
            logger.info(
                "filter_pharmacological (RDKit basic): %d -> %d molecules "
                "(%.1f%% pass rate) in %.1fs",
                total, len(result),
                len(result) / max(total, 1) * 100,
                elapsed,
            )
            return result
        except Exception as exc:
            logger.warning(
                "RDKit basic filter also failed: %s; falling back to mock",
                exc,
            )

    # --- Strategy 3: Hash-based mock (~10% pass rate) ---
    result = _filter_mock(smiles_list, progress_callback)
    elapsed = time.monotonic() - start_time
    logger.info(
        "filter_pharmacological (mock): %d -> %d molecules "
        "(%.1f%% pass rate) in %.1fs",
        total, len(result),
        len(result) / max(total, 1) * 100,
        elapsed,
    )
    return result


# =========================================================================
# STRATEGY 1: RDKit full filter
# =========================================================================

def _filter_rdkit_full(
    smiles_list: list[str],
    mw_ceiling: float,
    progress_callback: Optional[Callable[[dict], None]],
) -> list[str]:
    """Full RDKit filter: Lipinski + QED + PAINS + pocket MW constraint.

    Parameters
    ----------
    smiles_list : list[str]
        Input SMILES.
    mw_ceiling : float
        Maximum molecular weight derived from pocket volume.
    progress_callback : callable, optional
        Progress reporting callback.

    Returns
    -------
    list[str]
        SMILES passing all criteria.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, rdMolDescriptors

    # Conditionally import QED
    has_qed = _check_rdkit_qed()
    if has_qed:
        from rdkit.Chem import QED as QED_module
    else:
        QED_module = None  # type: ignore[assignment]

    # Suppress noisy RDKit warnings during batch processing
    RDLogger.DisableLog("rdApp.*")

    # Pre-compile PAINS patterns
    pains_patterns = _compile_pains_patterns(Chem)

    # Rejection reason counters
    stats: dict[str, int] = {
        "invalid_smiles": 0,
        "mw_over_500": 0,
        "logp_over_5": 0,
        "hbd_over_5": 0,
        "hba_over_10": 0,
        "lipinski_multi_violation": 0,
        "qed_low": 0,
        "pains_hit": 0,
        "mw_pocket_incompatible": 0,
        "mw_under_100": 0,
    }

    passed: list[str] = []
    batch_size = 10000
    total = len(smiles_list)

    for batch_start in range(0, total, batch_size):
        batch_end = min(batch_start + batch_size, total)
        batch = smiles_list[batch_start:batch_end]

        for smi in batch:
            if not smi or not smi.strip():
                stats["invalid_smiles"] += 1
                continue

            smi = smi.strip()

            # Parse SMILES
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                stats["invalid_smiles"] += 1
                continue

            # --- Lipinski Rule of 5 ---
            try:
                mw = Descriptors.ExactMolWt(mol)
            except Exception:
                stats["invalid_smiles"] += 1
                continue

            # Quick MW checks first (cheapest filter)
            if mw < 100.0:
                stats["mw_under_100"] += 1
                continue
            if mw > mw_ceiling:
                stats["mw_pocket_incompatible"] += 1
                continue

            try:
                logp = Descriptors.MolLogP(mol)
                hbd = rdMolDescriptors.CalcNumHBD(mol)
                hba = rdMolDescriptors.CalcNumHBA(mol)
            except Exception:
                stats["invalid_smiles"] += 1
                continue

            # Count Lipinski violations (max 1 allowed)
            violations = 0
            if mw > 500.0:
                violations += 1
                stats["mw_over_500"] += 1
            if logp > 5.0:
                violations += 1
                stats["logp_over_5"] += 1
            if hbd > 5:
                violations += 1
                stats["hbd_over_5"] += 1
            if hba > 10:
                violations += 1
                stats["hba_over_10"] += 1

            if violations > 1:
                stats["lipinski_multi_violation"] += 1
                continue

            # --- QED > 0.3 ---
            if has_qed and QED_module is not None:
                try:
                    qed_val = QED_module.qed(mol)
                    if qed_val < 0.3:
                        stats["qed_low"] += 1
                        continue
                except Exception:
                    # If QED fails, skip this filter rather than reject
                    pass

            # --- PAINS check ---
            if _has_pains_match(mol, pains_patterns):
                stats["pains_hit"] += 1
                continue

            # Molecule passed all filters
            canonical = Chem.MolToSmiles(mol)
            passed.append(canonical)

        # Report progress
        if progress_callback is not None:
            try:
                progress_callback({
                    "pass_name": "pharmacological_filter",
                    "processed": batch_end,
                    "total": total,
                    "kept": len(passed),
                })
            except Exception:
                pass

    # Restore RDKit logging
    try:
        RDLogger.EnableLog("rdApp.*")
    except Exception:
        pass

    # Log rejection stats
    _log_rejection_stats(stats, total, len(passed))

    return passed


# =========================================================================
# STRATEGY 2: RDKit basic filter (Lipinski only)
# =========================================================================

def _filter_rdkit_basic(
    smiles_list: list[str],
    mw_ceiling: float,
    progress_callback: Optional[Callable[[dict], None]],
) -> list[str]:
    """Basic RDKit filter: Lipinski Rule of 5 only, no QED or PAINS.

    This is the fallback when QED or PAINS-related imports fail.

    Parameters
    ----------
    smiles_list : list[str]
        Input SMILES.
    mw_ceiling : float
        Maximum molecular weight.
    progress_callback : callable, optional
        Progress reporting callback.

    Returns
    -------
    list[str]
        SMILES passing Lipinski with max 1 violation and MW < ceiling.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, rdMolDescriptors

    RDLogger.DisableLog("rdApp.*")

    stats: dict[str, int] = {
        "invalid_smiles": 0,
        "lipinski_multi_violation": 0,
        "mw_pocket_incompatible": 0,
        "mw_under_100": 0,
    }

    passed: list[str] = []
    batch_size = 10000
    total = len(smiles_list)

    for batch_start in range(0, total, batch_size):
        batch_end = min(batch_start + batch_size, total)
        batch = smiles_list[batch_start:batch_end]

        for smi in batch:
            if not smi or not smi.strip():
                stats["invalid_smiles"] += 1
                continue

            smi = smi.strip()
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                stats["invalid_smiles"] += 1
                continue

            try:
                mw = Descriptors.ExactMolWt(mol)
                if mw < 100.0:
                    stats["mw_under_100"] += 1
                    continue
                if mw > mw_ceiling:
                    stats["mw_pocket_incompatible"] += 1
                    continue

                logp = Descriptors.MolLogP(mol)
                hbd = rdMolDescriptors.CalcNumHBD(mol)
                hba = rdMolDescriptors.CalcNumHBA(mol)
            except Exception:
                stats["invalid_smiles"] += 1
                continue

            violations = 0
            if mw > 500.0:
                violations += 1
            if logp > 5.0:
                violations += 1
            if hbd > 5:
                violations += 1
            if hba > 10:
                violations += 1

            if violations > 1:
                stats["lipinski_multi_violation"] += 1
                continue

            canonical = Chem.MolToSmiles(mol)
            passed.append(canonical)

        if progress_callback is not None:
            try:
                progress_callback({
                    "pass_name": "pharmacological_filter_basic",
                    "processed": batch_end,
                    "total": total,
                    "kept": len(passed),
                })
            except Exception:
                pass

    try:
        RDLogger.EnableLog("rdApp.*")
    except Exception:
        pass

    _log_rejection_stats(stats, total, len(passed))
    return passed


# =========================================================================
# STRATEGY 3: Hash-based mock (~10% deterministic pass rate)
# =========================================================================

def _filter_mock(
    smiles_list: list[str],
    progress_callback: Optional[Callable[[dict], None]],
) -> list[str]:
    """Deterministic hash-based mock filter.

    Keeps approximately 10% of input molecules using a hash of the
    SMILES string. The same molecule always produces the same pass/fail
    decision, ensuring reproducibility.

    Parameters
    ----------
    smiles_list : list[str]
        Input SMILES.
    progress_callback : callable, optional
        Progress reporting callback.

    Returns
    -------
    list[str]
        Approximately 10% of input SMILES, deterministically selected.
    """
    logger.info(
        "filter_pharmacological mock: processing %d molecules "
        "(deterministic ~10%% pass rate)",
        len(smiles_list),
    )

    passed: list[str] = []
    batch_size = 10000
    total = len(smiles_list)

    for batch_start in range(0, total, batch_size):
        batch_end = min(batch_start + batch_size, total)
        batch = smiles_list[batch_start:batch_end]

        for smi in batch:
            if not smi or not smi.strip():
                continue

            smi = smi.strip()

            # Deterministic hash-based decision
            digest = hashlib.sha256(smi.encode("utf-8")).hexdigest()
            hash_val = int(digest[:8], 16)

            # Keep ~10% (hash modulo 10 == 0)
            if hash_val % 10 == 0:
                passed.append(smi)

        if progress_callback is not None:
            try:
                progress_callback({
                    "pass_name": "pharmacological_filter_mock",
                    "processed": batch_end,
                    "total": total,
                    "kept": len(passed),
                })
            except Exception:
                pass

    logger.info(
        "Mock pharma filter: %d -> %d (%.1f%%)",
        total, len(passed),
        len(passed) / max(total, 1) * 100,
    )
    return passed


# =========================================================================
# PAINS helpers
# =========================================================================

def _compile_pains_patterns(Chem: object) -> list:
    """Compile PAINS SMARTS patterns into RDKit Mol objects.

    Parameters
    ----------
    Chem : module
        The rdkit.Chem module.

    Returns
    -------
    list
        Compiled SMARTS pattern Mol objects. Empty list on failure,
        which means the PAINS filter is simply skipped.
    """
    compiled: list = []
    for smarts in _PAINS_SMARTS:
        try:
            pat = Chem.MolFromSmarts(smarts)  # type: ignore[union-attr]
            if pat is not None:
                compiled.append(pat)
        except Exception:
            continue

    if compiled:
        logger.debug("Compiled %d/%d PAINS SMARTS patterns", len(compiled), len(_PAINS_SMARTS))
    else:
        logger.warning("No PAINS patterns compiled; PAINS filter will be skipped")

    return compiled


def _has_pains_match(mol: object, patterns: list) -> bool:
    """Check if a molecule matches any pre-compiled PAINS pattern.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        The molecule to check.
    patterns : list
        Pre-compiled PAINS SMARTS Mol objects.

    Returns
    -------
    bool
        True if any PAINS pattern matches.
    """
    if not patterns:
        return False
    try:
        for pat in patterns:
            if mol.HasSubstructMatch(pat):  # type: ignore[union-attr]
                return True
        return False
    except Exception:
        return False


# =========================================================================
# Logging helpers
# =========================================================================

def _log_rejection_stats(
    stats: dict[str, int],
    total_input: int,
    total_passed: int,
) -> None:
    """Log a summary of rejection reasons.

    Parameters
    ----------
    stats : dict[str, int]
        Counter dict with rejection reason keys.
    total_input : int
        Total molecules submitted.
    total_passed : int
        Total molecules that passed.
    """
    total_rejected = total_input - total_passed
    logger.info(
        "Pharmacological filter stats: %d input, %d passed (%.1f%%), "
        "%d rejected",
        total_input, total_passed,
        total_passed / max(total_input, 1) * 100,
        total_rejected,
    )

    # Log individual rejection reasons (only non-zero)
    for reason, count in sorted(stats.items(), key=lambda x: -x[1]):
        if count > 0:
            logger.info(
                "  Rejection: %-30s  %7d  (%5.1f%%)",
                reason, count,
                count / max(total_input, 1) * 100,
            )


# =========================================================================
# CLI / Self-test
# =========================================================================

if __name__ == "__main__":
    import sys

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    print("=" * 70)
    print("DockIt -- filter_pharma.py self-test")
    print("=" * 70)

    # Test SMILES: mix of drug-like and non-drug-like
    test_smiles = [
        # Drug-like (should pass)
        "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OC)cc23)c1",     # Erlotinib-like
        "CC(=O)Oc1ccccc1C(=O)O",                         # Aspirin
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",                    # Ibuprofen
        "Cn1c(=O)c2c(ncn2C)n(C)c1=O",                    # Caffeine
        "CC(=O)Nc1ccc(O)cc1",                             # Acetaminophen
        # Non-drug-like (should fail)
        "C",                                               # Too small (methane)
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",            # Very long chain (MW > 500)
        "OC(=O)c1cc(O)c(O)c(O)c1",                       # Catechol (PAINS)
        "",                                                # Empty
        "invalid_smiles_xyz",                              # Invalid
    ]

    print(f"\nTesting with {len(test_smiles)} molecules...")

    progress_log: list[dict] = []

    def on_progress(info: dict) -> None:
        progress_log.append(info)
        print(f"  Progress: {info}")

    result = filter_pharmacological(
        test_smiles,
        pocket_volume=800.0,
        progress_callback=on_progress,
    )

    print(f"\nResult: {len(result)} molecules passed")
    for smi in result:
        print(f"  PASS: {smi}")

    # Basic assertions
    errors = 0

    if len(result) == 0:
        print("WARNING: No molecules passed (may be expected if RDKit not available)")
    elif len(result) > len(test_smiles):
        print("FAIL: More molecules passed than were input")
        errors += 1

    # Empty string should never pass
    if "" in result:
        print("FAIL: Empty string passed the filter")
        errors += 1

    # Check progress callback was called
    if progress_log:
        print(f"\nProgress callbacks received: {len(progress_log)}")
    else:
        print("WARNING: No progress callbacks received")

    if errors == 0:
        print("\nALL CHECKS PASSED")
    else:
        print(f"\n{errors} CHECK(S) FAILED")

    sys.exit(1 if errors > 0 else 0)
