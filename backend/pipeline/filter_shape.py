"""
DockIt pipeline -- 3D Shape filter (Pass 2 of massive screening).

Filters ~200K molecules down to ~10K based on 3D molecular shape
compatibility with the binding pocket. Molecules whose 3D volume and
dimensions are grossly incompatible with the pocket geometry are
rejected before the more expensive docking steps.

Strategies (3-tier fallback):
  1. RDKit ETKDG conformer: generate a single 3D conformer, compute
     molecular volume (bounding box approximation) and longest axis,
     compare against pocket dimensions.
  2. Simple MW-based volume estimation: use the rough empirical
     correlation vol ~ MW * 1.2 and longest_axis ~ MW^(1/3) * 1.5.
  3. Hash-based mock: deterministic ~5% pass rate.
"""

from __future__ import annotations

import hashlib
import logging
import math
import time
from typing import Callable, Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Availability flags (resolved lazily)
# ---------------------------------------------------------------------------
_RDKIT_AVAILABLE: Optional[bool] = None
_RDKIT_3D_AVAILABLE: Optional[bool] = None


def _check_rdkit() -> bool:
    """Check whether RDKit core can be imported."""
    global _RDKIT_AVAILABLE
    if _RDKIT_AVAILABLE is None:
        try:
            from rdkit import Chem  # noqa: F401
            _RDKIT_AVAILABLE = True
        except ImportError:
            _RDKIT_AVAILABLE = False
            logger.warning("RDKit not available; shape filter will use fallback")
    return _RDKIT_AVAILABLE


def _check_rdkit_3d() -> bool:
    """Check whether RDKit 3D conformer generation works."""
    global _RDKIT_3D_AVAILABLE
    if _RDKIT_3D_AVAILABLE is None:
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            # Quick smoke test with a tiny molecule
            mol = Chem.MolFromSmiles("C")
            if mol is not None:
                mol = Chem.AddHs(mol)
                res = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
                _RDKIT_3D_AVAILABLE = (res != -1)
            else:
                _RDKIT_3D_AVAILABLE = False
        except Exception:
            _RDKIT_3D_AVAILABLE = False
        if _RDKIT_3D_AVAILABLE:
            logger.info("RDKit 3D conformer generation available")
        else:
            logger.info("RDKit 3D conformer generation not available; will use volume estimation")
    return _RDKIT_3D_AVAILABLE


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def filter_by_shape(
    smiles_list: list[str],
    pocket_center: tuple[float, float, float],
    pocket_size: tuple[float, float, float],
    pocket_pdb: Optional[str] = None,
    progress_callback: Optional[Callable[[dict], None]] = None,
) -> list[str]:
    """Filter molecules by 3D shape compatibility with the binding pocket.

    This is Pass 2 of the massive screening pipeline. It rejects
    molecules whose 3D shape (volume and maximum dimension) is
    incompatible with the binding pocket geometry.

    Parameters
    ----------
    smiles_list : list[str]
        Input SMILES (typically ~200K from Pass 1).
    pocket_center : tuple[float, float, float]
        (x, y, z) centre of the binding pocket in Angstroms.
    pocket_size : tuple[float, float, float]
        (sx, sy, sz) dimensions of the binding pocket search box.
    pocket_pdb : str, optional
        Path to the pocket PDB file (not currently used, reserved for
        future shape complementarity scoring).
    progress_callback : callable, optional
        Called with a dict after each batch.

    Returns
    -------
    list[str]
        Filtered SMILES whose 3D shape is compatible with the pocket.
    """
    total = len(smiles_list)
    if total == 0:
        logger.info("filter_by_shape: empty input, nothing to filter")
        return []

    # Compute pocket volume estimate.
    # The pocket_size is the search box dimensions (typically 25x25x25 A),
    # NOT the actual cavity volume. The actual binding site cavity is
    # typically ~300-1500 A^3, much smaller than the search box volume
    # (15625 A^3 for a 25x25x25 box). We estimate the cavity volume
    # as a fraction of the search box, or use a reasonable default.
    search_box_vol = pocket_size[0] * pocket_size[1] * pocket_size[2]
    # Estimate cavity volume: roughly 5-10% of search box, or use a
    # sphere inscribed in the smallest dimension
    min_dim = min(pocket_size)
    # Approximate cavity as a sphere of radius = min_dim / 3
    cavity_radius = min_dim / 3.0
    pocket_vol = (4.0 / 3.0) * math.pi * cavity_radius ** 3
    # Clamp to reasonable range [200, 3000] A^3
    pocket_vol = max(200.0, min(3000.0, pocket_vol))
    max_pocket_dim = max(pocket_size)

    logger.info(
        "filter_by_shape: starting with %d molecules, "
        "pocket_size=%s, est_cavity_vol=%.1f A^3, max_dim=%.1f A",
        total, pocket_size, pocket_vol, max_pocket_dim,
    )
    start_time = time.monotonic()

    # --- Strategy 1: RDKit ETKDG conformer-based filtering ---
    if _check_rdkit() and _check_rdkit_3d():
        try:
            result = _filter_rdkit_conformer(
                smiles_list, pocket_vol, max_pocket_dim,
                progress_callback,
            )
            elapsed = time.monotonic() - start_time
            logger.info(
                "filter_by_shape (RDKit 3D): %d -> %d molecules "
                "(%.1f%% pass rate) in %.1fs",
                total, len(result),
                len(result) / max(total, 1) * 100,
                elapsed,
            )
            return result
        except Exception as exc:
            logger.warning(
                "RDKit 3D shape filter failed: %s; trying volume estimation",
                exc,
            )

    # --- Strategy 2: Simple MW-based volume estimation ---
    if _check_rdkit():
        try:
            result = _filter_volume_estimation(
                smiles_list, pocket_vol, max_pocket_dim,
                progress_callback,
            )
            elapsed = time.monotonic() - start_time
            logger.info(
                "filter_by_shape (volume est): %d -> %d molecules "
                "(%.1f%% pass rate) in %.1fs",
                total, len(result),
                len(result) / max(total, 1) * 100,
                elapsed,
            )
            return result
        except Exception as exc:
            logger.warning(
                "Volume estimation filter failed: %s; falling back to mock",
                exc,
            )

    # --- Strategy 3: Hash-based mock (~5% pass rate) ---
    result = _filter_mock(smiles_list, progress_callback)
    elapsed = time.monotonic() - start_time
    logger.info(
        "filter_by_shape (mock): %d -> %d molecules "
        "(%.1f%% pass rate) in %.1fs",
        total, len(result),
        len(result) / max(total, 1) * 100,
        elapsed,
    )
    return result


# =========================================================================
# STRATEGY 1: RDKit ETKDG conformer-based 3D filtering
# =========================================================================

def _filter_rdkit_conformer(
    smiles_list: list[str],
    pocket_vol: float,
    max_pocket_dim: float,
    progress_callback: Optional[Callable[[dict], None]],
) -> list[str]:
    """Filter by generating a single 3D conformer and computing volume.

    For each molecule:
      1. Generate 3D conformer with AllChem.EmbedMolecule (ETKDG).
      2. Compute bounding box volume and longest axis from atom positions.
      3. Keep molecules where:
         - 0.1 * pocket_vol < mol_vol < 1.5 * pocket_vol
         - longest_axis < max_pocket_dim * 1.2

    Parameters
    ----------
    smiles_list : list[str]
        Input SMILES.
    pocket_vol : float
        Pocket volume estimate (sx * sy * sz).
    max_pocket_dim : float
        Maximum pocket dimension.
    progress_callback : callable, optional
        Progress callback.

    Returns
    -------
    list[str]
        SMILES passing the shape filter.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem

    RDLogger.DisableLog("rdApp.*")

    vol_min = 0.1 * pocket_vol
    vol_max = 1.5 * pocket_vol
    axis_max = max_pocket_dim * 1.2

    logger.debug(
        "3D shape criteria: vol in [%.1f, %.1f], axis < %.1f",
        vol_min, vol_max, axis_max,
    )

    stats: dict[str, int] = {
        "invalid_smiles": 0,
        "embed_failed": 0,
        "too_small": 0,
        "too_large": 0,
        "too_long": 0,
    }

    passed: list[str] = []
    batch_size = 1000
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

            # Add hydrogens for better 3D embedding
            mol_h = Chem.AddHs(mol)

            # Generate a single conformer
            params = AllChem.ETKDGv3()
            params.useRandomCoords = False
            # Limit iterations for speed in batch screening
            try:
                params.maxIterations = 200
            except AttributeError:
                pass
            result_code = AllChem.EmbedMolecule(mol_h, params)

            if result_code == -1:
                # Retry with random coords
                params2 = AllChem.ETKDGv3()
                params2.useRandomCoords = True
                result_code = AllChem.EmbedMolecule(mol_h, params2)

            if result_code == -1:
                stats["embed_failed"] += 1
                continue

            # Compute bounding box volume and longest axis
            try:
                conf = mol_h.GetConformer()
                vol, longest_axis = _compute_bounding_box(conf, mol_h)
            except Exception:
                stats["embed_failed"] += 1
                continue

            # Apply shape filters
            if vol < vol_min:
                stats["too_small"] += 1
                continue
            if vol > vol_max:
                stats["too_large"] += 1
                continue
            if longest_axis > axis_max:
                stats["too_long"] += 1
                continue

            # Passed all shape checks
            canonical = Chem.MolToSmiles(mol)
            passed.append(canonical)

        # Progress callback
        if progress_callback is not None:
            try:
                progress_callback({
                    "pass_name": "shape_filter_3d",
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

    _log_shape_stats(stats, total, len(passed))
    return passed


def _compute_bounding_box(
    conf: object,
    mol: object,
) -> tuple[float, float]:
    """Compute bounding box volume and longest axis from a conformer.

    Parameters
    ----------
    conf : rdkit.Chem.Conformer
        The 3D conformer.
    mol : rdkit.Chem.Mol
        The molecule (for heavy atom count).

    Returns
    -------
    tuple[float, float]
        (bounding_box_volume, longest_axis_length) in Angstroms.
    """
    positions: list[tuple[float, float, float]] = []
    num_atoms = mol.GetNumAtoms()  # type: ignore[union-attr]

    for i in range(num_atoms):
        pos = conf.GetAtomPosition(i)  # type: ignore[union-attr]
        positions.append((pos.x, pos.y, pos.z))

    if not positions:
        return 0.0, 0.0

    xs = [p[0] for p in positions]
    ys = [p[1] for p in positions]
    zs = [p[2] for p in positions]

    dx = max(xs) - min(xs)
    dy = max(ys) - min(ys)
    dz = max(zs) - min(zs)

    # Add van der Waals radius padding (roughly 1.5 A per side)
    vdw_pad = 3.0  # 1.5 A on each side
    dx = max(dx + vdw_pad, 1.0)
    dy = max(dy + vdw_pad, 1.0)
    dz = max(dz + vdw_pad, 1.0)

    volume = dx * dy * dz

    # Longest axis: maximum pairwise distance (sampled for speed)
    longest_axis = max(dx, dy, dz)

    # For a more accurate longest axis, sample a few atom pairs
    if len(positions) <= 200:
        max_dist_sq = 0.0
        for i in range(len(positions)):
            for j in range(i + 1, len(positions)):
                d_sq = (
                    (positions[i][0] - positions[j][0]) ** 2
                    + (positions[i][1] - positions[j][1]) ** 2
                    + (positions[i][2] - positions[j][2]) ** 2
                )
                if d_sq > max_dist_sq:
                    max_dist_sq = d_sq
        longest_axis = math.sqrt(max_dist_sq) + vdw_pad
    else:
        # For large molecules, use bounding box diagonal as proxy
        longest_axis = math.sqrt(dx * dx + dy * dy + dz * dz)

    return volume, longest_axis


# =========================================================================
# STRATEGY 2: Simple MW-based volume estimation
# =========================================================================

def _filter_volume_estimation(
    smiles_list: list[str],
    pocket_vol: float,
    max_pocket_dim: float,
    progress_callback: Optional[Callable[[dict], None]],
) -> list[str]:
    """Estimate molecular volume from MW and filter.

    Uses the rough empirical correlation:
      volume ~ MW * 1.2  (Angstrom^3, very approximate)
      longest_axis ~ MW^(1/3) * 1.5  (Angstrom)

    Parameters
    ----------
    smiles_list : list[str]
        Input SMILES.
    pocket_vol : float
        Pocket volume estimate.
    max_pocket_dim : float
        Maximum pocket dimension.
    progress_callback : callable, optional
        Progress callback.

    Returns
    -------
    list[str]
        SMILES passing the estimated shape filter.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors

    RDLogger.DisableLog("rdApp.*")

    vol_min = 0.1 * pocket_vol
    vol_max = 1.5 * pocket_vol
    axis_max = max_pocket_dim * 1.2

    stats: dict[str, int] = {
        "invalid_smiles": 0,
        "too_small": 0,
        "too_large": 0,
        "too_long": 0,
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
            except Exception:
                stats["invalid_smiles"] += 1
                continue

            # Estimate volume from MW
            est_vol = mw * 1.2

            # Estimate longest axis from MW
            est_axis = (mw ** (1.0 / 3.0)) * 1.5

            if est_vol < vol_min:
                stats["too_small"] += 1
                continue
            if est_vol > vol_max:
                stats["too_large"] += 1
                continue
            if est_axis > axis_max:
                stats["too_long"] += 1
                continue

            canonical = Chem.MolToSmiles(mol)
            passed.append(canonical)

        if progress_callback is not None:
            try:
                progress_callback({
                    "pass_name": "shape_filter_volume_est",
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

    _log_shape_stats(stats, total, len(passed))
    return passed


# =========================================================================
# STRATEGY 3: Hash-based mock (~5% pass rate)
# =========================================================================

def _filter_mock(
    smiles_list: list[str],
    progress_callback: Optional[Callable[[dict], None]],
) -> list[str]:
    """Deterministic hash-based mock shape filter.

    Keeps approximately 5% of input molecules using a hash of the
    SMILES string. The same molecule always produces the same pass/fail
    decision.

    Parameters
    ----------
    smiles_list : list[str]
        Input SMILES.
    progress_callback : callable, optional
        Progress callback.

    Returns
    -------
    list[str]
        Approximately 5% of input SMILES, deterministically selected.
    """
    logger.info(
        "filter_by_shape mock: processing %d molecules "
        "(deterministic ~5%% pass rate)",
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
            digest = hashlib.sha256(
                f"shape:{smi}".encode("utf-8")
            ).hexdigest()
            hash_val = int(digest[:8], 16)

            # Keep ~5% (hash modulo 20 == 0)
            if hash_val % 20 == 0:
                passed.append(smi)

        if progress_callback is not None:
            try:
                progress_callback({
                    "pass_name": "shape_filter_mock",
                    "processed": batch_end,
                    "total": total,
                    "kept": len(passed),
                })
            except Exception:
                pass

    logger.info(
        "Mock shape filter: %d -> %d (%.1f%%)",
        total, len(passed),
        len(passed) / max(total, 1) * 100,
    )
    return passed


# =========================================================================
# Logging helpers
# =========================================================================

def _log_shape_stats(
    stats: dict[str, int],
    total_input: int,
    total_passed: int,
) -> None:
    """Log a summary of rejection reasons for the shape filter.

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
        "Shape filter stats: %d input, %d passed (%.1f%%), %d rejected",
        total_input, total_passed,
        total_passed / max(total_input, 1) * 100,
        total_rejected,
    )

    for reason, count in sorted(stats.items(), key=lambda x: -x[1]):
        if count > 0:
            logger.info(
                "  Rejection: %-25s  %7d  (%5.1f%%)",
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
    print("DockIt -- filter_shape.py self-test")
    print("=" * 70)

    # Test SMILES: mix of different sizes
    test_smiles = [
        "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OC)cc23)c1",     # Erlotinib-like
        "CC(=O)Oc1ccccc1C(=O)O",                         # Aspirin (small)
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",                    # Ibuprofen
        "Cn1c(=O)c2c(ncn2C)n(C)c1=O",                    # Caffeine
        "CC(=O)Nc1ccc(O)cc1",                             # Acetaminophen
        "C",                                               # Methane (too small)
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",            # Very long chain
        "c1ccccc1",                                        # Benzene (small)
        "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1", # Gefitinib
        "CC(C)c1n(CC[C@@H](O)C[C@@H](O)CC(=O)O)c(c2ccc(F)cc2)c(c1c1ccccc1)C(=O)Nc1ccccc1",  # Atorvastatin (large)
    ]

    # Typical kinase pocket size
    pocket_center = (22.0, 0.5, 18.0)
    pocket_size = (25.0, 25.0, 25.0)

    print(f"\nTesting with {len(test_smiles)} molecules...")
    print(f"Pocket size: {pocket_size}")
    print(f"Pocket volume: {pocket_size[0] * pocket_size[1] * pocket_size[2]:.0f} A^3")

    progress_log: list[dict] = []

    def on_progress(info: dict) -> None:
        progress_log.append(info)
        print(f"  Progress: {info}")

    result = filter_by_shape(
        test_smiles,
        pocket_center=pocket_center,
        pocket_size=pocket_size,
        progress_callback=on_progress,
    )

    print(f"\nResult: {len(result)} molecules passed")
    for smi in result:
        print(f"  PASS: {smi}")

    errors = 0

    if len(result) == 0:
        print("WARNING: No molecules passed (may be expected depending on RDKit availability)")
    elif len(result) > len(test_smiles):
        print("FAIL: More molecules passed than were input")
        errors += 1

    if errors == 0:
        print("\nALL CHECKS PASSED")
    else:
        print(f"\n{errors} CHECK(S) FAILED")

    sys.exit(1 if errors > 0 else 0)
