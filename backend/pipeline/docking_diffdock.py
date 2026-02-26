"""
DockIt pipeline -- AI-based molecular docking with DiffDock.

DiffDock uses a diffusion generative model for blind molecular docking,
often achieving higher accuracy than traditional scoring-function-based
methods like AutoDock Vina, especially on difficult or flexible targets.

This module provides an opt-in "precise mode" for DockIt V2.  It follows
the same pattern as the other pipeline modules: attempt the real tool
first, then fall back gracefully to a deterministic mock when DiffDock
(or a suitable GPU) is not available.

Real implementation:
    Calls the DiffDock inference CLI (``inference.py``) or its importable
    Python API, parses output poses and confidence scores, and converts
    the confidence metric to an estimated affinity scale comparable to
    Vina's kcal/mol output.

Mock implementation:
    Uses a hash-based deterministic RNG (same strategy as ``docking.py``)
    but produces tighter, stronger affinities (-8 to -14 kcal/mol) and a
    confidence score in [0, 1] derived from a pseudo-drug-likeness
    heuristic.  This makes results visually distinguishable as
    "AI docking" in the frontend.
"""

from __future__ import annotations

import hashlib
import logging
import math
import os
import random
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Callable, Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Availability check (cached after first call)
# ---------------------------------------------------------------------------

_DIFFDOCK_AVAILABLE: Optional[bool] = None
_DIFFDOCK_METHOD: Optional[str] = None  # "cli" or "python"


def check_diffdock_available() -> bool:
    """Check if DiffDock is installed and usable.

    The function probes two strategies in order:

    1. **CLI** -- look for the ``diffdock`` or ``DiffDock/inference.py``
       script on PATH or at a well-known location.
    2. **Python import** -- try to import the ``diffdock`` package
       (some installations expose a Python module).

    A GPU is *not* strictly required (DiffDock can run on CPU, albeit
    slowly), but this function logs a warning when CUDA is unavailable.

    Returns
    -------
    bool
        ``True`` if DiffDock can be invoked, ``False`` otherwise.
    """
    global _DIFFDOCK_AVAILABLE, _DIFFDOCK_METHOD

    if _DIFFDOCK_AVAILABLE is not None:
        return _DIFFDOCK_AVAILABLE

    # Strategy 1: CLI binary / script
    diffdock_bin = shutil.which("diffdock")
    diffdock_script = os.environ.get("DIFFDOCK_INFERENCE_SCRIPT")

    if diffdock_bin is not None:
        logger.info("DiffDock CLI found at %s", diffdock_bin)
        _DIFFDOCK_AVAILABLE = True
        _DIFFDOCK_METHOD = "cli"
    elif diffdock_script and Path(diffdock_script).is_file():
        logger.info("DiffDock inference script found at %s", diffdock_script)
        _DIFFDOCK_AVAILABLE = True
        _DIFFDOCK_METHOD = "cli"
    else:
        # Strategy 2: Python import
        try:
            import diffdock  # type: ignore[import-untyped]
            logger.info("DiffDock Python package importable")
            _DIFFDOCK_AVAILABLE = True
            _DIFFDOCK_METHOD = "python"
        except ImportError:
            _DIFFDOCK_AVAILABLE = False
            logger.warning(
                "DiffDock NOT found (no CLI binary, no inference script, "
                "no Python package); AI docking will be mocked"
            )

    # Check GPU availability (informational)
    if _DIFFDOCK_AVAILABLE:
        try:
            import torch
            if torch.cuda.is_available():
                logger.info("CUDA GPU available for DiffDock")
            else:
                logger.warning(
                    "No CUDA GPU detected; DiffDock will run on CPU (slow)"
                )
        except ImportError:
            logger.warning("PyTorch not importable; GPU check skipped")

    return _DIFFDOCK_AVAILABLE


# ---------------------------------------------------------------------------
# Single-ligand docking
# ---------------------------------------------------------------------------

def dock_diffdock(
    protein_pdb: Path,
    ligand_smiles: str,
    ligand_name: str,
    work_dir: Path,
    n_poses: int = 5,
    timeout_sec: int = 120,
) -> dict:
    """Dock a single ligand against a protein using DiffDock.

    Parameters
    ----------
    protein_pdb : Path
        Path to the protein PDB file.
    ligand_smiles : str
        SMILES string of the ligand to dock.
    ligand_name : str
        Human-readable ligand identifier (used for logging and file naming).
    work_dir : Path
        Working directory for intermediate and output files.
    n_poses : int
        Number of poses to generate (default 5).  The best pose (highest
        confidence) is selected as the primary result.
    timeout_sec : int
        Maximum wall-clock time in seconds before aborting (default 120).

    Returns
    -------
    dict
        Keys:

        - ``affinity`` (float): estimated binding score in kcal/mol
          (negative = better).
        - ``confidence`` (float): DiffDock confidence score in [0, 1].
        - ``pose_pdb_path`` (Path or None): path to the best-pose PDB
          file, or ``None`` if pose generation failed.
        - ``method`` (str): ``"diffdock"`` for real runs, or
          ``"diffdock_mock"`` when falling back.
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    if check_diffdock_available():
        result = _dock_diffdock_real(
            protein_pdb, ligand_smiles, ligand_name, work_dir,
            n_poses, timeout_sec,
        )
        if result is not None:
            return result
        logger.warning(
            "Real DiffDock run failed for %s; falling back to mock",
            ligand_name,
        )

    return _dock_diffdock_mock(
        protein_pdb, ligand_smiles, ligand_name, work_dir,
    )


# ---------------------------------------------------------------------------
# Real DiffDock execution
# ---------------------------------------------------------------------------

def _dock_diffdock_real(
    protein_pdb: Path,
    ligand_smiles: str,
    ligand_name: str,
    work_dir: Path,
    n_poses: int,
    timeout_sec: int,
) -> Optional[dict]:
    """Invoke the actual DiffDock tool and parse results.

    Supports two modes:

    * **CLI mode** -- writes a CSV input file, calls the inference script,
      and parses the output directory for ranked pose SDF files and a
      ``confidence`` field.
    * **Python mode** -- imports the diffdock module and calls its API
      directly (less common but supported by some forks).
    """
    safe_name = _sanitize_name(ligand_name)
    run_dir = work_dir / f"diffdock_{safe_name}"
    run_dir.mkdir(parents=True, exist_ok=True)

    if _DIFFDOCK_METHOD == "cli":
        return _run_diffdock_cli(
            protein_pdb, ligand_smiles, ligand_name,
            run_dir, n_poses, timeout_sec,
        )

    if _DIFFDOCK_METHOD == "python":
        return _run_diffdock_python(
            protein_pdb, ligand_smiles, ligand_name,
            run_dir, n_poses, timeout_sec,
        )

    return None


def _run_diffdock_cli(
    protein_pdb: Path,
    ligand_smiles: str,
    ligand_name: str,
    run_dir: Path,
    n_poses: int,
    timeout_sec: int,
) -> Optional[dict]:
    """Execute DiffDock via its CLI inference script."""
    # Prepare the input CSV that DiffDock expects:
    #   protein_path, ligand_description, complex_name
    csv_path = run_dir / "input.csv"
    csv_path.write_text(
        "complex_name,protein_path,ligand_description,protein_sequence\n"
        f"{ligand_name},{protein_pdb},{ligand_smiles},\n"
    )

    out_dir = run_dir / "output"
    out_dir.mkdir(exist_ok=True)

    # Build command
    diffdock_bin = shutil.which("diffdock")
    diffdock_script = os.environ.get("DIFFDOCK_INFERENCE_SCRIPT")

    if diffdock_bin:
        cmd = [diffdock_bin]
    elif diffdock_script:
        cmd = ["python", diffdock_script]
    else:
        return None

    cmd += [
        "--protein_ligand_csv", str(csv_path),
        "--out_dir", str(out_dir),
        "--samples_per_complex", str(n_poses),
        "--no_final_step_noise",
    ]

    try:
        logger.debug("DiffDock CLI command: %s", " ".join(cmd))
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout_sec,
            cwd=str(run_dir),
        )

        if proc.returncode != 0:
            logger.warning(
                "DiffDock CLI exited %d: %s",
                proc.returncode,
                proc.stderr[:500] if proc.stderr else "(no stderr)",
            )
            return None

        return _parse_diffdock_output(out_dir, ligand_name)

    except subprocess.TimeoutExpired:
        logger.warning(
            "DiffDock timed out after %ds for %s", timeout_sec, ligand_name,
        )
        return None
    except Exception as exc:
        logger.warning("DiffDock CLI execution error: %s", exc)
        return None


def _run_diffdock_python(
    protein_pdb: Path,
    ligand_smiles: str,
    ligand_name: str,
    run_dir: Path,
    n_poses: int,
    timeout_sec: int,
) -> Optional[dict]:
    """Execute DiffDock via its Python API (if available)."""
    try:
        import diffdock  # type: ignore[import-untyped]

        # The exact API varies by fork; try the most common pattern.
        result = diffdock.dock(
            protein_path=str(protein_pdb),
            ligand=ligand_smiles,
            out_dir=str(run_dir),
            samples_per_complex=n_poses,
        )

        if result is None:
            return None

        # Expect result to be a dict or list of dicts with confidence/pose
        if isinstance(result, list):
            result = max(result, key=lambda r: r.get("confidence", 0.0))

        confidence = float(result.get("confidence", 0.5))
        pose_path_str = result.get("pose_path") or result.get("sdf_path")
        pose_pdb_path: Optional[Path] = None

        if pose_path_str:
            pose_file = Path(pose_path_str)
            if pose_file.exists():
                pose_pdb_path = _convert_sdf_to_pdb(pose_file, run_dir)

        affinity = _confidence_to_affinity(confidence)

        return {
            "affinity": affinity,
            "confidence": round(confidence, 4),
            "pose_pdb_path": pose_pdb_path,
            "method": "diffdock",
        }

    except ImportError:
        logger.warning("diffdock Python import failed at runtime")
        return None
    except Exception as exc:
        logger.warning("DiffDock Python API error: %s", exc)
        return None


def _parse_diffdock_output(
    out_dir: Path,
    ligand_name: str,
) -> Optional[dict]:
    """Parse DiffDock output directory for the best pose and confidence.

    DiffDock typically outputs ranked SDF files named like:
        ``rank1_confidence-1.23.sdf``
    where the number after "confidence" is the raw confidence score
    (higher is better, can be negative in log-space).
    """
    import re

    complex_dirs = list(out_dir.iterdir()) if out_dir.exists() else []

    # DiffDock nests output under a subdirectory named after the complex
    search_dirs = complex_dirs if complex_dirs else [out_dir]

    best_confidence: float = -999.0
    best_pose: Optional[Path] = None

    for d in search_dirs:
        search_path = d if d.is_dir() else out_dir
        for sdf_file in sorted(search_path.glob("rank*_confidence*.sdf")):
            # Extract confidence from filename
            match = re.search(r"confidence([+-]?[\d.]+)", sdf_file.name)
            if match:
                try:
                    conf = float(match.group(1))
                    if conf > best_confidence:
                        best_confidence = conf
                        best_pose = sdf_file
                except ValueError:
                    pass

    if best_pose is None:
        # Fallback: look for any SDF file
        for d in search_dirs:
            search_path = d if d.is_dir() else out_dir
            sdf_files = sorted(search_path.glob("*.sdf"))
            if sdf_files:
                best_pose = sdf_files[0]
                best_confidence = 0.0
                break

    if best_pose is None:
        logger.warning("No DiffDock output poses found in %s", out_dir)
        return None

    # Normalise confidence to [0, 1] range
    # DiffDock raw confidence is typically in [-10, +2]; we use a sigmoid
    normalised_confidence = _sigmoid(best_confidence)

    # Convert to affinity estimate
    affinity = _confidence_to_affinity(normalised_confidence)

    # Convert best SDF pose to PDB for the 3D viewer
    pose_pdb_path = _convert_sdf_to_pdb(best_pose, best_pose.parent)

    return {
        "affinity": affinity,
        "confidence": round(normalised_confidence, 4),
        "pose_pdb_path": pose_pdb_path,
        "method": "diffdock",
    }


# ---------------------------------------------------------------------------
# Mock DiffDock (when not installed)
# ---------------------------------------------------------------------------

def _dock_diffdock_mock(
    protein_pdb: Path,
    ligand_smiles: str,
    ligand_name: str,
    work_dir: Path,
) -> dict:
    """Generate deterministic mock DiffDock results.

    Uses the same hash-seeded RNG approach as ``docking._run_vina_mock``
    but with a shifted distribution to simulate DiffDock's typically
    stronger and more consistent binding predictions:

    - Affinity range: -8.0 to -14.0 kcal/mol (vs Vina's -4.0 to -12.0)
    - Lower variance (tighter distribution around -11 kcal/mol)
    - Confidence score derived from SMILES-based pseudo-drug-likeness
    - Mock pose PDB files generated with offset coordinates

    Parameters
    ----------
    protein_pdb : Path
        Protein PDB (used for coordinate reference in mock poses).
    ligand_smiles : str
        Ligand SMILES string.
    ligand_name : str
        Ligand identifier.
    work_dir : Path
        Output directory.

    Returns
    -------
    dict
        Same schema as ``dock_diffdock`` return value, with
        ``method="diffdock_mock"``.
    """
    # Deterministic seed from ligand name + a DiffDock-specific salt
    seed_input = f"diffdock_v2_{ligand_name}_{ligand_smiles}"
    seed = int(hashlib.sha256(seed_input.encode()).hexdigest()[:8], 16)
    rng = random.Random(seed)

    # Compute pseudo-drug-likeness from SMILES to influence confidence
    drug_likeness = _estimate_drug_likeness(ligand_smiles)

    # Confidence: base from drug-likeness, add small random perturbation
    confidence = max(0.05, min(0.99, drug_likeness + rng.gauss(0.0, 0.08)))
    confidence = round(confidence, 4)

    # Affinity: DiffDock mock produces stronger affinities than Vina mock.
    # Mean around -11.0 kcal/mol with std ~1.5 (Vina mock: uniform -4 to -12).
    # Confidence also influences affinity (higher confidence -> stronger binding).
    base_affinity = rng.gauss(-11.0, 1.5)
    confidence_bonus = (confidence - 0.5) * 2.0  # [-1, +1] range
    affinity = base_affinity + confidence_bonus
    affinity = round(max(-14.0, min(-8.0, affinity)), 1)

    # Generate a mock pose PDB
    safe_name = _sanitize_name(ligand_name)
    pose_path = work_dir / f"{safe_name}_diffdock_pose.pdb"
    _generate_mock_pose_pdb(pose_path, ligand_name, ligand_smiles, affinity, confidence)

    logger.debug(
        "Mock DiffDock for %s: affinity=%.1f confidence=%.4f",
        ligand_name, affinity, confidence,
    )

    return {
        "affinity": affinity,
        "confidence": confidence,
        "pose_pdb_path": pose_path if pose_path.exists() else None,
        "method": "diffdock_mock",
    }


def _estimate_drug_likeness(smiles: str) -> float:
    """Estimate a rough drug-likeness score from a SMILES string.

    This is a lightweight heuristic used for mock confidence scoring.
    It does NOT replace proper QED computation (which requires RDKit).

    Criteria (Lipinski-inspired):

    - Molecular weight proxy: SMILES length (shorter = smaller = better)
    - Heteroatom ratio: N and O count vs total heavy atoms
    - Aromatic character: ring indicators present
    - Charged groups penalty: +/- characters

    Returns
    -------
    float
        Score in approximately [0.2, 0.9].
    """
    if not smiles:
        return 0.5

    length = len(smiles)
    n_count = smiles.upper().count("N")
    o_count = smiles.upper().count("O")
    ring_count = smiles.count("1") + smiles.count("2") + smiles.count("3")
    charged = smiles.count("+") + smiles.count("-")

    # Length penalty: optimal around 30-60 chars
    if 20 <= length <= 80:
        length_score = 0.8
    elif length < 20:
        length_score = 0.6  # too small, likely fragment
    else:
        length_score = max(0.3, 0.8 - (length - 80) * 0.005)

    # Heteroatom score: want some N/O but not too many
    hetero_count = n_count + o_count
    if 2 <= hetero_count <= 8:
        hetero_score = 0.8
    elif hetero_count < 2:
        hetero_score = 0.5
    else:
        hetero_score = max(0.3, 0.8 - (hetero_count - 8) * 0.05)

    # Aromatic/ring bonus
    ring_score = min(0.9, 0.5 + ring_count * 0.1)

    # Charge penalty
    charge_penalty = charged * 0.1

    score = (0.4 * length_score + 0.3 * hetero_score + 0.3 * ring_score
             - charge_penalty)
    return max(0.2, min(0.9, score))


def _generate_mock_pose_pdb(
    pose_path: Path,
    ligand_name: str,
    ligand_smiles: str,
    affinity: float,
    confidence: float,
) -> None:
    """Write a minimal mock PDB file representing a docked pose.

    The mock places pseudo-atoms in a small cluster to give the 3D
    viewer something to render.  Coordinates are derived deterministically
    from the ligand name so repeated runs are reproducible.

    Parameters
    ----------
    pose_path : Path
        Output PDB file path.
    ligand_name : str
        Used for REMARK headers and coordinate seeding.
    ligand_smiles : str
        Used to determine approximate atom count.
    affinity : float
        Written into REMARK header.
    confidence : float
        Written into REMARK header.
    """
    try:
        # Seed from name for reproducible coordinates
        seed = int(hashlib.md5(ligand_name.encode()).hexdigest()[:8], 16)
        rng = random.Random(seed)

        # Estimate heavy atom count from SMILES
        heavy_atoms = sum(
            1 for ch in ligand_smiles
            if ch.isalpha() and ch.upper() not in ("H",)
        )
        heavy_atoms = max(5, min(heavy_atoms, 50))

        lines: list[str] = []
        lines.append(
            f"REMARK   DiffDock AI Docking Result (mock)\n"
        )
        lines.append(
            f"REMARK   Ligand: {ligand_name}\n"
        )
        lines.append(
            f"REMARK   Affinity: {affinity:.1f} kcal/mol\n"
        )
        lines.append(
            f"REMARK   Confidence: {confidence:.4f}\n"
        )
        lines.append(
            f"REMARK   Method: diffdock_mock\n"
        )

        # Place atoms in a compact cluster around a random center
        cx = rng.uniform(10.0, 30.0)
        cy = rng.uniform(10.0, 30.0)
        cz = rng.uniform(10.0, 30.0)

        elements = ["C", "C", "C", "N", "O", "C", "C", "C", "N", "C"]

        for i in range(heavy_atoms):
            x = cx + rng.gauss(0, 2.0)
            y = cy + rng.gauss(0, 2.0)
            z = cz + rng.gauss(0, 2.0)
            elem = elements[i % len(elements)]
            atom_name = f"{elem}{i + 1}"[:4].ljust(4)
            line = (
                f"HETATM{i + 1:5d} {atom_name} LIG A   1    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
                f"  1.00  0.00           {elem:>2s}\n"
            )
            lines.append(line)

        lines.append("END\n")
        pose_path.write_text("".join(lines))

    except Exception as exc:
        logger.warning("Failed to create mock DiffDock pose file: %s", exc)


# ---------------------------------------------------------------------------
# Batch docking
# ---------------------------------------------------------------------------

def dock_all_diffdock(
    protein_pdb: Path,
    ligands: list[dict],
    work_dir: Path,
    n_poses: int = 5,
    progress_callback: Optional[Callable[[int, str], None]] = None,
) -> list[dict]:
    """Dock multiple ligands using DiffDock.

    Each ligand dict must have ``name`` and ``smiles`` keys.  Additional
    keys (e.g. ``source``, ``chembl_id``) are preserved in the output.

    Parameters
    ----------
    protein_pdb : Path
        Path to the protein PDB file.
    ligands : list[dict]
        List of ligand descriptors.  Required keys: ``name``, ``smiles``.
    work_dir : Path
        Working directory for intermediate files.
    n_poses : int
        Number of poses per ligand (default 5).
    progress_callback : callable or None
        Called with ``(percent: int, message: str)`` after each ligand.
        Percent is mapped to the 40%--90% range (matching the overall
        pipeline progress convention from ``tasks.py``).

    Returns
    -------
    list[dict]
        Each entry contains all original ligand keys plus:

        - ``affinity`` (float)
        - ``confidence`` (float)
        - ``pose_pdb_path`` (str or None)
        - ``method`` (str)

        Sorted by affinity ascending (most negative = strongest binding
        first).
    """
    ligand_dir = work_dir / "diffdock_ligands"
    ligand_dir.mkdir(parents=True, exist_ok=True)

    results: list[dict] = []
    total = len(ligands)

    for i, lig in enumerate(ligands):
        name = lig.get("name", f"ligand_{i}")
        smiles = lig.get("smiles", "")

        if not smiles:
            logger.warning("Skipping ligand %s: no SMILES", name)
            continue

        # Report progress (40% -> 90% range)
        if progress_callback is not None:
            pct = int(40 + (i / max(total, 1)) * 50)
            progress_callback(pct, f"AI docking {name} ({i + 1}/{total})")

        # Dock single ligand
        try:
            dock_result = dock_diffdock(
                protein_pdb=protein_pdb,
                ligand_smiles=smiles,
                ligand_name=name,
                work_dir=ligand_dir,
                n_poses=n_poses,
            )
        except Exception as exc:
            logger.warning("DiffDock failed for %s: %s", name, exc)
            continue

        results.append({
            **lig,  # preserve all original keys
            "affinity": dock_result["affinity"],
            "confidence": dock_result["confidence"],
            "pose_pdb_path": str(dock_result["pose_pdb_path"])
            if dock_result.get("pose_pdb_path") else None,
            "method": dock_result["method"],
        })

    # Sort by affinity ascending (most negative = strongest binding = best)
    results.sort(key=lambda r: r.get("affinity", 0.0))

    logger.info(
        "DiffDock docking complete: %d/%d ligands succeeded",
        len(results), total,
    )
    return results


# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

def _sanitize_name(name: str) -> str:
    """Convert a ligand name to a filesystem-safe string."""
    return "".join(
        c if c.isalnum() or c in "-_" else "_" for c in name
    )[:60]


def _sigmoid(x: float) -> float:
    """Standard sigmoid function, clamped to avoid overflow."""
    x_clamped = max(-20.0, min(20.0, x))
    return 1.0 / (1.0 + math.exp(-x_clamped))


def _confidence_to_affinity(confidence: float) -> float:
    """Convert a normalised DiffDock confidence [0, 1] to an estimated
    binding affinity in kcal/mol.

    The mapping is linear:

    - confidence 1.0  -->  -14.0 kcal/mol  (very strong binding)
    - confidence 0.0  -->   -6.0 kcal/mol  (weak binding)

    This produces values in a range comparable to Vina output.

    Parameters
    ----------
    confidence : float
        Normalised confidence score in [0, 1].

    Returns
    -------
    float
        Estimated affinity in kcal/mol (negative = better).
    """
    # Linear mapping: affinity = -6.0 - 8.0 * confidence
    affinity = -6.0 - 8.0 * confidence
    return round(affinity, 1)


def _convert_sdf_to_pdb(sdf_path: Path, out_dir: Path) -> Optional[Path]:
    """Convert an SDF file to PDB format for the 3D viewer.

    Tries RDKit first, then Open Babel, then returns None.

    Parameters
    ----------
    sdf_path : Path
        Input SDF file.
    out_dir : Path
        Output directory for the PDB file.

    Returns
    -------
    Path or None
        Path to the converted PDB file, or ``None`` if conversion fails.
    """
    pdb_path = out_dir / f"{sdf_path.stem}.pdb"

    # Strategy 1: RDKit
    try:
        from rdkit import Chem

        supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        mol = next(iter(supplier), None)
        if mol is not None:
            Chem.MolToPDBFile(mol, str(pdb_path))
            if pdb_path.exists() and pdb_path.stat().st_size > 10:
                logger.debug("RDKit SDF->PDB conversion: %s", pdb_path)
                return pdb_path
    except ImportError:
        pass
    except Exception as exc:
        logger.debug("RDKit SDF->PDB failed: %s", exc)

    # Strategy 2: Open Babel
    if shutil.which("obabel"):
        try:
            result = subprocess.run(
                ["obabel", str(sdf_path), "-O", str(pdb_path)],
                capture_output=True, text=True, timeout=30,
            )
            if (result.returncode == 0 and pdb_path.exists()
                    and pdb_path.stat().st_size > 10):
                logger.debug("obabel SDF->PDB conversion: %s", pdb_path)
                return pdb_path
        except Exception as exc:
            logger.debug("obabel SDF->PDB failed: %s", exc)

    logger.warning("Could not convert SDF to PDB: %s", sdf_path)
    return None


# ---------------------------------------------------------------------------
# Standalone test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    """Quick self-test: dock a few molecules and print results.

    Run with:
        cd /home/antho/dockit/backend
        python -m pipeline.docking_diffdock
    """
    import sys

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s %(levelname)-8s %(name)s: %(message)s",
    )

    print("=" * 70)
    print("DockIt -- DiffDock module self-test")
    print("=" * 70)

    # Check availability
    available = check_diffdock_available()
    print(f"\nDiffDock available: {available}")
    print(f"Method: {_DIFFDOCK_METHOD or 'mock'}")

    # Create a temporary working directory
    test_dir = Path(tempfile.mkdtemp(prefix="dockit_diffdock_test_"))
    print(f"Working directory: {test_dir}")

    # Create a minimal mock protein PDB
    mock_pdb = test_dir / "protein.pdb"
    mock_pdb.write_text(
        "HEADER    TEST PROTEIN\n"
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2  CA  GLY A   2       3.800   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  CA  LEU A   3       7.600   0.000   0.000  1.00  0.00           C\n"
        "END\n"
    )

    # Test ligands (including Erlotinib for the EGFR validation)
    test_ligands = [
        {
            "name": "Erlotinib",
            "smiles": "C=Cc1cc2c(cc1OC)nc(nc2/N=C/c3ccccc3)Nc4ccc(cc4)F",
            "source": "test",
        },
        {
            "name": "Aspirin",
            "smiles": "CC(=O)Oc1ccccc1C(=O)O",
            "source": "test",
        },
        {
            "name": "Ibuprofen",
            "smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
            "source": "test",
        },
        {
            "name": "Caffeine",
            "smiles": "Cn1c(=O)c2c(ncn2C)n(c1=O)C",
            "source": "test",
        },
        {
            "name": "Gefitinib",
            "smiles": "COc1cc2ncnc(c2cc1OCCCN3CCOCC3)Nc4ccc(cc4F)Cl",
            "source": "test",
        },
    ]

    # Test single docking
    print("\n--- Single ligand docking ---")
    single_result = dock_diffdock(
        protein_pdb=mock_pdb,
        ligand_smiles=test_ligands[0]["smiles"],
        ligand_name=test_ligands[0]["name"],
        work_dir=test_dir,
    )
    print(f"  Ligand:     {test_ligands[0]['name']}")
    print(f"  Affinity:   {single_result['affinity']:.1f} kcal/mol")
    print(f"  Confidence: {single_result['confidence']:.4f}")
    print(f"  Method:     {single_result['method']}")
    print(f"  Pose file:  {single_result['pose_pdb_path']}")

    # Test batch docking
    print("\n--- Batch docking ---")

    def test_progress(pct: int, msg: str) -> None:
        print(f"  [{pct:3d}%] {msg}")

    batch_results = dock_all_diffdock(
        protein_pdb=mock_pdb,
        ligands=test_ligands,
        work_dir=test_dir,
        progress_callback=test_progress,
    )

    print(f"\n--- Results (sorted by affinity) ---")
    print(f"{'Rank':<5} {'Name':<15} {'Affinity':>10} {'Confidence':>12} {'Method':<16}")
    print("-" * 60)
    for rank, r in enumerate(batch_results, 1):
        print(
            f"{rank:<5} {r['name']:<15} {r['affinity']:>10.1f} "
            f"{r['confidence']:>12.4f} {r['method']:<16}"
        )

    # Validate: affinities should be in -8 to -14 range for mock
    affinities = [r["affinity"] for r in batch_results]
    all_in_range = all(-14.0 <= a <= -8.0 for a in affinities)
    print(f"\nAll affinities in [-14, -8] range: {all_in_range}")

    # Validate: confidence scores should be in [0, 1]
    confidences = [r["confidence"] for r in batch_results]
    all_valid_conf = all(0.0 <= c <= 1.0 for c in confidences)
    print(f"All confidences in [0, 1] range:  {all_valid_conf}")

    # Validate: all results have required keys
    required_keys = {"affinity", "confidence", "pose_pdb_path", "method", "name", "smiles"}
    all_have_keys = all(required_keys.issubset(r.keys()) for r in batch_results)
    print(f"All results have required keys:   {all_have_keys}")

    # Validate: deterministic (run again and compare)
    batch_results_2 = dock_all_diffdock(
        protein_pdb=mock_pdb,
        ligands=test_ligands,
        work_dir=test_dir,
    )
    deterministic = all(
        r1["affinity"] == r2["affinity"] and r1["confidence"] == r2["confidence"]
        for r1, r2 in zip(batch_results, batch_results_2)
    )
    print(f"Results are deterministic:        {deterministic}")

    # Overall pass/fail
    all_pass = all_in_range and all_valid_conf and all_have_keys and deterministic
    print(f"\n{'PASS' if all_pass else 'FAIL'}: DiffDock module self-test")

    # Cleanup
    import shutil as _shutil
    _shutil.rmtree(test_dir, ignore_errors=True)

    sys.exit(0 if all_pass else 1)
