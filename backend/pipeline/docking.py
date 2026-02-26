"""
DockIt pipeline -- Molecular docking with GNINA / AutoDock Vina.

Wraps GNINA (CNN-scored docking, preferred) and ``vina`` CLI, with a
mock fallback for environments where neither is installed.

V5bis: GNINA support (CNN scoring), consensus ranking, 3-score system.
"""

from __future__ import annotations

import hashlib
import logging
import random
import re
import shutil
import subprocess
from pathlib import Path
from typing import Callable, Optional

from pipeline.prepare import prepare_ligand

logger = logging.getLogger(__name__)

# Whether tools are available (checked once)
_GNINA_AVAILABLE: Optional[bool] = None
_VINA_AVAILABLE: Optional[bool] = None


def _check_gnina() -> bool:
    """Check if the gnina binary is available on PATH or at known location."""
    global _GNINA_AVAILABLE
    if _GNINA_AVAILABLE is None:
        gnina_bin = shutil.which("gnina") or "/usr/local/bin/gnina"
        _GNINA_AVAILABLE = Path(gnina_bin).exists()
        if _GNINA_AVAILABLE:
            logger.info("GNINA found: %s", gnina_bin)
        else:
            logger.info("GNINA not found; will try Vina or mock")
    return _GNINA_AVAILABLE


def _check_vina() -> bool:
    """Check if the vina binary is available on PATH."""
    global _VINA_AVAILABLE
    if _VINA_AVAILABLE is None:
        _VINA_AVAILABLE = shutil.which("vina") is not None
        if _VINA_AVAILABLE:
            logger.info("AutoDock Vina found on PATH")
        else:
            logger.info("AutoDock Vina NOT found")
    return _VINA_AVAILABLE


# ---------------------------------------------------------------------------
# Initial ligand placement (standard docking preprocessing)
# ---------------------------------------------------------------------------

def _place_ligand_at_box_center(ligand_path: str, center: list[float]) -> str:
    """Place initial ligand conformer at the search box center.

    Standard docking preprocessing: GNINA/Vina require the initial ligand
    pose to be within the search space for the stochastic search to converge.
    RDKit-generated conformers are centered at the origin (0,0,0), which may
    be far from the binding pocket.

    This translates the ligand centroid to the box center as a starting
    point for docking. The OUTPUT coordinates from the docking engine are
    the final binding pose and are never modified.

    Parameters
    ----------
    ligand_path : str
        Path to the input ligand (SDF or PDBQT).
    center : list[float]
        [x, y, z] centre of the docking search box.

    Returns
    -------
    str
        Path to the placed ligand file (same format as input).
    """
    try:
        path = Path(ligand_path)

        if path.suffix.lower() in ('.sdf', '.mol'):
            return _place_sdf_at_center(ligand_path, center)
        elif path.suffix.lower() == '.pdbqt':
            return _place_pdbqt_at_center(ligand_path, center)
        else:
            return ligand_path

    except Exception as exc:
        logger.debug("Ligand placement failed for %s: %s", ligand_path, exc)
        return ligand_path


def _place_sdf_at_center(sdf_path: str, center: list[float]) -> str:
    """Place SDF ligand at box center using RDKit."""
    from rdkit import Chem
    from rdkit.Geometry import Point3D
    import numpy as np

    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
    mol = next(iter(suppl), None)
    if mol is None or mol.GetNumConformers() == 0:
        return sdf_path

    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    coords = np.array([[conf.GetAtomPosition(i).x,
                         conf.GetAtomPosition(i).y,
                         conf.GetAtomPosition(i).z] for i in range(n)])

    centroid = coords.mean(axis=0)
    shift = np.array(center) - centroid

    if np.linalg.norm(shift) < 3.0:
        return sdf_path  # Already at center

    for i in range(n):
        new = coords[i] + shift
        conf.SetAtomPosition(i, Point3D(float(new[0]), float(new[1]), float(new[2])))

    out_path = str(Path(sdf_path).with_suffix(".placed.sdf"))
    writer = Chem.SDWriter(out_path)
    writer.write(mol)
    writer.close()
    return out_path


def _place_pdbqt_at_center(pdbqt_path: str, center: list[float]) -> str:
    """Place PDBQT ligand at box center by shifting ATOM coordinates."""
    import numpy as np

    path = Path(pdbqt_path)
    lines = path.read_text().split("\n")

    coords, atom_indices = [], []
    for i, line in enumerate(lines):
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
            try:
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                atom_indices.append(i)
            except ValueError:
                continue

    if not coords:
        return pdbqt_path

    coords = np.array(coords)
    centroid = coords.mean(axis=0)
    shift = np.array(center) - centroid

    if np.linalg.norm(shift) < 3.0:
        return pdbqt_path

    new_lines = list(lines)
    for idx, line_i in enumerate(atom_indices):
        line = lines[line_i]
        new_c = coords[idx] + shift
        new_lines[line_i] = (
            line[:30] + f"{new_c[0]:8.3f}{new_c[1]:8.3f}{new_c[2]:8.3f}" + line[54:]
        )

    out_path = str(path.with_suffix(".placed.pdbqt"))
    Path(out_path).write_text("\n".join(new_lines))
    return out_path


# ---------------------------------------------------------------------------
# GNINA docking (V5bis)
# ---------------------------------------------------------------------------

def dock_gnina(
    ligand_sdf: str,
    receptor_pdb: str,
    center: list[float],
    size: list[float] | None = None,
    exhaustiveness: int = 8,
) -> Optional[list[dict]]:
    """Dock using GNINA (CNN-scored docking).

    GNINA extends AutoDock Vina with convolutional neural network scoring
    functions that provide pose confidence and CNN-predicted affinities
    in addition to classical Vina scores.

    Parameters
    ----------
    ligand_sdf : str
        Path to the ligand SDF file.
    receptor_pdb : str
        Path to the receptor PDB file.
    center : list[float]
        [x, y, z] centre of the search box.
    size : list[float] or None
        [sx, sy, sz] dimensions of the search box (default [22, 22, 22]).
    exhaustiveness : int
        Search exhaustiveness parameter.

    Returns
    -------
    list[dict] or None
        List of poses, each with: ``vina_score``, ``cnn_score``, ``cnn_affinity``,
        ``pose_path``. Returns None if GNINA is not available or fails.
    """
    if size is None:
        size = [22.0, 22.0, 22.0]

    gnina_bin = shutil.which("gnina") or "/usr/local/bin/gnina"

    if not Path(gnina_bin).exists():
        logger.debug("GNINA not found at %s, returning None", gnina_bin)
        return None

    output_sdf = str(Path(ligand_sdf).with_suffix(".docked.sdf"))

    # GNINA MCMC sampling explores the full search box — no ligand
    # pre-centering needed (CACHE Challenge #1; Sunseri & Koes 2021).
    # Passing the original 3D conformer directly.

    # Prefer PDB receptor over PDBQT: GNINA handles PDB natively via
    # OpenBabel (CACHE Challenge recommendation). If receptor_pdb is a
    # .pdbqt, try to find the .pdb sibling.
    receptor_file = receptor_pdb
    if receptor_pdb.endswith(".pdbqt"):
        pdb_sibling = receptor_pdb.rsplit(".pdbqt", 1)[0] + ".pdb"
        clean_pdb = receptor_pdb.rsplit(".pdbqt", 1)[0].replace("_receptor", "_clean") + ".pdb"
        for candidate in [pdb_sibling, clean_pdb]:
            if Path(candidate).exists():
                receptor_file = candidate
                break

    cmd = [
        gnina_bin,
        "--receptor", receptor_file,
        "--ligand", ligand_sdf,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--exhaustiveness", str(exhaustiveness),
        "--cnn_scoring", "rescore",
        # Default CNN ensemble (5 models) per McNutt et al. 2021 —
        # outperforms any single model on DUD-E and LIT-PCBA.
        "--num_modes", "9",
        "--seed", "0",
        "--min_rmsd_filter", "1.0",
        "--out", output_sdf,
    ]

    try:
        logger.debug("GNINA command: %s", " ".join(cmd))
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
        )

        if result.returncode != 0:
            logger.warning(
                "GNINA exited %d: %s", result.returncode, result.stderr[:500]
            )
            return None

        # Parse GNINA output: scores are in stdout or SDF properties
        poses = _parse_gnina_output(result.stdout, output_sdf)
        if poses:
            logger.info("GNINA produced %d poses for %s", len(poses), ligand_sdf)
            return poses

        logger.warning("Could not parse GNINA output for %s", ligand_sdf)
        return None

    except subprocess.TimeoutExpired:
        logger.warning("GNINA timed out for %s", ligand_sdf)
        return None
    except Exception as exc:
        logger.warning("GNINA execution error: %s", exc)
        return None


def _parse_gnina_output(stdout: str, output_sdf: str) -> list[dict]:
    """Parse GNINA stdout and output SDF for scores.

    GNINA v1.1 with ``--cnn_scoring rescore`` outputs 5 columns:
        mode |  affinity  |  intramol  |    CNN     |   CNN
             | (kcal/mol) | (kcal/mol) | pose score | affinity
        -----+------------+------------+------------+----------
          1       -8.3         -0.32       0.82        7.1

    Without rescore, outputs 4 columns (no intramol):
        mode |  affinity  |    CNN     |   CNN
             | (kcal/mol) | pose score | affinity
          1       -8.3        0.82        7.1

    Parameters
    ----------
    stdout : str
        GNINA standard output text.
    output_sdf : str
        Path to the docked SDF output.

    Returns
    -------
    list[dict]
        Parsed poses with ``vina_score``, ``cnn_score``, ``cnn_affinity``, ``pose_path``.
    """
    poses: list[dict] = []

    # Detect whether output has 5 columns (with intramol) or 4 columns
    has_intramol = "intramol" in stdout

    # Parse table from stdout
    for line in stdout.split("\n"):
        line = line.strip()
        if has_intramol:
            # 5-column format: mode, affinity, intramol, CNN score, CNN affinity
            match = re.match(
                r"^\s*(\d+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", line
            )
            if match:
                mode = int(match.group(1))
                vina_score = float(match.group(2))
                # group(3) is intramol energy -- skip
                cnn_score = float(match.group(4))
                cnn_affinity = float(match.group(5))
                poses.append({
                    "mode": mode,
                    "vina_score": vina_score,
                    "cnn_score": cnn_score,
                    "cnn_affinity": cnn_affinity,
                    "pose_path": output_sdf if Path(output_sdf).exists() else None,
                })
        else:
            # 4-column format: mode, affinity, CNN score, CNN affinity
            match = re.match(
                r"^\s*(\d+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", line
            )
            if match:
                mode = int(match.group(1))
                vina_score = float(match.group(2))
                cnn_score = float(match.group(3))
                cnn_affinity = float(match.group(4))
                poses.append({
                    "mode": mode,
                    "vina_score": vina_score,
                    "cnn_score": cnn_score,
                    "cnn_affinity": cnn_affinity,
                    "pose_path": output_sdf if Path(output_sdf).exists() else None,
                })

    return poses


# ---------------------------------------------------------------------------
# Single-ligand docking (GNINA -> Vina -> Mock)
# ---------------------------------------------------------------------------

def run_docking(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    center: tuple[float, float, float],
    size: tuple[float, float, float] = (20.0, 20.0, 20.0),
    exhaustiveness: int = 8,
    docking_engine: str = "auto",
) -> dict:
    """Dock a single ligand against a receptor.

    Engine selection:
      - ``"gnina"``: Only try GNINA, fall back to mock if unavailable.
      - ``"vina"``: Only try Vina, fall back to mock if unavailable.
      - ``"auto"`` (default): Try GNINA -> Vina -> mock.

    Parameters
    ----------
    receptor_pdbqt : Path
        Receptor PDBQT file.
    ligand_pdbqt : Path
        Ligand PDBQT file.
    center : tuple
        (x, y, z) centre of the search box.
    size : tuple
        (sx, sy, sz) dimensions of the search box in Angstroms.
    exhaustiveness : int
        Search exhaustiveness parameter.
    docking_engine : str
        Engine choice: ``"gnina"``, ``"vina"``, or ``"auto"`` (default).

    Returns
    -------
    dict
        Keys: ``affinity`` (float, kcal/mol), ``pose_pdbqt_path`` (Path or None),
        ``vina_score`` (float), ``cnn_score`` (float), ``cnn_affinity`` (float),
        ``docking_engine`` (str: "gnina", "vina", or "mock").
    """
    engine = docking_engine.lower().strip() if docking_engine else "auto"

    # ----- Engine: gnina_gpu (single molecule: same as gnina local) -----
    # GPU batch docking is handled at dock_all_ligands() level.
    # For single-molecule calls, gnina_gpu falls back to local GNINA.
    if engine == "gnina_gpu":
        engine = "gnina"

    # ----- Engine: gnina (only GNINA, then mock) -----
    if engine == "gnina":
        if _check_gnina():
            gnina_result = _try_gnina_single(
                receptor_pdbqt, ligand_pdbqt, center, size, exhaustiveness
            )
            if gnina_result is not None:
                return gnina_result
        logger.warning("GNINA requested but unavailable or failed; falling back to mock")
        return _run_mock_docking(ligand_pdbqt)

    # ----- Engine: vina (only Vina, then mock) -----
    if engine == "vina":
        if _check_vina():
            result = _run_vina_real(
                receptor_pdbqt, ligand_pdbqt, center, size, exhaustiveness
            )
            if result is not None:
                return result
        logger.warning("Vina requested but unavailable or failed; falling back to mock")
        return _run_mock_docking(ligand_pdbqt)

    # ----- Engine: auto (GNINA -> Vina -> mock) -----
    if _check_gnina():
        # GNINA works with PDB receptor and SDF ligand, but can also use PDBQT
        gnina_result = _try_gnina_single(
            receptor_pdbqt, ligand_pdbqt, center, size, exhaustiveness
        )
        if gnina_result is not None:
            return gnina_result

    if _check_vina():
        result = _run_vina_real(
            receptor_pdbqt, ligand_pdbqt, center, size, exhaustiveness
        )
        if result is not None:
            return result
        logger.warning("Real Vina run failed; falling back to mock")

    # ----- Mock fallback -----
    return _run_mock_docking(ligand_pdbqt, center=center)


def _extract_first_model_sdf(path: Path) -> None:
    """Keep only the first molecule from a multi-model SDF file.

    GNINA outputs multi-model SDF files where all poses are concatenated.
    3Dmol.js may render an incomplete fragment if the full file is loaded.
    This function truncates the file to contain only the first model.

    Parameters
    ----------
    path : Path
        Path to the SDF file to truncate in-place.
    """
    try:
        content = path.read_text()
        if "$$$$" in content:
            first_model = content.split("$$$$")[0] + "$$$$\n"
            path.write_text(first_model)
    except Exception as e:
        logger.warning("Could not extract first SDF model: %s", e)


def _try_gnina_single(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    exhaustiveness: int,
) -> Optional[dict]:
    """Attempt GNINA docking for a single ligand, return result dict or None.

    GNINA accepts both SDF and PDBQT input. SDF is preferred because
    PDBQT torsion trees from Meeko/obabel can cause "ligand outside box"
    errors even when coordinates are correct. SDF avoids this issue.
    """
    try:
        # Prefer SDF input: GNINA handles SDF natively and avoids PDBQT
        # torsion tree issues that cause "ligand outside box" errors.
        ligand_input = str(ligand_pdbqt)
        sdf_sibling = ligand_pdbqt.with_suffix(".sdf")
        if sdf_sibling.exists() and sdf_sibling.stat().st_size > 50:
            ligand_input = str(sdf_sibling)
            logger.debug("Using SDF input for GNINA: %s", sdf_sibling.name)

        poses = dock_gnina(
            ligand_sdf=ligand_input,
            receptor_pdb=str(receptor_pdbqt),
            center=list(center),
            size=list(size),
            exhaustiveness=exhaustiveness,
        )
        if poses and len(poses) > 0:
            best = poses[0]

            # Sanity check: positive Vina score = repulsive interaction = failed docking
            if best["vina_score"] > 0:
                logger.warning(
                    "GNINA returned positive Vina score (%.2f) for %s — "
                    "rejecting (likely steric clash or bad receptor prep)",
                    best["vina_score"], ligand_pdbqt.name,
                )
                return None

            # GNINA outputs docked poses as SDF — use directly for 3D viewing
            pose_sdf = ligand_pdbqt.parent / f"{ligand_pdbqt.stem}.docked.sdf"
            pose_path = pose_sdf  # Prefer SDF for 3D rendering
            if best.get("pose_path") and Path(best["pose_path"]).exists():
                try:
                    # Keep only first model from multi-model SDF
                    src = Path(best["pose_path"])
                    if src != pose_sdf:
                        shutil.copy2(src, pose_sdf)
                    _extract_first_model_sdf(pose_sdf)
                    pose_path = pose_sdf
                except Exception:
                    pass

            # Clamp CNN score to valid [0, 1] range
            raw_cnn = best.get("cnn_score", 0.0) or 0.0
            clamped_cnn = max(0.0, min(1.0, raw_cnn))
            cnn_aff = best.get("cnn_affinity", 0.0) or 0.0
            # CNN_VS composite (CNNscore × CNNaffinity) per CACHE Challenge / Sunseri 2021
            cnn_vs = round(clamped_cnn * cnn_aff, 4)
            return {
                "affinity": best["vina_score"],
                "pose_pdbqt_path": pose_path if pose_path.exists() else None,
                "vina_score": best["vina_score"],
                "cnn_score": round(clamped_cnn, 3),
                "cnn_affinity": cnn_aff,
                "cnn_vs": cnn_vs,
                "docking_engine": "gnina",
            }
    except Exception as exc:
        logger.warning("GNINA single-ligand docking failed: %s", exc)
    return None


# ---------------------------------------------------------------------------
# Vina real execution
# ---------------------------------------------------------------------------

def _run_vina_real(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    exhaustiveness: int,
) -> Optional[dict]:
    """Execute the actual Vina binary."""
    out_path = ligand_pdbqt.parent / f"{ligand_pdbqt.stem}_out.pdbqt"

    # NOTE: Vina 1.2+ does NOT support --log. Output goes to stdout.
    cmd = [
        "vina",
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", f"{center[0]:.3f}",
        "--center_y", f"{center[1]:.3f}",
        "--center_z", f"{center[2]:.3f}",
        "--size_x", f"{size[0]:.1f}",
        "--size_y", f"{size[1]:.1f}",
        "--size_z", f"{size[2]:.1f}",
        "--exhaustiveness", str(exhaustiveness),
        "--out", str(out_path),
    ]

    try:
        logger.debug("Vina command: %s", " ".join(cmd))
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
        )

        if proc.returncode != 0:
            logger.warning(
                "Vina exited %d: %s", proc.returncode, proc.stderr[:500]
            )
            return None

        # Parse affinity from stdout (Vina 1.2+ prints results table there)
        affinity = _parse_vina_affinity(proc.stdout)
        if affinity is None:
            affinity = _parse_vina_affinity_from_pdbqt(out_path)

        if affinity is None:
            logger.warning("Could not parse Vina affinity from stdout or PDBQT")
            return None

        return {
            "affinity": affinity,
            "pose_pdbqt_path": out_path if out_path.exists() else None,
            "vina_score": affinity,
            "cnn_score": 0.0,  # Not available from Vina
            "cnn_affinity": 0.0,  # Not available from Vina
            "docking_engine": "vina",
        }

    except subprocess.TimeoutExpired:
        logger.warning("Vina timed out for %s", ligand_pdbqt.name)
        return None
    except Exception as exc:
        logger.warning("Vina execution error: %s", exc)
        return None


def _parse_vina_affinity(stdout: str) -> Optional[float]:
    """Extract the best binding affinity from Vina stdout.

    Vina 1.2+ prints a results table to stdout like:
        mode |   affinity | dist from best mode
                          | rmsd l.b.| rmsd u.b.
        -----+------------+----------+----------
           1       -8.3          0.000      0.000
           2       -7.9          1.234      2.345
    """
    # Match lines like: "   1       -8.3      0.000      0.000"
    # The key change: allow multiple spaces between mode number and affinity
    matches = re.findall(
        r"^\s*(\d+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)",
        stdout,
        re.MULTILINE,
    )
    if matches:
        # First match is mode 1 (best pose)
        try:
            return float(matches[0][1])
        except (ValueError, IndexError):
            pass

    # Fallback: "Affinity: -8.3 (kcal/mol)"
    m = re.search(r"Affinity:\s*([-\d.]+)", stdout)
    if m:
        return float(m.group(1))

    return None


def _parse_vina_affinity_from_pdbqt(pdbqt_path: Path) -> Optional[float]:
    """Parse affinity from REMARK lines in the output PDBQT."""
    if not pdbqt_path.exists():
        return None
    try:
        for line in pdbqt_path.open():
            # REMARK VINA RESULT:   -8.3      0.000      0.000
            if "VINA RESULT" in line:
                parts = line.split()
                for p in parts:
                    try:
                        val = float(p)
                        if -20 < val < 0:
                            return val
                    except ValueError:
                        continue
        return None
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Mock docking with 3-score system (V5bis)
# ---------------------------------------------------------------------------

def _run_mock_docking(ligand_pdbqt: Path, center: tuple = None) -> dict:
    """Generate a mock docking result with 3 GNINA-style scores.

    Uses deterministic hash-based generation so repeated runs produce
    the same results for the same input.

    Mock docking does NOT produce a valid 3D pose — no pose file is
    generated because coordinates would be scientifically meaningless.
    Frontend shows 2D structure only when pose_pdbqt_path is None.

    Parameters
    ----------
    ligand_pdbqt : Path
        Ligand file used for hash-based score generation.
    center : tuple, optional
        Ignored (kept for API compat).

    Returns
    -------
    dict
        Mock result with ``affinity``, ``vina_score``, ``cnn_score``,
        ``cnn_affinity``, ``pose_pdbqt_path`` (None), ``docking_engine``.
    """
    # Use hash of ligand stem for deterministic scores
    stem = ligand_pdbqt.stem
    seed = int(hashlib.md5(stem.encode()).hexdigest()[:8], 16)
    rng = random.Random(seed)

    # Vina-style affinity
    vina_score = round(rng.uniform(-12.0, -4.0), 1)

    # CNN score (pose confidence 0-1)
    seed2 = int(hashlib.md5(stem.encode()).hexdigest()[8:16], 16)
    rng2 = random.Random(seed2)
    cnn_score = round(rng2.uniform(0.3, 0.95), 3)

    # CNN affinity (pK units)
    seed3 = int(hashlib.md5(stem.encode()).hexdigest()[16:24], 16)
    rng3 = random.Random(seed3)
    cnn_affinity = round(rng3.uniform(5.0, 9.5), 2)

    # Mock docking: NO pose file — coordinates would be meaningless
    pose_path = None

    logger.debug(
        "Mock docking for %s: vina=%.1f cnn_score=%.3f cnn_aff=%.2f",
        stem, vina_score, cnn_score, cnn_affinity,
    )
    return {
        "affinity": vina_score,
        "pose_pdbqt_path": pose_path,
        "vina_score": vina_score,
        "cnn_score": cnn_score,
        "cnn_affinity": cnn_affinity,
        "docking_engine": "mock",
    }


# ---------------------------------------------------------------------------
# Consensus ranking (V5bis)
# ---------------------------------------------------------------------------

def consensus_rank(molecules: list[dict]) -> list[dict]:
    """Compute consensus ranking from 3 GNINA scores.

    For each molecule, three independent ranks are computed:
      1. Rank by vina_score (lower = better)
      2. Rank by cnn_score (higher = better)
      3. Rank by cnn_affinity (higher = better)

    The consensus rank is the mean of the three ranks. A robustness filter
    flags molecules that appear in the top 50 of at least 2 of 3 methods.

    Parameters
    ----------
    molecules : list[dict]
        Must have ``vina_score``, ``cnn_score``, ``cnn_affinity`` keys
        (or defaults to 0).

    Returns
    -------
    list[dict]
        Same list, modified in-place with added keys:
        ``rank_vina``, ``rank_cnn_score``, ``rank_cnn_affinity``,
        ``consensus_rank``, ``consensus_robust``.
        Sorted by consensus_rank ascending (best first).
    """
    n = len(molecules)
    if n == 0:
        return molecules

    # Sort indices for each score
    by_vina = sorted(range(n), key=lambda i: molecules[i].get("vina_score", 0))
    by_cnn = sorted(range(n), key=lambda i: -molecules[i].get("cnn_score", 0))
    by_cnn_aff = sorted(
        range(n), key=lambda i: -molecules[i].get("cnn_affinity", 0)
    )

    # Build rank lookup dicts for O(1) access
    vina_rank: dict[int, int] = {idx: rank + 1 for rank, idx in enumerate(by_vina)}
    cnn_rank: dict[int, int] = {idx: rank + 1 for rank, idx in enumerate(by_cnn)}
    aff_rank: dict[int, int] = {idx: rank + 1 for rank, idx in enumerate(by_cnn_aff)}

    for i, mol in enumerate(molecules):
        r1 = vina_rank[i]
        r2 = cnn_rank[i]
        r3 = aff_rank[i]

        mol["rank_vina"] = r1
        mol["rank_cnn_score"] = r2
        mol["rank_cnn_affinity"] = r3
        mol["consensus_rank"] = round((r1 + r2 + r3) / 3, 2)

        # Robustness filter: in top 50 of at least 2/3 methods
        top50 = min(50, n)
        in_top = sum(1 for r in [r1, r2, r3] if r <= top50)
        mol["consensus_robust"] = in_top >= 2

    # Sort by consensus_rank (lower = better)
    molecules.sort(key=lambda m: m.get("consensus_rank", float("inf")))
    return molecules


# ---------------------------------------------------------------------------
# Batch docking
# ---------------------------------------------------------------------------

def dock_all_ligands(
    receptor_pdbqt: Path,
    ligands: list[dict],
    center: tuple[float, float, float],
    work_dir: Path,
    progress_callback: Optional[Callable[[int, str], None]] = None,
    size: tuple[float, float, float] = (20.0, 20.0, 20.0),
    exhaustiveness: int = 8,
    docking_engine: str = "auto",
) -> list[dict]:
    """Prepare, dock, and consensus-rank all ligands.

    Engine selection is delegated to :func:`run_docking`:
      - ``"gnina"``: Only try GNINA, fall back to mock if unavailable.
      - ``"vina"``: Only try Vina, fall back to mock if unavailable.
      - ``"auto"`` (default): Try GNINA -> Vina -> mock.

    Parameters
    ----------
    receptor_pdbqt : Path
        Prepared receptor.
    ligands : list[dict]
        Each entry must have ``name`` and ``smiles`` keys.
    center : tuple
        Pocket centre (x, y, z).
    work_dir : Path
        Working directory for intermediate files.
    progress_callback : callable, optional
        Called with (percent: int, message: str) after each ligand.
    size : tuple
        Search box size.
    exhaustiveness : int
        Search exhaustiveness.
    docking_engine : str
        Engine choice: ``"gnina"``, ``"vina"``, or ``"auto"`` (default).

    Returns
    -------
    list[dict]
        Each entry has: ``name``, ``smiles``, ``affinity``, ``vina_score``,
        ``cnn_score``, ``cnn_affinity``, ``consensus_rank``, ``consensus_robust``,
        ``pose_pdbqt_path``, ``source``, ``docking_engine``, plus all input keys.
        Sorted by consensus_rank (best first), then by affinity.
    """
    ligand_dir = work_dir / "ligands"
    ligand_dir.mkdir(parents=True, exist_ok=True)

    # --- GPU batch docking (gnina_gpu or auto with GPU available) ---
    engine = docking_engine.lower().strip() if docking_engine else "auto"
    if engine in ("gnina_gpu", "auto"):
        try:
            from pipeline.docking_gpu import is_gpu_available, dock_all_runpod_batch
            if engine == "gnina_gpu" or (engine == "auto" and is_gpu_available()):
                if is_gpu_available():
                    logger.info("Attempting GPU batch docking (%d ligands)", len(ligands))
                    gpu_results = dock_all_runpod_batch(
                        receptor_pdbqt=receptor_pdbqt,
                        ligands=ligands,
                        center=center,
                        work_dir=work_dir,
                        progress_callback=progress_callback,
                        size=size,
                        exhaustiveness=exhaustiveness,
                    )
                    if gpu_results:
                        logger.info("GPU batch docking succeeded: %d results", len(gpu_results))
                        return gpu_results
                    logger.warning("GPU batch returned empty results, falling back to local")
                else:
                    if engine == "gnina_gpu":
                        logger.warning("gnina_gpu requested but RUNPOD_API_KEY not set, falling back to local")
        except Exception as exc:
            logger.warning("GPU batch docking failed: %s, falling back to local", exc)

    results: list[dict] = []
    total = len(ligands)

    for i, lig in enumerate(ligands):
        name = lig.get("name", f"ligand_{i}")
        smiles = lig.get("smiles", "")

        if not smiles:
            logger.warning("Skipping ligand %s: no SMILES", name)
            continue

        # Report progress
        if progress_callback is not None:
            pct = int(40 + (i / max(total, 1)) * 50)  # 40% -> 90%
            progress_callback(pct, f"Docking {name} ({i+1}/{total})")

        # Prepare ligand PDBQT
        pdbqt_path = prepare_ligand(smiles, name, ligand_dir)
        if pdbqt_path is None:
            logger.warning("Skipping %s: preparation failed", name)
            continue

        # Dock using the selected engine
        try:
            dock_result = run_docking(
                receptor_pdbqt, pdbqt_path, center, size, exhaustiveness,
                docking_engine=docking_engine,
            )
        except Exception as exc:
            logger.warning("Docking failed for %s: %s", name, exc)
            continue

        # Clamp CNN score to valid [0, 1] range
        raw_cnn = dock_result.get("cnn_score", 0.0) or 0.0
        clamped_cnn = max(0.0, min(1.0, raw_cnn))

        results.append({
            **lig,  # preserve original keys
            "affinity": dock_result["affinity"],
            "vina_score": dock_result.get("vina_score", dock_result["affinity"]),
            "cnn_score": round(clamped_cnn, 3),
            "cnn_affinity": dock_result.get("cnn_affinity", 0.0),
            "docking_engine": dock_result.get("docking_engine", "unknown"),
            "pose_pdbqt_path": str(dock_result["pose_pdbqt_path"])
            if dock_result.get("pose_pdbqt_path") else None,
        })

    # Apply consensus ranking across all results
    if results:
        results = consensus_rank(results)

    logger.info(
        "Docking complete: %d/%d ligands succeeded (engine=%s)",
        len(results),
        total,
        docking_engine,
    )
    return results
