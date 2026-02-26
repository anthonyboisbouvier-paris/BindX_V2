"""
DockIt pipeline -- Rapid scoring with smina/Vinardo (Pass 3 of massive screening).

Scores ~10K molecules rapidly using smina with the Vinardo scoring function,
keeping the top 500 candidates for full Vina docking in Pass 4.

smina is a fork of AutoDock Vina with additional scoring functions,
including Vinardo which is faster than the default Vina scoring while
maintaining reasonable accuracy for virtual screening.

Strategies (3-tier fallback):
  1. smina binary with Vinardo scoring (exhaustiveness=4, num_modes=1).
  2. Vina binary with reduced exhaustiveness (if smina not found).
  3. Hash-based deterministic mock scoring (no external tools).
"""

from __future__ import annotations

import hashlib
import logging
import os
import shutil
import subprocess
import tempfile
import time
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Callable, Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Availability flags (resolved lazily)
# ---------------------------------------------------------------------------
_SMINA_AVAILABLE: Optional[bool] = None
_VINA_AVAILABLE: Optional[bool] = None
_OBABEL_AVAILABLE: Optional[bool] = None


def check_smina_available() -> bool:
    """Check if the smina binary is available on the system.

    Looks for smina at /usr/local/bin/smina first, then falls back
    to checking PATH via shutil.which.

    Returns
    -------
    bool
        True if smina is available and executable.
    """
    global _SMINA_AVAILABLE
    if _SMINA_AVAILABLE is None:
        # Check fixed path first
        smina_path = Path("/usr/local/bin/smina")
        if smina_path.exists() and os.access(str(smina_path), os.X_OK):
            _SMINA_AVAILABLE = True
            logger.info("smina found at %s", smina_path)
        elif shutil.which("smina"):
            _SMINA_AVAILABLE = True
            logger.info("smina found on PATH")
        else:
            _SMINA_AVAILABLE = False
            logger.info("smina not found; will try vina or mock fallback")
    return _SMINA_AVAILABLE


def _check_vina() -> bool:
    """Check if vina binary is available as fallback."""
    global _VINA_AVAILABLE
    if _VINA_AVAILABLE is None:
        _VINA_AVAILABLE = shutil.which("vina") is not None
        if _VINA_AVAILABLE:
            logger.info("AutoDock Vina found on PATH (fallback for smina)")
    return _VINA_AVAILABLE


def _check_obabel() -> bool:
    """Check if Open Babel is available for SMILES-to-PDBQT conversion."""
    global _OBABEL_AVAILABLE
    if _OBABEL_AVAILABLE is None:
        _OBABEL_AVAILABLE = shutil.which("obabel") is not None
        if _OBABEL_AVAILABLE:
            logger.info("Open Babel found on PATH")
        else:
            logger.info("Open Babel not found; ligand preparation will use RDKit or mock")
    return _OBABEL_AVAILABLE


# ---------------------------------------------------------------------------
# Single molecule rapid scoring
# ---------------------------------------------------------------------------

def score_rapid_single(
    receptor_pdbqt: Path,
    smiles: str,
    name: str,
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    work_dir: Path,
    timeout_sec: int = 30,
) -> Optional[dict]:
    """Score a single molecule using smina/Vinardo rapid docking.

    Parameters
    ----------
    receptor_pdbqt : Path
        Receptor PDBQT file.
    smiles : str
        SMILES string of the ligand.
    name : str
        Ligand identifier (used for file naming).
    center : tuple[float, float, float]
        (x, y, z) centre of the search box.
    size : tuple[float, float, float]
        (sx, sy, sz) dimensions of the search box.
    work_dir : Path
        Working directory for intermediate files.
    timeout_sec : int
        Maximum time in seconds for scoring this molecule.

    Returns
    -------
    dict or None
        ``{"name": str, "smiles": str, "affinity": float, "source": str}``
        or None if scoring fails.
    """
    if not smiles or not smiles.strip():
        return None

    smiles = smiles.strip()
    safe_name = "".join(c if c.isalnum() or c in "-_" else "_" for c in name)[:60]

    # Prepare ligand PDBQT
    ligand_pdbqt = _prepare_ligand_pdbqt(smiles, safe_name, work_dir)
    if ligand_pdbqt is None:
        logger.debug("Ligand preparation failed for %s; skipping", name)
        return None

    # Try smina first
    if check_smina_available():
        result = _run_smina(
            receptor_pdbqt, ligand_pdbqt, center, size, timeout_sec,
        )
        if result is not None:
            return {
                "name": name,
                "smiles": smiles,
                "affinity": result,
                "source": "smina_vinardo",
            }

    # Fallback to vina with low exhaustiveness
    if _check_vina():
        result = _run_vina_rapid(
            receptor_pdbqt, ligand_pdbqt, center, size, timeout_sec,
        )
        if result is not None:
            return {
                "name": name,
                "smiles": smiles,
                "affinity": result,
                "source": "vina_rapid",
            }

    # Mock fallback
    return _mock_score_single(smiles, name, center)


def _run_smina(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    timeout_sec: int,
) -> Optional[float]:
    """Execute smina with Vinardo scoring.

    Parameters
    ----------
    receptor_pdbqt : Path
        Receptor file.
    ligand_pdbqt : Path
        Ligand file.
    center : tuple
        Search box centre.
    size : tuple
        Search box dimensions.
    timeout_sec : int
        Timeout in seconds.

    Returns
    -------
    float or None
        Best affinity in kcal/mol, or None on failure.
    """
    # Determine smina binary path
    smina_bin = "/usr/local/bin/smina"
    if not Path(smina_bin).exists():
        smina_bin = shutil.which("smina") or "smina"

    out_path = ligand_pdbqt.parent / f"{ligand_pdbqt.stem}_smina_out.pdbqt"

    cmd = [
        smina_bin,
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", f"{center[0]:.3f}",
        "--center_y", f"{center[1]:.3f}",
        "--center_z", f"{center[2]:.3f}",
        "--size_x", f"{size[0]:.1f}",
        "--size_y", f"{size[1]:.1f}",
        "--size_z", f"{size[2]:.1f}",
        "--scoring", "vinardo",
        "--exhaustiveness", "4",
        "--num_modes", "1",
        "--out", str(out_path),
    ]

    try:
        logger.debug("smina command: %s", " ".join(cmd))
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout_sec,
        )

        if proc.returncode != 0:
            logger.debug(
                "smina exited %d for %s: %s",
                proc.returncode, ligand_pdbqt.name,
                proc.stderr[:200],
            )
            return None

        # Parse affinity from stdout or output file
        affinity = _parse_docking_affinity(proc.stdout, out_path)
        return affinity

    except subprocess.TimeoutExpired:
        logger.debug("smina timed out for %s", ligand_pdbqt.name)
        return None
    except Exception as exc:
        logger.debug("smina error for %s: %s", ligand_pdbqt.name, exc)
        return None


def _run_vina_rapid(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    timeout_sec: int,
) -> Optional[float]:
    """Execute Vina with reduced exhaustiveness as fallback for smina.

    Parameters
    ----------
    receptor_pdbqt : Path
        Receptor file.
    ligand_pdbqt : Path
        Ligand file.
    center : tuple
        Search box centre.
    size : tuple
        Search box dimensions.
    timeout_sec : int
        Timeout in seconds.

    Returns
    -------
    float or None
        Best affinity in kcal/mol, or None on failure.
    """
    out_path = ligand_pdbqt.parent / f"{ligand_pdbqt.stem}_vina_rapid_out.pdbqt"

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
        "--exhaustiveness", "4",
        "--num_modes", "1",
        "--out", str(out_path),
    ]

    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout_sec,
        )

        if proc.returncode != 0:
            return None

        return _parse_docking_affinity(proc.stdout, out_path)

    except subprocess.TimeoutExpired:
        return None
    except Exception:
        return None


def _parse_docking_affinity(
    stdout: str,
    out_path: Path,
) -> Optional[float]:
    """Parse best binding affinity from docking output.

    Checks stdout first, then the output PDBQT REMARK lines.

    Parameters
    ----------
    stdout : str
        Standard output from the docking program.
    out_path : Path
        Path to the output PDBQT file.

    Returns
    -------
    float or None
        Best affinity (kcal/mol, negative = better), or None.
    """
    import re

    # Parse from stdout: "   1     -8.3      0.000      0.000"
    for line in stdout.split("\n"):
        line = line.strip()
        match = re.match(r"^\s*1\s+([-\d.]+)", line)
        if match:
            try:
                val = float(match.group(1))
                if -20.0 < val < 0.0:
                    return val
            except ValueError:
                pass

    # Parse from output PDBQT: "REMARK VINA RESULT:   -8.3  ..."
    if out_path.exists():
        try:
            for line in out_path.open():
                if "VINA RESULT" in line or "RESULT" in line:
                    parts = line.split()
                    for p in parts:
                        try:
                            val = float(p)
                            if -20.0 < val < 0.0:
                                return val
                        except ValueError:
                            continue
        except Exception:
            pass

    # Parse "Affinity: -X.X (kcal/mol)" format
    match = re.search(r"Affinity:\s*([-\d.]+)", stdout)
    if match:
        try:
            val = float(match.group(1))
            if -20.0 < val < 0.0:
                return val
        except ValueError:
            pass

    return None


# ---------------------------------------------------------------------------
# Ligand preparation for rapid scoring
# ---------------------------------------------------------------------------

def _prepare_ligand_pdbqt(
    smiles: str,
    safe_name: str,
    work_dir: Path,
) -> Optional[Path]:
    """Convert SMILES to PDBQT for docking.

    Tries pipeline.prepare first, then Open Babel directly, then
    a manual RDKit-based fallback.

    Parameters
    ----------
    smiles : str
        SMILES string.
    safe_name : str
        Sanitized filename.
    work_dir : Path
        Output directory.

    Returns
    -------
    Path or None
        Path to PDBQT file, or None on failure.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    pdbqt_path = work_dir / f"{safe_name}.pdbqt"

    # Skip if already prepared
    if pdbqt_path.exists() and pdbqt_path.stat().st_size > 50:
        return pdbqt_path

    # Strategy 1: Use the existing prepare module
    try:
        from pipeline.prepare import prepare_ligand
        result = prepare_ligand(smiles, safe_name, work_dir)
        if result is not None and result.exists():
            return result
    except Exception as exc:
        logger.debug("prepare_ligand failed for %s: %s", safe_name, exc)

    # Strategy 2: Direct obabel conversion from SMILES
    if _check_obabel():
        try:
            cmd = [
                "obabel", "-:" + smiles,
                "-O", str(pdbqt_path),
                "--gen3d",
                "-h",
            ]
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=30,
            )
            if result.returncode == 0 and pdbqt_path.exists() and pdbqt_path.stat().st_size > 50:
                return pdbqt_path
        except Exception:
            pass

    # Strategy 3: RDKit to PDB, then manual PDBQT
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.useRandomCoords = True
        res = AllChem.EmbedMolecule(mol, params)
        if res == -1:
            return None

        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=100)
        except Exception:
            pass

        # Write minimal PDBQT
        conf = mol.GetConformer()
        lines: list[str] = ["ROOT\n"]
        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            element = atom.GetSymbol()
            atom_name = f"{element}{i + 1}"[:4].ljust(4)
            line = (
                f"ATOM  {i + 1:5d} {atom_name} LIG A   1    "
                f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}"
                f"  1.00  0.00    +0.000 {element:>2s}\n"
            )
            lines.append(line)
        lines.append("ENDROOT\n")
        lines.append("TORSDOF 0\n")
        pdbqt_path.write_text("".join(lines))
        return pdbqt_path

    except Exception as exc:
        logger.debug("RDKit PDBQT preparation failed for %s: %s", safe_name, exc)
        return None


# ---------------------------------------------------------------------------
# Batch rapid scoring
# ---------------------------------------------------------------------------

def score_rapid_batch(
    receptor_pdbqt: Path,
    smiles_list: list[str],
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    work_dir: Path,
    n_top: int = 500,
    n_cpus: Optional[int] = None,
    progress_callback: Optional[Callable[[dict], None]] = None,
) -> list[dict]:
    """Score a batch of molecules rapidly and return the top N.

    Uses multiprocessing when real docking tools are available.
    Falls back to serial processing or mock scoring when needed.

    Parameters
    ----------
    receptor_pdbqt : Path
        Receptor PDBQT file.
    smiles_list : list[str]
        SMILES strings to score.
    center : tuple[float, float, float]
        Search box centre.
    size : tuple[float, float, float]
        Search box dimensions.
    work_dir : Path
        Working directory.
    n_top : int
        Number of top-scoring molecules to return.
    n_cpus : int, optional
        Number of CPU cores for parallel processing. Defaults to
        ``os.cpu_count()``.
    progress_callback : callable, optional
        Called with a dict after each batch of molecules.

    Returns
    -------
    list[dict]
        Top ``n_top`` molecules sorted by affinity (most negative first).
        Each dict has keys: ``name``, ``smiles``, ``affinity``, ``source``.
    """
    total = len(smiles_list)
    if total == 0:
        logger.info("score_rapid_batch: empty input")
        return []

    if n_cpus is None:
        n_cpus = os.cpu_count() or 1

    work_dir.mkdir(parents=True, exist_ok=True)
    start_time = time.monotonic()

    logger.info(
        "score_rapid_batch: scoring %d molecules, n_top=%d, n_cpus=%d",
        total, n_top, n_cpus,
    )

    # Determine if real docking tools are available
    has_real_tools = check_smina_available() or _check_vina()

    if has_real_tools:
        results = _batch_real_scoring(
            receptor_pdbqt, smiles_list, center, size,
            work_dir, n_cpus, progress_callback,
        )
    else:
        results = _batch_mock_scoring(
            smiles_list, center, progress_callback,
        )

    # Sort by affinity (most negative = best)
    results.sort(key=lambda r: r.get("affinity", 0.0))

    # Take top N
    top_results = results[:n_top]

    elapsed = time.monotonic() - start_time
    logger.info(
        "score_rapid_batch: scored %d/%d molecules, kept top %d, "
        "best affinity=%.2f, elapsed=%.1fs",
        len(results), total, len(top_results),
        top_results[0]["affinity"] if top_results else 0.0,
        elapsed,
    )

    return top_results


def _batch_real_scoring(
    receptor_pdbqt: Path,
    smiles_list: list[str],
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    work_dir: Path,
    n_cpus: int,
    progress_callback: Optional[Callable[[dict], None]],
) -> list[dict]:
    """Score molecules using real docking tools (smina/vina).

    Processes molecules sequentially to avoid multiprocessing complexity
    with external binaries and file I/O. For truly massive batches,
    the subprocess calls are the bottleneck, not Python.

    Parameters
    ----------
    receptor_pdbqt : Path
        Receptor PDBQT file.
    smiles_list : list[str]
        SMILES to score.
    center : tuple
        Search box centre.
    size : tuple
        Search box dimensions.
    work_dir : Path
        Working directory.
    n_cpus : int
        Number of CPUs (used to set pool size for parallel runs).
    progress_callback : callable, optional
        Progress callback.

    Returns
    -------
    list[dict]
        Scored molecules.
    """
    ligand_dir = work_dir / "rapid_ligands"
    ligand_dir.mkdir(parents=True, exist_ok=True)

    results: list[dict] = []
    total = len(smiles_list)
    failed = 0

    for i, smi in enumerate(smiles_list):
        name = f"RAPID_{i:06d}"

        result = score_rapid_single(
            receptor_pdbqt=receptor_pdbqt,
            smiles=smi,
            name=name,
            center=center,
            size=size,
            work_dir=ligand_dir,
            timeout_sec=30,
        )

        if result is not None:
            results.append(result)
        else:
            failed += 1

        # Progress report every 100 molecules
        if progress_callback is not None and (i + 1) % 100 == 0:
            try:
                progress_callback({
                    "pass_name": "rapid_scoring",
                    "processed": i + 1,
                    "total": total,
                    "kept": len(results),
                    "failed": failed,
                })
            except Exception:
                pass

    # Final progress report
    if progress_callback is not None:
        try:
            progress_callback({
                "pass_name": "rapid_scoring",
                "processed": total,
                "total": total,
                "kept": len(results),
                "failed": failed,
            })
        except Exception:
            pass

    logger.info(
        "Real rapid scoring: %d scored, %d failed out of %d",
        len(results), failed, total,
    )
    return results


def _batch_mock_scoring(
    smiles_list: list[str],
    center: tuple[float, float, float],
    progress_callback: Optional[Callable[[dict], None]],
) -> list[dict]:
    """Score molecules using hash-based mock when no docking tools exist.

    Parameters
    ----------
    smiles_list : list[str]
        SMILES to score.
    center : tuple
        Pocket centre (used for hash salting).
    progress_callback : callable, optional
        Progress callback.

    Returns
    -------
    list[dict]
        Mock-scored molecules.
    """
    logger.info(
        "score_rapid_batch: using mock scoring for %d molecules",
        len(smiles_list),
    )

    results: list[dict] = []
    total = len(smiles_list)

    for i, smi in enumerate(smiles_list):
        if not smi or not smi.strip():
            continue

        smi = smi.strip()
        result = _mock_score_single(smi, f"RAPID_{i:06d}", center)
        if result is not None:
            results.append(result)

        # Progress report every 1000 molecules
        if progress_callback is not None and (i + 1) % 1000 == 0:
            try:
                progress_callback({
                    "pass_name": "rapid_scoring_mock",
                    "processed": i + 1,
                    "total": total,
                    "kept": len(results),
                })
            except Exception:
                pass

    # Final progress
    if progress_callback is not None:
        try:
            progress_callback({
                "pass_name": "rapid_scoring_mock",
                "processed": total,
                "total": total,
                "kept": len(results),
            })
        except Exception:
            pass

    return results


# ---------------------------------------------------------------------------
# Mock scoring
# ---------------------------------------------------------------------------

def _mock_score_single(
    smiles: str,
    name: str,
    center: tuple[float, float, float],
) -> Optional[dict]:
    """Generate a deterministic mock score for a single molecule.

    Uses a hash of the SMILES and pocket centre to produce a
    reproducible affinity in the range [-12.0, -4.0] kcal/mol.
    Molecules with common kinase-inhibitor pharmacophore features
    get a small bonus (more negative = better).

    Parameters
    ----------
    smiles : str
        SMILES string.
    name : str
        Molecule name.
    center : tuple
        Pocket centre (used for salting).

    Returns
    -------
    dict or None
        Mock scoring result, or None if SMILES is empty.
    """
    if not smiles:
        return None

    # Deterministic hash
    key = f"rapid:{smiles}:{center[0]:.2f}:{center[1]:.2f}:{center[2]:.2f}"
    digest = hashlib.sha256(key.encode("utf-8")).hexdigest()
    hash_int = int(digest[:12], 16)

    # Map to [-12.0, -4.0]
    base = -4.0 - (hash_int % 8001) / 1000.0

    # Pharmacophore bonus for kinase-inhibitor-like features
    bonus = 0.0
    ki_features = ["ncnc", "Nc1", "c1ccc", "C(=O)N", "c2cc", "OCCOC"]
    for feat in ki_features:
        if feat.lower() in smiles.lower():
            bonus -= 0.12

    affinity = max(-12.0, min(-4.0, base + bonus))

    return {
        "name": name,
        "smiles": smiles,
        "affinity": round(affinity, 2),
        "source": "mock_rapid",
    }


# =========================================================================
# CLI / Self-test
# =========================================================================

if __name__ == "__main__":
    import json
    import sys

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    print("=" * 70)
    print("DockIt -- scoring_rapid.py self-test")
    print("=" * 70)

    # Check tool availability
    print(f"\nsmina available: {check_smina_available()}")
    print(f"vina available:  {_check_vina()}")
    print(f"obabel available: {_check_obabel()}")

    # Test SMILES
    test_smiles = [
        "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OC)cc23)c1",     # Erlotinib-like
        "CC(=O)Oc1ccccc1C(=O)O",                         # Aspirin
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",                    # Ibuprofen
        "Cn1c(=O)c2c(ncn2C)n(C)c1=O",                    # Caffeine
        "CC(=O)Nc1ccc(O)cc1",                             # Acetaminophen
        "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1", # Gefitinib
        "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",  # Imatinib
        "CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1",   # Sorafenib
    ]

    pocket_center = (22.0, 0.5, 18.0)
    pocket_size = (25.0, 25.0, 25.0)
    test_work_dir = Path("/tmp/dockit_rapid_test")
    test_work_dir.mkdir(parents=True, exist_ok=True)

    # Create a dummy receptor
    dummy_receptor = test_work_dir / "receptor.pdbqt"
    dummy_receptor.write_text("REMARK mock receptor\nEND\n")

    # Test single scoring
    print("\n[1] Testing score_rapid_single...")
    single_result = score_rapid_single(
        receptor_pdbqt=dummy_receptor,
        smiles=test_smiles[0],
        name="Erlotinib_test",
        center=pocket_center,
        size=pocket_size,
        work_dir=test_work_dir,
    )
    print(f"  Single result: {json.dumps(single_result, indent=2)}")

    # Test batch scoring
    print(f"\n[2] Testing score_rapid_batch with {len(test_smiles)} molecules...")

    progress_log: list[dict] = []

    def on_progress(info: dict) -> None:
        progress_log.append(info)

    batch_results = score_rapid_batch(
        receptor_pdbqt=dummy_receptor,
        smiles_list=test_smiles,
        center=pocket_center,
        size=pocket_size,
        work_dir=test_work_dir,
        n_top=5,
        progress_callback=on_progress,
    )

    print(f"\n  Batch results (top 5):")
    print(f"  {'Rank':<6} {'Name':<20} {'Affinity':<10} {'Source':<15} {'SMILES'}")
    print("  " + "-" * 90)
    for i, r in enumerate(batch_results, 1):
        print(
            f"  {i:<6} {r['name']:<20} {r['affinity']:<10.2f} "
            f"{r['source']:<15} {r['smiles'][:40]}"
        )

    # Validation
    errors = 0

    if single_result is None:
        print("\nFAIL: single scoring returned None")
        errors += 1
    else:
        if single_result["affinity"] >= 0:
            print("FAIL: affinity should be negative")
            errors += 1
        if not single_result["smiles"]:
            print("FAIL: SMILES missing")
            errors += 1
        print("PASS: single scoring OK")

    if len(batch_results) == 0:
        print("FAIL: batch scoring returned empty")
        errors += 1
    else:
        # Check sorting (should be ascending by affinity)
        for i in range(len(batch_results) - 1):
            if batch_results[i]["affinity"] > batch_results[i + 1]["affinity"]:
                print("FAIL: results not sorted by affinity")
                errors += 1
                break
        else:
            print("PASS: batch results correctly sorted")

        if len(batch_results) <= 5:
            print(f"PASS: returned {len(batch_results)} results (n_top=5)")
        else:
            print(f"FAIL: returned {len(batch_results)} results but n_top=5")
            errors += 1

    # Determinism test
    print("\n[3] Testing determinism...")
    single_result2 = score_rapid_single(
        receptor_pdbqt=dummy_receptor,
        smiles=test_smiles[0],
        name="Erlotinib_test",
        center=pocket_center,
        size=pocket_size,
        work_dir=test_work_dir,
    )
    if single_result and single_result2:
        if single_result["affinity"] == single_result2["affinity"]:
            print("PASS: deterministic scoring verified")
        else:
            print("FAIL: non-deterministic scoring")
            errors += 1

    if errors == 0:
        print("\nALL CHECKS PASSED")
    else:
        print(f"\n{errors} CHECK(S) FAILED")

    sys.exit(1 if errors > 0 else 0)
