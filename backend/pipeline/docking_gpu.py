"""
DockIt pipeline -- GNINA GPU docking via RunPod serverless.

Sends batches of ligands as a multi-molecule SDF to a RunPod GPU endpoint
running GNINA, instead of docking one-by-one on CPU.

Optimized for throughput:
- Parallel batch submission (concurrent.futures)
- Parallel SDF preparation (ProcessPoolExecutor)
- Adaptive polling (1s → 5s ramp)
- Configurable batch size (default 250)
"""

from __future__ import annotations

import logging
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Callable, Optional

import requests

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Configuration from environment
# ---------------------------------------------------------------------------

RUNPOD_API_KEY = os.environ.get("RUNPOD_API_KEY", "")
GNINA_ENDPOINT_ID = os.environ.get("GNINA_ENDPOINT_ID", "weeeiy6z4jdsv3")
RUNPOD_TIMEOUT = int(os.environ.get("RUNPOD_TIMEOUT", "600"))
BATCH_SIZE = int(os.environ.get("GPU_BATCH_SIZE", "250"))
MAX_PARALLEL_BATCHES = int(os.environ.get("GPU_MAX_PARALLEL", "8"))


def is_gpu_available() -> bool:
    """Check if RunPod GPU docking is available (API key + endpoint configured)."""
    key = os.environ.get("RUNPOD_API_KEY", "")
    endpoint = os.environ.get("GNINA_ENDPOINT_ID", "weeeiy6z4jdsv3")
    available = bool(key and endpoint)
    if available:
        logger.info("GPU docking available (endpoint=%s)", endpoint)
    else:
        logger.debug("GPU docking not available (RUNPOD_API_KEY=%s)", "set" if key else "missing")
    return available


# ---------------------------------------------------------------------------
# SDF preparation (runs in worker processes for parallelism)
# ---------------------------------------------------------------------------

def _smiles_to_sdf_block(args: tuple) -> Optional[tuple[int, str, dict]]:
    """Convert a single SMILES to SDF block. Designed for ProcessPoolExecutor.

    Args: (index, ligand_dict)
    Returns: (index, sdf_block, ligand_dict) or None on failure.
    """
    idx, lig = args
    smiles = lig.get("smiles", "")
    name = lig.get("name", "unknown")
    if not smiles:
        return None

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        status = AllChem.EmbedMolecule(mol, params)
        if status != 0:
            params2 = AllChem.ETKDGv3()
            params2.useRandomCoords = True
            AllChem.EmbedMolecule(mol, params2)
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            pass
        mol.SetProp("_Name", name)
        block = Chem.MolToMolBlock(mol)
        return (idx, block + "$$$$\n", lig)
    except Exception:
        return None


def _prepare_sdf_parallel(ligands: list[dict]) -> tuple[list[str], list[dict]]:
    """Prepare SDF blocks from SMILES using parallel processes.

    For small batches (<50), runs sequentially to avoid thread overhead.
    For larger batches, uses ThreadPoolExecutor (RDKit releases the GIL
    during conformer generation, so threads give real parallelism).
    Note: ProcessPoolExecutor cannot be used inside Celery daemon workers.
    """
    n = len(ligands)
    args = [(i, lig) for i, lig in enumerate(ligands)]

    if n < 50:
        # Sequential — faster for small batches
        results = [_smiles_to_sdf_block(a) for a in args]
    else:
        # Threaded — works inside Celery daemons, RDKit releases GIL
        n_workers = min(os.cpu_count() or 4, 8)
        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            results = list(pool.map(_smiles_to_sdf_block, args))

    # Collect successful results, preserving order
    valid = [r for r in results if r is not None]
    valid.sort(key=lambda x: x[0])
    sdf_blocks = [r[1] for r in valid]
    valid_ligands = [r[2] for r in valid]
    return sdf_blocks, valid_ligands


# ---------------------------------------------------------------------------
# Low-level RunPod API: submit + poll + cancel
# ---------------------------------------------------------------------------

def _cancel_job(job_id: str) -> bool:
    """Cancel a RunPod job. Returns True if cancelled successfully."""
    api_key = os.environ.get("RUNPOD_API_KEY", RUNPOD_API_KEY)
    endpoint_id = os.environ.get("GNINA_ENDPOINT_ID", GNINA_ENDPOINT_ID)
    url = f"https://api.runpod.ai/v2/{endpoint_id}/cancel/{job_id}"
    headers = {"Authorization": f"Bearer {api_key}"}
    try:
        resp = requests.post(url, headers=headers, timeout=10)
        ok = resp.status_code in (200, 204)
        if ok:
            logger.info("Cancelled RunPod job %s", job_id)
        else:
            logger.warning("Failed to cancel RunPod job %s: %s", job_id, resp.status_code)
        return ok
    except Exception as e:
        logger.warning("Error cancelling RunPod job %s: %s", job_id, e)
        return False

def _submit_job(
    receptor_pdb_content: str,
    ligands_sdf_content: str,
    center: list[float],
    size: list[float],
    exhaustiveness: int = 8,
    num_modes: int = 9,
    energy_range: float = 3.0,
    cnn_scoring: str = "rescore",
    seed: int = 0,
    autobox_add: float = 4.0,
) -> tuple[Optional[str], Optional[str]]:
    """Submit a docking job to RunPod. Returns (job_id, error)."""
    api_key = os.environ.get("RUNPOD_API_KEY", RUNPOD_API_KEY)
    endpoint_id = os.environ.get("GNINA_ENDPOINT_ID", GNINA_ENDPOINT_ID)

    url = f"https://api.runpod.ai/v2/{endpoint_id}/run"
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
    }
    payload = {
        "input": {
            "receptor_pdb": receptor_pdb_content,
            "ligand_sdf": ligands_sdf_content,
            "center_x": center[0],
            "center_y": center[1],
            "center_z": center[2],
            "size_x": size[0],
            "size_y": size[1],
            "size_z": size[2],
            "exhaustiveness": exhaustiveness,
            "num_modes": num_modes,
            "energy_range": energy_range,
            "cnn_scoring": cnn_scoring,
            "autobox_add": autobox_add,
            "seed": seed,
            "min_rmsd_filter": 1.0,
        }
    }

    try:
        resp = requests.post(url, json=payload, headers=headers, timeout=30)
        if resp.status_code == 402:
            return None, "no_credits"
        resp.raise_for_status()
        data = resp.json()
        job_id = data.get("id")
        if not job_id:
            return None, "no_job_id"
        return job_id, None
    except requests.exceptions.RequestException as e:
        return None, f"network: {e}"


def _poll_job(job_id: str, timeout: int = 600) -> dict:
    """Poll a RunPod job until completion. Adaptive polling: 1s → 3s → 5s."""
    api_key = os.environ.get("RUNPOD_API_KEY", RUNPOD_API_KEY)
    endpoint_id = os.environ.get("GNINA_ENDPOINT_ID", GNINA_ENDPOINT_ID)
    headers = {"Authorization": f"Bearer {api_key}"}
    status_url = f"https://api.runpod.ai/v2/{endpoint_id}/status/{job_id}"

    start_time = time.time()
    poll_count = 0

    while True:
        elapsed = time.time() - start_time
        if elapsed > timeout:
            return {"error": f"timeout ({timeout}s)", "results": []}

        # Adaptive polling: fast at first, slower as job runs longer
        if poll_count < 5:
            time.sleep(1)
        elif poll_count < 15:
            time.sleep(2)
        else:
            time.sleep(4)
        poll_count += 1

        try:
            poll = requests.get(status_url, headers=headers, timeout=15)
            poll.raise_for_status()
            status_data = poll.json()
        except (requests.exceptions.ConnectionError, requests.exceptions.Timeout):
            continue
        except Exception:
            continue

        status = status_data.get("status", "UNKNOWN")

        if status == "COMPLETED":
            output = status_data.get("output", {})
            return {
                "results": output.get("results", []),
                "docked_sdf": output.get("docked_sdf", ""),
                "n_molecules": output.get("n_molecules", 0),
                "elapsed_seconds": output.get("elapsed_seconds", elapsed),
                "gpu_used": output.get("gpu_used", "unknown"),
                "error": None,
            }

        if status == "FAILED":
            error_msg = status_data.get("error", "unknown error")
            return {"error": f"GPU handler error: {error_msg}", "results": []}

    return {"error": "poll loop exited unexpectedly", "results": []}


def dock_batch_gpu(
    receptor_pdb_content: str,
    ligands_sdf_content: str,
    center: list[float],
    size: list[float],
    exhaustiveness: int = 8,
    num_modes: int = 9,
    energy_range: float = 3.0,
    cnn_scoring: str = "rescore",
    seed: int = 0,
    autobox_add: float = 4.0,
    progress_callback: Optional[Callable[[int, str], None]] = None,
) -> dict:
    """Send a batch of ligands to RunPod GPU for GNINA docking.

    Parameters
    ----------
    receptor_pdb_content : str
        Raw PDB file content of the receptor.
    ligands_sdf_content : str
        Multi-molecule SDF content (all ligands concatenated).
    center : list[float]
        [x, y, z] centre of the search box.
    size : list[float]
        [sx, sy, sz] dimensions of the search box.
    exhaustiveness : int
        GNINA search exhaustiveness (8-256).
    num_modes : int
        Number of docking poses per molecule.
    energy_range : float
        Max energy difference from best pose (kcal/mol).
    cnn_scoring : str
        CNN scoring mode: "rescore" (fast) or "refinement" (better poses).
    seed : int
        Random seed for reproducibility (0 = random).
    autobox_add : float
        Padding around pocket center (Angstroms).

    Returns
    -------
    dict
        Keys: results, docked_sdf, n_molecules, elapsed_seconds, error.
    """
    try:
        job_id, error = _submit_job(
            receptor_pdb_content, ligands_sdf_content,
            center, size, exhaustiveness, num_modes,
            energy_range, cnn_scoring, seed, autobox_add,
        )
        if error:
            logger.error("RunPod submit failed: %s", error)
            return {"error": error, "results": []}

        logger.info("RunPod job submitted: %s", job_id)

        if progress_callback:
            progress_callback(42, "GPU docking in progress...")

        timeout = int(os.environ.get("RUNPOD_TIMEOUT", str(RUNPOD_TIMEOUT)))
        result = _poll_job(job_id, timeout=timeout)

        if result.get("error"):
            logger.error("RunPod job %s failed: %s", job_id, result["error"])
        else:
            logger.info("RunPod job %s completed: %d molecules", job_id, result.get("n_molecules", 0))

        return result

    except Exception as e:
        logger.error("RunPod unexpected error: %s", e)
        return {"error": str(e), "results": []}


# ---------------------------------------------------------------------------
# Parse GPU batch results into standardized dicts
# ---------------------------------------------------------------------------

def _parse_batch_results(
    gpu_result: dict,
    batch_ligands: list[dict],
) -> list[dict]:
    """Parse RunPod GPU results into standardized molecule dicts."""
    batch_results_raw = gpu_result.get("results", [])

    # Group by molecule index
    mol_results: dict[int, list] = {}
    for r in batch_results_raw:
        mol_idx = r.get("mol_index", r.get("molecule_index", 0))
        mol_results.setdefault(mol_idx, []).append(r)

    parsed = []
    for local_idx, lig in enumerate(batch_ligands):
        poses = mol_results.get(local_idx, [])
        if not poses:
            # Fallback: flat list ordered by molecule
            if local_idx < len(batch_results_raw):
                poses = [batch_results_raw[local_idx]]

        if not poses:
            logger.warning("GPU: no results for %s (idx %d)", lig.get("name"), local_idx)
            continue

        best = min(poses, key=lambda p: p.get("minimizedAffinity", p.get("affinity", 0)))

        vina_score = best.get("minimizedAffinity", best.get("affinity", 0.0))
        cnn_score_raw = best.get("CNNscore", best.get("cnn_score", 0.0))
        cnn_affinity = best.get("CNNaffinity", best.get("cnn_affinity", 0.0))

        if vina_score > 0:
            logger.warning("GPU: positive Vina score (%.2f) for %s — rejecting", vina_score, lig.get("name"))
            continue

        cnn_score = max(0.0, min(1.0, cnn_score_raw or 0.0))
        cnn_vs = round(cnn_score * (cnn_affinity or 0.0), 4)

        parsed.append({
            **lig,
            "affinity": vina_score,
            "vina_score": vina_score,
            "cnn_score": round(cnn_score, 3),
            "cnn_affinity": cnn_affinity,
            "cnn_vs": cnn_vs,
            "docking_engine": "gnina_gpu",
            "pose_pdbqt_path": None,
        })

    return parsed


def _attach_pose_molblocks(
    all_results: list[dict],
    docked_sdf_contents: list[str],
    work_dir: Path,
) -> None:
    """Parse docked SDFs and attach pose molblocks to results (in-place)."""
    name_to_molblock: dict[str, str] = {}

    for sdf_idx, sdf_content in enumerate(docked_sdf_contents):
        if not sdf_content:
            continue
        try:
            # Save to disk
            sdf_path = work_dir / f"docked_gpu_{sdf_idx}.sdf"
            sdf_path.write_text(sdf_content)

            entries = sdf_content.split("$$$$")
            for entry in entries:
                entry = entry.strip()
                if not entry:
                    continue
                first_line = entry.split("\n", 1)[0].strip()
                mol_name = first_line if first_line else ""
                if mol_name and mol_name not in name_to_molblock:
                    name_to_molblock[mol_name] = entry + "\n$$$$\n"
        except Exception as e:
            logger.warning("GPU: could not parse docked SDF %d: %s", sdf_idx, e)

    attached = 0
    for result in all_results:
        name = result.get("name", "")
        if name in name_to_molblock:
            result["pose_molblock"] = name_to_molblock[name]
            attached += 1

    if name_to_molblock:
        logger.info("GPU: attached %d/%d pose molblocks", attached, len(all_results))


# ---------------------------------------------------------------------------
# High-level batch docking with parallel submission
# ---------------------------------------------------------------------------

def dock_all_runpod_batch(
    receptor_pdbqt: Path,
    ligands: list[dict],
    center: tuple[float, float, float],
    work_dir: Path,
    progress_callback: Optional[Callable[[int, str], None]] = None,
    size: tuple[float, float, float] = (20.0, 20.0, 20.0),
    exhaustiveness: int = 8,
    num_modes: int = 9,
    energy_range: float = 3.0,
    cnn_scoring: str = "rescore",
    seed: int = 0,
    autobox_add: float = 4.0,
) -> list[dict]:
    """Dock all ligands via RunPod GPU with parallel batch submission.

    Optimizations vs sequential version:
    - SDF preparation parallelized across CPU cores
    - Multiple batches submitted to RunPod simultaneously
    - Adaptive polling (1s→4s ramp)
    - Larger default batch size (250 vs 100)

    Parameters
    ----------
    receptor_pdbqt : Path
        Receptor file (PDB or PDBQT -- reads .pdb sibling if available).
    ligands : list[dict]
        Each must have 'name' and 'smiles' keys.
    center : tuple
        Pocket centre (x, y, z).
    work_dir : Path
        Working directory.
    progress_callback : callable, optional
        Called with (percent, message).
    size : tuple
        Search box size.
    exhaustiveness : int
        Search exhaustiveness.

    Returns
    -------
    list[dict]
        Same format as dock_all_ligands: name, smiles, affinity, vina_score,
        cnn_score, cnn_affinity, consensus_rank, docking_engine, etc.
    """
    from pipeline.docking import consensus_rank

    overall_start = time.time()

    # --- 1. Read receptor PDB content ---
    receptor_path = Path(receptor_pdbqt)
    pdb_path = receptor_path.with_suffix(".pdb")
    if not pdb_path.exists():
        pdb_candidates = list(receptor_path.parent.glob("*.pdb"))
        if pdb_candidates:
            pdb_path = pdb_candidates[0]
        else:
            pdb_path = receptor_path

    try:
        receptor_content = pdb_path.read_text()
    except Exception as e:
        logger.error("Cannot read receptor file %s: %s", pdb_path, e)
        return []

    logger.info("GPU batch docking: receptor=%s, %d ligands", pdb_path.name, len(ligands))

    # --- 2. Parallel SDF preparation ---
    if progress_callback:
        progress_callback(30, f"Preparing {len(ligands)} molecules (3D conformers)...")

    sdf_blocks, valid_ligands = _prepare_sdf_parallel(ligands)

    if not sdf_blocks:
        logger.warning("GPU: no valid SDF blocks generated")
        return []

    total = len(sdf_blocks)
    prep_elapsed = time.time() - overall_start
    logger.info("GPU: prepared %d/%d ligands as SDF in %.1fs", total, len(ligands), prep_elapsed)

    if progress_callback:
        progress_callback(38, f"Prepared {total} molecules, submitting to GPU...")

    # --- 3. Split into batches ---
    batch_size = BATCH_SIZE
    batches = []
    for i in range(0, total, batch_size):
        end = min(i + batch_size, total)
        batches.append({
            "sdf": "".join(sdf_blocks[i:end]),
            "ligands": valid_ligands[i:end],
            "idx": len(batches),
        })

    n_batches = len(batches)
    logger.info("GPU: %d batches of ~%d molecules (max %d parallel)",
                n_batches, batch_size, MAX_PARALLEL_BATCHES)

    # --- 4. Submit all batches in parallel + poll ---
    timeout = int(os.environ.get("RUNPOD_TIMEOUT", str(RUNPOD_TIMEOUT)))

    # Shared docking params
    dock_params = dict(
        center=list(center), size=list(size),
        exhaustiveness=exhaustiveness, num_modes=num_modes,
        energy_range=energy_range, cnn_scoring=cnn_scoring,
        seed=seed, autobox_add=autobox_add,
    )

    # Track all submitted job IDs for cleanup on failure
    submitted_jobs: list[str] = []
    import threading
    jobs_lock = threading.Lock()

    def _run_single_batch(batch: dict) -> tuple[int, dict, list[dict]]:
        """Submit and poll a single batch. Returns (batch_idx, gpu_result, batch_ligands)."""
        job_id, error = _submit_job(
            receptor_pdb_content=receptor_content,
            ligands_sdf_content=batch["sdf"],
            **dock_params,
        )
        if error:
            return (batch["idx"], {"error": error, "results": []}, batch["ligands"])

        with jobs_lock:
            submitted_jobs.append(job_id)

        logger.info("GPU batch %d/%d submitted: job=%s (%d mols)",
                     batch["idx"] + 1, n_batches, job_id, len(batch["ligands"]))

        result = _poll_job(job_id, timeout=timeout)

        # Remove from tracked jobs on completion
        with jobs_lock:
            if job_id in submitted_jobs:
                submitted_jobs.remove(job_id)

        return (batch["idx"], result, batch["ligands"])

    def _cancel_all_pending():
        """Cancel all tracked RunPod jobs that haven't completed."""
        with jobs_lock:
            to_cancel = list(submitted_jobs)
        if to_cancel:
            logger.warning("Cancelling %d pending RunPod jobs...", len(to_cancel))
            for jid in to_cancel:
                _cancel_job(jid)

    all_results = []
    docked_sdfs = []
    completed_batches = 0

    try:
        if n_batches == 1:
            # Single batch — no threading overhead
            batch_idx, gpu_result, batch_ligs = _run_single_batch(batches[0])
            if gpu_result.get("error"):
                raise RuntimeError(f"GPU docking failed: {gpu_result['error']}")
            all_results.extend(_parse_batch_results(gpu_result, batch_ligs))
            docked_sdfs.append(gpu_result.get("docked_sdf", ""))
            if progress_callback:
                progress_callback(85, f"GPU docking done: {len(all_results)} results")
        else:
            # Parallel batch submission
            max_workers = min(n_batches, MAX_PARALLEL_BATCHES)
            with ThreadPoolExecutor(max_workers=max_workers) as pool:
                futures = {pool.submit(_run_single_batch, b): b["idx"] for b in batches}

                for future in as_completed(futures):
                    batch_idx, gpu_result, batch_ligs = future.result()
                    completed_batches += 1

                    if gpu_result.get("error"):
                        error = gpu_result["error"]
                        logger.error("GPU batch %d failed: %s", batch_idx + 1, error)
                        # Cancel remaining jobs before raising
                        _cancel_all_pending()
                        raise RuntimeError(f"GPU batch {batch_idx + 1}/{n_batches} failed: {error}")

                    parsed = _parse_batch_results(gpu_result, batch_ligs)
                    all_results.extend(parsed)
                    docked_sdfs.append(gpu_result.get("docked_sdf", ""))

                    logger.info("GPU batch %d/%d done: %d results",
                               batch_idx + 1, n_batches, len(parsed))

                    if progress_callback:
                        pct = int(40 + (completed_batches / n_batches) * 45)
                        progress_callback(
                            pct,
                            f"GPU: {completed_batches}/{n_batches} batches done "
                            f"({len(all_results)} molecules)",
                        )
    except Exception:
        # Ensure all RunPod jobs are cancelled on any failure
        _cancel_all_pending()
        raise

    # --- 5. Attach pose molblocks from all docked SDFs ---
    _attach_pose_molblocks(all_results, docked_sdfs, work_dir)

    # --- 6. Consensus ranking ---
    if all_results:
        all_results = consensus_rank(all_results)

    total_elapsed = time.time() - overall_start
    throughput = len(all_results) / total_elapsed if total_elapsed > 0 else 0
    logger.info(
        "GPU docking complete: %d/%d ligands in %.1fs (%.0f mol/min, %d batches)",
        len(all_results), len(ligands), total_elapsed, throughput * 60, n_batches,
    )

    if progress_callback:
        progress_callback(
            88,
            f"GPU docking done: {len(all_results)} results in {total_elapsed:.0f}s "
            f"({throughput * 60:.0f} mol/min)",
        )

    return all_results
