"""
DockIt pipeline -- GNINA GPU docking via RunPod serverless.

Sends batches of ligands as a multi-molecule SDF to a RunPod GPU endpoint
running GNINA, instead of docking one-by-one on CPU.

Performance: ~50 molecules in 2-5 min (GPU) vs 25-50 min (CPU).
"""

from __future__ import annotations

import logging
import os
import time
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
BATCH_SIZE = int(os.environ.get("GPU_BATCH_SIZE", "100"))


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
# Low-level RunPod API call
# ---------------------------------------------------------------------------

def dock_batch_gpu(
    receptor_pdb_content: str,
    ligands_sdf_content: str,
    center: list[float],
    size: list[float],
    exhaustiveness: int = 8,
    num_modes: int = 9,
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
        GNINA search exhaustiveness.
    num_modes : int
        Number of docking poses per molecule.

    Returns
    -------
    dict
        Keys: results, docked_sdf, n_molecules, elapsed_seconds, error.
    """
    api_key = os.environ.get("RUNPOD_API_KEY", RUNPOD_API_KEY)
    endpoint_id = os.environ.get("GNINA_ENDPOINT_ID", GNINA_ENDPOINT_ID)
    timeout = int(os.environ.get("RUNPOD_TIMEOUT", str(RUNPOD_TIMEOUT)))

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
            "cnn_scoring": "rescore",
            # Use default CNN ensemble (5 models) — do NOT specify a single
            # model. Matches CPU config per McNutt et al. 2021.
            "seed": 0,
            "min_rmsd_filter": 1.0,
        }
    }

    try:
        # Submit job
        resp = requests.post(url, json=payload, headers=headers, timeout=30)
        if resp.status_code == 402:
            logger.warning("RunPod: no credits (402)")
            return {"error": "no_credits", "results": []}
        resp.raise_for_status()
        data = resp.json()
        job_id = data.get("id")
        if not job_id:
            logger.error("RunPod: no job_id in response: %s", data)
            return {"error": "no_job_id", "results": []}

        logger.info("RunPod job submitted: %s", job_id)

        # Poll for results
        status_url = f"https://api.runpod.ai/v2/{endpoint_id}/status/{job_id}"
        start_time = time.time()
        cold_start_warned = False

        while True:
            elapsed = time.time() - start_time
            if elapsed > timeout:
                logger.warning("RunPod job %s timed out after %ds", job_id, timeout)
                return {"error": "timeout", "results": []}

            time.sleep(5)

            try:
                poll = requests.get(status_url, headers=headers, timeout=15)
                poll.raise_for_status()
                status_data = poll.json()
            except Exception as poll_err:
                logger.warning("RunPod poll error: %s", poll_err)
                continue

            status = status_data.get("status", "UNKNOWN")

            if status == "COMPLETED":
                output = status_data.get("output", {})
                logger.info(
                    "RunPod job %s completed in %.1fs: %d molecules",
                    job_id, elapsed, output.get("n_molecules", 0),
                )
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
                logger.error("RunPod job %s FAILED: %s", job_id, error_msg)
                return {"error": f"handler_failed: {error_msg}", "results": []}

            # Cold start warning
            if status == "IN_QUEUE" and elapsed > 30 and not cold_start_warned:
                logger.info("RunPod job %s: GPU warming up (cold start)...", job_id)
                cold_start_warned = True

    except requests.exceptions.RequestException as e:
        logger.error("RunPod network error: %s", e)
        return {"error": f"network: {e}", "results": []}
    except Exception as e:
        logger.error("RunPod unexpected error: %s", e)
        return {"error": str(e), "results": []}


# ---------------------------------------------------------------------------
# High-level batch docking: replaces dock_all_ligands() when GPU available
# ---------------------------------------------------------------------------

def dock_all_runpod_batch(
    receptor_pdbqt: Path,
    ligands: list[dict],
    center: tuple[float, float, float],
    work_dir: Path,
    progress_callback: Optional[Callable[[int, str], None]] = None,
    size: tuple[float, float, float] = (20.0, 20.0, 20.0),
    exhaustiveness: int = 8,
) -> list[dict]:
    """Dock all ligands via RunPod GPU in batch mode.

    Prepares a multi-molecule SDF from SMILES, sends to RunPod in batches
    of BATCH_SIZE, parses results, applies consensus ranking.

    Returns results in the same format as dock_all_ligands().

    Parameters
    ----------
    receptor_pdbqt : Path
        Receptor file (PDB or PDBQT -- we read the .pdb sibling if available).
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

    # --- 1. Read receptor PDB content ---
    # Prefer .pdb over .pdbqt (RunPod handler expects PDB)
    receptor_path = Path(receptor_pdbqt)
    pdb_path = receptor_path.with_suffix(".pdb")
    if not pdb_path.exists():
        # Try broader search: any .pdb in the same directory
        pdb_candidates = list(receptor_path.parent.glob("*.pdb"))
        if pdb_candidates:
            pdb_path = pdb_candidates[0]
        else:
            pdb_path = receptor_path  # Use PDBQT as fallback

    try:
        receptor_content = pdb_path.read_text()
    except Exception as e:
        logger.error("Cannot read receptor file %s: %s", pdb_path, e)
        return []

    logger.info("GPU batch docking: receptor=%s, %d ligands", pdb_path.name, len(ligands))

    # --- 2. Generate 3D conformers and build multi-molecule SDF ---
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        logger.error("RDKit not available for GPU batch SDF generation")
        return []

    sdf_blocks = []
    valid_ligands = []  # Track which ligands succeeded

    for lig in ligands:
        smiles = lig.get("smiles", "")
        name = lig.get("name", "unknown")
        if not smiles:
            continue

        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.warning("GPU: invalid SMILES for %s, skipping", name)
                continue
            mol = Chem.AddHs(mol)
            # Generate 3D coordinates
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            status = AllChem.EmbedMolecule(mol, params)
            if status != 0:
                # Fallback: random coords
                AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            try:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            except Exception:
                pass
            mol.SetProp("_Name", name)
            block = Chem.MolToMolBlock(mol)
            sdf_blocks.append(block + "$$$$\n")
            valid_ligands.append(lig)
        except Exception as e:
            logger.warning("GPU: failed to generate 3D for %s: %s", name, e)
            continue

    if not sdf_blocks:
        logger.warning("GPU: no valid SDF blocks generated")
        return []

    total = len(sdf_blocks)
    logger.info("GPU: prepared %d/%d ligands as SDF", total, len(ligands))

    # --- 3. Split into batches and submit ---
    batch_size = BATCH_SIZE
    all_results = []
    n_batches = (total + batch_size - 1) // batch_size

    for batch_idx in range(n_batches):
        start = batch_idx * batch_size
        end = min(start + batch_size, total)
        batch_sdf = "".join(sdf_blocks[start:end])
        batch_ligands = valid_ligands[start:end]

        if progress_callback:
            pct = int(40 + (batch_idx / max(n_batches, 1)) * 45)
            msg = f"GPU docking batch {batch_idx + 1}/{n_batches} ({end - start} molecules)"
            if batch_idx == 0 and n_batches == 1:
                msg = f"GPU docking {total} molecules..."
            progress_callback(pct, msg)

        logger.info("GPU batch %d/%d: %d molecules", batch_idx + 1, n_batches, end - start)

        gpu_result = dock_batch_gpu(
            receptor_pdb_content=receptor_content,
            ligands_sdf_content=batch_sdf,
            center=list(center),
            size=list(size),
            exhaustiveness=exhaustiveness,
        )

        if gpu_result.get("error"):
            logger.error("GPU batch %d failed: %s", batch_idx + 1, gpu_result["error"])
            return []  # Signal failure so caller falls back to CPU

        # --- 4. Parse results ---
        batch_results_raw = gpu_result.get("results", [])

        # Results are ordered by molecule, then by rank
        # Group by molecule index
        mol_results: dict[int, list] = {}
        for r in batch_results_raw:
            mol_idx = r.get("mol_index", r.get("molecule_index", 0))
            mol_results.setdefault(mol_idx, []).append(r)

        for local_idx, lig in enumerate(batch_ligands):
            poses = mol_results.get(local_idx, [])
            if not poses:
                # Try matching by rank 1 at position local_idx
                # Some handlers return flat list ordered by molecule
                flat_idx = local_idx
                if flat_idx < len(batch_results_raw):
                    poses = [batch_results_raw[flat_idx]]

            if not poses:
                logger.warning("GPU: no results for %s (idx %d)", lig.get("name"), local_idx)
                continue

            # Take best pose (rank 1 or lowest minimizedAffinity)
            best = min(poses, key=lambda p: p.get("minimizedAffinity", p.get("affinity", 0)))

            # Map fields: minimizedAffinity->vina_score, CNNscore->cnn_score, CNNaffinity->cnn_affinity
            vina_score = best.get("minimizedAffinity", best.get("affinity", 0.0))
            cnn_score_raw = best.get("CNNscore", best.get("cnn_score", 0.0))
            cnn_affinity = best.get("CNNaffinity", best.get("cnn_affinity", 0.0))

            # Reject positive Vina scores (steric clash artifacts) — same as CPU
            if vina_score > 0:
                logger.warning(
                    "GPU: positive Vina score (%.2f) for %s — rejecting",
                    vina_score, lig.get("name"),
                )
                continue

            # Clamp CNN score to [0, 1]
            cnn_score = max(0.0, min(1.0, cnn_score_raw or 0.0))

            # CNN_VS composite (CNNscore × CNNaffinity) per CACHE Challenge / Sunseri 2021
            cnn_vs = round(cnn_score * (cnn_affinity or 0.0), 4)

            all_results.append({
                **lig,  # preserve original keys (name, smiles, source, etc.)
                "affinity": vina_score,
                "vina_score": vina_score,
                "cnn_score": round(cnn_score, 3),
                "cnn_affinity": cnn_affinity,
                "cnn_vs": cnn_vs,
                "docking_engine": "gnina_gpu",
                "pose_pdbqt_path": None,  # GPU doesn't return local file paths
            })

    # --- 5. Save docked SDF if available ---
    if gpu_result.get("docked_sdf"):
        try:
            sdf_path = work_dir / "docked_gpu.sdf"
            sdf_path.write_text(gpu_result["docked_sdf"])
            logger.info("GPU: saved docked SDF to %s", sdf_path)
        except Exception as e:
            logger.warning("GPU: could not save docked SDF: %s", e)

    # --- 6. Consensus ranking ---
    if all_results:
        all_results = consensus_rank(all_results)

    elapsed = gpu_result.get("elapsed_seconds", 0)
    gpu_name = gpu_result.get("gpu_used", "unknown")
    logger.info(
        "GPU docking complete: %d/%d ligands, %.1fs on %s",
        len(all_results), len(ligands), elapsed, gpu_name,
    )

    if progress_callback:
        progress_callback(
            88,
            f"GPU docking done: {len(all_results)} results in {elapsed:.0f}s",
        )

    return all_results
