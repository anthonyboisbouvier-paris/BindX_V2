"""
DockIt — Celery tasks that orchestrate the full docking pipeline.

V1 (basic/rapid mode): 8 steps — structure, pockets, prepare, ligands, dock, score, report, persist.
V2 (advanced/standard mode): 11 steps — adds generation, ADMET, retrosynthesis.
V3: auto ligand strategy, sequence input, structure_source tracking,
    pipeline_summary, step_details, deep screening task (4h), email notifications.
V5: audit log, off-target screening, confidence scoring, lead optimization task.
V6: disorder prediction (IDR), disorder_warning for pockets in IDRs.

Each step updates the job progress in the database so the frontend
can poll for real-time status.
"""

from __future__ import annotations

import json
import logging
import os
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import requests as http_requests

from celery_app import celery_app
from database import get_db, update_job

logger = logging.getLogger(__name__)

DATA_DIR = Path(os.environ.get("DOCKIT_DATA_DIR", "/data")) / "jobs"


def _update_progress(job_id: str, progress: int, step: str, step_details: Optional[str] = None) -> None:
    """Helper to update job progress in the database.

    Parameters
    ----------
    job_id : str
        Job UUID.
    progress : int
        Progress percentage (0-100).
    step : str
        Current step name for display.
    step_details : str, optional
        Additional detail string (e.g. counts).
    """
    detail_msg = f" | {step_details}" if step_details else ""
    logger.info("[%s] %d%% -- %s%s", job_id, progress, step, detail_msg)
    update_kwargs = {"progress": progress, "current_step": step}
    if step_details:
        update_kwargs["current_step"] = f"{step} - {step_details}"
    update_job(job_id, **update_kwargs)


def _store_pipeline_summary(job_id: str, summary: dict) -> None:
    """Store the pipeline summary as JSON in the database.

    Parameters
    ----------
    job_id : str
        Job UUID.
    summary : dict
        Pipeline summary dictionary.
    """
    try:
        update_job(job_id, pipeline_summary_json=json.dumps(summary))
    except Exception as exc:
        logger.warning("[%s] Failed to store pipeline summary: %s", job_id, exc)


def _run_single_interaction(protein_path: str, ligand_path: str, uniprot_id: str, smiles: str) -> Optional[dict]:
    """Run interaction analysis for a single molecule in an isolated subprocess.

    ProLIF/MDAnalysis can SIGSEGV on certain PDB structures. Running in a
    subprocess protects the Celery worker process from being killed.
    """
    import multiprocessing as mp

    def _worker(conn, prot, lig, uid, smi):
        try:
            from pipeline.interaction_analysis import analyze_interactions
            result = analyze_interactions(
                protein_path=prot, ligand_path=lig,
                uniprot_id=uid, smiles=smi,
            )
            conn.send(result)
        except Exception as exc:
            conn.send({"error": str(exc)})
        finally:
            conn.close()

    parent_conn, child_conn = mp.Pipe()
    proc = mp.Process(target=_worker, args=(child_conn, protein_path, ligand_path, uniprot_id, smiles))
    proc.start()
    proc.join(timeout=45)  # 45s max per molecule

    if proc.is_alive():
        proc.kill()
        proc.join(timeout=5)
        return None

    if proc.exitcode != 0:
        return None

    if parent_conn.poll():
        result = parent_conn.recv()
        if isinstance(result, dict) and "error" not in result:
            return result
    return None


def _safe_interaction_analysis(job_id: str, molecules: list[dict], pdb_path: str, uniprot_id: str) -> int:
    """Run interaction analysis for a list of molecules with crash isolation.

    Also computes structured pose quality metrics (V12).
    Returns the number of molecules successfully analyzed.
    """
    n_analyzed = 0
    for mol in molecules:
        smi = mol.get("smiles", "")
        try:
            result = _run_single_interaction(
                protein_path=pdb_path,
                ligand_path=mol.get("pose_pdbqt_path", ""),
                uniprot_id=uniprot_id,
                smiles=smi,
            )
            if result:
                mol["interactions"] = result
                mol["interaction_quality"] = result.get("interaction_quality", 0.5)
                n_analyzed += 1

                # V12: compute structured pose quality
                try:
                    from pipeline.interaction_analysis import compute_pose_quality
                    ligand_path = mol.get("pose_pdbqt_path", "")
                    if ligand_path:
                        pq = compute_pose_quality(pdb_path, ligand_path, uniprot_id, smi)
                        mol["pose_quality"] = pq
                except Exception as pq_exc:
                    logger.debug("[%s] Pose quality failed for %s: %s",
                                 job_id, mol.get("name", "?"), pq_exc)
            else:
                logger.warning("[%s] Interaction analysis returned no result for %s",
                               job_id, mol.get("name", "?"))
        except Exception as mol_exc:
            logger.warning("[%s] Interaction analysis skipped for %s: %s",
                           job_id, mol.get("name", "?"), mol_exc)
    return n_analyzed


def _fetch_protein_name(uniprot_id: str) -> Optional[str]:
    """Fetch the recommended protein name from the UniProt REST API.

    Parameters
    ----------
    uniprot_id : str
        UniProt accession (e.g. ``"P00533"``).

    Returns
    -------
    str or None
        The recommended full name, or None on failure.
    """
    try:
        resp = http_requests.get(
            f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json",
            timeout=5,
        )
        resp.raise_for_status()
        data = resp.json()
        protein_desc = data.get("proteinDescription", {})
        rec_name = protein_desc.get("recommendedName", {})
        full_name = rec_name.get("fullName", {})
        name = full_name.get("value")
        if name:
            return name
        # Fallback: submittedName
        submitted = protein_desc.get("submittedName", [])
        if submitted and isinstance(submitted, list):
            return submitted[0].get("fullName", {}).get("value")
        return None
    except Exception as exc:
        logger.warning("Failed to fetch protein name for %s: %s", uniprot_id, exc)
        return None


@celery_app.task(name="run_pipeline", bind=True, max_retries=0)
def run_pipeline(self, job_id: str, params: dict) -> dict:
    """Execute the virtual-screening pipeline (rapid or standard mode).

    Parameters
    ----------
    job_id : str
        UUID of the job.
    params : dict
        Keys: ``uniprot_id``, ``sequence``, ``use_chembl``, ``use_zinc``,
        ``max_ligands``, ``smiles_list``, ``mode``, ``docking_engine``,
        ``notification_email``, ``auto_strategy``,
        ``enable_generation``, ``enable_diffdock``, ``enable_retrosynthesis``,
        ``n_generated_molecules``.

    Returns
    -------
    dict
        Summary with ``status`` and ``result_count``.
    """
    work_dir = DATA_DIR / job_id
    work_dir.mkdir(parents=True, exist_ok=True)

    # Pipeline summary accumulator
    pipeline_summary: dict = {
        "steps_completed": [],
        "timings": {},
        "counts": {},
    }

    try:
        # Mark as running
        update_job(job_id, status="running", progress=0, current_step="Starting pipeline")

        # V5: Initialize audit log
        from pipeline.audit_log import AuditLog
        audit = AuditLog(job_id)
        audit.log("pipeline", "Pipeline started", {"params": {
            k: v for k, v in params.items()
            if k not in ("smiles_list", "target_config")  # exclude potentially large data
        }})

        uniprot_id: str = params.get("uniprot_id", "") or ""
        sequence: Optional[str] = params.get("sequence")
        use_chembl: bool = params.get("use_chembl", True)
        use_zinc: bool = params.get("use_zinc", False)
        max_ligands: int = params.get("max_ligands", 50)
        smiles_list: list[str] | None = params.get("smiles_list")

        # V2 params
        mode: str = params.get("mode", "rapid")
        enable_generation: bool = params.get("enable_generation", False)
        enable_diffdock: bool = params.get("enable_diffdock", False)
        enable_retrosynthesis: bool = params.get("enable_retrosynthesis", False)
        n_generated_molecules: int = params.get("n_generated_molecules", 100)

        # V3 params
        docking_engine: str = params.get("docking_engine", "auto")
        notification_email: Optional[str] = params.get("notification_email")
        auto_strategy: bool = params.get("auto_strategy", False)

        # V8: pre-computed target config from Target Setup (skips structure+pockets)
        target_config: Optional[dict] = params.get("target_config")

        enable_dmpk = params.get("enable_dmpk", True)

        # V12: granular analysis flags (backward compat: enable_dmpk=False disables all)
        master = enable_dmpk
        enable_admet = params.get("enable_admet", master)
        enable_synthesis = params.get("enable_synthesis", master)
        enable_selectivity = params.get("enable_selectivity", master)
        enable_herg = params.get("enable_herg", master)
        enable_safety = params.get("enable_safety", master)
        if not master:
            enable_admet = enable_synthesis = enable_selectivity = enable_herg = enable_safety = False

        # V12: box_size override (None = auto from pocket)
        box_size_override = params.get("box_size")

        # ------------------------------------------------------------------
        # Step 1: Fetch / restore protein structure (0% -> 8%)
        # ------------------------------------------------------------------
        structure_source = "unknown"
        pdb_info: Optional[dict] = None
        ligand_id: Optional[str] = None

        # V8: reuse pre-computed structure from Target Setup when available
        _target_config_used = False
        if target_config:
            try:
                structures_list = target_config.get("structures") or []
                sel_idx = int(target_config.get("selected_structure_idx", 0))
                sel_structure = (
                    structures_list[sel_idx] if structures_list and sel_idx < len(structures_list)
                    else structures_list[0] if structures_list
                    else target_config.get("structure") or {}
                )
                pdb_data_inline: Optional[str] = sel_structure.get("pdb_data")
                download_url: Optional[str] = sel_structure.get("download_url")
                structure_source = sel_structure.get("source", "cached")

                pdb_path = work_dir / "protein.pdb"
                if pdb_data_inline:
                    pdb_path.write_text(pdb_data_inline)
                    _target_config_used = True
                    logger.info("[%s] Step 1 skipped — using inline pdb_data from Target Setup", job_id)
                elif download_url:
                    _update_progress(job_id, 3, "Downloading structure from Target Setup")
                    resp = http_requests.get(download_url, timeout=30)
                    resp.raise_for_status()
                    pdb_path.write_bytes(resp.content)
                    _target_config_used = True
                    logger.info("[%s] Step 1 skipped — downloaded structure from %s", job_id, download_url)
                else:
                    logger.warning("[%s] target_config has no pdb_data or download_url — re-fetching", job_id)
            except Exception as exc:
                logger.warning("[%s] Failed to restore structure from target_config: %s — re-fetching", job_id, exc)
                _target_config_used = False

        if not _target_config_used:
            _update_progress(job_id, 3, "Fetching protein structure")
            if sequence and not uniprot_id:
                # V3: direct sequence input
                from pipeline.structure import fetch_structure_from_sequence
                pdb_path, structure_source = fetch_structure_from_sequence(sequence, work_dir)
                _update_progress(
                    job_id, 8, "Structure retrieved",
                    step_details=f"Folded from sequence ({len(sequence)} residues) via {structure_source}",
                )
            else:
                from pipeline.structure import fetch_structure
                pdb_path, structure_source = fetch_structure(uniprot_id, work_dir)
                _update_progress(
                    job_id, 8, "Structure retrieved",
                    step_details=f"Source: {structure_source}",
                )

            # V5bis: Retrieve PDB experimental metadata (ligand_id, resolution, etc.)
            if uniprot_id and structure_source == "pdb_experimental":
                try:
                    from pipeline.structure import get_pdb_info
                    pdb_info = get_pdb_info(work_dir, uniprot_id)
                    if pdb_info:
                        ligand_id = pdb_info.get("ligand_id")
                        pipeline_summary["pdb_info"] = pdb_info
                        logger.info("[%s] PDB info: %s (ligand=%s)", job_id,
                                    pdb_info.get("pdb_id", "?"), ligand_id or "none")
                except Exception as exc:
                    logger.warning("[%s] Failed to retrieve PDB info: %s", job_id, exc)

        _update_progress(job_id, 8, "Structure ready", step_details=f"Source: {structure_source}")

        # Store structure source and PDB path
        update_job(job_id, pdb_path=str(pdb_path), structure_source=structure_source)
        pipeline_summary["structure_source"] = structure_source
        pipeline_summary["steps_completed"].append("structure")

        # --- A1: Compute effective_structure_source early ---
        effective_structure_source = structure_source
        if structure_source == "pdb_experimental" and pdb_info and pdb_info.get("ligand_id"):
            effective_structure_source = "pdb_holo"
        pipeline_summary["effective_structure_source"] = effective_structure_source

        audit.log("structure", f"Structure ready from {structure_source}", {
            "source": structure_source,
            "pdb_path": str(pdb_path),
            "reused_from_target_setup": _target_config_used,
        })

        # ------------------------------------------------------------------
        # Step 1b (V6): Disorder prediction
        # ------------------------------------------------------------------
        disorder_info: Optional[dict] = None
        try:
            # Get the sequence: use provided sequence, or fetch from UniProt
            disorder_sequence: Optional[str] = sequence
            if not disorder_sequence and uniprot_id:
                from pipeline.structure import _fetch_uniprot_sequence
                disorder_sequence = _fetch_uniprot_sequence(uniprot_id)

            if disorder_sequence:
                from pipeline.structure import predict_disorder
                disorder_info = predict_disorder(disorder_sequence)
                pipeline_summary["disorder_info"] = {
                    "fraction_disordered": disorder_info.get("fraction_disordered", 0.0),
                    "n_idr_regions": len(disorder_info.get("idr_regions", [])),
                    "method": disorder_info.get("method", "unknown"),
                }
                logger.info(
                    "[%s] Disorder prediction: %.1f%% disordered, %d IDR(s) (method=%s)",
                    job_id,
                    disorder_info["fraction_disordered"] * 100,
                    len(disorder_info["idr_regions"]),
                    disorder_info["method"],
                )
                audit.log("disorder", "Disorder prediction complete", {
                    "fraction_disordered": disorder_info["fraction_disordered"],
                    "n_idr_regions": len(disorder_info["idr_regions"]),
                    "method": disorder_info["method"],
                })
            else:
                logger.info("[%s] No sequence available for disorder prediction", job_id)
        except Exception as exc:
            logger.warning("[%s] Disorder prediction failed: %s", job_id, exc)

        # ------------------------------------------------------------------
        # Step 2: Detect / restore binding pockets (8% -> 15%)
        # ------------------------------------------------------------------
        _update_progress(job_id, 10, "Restoring binding site from Target Setup" if target_config else "Detecting binding pockets")

        pockets = []
        pocket_size = (25.0, 25.0, 25.0)
        _pocket_reused = False

        if target_config:
            try:
                pockets_raw = target_config.get("pockets") or []
                sel_pocket_idx = int(target_config.get("selected_pocket_idx", 0))
                sel_pocket = (
                    pockets_raw[sel_pocket_idx] if pockets_raw and sel_pocket_idx < len(pockets_raw)
                    else pockets_raw[0] if pockets_raw
                    else None
                )
                if sel_pocket and sel_pocket.get("center"):
                    center_raw = sel_pocket["center"]
                    center = tuple(float(c) for c in center_raw[:3])
                    pocket_method = sel_pocket.get("method", "cached")
                    pockets = pockets_raw
                    best_pocket = sel_pocket
                    _pocket_reused = True
                    logger.info("[%s] Step 2 skipped — using pocket from Target Setup center=%s", job_id, center)
                else:
                    logger.warning("[%s] target_config pocket has no center — re-detecting", job_id)
            except Exception as exc:
                logger.warning("[%s] Failed to restore pocket from target_config: %s — re-detecting", job_id, exc)

        if not _pocket_reused:
            from pipeline.pockets import detect_pockets
            pockets = detect_pockets(pdb_path, work_dir, ligand_id=ligand_id)
            best_pocket = pockets[0]
            center = tuple(best_pocket["center"])
            pocket_method = best_pocket.get("method", "unknown")

        _update_progress(
            job_id, 15,
            f"Binding site: ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})",
            step_details=f"{'Reused from Target Setup' if _pocket_reused else str(len(pockets)) + ' pocket(s) detected'} via {pocket_method}",
        )
        pipeline_summary["n_pockets"] = len(pockets)
        pipeline_summary["best_pocket_center"] = list(center)
        pipeline_summary["best_pocket_method"] = pocket_method
        pipeline_summary["steps_completed"].append("pockets")
        audit.log("pockets", f"Pocket ready via {pocket_method}", {
            "n_pockets": len(pockets),
            "best_center": list(center),
            "reused_from_target_setup": _pocket_reused,
        })

        # V6: Check if pocket center falls within an IDR
        if disorder_info and disorder_info.get("idr_regions"):
            try:
                # The pocket center is in 3D space; we approximate the residue
                # index from the PDB. For a rough check, we use the z-coordinate
                # divided by ~3.8 A (rise per residue in extended chain) as a
                # proxy for residue index. This is a heuristic -- real mapping
                # would require parsing ATOM records.
                # A simpler approach: check if the fraction_disordered is high
                # and the pocket overlaps with any IDR based on residue numbering.
                # Since we do not have a direct residue->3D mapping here, we check
                # if ANY detected pocket residues (if available) fall within IDRs.
                pocket_residues = best_pocket.get("residues", [])
                if pocket_residues:
                    idr_set: set[int] = set()
                    for idr_start, idr_end in disorder_info["idr_regions"]:
                        idr_set.update(range(idr_start, idr_end))
                    overlap = [r for r in pocket_residues if r in idr_set]
                    if overlap:
                        pipeline_summary["disorder_warning"] = True
                        pipeline_summary["disorder_warning_detail"] = (
                            f"Pocket overlaps with IDR: {len(overlap)} residue(s) "
                            f"in disordered region(s)"
                        )
                        logger.warning(
                            "[%s] Pocket center overlaps with IDR: %d residue(s) in disordered region",
                            job_id, len(overlap),
                        )
                        audit.log("disorder_warning", "Pocket overlaps with IDR", {
                            "overlapping_residues": len(overlap),
                        })
                elif disorder_info.get("fraction_disordered", 0) > 0.5:
                    # If no residue info but protein is highly disordered, warn anyway
                    pipeline_summary["disorder_warning"] = True
                    pipeline_summary["disorder_warning_detail"] = (
                        f"Protein is highly disordered ({disorder_info['fraction_disordered']:.0%}), "
                        f"pocket reliability may be reduced"
                    )
                    logger.warning(
                        "[%s] Highly disordered protein (%.0f%%), pocket may be unreliable",
                        job_id, disorder_info["fraction_disordered"] * 100,
                    )
                    audit.log("disorder_warning", "Highly disordered protein", {
                        "fraction_disordered": disorder_info["fraction_disordered"],
                    })
            except Exception as exc:
                logger.debug("[%s] Disorder-pocket overlap check failed: %s", job_id, exc)

        # ------------------------------------------------------------------
        # Step 3: Prepare receptor (15% -> 20%)
        # ------------------------------------------------------------------
        _update_progress(job_id, 17, "Preparing receptor (PDB -> PDBQT)")
        from pipeline.prepare import prepare_receptor
        receptor_pdbqt = prepare_receptor(pdb_path, work_dir)
        _update_progress(job_id, 20, "Receptor prepared")
        pipeline_summary["steps_completed"].append("prepare")
        audit.log("prepare", "Receptor prepared (PDB -> PDBQT)", {
            "receptor_path": str(receptor_pdbqt),
        })

        # ------------------------------------------------------------------
        # Step 4: Fetch ligands (20% -> 30%) — V3: auto strategy
        # ------------------------------------------------------------------
        _update_progress(job_id, 22, "Fetching ligands")
        all_ligands: list[dict] = []
        strategy_message: Optional[str] = None

        # V3: Auto ligand strategy for standard mode
        if auto_strategy and uniprot_id and mode != "deep":
            _update_progress(job_id, 23, "Analyzing target in ChEMBL")
            from pipeline.ligands import auto_select_ligand_strategy
            strategy = auto_select_ligand_strategy(uniprot_id)
            strategy_message = strategy["message"]
            update_job(job_id, strategy_message=strategy_message)

            # Override flags based on strategy
            use_chembl = "chembl" in strategy["source"]
            use_zinc = strategy["use_zinc"] or "zinc" in strategy["source"]
            # V3: always enable generation in standard mode
            # auto_strategy only adjusts n_generated based on ChEMBL coverage
            if master:
                enable_generation = True
                if not strategy["use_generation"]:
                    # Well-documented target: generate fewer molecules (exploratory)
                    n_generated_molecules = min(n_generated_molecules, 20)
                    logger.info("Well-documented target: generation enabled with %d molecules", n_generated_molecules)

            pipeline_summary["ligand_strategy"] = strategy
            _update_progress(
                job_id, 24, "Ligand strategy determined",
                step_details=strategy_message,
            )

        if smiles_list:
            from pipeline.ligands import parse_user_smiles
            user_ligs = parse_user_smiles(smiles_list)
            all_ligands.extend(user_ligs)
            logger.info("Added %d user ligands", len(user_ligs))

        if use_chembl and uniprot_id:
            _update_progress(job_id, 25, "Querying ChEMBL database")
            from pipeline.ligands import fetch_chembl_ligands
            chembl_ligs = fetch_chembl_ligands(uniprot_id, max_count=max_ligands)
            all_ligands.extend(chembl_ligs)
            _update_progress(
                job_id, 27, "ChEMBL query complete",
                step_details=f"{len(chembl_ligs)} molecules found in ChEMBL",
            )
            logger.info("Added %d ChEMBL ligands", len(chembl_ligs))

        # V9: PubChem ligands (when gene name is available from strategy)
        if auto_strategy and strategy_message and uniprot_id:
            try:
                strategy_data = pipeline_summary.get("ligand_strategy", {})
                gene_name = strategy_data.get("gene_name")
                if gene_name and strategy_data.get("use_pubchem", False):
                    _update_progress(job_id, 27, f"Querying PubChem for {gene_name}")
                    from pipeline.ligands import fetch_pubchem_ligands
                    pubchem_ligs = fetch_pubchem_ligands(gene_name, max_count=min(30, max_ligands))
                    all_ligands.extend(pubchem_ligs)
                    logger.info("Added %d PubChem ligands for %s", len(pubchem_ligs), gene_name)
                    audit.log("ligands_pubchem", f"PubChem: {len(pubchem_ligs)} ligands for {gene_name}", {
                        "n_pubchem": len(pubchem_ligs),
                        "gene_name": gene_name,
                    })
            except Exception as exc:
                logger.warning("[%s] PubChem ligand fetch failed: %s", job_id, exc)

        if use_zinc:
            _update_progress(job_id, 28, "Loading ZINC molecules")
            from pipeline.ligands import fetch_zinc_ligands
            zinc_ligs = fetch_zinc_ligands(max_count=min(20, max_ligands))
            all_ligands.extend(zinc_ligs)
            logger.info("Added %d ZINC ligands", len(zinc_ligs))

        # V3: Auto-enable ZINC as fallback when sequence-only input has no ligands
        if not all_ligands and not uniprot_id:
            logger.info("No ligands from ChEMBL (no UniProt ID). Auto-enabling ZINC fallback.")
            _update_progress(job_id, 28, "Loading ZINC molecules (auto-fallback)")
            from pipeline.ligands import fetch_zinc_ligands
            zinc_ligs = fetch_zinc_ligands(max_count=min(50, max_ligands))
            all_ligands.extend(zinc_ligs)
            logger.info("Added %d ZINC ligands (auto-fallback)", len(zinc_ligs))

        # Trim to max_ligands
        if len(all_ligands) > max_ligands:
            all_ligands = all_ligands[:max_ligands]

        if not all_ligands:
            raise RuntimeError(
                f"No known ligands found for target {uniprot_id or 'sequence'} in any database. "
                "This is not a docking engine issue. "
                "Try providing custom SMILES molecules, or use a different ligand source."
            )

        _update_progress(
            job_id, 30, f"Collected {len(all_ligands)} ligands",
            step_details=f"{len(all_ligands)} molecules ready for docking",
        )
        pipeline_summary["n_ligands_total"] = len(all_ligands)
        pipeline_summary["steps_completed"].append("ligands")
        _store_pipeline_summary(job_id, pipeline_summary)
        audit.log("ligands", f"Collected {len(all_ligands)} ligands", {
            "n_ligands": len(all_ligands),
        })

        # ------------------------------------------------------------------
        # Step 5: Dock known ligands (30% -> 55%)
        # ------------------------------------------------------------------
        def progress_cb_dock(pct: int, msg: str) -> None:
            # Map 40-90% from V1 to 30-55% in V2/V3
            mapped = int(30 + (pct - 40) / 50 * 25) if pct >= 40 else 30
            _update_progress(job_id, min(mapped, 55), msg)

        if enable_diffdock and master:
            _update_progress(job_id, 31, "Docking with DiffDock (AI mode)")
            from pipeline.docking_diffdock import dock_all_diffdock
            docking_results = dock_all_diffdock(
                protein_pdb=pdb_path,
                ligands=all_ligands,
                work_dir=work_dir,
                progress_callback=progress_cb_dock,
            )
            for r in docking_results:
                r["docking_method"] = "diffdock"
        else:
            engine_label = {"auto": "Auto (GPU/GNINA/Vina)", "gnina": "GNINA", "gnina_gpu": "GNINA GPU (RunPod)", "vina": "Vina"}.get(docking_engine, docking_engine)
            _update_progress(job_id, 31, f"Docking with {engine_label}")
            from pipeline.docking import dock_all_ligands
            docking_results = dock_all_ligands(
                receptor_pdbqt=receptor_pdbqt,
                ligands=all_ligands,
                center=center,
                work_dir=work_dir,
                progress_callback=progress_cb_dock,
                docking_engine=docking_engine,
            )
            for r in docking_results:
                # Use actual engine from dock result, fall back to param
                r["docking_method"] = r.get("docking_engine", docking_engine)

        # Detect actual engine used (from first result)
        actual_engine = docking_engine
        if docking_results:
            actual_engine = docking_results[0].get("docking_engine", docking_engine)

        _update_progress(
            job_id, 55, f"Docking complete: {len(docking_results)} results",
            step_details=f"{len(docking_results)} poses computed ({actual_engine})",
        )
        pipeline_summary["n_docking_results"] = len(docking_results)
        pipeline_summary["docking_engine"] = actual_engine
        pipeline_summary["steps_completed"].append("docking")
        audit.log("docking", f"Docking complete: {len(docking_results)} results", {
            "n_results": len(docking_results),
            "engine": docking_engine,
        })

        # ------------------------------------------------------------------
        # Step 6: Generate new molecules (55% -> 68%) --- standard mode
        # ------------------------------------------------------------------
        generated_results: list[dict] = []

        if False and enable_dmpk and enable_generation:  # generation disabled in main pipeline; only in optimization
            _update_progress(job_id, 56, "Generating novel molecules (AI)")
            from pipeline.generation import generate_molecules

            # A2: Extract top 5 hit SMILES from docking results to seed generation
            top_hit_smiles = [r["smiles"] for r in sorted(docking_results, key=lambda x: x.get("affinity", 0))[:5] if r.get("smiles")]

            gen_mols = generate_molecules(
                pocket_center=center,
                pocket_size=pocket_size,
                receptor_pdbqt=receptor_pdbqt,
                work_dir=work_dir,
                n_molecules=n_generated_molecules,
                n_top=20,
                seed_smiles=top_hit_smiles,
            )
            logger.info("Generated %d novel molecules", len(gen_mols))

            if gen_mols:
                _update_progress(
                    job_id, 62, f"Docking {len(gen_mols)} generated molecules",
                    step_details=f"{len(gen_mols)} molecules generated by AI",
                )
                # Dock the generated molecules
                gen_ligands = [
                    {"name": m["name"], "smiles": m["smiles"], "source": m["source"]}
                    for m in gen_mols
                ]

                if enable_diffdock:
                    from pipeline.docking_diffdock import dock_all_diffdock
                    gen_dock_results = dock_all_diffdock(
                        protein_pdb=pdb_path,
                        ligands=gen_ligands,
                        work_dir=work_dir / "generated",
                    )
                    for r in gen_dock_results:
                        r["docking_method"] = "diffdock"
                else:
                    from pipeline.docking import dock_all_ligands
                    gen_dock_results = dock_all_ligands(
                        receptor_pdbqt=receptor_pdbqt,
                        ligands=gen_ligands,
                        center=center,
                        work_dir=work_dir / "generated",
                        docking_engine=docking_engine,
                    )
                    for r in gen_dock_results:
                        r["docking_method"] = r.get("docking_engine", docking_engine)

                # Merge generation metadata
                gen_meta_by_name = {m["name"]: m for m in gen_mols}
                for r in gen_dock_results:
                    meta = gen_meta_by_name.get(r.get("name", ""), {})
                    r["novelty_score"] = meta.get("novelty_score", 0.5)
                    r["source"] = meta.get("source", "mock_generation")

                generated_results = gen_dock_results

            _update_progress(
                job_id, 68, f"Generated {len(generated_results)} docked molecules",
                step_details=f"{len(generated_results)} molecules generated and docked",
            )
            pipeline_summary["n_generated"] = len(generated_results)
            pipeline_summary["steps_completed"].append("generation")
            audit.log("generation", f"Generated {len(generated_results)} molecules", {
                "n_generated": len(generated_results),
            })
        else:
            _update_progress(job_id, 68, "Skipping molecule generation")

        # ------------------------------------------------------------------
        # Step 7/8: ADMET predictions (68% -> 78%) — standard mode
        # ------------------------------------------------------------------
        admet_results: list[dict] = []
        all_for_admet = docking_results + generated_results

        if enable_admet:
            _update_progress(job_id, 70, "Predicting ADMET properties")
            from pipeline.admet import predict_admet

            all_smiles = list({r.get("smiles", "") for r in all_for_admet if r.get("smiles")})
            if all_smiles:
                admet_results = predict_admet(all_smiles)
                logger.info("ADMET predictions for %d molecules", len(admet_results))

            _update_progress(
                job_id, 78, f"ADMET complete for {len(admet_results)} molecules",
                step_details=f"Toxicity evaluated for {len(admet_results)} molecules",
            )
            pipeline_summary["n_admet"] = len(admet_results)
            pipeline_summary["steps_completed"].append("admet")
            audit.log("admet", f"ADMET predictions for {len(admet_results)} molecules", {
                "n_admet": len(admet_results),
            })
        else:
            _update_progress(job_id, 78, "Skipping ADMET")

        # ------------------------------------------------------------------
        # Step 9: Score and rank (78% -> 82%)
        # ------------------------------------------------------------------
        _update_progress(job_id, 80, "Computing scores")

        if enable_admet and admet_results:
            from pipeline.scoring import score_results_v2
            scored_known = score_results_v2(docking_results, admet_results)
            scored_generated = score_results_v2(generated_results, admet_results) if generated_results else []
        else:
            from pipeline.scoring import score_results
            scored_known = score_results(docking_results)
            scored_generated = score_results(generated_results) if generated_results else []

        _update_progress(job_id, 82, "Scoring complete")
        pipeline_summary["steps_completed"].append("scoring")
        audit.log("scoring", "Scoring and ranking complete", {
            "n_known_scored": len(scored_known),
            "n_generated_scored": len(scored_generated),
        })

        # ------------------------------------------------------------------
        # Step 9a-bis: V6.1 Consensus detail enrichment (z-scores, agreement)
        # ------------------------------------------------------------------
        try:
            from pipeline.scoring import enrich_consensus_detail
            if scored_known:
                enrich_consensus_detail(scored_known)
            if scored_generated:
                enrich_consensus_detail(scored_generated)
            pipeline_summary["steps_completed"].append("consensus_detail")
            audit.log("consensus_detail", "Consensus detail enrichment complete", {
                "n_known": len(scored_known),
                "n_generated": len(scored_generated),
            })
        except Exception as exc:
            logger.warning("[%s] Consensus detail enrichment failed: %s", job_id, exc)

        # ------------------------------------------------------------------
        # Step 9b: V5bis Hard cutoffs (82% -> 83%)
        # ------------------------------------------------------------------
        eliminated_known: list[dict] = []
        eliminated_generated: list[dict] = []
        try:
            _update_progress(job_id, 82, "Applying hard cutoffs")
            from pipeline.scoring import apply_hard_cutoffs
            scored_known, eliminated_known = apply_hard_cutoffs(scored_known)
            if scored_generated:
                scored_generated, eliminated_generated = apply_hard_cutoffs(scored_generated)
            total_eliminated = len(eliminated_known) + len(eliminated_generated)
            total_passed = len(scored_known) + len(scored_generated)
            # Count elimination reasons for display
            from collections import Counter
            all_eliminated = eliminated_known + eliminated_generated
            reason_counts = dict(Counter(
                mol.get("elimination_reason", "unknown") for mol in all_eliminated
            ))
            pipeline_summary["hard_cutoffs"] = {
                "passed": total_passed,
                "eliminated": total_eliminated,
                "reasons": reason_counts,
            }
            pipeline_summary["steps_completed"].append("hard_cutoffs")
            audit.log("hard_cutoffs", f"{total_passed} passed, {total_eliminated} eliminated", {
                "passed": total_passed,
                "eliminated": total_eliminated,
            })
            logger.info("[%s] Hard cutoffs: %d passed, %d eliminated", job_id, total_passed, total_eliminated)
        except Exception as exc:
            logger.warning("[%s] Hard cutoffs failed: %s", job_id, exc)

        # ------------------------------------------------------------------
        # Step 9c: V5bis Interaction analysis (83% -> 84%)
        # ------------------------------------------------------------------
        try:
            _update_progress(job_id, 83, "Analyzing protein-ligand interactions")
            all_scored = sorted(
                scored_known + scored_generated,
                key=lambda x: x.get("composite_score", 0), reverse=True,
            )
            top_for_interactions = all_scored[:10]
            n_analyzed = _safe_interaction_analysis(
                job_id, top_for_interactions, str(pdb_path), uniprot_id
            )
            if n_analyzed > 0:
                pipeline_summary["steps_completed"].append("interactions")
            audit.log("interactions", f"Interaction analysis for {n_analyzed}/{len(top_for_interactions)} candidates", {
                "n_analyzed": n_analyzed,
            })
        except Exception as exc:
            logger.warning("[%s] Interaction analysis failed: %s", job_id, exc)

        # ------------------------------------------------------------------
        # Step 9d: V5bis Butina clustering (84% -> 84%)
        # ------------------------------------------------------------------
        try:
            from pipeline.scoring import cluster_results
            all_for_clustering = scored_known + scored_generated
            cluster_results(all_for_clustering)
            n_clusters = len(set(m.get("cluster_id", 0) for m in all_for_clustering))
            pipeline_summary["n_clusters"] = n_clusters
            pipeline_summary["steps_completed"].append("clustering")
            audit.log("clustering", f"{len(all_for_clustering)} molecules -> {n_clusters} chemical families", {
                "n_clusters": n_clusters,
            })
            logger.info("[%s] Butina clustering: %d families", job_id, n_clusters)
        except Exception as exc:
            logger.warning("[%s] Butina clustering failed: %s", job_id, exc)

        # ------------------------------------------------------------------
        # Step 9e: V6.2 Pareto multi-objective ranking
        # ------------------------------------------------------------------
        try:
            from pipeline.scoring import pareto_ranking
            all_for_pareto = scored_known + scored_generated
            if all_for_pareto:
                pareto_ranking(all_for_pareto)
                front_size = sum(1 for m in all_for_pareto if m.get("pareto_front", False))
                pipeline_summary["pareto_front_size"] = front_size
                pipeline_summary["steps_completed"].append("pareto_ranking")
                audit.log("pareto_ranking", f"Pareto ranking: {front_size} molecules on front", {
                    "n_molecules": len(all_for_pareto),
                    "pareto_front_size": front_size,
                })
                logger.info("[%s] Pareto ranking: %d on front out of %d", job_id, front_size, len(all_for_pareto))
                # Re-split into scored_known and scored_generated after sorting
                # The pareto_ranking function sorts the combined list; we need to
                # preserve the known vs generated distinction for downstream steps.
                known_smiles_set = {m.get("smiles", "") for m in docking_results}
                scored_known = [m for m in all_for_pareto if m.get("smiles", "") in known_smiles_set]
                scored_generated = [m for m in all_for_pareto if m.get("smiles", "") not in known_smiles_set]
        except Exception as exc:
            logger.warning("[%s] Pareto ranking failed: %s", job_id, exc)

        # ------------------------------------------------------------------
        # Step 10: Off-target screening (84% -> 86%) --- standard mode, top 5
        # ------------------------------------------------------------------
        if enable_selectivity:
            _update_progress(job_id, 84, "Off-target safety screening")
            try:
                from pipeline.off_target import screen_candidates
                all_scored = sorted(
                    scored_known + scored_generated,
                    key=lambda x: x.get("composite_score", 0),
                    reverse=True,
                )
                top5_for_ot = all_scored[:5]
                top5_screened = screen_candidates(top5_for_ot, work_dir)
                # Merge off-target results back into the original lists
                ot_by_smiles = {
                    m.get("smiles", ""): m.get("off_target_results")
                    for m in top5_screened
                    if m.get("off_target_results")
                }
                for mol in scored_known + scored_generated:
                    smi = mol.get("smiles", "")
                    if smi in ot_by_smiles:
                        mol["off_target_results"] = ot_by_smiles[smi]
                pipeline_summary["steps_completed"].append("off_target")
                audit.log("off_target", f"Off-target screening for top 5 candidates", {
                    "n_screened": len(top5_screened),
                })
            except Exception as exc:
                logger.warning("[%s] Off-target screening failed: %s", job_id, exc)

        # ------------------------------------------------------------------
        # Step 10a (V6.3): Combined off-target screening for top 5
        # ------------------------------------------------------------------
        if enable_selectivity:
            try:
                _update_progress(job_id, 85, "V6.3: Combined off-target screening (SEA + docking)")
                from pipeline.off_target import combined_off_target_screening
                all_scored_for_combined = sorted(
                    scored_known + scored_generated,
                    key=lambda x: x.get("composite_score", 0),
                    reverse=True,
                )
                top5_for_combined = all_scored_for_combined[:5]
                combined_ot_summary: list[dict] = []
                for mol in top5_for_combined:
                    smi = mol.get("smiles", "")
                    if not smi:
                        continue
                    combined_result = combined_off_target_screening(smi, work_dir)
                    mol["combined_off_target"] = combined_result
                    combined_ot_summary.append({
                        "name": mol.get("name", "unknown"),
                        "combined_selectivity": combined_result.get("combined_selectivity", 0.0),
                        "tier1_hits": len(combined_result.get("tier1_hits", [])),
                        "tier2_safe": combined_result.get("tier2_safe_count", 0),
                    })
                # Merge back into scored lists by SMILES
                combined_by_smiles = {
                    m.get("smiles", ""): m.get("combined_off_target")
                    for m in top5_for_combined
                    if m.get("combined_off_target")
                }
                for mol in scored_known + scored_generated:
                    smi = mol.get("smiles", "")
                    if smi in combined_by_smiles:
                        mol["combined_off_target"] = combined_by_smiles[smi]
                pipeline_summary["combined_off_target_summary"] = combined_ot_summary
                pipeline_summary["steps_completed"].append("combined_off_target")
                audit.log("combined_off_target", "V6.3 combined off-target screening for top 5", {
                    "n_screened": len(top5_for_combined),
                    "summary": combined_ot_summary,
                })
                logger.info("[%s] V6.3 combined off-target: screened %d molecules", job_id, len(top5_for_combined))
            except Exception as exc:
                logger.warning("[%s] V6.3 combined off-target screening failed: %s", job_id, exc)

        # ------------------------------------------------------------------
        # Step 10a2 (V6.3): Specialized hERG for top 20 molecules
        # ------------------------------------------------------------------
        if enable_herg:
            try:
                _update_progress(job_id, 86, "V6.3: Specialized hERG screening (top 20)")
                from pipeline.admet import predict_herg_specialized
                all_scored_for_herg = sorted(
                    scored_known + scored_generated,
                    key=lambda x: x.get("composite_score", 0),
                    reverse=True,
                )
                top20_for_herg = all_scored_for_herg[:20]
                herg_summary_counts = {"LOW": 0, "MODERATE": 0, "HIGH": 0}
                for mol in top20_for_herg:
                    smi = mol.get("smiles", "")
                    if not smi:
                        continue
                    herg_result = predict_herg_specialized(smi)
                    mol["herg_specialized"] = herg_result
                    herg_summary_counts[herg_result.get("risk_level", "LOW")] += 1
                # Merge back by SMILES
                herg_by_smiles = {
                    m.get("smiles", ""): m.get("herg_specialized")
                    for m in top20_for_herg
                    if m.get("herg_specialized")
                }
                for mol in scored_known + scored_generated:
                    smi = mol.get("smiles", "")
                    if smi in herg_by_smiles:
                        mol["herg_specialized"] = herg_by_smiles[smi]
                pipeline_summary["herg_specialized_summary"] = herg_summary_counts
                pipeline_summary["steps_completed"].append("herg_specialized")
                audit.log("herg_specialized", "V6.3 specialized hERG for top 20", {
                    "n_screened": len(top20_for_herg),
                    "risk_counts": herg_summary_counts,
                })
                logger.info(
                    "[%s] V6.3 specialized hERG: %d screened (LOW=%d, MOD=%d, HIGH=%d)",
                    job_id, len(top20_for_herg),
                    herg_summary_counts["LOW"],
                    herg_summary_counts["MODERATE"],
                    herg_summary_counts["HIGH"],
                )
            except Exception as exc:
                logger.warning("[%s] V6.3 specialized hERG screening failed: %s", job_id, exc)

        # ------------------------------------------------------------------
        # Step 10b: Retrosynthesis of top 5 (86% -> 92%) --- standard mode
        # ------------------------------------------------------------------
        if enable_synthesis and enable_retrosynthesis:
            _update_progress(job_id, 87, "Planning retrosynthesis for top molecules")
            from pipeline.retrosynthesis import plan_synthesis

            # Take top 5 from both known and generated
            top_for_synthesis = (scored_known[:3] + scored_generated[:2])[:5]

            for i, mol in enumerate(top_for_synthesis):
                smiles = mol.get("smiles", "")
                name = mol.get("name", f"mol_{i}")
                if not smiles:
                    continue
                _update_progress(
                    job_id, 87 + i, f"Retrosynthesis: {name}",
                    step_details=f"Planning synthesis of {name}",
                )
                try:
                    route = plan_synthesis(smiles, max_depth=6, timeout_sec=120)
                    mol["synthesis_route"] = route
                except Exception as exc:
                    logger.warning("Retrosynthesis failed for %s: %s", name, exc)
                    mol["synthesis_route"] = None

            _update_progress(job_id, 92, "Retrosynthesis complete")
            pipeline_summary["steps_completed"].append("retrosynthesis")
            audit.log("retrosynthesis", "Retrosynthesis planning complete", {
                "n_planned": len(top_for_synthesis),
            })
        else:
            _update_progress(job_id, 92, "Skipping retrosynthesis")

        # ------------------------------------------------------------------
        # Step 10c: Confidence scoring (all results) — guarded by enable_safety
        # ------------------------------------------------------------------
        if enable_safety:
            try:
                from pipeline.confidence import calculate_confidence

                for mol in scored_known + scored_generated:
                    # V5bis: attach pocket method for confidence scoring
                    mol.setdefault("pocket_method", pocket_method)
                    # V6: pass disorder_info for potential penalty
                    mol["confidence"] = calculate_confidence(
                        mol, effective_structure_source, disorder_info=disorder_info
                    )
                pipeline_summary["steps_completed"].append("confidence")
                audit.log("confidence", "Confidence scores calculated for all results", {
                    "n_scored": len(scored_known) + len(scored_generated),
                    "effective_structure_source": effective_structure_source,
                    "disorder_penalty_applied": bool(
                        disorder_info and disorder_info.get("fraction_disordered", 0) > 0.3
                    ),
                })
            except Exception as exc:
                logger.warning("[%s] Confidence scoring failed: %s", job_id, exc)
        else:
            _update_progress(job_id, 93, "Skipping confidence scoring")

        # ------------------------------------------------------------------
        # Step 11: Generate report + ZIP (92% -> 98%)
        # ------------------------------------------------------------------
        _update_progress(job_id, 93, "Generating PDF report")
        # Resolve protein name from UniProt API (fallback to uniprot_id)
        protein_name: str = "Sequence-based target"
        if uniprot_id:
            protein_name = _fetch_protein_name(uniprot_id) or uniprot_id
        update_job(job_id, protein_name=protein_name)
        job_meta = {
            "job_id": job_id,
            "uniprot_id": uniprot_id or "N/A",
            "protein_name": protein_name,
            "mode": mode,
            "structure_source": structure_source,
            "strategy_message": strategy_message,
            # V5bis metadata
            "pdb_info": pdb_info,
            "pocket_method": pocket_method,
            "n_eliminated": len(eliminated_known) + len(eliminated_generated),
            "n_clusters": pipeline_summary.get("n_clusters", 0),
        }

        from pipeline.report import generate_csv, generate_pdf_report, generate_zip_archive

        all_results = scored_known + scored_generated

        # CSV
        csv_path = generate_csv(all_results, work_dir / "results.csv")

        # PDF
        pdf_path = generate_pdf_report(
            job_meta, scored_known, work_dir / "report.pdf",
            generated_results=scored_generated if master else None,
        )

        _update_progress(job_id, 97, "Creating ZIP archive")
        zip_path = generate_zip_archive(work_dir, work_dir / "results.zip")
        pipeline_summary["steps_completed"].append("report")
        audit.log("report", "Report and ZIP archive generated", {
            "pdf_path": str(pdf_path),
            "zip_path": str(zip_path),
        })

        # V5: Save audit log
        try:
            audit.log("pipeline", "Pipeline completed successfully")
            audit.save(work_dir / "audit_log.json")
        except Exception as exc:
            logger.warning("[%s] Failed to save audit log: %s", job_id, exc)

        # ------------------------------------------------------------------
        # Final: Persist results and mark complete (100%)
        # ------------------------------------------------------------------
        # V5bis: Include eliminated molecules at the end of results
        # (passed first, eliminated last) so the frontend can display them with red badges
        all_known = scored_known + eliminated_known
        all_generated = scored_generated + eliminated_generated

        serializable_known = _make_serializable(all_known)
        serializable_generated = _make_serializable(all_generated)

        results_json = json.dumps(serializable_known)
        generated_json = json.dumps(serializable_generated) if serializable_generated else None

        # Finalize pipeline summary
        pipeline_summary["total_results"] = len(all_known) + len(all_generated)
        pipeline_summary["n_passed"] = len(scored_known) + len(scored_generated)
        pipeline_summary["mode"] = mode
        if scored_known:
            pipeline_summary["best_affinity"] = min(r.get("affinity", 0.0) for r in scored_known)
            pipeline_summary["best_composite"] = max(r.get("composite_score", 0.0) for r in scored_known)
        _store_pipeline_summary(job_id, pipeline_summary)

        update_job(
            job_id,
            status="completed",
            progress=100,
            current_step="Complete",
            results_json=results_json,
            generated_json=generated_json,
            report_path=str(pdf_path),
            zip_path=str(zip_path),
            completed_at=datetime.now(timezone.utc),
        )
        total_results = len(scored_known) + len(scored_generated)
        logger.info("[%s] Pipeline completed (%s mode): %d results", job_id, mode, total_results)

        # V3: Send notification email if configured
        if notification_email:
            try:
                from notifications import send_notification_email
                send_notification_email(notification_email, job_id)
            except Exception as exc:
                logger.warning("[%s] Failed to send notification email: %s", job_id, exc)

        return {
            "status": "completed",
            "result_count": total_results,
            "mode": mode,
        }

    except Exception as exc:
        error_msg = f"{type(exc).__name__}: {exc}"
        logger.error("[%s] Pipeline failed: %s\n%s", job_id, error_msg, traceback.format_exc())
        update_job(
            job_id,
            status="failed",
            current_step="Failed",
            error_message=error_msg,
            completed_at=datetime.now(timezone.utc),
        )
        return {
            "status": "failed",
            "error": error_msg,
        }


@celery_app.task(name="run_deep_screening", bind=True, max_retries=0, time_limit=14400, soft_time_limit=14000)
def run_deep_screening(self, job_id: str, params: dict) -> dict:
    """Execute a deep (massive) screening pipeline, up to 4 hours.

    This task is dispatched for mode="deep". It orchestrates:
      structure -> pockets -> prepare -> auto_strategy -> screening (5 passes) ->
      ADMET -> retrosynthesis -> report

    Parameters
    ----------
    job_id : str
        UUID of the job.
    params : dict
        Same keys as run_pipeline params.

    Returns
    -------
    dict
        Summary with ``status`` and ``result_count``.
    """
    work_dir = DATA_DIR / job_id
    work_dir.mkdir(parents=True, exist_ok=True)

    pipeline_summary: dict = {
        "steps_completed": [],
        "mode": "deep",
        "timings": {},
        "counts": {},
    }

    try:
        update_job(job_id, status="running", progress=0, current_step="Starting deep screening pipeline")

        uniprot_id: str = params.get("uniprot_id", "") or ""
        sequence: Optional[str] = params.get("sequence")
        max_ligands: int = params.get("max_ligands", 5000)
        smiles_list: list[str] | None = params.get("smiles_list")
        docking_engine: str = params.get("docking_engine", "auto")
        notification_email: Optional[str] = params.get("notification_email")
        n_generated_molecules: int = params.get("n_generated_molecules", 200)

        # ------------------------------------------------------------------
        # Step 1: Fetch structure (0% -> 5%)
        # ------------------------------------------------------------------
        _update_progress(job_id, 2, "Deep screening: Fetching protein structure")

        structure_source = "unknown"
        if sequence and not uniprot_id:
            from pipeline.structure import fetch_structure_from_sequence
            pdb_path, structure_source = fetch_structure_from_sequence(sequence, work_dir)
        else:
            from pipeline.structure import fetch_structure
            pdb_path, structure_source = fetch_structure(uniprot_id, work_dir)

        update_job(job_id, pdb_path=str(pdb_path), structure_source=structure_source)
        pipeline_summary["structure_source"] = structure_source
        pipeline_summary["steps_completed"].append("structure")
        _update_progress(job_id, 5, "Structure retrieved", step_details=f"Source: {structure_source}")

        # V5bis: Retrieve PDB info for ligand-aware pocket detection
        pdb_info_deep: Optional[dict] = None
        ligand_id_deep: Optional[str] = None
        if uniprot_id and structure_source == "pdb_experimental":
            try:
                from pipeline.structure import get_pdb_info
                pdb_info_deep = get_pdb_info(work_dir, uniprot_id)
                if pdb_info_deep:
                    ligand_id_deep = pdb_info_deep.get("ligand_id")
                    pipeline_summary["pdb_info"] = pdb_info_deep
            except Exception as exc:
                logger.warning("[%s] Failed to retrieve PDB info in deep: %s", job_id, exc)

        # --- A1: Compute effective_structure_source early (deep) ---
        effective_structure_source_deep = structure_source
        if structure_source == "pdb_experimental" and pdb_info_deep and pdb_info_deep.get("ligand_id"):
            effective_structure_source_deep = "pdb_holo"
        pipeline_summary["effective_structure_source"] = effective_structure_source_deep

        # ------------------------------------------------------------------
        # Step 2: Detect pockets (5% -> 8%)
        # ------------------------------------------------------------------
        _update_progress(job_id, 6, "Deep screening: Detecting binding pockets")
        from pipeline.pockets import detect_pockets
        pockets = detect_pockets(pdb_path, work_dir, ligand_id=ligand_id_deep)
        best_pocket = pockets[0]
        center = tuple(best_pocket["center"])
        pocket_size = (25.0, 25.0, 25.0)
        pocket_method_deep = best_pocket.get("method", "unknown")
        pipeline_summary["n_pockets"] = len(pockets)
        pipeline_summary["best_pocket_method"] = pocket_method_deep
        pipeline_summary["steps_completed"].append("pockets")
        _update_progress(job_id, 8, "Pockets detected", step_details=f"{len(pockets)} pocket(s) via {pocket_method_deep}")

        # ------------------------------------------------------------------
        # Step 3: Prepare receptor (8% -> 10%)
        # ------------------------------------------------------------------
        _update_progress(job_id, 9, "Deep screening: Preparing receptor")
        from pipeline.prepare import prepare_receptor
        receptor_pdbqt = prepare_receptor(pdb_path, work_dir)
        pipeline_summary["steps_completed"].append("prepare")
        _update_progress(job_id, 10, "Receptor prepared")

        # ------------------------------------------------------------------
        # Step 4: Auto strategy + massive ligand collection (10% -> 20%)
        # ------------------------------------------------------------------
        _update_progress(job_id, 11, "Deep screening: Analyzing target and collecting ligands")
        all_ligands: list[dict] = []
        strategy_message: Optional[str] = None

        if uniprot_id:
            from pipeline.ligands import auto_select_ligand_strategy
            strategy = auto_select_ligand_strategy(uniprot_id)
            strategy_message = strategy["message"]
            update_job(job_id, strategy_message=strategy_message)
            pipeline_summary["ligand_strategy"] = strategy

        # Collect from all sources for deep screening
        if smiles_list:
            from pipeline.ligands import parse_user_smiles
            user_ligs = parse_user_smiles(smiles_list)
            all_ligands.extend(user_ligs)

        if uniprot_id:
            from pipeline.ligands import fetch_chembl_ligands
            chembl_ligs = fetch_chembl_ligands(uniprot_id, max_count=min(max_ligands, 500))
            all_ligands.extend(chembl_ligs)
            _update_progress(
                job_id, 14, "ChEMBL ligands collected",
                step_details=f"{len(chembl_ligs)} molecules from ChEMBL",
            )

        from pipeline.ligands import fetch_zinc_ligands
        zinc_ligs = fetch_zinc_ligands(max_count=min(100, max_ligands))
        all_ligands.extend(zinc_ligs)

        # V9: PubChem for deep screening
        try:
            from pipeline.ligands import resolve_gene_name, fetch_pubchem_ligands
            gene_name_deep = resolve_gene_name(uniprot_id) if uniprot_id else None
            if gene_name_deep:
                _update_progress(job_id, 16, f"Querying PubChem for {gene_name_deep}")
                pubchem_ligs = fetch_pubchem_ligands(gene_name_deep, max_count=min(50, max_ligands))
                all_ligands.extend(pubchem_ligs)
                logger.info("Deep: added %d PubChem ligands for %s", len(pubchem_ligs), gene_name_deep)
        except Exception as exc:
            logger.warning("[%s] Deep PubChem fetch failed: %s", job_id, exc)

        if not all_ligands:
            raise RuntimeError("No ligands found for deep screening.")

        # Trim to max_ligands
        if len(all_ligands) > max_ligands:
            all_ligands = all_ligands[:max_ligands]

        pipeline_summary["n_ligands_total"] = len(all_ligands)
        pipeline_summary["steps_completed"].append("ligands")
        _update_progress(
            job_id, 20, f"Collected {len(all_ligands)} ligands for deep screening",
            step_details=f"{len(all_ligands)} molecules ready",
        )
        _store_pipeline_summary(job_id, pipeline_summary)

        # ------------------------------------------------------------------
        # Step 5: Massive screening — 5 passes (20% -> 65%)
        # ------------------------------------------------------------------
        _update_progress(job_id, 21, "Deep screening: Starting massive docking (5 passes)")
        from pipeline.docking import dock_all_ligands

        # Split ligands into 5 batches for progress tracking
        batch_size = max(1, len(all_ligands) // 5)
        batches = [
            all_ligands[i:i + batch_size]
            for i in range(0, len(all_ligands), batch_size)
        ]

        all_docking_results: list[dict] = []
        for pass_idx, batch in enumerate(batches):
            pass_num = pass_idx + 1
            base_pct = 20 + (pass_idx * 9)  # 20, 29, 38, 47, 56
            _update_progress(
                job_id, base_pct,
                f"Deep screening: Pass {pass_num}/5 ({len(batch)} molecules)",
                step_details=f"Pass {pass_num}/5 - {len(batch)} molecules",
            )

            batch_results = dock_all_ligands(
                receptor_pdbqt=receptor_pdbqt,
                ligands=batch,
                center=center,
                work_dir=work_dir / f"pass_{pass_num}",
                docking_engine=docking_engine,
            )
            for r in batch_results:
                r["docking_method"] = r.get("docking_engine", docking_engine)
            all_docking_results.extend(batch_results)

        # Detect actual engine used
        if all_docking_results:
            pipeline_summary["docking_engine"] = all_docking_results[0].get("docking_engine", docking_engine)
        pipeline_summary["n_docking_results"] = len(all_docking_results)
        pipeline_summary["steps_completed"].append("massive_docking")
        _update_progress(
            job_id, 65, f"Massive docking complete: {len(all_docking_results)} results",
            step_details=f"{len(all_docking_results)} poses computed in 5 passes",
        )

        # ------------------------------------------------------------------
        # Step 6: AI molecule generation (65% -> 72%)
        # ------------------------------------------------------------------
        generated_results: list[dict] = []
        _update_progress(job_id, 66, "Deep screening: Generating novel molecules (AI)")
        try:
            from pipeline.generation import generate_molecules
            # A2: Extract top 5 hit SMILES from docking results to seed generation
            top_hit_smiles_deep = [r["smiles"] for r in sorted(all_docking_results, key=lambda x: x.get("affinity", 0))[:5] if r.get("smiles")]
            gen_mols = generate_molecules(
                pocket_center=center,
                pocket_size=pocket_size,
                receptor_pdbqt=receptor_pdbqt,
                work_dir=work_dir,
                n_molecules=n_generated_molecules,
                n_top=30,
                seed_smiles=top_hit_smiles_deep,
            )
            if gen_mols:
                gen_ligands = [
                    {"name": m["name"], "smiles": m["smiles"], "source": m["source"]}
                    for m in gen_mols
                ]
                gen_dock_results = dock_all_ligands(
                    receptor_pdbqt=receptor_pdbqt,
                    ligands=gen_ligands,
                    center=center,
                    work_dir=work_dir / "generated",
                    docking_engine=docking_engine,
                )
                gen_meta_by_name = {m["name"]: m for m in gen_mols}
                for r in gen_dock_results:
                    meta = gen_meta_by_name.get(r.get("name", ""), {})
                    r["novelty_score"] = meta.get("novelty_score", 0.5)
                    r["source"] = meta.get("source", "mock_generation")
                    r["docking_method"] = r.get("docking_engine", docking_engine)
                generated_results = gen_dock_results
        except Exception as exc:
            logger.warning("[%s] Generation failed in deep screening: %s", job_id, exc)

        pipeline_summary["n_generated"] = len(generated_results)
        pipeline_summary["steps_completed"].append("generation")
        _update_progress(job_id, 72, f"Generated {len(generated_results)} molecules")

        # ------------------------------------------------------------------
        # Step 7: ADMET (72% -> 80%)
        # ------------------------------------------------------------------
        _update_progress(job_id, 73, "Deep screening: ADMET predictions")
        admet_results: list[dict] = []
        all_for_admet = all_docking_results + generated_results
        try:
            from pipeline.admet import predict_admet
            all_smiles = list({r.get("smiles", "") for r in all_for_admet if r.get("smiles")})
            if all_smiles:
                admet_results = predict_admet(all_smiles)
        except Exception as exc:
            logger.warning("[%s] ADMET failed in deep screening: %s", job_id, exc)

        pipeline_summary["n_admet"] = len(admet_results)
        pipeline_summary["steps_completed"].append("admet")
        _update_progress(job_id, 80, f"ADMET complete for {len(admet_results)} molecules")

        # ------------------------------------------------------------------
        # Step 8: Score and rank (80% -> 83%)
        # ------------------------------------------------------------------
        _update_progress(job_id, 81, "Deep screening: Scoring and ranking")
        if admet_results:
            from pipeline.scoring import score_results_v2
            scored_known = score_results_v2(all_docking_results, admet_results)
            scored_generated = score_results_v2(generated_results, admet_results) if generated_results else []
        else:
            from pipeline.scoring import score_results
            scored_known = score_results(all_docking_results)
            scored_generated = score_results(generated_results) if generated_results else []

        pipeline_summary["steps_completed"].append("scoring")
        _update_progress(job_id, 83, "Scoring complete")

        # ------------------------------------------------------------------
        # Step 8a-bis: V6.1 Consensus detail enrichment (z-scores, agreement)
        # ------------------------------------------------------------------
        try:
            from pipeline.scoring import enrich_consensus_detail
            if scored_known:
                enrich_consensus_detail(scored_known)
            if scored_generated:
                enrich_consensus_detail(scored_generated)
            pipeline_summary["steps_completed"].append("consensus_detail")
        except Exception as exc:
            logger.warning("[%s] Deep consensus detail enrichment failed: %s", job_id, exc)

        # ------------------------------------------------------------------
        # Step 8b: V5bis Hard cutoffs + clustering (83% -> 85%)
        # ------------------------------------------------------------------
        eliminated_known_deep: list[dict] = []
        eliminated_generated_deep: list[dict] = []
        try:
            from pipeline.scoring import apply_hard_cutoffs, cluster_results
            scored_known, eliminated_known_deep = apply_hard_cutoffs(scored_known)
            if scored_generated:
                scored_generated, eliminated_generated_deep = apply_hard_cutoffs(scored_generated)
            cluster_results(scored_known + scored_generated)
            total_elim = len(eliminated_known_deep) + len(eliminated_generated_deep)
            total_pass = len(scored_known) + len(scored_generated)
            from collections import Counter
            reason_counts = dict(Counter(
                mol.get("elimination_reason", "unknown")
                for mol in eliminated_known_deep + eliminated_generated_deep
            ))
            pipeline_summary["hard_cutoffs"] = {
                "passed": total_pass,
                "eliminated": total_elim,
                "reasons": reason_counts,
            }
            pipeline_summary["steps_completed"].append("hard_cutoffs")
            pipeline_summary["steps_completed"].append("clustering")
        except Exception as exc:
            logger.warning("[%s] Deep V5bis cutoffs/clustering failed: %s", job_id, exc)

        _update_progress(job_id, 85, "Cutoffs and clustering complete")

        # ------------------------------------------------------------------
        # Step 8c: V6.2 Pareto multi-objective ranking (deep)
        # ------------------------------------------------------------------
        try:
            from pipeline.scoring import pareto_ranking
            all_for_pareto_deep = scored_known + scored_generated
            if all_for_pareto_deep:
                pareto_ranking(all_for_pareto_deep)
                front_size_deep = sum(1 for m in all_for_pareto_deep if m.get("pareto_front", False))
                pipeline_summary["pareto_front_size"] = front_size_deep
                pipeline_summary["steps_completed"].append("pareto_ranking")
                logger.info("[%s] Deep Pareto ranking: %d on front out of %d", job_id, front_size_deep, len(all_for_pareto_deep))
                # Re-split after Pareto sorting
                known_smiles_set_deep = {m.get("smiles", "") for m in all_docking_results}
                scored_known = [m for m in all_for_pareto_deep if m.get("smiles", "") in known_smiles_set_deep]
                scored_generated = [m for m in all_for_pareto_deep if m.get("smiles", "") not in known_smiles_set_deep]
        except Exception as exc:
            logger.warning("[%s] Deep Pareto ranking failed: %s", job_id, exc)

        # ------------------------------------------------------------------
        # Step 8d: Off-target screening (deep) — top 5
        # ------------------------------------------------------------------
        try:
            _update_progress(job_id, 86, "Deep screening: Off-target safety screening")
            from pipeline.off_target import screen_candidates
            all_scored_deep_ot = sorted(
                scored_known + scored_generated,
                key=lambda x: x.get("composite_score", 0),
                reverse=True,
            )
            top5_deep_ot = all_scored_deep_ot[:5]
            top5_deep_screened = screen_candidates(top5_deep_ot, work_dir)
            ot_by_smiles_deep = {
                m.get("smiles", ""): m.get("off_target_results")
                for m in top5_deep_screened
                if m.get("off_target_results")
            }
            for mol in scored_known + scored_generated:
                smi = mol.get("smiles", "")
                if smi in ot_by_smiles_deep:
                    mol["off_target_results"] = ot_by_smiles_deep[smi]
            pipeline_summary["steps_completed"].append("off_target")
            logger.info("[%s] Deep off-target screening: screened %d molecules", job_id, len(top5_deep_screened))
        except Exception as exc:
            logger.warning("[%s] Deep off-target screening failed: %s", job_id, exc)

        # ------------------------------------------------------------------
        # Step 8e: Confidence scoring (deep) — all results
        # ------------------------------------------------------------------
        try:
            from pipeline.confidence import calculate_confidence
            for mol in scored_known + scored_generated:
                mol.setdefault("pocket_method", pocket_method_deep)
                mol["confidence"] = calculate_confidence(mol, effective_structure_source_deep)
            pipeline_summary["steps_completed"].append("confidence")
            logger.info("[%s] Deep confidence scoring complete", job_id)
        except Exception as exc:
            logger.warning("[%s] Deep confidence scoring failed: %s", job_id, exc)

        # ------------------------------------------------------------------
        # Step 8f: Interaction analysis (deep) — top 10
        # ------------------------------------------------------------------
        try:
            _update_progress(job_id, 87, "Deep screening: Analyzing protein-ligand interactions")
            all_scored_deep_ia = sorted(
                scored_known + scored_generated,
                key=lambda x: x.get("composite_score", 0), reverse=True,
            )
            top_for_interactions_deep = all_scored_deep_ia[:10]
            n_analyzed_deep = _safe_interaction_analysis(
                job_id, top_for_interactions_deep, str(pdb_path), uniprot_id
            )
            if n_analyzed_deep > 0:
                pipeline_summary["steps_completed"].append("interactions")
            logger.info("[%s] Deep interaction analysis: analyzed %d/%d molecules",
                        job_id, n_analyzed_deep, len(top_for_interactions_deep))
        except Exception as exc:
            logger.warning("[%s] Deep interaction analysis failed: %s", job_id, exc)

        # ------------------------------------------------------------------
        # Step 9: Retrosynthesis (85% -> 92%)
        # ------------------------------------------------------------------
        _update_progress(job_id, 86, "Deep screening: Retrosynthesis planning")
        try:
            from pipeline.retrosynthesis import plan_synthesis
            top_for_synthesis = (scored_known[:5] + scored_generated[:3])[:8]
            for i, mol in enumerate(top_for_synthesis):
                smiles = mol.get("smiles", "")
                name = mol.get("name", f"mol_{i}")
                if not smiles:
                    continue
                _update_progress(job_id, 86 + i, f"Retrosynthesis: {name}")
                try:
                    route = plan_synthesis(smiles, max_depth=6, timeout_sec=120)
                    mol["synthesis_route"] = route
                except Exception as exc:
                    logger.warning("Retrosynthesis failed for %s: %s", name, exc)
                    mol["synthesis_route"] = None
        except Exception as exc:
            logger.warning("[%s] Retrosynthesis step failed: %s", job_id, exc)

        pipeline_summary["steps_completed"].append("retrosynthesis")
        _update_progress(job_id, 92, "Retrosynthesis complete")

        # ------------------------------------------------------------------
        # Step 10: Report + ZIP (92% -> 98%)
        # ------------------------------------------------------------------
        _update_progress(job_id, 93, "Deep screening: Generating report")
        # Resolve protein name from UniProt API (fallback to uniprot_id)
        protein_name_deep: str = "Sequence-based target"
        if uniprot_id:
            protein_name_deep = _fetch_protein_name(uniprot_id) or uniprot_id
        update_job(job_id, protein_name=protein_name_deep)
        job_meta = {
            "job_id": job_id,
            "uniprot_id": uniprot_id or "N/A",
            "protein_name": protein_name_deep,
            "mode": "deep",
            "structure_source": structure_source,
            "strategy_message": strategy_message,
        }

        from pipeline.report import generate_csv, generate_pdf_report, generate_zip_archive

        all_results = scored_known + scored_generated
        csv_path = generate_csv(all_results, work_dir / "results.csv")
        pdf_path = generate_pdf_report(
            job_meta, scored_known, work_dir / "report.pdf",
            generated_results=scored_generated,
        )
        _update_progress(job_id, 97, "Creating ZIP archive")
        zip_path = generate_zip_archive(work_dir, work_dir / "results.zip")
        pipeline_summary["steps_completed"].append("report")

        # ------------------------------------------------------------------
        # Final: Persist and mark complete (100%)
        # ------------------------------------------------------------------
        # V5bis: Include eliminated molecules at the end of results
        all_known_deep = scored_known + eliminated_known_deep
        all_generated_deep = scored_generated + eliminated_generated_deep

        serializable_known = _make_serializable(all_known_deep)
        serializable_generated = _make_serializable(all_generated_deep)

        results_json = json.dumps(serializable_known)
        generated_json = json.dumps(serializable_generated) if serializable_generated else None

        pipeline_summary["total_results"] = len(all_known_deep) + len(all_generated_deep)
        pipeline_summary["n_passed"] = len(scored_known) + len(scored_generated)
        if scored_known:
            pipeline_summary["best_affinity"] = min(r.get("affinity", 0.0) for r in scored_known)
            pipeline_summary["best_composite"] = max(r.get("composite_score", 0.0) for r in scored_known)
        _store_pipeline_summary(job_id, pipeline_summary)

        update_job(
            job_id,
            status="completed",
            progress=100,
            current_step="Complete",
            results_json=results_json,
            generated_json=generated_json,
            report_path=str(pdf_path),
            zip_path=str(zip_path),
            completed_at=datetime.now(timezone.utc),
        )
        total_results = len(scored_known) + len(scored_generated)
        logger.info("[%s] Deep screening completed: %d results", job_id, total_results)

        # V3: Send notification email
        if notification_email:
            try:
                from notifications import send_notification_email
                send_notification_email(notification_email, job_id)
            except Exception as exc:
                logger.warning("[%s] Failed to send notification email: %s", job_id, exc)

        return {
            "status": "completed",
            "result_count": total_results,
            "mode": "deep",
        }

    except Exception as exc:
        error_msg = f"{type(exc).__name__}: {exc}"
        logger.error("[%s] Deep screening failed: %s\n%s", job_id, error_msg, traceback.format_exc())
        update_job(
            job_id,
            status="failed",
            current_step="Failed",
            error_message=error_msg,
            completed_at=datetime.now(timezone.utc),
        )
        return {
            "status": "failed",
            "error": error_msg,
        }


@celery_app.task(name="run_lead_optimization", bind=True, max_retries=0)
def run_lead_optimization(self, opt_id: str, params: dict) -> dict:
    """Execute lead optimization as a Celery task.

    Parameters
    ----------
    opt_id : str
        Optimization UUID.
    params : dict
        Keys: ``job_id``, ``optimization_id``, ``smiles``, ``molecule_name``,
        ``weights``, ``n_iterations``, ``variants_per_iter``.

    Returns
    -------
    dict
        Summary with ``status`` and optimization ``result``.
    """
    # Import here to get access to the in-memory results store in main.py
    # For Celery workers running in a separate process, we store results
    # in a JSON file under the job's work directory.
    job_id = params["job_id"]
    work_dir = DATA_DIR / job_id
    opt_work_dir = work_dir / "optimization" / opt_id
    opt_work_dir.mkdir(parents=True, exist_ok=True)

    try:
        from pipeline.lead_optimization import run_optimization

        # Get pocket center and size from pipeline summary
        pocket_center = [22.0, 0.5, 18.0]  # fallback
        pocket_size = params.get("pocket_size", [25.0, 25.0, 25.0])
        try:
            from database import get_job
            job = get_job(job_id)
            if job and job.pipeline_summary_json:
                summary = json.loads(job.pipeline_summary_json)
                pocket_center = summary.get("best_pocket_center", pocket_center)
                # Extract pocket size from pipeline summary if available
                if "pocket_size" in summary:
                    pocket_size = summary["pocket_size"]
        except Exception as exc:
            logger.warning("[%s] Could not load pocket center: %s", opt_id, exc)

        # Resolve receptor path from pipeline summary (the real file is
        # named e.g. P00533_receptor.pdbqt, stored at prepare step).
        target_pdbqt = None
        try:
            from database import get_job as _gj
            parent = _gj(job_id)
            if parent and parent.pipeline_summary_json:
                ps = json.loads(parent.pipeline_summary_json)
                candidate = ps.get("receptor_path", "")
                if candidate and Path(candidate).exists():
                    target_pdbqt = candidate
                    logger.info("[opt:%s] Receptor from pipeline_summary: %s", opt_id[:8], target_pdbqt)
        except Exception as rp_exc:
            logger.warning("[opt:%s] Could not read receptor_path from summary: %s", opt_id[:8], rp_exc)

        # Fallback: glob for *receptor*.pdbqt in work_dir
        if not target_pdbqt:
            receptor_candidates = sorted(work_dir.glob("*receptor*.pdbqt"))
            if receptor_candidates:
                target_pdbqt = str(receptor_candidates[0])
                logger.info("[opt:%s] Receptor from glob: %s", opt_id[:8], target_pdbqt)
            else:
                target_pdbqt = str(work_dir / "receptor.pdbqt")
                logger.warning("[opt:%s] No receptor found, using default: %s", opt_id[:8], target_pdbqt)

        n_iterations = params.get("n_iterations", 10)

        def progress_cb(iteration: int, score: float, message: str) -> None:
            logger.info("[opt:%s] %s", opt_id[:8], message)

        result = run_optimization(
            starting_smiles=params["smiles"],
            starting_name=params.get("molecule_name", "molecule"),
            target_pdbqt=target_pdbqt,
            pocket_center=pocket_center,
            weights=params.get("weights"),
            n_iterations=n_iterations,
            variants_per_iter=params.get("variants_per_iter", 50),
            work_dir=opt_work_dir,
            progress_callback=progress_cb,
            structural_rules=params.get("structural_rules"),
            docking_engine=params.get("docking_engine", "auto"),
            pocket_size=pocket_size,
            exhaustiveness=params.get("exhaustiveness", 8),
            dock_top_n=params.get("dock_top_n", 20),
        )

        # Store result as JSON file for retrieval
        result_path = opt_work_dir / "result.json"
        with open(result_path, "w") as f:
            json.dump(result, f, indent=2)

        logger.info(
            "[opt:%s] Optimization completed: score %.3f -> %.3f",
            opt_id[:8],
            result["starting_molecule"]["score"],
            result["final_lead"]["score"],
        )

        # Auto-create a new job in the project with optimization results
        new_job_id = _create_job_from_optimization(opt_id, params, result)

        return {
            "status": "completed",
            "optimization_id": opt_id,
            "result": result,
            "created_job_id": new_job_id,
        }

    except Exception as exc:
        error_msg = f"{type(exc).__name__}: {exc}"
        logger.error("[opt:%s] Optimization failed: %s\n%s", opt_id[:8], error_msg, traceback.format_exc())
        return {
            "status": "failed",
            "optimization_id": opt_id,
            "error": error_msg,
        }


def _create_job_from_optimization(opt_id: str, params: dict, result: dict) -> Optional[str]:
    """Create a completed job in the project from optimization results.

    Called after optimization finishes. Builds molecule dicts from the
    optimization output, then runs the full enrichment pipeline (ADMET,
    scoring, clustering, Pareto, off-target, hERG, confidence) so that
    optimization results have the same quality of data as HTS runs.

    Parameters
    ----------
    opt_id : str
        Optimization UUID.
    params : dict
        Original optimization parameters (contains ``job_id``).
    result : dict
        Optimization result dict with ``final_lead``, ``best_molecule``, etc.

    Returns
    -------
    str or None
        The new job_id, or None on failure.
    """
    from uuid import uuid4 as _uuid4

    try:
        parent_job_id = params["job_id"]
        from database import get_job as _get_job, create_job as _create_job, update_job as _update_job
        parent_job = _get_job(parent_job_id)
        if parent_job is None:
            logger.warning("Cannot create opt job: parent job %s not found", parent_job_id)
            return None

        new_job_id = str(_uuid4())
        project_id = getattr(parent_job, "project_id", None)
        uniprot_id = getattr(parent_job, "uniprot_id", "") or ""
        user_id = getattr(parent_job, "user_id", None)

        from pipeline.scoring import generate_2d_svg, compute_properties

        def _opt_mol_entry(smiles, name, affinity, score, admet_data, source="optimization", dock_meta=None):
            """Build a result entry for an optimization molecule with RDKit properties."""
            props = compute_properties(smiles) or {}
            dm = dock_meta or {}
            return {
                "name": name,
                "smiles": smiles,
                "affinity": affinity,
                "composite_score": score,
                "source": source,
                "vina_score": dm.get("vina_score", affinity),
                "cnn_score": dm.get("cnn_score", 0.0),
                "cnn_affinity": dm.get("cnn_affinity", 0.0),
                "docking_engine": dm.get("docking_engine", "optimization"),
                "pose_pdbqt_path": dm.get("pose_pdbqt_path", ""),
                "mw": props.get("MW") or admet_data.get("mw", 0),
                "logp": props.get("logP") or admet_data.get("logp", 0),
                "qed": props.get("qed") or admet_data.get("qed", 0),
                "tpsa": props.get("tpsa") or admet_data.get("tpsa", 0),
                "hbd": props.get("hbd"),
                "hba": props.get("hba"),
                "rotatable_bonds": props.get("rotatable_bonds"),
                "toxicity_score": admet_data.get("toxicity_score", 0),
                "bioavailability_score": admet_data.get("bioavailability_score", 0),
                "synthesis_score": admet_data.get("synthesis_score", 0),
                "svg": generate_2d_svg(smiles),
            }

        # Collect all molecules from optimization
        opt_results = []

        # 1. Add best/final lead
        best_mol = result.get("best_molecule") or result.get("final_lead") or {}
        starting_mol = result.get("starting_molecule") or {}
        mol_name = best_mol.get("name", params.get("molecule_name", "optimized"))
        admet = best_mol.get("admet") or {}

        best_smiles = best_mol.get("smiles", params.get("smiles", ""))

        # Extract docking metadata from best molecule
        best_dock_meta = {k: best_mol.get(k) for k in
            ("vina_score", "cnn_score", "cnn_affinity", "docking_engine", "pose_pdbqt_path")
            if best_mol.get(k) is not None}

        entry = _opt_mol_entry(best_smiles, mol_name, best_mol.get("affinity", -8.0), best_mol.get("score", 0.5), admet,
                               dock_meta=best_dock_meta)
        entry["optimization_origin"] = {
            "optimization_id": opt_id,
            "parent_job_id": parent_job_id,
            "starting_molecule": starting_mol.get("name", ""),
            "starting_score": starting_mol.get("score", 0),
            "final_score": best_mol.get("score", 0),
            "iterations": len(result.get("iterations", [])),
        }
        opt_results.append(entry)

        # 2. Add top molecules from iterations (deduplicated)
        seen = {best_smiles}
        for tm in result.get("top_molecules", []):
            smi = tm.get("smiles", "")
            if smi and smi not in seen:
                seen.add(smi)
                tm_admet = tm.get("admet") or {}
                tm_dock_meta = {k: tm.get(k) for k in
                    ("vina_score", "cnn_score", "cnn_affinity", "docking_engine", "pose_pdbqt_path")
                    if tm.get(k) is not None}
                tm_entry = _opt_mol_entry(smi, tm.get("name", f"opt_variant_{len(opt_results)}"),
                                          tm.get("affinity", -7.0), tm.get("score", 0.4), tm_admet,
                                          dock_meta=tm_dock_meta)
                tm_entry["optimization_origin"] = {
                    "optimization_id": opt_id,
                    "parent_job_id": parent_job_id,
                    "iteration": tm.get("iteration", 0),
                }
                opt_results.append(tm_entry)

        # 3. Add starting molecule for comparison
        start_smiles = starting_mol.get("smiles", params.get("smiles", ""))
        start_entry = _opt_mol_entry(start_smiles, f"{starting_mol.get('name', 'starting')} (original)",
                                     starting_mol.get("affinity", -5.0), starting_mol.get("score", 0.3), {},
                                     source="optimization_reference")
        start_entry["docking_engine"] = "reference"
        opt_results.append(start_entry)

        # -----------------------------------------------------------
        # Run full enrichment pipeline (ADMET, scoring, clustering,
        # Pareto, off-target, hERG, confidence)
        # -----------------------------------------------------------
        try:
            from pipeline.enrich import enrich_results, generate_3d_conformer

            enrich_work_dir = DATA_DIR / parent_job_id / "optimization" / opt_id / "enrich"
            enrich_work_dir.mkdir(parents=True, exist_ok=True)

            structure_source = getattr(parent_job, "structure_source", "cached") or "cached"

            # V12: forward individual analysis flags from optimization params
            opt_enable_admet = params.get("enable_admet", True)
            opt_enable_synthesis = params.get("enable_synthesis", True)
            opt_enable_selectivity = params.get("enable_selectivity", True)
            opt_enable_herg = params.get("enable_herg", True)
            opt_enable_safety = params.get("enable_safety", True)

            passed, eliminated, enrich_summary = enrich_results(
                molecules=opt_results,
                work_dir=enrich_work_dir,
                structure_source=structure_source,
                enable_retrosynthesis=opt_enable_synthesis,
                enable_admet=opt_enable_admet,
                enable_selectivity=opt_enable_selectivity,
                enable_herg=opt_enable_herg,
                enable_safety=opt_enable_safety,
            )

            # Use enriched results (passed + eliminated for completeness)
            opt_results = passed + eliminated
            logger.info(
                "[opt:%s] Enrichment: %d passed, %d eliminated, steps=%s",
                opt_id[:8], len(passed), len(eliminated),
                enrich_summary.get("steps_completed", []),
            )

            # Generate 3D conformers for pose-less molecules
            for i, mol in enumerate(opt_results):
                smi = mol.get("smiles", "")
                if smi and not mol.get("pose_sdf_path"):
                    sdf_path = enrich_work_dir / f"conformer_{i}.sdf"
                    result_path = generate_3d_conformer(smi, sdf_path)
                    if result_path:
                        mol["pose_sdf_path"] = str(result_path)

            # Run interaction analysis for molecules with real pose files
            pdb_path = getattr(parent_job, "pdb_path", None)
            if pdb_path:
                mols_with_pose = [m for m in opt_results
                                  if m.get("pose_pdbqt_path") and os.path.exists(str(m["pose_pdbqt_path"]))]
                if mols_with_pose:
                    try:
                        n_analyzed = _safe_interaction_analysis(
                            new_job_id, mols_with_pose, pdb_path, uniprot_id)
                        logger.info("[opt:%s] Interaction analysis: %d/%d molecules",
                                    opt_id[:8], n_analyzed, len(mols_with_pose))
                    except Exception as ia_exc:
                        logger.warning("[opt:%s] Interaction analysis failed: %s", opt_id[:8], ia_exc)
        except Exception as enrich_exc:
            logger.warning("[opt:%s] Enrichment failed (using raw results): %s", opt_id[:8], enrich_exc)
            enrich_summary = {"steps_completed": [], "steps_failed": ["all"]}

        results_json = json.dumps(opt_results)
        pipeline_summary = {
            "mode": "optimization",
            "steps_completed": ["optimization"] + enrich_summary.get("steps_completed", []),
            "optimization_id": opt_id,
            "parent_job_id": parent_job_id,
            "total_results": len(opt_results),
            "n_passed": len([m for m in opt_results if not m.get("eliminated")]),
            "best_affinity": best_mol.get("affinity", -8.0),
            "best_composite": best_mol.get("score", 0.5),
            "comparison": result.get("comparison", {}),
            "enrichment": enrich_summary,
        }

        _create_job(
            job_id=new_job_id,
            uniprot_id=uniprot_id,
            use_chembl=False,
            use_zinc=False,
            max_ligands=1,
            smiles_list=[best_mol.get("smiles", "")],
            mode="optimization",
            project_id=project_id,
            user_id=user_id,
        )

        _update_job(
            new_job_id,
            status="completed",
            progress=100,
            current_step="Optimization complete",
            results_json=results_json,
            pipeline_summary_json=json.dumps(pipeline_summary),
            protein_name=getattr(parent_job, "protein_name", uniprot_id) or uniprot_id,
            pdb_path=getattr(parent_job, "pdb_path", None),
            structure_source=getattr(parent_job, "structure_source", "cached"),
            completed_at=datetime.now(timezone.utc),
        )

        # Write meta.json for disk recovery
        meta_dir = DATA_DIR / parent_job_id / "optimization" / opt_id
        meta_dir.mkdir(parents=True, exist_ok=True)
        try:
            (meta_dir / "meta.json").write_text(
                json.dumps({"created_job_id": new_job_id})
            )
        except Exception:
            pass

        logger.info(
            "[opt:%s] Created optimization job %s (project=%s)",
            opt_id[:8], new_job_id, project_id,
        )
        return new_job_id

    except Exception as exc:
        logger.warning("[opt:%s] Failed to create job from optimization: %s", opt_id[:8], exc)
        return None


def _make_serializable(results: list[dict]) -> list[dict]:
    """Strip or convert non-JSON-serializable values from result dicts."""
    clean: list[dict] = []
    for r in results:
        entry: dict = {}
        for k, v in r.items():
            if v is None:
                entry[k] = None
            elif isinstance(v, (str, int, float, bool)):
                entry[k] = v
            elif isinstance(v, list):
                entry[k] = v
            elif isinstance(v, dict):
                entry[k] = v
            elif isinstance(v, Path):
                entry[k] = str(v)
            else:
                entry[k] = str(v)
        clean.append(entry)
    return clean
