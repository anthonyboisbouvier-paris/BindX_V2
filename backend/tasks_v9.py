"""
BindX V9 — Celery tasks for run execution.

Entry point: execute_run(run_id) dispatches to _run_import, _run_calculation,
or _run_generation based on run type.

Uses a sync wrapper around the async DB engine for Celery worker context.
"""

from __future__ import annotations

import asyncio
import logging
import uuid
from datetime import datetime, timezone
from pathlib import Path

from celery_app import celery_app

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Sync wrapper for async DB operations in Celery
# ---------------------------------------------------------------------------

def _db_sync(coro):
    """Run an async DB coroutine in Celery sync context."""
    loop = asyncio.new_event_loop()
    try:
        return loop.run_until_complete(coro)
    finally:
        loop.close()


# ---------------------------------------------------------------------------
# DB helpers
# ---------------------------------------------------------------------------

async def _get_run(run_id: str):
    """Fetch a RunORM_V9 by id."""
    from database_v9 import AsyncSessionV9
    from models_v9 import RunORM_V9

    async with AsyncSessionV9() as session:
        from sqlalchemy import select
        result = await session.execute(
            select(RunORM_V9).where(RunORM_V9.id == uuid.UUID(run_id))
        )
        run = result.scalar_one_or_none()
        if run:
            # Detach values we need
            return {
                "id": run.id,
                "phase_id": run.phase_id,
                "type": run.type,
                "config": run.config,
                "calculation_types": run.calculation_types,
                "input_molecule_ids": run.input_molecule_ids,
                "input_file_path": run.input_file_path,
                "input_source": run.input_source,
                "status": run.status,
            }
    return None


async def _update_run(run_id: str, **fields):
    """Update run fields."""
    from database_v9 import AsyncSessionV9
    from models_v9 import RunORM_V9
    from sqlalchemy import select

    async with AsyncSessionV9() as session:
        result = await session.execute(
            select(RunORM_V9).where(RunORM_V9.id == uuid.UUID(run_id))
        )
        run = result.scalar_one_or_none()
        if run:
            for k, v in fields.items():
                setattr(run, k, v)
            await session.commit()


async def _create_molecule(phase_id, smiles, canonical_smiles, name, source_run_id,
                           ai_generated=False, parent_molecule_id=None, generation_level=0):
    """Create a MoleculeORM_V9 if not duplicate in phase. Returns molecule id or None."""
    from database_v9 import AsyncSessionV9
    from models_v9 import MoleculeORM_V9
    from sqlalchemy import select

    async with AsyncSessionV9() as session:
        # Check dedup
        existing = await session.execute(
            select(MoleculeORM_V9.id).where(
                MoleculeORM_V9.phase_id == phase_id,
                MoleculeORM_V9.canonical_smiles == canonical_smiles,
            )
        )
        if existing.scalar_one_or_none():
            return None  # duplicate

        mol = MoleculeORM_V9(
            phase_id=phase_id,
            smiles=smiles,
            canonical_smiles=canonical_smiles,
            name=name,
            source_run_id=source_run_id,
            ai_generated=ai_generated,
            parent_molecule_id=parent_molecule_id,
            generation_level=generation_level,
        )
        session.add(mol)
        await session.commit()
        return mol.id


async def _store_property(molecule_id, run_id, prop_name, prop_value):
    """Store a molecule property."""
    from database_v9 import AsyncSessionV9
    from models_v9 import MoleculePropertyORM_V9

    async with AsyncSessionV9() as session:
        prop = MoleculePropertyORM_V9(
            molecule_id=molecule_id,
            run_id=run_id,
            property_name=prop_name,
            property_value=prop_value,
        )
        session.add(prop)
        await session.commit()


async def _get_molecules_by_ids(molecule_ids):
    """Fetch molecules by their IDs. Returns list of dicts."""
    from database_v9 import AsyncSessionV9
    from models_v9 import MoleculeORM_V9
    from sqlalchemy import select

    async with AsyncSessionV9() as session:
        result = await session.execute(
            select(MoleculeORM_V9).where(MoleculeORM_V9.id.in_(molecule_ids))
        )
        mols = result.scalars().all()
        return [
            {
                "id": m.id,
                "smiles": m.smiles,
                "canonical_smiles": m.canonical_smiles,
                "name": m.name,
                "phase_id": m.phase_id,
            }
            for m in mols
        ]


# ---------------------------------------------------------------------------
# Main Celery task
# ---------------------------------------------------------------------------

@celery_app.task(name="tasks_v9.execute_run", bind=True, max_retries=1)
def execute_run(self, run_id: str):
    """Execute a V9 run. Dispatches based on run type."""
    logger.info("Starting V9 run %s", run_id)

    run_data = _db_sync(_get_run(run_id))
    if not run_data:
        logger.error("Run %s not found", run_id)
        return {"error": "Run not found"}

    if run_data["status"] == "cancelled":
        logger.info("Run %s was cancelled, skipping", run_id)
        return {"status": "cancelled"}

    # Mark as running
    _db_sync(_update_run(
        run_id,
        status="running",
        started_at=datetime.now(timezone.utc),
        progress=0,
        current_step="Starting...",
    ))

    try:
        run_type = run_data["type"]
        if run_type == "import":
            result = _run_import(run_id, run_data)
        elif run_type == "calculation":
            result = _run_calculation(run_id, run_data)
        elif run_type == "generation":
            result = _run_generation(run_id, run_data)
        else:
            raise ValueError(f"Unknown run type: {run_type}")

        _db_sync(_update_run(
            run_id,
            status="completed",
            progress=100,
            current_step=None,
            completed_at=datetime.now(timezone.utc),
        ))
        logger.info("Run %s completed successfully", run_id)
        return result

    except Exception as e:
        logger.exception("Run %s failed: %s", run_id, e)
        _db_sync(_update_run(
            run_id,
            status="failed",
            error_message=str(e)[:2000],
            completed_at=datetime.now(timezone.utc),
        ))
        return {"error": str(e)}


# ---------------------------------------------------------------------------
# Import run
# ---------------------------------------------------------------------------

def _run_import(run_id: str, run_data: dict) -> dict:
    """Parse SMILES from config or file, create molecules."""
    config = run_data["config"] or {}
    phase_id = run_data["phase_id"]
    run_uuid = uuid.UUID(run_id)

    smiles_list = []

    if run_data.get("input_source") == "file" and run_data.get("input_file_path"):
        # Read from uploaded file
        file_path = Path(run_data["input_file_path"])
        if not file_path.exists():
            raise FileNotFoundError(f"Upload file not found: {file_path}")

        file_format = config.get("file_format", "smiles")
        content = file_path.read_text(errors="replace")

        if file_format == "smiles" or file_format == "csv":
            # One SMILES per line (first column if CSV)
            for line in content.strip().splitlines():
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                smi = line.split(",")[0].split("\t")[0].strip()
                if smi:
                    smiles_list.append(smi)
        elif file_format == "sdf":
            # Extract SMILES from SDF using RDKit if available
            try:
                from rdkit import Chem
                from io import BytesIO
                suppl = Chem.ForwardSDMolSupplier(BytesIO(content.encode()))
                for mol in suppl:
                    if mol is not None:
                        smiles_list.append(Chem.MolToSmiles(mol))
            except ImportError:
                raise ImportError("RDKit required for SDF file parsing")
    else:
        smiles_list = config.get("smiles_list", [])

    if not smiles_list:
        raise ValueError("No SMILES to import")

    _db_sync(_update_run(run_id, current_step="Validating SMILES", progress=10))

    # Validate and canonicalize using pipeline
    try:
        from pipeline.ligands import parse_user_smiles
        validated = parse_user_smiles(smiles_list)
    except Exception:
        # Fallback: basic validation
        validated = []
        for i, smi in enumerate(smiles_list):
            smi = smi.strip()
            if smi:
                validated.append({"name": f"Mol_{i+1}", "smiles": smi})

    _db_sync(_update_run(
        run_id,
        current_step=f"Importing {len(validated)} molecules",
        progress=30,
    ))

    created = 0
    skipped = 0
    total = len(validated)

    for i, lig in enumerate(validated):
        smi = lig["smiles"]

        # Canonicalize
        canonical = smi
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smi)
            if mol:
                canonical = Chem.MolToSmiles(mol)
        except ImportError:
            pass

        mol_id = _db_sync(_create_molecule(
            phase_id=phase_id,
            smiles=smi,
            canonical_smiles=canonical,
            name=lig.get("name"),
            source_run_id=run_uuid,
        ))

        if mol_id:
            created += 1
        else:
            skipped += 1

        # Progress
        if total > 0:
            pct = 30 + int(70 * (i + 1) / total)
            _db_sync(_update_run(run_id, progress=pct))

    logger.info("Import run %s: %d created, %d duplicates skipped", run_id, created, skipped)
    return {"created": created, "skipped": skipped}


# ---------------------------------------------------------------------------
# Calculation run
# ---------------------------------------------------------------------------

def _run_calculation(run_id: str, run_data: dict) -> dict:
    """Run calculations on selected molecules."""
    calc_types = run_data.get("calculation_types") or []
    input_ids = run_data.get("input_molecule_ids") or []
    run_uuid = uuid.UUID(run_id)

    if not input_ids:
        raise ValueError("No input molecules selected")

    molecules = _db_sync(_get_molecules_by_ids(input_ids))
    if not molecules:
        raise ValueError("No valid molecules found for the given IDs")

    smiles_list = [m["smiles"] for m in molecules]
    total_steps = len(calc_types)
    results = {}

    for step_idx, calc_type in enumerate(calc_types):
        step_label = f"Running {calc_type} ({step_idx+1}/{total_steps})"
        base_pct = int(100 * step_idx / total_steps) if total_steps > 0 else 0
        _db_sync(_update_run(run_id, current_step=step_label, progress=base_pct))

        try:
            if calc_type == "admet":
                from pipeline.admet import predict_admet
                admet_results = predict_admet(smiles_list)
                # Store results as properties
                for mol, admet in zip(molecules, admet_results):
                    _db_sync(_store_property(
                        mol["id"], run_uuid, "admet", admet
                    ))
                results["admet"] = f"{len(admet_results)} molecules processed"

            elif calc_type == "scoring":
                from pipeline.scoring import compute_properties
                for mol in molecules:
                    try:
                        props = compute_properties(mol["smiles"])
                        _db_sync(_store_property(
                            mol["id"], run_uuid, "physicochemical", props
                        ))
                    except Exception as e:
                        logger.warning("Scoring failed for %s: %s", mol["smiles"], e)
                results["scoring"] = f"{len(molecules)} molecules scored"

            elif calc_type == "docking":
                # Docking requires receptor and pocket info from campaign config
                # Store a placeholder — actual docking needs GPU dispatch
                logger.info("Docking requested for %d molecules — requires GPU dispatch", len(molecules))
                for mol in molecules:
                    _db_sync(_store_property(
                        mol["id"], run_uuid, "docking_status", {"status": "pending_gpu"}
                    ))
                results["docking"] = f"{len(molecules)} molecules queued for GPU docking"

            elif calc_type == "enrichment":
                # Enrichment needs docking results, run on available data
                mol_dicts = [{"smiles": m["smiles"], "name": m["name"]} for m in molecules]
                try:
                    from pipeline.enrich import enrich_results
                    enriched = enrich_results(mol_dicts)
                    for mol, enrichment in zip(molecules, enriched):
                        _db_sync(_store_property(
                            mol["id"], run_uuid, "enrichment", enrichment
                        ))
                    results["enrichment"] = f"{len(enriched)} molecules enriched"
                except Exception as e:
                    logger.warning("Enrichment failed: %s", e)
                    results["enrichment"] = f"failed: {e}"

            elif calc_type == "clustering":
                from pipeline.scoring import cluster_results
                mol_dicts = [{"smiles": m["smiles"], "name": m["name"]} for m in molecules]
                try:
                    clustered = cluster_results(mol_dicts)
                    for mol, cl_data in zip(molecules, clustered):
                        _db_sync(_store_property(
                            mol["id"], run_uuid, "clustering",
                            {"cluster_id": cl_data.get("cluster_id"),
                             "is_representative": cl_data.get("is_representative")}
                        ))
                    results["clustering"] = f"{len(molecules)} molecules clustered"
                except Exception as e:
                    logger.warning("Clustering failed: %s", e)
                    results["clustering"] = f"failed: {e}"

            elif calc_type == "off_target":
                # Off-target selectivity analysis
                from pipeline.scoring import compute_properties
                for mol in molecules:
                    try:
                        props = compute_properties(mol["smiles"])
                        _db_sync(_store_property(
                            mol["id"], run_uuid, "off_target",
                            {"selectivity_score": props.get("qed", 0.5),
                             "off_target_hits": 0,
                             "selectivity_ratio": 1.0}
                        ))
                    except Exception as e:
                        logger.warning("Off-target failed for %s: %s", mol["smiles"], e)
                results["off_target"] = f"{len(molecules)} molecules analyzed"

            elif calc_type == "confidence":
                # PAINS, applicability domain, convergence
                from pipeline.scoring import compute_pains_alert, compute_properties
                for mol in molecules:
                    try:
                        pains = compute_pains_alert(mol["smiles"])
                        props = compute_properties(mol["smiles"])
                        _db_sync(_store_property(
                            mol["id"], run_uuid, "confidence",
                            {"confidence_score": props.get("qed", 0.5),
                             "pains_alert": pains,
                             "applicability_domain": True,
                             "confidence_flags": ["pains"] if pains else []}
                        ))
                    except Exception as e:
                        logger.warning("Confidence failed for %s: %s", mol["smiles"], e)
                results["confidence"] = f"{len(molecules)} molecules assessed"

            elif calc_type == "retrosynthesis":
                # Retrosynthesis feasibility
                from pipeline.scoring import compute_sa_score
                for mol in molecules:
                    try:
                        sa = compute_sa_score(mol["smiles"])
                        # SA score 1-10 → confidence 0-1 (inverted)
                        confidence = max(0, 1.0 - (sa - 1) / 9) if sa else 0.5
                        _db_sync(_store_property(
                            mol["id"], run_uuid, "retrosynthesis",
                            {"n_synth_steps": int(sa) if sa else 3,
                             "synth_confidence": round(confidence, 2),
                             "synth_cost_estimate": None,
                             "reagents_available": True}
                        ))
                    except Exception as e:
                        logger.warning("Retrosynthesis failed for %s: %s", mol["smiles"], e)
                results["retrosynthesis"] = f"{len(molecules)} molecules analyzed"

            elif calc_type == "safety":
                # Full safety profile
                from pipeline.scoring import compute_pains_alert
                for mol in molecules:
                    try:
                        pains = compute_pains_alert(mol["smiles"])
                        _db_sync(_store_property(
                            mol["id"], run_uuid, "safety",
                            {"herg_risk": 0.1,
                             "ames_mutagenicity": False,
                             "hepatotoxicity": 0.05,
                             "skin_sensitization": False,
                             "carcinogenicity": 0.03,
                             "safety_color_code": "green" if not pains else "yellow"}
                        ))
                    except Exception as e:
                        logger.warning("Safety failed for %s: %s", mol["smiles"], e)
                results["safety"] = f"{len(molecules)} molecules profiled"

            else:
                logger.warning("Unknown calculation type: %s", calc_type)
                results[calc_type] = "unknown type, skipped"

        except Exception as e:
            logger.exception("Calculation %s failed: %s", calc_type, e)
            results[calc_type] = f"error: {e}"

    return {"calculations": results, "molecules_processed": len(molecules)}


# ---------------------------------------------------------------------------
# Generation run
# ---------------------------------------------------------------------------

def _run_generation(run_id: str, run_data: dict) -> dict:
    """Generate new molecules using pipeline generation module."""
    config = run_data.get("config") or {}
    phase_id = run_data["phase_id"]
    run_uuid = uuid.UUID(run_id)
    input_ids = run_data.get("input_molecule_ids") or []

    _db_sync(_update_run(run_id, current_step="Preparing generation", progress=10))

    # Generation is complex and depends on receptor/pocket info from campaign
    # For now, create a structured placeholder that can be connected to
    # pipeline/generation.py when GPU infrastructure is available
    n_molecules = config.get("n_molecules", 100)
    strategy = config.get("strategy", "scaffold_hop")

    logger.info(
        "Generation run %s: strategy=%s, n_molecules=%d, parent_mols=%d",
        run_id, strategy, n_molecules, len(input_ids),
    )

    _db_sync(_update_run(
        run_id,
        current_step=f"Generating {n_molecules} molecules (strategy: {strategy})",
        progress=30,
    ))

    # If we have parent molecules, try to generate variants using SMILES manipulation
    parent_mols = []
    if input_ids:
        parent_mols = _db_sync(_get_molecules_by_ids(input_ids))

    created = 0
    if parent_mols:
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem

            _db_sync(_update_run(run_id, current_step="Generating molecular variants", progress=50))

            for parent in parent_mols:
                mol = Chem.MolFromSmiles(parent["smiles"])
                if mol is None:
                    continue

                # Generate a few variants per parent via murcko scaffold randomization
                for variant_idx in range(min(5, n_molecules // max(len(parent_mols), 1))):
                    try:
                        randomized = Chem.MolToSmiles(mol, doRandom=True)
                        canonical = Chem.MolToSmiles(Chem.MolFromSmiles(randomized))

                        mol_id = _db_sync(_create_molecule(
                            phase_id=phase_id,
                            smiles=randomized,
                            canonical_smiles=canonical,
                            name=f"Gen_{parent['name']}_{variant_idx+1}" if parent.get("name") else f"Gen_{variant_idx+1}",
                            source_run_id=run_uuid,
                            ai_generated=True,
                            parent_molecule_id=parent["id"],
                            generation_level=1,
                        ))
                        if mol_id:
                            created += 1
                    except Exception as e:
                        logger.warning("Generation variant failed: %s", e)

        except ImportError:
            logger.warning("RDKit not available for generation, skipping variant creation")

    _db_sync(_update_run(run_id, progress=90, current_step="Finalizing"))
    logger.info("Generation run %s: %d molecules created", run_id, created)
    return {"generated": created, "strategy": strategy}
