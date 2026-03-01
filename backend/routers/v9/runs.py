"""
BindX V9 — Run CRUD router.

6 endpoints: create run, list runs, get run, cancel run, archive run, import file.
Dispatches Celery tasks for async execution.
"""

from __future__ import annotations

import logging
import os
import tempfile
from pathlib import Path
from uuid import UUID

from fastapi import APIRouter, Body, Depends, File, HTTPException, UploadFile
from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from auth_v9 import require_v9_user
from database_v9 import get_v9_db
from models_v9 import RunORM_V9, RunLogORM_V9
from sqlalchemy.orm import selectinload
from routers.v9.deps import get_phase_owned, get_run_owned
from schemas_v9 import RunCreate, RunListItem, RunResponse

logger = logging.getLogger(__name__)

router = APIRouter()

VALID_RUN_TYPES = {"import", "calculation", "generation"}
VALID_CALCULATION_TYPES = {
    "docking", "admet", "scoring", "enrichment", "clustering",
    "off_target", "confidence", "retrosynthesis", "safety",
}


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------

@router.post(
    "/phases/{phase_id}/runs",
    response_model=RunResponse,
    status_code=201,
)
async def create_run(
    phase_id: UUID,
    body: RunCreate,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Create a new run and dispatch Celery task."""
    if body.type not in VALID_RUN_TYPES:
        raise HTTPException(
            status_code=422,
            detail=f"Invalid run type. Must be one of: {', '.join(VALID_RUN_TYPES)}",
        )

    phase = await get_phase_owned(phase_id, user_id, db)

    if phase.status == "frozen":
        raise HTTPException(status_code=400, detail="Cannot create runs on a frozen phase")

    # Validate calculation types
    if body.type == "calculation":
        if not body.calculation_types:
            raise HTTPException(
                status_code=422, detail="calculation_types required for calculation runs"
            )
        invalid = set(body.calculation_types) - VALID_CALCULATION_TYPES
        if invalid:
            raise HTTPException(
                status_code=422,
                detail=f"Invalid calculation types: {', '.join(invalid)}",
            )
        if not body.input_molecule_ids:
            raise HTTPException(
                status_code=422, detail="input_molecule_ids required for calculation runs"
            )

    # Validate import has data
    if body.type == "import":
        config = body.config or {}
        has_smiles = bool(config.get("smiles_list"))
        has_phase_source = config.get("source") == "phase_selection" and config.get("source_phase_id")
        if not has_smiles and not has_phase_source:
            raise HTTPException(
                status_code=422,
                detail="config.smiles_list or config.source='phase_selection' required for import runs",
            )

    run = RunORM_V9(
        phase_id=phase_id,
        type=body.type,
        config=body.config or {},
        calculation_types=body.calculation_types,
        input_molecule_ids=body.input_molecule_ids,
        input_source="api",
        status="created",
    )
    db.add(run)
    await db.flush()

    # Dispatch Celery task
    try:
        from celery_app import celery_app
        celery_app.send_task("tasks_v9.execute_run", args=[str(run.id)])
        logger.info("Dispatched Celery task for run %s", run.id)
    except Exception as e:
        logger.warning("Failed to dispatch Celery task for run %s: %s", run.id, e)
        run.status = "failed"
        run.error_message = f"Failed to dispatch task: {e}"
        await db.flush()

    # Reload with logs eagerly loaded (avoid lazy-load crash in async)
    stmt = select(RunORM_V9).where(RunORM_V9.id == run.id).options(selectinload(RunORM_V9.logs))
    result = await db.execute(stmt)
    return result.scalar_one()


@router.get(
    "/phases/{phase_id}/runs",
    response_model=list[RunListItem],
)
async def list_runs(
    phase_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """List all runs for a phase."""
    await get_phase_owned(phase_id, user_id, db)

    stmt = (
        select(RunORM_V9)
        .where(RunORM_V9.phase_id == phase_id)
        .options(selectinload(RunORM_V9.logs))
        .order_by(RunORM_V9.created_at.desc())
    )
    result = await db.execute(stmt)
    return result.scalars().all()


@router.get("/runs/{run_id}", response_model=RunResponse)
async def get_run(
    run_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Get run detail with progress."""
    return await get_run_owned(run_id, user_id, db)


@router.post("/runs/{run_id}/cancel", response_model=RunResponse)
async def cancel_run(
    run_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Cancel a running run."""
    run = await get_run_owned(run_id, user_id, db)
    if run.status not in ("created", "running"):
        raise HTTPException(
            status_code=400, detail=f"Cannot cancel run with status '{run.status}'"
        )
    run.status = "cancelled"
    await db.flush()

    # Try to revoke Celery task
    try:
        from celery_app import celery_app
        celery_app.control.revoke(str(run.id), terminate=True)
    except Exception as e:
        logger.warning("Failed to revoke Celery task for run %s: %s", run.id, e)

    return run  # logs already eager-loaded by get_run_owned


@router.post("/runs/{run_id}/archive", response_model=RunResponse)
async def archive_run(
    run_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Archive a run (soft delete)."""
    run = await get_run_owned(run_id, user_id, db)
    if run.status in ("created", "running"):
        raise HTTPException(
            status_code=400, detail="Cannot archive a run that is still running"
        )
    run.archived = True
    await db.flush()
    return run


@router.post(
    "/phases/{phase_id}/runs/import-file",
    response_model=RunResponse,
    status_code=201,
)
async def import_file(
    phase_id: UUID,
    file: UploadFile = File(...),
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Upload an SDF/SMILES/CSV file and create an import run."""
    phase = await get_phase_owned(phase_id, user_id, db)

    if phase.status == "frozen":
        raise HTTPException(status_code=400, detail="Cannot import to a frozen phase")

    # Determine format from extension
    filename = file.filename or "unknown"
    suffix = Path(filename).suffix.lower()
    format_map = {".sdf": "sdf", ".smi": "smiles", ".smiles": "smiles", ".csv": "csv"}
    file_format = format_map.get(suffix)
    if not file_format:
        raise HTTPException(
            status_code=422,
            detail=f"Unsupported file format '{suffix}'. Accepted: .sdf, .smi, .smiles, .csv",
        )

    # Save uploaded file to shared volume (accessible by both backend & celery worker)
    content = await file.read()
    upload_dir = Path(os.environ.get("DATA_DIR", "/data")) / "uploads"
    upload_dir.mkdir(parents=True, exist_ok=True)
    file_path = upload_dir / f"{phase_id}_{filename}"
    file_path.write_bytes(content)

    run = RunORM_V9(
        phase_id=phase_id,
        type="import",
        config={"file_format": file_format, "original_filename": filename},
        input_source="file",
        input_file_path=str(file_path),
        status="created",
    )
    db.add(run)
    await db.flush()

    # Dispatch Celery task
    try:
        from celery_app import celery_app
        celery_app.send_task("tasks_v9.execute_run", args=[str(run.id)])
        logger.info("Dispatched Celery task for file import run %s", run.id)
    except Exception as e:
        logger.warning("Failed to dispatch Celery task for file run %s: %s", run.id, e)
        run.status = "failed"
        run.error_message = f"Failed to dispatch task: {e}"
        await db.flush()

    # Reload with logs eagerly loaded
    stmt = select(RunORM_V9).where(RunORM_V9.id == run.id).options(selectinload(RunORM_V9.logs))
    result = await db.execute(stmt)
    return result.scalar_one()


# ---------------------------------------------------------------------------
# Database import
# ---------------------------------------------------------------------------

VALID_DATABASES = {"chembl", "pubchem", "zinc", "enamine", "fragments"}


@router.post(
    "/phases/{phase_id}/runs/import-database",
    response_model=RunResponse,
    status_code=201,
)
async def import_from_database(
    phase_id: UUID,
    body: dict = Body(...),
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Fetch compounds from public databases and create an import run.

    Body: {databases: ["chembl","pubchem",...], uniprot_id, max_per_source}
    """
    phase = await get_phase_owned(phase_id, user_id, db)

    if phase.status == "frozen":
        raise HTTPException(status_code=400, detail="Cannot import to a frozen phase")

    databases = body.get("databases", [])
    if not databases:
        raise HTTPException(status_code=422, detail="No databases selected")

    invalid = set(databases) - VALID_DATABASES
    if invalid:
        raise HTTPException(
            status_code=422,
            detail=f"Invalid databases: {', '.join(invalid)}. "
            f"Valid: {', '.join(sorted(VALID_DATABASES))}",
        )

    uniprot_id = body.get("uniprot_id", "")
    max_per_source = body.get("max_per_source", 50)

    import asyncio as _asyncio

    def _fetch():
        all_smiles = []
        sources_summary = {}

        for db_name in databases:
            try:
                if db_name == "chembl" and uniprot_id:
                    from pipeline.ligands import fetch_chembl_ligands
                    ligands = fetch_chembl_ligands(uniprot_id, max_count=max_per_source)
                    all_smiles.extend(ligands)
                    sources_summary["chembl"] = len(ligands)

                elif db_name == "pubchem" and uniprot_id:
                    from pipeline.ligands import fetch_pubchem_ligands, resolve_gene_name
                    gene = resolve_gene_name(uniprot_id)
                    if gene:
                        ligands = fetch_pubchem_ligands(gene, max_count=max_per_source)
                        all_smiles.extend(ligands)
                        sources_summary["pubchem"] = len(ligands)
                    else:
                        sources_summary["pubchem"] = 0

                elif db_name == "zinc":
                    from pipeline.ligands import fetch_zinc_ligands
                    ligands = fetch_zinc_ligands(max_count=max_per_source)
                    all_smiles.extend(ligands)
                    sources_summary["zinc"] = len(ligands)

                elif db_name == "enamine":
                    from pipeline.ligands import sample_enamine_real
                    ligands_raw = sample_enamine_real(n=max_per_source)
                    for smi in ligands_raw:
                        if isinstance(smi, dict):
                            all_smiles.append(smi)
                        else:
                            all_smiles.append({"smiles": smi, "name": None, "source": "enamine_real"})
                    sources_summary["enamine"] = len(ligands_raw)

                elif db_name == "fragments":
                    from pipeline.ligands import load_fragment_library
                    frags = load_fragment_library()
                    all_smiles.extend(frags)
                    sources_summary["fragments"] = len(frags)

            except Exception as exc:
                logger.warning("Failed to fetch from %s: %s", db_name, exc)
                sources_summary[db_name] = 0

        return all_smiles, sources_summary

    loop = _asyncio.get_event_loop()
    all_smiles, sources_summary = await loop.run_in_executor(None, _fetch)

    if not all_smiles:
        raise HTTPException(
            status_code=422, detail="No compounds found from selected databases",
        )

    smiles_list = [lig["smiles"] for lig in all_smiles if lig.get("smiles")]

    run = RunORM_V9(
        phase_id=phase_id,
        type="import",
        config={
            "smiles_list": smiles_list,
            "databases": databases,
            "sources_summary": sources_summary,
            "uniprot_id": uniprot_id,
            "max_per_source": max_per_source,
        },
        input_source="connected_library",
        status="created",
    )
    db.add(run)
    await db.flush()

    try:
        from celery_app import celery_app
        celery_app.send_task("tasks_v9.execute_run", args=[str(run.id)])
        logger.info(
            "Dispatched database import run %s: %s (%d compounds)",
            run.id, sources_summary, len(smiles_list),
        )
    except Exception as e:
        logger.warning("Failed to dispatch Celery task for db run %s: %s", run.id, e)
        run.status = "failed"
        run.error_message = f"Failed to dispatch task: {e}"
        await db.flush()

    # Reload with logs eagerly loaded
    stmt = select(RunORM_V9).where(RunORM_V9.id == run.id).options(selectinload(RunORM_V9.logs))
    result = await db.execute(stmt)
    return result.scalar_one()
