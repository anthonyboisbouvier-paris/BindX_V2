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
    "pharmacophore", "activity_cliffs",
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
        run_all = (body.config or {}).get("run_all_molecules", False)
        if not body.input_molecule_ids and not run_all:
            raise HTTPException(
                status_code=422, detail="input_molecule_ids required for calculation runs"
            )

    # Validate import has data
    if body.type == "import":
        config = body.config or {}
        has_smiles = bool(config.get("smiles_list"))
        has_phase_source = config.get("source") == "phase_selection" and config.get("source_phase_id")
        has_db_source = config.get("source") == "database" and config.get("databases")
        has_ligand_list = bool(config.get("ligand_list"))
        if not has_smiles and not has_phase_source and not has_db_source and not has_ligand_list:
            raise HTTPException(
                status_code=422,
                detail="config.smiles_list, config.source='phase_selection', or config.source='database' required for import runs",
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
    await db.commit()

    # Dispatch Celery task (after commit so the worker can find the run)
    try:
        from celery_app import celery_app
        celery_app.send_task("tasks_v9.execute_run", args=[str(run.id)])
        logger.info("Dispatched Celery task for run %s", run.id)
    except Exception as e:
        logger.warning("Failed to dispatch Celery task for run %s: %s", run.id, e)
        run.status = "failed"
        run.error_message = f"Failed to dispatch task: {e}"
        await db.commit()

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


@router.get("/runs/{run_id}/progress")
async def get_run_progress(
    run_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Lightweight progress endpoint for polling during run execution."""
    run = await get_run_owned(run_id, user_id, db)
    return {
        "id": run.id,
        "status": run.status,
        "progress": run.progress,
        "current_step": run.current_step,
        "estimated_duration_seconds": run.estimated_duration_seconds,
        "started_at": run.started_at,
        "completed_at": run.completed_at,
        "error_message": run.error_message,
    }


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
    await db.commit()

    # Dispatch Celery task (after commit so the worker can find the run)
    try:
        from celery_app import celery_app
        celery_app.send_task("tasks_v9.execute_run", args=[str(run.id)])
        logger.info("Dispatched Celery task for file import run %s", run.id)
    except Exception as e:
        logger.warning("Failed to dispatch Celery task for file run %s: %s", run.id, e)
        run.status = "failed"
        run.error_message = f"Failed to dispatch task: {e}"
        await db.commit()

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
    """Create a database import run. Fetch happens async in Celery worker.

    Body: {databases: ["chembl","pubchem",...], uniprot_id, max_per_source}
    Response is instant — no ligand_list stored in config.
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

    # Store only parameters — fetch happens in Celery worker
    run = RunORM_V9(
        phase_id=phase_id,
        type="import",
        config={
            "source": "database",
            "databases": databases,
            "uniprot_id": body.get("uniprot_id", ""),
            "max_per_source": body.get("max_per_source", 50),
        },
        input_source="connected_library",
        status="created",
    )
    db.add(run)
    await db.commit()

    # Dispatch Celery task (after commit so the worker can find the run)
    try:
        from celery_app import celery_app
        celery_app.send_task("tasks_v9.execute_run", args=[str(run.id)])
        logger.info("Dispatched database import run %s: %s", run.id, databases)
    except Exception as e:
        logger.warning("Failed to dispatch Celery task for db run %s: %s", run.id, e)
        run.status = "failed"
        run.error_message = f"Failed to dispatch task: {e}"
        await db.commit()

    # Reload with logs eagerly loaded
    stmt = select(RunORM_V9).where(RunORM_V9.id == run.id).options(selectinload(RunORM_V9.logs))
    result = await db.execute(stmt)
    return result.scalar_one()
