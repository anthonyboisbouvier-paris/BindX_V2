"""
BindX V9 — Run CRUD router.

6 endpoints: create run, list runs, get run, cancel run, archive run, import file.
Dispatches Celery tasks for async execution.
"""

from __future__ import annotations

import logging
import tempfile
from pathlib import Path
from uuid import UUID

from fastapi import APIRouter, Depends, File, HTTPException, UploadFile
from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.orm import selectinload

from auth_v9 import require_v9_user
from database_v9 import get_v9_db
from models_v9 import (
    CampaignORM_V9,
    PhaseORM_V9,
    ProjectORM_V9,
    RunORM_V9,
)
from schemas_v9 import RunCreate, RunListItem, RunResponse

logger = logging.getLogger(__name__)

router = APIRouter()

VALID_RUN_TYPES = {"import", "calculation", "generation"}
VALID_CALCULATION_TYPES = {
    "docking", "admet", "scoring", "enrichment", "clustering",
    "off_target", "confidence", "retrosynthesis", "safety",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

async def _get_phase_for_run(
    phase_id: UUID,
    user_id: str,
    db: AsyncSession,
) -> PhaseORM_V9:
    """Fetch phase and verify ownership via campaign→project."""
    stmt = (
        select(PhaseORM_V9)
        .where(PhaseORM_V9.id == phase_id)
        .options(
            selectinload(PhaseORM_V9.campaign)
            .selectinload(CampaignORM_V9.project)
        )
    )
    result = await db.execute(stmt)
    phase = result.scalar_one_or_none()
    if not phase:
        raise HTTPException(status_code=404, detail="Phase not found")
    if str(phase.campaign.project.user_id) != user_id:
        raise HTTPException(status_code=403, detail="Not your phase")
    return phase


async def _get_run_owned(
    run_id: UUID,
    user_id: str,
    db: AsyncSession,
) -> RunORM_V9:
    """Fetch run and verify ownership via phase→campaign→project."""
    stmt = (
        select(RunORM_V9)
        .where(RunORM_V9.id == run_id)
        .options(
            selectinload(RunORM_V9.phase)
            .selectinload(PhaseORM_V9.campaign)
            .selectinload(CampaignORM_V9.project)
        )
    )
    result = await db.execute(stmt)
    run = result.scalar_one_or_none()
    if not run:
        raise HTTPException(status_code=404, detail="Run not found")
    if str(run.phase.campaign.project.user_id) != user_id:
        raise HTTPException(status_code=403, detail="Not your run")
    return run


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

    phase = await _get_phase_for_run(phase_id, user_id, db)

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
        if not config.get("smiles_list"):
            raise HTTPException(
                status_code=422,
                detail="config.smiles_list required for import runs (or use import-file endpoint)",
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

    return run


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
    await _get_phase_for_run(phase_id, user_id, db)

    stmt = (
        select(RunORM_V9)
        .where(RunORM_V9.phase_id == phase_id)
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
    return await _get_run_owned(run_id, user_id, db)


@router.post("/runs/{run_id}/cancel", response_model=RunResponse)
async def cancel_run(
    run_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Cancel a running run."""
    run = await _get_run_owned(run_id, user_id, db)
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

    return run


@router.post("/runs/{run_id}/archive", response_model=RunResponse)
async def archive_run(
    run_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Archive a run (soft delete)."""
    run = await _get_run_owned(run_id, user_id, db)
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
    phase = await _get_phase_for_run(phase_id, user_id, db)

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

    # Save uploaded file
    content = await file.read()
    tmp_dir = Path(tempfile.gettempdir()) / "bindx_uploads"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    file_path = tmp_dir / f"{phase_id}_{filename}"
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
        logger.warning("Failed to dispatch Celery task for run %s: %s", run.id, e)
        run.status = "failed"
        run.error_message = f"Failed to dispatch task: {e}"
        await db.flush()

    return run
