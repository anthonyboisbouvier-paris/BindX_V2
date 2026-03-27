"""
BindX V9 — AFVS Router (BDX-41)

3 endpoints for Adaptive Focused Virtual Screening:
  POST   /phases/{phase_id}/runs/afvs   — launch AFVS run
  DELETE /runs/{run_id}/afvs            — cancel AFVS run
  GET    /runs/{run_id}/afvs/status     — detailed status
"""

from __future__ import annotations

import logging
import os
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.orm import selectinload

from auth_v9 import require_v9_user
from database_v9 import get_v9_db
from models_v9 import (
    AFVSJobORM_V9,
    CampaignORM_V9,
    PhaseORM_V9,
    ProjectORM_V9,
    RunORM_V9,
)
from routers.v9.deps import get_phase_owned, get_run_owned
from schemas_v9 import AFVSRunCreate, AFVSStatusResponse

logger = logging.getLogger(__name__)

router = APIRouter()


# ---------------------------------------------------------------------------
# POST /phases/{phase_id}/runs/afvs — Launch AFVS run
# ---------------------------------------------------------------------------

@router.post(
    "/phases/{phase_id}/runs/afvs",
    response_model=AFVSStatusResponse,
    status_code=201,
)
async def launch_afvs(
    phase_id: UUID,
    body: AFVSRunCreate,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Launch an AFVS ultra-large screening run."""

    # Check EC2 is configured
    if not os.environ.get("AFVS_EC2_HOST"):
        raise HTTPException(
            status_code=503,
            detail="AFVS infrastructure not configured. Set AFVS_EC2_HOST environment variable.",
        )

    # Fetch phase with campaign + project
    phase = await get_phase_owned(phase_id, user_id, db)

    if phase.status == "frozen":
        raise HTTPException(status_code=400, detail="Cannot create runs on a frozen phase")

    # AFVS only in Phase A (hit_discovery)
    if phase.type != "hit_discovery":
        raise HTTPException(
            status_code=400,
            detail="AFVS runs are only available in Phase A (Hit Discovery).",
        )

    # Get project to check receptor
    campaign = phase.campaign
    stmt_proj = select(ProjectORM_V9).where(ProjectORM_V9.id == campaign.project_id)
    result = await db.execute(stmt_proj)
    project = result.scalar_one_or_none()

    if not project or not project.receptor_prep_report:
        raise HTTPException(
            status_code=400,
            detail="Receptor not prepared. Go to Target Setup and prepare the receptor before launching AFVS.",
        )

    # Check pocket is configured
    pocket_config = campaign.pocket_config or {}
    pocket_center = pocket_config.get("center")
    if not pocket_center:
        pockets = project.pockets_detected or []
        tp = project.target_preview or {}
        selected_idx = tp.get("selected_pocket_index", 0)
        if pockets and selected_idx < len(pockets):
            pocket_center = pockets[selected_idx].get("center")
        elif pockets:
            pocket_center = pockets[0].get("center")

    if not pocket_center:
        raise HTTPException(
            status_code=400,
            detail="No pocket configured. Configure a pocket in Campaign settings.",
        )

    # Build docking box from pocket
    pocket_size = pocket_config.get("size", [22, 22, 22])
    docking_box = {
        "center_x": pocket_center[0],
        "center_y": pocket_center[1],
        "center_z": pocket_center[2],
        "size_x": pocket_size[0] if isinstance(pocket_size, list) else 22,
        "size_y": pocket_size[1] if isinstance(pocket_size, list) else 22,
        "size_z": pocket_size[2] if isinstance(pocket_size, list) else 22,
    }

    # Create Run (type=afvs)
    config = {
        "strategy": body.strategy,
        "docking_scenario": body.docking_scenario,
        "exhaustiveness_prescreen": body.exhaustiveness_prescreen,
        "exhaustiveness_primary": body.exhaustiveness_primary,
        "reps_per_tranche": body.reps_per_tranche,
        "tranche_pct": body.tranche_pct,
        "mw_max": body.mw_max,
        "logp_max": body.logp_max,
        "hbd_max": body.hbd_max,
        "hba_max": body.hba_max,
        "top_n_results": body.top_n_results,
        "budget_cap_usd": body.budget_cap_usd,
    }

    run = RunORM_V9(
        phase_id=phase_id,
        type="afvs",
        config=config,
        input_source="afvs",
        status="created",
    )
    db.add(run)
    await db.flush()  # get run.id

    # Create AFVSJob
    afvs_job = AFVSJobORM_V9(
        run_id=run.id,
        receptor_s3_path="pending",  # will be set by Celery task
        docking_box=docking_box,
        docking_scenario=body.docking_scenario,
        top_n_results=body.top_n_results,
        budget_cap_usd=body.budget_cap_usd,
        s3_output_prefix=f"runs/{run.id}",
    )
    db.add(afvs_job)
    await db.commit()

    # Dispatch Celery task
    try:
        from celery_app import celery_app
        task = celery_app.send_task("tasks_v9.execute_run", args=[str(run.id)])
        run.config = {**config, "_celery_task_id": task.id}
        await db.commit()
        logger.info("Dispatched AFVS task %s for run %s", task.id, run.id)
    except Exception as e:
        logger.warning("Failed to dispatch AFVS task for run %s: %s", run.id, e)
        run.status = "failed"
        run.error_message = f"Failed to dispatch task: {e}"
        await db.commit()
        raise HTTPException(status_code=503, detail="Task dispatch failed. Please try again.")

    return AFVSStatusResponse(
        run_id=run.id,
        status=run.status,
        afvs_phase=afvs_job.afvs_phase,
        progress=0,
        current_step="Queued",
        budget_cap_usd=body.budget_cap_usd,
    )


# ---------------------------------------------------------------------------
# DELETE /runs/{run_id}/afvs — Cancel AFVS run
# ---------------------------------------------------------------------------

@router.delete("/runs/{run_id}/afvs")
async def cancel_afvs(
    run_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Cancel a running AFVS screening job."""
    run = await get_run_owned(run_id, user_id, db)

    if run.type != "afvs":
        raise HTTPException(status_code=400, detail="Not an AFVS run")

    if run.status not in ("created", "running"):
        raise HTTPException(
            status_code=400, detail=f"Cannot cancel AFVS run with status '{run.status}'"
        )

    run.status = "cancelled"

    # Update AFVS job phase
    stmt = select(AFVSJobORM_V9).where(AFVSJobORM_V9.run_id == run_id)
    result = await db.execute(stmt)
    afvs_job = result.scalar_one_or_none()
    if afvs_job:
        afvs_job.afvs_phase = "failed"

    # Try to revoke Celery task
    celery_task_id = (run.config or {}).get("_celery_task_id")
    if celery_task_id:
        try:
            from celery_app import celery_app
            celery_app.control.revoke(celery_task_id, terminate=True, signal="SIGTERM")
        except Exception as e:
            logger.warning("Failed to revoke AFVS Celery task %s: %s", celery_task_id, e)

    await db.commit()
    return {"status": "cancelled", "run_id": str(run_id)}


# ---------------------------------------------------------------------------
# GET /runs/{run_id}/afvs/status — Detailed AFVS status
# ---------------------------------------------------------------------------

@router.get("/runs/{run_id}/afvs/status", response_model=AFVSStatusResponse)
async def get_afvs_status(
    run_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Get detailed status of an AFVS run including job metadata."""
    run = await get_run_owned(run_id, user_id, db)

    if run.type != "afvs":
        raise HTTPException(status_code=400, detail="Not an AFVS run")

    # Fetch AFVS job details
    stmt = select(AFVSJobORM_V9).where(AFVSJobORM_V9.run_id == run_id)
    result = await db.execute(stmt)
    afvs_job = result.scalar_one_or_none()

    return AFVSStatusResponse(
        run_id=run.id,
        status=run.status,
        afvs_phase=afvs_job.afvs_phase if afvs_job else None,
        progress=run.progress,
        current_step=run.current_step,
        molecules_screened=afvs_job.molecules_screened if afvs_job else None,
        molecules_imported=afvs_job.molecules_imported if afvs_job else None,
        actual_cost_usd=afvs_job.actual_cost_usd if afvs_job else None,
        budget_cap_usd=afvs_job.budget_cap_usd if afvs_job else None,
        started_at=run.started_at,
        completed_at=run.completed_at,
        error_message=run.error_message,
    )
