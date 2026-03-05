"""
BindX V9 — Project CRUD router.

5 endpoints: create, list, get, update, delete.
Auto-creates a "Default Campaign" on project creation.
"""

from __future__ import annotations

import asyncio
import logging
from uuid import UUID

from fastapi import APIRouter, Body, Depends, HTTPException
from sqlalchemy import func, select
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.orm import selectinload

from auth_v9 import require_v9_user
from database_v9 import get_v9_db
from models_v9 import CampaignORM_V9, MoleculeORM_V9, PhaseORM_V9, ProjectORM_V9, RunORM_V9
from routers.v9.deps import get_project_owned
from schemas_v9 import (
    ProjectCreate,
    ProjectListItem,
    ProjectResponse,
    ProjectUpdate,
)

logger = logging.getLogger(__name__)

router = APIRouter()


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------

@router.post("/projects", response_model=ProjectResponse, status_code=201)
async def create_project(
    body: ProjectCreate,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Create a new project with a default campaign."""
    project = ProjectORM_V9(
        user_id=user_id,
        name=body.name,
        description=body.description,
        target_input_type=body.target_input_type,
        target_input_value=body.target_input_value,
        target_name=body.target_name,
        target_pdb_id=body.target_pdb_id,
    )
    db.add(project)
    await db.flush()

    # Auto-create default campaign
    default_campaign = CampaignORM_V9(
        project_id=project.id,
        name="Default Campaign",
    )
    db.add(default_campaign)
    await db.flush()

    # Reload with relationships
    proj = await get_project_owned(project.id, user_id, db)
    return proj


@router.get("/projects", response_model=list[ProjectListItem])
async def list_projects(
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """List all projects for the authenticated user."""
    stmt = (
        select(ProjectORM_V9)
        .where(ProjectORM_V9.user_id == user_id)
        .options(
            selectinload(ProjectORM_V9.campaigns)
            .selectinload(CampaignORM_V9.phases)
        )
        .order_by(ProjectORM_V9.updated_at.desc())
    )
    result = await db.execute(stmt)
    projects = result.scalars().all()

    # Compute per-phase stats (molecules count, bookmarked, runs)
    all_phase_ids = []
    for p in projects:
        for c in p.campaigns:
            for ph in c.phases:
                all_phase_ids.append(ph.id)

    if all_phase_ids:
        # Molecule counts per phase
        mol_stats_stmt = (
            select(
                MoleculeORM_V9.phase_id,
                func.count().label("total"),
                func.count().filter(MoleculeORM_V9.bookmarked == True).label("bookmarked"),  # noqa: E712
            )
            .where(MoleculeORM_V9.phase_id.in_(all_phase_ids))
            .group_by(MoleculeORM_V9.phase_id)
        )
        mol_result = await db.execute(mol_stats_stmt)
        mol_stats = {row.phase_id: {"total_molecules": row.total, "bookmarked": row.bookmarked} for row in mol_result.all()}

        # Run counts per phase
        run_stats_stmt = (
            select(
                RunORM_V9.phase_id,
                func.count().filter(RunORM_V9.status == "completed").label("runs_completed"),
                func.count().filter(RunORM_V9.status.in_(["running", "created"])).label("runs_running"),
            )
            .where(RunORM_V9.phase_id.in_(all_phase_ids))
            .group_by(RunORM_V9.phase_id)
        )
        run_result = await db.execute(run_stats_stmt)
        run_stats = {row.phase_id: {"runs_completed": row.runs_completed, "runs_running": row.runs_running} for row in run_result.all()}

        # Inject stats into phases
        for p in projects:
            for c in p.campaigns:
                for ph in c.phases:
                    ms = mol_stats.get(ph.id, {})
                    rs = run_stats.get(ph.id, {})
                    ph.stats = {
                        "total_molecules": ms.get("total_molecules", 0),
                        "bookmarked": ms.get("bookmarked", 0),
                        "runs_completed": rs.get("runs_completed", 0),
                        "runs_running": rs.get("runs_running", 0),
                    }

    return projects


@router.get("/projects/{project_id}", response_model=ProjectResponse)
async def get_project(
    project_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Get project detail with campaigns and phases."""
    return await get_project_owned(project_id, user_id, db)


@router.put("/projects/{project_id}", response_model=ProjectResponse)
async def update_project(
    project_id: UUID,
    body: ProjectUpdate,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Update project fields (partial update)."""
    project = await get_project_owned(project_id, user_id, db)
    update_data = body.model_dump(exclude_unset=True)
    for field, value in update_data.items():
        setattr(project, field, value)
    await db.flush()
    # Reload
    return await get_project_owned(project_id, user_id, db)


@router.delete("/projects/{project_id}", status_code=204)
async def delete_project(
    project_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Delete a project and all its children (cascade)."""
    project = await get_project_owned(project_id, user_id, db)
    await db.delete(project)


@router.post("/projects/{project_id}/preview-target")
async def preview_target_v9(
    project_id: UUID,
    body: dict = Body(...),
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """V9-native target preview. Delegates to the V8 legacy logic but scoped to project ownership."""
    await get_project_owned(project_id, user_id, db)

    uniprot_id = body.get("uniprot_id", "").strip().upper()
    if not uniprot_id:
        raise HTTPException(status_code=400, detail="uniprot_id is required")

    # Reuse the V8 preview-target logic
    from routers.v8_legacy import preview_target
    return await preview_target({"uniprot_id": uniprot_id})
