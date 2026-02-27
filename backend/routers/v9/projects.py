"""
BindX V9 — Project CRUD router.

5 endpoints: create, list, get, update, delete.
Auto-creates a "Default Campaign" on project creation.
"""

from __future__ import annotations

import logging
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.orm import selectinload

from auth_v9 import require_v9_user
from database_v9 import get_v9_db
from models_v9 import CampaignORM_V9, ProjectORM_V9
from schemas_v9 import (
    ProjectCreate,
    ProjectListItem,
    ProjectResponse,
    ProjectUpdate,
)

logger = logging.getLogger(__name__)

router = APIRouter()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

async def _get_project_owned(
    project_id: UUID,
    user_id: str,
    db: AsyncSession,
) -> ProjectORM_V9:
    """Fetch a project with campaigns→phases, verify ownership."""
    stmt = (
        select(ProjectORM_V9)
        .where(ProjectORM_V9.id == project_id)
        .options(
            selectinload(ProjectORM_V9.campaigns)
            .selectinload(CampaignORM_V9.phases)
        )
    )
    result = await db.execute(stmt)
    project = result.scalar_one_or_none()
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    if str(project.user_id) != user_id:
        raise HTTPException(status_code=403, detail="Not your project")
    return project


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
    proj = await _get_project_owned(project.id, user_id, db)
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
    return result.scalars().all()


@router.get("/projects/{project_id}", response_model=ProjectResponse)
async def get_project(
    project_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Get project detail with campaigns and phases."""
    return await _get_project_owned(project_id, user_id, db)


@router.put("/projects/{project_id}", response_model=ProjectResponse)
async def update_project(
    project_id: UUID,
    body: ProjectUpdate,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Update project fields (partial update)."""
    project = await _get_project_owned(project_id, user_id, db)
    update_data = body.model_dump(exclude_unset=True)
    for field, value in update_data.items():
        setattr(project, field, value)
    await db.flush()
    # Reload
    return await _get_project_owned(project_id, user_id, db)


@router.delete("/projects/{project_id}", status_code=204)
async def delete_project(
    project_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Delete a project and all its children (cascade)."""
    project = await _get_project_owned(project_id, user_id, db)
    await db.delete(project)
