"""
BindX V9 — Campaign CRUD router.

4 endpoints: create (under project), list (under project), get, update.
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
from models_v9 import CampaignORM_V9, PhaseORM_V9, ProjectORM_V9
from routers.v9.deps import get_campaign_owned, verify_project_ownership
from schemas_v9 import CampaignCreate, CampaignResponse, CampaignUpdate

logger = logging.getLogger(__name__)

router = APIRouter()


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------

@router.post(
    "/projects/{project_id}/campaigns",
    response_model=CampaignResponse,
    status_code=201,
)
async def create_campaign(
    project_id: UUID,
    body: CampaignCreate,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Create a new campaign under a project."""
    await verify_project_ownership(project_id, user_id, db)

    campaign = CampaignORM_V9(
        project_id=project_id,
        name=body.name,
        pocket_config=body.pocket_config,
    )
    db.add(campaign)
    await db.flush()

    return await get_campaign_owned(campaign.id, user_id, db)


@router.get(
    "/projects/{project_id}/campaigns",
    response_model=list[CampaignResponse],
)
async def list_campaigns(
    project_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """List all campaigns for a project."""
    await verify_project_ownership(project_id, user_id, db)

    stmt = (
        select(CampaignORM_V9)
        .where(CampaignORM_V9.project_id == project_id)
        .options(selectinload(CampaignORM_V9.phases))
        .order_by(CampaignORM_V9.created_at)
    )
    result = await db.execute(stmt)
    return result.scalars().all()


@router.get("/campaigns/{campaign_id}", response_model=CampaignResponse)
async def get_campaign(
    campaign_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Get campaign detail with phases."""
    return await get_campaign_owned(campaign_id, user_id, db)


@router.put("/campaigns/{campaign_id}", response_model=CampaignResponse)
async def update_campaign(
    campaign_id: UUID,
    body: CampaignUpdate,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Update campaign fields (partial update)."""
    campaign = await get_campaign_owned(campaign_id, user_id, db)
    update_data = body.model_dump(exclude_unset=True)
    for field, value in update_data.items():
        setattr(campaign, field, value)
    await db.flush()
    return await get_campaign_owned(campaign_id, user_id, db)
