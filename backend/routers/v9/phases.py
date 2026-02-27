"""
BindX V9 — Phase CRUD router.

6 endpoints: create (under campaign), list, get (with stats), update,
freeze, unfreeze.

Business rules:
- Duplicate phase type per campaign → 409
- Update on frozen phase → 400
- Unfreeze warns if downstream phases have runs
"""

from __future__ import annotations

import logging
from datetime import datetime, timezone
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy import func, select
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.orm import selectinload

from auth_v9 import require_v9_user
from database_v9 import get_v9_db
from models_v9 import (
    CampaignORM_V9,
    MoleculeORM_V9,
    PhaseORM_V9,
    ProjectORM_V9,
    RunORM_V9,
)
from schemas_v9 import (
    FreezeResponse,
    PhaseCreate,
    PhaseResponse,
    PhaseSummary,
    PhaseUpdate,
    UnfreezeResponse,
)

logger = logging.getLogger(__name__)

router = APIRouter()

# Phase type ordering for downstream checks
PHASE_ORDER = {
    "hit_discovery": 0,
    "hit_to_lead": 1,
    "lead_optimization": 2,
}

VALID_PHASE_TYPES = set(PHASE_ORDER.keys())


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

async def _verify_campaign_ownership(
    campaign_id: UUID,
    user_id: str,
    db: AsyncSession,
) -> CampaignORM_V9:
    """Fetch campaign and verify ownership via parent project."""
    stmt = (
        select(CampaignORM_V9)
        .where(CampaignORM_V9.id == campaign_id)
        .options(selectinload(CampaignORM_V9.project))
    )
    result = await db.execute(stmt)
    campaign = result.scalar_one_or_none()
    if not campaign:
        raise HTTPException(status_code=404, detail="Campaign not found")
    if str(campaign.project.user_id) != user_id:
        raise HTTPException(status_code=403, detail="Not your campaign")
    return campaign


async def _get_phase_owned(
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


async def _phase_stats(phase_id: UUID, db: AsyncSession) -> dict:
    """Compute molecule_count, bookmarked_count, run_count for a phase."""
    mol_q = await db.execute(
        select(func.count()).select_from(MoleculeORM_V9).where(
            MoleculeORM_V9.phase_id == phase_id
        )
    )
    mol_count = mol_q.scalar() or 0

    bm_q = await db.execute(
        select(func.count()).select_from(MoleculeORM_V9).where(
            MoleculeORM_V9.phase_id == phase_id,
            MoleculeORM_V9.bookmarked == True,  # noqa: E712
        )
    )
    bm_count = bm_q.scalar() or 0

    run_q = await db.execute(
        select(func.count()).select_from(RunORM_V9).where(
            RunORM_V9.phase_id == phase_id
        )
    )
    run_count = run_q.scalar() or 0

    return {
        "molecule_count": mol_count,
        "bookmarked_count": bm_count,
        "run_count": run_count,
    }


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------

@router.post(
    "/campaigns/{campaign_id}/phases",
    response_model=PhaseResponse,
    status_code=201,
)
async def create_phase(
    campaign_id: UUID,
    body: PhaseCreate,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Create a new phase. Rejects duplicate type per campaign (409)."""
    if body.type not in VALID_PHASE_TYPES:
        raise HTTPException(
            status_code=422,
            detail=f"Invalid phase type. Must be one of: {', '.join(VALID_PHASE_TYPES)}",
        )

    campaign = await _verify_campaign_ownership(campaign_id, user_id, db)

    # Check for duplicate type
    dup = await db.execute(
        select(PhaseORM_V9).where(
            PhaseORM_V9.campaign_id == campaign_id,
            PhaseORM_V9.type == body.type,
        )
    )
    if dup.scalar_one_or_none():
        raise HTTPException(
            status_code=409,
            detail=f"Phase of type '{body.type}' already exists in this campaign",
        )

    phase = PhaseORM_V9(
        campaign_id=campaign_id,
        type=body.type,
        status="active",
    )
    db.add(phase)
    await db.flush()

    return PhaseResponse(
        id=phase.id,
        campaign_id=phase.campaign_id,
        type=phase.type,
        status=phase.status,
        frozen_at=phase.frozen_at,
        created_at=phase.created_at,
        molecule_count=0,
        bookmarked_count=0,
        run_count=0,
    )


@router.get(
    "/campaigns/{campaign_id}/phases",
    response_model=list[PhaseSummary],
)
async def list_phases(
    campaign_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """List all phases for a campaign."""
    await _verify_campaign_ownership(campaign_id, user_id, db)

    stmt = (
        select(PhaseORM_V9)
        .where(PhaseORM_V9.campaign_id == campaign_id)
        .order_by(PhaseORM_V9.created_at)
    )
    result = await db.execute(stmt)
    return result.scalars().all()


@router.get("/phases/{phase_id}", response_model=PhaseResponse)
async def get_phase(
    phase_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Get phase detail with computed stats."""
    phase = await _get_phase_owned(phase_id, user_id, db)
    stats = await _phase_stats(phase_id, db)
    return PhaseResponse(
        id=phase.id,
        campaign_id=phase.campaign_id,
        type=phase.type,
        status=phase.status,
        frozen_at=phase.frozen_at,
        created_at=phase.created_at,
        **stats,
    )


@router.put("/phases/{phase_id}", response_model=PhaseSummary)
async def update_phase(
    phase_id: UUID,
    body: PhaseUpdate,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Update phase. Rejects if frozen (400)."""
    phase = await _get_phase_owned(phase_id, user_id, db)
    if phase.status == "frozen":
        raise HTTPException(status_code=400, detail="Cannot update a frozen phase")

    update_data = body.model_dump(exclude_unset=True)
    for field, value in update_data.items():
        setattr(phase, field, value)
    await db.flush()
    return phase


@router.post("/phases/{phase_id}/freeze", response_model=FreezeResponse)
async def freeze_phase(
    phase_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Freeze a phase: set status=frozen, frozen_at=now(). Returns bookmarked count."""
    phase = await _get_phase_owned(phase_id, user_id, db)
    if phase.status == "frozen":
        raise HTTPException(status_code=400, detail="Phase is already frozen")

    phase.status = "frozen"
    phase.frozen_at = datetime.now(timezone.utc)
    await db.flush()

    # Count bookmarked molecules
    bm_q = await db.execute(
        select(func.count()).select_from(MoleculeORM_V9).where(
            MoleculeORM_V9.phase_id == phase_id,
            MoleculeORM_V9.bookmarked == True,  # noqa: E712
        )
    )
    bm_count = bm_q.scalar() or 0

    return FreezeResponse(
        id=phase.id,
        status=phase.status,
        frozen_at=phase.frozen_at,
        bookmarked_molecule_count=bm_count,
    )


@router.post("/phases/{phase_id}/unfreeze", response_model=UnfreezeResponse)
async def unfreeze_phase(
    phase_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Unfreeze a phase. Warns if downstream phases have runs."""
    phase = await _get_phase_owned(phase_id, user_id, db)
    if phase.status != "frozen":
        raise HTTPException(status_code=400, detail="Phase is not frozen")

    # Check downstream phases for runs
    warning = None
    current_order = PHASE_ORDER.get(phase.type, -1)
    downstream_types = [t for t, o in PHASE_ORDER.items() if o > current_order]

    if downstream_types:
        # Find sibling phases of higher type in the same campaign
        sibling_q = await db.execute(
            select(PhaseORM_V9.id).where(
                PhaseORM_V9.campaign_id == phase.campaign_id,
                PhaseORM_V9.type.in_(downstream_types),
            )
        )
        sibling_ids = [row[0] for row in sibling_q.fetchall()]

        if sibling_ids:
            run_count_q = await db.execute(
                select(func.count()).select_from(RunORM_V9).where(
                    RunORM_V9.phase_id.in_(sibling_ids)
                )
            )
            downstream_runs = run_count_q.scalar() or 0
            if downstream_runs > 0:
                warning = (
                    f"Warning: {downstream_runs} run(s) exist in downstream phases. "
                    "Unfreezing may invalidate their results."
                )

    phase.status = "active"
    phase.frozen_at = None
    await db.flush()

    return UnfreezeResponse(
        id=phase.id,
        status=phase.status,
        warning=warning,
    )
