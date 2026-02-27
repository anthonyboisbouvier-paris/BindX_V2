"""
BindX V9 — Molecule router.

4 endpoints: list molecules (paginated), get molecule, bookmark, bookmark batch.
"""

from __future__ import annotations

import logging
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy import func, select
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.orm import selectinload

from auth_v9 import require_v9_user
from database_v9 import get_v9_db
from models_v9 import (
    CampaignORM_V9,
    MoleculeORM_V9,
    MoleculePropertyORM_V9,
    PhaseORM_V9,
    ProjectORM_V9,
)
from schemas_v9 import BookmarkRequest, MoleculeListResponse, MoleculeResponse

logger = logging.getLogger(__name__)

router = APIRouter()

VALID_SORT_FIELDS = {
    "created_at", "smiles", "name", "bookmarked", "canonical_smiles",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

async def _verify_phase_ownership(
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


async def _get_molecule_owned(
    molecule_id: UUID,
    user_id: str,
    db: AsyncSession,
) -> MoleculeORM_V9:
    """Fetch molecule and verify ownership via phase→campaign→project."""
    stmt = (
        select(MoleculeORM_V9)
        .where(MoleculeORM_V9.id == molecule_id)
        .options(
            selectinload(MoleculeORM_V9.properties),
            selectinload(MoleculeORM_V9.phase)
            .selectinload(PhaseORM_V9.campaign)
            .selectinload(CampaignORM_V9.project),
        )
    )
    result = await db.execute(stmt)
    mol = result.scalar_one_or_none()
    if not mol:
        raise HTTPException(status_code=404, detail="Molecule not found")
    if str(mol.phase.campaign.project.user_id) != user_id:
        raise HTTPException(status_code=403, detail="Not your molecule")
    return mol


def _mol_to_response(mol: MoleculeORM_V9) -> MoleculeResponse:
    """Convert ORM molecule + loaded properties to response schema."""
    props = {}
    if mol.properties:
        for p in mol.properties:
            props[p.property_name] = p.property_value
    return MoleculeResponse(
        id=mol.id,
        phase_id=mol.phase_id,
        smiles=mol.smiles,
        canonical_smiles=mol.canonical_smiles,
        name=mol.name,
        source_run_id=mol.source_run_id,
        bookmarked=mol.bookmarked,
        generation_level=mol.generation_level,
        parent_molecule_id=mol.parent_molecule_id,
        ai_generated=mol.ai_generated,
        created_at=mol.created_at,
        properties=props if props else None,
    )


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------

@router.get(
    "/phases/{phase_id}/molecules",
    response_model=MoleculeListResponse,
)
async def list_molecules(
    phase_id: UUID,
    sort_by: str = Query("created_at", description="Sort field"),
    sort_dir: str = Query("desc", regex="^(asc|desc)$"),
    offset: int = Query(0, ge=0),
    limit: int = Query(50, ge=1, le=500),
    bookmarked_only: bool = Query(False),
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """List molecules for a phase with properties, sorted and paginated."""
    await _verify_phase_ownership(phase_id, user_id, db)

    # Sanitize sort field
    if sort_by not in VALID_SORT_FIELDS:
        sort_by = "created_at"

    # Total count
    count_stmt = (
        select(func.count())
        .select_from(MoleculeORM_V9)
        .where(MoleculeORM_V9.phase_id == phase_id)
    )
    if bookmarked_only:
        count_stmt = count_stmt.where(MoleculeORM_V9.bookmarked == True)  # noqa: E712
    total_result = await db.execute(count_stmt)
    total = total_result.scalar() or 0

    # Query with properties
    sort_col = getattr(MoleculeORM_V9, sort_by, MoleculeORM_V9.created_at)
    order = sort_col.desc() if sort_dir == "desc" else sort_col.asc()

    stmt = (
        select(MoleculeORM_V9)
        .where(MoleculeORM_V9.phase_id == phase_id)
        .options(selectinload(MoleculeORM_V9.properties))
        .order_by(order)
        .offset(offset)
        .limit(limit)
    )
    if bookmarked_only:
        stmt = stmt.where(MoleculeORM_V9.bookmarked == True)  # noqa: E712

    result = await db.execute(stmt)
    mols = result.scalars().all()

    return MoleculeListResponse(
        molecules=[_mol_to_response(m) for m in mols],
        total=total,
        offset=offset,
        limit=limit,
    )


@router.get("/molecules/{molecule_id}", response_model=MoleculeResponse)
async def get_molecule(
    molecule_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Get molecule detail with all properties."""
    mol = await _get_molecule_owned(molecule_id, user_id, db)
    return _mol_to_response(mol)


@router.put("/molecules/{molecule_id}/bookmark", response_model=MoleculeResponse)
async def bookmark_molecule(
    molecule_id: UUID,
    bookmarked: bool = Query(..., description="Set bookmark state"),
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Toggle bookmark on a single molecule."""
    mol = await _get_molecule_owned(molecule_id, user_id, db)
    mol.bookmarked = bookmarked
    await db.flush()
    return _mol_to_response(mol)


@router.post("/phases/{phase_id}/molecules/bookmark-batch")
async def bookmark_batch(
    phase_id: UUID,
    body: BookmarkRequest,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Bookmark or unbookmark multiple molecules at once."""
    await _verify_phase_ownership(phase_id, user_id, db)

    # Fetch all target molecules in this phase
    stmt = (
        select(MoleculeORM_V9)
        .where(
            MoleculeORM_V9.phase_id == phase_id,
            MoleculeORM_V9.id.in_(body.molecule_ids),
        )
    )
    result = await db.execute(stmt)
    mols = result.scalars().all()

    updated = 0
    for mol in mols:
        mol.bookmarked = body.bookmarked
        updated += 1

    await db.flush()
    return {"updated": updated, "bookmarked": body.bookmarked}
