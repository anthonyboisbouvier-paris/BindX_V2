"""
BindX V9 — Shared ownership verification helpers.

All V9 routers use these functions to check entity ownership
via the chain: entity → ... → project.user_id == current_user.
"""

from __future__ import annotations

from uuid import UUID

from fastapi import HTTPException
from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.orm import selectinload

from models_v9 import (
    CampaignORM_V9,
    MoleculeORM_V9,
    MoleculePropertyORM_V9,
    PhaseORM_V9,
    ProjectORM_V9,
    RunORM_V9,
)


async def get_project_owned(
    project_id: UUID, user_id: str, db: AsyncSession,
) -> ProjectORM_V9:
    """Fetch project with campaigns→phases, verify ownership."""
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


async def verify_project_ownership(
    project_id: UUID, user_id: str, db: AsyncSession,
) -> ProjectORM_V9:
    """Lightweight project ownership check (no eager loading)."""
    stmt = select(ProjectORM_V9).where(ProjectORM_V9.id == project_id)
    result = await db.execute(stmt)
    project = result.scalar_one_or_none()
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    if str(project.user_id) != user_id:
        raise HTTPException(status_code=403, detail="Not your project")
    return project


async def get_campaign_owned(
    campaign_id: UUID, user_id: str, db: AsyncSession,
) -> CampaignORM_V9:
    """Fetch campaign with phases, verify ownership via project."""
    stmt = (
        select(CampaignORM_V9)
        .where(CampaignORM_V9.id == campaign_id)
        .options(
            selectinload(CampaignORM_V9.project),
            selectinload(CampaignORM_V9.phases),
        )
    )
    result = await db.execute(stmt)
    campaign = result.scalar_one_or_none()
    if not campaign:
        raise HTTPException(status_code=404, detail="Campaign not found")
    if str(campaign.project.user_id) != user_id:
        raise HTTPException(status_code=403, detail="Not your campaign")
    return campaign


async def get_phase_owned(
    phase_id: UUID, user_id: str, db: AsyncSession,
) -> PhaseORM_V9:
    """Fetch phase, verify ownership via campaign→project."""
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


async def get_run_owned(
    run_id: UUID, user_id: str, db: AsyncSession,
) -> RunORM_V9:
    """Fetch run, verify ownership via phase→campaign→project."""
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


async def get_molecule_owned(
    molecule_id: UUID, user_id: str, db: AsyncSession,
) -> MoleculeORM_V9:
    """Fetch molecule with properties, verify ownership via phase→campaign→project."""
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
