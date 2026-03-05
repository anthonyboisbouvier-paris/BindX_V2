"""
BindX V9 — Phase CRUD router.

7 endpoints: create (under campaign), list, get (with stats), update,
freeze, unfreeze, delete.

Business rules:
- Duplicate phase type per campaign → 409
- Update on frozen phase → 400
- Unfreeze warns if downstream phases have runs
- Delete blocked if active runs exist (created/running)
"""

from __future__ import annotations

import csv
import io
import logging
from datetime import datetime, timezone
from uuid import UUID

from fastapi import APIRouter, Body, Depends, HTTPException
from fastapi.responses import StreamingResponse
from sqlalchemy import func, select
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.orm import selectinload

from auth_v9 import require_v9_user
from database_v9 import get_v9_db
from models_v9 import (
    MoleculeORM_V9,
    MoleculePropertyORM_V9,
    PhaseORM_V9,
    RunORM_V9,
)
from routers.v9.deps import get_campaign_owned, get_phase_owned
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

    await get_campaign_owned(campaign_id, user_id, db)

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
    await get_campaign_owned(campaign_id, user_id, db)

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
    phase = await get_phase_owned(phase_id, user_id, db)
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
    phase = await get_phase_owned(phase_id, user_id, db)
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
    phase = await get_phase_owned(phase_id, user_id, db)
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
    phase = await get_phase_owned(phase_id, user_id, db)
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


@router.delete("/phases/{phase_id}", status_code=204)
async def delete_phase(
    phase_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Delete a phase and all its runs/molecules (cascade).

    Blocked if any run is still active (created/running).
    """
    phase = await get_phase_owned(phase_id, user_id, db)

    # Block deletion if there are active runs
    active_q = await db.execute(
        select(func.count()).select_from(RunORM_V9).where(
            RunORM_V9.phase_id == phase_id,
            RunORM_V9.status.in_(("created", "running")),
        )
    )
    active_count = active_q.scalar() or 0
    if active_count > 0:
        raise HTTPException(
            status_code=400,
            detail=f"Cannot delete phase with {active_count} active run(s). Cancel them first.",
        )

    await db.delete(phase)
    await db.commit()


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------

@router.post("/phases/{phase_id}/export")
async def export_phase(
    phase_id: UUID,
    body: dict = Body(...),
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Export phase molecules as CSV or SDF.

    Body: {format: "csv"|"sdf", bookmarked_only: bool}
    """
    phase = await get_phase_owned(phase_id, user_id, db)
    fmt = body.get("format", "csv")
    bookmarked_only = body.get("bookmarked_only", False)

    if fmt not in ("csv", "sdf"):
        raise HTTPException(status_code=422, detail="format must be 'csv' or 'sdf'")

    # Fetch molecules with properties
    where = [MoleculeORM_V9.phase_id == phase_id]
    if bookmarked_only:
        where.append(MoleculeORM_V9.bookmarked == True)  # noqa: E712

    stmt = (
        select(MoleculeORM_V9)
        .where(*where)
        .options(selectinload(MoleculeORM_V9.properties))
        .order_by(MoleculeORM_V9.created_at)
    )
    result = await db.execute(stmt)
    mols = result.scalars().all()

    if fmt == "csv":
        return _export_csv(mols, phase)
    else:
        return _export_sdf(mols, phase)


def _export_csv(mols, phase) -> StreamingResponse:
    """Build CSV export with flattened properties."""
    # Collect all property names across molecules
    all_prop_names: set[str] = set()
    mol_props_map = {}
    for mol in mols:
        flat = {}
        for p in (mol.properties or []):
            val = p.property_value
            if isinstance(val, dict):
                # Flatten dict: "admet.logP", "admet.MW", etc.
                for k, v in val.items():
                    key = f"{p.property_name}.{k}" if p.property_name not in ("physicochemical",) else k
                    flat[key] = v
                    all_prop_names.add(key)
            else:
                flat[p.property_name] = val
                all_prop_names.add(p.property_name)
        mol_props_map[mol.id] = flat

    prop_cols = sorted(all_prop_names)
    header = ["smiles", "canonical_smiles", "name", "bookmarked"] + prop_cols

    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow(header)
    for mol in mols:
        flat = mol_props_map.get(mol.id, {})
        row = [mol.smiles, mol.canonical_smiles, mol.name or "", mol.bookmarked]
        for col in prop_cols:
            row.append(flat.get(col, ""))
        writer.writerow(row)

    buf.seek(0)
    filename = f"phase_{phase.type}_export.csv"
    return StreamingResponse(
        iter([buf.getvalue()]),
        media_type="text/csv",
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


def _export_sdf(mols, phase) -> StreamingResponse:
    """Build SDF export using RDKit if available, fallback to SMILES list."""
    buf = io.StringIO()

    try:
        from rdkit import Chem

        for mol_orm in mols:
            mol = Chem.MolFromSmiles(mol_orm.smiles)
            if mol is None:
                continue
            block = Chem.MolToMolBlock(mol)
            buf.write(block)
            buf.write(f">  <Name>\n{mol_orm.name or ''}\n\n")
            buf.write(f">  <SMILES>\n{mol_orm.smiles}\n\n")
            buf.write(f">  <Bookmarked>\n{mol_orm.bookmarked}\n\n")
            # Add properties
            for p in (mol_orm.properties or []):
                val = p.property_value
                if isinstance(val, dict):
                    for k, v in val.items():
                        buf.write(f">  <{p.property_name}_{k}>\n{v}\n\n")
                else:
                    buf.write(f">  <{p.property_name}>\n{val}\n\n")
            buf.write("$$$$\n")
    except ImportError:
        # Fallback: just write SMILES
        for mol_orm in mols:
            buf.write(f"{mol_orm.smiles} {mol_orm.name or ''}\n")

    buf.seek(0)
    filename = f"phase_{phase.type}_export.sdf"
    return StreamingResponse(
        iter([buf.getvalue()]),
        media_type="chemical/x-mdl-sdfile",
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )
