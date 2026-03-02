"""
BindX V9 — Molecule router.

Endpoints: list (cursor+offset), stats, get, bookmark, bookmark batch.
"""

from __future__ import annotations

import base64
import json
import logging
from datetime import datetime
from typing import Optional
from uuid import UUID

from fastapi import APIRouter, Body, Depends, HTTPException, Query
from sqlalchemy import case, func, literal, or_, select, text, tuple_
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.orm import selectinload

from auth_v9 import require_v9_user
from database_v9 import get_v9_db
from models_v9 import (
    MoleculeORM_V9,
    MoleculePropertyORM_V9,
)
from routers.v9.deps import get_molecule_owned, get_phase_owned
from schemas_v9 import (
    BookmarkRequest,
    MoleculeListResponse,
    MoleculeResponse,
    MoleculeStatsResponse,
)

logger = logging.getLogger(__name__)

router = APIRouter()

# Direct columns on the molecules table that can be sorted
VALID_SORT_FIELDS = {
    "created_at", "smiles", "name", "bookmarked", "canonical_smiles",
}

# Property-based sort fields — resolved via lateral join on molecule_properties
PROPERTY_SORT_FIELDS = {
    "docking_score", "QED", "logP", "MW", "TPSA", "HBD", "HBA",
    "cnn_score", "cnn_affinity", "composite_score",
}


# ---------------------------------------------------------------------------
# Cursor helpers
# ---------------------------------------------------------------------------

def _encode_cursor(sort_value, mol_id: UUID) -> str:
    """Encode (sort_value, id) into an opaque base64 cursor."""
    payload = json.dumps({"v": sort_value, "id": str(mol_id)})
    return base64.urlsafe_b64encode(payload.encode()).decode()


def _decode_cursor(cursor: str) -> tuple:
    """Decode cursor → (sort_value, UUID)."""
    try:
        payload = json.loads(base64.urlsafe_b64decode(cursor))
        return payload["v"], UUID(payload["id"])
    except Exception:
        raise HTTPException(status_code=400, detail="Invalid cursor")


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
    ai_generated_only: bool = Query(False),
    search: Optional[str] = Query(None, description="Search name or SMILES"),
    cursor: Optional[str] = Query(None, description="Keyset pagination cursor"),
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """List molecules for a phase with properties, sorted and paginated.

    Supports keyset cursor pagination (preferred for large datasets) and
    offset-based pagination (fallback).
    """
    await get_phase_owned(phase_id, user_id, db)

    # --- Base filter ---
    base_where = [MoleculeORM_V9.phase_id == phase_id]
    if bookmarked_only:
        base_where.append(MoleculeORM_V9.bookmarked == True)  # noqa: E712
    if ai_generated_only:
        base_where.append(MoleculeORM_V9.ai_generated == True)  # noqa: E712
    if search:
        pattern = f"%{search}%"
        base_where.append(
            or_(
                MoleculeORM_V9.name.ilike(pattern),
                MoleculeORM_V9.canonical_smiles.ilike(pattern),
            )
        )

    # --- Total count (with filters applied) ---
    count_stmt = (
        select(func.count())
        .select_from(MoleculeORM_V9)
        .where(*base_where)
    )
    total_result = await db.execute(count_stmt)
    total = total_result.scalar() or 0

    # --- Determine sort column ---
    is_property_sort = sort_by in PROPERTY_SORT_FIELDS
    if not is_property_sort and sort_by not in VALID_SORT_FIELDS:
        sort_by = "created_at"

    if is_property_sort:
        # Lateral join to get the property value for sorting
        prop_sub = (
            select(
                MoleculePropertyORM_V9.property_value["value"].as_float().label("prop_val")
            )
            .where(
                MoleculePropertyORM_V9.molecule_id == MoleculeORM_V9.id,
                MoleculePropertyORM_V9.property_name == sort_by,
            )
            .limit(1)
            .correlate(MoleculeORM_V9)
            .scalar_subquery()
        )
        sort_expr = prop_sub
        if sort_dir == "desc":
            order = sort_expr.desc().nullslast()
        else:
            order = sort_expr.asc().nullslast()

        # For property sorts, fall back to offset pagination (cursor is complex with subquery)
        stmt = (
            select(MoleculeORM_V9)
            .where(*base_where)
            .options(selectinload(MoleculeORM_V9.properties))
            .order_by(order, MoleculeORM_V9.id)
            .offset(offset)
            .limit(limit + 1)  # fetch one extra to detect has_more
        )
    else:
        # Direct column sort — supports cursor pagination
        sort_col = getattr(MoleculeORM_V9, sort_by, MoleculeORM_V9.created_at)

        if sort_dir == "desc":
            order = sort_col.desc()
        else:
            order = sort_col.asc()

        stmt = (
            select(MoleculeORM_V9)
            .where(*base_where)
            .options(selectinload(MoleculeORM_V9.properties))
        )

        # Apply cursor if provided
        if cursor:
            cursor_val, cursor_id = _decode_cursor(cursor)
            # If the sort column is a datetime, parse the cursor value back
            if hasattr(sort_col.type, 'python_type') and issubclass(sort_col.type.python_type, datetime):
                try:
                    cursor_val = datetime.fromisoformat(cursor_val)
                except (ValueError, TypeError):
                    raise HTTPException(status_code=400, detail="Invalid cursor value for datetime column")
            if sort_dir == "desc":
                stmt = stmt.where(
                    or_(
                        sort_col < cursor_val,
                        (sort_col == cursor_val) & (MoleculeORM_V9.id < cursor_id),
                    )
                )
            else:
                stmt = stmt.where(
                    or_(
                        sort_col > cursor_val,
                        (sort_col == cursor_val) & (MoleculeORM_V9.id > cursor_id),
                    )
                )
        else:
            stmt = stmt.offset(offset)

        stmt = stmt.order_by(order, MoleculeORM_V9.id).limit(limit + 1)

    result = await db.execute(stmt)
    mols = list(result.scalars().all())

    # Detect has_more and compute next_cursor
    has_more = len(mols) > limit
    if has_more:
        mols = mols[:limit]

    next_cursor = None
    if has_more and mols and not is_property_sort:
        last = mols[-1]
        last_val = getattr(last, sort_by, last.created_at)
        # Convert datetime to ISO string for JSON serialization
        if hasattr(last_val, 'isoformat'):
            last_val = last_val.isoformat()
        next_cursor = _encode_cursor(last_val, last.id)

    return MoleculeListResponse(
        molecules=[_mol_to_response(m) for m in mols],
        total=total,
        offset=offset,
        limit=limit,
        next_cursor=next_cursor,
        has_more=has_more,
    )


@router.get(
    "/phases/{phase_id}/molecules/stats",
    response_model=MoleculeStatsResponse,
)
async def molecule_stats(
    phase_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Lightweight stats for a phase: total, bookmarked, ai_generated."""
    await get_phase_owned(phase_id, user_id, db)

    stmt = select(
        func.count().label("total"),
        func.count().filter(MoleculeORM_V9.bookmarked == True).label("bookmarked"),  # noqa: E712
        func.count().filter(MoleculeORM_V9.ai_generated == True).label("ai_generated"),  # noqa: E712
    ).where(MoleculeORM_V9.phase_id == phase_id)

    result = await db.execute(stmt)
    row = result.one()

    return MoleculeStatsResponse(
        total=row.total,
        bookmarked=row.bookmarked,
        ai_generated=row.ai_generated,
    )


@router.get("/molecules/{molecule_id}", response_model=MoleculeResponse)
async def get_molecule(
    molecule_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Get molecule detail with all properties."""
    mol = await get_molecule_owned(molecule_id, user_id, db)
    return _mol_to_response(mol)


@router.put("/molecules/{molecule_id}/bookmark", response_model=MoleculeResponse)
async def bookmark_molecule(
    molecule_id: UUID,
    bookmarked: bool = Query(..., description="Set bookmark state"),
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Toggle bookmark on a single molecule."""
    mol = await get_molecule_owned(molecule_id, user_id, db)
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
    await get_phase_owned(phase_id, user_id, db)

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


# ---------------------------------------------------------------------------
# SMILES → 3D molblock
# ---------------------------------------------------------------------------

@router.post("/molecule/smiles-to-3d")
async def smiles_to_3d(
    body: dict = Body(...),
    user_id: str = Depends(require_v9_user),
):
    """Convert a SMILES string to a 3D molblock using RDKit."""
    smiles = body.get("smiles")
    if not smiles or not isinstance(smiles, str):
        raise HTTPException(status_code=422, detail="Missing or invalid 'smiles' field")

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        raise HTTPException(
            status_code=422,
            detail="RDKit is not available on this server",
        )

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise HTTPException(status_code=422, detail="Invalid SMILES string")

    try:
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        result = AllChem.EmbedMolecule(mol, params)
        if result != 0:
            raise ValueError("Embedding failed")
        AllChem.MMFFOptimizeMolecule(mol)
        mol = Chem.RemoveHs(mol)
        molblock = Chem.MolToMolBlock(mol)
    except Exception as exc:
        logger.warning("SMILES→3D failed for %s: %s", smiles, exc)
        raise HTTPException(
            status_code=422,
            detail=f"3D coordinate generation failed: {exc}",
        )

    return {"molblock": molblock}
