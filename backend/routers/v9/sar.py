"""
BindX V9 — SAR (Structure-Activity Relationship) Explorer router.

Endpoints:
- GET  /phases/{phase_id}/sar              — full SAR analysis
- POST /phases/{phase_id}/sar/apply-mmp    — preview MMP-guided suggestions
- POST /phases/{phase_id}/sar/import-suggestions — import suggestions into phase
"""

from __future__ import annotations

import logging
from datetime import datetime, timezone
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.orm import selectinload

from auth_v9 import require_v9_user
from database_v9 import get_v9_db
from models_v9 import MoleculeORM_V9, MoleculePropertyORM_V9, RunORM_V9
from pipeline.sar_analysis import (
    analyze_sar,
    apply_mmp_transform,
    get_cached_sar,
    set_cached_sar,
)
from routers.v9.deps import get_phase_owned
from routers.v9.reports import flatten_molecule
from schemas_v9 import MMPApplyRequest, MMPApplyResponse, MMPImportRequest, MMPSuggestion

logger = logging.getLogger(__name__)

router = APIRouter()


# ---------------------------------------------------------------------------
# GET — Full SAR analysis
# ---------------------------------------------------------------------------

@router.get("/phases/{phase_id}/sar")
async def get_sar_analysis(
    phase_id: UUID,
    property_key: str = Query("auto", description="Property to analyze (auto = best available)"),
    scaffold: str | None = Query(None, description="Scaffold SMILES to scope to a chemical series"),
    matrix_r1: str | None = Query(None, description="R-group key for matrix rows (e.g. R1)"),
    matrix_r2: str | None = Query(None, description="R-group key for matrix columns (e.g. R2)"),
    r_filters: str | None = Query(None, description="JSON dict of R-group filters (e.g. {\"R3\": \"F\"})"),
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Run SAR analysis on all molecules in a phase, optionally scoped to a scaffold series."""

    phase = await get_phase_owned(phase_id, user_id, db)

    # Load molecules with properties
    stmt = (
        select(MoleculeORM_V9)
        .where(MoleculeORM_V9.phase_id == phase.id)
        .options(selectinload(MoleculeORM_V9.properties))
    )
    result = await db.execute(stmt)
    molecules_orm = result.scalars().all()

    if not molecules_orm:
        raise HTTPException(status_code=404, detail="No molecules in this phase")

    n_mols = len(molecules_orm)
    scaffold_key = scaffold or "__auto__"

    # Check cache (include matrix axes in cache key)
    cache_suffix = f"{scaffold_key}|{matrix_r1}|{matrix_r2}|{r_filters or ''}"
    cache_key = str(phase.id)
    cached = get_cached_sar(cache_key, property_key, n_mols, cache_suffix)
    if cached is not None:
        logger.info("SAR cache hit for phase %s (%d mols, prop=%s)", phase_id, n_mols, property_key)
        return cached

    # Parse r_filters JSON
    parsed_filters = None
    if r_filters:
        try:
            import json
            parsed_filters = json.loads(r_filters)
        except (json.JSONDecodeError, TypeError):
            logger.warning("Invalid r_filters JSON: %s", r_filters)

    # Flatten molecules
    molecules = [flatten_molecule(m) for m in molecules_orm]

    # Run analysis
    try:
        sar_result = analyze_sar(
            molecules,
            property_key=property_key,
            scaffold_filter=scaffold,
            matrix_r1=matrix_r1,
            matrix_r2=matrix_r2,
            r_filters=parsed_filters,
        )
    except Exception as exc:
        logger.error("SAR analysis failed: %s", exc, exc_info=True)
        raise HTTPException(status_code=500, detail=f"SAR analysis failed: {exc}")

    # Cache result
    set_cached_sar(cache_key, property_key, n_mols, sar_result, cache_suffix)

    return sar_result


# ---------------------------------------------------------------------------
# POST — Apply MMP transformation (preview, nothing saved)
# ---------------------------------------------------------------------------

@router.post("/phases/{phase_id}/sar/apply-mmp", response_model=MMPApplyResponse)
async def apply_mmp(
    phase_id: UUID,
    body: MMPApplyRequest,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Apply an MMP transformation to generate molecule suggestions (preview only)."""

    phase = await get_phase_owned(phase_id, user_id, db)

    # Load molecules
    stmt = (
        select(MoleculeORM_V9)
        .where(MoleculeORM_V9.phase_id == phase.id)
        .options(selectinload(MoleculeORM_V9.properties))
    )
    if body.molecule_ids:
        stmt = stmt.where(MoleculeORM_V9.id.in_(body.molecule_ids))

    result = await db.execute(stmt)
    molecules_orm = result.scalars().all()

    if not molecules_orm:
        raise HTTPException(status_code=404, detail="No molecules found")

    molecules = [flatten_molecule(m) for m in molecules_orm]

    # Apply transform
    suggestions = apply_mmp_transform(
        molecules,
        from_smiles=body.from_smiles,
        to_smiles=body.to_smiles,
        max_results=min(body.max_results, 100),
    )

    transform_label = f"{body.from_smiles} → {body.to_smiles}"

    return MMPApplyResponse(
        suggestions=[MMPSuggestion(**s) for s in suggestions],
        transform=transform_label,
        n_input_molecules=len(molecules),
    )


# ---------------------------------------------------------------------------
# POST — Import MMP suggestions into phase
# ---------------------------------------------------------------------------

@router.post("/phases/{phase_id}/sar/import-suggestions")
async def import_mmp_suggestions(
    phase_id: UUID,
    body: MMPImportRequest,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Import MMP-suggested molecules into the phase as a new import run."""

    phase = await get_phase_owned(phase_id, user_id, db)

    if phase.status == "frozen":
        raise HTTPException(status_code=400, detail="Cannot import into a frozen phase")

    if not body.smiles_list:
        raise HTTPException(status_code=422, detail="No SMILES provided")

    if len(body.smiles_list) > 200:
        raise HTTPException(status_code=422, detail="Maximum 200 molecules per import")

    try:
        from rdkit import Chem
    except ImportError:
        raise HTTPException(status_code=500, detail="RDKit not available")

    now = datetime.now(timezone.utc)

    # Create import run
    run = RunORM_V9(
        phase_id=phase.id,
        type="import",
        status="completed",
        config={
            "source": "mmp_suggestion",
            "transform": body.transform or "",
            "n_imported": 0,  # updated below
        },
        progress=100,
        current_step="Import complete",
        started_at=now,
        completed_at=now,
    )
    db.add(run)
    await db.flush()  # get run.id

    # Create molecules
    imported = 0
    for i, smiles in enumerate(body.smiles_list):
        rdmol = Chem.MolFromSmiles(smiles)
        if rdmol is None:
            continue

        canonical = Chem.MolToSmiles(rdmol)
        name = body.names[i] if body.names and i < len(body.names) else f"MMP-{imported + 1}"

        mol_orm = MoleculeORM_V9(
            phase_id=phase.id,
            smiles=smiles,
            canonical_smiles=canonical,
            name=name,
            source_run_id=run.id,
            ai_generated=True,
            tags=["mmp_suggestion"],
        )
        db.add(mol_orm)
        imported += 1

    # Update run config with actual count
    run.config = {**run.config, "n_imported": imported}

    await db.commit()

    return {
        "imported": imported,
        "run_id": str(run.id),
        "skipped": len(body.smiles_list) - imported,
    }
