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

from fastapi import APIRouter, Body, Depends, HTTPException, Query, Response
from sqlalchemy import Float, case, func, literal, or_, select, text, tuple_
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
    MoleculeAnnotationUpdate,
    MoleculeListResponse,
    MoleculeResponse,
    MoleculeStatsResponse,
)

logger = logging.getLogger(__name__)

router = APIRouter()

# Direct columns on the molecules table that can be sorted
VALID_SORT_FIELDS = {
    "created_at", "smiles", "name", "bookmarked", "canonical_smiles",
    "invalidated", "source_run_id", "generation_level",
}

# Legacy property_name aliases for backward compat with old "admet" runs.
# When sorting/filtering, we search both the new name and the legacy name.
_PROP_NAME_ALIASES = {
    "adme": "admet",      # new ADME runs store as "adme", old ones as "admet"
    "toxicity": "admet",  # new Toxicity runs store as "toxicity", old ones as "admet"
}

# Maps frontend sort key → (property_name, jsonb_path_keys)
# Used to resolve nested JSON properties for server-side sorting.
# json_path=[] means flat EAV (property_value is a bare scalar).
SORT_FIELD_MAP = {
    # --- physicochemical ---
    "MW": ("physicochemical", ["MW"]),
    "HBD": ("physicochemical", ["hbd"]),
    "HBA": ("physicochemical", ["hba"]),
    "QED": ("physicochemical", ["qed"]),
    "logP": ("physicochemical", ["logP"]),
    "TPSA": ("physicochemical", ["tpsa"]),
    "rotatable_bonds": ("physicochemical", ["rotatable_bonds"]),
    "lipinski_pass": ("physicochemical", ["lipinski_pass"]),
    "heavy_atom_count": ("physicochemical", ["heavy_atom_count"]),
    "ro3_pass": ("physicochemical", ["ro3_pass"]),
    "inchikey": ("physicochemical", ["inchikey"]),
    # --- sa_score ---
    "sa_score": ("sa_score", ["sa_score"]),
    # --- ligand_efficiency ---
    "ligand_efficiency": ("ligand_efficiency", ["ligand_efficiency"]),
    # --- docking ---
    "docking_score": ("docking", ["affinity"]),
    "cnn_score": ("docking", ["cnn_score"]),
    "cnn_affinity": ("docking", ["cnn_affinity"]),
    "cnn_vs": ("docking", ["cnn_vs"]),
    "consensus_ecr": ("docking", ["consensus_ecr"]),
    "pocket_distance": ("docking", ["pocket_distance"]),
    # --- adme — absorption ---
    "solubility": ("adme", ["absorption", "solubility"]),
    "oral_bioavailability": ("adme", ["absorption", "oral_bioavailability"]),
    "pgp_substrate": ("adme", ["absorption", "pgp_substrate"]),
    "intestinal_permeability": ("adme", ["absorption", "intestinal_permeability"]),
    # --- adme — distribution ---
    "BBB": ("adme", ["distribution", "bbb_permeability"]),
    "vd": ("adme", ["distribution", "vd"]),
    "plasma_protein_binding": ("adme", ["distribution", "plasma_protein_binding"]),
    # --- adme — metabolism ---
    "cyp1a2_inhibitor": ("adme", ["metabolism", "cyp1a2_inhibitor"]),
    "cyp2c9_inhibitor": ("adme", ["metabolism", "cyp2c9_inhibitor"]),
    "cyp2d6_inhibitor": ("adme", ["metabolism", "cyp2d6_inhibitor"]),
    "cyp3a4_inhibitor": ("adme", ["metabolism", "cyp3a4_inhibitor"]),
    "cyp2c19_inhibitor": ("adme", ["metabolism", "cyp2c19_inhibitor"]),
    # --- adme — excretion ---
    "half_life": ("adme", ["excretion", "half_life"]),
    # --- toxicity ---
    "composite_score": ("toxicity", ["composite_score"]),
    "safety_color_code": ("toxicity", ["color_code"]),
    "hERG": ("toxicity", ["toxicity", "herg_inhibition"]),
    "hepatotoxicity": ("toxicity", ["toxicity", "hepatotoxicity"]),
    "carcinogenicity": ("toxicity", ["toxicity", "carcinogenicity"]),
    "ames_mutagenicity": ("toxicity", ["toxicity", "ames_mutagenicity"]),
    "skin_sensitization": ("toxicity", ["toxicity", "skin_sensitization"]),
    # --- druglikeness_rules (piggyback on ADME run) ---
    "cns_mpo": ("druglikeness_rules", ["cns_mpo"]),
    "pfizer_alert": ("druglikeness_rules", ["pfizer_alert"]),
    "gsk_alert": ("druglikeness_rules", ["gsk_alert"]),
    "brenk_alert": ("druglikeness_rules", ["brenk_alert"]),
    # --- confidence ---
    "confidence_score": ("confidence", ["confidence_score"]),
    "pains_alert": ("confidence", ["pains_alert"]),
    "applicability_domain": ("confidence", ["applicability_domain"]),
    # --- clustering ---
    "cluster_id": ("clustering", ["cluster_id"]),
    "tanimoto_to_centroid": ("clustering", ["tanimoto_to_centroid"]),
    "is_representative": ("clustering", ["is_representative"]),
    # --- off_target ---
    "selectivity_score": ("off_target", ["selectivity_score"]),
    "off_target_hits": ("off_target", ["off_target_hits"]),
    "selectivity_ratio": ("off_target", ["selectivity_ratio"]),
    # --- interactions ---
    "n_interactions": ("interactions", ["n_interactions"]),
    "functional_contacts": ("interactions", ["functional_contacts"]),
    "interaction_quality": ("interactions", ["interaction_quality"]),
    "key_hbonds": ("interactions", ["key_hbonds"]),
    # --- scaffold ---
    "scaffold_smiles": ("scaffold", ["scaffold_smiles"]),
    "n_modifiable_positions": ("scaffold", ["n_modifiable_positions"]),
    "brics_bond_count": ("scaffold", ["brics_bond_count"]),
    "scaffold_n_rings": ("scaffold", ["scaffold_n_rings"]),
    "scaffold_mw": ("scaffold", ["scaffold_mw"]),
    "scaffold_group": ("scaffold", ["scaffold_group"]),
    "scaffold_group_size": ("scaffold", ["scaffold_group_size"]),
    # --- retrosynthesis ---
    "n_synth_steps": ("retrosynthesis", ["n_synth_steps"]),
    "synth_confidence": ("retrosynthesis", ["synth_confidence"]),
    "synth_cost_estimate": ("retrosynthesis", ["synth_cost_estimate"]),
    "reagents_available": ("retrosynthesis", ["reagents_available"]),
    # --- pharmacophore ---
    "pharmacophore_features": ("pharmacophore", ["pharmacophore_features"]),
    "pharmacophore_similarity": ("pharmacophore", ["pharmacophore_similarity"]),
    # --- Import metadata (flat EAV — property_value is a bare scalar) ---
    "pchembl_value": ("pchembl_value", []),
    "activity_value_nM": ("activity_value_nM", []),
    "activity_type": ("activity_type", []),
    "import_mwt": ("mwt", []),
    "import_logp": ("logp", []),
    "import_hbd": ("import_hbd", []),
    "import_hba": ("import_hba", []),
    "import_tpsa": ("import_tpsa", []),
    "source": ("source", []),
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
        user_comment=mol.user_comment,
        ai_comment=mol.ai_comment,
        tags=mol.tags or [],
        invalidated=mol.invalidated,
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
        escaped = search.replace("%", "\\%").replace("_", "\\_")
        pattern = f"%{escaped}%"
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
    is_property_sort = sort_by not in VALID_SORT_FIELDS

    if is_property_sort:
        # Resolve property_name and JSON path from SORT_FIELD_MAP
        mapping = SORT_FIELD_MAP.get(sort_by)
        if mapping:
            prop_name, json_path = mapping
        else:
            # Fallback: treat sort_by as both property_name and single-key path
            prop_name = sort_by
            json_path = ["value"]

        # Build sort expression from JSONB path
        json_col = MoleculePropertyORM_V9.property_value
        if not json_path:
            # Flat EAV: property_value is a bare scalar (e.g. "6.5"), cast directly
            json_col = json_col.astext
        else:
            # Nested JSONB: ->'key1'->...->>last_key  (astext on final key for castability)
            for key in json_path[:-1]:
                json_col = json_col[key]          # -> returns JSONB
            json_col = json_col[json_path[-1]].astext   # ->> returns text

        # Support legacy property names (e.g. "admet" → "adme"/"toxicity")
        prop_names = [prop_name]
        legacy = _PROP_NAME_ALIASES.get(prop_name)
        if legacy:
            prop_names.append(legacy)

        prop_sub = (
            select(
                json_col.cast(Float).label("prop_val")
            )
            .where(
                MoleculePropertyORM_V9.molecule_id == MoleculeORM_V9.id,
                MoleculePropertyORM_V9.property_name.in_(prop_names),
            )
            .order_by(MoleculePropertyORM_V9.created_at.desc())
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

    responses = [_mol_to_response(m) for m in mols]

    return MoleculeListResponse(
        molecules=responses,
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
    """Stats for a phase: counts + global min/max per numeric column."""
    await get_phase_owned(phase_id, user_id, db)

    stmt = select(
        func.count().label("total"),
        func.count().filter(MoleculeORM_V9.bookmarked == True).label("bookmarked"),  # noqa: E712
        func.count().filter(MoleculeORM_V9.ai_generated == True).label("ai_generated"),  # noqa: E712
    ).where(MoleculeORM_V9.phase_id == phase_id)

    result = await db.execute(stmt)
    row = result.one()

    # Compute global min/max per numeric column from properties
    column_ranges = await _compute_column_ranges(phase_id, db)

    return MoleculeStatsResponse(
        total=row.total,
        bookmarked=row.bookmarked,
        ai_generated=row.ai_generated,
        column_ranges=column_ranges,
    )


async def _compute_column_ranges(phase_id: UUID, db: AsyncSession) -> dict:
    """Compute global min/max for each numeric property column using SQL aggregation.

    Uses JSONB path extraction + SQL MIN/MAX to avoid loading all properties into Python.
    Groups queries by property_name to minimize round-trips.
    """
    # Group SORT_FIELD_MAP entries by property_name for batch queries
    by_prop: dict[str, list[tuple[str, list[str]]]] = {}
    for col_key, (prop_name, json_path) in SORT_FIELD_MAP.items():
        by_prop.setdefault(prop_name, []).append((col_key, json_path))

    mol_subq = select(MoleculeORM_V9.id).where(MoleculeORM_V9.phase_id == phase_id)
    ranges: dict = {}

    for prop_name, columns in by_prop.items():
        # Build list of property names to match (including legacy aliases)
        prop_names = [prop_name]
        legacy = _PROP_NAME_ALIASES.get(prop_name)
        if legacy:
            prop_names.append(legacy)

        for col_key, json_path in columns:
            # Build JSONB path extraction expression
            if not json_path:
                # Flat EAV — property_value is a scalar, cast to float
                val_expr = MoleculePropertyORM_V9.property_value.cast(Float)
            elif len(json_path) == 1:
                val_expr = (
                    MoleculePropertyORM_V9.property_value[json_path[0]].as_float()
                )
            elif len(json_path) == 2:
                val_expr = (
                    MoleculePropertyORM_V9.property_value[json_path[0]][json_path[1]].as_float()
                )
            elif len(json_path) == 3:
                val_expr = (
                    MoleculePropertyORM_V9.property_value[json_path[0]][json_path[1]][json_path[2]].as_float()
                )
            else:
                continue  # Skip deeper paths

            stmt = (
                select(
                    func.min(val_expr).label("min_val"),
                    func.max(val_expr).label("max_val"),
                )
                .where(
                    MoleculePropertyORM_V9.molecule_id.in_(mol_subq),
                    MoleculePropertyORM_V9.property_name.in_(prop_names),
                )
            )

            try:
                result = await db.execute(stmt)
                row = result.one()
                if row.min_val is not None and row.max_val is not None:
                    ranges[col_key] = {"min": float(row.min_val), "max": float(row.max_val)}
            except Exception:
                # Skip columns where JSONB extraction fails (non-numeric values)
                pass

    return ranges


@router.get("/molecules/{molecule_id}", response_model=MoleculeResponse)
async def get_molecule(
    molecule_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Get molecule detail with all properties."""
    mol = await get_molecule_owned(molecule_id, user_id, db)
    return _mol_to_response(mol)


@router.get("/molecules/{molecule_id}/feature-map")
async def get_feature_map(
    molecule_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Generate a 2D SVG with pharmacophore features highlighted."""
    mol = await get_molecule_owned(molecule_id, user_id, db)
    smiles = mol.canonical_smiles or mol.smiles
    if not smiles:
        raise HTTPException(status_code=422, detail="Molecule has no SMILES")

    try:
        from pipeline.pharmacophore import generate_feature_map_svg
        svg = generate_feature_map_svg(smiles)
    except Exception as e:
        logger.error("Feature map generation failed: %s", e)
        raise HTTPException(status_code=500, detail=f"Feature map generation failed: {e}")

    return Response(content=svg, media_type="image/svg+xml")


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


@router.patch("/molecules/{molecule_id}/annotations", response_model=MoleculeResponse)
async def update_annotations(
    molecule_id: UUID,
    body: MoleculeAnnotationUpdate,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Update user annotations (comment, tags, invalidated) on a molecule."""
    mol = await get_molecule_owned(molecule_id, user_id, db)
    if body.user_comment is not None:
        mol.user_comment = body.user_comment
    if body.tags is not None:
        mol.tags = body.tags
    if body.invalidated is not None:
        mol.invalidated = body.invalidated
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
# SMILES → feature map SVG
# ---------------------------------------------------------------------------

@router.post("/feature-map")
async def feature_map_from_smiles(
    body: dict = Body(...),
    user_id: str = Depends(require_v9_user),
):
    """Generate a 2D SVG with pharmacophore features from a SMILES string."""
    smiles = body.get("smiles")
    if not smiles or not isinstance(smiles, str):
        raise HTTPException(status_code=422, detail="Missing or invalid 'smiles' field")

    try:
        from pipeline.pharmacophore import generate_feature_map_svg
        svg = generate_feature_map_svg(smiles, width=300, height=200)
    except Exception as e:
        logger.error("Feature map generation failed for SMILES: %s", e)
        raise HTTPException(status_code=500, detail=f"Feature map generation failed: {e}")

    return Response(content=svg, media_type="image/svg+xml")


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


# ---------------------------------------------------------------------------
# Property detail endpoints (safety, synthesis, confidence)
# ---------------------------------------------------------------------------

async def _get_molecule_property(
    molecule_id: UUID, property_name: str, user_id: str, db: AsyncSession
) -> dict:
    """Fetch a specific property for a molecule, return its value or 404."""
    mol = await get_molecule_owned(molecule_id, user_id, db)
    stmt = (
        select(MoleculePropertyORM_V9)
        .where(
            MoleculePropertyORM_V9.molecule_id == mol.id,
            MoleculePropertyORM_V9.property_name == property_name,
        )
        .order_by(MoleculePropertyORM_V9.created_at.desc())
        .limit(1)
    )
    result = await db.execute(stmt)
    prop = result.scalar_one_or_none()
    if not prop:
        raise HTTPException(
            status_code=404,
            detail=f"No {property_name} data for this molecule",
        )
    return prop.property_value


@router.get("/molecules/{molecule_id}/safety")
async def get_molecule_safety(
    molecule_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Get safety profile for a molecule."""
    return await _get_molecule_property(molecule_id, "safety", user_id, db)


@router.get("/molecules/{molecule_id}/synthesis")
async def get_molecule_synthesis(
    molecule_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Get retrosynthesis route for a molecule."""
    return await _get_molecule_property(molecule_id, "retrosynthesis", user_id, db)


@router.get("/molecules/{molecule_id}/confidence")
async def get_molecule_confidence(
    molecule_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Get confidence breakdown for a molecule."""
    return await _get_molecule_property(molecule_id, "confidence", user_id, db)


# ---------------------------------------------------------------------------
# Scaffold analysis
# ---------------------------------------------------------------------------

@router.post("/molecule/smiles-to-svg")
async def smiles_to_svg(
    body: dict = Body(...),
    user_id: str = Depends(require_v9_user),
):
    """Render a SMILES string as a 2D SVG image using RDKit."""
    smiles = body.get("smiles")
    if not smiles or not isinstance(smiles, str):
        raise HTTPException(status_code=422, detail="Missing or invalid 'smiles' field")

    width = int(body.get("width", 120))
    height = int(body.get("height", 80))

    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError:
        raise HTTPException(status_code=500, detail="RDKit not available")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise HTTPException(status_code=422, detail="Invalid SMILES")

    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.drawOptions().addStereoAnnotation = False
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return Response(content=drawer.GetDrawingText(), media_type="image/svg+xml")


@router.post("/molecule/smiles-to-svg-diff")
async def smiles_to_svg_diff(
    body: dict = Body(...),
    user_id: str = Depends(require_v9_user),
):
    """Render two SMILES as SVGs with differing atoms highlighted in yellow."""
    smiles_a = body.get("smiles_a")
    smiles_b = body.get("smiles_b")
    width = int(body.get("width", 150))
    height = int(body.get("height", 100))

    if not smiles_a or not smiles_b:
        raise HTTPException(status_code=422, detail="Missing smiles_a or smiles_b")

    try:
        from rdkit import Chem
        from rdkit.Chem import rdFMCS
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError:
        raise HTTPException(status_code=500, detail="RDKit not available")

    mol_a = Chem.MolFromSmiles(smiles_a)
    mol_b = Chem.MolFromSmiles(smiles_b)
    if not mol_a or not mol_b:
        raise HTTPException(status_code=422, detail="Invalid SMILES")

    def _mol_to_svg_highlighted(mol, highlight_atoms, highlight_bonds, w, h):
        drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
        drawer.drawOptions().addStereoAnnotation = False
        atom_colors = {i: (1.0, 0.84, 0.0) for i in highlight_atoms}
        bond_colors = {i: (1.0, 0.84, 0.0) for i in highlight_bonds}
        atom_radii = {i: 0.4 for i in highlight_atoms}
        drawer.DrawMolecule(
            mol,
            highlightAtoms=list(highlight_atoms),
            highlightAtomColors=atom_colors,
            highlightBonds=list(highlight_bonds),
            highlightBondColors=bond_colors,
            highlightAtomRadii=atom_radii,
        )
        drawer.FinishDrawing()
        return drawer.GetDrawingText()

    def _plain_svg(mol, w, h):
        drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
        drawer.drawOptions().addStereoAnnotation = False
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()

    def _get_diff(mol, mcs_mol):
        if not mcs_mol:
            return set(range(mol.GetNumAtoms())), set()
        match = mol.GetSubstructMatch(mcs_mol)
        if not match:
            # Try all matches and pick the one with most atoms
            matches = mol.GetSubstructMatches(mcs_mol, uniquify=True)
            if not matches:
                return set(range(mol.GetNumAtoms())), set()
            match = max(matches, key=len)
        matched_atoms = set(match)
        # Sanity: if MCS matched < 40% of atoms, diff is not meaningful
        if len(matched_atoms) < 0.4 * mol.GetNumAtoms():
            return set(), set()
        diff_atoms = set(range(mol.GetNumAtoms())) - matched_atoms
        diff_bonds = {
            b.GetIdx() for b in mol.GetBonds()
            if b.GetBeginAtomIdx() in diff_atoms
            or b.GetEndAtomIdx() in diff_atoms
        }
        return diff_atoms, diff_bonds

    try:
        # Try progressively relaxed MCS parameters
        mcs_mol = None
        for params in [
            # 1. Relaxed: any bond, ring matches ring, generous timeout
            dict(atomCompare=rdFMCS.AtomCompare.CompareElements,
                 bondCompare=rdFMCS.BondCompare.CompareAny,
                 ringMatchesRingOnly=True, timeout=5),
            # 2. Even more relaxed: no ring constraint
            dict(atomCompare=rdFMCS.AtomCompare.CompareElements,
                 bondCompare=rdFMCS.BondCompare.CompareAny,
                 timeout=3),
        ]:
            mcs = rdFMCS.FindMCS([mol_a, mol_b], **params)
            if mcs.smartsString:
                candidate = Chem.MolFromSmarts(mcs.smartsString)
                if candidate and candidate.GetNumAtoms() >= 0.4 * min(
                    mol_a.GetNumAtoms(), mol_b.GetNumAtoms()
                ):
                    mcs_mol = candidate
                    break

        diff_a_atoms, diff_a_bonds = _get_diff(mol_a, mcs_mol)
        diff_b_atoms, diff_b_bonds = _get_diff(mol_b, mcs_mol)

        # If both diffs are empty (MCS too small → fallback triggered), render plain
        if not diff_a_atoms and not diff_b_atoms:
            svg_a = _plain_svg(mol_a, width, height)
            svg_b = _plain_svg(mol_b, width, height)
        else:
            svg_a = _mol_to_svg_highlighted(mol_a, diff_a_atoms, diff_a_bonds, width, height)
            svg_b = _mol_to_svg_highlighted(mol_b, diff_b_atoms, diff_b_bonds, width, height)
    except Exception:
        logger.warning("MCS diff failed for %s vs %s, falling back to plain SVG", smiles_a, smiles_b)
        svg_a = _plain_svg(mol_a, width, height)
        svg_b = _plain_svg(mol_b, width, height)

    return {"svg_a": svg_a, "svg_b": svg_b}


@router.post("/scaffold-analysis")
async def scaffold_analysis(
    body: dict = Body(...),
    user_id: str = Depends(require_v9_user),
):
    """Analyze a molecule's scaffold and R-group positions."""
    smiles = body.get("smiles")
    if not smiles or not isinstance(smiles, str):
        raise HTTPException(status_code=422, detail="Missing or invalid 'smiles' field")

    try:
        from pipeline.scaffold_analysis import analyze_scaffold
        result = analyze_scaffold(smiles)
        return result
    except Exception as e:
        logger.error("Scaffold analysis failed: %s", e)
        raise HTTPException(status_code=500, detail=f"Scaffold analysis failed: {e}")


# ---------------------------------------------------------------------------
# ESMFold structure prediction
# ---------------------------------------------------------------------------

@router.post("/predict-structure")
async def predict_structure(
    body: dict = Body(...),
    user_id: str = Depends(require_v9_user),
):
    """Predict protein structure from FASTA sequence using ESMFold API."""
    import httpx

    sequence = body.get("sequence", "").strip()
    if not sequence:
        raise HTTPException(status_code=422, detail="Missing 'sequence' field")

    # Clean FASTA: remove header lines and whitespace
    lines = sequence.split("\n")
    clean_seq = "".join(l.strip() for l in lines if not l.startswith(">"))

    if len(clean_seq) < 10:
        raise HTTPException(status_code=422, detail="Sequence too short (min 10 residues)")
    if len(clean_seq) > 2000:
        raise HTTPException(status_code=422, detail="Sequence too long (max 2000 residues)")

    # Validate amino acid characters
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    invalid = set(clean_seq.upper()) - valid_aa
    if invalid:
        raise HTTPException(status_code=422, detail=f"Invalid amino acid characters: {invalid}")

    try:
        async with httpx.AsyncClient(timeout=300.0) as client:
            resp = await client.post(
                "https://api.esmatlas.com/foldSequence/v1/pdb/",
                content=clean_seq,
                headers={"Content-Type": "text/plain"},
            )
            if resp.status_code != 200:
                raise HTTPException(
                    status_code=502,
                    detail=f"ESMFold API returned {resp.status_code}: {resp.text[:200]}"
                )
            pdb_text = resp.text
            return {
                "pdb_text": pdb_text,
                "sequence_length": len(clean_seq),
                "source": "esmfold",
                "method": "ESMFold v1",
            }
    except httpx.TimeoutException:
        raise HTTPException(status_code=504, detail="ESMFold prediction timed out (>5 min)")
    except HTTPException:
        raise
    except Exception as e:
        logger.error("ESMFold prediction failed: %s", e)
        raise HTTPException(status_code=500, detail=f"Structure prediction failed: {e}")
