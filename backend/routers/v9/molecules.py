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
    # --- enrichment ---
    "interactions_count": ("enrichment", ["interactions_count"]),
    "scaffold": ("enrichment", ["scaffold"]),
    # --- retrosynthesis ---
    "n_synth_steps": ("retrosynthesis", ["n_synth_steps"]),
    "synth_confidence": ("retrosynthesis", ["synth_confidence"]),
    "synth_cost_estimate": ("retrosynthesis", ["synth_cost_estimate"]),
    "reagents_available": ("retrosynthesis", ["reagents_available"]),
    # --- activity_cliffs ---
    "is_cliff": ("activity_cliffs", ["is_cliff"]),
    "sali_max": ("activity_cliffs", ["sali_max"]),
    "n_cliffs": ("activity_cliffs", ["n_cliffs"]),
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
    """Compute global min/max for each numeric property column in a phase."""
    # Load all property rows for this phase
    prop_stmt = (
        select(MoleculePropertyORM_V9)
        .where(
            MoleculePropertyORM_V9.molecule_id.in_(
                select(MoleculeORM_V9.id).where(MoleculeORM_V9.phase_id == phase_id)
            )
        )
    )
    prop_result = await db.execute(prop_stmt)
    all_props = prop_result.scalars().all()

    ranges: dict = {}
    for p in all_props:
        for col_key, (prop_name, json_path) in SORT_FIELD_MAP.items():
            # Match current name or legacy alias
            legacy = _PROP_NAME_ALIASES.get(prop_name)
            if p.property_name != prop_name and p.property_name != legacy:
                continue
            val = p.property_value
            for k in json_path:
                if isinstance(val, dict):
                    val = val.get(k)
                else:
                    val = None
                    break
            if val is not None and isinstance(val, (int, float)):
                if col_key not in ranges:
                    ranges[col_key] = {"min": val, "max": val}
                else:
                    ranges[col_key]["min"] = min(ranges[col_key]["min"], val)
                    ranges[col_key]["max"] = max(ranges[col_key]["max"], val)
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
