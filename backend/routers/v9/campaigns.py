"""
BindX V9 — Campaign CRUD router.

Endpoints: create, list, get, update, set reference ligand, suggest ligands.
"""

from __future__ import annotations

import asyncio
import logging
from typing import Optional
from uuid import UUID

import httpx
from fastapi import APIRouter, Depends, HTTPException, Query
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


@router.post("/campaigns/{campaign_id}/reference-ligands")
async def add_reference_ligand(
    campaign_id: UUID,
    body: dict,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Add a reference ligand to the campaign (max 5)."""
    campaign = await get_campaign_owned(campaign_id, user_id, db)

    current = list(campaign.reference_ligands or [])
    if len(current) >= 5:
        raise HTTPException(400, "Maximum 5 reference ligands allowed.")

    smiles = body.get("smiles", "").strip()
    name = body.get("name", "").strip() or "Reference"
    source = body.get("source", "").strip() or None

    if not smiles:
        raise HTTPException(422, "Missing 'smiles' field.")

    # Validate SMILES with RDKit
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES")
    except Exception:
        raise HTTPException(422, f"Invalid SMILES: {smiles}")

    # Generate feature map SVG and count features
    try:
        from pipeline.pharmacophore import generate_feature_map_svg, _count_features
        feature_map_svg = generate_feature_map_svg(smiles, 300, 200)
        ref_feature_counts = _count_features(mol)
    except Exception as e:
        logger.error("Feature map generation failed: %s", e)
        raise HTTPException(500, f"Failed to generate feature map: {e}")

    ref_data = {
        "index": len(current),
        "smiles": smiles,
        "name": name,
        "feature_map_svg": feature_map_svg,
        "ref_feature_counts": ref_feature_counts,
    }
    if source:
        ref_data["source"] = source

    current.append(ref_data)
    campaign.reference_ligands = current
    await db.flush()

    return ref_data


@router.delete("/campaigns/{campaign_id}/reference-ligands/{index}")
async def remove_reference_ligand(
    campaign_id: UUID,
    index: int,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Remove a reference ligand by index and re-index the rest."""
    campaign = await get_campaign_owned(campaign_id, user_id, db)

    current = list(campaign.reference_ligands or [])
    if index < 0 or index >= len(current):
        raise HTTPException(404, f"Reference ligand index {index} not found.")

    current.pop(index)
    # Re-index
    for i, ref in enumerate(current):
        ref["index"] = i

    campaign.reference_ligands = current
    await db.flush()

    return {"reference_ligands": current}


# ---------------------------------------------------------------------------
# Suggest reference ligands (UniProt + RCSB + ChEMBL)
# ---------------------------------------------------------------------------

# Crystallization buffer/ion artifacts to exclude from PDB ligands
_PDB_EXCLUDE = frozenset({
    "SO4", "PO4", "GOL", "EDO", "ACT", "MG", "ZN", "CA", "NA", "CL",
    "DMS", "PEG", "MPD", "BME", "FMT", "NO3", "MES", "TRS", "EPE",
    "IMD", "ACE", "ACY", "NH4", "IOD", "BR", "12P", "BU1", "ACM",
    "P6G", "1PE", "PGE", "NAG", "FUC", "MAN", "GAL", "BGC", "GLC",
    "SUC", "BOG", "CLR", "PLM", "OLC", "OLA", "MYR", "STE", "CIT",
    "HEX", "UNX", "UNL", "HOH", "DOD", "MLI", "TAR", "FLC",
})


async def _fetch_uniprot_ligands(client: httpx.AsyncClient, uniprot_id: str) -> list[dict]:
    """Fetch endogenous/natural ligands from UniProt binding site annotations."""
    suggestions = []
    try:
        resp = await client.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json", timeout=15)
        if resp.status_code != 200:
            return []
        data = resp.json()
        seen = set()
        for feat in data.get("features", []):
            if feat.get("type") != "Binding site":
                continue
            ligand = feat.get("ligand", {})
            name = ligand.get("name", "")
            chebi_id = ligand.get("id", "")  # e.g. "ChEBI:CHEBI:18243"
            if not name or name in seen:
                continue
            seen.add(name)
            # Try to get SMILES from ChEBI via OLS4
            smiles = None
            if chebi_id and chebi_id.startswith("ChEBI:"):
                obo_id = chebi_id.replace("ChEBI:", "")
                try:
                    ols_resp = await client.get(
                        f"https://www.ebi.ac.uk/ols4/api/ontologies/chebi/terms",
                        params={"obo_id": obo_id},
                        timeout=10,
                    )
                    if ols_resp.status_code == 200:
                        terms = ols_resp.json().get("_embedded", {}).get("terms", [])
                        if terms:
                            smiles = (terms[0].get("annotation", {}).get("smiles_string") or [None])[0]
                except Exception:
                    pass
            if smiles:
                suggestions.append({
                    "smiles": smiles,
                    "name": name,
                    "source": "uniprot",
                    "source_label": "Endogenous ligand",
                })
    except Exception as e:
        logger.warning("UniProt ligand fetch failed for %s: %s", uniprot_id, e)
    return suggestions


async def _fetch_cocrystal_ligands(client: httpx.AsyncClient, pdb_id: str, known_ligand: str | None) -> list[dict]:
    """Fetch co-crystallized ligands from a specific PDB entry."""
    suggestions = []
    try:
        # Get all non-polymer entities
        resp = await client.get(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}", timeout=10)
        if resp.status_code != 200:
            return []
        entry = resp.json()
        np_ids = entry.get("rcsb_entry_container_identifiers", {}).get("non_polymer_entity_ids", [])

        for np_id in np_ids:
            try:
                np_resp = await client.get(
                    f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id}/{np_id}",
                    timeout=10,
                )
                if np_resp.status_code != 200:
                    continue
                np_data = np_resp.json()
                comp_id = np_data.get("pdbx_entity_nonpoly", {}).get("comp_id", "")
                if not comp_id or comp_id in _PDB_EXCLUDE:
                    continue
                # Get SMILES
                chem_resp = await client.get(
                    f"https://data.rcsb.org/rest/v1/core/chemcomp/{comp_id}",
                    timeout=10,
                )
                if chem_resp.status_code != 200:
                    continue
                chem_data = chem_resp.json()
                descriptors = chem_data.get("rcsb_chem_comp_descriptor", {})
                smiles = None
                if isinstance(descriptors, dict):
                    # Flat object format: {smiles, smilesstereo, ...}
                    smiles = descriptors.get("smilesstereo") or descriptors.get("smiles")
                elif isinstance(descriptors, list):
                    # Legacy list format: [{type, descriptor}, ...]
                    for d in descriptors:
                        if d.get("type") == "SMILES_CANONICAL":
                            smiles = d.get("descriptor")
                            break
                    if not smiles:
                        for d in descriptors:
                            if d.get("type") == "SMILES":
                                smiles = d.get("descriptor")
                                break
                if smiles:
                    chem_name = (chem_data.get("chem_comp", {}).get("name") or "").strip()
                    is_primary = (comp_id == known_ligand)
                    # Build a human-readable name: "Imatinib (STI)" or fallback to comp_id
                    display_name = f"{chem_name.capitalize()} ({comp_id})" if chem_name else comp_id
                    suggestions.append({
                        "smiles": smiles,
                        "name": display_name,
                        "source": "cocrystal",
                        "source_label": "Co-crystallized ligand" if is_primary else "PDB ligand",
                        "comp_id": comp_id,
                        "priority": 0 if is_primary else 1,
                    })
            except Exception:
                continue
    except Exception as e:
        logger.warning("RCSB ligand fetch failed for %s: %s", pdb_id, e)
    # Sort: primary co-crystal first
    suggestions.sort(key=lambda x: x.get("priority", 1))
    return suggestions


async def _fetch_chembl_actives(client: httpx.AsyncClient, uniprot_id: str, limit: int = 3) -> list[dict]:
    """Fetch top active compounds from ChEMBL by pChEMBL value."""
    suggestions = []
    try:
        # Step 1: get ChEMBL target ID
        target_resp = await client.get(
            "https://www.ebi.ac.uk/chembl/api/data/target",
            params={"target_components__accession": uniprot_id, "format": "json"},
            timeout=15,
        )
        if target_resp.status_code != 200:
            return []
        targets = target_resp.json().get("targets", [])
        chembl_target_id = None
        for t in targets:
            if t.get("target_type") == "SINGLE PROTEIN":
                chembl_target_id = t["target_chembl_id"]
                break
        if not chembl_target_id:
            return []

        # Step 2: get top activities
        act_resp = await client.get(
            "https://www.ebi.ac.uk/chembl/api/data/activity",
            params={
                "target_chembl_id": chembl_target_id,
                "standard_type__in": "IC50,Ki",
                "pchembl_value__isnull": "false",
                "order_by": "-pchembl_value",
                "limit": 20,
                "format": "json",
            },
            timeout=15,
        )
        if act_resp.status_code != 200:
            return []

        activities = act_resp.json().get("activities", [])
        seen_smiles = set()
        for act in activities:
            smi = act.get("canonical_smiles")
            if not smi or smi in seen_smiles:
                continue
            seen_smiles.add(smi)
            mol_id = act.get("molecule_chembl_id", "")
            pchembl = act.get("pchembl_value", "")
            act_type = act.get("standard_type", "IC50")
            act_val = act.get("standard_value", "")
            act_units = act.get("standard_units", "nM")
            suggestions.append({
                "smiles": smi,
                "name": f"{mol_id} ({act_type}={act_val} {act_units})",
                "source": "chembl",
                "source_label": f"ChEMBL active (pChEMBL {pchembl})",
            })
            if len(suggestions) >= limit:
                break
    except Exception as e:
        logger.warning("ChEMBL fetch failed for %s: %s", uniprot_id, e)
    return suggestions


@router.get("/projects/{project_id}/suggest-ligands")
async def suggest_ligands(
    project_id: UUID,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Suggest reference ligands from UniProt (endogenous), RCSB (co-crystal), and ChEMBL (actives)."""
    await verify_project_ownership(project_id, user_id, db)

    stmt = select(ProjectORM_V9).where(ProjectORM_V9.id == project_id)
    result = await db.execute(stmt)
    project = result.scalar_one_or_none()
    if not project:
        raise HTTPException(404, "Project not found")

    pdb_id = project.target_pdb_id
    cocrystal = project.cocrystal_ligand
    uniprot_id = project.target_input_value if project.target_input_type == "uniprot" else None

    suggestions = []

    async def _noop():
        return []

    async with httpx.AsyncClient() as client:
        results = await asyncio.gather(
            _fetch_cocrystal_ligands(client, pdb_id, cocrystal) if pdb_id else _noop(),
            _fetch_uniprot_ligands(client, uniprot_id) if uniprot_id else _noop(),
            _fetch_chembl_actives(client, uniprot_id, limit=3) if uniprot_id else _noop(),
            return_exceptions=True,
        )

    for r in results:
        if isinstance(r, list):
            suggestions.extend(r)

    # Deduplicate by SMILES
    seen = set()
    unique = []
    for s in suggestions:
        if s["smiles"] not in seen:
            seen.add(s["smiles"])
            unique.append(s)

    return {"suggestions": unique}
