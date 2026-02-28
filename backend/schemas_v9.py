"""
BindX V9 — Pydantic v2 request/response schemas.

Covers Project, Campaign, Phase CRUD (Epic 2).
"""

from __future__ import annotations

import uuid
from datetime import datetime
from typing import List, Optional

from pydantic import BaseModel, ConfigDict


# ---------------------------------------------------------------------------
# Campaign schemas (declared first to avoid forward-ref issues)
# ---------------------------------------------------------------------------

class PhaseSummary(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: uuid.UUID
    campaign_id: uuid.UUID
    type: str
    status: str
    frozen_at: Optional[datetime] = None
    created_at: datetime


class CampaignSummary(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: uuid.UUID
    project_id: uuid.UUID
    name: str
    pocket_config: Optional[dict] = None
    strategy_notes: Optional[str] = None
    created_at: datetime
    phases: List[PhaseSummary] = []


# ---------------------------------------------------------------------------
# Project schemas
# ---------------------------------------------------------------------------

class ProjectCreate(BaseModel):
    name: str
    target_input_type: Optional[str] = "uniprot"
    target_input_value: Optional[str] = ""
    description: Optional[str] = None
    target_name: Optional[str] = None
    target_pdb_id: Optional[str] = None


class ProjectUpdate(BaseModel):
    name: Optional[str] = None
    description: Optional[str] = None
    target_name: Optional[str] = None
    target_input_type: Optional[str] = None
    target_input_value: Optional[str] = None
    target_pdb_id: Optional[str] = None
    target_preview: Optional[dict] = None
    structure_source: Optional[str] = None
    structure_resolution: Optional[float] = None
    structure_method: Optional[str] = None
    cocrystal_ligand: Optional[str] = None
    pockets_detected: Optional[list] = None
    chembl_actives_count: Optional[int] = None
    chembl_median_ic50: Optional[float] = None
    status: Optional[str] = None
    notification_email: Optional[str] = None


class ProjectResponse(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: uuid.UUID
    user_id: uuid.UUID
    name: str
    description: Optional[str] = None
    target_input_type: Optional[str] = None
    target_input_value: Optional[str] = None
    target_name: Optional[str] = None
    target_pdb_id: Optional[str] = None
    target_preview: Optional[dict] = None
    structure_source: Optional[str] = None
    structure_resolution: Optional[float] = None
    structure_method: Optional[str] = None
    cocrystal_ligand: Optional[str] = None
    pockets_detected: Optional[list] = None
    chembl_actives_count: Optional[int] = None
    chembl_median_ic50: Optional[float] = None
    status: str
    notification_email: Optional[str] = None
    created_at: datetime
    updated_at: datetime
    campaigns: List[CampaignSummary] = []


class ProjectListItem(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: uuid.UUID
    name: str
    description: Optional[str] = None
    target_input_type: Optional[str] = None
    target_input_value: Optional[str] = None
    target_name: Optional[str] = None
    target_pdb_id: Optional[str] = None
    target_preview: Optional[dict] = None
    structure_source: Optional[str] = None
    structure_resolution: Optional[float] = None
    structure_method: Optional[str] = None
    pockets_detected: Optional[list] = None
    chembl_actives_count: Optional[int] = None
    status: str
    created_at: datetime
    updated_at: datetime
    campaigns: List[CampaignSummary] = []


# ---------------------------------------------------------------------------
# Campaign CRUD schemas
# ---------------------------------------------------------------------------

class CampaignCreate(BaseModel):
    name: str
    pocket_config: Optional[dict] = None
    strategy_notes: Optional[str] = None


class CampaignUpdate(BaseModel):
    name: Optional[str] = None
    pocket_config: Optional[dict] = None
    strategy_notes: Optional[str] = None


class CampaignResponse(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: uuid.UUID
    project_id: uuid.UUID
    name: str
    pocket_config: Optional[dict] = None
    strategy_notes: Optional[str] = None
    created_at: datetime
    phases: List[PhaseSummary] = []


# ---------------------------------------------------------------------------
# Phase CRUD schemas
# ---------------------------------------------------------------------------

class PhaseCreate(BaseModel):
    type: str  # hit_discovery | hit_to_lead | lead_optimization


class PhaseUpdate(BaseModel):
    status: Optional[str] = None  # active | completed


class PhaseResponse(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: uuid.UUID
    campaign_id: uuid.UUID
    type: str
    status: str
    frozen_at: Optional[datetime] = None
    created_at: datetime
    molecule_count: int = 0
    bookmarked_count: int = 0
    run_count: int = 0


class FreezeResponse(BaseModel):
    id: uuid.UUID
    status: str
    frozen_at: Optional[datetime] = None
    bookmarked_molecule_count: int = 0


class UnfreezeResponse(BaseModel):
    id: uuid.UUID
    status: str
    warning: Optional[str] = None


# ---------------------------------------------------------------------------
# Run config schemas
# ---------------------------------------------------------------------------

class ImportRunConfig(BaseModel):
    """Config for import runs."""
    smiles_list: Optional[List[str]] = None
    file_format: Optional[str] = None  # "sdf" | "smiles" | "csv"
    pre_filters: Optional[dict] = None


class CalculationRunConfig(BaseModel):
    """Config for calculation runs."""
    docking_engine: Optional[str] = "gnina"
    scoring_weights: Optional[dict] = None
    exhaustiveness: Optional[int] = 8


class GenerationRunConfig(BaseModel):
    """Config for generation runs."""
    n_molecules: Optional[int] = 100
    strategy: Optional[str] = "scaffold_hop"


# ---------------------------------------------------------------------------
# Run schemas
# ---------------------------------------------------------------------------

class RunCreate(BaseModel):
    type: str  # "import" | "calculation" | "generation"
    config: Optional[dict] = {}
    input_molecule_ids: Optional[List[uuid.UUID]] = None
    calculation_types: Optional[List[str]] = None  # ["docking", "admet", "scoring", ...]


class RunResponse(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: uuid.UUID
    phase_id: uuid.UUID
    type: str
    calculation_types: Optional[List[str]] = None
    status: str
    config: dict
    input_molecule_ids: Optional[List[uuid.UUID]] = None
    input_source: Optional[str] = None
    input_file_path: Optional[str] = None
    progress: int = 0
    current_step: Optional[str] = None
    estimated_duration_seconds: Optional[int] = None
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    error_message: Optional[str] = None
    archived: bool = False
    created_at: datetime


class RunListItem(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: uuid.UUID
    phase_id: uuid.UUID
    type: str
    status: str
    progress: int = 0
    current_step: Optional[str] = None
    archived: bool = False
    created_at: datetime


# ---------------------------------------------------------------------------
# Molecule schemas
# ---------------------------------------------------------------------------

class MoleculeResponse(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: uuid.UUID
    phase_id: uuid.UUID
    smiles: str
    canonical_smiles: str
    name: Optional[str] = None
    source_run_id: Optional[uuid.UUID] = None
    bookmarked: bool = False
    generation_level: int = 0
    parent_molecule_id: Optional[uuid.UUID] = None
    ai_generated: bool = False
    created_at: datetime
    properties: Optional[dict] = None  # {prop_name: value}


class MoleculeListResponse(BaseModel):
    molecules: List[MoleculeResponse]
    total: int
    offset: int = 0
    limit: int = 50


class BookmarkRequest(BaseModel):
    molecule_ids: List[uuid.UUID]
    bookmarked: bool
