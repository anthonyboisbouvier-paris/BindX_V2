"""
BindX V9 — Pydantic v2 request/response schemas.

Covers Project, Campaign, Phase CRUD (Epic 2).
"""

from __future__ import annotations

import uuid
from datetime import datetime
from typing import List, Optional

from pydantic import BaseModel, ConfigDict, field_validator


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
    stats: Optional[dict] = None


class CampaignSummary(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: uuid.UUID
    project_id: uuid.UUID
    name: str
    pocket_config: Optional[dict] = None
    strategy_notes: Optional[str] = None
    reference_ligands: Optional[list] = []
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
    pubchem_compounds_count: Optional[int] = None
    receptor_prep_config: Optional[dict] = None
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
    pubchem_compounds_count: Optional[int] = None
    receptor_prep_config: Optional[dict] = None
    receptor_prep_report: Optional[dict] = None
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
    pubchem_compounds_count: Optional[int] = None
    receptor_prep_config: Optional[dict] = None
    receptor_prep_report: Optional[dict] = None
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
    reference_ligands: Optional[list] = []
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

class RunLogEntryOut(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: uuid.UUID
    level: str = "info"
    message: str
    created_at: datetime


class RunCreate(BaseModel):
    type: str  # "import" | "calculation" | "generation"
    config: Optional[dict] = {}
    input_molecule_ids: Optional[List[uuid.UUID]] = None
    calculation_types: Optional[List[str]] = None  # ["docking", "admet", "scoring", ...]


def _clean_config(v):
    """Strip internal keys (prefixed with _) from run config."""
    if isinstance(v, dict):
        return {k: val for k, val in v.items() if not k.startswith("_")}
    return v


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
    progress: int = 0
    current_step: Optional[str] = None
    estimated_duration_seconds: Optional[int] = None
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    error_message: Optional[str] = None
    archived: bool = False
    created_at: datetime
    logs: List[RunLogEntryOut] = []

    @field_validator("config", mode="before")
    @classmethod
    def _clean(cls, v):
        return _clean_config(v)


class RunListItem(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: uuid.UUID
    phase_id: uuid.UUID
    type: str
    calculation_types: Optional[List[str]] = None
    status: str
    config: dict = {}
    progress: int = 0
    current_step: Optional[str] = None
    error_message: Optional[str] = None
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    archived: bool = False
    created_at: datetime
    logs: List[RunLogEntryOut] = []

    @field_validator("config", mode="before")
    @classmethod
    def _clean(cls, v):
        return _clean_config(v)


# ---------------------------------------------------------------------------
# Molecule schemas
# ---------------------------------------------------------------------------

class MoleculeAnnotationUpdate(BaseModel):
    """Partial update for molecule annotations (user_comment, tags, invalidated)."""
    user_comment: Optional[str] = None
    tags: Optional[List[str]] = None
    invalidated: Optional[bool] = None


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
    user_comment: Optional[str] = None
    ai_comment: Optional[str] = None
    tags: Optional[List[str]] = None
    invalidated: bool = False
    created_at: datetime
    properties: Optional[dict] = None  # {prop_name: value}


class MoleculeListResponse(BaseModel):
    molecules: List[MoleculeResponse]
    total: int
    offset: int = 0
    limit: int = 50
    next_cursor: Optional[str] = None
    has_more: bool = False


class MoleculeStatsResponse(BaseModel):
    total: int = 0
    bookmarked: int = 0
    ai_generated: int = 0
    column_ranges: dict = {}  # { col_key: { min, max } }


class BookmarkRequest(BaseModel):
    molecule_ids: List[uuid.UUID]
    bookmarked: bool


# ---------------------------------------------------------------------------
# SAR / MMP schemas
# ---------------------------------------------------------------------------

class MMPApplyRequest(BaseModel):
    from_smiles: str
    to_smiles: str
    molecule_ids: Optional[List[uuid.UUID]] = None
    max_results: int = 50


class MMPSuggestion(BaseModel):
    smiles: str
    source_mol_id: Optional[str] = None
    source_mol_name: Optional[str] = None
    qed: Optional[float] = None
    mw: Optional[float] = None
    logp: Optional[float] = None
    hbd: Optional[int] = None
    hba: Optional[int] = None
    valid: bool = True


class MMPApplyResponse(BaseModel):
    suggestions: List[MMPSuggestion]
    transform: str
    n_input_molecules: int


class MMPImportRequest(BaseModel):
    smiles_list: List[str]
    names: Optional[List[str]] = None
    transform: Optional[str] = None


# ---------------------------------------------------------------------------
# AFVS schemas (BDX-41)
# ---------------------------------------------------------------------------

class AFVSRunCreate(BaseModel):
    """Config for launching an AFVS ultra-large screening run."""
    strategy: str = "balanced"  # lean | balanced | aggressive
    docking_scenario: str = "smina_vinardo"  # smina_vinardo | qvina2
    exhaustiveness_prescreen: int = 1
    exhaustiveness_primary: int = 8
    reps_per_tranche: int = 1
    tranche_pct: float = 0.5  # percentage of tranches for primary screen
    mw_max: float = 500
    logp_max: float = 5
    hbd_max: int = 5
    hba_max: int = 10
    top_n_results: int = 10000
    budget_cap_usd: float = 500

    @field_validator("top_n_results")
    @classmethod
    def _cap_top_n(cls, v):
        if v > 100_000:
            raise ValueError("top_n_results cannot exceed 100,000")
        return v

    @field_validator("budget_cap_usd")
    @classmethod
    def _budget_positive(cls, v):
        if v <= 0:
            raise ValueError("budget_cap_usd must be positive")
        return v


class AFVSStatusResponse(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    run_id: uuid.UUID
    status: str
    afvs_phase: Optional[str] = None
    progress: int = 0
    current_step: Optional[str] = None
    molecules_screened: Optional[int] = None
    molecules_imported: Optional[int] = None
    actual_cost_usd: Optional[float] = None
    budget_cap_usd: Optional[float] = None
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    error_message: Optional[str] = None
