"""
BindX V9 — ORM Models (SQLAlchemy 2.0 Mapped style).

9 models matching the V9 Supabase schema.
user_id stored as UUID without ORM FK (auth.users is managed by Supabase).
"""

from __future__ import annotations

import uuid
from datetime import datetime
from typing import List, Optional

from sqlalchemy import (
    Boolean,
    Float,
    ForeignKey,
    Integer,
    Text,
    UniqueConstraint,
)
from sqlalchemy.dialects.postgresql import ARRAY, JSONB, UUID, TIMESTAMP
from sqlalchemy.orm import Mapped, mapped_column, relationship

from database_v9 import BaseV9


# ---------------------------------------------------------------------------
# 1. projects
# ---------------------------------------------------------------------------

class ProjectORM_V9(BaseV9):
    __tablename__ = "projects"

    id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    user_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), nullable=False)
    name: Mapped[str] = mapped_column(Text, nullable=False)
    description: Mapped[Optional[str]] = mapped_column(Text)
    # Target info
    target_input_type: Mapped[Optional[str]] = mapped_column(Text)
    target_input_value: Mapped[Optional[str]] = mapped_column(Text, default="")
    target_name: Mapped[Optional[str]] = mapped_column(Text)
    target_pdb_id: Mapped[Optional[str]] = mapped_column(Text)
    target_preview: Mapped[Optional[dict]] = mapped_column(JSONB)
    # Structure source info
    structure_source: Mapped[Optional[str]] = mapped_column(Text)
    structure_resolution: Mapped[Optional[float]] = mapped_column(Float)
    structure_method: Mapped[Optional[str]] = mapped_column(Text)
    cocrystal_ligand: Mapped[Optional[str]] = mapped_column(Text)
    # Pockets
    pockets_detected: Mapped[Optional[list]] = mapped_column(JSONB)
    # ChEMBL
    chembl_actives_count: Mapped[Optional[int]] = mapped_column(Integer)
    chembl_median_ic50: Mapped[Optional[float]] = mapped_column(Float)
    #
    status: Mapped[str] = mapped_column(Text, default="active")
    notification_email: Mapped[Optional[str]] = mapped_column(Text)
    created_at: Mapped[datetime] = mapped_column(TIMESTAMP(timezone=True), default=datetime.utcnow)
    updated_at: Mapped[datetime] = mapped_column(TIMESTAMP(timezone=True), default=datetime.utcnow, onupdate=datetime.utcnow)

    # Relationships
    campaigns: Mapped[List[CampaignORM_V9]] = relationship(back_populates="project", cascade="all, delete-orphan")


# ---------------------------------------------------------------------------
# 2. campaigns
# ---------------------------------------------------------------------------

class CampaignORM_V9(BaseV9):
    __tablename__ = "campaigns"

    id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    project_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), ForeignKey("projects.id", ondelete="CASCADE"), nullable=False)
    name: Mapped[str] = mapped_column(Text, nullable=False)
    pocket_config: Mapped[Optional[dict]] = mapped_column(JSONB)
    strategy_notes: Mapped[Optional[str]] = mapped_column(Text)
    created_at: Mapped[datetime] = mapped_column(TIMESTAMP(timezone=True), default=datetime.utcnow)

    # Relationships
    project: Mapped[ProjectORM_V9] = relationship(back_populates="campaigns")
    phases: Mapped[List[PhaseORM_V9]] = relationship(back_populates="campaign", cascade="all, delete-orphan")


# ---------------------------------------------------------------------------
# 3. phases
# ---------------------------------------------------------------------------

class PhaseORM_V9(BaseV9):
    __tablename__ = "phases"

    id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    campaign_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), ForeignKey("campaigns.id", ondelete="CASCADE"), nullable=False)
    type: Mapped[str] = mapped_column(Text, nullable=False)
    status: Mapped[str] = mapped_column(Text, default="active")
    frozen_at: Mapped[Optional[datetime]] = mapped_column(TIMESTAMP(timezone=True))
    created_at: Mapped[datetime] = mapped_column(TIMESTAMP(timezone=True), default=datetime.utcnow)

    # Relationships
    campaign: Mapped[CampaignORM_V9] = relationship(back_populates="phases")
    runs: Mapped[List[RunORM_V9]] = relationship(back_populates="phase", cascade="all, delete-orphan")
    molecules: Mapped[List[MoleculeORM_V9]] = relationship(back_populates="phase", cascade="all, delete-orphan")


# ---------------------------------------------------------------------------
# 4. runs
# ---------------------------------------------------------------------------

class RunORM_V9(BaseV9):
    __tablename__ = "runs"

    id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    phase_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), ForeignKey("phases.id", ondelete="CASCADE"), nullable=False)
    type: Mapped[str] = mapped_column(Text, nullable=False)
    calculation_types: Mapped[Optional[list]] = mapped_column(ARRAY(Text))
    status: Mapped[str] = mapped_column(Text, default="created")
    config: Mapped[dict] = mapped_column(JSONB, nullable=False)
    input_molecule_ids: Mapped[Optional[list]] = mapped_column(ARRAY(UUID(as_uuid=True)))
    input_source: Mapped[Optional[str]] = mapped_column(Text)
    input_file_path: Mapped[Optional[str]] = mapped_column(Text)
    progress: Mapped[int] = mapped_column(Integer, default=0)
    current_step: Mapped[Optional[str]] = mapped_column(Text)
    estimated_duration_seconds: Mapped[Optional[int]] = mapped_column(Integer)
    started_at: Mapped[Optional[datetime]] = mapped_column(TIMESTAMP(timezone=True))
    completed_at: Mapped[Optional[datetime]] = mapped_column(TIMESTAMP(timezone=True))
    error_message: Mapped[Optional[str]] = mapped_column(Text)
    archived: Mapped[bool] = mapped_column(Boolean, default=False)
    created_at: Mapped[datetime] = mapped_column(TIMESTAMP(timezone=True), default=datetime.utcnow)

    # Relationships
    phase: Mapped[PhaseORM_V9] = relationship(back_populates="runs")
    artifacts: Mapped[List[ArtifactORM_V9]] = relationship(back_populates="run", cascade="all, delete-orphan")


# ---------------------------------------------------------------------------
# 5. molecules
# ---------------------------------------------------------------------------

class MoleculeORM_V9(BaseV9):
    __tablename__ = "molecules"
    __table_args__ = (
        UniqueConstraint("phase_id", "canonical_smiles", name="uq_molecules_phase_canonical"),
    )

    id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    phase_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), ForeignKey("phases.id", ondelete="CASCADE"), nullable=False)
    smiles: Mapped[str] = mapped_column(Text, nullable=False)
    canonical_smiles: Mapped[str] = mapped_column(Text, nullable=False)
    name: Mapped[Optional[str]] = mapped_column(Text)
    source_run_id: Mapped[Optional[uuid.UUID]] = mapped_column(UUID(as_uuid=True), ForeignKey("runs.id"))
    bookmarked: Mapped[bool] = mapped_column(Boolean, default=False)
    generation_level: Mapped[int] = mapped_column(Integer, default=0)
    parent_molecule_id: Mapped[Optional[uuid.UUID]] = mapped_column(UUID(as_uuid=True), ForeignKey("molecules.id"))
    ai_generated: Mapped[bool] = mapped_column(Boolean, default=False)
    created_at: Mapped[datetime] = mapped_column(TIMESTAMP(timezone=True), default=datetime.utcnow)

    # Relationships
    phase: Mapped[PhaseORM_V9] = relationship(back_populates="molecules")
    properties: Mapped[List[MoleculePropertyORM_V9]] = relationship(back_populates="molecule", cascade="all, delete-orphan")


# ---------------------------------------------------------------------------
# 6. molecule_properties
# ---------------------------------------------------------------------------

class MoleculePropertyORM_V9(BaseV9):
    __tablename__ = "molecule_properties"
    __table_args__ = (
        UniqueConstraint("molecule_id", "property_name", "run_id", name="uq_mol_props_mol_name_run"),
    )

    id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    molecule_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), ForeignKey("molecules.id", ondelete="CASCADE"), nullable=False)
    run_id: Mapped[Optional[uuid.UUID]] = mapped_column(UUID(as_uuid=True), ForeignKey("runs.id"))
    property_name: Mapped[str] = mapped_column(Text, nullable=False)
    property_value: Mapped[dict] = mapped_column(JSONB, nullable=False)
    created_at: Mapped[datetime] = mapped_column(TIMESTAMP(timezone=True), default=datetime.utcnow)

    # Relationships
    molecule: Mapped[MoleculeORM_V9] = relationship(back_populates="properties")


# ---------------------------------------------------------------------------
# 7. calculation_cache
# ---------------------------------------------------------------------------

class CalculationCacheORM_V9(BaseV9):
    __tablename__ = "calculation_cache"

    id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    cache_key: Mapped[str] = mapped_column(Text, unique=True, nullable=False)
    result: Mapped[dict] = mapped_column(JSONB, nullable=False)
    created_at: Mapped[datetime] = mapped_column(TIMESTAMP(timezone=True), default=datetime.utcnow)


# ---------------------------------------------------------------------------
# 8. artifacts
# ---------------------------------------------------------------------------

class ArtifactORM_V9(BaseV9):
    __tablename__ = "artifacts"

    id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    run_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), ForeignKey("runs.id", ondelete="CASCADE"), nullable=False)
    type: Mapped[str] = mapped_column(Text, nullable=False)
    storage_path: Mapped[str] = mapped_column(Text, nullable=False)
    created_at: Mapped[datetime] = mapped_column(TIMESTAMP(timezone=True), default=datetime.utcnow)

    # Relationships
    run: Mapped[RunORM_V9] = relationship(back_populates="artifacts")


# ---------------------------------------------------------------------------
# 9. audit_log
# ---------------------------------------------------------------------------

class AuditLogORM_V9(BaseV9):
    __tablename__ = "audit_log"

    id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    user_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), nullable=False)
    action: Mapped[str] = mapped_column(Text, nullable=False)
    entity_type: Mapped[Optional[str]] = mapped_column(Text)
    entity_id: Mapped[Optional[uuid.UUID]] = mapped_column(UUID(as_uuid=True))
    details: Mapped[Optional[dict]] = mapped_column(JSONB)
    created_at: Mapped[datetime] = mapped_column(TIMESTAMP(timezone=True), default=datetime.utcnow)
