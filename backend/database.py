"""
DockIt — Database setup (SQLite via SQLAlchemy).

Provides engine, session factory, and helper functions for job persistence.
V3: supports sequence, notification_email, docking_engine fields.
"""

from __future__ import annotations

import json
import logging
import os
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from typing import Generator, Optional

from sqlalchemy import create_engine
from sqlalchemy.orm import Session, sessionmaker

from models import Base, JobORM, ProjectORM, TargetAssessmentORM, UserORM

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Engine / Session
# ---------------------------------------------------------------------------

DATA_DIR = Path(os.environ.get("DOCKIT_DATA_DIR", "/data"))
DATA_DIR.mkdir(parents=True, exist_ok=True)

DATABASE_URL = f"sqlite:///{DATA_DIR / 'dockit.db'}"

engine = create_engine(
    DATABASE_URL,
    connect_args={"check_same_thread": False},
    echo=False,
)

SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)


def init_db() -> None:
    """Create all tables if they do not exist, and migrate existing tables."""
    Base.metadata.create_all(bind=engine)

    # V7 migration: add new columns to existing jobs table if missing
    from sqlalchemy import text, inspect
    insp = inspect(engine)
    if "jobs" in insp.get_table_names():
        existing_cols = {c["name"] for c in insp.get_columns("jobs")}
        with engine.begin() as conn:
            if "project_id" not in existing_cols:
                conn.execute(text("ALTER TABLE jobs ADD COLUMN project_id VARCHAR(36)"))
                logger.info("Migrated: added project_id column to jobs table")
            if "user_id" not in existing_cols:
                conn.execute(text("ALTER TABLE jobs ADD COLUMN user_id VARCHAR(36)"))
                logger.info("Migrated: added user_id column to jobs table")

    logger.info("Database initialised at %s", DATABASE_URL)


@contextmanager
def get_db() -> Generator[Session, None, None]:
    """Context manager that yields a SQLAlchemy session and handles cleanup."""
    db = SessionLocal()
    try:
        yield db
        db.commit()
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()


# ---------------------------------------------------------------------------
# Convenience helpers
# ---------------------------------------------------------------------------

def create_job(
    job_id: str,
    uniprot_id: Optional[str] = None,
    use_chembl: bool = True,
    use_zinc: bool = False,
    max_ligands: int = 50,
    smiles_list: Optional[list[str]] = None,
    mode: str = "rapid",
    enable_generation: bool = False,
    enable_diffdock: bool = False,
    enable_retrosynthesis: bool = False,
    n_generated_molecules: int = 100,
    # V3 fields
    sequence: Optional[str] = None,
    notification_email: Optional[str] = None,
    docking_engine: str = "vina",
    # V7 fields
    project_id: Optional[str] = None,
    user_id: Optional[str] = None,
) -> JobORM:
    """Insert a new job row and return it.

    Parameters
    ----------
    job_id : str
        UUID4 identifier for the job.
    uniprot_id : str, optional
        UniProt accession. Can be None if sequence is provided.
    use_chembl : bool
        Whether to query ChEMBL for known ligands.
    use_zinc : bool
        Whether to include ZINC drug-like molecules.
    max_ligands : int
        Maximum number of ligands to screen.
    smiles_list : list[str], optional
        User-supplied SMILES strings.
    mode : str
        Pipeline mode: rapid, standard, or deep.
    enable_generation : bool
        Enable AI molecule generation (V2 compat).
    enable_diffdock : bool
        Enable DiffDock (V2 compat).
    enable_retrosynthesis : bool
        Enable retrosynthesis planning (V2 compat).
    n_generated_molecules : int
        Number of molecules to generate (V2 compat).
    sequence : str, optional
        Raw amino acid sequence (V3).
    notification_email : str, optional
        Email for completion notification (V3).
    docking_engine : str
        Docking engine to use, default "vina" (V3).
    project_id : str, optional
        Project to associate with (V7).
    user_id : str, optional
        User who owns this job (V7).

    Returns
    -------
    JobORM
        The created job record.
    """
    with get_db() as db:
        job = JobORM(
            id=job_id,
            uniprot_id=uniprot_id or "",
            status="queued",
            progress=0,
            current_step="Queued",
            use_chembl=int(use_chembl),
            use_zinc=int(use_zinc),
            max_ligands=max_ligands,
            smiles_list=json.dumps(smiles_list) if smiles_list else None,
            mode=mode,
            enable_generation=int(enable_generation),
            enable_diffdock=int(enable_diffdock),
            enable_retrosynthesis=int(enable_retrosynthesis),
            n_generated_molecules=n_generated_molecules,
            # V3 fields
            sequence=sequence,
            notification_email=notification_email,
            docking_engine=docking_engine,
            # V7 fields
            project_id=project_id,
            user_id=user_id,
            created_at=datetime.utcnow(),
        )
        db.add(job)
        db.flush()
        db.expunge(job)
        display_id = uniprot_id or f"seq[{len(sequence or '')}aa]"
        logger.info("Job %s created for %s (mode=%s, engine=%s)", job_id, display_id, mode, docking_engine)
    return job


def get_job(job_id: str) -> Optional[JobORM]:
    """Fetch a job by ID, or return None."""
    with get_db() as db:
        job = db.get(JobORM, job_id)
        if job is not None:
            # Detach from session so caller can use fields freely
            db.expunge(job)
        return job


def list_jobs(
    limit: int = 50,
    project_id: Optional[str] = None,
    user_id: Optional[str] = None,
) -> list[JobORM]:
    """Return the most recent jobs, ordered by created_at descending.

    Parameters
    ----------
    limit : int
        Maximum number of jobs to return (default 50).
    project_id : str, optional
        Filter by project ID (V7).
    user_id : str, optional
        Filter by user ID (V7).

    Returns
    -------
    list[JobORM]
        List of detached JobORM instances.
    """
    with get_db() as db:
        query = db.query(JobORM)
        if project_id is not None:
            query = query.filter(JobORM.project_id == project_id)
        if user_id is not None:
            query = query.filter(JobORM.user_id == user_id)
        jobs = (
            query
            .order_by(JobORM.created_at.desc())
            .limit(limit)
            .all()
        )
        for job in jobs:
            db.expunge(job)
        return jobs


def update_job(job_id: str, **kwargs) -> None:
    """Update arbitrary columns on a job row."""
    with get_db() as db:
        job = db.get(JobORM, job_id)
        if job is None:
            logger.warning("update_job: job %s not found", job_id)
            return
        for key, value in kwargs.items():
            setattr(job, key, value)
        db.flush()


# ---------------------------------------------------------------------------
# V7: User helpers
# ---------------------------------------------------------------------------

def create_user(
    user_id: str,
    email: str,
    username: str,
    password_hash: str,
) -> UserORM:
    """Insert a new user row and return it.

    Parameters
    ----------
    user_id : str
        UUID4 identifier for the user.
    email : str
        User email (unique).
    username : str
        Display name.
    password_hash : str
        Bcrypt hash of the password.

    Returns
    -------
    UserORM
        The created user record.
    """
    with get_db() as db:
        user = UserORM(
            id=user_id,
            email=email,
            username=username,
            password_hash=password_hash,
            created_at=datetime.utcnow(),
        )
        db.add(user)
        db.flush()
        db.expunge(user)
        logger.info("User %s created (email=%s)", user_id, email)
    return user


def get_user(user_id: str) -> Optional[UserORM]:
    """Fetch a user by ID, or return None."""
    with get_db() as db:
        user = db.get(UserORM, user_id)
        if user is not None:
            db.expunge(user)
        return user


def get_user_by_email(email: str) -> Optional[UserORM]:
    """Fetch a user by email address, or return None.

    Parameters
    ----------
    email : str
        The email address to look up.

    Returns
    -------
    UserORM or None
        The user record, or None if not found.
    """
    with get_db() as db:
        user = db.query(UserORM).filter(UserORM.email == email).first()
        if user is not None:
            db.expunge(user)
        return user


# ---------------------------------------------------------------------------
# V7: Project helpers
# ---------------------------------------------------------------------------

def create_project(
    project_id: str,
    user_id: str,
    name: str,
    uniprot_id: Optional[str] = None,
    sequence: Optional[str] = None,
    description: Optional[str] = None,
    target_preview_json: Optional[str] = None,
) -> ProjectORM:
    """Insert a new project row and return it.

    Parameters
    ----------
    project_id : str
        UUID4 identifier for the project.
    user_id : str
        UUID of the owning user.
    name : str
        Project display name.
    uniprot_id : str, optional
        Target UniProt accession.
    sequence : str, optional
        Target amino acid sequence.
    description : str, optional
        Free-text description.
    target_preview_json : str, optional
        JSON-encoded preview data.

    Returns
    -------
    ProjectORM
        The created project record.
    """
    with get_db() as db:
        project = ProjectORM(
            id=project_id,
            user_id=user_id,
            name=name,
            uniprot_id=uniprot_id,
            sequence=sequence,
            description=description,
            target_preview_json=target_preview_json,
            created_at=datetime.utcnow(),
        )
        db.add(project)
        db.flush()
        db.expunge(project)
        logger.info("Project %s created by user %s", project_id, user_id)
    return project


def get_project(project_id: str) -> Optional[ProjectORM]:
    """Fetch a project by ID, or return None."""
    with get_db() as db:
        project = db.get(ProjectORM, project_id)
        if project is not None:
            db.expunge(project)
        return project


def list_projects(user_id: str, limit: int = 50) -> list[ProjectORM]:
    """Return all projects for a user, ordered by created_at descending.

    Parameters
    ----------
    user_id : str
        UUID of the user.
    limit : int
        Maximum number of projects to return (default 50).

    Returns
    -------
    list[ProjectORM]
        List of detached ProjectORM instances.
    """
    with get_db() as db:
        projects = (
            db.query(ProjectORM)
            .filter(ProjectORM.user_id == user_id)
            .order_by(ProjectORM.created_at.desc())
            .limit(limit)
            .all()
        )
        for p in projects:
            db.expunge(p)
        return projects


def update_project(project_id: str, **kwargs) -> None:
    """Update arbitrary columns on a project row.

    Parameters
    ----------
    project_id : str
        UUID of the project.
    **kwargs
        Column name/value pairs to update.
    """
    with get_db() as db:
        project = db.get(ProjectORM, project_id)
        if project is None:
            logger.warning("update_project: project %s not found", project_id)
            return
        for key, value in kwargs.items():
            setattr(project, key, value)
        db.flush()


def count_project_jobs(project_id: str) -> int:
    """Count the number of jobs associated with a project.

    Parameters
    ----------
    project_id : str
        UUID of the project.

    Returns
    -------
    int
        Number of jobs.
    """
    with get_db() as db:
        count = db.query(JobORM).filter(JobORM.project_id == project_id).count()
        return count


# ---------------------------------------------------------------------------
# BindX: Target Assessment helpers
# ---------------------------------------------------------------------------

def create_assessment(
    assessment_id: str,
    uniprot_id: str,
    assessment_json: str,
    project_id: Optional[str] = None,
    disease_context: Optional[str] = None,
    agent_responses_json: Optional[str] = None,
    input_hash: Optional[str] = None,
) -> TargetAssessmentORM:
    """Insert a new target assessment row and return it."""
    with get_db() as db:
        row = TargetAssessmentORM(
            id=assessment_id,
            project_id=project_id,
            uniprot_id=uniprot_id,
            disease_context=disease_context,
            assessment_json=assessment_json,
            agent_responses_json=agent_responses_json,
            input_hash=input_hash,
            created_at=datetime.utcnow(),
        )
        db.add(row)
        db.flush()
        db.expunge(row)
        logger.info("Assessment %s created for %s", assessment_id, uniprot_id)
    return row


def get_assessment(assessment_id: str) -> Optional[TargetAssessmentORM]:
    """Fetch a target assessment by ID."""
    with get_db() as db:
        row = db.get(TargetAssessmentORM, assessment_id)
        if row is not None:
            db.expunge(row)
        return row


def get_assessment_by_hash(input_hash: str) -> Optional[TargetAssessmentORM]:
    """Fetch a cached assessment by input hash (idempotence)."""
    with get_db() as db:
        row = (
            db.query(TargetAssessmentORM)
            .filter(TargetAssessmentORM.input_hash == input_hash)
            .order_by(TargetAssessmentORM.created_at.desc())
            .first()
        )
        if row is not None:
            db.expunge(row)
        return row


def update_assessment_agent(assessment_id: str, agent_responses_json: str) -> None:
    """Update the agent_responses_json for an existing assessment."""
    with get_db() as db:
        row = db.query(TargetAssessmentORM).filter(
            TargetAssessmentORM.id == assessment_id
        ).first()
        if row:
            row.agent_responses_json = agent_responses_json


def list_project_jobs(project_id: str, limit: int = 50) -> list[JobORM]:
    """Return all jobs for a project, ordered by created_at descending.

    Parameters
    ----------
    project_id : str
        UUID of the project.
    limit : int
        Maximum number of jobs to return (default 50).

    Returns
    -------
    list[JobORM]
        List of detached JobORM instances.
    """
    with get_db() as db:
        jobs = (
            db.query(JobORM)
            .filter(JobORM.project_id == project_id)
            .order_by(JobORM.created_at.desc())
            .limit(limit)
            .all()
        )
        for job in jobs:
            db.expunge(job)
        return jobs
