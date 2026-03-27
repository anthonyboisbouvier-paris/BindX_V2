"""
BindX V9 — Celery tasks for run execution.

Entry point: execute_run(run_id) dispatches to _run_import, _run_calculation,
or _run_generation based on run type.

Uses a sync wrapper around the async DB engine for Celery worker context.
"""

from __future__ import annotations

import asyncio
import logging
import uuid
from datetime import datetime, timezone
from pathlib import Path

from celery_app import celery_app

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Sync wrapper for async DB operations in Celery
# ---------------------------------------------------------------------------
# IMPORTANT: asyncpg connections are bound to the event loop that created them.
# We must reuse a single loop per worker process so the connection pool works
# across multiple _db_sync() calls within a task.

_worker_loop = None


def _db_sync(coro):
    """Run an async DB coroutine in Celery sync context.

    Reuses a single event loop per worker process to avoid asyncpg
    'Future attached to a different loop' errors.
    """
    global _worker_loop
    if _worker_loop is None or _worker_loop.is_closed():
        _worker_loop = asyncio.new_event_loop()
    return _worker_loop.run_until_complete(coro)


# ---------------------------------------------------------------------------
# DB helpers
# ---------------------------------------------------------------------------

async def _get_run(run_id: str):
    """Fetch a RunORM_V9 by id."""
    from database_v9 import AsyncSessionV9
    from models_v9 import RunORM_V9

    async with AsyncSessionV9() as session:
        from sqlalchemy import select
        result = await session.execute(
            select(RunORM_V9).where(RunORM_V9.id == uuid.UUID(run_id))
        )
        run = result.scalar_one_or_none()
        if run:
            # Detach values we need
            return {
                "id": run.id,
                "phase_id": run.phase_id,
                "type": run.type,
                "config": run.config,
                "calculation_types": run.calculation_types,
                "input_molecule_ids": run.input_molecule_ids,
                "input_file_path": run.input_file_path,
                "input_source": run.input_source,
                "status": run.status,
            }
    return None


async def _update_run(run_id: str, **fields):
    """Update run fields."""
    from database_v9 import AsyncSessionV9
    from models_v9 import RunORM_V9
    from sqlalchemy import select

    async with AsyncSessionV9() as session:
        result = await session.execute(
            select(RunORM_V9).where(RunORM_V9.id == uuid.UUID(run_id))
        )
        run = result.scalar_one_or_none()
        if run:
            for k, v in fields.items():
                setattr(run, k, v)
            await session.commit()


async def _create_molecule(phase_id, smiles, canonical_smiles, name, source_run_id,
                           ai_generated=False, parent_molecule_id=None, generation_level=0):
    """Create a MoleculeORM_V9 if not duplicate in phase. Returns molecule id or None."""
    from database_v9 import AsyncSessionV9
    from models_v9 import MoleculeORM_V9
    from sqlalchemy import select

    async with AsyncSessionV9() as session:
        # Check dedup
        existing = await session.execute(
            select(MoleculeORM_V9.id).where(
                MoleculeORM_V9.phase_id == phase_id,
                MoleculeORM_V9.canonical_smiles == canonical_smiles,
            )
        )
        if existing.scalar_one_or_none():
            return None  # duplicate

        mol = MoleculeORM_V9(
            phase_id=phase_id,
            smiles=smiles,
            canonical_smiles=canonical_smiles,
            name=name,
            source_run_id=source_run_id,
            ai_generated=ai_generated,
            parent_molecule_id=parent_molecule_id,
            generation_level=generation_level,
        )
        session.add(mol)
        await session.commit()
        return mol.id


async def _create_molecules_batch(phase_id, molecules: list[dict], source_run_id):
    """Batch INSERT ... ON CONFLICT DO NOTHING.

    Returns (created_count, skipped_count, created_rows) where created_rows
    is a list of (id, canonical_smiles) for newly inserted molecules.

    Each molecule dict: {"smiles": str, "canonical_smiles": str, "name": str|None,
    plus optional metadata keys (mwt, logp, etc.)}
    Uses a single INSERT with N rows per call — dedup via DB constraint.
    """
    from database_v9 import AsyncSessionV9
    from sqlalchemy import text

    if not molecules:
        return 0, 0, []

    async with AsyncSessionV9() as session:
        # Build parameterized VALUES
        values_parts = []
        params = {"phase_id": phase_id, "source_run_id": source_run_id}
        for i, mol in enumerate(molecules):
            mol_id = uuid.uuid4()
            values_parts.append(
                f"(:id_{i}, :phase_id, :smi_{i}, :can_{i}, :name_{i}, :source_run_id)"
            )
            params[f"id_{i}"] = mol_id
            params[f"smi_{i}"] = mol["smiles"]
            params[f"can_{i}"] = mol["canonical_smiles"]
            params[f"name_{i}"] = mol.get("name")

        sql = text(
            "INSERT INTO molecules (id, phase_id, smiles, canonical_smiles, name, source_run_id) "
            f"VALUES {', '.join(values_parts)} "
            "ON CONFLICT (phase_id, canonical_smiles) DO NOTHING "
            "RETURNING id, canonical_smiles"
        )
        result = await session.execute(sql, params)
        rows = result.fetchall()
        created = len(rows)
        await session.commit()
        skipped = len(molecules) - created
        return created, skipped, rows


async def _store_property(molecule_id, run_id, prop_name, prop_value):
    """Store a molecule property."""
    from database_v9 import AsyncSessionV9
    from models_v9 import MoleculePropertyORM_V9

    async with AsyncSessionV9() as session:
        prop = MoleculePropertyORM_V9(
            molecule_id=molecule_id,
            run_id=run_id,
            property_name=prop_name,
            property_value=prop_value,
        )
        session.add(prop)
        await session.commit()


async def _store_properties_batch(entries: list[tuple], batch_size: int = 500):
    """Bulk-upsert molecule properties in batches.

    entries: list of (molecule_id, run_id, prop_name, prop_value) tuples.
    Uses PostgreSQL ON CONFLICT DO UPDATE to handle re-runs gracefully.
    Commits in chunks of batch_size to avoid oversized transactions.
    """
    from database_v9 import AsyncSessionV9
    from sqlalchemy.dialects.postgresql import insert as pg_insert
    from models_v9 import MoleculePropertyORM_V9

    if not entries:
        return

    table = MoleculePropertyORM_V9.__table__

    for i in range(0, len(entries), batch_size):
        chunk = entries[i : i + batch_size]
        rows = [
            {
                "molecule_id": mol_id,
                "run_id": run_id,
                "property_name": prop_name,
                "property_value": prop_value,
            }
            for mol_id, run_id, prop_name, prop_value in chunk
        ]
        stmt = pg_insert(table).values(rows)
        stmt = stmt.on_conflict_do_update(
            index_elements=["molecule_id", "property_name", "run_id"],
            set_={"property_value": stmt.excluded.property_value},
        )
        async with AsyncSessionV9() as session:
            await session.execute(stmt)
            await session.commit()

    logger.info("Batch upserted %d properties (%d chunks of %d)",
                len(entries), (len(entries) + batch_size - 1) // batch_size, batch_size)


async def _cache_get(cache_key: str):
    """Look up a cached calculation result. Returns dict or None."""
    from database_v9 import AsyncSessionV9
    from models_v9 import CalculationCacheORM_V9
    from sqlalchemy import select

    async with AsyncSessionV9() as session:
        result = await session.execute(
            select(CalculationCacheORM_V9.result).where(
                CalculationCacheORM_V9.cache_key == cache_key
            )
        )
        row = result.scalar_one_or_none()
        return row if row else None


async def _cache_get_batch(cache_keys: list[str]) -> dict[str, dict]:
    """Batch lookup cached results. Returns {cache_key: result} for hits."""
    from database_v9 import AsyncSessionV9
    from models_v9 import CalculationCacheORM_V9
    from sqlalchemy import select

    if not cache_keys:
        return {}

    hits = {}
    # Query in chunks to avoid oversized IN clauses
    chunk_size = 500
    for i in range(0, len(cache_keys), chunk_size):
        chunk = cache_keys[i : i + chunk_size]
        async with AsyncSessionV9() as session:
            result = await session.execute(
                select(
                    CalculationCacheORM_V9.cache_key,
                    CalculationCacheORM_V9.result,
                ).where(CalculationCacheORM_V9.cache_key.in_(chunk))
            )
            for row in result.fetchall():
                hits[row[0]] = row[1]
    return hits


async def _cache_set(cache_key: str, result: dict):
    """Store a calculation result in cache."""
    from database_v9 import AsyncSessionV9
    from models_v9 import CalculationCacheORM_V9

    async with AsyncSessionV9() as session:
        entry = CalculationCacheORM_V9(cache_key=cache_key, result=result)
        session.add(entry)
        try:
            await session.commit()
        except Exception:
            await session.rollback()  # duplicate key — ignore


async def _cache_set_batch(entries: list[tuple[str, dict]]):
    """Batch store cache entries. entries: [(cache_key, result), ...]."""
    from database_v9 import AsyncSessionV9
    from models_v9 import CalculationCacheORM_V9
    from sqlalchemy.dialects.postgresql import insert as pg_insert

    if not entries:
        return

    table = CalculationCacheORM_V9.__table__
    chunk_size = 500
    for i in range(0, len(entries), chunk_size):
        chunk = entries[i : i + chunk_size]
        rows = [{"cache_key": k, "result": v} for k, v in chunk]
        stmt = pg_insert(table).values(rows)
        stmt = stmt.on_conflict_do_nothing(index_elements=["cache_key"])
        async with AsyncSessionV9() as session:
            await session.execute(stmt)
            await session.commit()


async def _get_reference_ligands(phase_id):
    """Fetch all reference ligands from phase→campaign (list of dicts)."""
    from database_v9 import AsyncSessionV9
    from models_v9 import PhaseORM_V9, CampaignORM_V9
    from sqlalchemy import select

    async with AsyncSessionV9() as session:
        result = await session.execute(
            select(PhaseORM_V9).where(PhaseORM_V9.id == phase_id)
        )
        phase = result.scalar_one_or_none()
        if not phase:
            return []

        result = await session.execute(
            select(CampaignORM_V9).where(CampaignORM_V9.id == phase.campaign_id)
        )
        campaign = result.scalar_one_or_none()
        if not campaign:
            return []

        return list(campaign.reference_ligands or [])


async def _get_docking_context(phase_id):
    """Fetch receptor URL + pocket center for docking, navigating phase→campaign→project."""
    from database_v9 import AsyncSessionV9
    from models_v9 import PhaseORM_V9, CampaignORM_V9, ProjectORM_V9
    from sqlalchemy import select
    from sqlalchemy.orm import selectinload

    async with AsyncSessionV9() as session:
        result = await session.execute(
            select(PhaseORM_V9).where(PhaseORM_V9.id == phase_id)
        )
        phase = result.scalar_one_or_none()
        if not phase:
            return None

        result = await session.execute(
            select(CampaignORM_V9).where(CampaignORM_V9.id == phase.campaign_id)
        )
        campaign = result.scalar_one_or_none()
        if not campaign:
            return None

        result = await session.execute(
            select(ProjectORM_V9).where(ProjectORM_V9.id == campaign.project_id)
        )
        project = result.scalar_one_or_none()
        if not project:
            return None

        # Get pocket center — from campaign pocket_config or project's selected pocket
        pocket_center = None
        pocket_config = campaign.pocket_config or {}
        if pocket_config.get("center"):
            pocket_center = tuple(pocket_config["center"])

        if not pocket_center:
            # Fall back to project's detected pockets
            pockets = project.pockets_detected or []
            tp = project.target_preview or {}
            selected_idx = tp.get("selected_pocket_index", 0)
            if pockets and selected_idx < len(pockets):
                c = pockets[selected_idx].get("center")
                if c:
                    pocket_center = tuple(c)
            elif pockets:
                c = pockets[0].get("center")
                if c:
                    pocket_center = tuple(c)

        # Get receptor PDB download URL
        pdb_url = None
        tp = project.target_preview or {}
        structure = tp.get("structure") or {}
        if structure.get("download_url"):
            pdb_url = structure["download_url"]
        elif project.target_pdb_id:
            pdb_url = f"https://files.rcsb.org/download/{project.target_pdb_id}.pdb"

        # Check for project-level prepared receptor files
        prep_report = project.receptor_prep_report or {}
        project_dir = Path(f"/data/projects/{str(project.id)}")
        prepared_pdbqt = None
        prepared_pdb = None

        if prep_report and not prep_report.get("error"):
            pdbqt_candidate = Path(prep_report.get("pdbqt_path", ""))
            pdb_candidate = Path(prep_report.get("prepared_pdb_path", ""))
            if pdbqt_candidate.exists() and pdbqt_candidate.stat().st_size > 100:
                prepared_pdbqt = pdbqt_candidate
            if pdb_candidate.exists() and pdb_candidate.stat().st_size > 100:
                prepared_pdb = pdb_candidate

        # Resolve UniProt ID from input type or target_preview
        uniprot_id = ""
        if project.target_input_type == "uniprot" and project.target_input_value:
            uniprot_id = project.target_input_value
        elif tp.get("uniprot", {}).get("uniprotId"):
            uniprot_id = tp["uniprot"]["uniprotId"]

        # Functional residues from target_preview (user-defined + auto-detected)
        functional_residues = tp.get("functional_residues", None)

        return {
            "pdb_url": pdb_url,
            "pocket_center": pocket_center,
            "pocket_config": pocket_config,
            "pdb_id": project.target_pdb_id,
            "structure_source": project.structure_source,
            "prepared_pdbqt": prepared_pdbqt,
            "prepared_pdb": prepared_pdb,
            "uniprot_id": uniprot_id,
            "functional_residues": functional_residues,
        }


async def _get_molecules_by_ids(molecule_ids):
    """Fetch molecules by their IDs. Returns list of dicts."""
    from database_v9 import AsyncSessionV9
    from models_v9 import MoleculeORM_V9
    from sqlalchemy import select

    async with AsyncSessionV9() as session:
        result = await session.execute(
            select(MoleculeORM_V9).where(MoleculeORM_V9.id.in_(molecule_ids))
        )
        mols = result.scalars().all()
        return [
            {
                "id": m.id,
                "smiles": m.smiles,
                "canonical_smiles": m.canonical_smiles,
                "name": m.name,
                "phase_id": m.phase_id,
                "generation_level": m.generation_level,
            }
            for m in mols
        ]


async def _get_all_phase_molecules(phase_id):
    """Fetch ALL molecules in a phase. Returns list of dicts."""
    from database_v9 import AsyncSessionV9
    from models_v9 import MoleculeORM_V9
    from sqlalchemy import select

    pid = uuid.UUID(phase_id) if isinstance(phase_id, str) else phase_id
    async with AsyncSessionV9() as session:
        result = await session.execute(
            select(MoleculeORM_V9).where(MoleculeORM_V9.phase_id == pid)
        )
        mols = result.scalars().all()
        return [
            {
                "id": m.id,
                "smiles": m.smiles,
                "canonical_smiles": m.canonical_smiles,
                "name": m.name,
                "phase_id": m.phase_id,
                "generation_level": m.generation_level,
            }
            for m in mols
        ]


async def _get_molecules_by_run(run_id: str, phase_id: str):
    """Fetch molecules created by a specific run. Returns list of dicts."""
    from database_v9 import AsyncSessionV9
    from models_v9 import MoleculeORM_V9
    from sqlalchemy import select

    async with AsyncSessionV9() as session:
        result = await session.execute(
            select(MoleculeORM_V9).where(
                MoleculeORM_V9.source_run_id == uuid.UUID(run_id),
                MoleculeORM_V9.phase_id == uuid.UUID(phase_id) if isinstance(phase_id, str) else MoleculeORM_V9.phase_id == phase_id,
            )
        )
        mols = result.scalars().all()
        return [
            {
                "id": m.id,
                "smiles": m.smiles,
                "canonical_smiles": m.canonical_smiles,
                "name": m.name,
                "phase_id": m.phase_id,
                "generation_level": m.generation_level,
            }
            for m in mols
        ]


async def _get_existing_properties(molecule_ids: list) -> dict:
    """Fetch all existing properties for a list of molecules.

    Returns {molecule_id: {flat_prop_key: value, ...}} where nested
    property dicts are flattened one level deep.
    """
    from database_v9 import AsyncSessionV9
    from models_v9 import MoleculePropertyORM_V9
    from sqlalchemy import select

    result_map: dict = {}
    if not molecule_ids:
        return result_map

    chunk_size = 500
    for i in range(0, len(molecule_ids), chunk_size):
        chunk = molecule_ids[i : i + chunk_size]
        async with AsyncSessionV9() as session:
            stmt = select(MoleculePropertyORM_V9).where(
                MoleculePropertyORM_V9.molecule_id.in_(chunk)
            )
            result = await session.execute(stmt)
            rows = result.scalars().all()
            for row in rows:
                mid = row.molecule_id
                if mid not in result_map:
                    result_map[mid] = {}
                val = row.property_value
                if isinstance(val, dict):
                    # Flatten one level (e.g. physicochemical → logP, MW, etc.)
                    for k, v in val.items():
                        if isinstance(v, dict):
                            for kk, vv in v.items():
                                result_map[mid][kk] = vv
                        else:
                            result_map[mid][k] = v
                else:
                    result_map[mid][row.property_name] = val
    return result_map


# ---------------------------------------------------------------------------
# Run log helper
# ---------------------------------------------------------------------------

async def _write_log_async(run_id: str, level: str, message: str):
    """Write a log entry for a run."""
    from database_v9 import AsyncSessionV9
    from models_v9 import RunLogORM_V9

    async with AsyncSessionV9() as session:
        log = RunLogORM_V9(
            run_id=uuid.UUID(run_id),
            level=level,
            message=message,
        )
        session.add(log)
        await session.commit()


def _write_log(run_id: str, level: str, message: str):
    """Write a run log entry (sync wrapper)."""
    try:
        _db_sync(_write_log_async(run_id, level, message))
    except Exception:
        logger.warning("Failed to write run log: %s", message)


# ---------------------------------------------------------------------------
# Preparation report helpers
# ---------------------------------------------------------------------------

def _aggregate_prep_stats(dock_results: list[dict], failure_metas: list[dict] | None = None) -> dict:
    """Aggregate preparation metadata from docking results + failures.

    Returns stats dict with counts for each step.
    """
    stats = {
        "total_input": 0,
        "passed": 0,
        "failed_parse": 0,
        "failed_protonation": 0,
        "failed_conformer": 0,
        "failed_rotatable_bonds": 0,
        "failed_other": 0,
        "tautomer_changed": 0,
        "standardized": 0,
        "protonated": 0,
        "minimization_mmff94": 0,
        "minimization_uff": 0,
        "minimization_unoptimized": 0,
    }

    # Count successes from dock_results
    for dr in dock_results:
        meta = dr.get("_prep_meta", {})
        if not meta:
            stats["passed"] += 1
            continue
        stats["passed"] += 1
        if meta.get("tautomer_changed"):
            stats["tautomer_changed"] += 1
        if meta.get("standardized"):
            stats["standardized"] += 1
        if meta.get("protonated"):
            stats["protonated"] += 1
        mini = meta.get("minimization", "")
        if mini == "MMFF94":
            stats["minimization_mmff94"] += 1
        elif mini == "UFF":
            stats["minimization_uff"] += 1
        elif mini == "unoptimized":
            stats["minimization_unoptimized"] += 1

    # Count failures
    for fm in (failure_metas or []):
        meta = fm.get("_prep_meta", fm)
        failure = meta.get("failure", "")
        if failure == "parse":
            stats["failed_parse"] += 1
        elif failure == "protonation":
            stats["failed_protonation"] += 1
        elif failure == "conformer":
            stats["failed_conformer"] += 1
        elif failure == "rotatable_bonds":
            stats["failed_rotatable_bonds"] += 1
        elif failure:
            stats["failed_other"] += 1

    stats["total_input"] = stats["passed"] + stats["failed_parse"] + stats["failed_protonation"] + \
        stats["failed_conformer"] + stats["failed_rotatable_bonds"] + stats["failed_other"]
    return stats


def _format_prep_summary(stats: dict) -> str:
    """Format preparation stats into a human-readable log message."""
    total = stats["total_input"]
    passed = stats["passed"]
    lines = [f"Preparation summary: {passed}/{total} passed"]

    failures = []
    if stats["failed_parse"]:
        failures.append(f"  - {stats['failed_parse']} failed parse")
    if stats["failed_rotatable_bonds"]:
        failures.append(f"  - {stats['failed_rotatable_bonds']} filtered (rotatable bonds > 15)")
    if stats["failed_protonation"]:
        failures.append(f"  - {stats['failed_protonation']} failed protonation")
    if stats["failed_conformer"]:
        failures.append(f"  - {stats['failed_conformer']} failed conformer generation")
    if stats["failed_other"]:
        failures.append(f"  - {stats['failed_other']} failed (other)")
    if failures:
        lines.extend(failures)

    details = []
    if stats["tautomer_changed"]:
        details.append(f"  - {stats['tautomer_changed']} had tautomer change")
    if stats["standardized"]:
        details.append(f"  - {stats['standardized']} were standardized (salts/fragments)")
    mini_parts = []
    if stats["minimization_mmff94"]:
        mini_parts.append(f"{stats['minimization_mmff94']} MMFF94")
    if stats["minimization_uff"]:
        mini_parts.append(f"{stats['minimization_uff']} UFF fallback")
    if stats["minimization_unoptimized"]:
        mini_parts.append(f"{stats['minimization_unoptimized']} unoptimized")
    if mini_parts:
        details.append(f"  - Minimization: {', '.join(mini_parts)}")
    if details:
        lines.extend(details)

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Receptor preparation Celery task (project-level, one-time)
# ---------------------------------------------------------------------------

@celery_app.task(name="tasks_v9.prepare_receptor_task", bind=True, max_retries=1)
def prepare_receptor_task(self, project_id: str, force: bool = False):
    """Prepare receptor at project level. Downloads PDB, runs 7-step pipeline,
    stores prepared files in /data/projects/{project_id}/."""
    import requests as http_requests

    logger.info("Starting receptor preparation for project %s (force=%s)", project_id, force)

    # Fetch project data
    async def _get_project_for_prep(pid):
        from database_v9 import AsyncSessionV9
        from models_v9 import ProjectORM_V9
        from sqlalchemy import select

        async with AsyncSessionV9() as session:
            result = await session.execute(
                select(ProjectORM_V9).where(ProjectORM_V9.id == uuid.UUID(pid))
            )
            p = result.scalar_one_or_none()
            if p:
                return {
                    "id": str(p.id),
                    "target_pdb_id": p.target_pdb_id,
                    "target_preview": p.target_preview,
                    "structure_source": p.structure_source,
                    "receptor_prep_config": p.receptor_prep_config,
                }
        return None

    async def _save_prep_result(pid, report):
        from database_v9 import AsyncSessionV9
        from models_v9 import ProjectORM_V9
        from sqlalchemy import select

        async with AsyncSessionV9() as session:
            result = await session.execute(
                select(ProjectORM_V9).where(ProjectORM_V9.id == uuid.UUID(pid))
            )
            p = result.scalar_one_or_none()
            if p:
                p.receptor_prep_report = report
                await session.commit()

    try:
        project = _db_sync(_get_project_for_prep(project_id))
        if not project:
            logger.error("Project %s not found", project_id)
            return {"error": "Project not found"}

        # Get PDB URL
        tp = project.get("target_preview") or {}
        structure = tp.get("structure") or {}
        pdb_url = structure.get("download_url")
        if not pdb_url and project.get("target_pdb_id"):
            pdb_url = f"https://files.rcsb.org/download/{project['target_pdb_id']}.pdb"
        if not pdb_url:
            return {"error": "No PDB URL configured"}

        # Create persistent project directory
        project_dir = Path(f"/data/projects/{project_id}")
        project_dir.mkdir(parents=True, exist_ok=True)

        pdb_filename = project.get("target_pdb_id") or "receptor"
        pdb_path = project_dir / f"{pdb_filename}.pdb"

        # Clear cached files when force re-preparing (before download)
        if force:
            for old_file in project_dir.glob(f"{pdb_filename}*"):
                logger.info("Force: removing cached file %s", old_file)
                old_file.unlink(missing_ok=True)

        # Download PDB
        logger.info("Downloading receptor PDB from %s", pdb_url)
        resp = http_requests.get(pdb_url, timeout=(10, 60))
        resp.raise_for_status()
        pdb_path.write_text(resp.text)

        # Get prep config (defaults = all checked)
        prep_config = project.get("receptor_prep_config") or {}
        structure_source = project.get("structure_source") or "pdb"

        # Run 7-step preparation pipeline
        from pipeline.prepare import prepare_receptor

        pdbqt_path, prep_report = prepare_receptor(
            pdb_path=pdb_path,
            work_dir=project_dir,
            structure_source=structure_source,
            keep_metals=prep_config.get("keep_metals", True),
            keep_cofactors=prep_config.get("keep_cofactors", True),
            minimize=prep_config.get("minimize"),
            ph=prep_config.get("ph", 7.4),
        )

        # Add file paths to report
        prep_report["pdbqt_path"] = str(pdbqt_path)
        prep_report["prepared_pdb_path"] = str(project_dir / f"{pdb_path.stem}_prepared.pdb")
        prep_report["project_dir"] = str(project_dir)

        # Save report to DB
        _db_sync(_save_prep_result(project_id, prep_report))

        logger.info("Receptor preparation complete for project %s: %s",
                     project_id, prep_report.get("conversion_method"))
        return {"status": "completed", "prep_report": prep_report}

    except Exception as exc:
        logger.error("Receptor preparation failed for project %s: %s", project_id, exc)
        error_report = {"error": str(exc), "steps_completed": [], "steps_skipped": []}
        _db_sync(_save_prep_result(project_id, error_report))
        return {"error": str(exc)}


# ---------------------------------------------------------------------------
# Main Celery task
# ---------------------------------------------------------------------------

@celery_app.task(name="tasks_v9.execute_run", bind=True, max_retries=1)
def execute_run(self, run_id: str):
    """Execute a V9 run. Dispatches based on run type."""
    logger.info("Starting V9 run %s", run_id)

    run_data = _db_sync(_get_run(run_id))
    if not run_data:
        logger.error("Run %s not found", run_id)
        return {"error": "Run not found"}

    if run_data["status"] == "cancelled":
        logger.info("Run %s was cancelled, skipping", run_id)
        return {"status": "cancelled"}

    # Mark as running
    _db_sync(_update_run(
        run_id,
        status="running",
        started_at=datetime.now(timezone.utc),
        progress=0,
        current_step="Starting...",
    ))
    _write_log(run_id, "info", f"Run started — type: {run_data['type']}")

    try:
        run_type = run_data["type"]
        if run_type == "import":
            result = _run_import(run_id, run_data)
        elif run_type == "calculation":
            result = _run_calculation(run_id, run_data)
        elif run_type == "generation":
            result = _run_generation(run_id, run_data)
        elif run_type == "afvs":
            result = _run_afvs(run_id, run_data)
        else:
            raise ValueError(f"Unknown run type: {run_type}")

        # Re-check status before marking completed (user may have cancelled)
        fresh = _db_sync(_get_run(run_id))
        if fresh and fresh["status"] == "cancelled":
            logger.info("Run %s was cancelled during execution, not overwriting", run_id)
            _write_log(run_id, "info", "Run was cancelled during execution")
            return {"status": "cancelled"}

        _db_sync(_update_run(
            run_id,
            status="completed",
            progress=100,
            current_step=None,
            completed_at=datetime.now(timezone.utc),
        ))
        _write_log(run_id, "success", f"Run completed — {result}")
        logger.info("Run %s completed successfully", run_id)
        return result

    except Exception as e:
        logger.exception("Run %s failed: %s", run_id, e)
        _write_log(run_id, "error", f"Run failed: {str(e)[:500]}")
        _db_sync(_update_run(
            run_id,
            status="failed",
            error_message=str(e)[:2000],
            completed_at=datetime.now(timezone.utc),
        ))
        return {"error": str(e)}


# ---------------------------------------------------------------------------
# Import run — streaming & batch helpers
# ---------------------------------------------------------------------------

IMPORT_BATCH_SIZE = 500


def _canonicalize_stream(raw_iter):
    """Generator: canonicalize + validate SMILES in a single RDKit pass.

    Accepts str (raw SMILES) or dict ({"smiles", "name", ...}).
    Yields {"smiles": str, "canonical_smiles": str, "name": str|None}.
    """
    try:
        from rdkit import Chem
        has_rdkit = True
    except ImportError:
        has_rdkit = False

    counter = 0
    for item in raw_iter:
        if isinstance(item, dict):
            smi = (item.get("smiles") or "").strip()
            name = item.get("name") or item.get("chembl_id") or item.get("zinc_id")
        else:
            smi = str(item).strip()
            name = None

        if not smi:
            continue

        # Shortcut: if source already provides canonical SMILES, skip RDKit re-parse
        pre_canonical = item.get("canonical_smiles") if isinstance(item, dict) else None
        if pre_canonical:
            canonical = pre_canonical
        elif has_rdkit:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            canonical = Chem.MolToSmiles(mol)
        else:
            canonical = smi

        counter += 1
        if not name:
            name = f"Mol_{counter}"

        result = {"smiles": smi, "canonical_smiles": canonical, "name": name}
        # Preserve source metadata (mwt, logp, hbd, hba, tpsa, pchembl, etc.)
        if isinstance(item, dict):
            for meta_key in ("mwt", "logp", "import_hbd", "import_hba", "import_tpsa", "pchembl_value", "activity_value_nM", "activity_type", "assay_name", "source", "mw"):
                if item.get(meta_key) is not None:
                    result[meta_key] = item[meta_key]
        yield result


def _stream_smiles_file(file_path: Path, file_format: str):
    """Generator: stream SMILES from a file line-by-line (no full load in memory).

    Supports smiles/csv (text line) and sdf (RDKit ForwardSDMolSupplier).
    """
    if file_format in ("smiles", "csv"):
        with open(file_path, "r", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                smi = line.split(",")[0].split("\t")[0].strip()
                if smi:
                    yield smi
    elif file_format == "sdf":
        try:
            from rdkit import Chem
            suppl = Chem.ForwardSDMolSupplier(open(str(file_path), 'rb'))
            for mol in suppl:
                if mol is not None:
                    yield Chem.MolToSmiles(mol)
        except ImportError:
            raise ImportError("RDKit required for SDF file parsing")


def _count_lines_fast(file_path: Path) -> int:
    """Fast line count for progress estimation."""
    count = 0
    with open(file_path, "rb") as f:
        for _ in f:
            count += 1
    return count


_IMPORT_META_KEYS = ("mwt", "logp", "import_hbd", "import_hba", "import_tpsa", "pchembl_value", "activity_value_nM", "activity_type", "assay_name", "mw", "source")


def _store_batch_metadata(batch: list[dict], created_rows: list, source_run_id):
    """Store source metadata as molecule properties for newly created molecules."""
    if not created_rows:
        return
    # Build lookup: canonical_smiles → metadata dict
    meta_by_smiles: dict[str, dict] = {}
    for mol in batch:
        can = mol.get("canonical_smiles", "")
        meta = {k: mol[k] for k in _IMPORT_META_KEYS if mol.get(k) is not None}
        if meta and can:
            meta_by_smiles[can] = meta

    if not meta_by_smiles:
        return

    entries = []
    for row in created_rows:
        mol_id, can_smi = row[0], row[1]
        meta = meta_by_smiles.get(can_smi)
        if meta:
            for prop_name, prop_value in meta.items():
                entries.append((mol_id, source_run_id, prop_name, prop_value))

    if entries:
        _db_sync(_store_properties_batch(entries))


def _import_molecules_chunked(run_id: str, phase_id, source_run_id, mol_iter, total_estimate: int = 0) -> dict:
    """Consume a molecule iterator in batches, insert via _create_molecules_batch.

    Updates progress once per batch (not per molecule).
    Stores source metadata (mwt, logp, etc.) as molecule properties.
    Returns {"created": N, "skipped": N}.
    """
    created_total = 0
    skipped_total = 0
    processed = 0
    batch = []

    for mol in mol_iter:
        batch.append(mol)
        if len(batch) >= IMPORT_BATCH_SIZE:
            c, s, rows = _db_sync(_create_molecules_batch(phase_id, batch, source_run_id))
            created_total += c
            skipped_total += s
            _store_batch_metadata(batch, rows, source_run_id)
            processed += len(batch)
            batch = []

            # Progress: 20% → 95% proportional to processed/total
            if total_estimate > 0:
                pct = 20 + int(75 * min(processed, total_estimate) / total_estimate)
            else:
                pct = min(20 + processed // 100, 95)
            _db_sync(_update_run(
                run_id,
                progress=pct,
                current_step=f"Importing... {processed:,} processed ({created_total:,} new)",
            ))

    # Flush remaining
    if batch:
        c, s, rows = _db_sync(_create_molecules_batch(phase_id, batch, source_run_id))
        created_total += c
        skipped_total += s
        _store_batch_metadata(batch, rows, source_run_id)

    return {"created": created_total, "skipped": skipped_total}


async def _fetch_bookmarked_smiles(src_phase_id: str):
    """Fetch bookmarked molecules from a source phase."""
    from database_v9 import AsyncSessionV9
    from models_v9 import MoleculeORM_V9
    from sqlalchemy import select as sa_select

    async with AsyncSessionV9() as session:
        result = await session.execute(
            sa_select(
                MoleculeORM_V9.id,
                MoleculeORM_V9.smiles,
                MoleculeORM_V9.canonical_smiles,
                MoleculeORM_V9.name,
            ).where(
                MoleculeORM_V9.phase_id == uuid.UUID(src_phase_id),
                MoleculeORM_V9.bookmarked == True,
            )
        )
        return [{"smiles": row.smiles, "name": row.name} for row in result.all()]


async def _copy_phase_properties(source_phase_id: str, target_phase_id: str, import_run_id: uuid.UUID) -> int:
    """Copy all calculated properties from bookmarked source molecules to target phase molecules.

    Matches molecules by canonical_smiles. For each property, takes the most recent value.
    Uses ON CONFLICT DO NOTHING so re-imports are safe.
    Returns count of properties copied.
    """
    from database_v9 import AsyncSessionV9
    from sqlalchemy import text

    sql = text("""
        INSERT INTO molecule_properties (id, molecule_id, run_id, property_name, property_value, created_at)
        SELECT gen_random_uuid(), t.id, :import_run_id, sub.property_name, sub.property_value, now()
        FROM molecules s
        JOIN molecules t
          ON t.canonical_smiles = s.canonical_smiles
          AND t.phase_id = :target_phase_id
        JOIN LATERAL (
            SELECT DISTINCT ON (property_name) property_name, property_value
            FROM molecule_properties WHERE molecule_id = s.id
            ORDER BY property_name, created_at DESC
        ) sub ON true
        WHERE s.phase_id = :source_phase_id AND s.bookmarked = true
        ON CONFLICT (molecule_id, property_name, run_id) DO NOTHING
    """)

    async with AsyncSessionV9() as session:
        result = await session.execute(sql, {
            "import_run_id": import_run_id,
            "source_phase_id": uuid.UUID(str(source_phase_id)),
            "target_phase_id": uuid.UUID(str(target_phase_id)),
        })
        await session.commit()
        return result.rowcount


def _fetch_database_ligands(config: dict):
    """Fetch ligands from a public database.

    Supports single source (new) and multi-source (legacy).
    Runs in Celery worker context. Returns list of dicts.
    """
    databases = config.get("databases", [])
    uniprot_id = config.get("uniprot_id", "")
    max_per_source = config.get("max_compounds", config.get("max_per_source", 50))
    filters = config.get("filters", {})
    all_ligands = []

    for db_name in databases:
        try:
            if db_name == "chembl":
                if not uniprot_id:
                    raise ValueError(
                        "ChEMBL import requires a UniProt ID. "
                        "Please configure a target with a UniProt accession in Target Setup."
                    )
                from pipeline.ligands import fetch_chembl_ligands
                all_ligands.extend(fetch_chembl_ligands(
                    uniprot_id, max_count=max_per_source, filters=filters if filters else None,
                ))

            elif db_name == "pubchem":
                if not uniprot_id:
                    raise ValueError(
                        "PubChem import requires a UniProt ID. "
                        "Please configure a target with a UniProt accession in Target Setup."
                    )
                from pipeline.ligands import fetch_pubchem_ligands, resolve_gene_name
                gene = resolve_gene_name(uniprot_id)
                if gene:
                    all_ligands.extend(fetch_pubchem_ligands(gene, max_count=max_per_source))
                else:
                    raise ValueError(
                        f"Could not resolve gene name for UniProt ID '{uniprot_id}'. "
                        "PubChem requires a valid gene name to search."
                    )

            elif db_name == "zinc":
                from pipeline.ligands import fetch_zinc_ligands
                subset = config.get("zinc_subset")
                all_ligands.extend(fetch_zinc_ligands(
                    max_count=max_per_source,
                    subset=subset,
                    filters=filters if filters else None,
                ))

            elif db_name == "enamine":
                from pipeline.ligands import sample_enamine_real
                for smi in sample_enamine_real(n=max_per_source):
                    if isinstance(smi, dict):
                        all_ligands.append(smi)
                    else:
                        all_ligands.append({"smiles": smi, "name": None, "source": "enamine_real"})

            elif db_name == "fragments":
                from pipeline.ligands import load_fragment_library
                all_ligands.extend(load_fragment_library())

        except Exception as exc:
            logger.warning("Failed to fetch from %s: %s", db_name, exc)
            raise  # Propagate to mark run as failed with clear error

    return all_ligands


# ---------------------------------------------------------------------------
# Import run — main entry
# ---------------------------------------------------------------------------

def _run_import(run_id: str, run_data: dict) -> dict:
    """Import molecules from various sources using streaming + batch inserts."""
    config = run_data["config"] or {}
    phase_id = run_data["phase_id"]
    run_uuid = uuid.UUID(run_id)

    # ── Branch 1: Database import (fetch in worker) ──────────────────────
    if config.get("source") == "database":
        db_name = config.get("database", config.get("databases", ["unknown"])[0])
        max_compounds = config.get("max_compounds", config.get("max_per_source", 50))
        _db_sync(_update_run(run_id, current_step=f"Fetching from {db_name}...", progress=5))
        _write_log(run_id, "info", f"Fetching up to {max_compounds} compounds from {db_name}")

        ligands = _fetch_database_ligands(config)
        if not ligands:
            raise ValueError("No compounds found from selected databases")

        _write_log(run_id, "info", f"Fetched {len(ligands)} compounds — importing...")
        mol_iter = _canonicalize_stream(iter(ligands))
        result = _import_molecules_chunked(run_id, phase_id, run_uuid, mol_iter, len(ligands))

    # ── Branch 2: File upload (streaming) ────────────────────────────────
    elif run_data.get("input_source") == "file" and run_data.get("input_file_path"):
        file_path = Path(run_data["input_file_path"])
        if not file_path.exists():
            raise FileNotFoundError(f"Upload file not found: {file_path}")

        file_format = config.get("file_format", "smiles")
        _db_sync(_update_run(run_id, current_step="Counting molecules...", progress=5))

        total_estimate = _count_lines_fast(file_path) if file_format != "sdf" else 0
        _write_log(run_id, "info", f"File: {config.get('original_filename', file_path.name)} (~{total_estimate} lines)")

        raw_iter = _stream_smiles_file(file_path, file_format)
        mol_iter = _canonicalize_stream(raw_iter)
        result = _import_molecules_chunked(run_id, phase_id, run_uuid, mol_iter, total_estimate)

    # ── Branch 3: Phase selection (bookmarked) ───────────────────────────
    elif config.get("source") == "phase_selection":
        source_phase_id = config.get("source_phase_id")
        if not source_phase_id:
            raise ValueError("source_phase_id required for phase_selection import")

        _db_sync(_update_run(run_id, current_step="Fetching bookmarked molecules...", progress=5))
        bookmarked = _db_sync(_fetch_bookmarked_smiles(source_phase_id))
        logger.info("Phase selection: %d bookmarked from phase %s", len(bookmarked), source_phase_id)

        if not bookmarked:
            raise ValueError("No bookmarked molecules in source phase")

        mol_iter = _canonicalize_stream(iter(bookmarked))
        result = _import_molecules_chunked(run_id, phase_id, run_uuid, mol_iter, len(bookmarked))

        # Copy all calculated properties from source phase to target phase
        if result["created"] > 0:
            _db_sync(_update_run(run_id, current_step="Copying properties from previous phase...", progress=92))
            copied = _db_sync(_copy_phase_properties(source_phase_id, phase_id, run_uuid))
            _write_log(run_id, "info", f"Copied {copied} properties from source phase")

    # ── Branch 4: Inline smiles_list (manual / small) ────────────────────
    elif config.get("smiles_list"):
        smiles_list = config["smiles_list"]
        _write_log(run_id, "info", f"Inline import: {len(smiles_list)} SMILES")
        _db_sync(_update_run(run_id, current_step="Validating SMILES...", progress=10))

        mol_iter = _canonicalize_stream(iter(smiles_list))
        result = _import_molecules_chunked(run_id, phase_id, run_uuid, mol_iter, len(smiles_list))

    # ── Branch 5: Legacy ligand_list in config (backward compat) ─────────
    elif config.get("ligand_list"):
        ligand_list = config["ligand_list"]
        _write_log(run_id, "info", f"Legacy import: {len(ligand_list)} ligands from config")
        _db_sync(_update_run(run_id, current_step="Importing ligands...", progress=10))

        mol_iter = _canonicalize_stream(iter(ligand_list))
        result = _import_molecules_chunked(run_id, phase_id, run_uuid, mol_iter, len(ligand_list))

    else:
        raise ValueError("No SMILES source found in run config")

    _write_log(run_id, "success", f"Import done: {result['created']} created, {result['skipped']} duplicates skipped")
    logger.info("Import run %s: %d created, %d skipped", run_id, result["created"], result["skipped"])
    return result


# ---------------------------------------------------------------------------
# Calculation run
# ---------------------------------------------------------------------------

# Frontend column key → backend property key(s) mapping.
# The frontend uses display keys (e.g. "BBB", "HBD") while the backend pipeline
# stores different keys (e.g. "bbb_permeability", "hbd"). This map allows
# _filter_cols to match backend keys against frontend included_columns.
_FRONTEND_TO_BACKEND = {
    "HBD": "hbd", "HBA": "hba", "QED": "qed", "TPSA": "tpsa",
    "BBB": "bbb_permeability", "hERG": "herg_inhibition",
    "safety_color_code": "color_code",
    "docking_score": "affinity",
}
# Build reverse: backend key → frontend key
_BACKEND_TO_FRONTEND = {v: k for k, v in _FRONTEND_TO_BACKEND.items()}


def _expand_included(included: set | None) -> set | None:
    """Expand frontend column keys to include their backend aliases."""
    if included is None:
        return None
    expanded = set(included)
    for frontend_key in list(included):
        backend_key = _FRONTEND_TO_BACKEND.get(frontend_key)
        if backend_key:
            expanded.add(backend_key)
    return expanded


_ALWAYS_KEEP = {
    "interactions_method", "interactions_detail", "mapped_functional",
    "pose_molblock", "method", "summary", "herg_specialized",
    "cluster_size",
    "scaffold_group", "scaffold_group_size",
    "pains_alert",
    "matched_types", "missing_types", "feature_counts", "ref_feature_counts",
    "pharmacophore_fit_1", "pharmacophore_fit_2", "pharmacophore_fit_3",
    "pharmacophore_fit_4", "pharmacophore_fit_5", "pharmacophore_fit_best",
    "pharmacophore_fit_avg",
}


def _filter_cols(props: dict, included: set | None) -> dict:
    """Filter a properties dict keeping only keys whose frontend column is included.

    Handles nested dicts (e.g. ADMET absorption/distribution/toxicity sub-dicts)
    by recursively filtering. Pass-through if included is None.
    Keys in _ALWAYS_KEEP are never filtered out (metadata fields).
    """
    if included is None:
        return props
    result = {}
    for k, v in props.items():
        # Always keep metadata keys
        if k in _ALWAYS_KEEP:
            result[k] = v
        # Check if this key (or its frontend alias) is in the included set
        elif k in included or _BACKEND_TO_FRONTEND.get(k) in included:
            result[k] = v
        elif isinstance(v, dict) and k not in included:
            # Recurse into nested dicts (e.g. admet.absorption, admet.toxicity)
            filtered_sub = _filter_cols(v, included)
            if filtered_sub:
                result[k] = filtered_sub
    return result


def _run_calculation(run_id: str, run_data: dict) -> dict:
    """Run calculations on selected molecules."""
    calc_types = run_data.get("calculation_types") or []
    input_ids = run_data.get("input_molecule_ids") or []
    run_uuid = uuid.UUID(run_id)
    # Column filtering from RunCreator checklist
    _run_config = run_data.get("config") or {}
    _inc_raw = set(_run_config["included_columns"]) if _run_config.get("included_columns") else None
    _inc = _expand_included(_inc_raw)

    # Run on all molecules in phase (ignore selection)
    if _run_config.get("run_all_molecules"):
        phase_id = run_data.get("phase_id")
        if not phase_id:
            raise ValueError("Phase ID required for run_all_molecules")
        molecules = _db_sync(_get_all_phase_molecules(phase_id))
    elif input_ids:
        molecules = _db_sync(_get_molecules_by_ids(input_ids))
    else:
        raise ValueError("No input molecules selected")

    if not molecules:
        raise ValueError("No valid molecules found")

    smiles_list = [m["smiles"] for m in molecules]
    mol_ids = [m["id"] for m in molecules]
    total_steps = len(calc_types)
    results = {}

    def _cached_compute(calc_type, smiles, compute_fn):
        """Try cache first, then compute and cache the result."""
        cache_key = f"{calc_type}:v2:{smiles}"
        cached = _db_sync(_cache_get(cache_key))
        if cached is not None:
            return cached
        result = compute_fn(smiles)
        _db_sync(_cache_set(cache_key, result))
        return result

    for step_idx, calc_type in enumerate(calc_types):
        # Check if user cancelled between steps
        _check = _db_sync(_get_run(run_id))
        if _check and _check["status"] == "cancelled":
            _write_log(run_id, "info", f"Run cancelled before {calc_type}")
            return {"calculations": results, "cancelled": True}

        step_label = f"Running {calc_type} ({step_idx+1}/{total_steps})"
        base_pct = int(100 * step_idx / total_steps) if total_steps > 0 else 0
        _db_sync(_update_run(run_id, current_step=step_label, progress=base_pct))

        try:
            _write_log(run_id, "info", f"Starting {calc_type} on {len(molecules)} molecules")

            if calc_type in ("adme", "toxicity", "admet"):
                from pipeline.admet import predict_admet
                from pipeline.scoring import (
                    compute_properties, compute_cns_mpo,
                    compute_druglikeness_rules, compute_brenk_alert,
                )
                # predict_admet() returns full ADME+tox in one call
                admet_results = predict_admet(smiles_list)

                # Keys that belong to each subset
                _ADME_SECTIONS = {"absorption", "distribution", "metabolism", "excretion"}
                _TOX_KEYS = {"toxicity", "composite_score", "flags", "color_code"}

                adme_entries = []
                tox_entries = []
                druglike_entries = []

                for mol, admet in zip(molecules, admet_results):
                    # ── ADME subset ──
                    if calc_type in ("adme", "admet"):
                        adme_props = {k: v for k, v in admet.items() if k in _ADME_SECTIONS}
                        adme_filtered = _filter_cols(adme_props, _inc)
                        # Always preserve metabolism (CYP inhibitions) — detail-level data
                        # not in column checklist but shown in ADME detail panel
                        if "metabolism" in adme_props and "metabolism" not in adme_filtered:
                            adme_filtered["metabolism"] = adme_props["metabolism"]
                        adme_entries.append((mol["id"], run_uuid, "adme", adme_filtered))
                        # CNS MPO + Pfizer/GSK + Brenk (pharmacokinetics rules)
                        _piggyback_keys = {"cns_mpo", "pfizer_alert", "gsk_alert", "brenk_alert"}
                        if _inc is None or _piggyback_keys & _inc:
                            try:
                                props = compute_properties(mol["smiles"])
                                cns_mpo = compute_cns_mpo(props)
                                dl_rules = compute_druglikeness_rules(props)
                                brenk = compute_brenk_alert(mol["smiles"])
                                extra = {
                                    "cns_mpo": cns_mpo,
                                    "pfizer_alert": dl_rules["pfizer_alert"],
                                    "gsk_alert": dl_rules["gsk_alert"],
                                    "brenk_alert": brenk,
                                }
                                druglike_entries.append((mol["id"], run_uuid, "druglikeness_rules", _filter_cols(extra, _inc)))
                            except Exception as e:
                                logger.debug("CNS MPO/Brenk failed for %s: %s", mol["smiles"][:40], e)

                    # ── Toxicity subset ──
                    if calc_type in ("toxicity", "admet"):
                        tox_props = {k: v for k, v in admet.items() if k in _TOX_KEYS}
                        # Enrich with hERG specialized (IC50 + risk level) — toxicity only
                        if calc_type == "toxicity":
                            try:
                                from pipeline.admet import predict_herg_specialized
                                herg_detail = predict_herg_specialized(mol["smiles"])
                                tox_props["herg_specialized"] = {
                                    "ic50_um": herg_detail["ic50_um"],
                                    "risk_level": herg_detail["risk_level"],
                                    "method": herg_detail["method"],
                                }
                            except Exception as e:
                                logger.debug("hERG specialized failed for %s: %s", mol["smiles"][:40], e)
                            # PAINS alert (480 patterns — structural interference)
                            try:
                                from pipeline.scoring import compute_pains_alert
                                tox_props["pains_alert"] = compute_pains_alert(mol["smiles"])
                            except Exception as e:
                                logger.debug("PAINS check failed for %s: %s", mol["smiles"][:40], e)
                        tox_entries.append((mol["id"], run_uuid, "toxicity", _filter_cols(tox_props, _inc)))

                if adme_entries:
                    _db_sync(_store_properties_batch(adme_entries))
                if tox_entries:
                    _db_sync(_store_properties_batch(tox_entries))
                if druglike_entries:
                    _db_sync(_store_properties_batch(druglike_entries))

                label = {"adme": "ADME", "toxicity": "Toxicity", "admet": "ADMET (ADME + Toxicity)"}[calc_type]
                _write_log(run_id, "info", f"{label} complete: {len(admet_results)} molecules profiled")
                results[calc_type] = f"{len(admet_results)} molecules processed"

            elif calc_type == "scoring":
                from pipeline.scoring import (
                    compute_properties, compute_sa_score,
                )
                from pipeline.pharmacophore import _count_features
                from rdkit import Chem
                # Batch cache lookup (1 query instead of N)
                cache_keys = [f"scoring:{m['smiles']}" for m in molecules]
                cached = _db_sync(_cache_get_batch(cache_keys))

                physchem_entries = []
                sa_entries = []
                new_cache = []  # (key, result) pairs to batch-store
                for mol, ck in zip(molecules, cache_keys):
                    try:
                        props = cached.get(ck)
                        if props is None:
                            props = compute_properties(mol["smiles"])
                            new_cache.append((ck, props))
                        # Add pharmacophore feature counts into physicochemical props
                        try:
                            rdmol = Chem.MolFromSmiles(mol["smiles"])
                            if rdmol is not None:
                                fc = _count_features(rdmol)
                                props["n_hydrophobic"] = fc.get("hydrophobic", 0)
                                props["n_aromatic"] = fc.get("aromatic", 0)
                                props["n_pos_ionizable"] = fc.get("positive", 0)
                                props["n_neg_ionizable"] = fc.get("negative", 0)
                        except Exception as e_feat:
                            logger.warning("Feature count failed for %s: %s", mol["smiles"], e_feat)
                        physchem_entries.append((mol["id"], run_uuid, "physicochemical", _filter_cols(props, _inc)))
                        if _inc is None or "sa_score" in _inc:
                            sa = compute_sa_score(mol["smiles"])
                            if sa is not None:
                                sa_entries.append((mol["id"], run_uuid, "sa_score", {"sa_score": sa}))
                    except Exception as e:
                        logger.warning("Scoring failed for %s: %s", mol["smiles"], e)
                # Batch store: cache + properties
                if new_cache:
                    _db_sync(_cache_set_batch(new_cache))
                _db_sync(_store_properties_batch(physchem_entries))
                if sa_entries:
                    _db_sync(_store_properties_batch(sa_entries))
                _write_log(run_id, "info", f"Scoring complete: {len(physchem_entries)}/{len(molecules)} molecules scored")
                results["scoring"] = f"{len(molecules)} molecules scored"

            elif calc_type == "docking":
                import tempfile
                import requests as http_requests

                # Fetch receptor + pocket context from campaign/project
                phase_id = molecules[0]["phase_id"] if molecules else None
                docking_ctx = _db_sync(_get_docking_context(phase_id)) if phase_id else None

                if not docking_ctx or not docking_ctx.get("pdb_url"):
                    raise ValueError(
                        "Docking requires a target with a PDB structure. "
                        "Please configure the target in Target Setup first."
                    )
                if not docking_ctx.get("pocket_center"):
                    raise ValueError(
                        "No binding pocket center found. "
                        "Please detect pockets in Target Setup first."
                    )

                pdb_url = docking_ctx["pdb_url"]
                pocket_center = docking_ctx["pocket_center"]
                pocket_size = tuple(docking_ctx["pocket_config"].get("size", [22, 22, 22]))

                # Check for project-level prepared receptor (avoids re-preparing)
                prepared_pdbqt = docking_ctx.get("prepared_pdbqt")
                prepared_pdb = docking_ctx.get("prepared_pdb")

                if prepared_pdbqt:
                    _write_log(run_id, "info", "Using project-level prepared receptor (skipping re-preparation)")
                    _db_sync(_update_run(run_id, current_step="Using prepared receptor"))
                    receptor_pdbqt = prepared_pdbqt
                    # Use the prepared PDB's parent dir as work_dir for docking temp files
                    tmp_context = tempfile.TemporaryDirectory(prefix="docking_v9_")
                    tmp_path = Path(tmp_context.__enter__())
                else:
                    raise ValueError(
                        "Receptor not prepared at project level. "
                        "Please prepare the receptor in Target Setup before running docking."
                    )

                try:
                    _write_log(run_id, "info",
                        f"Docking {len(molecules)} molecules — center={pocket_center}, size={pocket_size}")

                    # Prepare ligand list — use molecule UUID as name for collision-free mapping
                    ligands = [{"name": str(m["id"]), "smiles": m["smiles"]}
                               for m in molecules]
                    # Build id→molecule mapping BEFORE dock_all_ligands so
                    # on_batch_done callback can resolve names
                    name_to_mol = {str(m["id"]): m for m in molecules}

                    # Read advanced docking config from run config
                    run_config = run_data.get("config") or {}
                    dock_cfg = run_config.get("docking", {}) if isinstance(run_config.get("docking"), dict) else run_config
                    dock_exhaustiveness = int(dock_cfg.get("exhaustiveness", 8))
                    dock_num_modes = int(dock_cfg.get("num_modes", 9))
                    dock_engine = dock_cfg.get("engine") or "gnina_gpu"
                    dock_num_cpu = int(dock_cfg.get("num_cpu") or 0)
                    dock_cnn_scoring = dock_cfg.get("cnn_scoring") or "rescore"

                    # Progress callback
                    def progress_cb(pct, msg):
                        _db_sync(_update_run(run_id, current_step=msg, progress=base_pct + int(pct / total_steps)))

                    # Partial results callback — save each batch immediately
                    def on_batch_done(batch_results):
                        partial_entries = []
                        for dr in batch_results:
                            mol = name_to_mol.get(dr.get("name"))
                            if not mol:
                                continue
                            docking_props = {
                                "affinity": dr.get("affinity"),
                                "vina_score": dr.get("vina_score"),
                                "docking_engine": dr.get("docking_engine"),
                            }
                            if "cnn_score" in dr:
                                docking_props["cnn_score"] = dr["cnn_score"]
                                docking_props["cnn_affinity"] = dr.get("cnn_affinity")
                                docking_props["cnn_vs"] = dr.get("cnn_vs")
                            if dr.get("pocket_distance") is not None:
                                docking_props["pocket_distance"] = dr["pocket_distance"]
                            if dr.get("pose_molblock"):
                                docking_props["pose_molblock"] = dr["pose_molblock"][:50000]
                            filtered = _filter_cols(docking_props, _inc)
                            for always_keep in ("pose_molblock", "docking_engine", "pocket_distance"):
                                if always_keep in docking_props:
                                    filtered[always_keep] = docking_props[always_keep]
                            partial_entries.append((mol["id"], run_uuid, "docking", filtered))
                        if partial_entries:
                            _db_sync(_store_properties_batch(partial_entries))
                            _write_log(run_id, "info",
                                f"Partial results saved: {len(partial_entries)} molecules from batch")

                    from pipeline.docking import dock_all_ligands
                    dock_results = dock_all_ligands(
                        receptor_pdbqt=receptor_pdbqt,
                        ligands=ligands,
                        center=pocket_center,
                        work_dir=tmp_path,
                        progress_callback=progress_cb,
                        size=pocket_size,
                        exhaustiveness=dock_exhaustiveness,
                        docking_engine=dock_engine,
                        num_cpu=dock_num_cpu,
                        cnn_scoring=dock_cnn_scoring,
                        batch_callback=on_batch_done,
                    )

                    # Store results as molecule properties (batch insert)
                    dock_entries = []  # (mol_id, run_id, prop_name, prop_value)
                    for dr in dock_results:
                        mol = name_to_mol.get(dr.get("name"))
                        if not mol:
                            continue
                        docking_props = {
                            "affinity": dr.get("affinity"),
                            "vina_score": dr.get("vina_score"),
                            "consensus_rank": dr.get("consensus_rank"),
                            "consensus_ecr": dr.get("consensus_ecr"),
                            "docking_engine": dr.get("docking_engine"),
                        }
                        # CNN metrics only present for GNINA engine
                        if "cnn_score" in dr:
                            docking_props["cnn_score"] = dr["cnn_score"]
                            docking_props["cnn_affinity"] = dr.get("cnn_affinity")
                            docking_props["cnn_vs"] = dr.get("cnn_vs")
                        # Read pose 3D coordinates:
                        # 1. From file path (CPU docking — GNINA/Vina)
                        # 2. From inline molblock (GPU docking — RunPod)
                        pose_path = dr.get("pose_pdbqt_path")
                        if pose_path:
                            try:
                                pose_file = Path(pose_path)
                                if pose_file.exists():
                                    docking_props["pose_molblock"] = pose_file.read_text()[:50000]
                            except Exception as pe:
                                logger.debug("Could not read pose file for %s: %s", dr.get("name"), pe)
                        elif dr.get("pose_molblock"):
                            docking_props["pose_molblock"] = dr["pose_molblock"][:50000]
                        # Pocket distance (computed in docking_gpu from pose centroid)
                        if dr.get("pocket_distance") is not None:
                            docking_props["pocket_distance"] = dr["pocket_distance"]
                        # Filter but always keep pose_molblock + docking_engine
                        filtered_dock = _filter_cols(docking_props, _inc)
                        for always_keep in ("pose_molblock", "docking_engine", "pocket_distance"):
                            if always_keep in docking_props:
                                filtered_dock[always_keep] = docking_props[always_keep]
                        dock_entries.append((mol["id"], run_uuid, "docking", filtered_dock))

                    _db_sync(_store_properties_batch(dock_entries))
                    stored = len(dock_entries)

                    # --- Preparation report: aggregate stats + structured log ---
                    prep_stats = _aggregate_prep_stats(dock_results)
                    prep_summary = _format_prep_summary(prep_stats)
                    _write_log(run_id, "info", prep_summary)

                    # Store preparation metadata per molecule
                    prep_entries = []
                    for dr in dock_results:
                        meta = dr.get("_prep_meta")
                        if not meta:
                            continue
                        mol = name_to_mol.get(dr.get("name"))
                        if not mol:
                            continue
                        # Store clean prep metadata (no internal keys)
                        prep_data = {
                            k: v for k, v in meta.items()
                            if k not in ("failure",) and v is not None
                        }
                        if prep_data:
                            prep_entries.append((mol["id"], run_uuid, "preparation", prep_data))
                    if prep_entries:
                        _db_sync(_store_properties_batch(prep_entries))

                    # Ligand Efficiency (computed after docking scores are available)
                    if _inc is None or "ligand_efficiency" in _inc:
                        from pipeline.scoring import compute_ligand_efficiency, compute_properties
                        le_entries = []
                        for dr in dock_results:
                            mol = name_to_mol.get(dr.get("name"))
                            if not mol:
                                continue
                            docking_score_val = dr.get("affinity")
                            if docking_score_val is None:
                                continue
                            try:
                                props = compute_properties(mol["smiles"])
                                ha_count = props.get("heavy_atom_count")
                            except Exception:
                                ha_count = None
                            le = compute_ligand_efficiency(docking_score_val, ha_count)
                            if le is not None:
                                le_entries.append((mol["id"], run_uuid, "ligand_efficiency", {"ligand_efficiency": le}))
                        if le_entries:
                            _db_sync(_store_properties_batch(le_entries))
                            _write_log(run_id, "info", f"Ligand efficiency computed for {len(le_entries)} molecules")

                    _write_log(run_id, "info",
                        f"Docking complete: {stored}/{len(molecules)} molecules docked")
                    results["docking"] = f"{stored}/{len(molecules)} molecules docked"
                finally:
                    # Clean up temp directory
                    try:
                        tmp_context.__exit__(None, None, None)
                    except Exception:
                        pass

            elif calc_type == "interactions":
                import tempfile as _tmpfile
                from pipeline.interaction_analysis import analyze_interactions

                # Get protein path from project (same navigation as docking)
                phase_id = molecules[0]["phase_id"] if molecules else None
                docking_ctx = _db_sync(_get_docking_context(phase_id)) if phase_id else None
                protein_path = None
                if docking_ctx:
                    pp = docking_ctx.get("prepared_pdb") or docking_ctx.get("prepared_pdbqt")
                    if pp:
                        protein_path = str(pp)

                # Load existing docking poses from DB
                existing_props = _db_sync(_get_existing_properties(mol_ids))

                inter_entries = []
                skipped = 0
                for mol in molecules:
                    props = existing_props.get(mol["id"], {})
                    pose_molblock = props.get("pose_molblock")

                    if not pose_molblock:
                        logger.info("Skipping %s — no docked pose available", mol["smiles"][:40])
                        skipped += 1
                        continue

                    # Write pose to temp file
                    ligand_path = None
                    if pose_molblock:
                        tmp = _tmpfile.NamedTemporaryFile(suffix=".sdf", delete=False, mode="w")
                        tmp.write(pose_molblock)
                        tmp.close()
                        ligand_path = tmp.name

                    try:
                        _uniprot_id = docking_ctx.get("uniprot_id", "") if docking_ctx else ""
                        _func_residues = docking_ctx.get("functional_residues") if docking_ctx else None

                        # Apply excluded_residues from run config
                        _excluded = _run_config.get("excluded_residues", [])
                        if _excluded and _func_residues and _func_residues.get("residues"):
                            _func_residues = dict(_func_residues)
                            _func_residues["residues"] = [
                                r for r in _func_residues["residues"]
                                if (r["number"] if isinstance(r, dict) else r) not in _excluded
                            ]
                            if _func_residues.get("key_hbond_residues"):
                                _func_residues["key_hbond_residues"] = [
                                    n for n in _func_residues["key_hbond_residues"] if n not in _excluded
                                ]

                        result = analyze_interactions(
                            protein_path=protein_path or "",
                            ligand_path=ligand_path or "",
                            uniprot_id=_uniprot_id,
                            smiles=mol["smiles"],
                            functional_residues=_func_residues,
                        )
                        _inters = result.get("interactions", [])
                        inter_data = {
                            "n_interactions": len({i.get("residue_number") for i in _inters}),
                            "functional_contacts": result.get("functional_contacts", 0),
                            "total_functional": result.get("total_functional", 0),
                            "interaction_quality": result.get("interaction_quality"),
                            "key_hbonds": result.get("key_hbonds", 0),
                            "interactions_method": result.get("method", "unknown"),
                            "interactions_detail": result.get("interactions", []),
                            "mapped_functional": result.get("mapped_functional", []),
                        }
                        inter_entries.append((mol["id"], run_uuid, "interactions", _filter_cols(inter_data, _inc)))
                    except Exception as e:
                        logger.warning("Interactions failed for %s: %s", mol["smiles"], e)
                    finally:
                        if ligand_path:
                            Path(ligand_path).unlink(missing_ok=True)

                _db_sync(_store_properties_batch(inter_entries))
                _write_log(run_id, "info", f"Interactions analysis complete: {len(inter_entries)}/{len(molecules)} molecules ({skipped} skipped — no pose)")
                results["interactions"] = f"{len(inter_entries)} molecules analyzed"

            elif calc_type == "scaffold":
                from pipeline.scaffold_analysis import analyze_scaffold

                scaffold_entries = []
                for mol in molecules:
                    try:
                        result = analyze_scaffold(mol["smiles"])
                        stats = result.get("stats", {})
                        scaffold_data = {
                            "scaffold_smiles": result.get("scaffold_smiles"),
                            "n_modifiable_positions": stats.get("n_positions", 0),
                            "brics_bond_count": stats.get("n_brics_bonds", 0),
                            "scaffold_n_rings": stats.get("n_rings", 0),
                            "scaffold_mw": stats.get("scaffold_mw"),
                            "scaffold_positions": result.get("positions", []),
                            "scaffold_svg": result.get("annotated_svg", ""),
                        }
                        scaffold_entries.append((mol["id"], run_uuid, "scaffold", _filter_cols(scaffold_data, _inc)))
                    except Exception as e:
                        logger.warning("Scaffold analysis failed for %s: %s", mol["smiles"], e)

                # Post-processing: group by shared scaffold (chemical series)
                from collections import Counter as _ScaffoldCounter
                scaffold_map = {}
                group_counter = 0
                for entry in scaffold_entries:
                    smi = entry[3].get("scaffold_smiles")
                    if smi and smi not in scaffold_map:
                        scaffold_map[smi] = group_counter
                        group_counter += 1

                group_sizes = _ScaffoldCounter()
                for entry in scaffold_entries:
                    smi = entry[3].get("scaffold_smiles")
                    if smi:
                        group_sizes[scaffold_map[smi]] += 1

                for i, (mol_id, ruuid, ptype, props) in enumerate(scaffold_entries):
                    smi = props.get("scaffold_smiles")
                    gid = scaffold_map.get(smi, -1) if smi else -1
                    props["scaffold_group"] = gid
                    props["scaffold_group_size"] = group_sizes.get(gid, 1)
                    scaffold_entries[i] = (mol_id, ruuid, ptype, props)

                _db_sync(_store_properties_batch(scaffold_entries))
                _write_log(run_id, "info", f"Scaffold analysis complete: {len(scaffold_entries)}/{len(molecules)} molecules, {group_counter} series")
                results["scaffold"] = f"{len(scaffold_entries)} molecules analyzed, {group_counter} series"

                # --- Diversity clustering (Butina on Morgan FP) — merged into scaffold ---
                from pipeline.scoring import cluster_results as _cluster_results
                _scaffold_cfg = _run_config.get("scaffold", {}) if isinstance(_run_config.get("scaffold"), dict) else {}
                _cutoff = float(_scaffold_cfg.get("cutoff", _run_config.get("cutoff", 0.65)))
                mol_dicts = [{"smiles": m["smiles"], "name": m["name"],
                              "composite_score": m.get("composite_score", 0)} for m in molecules]
                try:
                    clustered = _cluster_results(mol_dicts, cutoff=_cutoff)
                    from collections import Counter as _ClusterCounter
                    _cluster_counts = _ClusterCounter(c.get("cluster_id", 0) for c in clustered)
                    n_clusters = len(_cluster_counts)
                    cluster_entries = [
                        (mol["id"], run_uuid, "clustering", _filter_cols(
                            {"cluster_id": cl.get("cluster_id"),
                             "cluster_size": _cluster_counts.get(cl.get("cluster_id", 0), 0),
                             "is_representative": cl.get("is_representative"),
                             "tanimoto_to_centroid": cl.get("tanimoto_to_centroid")}, _inc))
                        for mol, cl in zip(molecules, clustered)
                    ]
                    _db_sync(_store_properties_batch(cluster_entries))
                    _write_log(run_id, "info", f"Clustering: {n_clusters} clusters (cutoff={_cutoff})")
                except Exception as e:
                    logger.warning("Clustering within scaffold failed: %s", e)
                    _write_log(run_id, "warning", f"Clustering step failed: {e}")

            elif calc_type == "clustering":
                from pipeline.scoring import cluster_results
                mol_dicts = [{"smiles": m["smiles"], "name": m["name"], "composite_score": m.get("composite_score", 0)} for m in molecules]
                try:
                    _cutoff = float(_run_config.get("cutoff", 0.4))
                    clustered = cluster_results(mol_dicts, cutoff=_cutoff)
                    # Compute cluster sizes
                    from collections import Counter
                    _cluster_counts = Counter(c.get("cluster_id", 0) for c in clustered)
                    n_clusters = len(_cluster_counts)
                    cluster_entries = [
                        (mol["id"], run_uuid, "clustering", _filter_cols(
                            {"cluster_id": cl_data.get("cluster_id"),
                             "cluster_size": _cluster_counts.get(cl_data.get("cluster_id", 0), 0),
                             "is_representative": cl_data.get("is_representative"),
                             "scaffold_smiles": cl_data.get("scaffold_smiles", cl_data.get("scaffold")),
                             "tanimoto_to_centroid": cl_data.get("tanimoto_to_centroid")}, _inc))
                        for mol, cl_data in zip(molecules, clustered)
                    ]
                    _db_sync(_store_properties_batch(cluster_entries))
                    _write_log(run_id, "info", f"Clustering complete: {n_clusters} clusters from {len(molecules)} molecules")
                    results["clustering"] = f"{len(molecules)} molecules clustered"
                except Exception as e:
                    logger.warning("Clustering failed: %s", e)
                    results["clustering"] = f"failed: {e}"

            elif calc_type == "off_target":
                # Off-target selectivity analysis against 10 anti-target panel
                from pipeline.off_target import screen_off_targets
                import tempfile
                tmp_dir = Path(tempfile.mkdtemp(prefix="off_target_"))
                ot_entries = []
                for mol in molecules:
                    try:
                        ot = screen_off_targets(mol["smiles"], tmp_dir)
                        off_target_hits = ot["n_total"] - ot["n_safe"]
                        selectivity_ratio = ot["n_safe"] / ot["n_total"] if ot["n_total"] > 0 else 0.0
                        ot_cols = {"selectivity_score": ot["selectivity_score"],
                                   "off_target_hits": off_target_hits,
                                   "selectivity_ratio": round(selectivity_ratio, 3)}
                        # Metadata needed by detail panel — never filtered
                        ot_meta = {"n_safe": ot["n_safe"],
                                   "n_total": ot["n_total"],
                                   "results": ot["results"],
                                   "warnings": ot["warnings"]}
                        ot_filtered = _filter_cols(ot_cols, _inc)
                        ot_filtered.update(ot_meta)
                        ot_entries.append((mol["id"], run_uuid, "off_target", ot_filtered))
                    except Exception as e:
                        logger.warning("Off-target failed for %s: %s", mol["smiles"], e)
                _db_sync(_store_properties_batch(ot_entries))
                _write_log(run_id, "info", f"Off-target screening complete: {len(molecules)} molecules vs 10 anti-targets")
                results["off_target"] = f"{len(molecules)} molecules analyzed"

            elif calc_type == "confidence":
                # Multi-component confidence scoring
                from pipeline.confidence import calculate_confidence
                from pipeline.scoring import compute_pains_alert
                existing_props = _db_sync(_get_existing_properties(mol_ids))
                conf_entries = []
                for mol in molecules:
                    try:
                        mol_dict = {"smiles": mol["smiles"], "name": mol["name"]}
                        props = existing_props.get(mol["id"], {})
                        for prop_key in ("affinity", "vina_score", "cnn_score", "admet",
                                         "synthesis_route", "docking_method", "source"):
                            if props.get(prop_key) is not None:
                                mol_dict[prop_key] = props[prop_key]
                        conf = calculate_confidence(mol_dict)
                        pains = compute_pains_alert(mol["smiles"])
                        conf["pains_alert"] = pains
                        conf["confidence_score"] = conf["overall"]
                        conf["confidence_flags"] = ["pains"] if pains else []
                        conf_entries.append((mol["id"], run_uuid, "confidence", _filter_cols(conf, _inc)))
                    except Exception as e:
                        logger.warning("Confidence failed for %s: %s", mol["smiles"], e)
                _db_sync(_store_properties_batch(conf_entries))
                _write_log(run_id, "info", f"Confidence scoring complete: {len(molecules)} molecules assessed")
                results["confidence"] = f"{len(molecules)} molecules assessed"

            elif calc_type == "retrosynthesis":
                # Full retrosynthesis planning via pipeline/retrosynthesis.py
                from pipeline.retrosynthesis import plan_synthesis
                retro_cfg = _run_config.get("retrosynthesis", {}) if isinstance(_run_config.get("retrosynthesis"), dict) else {}
                retro_method = retro_cfg.get("method", "auto")
                _write_log(run_id, "info", f"Retrosynthesis method: {retro_method}")
                retro_entries = []
                for mol in molecules:
                    try:
                        route = plan_synthesis(mol["smiles"], method=retro_method)
                        retro_cols = {
                            "n_synth_steps": route.get("n_steps", 0),
                            "synth_confidence": route.get("confidence", 0),
                            "synth_cost_estimate": route.get("estimated_cost"),
                            "reagents_available": route.get("all_reagents_available", True),
                        }
                        # Metadata needed by detail panel — never filtered
                        retro_meta = {
                            "n_steps": route.get("n_steps", 0),
                            "confidence": route.get("confidence", 0),
                            "estimated_cost": route.get("estimated_cost"),
                            "all_reagents_available": route.get("all_reagents_available"),
                            "steps": route.get("steps", []),
                            "tree": route.get("tree"),
                            "cost_estimate": route.get("cost_estimate"),
                            "reagent_availability": route.get("reagent_availability", []),
                            "method": route.get("method", "unknown"),
                        }
                        retro_filtered = _filter_cols(retro_cols, _inc)
                        retro_filtered.update(retro_meta)
                        retro_entries.append((mol["id"], run_uuid, "retrosynthesis", retro_filtered))
                    except Exception as e:
                        logger.warning("Retrosynthesis failed for %s: %s", mol["smiles"], e)
                        _write_log(run_id, "warning", f"No route found for {mol.get('name', mol['smiles'][:30])}: {e}")
                        # Store a "no route" result so the user sees something
                        from pipeline.retrosynthesis import _empty_result
                        fallback = _empty_result(mol["smiles"])
                        retro_entries.append((mol["id"], run_uuid, "retrosynthesis", {
                            "n_synth_steps": 0,
                            "synth_confidence": 0.0,
                            "synth_cost_estimate": None,
                            "reagents_available": None,
                            "n_steps": 0,
                            "confidence": 0.0,
                            "estimated_cost": None,
                            "all_reagents_available": None,
                            "steps": [],
                            "tree": fallback.get("tree"),
                            "method": "failed",
                        }))
                _db_sync(_store_properties_batch(retro_entries))
                n_success = sum(1 for _, _, _, d in retro_entries if d.get("n_steps", 0) > 0)
                _write_log(run_id, "info", f"Retrosynthesis complete: {n_success}/{len(molecules)} routes found")
                results["retrosynthesis"] = f"{len(molecules)} molecules analyzed"

            elif calc_type == "pharmacophore":
                from pipeline.pharmacophore import score_against_reference

                # Fetch all reference ligands from campaign
                ref_ligands = _db_sync(_get_reference_ligands(phase_id))
                if not ref_ligands:
                    _write_log(run_id, "error", "No reference ligands set. Go to Project Home → Campaign to add reference ligands.")
                    results["pharmacophore"] = "FAILED — no reference ligands"
                    continue

                pharma_entries = []
                ref_names = []
                for mol in molecules:
                    props = {"feature_counts": None}
                    fit_values = []
                    for ref_idx, ref_lig in enumerate(ref_ligands):
                        ref_smiles = ref_lig.get("smiles", "")
                        if not ref_smiles:
                            continue
                        try:
                            result = score_against_reference(mol["smiles"], ref_smiles)
                            fit_key = f"pharmacophore_fit_{ref_idx + 1}"
                            props[fit_key] = result["pharmacophore_fit"]
                            fit_values.append(result["pharmacophore_fit"])
                            # Store per-ref metadata (matched/missing types from last ref for feature_counts)
                            if props["feature_counts"] is None:
                                props["feature_counts"] = result["feature_counts"]
                        except Exception as e:
                            logger.warning("Pharmacophore scoring failed for %s vs ref %d: %s", mol["smiles"], ref_idx, e)
                    # Best + average fit across all references
                    if fit_values:
                        props["pharmacophore_fit_best"] = max(fit_values)
                        props["pharmacophore_fit_avg"] = round(sum(fit_values) / len(fit_values), 4)
                    # Store ref_feature_counts from all references for detail panel
                    props["ref_feature_counts"] = {
                        f"ref_{i+1}": rl.get("ref_feature_counts", {})
                        for i, rl in enumerate(ref_ligands)
                    }
                    pharma_entries.append((mol["id"], run_uuid, "pharmacophore", props))

                _db_sync(_store_properties_batch(pharma_entries))

                ref_names = [rl.get("name", f"ref {i+1}") for i, rl in enumerate(ref_ligands)]
                _write_log(run_id, "info", f"Reference similarity complete: {len(molecules)} molecules scored against {len(ref_ligands)} references ({', '.join(ref_names)})")
                results["pharmacophore"] = f"{len(molecules)} molecules scored against {len(ref_ligands)} references"

            elif calc_type == "composite":
                from pipeline.scoring import compute_weighted_composite, DEFAULT_COMPOSITE_WEIGHTS

                weights = _run_config.get("weights", DEFAULT_COMPOSITE_WEIGHTS)

                # Fetch existing properties for all input molecules
                mol_ids = [m["id"] for m in molecules]
                existing_props = _db_sync(_get_existing_properties(mol_ids))

                composite_entries = []
                scored_count = 0
                for mol in molecules:
                    props = existing_props.get(mol["id"], {})
                    # Map property keys to weight keys
                    mapped = {}
                    if "affinity" in props:
                        mapped["docking_score"] = props["affinity"]
                    elif "docking_score" in props:
                        mapped["docking_score"] = props["docking_score"]
                    if "cnn_score" in props:
                        mapped["cnn_score"] = props["cnn_score"]
                    if "logP" in props or "logp" in props:
                        mapped["logP"] = props.get("logP") or props.get("logp")
                    if "solubility" in props:
                        mapped["solubility"] = props["solubility"]
                    if "selectivity_score" in props:
                        mapped["selectivity"] = props["selectivity_score"]
                    if "qed" in props:
                        mapped["qed"] = props["qed"]
                    if "composite_score" in props:
                        mapped["safety"] = props["composite_score"]
                    if "novelty_score" in props or "novelty" in props:
                        mapped["novelty"] = props.get("novelty_score") or props.get("novelty")

                    final, breakdown = compute_weighted_composite(mapped, weights)

                    composite_entries.append((mol["id"], run_uuid, "composite", {
                        "weighted_score": final,
                        "breakdown": breakdown,
                        "weights_used": weights,
                    }))
                    if final is not None:
                        scored_count += 1

                _db_sync(_store_properties_batch(composite_entries))
                _write_log(run_id, "info",
                    f"Composite score complete: {scored_count}/{len(molecules)} molecules scored")
                results["composite"] = f"{scored_count} molecules scored"

            else:
                logger.warning("Unknown calculation type: %s", calc_type)
                results[calc_type] = "unknown type, skipped"

        except Exception as e:
            logger.exception("Calculation %s failed: %s", calc_type, e)
            _write_log(run_id, "error", f"{calc_type} failed: {str(e)[:500]}")
            results[calc_type] = f"error: {e}"

    return {"calculations": results, "molecules_processed": len(molecules)}


# ---------------------------------------------------------------------------
# Generation run
# ---------------------------------------------------------------------------

def _run_generation(run_id: str, run_data: dict) -> dict:
    """Generate new molecules using pipeline/generation.py (REINVENT4 → RDKit mock → hardcoded)."""
    import tempfile
    import requests as http_requests

    config = run_data.get("config") or {}
    phase_id = run_data["phase_id"]
    run_uuid = uuid.UUID(run_id)
    input_ids = run_data.get("input_molecule_ids") or []

    # Dispatch to optimization handler if mode is optimization
    if config.get("mode") == "optimization":
        return _run_lead_optimization(run_id, run_data)

    # Per-molecule mode with scaffold rules → uses lead optimization with structural rules
    if config.get("mode") == "molecule" and config.get("scaffold_rules"):
        return _run_lead_optimization(run_id, run_data)

    n_molecules = config.get("n_molecules", 100)
    n_top = config.get("n_top", 20)
    strategy = config.get("strategy", "scaffold_hop")

    _db_sync(_update_run(run_id, current_step="Preparing generation", progress=5))
    _write_log(run_id, "info", f"Generation: strategy={strategy}, n_molecules={n_molecules}, n_top={n_top}")

    # Get parent molecules as seed SMILES
    parent_mols = []
    seed_smiles = []
    if input_ids:
        parent_mols = _db_sync(_get_molecules_by_ids(input_ids))
        seed_smiles = [m["smiles"] for m in parent_mols if m.get("smiles")]

    logger.info(
        "Generation run %s: strategy=%s, n_molecules=%d, seeds=%d",
        run_id, strategy, n_molecules, len(seed_smiles),
    )

    # Get receptor + pocket context (same as docking)
    docking_ctx = _db_sync(_get_docking_context(phase_id)) if phase_id else None

    if not docking_ctx or not docking_ctx.get("pdb_url"):
        raise ValueError(
            "Generation requires a target with a PDB structure. "
            "Please configure the target in Target Setup first."
        )
    if not docking_ctx.get("pocket_center"):
        raise ValueError(
            "No binding pocket center found. "
            "Please detect pockets in Target Setup first."
        )

    pdb_url = docking_ctx["pdb_url"]
    pocket_center = docking_ctx["pocket_center"]
    pocket_size = tuple(docking_ctx["pocket_config"].get("size", [22, 22, 22]))

    _db_sync(_update_run(run_id, current_step="Downloading receptor", progress=10))

    with tempfile.TemporaryDirectory(prefix="generation_v9_") as tmp_dir:
        tmp_path = Path(tmp_dir)
        pdb_path = tmp_path / f"{docking_ctx.get('pdb_id', 'receptor')}.pdb"

        # Download PDB
        resp = http_requests.get(pdb_url, timeout=(10, 60))
        resp.raise_for_status()
        pdb_path.write_text(resp.text)

        _db_sync(_update_run(
            run_id,
            current_step=f"Generating {n_molecules} molecules",
            progress=20,
        ))
        method = config.get("method", "scaffold_hopping")
        _write_log(run_id, "info",
            f"Calling generation: method={method}, center={pocket_center}, seeds={len(seed_smiles)}")

        # Dispatch by generation method
        if method == "fragment_growing":
            from pipeline.generation import generate_fragment_growing
            gen_results = generate_fragment_growing(
                parent_smiles=seed_smiles,
                pocket_center=pocket_center,
                n_molecules=n_molecules,
                n_top=n_top,
            )
        elif method == "de_novo":
            from pipeline.generation import generate_molecules
            gen_results = generate_molecules(
                pocket_center=pocket_center,
                pocket_size=pocket_size,
                receptor_pdbqt=pdb_path,
                work_dir=tmp_path,
                n_molecules=max(n_molecules, 200),
                n_top=min(n_top, 50),
                seed_smiles=None,  # unconstrained de novo
            )
        elif method == "bioisosteric":
            from pipeline.generation import generate_bioisosteric
            gen_results = generate_bioisosteric(
                parent_smiles=seed_smiles,
                pocket_center=pocket_center,
                n_molecules=n_molecules,
                n_top=n_top,
            )
        elif method == "fragment_linking":
            from pipeline.generation import generate_fragment_linking
            gen_results = generate_fragment_linking(
                parent_smiles=seed_smiles,
                pocket_center=pocket_center,
                n_molecules=n_molecules,
                n_top=n_top,
            )
        elif method == "matched_pairs":
            from pipeline.generation import generate_matched_pairs
            gen_results = generate_matched_pairs(
                parent_smiles=seed_smiles,
                pocket_center=pocket_center,
                n_molecules=n_molecules,
                n_top=n_top,
            )
        else:
            # Default: scaffold hopping
            from pipeline.generation import generate_molecules
            gen_results = generate_molecules(
                pocket_center=pocket_center,
                pocket_size=pocket_size,
                receptor_pdbqt=pdb_path,
                work_dir=tmp_path,
                n_molecules=n_molecules,
                n_top=n_top,
                seed_smiles=seed_smiles if seed_smiles else None,
            )

        _db_sync(_update_run(
            run_id,
            current_step=f"Storing {len(gen_results)} generated molecules",
            progress=80,
        ))

        # Store generated molecules in the phase
        created = 0
        for i, gen in enumerate(gen_results):
            smiles = gen.get("smiles", "")
            if not smiles:
                continue

            try:
                from rdkit import Chem
                canonical = Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) if Chem.MolFromSmiles(smiles) else smiles
            except Exception:
                canonical = smiles

            name = gen.get("name", f"GEN_{i+1:03d}")

            # Find parent molecule and inherit generation level
            parent_id = None
            parent_gen_level = 0
            if parent_mols:
                parent = parent_mols[i % len(parent_mols)]
                parent_id = parent["id"]
                parent_gen_level = parent.get("generation_level", 0)

            mol_id = _db_sync(_create_molecule(
                phase_id=phase_id,
                smiles=smiles,
                canonical_smiles=canonical,
                name=name,
                source_run_id=run_uuid,
                ai_generated=True,
                parent_molecule_id=parent_id,
                generation_level=parent_gen_level + 1,
            ))

            if mol_id:
                created += 1
                # Store generation properties on the new molecule
                gen_props = {
                    "estimated_affinity": gen.get("estimated_affinity"),
                    "qed": gen.get("qed"),
                    "novelty_score": gen.get("novelty_score"),
                    "lipinski_violations": gen.get("lipinski_violations"),
                    "source": gen.get("source", "generation"),
                }
                _db_sync(_store_property(mol_id, run_uuid, "generation", gen_props))

    # --- Post-generation analyses (chained) ---
    post_analyses = []
    if config.get("include_docking"):
        post_analyses.append("docking")
    if config.get("include_admet"):
        post_analyses.append("admet")
    if config.get("include_scoring"):
        post_analyses.append("scoring")

    if post_analyses and created > 0:
        _db_sync(_update_run(run_id, current_step=f"Running post-analyses: {', '.join(post_analyses)}", progress=85))
        _write_log(run_id, "info", f"Running post-generation analyses: {post_analyses}")

        # Get the generated molecule IDs
        gen_mol_ids = []
        for i, gen in enumerate(gen_results):
            smi = gen.get("smiles", "")
            if smi:
                gen_mol_ids.append(smi)

        # Build a minimal run_data for calculation
        post_run_data = {
            "phase_id": phase_id,
            "config": {
                "calculation_types": post_analyses,
            },
        }

        try:
            # Get all generated molecules from the phase for this run
            gen_molecules = _db_sync(_get_molecules_by_run(run_id, phase_id))
            if gen_molecules:
                post_run_data["input_molecule_ids"] = [m["id"] for m in gen_molecules]
                _run_calculation(run_id, post_run_data)
                _write_log(run_id, "info", f"Post-generation analyses complete: {post_analyses}")
        except Exception as e:
            logger.warning("Post-generation analyses failed: %s", e)
            _write_log(run_id, "warning", f"Post-generation analyses failed: {e}")

    _db_sync(_update_run(run_id, progress=95, current_step="Finalizing"))
    _write_log(run_id, "info", f"Generation complete: {created}/{len(gen_results)} molecules created")
    logger.info("Generation run %s: %d molecules created", run_id, created)
    return {"generated": created, "strategy": strategy, "total_candidates": len(gen_results)}


# ---------------------------------------------------------------------------
# Lead Optimization run
# ---------------------------------------------------------------------------

def _run_lead_optimization(run_id: str, run_data: dict) -> dict:
    """Multi-iteration lead optimization: generate variants → dock → score → iterate."""
    import tempfile
    import requests as http_requests

    config = run_data.get("config") or {}
    phase_id = run_data["phase_id"]
    run_uuid = uuid.UUID(run_id)
    input_ids = run_data.get("input_molecule_ids") or []

    _db_sync(_update_run(run_id, current_step="Preparing optimization", progress=5))
    _write_log(run_id, "info", "Lead optimization started")

    # Get parent molecule (optimization works on a single seed)
    parent_mols = []
    seed_smiles = []
    if input_ids:
        parent_mols = _db_sync(_get_molecules_by_ids(input_ids))
        seed_smiles = [m["smiles"] for m in parent_mols if m.get("smiles")]

    if not seed_smiles:
        raise ValueError("Lead optimization requires at least one seed molecule.")

    starting_smiles = seed_smiles[0]
    starting_name = parent_mols[0].get("name", "MOL") if parent_mols else "MOL"

    # Get receptor + pocket context
    docking_ctx = _db_sync(_get_docking_context(phase_id)) if phase_id else None

    if not docking_ctx or not docking_ctx.get("pdb_url"):
        raise ValueError("Optimization requires a target with a PDB structure.")
    if not docking_ctx.get("pocket_center"):
        raise ValueError("No binding pocket center found.")

    pdb_url = docking_ctx["pdb_url"]
    pocket_center = docking_ctx["pocket_center"]
    pocket_size = tuple(docking_ctx["pocket_config"].get("size", [22, 22, 22]))

    # Optimization parameters from config
    n_iterations = config.get("iterations", 5)
    variants_per_iter = config.get("variants_per_iteration", 50)
    weights = config.get("weights", {
        "binding_affinity": 0.35,
        "toxicity": 0.25,
        "bioavailability": 0.20,
        "synthesis_ease": 0.20,
    })
    structural_rules = config.get("scaffold_rules")

    _db_sync(_update_run(run_id, current_step="Downloading receptor", progress=10))

    with tempfile.TemporaryDirectory(prefix="optimization_v9_") as tmp_dir:
        tmp_path = Path(tmp_dir)
        pdb_path = tmp_path / f"{docking_ctx.get('pdb_id', 'receptor')}.pdb"

        # Download PDB
        resp = http_requests.get(pdb_url, timeout=(10, 60))
        resp.raise_for_status()
        pdb_path.write_text(resp.text)

        _db_sync(_update_run(
            run_id,
            current_step=f"Optimizing: iteration 0/{n_iterations}",
            progress=15,
        ))
        _write_log(run_id, "info",
            f"Starting optimization: {n_iterations} iterations, {variants_per_iter} variants/iter, "
            f"seed={starting_name}, weights={weights}")

        # Progress callback: maps iteration progress to 15-85% of run progress
        def progress_cb(iteration, best_score, message):
            pct = 15 + int(70 * iteration / n_iterations)
            _db_sync(_update_run(
                run_id,
                current_step=f"Optimizing: iteration {iteration}/{n_iterations} (score={best_score:.3f})",
                progress=pct,
            ))
            _write_log(run_id, "info", message)

        from pipeline.lead_optimization import run_optimization
        opt_result = run_optimization(
            starting_smiles=starting_smiles,
            starting_name=starting_name,
            target_pdbqt=str(pdb_path),
            pocket_center=pocket_center,
            weights=weights,
            n_iterations=n_iterations,
            variants_per_iter=variants_per_iter,
            work_dir=tmp_path,
            progress_callback=progress_cb,
            structural_rules=structural_rules,
            pocket_size=list(pocket_size),
        )

        _db_sync(_update_run(
            run_id,
            current_step="Storing optimized molecules",
            progress=85,
        ))

        # Store top molecules + final lead in the phase
        created = 0
        top_molecules = opt_result.get("top_molecules", [])
        final_lead = opt_result.get("final_lead", {})

        # Ensure final lead is in top_molecules
        final_smiles = final_lead.get("smiles", "")
        top_smiles_set = {m.get("smiles") for m in top_molecules}
        if final_smiles and final_smiles not in top_smiles_set:
            top_molecules.insert(0, final_lead)

        parent_id = parent_mols[0]["id"] if parent_mols else None
        parent_gen_level = parent_mols[0].get("generation_level", 0) if parent_mols else 0

        for i, mol_data in enumerate(top_molecules):
            smiles = mol_data.get("smiles", "")
            if not smiles:
                continue

            try:
                from rdkit import Chem
                canonical = Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) if Chem.MolFromSmiles(smiles) else smiles
            except Exception:
                canonical = smiles

            name = mol_data.get("name", f"OPT_{i+1:03d}")

            mol_id = _db_sync(_create_molecule(
                phase_id=phase_id,
                smiles=smiles,
                canonical_smiles=canonical,
                name=name,
                source_run_id=run_uuid,
                ai_generated=True,
                parent_molecule_id=parent_id,
                generation_level=parent_gen_level + 1,
            ))

            if mol_id:
                created += 1
                # Store optimization properties
                opt_props = {
                    "optimization_score": mol_data.get("score"),
                    "estimated_affinity": mol_data.get("affinity"),
                    "iteration": mol_data.get("iteration"),
                    "is_final_lead": (smiles == final_smiles),
                    "source": "lead_optimization",
                }
                # Include ADMET data if present
                admet = mol_data.get("admet", {})
                if admet:
                    opt_props["toxicity_score"] = admet.get("toxicity_score")
                    opt_props["bioavailability_score"] = admet.get("bioavailability_score")
                    opt_props["synthesis_score"] = admet.get("synthesis_score")
                # Include docking metadata
                for dk in ("vina_score", "cnn_score", "cnn_affinity", "docking_engine"):
                    if mol_data.get(dk) is not None:
                        opt_props[dk] = mol_data[dk]

                _db_sync(_store_property(mol_id, run_uuid, "generation", opt_props))

    # Store optimization summary in run results
    summary = {
        "mode": "optimization",
        "generated": created,
        "total_tested": opt_result.get("total_tested", 0),
        "total_improved": opt_result.get("total_improved", 0),
        "iterations_completed": len(opt_result.get("iterations", [])),
        "starting_score": opt_result.get("starting_molecule", {}).get("score"),
        "final_score": final_lead.get("score"),
        "docking_mode": opt_result.get("docking_mode", "mock"),
        "objectives": opt_result.get("objectives", {}),
        "comparison": opt_result.get("comparison", {}),
    }

    _db_sync(_update_run(run_id, progress=95, current_step="Finalizing"))
    _write_log(run_id, "info",
        f"Optimization complete: {created} molecules stored, "
        f"score {summary.get('starting_score', '?')} → {summary.get('final_score', '?')}, "
        f"{summary['total_tested']} variants tested across {summary['iterations_completed']} iterations")
    logger.info("Optimization run %s: %d molecules created", run_id, created)
    return summary


# ---------------------------------------------------------------------------
# AFVS run (BDX-41) — Adaptive Focused Virtual Screening
# ---------------------------------------------------------------------------

async def _get_afvs_job(run_id: str) -> dict | None:
    """Fetch AFVS job details for a run."""
    from database_v9 import AsyncSessionV9
    from models_v9 import AFVSJobORM_V9
    from sqlalchemy import select

    async with AsyncSessionV9() as session:
        result = await session.execute(
            select(AFVSJobORM_V9).where(
                AFVSJobORM_V9.run_id == uuid.UUID(run_id)
            )
        )
        job = result.scalar_one_or_none()
        if job:
            return {
                "id": job.id,
                "run_id": job.run_id,
                "receptor_s3_path": job.receptor_s3_path,
                "docking_box": job.docking_box,
                "docking_scenario": job.docking_scenario,
                "top_n_results": job.top_n_results,
                "budget_cap_usd": job.budget_cap_usd,
                "s3_output_prefix": job.s3_output_prefix,
                "afvs_phase": job.afvs_phase,
                "molecules_screened": job.molecules_screened,
                "molecules_imported": job.molecules_imported,
                "actual_cost_usd": job.actual_cost_usd,
            }
    return None


async def _update_afvs_job(run_id: str, **fields):
    """Update AFVS job fields."""
    from database_v9 import AsyncSessionV9
    from models_v9 import AFVSJobORM_V9
    from sqlalchemy import select

    async with AsyncSessionV9() as session:
        result = await session.execute(
            select(AFVSJobORM_V9).where(
                AFVSJobORM_V9.run_id == uuid.UUID(run_id)
            )
        )
        job = result.scalar_one_or_none()
        if job:
            for k, v in fields.items():
                setattr(job, k, v)
            await session.commit()


def _run_afvs(run_id: str, run_data: dict) -> dict:
    """Execute AFVS screening: SSH to EC2, submit AWS Batch jobs, poll, import results.

    Single-phase flow (AdaptiveFlow v1):
      1. Load docking context (receptor + pocket box)
      2. SSH connect → cleanup previous → setup_job (upload receptor, write configs, prepare workunits)
      3. Submit all workunits to AWS Batch
      4. Poll job status every 60s (72h timeout, cancel/budget checks)
      5. Collect top N results from S3 CSV.gz files
      6. Download + parse results CSV
      7. Import molecules + docking scores into Phase A
    """
    import time as _time
    from pipeline.afvs import (
        AFVSController, parse_afvs_results_csv,
        download_s3_csv, download_s3_csv_via_ec2,
    )

    config = run_data.get("config", {})
    phase_id = run_data["phase_id"]
    str_run_id = str(run_id)
    job_name = None
    scenario_name = None

    # Step 1: Load docking context
    _db_sync(_update_run(str_run_id, progress=5, current_step="Loading receptor"))
    _write_log(str_run_id, "info", "Fetching docking context (receptor + pocket)")

    ctx = _db_sync(_get_docking_context(phase_id))
    if not ctx:
        raise RuntimeError("Docking context not found — configure target + pocket first")

    afvs_job = _db_sync(_get_afvs_job(str_run_id))
    if not afvs_job:
        raise RuntimeError("AFVS job record not found")

    prepared_pdb = ctx.get("prepared_pdb")
    if not prepared_pdb or not Path(str(prepared_pdb)).exists():
        raise RuntimeError("Prepared receptor PDB not found on disk")

    # Build docking box from pocket config
    pocket = ctx.get("pocket_config") or {}
    docking_box = afvs_job.get("docking_box") or {
        "center_x": pocket.get("center", [0, 0, 0])[0],
        "center_y": pocket.get("center", [0, 0, 0])[1],
        "center_z": pocket.get("center", [0, 0, 0])[2],
        "size_x": pocket.get("size", [22, 22, 22])[0],
        "size_y": pocket.get("size", [22, 22, 22])[1],
        "size_z": pocket.get("size", [22, 22, 22])[2],
    }

    controller = AFVSController()

    try:
        # Step 2: Connect + setup job on EC2
        _db_sync(_update_run(str_run_id, progress=10, current_step="Connecting to EC2"))
        _write_log(str_run_id, "info", "Connecting to AFVS EC2 controller")
        controller.connect()

        _db_sync(_update_run(str_run_id, progress=15, current_step="Setting up AFVS job"))
        _write_log(str_run_id, "info", "Cleaning previous workflow + preparing new job")
        controller.cleanup_previous()

        job_name, scenario_name, num_workunits = controller.setup_job(
            run_id=str_run_id,
            receptor_pdb_path=str(prepared_pdb),
            docking_box=docking_box,
            config=config,
        )
        _write_log(str_run_id, "info",
            f"Job {job_name} prepared: {num_workunits} workunits, scenario {scenario_name}")
        _db_sync(_update_afvs_job(
            str_run_id,
            s3_output_prefix=f"AF/AFVS/jobs/{job_name}",
            afvs_work_units=num_workunits,
        ))

        # Step 3: Submit to AWS Batch
        _db_sync(_update_run(str_run_id, progress=25, current_step="Submitting to AWS Batch"))
        _write_log(str_run_id, "info", f"Submitting {num_workunits} workunits to AWS Batch")
        controller.submit_jobs(num_workunits)
        _db_sync(_update_afvs_job(str_run_id, afvs_phase="primary_screen"))

        # Step 4: Poll completion (72h timeout, 60s interval)
        _db_sync(_update_run(str_run_id, progress=30, current_step="Screening in progress"))
        _write_log(str_run_id, "info", "Polling AWS Batch job status (60s interval, 72h timeout)")

        JOB_TIMEOUT = 72 * 3600  # 72h
        poll_start = _time.time()
        last_pct = 0

        while True:
            # Check cancellation
            fresh = _db_sync(_get_run(str_run_id))
            if fresh and fresh["status"] == "cancelled":
                _write_log(str_run_id, "info", "Cancellation requested — killing all Batch jobs")
                controller.cancel_all_jobs()
                _db_sync(_update_afvs_job(str_run_id, afvs_phase="failed"))
                return {"status": "cancelled", "type": "afvs"}

            # Query batch status
            status = controller.get_job_status()
            if status.get("error"):
                _write_log(str_run_id, "warning", f"Status query issue: {status['error']}")
            else:
                pct = status.get("pct", 0)
                succeeded = status.get("succeeded", 0)
                failed = status.get("failed", 0)
                total = status.get("total", 0)
                running = status.get("running", 0)

                # Update progress (30-75 range for polling phase)
                progress = 30 + int(pct * 0.45)
                if pct != last_pct:
                    _db_sync(_update_run(str_run_id, progress=progress,
                        current_step=f"Screening: {pct}% ({succeeded}/{total} done, {running} running)"))
                    last_pct = pct

                if status.get("done"):
                    _write_log(str_run_id, "info",
                        f"All jobs complete: {succeeded} succeeded, {failed} failed out of {total}")
                    _db_sync(_update_afvs_job(str_run_id,
                        molecules_screened=succeeded * 1000))  # ~1000 ligands/workunit
                    break

            # Timeout check
            elapsed = _time.time() - poll_start
            if elapsed > JOB_TIMEOUT:
                _write_log(str_run_id, "warning", "72h timeout reached — collecting partial results")
                break

            _time.sleep(60)

        # Step 5: Collect top results from S3
        top_n = config.get("top_n_results", 10_000)
        _db_sync(_update_run(str_run_id, progress=78, current_step=f"Collecting top {top_n} results"))
        _write_log(str_run_id, "info", f"Collecting top {top_n} results from S3 CSV.gz files")

        s3_csv_path = controller.collect_top_results(job_name, scenario_name, top_n)
        _write_log(str_run_id, "info", f"Results collected: {s3_csv_path}")

        # Step 6: Download + parse CSV
        _db_sync(_update_run(str_run_id, progress=82, current_step="Downloading results"))
        _write_log(str_run_id, "info", "Downloading results CSV from S3")

        try:
            csv_content = download_s3_csv(s3_csv_path)
        except Exception as e:
            _write_log(str_run_id, "warning",
                f"Direct S3 download failed ({e}), using EC2 fallback")
            csv_content = download_s3_csv_via_ec2(controller, s3_csv_path)

        molecules = parse_afvs_results_csv(csv_content)

        if not molecules:
            _write_log(str_run_id, "warning", "No molecules parsed from AFVS results")
            _db_sync(_update_afvs_job(str_run_id, afvs_phase="completed", molecules_imported=0))
            return {"type": "afvs", "molecules_imported": 0, "molecules_skipped": 0}

        _write_log(str_run_id, "info", f"Parsed {len(molecules)} molecules from AFVS results")

        # Step 7: Import molecules + store docking scores
        _db_sync(_update_run(str_run_id, progress=87,
            current_step=f"Importing {len(molecules)} molecules"))
        _write_log(str_run_id, "info", f"Importing {len(molecules)} molecules into Phase A")

        created, skipped, created_rows = _db_sync(
            _create_molecules_batch(phase_id, molecules, uuid.UUID(str_run_id))
        )
        _write_log(str_run_id, "info", f"Created {created} molecules, {skipped} duplicates skipped")

        # Store docking scores as properties
        if created_rows:
            mol_map = {row[1]: row[0] for row in created_rows}  # canonical_smiles → mol_id
            score_entries = []
            for mol in molecules:
                mol_id = mol_map.get(mol["canonical_smiles"])
                if mol_id and mol.get("docking_score") is not None:
                    score_entries.append((
                        mol_id,
                        uuid.UUID(str_run_id),
                        "docking",
                        {"docking_score": mol["docking_score"], "source": "afvs"},
                    ))
            if score_entries:
                _db_sync(_store_properties_batch(score_entries))
                _write_log(str_run_id, "info", f"Stored {len(score_entries)} docking scores")

        # Step 8: Finalize
        _db_sync(_update_run(str_run_id, progress=95, current_step="Finalizing"))
        _db_sync(_update_afvs_job(
            str_run_id,
            afvs_phase="completed",
            molecules_imported=created,
        ))

        summary = {
            "type": "afvs",
            "molecules_imported": created,
            "molecules_skipped": skipped,
            "top_n_requested": top_n,
            "job_name": job_name,
        }
        _write_log(str_run_id, "success",
            f"AFVS complete: {created} molecules imported from {len(molecules)} results")

        return summary

    except Exception:
        _db_sync(_update_afvs_job(str_run_id, afvs_phase="failed"))
        # Try to cancel any running jobs on failure
        try:
            controller.cancel_all_jobs()
        except Exception:
            pass
        raise
    finally:
        # Cleanup EC2 temporary files (keep scenario for debugging)
        try:
            controller.run_command("rm -f /tmp/bindx_status.py /tmp/bindx_collect.py")
        except Exception:
            pass
