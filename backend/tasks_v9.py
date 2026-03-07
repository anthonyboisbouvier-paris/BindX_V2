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
    """Batch INSERT ... ON CONFLICT DO NOTHING. Returns (created_count, skipped_count).

    Each molecule dict: {"smiles": str, "canonical_smiles": str, "name": str|None}
    Uses a single INSERT with N rows per call — dedup via DB constraint.
    """
    from database_v9 import AsyncSessionV9
    from sqlalchemy import text

    if not molecules:
        return 0, 0

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
            "RETURNING id"
        )
        result = await session.execute(sql, params)
        created = len(result.fetchall())
        await session.commit()
        skipped = len(molecules) - created
        return created, skipped


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
            constraint="uq_mol_props_mol_name_run",
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

        return {
            "pdb_url": pdb_url,
            "pocket_center": pocket_center,
            "pocket_config": pocket_config,
            "pdb_id": project.target_pdb_id,
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
        else:
            raise ValueError(f"Unknown run type: {run_type}")

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

        if has_rdkit:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            canonical = Chem.MolToSmiles(mol)
        else:
            canonical = smi

        counter += 1
        if not name:
            name = f"Mol_{counter}"

        yield {"smiles": smi, "canonical_smiles": canonical, "name": name}


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
            suppl = Chem.ForwardSDMolSupplier(str(file_path))
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


def _import_molecules_chunked(run_id: str, phase_id, source_run_id, mol_iter, total_estimate: int = 0) -> dict:
    """Consume a molecule iterator in batches, insert via _create_molecules_batch.

    Updates progress once per batch (not per molecule).
    Returns {"created": N, "skipped": N}.
    """
    created_total = 0
    skipped_total = 0
    processed = 0
    batch = []

    for mol in mol_iter:
        batch.append(mol)
        if len(batch) >= IMPORT_BATCH_SIZE:
            c, s = _db_sync(_create_molecules_batch(phase_id, batch, source_run_id))
            created_total += c
            skipped_total += s
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
        c, s = _db_sync(_create_molecules_batch(phase_id, batch, source_run_id))
        created_total += c
        skipped_total += s

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
    """Fetch ligands from public databases (ChEMBL, PubChem, ZINC, Enamine, Fragments).

    Runs in Celery worker context. Returns iterator of dicts.
    """
    databases = config.get("databases", [])
    uniprot_id = config.get("uniprot_id", "")
    max_per_source = config.get("max_per_source", 50)
    all_ligands = []

    for db_name in databases:
        try:
            if db_name == "chembl" and uniprot_id:
                from pipeline.ligands import fetch_chembl_ligands
                all_ligands.extend(fetch_chembl_ligands(uniprot_id, max_count=max_per_source))

            elif db_name == "pubchem" and uniprot_id:
                from pipeline.ligands import fetch_pubchem_ligands, resolve_gene_name
                gene = resolve_gene_name(uniprot_id)
                if gene:
                    all_ligands.extend(fetch_pubchem_ligands(gene, max_count=max_per_source))

            elif db_name == "zinc":
                from pipeline.ligands import fetch_zinc_ligands
                all_ligands.extend(fetch_zinc_ligands(max_count=max_per_source))

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
        _db_sync(_update_run(run_id, current_step="Fetching from databases...", progress=5))
        _write_log(run_id, "info", f"Fetching from: {config.get('databases', [])}")

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


def _filter_cols(props: dict, included: set | None) -> dict:
    """Filter a properties dict keeping only keys whose frontend column is included.

    Handles nested dicts (e.g. ADMET absorption/distribution/toxicity sub-dicts)
    by recursively filtering. Pass-through if included is None.
    """
    if included is None:
        return props
    result = {}
    for k, v in props.items():
        # Check if this key (or its frontend alias) is in the included set
        if k in included or _BACKEND_TO_FRONTEND.get(k) in included:
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
    total_steps = len(calc_types)
    results = {}

    def _cached_compute(calc_type, smiles, compute_fn):
        """Try cache first, then compute and cache the result."""
        cache_key = f"{calc_type}:{smiles}"
        cached = _db_sync(_cache_get(cache_key))
        if cached is not None:
            return cached
        result = compute_fn(smiles)
        _db_sync(_cache_set(cache_key, result))
        return result

    for step_idx, calc_type in enumerate(calc_types):
        step_label = f"Running {calc_type} ({step_idx+1}/{total_steps})"
        base_pct = int(100 * step_idx / total_steps) if total_steps > 0 else 0
        _db_sync(_update_run(run_id, current_step=step_label, progress=base_pct))

        try:
            _write_log(run_id, "info", f"Starting {calc_type} on {len(molecules)} molecules")

            if calc_type == "admet":
                from pipeline.admet import predict_admet
                from pipeline.scoring import (
                    compute_properties, compute_cns_mpo,
                    compute_druglikeness_rules, compute_brenk_alert,
                )
                admet_results = predict_admet(smiles_list)
                admet_entries = []
                druglike_entries = []
                for mol, admet in zip(molecules, admet_results):
                    admet_entries.append((mol["id"], run_uuid, "admet", _filter_cols(admet, _inc)))
                    # CNS MPO + Pfizer/GSK + Brenk (piggyback on ADMET run)
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
                _db_sync(_store_properties_batch(admet_entries))
                if druglike_entries:
                    _db_sync(_store_properties_batch(druglike_entries))
                _write_log(run_id, "info", f"ADMET complete: {len(admet_results)} molecules profiled (+ CNS MPO, Pfizer/GSK, Brenk)")
                results["admet"] = f"{len(admet_results)} molecules processed"

            elif calc_type == "scoring":
                from pipeline.scoring import (
                    compute_properties, compute_sa_score,
                )
                physchem_entries = []
                sa_entries = []
                for mol in molecules:
                    try:
                        props = _cached_compute("scoring", mol["smiles"], compute_properties)
                        physchem_entries.append((mol["id"], run_uuid, "physicochemical", _filter_cols(props, _inc)))
                        if _inc is None or "sa_score" in _inc:
                            sa = compute_sa_score(mol["smiles"])
                            if sa is not None:
                                sa_entries.append((mol["id"], run_uuid, "sa_score", {"sa_score": sa}))
                    except Exception as e:
                        logger.warning("Scoring failed for %s: %s", mol["smiles"], e)
                _db_sync(_store_properties_batch(physchem_entries))
                if sa_entries:
                    _db_sync(_store_properties_batch(sa_entries))
                _write_log(run_id, "info", f"Scoring complete: {len(physchem_entries)}/{len(molecules)} molecules scored (+ SA, InChIKey, Ro3)")
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

                _write_log(run_id, "info", f"Downloading receptor from {pdb_url}")
                _db_sync(_update_run(run_id, current_step="Downloading receptor PDB"))

                with tempfile.TemporaryDirectory(prefix="docking_v9_") as tmp_dir:
                    tmp_path = Path(tmp_dir)
                    pdb_path = tmp_path / f"{docking_ctx.get('pdb_id', 'receptor')}.pdb"

                    # Download PDB
                    resp = http_requests.get(pdb_url, timeout=(10, 60))
                    resp.raise_for_status()
                    pdb_path.write_text(resp.text)

                    _write_log(run_id, "info",
                        f"Docking {len(molecules)} molecules — center={pocket_center}, size={pocket_size}")

                    # Prepare ligand list — use molecule UUID as name for collision-free mapping
                    ligands = [{"name": str(m["id"]), "smiles": m["smiles"]}
                               for m in molecules]

                    # Read advanced docking config from run config
                    run_config = run_data.get("config") or {}
                    dock_cfg = run_config.get("docking", {}) if isinstance(run_config.get("docking"), dict) else run_config
                    dock_exhaustiveness = int(dock_cfg.get("exhaustiveness", 8))
                    dock_num_modes = int(dock_cfg.get("num_modes", 9))
                    dock_engine = dock_cfg.get("engine", "gnina_gpu")
                    # Frontend engine names pass through directly to dock_all_ligands
                    # gnina_gpu → GPU only (error if unavailable)
                    # gnina_cpu → local GNINA only
                    # vina      → local Vina only
                    # auto      → GPU → GNINA → Vina → mock

                    # Prepare receptor: convert PDB → PDBQT (required for Vina, accepted by GNINA)
                    from pipeline.prepare import prepare_receptor
                    receptor_pdbqt = prepare_receptor(pdb_path, tmp_path)
                    _write_log(run_id, "info", f"Receptor prepared: {receptor_pdbqt.name}")

                    # Progress callback
                    def progress_cb(pct, msg):
                        _db_sync(_update_run(run_id, current_step=msg, progress=base_pct + int(pct / total_steps)))

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
                    )

                    # Store results as molecule properties (batch insert)
                    # Build id→molecule mapping (str keys to match docking result names)
                    name_to_mol = {str(m["id"]): m for m in molecules}

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
                        # Filter but always keep pose_molblock + docking_engine
                        filtered_dock = _filter_cols(docking_props, _inc)
                        for always_keep in ("pose_molblock", "docking_engine"):
                            if always_keep in docking_props:
                                filtered_dock[always_keep] = docking_props[always_keep]
                        dock_entries.append((mol["id"], run_uuid, "docking", filtered_dock))

                    _db_sync(_store_properties_batch(dock_entries))
                    stored = len(dock_entries)

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

            elif calc_type == "enrichment":
                # Enrichment needs docking results — preserve existing properties
                mol_dicts = []
                for m in molecules:
                    d = {"smiles": m["smiles"], "name": m["name"]}
                    # Carry forward scored properties for composite scoring/Pareto
                    for prop_key in ("affinity", "vina_score", "cnn_score", "cnn_affinity",
                                     "composite_score", "qed", "logP", "source"):
                        if m.get(prop_key) is not None:
                            d[prop_key] = m[prop_key]
                    mol_dicts.append(d)
                try:
                    from pipeline.enrich import enrich_results
                    scored, eliminated, summary = enrich_results(mol_dicts)
                    enrich_entries = [
                        (mol["id"], run_uuid, "enrichment", _filter_cols(enrichment, _inc))
                        for mol, enrichment in zip(molecules, scored)
                    ]
                    _db_sync(_store_properties_batch(enrich_entries))
                    _write_log(run_id, "info", f"Enrichment complete: {len(scored)} enriched, {len(eliminated)} eliminated")
                    results["enrichment"] = f"{len(scored)} molecules enriched ({len(eliminated)} eliminated)"
                except Exception as e:
                    logger.warning("Enrichment failed: %s", e)
                    results["enrichment"] = f"failed: {e}"

            elif calc_type == "clustering":
                from pipeline.scoring import cluster_results
                mol_dicts = [{"smiles": m["smiles"], "name": m["name"]} for m in molecules]
                try:
                    clustered = cluster_results(mol_dicts)
                    n_clusters = len(set(c.get("cluster_id", 0) for c in clustered))
                    cluster_entries = [
                        (mol["id"], run_uuid, "clustering", _filter_cols(
                            {"cluster_id": cl_data.get("cluster_id"),
                             "is_representative": cl_data.get("is_representative")}, _inc))
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
                        ot_props = {"selectivity_score": ot["selectivity_score"],
                                    "off_target_hits": off_target_hits,
                                    "selectivity_ratio": round(selectivity_ratio, 3),
                                    "n_safe": ot["n_safe"],
                                    "n_total": ot["n_total"],
                                    "results": ot["results"],
                                    "warnings": ot["warnings"]}
                        ot_entries.append((mol["id"], run_uuid, "off_target", _filter_cols(ot_props, _inc)))
                    except Exception as e:
                        logger.warning("Off-target failed for %s: %s", mol["smiles"], e)
                _db_sync(_store_properties_batch(ot_entries))
                _write_log(run_id, "info", f"Off-target screening complete: {len(molecules)} molecules vs 10 anti-targets")
                results["off_target"] = f"{len(molecules)} molecules analyzed"

            elif calc_type == "confidence":
                # Multi-component confidence scoring
                from pipeline.confidence import calculate_confidence
                from pipeline.scoring import compute_pains_alert
                conf_entries = []
                for mol in molecules:
                    try:
                        mol_dict = {"smiles": mol["smiles"], "name": mol["name"]}
                        for prop_key in ("affinity", "vina_score", "cnn_score", "admet",
                                         "synthesis_route", "docking_method", "source"):
                            if mol.get(prop_key) is not None:
                                mol_dict[prop_key] = mol[prop_key]
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
                retro_entries = []
                for mol in molecules:
                    try:
                        route = plan_synthesis(mol["smiles"])
                        stored_val = {
                            "n_synth_steps": route.get("n_steps", 0),
                            "synth_confidence": route.get("confidence", 0),
                            "synth_cost_estimate": route.get("estimated_cost"),
                            "reagents_available": route.get("all_reagents_available", True),
                            "n_steps": route.get("n_steps", 0),
                            "confidence": route.get("confidence", 0),
                            "estimated_cost": route.get("estimated_cost"),
                            "all_reagents_available": route.get("all_reagents_available"),
                            "steps": route.get("steps", []),
                            "tree": route.get("tree"),
                            "cost_estimate": route.get("cost_estimate"),
                            "reagent_availability": route.get("reagent_availability", []),
                        }
                        retro_entries.append((mol["id"], run_uuid, "retrosynthesis", _filter_cols(stored_val, _inc)))
                    except Exception as e:
                        logger.warning("Retrosynthesis failed for %s: %s", mol["smiles"], e)
                _db_sync(_store_properties_batch(retro_entries))
                _write_log(run_id, "info", f"Retrosynthesis complete: {len(molecules)} molecules — routes planned")
                results["retrosynthesis"] = f"{len(molecules)} molecules analyzed"

            elif calc_type == "pharmacophore":
                from pipeline.pharmacophore import extract_pharmacophore, compute_pharmacophore_similarity
                pharma_entries = []
                for mol in molecules:
                    try:
                        pharma = extract_pharmacophore(mol["smiles"])
                        pharma_entries.append((mol["id"], run_uuid, "pharmacophore", {
                            "feature_counts": pharma.get("feature_counts", {}),
                            "n_features": pharma.get("n_features", 0),
                            "fingerprint_bits": len(pharma.get("fingerprint", [])),
                        }))
                    except Exception as e:
                        logger.warning("Pharmacophore failed for %s: %s", mol["smiles"], e)
                _db_sync(_store_properties_batch(pharma_entries))

                # Compute pairwise similarity
                try:
                    all_smiles = [m["smiles"] for m in molecules]
                    sim_result = compute_pharmacophore_similarity(all_smiles)
                    if molecules:
                        _db_sync(_store_property(
                            molecules[0]["id"], run_uuid, "pharmacophore_similarity", sim_result
                        ))
                except Exception as e:
                    logger.warning("Pharmacophore similarity failed: %s", e)

                _write_log(run_id, "info", f"Pharmacophore mapping complete: {len(molecules)} molecules analyzed")
                results["pharmacophore"] = f"{len(molecules)} molecules mapped"

            elif calc_type == "activity_cliffs":
                from pipeline.activity_cliffs import detect_activity_cliffs
                # Use docking_score as primary activity metric
                mol_dicts = []
                for m in molecules:
                    d = {"smiles": m["smiles"], "name": m["name"]}
                    for prop_key in ("docking_score", "affinity", "vina_score",
                                     "composite_score", "cnn_score"):
                        if m.get(prop_key) is not None:
                            d[prop_key] = m[prop_key]
                    # Use best available activity metric
                    if "docking_score" not in d and "affinity" in d:
                        d["docking_score"] = d["affinity"]
                    mol_dicts.append(d)
                cliff_results = detect_activity_cliffs(mol_dicts, activity_key="docking_score")
                cliff_entries = []
                n_cliffs = 0
                for mol, cliff in zip(molecules, cliff_results):
                    cliff_entries.append((mol["id"], run_uuid, "activity_cliffs", _filter_cols(cliff, _inc)))
                    if cliff.get("is_cliff"):
                        n_cliffs += 1
                _db_sync(_store_properties_batch(cliff_entries))
                _write_log(run_id, "info",
                    f"Activity cliff detection complete: {n_cliffs}/{len(molecules)} molecules are cliffs")
                results["activity_cliffs"] = f"{n_cliffs} cliffs detected in {len(molecules)} molecules"

            else:
                logger.warning("Unknown calculation type: %s", calc_type)
                results[calc_type] = "unknown type, skipped"

        except Exception as e:
            logger.exception("Calculation %s failed: %s", calc_type, e)
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
