#!/usr/bin/env python3
"""
PharmacoDB — TURBO Ingestion (10x faster)
Parallel API fetching + bulk DB inserts for compounds, assays, bioactivities.
Uses ThreadPoolExecutor with 8 workers to fetch pages in parallel.
"""

import json
import time
import sys
import os
import urllib.request
import urllib.error
import psycopg2
import psycopg2.extras
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
DB_CONFIG = {
    "host": "localhost", "port": 5433,
    "dbname": "pharmaco", "user": "postgres", "password": "pharmaco_secret"
}

# Turbo settings
FETCH_WORKERS = 8       # Parallel API fetchers
PAGE_SIZE = 1000        # ChEMBL supports up to 1000
INSERT_BATCH = 5000     # Bulk insert size
MAX_RETRIES = 5
RETRY_DELAY = 3

def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

def _float(v):
    try: return float(v) if v is not None else None
    except: return None

def _int(v):
    try: return int(v) if v is not None else None
    except: return None

# ============================================================
# Fast parallel API fetcher
# ============================================================
def api_get(url, retries=MAX_RETRIES):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={
                "Accept": "application/json",
                "User-Agent": "PharmacoDB-Turbo/2.0"
            })
            with urllib.request.urlopen(req, timeout=90) as resp:
                return json.loads(resp.read().decode())
        except urllib.error.HTTPError as e:
            if e.code == 429:
                wait = RETRY_DELAY * (attempt + 1) * 2
                log(f"  Rate limited, waiting {wait}s...")
                time.sleep(wait)
            elif e.code >= 500:
                time.sleep(RETRY_DELAY * (attempt + 1))
            else:
                log(f"  HTTP {e.code} for {url}")
                return None
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(RETRY_DELAY)
    return None

def fetch_page(base_url, offset):
    """Fetch a single page. Returns (offset, data)."""
    url = f"{base_url}&limit={PAGE_SIZE}&offset={offset}"
    data = api_get(url)
    return (offset, data)

def parallel_fetch(base_url, start_offset, count_key, callback, max_pages=None):
    """
    Fetch pages in parallel using a sliding window.
    callback(items) is called with each page's items.
    Returns total items fetched.
    """
    total = 0
    offset = start_offset
    done = False
    pages_fetched = 0

    while not done:
        # Submit FETCH_WORKERS pages in parallel
        futures = {}
        with ThreadPoolExecutor(max_workers=FETCH_WORKERS) as pool:
            for i in range(FETCH_WORKERS):
                f = pool.submit(fetch_page, base_url, offset + i * PAGE_SIZE)
                futures[f] = offset + i * PAGE_SIZE

            # Collect results in order
            results = {}
            for f in as_completed(futures):
                page_offset = futures[f]
                try:
                    _, data = f.result()
                    results[page_offset] = data
                except Exception as e:
                    log(f"  Error at offset {page_offset}: {e}")
                    results[page_offset] = None

        # Process in order
        for i in range(FETCH_WORKERS):
            page_offset = offset + i * PAGE_SIZE
            data = results.get(page_offset)

            if not data or not data.get(count_key):
                done = True
                break

            items = data[count_key]
            callback(items)
            total += len(items)
            pages_fetched += 1

            if max_pages and pages_fetched >= max_pages:
                done = True
                break

            if not data.get("page_meta", {}).get("next"):
                done = True
                break

        offset += FETCH_WORKERS * PAGE_SIZE

        # Brief pause between batches to be nice
        time.sleep(0.1)

    return total

# ============================================================
# DB helpers
# ============================================================
def get_conn():
    return psycopg2.connect(**DB_CONFIG)

def start_log(conn, source, step):
    with conn.cursor() as cur:
        cur.execute("""
            INSERT INTO pharmaco_db.ingestion_log (source, step, status, started_at)
            VALUES (%s, %s, 'running', NOW()) RETURNING id
        """, (source, step))
        log_id = cur.fetchone()[0]
    conn.commit()
    return log_id

def update_log(conn, log_id, status, rows):
    with conn.cursor() as cur:
        cur.execute("""
            UPDATE pharmaco_db.ingestion_log
            SET status=%s, rows_inserted=%s, completed_at=NOW()
            WHERE id=%s
        """, (status, rows, log_id))
    conn.commit()

def bulk_insert(conn, table, columns, rows):
    """Fast bulk insert using execute_values."""
    if not rows:
        return 0
    col_list = ", ".join(columns)
    template = "(" + ", ".join(["%s"] * len(columns)) + ")"
    sql = f"""
        INSERT INTO pharmaco_db.{table} ({col_list})
        VALUES %s
        ON CONFLICT DO NOTHING
    """
    try:
        with conn.cursor() as cur:
            psycopg2.extras.execute_values(cur, sql, rows, template=template, page_size=2000)
        conn.commit()
        return len(rows)
    except Exception as e:
        conn.rollback()
        log(f"  Insert error on {table}: {e}")
        return 0

# ============================================================
# COMPOUNDS (turbo)
# ============================================================
def turbo_compounds(conn):
    log("=" * 60)
    log("TURBO: Compounds")
    log("=" * 60)

    with conn.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.compounds")
        existing = cur.fetchone()[0]
    log(f"  Existing: {existing:,}")

    log_id = start_log(conn, 'chembl', 'compounds_turbo')
    base_url = (f"{CHEMBL_API}/molecule.json?"
                f"molecule_properties__mw_freebase__gte=100"
                f"&molecule_properties__mw_freebase__lte=900"
                f"&molecule_structures__canonical_smiles__isnull=false")

    batch = []
    total_inserted = 0
    t0 = time.time()

    def process_molecules(molecules):
        nonlocal batch, total_inserted
        for m in molecules:
            structs = m.get("molecule_structures") or {}
            props = m.get("molecule_properties") or {}
            smiles = structs.get("canonical_smiles")
            if not smiles:
                continue
            batch.append((
                m.get("molecule_chembl_id"),
                smiles,
                structs.get("standard_inchi"),
                structs.get("standard_inchi_key"),
                m.get("pref_name"),
                _float(props.get("mw_freebase")),
                _float(props.get("alogp")),
                _int(props.get("hba")),
                _int(props.get("hbd")),
                _float(props.get("psa")),
                _int(props.get("rtb")),
                _int(props.get("num_ro5_violations")),
                _int(props.get("aromatic_rings")),
                _int(props.get("heavy_atom_count")),
                _float(props.get("qed_weighted")),
                props.get("full_molformula"),
                _int(m.get("max_phase")),
                m.get("natural_product") == 1 if m.get("natural_product") is not None else None,
            ))

        if len(batch) >= INSERT_BATCH:
            n = bulk_insert(conn, "compounds", [
                "chembl_id", "canonical_smiles", "inchi", "inchi_key", "pref_name",
                "molecular_weight", "alogp", "hba", "hbd", "psa", "rtb",
                "num_ro5_violations", "aromatic_rings", "heavy_atoms", "qed_weighted",
                "molecular_formula", "max_phase", "is_natural_product",
            ], batch)
            total_inserted += n
            batch.clear()
            elapsed = time.time() - t0
            rate = total_inserted / elapsed if elapsed > 0 else 0
            log(f"  Inserted {total_inserted:,} new ({rate:.0f}/s)")

    total_fetched = parallel_fetch(base_url, existing, "molecules", process_molecules)

    # Flush remaining
    if batch:
        n = bulk_insert(conn, "compounds", [
            "chembl_id", "canonical_smiles", "inchi", "inchi_key", "pref_name",
            "molecular_weight", "alogp", "hba", "hbd", "psa", "rtb",
            "num_ro5_violations", "aromatic_rings", "heavy_atoms", "qed_weighted",
            "molecular_formula", "max_phase", "is_natural_product",
        ], batch)
        total_inserted += n

    elapsed = time.time() - t0
    log(f"  DONE: {total_inserted:,} compounds in {int(elapsed)}s ({total_inserted/elapsed:.0f}/s)")
    update_log(conn, log_id, 'completed', total_inserted)
    return total_inserted

# ============================================================
# ASSAYS (turbo)
# ============================================================
def turbo_assays(conn):
    log("=" * 60)
    log("TURBO: Assays")
    log("=" * 60)

    with conn.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.assays")
        existing = cur.fetchone()[0]
    if existing > 0:
        log(f"  Already have {existing:,} assays, skipping.")
        return 0

    # Target map
    target_map = {}
    with conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.targets")
        target_map = dict(cur.fetchall())
    log(f"  Target map: {len(target_map):,} entries")

    log_id = start_log(conn, 'chembl', 'assays_turbo')
    total_inserted = 0
    t0 = time.time()

    for assay_type in ["B", "F"]:
        log(f"  Fetching {assay_type}-type assays...")
        base_url = (f"{CHEMBL_API}/assay.json?"
                    f"assay_type={assay_type}"
                    f"&confidence_score__gte=7"
                    f"&assay_organism=Homo%20sapiens")

        batch = []

        def process_assays(assays):
            nonlocal batch, total_inserted
            for a in assays:
                target_chembl = a.get("target_chembl_id")
                src = a.get("src_id")
                batch.append((
                    a.get("assay_chembl_id"),
                    (a.get("description") or "")[:500],
                    a.get("assay_type"),
                    a.get("assay_category"),
                    _int(src) if src else None,
                    None,  # src_description
                    None,  # journal
                    None,  # year
                    None,  # doi
                    target_map.get(target_chembl),
                    target_chembl,
                    a.get("assay_organism"),
                    a.get("assay_cell_type"),
                    _int(a.get("confidence_score")),
                ))

            if len(batch) >= INSERT_BATCH:
                n = bulk_insert(conn, "assays", [
                    "chembl_id", "description", "assay_type", "assay_category",
                    "src_id", "src_description", "journal", "year", "doi",
                    "target_id", "target_chembl_id",
                    "assay_organism", "assay_cell_type", "confidence_score",
                ], batch)
                total_inserted += n
                batch.clear()
                elapsed = time.time() - t0
                log(f"  Assays inserted: {total_inserted:,} ({total_inserted/elapsed:.0f}/s)")

        parallel_fetch(base_url, 0, "assays", process_assays)

        # Flush
        if batch:
            n = bulk_insert(conn, "assays", [
                "chembl_id", "description", "assay_type", "assay_category",
                "src_id", "src_description", "journal", "year", "doi",
                "target_id", "target_chembl_id",
                "assay_organism", "assay_cell_type", "confidence_score",
            ], batch)
            total_inserted += n
            batch = []

    elapsed = time.time() - t0
    log(f"  DONE: {total_inserted:,} assays in {int(elapsed)}s")
    update_log(conn, log_id, 'completed', total_inserted)
    return total_inserted

# ============================================================
# BIOACTIVITIES (turbo) — the big one
# ============================================================
def turbo_bioactivities(conn):
    log("=" * 60)
    log("TURBO: Bioactivities (the big one)")
    log("=" * 60)

    with conn.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.bioactivities")
        existing = cur.fetchone()[0]
    if existing > 0:
        log(f"  Already have {existing:,}, resuming from offset {existing}")

    # Build lookup maps
    log("  Building lookup maps...")
    compound_map = {}
    with conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.compounds")
        compound_map = dict(cur.fetchall())
    log(f"  Compound map: {len(compound_map):,}")

    target_map = {}
    with conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.targets")
        target_map = dict(cur.fetchall())
    log(f"  Target map: {len(target_map):,}")

    assay_map = {}
    with conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.assays")
        assay_map = dict(cur.fetchall())
    log(f"  Assay map: {len(assay_map):,}")

    log_id = start_log(conn, 'chembl', 'bioactivities_turbo')
    total_inserted = 0
    total_skipped = 0
    t0 = time.time()

    base_url = (f"{CHEMBL_API}/activity.json?"
                f"pchembl_value__isnull=false"
                f"&target_organism=Homo%20sapiens")

    batch = []

    def process_activities(activities):
        nonlocal batch, total_inserted, total_skipped
        for act in activities:
            mol_chembl = act.get("molecule_chembl_id")
            compound_id = compound_map.get(mol_chembl)
            if not compound_id:
                total_skipped += 1
                continue

            target_chembl = act.get("target_chembl_id")
            assay_chembl = act.get("assay_chembl_id")

            pchembl = _float(act.get("pchembl_value"))
            value = _float(act.get("value"))
            units = act.get("units")

            # Classify activity
            activity_class = None
            if pchembl is not None:
                if pchembl >= 7:
                    activity_class = "active"
                elif pchembl >= 5:
                    activity_class = "intermediate"
                else:
                    activity_class = "inactive"

            batch.append((
                compound_id,
                target_map.get(target_chembl),
                assay_map.get(assay_chembl),
                act.get("standard_type"),
                act.get("standard_relation"),
                value,
                units,
                pchembl,
                activity_class,
                "ChEMBL",
                act.get("activity_id"),
                act.get("data_validity_comment"),
                act.get("potential_duplicate"),
            ))

        if len(batch) >= INSERT_BATCH:
            n = bulk_insert(conn, "bioactivities", [
                "compound_id", "target_id", "assay_id",
                "activity_type", "relation", "value", "units", "pchembl_value",
                "activity_class", "source", "source_id",
                "data_validity", "potential_duplicate",
            ], batch)
            total_inserted += n
            batch.clear()
            elapsed = time.time() - t0
            rate = total_inserted / elapsed if elapsed > 0 else 0
            eta_s = int((20_000_000 - total_inserted) / rate) if rate > 0 else 0
            log(f"  Bio: {total_inserted:,} inserted, {total_skipped:,} skipped ({rate:.0f}/s, ETA ~{eta_s//3600}h{(eta_s%3600)//60:02d}m)")

    total_fetched = parallel_fetch(base_url, existing, "activities", process_activities)

    # Flush
    if batch:
        n = bulk_insert(conn, "bioactivities", [
            "compound_id", "target_id", "assay_id",
            "activity_type", "relation", "value", "units", "pchembl_value",
            "activity_class", "source", "source_id",
            "data_validity", "potential_duplicate",
        ], batch)
        total_inserted += n

    elapsed = time.time() - t0
    log(f"  DONE: {total_inserted:,} bioactivities in {int(elapsed)}s ({total_inserted/elapsed:.0f}/s)")
    update_log(conn, log_id, 'completed', total_inserted)
    return total_inserted

# ============================================================
# MECHANISMS (turbo)
# ============================================================
def turbo_mechanisms(conn):
    log("=" * 60)
    log("TURBO: Drug Mechanisms")
    log("=" * 60)

    with conn.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.drug_mechanisms")
        if cur.fetchone()[0] > 0:
            log("  Already populated, skipping.")
            return 0

    import ingest_chembl
    compound_map = {}
    with conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.compounds")
        compound_map = dict(cur.fetchall())

    target_map = {}
    with conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.targets")
        target_map = dict(cur.fetchall())

    log_id = start_log(conn, 'chembl', 'mechanisms_turbo')
    total_inserted = 0
    t0 = time.time()
    batch = []

    base_url = f"{CHEMBL_API}/mechanism.json?"

    def process_mechs(mechs):
        nonlocal batch, total_inserted
        for m in mechs:
            mol_chembl = m.get("molecule_chembl_id")
            target_chembl = m.get("target_chembl_id")
            batch.append((
                compound_map.get(mol_chembl),
                target_map.get(target_chembl),
                m.get("mechanism_of_action"),
                m.get("action_type"),
                m.get("direct_interaction"),
                m.get("molecular_mechanism"),
                m.get("selectivity_comment"),
            ))

        if len(batch) >= INSERT_BATCH:
            n = bulk_insert(conn, "drug_mechanisms", [
                "compound_id", "target_id", "mechanism_of_action",
                "action_type", "direct_interaction", "molecular_mechanism",
                "selectivity_comment",
            ], batch)
            total_inserted += n
            batch.clear()

    parallel_fetch(base_url, 0, "mechanisms", process_mechs)

    if batch:
        n = bulk_insert(conn, "drug_mechanisms", [
            "compound_id", "target_id", "mechanism_of_action",
            "action_type", "direct_interaction", "molecular_mechanism",
            "selectivity_comment",
        ], batch)
        total_inserted += n

    elapsed = time.time() - t0
    log(f"  DONE: {total_inserted:,} mechanisms in {int(elapsed)}s")
    update_log(conn, log_id, 'completed', total_inserted)
    return total_inserted

# ============================================================
# INDICATIONS (turbo)
# ============================================================
def turbo_indications(conn):
    log("=" * 60)
    log("TURBO: Drug Indications")
    log("=" * 60)

    with conn.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.drug_indications")
        if cur.fetchone()[0] > 0:
            log("  Already populated, skipping.")
            return 0

    compound_map = {}
    with conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.compounds")
        compound_map = dict(cur.fetchall())

    log_id = start_log(conn, 'chembl', 'indications_turbo')
    total_inserted = 0
    t0 = time.time()
    batch = []

    base_url = f"{CHEMBL_API}/drug_indication.json?"

    def process_indications(indications):
        nonlocal batch, total_inserted
        for ind in indications:
            mol_chembl = ind.get("molecule_chembl_id")
            refs = ind.get("indication_refs") or []
            ref_str = "; ".join([f"{r.get('ref_type','')}: {r.get('ref_id','')}" for r in refs]) if refs else None

            batch.append((
                compound_map.get(mol_chembl),
                ind.get("mesh_id"),
                ind.get("mesh_heading"),
                ind.get("efo_id"),
                ind.get("efo_term"),
                _int(ind.get("max_phase_for_ind")),
                ref_str,
            ))

        if len(batch) >= INSERT_BATCH:
            n = bulk_insert(conn, "drug_indications", [
                "compound_id", "mesh_id", "mesh_heading",
                "efo_id", "efo_term", "max_phase", "indication_refs",
            ], batch)
            total_inserted += n
            batch.clear()

    parallel_fetch(base_url, 0, "drug_indications", process_indications)

    if batch:
        n = bulk_insert(conn, "drug_indications", [
            "compound_id", "mesh_id", "mesh_heading",
            "efo_id", "efo_term", "max_phase", "indication_refs",
        ], batch)
        total_inserted += n

    elapsed = time.time() - t0
    log(f"  DONE: {total_inserted:,} indications in {int(elapsed)}s")
    update_log(conn, log_id, 'completed', total_inserted)
    return total_inserted

# ============================================================
# MAIN
# ============================================================
def main():
    log("=" * 60)
    log("PharmacoDB — TURBO INGESTION (8 parallel workers)")
    log(f"  Page size: {PAGE_SIZE}, Insert batch: {INSERT_BATCH}")
    log("=" * 60)

    conn = get_conn()
    t0 = time.time()

    # 1. Finish compounds
    turbo_compounds(conn)

    # 2. Assays
    turbo_assays(conn)

    # 3. Bioactivities (the big one)
    turbo_bioactivities(conn)

    # 4. Mechanisms
    turbo_mechanisms(conn)

    # 5. Indications
    turbo_indications(conn)

    # Final stats
    log("\n" + "=" * 60)
    log("FINAL STATS:")
    with conn.cursor() as cur:
        for t in ['compounds', 'targets', 'assays', 'bioactivities', 'drug_mechanisms', 'drug_indications']:
            cur.execute(f"SELECT COUNT(*) FROM pharmaco_db.{t}")
            log(f"  {t:>20}: {cur.fetchone()[0]:>12,}")
        cur.execute("SELECT pg_size_pretty(pg_database_size('pharmaco'))")
        log(f"  {'DB size':>20}: {cur.fetchone()[0]}")

    elapsed = time.time() - t0
    log(f"\nTOTAL TIME: {int(elapsed)}s ({elapsed/3600:.1f}h)")
    log("=" * 60)
    conn.close()

if __name__ == "__main__":
    main()
