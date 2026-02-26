#!/usr/bin/env python3
"""
PharmacoDB — Resume Ingestion from last checkpoint.
Skips completed phases, resumes compounds at correct offset,
then runs remaining phases (assays, bioactivities, mechanisms, indications).
"""

import sys
import os
import time
import psycopg2
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

DB_CONFIG = {
    "host": "localhost",
    "port": 5433,
    "dbname": "pharmaco",
    "user": "postgres",
    "password": "pharmaco_secret"
}

def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

def main():
    log("=" * 60)
    log("PharmacoDB — RESUME INGESTION")
    log("=" * 60)

    conn = psycopg2.connect(**DB_CONFIG)

    # Check current state
    with conn.cursor() as cur:
        cur.execute("SELECT step, status, rows_inserted FROM pharmaco_db.ingestion_log ORDER BY started_at")
        steps = cur.fetchall()
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.compounds")
        n_compounds = cur.fetchone()[0]
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.targets")
        n_targets = cur.fetchone()[0]

    log(f"Current state: {n_targets} targets, {n_compounds} compounds")
    for step, status, rows in steps:
        log(f"  {step}: {status} ({rows} rows)")

    # Import the ingestion module
    import ingest_chembl

    # --- RESUME COMPOUNDS ---
    # Mark old running entry as failed
    with conn.cursor() as cur:
        cur.execute("""
            UPDATE pharmaco_db.ingestion_log
            SET status='failed', error_message='Interrupted — resuming', completed_at=NOW()
            WHERE step='compounds' AND status='running'
        """)
    conn.commit()

    log(f"\n>>> RESUMING: Compounds from offset {n_compounds}")
    log_id = ingest_chembl.start_log(conn, 'chembl', 'compounds_resume')

    total_inserted = 0
    offset = n_compounds  # Resume from where we left off
    batch = []

    while True:
        url = (f"{ingest_chembl.CHEMBL_API}/molecule.json?"
               f"molecule_properties__mw_freebase__gte=100"
               f"&molecule_properties__mw_freebase__lte=900"
               f"&molecule_structures__canonical_smiles__isnull=false"
               f"&limit={ingest_chembl.BATCH_SIZE}&offset={offset}")
        data = ingest_chembl.api_get(url)
        if not data or not data.get("molecules"):
            break

        molecules = data["molecules"]
        log(f"  Fetched {len(molecules)} compounds (offset={offset}, new={total_inserted})")

        for m in molecules:
            structs = m.get("molecule_structures") or {}
            props = m.get("molecule_properties") or {}
            smiles = structs.get("canonical_smiles")
            if not smiles:
                continue

            batch.append({
                "chembl_id": m.get("molecule_chembl_id"),
                "canonical_smiles": smiles,
                "inchi": structs.get("standard_inchi"),
                "inchi_key": structs.get("standard_inchi_key"),
                "pref_name": m.get("pref_name"),
                "molecular_weight": ingest_chembl._float(props.get("mw_freebase")),
                "alogp": ingest_chembl._float(props.get("alogp")),
                "hba": ingest_chembl._int(props.get("hba")),
                "hbd": ingest_chembl._int(props.get("hbd")),
                "psa": ingest_chembl._float(props.get("psa")),
                "rtb": ingest_chembl._int(props.get("rtb")),
                "num_ro5_violations": ingest_chembl._int(props.get("num_ro5_violations")),
                "aromatic_rings": ingest_chembl._int(props.get("aromatic_rings")),
                "heavy_atoms": ingest_chembl._int(props.get("heavy_atom_count")),
                "qed_weighted": ingest_chembl._float(props.get("qed_weighted")),
                "molecular_formula": props.get("full_molformula"),
                "max_phase": ingest_chembl._int(m.get("max_phase")),
                "is_natural_product": m.get("natural_product") == 1 if m.get("natural_product") is not None else None,
            })

        if len(batch) >= ingest_chembl.INSERT_BATCH:
            n = ingest_chembl._insert_compounds(conn, batch)
            total_inserted += n
            batch = []
            if total_inserted % 10000 == 0:
                log(f"  Total new compounds: {total_inserted}")

        if not data.get("page_meta", {}).get("next"):
            break
        offset += ingest_chembl.BATCH_SIZE
        time.sleep(0.2)

    if batch:
        n = ingest_chembl._insert_compounds(conn, batch)
        total_inserted += n

    log(f"  DONE: {total_inserted} new compounds inserted (total now: {n_compounds + total_inserted})")
    ingest_chembl.update_log(conn, log_id, 'completed', total_inserted)

    # --- ASSAYS ---
    log("\n>>> STARTING: Assays")
    ingest_chembl.ingest_assays(conn)

    # --- BIOACTIVITIES ---
    log("\n>>> STARTING: Bioactivities (this is the big one)")
    ingest_chembl.ingest_bioactivities(conn)

    # --- MECHANISMS ---
    log("\n>>> STARTING: Drug Mechanisms")
    ingest_chembl.ingest_mechanisms(conn)

    # --- INDICATIONS ---
    log("\n>>> STARTING: Drug Indications")
    ingest_chembl.ingest_indications(conn)

    # --- FINAL STATS ---
    with conn.cursor() as cur:
        for table in ['compounds', 'targets', 'assays', 'bioactivities', 'drug_mechanisms', 'drug_indications']:
            cur.execute(f"SELECT COUNT(*) FROM pharmaco_db.{table}")
            count = cur.fetchone()[0]
            log(f"  {table}: {count:,}")

        cur.execute("SELECT pg_size_pretty(pg_database_size('pharmaco'))")
        size = cur.fetchone()[0]
        log(f"  DB size: {size}")

    conn.close()
    log("\n" + "=" * 60)
    log("CHEMBL INGESTION COMPLETE — Ready for Phase 2 (enrichment)")
    log("=" * 60)

if __name__ == "__main__":
    main()
