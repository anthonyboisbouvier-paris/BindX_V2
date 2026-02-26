#!/usr/bin/env python3
"""
PharmacoDB — Master Ingestion Orchestrator
Runs all ingestion pipelines in order:
1. ChEMBL (targets, compounds, assays, activities, mechanisms, indications)
2. Target classification enrichment (ChEMBL + NCBI)
3. UniProt enrichment (function, GO, pathways, PDB)
4. PubChem (CID resolution + bioassays)
5. Final statistics
"""

import sys
import os
import time
import psycopg2
from datetime import datetime

# Add current dir to path
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

def get_stats(conn):
    """Print comprehensive database statistics."""
    stats = {}
    with conn.cursor() as cur:
        tables = [
            ("compounds", "pharmaco_db.compounds"),
            ("targets", "pharmaco_db.targets"),
            ("assays", "pharmaco_db.assays"),
            ("bioactivities", "pharmaco_db.bioactivities"),
            ("drug_mechanisms", "pharmaco_db.drug_mechanisms"),
            ("drug_indications", "pharmaco_db.drug_indications"),
            ("cross_references", "pharmaco_db.cross_references"),
        ]
        for name, table in tables:
            cur.execute(f"SELECT COUNT(*) FROM {table}")
            stats[name] = cur.fetchone()[0]

        # Extra stats
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.targets WHERE uniprot_function IS NOT NULL")
        stats["targets_with_function"] = cur.fetchone()[0]

        cur.execute("SELECT COUNT(*) FROM pharmaco_db.targets WHERE protein_family IS NOT NULL")
        stats["targets_with_family"] = cur.fetchone()[0]

        cur.execute("SELECT COUNT(*) FROM pharmaco_db.targets WHERE is_druggable = TRUE")
        stats["druggable_targets"] = cur.fetchone()[0]

        cur.execute("SELECT COUNT(*) FROM pharmaco_db.compounds WHERE pubchem_cid IS NOT NULL")
        stats["compounds_with_cid"] = cur.fetchone()[0]

        cur.execute("SELECT COUNT(*) FROM pharmaco_db.compounds WHERE max_phase >= 4")
        stats["approved_drugs"] = cur.fetchone()[0]

        cur.execute("SELECT COUNT(DISTINCT compound_id) FROM pharmaco_db.bioactivities")
        stats["compounds_with_activity"] = cur.fetchone()[0]

        cur.execute("SELECT COUNT(DISTINCT target_id) FROM pharmaco_db.bioactivities WHERE target_id IS NOT NULL")
        stats["targets_with_activity"] = cur.fetchone()[0]

        cur.execute("""
            SELECT activity_type, COUNT(*) FROM pharmaco_db.bioactivities
            GROUP BY activity_type ORDER BY COUNT(*) DESC LIMIT 10
        """)
        stats["top_activity_types"] = cur.fetchall()

        cur.execute("""
            SELECT protein_family, COUNT(*) FROM pharmaco_db.targets
            WHERE protein_family IS NOT NULL
            GROUP BY protein_family ORDER BY COUNT(*) DESC LIMIT 10
        """)
        stats["top_families"] = cur.fetchall()

        # DB size
        cur.execute("""
            SELECT pg_size_pretty(pg_total_relation_size('pharmaco_db.compounds')) as compounds_size,
                   pg_size_pretty(pg_total_relation_size('pharmaco_db.bioactivities')) as activities_size,
                   pg_size_pretty(pg_database_size('pharmaco')) as total_size
        """)
        sizes = cur.fetchone()
        stats["size_compounds"] = sizes[0]
        stats["size_activities"] = sizes[1]
        stats["size_total"] = sizes[2]

    return stats

def print_stats(stats):
    log("=" * 60)
    log("PHARMACO DB — FINAL STATISTICS")
    log("=" * 60)
    log(f"  Compounds:          {stats['compounds']:>10,}")
    log(f"    with PubChem CID: {stats['compounds_with_cid']:>10,}")
    log(f"    approved drugs:   {stats['approved_drugs']:>10,}")
    log(f"    with activity:    {stats['compounds_with_activity']:>10,}")
    log(f"  Targets:            {stats['targets']:>10,}")
    log(f"    with function:    {stats['targets_with_function']:>10,}")
    log(f"    with family:      {stats['targets_with_family']:>10,}")
    log(f"    druggable:        {stats['druggable_targets']:>10,}")
    log(f"    with activity:    {stats['targets_with_activity']:>10,}")
    log(f"  Assays:             {stats['assays']:>10,}")
    log(f"  Bioactivities:      {stats['bioactivities']:>10,}")
    log(f"  Drug Mechanisms:    {stats['drug_mechanisms']:>10,}")
    log(f"  Drug Indications:   {stats['drug_indications']:>10,}")
    log(f"  Cross References:   {stats['cross_references']:>10,}")
    log("")
    log("  Top activity types:")
    for atype, count in stats.get("top_activity_types", []):
        log(f"    {atype:>15}: {count:>10,}")
    log("")
    log("  Top target families:")
    for family, count in stats.get("top_families", []):
        log(f"    {family:>25}: {count:>6,}")
    log("")
    log(f"  DB Size: {stats['size_total']}")
    log(f"    Compounds table:    {stats['size_compounds']}")
    log(f"    Activities table:   {stats['size_activities']}")
    log("=" * 60)

def main():
    start_time = time.time()
    log("=" * 60)
    log("PharmacoDB — MASTER INGESTION PIPELINE")
    log(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log("=" * 60)

    # ---- PHASE 1: ChEMBL Core Data ----
    log("\n>>> STARTING: ChEMBL Core Data Ingestion")
    import ingest_chembl
    ingest_chembl.main()

    # ---- PHASE 2: Target Classification ----
    log("\n>>> STARTING: Target Classification Enrichment")
    import ingest_target_classes
    ingest_target_classes.main()

    # ---- PHASE 3: UniProt Enrichment ----
    log("\n>>> STARTING: UniProt Target Enrichment")
    import ingest_uniprot
    ingest_uniprot.main()

    # ---- PHASE 4: PubChem Data ----
    log("\n>>> STARTING: PubChem Data Ingestion")
    import ingest_pubchem
    ingest_pubchem.main()

    # ---- FINAL: Statistics ----
    conn = psycopg2.connect(**DB_CONFIG)
    try:
        stats = get_stats(conn)
        print_stats(stats)
    finally:
        conn.close()

    elapsed = time.time() - start_time
    hours = int(elapsed // 3600)
    minutes = int((elapsed % 3600) // 60)
    log(f"\nTotal time: {hours}h {minutes}m")
    log("INGESTION COMPLETE")

if __name__ == "__main__":
    main()
