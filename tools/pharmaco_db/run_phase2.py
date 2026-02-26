#!/usr/bin/env python3
"""
PharmacoDB — Phase 2 Pipeline (runs AFTER ChEMBL core finishes)
1. Apply extended tag schema
2. Target classification + NCBI enrichment
3. UniProt enrichment
4. PubChem CIDs + bioassays
5. Extra sources (ATC, GtoPdb, Open Targets, BindingDB)
6. Compute ALL tags
7. Final statistics
"""

import sys
import os
import time
import subprocess
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

def get_stats(conn):
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
            ("compound_tags", "pharmaco_db.compound_tags"),
            ("target_tags", "pharmaco_db.target_tags"),
            ("activity_tags", "pharmaco_db.activity_tags"),
            ("atc_classification", "pharmaco_db.atc_classification"),
            ("gtop_interactions", "pharmaco_db.gtop_interactions"),
            ("disease_associations", "pharmaco_db.disease_target_associations"),
            ("selectivity_profiles", "pharmaco_db.selectivity_profiles"),
        ]
        for name, table in tables:
            try:
                cur.execute(f"SELECT COUNT(*) FROM {table}")
                stats[name] = cur.fetchone()[0]
            except:
                conn.rollback()
                stats[name] = 0

        # Enrichment stats
        try:
            cur.execute("SELECT COUNT(*) FROM pharmaco_db.targets WHERE uniprot_function IS NOT NULL")
            stats["targets_enriched_uniprot"] = cur.fetchone()[0]
        except:
            conn.rollback()
            stats["targets_enriched_uniprot"] = 0

        try:
            cur.execute("SELECT COUNT(*) FROM pharmaco_db.targets WHERE protein_family IS NOT NULL")
            stats["targets_with_family"] = cur.fetchone()[0]
        except:
            conn.rollback()
            stats["targets_with_family"] = 0

        try:
            cur.execute("SELECT COUNT(*) FROM pharmaco_db.compounds WHERE pubchem_cid IS NOT NULL")
            stats["compounds_with_cid"] = cur.fetchone()[0]
        except:
            conn.rollback()
            stats["compounds_with_cid"] = 0

        try:
            cur.execute("SELECT COUNT(*) FROM pharmaco_db.compounds WHERE max_phase >= 4")
            stats["approved_drugs"] = cur.fetchone()[0]
        except:
            conn.rollback()
            stats["approved_drugs"] = 0

        try:
            cur.execute("SELECT COUNT(*) FROM pharmaco_db.compound_tags WHERE lipinski_pass = TRUE")
            stats["lipinski_pass"] = cur.fetchone()[0]
        except:
            conn.rollback()
            stats["lipinski_pass"] = 0

        try:
            cur.execute("SELECT COUNT(*) FROM pharmaco_db.compound_tags WHERE drug_like = TRUE")
            stats["drug_like"] = cur.fetchone()[0]
        except:
            conn.rollback()
            stats["drug_like"] = 0

        try:
            cur.execute("SELECT COUNT(*) FROM pharmaco_db.compound_tags WHERE cns_penetrant = TRUE")
            stats["cns_penetrant"] = cur.fetchone()[0]
        except:
            conn.rollback()
            stats["cns_penetrant"] = 0

        try:
            cur.execute("""
                SELECT source, COUNT(*) FROM pharmaco_db.bioactivities
                GROUP BY source ORDER BY COUNT(*) DESC
            """)
            stats["activities_by_source"] = cur.fetchall()
        except:
            conn.rollback()
            stats["activities_by_source"] = []

        try:
            cur.execute("""
                SELECT target_class, COUNT(*) FROM pharmaco_db.target_tags
                WHERE target_class IS NOT NULL
                GROUP BY target_class ORDER BY COUNT(*) DESC LIMIT 15
            """)
            stats["target_classes"] = cur.fetchall()
        except:
            conn.rollback()
            stats["target_classes"] = []

        try:
            cur.execute("""
                SELECT pchembl_bin, COUNT(*) FROM pharmaco_db.activity_tags
                WHERE pchembl_bin IS NOT NULL
                GROUP BY pchembl_bin ORDER BY pchembl_bin
            """)
            stats["activity_bins"] = cur.fetchall()
        except:
            conn.rollback()
            stats["activity_bins"] = []

        try:
            cur.execute("SELECT pg_size_pretty(pg_database_size('pharmaco'))")
            stats["db_size"] = cur.fetchone()[0]
        except:
            conn.rollback()
            stats["db_size"] = "unknown"

    return stats

def print_stats(stats):
    log("=" * 70)
    log("  PHARMACO DB — FINAL COMPREHENSIVE STATISTICS")
    log("=" * 70)
    log("")
    log("  === CORE DATA ===")
    log(f"    Compounds:              {stats.get('compounds', 0):>12,}")
    log(f"      Approved drugs:       {stats.get('approved_drugs', 0):>12,}")
    log(f"      With PubChem CID:     {stats.get('compounds_with_cid', 0):>12,}")
    log(f"    Targets:                {stats.get('targets', 0):>12,}")
    log(f"      Enriched (UniProt):   {stats.get('targets_enriched_uniprot', 0):>12,}")
    log(f"      With protein family:  {stats.get('targets_with_family', 0):>12,}")
    log(f"    Assays:                 {stats.get('assays', 0):>12,}")
    log(f"    Bioactivities:          {stats.get('bioactivities', 0):>12,}")
    log(f"    Drug Mechanisms:        {stats.get('drug_mechanisms', 0):>12,}")
    log(f"    Drug Indications:       {stats.get('drug_indications', 0):>12,}")
    log(f"    Cross References:       {stats.get('cross_references', 0):>12,}")
    log("")
    log("  === ADDITIONAL SOURCES ===")
    log(f"    ATC Classifications:    {stats.get('atc_classification', 0):>12,}")
    log(f"    GtoPdb Interactions:    {stats.get('gtop_interactions', 0):>12,}")
    log(f"    Disease Associations:   {stats.get('disease_associations', 0):>12,}")
    log(f"    Selectivity Profiles:   {stats.get('selectivity_profiles', 0):>12,}")
    log("")
    log("  === ML TAGS ===")
    log(f"    Compound tags:          {stats.get('compound_tags', 0):>12,}")
    log(f"      Lipinski pass:        {stats.get('lipinski_pass', 0):>12,}")
    log(f"      Drug-like (QED>0.5):  {stats.get('drug_like', 0):>12,}")
    log(f"      CNS penetrant:        {stats.get('cns_penetrant', 0):>12,}")
    log(f"    Target tags:            {stats.get('target_tags', 0):>12,}")
    log(f"    Activity tags:          {stats.get('activity_tags', 0):>12,}")
    log("")

    if stats.get("activities_by_source"):
        log("  === ACTIVITIES BY SOURCE ===")
        for src, cnt in stats["activities_by_source"]:
            log(f"    {src:>20}: {cnt:>12,}")
        log("")

    if stats.get("target_classes"):
        log("  === TARGET CLASSES ===")
        for cls, cnt in stats["target_classes"]:
            log(f"    {cls:>25}: {cnt:>8,}")
        log("")

    if stats.get("activity_bins"):
        log("  === ACTIVITY POTENCY DISTRIBUTION ===")
        for bin_name, cnt in stats["activity_bins"]:
            log(f"    {bin_name:>20}: {cnt:>12,}")
        log("")

    log(f"  Database size: {stats.get('db_size', 'unknown')}")
    log("=" * 70)

def main():
    start_time = time.time()
    log("=" * 70)
    log("PharmacoDB — PHASE 2: Enrichment + Tags + Extra Sources")
    log(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log("=" * 70)

    # ---- Step 1: Target Classification ----
    log("\n>>> STEP 1: Target Classification")
    import ingest_target_classes
    ingest_target_classes.main()

    # ---- Step 2: UniProt Enrichment ----
    log("\n>>> STEP 2: UniProt Enrichment")
    import ingest_uniprot
    ingest_uniprot.main()

    # ---- Step 3: PubChem Data ----
    log("\n>>> STEP 3: PubChem Data")
    import ingest_pubchem
    ingest_pubchem.main()

    # ---- Step 4: Extra Sources ----
    log("\n>>> STEP 4: Extra Sources (ATC, GtoPdb, Open Targets, BindingDB)")
    import ingest_extra_sources
    ingest_extra_sources.main()

    # ---- Step 5: Compute ALL Tags ----
    log("\n>>> STEP 5: Computing Comprehensive Tags")
    import compute_tags
    compute_tags.main()

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
    log(f"\nPhase 2 total time: {hours}h {minutes}m")
    log("PHASE 2 COMPLETE — Database ready for ML/SAR analysis")

if __name__ == "__main__":
    main()
