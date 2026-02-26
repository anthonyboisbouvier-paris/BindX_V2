#!/usr/bin/env python3
"""
Fast push to Supabase: DROP indexes → bulk INSERT → RECREATE indexes.
Uses COPY-like speed via execute_values with large pages.
"""
import psycopg2
import psycopg2.extras
import time
import sys
from datetime import datetime

LOCAL_DB = {
    "host": "localhost", "port": 5433,
    "dbname": "pharmaco", "user": "postgres", "password": "postgres"
}
SUPABASE_CONN = "postgresql://postgres.webkntghfzscrnuixfba:n1qsKxvAWJV954vR@aws-0-us-west-2.pooler.supabase.com:6543/postgres"

BATCH = 5000  # Big batches now that indexes are dropped

def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

def get_local():
    return psycopg2.connect(**LOCAL_DB)

def get_supa(timeout_ms=600000):
    return psycopg2.connect(SUPABASE_CONN, connect_timeout=30,
                            options=f"-c statement_timeout={timeout_ms}")

def get_columns(conn, table):
    with conn.cursor() as cur:
        cur.execute("""SELECT column_name FROM information_schema.columns
                       WHERE table_schema='pharmaco_db' AND table_name=%s
                       ORDER BY ordinal_position""", (table,))
        return [r[0] for r in cur.fetchall()]

# ---- INDEX DEFINITIONS (from Supabase) ----
# We keep PKEYs, drop everything else
SECONDARY_INDEXES = {
    "bioactivities": [
        "CREATE INDEX idx_bioact_assay ON pharmaco_db.bioactivities USING btree (assay_id)",
        "CREATE INDEX idx_bioact_compound ON pharmaco_db.bioactivities USING btree (compound_id)",
        "CREATE INDEX idx_bioact_compound_target ON pharmaco_db.bioactivities USING btree (compound_id, target_id)",
        "CREATE INDEX idx_bioact_pchembl ON pharmaco_db.bioactivities USING btree (pchembl_value)",
        "CREATE INDEX idx_bioact_source ON pharmaco_db.bioactivities USING btree (source)",
        "CREATE INDEX idx_bioact_target ON pharmaco_db.bioactivities USING btree (target_id)",
        "CREATE INDEX idx_bioact_type ON pharmaco_db.bioactivities USING btree (activity_type)",
    ],
    "activity_tags": [
        "CREATE INDEX idx_atags_cliff ON pharmaco_db.activity_tags USING btree (is_activity_cliff)",
        "CREATE INDEX idx_atags_confidence ON pharmaco_db.activity_tags USING btree (confidence_level)",
        "CREATE INDEX idx_atags_pchembl_bin ON pharmaco_db.activity_tags USING btree (pchembl_bin)",
    ],
    "clinical_trials": [
        "CREATE INDEX idx_ct_compound ON pharmaco_db.clinical_trials USING btree (compound_id)",
        "CREATE INDEX idx_ct_phase ON pharmaco_db.clinical_trials USING btree (phase)",
        "CREATE INDEX idx_ct_status ON pharmaco_db.clinical_trials USING btree (status)",
        "CREATE UNIQUE INDEX uq_ct_nct ON pharmaco_db.clinical_trials USING btree (nct_id)",
    ],
    "compounds": [
        "CREATE INDEX idx_compounds_alogp ON pharmaco_db.compounds USING btree (alogp)",
        "CREATE INDEX idx_compounds_chembl ON pharmaco_db.compounds USING btree (chembl_id)",
        "CREATE INDEX idx_compounds_inchi_key ON pharmaco_db.compounds USING btree (inchi_key)",
        "CREATE INDEX idx_compounds_max_phase ON pharmaco_db.compounds USING btree (max_phase)",
        "CREATE INDEX idx_compounds_mw ON pharmaco_db.compounds USING btree (molecular_weight)",
        "CREATE INDEX idx_compounds_name_trgm ON pharmaco_db.compounds USING gin (pref_name pharmaco_db.gin_trgm_ops)",
        "CREATE INDEX idx_compounds_pubchem ON pharmaco_db.compounds USING btree (pubchem_cid)",
        "CREATE INDEX idx_compounds_smiles ON pharmaco_db.compounds USING hash (canonical_smiles)",
    ],
    "compound_tags": [
        "CREATE INDEX idx_ctags_drug ON pharmaco_db.compound_tags USING btree (drug_like)",
        "CREATE INDEX idx_ctags_lead ON pharmaco_db.compound_tags USING btree (lead_like)",
        "CREATE INDEX idx_ctags_lipinski ON pharmaco_db.compound_tags USING btree (lipinski_pass)",
        "CREATE INDEX idx_ctags_logp_bin ON pharmaco_db.compound_tags USING btree (logp_bin)",
        "CREATE INDEX idx_ctags_mw_bin ON pharmaco_db.compound_tags USING btree (mw_bin)",
        "CREATE INDEX idx_ctags_phase ON pharmaco_db.compound_tags USING btree (clinical_phase)",
        "CREATE INDEX idx_ctags_scaffold ON pharmaco_db.compound_tags USING btree (murcko_scaffold)",
    ],
    "disease_target_associations": [
        "CREATE INDEX idx_dta_area ON pharmaco_db.disease_target_associations USING btree (therapeutic_area)",
        "CREATE INDEX idx_dta_disease ON pharmaco_db.disease_target_associations USING btree (disease_id)",
        "CREATE INDEX idx_dta_score ON pharmaco_db.disease_target_associations USING btree (overall_score)",
        "CREATE INDEX idx_dta_target ON pharmaco_db.disease_target_associations USING btree (target_id)",
    ],
    "compound_structures": [],  # Only pkey, nothing to drop
}

# Table -> (order_by_col, conflict_col)
TABLE_CONFIG = {
    "clinical_trials":            ("id",          "id"),
    "compound_structures":        ("compound_id", "compound_id"),
    "bioactivities":              ("id",          "id"),
    "activity_tags":              ("activity_id", "activity_id"),
    "compounds":                  ("id",          "id"),
    "disease_target_associations":("id",          "id"),
    "compound_tags":              ("compound_id", "compound_id"),
}


def drop_indexes(supa, table):
    """Drop secondary indexes (keep pkey)."""
    indexes = SECONDARY_INDEXES.get(table, [])
    if not indexes:
        return
    log(f"  {table}: dropping {len(indexes)} secondary indexes...")
    for idx_def in indexes:
        # Extract index name from CREATE INDEX idx_name ON ...
        idx_name = idx_def.split()[2]
        try:
            with supa.cursor() as cur:
                cur.execute(f"DROP INDEX IF EXISTS pharmaco_db.{idx_name}")
            supa.commit()
        except Exception as e:
            supa.rollback()
            log(f"    WARN drop {idx_name}: {str(e)[:80]}")


def recreate_indexes(supa, table):
    """Recreate secondary indexes with extended timeout."""
    indexes = SECONDARY_INDEXES.get(table, [])
    if not indexes:
        return
    log(f"  {table}: recreating {len(indexes)} indexes...")
    for idx_def in indexes:
        idx_name = idx_def.split()[2]
        # Use IF NOT EXISTS for idempotency
        safe_def = idx_def.replace("CREATE INDEX ", "CREATE INDEX IF NOT EXISTS ").replace("CREATE UNIQUE INDEX ", "CREATE UNIQUE INDEX IF NOT EXISTS ")
        try:
            with supa.cursor() as cur:
                cur.execute("SET statement_timeout = '1800000'")  # 30 min for index creation
                cur.execute(safe_def)
            supa.commit()
            log(f"    OK {idx_name}")
        except Exception as e:
            supa.rollback()
            log(f"    WARN recreate {idx_name}: {str(e)[:80]}")


def push_table(table):
    """Push one table: drop indexes, bulk insert, recreate indexes."""
    order_col, conflict_col = TABLE_CONFIG[table]

    local = get_local()
    supa = get_supa()

    # Get counts
    with local.cursor() as cur:
        cur.execute(f"SELECT count(*) FROM pharmaco_db.{table}")
        local_n = cur.fetchone()[0]

    with supa.cursor() as cur:
        try:
            cur.execute(f"""SELECT reltuples::bigint FROM pg_class
                           JOIN pg_namespace ON pg_namespace.oid=relnamespace
                           WHERE relname=%s AND nspname='pharmaco_db'""", (table,))
            supa_n = int(cur.fetchone()[0])
            if supa_n < 0:
                supa_n = 0
        except:
            supa.rollback()
            supa_n = 0

    if local_n == 0:
        log(f"  {table}: empty locally, skip")
        local.close(); supa.close()
        return

    if supa_n >= local_n * 0.99:  # within 1% = up to date
        log(f"  {table}: Supabase ~up to date ({supa_n:,} vs {local_n:,}), skip")
        local.close(); supa.close()
        return

    to_push = local_n - supa_n
    log(f"  {table}: need ~{to_push:,} rows (local={local_n:,}, supa~{supa_n:,})")

    # Step 1: Drop secondary indexes
    drop_indexes(supa, table)
    supa.close()

    # Step 2: Bulk insert
    supa = get_supa()
    supa_cols = get_columns(supa, table)
    local_cols = get_columns(local, table)
    cols = [c for c in local_cols if c in supa_cols]
    col_list = ", ".join(cols)

    sql = f"INSERT INTO pharmaco_db.{table} ({col_list}) VALUES %s ON CONFLICT ({conflict_col}) DO NOTHING"

    start = time.time()
    total = 0
    errors = 0
    offset = max(0, supa_n)

    log(f"  {table}: starting bulk insert from offset {offset:,}...")

    while True:
        try:
            with local.cursor() as cur:
                cur.execute(
                    f"SELECT {col_list} FROM pharmaco_db.{table} "
                    f"ORDER BY {order_col} LIMIT {BATCH} OFFSET {offset}"
                )
                rows = cur.fetchall()
        except Exception as e:
            log(f"  {table}: local read error at {offset}: {e}")
            break

        if not rows:
            break

        try:
            with supa.cursor() as cur:
                psycopg2.extras.execute_values(cur, sql, rows, page_size=1000)
            supa.commit()
            total += len(rows)
            errors = 0

            if total % 20000 == 0:
                rate = total / max(time.time() - start, 1)
                remaining = max(0, to_push - total) / max(rate, 1)
                log(f"  {table}: {total + supa_n:,}/{local_n:,} "
                    f"({int(rate)}/s, ~{int(remaining/60)}m left)")
        except Exception as e:
            errors += 1
            err_msg = str(e)[:100]
            log(f"  {table}: ERROR at {offset}: {err_msg}")

            if errors >= 5:
                log(f"  {table}: too many errors, stopping")
                break

            # Safe reconnect
            try:
                supa.rollback()
            except:
                pass
            try:
                supa.close()
            except:
                pass
            time.sleep(3)
            try:
                supa = get_supa()
            except Exception as e2:
                log(f"  {table}: FATAL reconnect failed: {e2}")
                break

        offset += BATCH

    elapsed = time.time() - start
    rate = total / max(elapsed, 1)
    log(f"  {table}: bulk insert done — +{total:,} in {int(elapsed)}s ({int(rate)}/s)")

    # Step 3: Recreate indexes
    try:
        supa.close()
    except:
        pass
    supa = get_supa()
    recreate_indexes(supa, table)

    supa.close()
    local.close()
    return total


def main():
    log("=" * 60)
    log("PharmacoDB — FAST PUSH (drop idx → insert → recreate idx)")
    log("=" * 60)

    # Test connections
    try:
        l = get_local(); l.close()
        log("  Local DB OK")
    except Exception as e:
        log(f"  FATAL local: {e}"); sys.exit(1)
    try:
        s = get_supa(); s.close()
        log("  Supabase OK")
    except Exception as e:
        log(f"  FATAL supabase: {e}"); sys.exit(1)

    # Allow single table arg
    if len(sys.argv) > 1 and sys.argv[1] != "all":
        table = sys.argv[1]
        if table in TABLE_CONFIG:
            push_table(table)
        else:
            log(f"Unknown table: {table}")
            log(f"Available: {list(TABLE_CONFIG.keys())}")
        return

    start = time.time()

    # Push in order (smallest first)
    order = [
        "clinical_trials",          # ~5K remaining
        "compound_structures",      # ~2.7M
        "disease_target_associations",  # check
        "compounds",                # check
        "compound_tags",            # check
        "bioactivities",            # ~2.6M
        "activity_tags",            # ~4.8M
    ]

    for table in order:
        try:
            push_table(table)
        except Exception as e:
            log(f"  {table}: FATAL {e}")
        log("")

    elapsed = time.time() - start
    log(f"{'='*60}")
    log(f"ALL DONE in {int(elapsed)}s ({int(elapsed/60)}min)")
    log(f"{'='*60}")


if __name__ == "__main__":
    main()
