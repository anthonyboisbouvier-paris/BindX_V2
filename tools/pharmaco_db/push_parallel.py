#!/usr/bin/env python3
"""
Push ALL PharmacoDB tables to Supabase in parallel.
Small tables first, big tables in parallel threads.
v2: fix ORDER BY for tables without 'id', smaller batches, longer timeouts.
"""
import psycopg2
import psycopg2.extras
import time
import sys
import threading
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

LOCAL_DB = {
    "host": "localhost", "port": 5433,
    "dbname": "pharmaco", "user": "postgres", "password": "postgres"
}
SUPABASE_DB = {
    "host": "aws-0-us-west-2.pooler.supabase.com", "port": 6543,
    "dbname": "postgres",
    "user": "postgres.webkntghfzscrnuixfba",
    "password": "n1qsKxvAWJV954vR",
    "options": "-c statement_timeout=300000",  # 5 min timeout
}

BATCH = 1000  # Smaller batches to avoid timeouts
lock = threading.Lock()

def log(msg):
    with lock:
        print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

def get_local():
    return psycopg2.connect(**LOCAL_DB)

def get_supa():
    return psycopg2.connect(**SUPABASE_DB, connect_timeout=30)

def count_fast(conn, table):
    """Fast approximate count, fallback to exact."""
    with conn.cursor() as cur:
        # Try fast estimate first
        cur.execute(f"""SELECT reltuples::bigint FROM pg_class
                       JOIN pg_namespace ON pg_namespace.oid = relnamespace
                       WHERE relname=%s AND nspname='pharmaco_db'""", (table,))
        row = cur.fetchone()
        if row and row[0] > 0:
            return int(row[0])
        # Fallback to exact
        cur.execute(f"SELECT COUNT(*) FROM pharmaco_db.{table}")
        return cur.fetchone()[0]

def count_exact(conn, table):
    with conn.cursor() as cur:
        cur.execute(f"SELECT COUNT(*) FROM pharmaco_db.{table}")
        return cur.fetchone()[0]

def get_columns(conn, table):
    with conn.cursor() as cur:
        cur.execute("""SELECT column_name FROM information_schema.columns
                       WHERE table_schema='pharmaco_db' AND table_name=%s
                       ORDER BY ordinal_position""", (table,))
        return [r[0] for r in cur.fetchall()]

def get_first_column(conn, table):
    """Get first column name (for ORDER BY when no 'id')."""
    cols = get_columns(conn, table)
    return cols[0] if cols else None


# Map table -> order_by column (for tables without 'id')
ORDER_BY_MAP = {
    "compound_structures": "compound_id",
    "activity_tags": "activity_id",
    "compounds": "id",
    "bioactivities": "id",
    "compound_tags": "id",
    "disease_target_associations": "id",
    "clinical_trials": "id",
}


def push_table(table, conflict_col=None):
    """Push one table from local to Supabase."""
    local = get_local()
    supa = get_supa()

    try:
        local_n = count_exact(local, table)
    except Exception as e:
        log(f"  {table}: ERROR counting local: {e}")
        local.close(); supa.close()
        return 0

    try:
        supa_n = count_fast(supa, table)
    except Exception as e:
        log(f"  {table}: ERROR counting Supabase (using 0): {e}")
        supa.rollback()
        supa_n = 0

    if local_n == 0:
        log(f"  {table}: empty locally, skip")
        local.close(); supa.close()
        return 0

    if supa_n >= local_n:
        log(f"  {table}: Supabase up to date ({supa_n:,}), skip")
        local.close(); supa.close()
        return 0

    # Get column list from Supabase (to match schema)
    supa_cols = get_columns(supa, table)
    local_cols = get_columns(local, table)
    cols = [c for c in local_cols if c in supa_cols]
    col_list = ", ".join(cols)

    # Determine ORDER BY column
    order_col = ORDER_BY_MAP.get(table)
    if not order_col or order_col not in cols:
        order_col = cols[0]  # fallback to first column

    # Determine conflict clause
    conflict_clause = ""
    if conflict_col and conflict_col in cols:
        conflict_clause = f"ON CONFLICT ({conflict_col}) DO NOTHING"
    elif 'id' in cols:
        conflict_clause = "ON CONFLICT (id) DO NOTHING"

    sql = f"INSERT INTO pharmaco_db.{table} ({col_list}) VALUES %s {conflict_clause}"

    to_push = local_n - supa_n
    log(f"  {table}: pushing ~{to_push:,} rows (local={local_n:,}, supa={supa_n:,})...")
    start = time.time()
    total = 0
    errors = 0
    offset = supa_n  # Start from where Supabase left off

    while True:
        try:
            with local.cursor() as cur:
                cur.execute(
                    f"SELECT {col_list} FROM pharmaco_db.{table} "
                    f"ORDER BY {order_col} LIMIT {BATCH} OFFSET {offset}"
                )
                rows = cur.fetchall()
        except Exception as e:
            log(f"  {table}: ERROR reading local at offset {offset}: {str(e)[:80]}")
            break

        if not rows:
            break

        try:
            with supa.cursor() as cur:
                psycopg2.extras.execute_values(cur, sql, rows, page_size=200)
            supa.commit()
            total += len(rows)
            errors = 0  # Reset consecutive error count

            if total % 20000 == 0:
                rate = total / max(time.time() - start, 1)
                remaining = max(0, to_push - total) / max(rate, 1)
                log(f"  {table}: {total + supa_n:,}/{local_n:,} "
                    f"({int(rate)}/s, ~{int(remaining/60)}m left)")
        except Exception as e:
            errors += 1
            supa.rollback()
            log(f"  {table}: ERROR at offset {offset}: {str(e)[:100]}")
            if errors >= 5:
                log(f"  {table}: too many consecutive errors, stopping")
                break
            # Reconnect
            try:
                supa.close()
            except:
                pass
            try:
                supa = get_supa()
            except Exception as e2:
                log(f"  {table}: FATAL reconnect failed: {e2}")
                break
            time.sleep(2)  # Back off

        offset += BATCH

    elapsed = time.time() - start
    try:
        final_n = count_fast(supa, table)
        log(f"  {table}: DONE ~{final_n:,} rows on Supabase in {int(elapsed)}s (+{total:,} pushed)")
    except:
        log(f"  {table}: DONE in {int(elapsed)}s (+{total:,} pushed)")

    local.close()
    try:
        supa.close()
    except:
        pass
    return total


# Tables to push: (table_name, conflict_column)
# conflict_column used in ON CONFLICT clause
ALL_TABLES = [
    # Small (< 100K)
    ("target_interactions", "id"),
    ("kegg_drug_info", "id"),
    ("clinical_trials", "id"),
    # Large — the real work
    ("compound_structures", "compound_id"),
    ("bioactivities", "id"),
    ("activity_tags", "activity_id"),
    ("compounds", "chembl_id"),
    ("disease_target_associations", "id"),
    ("compound_tags", "id"),
]


def main():
    log("=" * 60)
    log("PharmacoDB — PUSH TO SUPABASE (v2)")
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

    step = sys.argv[1] if len(sys.argv) > 1 else "all"

    if step != "all":
        # Push single table
        conflict = dict(ALL_TABLES).get(step, "id")
        push_table(step, conflict)
        return

    start = time.time()

    # Push tables sequentially (safer than parallel for Supabase pooler)
    for table, conflict in ALL_TABLES:
        try:
            push_table(table, conflict)
        except Exception as e:
            log(f"  {table}: ERROR {e}")

    elapsed = time.time() - start
    log(f"\n{'='*60}")
    log(f"PUSH COMPLETE in {int(elapsed)}s ({int(elapsed/60)}min)")
    log("=" * 60)

if __name__ == "__main__":
    main()
