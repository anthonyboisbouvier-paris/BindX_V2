#!/usr/bin/env python3
"""
PharmacoDB — UniProt Individual Enrichment for Secondary Accessions
Fetches remaining unenriched targets individually (with 303 redirect following).
Uses ThreadPoolExecutor for moderate parallelism (3 workers).
"""

import sys
import os
import time
import json
import urllib.request
import urllib.error
import psycopg2
import psycopg2.extras
import traceback
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

DB_CONFIG = {
    "host": "localhost", "port": 5433,
    "dbname": "pharmaco", "user": "postgres", "password": "pharmaco_secret"
}

UNIPROT_API = "https://rest.uniprot.org/uniprotkb"
WORKERS = 3
TIMEOUT = 15

def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

def fetch_entry(accession):
    """Fetch a single UniProt entry with redirect following."""
    url = f"{UNIPROT_API}/{accession}.json"
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={
                "Accept": "application/json",
                "User-Agent": "PharmacoDB/2.0"
            })
            opener = urllib.request.build_opener(urllib.request.HTTPRedirectHandler)
            resp = opener.open(req, timeout=TIMEOUT)
            data = json.loads(resp.read().decode())
            resp.close()
            data["_queried_accession"] = accession
            return data
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if e.code == 429:
                time.sleep(3 * (attempt + 1))
                continue
        except Exception:
            if attempt < 2:
                time.sleep(1)
    return None

def main():
    # Import parse/update logic from the full script
    from ingest_uniprot_full import (
        parse_entry, build_update_query, build_update_params,
        insert_cross_references, start_log, update_log
    )

    conn = psycopg2.connect(**DB_CONFIG)

    # Get unenriched targets
    with conn.cursor() as cur:
        cur.execute("""
            SELECT id, uniprot_id, chembl_id FROM pharmaco_db.targets
            WHERE uniprot_id IS NOT NULL AND uniprot_enriched_at IS NULL
            ORDER BY id
        """)
        targets = cur.fetchall()

    total = len(targets)
    log(f"Found {total} unenriched targets (individual mode, {WORKERS} workers)")

    if total == 0:
        log("Nothing to do.")
        conn.close()
        return

    log_id = start_log(conn, "uniprot_individual", "secondary_accessions")
    update_query = build_update_query()

    enriched = 0
    errors = 0
    skipped = 0
    start_time = time.time()

    # Process in chunks with ThreadPoolExecutor
    chunk_size = WORKERS * 10  # process 30 at a time

    for chunk_start in range(0, total, chunk_size):
        chunk = targets[chunk_start:chunk_start + chunk_size]

        # Fetch entries in parallel
        results = {}
        with ThreadPoolExecutor(max_workers=WORKERS) as pool:
            futures = {}
            for target_id, uniprot_id, chembl_id in chunk:
                f = pool.submit(fetch_entry, uniprot_id)
                futures[f] = (target_id, uniprot_id, chembl_id)

            for f in as_completed(futures):
                target_id, uniprot_id, chembl_id = futures[f]
                entry = f.result()
                results[(target_id, uniprot_id, chembl_id)] = entry

        # Process results sequentially (DB writes)
        for (target_id, uniprot_id, chembl_id), entry in results.items():
            if not entry:
                skipped += 1
                # Mark as enriched so we don't retry
                try:
                    with conn.cursor() as cur:
                        cur.execute(
                            "UPDATE pharmaco_db.targets SET uniprot_enriched_at = NOW() WHERE id = %s",
                            (target_id,)
                        )
                    conn.commit()
                except:
                    conn.rollback()
                continue

            try:
                parsed = parse_entry(entry)
                params = build_update_params(parsed, uniprot_id)

                with conn.cursor() as cur:
                    cur.execute(update_query, params)

                insert_cross_references(conn, target_id, parsed)
                conn.commit()
                enriched += 1

            except Exception as e:
                conn.rollback()
                errors += 1
                if errors <= 5:
                    log(f"  Error {uniprot_id}: {e}")
                    traceback.print_exc()

        # Progress
        done = chunk_start + len(chunk)
        if done % 100 == 0 or done >= total:
            elapsed = time.time() - start_time
            rate = (enriched + skipped) / elapsed if elapsed > 0 else 0
            remaining = total - done
            eta = int(remaining / rate) if rate > 0 else 0
            log(f"  Progress: {done}/{total} ({enriched} enriched, {skipped} skipped/404, {errors} errors) [{rate:.1f}/s, ETA {eta//60}m{eta%60:02d}s]")

        # Small sleep between chunks to be nice to UniProt
        time.sleep(0.2)

    elapsed = time.time() - start_time
    log(f"\nDONE in {int(elapsed)}s: {enriched} enriched, {skipped} skipped, {errors} errors")
    update_log(conn, log_id, "completed", enriched)
    conn.close()

if __name__ == "__main__":
    main()
