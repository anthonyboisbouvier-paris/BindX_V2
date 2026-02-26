#!/usr/bin/env python3
"""
Fix KEGG drug_info: re-fetch pathway_ids, pathway_names, therapeutic_category.
Uses ThreadPoolExecutor for 10x speedup.
"""
import psycopg2
import psycopg2.extras
import urllib.request
import time
import sys
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

LOCAL_DB = {
    "host": "localhost", "port": 5433,
    "dbname": "pharmaco", "user": "postgres", "password": "postgres"
}

WORKERS = 15  # Parallel threads (KEGG allows ~10 req/s, burst OK)


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def fetch_kegg_detail(kegg_id):
    """Fetch and parse a single KEGG drug entry. Returns (kegg_id, efficacy, pids, pnames)."""
    url = f"https://rest.kegg.jp/get/{kegg_id}"
    try:
        req = urllib.request.Request(url)
        resp = urllib.request.urlopen(req, timeout=15)
        text = resp.read().decode("utf-8")
    except Exception:
        return kegg_id, None, None, None

    efficacy = None
    pathway_ids = []
    pathway_names = []
    in_target = False

    for line in text.split("\n"):
        if line.startswith("EFFICACY"):
            parts = line.split(None, 1)
            if len(parts) > 1:
                efficacy = parts[1].strip()
        elif line.startswith("TARGET"):
            in_target = True
        elif line.lstrip().startswith("PATHWAY") and in_target:
            rest = line.strip()
            if rest.startswith("PATHWAY"):
                rest = rest[7:].strip()
            if rest:
                parts = rest.split(None, 1)
                pid = parts[0].split("(")[0]
                pathway_ids.append(pid)
                if len(parts) > 1:
                    pathway_names.append(parts[1])
        elif line.startswith("BRITE") or line.startswith("DBLINKS"):
            in_target = False
        elif not line.startswith(" ") and not line.startswith("\t") and not line.startswith("///"):
            in_target = False
        elif line.startswith("///"):
            break

    return (
        kegg_id,
        [efficacy] if efficacy else None,  # Array
        pathway_ids if pathway_ids else None,
        pathway_names if pathway_names else None,
    )


def main():
    log("=== KEGG Drug Info Fix (parallel) ===")

    conn = psycopg2.connect(**LOCAL_DB)

    with conn.cursor() as cur:
        cur.execute("""SELECT id, kegg_drug_id FROM pharmaco_db.kegg_drug_info
                       WHERE kegg_drug_id IS NOT NULL ORDER BY id""")
        entries = cur.fetchall()

    total = len(entries)
    log(f"Found {total:,} entries — fetching with {WORKERS} threads")

    id_map = {kegg_id: row_id for row_id, kegg_id in entries}
    kegg_ids = [kegg_id for _, kegg_id in entries]

    start = time.time()
    updated = 0
    errors = 0
    results = []

    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futures = {pool.submit(fetch_kegg_detail, kid): kid for kid in kegg_ids}

        for i, future in enumerate(as_completed(futures)):
            kegg_id, categories, pids, pnames = future.result()

            if categories or pids:
                results.append((categories, pids, pnames, id_map[kegg_id]))
                updated += 1
            elif categories is None and pids is None:
                errors += 1

            if (i + 1) % 500 == 0:
                # Batch update
                with conn.cursor() as cur:
                    for cat, pid, pname, rid in results:
                        cur.execute("""UPDATE pharmaco_db.kegg_drug_info
                                       SET therapeutic_category = %s,
                                           pathway_ids = %s,
                                           pathway_names = %s
                                       WHERE id = %s""", (cat, pid, pname, rid))
                conn.commit()
                results.clear()

                elapsed = time.time() - start
                rate = (i + 1) / elapsed
                remaining = (total - i - 1) / max(rate, 1)
                log(f"  {i+1:,}/{total:,} — updated={updated}, errors={errors} "
                    f"({rate:.0f}/s, ~{remaining/60:.0f}min left)")

    # Final batch
    if results:
        with conn.cursor() as cur:
            for cat, pid, pname, rid in results:
                cur.execute("""UPDATE pharmaco_db.kegg_drug_info
                               SET therapeutic_category = %s,
                                   pathway_ids = %s,
                                   pathway_names = %s
                               WHERE id = %s""", (cat, pid, pname, rid))
        conn.commit()

    elapsed = time.time() - start
    log(f"Done in {int(elapsed)}s: {updated:,} updated, {errors:,} errors")
    conn.close()


if __name__ == "__main__":
    main()
