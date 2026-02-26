#!/usr/bin/env python3
"""Fill adverse_events from OpenFDA drug/event API."""
import json
import time
import urllib.request
import urllib.parse
import urllib.error
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

import psycopg2
import psycopg2.extras

DB = {"host": "localhost", "port": 5433, "dbname": "pharmaco", "user": "postgres", "password": "postgres"}

def log(msg: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

def main():
    log("OpenFDA Adverse Events")
    conn = psycopg2.connect(**DB)

    # Get drugs with names - try multiple name formats
    with conn.cursor() as cur:
        cur.execute("""
            SELECT id, LOWER(pref_name) FROM pharmaco_db.compounds
            WHERE is_drug = true AND pref_name IS NOT NULL AND max_phase >= 4
            ORDER BY id
        """)
        drugs = cur.fetchall()
    log(f"  Approved drugs (phase 4): {len(drugs):,}")

    # Already done
    with conn.cursor() as cur:
        cur.execute("SELECT DISTINCT compound_id FROM pharmaco_db.adverse_events")
        done = {r[0] for r in cur.fetchall()}
    log(f"  Already done: {len(done):,}")

    batch = []
    queried = 0
    found = 0
    errors = 0

    for i, (cid, name) in enumerate(drugs):
        if cid in done:
            continue

        # Try different search strategies
        found_any = False
        for search_field in ["generic_name", "brand_name", "substance_name"]:
            encoded = urllib.parse.quote(name)
            url = (f"https://api.fda.gov/drug/event.json?"
                   f"search=patient.drug.openfda.{search_field}:%22{encoded}%22"
                   f"&count=patient.reaction.reactionmeddrapt.exact&limit=25")
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "PharmacoDB/1.0"})
                with urllib.request.urlopen(req, timeout=15) as resp:
                    data = json.loads(resp.read().decode())
                results = data.get("results", [])
                for entry in results:
                    reaction = entry.get("term", "")
                    count = entry.get("count", 0)
                    if reaction and count > 0:
                        batch.append((cid, name, reaction, count, "openfda"))
                if results:
                    found += 1
                    found_any = True
                    break
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    continue  # Try next search field
                elif e.code == 429:
                    time.sleep(60)
                    errors += 1
                    break
                else:
                    continue
            except Exception:
                errors += 1
                break

        queried += 1
        time.sleep(0.35)

        if (i + 1) % 100 == 0:
            if batch:
                with conn.cursor() as cur:
                    psycopg2.extras.execute_values(cur,
                        """INSERT INTO pharmaco_db.adverse_events
                           (compound_id, drug_name, reaction, count, source)
                           VALUES %s ON CONFLICT DO NOTHING""",
                        batch, page_size=2000)
                conn.commit()
                batch = []
            with conn.cursor() as cur:
                cur.execute("SELECT COUNT(*) FROM pharmaco_db.adverse_events")
                total = cur.fetchone()[0]
            log(f"  [{i+1}/{len(drugs)}] queried={queried} found={found} AEs={total:,} errors={errors}")

        if errors > 100:
            log("  Too many errors, stopping")
            break

    if batch:
        with conn.cursor() as cur:
            psycopg2.extras.execute_values(cur,
                """INSERT INTO pharmaco_db.adverse_events
                   (compound_id, drug_name, reaction, count, source)
                   VALUES %s ON CONFLICT DO NOTHING""",
                batch, page_size=2000)
        conn.commit()

    with conn.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.adverse_events")
        total = cur.fetchone()[0]
    conn.close()
    log(f"  DONE: {total:,} adverse events")

if __name__ == "__main__":
    main()
