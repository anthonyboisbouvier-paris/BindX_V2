#!/usr/bin/env python3
"""
PharmacoDB — PubChem Bioassay Ingestion Pipeline
Fetches bioassay data from PubChem PUG REST for compounds already in the DB.
Adds CID cross-references and additional bioactivity data.
"""

import json
import time
import sys
import urllib.request
import urllib.error
import psycopg2
import psycopg2.extras
from datetime import datetime

DB_CONFIG = {
    "host": "localhost",
    "port": 5433,
    "dbname": "pharmaco",
    "user": "postgres",
    "password": "pharmaco_secret"
}

PUBCHEM_API = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
MAX_RETRIES = 3
BATCH_CID_LOOKUP = 100  # PubChem allows batch queries

def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

def api_get(url, retries=MAX_RETRIES):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={
                "Accept": "application/json",
                "User-Agent": "PharmacoDB-Ingestion/1.0"
            })
            with urllib.request.urlopen(req, timeout=30) as resp:
                return json.loads(resp.read().decode())
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if e.code == 503 or e.code == 429:
                time.sleep(5 * (attempt + 1))
            elif attempt < retries - 1:
                time.sleep(3 * (attempt + 1))
            else:
                return None
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(3 * (attempt + 1))
            else:
                return None
    return None

def api_post(url, data_str, retries=MAX_RETRIES):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url,
                data=data_str.encode('utf-8'),
                headers={
                    "Accept": "application/json",
                    "Content-Type": "application/x-www-form-urlencoded",
                    "User-Agent": "PharmacoDB-Ingestion/1.0"
                })
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.loads(resp.read().decode())
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if attempt < retries - 1:
                time.sleep(5 * (attempt + 1))
            else:
                return None
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(3 * (attempt + 1))
            else:
                return None
    return None

def get_conn():
    return psycopg2.connect(**DB_CONFIG)

def start_log(conn, source, step):
    with conn.cursor() as cur:
        cur.execute("""
            INSERT INTO pharmaco_db.ingestion_log (source, step, status)
            VALUES (%s, %s, 'running') RETURNING id
        """, (source, step))
        log_id = cur.fetchone()[0]
    conn.commit()
    return log_id

def update_log(conn, log_id, status, rows=None, error=None):
    with conn.cursor() as cur:
        cur.execute("""
            UPDATE pharmaco_db.ingestion_log
            SET status=%s, rows_inserted=%s, error_message=%s, completed_at=NOW()
            WHERE id=%s
        """, (status, rows, error, log_id))
    conn.commit()

# ============================================================
# 1. Resolve PubChem CIDs for existing compounds via InChIKey
# ============================================================
def resolve_pubchem_cids(conn):
    log("=" * 60)
    log("PHASE 1: Resolving PubChem CIDs from InChIKey")
    log("=" * 60)
    log_id = start_log(conn, 'pubchem', 'cid_resolution')

    # Get compounds with InChIKey but no PubChem CID
    with conn.cursor() as cur:
        cur.execute("""
            SELECT id, inchi_key FROM pharmaco_db.compounds
            WHERE inchi_key IS NOT NULL
            AND pubchem_cid IS NULL
            ORDER BY id
            LIMIT 100000
        """)
        compounds = cur.fetchall()

    log(f"  Found {len(compounds)} compounds to resolve")
    resolved = 0
    errors = 0

    # Process in batches
    for i in range(0, len(compounds), BATCH_CID_LOOKUP):
        batch = compounds[i:i + BATCH_CID_LOOKUP]
        inchikeys = [c[1] for c in batch]
        ids = [c[0] for c in batch]

        # PubChem batch lookup by InChIKey
        inchikey_str = ",".join(inchikeys)
        url = f"{PUBCHEM_API}/compound/inchikey/property/CID/JSON"
        data = api_post(url, f"inchikey={inchikey_str}")

        if data and "PropertyTable" in data:
            props = data["PropertyTable"].get("Properties", [])
            cid_map = {}
            for p in props:
                cid = p.get("CID")
                if cid:
                    # Map back via position (PubChem returns in order)
                    cid_map[p.get("CID")] = cid

            # Update CIDs - match by InChIKey
            for compound_id, inchikey in batch:
                # Try individual lookup for each
                single_url = f"{PUBCHEM_API}/compound/inchikey/{inchikey}/property/CID/JSON"
                single_data = api_get(single_url)
                if single_data and "PropertyTable" in single_data:
                    cid_props = single_data["PropertyTable"].get("Properties", [])
                    if cid_props:
                        cid = cid_props[0].get("CID")
                        if cid:
                            with conn.cursor() as cur:
                                cur.execute("""
                                    UPDATE pharmaco_db.compounds
                                    SET pubchem_cid = %s, updated_at = NOW()
                                    WHERE id = %s
                                """, (cid, compound_id))
                            resolved += 1
                time.sleep(0.25)  # PubChem rate limit: ~5 req/s

            conn.commit()
        else:
            errors += len(batch)

        if (i // BATCH_CID_LOOKUP) % 10 == 0:
            log(f"  Progress: {i}/{len(compounds)} ({resolved} resolved)")

        time.sleep(0.5)

    log(f"  DONE: {resolved} CIDs resolved, {errors} errors")
    update_log(conn, log_id, 'completed', resolved)
    return resolved

# ============================================================
# 2. Fetch PubChem BioAssay data for top compounds
# ============================================================
def ingest_pubchem_bioassays(conn):
    log("=" * 60)
    log("PHASE 2: Ingesting PubChem BioAssay data")
    log("=" * 60)
    log_id = start_log(conn, 'pubchem', 'bioassays')

    # Get compounds with PubChem CIDs (prioritize approved drugs)
    with conn.cursor() as cur:
        cur.execute("""
            SELECT id, pubchem_cid, chembl_id FROM pharmaco_db.compounds
            WHERE pubchem_cid IS NOT NULL
            ORDER BY max_phase DESC NULLS LAST, id
            LIMIT 10000
        """)
        compounds = cur.fetchall()

    log(f"  Found {len(compounds)} compounds with PubChem CIDs")
    total_activities = 0
    target_map = {}

    # Build target name -> id map for fuzzy matching
    with conn.cursor() as cur:
        cur.execute("SELECT id, gene_name, pref_name FROM pharmaco_db.targets WHERE gene_name IS NOT NULL")
        for row in cur.fetchall():
            target_map[row[1].upper()] = row[0]
            if row[2]:
                target_map[row[2].upper()] = row[0]

    for i, (compound_id, cid, chembl_id) in enumerate(compounds):
        if i % 100 == 0 and i > 0:
            log(f"  Progress: {i}/{len(compounds)} ({total_activities} activities)")

        # Get bioassay results for this CID
        url = f"{PUBCHEM_API}/compound/cid/{cid}/assaysummary/JSON"
        data = api_get(url)

        if not data or "Table" not in data:
            time.sleep(0.3)
            continue

        table = data["Table"]
        columns = table.get("Columns", {}).get("Column", [])
        rows = table.get("Row", [])

        if not rows:
            time.sleep(0.3)
            continue

        # Parse column indices
        col_idx = {c: i for i, c in enumerate(columns)}
        aid_idx = col_idx.get("AID")
        activity_idx = col_idx.get("Activity Outcome")
        target_name_idx = col_idx.get("Target Name")
        target_gi_idx = col_idx.get("Target GI")

        batch = []
        for row in rows[:50]:  # Limit per compound
            cells = row.get("Cell", [])
            if not cells:
                continue

            activity_outcome = cells[activity_idx] if activity_idx is not None and activity_idx < len(cells) else None
            if activity_outcome not in ["Active", "Inactive"]:
                continue

            target_name = cells[target_name_idx] if target_name_idx is not None and target_name_idx < len(cells) else None
            target_id = None
            if target_name:
                target_id = target_map.get(target_name.upper())

            batch.append({
                "compound_id": compound_id,
                "target_id": target_id,
                "activity_type": "PubChem_BioAssay",
                "activity_class": "active" if activity_outcome == "Active" else "inactive",
                "source": "pubchem",
                "source_id": str(cells[aid_idx]) if aid_idx is not None and aid_idx < len(cells) else None,
            })

        if batch:
            with conn.cursor() as cur:
                psycopg2.extras.execute_values(cur, """
                    INSERT INTO pharmaco_db.bioactivities
                    (compound_id, target_id, activity_type, activity_class, source, source_id)
                    VALUES %s
                """, [(
                    a["compound_id"], a["target_id"], a["activity_type"],
                    a["activity_class"], a["source"], a["source_id"]
                ) for a in batch], page_size=200)
            conn.commit()
            total_activities += len(batch)

        # Add cross-reference
        with conn.cursor() as cur:
            cur.execute("""
                INSERT INTO pharmaco_db.cross_references (compound_id, target_id, db_name, db_id, url)
                VALUES (%s, NULL, 'PubChem', %s, %s)
                ON CONFLICT DO NOTHING
            """, (compound_id, str(cid), f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"))
        conn.commit()

        time.sleep(0.3)  # PubChem rate limit

    log(f"  DONE: {total_activities} PubChem activities inserted")
    update_log(conn, log_id, 'completed', total_activities)
    return total_activities

def main():
    log("=" * 60)
    log("PharmacoDB — PubChem Ingestion Pipeline")
    log("=" * 60)

    conn = get_conn()
    try:
        resolve_pubchem_cids(conn)
        ingest_pubchem_bioassays(conn)
    except Exception as e:
        log(f"FATAL ERROR: {e}")
        import traceback
        traceback.print_exc()
    finally:
        conn.close()

if __name__ == "__main__":
    main()
