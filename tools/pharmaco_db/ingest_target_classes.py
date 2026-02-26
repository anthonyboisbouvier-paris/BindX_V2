#!/usr/bin/env python3
"""
PharmacoDB — Target Classification Enrichment
Uses ChEMBL target classification API to add protein family hierarchy.
Also fetches NCBI Gene summaries for targets.
"""

import json
import time
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

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
NCBI_API = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

def api_get(url, retries=3):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={
                "Accept": "application/json",
                "User-Agent": "PharmacoDB-Ingestion/1.0"
            })
            with urllib.request.urlopen(req, timeout=30) as resp:
                return json.loads(resp.read().decode())
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
# 1. ChEMBL Target Classification
# ============================================================
def enrich_target_classifications(conn):
    log("=" * 60)
    log("Target Classification Enrichment (ChEMBL)")
    log("=" * 60)
    log_id = start_log(conn, 'chembl', 'target_classification')

    # Get all targets
    with conn.cursor() as cur:
        cur.execute("""
            SELECT id, chembl_id FROM pharmaco_db.targets
            WHERE chembl_id IS NOT NULL
            AND protein_class_l1 IS NULL
            ORDER BY id
        """)
        targets = cur.fetchall()

    log(f"  Found {len(targets)} targets to classify")
    enriched = 0

    for i, (target_id, chembl_id) in enumerate(targets):
        if i % 200 == 0 and i > 0:
            log(f"  Progress: {i}/{len(targets)} ({enriched} enriched)")

        url = f"{CHEMBL_API}/target_classification.json?target_chembl_id={chembl_id}&limit=10"
        data = api_get(url)

        if not data or not data.get("target_classifications"):
            time.sleep(0.15)
            continue

        classifications = data["target_classifications"]
        if classifications:
            cls = classifications[0]
            l1 = cls.get("l1")
            l2 = cls.get("l2")
            l3 = cls.get("l3")
            protein_class = cls.get("protein_class_desc")

            # Determine protein family from classification
            family = None
            if l1:
                family_map = {
                    "Enzyme": l2 or "Enzyme",
                    "Membrane receptor": l2 or "Receptor",
                    "Ion channel": "Ion channel",
                    "Transporter": "Transporter",
                    "Transcription factor": "Transcription factor",
                    "Epigenetic regulator": "Epigenetic regulator",
                    "Structural protein": "Structural protein",
                    "Surface antigen": "Surface antigen",
                    "Secreted protein": "Secreted protein",
                }
                family = family_map.get(l1, l1)

            with conn.cursor() as cur:
                cur.execute("""
                    UPDATE pharmaco_db.targets SET
                        protein_class_l1 = %s,
                        protein_class_l2 = %s,
                        protein_class_l3 = %s,
                        protein_family = COALESCE(protein_family, %s),
                        updated_at = NOW()
                    WHERE id = %s
                """, (l1, l2, l3, family, target_id))
            conn.commit()
            enriched += 1

        time.sleep(0.15)

    log(f"  DONE: {enriched} targets classified")
    update_log(conn, log_id, 'completed', enriched)
    return enriched

# ============================================================
# 2. NCBI Gene Summaries
# ============================================================
def enrich_ncbi_genes(conn):
    log("=" * 60)
    log("NCBI Gene Summary Enrichment")
    log("=" * 60)
    log_id = start_log(conn, 'ncbi', 'gene_summaries')

    # Get targets with NCBI gene IDs but no summary
    with conn.cursor() as cur:
        cur.execute("""
            SELECT id, ncbi_gene_id FROM pharmaco_db.targets
            WHERE ncbi_gene_id IS NOT NULL
            AND ncbi_gene_summary IS NULL
            ORDER BY id
        """)
        targets = cur.fetchall()

    log(f"  Found {len(targets)} targets with NCBI gene IDs")
    enriched = 0

    # Process in batches of 50
    for i in range(0, len(targets), 50):
        batch = targets[i:i + 50]
        gene_ids = [str(t[1]) for t in batch]
        id_str = ",".join(gene_ids)

        url = f"{NCBI_API}/esummary.fcgi?db=gene&id={id_str}&retmode=json"
        data = api_get(url)

        if not data or "result" not in data:
            time.sleep(1)
            continue

        results = data["result"]
        for target_id, ncbi_gene_id in batch:
            gene_data = results.get(str(ncbi_gene_id))
            if not gene_data or "error" in gene_data:
                continue

            summary = gene_data.get("summary", "")
            if summary:
                with conn.cursor() as cur:
                    cur.execute("""
                        UPDATE pharmaco_db.targets SET
                            ncbi_gene_summary = %s,
                            updated_at = NOW()
                        WHERE id = %s
                    """, (summary[:3000], target_id))
                enriched += 1

        conn.commit()

        if (i // 50) % 20 == 0 and i > 0:
            log(f"  Progress: {i}/{len(targets)} ({enriched} enriched)")

        time.sleep(0.4)  # NCBI rate limit: ~3 req/s without API key

    log(f"  DONE: {enriched} gene summaries added")
    update_log(conn, log_id, 'completed', enriched)
    return enriched

# ============================================================
# 3. Mark druggable targets & count approved drugs
# ============================================================
def compute_druggability(conn):
    log("=" * 60)
    log("Computing target druggability metrics")
    log("=" * 60)
    log_id = start_log(conn, 'computed', 'druggability')

    with conn.cursor() as cur:
        # Count approved drugs per target
        cur.execute("""
            UPDATE pharmaco_db.targets t SET
                num_approved_drugs = sub.cnt,
                is_druggable = TRUE,
                updated_at = NOW()
            FROM (
                SELECT dm.target_id, COUNT(DISTINCT dm.compound_id) as cnt
                FROM pharmaco_db.drug_mechanisms dm
                JOIN pharmaco_db.compounds c ON c.id = dm.compound_id
                WHERE c.max_phase >= 4
                GROUP BY dm.target_id
            ) sub
            WHERE t.id = sub.target_id
        """)
        n1 = cur.rowcount
        log(f"  Marked {n1} targets with approved drugs")

        # Mark targets with active compounds as potentially druggable
        cur.execute("""
            UPDATE pharmaco_db.targets t SET
                is_druggable = TRUE,
                updated_at = NOW()
            FROM (
                SELECT DISTINCT target_id
                FROM pharmaco_db.bioactivities
                WHERE pchembl_value >= 6.0
                AND target_id IS NOT NULL
            ) sub
            WHERE t.id = sub.target_id
            AND t.is_druggable IS NULL
        """)
        n2 = cur.rowcount
        log(f"  Marked {n2} additional targets as druggable (pChEMBL >= 6)")

    conn.commit()
    update_log(conn, log_id, 'completed', n1 + n2)
    return n1 + n2

def main():
    conn = get_conn()
    try:
        enrich_target_classifications(conn)
        enrich_ncbi_genes(conn)
        compute_druggability(conn)
    except Exception as e:
        log(f"FATAL ERROR: {e}")
        import traceback
        traceback.print_exc()
    finally:
        conn.close()

if __name__ == "__main__":
    main()
