#!/usr/bin/env python3
"""Fill disease_target_associations from Open Targets GraphQL API."""
import json
import time
import urllib.request
from datetime import datetime

import psycopg2
import psycopg2.extras

DB = {"host": "localhost", "port": 5433, "dbname": "pharmaco", "user": "postgres", "password": "postgres"}
API = "https://api.platform.opentargets.org/api/v4/graphql"

def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

def gql(query, retries=3):
    payload = json.dumps({"query": query}).encode()
    for attempt in range(retries):
        try:
            req = urllib.request.Request(API, data=payload,
                headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=20) as resp:
                return json.loads(resp.read().decode())
        except Exception:
            if attempt < retries - 1:
                time.sleep(2 * (attempt + 1))
            else:
                raise

def main():
    log("Open Targets Disease-Target Associations")
    conn = psycopg2.connect(**DB)

    with conn.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.disease_target_associations")
        existing = cur.fetchone()[0]
    log(f"  Existing: {existing:,}")

    done_tids = set()
    if existing > 0:
        with conn.cursor() as cur:
            cur.execute("SELECT DISTINCT target_id FROM pharmaco_db.disease_target_associations")
            done_tids = {r[0] for r in cur.fetchall()}

    with conn.cursor() as cur:
        cur.execute("""SELECT id, gene_name FROM pharmaco_db.targets
                       WHERE gene_name IS NOT NULL AND organism = 'Homo sapiens'
                       ORDER BY id""")
        all_genes = cur.fetchall()

    genes = [(r[0], r[1]) for r in all_genes if r[0] not in done_tids]
    log(f"  Targets to query: {len(genes):,} (skipping {len(done_tids):,} done)")

    batch = []
    queried = 0
    found = 0
    errors = 0

    for i, (tid, gene) in enumerate(genes):
        try:
            # Step 1: gene -> Ensembl ID
            q1 = '{ search(queryString: "%s", entityNames: ["target"], page: {index: 0, size: 1}) { hits { id name } } }' % gene.replace('"', '')
            d1 = gql(q1)
            hits = d1.get("data", {}).get("search", {}).get("hits", [])
            eid = hits[0]["id"] if hits else None

            if not eid:
                queried += 1
                time.sleep(0.2)
                continue

            # Step 2: diseases
            q2 = '{ target(ensemblId: "%s") { associatedDiseases(page: {index: 0, size: 50}) { rows { disease { id name therapeuticAreas { id name } } score datasourceScores { id score } } } } }' % eid
            d2 = gql(q2)
            target = d2.get("data", {}).get("target")
            if not target:
                queried += 1
                time.sleep(0.2)
                continue

            rows = target.get("associatedDiseases", {}).get("rows", [])
            for r in rows:
                disease = r.get("disease", {})
                ta_list = disease.get("therapeuticAreas", [])
                scores = {s["id"]: s["score"] for s in r.get("datasourceScores", [])}
                batch.append((
                    tid, disease.get("id", ""), disease.get("name", ""),
                    ta_list[0]["name"] if ta_list else "",
                    r.get("score", 0),
                    scores.get("ot_genetics_portal"),
                    scores.get("cancer_gene_census"),
                    scores.get("chembl"),
                    scores.get("europepmc"),
                    scores.get("expression_atlas"),
                    scores.get("impc"),
                    "open_targets",
                ))

            if rows:
                found += 1
            queried += 1
        except Exception as e:
            errors += 1
            if errors <= 5:
                log(f"  Error #{errors} on {gene}: {e}")
            if errors % 100 == 0:
                log(f"  {errors} errors, pausing 15s...")
                time.sleep(15)
            time.sleep(1)
            continue

        time.sleep(0.25)

        # Batch insert every 100 targets
        if (i + 1) % 100 == 0 and batch:
            with conn.cursor() as cur:
                psycopg2.extras.execute_values(cur,
                    """INSERT INTO pharmaco_db.disease_target_associations
                       (target_id, disease_id, disease_name, therapeutic_area,
                        overall_score, genetic_score, somatic_score,
                        known_drug_score, literature_score, rna_expression_score,
                        animal_model_score, source)
                       VALUES %s ON CONFLICT DO NOTHING""",
                    batch, page_size=5000)
            conn.commit()
            batch_len = len(batch)
            batch = []

            with conn.cursor() as cur:
                cur.execute("SELECT COUNT(*) FROM pharmaco_db.disease_target_associations")
                total = cur.fetchone()[0]
            log(f"  [{i+1}/{len(genes)}] queried={queried} found={found} assocs={total:,} errors={errors}")

    # Final batch
    if batch:
        with conn.cursor() as cur:
            psycopg2.extras.execute_values(cur,
                """INSERT INTO pharmaco_db.disease_target_associations
                   (target_id, disease_id, disease_name, therapeutic_area,
                    overall_score, genetic_score, somatic_score,
                    known_drug_score, literature_score, rna_expression_score,
                    animal_model_score, source)
                   VALUES %s ON CONFLICT DO NOTHING""",
                batch, page_size=5000)
        conn.commit()

    with conn.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.disease_target_associations")
        total = cur.fetchone()[0]
    conn.close()
    log(f"  DONE: {total:,} disease-target associations (queried={queried}, found={found}, errors={errors})")

if __name__ == "__main__":
    main()
