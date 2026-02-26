#!/usr/bin/env python3
"""Fill disease_target_associations from Open Targets BULK FTP files (200 parts, ~250MB total)."""
import json
import os
import time
import urllib.request
import re
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Set, Tuple

import psycopg2
import psycopg2.extras

DB = {"host": "localhost", "port": 5433, "dbname": "pharmaco", "user": "postgres", "password": "postgres"}
BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/24.09/output/etl/json/associationByOverallDirect/"
DL_DIR = Path("/tmp/opentargets_bulk")

def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

def main():
    log("Open Targets BULK Download")
    DL_DIR.mkdir(parents=True, exist_ok=True)

    conn = psycopg2.connect(**DB)

    # Build Ensembl ID -> target_id map via gene names
    # We need to map Open Targets target IDs (Ensembl) to our target_id
    # Open Targets uses Ensembl gene IDs (ENSG...) as target identifiers
    # We need: gene_name -> target_id, then we'll cross-reference
    with conn.cursor() as cur:
        cur.execute("""SELECT id, gene_name FROM pharmaco_db.targets
                       WHERE gene_name IS NOT NULL AND organism = 'Homo sapiens'""")
        gene_to_tid = {r[1]: r[0] for r in cur.fetchall()}
    log(f"  Gene map: {len(gene_to_tid):,} targets")

    # Get Ensembl -> gene mapping from Open Targets targets file
    targets_base = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/24.09/output/etl/json/targets/"
    log("  Downloading targets mapping (Ensembl -> gene symbol)...")

    # Get target file list
    req = urllib.request.Request(targets_base, headers={"User-Agent": "PharmacoDB/1.0"})
    html = urllib.request.urlopen(req, timeout=15).read().decode()
    target_files = re.findall(r'href="(part-\d+-[a-f0-9-]+-c000\.json)"', html)
    log(f"  Target files: {len(target_files)}")

    ensembl_to_gene = {}
    for fi, fname in enumerate(target_files):
        fpath = DL_DIR / f"targets_{fname}"
        if not fpath.exists():
            url = targets_base + fname
            urllib.request.urlretrieve(url, fpath)

        with open(fpath) as f:
            for line in f:
                try:
                    obj = json.loads(line)
                    eid = obj.get("id", "")
                    symbol = obj.get("approvedSymbol", "")
                    if eid and symbol and symbol in gene_to_tid:
                        ensembl_to_gene[eid] = symbol
                except:
                    continue

        if (fi + 1) % 50 == 0:
            log(f"    Targets: {fi+1}/{len(target_files)} files, {len(ensembl_to_gene):,} mapped")

    log(f"  Ensembl -> gene map: {len(ensembl_to_gene):,} targets")

    # Now download association files
    log("  Getting association file list...")
    req = urllib.request.Request(BASE_URL, headers={"User-Agent": "PharmacoDB/1.0"})
    html = urllib.request.urlopen(req, timeout=15).read().decode()
    assoc_files = re.findall(r'href="(part-\d+-[a-f0-9-]+-c000\.json)"', html)
    log(f"  Association files: {len(assoc_files)}")

    # Also need disease info to get names
    # Open Targets associations have targetId and diseaseId but no names inline
    # We'll need the diseases file too
    diseases_base = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/24.09/output/etl/json/diseases/"
    log("  Downloading disease names...")
    req = urllib.request.Request(diseases_base, headers={"User-Agent": "PharmacoDB/1.0"})
    html = urllib.request.urlopen(req, timeout=15).read().decode()
    disease_files = re.findall(r'href="(part-\d+-[a-f0-9-]+-c000\.json)"', html)

    disease_names = {}
    disease_areas = {}
    for fi, fname in enumerate(disease_files):
        fpath = DL_DIR / f"diseases_{fname}"
        if not fpath.exists():
            url = diseases_base + fname
            urllib.request.urlretrieve(url, fpath)

        with open(fpath) as f:
            for line in f:
                try:
                    obj = json.loads(line)
                    did = obj.get("id", "")
                    dname = obj.get("name", "")
                    tas = obj.get("therapeuticAreas", [])
                    if did and dname:
                        disease_names[did] = dname
                        if tas:
                            disease_areas[did] = tas[0].get("name", "") if isinstance(tas[0], dict) else str(tas[0])
                except:
                    continue

        if (fi + 1) % 20 == 0:
            log(f"    Diseases: {fi+1}/{len(disease_files)} files, {len(disease_names):,} names")

    log(f"  Disease names: {len(disease_names):,}")

    # Clear existing and process associations
    with conn.cursor() as cur:
        cur.execute("TRUNCATE pharmaco_db.disease_target_associations")
    conn.commit()

    total_inserted = 0
    total_skipped = 0
    batch = []

    for fi, fname in enumerate(assoc_files):
        fpath = DL_DIR / f"assoc_{fname}"
        if not fpath.exists():
            url = BASE_URL + fname
            urllib.request.urlretrieve(url, fpath)

        with open(fpath) as f:
            for line in f:
                try:
                    obj = json.loads(line)
                    target_eid = obj.get("targetId", "")
                    disease_id = obj.get("diseaseId", "")
                    overall_score = obj.get("score", 0)

                    gene = ensembl_to_gene.get(target_eid)
                    if not gene:
                        total_skipped += 1
                        continue

                    tid = gene_to_tid.get(gene)
                    if not tid:
                        total_skipped += 1
                        continue

                    dname = disease_names.get(disease_id, "")
                    ta = disease_areas.get(disease_id, "")

                    # Parse datatype scores
                    dtypes = obj.get("datatypeScores", [])
                    scores = {d.get("componentId", d.get("id", "")): d.get("score", 0) for d in dtypes}

                    batch.append((
                        tid, disease_id, dname, ta, overall_score,
                        scores.get("ot_genetics_portal"),
                        scores.get("cancer_gene_census"),
                        scores.get("chembl"),
                        scores.get("europepmc"),
                        scores.get("expression_atlas"),
                        scores.get("impc"),
                        "open_targets",
                    ))
                except:
                    continue

        if len(batch) >= 50000:
            with conn.cursor() as cur:
                psycopg2.extras.execute_values(cur,
                    """INSERT INTO pharmaco_db.disease_target_associations
                       (target_id, disease_id, disease_name, therapeutic_area,
                        overall_score, genetic_score, somatic_score,
                        known_drug_score, literature_score, rna_expression_score,
                        animal_model_score, source)
                       VALUES %s ON CONFLICT DO NOTHING""",
                    batch, page_size=10000)
            conn.commit()
            total_inserted += len(batch)
            batch = []

        if (fi + 1) % 20 == 0:
            log(f"  Associations: {fi+1}/{len(assoc_files)} files, {total_inserted:,} inserted, {total_skipped:,} skipped")

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
                batch, page_size=10000)
        conn.commit()
        total_inserted += len(batch)

    with conn.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.disease_target_associations")
        final = cur.fetchone()[0]
    conn.close()
    log(f"  DONE: {final:,} disease-target associations ({total_inserted:,} inserted, {total_skipped:,} skipped)")

if __name__ == "__main__":
    main()
