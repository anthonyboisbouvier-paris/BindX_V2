#!/usr/bin/env python3
"""
MegaPharmaDB — Comprehensive Consolidation Engine
Mines ChEMBL dump + computes tags + ingests free public sources.
Target: maximize data density in pharmaco_db as fast as possible.

Phase 1: Mine ChEMBL 36 dump (already in local container)
Phase 2: Compute all ML tags (compound, target, activity, selectivity)
Phase 3: Download & ingest free public sources (STRING, SIDER, DGIdb, etc.)
"""

import psycopg2
import psycopg2.extras
import os
import sys
import time
import gzip
import csv
import io
import json
import urllib.request
import urllib.error
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# ============================================================
# CONFIG
# ============================================================
DB_CONFIG = {
    "host": "localhost",
    "port": 5433,
    "dbname": "pharmaco",
    "user": "postgres",
    "password": "pharmaco_secret",
}

CHEMBL_DB = {
    "host": "localhost",
    "port": 5433,
    "dbname": "chembl_36",
    "user": "postgres",
    "password": "pharmaco_secret",
}

DOWNLOAD_DIR = Path("/tmp/megapharma_downloads")
BATCH = 10000

def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

def get_conn(config=None):
    return psycopg2.connect(**(config or DB_CONFIG))

def get_chembl_conn():
    return psycopg2.connect(**CHEMBL_DB)

def download_file(url, dest, desc=""):
    """Download a file with progress."""
    if dest.exists():
        log(f"  Already downloaded: {dest.name}")
        return dest
    log(f"  Downloading {desc or url}...")
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "MegaPharmaDB/1.0"})
        with urllib.request.urlopen(req, timeout=300) as resp:
            total = int(resp.headers.get("Content-Length", 0))
            downloaded = 0
            with open(dest, "wb") as f:
                while True:
                    chunk = resp.read(1024 * 1024)
                    if not chunk:
                        break
                    f.write(chunk)
                    downloaded += len(chunk)
            log(f"  Downloaded {downloaded / 1024 / 1024:.1f} MB")
        return dest
    except Exception as e:
        log(f"  DOWNLOAD FAILED: {e}")
        if dest.exists():
            dest.unlink()
        return None

# ============================================================
# PHASE 1: MINE CHEMBL DUMP
# ============================================================

def mine_ligand_efficiency(pharmaco_conn, chembl_conn):
    """Extract ligand efficiency metrics from ChEMBL dump (2.2M rows).
    Maps to compound_tags: adds LE, BEI, SEI columns."""
    log("  Mining ligand efficiency...")

    # First add columns if missing
    with pharmaco_conn.cursor() as cur:
        for col, typ in [
            ("le", "REAL"),
            ("bei", "REAL"),
            ("sei", "REAL"),
            ("lle", "REAL"),
        ]:
            try:
                cur.execute(f"ALTER TABLE pharmaco_db.compound_tags ADD COLUMN {col} {typ}")
                pharmaco_conn.commit()
            except psycopg2.errors.DuplicateColumn:
                pharmaco_conn.rollback()

    # Build compound chembl_id -> id map
    with pharmaco_conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.compounds WHERE chembl_id IS NOT NULL")
        cmap = dict(cur.fetchall())
    log(f"    Compound map: {len(cmap):,} entries")

    # Read from chembl_36 using server-side cursor
    # ligand_eff joins via activity_id → activities → molecule_dictionary
    read_cur = chembl_conn.cursor("le_cursor")
    read_cur.itersize = BATCH
    read_cur.execute("""
        SELECT DISTINCT ON (md.chembl_id)
               md.chembl_id, le.bei, le.sei, le.le, le.lle
        FROM ligand_eff le
        JOIN activities act ON act.activity_id = le.activity_id
        JOIN molecule_dictionary md ON md.molregno = act.molregno
        ORDER BY md.chembl_id, le.bei DESC NULLS LAST
    """)

    batch = []
    total = 0
    for row in read_cur:
        chembl_id = row[0]
        cid = cmap.get(chembl_id)
        if not cid:
            continue
        batch.append((row[2], row[1], row[3], row[4], cid))  # bei, sei, le, lle, compound_id
        if len(batch) >= BATCH:
            with pharmaco_conn.cursor() as wr:
                psycopg2.extras.execute_batch(wr, """
                    UPDATE pharmaco_db.compound_tags SET
                        bei = %s, sei = %s, le = %s, lle = %s
                    WHERE compound_id = %s
                """, batch, page_size=5000)
            pharmaco_conn.commit()
            total += len(batch)
            batch = []
    if batch:
        with pharmaco_conn.cursor() as wr:
            psycopg2.extras.execute_batch(wr, """
                UPDATE pharmaco_db.compound_tags SET
                    bei = %s, sei = %s, le = %s, lle = %s
                WHERE compound_id = %s
            """, batch, page_size=5000)
        pharmaco_conn.commit()
        total += len(batch)

    read_cur.close()
    log(f"    Ligand efficiency: {total:,} compounds updated")

def mine_structural_alerts(pharmaco_conn, chembl_conn):
    """Extract PAINS + Brenk alert counts from ChEMBL dump (4.9M rows).
    Updates compound_tags.pains_alerts and brenk_alerts."""
    log("  Mining structural alerts...")

    # Get alert set IDs
    with chembl_conn.cursor() as cur:
        cur.execute("SELECT set_name, alert_set_id FROM structural_alert_sets")
        sets = dict(cur.fetchall())
    log(f"    Alert sets: {list(sets.keys())}")
    pains_id = sets.get("PAINS")  # id=4
    brenk_id = sets.get("BMS")   # id=3 (closest to Brenk)
    log(f"    Using: PAINS={pains_id}, BMS/Brenk={brenk_id}")

    # Build compound map
    with pharmaco_conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.compounds WHERE chembl_id IS NOT NULL")
        cmap = dict(cur.fetchall())

    # Count alerts per compound, grouped by set
    read_cur = chembl_conn.cursor("alerts_cursor")
    read_cur.itersize = 50000
    read_cur.execute("""
        SELECT md.chembl_id, sa.alert_set_id, COUNT(*) as cnt
        FROM compound_structural_alerts csa
        JOIN structural_alerts sa ON sa.alert_id = csa.alert_id
        JOIN molecule_dictionary md ON md.molregno = csa.molregno
        GROUP BY md.chembl_id, sa.alert_set_id
    """)

    pains_counts = {}
    brenk_counts = {}
    for row in read_cur:
        chembl_id, set_id, cnt = row
        cid = cmap.get(chembl_id)
        if not cid:
            continue
        if set_id == pains_id:
            pains_counts[cid] = cnt
        elif set_id == brenk_id:
            brenk_counts[cid] = cnt
    read_cur.close()

    # Batch update
    batch = []
    for cid in set(list(pains_counts.keys()) + list(brenk_counts.keys())):
        batch.append((pains_counts.get(cid, 0), brenk_counts.get(cid, 0), cid))
        if len(batch) >= BATCH:
            with pharmaco_conn.cursor() as wr:
                psycopg2.extras.execute_batch(wr, """
                    UPDATE pharmaco_db.compound_tags SET
                        pains_alerts = %s, brenk_alerts = %s,
                        has_reactive_group = (%s > 0 OR %s > 0)
                    WHERE compound_id = %s
                """, [(p, b, p, b, c) for p, b, c in batch], page_size=5000)
            pharmaco_conn.commit()
            batch = []
    if batch:
        with pharmaco_conn.cursor() as wr:
            psycopg2.extras.execute_batch(wr, """
                UPDATE pharmaco_db.compound_tags SET
                    pains_alerts = %s, brenk_alerts = %s,
                    has_reactive_group = (%s > 0 OR %s > 0)
                WHERE compound_id = %s
            """, [(p, b, p, b, c) for p, b, c in batch], page_size=5000)
        pharmaco_conn.commit()

    log(f"    PAINS alerts: {len(pains_counts):,} compounds")
    log(f"    Brenk alerts: {len(brenk_counts):,} compounds")

def mine_atc_classification(pharmaco_conn, chembl_conn):
    """Extract ATC codes from ChEMBL dump."""
    log("  Mining ATC classification...")

    with pharmaco_conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.compounds WHERE chembl_id IS NOT NULL")
        cmap = dict(cur.fetchall())

    with chembl_conn.cursor() as cur:
        cur.execute("""
            SELECT md.chembl_id,
                   ac.level1, ac.level2, ac.level3, ac.level4, ac.level5,
                   ac.level1_description, ac.level2_description,
                   ac.level3_description, ac.level4_description,
                   ac.who_name
            FROM molecule_atc_classification mac
            JOIN molecule_dictionary md ON md.molregno = mac.molregno
            JOIN atc_classification ac ON ac.level5 = mac.level5
        """)
        rows = cur.fetchall()

    batch = []
    for row in rows:
        chembl_id = row[0]
        cid = cmap.get(chembl_id)
        if not cid:
            continue
        atc_code = row[5] or row[4] or row[3] or row[2] or row[1]
        batch.append((
            cid, atc_code,
            row[1], row[6],   # l1, l1_desc
            row[2], row[7],   # l2, l2_desc
            row[3], row[8],   # l3, l3_desc
            row[4], row[9],   # l4, l4_desc
            row[5], None,     # l5, l5_desc
            row[10],          # who_name
        ))

    if batch:
        with pharmaco_conn.cursor() as wr:
            # Truncate first
            wr.execute("TRUNCATE pharmaco_db.atc_classification")
            psycopg2.extras.execute_values(wr, """
                INSERT INTO pharmaco_db.atc_classification
                (compound_id, atc_code, level1, level1_desc, level2, level2_desc,
                 level3, level3_desc, level4, level4_desc, level5, level5_desc, who_name)
                VALUES %s
            """, batch, page_size=5000)
        pharmaco_conn.commit()
    log(f"    ATC classification: {len(batch):,} entries")

def mine_binding_sites(pharmaco_conn, chembl_conn):
    """Extract binding site information from ChEMBL dump."""
    log("  Mining binding sites...")

    # Add binding_sites table extension if needed
    with pharmaco_conn.cursor() as cur:
        cur.execute("""
            CREATE TABLE IF NOT EXISTS pharmaco_db.binding_sites (
                id BIGSERIAL PRIMARY KEY,
                target_id BIGINT REFERENCES pharmaco_db.targets(id),
                site_name TEXT,
                site_id INTEGER,
                domain_type TEXT,
                domain_name TEXT,
                domain_description TEXT,
                start_position INTEGER,
                end_position INTEGER,
                component_type TEXT,
                source TEXT DEFAULT 'chembl',
                created_at TIMESTAMPTZ DEFAULT NOW()
            )
        """)
        pharmaco_conn.commit()

    # Map targets
    with pharmaco_conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.targets WHERE chembl_id IS NOT NULL")
        tmap = dict(cur.fetchall())

    with chembl_conn.cursor() as cur:
        cur.execute("""
            SELECT td.chembl_id as target_chembl_id,
                   bs.site_name, bs.site_id,
                   d.domain_type, d.domain_name, d.domain_description,
                   cd.start_position, cd.end_position,
                   'protein' as component_type
            FROM binding_sites bs
            JOIN target_dictionary td ON td.tid = bs.tid
            LEFT JOIN site_components sc ON sc.site_id = bs.site_id
            LEFT JOIN domains d ON d.domain_id = sc.domain_id
            LEFT JOIN component_domains cd ON cd.domain_id = sc.domain_id
                AND cd.component_id = sc.component_id
        """)
        rows = cur.fetchall()

    batch = []
    for row in rows:
        tid = tmap.get(row[0])
        if not tid:
            continue
        batch.append((tid, row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8]))

    if batch:
        with pharmaco_conn.cursor() as wr:
            wr.execute("TRUNCATE pharmaco_db.binding_sites CASCADE")
            psycopg2.extras.execute_values(wr, """
                INSERT INTO pharmaco_db.binding_sites
                (target_id, site_name, site_id, domain_type, domain_name,
                 domain_description, start_position, end_position, component_type)
                VALUES %s
            """, batch, page_size=5000)
        pharmaco_conn.commit()
    log(f"    Binding sites: {len(batch):,} entries")

def mine_target_relations(pharmaco_conn, chembl_conn):
    """Extract target-target relationships from ChEMBL dump (98K rows)."""
    log("  Mining target relations...")

    with pharmaco_conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.targets WHERE chembl_id IS NOT NULL")
        tmap = dict(cur.fetchall())

    with chembl_conn.cursor() as cur:
        cur.execute("""
            SELECT td1.chembl_id, td2.chembl_id, tr.relationship
            FROM target_relations tr
            JOIN target_dictionary td1 ON td1.tid = tr.tid
            JOIN target_dictionary td2 ON td2.tid = tr.related_tid
        """)
        rows = cur.fetchall()

    batch = []
    for row in rows:
        tid1 = tmap.get(row[0])
        tid2 = tmap.get(row[1])
        if not tid1 or not tid2 or tid1 == tid2:
            continue
        batch.append((tid1, tid2, row[2], "chembl", 1.0))

    if batch:
        with pharmaco_conn.cursor() as wr:
            wr.execute("DELETE FROM pharmaco_db.target_interactions WHERE source = 'chembl'")
            psycopg2.extras.execute_values(wr, """
                INSERT INTO pharmaco_db.target_interactions
                (target_id_1, target_id_2, interaction_type, source, confidence)
                VALUES %s
                ON CONFLICT (target_id_1, target_id_2, interaction_type) DO NOTHING
            """, batch, page_size=5000)
        pharmaco_conn.commit()
    log(f"    Target relations: {len(batch):,} pairs")

def mine_predicted_binding_domains(pharmaco_conn, chembl_conn):
    """Extract predicted binding domains for bioactivities."""
    log("  Mining predicted binding domains...")

    # Add column to bioactivities for domain info
    with pharmaco_conn.cursor() as cur:
        for col, typ in [
            ("predicted_domain", "TEXT"),
            ("predicted_domain_type", "TEXT"),
        ]:
            try:
                cur.execute(f"ALTER TABLE pharmaco_db.bioactivities ADD COLUMN {col} {typ}")
                pharmaco_conn.commit()
            except psycopg2.errors.DuplicateColumn:
                pharmaco_conn.rollback()

    # Build assay_chembl_id -> pharmaco assay_id map
    with pharmaco_conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.assays WHERE chembl_id IS NOT NULL")
        amap = dict(cur.fetchall())

    # Read domain predictions (sample - too many for full join)
    read_cur = chembl_conn.cursor("domain_cursor")
    read_cur.itersize = 50000
    read_cur.execute("""
        SELECT DISTINCT ON (a.chembl_id)
               a.chembl_id as assay_chembl_id,
               cd.domain_name, cd.domain_type
        FROM predicted_binding_domains pbd
        JOIN binding_sites bs ON bs.site_id = pbd.site_id
        JOIN site_components sc ON sc.site_id = bs.site_id
        LEFT JOIN component_domains cd ON cd.compd_id = sc.domain_id
        JOIN assays a ON a.assay_id = pbd.activity_id  -- This is actually assay_id in chembl
        WHERE cd.domain_name IS NOT NULL
        LIMIT 500000
    """)

    # Just store the unique assay -> domain mapping
    domain_map = {}
    for row in read_cur:
        assay_chembl_id, domain_name, domain_type = row
        aid = amap.get(assay_chembl_id)
        if aid:
            domain_map[aid] = (domain_name, domain_type)
    read_cur.close()

    log(f"    Found {len(domain_map):,} assay-domain mappings (skipping bulk update for speed)")

def mine_go_terms_from_chembl(pharmaco_conn, chembl_conn):
    """Extract GO terms from ChEMBL component_go table for targets missing UniProt enrichment."""
    log("  Mining GO terms from ChEMBL...")

    with pharmaco_conn.cursor() as cur:
        cur.execute("""
            SELECT chembl_id, id FROM pharmaco_db.targets
            WHERE go_molecular_function IS NULL AND go_biological_process IS NULL
        """)
        targets_missing_go = dict(cur.fetchall())
    log(f"    Targets missing GO terms: {len(targets_missing_go):,}")

    if not targets_missing_go:
        return

    with chembl_conn.cursor() as cur:
        cur.execute("""
            SELECT td.chembl_id, gc.go_id, gc.aspect
            FROM component_go cg
            JOIN go_classification gc ON gc.go_id = cg.go_id
            JOIN target_components tc ON tc.component_id = cg.component_id
            JOIN target_dictionary td ON td.tid = tc.tid
            WHERE gc.aspect IS NOT NULL
        """)
        rows = cur.fetchall()

    # Group by target
    go_by_target = {}
    for chembl_id, go_id, go_type in rows:
        tid = targets_missing_go.get(chembl_id)
        if not tid:
            continue
        if tid not in go_by_target:
            go_by_target[tid] = {"F": [], "P": [], "C": []}
        if go_type in go_by_target[tid]:
            go_by_target[tid][go_type].append(go_id)

    batch = []
    for tid, terms in go_by_target.items():
        batch.append((terms.get("F", []) or None, terms.get("P", []) or None, terms.get("C", []) or None, tid))

    if batch:
        with pharmaco_conn.cursor() as wr:
            psycopg2.extras.execute_batch(wr, """
                UPDATE pharmaco_db.targets SET
                    go_molecular_function = COALESCE(go_molecular_function, %s),
                    go_biological_process = COALESCE(go_biological_process, %s),
                    go_cellular_component = COALESCE(go_cellular_component, %s)
                WHERE id = %s
            """, batch, page_size=5000)
        pharmaco_conn.commit()
    log(f"    Updated GO terms for {len(batch):,} targets")

def mine_molecule_synonyms(pharmaco_conn, chembl_conn):
    """Extract molecule synonyms from ChEMBL dump (132K names)."""
    log("  Mining molecule synonyms...")

    with pharmaco_conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.compounds WHERE chembl_id IS NOT NULL AND pref_name IS NULL")
        cmap = dict(cur.fetchall())
    log(f"    Compounds without names: {len(cmap):,}")

    if not cmap:
        return

    with chembl_conn.cursor() as cur:
        cur.execute("""
            SELECT md.chembl_id, ms.synonyms, ms.syn_type
            FROM molecule_synonyms ms
            JOIN molecule_dictionary md ON md.molregno = ms.molregno
            WHERE ms.syn_type IN ('INN', 'USAN', 'BAN', 'TRADE_NAME', 'OTHER')
            ORDER BY md.chembl_id,
                CASE ms.syn_type
                    WHEN 'INN' THEN 1 WHEN 'USAN' THEN 2
                    WHEN 'BAN' THEN 3 WHEN 'TRADE_NAME' THEN 4
                    ELSE 5 END
        """)
        rows = cur.fetchall()

    # First synonym per compound
    seen = set()
    batch = []
    for chembl_id, name, syn_type in rows:
        if chembl_id in seen:
            continue
        cid = cmap.get(chembl_id)
        if not cid:
            continue
        seen.add(chembl_id)
        batch.append((name, cid))

    if batch:
        with pharmaco_conn.cursor() as wr:
            psycopg2.extras.execute_batch(wr, """
                UPDATE pharmaco_db.compounds SET pref_name = %s WHERE id = %s
            """, batch, page_size=5000)
        pharmaco_conn.commit()
    log(f"    Updated names for {len(batch):,} compounds")

def mine_assay_publications(pharmaco_conn, chembl_conn):
    """Extract journal/year/DOI from ChEMBL docs table for assays."""
    log("  Mining assay publications...")

    with pharmaco_conn.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.assays WHERE journal IS NULL")
        missing = cur.fetchone()[0]
    log(f"    Assays without publications: {missing:,}")

    if missing == 0:
        return

    # Build assay map
    with pharmaco_conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.assays WHERE chembl_id IS NOT NULL AND journal IS NULL")
        amap = dict(cur.fetchall())

    read_cur = chembl_conn.cursor("docs_cursor")
    read_cur.itersize = 50000
    read_cur.execute("""
        SELECT a.chembl_id as assay_chembl_id, d.journal, d.year, d.doi
        FROM assays a
        JOIN docs d ON d.doc_id = a.doc_id
        WHERE d.journal IS NOT NULL
    """)

    batch = []
    total = 0
    for row in read_cur:
        assay_chembl_id, journal, year, doi = row
        aid = amap.get(assay_chembl_id)
        if not aid:
            continue
        batch.append((journal, year, doi, aid))
        if len(batch) >= BATCH:
            with pharmaco_conn.cursor() as wr:
                psycopg2.extras.execute_batch(wr, """
                    UPDATE pharmaco_db.assays SET journal = %s, year = %s, doi = %s WHERE id = %s
                """, batch, page_size=5000)
            pharmaco_conn.commit()
            total += len(batch)
            batch = []
    if batch:
        with pharmaco_conn.cursor() as wr:
            psycopg2.extras.execute_batch(wr, """
                UPDATE pharmaco_db.assays SET journal = %s, year = %s, doi = %s WHERE id = %s
            """, batch, page_size=5000)
        pharmaco_conn.commit()
        total += len(batch)
    read_cur.close()
    log(f"    Updated publications for {total:,} assays")


def run_phase1_chembl_mining():
    """Phase 1: Mine everything useful from the ChEMBL dump."""
    log("=" * 60)
    log("PHASE 1: Mining ChEMBL 36 dump")
    log("=" * 60)

    pharmaco_conn = get_conn()
    chembl_conn = get_chembl_conn()

    try:
        # First ensure compound_tags has records for all compounds
        with pharmaco_conn.cursor() as cur:
            cur.execute("""
                INSERT INTO pharmaco_db.compound_tags (compound_id)
                SELECT id FROM pharmaco_db.compounds
                ON CONFLICT (compound_id) DO NOTHING
            """)
            pharmaco_conn.commit()
            log(f"  Ensured compound_tags: {cur.rowcount} new records")

        mine_molecule_synonyms(pharmaco_conn, chembl_conn)
        mine_atc_classification(pharmaco_conn, chembl_conn)
        mine_structural_alerts(pharmaco_conn, chembl_conn)
        mine_ligand_efficiency(pharmaco_conn, chembl_conn)
        mine_binding_sites(pharmaco_conn, chembl_conn)
        mine_target_relations(pharmaco_conn, chembl_conn)
        mine_go_terms_from_chembl(pharmaco_conn, chembl_conn)
        mine_assay_publications(pharmaco_conn, chembl_conn)

    finally:
        pharmaco_conn.close()
        chembl_conn.close()

    log("Phase 1 COMPLETE")


# ============================================================
# PHASE 2: COMPUTE TAGS
# ============================================================

def run_phase2_compute_tags():
    """Phase 2: Compute all ML-relevant tags."""
    log("=" * 60)
    log("PHASE 2: Computing ML Tags")
    log("=" * 60)

    conn = get_conn()
    try:
        compute_compound_tags(conn)
        compute_target_tags(conn)
        compute_activity_tags(conn)
        compute_selectivity(conn)
    finally:
        conn.close()

    log("Phase 2 COMPLETE")


def compute_compound_tags(conn):
    log("  Computing compound tags...")
    with conn.cursor() as cur:
        # Ensure records exist
        cur.execute("""
            INSERT INTO pharmaco_db.compound_tags (compound_id)
            SELECT id FROM pharmaco_db.compounds
            ON CONFLICT (compound_id) DO NOTHING
        """)
        conn.commit()

        # Drug-likeness rules
        cur.execute("""
            UPDATE pharmaco_db.compound_tags ct SET
                lipinski_pass = (c.molecular_weight <= 500 AND c.alogp <= 5
                    AND c.hbd <= 5 AND c.hba <= 10),
                lipinski_violations = c.num_ro5_violations,
                veber_pass = (c.psa <= 140 AND c.rtb <= 10),
                lead_like = (c.molecular_weight BETWEEN 200 AND 350
                    AND c.alogp BETWEEN -1 AND 3
                    AND COALESCE(c.hbd, 0) <= 3 AND COALESCE(c.hba, 0) <= 6),
                fragment_like = (c.molecular_weight <= 300 AND c.alogp <= 3
                    AND COALESCE(c.hbd, 0) <= 3 AND COALESCE(c.hba, 0) <= 3
                    AND COALESCE(c.rtb, 0) <= 3),
                drug_like = (COALESCE(c.qed_weighted, 0) >= 0.5),
                ppi_inhibitor_like = (c.molecular_weight > 400 AND c.alogp > 4
                    AND COALESCE(c.aromatic_rings, 0) >= 4),
                clinical_phase = c.max_phase,
                is_approved = (c.max_phase >= 4),
                updated_at = NOW()
            FROM pharmaco_db.compounds c
            WHERE c.id = ct.compound_id
        """)
        conn.commit()
        log(f"    Drug-likeness rules: {cur.rowcount:,} rows")

        # Physicochemical bins
        cur.execute("""
            UPDATE pharmaco_db.compound_tags ct SET
                mw_bin = CASE
                    WHEN c.molecular_weight < 200 THEN '<200'
                    WHEN c.molecular_weight BETWEEN 200 AND 350 THEN '200-350'
                    WHEN c.molecular_weight BETWEEN 350 AND 500 THEN '350-500'
                    ELSE '>500' END,
                logp_bin = CASE
                    WHEN c.alogp < 0 THEN '<0'
                    WHEN c.alogp BETWEEN 0 AND 2 THEN '0-2'
                    WHEN c.alogp BETWEEN 2 AND 4 THEN '2-4'
                    ELSE '>4' END,
                tpsa_bin = CASE
                    WHEN c.psa < 60 THEN '<60'
                    WHEN c.psa BETWEEN 60 AND 90 THEN '60-90'
                    WHEN c.psa BETWEEN 90 AND 140 THEN '90-140'
                    ELSE '>140' END,
                complexity_bin = CASE
                    WHEN c.heavy_atoms < 15 THEN 'low'
                    WHEN c.heavy_atoms BETWEEN 15 AND 25 THEN 'medium'
                    WHEN c.heavy_atoms BETWEEN 25 AND 40 THEN 'high'
                    ELSE 'very_high' END
            FROM pharmaco_db.compounds c
            WHERE c.id = ct.compound_id
        """)
        conn.commit()
        log(f"    Physicochemical bins: {cur.rowcount:,} rows")

        # CNS MPO
        cur.execute("""
            UPDATE pharmaco_db.compound_tags ct SET
                cns_mpo_score = (
                    CASE WHEN c.molecular_weight <= 360 THEN 1.0
                         WHEN c.molecular_weight <= 500 THEN 1.0 - (c.molecular_weight - 360) / 140
                         ELSE 0 END
                    + CASE WHEN c.alogp BETWEEN 1 AND 3 THEN 1.0
                           WHEN c.alogp < 1 THEN GREATEST(0, 1.0 - (1 - c.alogp))
                           WHEN c.alogp > 3 THEN GREATEST(0, 1.0 - (c.alogp - 3) / 2)
                           ELSE 0 END
                    + CASE WHEN c.psa BETWEEN 40 AND 90 THEN 1.0
                           WHEN c.psa < 40 THEN GREATEST(0, c.psa / 40)
                           WHEN c.psa > 90 THEN GREATEST(0, 1.0 - (c.psa - 90) / 50)
                           ELSE 0 END
                    + CASE WHEN COALESCE(c.hbd, 0) <= 1 THEN 1.0
                           WHEN c.hbd = 2 THEN 0.75
                           WHEN c.hbd = 3 THEN 0.25
                           ELSE 0 END
                    + CASE WHEN COALESCE(c.hba, 0) <= 5 THEN 1.0
                           ELSE GREATEST(0, 1.0 - (c.hba - 5) / 5.0) END
                    + CASE WHEN COALESCE(c.aromatic_rings, 0) <= 2 THEN 1.0
                           WHEN c.aromatic_rings = 3 THEN 0.5
                           ELSE 0 END
                ),
                cns_penetrant = (
                    CASE WHEN c.molecular_weight <= 360 THEN 1.0
                         WHEN c.molecular_weight <= 500 THEN 1.0 - (c.molecular_weight - 360) / 140
                         ELSE 0 END
                    + CASE WHEN c.alogp BETWEEN 1 AND 3 THEN 1.0
                           WHEN c.alogp < 1 THEN GREATEST(0, 1.0 - (1 - c.alogp))
                           WHEN c.alogp > 3 THEN GREATEST(0, 1.0 - (c.alogp - 3) / 2)
                           ELSE 0 END
                    + CASE WHEN c.psa BETWEEN 40 AND 90 THEN 1.0
                           WHEN c.psa < 40 THEN GREATEST(0, c.psa / 40)
                           WHEN c.psa > 90 THEN GREATEST(0, 1.0 - (c.psa - 90) / 50)
                           ELSE 0 END
                    + CASE WHEN COALESCE(c.hbd, 0) <= 1 THEN 1.0
                           WHEN c.hbd = 2 THEN 0.75
                           WHEN c.hbd = 3 THEN 0.25
                           ELSE 0 END
                    + CASE WHEN COALESCE(c.hba, 0) <= 5 THEN 1.0
                           ELSE GREATEST(0, 1.0 - (c.hba - 5) / 5.0) END
                    + CASE WHEN COALESCE(c.aromatic_rings, 0) <= 2 THEN 1.0
                           WHEN c.aromatic_rings = 3 THEN 0.5
                           ELSE 0 END
                ) >= 4.0,
                updated_at = NOW()
            FROM pharmaco_db.compounds c
            WHERE c.id = ct.compound_id AND c.molecular_weight IS NOT NULL AND c.alogp IS NOT NULL
        """)
        conn.commit()
        log(f"    CNS MPO: {cur.rowcount:,} rows")

        # Data richness
        cur.execute("""
            UPDATE pharmaco_db.compound_tags ct SET
                num_bioactivities = sub.n_act,
                num_targets = sub.n_tgt,
                num_assays = sub.n_asy
            FROM (
                SELECT compound_id,
                       COUNT(*) as n_act,
                       COUNT(DISTINCT target_id) as n_tgt,
                       COUNT(DISTINCT assay_id) as n_asy
                FROM pharmaco_db.bioactivities
                GROUP BY compound_id
            ) sub
            WHERE sub.compound_id = ct.compound_id
        """)
        conn.commit()
        log(f"    Data richness: {cur.rowcount:,} rows")

        # Data completeness
        cur.execute("""
            UPDATE pharmaco_db.compound_tags ct SET
                data_completeness = (
                    0.1 + CASE WHEN c.molecular_weight IS NOT NULL THEN 0.2 ELSE 0 END
                    + CASE WHEN ct.num_bioactivities > 0 THEN 0.3 ELSE 0 END
                    + CASE WHEN ct.num_targets > 0 THEN 0.2 ELSE 0 END
                    + CASE WHEN c.max_phase > 0 THEN 0.2 ELSE 0 END
                )
            FROM pharmaco_db.compounds c
            WHERE c.id = ct.compound_id
        """)
        conn.commit()
        log(f"    Data completeness: {cur.rowcount:,} rows")

    log("  Compound tags DONE")

def compute_target_tags(conn):
    log("  Computing target tags...")
    with conn.cursor() as cur:
        cur.execute("""
            INSERT INTO pharmaco_db.target_tags (target_id)
            SELECT id FROM pharmaco_db.targets
            ON CONFLICT (target_id) DO NOTHING
        """)
        conn.commit()

        # Classification
        cur.execute("""
            UPDATE pharmaco_db.target_tags tt SET
                target_class = CASE
                    WHEN t.protein_class_l1 = 'Enzyme' AND t.protein_class_l2 LIKE '%%kinase%%' THEN 'kinase'
                    WHEN t.protein_class_l1 = 'Enzyme' AND t.protein_class_l2 LIKE '%%protease%%' THEN 'protease'
                    WHEN t.protein_class_l1 = 'Enzyme' AND t.protein_class_l2 LIKE '%%phosphodiesterase%%' THEN 'phosphodiesterase'
                    WHEN t.protein_class_l1 = 'Enzyme' THEN 'enzyme_other'
                    WHEN t.protein_family = 'G-protein coupled receptor' OR t.protein_class_l1 = 'Membrane receptor' THEN 'gpcr'
                    WHEN t.protein_family = 'Ion channel' OR t.protein_class_l1 = 'Ion channel' THEN 'ion_channel'
                    WHEN t.protein_family = 'Nuclear receptor' THEN 'nuclear_receptor'
                    WHEN t.protein_family = 'Transporter' OR t.protein_class_l1 = 'Transporter' THEN 'transporter'
                    WHEN t.protein_class_l1 = 'Epigenetic regulator' THEN 'epigenetic'
                    WHEN t.protein_class_l1 = 'Transcription factor' THEN 'transcription_factor'
                    WHEN t.protein_class_l1 IS NOT NULL THEN LOWER(REPLACE(t.protein_class_l1, ' ', '_'))
                    WHEN t.protein_family IS NOT NULL THEN LOWER(REPLACE(t.protein_family, ' ', '_'))
                    ELSE 'unclassified'
                END,
                target_subclass = COALESCE(t.protein_class_l2, t.protein_class_l3),
                enzyme_class_ec = t.ec_number,
                has_crystal_structure = (t.pdb_ids IS NOT NULL AND array_length(t.pdb_ids, 1) > 0),
                num_pdb_structures = COALESCE(array_length(t.pdb_ids, 1), 0),
                tissue_expression = t.tissue_expression,
                tissue_specificity = t.tissue_specificity,
                updated_at = NOW()
            FROM pharmaco_db.targets t
            WHERE t.id = tt.target_id
        """)
        conn.commit()
        log(f"    Target classification: {cur.rowcount:,} rows")

        # Activity stats
        cur.execute("""
            UPDATE pharmaco_db.target_tags tt SET
                num_active_compounds = sub.n_active,
                num_total_activities = sub.n_total,
                num_approved_drugs = COALESCE(t.num_approved_drugs, 0),
                median_pchembl = sub.med_pchembl,
                data_maturity = CASE
                    WHEN sub.n_active >= 100 THEN 'mature'
                    WHEN sub.n_active >= 10 THEN 'emerging'
                    ELSE 'sparse' END
            FROM (
                SELECT target_id,
                       COUNT(*) FILTER (WHERE activity_class = 'active') as n_active,
                       COUNT(*) as n_total,
                       PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY pchembl_value) as med_pchembl
                FROM pharmaco_db.bioactivities
                WHERE target_id IS NOT NULL AND pchembl_value IS NOT NULL
                GROUP BY target_id
            ) sub
            JOIN pharmaco_db.targets t ON t.id = sub.target_id
            WHERE sub.target_id = tt.target_id
        """)
        conn.commit()
        log(f"    Activity stats: {cur.rowcount:,} rows")

        # Druggability score
        cur.execute("""
            UPDATE pharmaco_db.target_tags tt SET
                druggability_score = (
                    CASE WHEN tt.has_crystal_structure THEN 0.2 ELSE 0 END
                    + CASE WHEN t.num_approved_drugs > 0 THEN 0.3 ELSE 0 END
                    + CASE WHEN tt.num_active_compounds > 10 THEN 0.2
                           WHEN tt.num_active_compounds > 0 THEN 0.1 ELSE 0 END
                    + CASE WHEN tt.target_class IN ('kinase','gpcr','ion_channel',
                            'nuclear_receptor','protease','phosphodiesterase') THEN 0.2 ELSE 0 END
                    + CASE WHEN tt.data_maturity = 'mature' THEN 0.1
                           WHEN tt.data_maturity = 'emerging' THEN 0.05 ELSE 0 END
                ),
                updated_at = NOW()
            FROM pharmaco_db.targets t
            WHERE t.id = tt.target_id
        """)
        conn.commit()
        log(f"    Druggability scores: {cur.rowcount:,} rows")

        # Family size
        cur.execute("""
            UPDATE pharmaco_db.target_tags tt SET
                family_size = sub.fam_size,
                selectivity_challenge = CASE
                    WHEN sub.fam_size > 50 THEN 'hard'
                    WHEN sub.fam_size > 10 THEN 'moderate'
                    ELSE 'easy' END
            FROM (
                SELECT t.id, COUNT(*) OVER (PARTITION BY t.protein_family) as fam_size
                FROM pharmaco_db.targets t WHERE t.protein_family IS NOT NULL
            ) sub
            WHERE sub.id = tt.target_id
        """)
        conn.commit()
        log(f"    Family sizes: {cur.rowcount:,} rows")

    log("  Target tags DONE")

def compute_activity_tags(conn):
    log("  Computing activity tags (batched)...")
    with conn.cursor() as cur:
        cur.execute("SELECT MAX(id) FROM pharmaco_db.bioactivities")
        max_id = cur.fetchone()[0] or 0
        log(f"    Max bioactivity ID: {max_id:,}")

    batch_size = 500000
    offset = 0
    total = 0

    while offset < max_id:
        with conn.cursor() as cur:
            cur.execute("""
                INSERT INTO pharmaco_db.activity_tags (activity_id, confidence_level, pchembl_bin,
                    is_direct_binding, is_cell_based)
                SELECT
                    b.id,
                    CASE
                        WHEN a.confidence_score >= 8 THEN 'high'
                        WHEN a.confidence_score >= 6 THEN 'medium'
                        ELSE 'low' END,
                    CASE
                        WHEN b.pchembl_value >= 8 THEN 'very_potent(>8)'
                        WHEN b.pchembl_value >= 7 THEN 'potent(7-8)'
                        WHEN b.pchembl_value >= 6 THEN 'moderate(6-7)'
                        WHEN b.pchembl_value >= 5 THEN 'weak(5-6)'
                        ELSE 'inactive(<5)' END,
                    (a.assay_type = 'B'),
                    (a.assay_type = 'F')
                FROM pharmaco_db.bioactivities b
                LEFT JOIN pharmaco_db.assays a ON a.id = b.assay_id
                WHERE b.id > %s AND b.id <= %s
                AND b.pchembl_value IS NOT NULL
                ON CONFLICT (activity_id) DO UPDATE SET
                    confidence_level = EXCLUDED.confidence_level,
                    pchembl_bin = EXCLUDED.pchembl_bin,
                    is_direct_binding = EXCLUDED.is_direct_binding,
                    is_cell_based = EXCLUDED.is_cell_based
            """, (offset, offset + batch_size))
            conn.commit()
            total += cur.rowcount
        offset += batch_size
        log(f"    Progress: {min(offset, max_id):,}/{max_id:,} ({total:,} tagged)")

    log(f"  Activity tags DONE: {total:,} total")

def compute_selectivity(conn):
    log("  Computing selectivity profiles...")
    with conn.cursor() as cur:
        cur.execute("TRUNCATE pharmaco_db.selectivity_profiles")
        cur.execute("""
            INSERT INTO pharmaco_db.selectivity_profiles
            (compound_id, primary_target_id, off_target_id, primary_pchembl, off_target_pchembl, selectivity_index)
            SELECT
                b1.compound_id,
                b1.target_id,
                b2.target_id,
                MAX(b1.pchembl_value),
                MAX(b2.pchembl_value),
                MAX(b1.pchembl_value) - MAX(b2.pchembl_value)
            FROM pharmaco_db.bioactivities b1
            JOIN pharmaco_db.bioactivities b2 ON b2.compound_id = b1.compound_id
                AND b2.target_id != b1.target_id AND b2.target_id IS NOT NULL
            WHERE b1.target_id IS NOT NULL
                AND b1.pchembl_value >= 6.0 AND b2.pchembl_value IS NOT NULL
            GROUP BY b1.compound_id, b1.target_id, b2.target_id
            LIMIT 500000
        """)
        conn.commit()
        log(f"    Selectivity pairs: {cur.rowcount:,}")

    log("  Selectivity profiles DONE")


# ============================================================
# PHASE 3: FREE PUBLIC DATA SOURCES
# ============================================================

def ingest_string_ppi(conn):
    """Download and ingest STRING protein-protein interactions.
    Schema: protein_interactions(target_id_1, target_id_2, gene_1, gene_2,
            combined_score, experimental_score, database_score, textmining_score,
            coexpression_score, source)
    """
    log("  Ingesting STRING PPI data...")

    # Use detailed links for sub-scores
    url = "https://stringdb-downloads.org/download/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.txt.gz"
    dest = DOWNLOAD_DIR / "string_9606_links_detailed.txt.gz"
    f = download_file(url, dest, "STRING 9606 detailed links")
    if not f:
        # Fallback to simple links
        url = "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz"
        dest = DOWNLOAD_DIR / "string_9606_links.txt.gz"
        f = download_file(url, dest, "STRING 9606 links")
    if not f:
        log("    SKIP: STRING download failed")
        return

    # Build gene_name -> target_id map
    with conn.cursor() as cur:
        cur.execute("SELECT gene_name, id FROM pharmaco_db.targets WHERE gene_name IS NOT NULL")
        gene_map = dict(cur.fetchall())

    # STRING ID -> gene_name mapping via info file
    url_info = "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz"
    dest_info = DOWNLOAD_DIR / "string_9606_info.txt.gz"
    f_info = download_file(url_info, dest_info, "STRING 9606 info")

    string_to_gene = {}
    if f_info and f_info.exists():
        with gzip.open(f_info, "rt") as fh:
            header = next(fh)
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    string_to_gene[parts[0]] = parts[1]

    log(f"    STRING->gene mapping: {len(string_to_gene):,}")

    # Parse links (detailed format has sub-scores)
    interactions = []
    detailed = "detailed" in str(dest)
    with gzip.open(f, "rt") as fh:
        header = next(fh)
        for line in fh:
            parts = line.strip().split()
            if len(parts) < 3:
                continue

            p1, p2 = parts[0], parts[1]
            if detailed and len(parts) >= 10:
                # Columns: protein1 protein2 neighborhood fusion cooccurence coexpression
                #          experimental database textmining combined_score
                combined = int(parts[-1])
                experimental = int(parts[6]) if len(parts) > 6 else 0
                database = int(parts[7]) if len(parts) > 7 else 0
                textmining = int(parts[8]) if len(parts) > 8 else 0
                coexpression = int(parts[5]) if len(parts) > 5 else 0
            else:
                combined = int(parts[2])
                experimental = database = textmining = coexpression = 0

            if combined < 700:
                continue

            g1 = string_to_gene.get(p1)
            g2 = string_to_gene.get(p2)
            if not g1 or not g2:
                continue

            tid1 = gene_map.get(g1)
            tid2 = gene_map.get(g2)
            if not tid1 or not tid2 or tid1 == tid2:
                continue

            # Ensure gene_1 < gene_2 for unique constraint
            if g1 > g2:
                g1, g2 = g2, g1
                tid1, tid2 = tid2, tid1

            interactions.append((tid1, tid2, g1, g2, combined, experimental, database, textmining, coexpression))

    if interactions:
        # Deduplicate
        seen = set()
        unique = []
        for row in interactions:
            key = (row[2], row[3])  # gene_1, gene_2
            if key not in seen:
                seen.add(key)
                unique.append(row)

        with conn.cursor() as cur:
            cur.execute("DELETE FROM pharmaco_db.protein_interactions WHERE source = 'string'")
            psycopg2.extras.execute_values(cur, """
                INSERT INTO pharmaco_db.protein_interactions
                (target_id_1, target_id_2, gene_1, gene_2,
                 combined_score, experimental_score, database_score,
                 textmining_score, coexpression_score, source)
                VALUES %s
                ON CONFLICT DO NOTHING
            """, [(t1, t2, g1, g2, cs, es, ds, ts, ces, "string")
                  for t1, t2, g1, g2, cs, es, ds, ts, ces in unique],
            page_size=5000)
        conn.commit()
        log(f"    STRING PPI: {len(unique):,} high-confidence interactions")
    else:
        log("    STRING PPI: 0 mapped")

def ingest_dgidb(conn):
    """Download and ingest DGIdb drug-gene interactions.
    Schema: dgidb_interactions(compound_id, target_id, gene_name, drug_name,
            drug_chembl_id, interaction_type, interaction_score, sources, pmids)
    """
    log("  Ingesting DGIdb interactions...")

    url = "https://www.dgidb.org/data/monthly_tsvs/2024-Feb/interactions.tsv"
    dest = DOWNLOAD_DIR / "dgidb_interactions.tsv"
    f = download_file(url, dest, "DGIdb interactions")
    if not f:
        url = "https://old.dgidb.org/data/monthly_tsvs/2022-Feb/interactions.tsv"
        f = download_file(url, dest, "DGIdb interactions (alt)")
    if not f:
        log("    SKIP: DGIdb download failed")
        return

    with conn.cursor() as cur:
        cur.execute("SELECT gene_name, id FROM pharmaco_db.targets WHERE gene_name IS NOT NULL")
        gene_map = dict(cur.fetchall())
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.compounds WHERE chembl_id IS NOT NULL")
        chembl_map = dict(cur.fetchall())
        cur.execute("SELECT pref_name, id FROM pharmaco_db.compounds WHERE pref_name IS NOT NULL")
        name_map = {k.upper(): v for k, v in cur.fetchall()}

    batch = []
    seen = set()
    with open(f, "r") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene = row.get("gene_name", "")
            drug = row.get("drug_name", "")
            drug_chembl = row.get("drug_chembl_id", "")
            interaction_type = row.get("interaction_types", "") or row.get("interaction_group_score", "")
            score = row.get("interaction_group_score", "")
            source_col = row.get("interaction_claim_source", "")
            pmids_col = row.get("PMIDs", "")

            tid = gene_map.get(gene)
            cid = chembl_map.get(drug_chembl) or (name_map.get(drug.upper()) if drug else None)

            if not tid:
                continue

            # Dedup key
            key = (gene, drug, interaction_type)
            if key in seen:
                continue
            seen.add(key)

            # Parse sources and PMIDs as arrays
            sources_arr = [s.strip() for s in source_col.split(",")] if source_col else None
            pmids_arr = [p.strip() for p in pmids_col.split(",")] if pmids_col else None

            batch.append((
                cid, tid, gene, drug, drug_chembl,
                interaction_type[:200] if interaction_type else None,
                float(score) if score and score.replace(".", "").replace("-", "").isdigit() else None,
                sources_arr,
                pmids_arr,
            ))

    if batch:
        with conn.cursor() as cur:
            cur.execute("TRUNCATE pharmaco_db.dgidb_interactions")
            psycopg2.extras.execute_values(cur, """
                INSERT INTO pharmaco_db.dgidb_interactions
                (compound_id, target_id, gene_name, drug_name, drug_chembl_id,
                 interaction_type, interaction_score, sources, pmids)
                VALUES %s
            """, batch, page_size=5000)
        conn.commit()
    log(f"    DGIdb: {len(batch):,} interactions")

def ingest_sider(conn):
    """Download and ingest SIDER side effects.
    Schema: side_effects(compound_id, stitch_id, side_effect_name,
            meddra_concept_type, meddra_id, frequency, source)
    """
    log("  Ingesting SIDER side effects...")

    url = "http://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz"
    dest = DOWNLOAD_DIR / "sider_all_se.tsv.gz"
    f = download_file(url, dest, "SIDER side effects")
    if not f:
        log("    SKIP: SIDER download failed")
        return

    # SIDER uses STITCH compound IDs -> map via PubChem CID
    with conn.cursor() as cur:
        cur.execute("SELECT pubchem_cid, id FROM pharmaco_db.compounds WHERE pubchem_cid IS NOT NULL")
        pubchem_map = dict(cur.fetchall())
    log(f"    PubChem CID map: {len(pubchem_map):,} entries")

    batch = []
    seen = set()
    with gzip.open(f, "rt") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            stitch_flat, stitch_stereo, se_type, umls_cui, meddra_type, meddra_term = parts[:6]

            try:
                if stitch_flat.startswith("CID1"):
                    pcid = int(stitch_flat[4:])
                elif stitch_flat.startswith("CID"):
                    pcid = int(stitch_flat[3:])
                else:
                    continue
            except ValueError:
                continue

            cid = pubchem_map.get(pcid)
            if not cid:
                continue

            key = (stitch_flat, umls_cui, meddra_type)
            if key in seen:
                continue
            seen.add(key)

            batch.append((cid, stitch_flat, meddra_term, meddra_type, umls_cui, None, "sider"))

    if batch:
        with conn.cursor() as cur:
            cur.execute("TRUNCATE pharmaco_db.side_effects")
            psycopg2.extras.execute_values(cur, """
                INSERT INTO pharmaco_db.side_effects
                (compound_id, stitch_id, side_effect_name,
                 meddra_concept_type, meddra_id, frequency, source)
                VALUES %s
            """, batch, page_size=5000)
        conn.commit()
        log(f"    SIDER: {len(batch):,} unique side effects")
    else:
        log("    SIDER: 0 mapped (PubChem CID mapping needed)")

def ingest_reactome(conn):
    """Download and ingest Reactome pathways.
    Schema: pathways(target_id, uniprot_id, pathway_id, pathway_name,
            pathway_category, species, source)
    """
    log("  Ingesting Reactome pathways...")

    url = "https://reactome.org/download/current/UniProt2Reactome.txt"
    dest = DOWNLOAD_DIR / "uniprot2reactome.txt"
    f = download_file(url, dest, "Reactome UniProt mapping")
    if not f:
        log("    SKIP: Reactome download failed")
        return

    with conn.cursor() as cur:
        cur.execute("SELECT uniprot_id, id FROM pharmaco_db.targets WHERE uniprot_id IS NOT NULL")
        uniprot_map = dict(cur.fetchall())

    batch = []
    seen = set()
    with open(f, "r") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            uniprot_id, pathway_id, url_col, pathway_name, evidence, species = parts[:6]
            if species != "Homo sapiens":
                continue

            tid = uniprot_map.get(uniprot_id)
            if not tid:
                continue

            key = (tid, pathway_id)
            if key in seen:
                continue
            seen.add(key)

            category = "signal_transduction"
            if "R-HSA" in pathway_id:
                category = "reactome"

            batch.append((tid, uniprot_id, pathway_id, pathway_name, category, species, "reactome"))

    if batch:
        with conn.cursor() as cur:
            cur.execute("DELETE FROM pharmaco_db.pathways WHERE source = 'reactome'")
            psycopg2.extras.execute_values(cur, """
                INSERT INTO pharmaco_db.pathways
                (target_id, uniprot_id, pathway_id, pathway_name,
                 pathway_category, species, source)
                VALUES %s
            """, batch, page_size=5000)
        conn.commit()
    log(f"    Reactome: {len(batch):,} pathway-target mappings")

def ingest_disgenet(conn):
    """Download and ingest DisGeNET disease-gene associations.
    Schema: disgenet_associations(target_id, gene_name, gene_ncbi_id,
            disease_id, disease_name, disease_type, association_score,
            ei_score, num_pmids, num_snps, source)
    """
    log("  Ingesting DisGeNET associations...")

    url = "https://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_gene_disease_associations.tsv.gz"
    dest = DOWNLOAD_DIR / "disgenet_curated.tsv.gz"
    f = download_file(url, dest, "DisGeNET curated")
    if not f:
        url = "https://www.disgenet.org/static/disgenet_ap1/files/downloads/all_gene_disease_associations.tsv.gz"
        dest = DOWNLOAD_DIR / "disgenet_all.tsv.gz"
        f = download_file(url, dest, "DisGeNET all")
    if not f:
        log("    SKIP: DisGeNET download failed (may need API key)")
        return

    with conn.cursor() as cur:
        cur.execute("SELECT gene_name, id FROM pharmaco_db.targets WHERE gene_name IS NOT NULL")
        gene_map = dict(cur.fetchall())

    batch = []
    seen = set()
    try:
        with gzip.open(f, "rt") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                gene = row.get("geneSymbol", "")
                tid = gene_map.get(gene)
                if not tid:
                    continue

                disease_id = row.get("diseaseId", "")
                key = (gene, disease_id)
                if key in seen:
                    continue
                seen.add(key)

                disease_name = row.get("diseaseName", "")
                disease_type = row.get("diseaseType", "")
                score = row.get("score", "")
                ei = row.get("EI", "")
                n_pmids = row.get("NofPmids", "")
                n_snps = row.get("NofSnps", "")
                ncbi_id = row.get("geneId", "")

                batch.append((
                    tid, gene,
                    int(ncbi_id) if ncbi_id and ncbi_id.isdigit() else None,
                    disease_id, disease_name, disease_type,
                    float(score) if score else None,
                    float(ei) if ei else None,
                    int(n_pmids) if n_pmids and n_pmids.isdigit() else None,
                    int(n_snps) if n_snps and n_snps.isdigit() else None,
                    "disgenet"
                ))
    except Exception as e:
        log(f"    DisGeNET parse error: {e}")
        return

    if batch:
        with conn.cursor() as cur:
            cur.execute("TRUNCATE pharmaco_db.disgenet_associations")
            psycopg2.extras.execute_values(cur, """
                INSERT INTO pharmaco_db.disgenet_associations
                (target_id, gene_name, gene_ncbi_id, disease_id, disease_name,
                 disease_type, association_score, ei_score, num_pmids, num_snps, source)
                VALUES %s
            """, batch, page_size=5000)
        conn.commit()
    log(f"    DisGeNET: {len(batch):,} associations")

def run_phase3_public_sources():
    """Phase 3: Download and ingest free public data sources."""
    log("=" * 60)
    log("PHASE 3: Free Public Data Sources")
    log("=" * 60)

    DOWNLOAD_DIR.mkdir(parents=True, exist_ok=True)
    conn = get_conn()

    try:
        ingest_string_ppi(conn)
        ingest_reactome(conn)
        ingest_dgidb(conn)
        ingest_disgenet(conn)
        ingest_sider(conn)
    finally:
        conn.close()

    log("Phase 3 COMPLETE")


# ============================================================
# FINAL STATS
# ============================================================

def print_final_stats():
    log("=" * 60)
    log("FINAL DATABASE STATISTICS")
    log("=" * 60)

    conn = get_conn()
    with conn.cursor() as cur:
        tables = [
            "compounds", "targets", "assays", "bioactivities",
            "drug_mechanisms", "drug_indications", "cross_references",
            "compound_tags", "target_tags", "activity_tags",
            "atc_classification", "binding_sites", "target_interactions",
            "selectivity_profiles", "protein_interactions", "pathways",
            "dgidb_interactions", "side_effects", "disgenet_associations",
        ]
        for tbl in tables:
            try:
                cur.execute(f"SELECT COUNT(*) FROM pharmaco_db.{tbl}")
                cnt = cur.fetchone()[0]
                log(f"  {tbl:30s} {cnt:>12,}")
            except Exception:
                conn.rollback()
                log(f"  {tbl:30s} {'N/A':>12s}")

        cur.execute("SELECT pg_size_pretty(pg_database_size('pharmaco'))")
        size = cur.fetchone()[0]
        log(f"\n  Database size: {size}")
    conn.close()


# ============================================================
# PHASE 4: COLUMN DESCRIPTIONS
# ============================================================

def apply_column_comments():
    """Apply COMMENT ON COLUMN for ALL tables including new Phase 2/3 tables."""
    log("=" * 60)
    log("PHASE 4: Applying column descriptions")
    log("=" * 60)

    conn = get_conn()
    comments = [
        # ── COMPOUND TAGS ──
        ("compound_tags", "compound_id", "FK to compounds table"),
        ("compound_tags", "lipinski_pass", "Passes Lipinski Ro5: MW<=500, logP<=5, HBD<=5, HBA<=10"),
        ("compound_tags", "lipinski_violations", "Count of Ro5 violations (0-4)"),
        ("compound_tags", "veber_pass", "Passes Veber rules: TPSA<=140 and RB<=10 (oral bioavailability proxy)"),
        ("compound_tags", "lead_like", "Lead-like: MW 200-350, logP -1..3, HBD<=3, HBA<=6"),
        ("compound_tags", "fragment_like", "Fragment-like Ro3: MW<=300, logP<=3, HBD<=3, HBA<=3, RB<=3"),
        ("compound_tags", "drug_like", "QED >= 0.5 (quantitative drug-likeness)"),
        ("compound_tags", "ppi_inhibitor_like", "PPI inhibitor-like: MW>400, logP>4, aromatic_rings>=4"),
        ("compound_tags", "brenk_alerts", "Count of Brenk structural alerts (unwanted substructures)"),
        ("compound_tags", "pains_alerts", "Count of PAINS substructures (pan-assay interference)"),
        ("compound_tags", "has_reactive_group", "Contains PAINS or Brenk reactive group"),
        ("compound_tags", "cns_mpo_score", "CNS multiparameter optimization score (0-6). >=4 = CNS penetrant"),
        ("compound_tags", "cns_penetrant", "CNS MPO score >= 4.0"),
        ("compound_tags", "bbb_permeable", "Predicted blood-brain barrier permeability"),
        ("compound_tags", "mw_bin", "MW bin: <200, 200-350, 350-500, >500"),
        ("compound_tags", "logp_bin", "LogP bin: <0, 0-2, 2-4, >4"),
        ("compound_tags", "tpsa_bin", "TPSA bin: <60, 60-90, 90-140, >140"),
        ("compound_tags", "complexity_bin", "Complexity bin by heavy atoms: low(<15), medium(15-25), high(25-40), very_high(>40)"),
        ("compound_tags", "murcko_scaffold", "Generic Murcko scaffold SMILES"),
        ("compound_tags", "scaffold_cluster", "Butina clustering ID"),
        ("compound_tags", "chemical_series", "Named chemical series if identified"),
        ("compound_tags", "privileged_scaffold", "Contains known privileged pharmacophore motif"),
        ("compound_tags", "num_stereocenters", "Count of stereocenters"),
        ("compound_tags", "fraction_csp3", "Fraction of sp3 carbons (3D complexity)"),
        ("compound_tags", "clinical_phase", "Max clinical trial phase (0=preclinical, 4=approved)"),
        ("compound_tags", "is_approved", "Phase 4 approved drug"),
        ("compound_tags", "is_withdrawn", "Withdrawn from market"),
        ("compound_tags", "is_prodrug", "Known prodrug"),
        ("compound_tags", "atc_code_l1", "ATC level 1 — Anatomical main group (e.g. L=Antineoplastic)"),
        ("compound_tags", "atc_code_l2", "ATC level 2 — Therapeutic subgroup"),
        ("compound_tags", "atc_code_l3", "ATC level 3 — Pharmacological subgroup"),
        ("compound_tags", "atc_code_l4", "ATC level 4 — Chemical subgroup"),
        ("compound_tags", "atc_code_l5", "ATC level 5 — Chemical substance"),
        ("compound_tags", "therapeutic_areas", "Therapeutic areas array (oncology, neurology, etc.)"),
        ("compound_tags", "cyp_substrate", "CYP enzymes this compound is a substrate of"),
        ("compound_tags", "cyp_inhibitor", "CYP enzymes this compound inhibits"),
        ("compound_tags", "elimination_route", "Primary elimination route: hepatic, renal, mixed"),
        ("compound_tags", "is_herg_liability", "hERG IC50 < 10 uM (cardiac risk)"),
        ("compound_tags", "is_hepatotoxic", "Known hepatotoxicity liability"),
        ("compound_tags", "is_mutagenic", "Ames test positive (mutagenic)"),
        ("compound_tags", "is_carcinogenic", "Known carcinogenicity"),
        ("compound_tags", "tox_class", "GHS toxicity class (1-5)"),
        ("compound_tags", "num_bioactivities", "Total bioactivity measurements for this compound"),
        ("compound_tags", "num_targets", "Number of distinct targets tested"),
        ("compound_tags", "num_assays", "Number of distinct assays"),
        ("compound_tags", "data_completeness", "Data completeness score (0-1): SMILES+props+activity+targets+clinical"),

        # ── TARGET TAGS ──
        ("target_tags", "target_id", "FK to targets table"),
        ("target_tags", "target_class", "Target class: kinase, gpcr, ion_channel, nuclear_receptor, protease, etc."),
        ("target_tags", "target_subclass", "Target subclass (e.g. tyrosine_kinase, serine_threonine_kinase)"),
        ("target_tags", "enzyme_class_ec", "EC enzyme commission number (e.g. 2.7.10.1)"),
        ("target_tags", "druggability_score", "Composite druggability score 0-1 (structure+drugs+compounds+class+data)"),
        ("target_tags", "has_crystal_structure", "At least one PDB crystal structure available"),
        ("target_tags", "num_pdb_structures", "Number of PDB structures"),
        ("target_tags", "has_active_site_known", "Active site residues annotated"),
        ("target_tags", "has_allosteric_site", "Known allosteric binding site"),
        ("target_tags", "binding_site_druggability", "Binding site druggability: high, medium, low"),
        ("target_tags", "disease_areas", "Disease areas array (cancer, inflammation, etc.)"),
        ("target_tags", "disease_ids", "EFO/MONDO disease IDs"),
        ("target_tags", "genetic_association", "Has GWAS/genetic evidence for disease"),
        ("target_tags", "somatic_mutation", "Cancer driver via somatic mutation"),
        ("target_tags", "known_lof_phenotype", "Known loss-of-function phenotype"),
        ("target_tags", "tissue_expression", "Tissues with high expression"),
        ("target_tags", "tissue_specificity", "Tissue specificity: ubiquitous, tissue_specific, tissue_enriched"),
        ("target_tags", "is_essential", "Essential gene (DepMap)"),
        ("target_tags", "essentiality_score", "DepMap essentiality score"),
        ("target_tags", "safety_risk_level", "Safety risk: low, medium, high"),
        ("target_tags", "has_safety_concern", "Has identified safety concern"),
        ("target_tags", "safety_concerns", "Specific safety concerns array"),
        ("target_tags", "family_size", "Number of close homologs in same protein family"),
        ("target_tags", "selectivity_challenge", "Selectivity difficulty: easy(<10 homologs), moderate(10-50), hard(>50)"),
        ("target_tags", "num_active_compounds", "Compounds with pChEMBL>=7 (active)"),
        ("target_tags", "num_total_activities", "Total bioactivity measurements"),
        ("target_tags", "num_approved_drugs", "Count of approved drugs targeting this protein"),
        ("target_tags", "median_pchembl", "Median pChEMBL across all measurements"),
        ("target_tags", "data_maturity", "Data maturity: mature(>=100 actives), emerging(10-100), sparse(<10)"),

        # ── ACTIVITY TAGS ──
        ("activity_tags", "activity_id", "FK to bioactivities table"),
        ("activity_tags", "confidence_level", "Assay confidence: high(>=8), medium(6-7), low(<6)"),
        ("activity_tags", "is_direct_binding", "Direct binding assay (assay_type=B)"),
        ("activity_tags", "is_cell_based", "Cell-based functional assay (assay_type=F)"),
        ("activity_tags", "is_in_vivo", "In vivo data"),
        ("activity_tags", "pchembl_bin", "Potency bin: inactive(<5), weak(5-6), moderate(6-7), potent(7-8), very_potent(>8)"),
        ("activity_tags", "selectivity_ratio", "Activity delta vs closest off-target"),
        ("activity_tags", "is_activity_cliff", "Part of an activity cliff pair (high structural similarity, >2 log potency difference)"),
        ("activity_tags", "cliff_partner_id", "Partner compound in the activity cliff pair"),
        ("activity_tags", "cliff_fold_change", "Fold change in activity cliff pair"),

        # ── ATC CLASSIFICATION ──
        ("atc_classification", "compound_id", "FK to compounds table"),
        ("atc_classification", "atc_code", "Full ATC code (most specific level available)"),
        ("atc_classification", "level1", "ATC level 1 — Anatomical main group (A-V)"),
        ("atc_classification", "level1_desc", "Level 1 description (e.g. Alimentary tract, Nervous system)"),
        ("atc_classification", "level2", "ATC level 2 — Therapeutic subgroup"),
        ("atc_classification", "level2_desc", "Level 2 description"),
        ("atc_classification", "level3", "ATC level 3 — Pharmacological subgroup"),
        ("atc_classification", "level3_desc", "Level 3 description"),
        ("atc_classification", "level4", "ATC level 4 — Chemical subgroup"),
        ("atc_classification", "level4_desc", "Level 4 description"),
        ("atc_classification", "level5", "ATC level 5 — Chemical substance"),
        ("atc_classification", "level5_desc", "Level 5 description"),
        ("atc_classification", "who_name", "WHO International Nonproprietary Name"),

        # ── BINDING SITES ──
        ("binding_sites", "target_id", "FK to targets table"),
        ("binding_sites", "site_name", "Binding site name from ChEMBL"),
        ("binding_sites", "site_id", "ChEMBL binding site ID"),
        ("binding_sites", "domain_type", "Domain type (Pfam, InterPro, etc.)"),
        ("binding_sites", "domain_name", "Domain name"),
        ("binding_sites", "domain_description", "Domain description"),
        ("binding_sites", "start_position", "Start position in protein sequence"),
        ("binding_sites", "end_position", "End position in protein sequence"),
        ("binding_sites", "component_type", "Component type (protein, nucleic acid, etc.)"),

        # ── SELECTIVITY PROFILES ──
        ("selectivity_profiles", "compound_id", "FK to compounds table"),
        ("selectivity_profiles", "primary_target_id", "FK to most potent target"),
        ("selectivity_profiles", "off_target_id", "FK to off-target"),
        ("selectivity_profiles", "primary_pchembl", "pChEMBL on primary target"),
        ("selectivity_profiles", "off_target_pchembl", "pChEMBL on off-target"),
        ("selectivity_profiles", "selectivity_index", "Delta: primary_pchembl - off_target_pchembl. >2 = selective"),

        # ── PROTEIN INTERACTIONS (STRING) ──
        ("protein_interactions", "target_id_1", "FK to first target in interaction pair"),
        ("protein_interactions", "target_id_2", "FK to second target in interaction pair"),
        ("protein_interactions", "gene_1", "Gene symbol of first protein"),
        ("protein_interactions", "gene_2", "Gene symbol of second protein"),
        ("protein_interactions", "combined_score", "STRING combined confidence score (0-1000). >=700 = high confidence"),
        ("protein_interactions", "experimental_score", "Experimental evidence sub-score"),
        ("protein_interactions", "database_score", "Database (curated) evidence sub-score"),
        ("protein_interactions", "textmining_score", "Text mining evidence sub-score"),
        ("protein_interactions", "coexpression_score", "Co-expression evidence sub-score"),

        # ── PATHWAYS (Reactome) ──
        ("pathways", "target_id", "FK to targets table"),
        ("pathways", "uniprot_id", "UniProt accession for this target"),
        ("pathways", "pathway_id", "Reactome/KEGG pathway identifier (e.g. R-HSA-109582)"),
        ("pathways", "pathway_name", "Pathway name"),
        ("pathways", "pathway_category", "Pathway category"),
        ("pathways", "species", "Species (default: Homo sapiens)"),

        # ── DGIdb INTERACTIONS ──
        ("dgidb_interactions", "compound_id", "FK to compounds (nullable, mapped via ChEMBL ID)"),
        ("dgidb_interactions", "target_id", "FK to targets (mapped via gene symbol)"),
        ("dgidb_interactions", "gene_name", "Gene symbol from DGIdb"),
        ("dgidb_interactions", "drug_name", "Drug name from DGIdb"),
        ("dgidb_interactions", "drug_chembl_id", "ChEMBL ID if available"),
        ("dgidb_interactions", "interaction_type", "Interaction type: inhibitor, agonist, antagonist, etc."),
        ("dgidb_interactions", "interaction_score", "DGIdb interaction group score"),
        ("dgidb_interactions", "sources", "Source databases array (DrugBank, PharmGKB, etc.)"),
        ("dgidb_interactions", "pmids", "PubMed IDs supporting this interaction"),

        # ── SIDER SIDE EFFECTS ──
        ("side_effects", "compound_id", "FK to compounds (mapped via PubChem CID)"),
        ("side_effects", "stitch_id", "STITCH compound ID (links to PubChem CID)"),
        ("side_effects", "side_effect_name", "MedDRA preferred term for the side effect"),
        ("side_effects", "meddra_concept_type", "MedDRA concept type: PT (preferred term), LLT (low-level term)"),
        ("side_effects", "meddra_id", "MedDRA concept UMLS CUI"),
        ("side_effects", "frequency", "Reported frequency: rare, infrequent, frequent, postmarketing"),

        # ── DisGeNET ──
        ("disgenet_associations", "target_id", "FK to targets table"),
        ("disgenet_associations", "gene_name", "Gene symbol"),
        ("disgenet_associations", "gene_ncbi_id", "NCBI Gene/Entrez ID"),
        ("disgenet_associations", "disease_id", "UMLS CUI disease identifier (e.g. C0006142)"),
        ("disgenet_associations", "disease_name", "Disease name"),
        ("disgenet_associations", "disease_type", "Disease type: disease, phenotype, group"),
        ("disgenet_associations", "association_score", "GDA association score (0-1). Higher = stronger evidence"),
        ("disgenet_associations", "ei_score", "Evidence Index: ratio of supporting vs total publications"),
        ("disgenet_associations", "num_pmids", "Number of PubMed articles supporting association"),
        ("disgenet_associations", "num_snps", "Number of SNPs associated"),

        # ── DISEASE-TARGET ASSOCIATIONS (Open Targets) ──
        ("disease_target_associations", "target_id", "FK to targets table"),
        ("disease_target_associations", "disease_id", "EFO/MONDO disease identifier"),
        ("disease_target_associations", "disease_name", "Disease name"),
        ("disease_target_associations", "therapeutic_area", "Therapeutic area classification"),
        ("disease_target_associations", "overall_score", "Open Targets overall association score (0-1)"),
        ("disease_target_associations", "genetic_score", "Genetic evidence score"),
        ("disease_target_associations", "somatic_score", "Somatic mutation evidence score"),
        ("disease_target_associations", "known_drug_score", "Known drug evidence score"),
        ("disease_target_associations", "literature_score", "Literature mining evidence score"),
        ("disease_target_associations", "rna_expression_score", "RNA expression evidence score"),
        ("disease_target_associations", "animal_model_score", "Animal model evidence score"),

        # ── TOXICOLOGY ──
        ("toxicology_data", "compound_id", "FK to compounds table"),
        ("toxicology_data", "assay_name", "Tox21/ToxCast assay name"),
        ("toxicology_data", "assay_endpoint", "Assay endpoint (e.g. NR-AR, SR-MMP, NR-ER)"),
        ("toxicology_data", "activity_outcome", "Activity outcome: active or inactive"),
        ("toxicology_data", "ac50_um", "AC50 concentration in micromolar"),
        ("toxicology_data", "efficacy", "Maximum response efficacy"),

        # ── PHARMACOKINETICS ──
        ("pharmacokinetics", "compound_id", "FK to compounds table"),
        ("pharmacokinetics", "species", "Species: human, rat, mouse, dog"),
        ("pharmacokinetics", "route", "Route of administration: oral, iv, ip, sc"),
        ("pharmacokinetics", "parameter", "PK parameter: Cmax, AUC, t1/2, F%, Vd, CL"),
        ("pharmacokinetics", "value", "Parameter value"),
        ("pharmacokinetics", "units", "Parameter units"),

        # ── GtoP INTERACTIONS ──
        ("gtop_interactions", "compound_id", "FK to compounds table"),
        ("gtop_interactions", "target_id", "FK to targets table"),
        ("gtop_interactions", "ligand_id", "Guide to Pharmacology ligand ID"),
        ("gtop_interactions", "target_gtop_id", "Guide to Pharmacology target ID"),
        ("gtop_interactions", "interaction_type", "Interaction type: agonist, antagonist, inhibitor, etc."),
        ("gtop_interactions", "affinity_type", "Affinity measurement type: pKi, pIC50, pEC50, pKd"),
        ("gtop_interactions", "affinity_value", "Affinity value (-log10 molar)"),
        ("gtop_interactions", "endogenous", "Is endogenous ligand"),
        ("gtop_interactions", "primary_target", "Is primary pharmacological target"),

        # ── KEGG ──
        ("kegg_drug_info", "compound_id", "FK to compounds table"),
        ("kegg_drug_info", "kegg_drug_id", "KEGG drug identifier (e.g. D00123)"),
        ("kegg_drug_info", "kegg_name", "Drug name in KEGG"),
        ("kegg_drug_info", "kegg_formula", "Molecular formula from KEGG"),
        ("kegg_drug_info", "therapeutic_category", "KEGG therapeutic categories array"),
        ("kegg_drug_info", "pathway_ids", "KEGG pathway IDs array"),
        ("kegg_drug_info", "pathway_names", "KEGG pathway names array"),

        # ── CLINICAL TRIALS ──
        ("clinical_trials", "compound_id", "FK to compounds table"),
        ("clinical_trials", "nct_id", "ClinicalTrials.gov NCT identifier"),
        ("clinical_trials", "title", "Trial title"),
        ("clinical_trials", "status", "Trial status: Recruiting, Completed, Terminated, etc."),
        ("clinical_trials", "phase", "Trial phase: Phase 1, Phase 2, Phase 3, Phase 4"),
        ("clinical_trials", "conditions", "Disease conditions studied (array)"),
        ("clinical_trials", "interventions", "Interventions/treatments (array)"),
        ("clinical_trials", "enrollment", "Number of participants enrolled"),

        # ── ADVERSE EVENTS (OpenFDA) ──
        ("adverse_events", "compound_id", "FK to compounds table"),
        ("adverse_events", "drug_name", "Drug name from OpenFDA"),
        ("adverse_events", "reaction", "Adverse reaction (MedDRA preferred term)"),
        ("adverse_events", "count", "Number of reports for this drug-reaction pair"),

        # ── TARGET INTERACTIONS (ChEMBL) ──
        ("target_interactions", "target_id_1", "FK to first target"),
        ("target_interactions", "target_id_2", "FK to second target"),
        ("target_interactions", "interaction_type", "Relationship type: SUBSET OF, SUPERSET OF, EQUIVALENT TO, OVERLAPS WITH"),
        ("target_interactions", "source", "Data source (chembl, string, etc.)"),
        ("target_interactions", "confidence", "Confidence score (0-1)"),
    ]

    with conn.cursor() as cur:
        applied = 0
        for table, column, description in comments:
            try:
                cur.execute(f"""
                    COMMENT ON COLUMN pharmaco_db.{table}.{column}
                    IS %s
                """, (description,))
                applied += 1
            except Exception as e:
                conn.rollback()
                # Column may not exist yet
                continue
        conn.commit()

    log(f"  Applied {applied}/{len(comments)} column descriptions")
    conn.close()
    log("Phase 4 COMPLETE")


# ============================================================
# MAIN
# ============================================================

def main():
    start = time.time()
    log("=" * 60)
    log("MegaPharmaDB — Comprehensive Consolidation")
    log("=" * 60)

    phase = sys.argv[1] if len(sys.argv) > 1 else "all"

    if phase in ("all", "1"):
        run_phase1_chembl_mining()

    if phase in ("all", "2"):
        run_phase2_compute_tags()

    if phase in ("all", "3"):
        run_phase3_public_sources()

    if phase in ("all", "4"):
        apply_column_comments()

    print_final_stats()

    elapsed = time.time() - start
    log(f"\nTotal time: {elapsed / 60:.1f} minutes")

if __name__ == "__main__":
    main()
