#!/usr/bin/env python3
"""
PharmacoDB — Fast migration from ChEMBL 36 dump.
Both databases are in the same PostgreSQL container.
Uses server-side cursors for memory efficiency on 24M+ rows.
"""

import psycopg2
import psycopg2.extras
import time
from datetime import datetime

PHARMACO_DB = {
    "host": "localhost", "port": 5433,
    "dbname": "pharmaco", "user": "postgres", "password": "pharmaco_secret"
}
CHEMBL_DB = {
    "host": "localhost", "port": 5433,
    "dbname": "chembl_36", "user": "postgres", "password": "pharmaco_secret"
}

BATCH = 10000

def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

def migrate_table(src_conn, dst_conn, query, insert_sql, label, conflict=""):
    """Generic migration: read from ChEMBL, write to pharmaco_db in batches."""
    t0 = time.time()
    total = 0

    # Server-side cursor for memory efficiency
    if src_conn.autocommit:
        src_conn.autocommit = False
    # Ensure clean transaction state
    src_conn.rollback()
    with src_conn.cursor(name='migrate_cursor') as src_cur:
        src_cur.itersize = BATCH
        src_cur.execute(query)

        batch = []
        for row in src_cur:
            batch.append(row)
            if len(batch) >= BATCH:
                with dst_conn.cursor() as dst_cur:
                    psycopg2.extras.execute_values(
                        dst_cur,
                        insert_sql + conflict,
                        batch,
                        page_size=5000
                    )
                dst_conn.commit()
                total += len(batch)
                batch = []
                elapsed = time.time() - t0
                rate = total / elapsed if elapsed > 0 else 0
                log(f"  {label}: {total:,} rows ({rate:,.0f}/s)")

        # Flush remainder
        if batch:
            with dst_conn.cursor() as dst_cur:
                psycopg2.extras.execute_values(
                    dst_cur,
                    insert_sql + conflict,
                    batch,
                    page_size=5000
                )
            dst_conn.commit()
            total += len(batch)

    elapsed = time.time() - t0
    rate = total / elapsed if elapsed > 0 else 0
    log(f"  {label} DONE: {total:,} rows in {int(elapsed)}s ({rate:,.0f}/s)")
    return total

def main():
    log("=" * 60)
    log("PharmacoDB — MIGRATE FROM CHEMBL 36 DUMP")
    log("=" * 60)

    src = psycopg2.connect(**CHEMBL_DB)
    dst = psycopg2.connect(**PHARMACO_DB)
    t_start = time.time()

    # ─── 0. Update existing compounds (max_phase, is_drug) ─────────
    log("\n>>> Updating existing compounds max_phase & is_drug...")
    with src.cursor() as cur:
        cur.execute("""
            SELECT chembl_id, max_phase
            FROM molecule_dictionary
            WHERE max_phase > 0
        """)
        updates = cur.fetchall()
    log(f"  Found {len(updates):,} compounds with max_phase > 0")

    batch = []
    updated = 0
    for chembl_id, max_phase in updates:
        batch.append((max_phase, max_phase >= 4, chembl_id))
        if len(batch) >= BATCH:
            with dst.cursor() as cur:
                psycopg2.extras.execute_batch(cur, """
                    UPDATE pharmaco_db.compounds
                    SET max_phase = %s, is_drug = %s
                    WHERE chembl_id = %s
                """, batch, page_size=2000)
            dst.commit()
            updated += len(batch)
            batch = []
    if batch:
        with dst.cursor() as cur:
            psycopg2.extras.execute_batch(cur, """
                UPDATE pharmaco_db.compounds
                SET max_phase = %s, is_drug = %s
                WHERE chembl_id = %s
            """, batch, page_size=2000)
        dst.commit()
        updated += len(batch)
    log(f"  Updated {updated:,} compounds")

    # ─── 1. Insert missing compounds ──────────────────────────────
    log("\n>>> Inserting missing compounds from dump...")

    with dst.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.compounds")
        before = cur.fetchone()[0]
    log(f"  Before: {before:,} compounds")

    migrate_table(src, dst,
        query="""
            SELECT md.chembl_id, cs.canonical_smiles, cs.standard_inchi, cs.standard_inchi_key,
                   md.pref_name,
                   cp.mw_freebase::double precision, cp.alogp::double precision,
                   cp.hba::int, cp.hbd::int, cp.psa::double precision, cp.rtb::int,
                   cp.num_ro5_violations::int, cp.aromatic_rings::int, cp.heavy_atoms::int,
                   cp.qed_weighted::double precision, cp.full_molformula,
                   md.max_phase,
                   CASE WHEN md.natural_product = 1 THEN true WHEN md.natural_product = 0 THEN false ELSE NULL END
            FROM molecule_dictionary md
            JOIN compound_structures cs ON cs.molregno = md.molregno
            JOIN compound_properties cp ON cp.molregno = md.molregno
            WHERE cs.canonical_smiles IS NOT NULL
              AND cp.mw_freebase >= 100 AND cp.mw_freebase <= 900
        """,
        insert_sql="""
            INSERT INTO pharmaco_db.compounds (
                chembl_id, canonical_smiles, inchi, inchi_key, pref_name,
                molecular_weight, alogp, hba, hbd, psa, rtb,
                num_ro5_violations, aromatic_rings, heavy_atoms, qed_weighted,
                molecular_formula, max_phase, is_natural_product
            ) VALUES %s
        """,
        label="compounds",
        conflict=" ON CONFLICT (chembl_id) DO NOTHING"
    )

    with dst.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.compounds")
        after = cur.fetchone()[0]
    log(f"  After: {after:,} compounds (+{after - before:,} new)")

    # ─── 2. Build lookup maps ─────────────────────────────────────
    log("\n>>> Building lookup maps...")

    compound_map = {}
    with dst.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.compounds")
        compound_map = dict(cur.fetchall())
    log(f"  Compound map: {len(compound_map):,}")

    target_map = {}
    with dst.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.targets")
        target_map = dict(cur.fetchall())
    log(f"  Target map: {len(target_map):,}")

    # ─── 3. Assays ────────────────────────────────────────────────
    log("\n>>> Inserting assays...")

    with dst.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.assays")
        existing = cur.fetchone()[0]
    # If we have partial assays from a previous crash, truncate and redo
    if existing > 0 and existing < 500000:
        log(f"  Partial assays ({existing:,}), truncating and redoing...")
        with dst.cursor() as cur:
            cur.execute("TRUNCATE pharmaco_db.assays CASCADE")
        dst.commit()
        existing = 0
    if existing >= 500000:
        log(f"  Already have {existing:,} assays, skipping.")
    else:
        # Read all assays, resolve target_id
        t0 = time.time()
        total = 0
        src.rollback()
        with src.cursor(name='assay_cursor') as cur:
            cur.itersize = BATCH
            cur.execute("""
                SELECT a.chembl_id, LEFT(a.description, 500), a.assay_type, a.assay_category,
                       a.src_id, td.chembl_id as target_chembl_id,
                       a.assay_organism, a.assay_cell_type, a.confidence_score
                FROM assays a
                LEFT JOIN target_dictionary td ON td.tid = a.tid
                WHERE a.assay_organism = 'Homo sapiens'
                  AND a.confidence_score >= 4
            """)
            batch = []
            for row in cur:
                target_id = target_map.get(row[5])
                # Clean description: replace any non-UTF8 safe chars
                desc = row[1]
                if desc:
                    desc = desc.encode('utf-8', errors='replace').decode('utf-8')
                batch.append((
                    row[0], desc, row[2], row[3],      # chembl_id, desc, type, category
                    row[4], None, None, None, None,     # src_id, src_desc, journal, year, doi
                    target_id, row[5],                   # target_id, target_chembl_id
                    row[6], row[7], row[8],              # organism, cell_type, confidence
                ))
                if len(batch) >= BATCH:
                    with dst.cursor() as dc:
                        psycopg2.extras.execute_values(dc, """
                            INSERT INTO pharmaco_db.assays (
                                chembl_id, description, assay_type, assay_category,
                                src_id, src_description, journal, year, doi,
                                target_id, target_chembl_id,
                                assay_organism, assay_cell_type, confidence_score
                            ) VALUES %s ON CONFLICT (chembl_id) DO NOTHING
                        """, batch, page_size=5000)
                    dst.commit()
                    total += len(batch)
                    batch = []
                    elapsed = time.time() - t0
                    log(f"  Assays: {total:,} ({total/elapsed:,.0f}/s)")

            if batch:
                with dst.cursor() as dc:
                    psycopg2.extras.execute_values(dc, """
                        INSERT INTO pharmaco_db.assays (
                            chembl_id, description, assay_type, assay_category,
                            src_id, src_description, journal, year, doi,
                            target_id, target_chembl_id,
                            assay_organism, assay_cell_type, confidence_score
                        ) VALUES %s ON CONFLICT (chembl_id) DO NOTHING
                    """, batch, page_size=5000)
                dst.commit()
                total += len(batch)

        elapsed = time.time() - t0
        log(f"  Assays DONE: {total:,} in {int(elapsed)}s")

    # Build assay map after insert
    assay_map = {}
    with dst.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.assays")
        assay_map = dict(cur.fetchall())
    log(f"  Assay map: {len(assay_map):,}")

    # ─── 4. BIOACTIVITIES (the big one) ──────────────────────────
    log("\n>>> Inserting bioactivities (24M source rows, filtered by pChEMBL)...")

    with dst.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.bioactivities")
        existing = cur.fetchone()[0]

    if existing > 0:
        log(f"  Already have {existing:,} bioactivities, skipping.")
    else:
        t0 = time.time()
        total = 0
        skipped = 0

        src.rollback()
        with src.cursor(name='bio_cursor') as cur:
            cur.itersize = BATCH
            cur.execute("""
                SELECT act.activity_id, act.standard_type, act.standard_relation,
                       act.standard_value, act.standard_units, act.pchembl_value,
                       act.data_validity_comment, act.potential_duplicate,
                       md.chembl_id as mol_chembl_id,
                       a.chembl_id as assay_chembl_id,
                       td.chembl_id as target_chembl_id
                FROM activities act
                JOIN molecule_dictionary md ON md.molregno = act.molregno
                JOIN assays a ON a.assay_id = act.assay_id
                LEFT JOIN target_dictionary td ON td.tid = a.tid
                WHERE act.pchembl_value IS NOT NULL
            """)

            batch = []
            for row in cur:
                compound_id = compound_map.get(row[8])
                if not compound_id:
                    skipped += 1
                    continue

                target_id = target_map.get(row[10])
                assay_id = assay_map.get(row[9])
                pchembl = row[5]

                activity_class = None
                if pchembl is not None:
                    if pchembl >= 7:
                        activity_class = 'active'
                    elif pchembl >= 5:
                        activity_class = 'intermediate'
                    else:
                        activity_class = 'inactive'

                batch.append((
                    compound_id, target_id, assay_id,
                    row[1], row[2], row[3], row[4], pchembl,
                    activity_class, 'ChEMBL', str(row[0]),
                    row[6],
                    True if row[7] == 1 else (False if row[7] == 0 else None),
                ))

                if len(batch) >= BATCH:
                    with dst.cursor() as dc:
                        psycopg2.extras.execute_values(dc, """
                            INSERT INTO pharmaco_db.bioactivities (
                                compound_id, target_id, assay_id,
                                activity_type, relation, value, units, pchembl_value,
                                activity_class, source, source_id,
                                data_validity, potential_duplicate
                            ) VALUES %s
                        """, batch, page_size=5000)
                    dst.commit()
                    total += len(batch)
                    batch = []
                    elapsed = time.time() - t0
                    rate = total / elapsed if elapsed > 0 else 0
                    log(f"  Bioactivities: {total:,} inserted, {skipped:,} skipped ({rate:,.0f}/s)")

            if batch:
                with dst.cursor() as dc:
                    psycopg2.extras.execute_values(dc, """
                        INSERT INTO pharmaco_db.bioactivities (
                            compound_id, target_id, assay_id,
                            activity_type, relation, value, units, pchembl_value,
                            activity_class, source, source_id,
                            data_validity, potential_duplicate
                        ) VALUES %s
                    """, batch, page_size=5000)
                dst.commit()
                total += len(batch)

        elapsed = time.time() - t0
        log(f"  Bioactivities DONE: {total:,} in {int(elapsed)}s ({total/elapsed:,.0f}/s)")

    # ─── 5. DRUG MECHANISMS ──────────────────────────────────────
    log("\n>>> Inserting drug mechanisms...")

    with dst.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.drug_mechanisms")
        if cur.fetchone()[0] > 0:
            log("  Already populated, skipping.")
        else:
            with src.cursor() as cur:
                cur.execute("""
                    SELECT md.chembl_id, td.chembl_id,
                           dm.mechanism_of_action, dm.action_type,
                           dm.direct_interaction, dm.molecular_mechanism,
                           dm.selectivity_comment
                    FROM drug_mechanism dm
                    JOIN molecule_dictionary md ON md.molregno = dm.molregno
                    LEFT JOIN target_dictionary td ON td.tid = dm.tid
                """)
                rows = []
                for r in cur:
                    cid = compound_map.get(r[0])
                    if not cid:
                        continue  # Skip mechanisms without matching compound
                    rows.append((
                        cid,
                        target_map.get(r[1]),
                        r[2], r[3],
                        True if r[4] == 1 else (False if r[4] == 0 else None),  # direct_interaction
                        str(r[5]) if r[5] is not None else None,  # molecular_mechanism (may be int)
                        r[6]
                    ))

            with dst.cursor() as dc:
                psycopg2.extras.execute_values(dc, """
                    INSERT INTO pharmaco_db.drug_mechanisms (
                        compound_id, target_id, mechanism_of_action,
                        action_type, direct_interaction, molecular_mechanism,
                        selectivity_comment
                    ) VALUES %s
                """, rows, page_size=2000)
            dst.commit()
            log(f"  Drug mechanisms DONE: {len(rows):,}")

    # ─── 6. DRUG INDICATIONS ─────────────────────────────────────
    log("\n>>> Inserting drug indications...")

    with dst.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.drug_indications")
        if cur.fetchone()[0] > 0:
            log("  Already populated, skipping.")
        else:
            with src.cursor() as cur:
                cur.execute("""
                    SELECT md.chembl_id,
                           di.mesh_id, di.mesh_heading,
                           di.efo_id, di.efo_term,
                           di.max_phase_for_ind
                    FROM drug_indication di
                    JOIN molecule_dictionary md ON md.molregno = di.molregno
                """)
                rows = []
                for r in cur:
                    cid = compound_map.get(r[0])
                    if cid:
                        rows.append((cid, r[1], r[2], r[3], r[4], r[5]))

            with dst.cursor() as dc:
                psycopg2.extras.execute_values(dc, """
                    INSERT INTO pharmaco_db.drug_indications (
                        compound_id, mesh_id, mesh_heading,
                        efo_id, efo_term, max_phase
                    ) VALUES %s
                """, rows, page_size=2000)
            dst.commit()
            log(f"  Drug indications DONE: {len(rows):,}")

    # ─── 7. Update target druggability ───────────────────────────
    log("\n>>> Updating target druggability...")
    with dst.cursor() as cur:
        cur.execute("""
            UPDATE pharmaco_db.targets t SET
                num_approved_drugs = sub.n,
                is_druggable = true
            FROM (
                SELECT target_id, COUNT(DISTINCT compound_id) AS n
                FROM pharmaco_db.drug_mechanisms
                WHERE target_id IS NOT NULL
                GROUP BY target_id
            ) sub
            WHERE t.id = sub.target_id
        """)
        log(f"  Updated {cur.rowcount:,} targets as druggable")
    dst.commit()

    # ─── Log ─────────────────────────────────────────────────────
    with dst.cursor() as cur:
        cur.execute("""
            INSERT INTO pharmaco_db.ingestion_log (source, step, status, rows_inserted, started_at, completed_at)
            VALUES ('chembl_36_dump', 'full_migration', 'completed', 0, %s, NOW())
        """, (datetime.utcnow(),))
    dst.commit()

    # ─── FINAL STATS ─────────────────────────────────────────────
    log("\n" + "=" * 60)
    log("MIGRATION COMPLETE — Final counts:")
    log("=" * 60)

    with dst.cursor() as cur:
        for t in ['compounds', 'targets', 'assays', 'bioactivities',
                   'drug_mechanisms', 'drug_indications', 'cross_references']:
            cur.execute(f"SELECT COUNT(*) FROM pharmaco_db.{t}")
            log(f"  {t:>20}: {cur.fetchone()[0]:>12,}")
        cur.execute("SELECT pg_size_pretty(pg_database_size('pharmaco'))")
        log(f"  {'DB size':>20}: {cur.fetchone()[0]}")

    total_elapsed = time.time() - t_start
    log(f"\nTOTAL TIME: {int(total_elapsed)}s ({total_elapsed/60:.1f}min)")
    log("=" * 60)

    src.close()
    dst.close()

if __name__ == "__main__":
    main()
