#!/usr/bin/env python3
"""
PharmacoDB — Comprehensive Tag Computation Engine
Computes ALL ML-relevant tags for compounds, targets, and activities.
"""

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

def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

def get_conn():
    return psycopg2.connect(**DB_CONFIG)

# ============================================================
# COMPOUND TAGS
# ============================================================
def compute_compound_tags(conn):
    log("=" * 60)
    log("Computing compound tags")
    log("=" * 60)

    with conn.cursor() as cur:
        # Step 1: Insert base records for all compounds
        log("  Step 1: Creating compound_tags records...")
        cur.execute("""
            INSERT INTO pharmaco_db.compound_tags (compound_id)
            SELECT id FROM pharmaco_db.compounds
            ON CONFLICT (compound_id) DO NOTHING
        """)
        n = cur.rowcount
        conn.commit()
        log(f"    Created {n} new tag records")

        # Step 2: Drug-likeness rules
        log("  Step 2: Computing drug-likeness rules...")
        cur.execute("""
            UPDATE pharmaco_db.compound_tags ct SET
                -- Lipinski Rule of 5
                lipinski_pass = (c.molecular_weight <= 500 AND c.alogp <= 5
                    AND c.hbd <= 5 AND c.hba <= 10),
                lipinski_violations = c.num_ro5_violations,
                -- Veber rules (oral bioavailability)
                veber_pass = (c.psa <= 140 AND c.rtb <= 10),
                -- Lead-likeness
                lead_like = (c.molecular_weight BETWEEN 200 AND 350
                    AND c.alogp BETWEEN -1 AND 3
                    AND COALESCE(c.hbd, 0) <= 3 AND COALESCE(c.hba, 0) <= 6),
                -- Fragment-likeness (Ro3)
                fragment_like = (c.molecular_weight <= 300 AND c.alogp <= 3
                    AND COALESCE(c.hbd, 0) <= 3 AND COALESCE(c.hba, 0) <= 3
                    AND COALESCE(c.rtb, 0) <= 3),
                -- Drug-likeness (QED)
                drug_like = (COALESCE(c.qed_weighted, 0) >= 0.5),
                -- PPI inhibitor-like
                ppi_inhibitor_like = (c.molecular_weight > 400 AND c.alogp > 4
                    AND COALESCE(c.aromatic_rings, 0) >= 4),
                -- Clinical status
                clinical_phase = c.max_phase,
                is_approved = (c.max_phase >= 4),
                updated_at = NOW()
            FROM pharmaco_db.compounds c
            WHERE c.id = ct.compound_id
        """)
        conn.commit()
        log(f"    Updated {cur.rowcount} rows")

        # Step 3: Physicochemical bins
        log("  Step 3: Computing physicochemical bins...")
        cur.execute("""
            UPDATE pharmaco_db.compound_tags ct SET
                mw_bin = CASE
                    WHEN c.molecular_weight < 200 THEN '<200'
                    WHEN c.molecular_weight BETWEEN 200 AND 350 THEN '200-350'
                    WHEN c.molecular_weight BETWEEN 350 AND 500 THEN '350-500'
                    ELSE '>500'
                END,
                logp_bin = CASE
                    WHEN c.alogp < 0 THEN '<0'
                    WHEN c.alogp BETWEEN 0 AND 2 THEN '0-2'
                    WHEN c.alogp BETWEEN 2 AND 4 THEN '2-4'
                    ELSE '>4'
                END,
                tpsa_bin = CASE
                    WHEN c.psa < 60 THEN '<60'
                    WHEN c.psa BETWEEN 60 AND 90 THEN '60-90'
                    WHEN c.psa BETWEEN 90 AND 140 THEN '90-140'
                    ELSE '>140'
                END,
                complexity_bin = CASE
                    WHEN c.heavy_atoms < 15 THEN 'low'
                    WHEN c.heavy_atoms BETWEEN 15 AND 25 THEN 'medium'
                    WHEN c.heavy_atoms BETWEEN 25 AND 40 THEN 'high'
                    ELSE 'very_high'
                END
            FROM pharmaco_db.compounds c
            WHERE c.id = ct.compound_id
        """)
        conn.commit()
        log(f"    Updated {cur.rowcount} rows")

        # Step 4: CNS MPO score (approximate)
        log("  Step 4: Computing CNS MPO scores...")
        cur.execute("""
            UPDATE pharmaco_db.compound_tags ct SET
                cns_mpo_score = (
                    -- MW component (0-1): ideal 250-360
                    CASE WHEN c.molecular_weight <= 360 THEN 1.0
                         WHEN c.molecular_weight <= 500 THEN 1.0 - (c.molecular_weight - 360) / 140
                         ELSE 0 END
                    +
                    -- logP component (0-1): ideal 1-3
                    CASE WHEN c.alogp BETWEEN 1 AND 3 THEN 1.0
                         WHEN c.alogp < 1 THEN GREATEST(0, 1.0 - (1 - c.alogp))
                         WHEN c.alogp > 3 THEN GREATEST(0, 1.0 - (c.alogp - 3) / 2)
                         ELSE 0 END
                    +
                    -- TPSA component (0-1): ideal 40-90
                    CASE WHEN c.psa BETWEEN 40 AND 90 THEN 1.0
                         WHEN c.psa < 40 THEN GREATEST(0, c.psa / 40)
                         WHEN c.psa > 90 THEN GREATEST(0, 1.0 - (c.psa - 90) / 50)
                         ELSE 0 END
                    +
                    -- HBD component (0-1): ideal 0-1
                    CASE WHEN COALESCE(c.hbd, 0) <= 1 THEN 1.0
                         WHEN c.hbd = 2 THEN 0.75
                         WHEN c.hbd = 3 THEN 0.25
                         ELSE 0 END
                    +
                    -- pKa proxy via HBA (0-1)
                    CASE WHEN COALESCE(c.hba, 0) <= 5 THEN 1.0
                         ELSE GREATEST(0, 1.0 - (c.hba - 5) / 5.0) END
                    +
                    -- Aromatic rings (0-1): ideal <=2
                    CASE WHEN COALESCE(c.aromatic_rings, 0) <= 2 THEN 1.0
                         WHEN c.aromatic_rings = 3 THEN 0.5
                         ELSE 0 END
                ),
                updated_at = NOW()
            FROM pharmaco_db.compounds c
            WHERE c.id = ct.compound_id
            AND c.molecular_weight IS NOT NULL
            AND c.alogp IS NOT NULL
        """)
        conn.commit()
        log(f"    Updated {cur.rowcount} rows")

        # Mark CNS penetrant
        cur.execute("""
            UPDATE pharmaco_db.compound_tags SET
                cns_penetrant = (cns_mpo_score >= 4.0)
            WHERE cns_mpo_score IS NOT NULL
        """)
        conn.commit()

        # Step 5: Data richness stats
        log("  Step 5: Computing data richness...")
        cur.execute("""
            UPDATE pharmaco_db.compound_tags ct SET
                num_bioactivities = sub.n_act,
                num_targets = sub.n_tgt,
                num_assays = sub.n_asy
            FROM (
                SELECT
                    compound_id,
                    COUNT(*) as n_act,
                    COUNT(DISTINCT target_id) as n_tgt,
                    COUNT(DISTINCT assay_id) as n_asy
                FROM pharmaco_db.bioactivities
                GROUP BY compound_id
            ) sub
            WHERE sub.compound_id = ct.compound_id
        """)
        conn.commit()
        log(f"    Updated {cur.rowcount} rows")

        # Step 6: Data completeness score
        log("  Step 6: Computing data completeness...")
        cur.execute("""
            UPDATE pharmaco_db.compound_tags ct SET
                data_completeness = (
                    -- Has SMILES: 0.1
                    0.1
                    -- Has properties: 0.2
                    + CASE WHEN c.molecular_weight IS NOT NULL THEN 0.2 ELSE 0 END
                    -- Has bioactivity: 0.3
                    + CASE WHEN ct.num_bioactivities > 0 THEN 0.3 ELSE 0 END
                    -- Has target data: 0.2
                    + CASE WHEN ct.num_targets > 0 THEN 0.2 ELSE 0 END
                    -- Has clinical data: 0.2
                    + CASE WHEN c.max_phase > 0 THEN 0.2 ELSE 0 END
                )
            FROM pharmaco_db.compounds c
            WHERE c.id = ct.compound_id
        """)
        conn.commit()
        log(f"    Updated {cur.rowcount} rows")

    log("  Compound tags DONE")

# ============================================================
# TARGET TAGS
# ============================================================
def compute_target_tags(conn):
    log("=" * 60)
    log("Computing target tags")
    log("=" * 60)

    with conn.cursor() as cur:
        # Step 1: Create records
        log("  Step 1: Creating target_tags records...")
        cur.execute("""
            INSERT INTO pharmaco_db.target_tags (target_id)
            SELECT id FROM pharmaco_db.targets
            ON CONFLICT (target_id) DO NOTHING
        """)
        conn.commit()
        log(f"    Created {cur.rowcount} new tag records")

        # Step 2: Target classification
        log("  Step 2: Target classification...")
        cur.execute("""
            UPDATE pharmaco_db.target_tags tt SET
                target_class = CASE
                    WHEN t.protein_class_l1 = 'Enzyme' AND t.protein_class_l2 LIKE '%kinase%' THEN 'kinase'
                    WHEN t.protein_class_l1 = 'Enzyme' AND t.protein_class_l2 LIKE '%protease%' THEN 'protease'
                    WHEN t.protein_class_l1 = 'Enzyme' AND t.protein_class_l2 LIKE '%phosphodiesterase%' THEN 'phosphodiesterase'
                    WHEN t.protein_class_l1 = 'Enzyme' THEN 'enzyme_other'
                    WHEN t.protein_family = 'G-protein coupled receptor' OR t.protein_class_l1 = 'Membrane receptor' THEN 'gpcr'
                    WHEN t.protein_family = 'Ion channel' OR t.protein_class_l1 = 'Ion channel' THEN 'ion_channel'
                    WHEN t.protein_family = 'Nuclear receptor' THEN 'nuclear_receptor'
                    WHEN t.protein_family = 'Transporter' OR t.protein_class_l1 = 'Transporter' THEN 'transporter'
                    WHEN t.protein_class_l1 = 'Epigenetic regulator' THEN 'epigenetic'
                    WHEN t.protein_class_l1 = 'Transcription factor' THEN 'transcription_factor'
                    WHEN t.protein_class_l1 = 'Structural protein' THEN 'structural'
                    WHEN t.protein_class_l1 IS NOT NULL THEN LOWER(REPLACE(t.protein_class_l1, ' ', '_'))
                    WHEN t.protein_family IS NOT NULL THEN LOWER(REPLACE(t.protein_family, ' ', '_'))
                    ELSE 'unclassified'
                END,
                target_subclass = COALESCE(t.protein_class_l2, t.protein_class_l3),
                has_crystal_structure = (t.pdb_ids IS NOT NULL AND array_length(t.pdb_ids, 1) > 0),
                num_pdb_structures = COALESCE(array_length(t.pdb_ids, 1), 0),
                updated_at = NOW()
            FROM pharmaco_db.targets t
            WHERE t.id = tt.target_id
        """)
        conn.commit()
        log(f"    Updated {cur.rowcount} rows")

        # Step 3: Activity statistics
        log("  Step 3: Computing activity statistics...")
        cur.execute("""
            UPDATE pharmaco_db.target_tags tt SET
                num_active_compounds = sub.n_active,
                num_total_activities = sub.n_total,
                median_pchembl = sub.med_pchembl,
                data_maturity = CASE
                    WHEN sub.n_active >= 100 THEN 'mature'
                    WHEN sub.n_active >= 10 THEN 'emerging'
                    ELSE 'sparse'
                END
            FROM (
                SELECT
                    target_id,
                    COUNT(*) FILTER (WHERE activity_class = 'active') as n_active,
                    COUNT(*) as n_total,
                    PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY pchembl_value) as med_pchembl
                FROM pharmaco_db.bioactivities
                WHERE target_id IS NOT NULL AND pchembl_value IS NOT NULL
                GROUP BY target_id
            ) sub
            WHERE sub.target_id = tt.target_id
        """)
        conn.commit()
        log(f"    Updated {cur.rowcount} rows")

        # Step 4: Druggability score (composite)
        log("  Step 4: Computing druggability scores...")
        cur.execute("""
            UPDATE pharmaco_db.target_tags tt SET
                druggability_score = (
                    -- Has crystal structure: +0.2
                    CASE WHEN tt.has_crystal_structure THEN 0.2 ELSE 0 END
                    -- Has approved drugs: +0.3
                    + CASE WHEN t.num_approved_drugs > 0 THEN 0.3 ELSE 0 END
                    -- Has active compounds: +0.2
                    + CASE WHEN tt.num_active_compounds > 10 THEN 0.2
                           WHEN tt.num_active_compounds > 0 THEN 0.1
                           ELSE 0 END
                    -- Druggable class: +0.2
                    + CASE WHEN tt.target_class IN ('kinase','gpcr','ion_channel',
                            'nuclear_receptor','protease','phosphodiesterase') THEN 0.2
                           ELSE 0 END
                    -- Data maturity: +0.1
                    + CASE WHEN tt.data_maturity = 'mature' THEN 0.1
                           WHEN tt.data_maturity = 'emerging' THEN 0.05
                           ELSE 0 END
                ),
                updated_at = NOW()
            FROM pharmaco_db.targets t
            WHERE t.id = tt.target_id
        """)
        conn.commit()
        log(f"    Updated {cur.rowcount} rows")

        # Step 5: Family size (for selectivity assessment)
        log("  Step 5: Computing family sizes...")
        cur.execute("""
            UPDATE pharmaco_db.target_tags tt SET
                family_size = sub.fam_size,
                selectivity_challenge = CASE
                    WHEN sub.fam_size > 50 THEN 'hard'
                    WHEN sub.fam_size > 10 THEN 'moderate'
                    ELSE 'easy'
                END
            FROM (
                SELECT t.id, COUNT(*) OVER (PARTITION BY t.protein_family) as fam_size
                FROM pharmaco_db.targets t
                WHERE t.protein_family IS NOT NULL
            ) sub
            WHERE sub.id = tt.target_id
        """)
        conn.commit()
        log(f"    Updated {cur.rowcount} rows")

    log("  Target tags DONE")

# ============================================================
# ACTIVITY TAGS
# ============================================================
def compute_activity_tags(conn):
    log("=" * 60)
    log("Computing activity tags (batch processing)")
    log("=" * 60)

    with conn.cursor() as cur:
        # Get total count
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.bioactivities")
        total = cur.fetchone()[0]
        log(f"  Total activities to tag: {total:,}")

        # Process in batches to avoid memory issues
        batch_size = 100000
        offset = 0

        while offset < total:
            # Step 1: Create tag records
            cur.execute("""
                INSERT INTO pharmaco_db.activity_tags (activity_id, confidence_level, pchembl_bin)
                SELECT
                    b.id,
                    CASE
                        WHEN a.confidence_score >= 8 THEN 'high'
                        WHEN a.confidence_score >= 6 THEN 'medium'
                        ELSE 'low'
                    END,
                    CASE
                        WHEN b.pchembl_value >= 8 THEN 'very_potent(>8)'
                        WHEN b.pchembl_value >= 7 THEN 'potent(7-8)'
                        WHEN b.pchembl_value >= 6 THEN 'moderate(6-7)'
                        WHEN b.pchembl_value >= 5 THEN 'weak(5-6)'
                        ELSE 'inactive(<5)'
                    END
                FROM pharmaco_db.bioactivities b
                LEFT JOIN pharmaco_db.assays a ON a.id = b.assay_id
                WHERE b.id > %s AND b.id <= %s + %s
                AND b.pchembl_value IS NOT NULL
                ON CONFLICT (activity_id) DO UPDATE SET
                    confidence_level = EXCLUDED.confidence_level,
                    pchembl_bin = EXCLUDED.pchembl_bin
            """, (offset, offset, batch_size))
            conn.commit()

            offset += batch_size
            if offset % 500000 == 0:
                log(f"    Progress: {offset:,}/{total:,}")

        # Step 2: Mark cell-based vs biochemical
        log("  Step 2: Marking assay types...")
        cur.execute("""
            UPDATE pharmaco_db.activity_tags at SET
                is_direct_binding = (a.assay_type = 'B'),
                is_cell_based = (a.assay_type = 'F')
            FROM pharmaco_db.bioactivities b
            JOIN pharmaco_db.assays a ON a.id = b.assay_id
            WHERE b.id = at.activity_id
        """)
        conn.commit()
        log(f"    Updated {cur.rowcount} rows")

    log("  Activity tags DONE")

# ============================================================
# SELECTIVITY PROFILES
# ============================================================
def compute_selectivity(conn):
    log("=" * 60)
    log("Computing selectivity profiles")
    log("=" * 60)

    with conn.cursor() as cur:
        # For compounds tested on multiple targets, compute selectivity
        log("  Finding multi-target compounds...")
        cur.execute("""
            INSERT INTO pharmaco_db.selectivity_profiles
            (compound_id, primary_target_id, off_target_id, primary_pchembl, off_target_pchembl, selectivity_index)
            SELECT
                b1.compound_id,
                b1.target_id as primary_target_id,
                b2.target_id as off_target_id,
                MAX(b1.pchembl_value) as primary_pchembl,
                MAX(b2.pchembl_value) as off_target_pchembl,
                MAX(b1.pchembl_value) - MAX(b2.pchembl_value) as selectivity_index
            FROM pharmaco_db.bioactivities b1
            JOIN pharmaco_db.bioactivities b2 ON b2.compound_id = b1.compound_id
                AND b2.target_id != b1.target_id
                AND b2.target_id IS NOT NULL
            WHERE b1.target_id IS NOT NULL
                AND b1.pchembl_value >= 6.0
                AND b2.pchembl_value IS NOT NULL
            GROUP BY b1.compound_id, b1.target_id, b2.target_id
            HAVING MAX(b1.pchembl_value) - MAX(b2.pchembl_value) IS NOT NULL
            LIMIT 500000
        """)
        conn.commit()
        log(f"    Inserted {cur.rowcount} selectivity pairs")

    log("  Selectivity profiles DONE")

def main():
    log("=" * 60)
    log("PharmacoDB — Tag Computation Engine")
    log("=" * 60)

    conn = get_conn()
    try:
        compute_compound_tags(conn)
        compute_target_tags(conn)
        compute_activity_tags(conn)
        compute_selectivity(conn)
    except Exception as e:
        log(f"FATAL ERROR: {e}")
        import traceback
        traceback.print_exc()
    finally:
        conn.close()

if __name__ == "__main__":
    main()
