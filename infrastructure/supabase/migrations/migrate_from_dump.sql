-- ============================================================
-- PharmacoDB — Migrate from ChEMBL 36 dump to pharmaco_db schema
-- Run inside the pharmaco container:
--   psql -U postgres -d pharmaco -f /tmp/migrate.sql
-- ============================================================

-- We need access to the chembl_36 database via dblink
CREATE EXTENSION IF NOT EXISTS dblink;

-- Create a persistent connection to chembl_36
SELECT dblink_connect('chembl', 'dbname=chembl_36 user=postgres');

-- ============================================================
-- 1. COMPOUNDS — upsert from ChEMBL dump (fills max_phase, is_drug, etc.)
-- ============================================================
\echo '>>> Updating compounds (max_phase, is_drug, drugbank, pubchem)...'

UPDATE pharmaco_db.compounds c SET
    max_phase = sub.max_phase,
    is_drug = sub.is_drug
FROM dblink('chembl', '
    SELECT molecule_chembl_id, max_phase,
           CASE WHEN max_phase >= 4 THEN true ELSE false END as is_drug
    FROM molecule_dictionary
    WHERE max_phase > 0
') AS sub(chembl_id text, max_phase int, is_drug boolean)
WHERE c.chembl_id = sub.chembl_id;

\echo '   compounds max_phase updated'

-- Insert missing compounds
\echo '>>> Inserting missing compounds...'

INSERT INTO pharmaco_db.compounds (
    chembl_id, canonical_smiles, inchi, inchi_key, pref_name,
    molecular_weight, alogp, hba, hbd, psa, rtb,
    num_ro5_violations, aromatic_rings, heavy_atoms, qed_weighted,
    molecular_formula, max_phase, is_natural_product
)
SELECT
    md.chembl_id,
    cs.canonical_smiles,
    cs.standard_inchi,
    cs.standard_inchi_key,
    md.pref_name,
    cp.mw_freebase::double precision,
    cp.alogp::double precision,
    cp.hba::int,
    cp.hbd::int,
    cp.psa::double precision,
    cp.rtb::int,
    cp.num_ro5_violations::int,
    cp.aromatic_rings::int,
    cp.heavy_atoms::int,
    cp.qed_weighted::double precision,
    cp.full_molformula,
    md.max_phase,
    CASE WHEN md.natural_product = 1 THEN true WHEN md.natural_product = 0 THEN false ELSE NULL END
FROM dblink('chembl', '
    SELECT md.chembl_id, md.pref_name, md.max_phase, md.natural_product,
           cs.canonical_smiles, cs.standard_inchi, cs.standard_inchi_key,
           cp.mw_freebase, cp.alogp, cp.hba, cp.hbd, cp.psa, cp.rtb,
           cp.num_ro5_violations, cp.aromatic_rings, cp.heavy_atoms,
           cp.qed_weighted, cp.full_molformula
    FROM molecule_dictionary md
    JOIN compound_structures cs ON cs.molregno = md.molregno
    JOIN compound_properties cp ON cp.molregno = md.molregno
    WHERE cs.canonical_smiles IS NOT NULL
      AND cp.mw_freebase >= 100 AND cp.mw_freebase <= 900
') AS sub(
    chembl_id text, pref_name text, max_phase int, natural_product int,
    canonical_smiles text, standard_inchi text, standard_inchi_key text,
    mw_freebase text, alogp text, hba text, hbd text, psa text, rtb text,
    num_ro5_violations text, aromatic_rings text, heavy_atoms text,
    qed_weighted text, full_molformula text
)
    md ON true
    cs ON true
    cp ON true
ON CONFLICT (chembl_id) DO NOTHING;

\echo '   compounds done'

-- ============================================================
-- 2. ASSAYS
-- ============================================================
\echo '>>> Inserting assays...'

-- Build target mapping
CREATE TEMP TABLE _target_map AS
SELECT chembl_id, id FROM pharmaco_db.targets;
CREATE INDEX ON _target_map(chembl_id);

INSERT INTO pharmaco_db.assays (
    chembl_id, description, assay_type, assay_category,
    src_id, target_id, target_chembl_id,
    assay_organism, assay_cell_type, confidence_score
)
SELECT
    a.chembl_id, LEFT(a.description, 500), a.assay_type, a.assay_category,
    a.src_id, tm.id, a.target_chembl_id,
    a.assay_organism, a.assay_cell_type, a.confidence_score
FROM dblink('chembl', '
    SELECT a.chembl_id, a.description, a.assay_type, a.assay_category,
           a.src_id, td.chembl_id as target_chembl_id,
           a.assay_organism, a.assay_cell_type, a.confidence_score
    FROM assays a
    LEFT JOIN target_dictionary td ON td.tid = a.tid
    WHERE a.assay_organism = ''Homo sapiens''
      AND a.confidence_score >= 4
') AS a(
    chembl_id text, description text, assay_type text, assay_category text,
    src_id int, target_chembl_id text,
    assay_organism text, assay_cell_type text, confidence_score int
)
LEFT JOIN _target_map tm ON tm.chembl_id = a.target_chembl_id
ON CONFLICT (chembl_id) DO NOTHING;

\echo '   assays done'

-- ============================================================
-- 3. BIOACTIVITIES (the big one — 24M rows, filtered to ~10M with pChEMBL)
-- ============================================================
\echo '>>> Inserting bioactivities (this will take a few minutes)...'

-- Build compound mapping
CREATE TEMP TABLE _compound_map AS
SELECT chembl_id, id FROM pharmaco_db.compounds;
CREATE INDEX ON _compound_map(chembl_id);

-- Build assay mapping
CREATE TEMP TABLE _assay_map AS
SELECT chembl_id, id FROM pharmaco_db.assays;
CREATE INDEX ON _assay_map(chembl_id);

INSERT INTO pharmaco_db.bioactivities (
    compound_id, target_id, assay_id,
    activity_type, relation, value, units, pchembl_value,
    activity_class, source, source_id,
    data_validity, potential_duplicate
)
SELECT
    cm.id,
    tm.id,
    am.id,
    act.standard_type,
    act.standard_relation,
    act.standard_value,
    act.standard_units,
    act.pchembl_value,
    CASE
        WHEN act.pchembl_value >= 7 THEN 'active'
        WHEN act.pchembl_value >= 5 THEN 'intermediate'
        WHEN act.pchembl_value IS NOT NULL THEN 'inactive'
        ELSE NULL
    END,
    'ChEMBL',
    act.activity_id::text,
    act.data_validity_comment,
    act.potential_duplicate
FROM dblink('chembl', '
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
') AS act(
    activity_id bigint, standard_type text, standard_relation text,
    standard_value double precision, standard_units text, pchembl_value double precision,
    data_validity_comment text, potential_duplicate boolean,
    mol_chembl_id text, assay_chembl_id text, target_chembl_id text
)
LEFT JOIN _compound_map cm ON cm.chembl_id = act.mol_chembl_id
LEFT JOIN _assay_map am ON am.chembl_id = act.assay_chembl_id
LEFT JOIN _target_map tm ON tm.chembl_id = act.target_chembl_id
WHERE cm.id IS NOT NULL;

\echo '   bioactivities done'

-- ============================================================
-- 4. DRUG MECHANISMS
-- ============================================================
\echo '>>> Inserting drug mechanisms...'

INSERT INTO pharmaco_db.drug_mechanisms (
    compound_id, target_id, mechanism_of_action,
    action_type, direct_interaction, molecular_mechanism,
    selectivity_comment
)
SELECT
    cm.id,
    tm.id,
    dm.mechanism_of_action,
    dm.action_type,
    dm.direct_interaction,
    dm.molecular_mechanism,
    dm.selectivity_comment
FROM dblink('chembl', '
    SELECT md.chembl_id as mol_chembl_id,
           td.chembl_id as target_chembl_id,
           dm.mechanism_of_action, dm.action_type,
           dm.direct_interaction, dm.molecular_mechanism,
           dm.selectivity_comment
    FROM drug_mechanism dm
    JOIN molecule_dictionary md ON md.molregno = dm.molregno
    LEFT JOIN target_dictionary td ON td.tid = dm.tid
') AS dm(
    mol_chembl_id text, target_chembl_id text,
    mechanism_of_action text, action_type text,
    direct_interaction boolean, molecular_mechanism text,
    selectivity_comment text
)
LEFT JOIN _compound_map cm ON cm.chembl_id = dm.mol_chembl_id
LEFT JOIN _target_map tm ON tm.chembl_id = dm.target_chembl_id;

\echo '   drug_mechanisms done'

-- ============================================================
-- 5. DRUG INDICATIONS
-- ============================================================
\echo '>>> Inserting drug indications...'

INSERT INTO pharmaco_db.drug_indications (
    compound_id, mesh_id, mesh_heading,
    efo_id, efo_term, max_phase
)
SELECT
    cm.id,
    di.mesh_id,
    di.mesh_heading,
    di.efo_id,
    di.efo_term,
    di.max_phase_for_ind
FROM dblink('chembl', '
    SELECT md.chembl_id as mol_chembl_id,
           di.mesh_id, di.mesh_heading,
           di.efo_id, di.efo_term,
           di.max_phase_for_ind
    FROM drug_indication di
    JOIN molecule_dictionary md ON md.molregno = di.molregno
') AS di(
    mol_chembl_id text,
    mesh_id text, mesh_heading text,
    efo_id text, efo_term text,
    max_phase_for_ind int
)
LEFT JOIN _compound_map cm ON cm.chembl_id = di.mol_chembl_id
WHERE cm.id IS NOT NULL;

\echo '   drug_indications done'

-- ============================================================
-- 6. Update druggability fields on targets
-- ============================================================
\echo '>>> Updating target druggability...'

UPDATE pharmaco_db.targets t SET
    num_approved_drugs = sub.n_approved,
    is_druggable = sub.n_approved > 0
FROM (
    SELECT target_id, COUNT(DISTINCT compound_id) as n_approved
    FROM pharmaco_db.drug_mechanisms
    WHERE target_id IS NOT NULL
    GROUP BY target_id
) sub
WHERE t.id = sub.target_id;

\echo '   targets druggability updated'

-- Cleanup
DROP TABLE IF EXISTS _target_map;
DROP TABLE IF EXISTS _compound_map;
DROP TABLE IF EXISTS _assay_map;

SELECT dblink_disconnect('chembl');

-- ============================================================
-- FINAL STATS
-- ============================================================
\echo ''
\echo '============================================================'
\echo 'MIGRATION COMPLETE — Final counts:'
\echo '============================================================'

SELECT 'compounds' AS table_name, COUNT(*) AS rows FROM pharmaco_db.compounds
UNION ALL SELECT 'targets', COUNT(*) FROM pharmaco_db.targets
UNION ALL SELECT 'assays', COUNT(*) FROM pharmaco_db.assays
UNION ALL SELECT 'bioactivities', COUNT(*) FROM pharmaco_db.bioactivities
UNION ALL SELECT 'drug_mechanisms', COUNT(*) FROM pharmaco_db.drug_mechanisms
UNION ALL SELECT 'drug_indications', COUNT(*) FROM pharmaco_db.drug_indications
UNION ALL SELECT 'cross_references', COUNT(*) FROM pharmaco_db.cross_references
ORDER BY 1;

SELECT pg_size_pretty(pg_database_size('pharmaco')) AS db_size;
