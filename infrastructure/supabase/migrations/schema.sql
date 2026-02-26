-- ============================================================
-- PharmacoDB — Comprehensive Pharmacological Database Schema
-- Designed for SAR analysis & Machine Learning
-- ============================================================

CREATE SCHEMA IF NOT EXISTS pharmaco_db;
SET search_path TO pharmaco_db, public;

-- Enable useful extensions
CREATE EXTENSION IF NOT EXISTS pg_trgm;

-- ============================================================
-- 1. COMPOUNDS — Chemical entities with SMILES + descriptors
-- ============================================================
CREATE TABLE pharmaco_db.compounds (
    id              BIGSERIAL PRIMARY KEY,
    canonical_smiles TEXT NOT NULL,
    inchi           TEXT,
    inchi_key       TEXT,
    -- Source IDs
    chembl_id       TEXT UNIQUE,
    pubchem_cid     BIGINT,
    drugbank_id     TEXT,
    bindingdb_id    TEXT,
    -- Names
    pref_name       TEXT,
    -- Molecular properties (precomputed for ML)
    molecular_weight REAL,
    alogp           REAL,
    hba             SMALLINT,  -- H-bond acceptors
    hbd             SMALLINT,  -- H-bond donors
    psa             REAL,      -- Polar surface area
    rtb             SMALLINT,  -- Rotatable bonds
    num_ro5_violations SMALLINT,
    aromatic_rings  SMALLINT,
    heavy_atoms     SMALLINT,
    qed_weighted    REAL,      -- Drug-likeness score
    -- Molecular formula
    molecular_formula TEXT,
    -- Flags
    is_drug         BOOLEAN DEFAULT FALSE,
    max_phase       SMALLINT DEFAULT 0,  -- Clinical phase (4 = approved)
    is_natural_product BOOLEAN,
    -- Morgan fingerprint as bit string for similarity search
    morgan_fp       TEXT,  -- Base64-encoded 2048-bit Morgan FP
    -- Timestamps
    created_at      TIMESTAMPTZ DEFAULT NOW(),
    updated_at      TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_compounds_smiles ON pharmaco_db.compounds USING hash(canonical_smiles);
CREATE INDEX idx_compounds_inchi_key ON pharmaco_db.compounds(inchi_key);
CREATE INDEX idx_compounds_chembl ON pharmaco_db.compounds(chembl_id);
CREATE INDEX idx_compounds_pubchem ON pharmaco_db.compounds(pubchem_cid);
CREATE INDEX idx_compounds_mw ON pharmaco_db.compounds(molecular_weight);
CREATE INDEX idx_compounds_alogp ON pharmaco_db.compounds(alogp);
CREATE INDEX idx_compounds_max_phase ON pharmaco_db.compounds(max_phase);
CREATE INDEX idx_compounds_name_trgm ON pharmaco_db.compounds USING gin(pref_name gin_trgm_ops);

-- ============================================================
-- 2. TARGETS — Biological targets (proteins, receptors, etc.)
-- ============================================================
CREATE TABLE pharmaco_db.targets (
    id              BIGSERIAL PRIMARY KEY,
    chembl_id       TEXT UNIQUE,
    uniprot_id      TEXT,
    -- Identity
    pref_name       TEXT NOT NULL,
    gene_name       TEXT,
    organism        TEXT DEFAULT 'Homo sapiens',
    tax_id          INTEGER,
    -- Classification
    target_type     TEXT,  -- SINGLE PROTEIN, PROTEIN COMPLEX, etc.
    protein_family  TEXT,  -- Kinase, GPCR, Ion channel, etc.
    protein_class_l1 TEXT, -- Level 1 classification
    protein_class_l2 TEXT, -- Level 2
    protein_class_l3 TEXT, -- Level 3
    -- UniProt enrichment
    uniprot_function TEXT,
    uniprot_subcellular TEXT,
    go_molecular_function TEXT[],
    go_biological_process TEXT[],
    go_cellular_component TEXT[],
    pathway_kegg    TEXT[],
    pathway_reactome TEXT[],
    -- PDB structures
    pdb_ids         TEXT[],
    -- NCBI
    ncbi_gene_id    INTEGER,
    ncbi_gene_summary TEXT,
    -- Drug target info
    is_druggable    BOOLEAN,
    num_approved_drugs INTEGER DEFAULT 0,
    -- Timestamps
    created_at      TIMESTAMPTZ DEFAULT NOW(),
    updated_at      TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_targets_chembl ON pharmaco_db.targets(chembl_id);
CREATE INDEX idx_targets_uniprot ON pharmaco_db.targets(uniprot_id);
CREATE INDEX idx_targets_gene ON pharmaco_db.targets(gene_name);
CREATE INDEX idx_targets_organism ON pharmaco_db.targets(organism);
CREATE INDEX idx_targets_family ON pharmaco_db.targets(protein_family);
CREATE INDEX idx_targets_type ON pharmaco_db.targets(target_type);
CREATE INDEX idx_targets_name_trgm ON pharmaco_db.targets USING gin(pref_name gin_trgm_ops);

-- ============================================================
-- 3. ASSAYS — Experimental protocols
-- ============================================================
CREATE TABLE pharmaco_db.assays (
    id              BIGSERIAL PRIMARY KEY,
    chembl_id       TEXT UNIQUE,
    -- Description
    description     TEXT,
    assay_type      TEXT,  -- B=Binding, F=Functional, A=ADMET, T=Toxicity, P=Physicochemical
    assay_category  TEXT,  -- confirmatory, screening, panel, etc.
    -- Source
    src_id          INTEGER,
    src_description TEXT,
    journal         TEXT,
    year            SMALLINT,
    doi             TEXT,
    -- Target link
    target_id       BIGINT REFERENCES pharmaco_db.targets(id),
    target_chembl_id TEXT,
    -- Organism / cell line
    assay_organism  TEXT,
    assay_cell_type TEXT,
    -- Confidence
    confidence_score SMALLINT,  -- 0-9, 9 = direct single protein target
    -- Timestamps
    created_at      TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_assays_chembl ON pharmaco_db.assays(chembl_id);
CREATE INDEX idx_assays_target ON pharmaco_db.assays(target_id);
CREATE INDEX idx_assays_type ON pharmaco_db.assays(assay_type);
CREATE INDEX idx_assays_confidence ON pharmaco_db.assays(confidence_score);

-- ============================================================
-- 4. BIOACTIVITIES — The core SAR data
-- ============================================================
CREATE TABLE pharmaco_db.bioactivities (
    id              BIGSERIAL PRIMARY KEY,
    -- Links
    compound_id     BIGINT NOT NULL REFERENCES pharmaco_db.compounds(id),
    target_id       BIGINT REFERENCES pharmaco_db.targets(id),
    assay_id        BIGINT REFERENCES pharmaco_db.assays(id),
    -- Activity
    activity_type   TEXT NOT NULL,  -- IC50, Ki, Kd, EC50, etc.
    relation        TEXT,           -- =, >, <, >=, <=, ~
    value           REAL,           -- In nM
    units           TEXT,           -- nM
    pchembl_value   REAL,           -- -log10(value/1e9), standardized
    -- Classification
    activity_class  TEXT,           -- active, inactive, intermediate
    -- Source
    source          TEXT NOT NULL DEFAULT 'chembl',  -- chembl, pubchem, bindingdb
    source_id       TEXT,           -- Original activity ID
    -- Data quality
    data_validity   TEXT,           -- valid, outside typical range, etc.
    potential_duplicate BOOLEAN DEFAULT FALSE,
    -- Timestamps
    created_at      TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_bioact_compound ON pharmaco_db.bioactivities(compound_id);
CREATE INDEX idx_bioact_target ON pharmaco_db.bioactivities(target_id);
CREATE INDEX idx_bioact_assay ON pharmaco_db.bioactivities(assay_id);
CREATE INDEX idx_bioact_type ON pharmaco_db.bioactivities(activity_type);
CREATE INDEX idx_bioact_pchembl ON pharmaco_db.bioactivities(pchembl_value);
CREATE INDEX idx_bioact_source ON pharmaco_db.bioactivities(source);
CREATE INDEX idx_bioact_compound_target ON pharmaco_db.bioactivities(compound_id, target_id);

-- ============================================================
-- 5. COMPOUND_STRUCTURES — Extended structural data for ML
-- ============================================================
CREATE TABLE pharmaco_db.compound_structures (
    compound_id     BIGINT PRIMARY KEY REFERENCES pharmaco_db.compounds(id),
    -- Standard representations
    molfile         TEXT,
    canonical_smiles TEXT NOT NULL,
    -- Fingerprints (for similarity & ML)
    morgan_fp_2048  TEXT,   -- Morgan radius=2, 2048 bits, base64
    maccs_fp        TEXT,   -- MACCS 166 keys, base64
    rdkit_fp        TEXT,   -- RDKit topological FP, base64
    -- Descriptors for ML (precomputed)
    num_atoms       SMALLINT,
    num_bonds       SMALLINT,
    num_rings       SMALLINT,
    num_heteroatoms SMALLINT,
    fraction_csp3   REAL,
    tpsa            REAL,
    molar_refractivity REAL,
    num_valence_electrons INTEGER,
    num_radical_electrons SMALLINT,
    formal_charge   SMALLINT,
    -- Scaffold
    murcko_scaffold TEXT,  -- Murcko generic scaffold SMILES
    -- Timestamps
    created_at      TIMESTAMPTZ DEFAULT NOW()
);

-- ============================================================
-- 6. TARGET_INTERACTIONS — Target-target relationships
-- ============================================================
CREATE TABLE pharmaco_db.target_interactions (
    id              BIGSERIAL PRIMARY KEY,
    target_id_1     BIGINT NOT NULL REFERENCES pharmaco_db.targets(id),
    target_id_2     BIGINT NOT NULL REFERENCES pharmaco_db.targets(id),
    interaction_type TEXT,  -- pathway, complex, paralog, functional
    source          TEXT,  -- string, kegg, reactome
    confidence      REAL,
    UNIQUE(target_id_1, target_id_2, interaction_type)
);

-- ============================================================
-- 7. DRUG_INDICATIONS — Therapeutic use
-- ============================================================
CREATE TABLE pharmaco_db.drug_indications (
    id              BIGSERIAL PRIMARY KEY,
    compound_id     BIGINT NOT NULL REFERENCES pharmaco_db.compounds(id),
    mesh_id         TEXT,
    mesh_heading    TEXT,
    efo_id          TEXT,
    efo_term        TEXT,
    max_phase       SMALLINT,
    indication_refs TEXT
);

CREATE INDEX idx_indications_compound ON pharmaco_db.drug_indications(compound_id);

-- ============================================================
-- 8. DRUG_MECHANISMS — Mechanism of action
-- ============================================================
CREATE TABLE pharmaco_db.drug_mechanisms (
    id              BIGSERIAL PRIMARY KEY,
    compound_id     BIGINT NOT NULL REFERENCES pharmaco_db.compounds(id),
    target_id       BIGINT REFERENCES pharmaco_db.targets(id),
    mechanism_of_action TEXT,
    action_type     TEXT,  -- INHIBITOR, AGONIST, ANTAGONIST, etc.
    direct_interaction BOOLEAN,
    molecular_mechanism TEXT,
    selectivity_comment TEXT
);

CREATE INDEX idx_mechanisms_compound ON pharmaco_db.drug_mechanisms(compound_id);
CREATE INDEX idx_mechanisms_target ON pharmaco_db.drug_mechanisms(target_id);
CREATE INDEX idx_mechanisms_action ON pharmaco_db.drug_mechanisms(action_type);

-- ============================================================
-- 9. CROSS_REFERENCES — External DB links
-- ============================================================
CREATE TABLE pharmaco_db.cross_references (
    id              BIGSERIAL PRIMARY KEY,
    compound_id     BIGINT REFERENCES pharmaco_db.compounds(id),
    target_id       BIGINT REFERENCES pharmaco_db.targets(id),
    db_name         TEXT NOT NULL,  -- PubChem, BindingDB, DrugBank, UniProt, KEGG, etc.
    db_id           TEXT NOT NULL,
    url             TEXT,
    UNIQUE(compound_id, target_id, db_name, db_id)
);

CREATE INDEX idx_xref_compound ON pharmaco_db.cross_references(compound_id);
CREATE INDEX idx_xref_target ON pharmaco_db.cross_references(target_id);
CREATE INDEX idx_xref_db ON pharmaco_db.cross_references(db_name, db_id);

-- ============================================================
-- 10. INGESTION_LOG — Track data loading progress
-- ============================================================
CREATE TABLE pharmaco_db.ingestion_log (
    id              BIGSERIAL PRIMARY KEY,
    source          TEXT NOT NULL,
    step            TEXT NOT NULL,
    status          TEXT NOT NULL DEFAULT 'running',  -- running, completed, failed
    rows_inserted   INTEGER DEFAULT 0,
    error_message   TEXT,
    started_at      TIMESTAMPTZ DEFAULT NOW(),
    completed_at    TIMESTAMPTZ
);

-- ============================================================
-- VIEWS for quick analytics
-- ============================================================

-- Compound-target activity matrix (for SAR)
CREATE VIEW pharmaco_db.v_sar_matrix AS
SELECT
    c.canonical_smiles,
    c.chembl_id AS compound_chembl_id,
    c.molecular_weight,
    c.alogp,
    t.gene_name,
    t.uniprot_id,
    t.protein_family,
    b.activity_type,
    b.pchembl_value,
    b.value AS activity_value_nm,
    b.activity_class
FROM pharmaco_db.bioactivities b
JOIN pharmaco_db.compounds c ON c.id = b.compound_id
JOIN pharmaco_db.targets t ON t.id = b.target_id
WHERE b.pchembl_value IS NOT NULL;

-- Target druggability summary
CREATE VIEW pharmaco_db.v_target_summary AS
SELECT
    t.id,
    t.pref_name,
    t.gene_name,
    t.uniprot_id,
    t.protein_family,
    t.organism,
    COUNT(DISTINCT b.compound_id) AS num_compounds,
    COUNT(b.id) AS num_activities,
    AVG(b.pchembl_value) AS avg_pchembl,
    MAX(b.pchembl_value) AS max_pchembl,
    t.num_approved_drugs
FROM pharmaco_db.targets t
LEFT JOIN pharmaco_db.bioactivities b ON b.target_id = t.id
GROUP BY t.id;

-- Compound selectivity profile
CREATE VIEW pharmaco_db.v_compound_selectivity AS
SELECT
    c.id AS compound_id,
    c.canonical_smiles,
    c.chembl_id,
    COUNT(DISTINCT b.target_id) AS num_targets,
    ARRAY_AGG(DISTINCT t.gene_name) AS target_genes,
    MIN(b.pchembl_value) AS min_pchembl,
    MAX(b.pchembl_value) AS max_pchembl,
    MAX(b.pchembl_value) - MIN(b.pchembl_value) AS selectivity_window
FROM pharmaco_db.compounds c
JOIN pharmaco_db.bioactivities b ON b.compound_id = c.id
JOIN pharmaco_db.targets t ON t.id = b.target_id
WHERE b.pchembl_value IS NOT NULL
GROUP BY c.id
HAVING COUNT(DISTINCT b.target_id) > 1;

-- ============================================================
-- FUNCTIONS for SAR queries
-- ============================================================

-- Get activity cliff pairs (similar structures, different activity)
CREATE OR REPLACE FUNCTION pharmaco_db.get_sar_pairs(
    p_target_gene TEXT,
    p_activity_type TEXT DEFAULT 'IC50',
    p_min_pchembl_diff REAL DEFAULT 2.0
)
RETURNS TABLE(
    compound_1 TEXT,
    smiles_1 TEXT,
    pchembl_1 REAL,
    compound_2 TEXT,
    smiles_2 TEXT,
    pchembl_2 REAL,
    pchembl_diff REAL
) AS $$
SELECT
    c1.chembl_id, c1.canonical_smiles, b1.pchembl_value,
    c2.chembl_id, c2.canonical_smiles, b2.pchembl_value,
    ABS(b1.pchembl_value - b2.pchembl_value) as pchembl_diff
FROM pharmaco_db.bioactivities b1
JOIN pharmaco_db.compounds c1 ON c1.id = b1.compound_id
JOIN pharmaco_db.targets t1 ON t1.id = b1.target_id
JOIN pharmaco_db.bioactivities b2 ON b2.target_id = b1.target_id
    AND b2.activity_type = b1.activity_type
    AND b2.compound_id > b1.compound_id
JOIN pharmaco_db.compounds c2 ON c2.id = b2.compound_id
WHERE t1.gene_name = p_target_gene
    AND b1.activity_type = p_activity_type
    AND b1.pchembl_value IS NOT NULL
    AND b2.pchembl_value IS NOT NULL
    AND ABS(b1.pchembl_value - b2.pchembl_value) >= p_min_pchembl_diff
ORDER BY pchembl_diff DESC
LIMIT 100;
$$ LANGUAGE SQL STABLE;
