-- ============================================================
-- PharmacoDB — Phase 3 Schema: Additional Public Databases
-- Sources: DGIdb, SIDER, OpenFDA, ClinicalTrials.gov, UniChem,
--          ChEBI, STRING, Reactome, DisGeNET, KEGG
-- ============================================================

SET search_path TO pharmaco_db, public;

-- ============================================================
-- 1. DGIdb — Drug-Gene Interactions
-- ============================================================
CREATE TABLE IF NOT EXISTS pharmaco_db.dgidb_interactions (
    id                  BIGSERIAL PRIMARY KEY,
    target_id           BIGINT REFERENCES pharmaco_db.targets(id),
    gene_name           TEXT NOT NULL,
    drug_name           TEXT,
    drug_chembl_id      TEXT,
    compound_id         BIGINT REFERENCES pharmaco_db.compounds(id),
    interaction_type    TEXT,
    interaction_score   REAL,
    sources             TEXT[],
    pmids               TEXT[],
    created_at          TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_dgidb_target ON pharmaco_db.dgidb_interactions(target_id);
CREATE INDEX IF NOT EXISTS idx_dgidb_compound ON pharmaco_db.dgidb_interactions(compound_id);
CREATE INDEX IF NOT EXISTS idx_dgidb_gene ON pharmaco_db.dgidb_interactions(gene_name);
CREATE INDEX IF NOT EXISTS idx_dgidb_type ON pharmaco_db.dgidb_interactions(interaction_type);
CREATE UNIQUE INDEX IF NOT EXISTS uq_dgidb_gene_drug ON pharmaco_db.dgidb_interactions(gene_name, drug_name, interaction_type)
    WHERE drug_name IS NOT NULL;

-- ============================================================
-- 2. SIDER — Side Effects
-- ============================================================
CREATE TABLE IF NOT EXISTS pharmaco_db.side_effects (
    id                  BIGSERIAL PRIMARY KEY,
    compound_id         BIGINT REFERENCES pharmaco_db.compounds(id),
    stitch_id           TEXT,
    side_effect_name    TEXT NOT NULL,
    meddra_concept_type TEXT,
    meddra_id           TEXT,
    frequency           TEXT,
    source              TEXT DEFAULT 'sider',
    created_at          TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_se_compound ON pharmaco_db.side_effects(compound_id);
CREATE INDEX IF NOT EXISTS idx_se_stitch ON pharmaco_db.side_effects(stitch_id);
CREATE INDEX IF NOT EXISTS idx_se_name ON pharmaco_db.side_effects(side_effect_name);
CREATE INDEX IF NOT EXISTS idx_se_meddra ON pharmaco_db.side_effects(meddra_id);
CREATE UNIQUE INDEX IF NOT EXISTS uq_se_stitch_meddra ON pharmaco_db.side_effects(stitch_id, meddra_id, meddra_concept_type)
    WHERE stitch_id IS NOT NULL AND meddra_id IS NOT NULL;

-- ============================================================
-- 3. Clinical Trials (ClinicalTrials.gov)
-- ============================================================
CREATE TABLE IF NOT EXISTS pharmaco_db.clinical_trials (
    id                  BIGSERIAL PRIMARY KEY,
    compound_id         BIGINT REFERENCES pharmaco_db.compounds(id),
    nct_id              TEXT NOT NULL,
    title               TEXT,
    status              TEXT,
    phase               TEXT,
    conditions          TEXT[],
    interventions       TEXT[],
    start_date          TEXT,
    completion_date     TEXT,
    enrollment          INTEGER,
    source              TEXT DEFAULT 'clinicaltrials.gov',
    created_at          TIMESTAMPTZ DEFAULT NOW()
);

CREATE UNIQUE INDEX IF NOT EXISTS uq_ct_nct ON pharmaco_db.clinical_trials(nct_id);
CREATE INDEX IF NOT EXISTS idx_ct_compound ON pharmaco_db.clinical_trials(compound_id);
CREATE INDEX IF NOT EXISTS idx_ct_status ON pharmaco_db.clinical_trials(status);
CREATE INDEX IF NOT EXISTS idx_ct_phase ON pharmaco_db.clinical_trials(phase);

-- ============================================================
-- 4. OpenFDA Adverse Events
-- ============================================================
CREATE TABLE IF NOT EXISTS pharmaco_db.adverse_events (
    id                  BIGSERIAL PRIMARY KEY,
    compound_id         BIGINT REFERENCES pharmaco_db.compounds(id),
    drug_name           TEXT,
    reaction            TEXT NOT NULL,
    count               INTEGER,
    source              TEXT DEFAULT 'openfda',
    created_at          TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_ae_compound ON pharmaco_db.adverse_events(compound_id);
CREATE INDEX IF NOT EXISTS idx_ae_drug ON pharmaco_db.adverse_events(drug_name);
CREATE INDEX IF NOT EXISTS idx_ae_reaction ON pharmaco_db.adverse_events(reaction);
CREATE UNIQUE INDEX IF NOT EXISTS uq_ae_drug_reaction ON pharmaco_db.adverse_events(drug_name, reaction)
    WHERE drug_name IS NOT NULL;

-- ============================================================
-- 5. STRING — Protein-Protein Interactions
-- ============================================================
CREATE TABLE IF NOT EXISTS pharmaco_db.protein_interactions (
    id                  BIGSERIAL PRIMARY KEY,
    target_id_1         BIGINT REFERENCES pharmaco_db.targets(id),
    target_id_2         BIGINT REFERENCES pharmaco_db.targets(id),
    gene_1              TEXT NOT NULL,
    gene_2              TEXT NOT NULL,
    combined_score      INTEGER,
    experimental_score  INTEGER,
    database_score      INTEGER,
    textmining_score    INTEGER,
    coexpression_score  INTEGER,
    source              TEXT DEFAULT 'string',
    created_at          TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_pi_target1 ON pharmaco_db.protein_interactions(target_id_1);
CREATE INDEX IF NOT EXISTS idx_pi_target2 ON pharmaco_db.protein_interactions(target_id_2);
CREATE INDEX IF NOT EXISTS idx_pi_gene1 ON pharmaco_db.protein_interactions(gene_1);
CREATE INDEX IF NOT EXISTS idx_pi_gene2 ON pharmaco_db.protein_interactions(gene_2);
CREATE INDEX IF NOT EXISTS idx_pi_score ON pharmaco_db.protein_interactions(combined_score);
CREATE UNIQUE INDEX IF NOT EXISTS uq_pi_genes ON pharmaco_db.protein_interactions(gene_1, gene_2)
    WHERE gene_1 < gene_2;

-- ============================================================
-- 6. Reactome — Biological Pathways
-- ============================================================
CREATE TABLE IF NOT EXISTS pharmaco_db.pathways (
    id                  BIGSERIAL PRIMARY KEY,
    target_id           BIGINT REFERENCES pharmaco_db.targets(id),
    uniprot_id          TEXT,
    pathway_id          TEXT NOT NULL,
    pathway_name        TEXT,
    pathway_category    TEXT,
    species             TEXT DEFAULT 'Homo sapiens',
    source              TEXT NOT NULL,
    created_at          TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_pw_target ON pharmaco_db.pathways(target_id);
CREATE INDEX IF NOT EXISTS idx_pw_uniprot ON pharmaco_db.pathways(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_pw_pathway ON pharmaco_db.pathways(pathway_id);
CREATE INDEX IF NOT EXISTS idx_pw_source ON pharmaco_db.pathways(source);
CREATE UNIQUE INDEX IF NOT EXISTS uq_pw_target_pathway ON pharmaco_db.pathways(target_id, pathway_id, source)
    WHERE target_id IS NOT NULL;

-- ============================================================
-- 7. KEGG — Drug & Pathway Info
-- ============================================================
CREATE TABLE IF NOT EXISTS pharmaco_db.kegg_drug_info (
    id                  BIGSERIAL PRIMARY KEY,
    compound_id         BIGINT REFERENCES pharmaco_db.compounds(id),
    kegg_drug_id        TEXT,
    kegg_name           TEXT,
    kegg_formula        TEXT,
    therapeutic_category TEXT[],
    pathway_ids         TEXT[],
    pathway_names       TEXT[],
    created_at          TIMESTAMPTZ DEFAULT NOW()
);

CREATE UNIQUE INDEX IF NOT EXISTS uq_kegg_drug ON pharmaco_db.kegg_drug_info(kegg_drug_id)
    WHERE kegg_drug_id IS NOT NULL;
CREATE INDEX IF NOT EXISTS idx_kegg_compound ON pharmaco_db.kegg_drug_info(compound_id);
CREATE INDEX IF NOT EXISTS idx_kegg_drug_id ON pharmaco_db.kegg_drug_info(kegg_drug_id);

-- ============================================================
-- 8. ChEBI — Chemical Ontology (stored in cross_references)
-- No new table needed, uses existing cross_references table
-- ============================================================

-- ============================================================
-- 9. DisGeNET — Disease-Gene Associations
-- (supplements existing disease_target_associations table)
-- We add a separate table for richer DisGeNET-specific data
-- ============================================================
CREATE TABLE IF NOT EXISTS pharmaco_db.disgenet_associations (
    id                  BIGSERIAL PRIMARY KEY,
    target_id           BIGINT REFERENCES pharmaco_db.targets(id),
    gene_name           TEXT NOT NULL,
    gene_ncbi_id        INTEGER,
    disease_id          TEXT NOT NULL,       -- UMLS CUI (e.g., C0006142)
    disease_name        TEXT,
    disease_type        TEXT,               -- disease, phenotype, group
    association_score   REAL,               -- GDA score 0-1
    ei_score            REAL,               -- Evidence Index
    num_pmids           INTEGER,
    num_snps            INTEGER,
    source              TEXT DEFAULT 'disgenet',
    created_at          TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_dgn_target ON pharmaco_db.disgenet_associations(target_id);
CREATE INDEX IF NOT EXISTS idx_dgn_gene ON pharmaco_db.disgenet_associations(gene_name);
CREATE INDEX IF NOT EXISTS idx_dgn_disease ON pharmaco_db.disgenet_associations(disease_id);
CREATE INDEX IF NOT EXISTS idx_dgn_score ON pharmaco_db.disgenet_associations(association_score);
CREATE UNIQUE INDEX IF NOT EXISTS uq_dgn_gene_disease ON pharmaco_db.disgenet_associations(gene_name, disease_id);

-- ============================================================
-- VIEWS for Phase 3 data
-- ============================================================

-- Drug safety profile: side effects + adverse events combined
CREATE OR REPLACE VIEW pharmaco_db.v_drug_safety_profile AS
SELECT
    c.id AS compound_id,
    c.chembl_id,
    c.pref_name AS drug_name,
    se.side_effect_name,
    se.frequency AS se_frequency,
    se.meddra_concept_type,
    ae.reaction AS adverse_reaction,
    ae.count AS adverse_event_count,
    ae.source AS ae_source
FROM pharmaco_db.compounds c
LEFT JOIN pharmaco_db.side_effects se ON se.compound_id = c.id
LEFT JOIN pharmaco_db.adverse_events ae ON ae.compound_id = c.id
WHERE se.id IS NOT NULL OR ae.id IS NOT NULL;

-- Target network: PPI + pathways combined
CREATE OR REPLACE VIEW pharmaco_db.v_target_network AS
SELECT
    t.id AS target_id,
    t.gene_name,
    t.protein_family,
    pi.gene_2 AS interactor_gene,
    pi.combined_score AS ppi_score,
    pw.pathway_name,
    pw.source AS pathway_source,
    dg.disease_name,
    dg.association_score AS disease_score
FROM pharmaco_db.targets t
LEFT JOIN pharmaco_db.protein_interactions pi ON pi.target_id_1 = t.id
LEFT JOIN pharmaco_db.pathways pw ON pw.target_id = t.id
LEFT JOIN pharmaco_db.disgenet_associations dg ON dg.target_id = t.id;
