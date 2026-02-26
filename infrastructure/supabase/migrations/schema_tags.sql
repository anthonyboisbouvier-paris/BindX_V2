-- ============================================================
-- PharmacoDB — Extended Schema: ML Tags & Additional Sources
-- ============================================================

SET search_path TO pharmaco_db, public;

-- ============================================================
-- COMPOUND TAGS — Comprehensive ML-ready annotations
-- ============================================================
CREATE TABLE IF NOT EXISTS pharmaco_db.compound_tags (
    compound_id     BIGINT PRIMARY KEY REFERENCES pharmaco_db.compounds(id),

    -- ===== DRUG-LIKENESS RULES =====
    lipinski_pass   BOOLEAN,          -- Ro5: MW<=500, logP<=5, HBD<=5, HBA<=10
    lipinski_violations SMALLINT,
    veber_pass      BOOLEAN,          -- TPSA<=140, RB<=10
    lead_like       BOOLEAN,          -- MW 200-350, logP -1..3, HBD<=3, HBA<=6
    fragment_like   BOOLEAN,          -- Ro3: MW<=300, logP<=3, HBD<=3, HBA<=3, RB<=3
    drug_like       BOOLEAN,          -- QED >= 0.5
    ppi_inhibitor_like BOOLEAN,       -- MW>400, logP>4, aromatic_rings>=4
    brenk_alerts    SMALLINT DEFAULT 0,  -- Count of Brenk structural alerts
    pains_alerts    SMALLINT DEFAULT 0,  -- Count of PAINS substructures
    has_reactive_group BOOLEAN DEFAULT FALSE,

    -- ===== CNS PENETRATION =====
    cns_mpo_score   REAL,             -- Multiparameter optimization score (0-6)
    cns_penetrant   BOOLEAN,          -- cns_mpo >= 4
    bbb_permeable   BOOLEAN,          -- Predicted BBB permeability

    -- ===== PHYSICOCHEMICAL BINS (for stratified ML) =====
    mw_bin          TEXT,             -- '<200','200-350','350-500','>500'
    logp_bin        TEXT,             -- '<0','0-2','2-4','>4'
    tpsa_bin        TEXT,             -- '<60','60-90','90-140','>140'
    complexity_bin  TEXT,             -- 'low','medium','high','very_high'

    -- ===== STRUCTURAL CLASSIFICATION =====
    murcko_scaffold TEXT,             -- Generic Murcko scaffold SMILES
    scaffold_cluster INTEGER,         -- Butina clustering ID
    chemical_series TEXT,             -- Named series if identified
    privileged_scaffold BOOLEAN,      -- Contains known privileged motif
    num_stereocenters SMALLINT,
    fraction_csp3   REAL,             -- sp3 fraction (3D complexity)

    -- ===== CLINICAL STATUS =====
    clinical_phase  SMALLINT DEFAULT 0,  -- 0-4
    is_approved     BOOLEAN DEFAULT FALSE,
    is_withdrawn    BOOLEAN DEFAULT FALSE,
    is_prodrug      BOOLEAN DEFAULT FALSE,

    -- ===== ATC CLASSIFICATION =====
    atc_code_l1     TEXT,             -- Anatomical main group (e.g., L=Antineoplastic)
    atc_code_l2     TEXT,             -- Therapeutic subgroup
    atc_code_l3     TEXT,             -- Pharmacological subgroup
    atc_code_l4     TEXT,             -- Chemical subgroup
    atc_code_l5     TEXT,             -- Chemical substance

    -- ===== THERAPEUTIC AREA =====
    therapeutic_areas TEXT[],          -- ['oncology','neurology',...]

    -- ===== METABOLISM =====
    cyp_substrate   TEXT[],           -- ['CYP3A4','CYP2D6',...]
    cyp_inhibitor   TEXT[],           -- CYPs this compound inhibits
    elimination_route TEXT,           -- hepatic, renal, mixed

    -- ===== TOXICOLOGY FLAGS =====
    is_herg_liability BOOLEAN,        -- hERG IC50 < 10 uM
    is_hepatotoxic    BOOLEAN,        -- Known hepatotoxicity
    is_mutagenic      BOOLEAN,        -- Ames test positive
    is_carcinogenic   BOOLEAN,
    tox_class         TEXT,           -- GHS toxicity class (1-5)

    -- ===== DATA RICHNESS =====
    num_bioactivities INTEGER DEFAULT 0,
    num_targets       INTEGER DEFAULT 0,
    num_assays        INTEGER DEFAULT 0,
    data_completeness REAL,            -- 0-1 score of how complete the data is

    created_at      TIMESTAMPTZ DEFAULT NOW(),
    updated_at      TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_ctags_lipinski ON pharmaco_db.compound_tags(lipinski_pass);
CREATE INDEX idx_ctags_lead ON pharmaco_db.compound_tags(lead_like);
CREATE INDEX idx_ctags_drug ON pharmaco_db.compound_tags(drug_like);
CREATE INDEX idx_ctags_scaffold ON pharmaco_db.compound_tags(murcko_scaffold);
CREATE INDEX idx_ctags_phase ON pharmaco_db.compound_tags(clinical_phase);
CREATE INDEX idx_ctags_mw_bin ON pharmaco_db.compound_tags(mw_bin);
CREATE INDEX idx_ctags_logp_bin ON pharmaco_db.compound_tags(logp_bin);

-- ============================================================
-- TARGET TAGS — ML-ready target annotations
-- ============================================================
CREATE TABLE IF NOT EXISTS pharmaco_db.target_tags (
    target_id       BIGINT PRIMARY KEY REFERENCES pharmaco_db.targets(id),

    -- ===== CLASSIFICATION =====
    target_class    TEXT,             -- kinase, gpcr, ion_channel, nuclear_receptor, etc.
    target_subclass TEXT,             -- tyrosine_kinase, serine_threonine_kinase, etc.
    enzyme_class_ec TEXT,             -- EC number (e.g., 2.7.10.1)

    -- ===== DRUGGABILITY =====
    druggability_score REAL,          -- 0-1 composite score
    has_crystal_structure BOOLEAN,
    num_pdb_structures INTEGER DEFAULT 0,
    has_active_site_known BOOLEAN,
    has_allosteric_site BOOLEAN,
    binding_site_druggability TEXT,   -- 'high','medium','low'

    -- ===== DISEASE ASSOCIATION =====
    disease_areas   TEXT[],           -- ['cancer','inflammation',...]
    disease_ids     TEXT[],           -- EFO/MONDO IDs
    genetic_association BOOLEAN,      -- GWAS/genetic evidence
    somatic_mutation BOOLEAN,         -- Cancer driver
    known_lof_phenotype BOOLEAN,

    -- ===== EXPRESSION =====
    tissue_expression TEXT[],         -- Tissues with high expression
    tissue_specificity TEXT,          -- 'ubiquitous','tissue_specific','tissue_enriched'

    -- ===== ESSENTIALITY =====
    is_essential    BOOLEAN,          -- DepMap essentiality
    essentiality_score REAL,

    -- ===== SAFETY =====
    safety_risk_level TEXT,           -- 'low','medium','high'
    has_safety_concern BOOLEAN,
    safety_concerns TEXT[],           -- Specific concerns

    -- ===== SELECTIVITY =====
    family_size     INTEGER,          -- Number of close homologs
    selectivity_challenge TEXT,       -- 'easy','moderate','hard'

    -- ===== DATA RICHNESS =====
    num_active_compounds INTEGER DEFAULT 0,
    num_total_activities INTEGER DEFAULT 0,
    num_approved_drugs INTEGER DEFAULT 0,
    median_pchembl REAL,
    data_maturity   TEXT,             -- 'mature','emerging','sparse'

    created_at      TIMESTAMPTZ DEFAULT NOW(),
    updated_at      TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_ttags_class ON pharmaco_db.target_tags(target_class);
CREATE INDEX idx_ttags_druggability ON pharmaco_db.target_tags(druggability_score);
CREATE INDEX idx_ttags_safety ON pharmaco_db.target_tags(safety_risk_level);

-- ============================================================
-- ACTIVITY TAGS — Per-activity annotations for ML
-- ============================================================
CREATE TABLE IF NOT EXISTS pharmaco_db.activity_tags (
    activity_id     BIGINT PRIMARY KEY REFERENCES pharmaco_db.bioactivities(id),

    -- ===== QUALITY =====
    confidence_level TEXT,            -- 'high','medium','low'
    is_direct_binding BOOLEAN,        -- Direct binding assay
    is_cell_based   BOOLEAN,          -- Cell-based vs biochemical
    is_in_vivo      BOOLEAN,          -- In vivo data

    -- ===== CLASSIFICATION =====
    pchembl_bin     TEXT,             -- 'inactive(<5)','weak(5-6)','moderate(6-7)','potent(7-8)','very_potent(>8)'
    selectivity_ratio REAL,           -- Activity vs closest off-target

    -- ===== SAR FLAGS =====
    is_activity_cliff BOOLEAN,        -- Part of an activity cliff pair
    cliff_partner_id BIGINT,          -- The other compound in the cliff
    cliff_fold_change REAL,

    created_at      TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_atags_confidence ON pharmaco_db.activity_tags(confidence_level);
CREATE INDEX idx_atags_pchembl_bin ON pharmaco_db.activity_tags(pchembl_bin);
CREATE INDEX idx_atags_cliff ON pharmaco_db.activity_tags(is_activity_cliff);

-- ============================================================
-- DISEASE-TARGET ASSOCIATIONS (Open Targets)
-- ============================================================
CREATE TABLE IF NOT EXISTS pharmaco_db.disease_target_associations (
    id              BIGSERIAL PRIMARY KEY,
    target_id       BIGINT REFERENCES pharmaco_db.targets(id),
    disease_id      TEXT NOT NULL,     -- EFO/MONDO ID
    disease_name    TEXT,
    therapeutic_area TEXT,
    overall_score   REAL,              -- Open Targets overall association score
    genetic_score   REAL,
    somatic_score   REAL,
    known_drug_score REAL,
    literature_score REAL,
    rna_expression_score REAL,
    animal_model_score REAL,
    source          TEXT DEFAULT 'open_targets',
    created_at      TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_dta_target ON pharmaco_db.disease_target_associations(target_id);
CREATE INDEX idx_dta_disease ON pharmaco_db.disease_target_associations(disease_id);
CREATE INDEX idx_dta_score ON pharmaco_db.disease_target_associations(overall_score);
CREATE INDEX idx_dta_area ON pharmaco_db.disease_target_associations(therapeutic_area);

-- ============================================================
-- TOXICOLOGY DATA (ToxCast/Tox21)
-- ============================================================
CREATE TABLE IF NOT EXISTS pharmaco_db.toxicology_data (
    id              BIGSERIAL PRIMARY KEY,
    compound_id     BIGINT REFERENCES pharmaco_db.compounds(id),
    assay_name      TEXT,
    assay_endpoint  TEXT,              -- e.g., 'NR-AR', 'SR-MMP', 'NR-ER'
    activity_outcome TEXT,             -- active/inactive
    ac50_um         REAL,              -- AC50 in uM
    efficacy        REAL,
    source          TEXT DEFAULT 'tox21',
    created_at      TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_tox_compound ON pharmaco_db.toxicology_data(compound_id);
CREATE INDEX idx_tox_endpoint ON pharmaco_db.toxicology_data(assay_endpoint);

-- ============================================================
-- PHARMACOKINETICS DATA
-- ============================================================
CREATE TABLE IF NOT EXISTS pharmaco_db.pharmacokinetics (
    id              BIGSERIAL PRIMARY KEY,
    compound_id     BIGINT REFERENCES pharmaco_db.compounds(id),
    species         TEXT,              -- human, rat, mouse, dog
    route           TEXT,              -- oral, iv, ip, sc
    parameter       TEXT,              -- Cmax, AUC, t1/2, F%, Vd, CL
    value           REAL,
    units           TEXT,
    source          TEXT,
    reference       TEXT,
    created_at      TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_pk_compound ON pharmaco_db.pharmacokinetics(compound_id);
CREATE INDEX idx_pk_species ON pharmaco_db.pharmacokinetics(species);
CREATE INDEX idx_pk_param ON pharmaco_db.pharmacokinetics(parameter);

-- ============================================================
-- GUIDE TO PHARMACOLOGY DATA
-- ============================================================
CREATE TABLE IF NOT EXISTS pharmaco_db.gtop_interactions (
    id              BIGSERIAL PRIMARY KEY,
    compound_id     BIGINT REFERENCES pharmaco_db.compounds(id),
    target_id       BIGINT REFERENCES pharmaco_db.targets(id),
    ligand_id       INTEGER,           -- GtoP ligand ID
    target_gtop_id  INTEGER,           -- GtoP target ID
    interaction_type TEXT,              -- agonist, antagonist, inhibitor, etc.
    affinity_type   TEXT,              -- pKi, pIC50, pEC50, pKd
    affinity_value  REAL,              -- -log10(M) value
    affinity_high   REAL,              -- Range high
    affinity_low    REAL,              -- Range low
    species         TEXT,
    endogenous      BOOLEAN DEFAULT FALSE,
    primary_target  BOOLEAN DEFAULT FALSE,
    source          TEXT DEFAULT 'guide_to_pharmacology',
    created_at      TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_gtop_compound ON pharmaco_db.gtop_interactions(compound_id);
CREATE INDEX idx_gtop_target ON pharmaco_db.gtop_interactions(target_id);
CREATE INDEX idx_gtop_type ON pharmaco_db.gtop_interactions(interaction_type);

-- ============================================================
-- CHEMBL ATC CLASSIFICATION
-- ============================================================
CREATE TABLE IF NOT EXISTS pharmaco_db.atc_classification (
    id              BIGSERIAL PRIMARY KEY,
    compound_id     BIGINT REFERENCES pharmaco_db.compounds(id),
    atc_code        TEXT NOT NULL,
    level1          TEXT,              -- Anatomical (A, B, C, D, G, H, J, L, M, N, P, R, S, V)
    level1_desc     TEXT,
    level2          TEXT,
    level2_desc     TEXT,
    level3          TEXT,
    level3_desc     TEXT,
    level4          TEXT,
    level4_desc     TEXT,
    level5          TEXT,
    level5_desc     TEXT,
    who_name        TEXT,
    source          TEXT DEFAULT 'chembl',
    created_at      TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_atc_compound ON pharmaco_db.atc_classification(compound_id);
CREATE INDEX idx_atc_code ON pharmaco_db.atc_classification(atc_code);
CREATE INDEX idx_atc_l1 ON pharmaco_db.atc_classification(level1);

-- ============================================================
-- SELECTIVITY DATA — Compound selectivity panels
-- ============================================================
CREATE TABLE IF NOT EXISTS pharmaco_db.selectivity_profiles (
    id              BIGSERIAL PRIMARY KEY,
    compound_id     BIGINT NOT NULL REFERENCES pharmaco_db.compounds(id),
    primary_target_id BIGINT REFERENCES pharmaco_db.targets(id),
    off_target_id   BIGINT REFERENCES pharmaco_db.targets(id),
    primary_pchembl REAL,
    off_target_pchembl REAL,
    selectivity_index REAL,            -- primary_pchembl - off_target_pchembl
    created_at      TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_sel_compound ON pharmaco_db.selectivity_profiles(compound_id);
CREATE INDEX idx_sel_primary ON pharmaco_db.selectivity_profiles(primary_target_id);

-- ============================================================
-- ML FEATURE VIEWS
-- ============================================================

-- Complete compound feature vector for ML
CREATE OR REPLACE VIEW pharmaco_db.v_ml_compound_features AS
SELECT
    c.id, c.chembl_id, c.canonical_smiles,
    c.molecular_weight, c.alogp, c.hba, c.hbd, c.psa, c.rtb,
    c.num_ro5_violations, c.aromatic_rings, c.heavy_atoms, c.qed_weighted,
    ct.lipinski_pass, ct.veber_pass, ct.lead_like, ct.fragment_like,
    ct.drug_like, ct.ppi_inhibitor_like, ct.cns_mpo_score, ct.cns_penetrant,
    ct.mw_bin, ct.logp_bin, ct.tpsa_bin, ct.complexity_bin,
    ct.fraction_csp3, ct.num_stereocenters,
    ct.pains_alerts, ct.brenk_alerts, ct.has_reactive_group,
    ct.clinical_phase, ct.is_approved,
    ct.num_bioactivities, ct.num_targets, ct.data_completeness
FROM pharmaco_db.compounds c
LEFT JOIN pharmaco_db.compound_tags ct ON ct.compound_id = c.id;

-- Complete target feature vector for ML
CREATE OR REPLACE VIEW pharmaco_db.v_ml_target_features AS
SELECT
    t.id, t.chembl_id, t.uniprot_id, t.gene_name, t.pref_name,
    t.protein_family, t.protein_class_l1, t.protein_class_l2,
    t.organism, t.is_druggable, t.num_approved_drugs,
    tt.target_class, tt.target_subclass,
    tt.druggability_score, tt.has_crystal_structure, tt.num_pdb_structures,
    tt.tissue_specificity, tt.is_essential,
    tt.safety_risk_level, tt.family_size, tt.selectivity_challenge,
    tt.num_active_compounds, tt.num_total_activities, tt.median_pchembl,
    tt.data_maturity
FROM pharmaco_db.targets t
LEFT JOIN pharmaco_db.target_tags tt ON tt.target_id = t.id;

-- SAR-ready dataset: compound-target pairs with all features
CREATE OR REPLACE VIEW pharmaco_db.v_ml_sar_dataset AS
SELECT
    b.id AS activity_id,
    -- Compound features
    c.canonical_smiles, c.molecular_weight, c.alogp, c.hba, c.hbd,
    c.psa, c.rtb, c.aromatic_rings, c.heavy_atoms, c.qed_weighted,
    ct.fraction_csp3, ct.mw_bin, ct.logp_bin,
    -- Target features
    t.gene_name, t.protein_family,
    tt.target_class, tt.druggability_score,
    -- Activity
    b.activity_type, b.pchembl_value, b.activity_class,
    -- Tags
    at.confidence_level, at.pchembl_bin, at.is_activity_cliff
FROM pharmaco_db.bioactivities b
JOIN pharmaco_db.compounds c ON c.id = b.compound_id
LEFT JOIN pharmaco_db.compound_tags ct ON ct.compound_id = c.id
LEFT JOIN pharmaco_db.targets t ON t.id = b.target_id
LEFT JOIN pharmaco_db.target_tags tt ON tt.target_id = t.id
LEFT JOIN pharmaco_db.activity_tags at ON at.activity_id = b.id
WHERE b.pchembl_value IS NOT NULL;
