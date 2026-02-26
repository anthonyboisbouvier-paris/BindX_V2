-- ============================================================
-- PharmacoDB — Extended UniProt Enrichment Columns
-- Run this BEFORE ingest_uniprot_full.py
-- ============================================================

SET search_path TO pharmaco_db, public;

-- Sequence & protein metadata
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS sequence_length INTEGER;
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS mass INTEGER;
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS ec_number TEXT;

-- Structural features
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS signal_peptide TEXT;
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS transmembrane_regions TEXT[];
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS active_site TEXT;
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS binding_sites TEXT[];
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS disulfide_bonds INTEGER;
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS glycosylation_sites INTEGER;
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS phosphorylation_sites INTEGER;

-- Domain & family annotations
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS domains TEXT[];
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS interpro_ids TEXT[];
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS pfam_ids TEXT[];

-- Expression & tissue
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS tissue_expression TEXT[];
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS tissue_specificity TEXT;

-- Disease associations
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS disease_associations TEXT[];
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS involvement_in_disease TEXT;

-- Keywords & evidence
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS keywords TEXT[];
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS protein_existence TEXT;

-- Cross-references to external databases
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS cross_refs_omim TEXT[];
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS cross_refs_orphanet TEXT[];
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS cross_refs_pharmgkb TEXT[];
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS cross_refs_reactome TEXT[];
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS cross_refs_string TEXT;
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS cross_refs_intact TEXT[];

-- Isoforms & naming
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS isoforms INTEGER;
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS alternative_names TEXT[];

-- Genomic location
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS chromosome TEXT;
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS gene_location TEXT;

-- Enrichment tracking (prevents re-processing sparse entries)
ALTER TABLE pharmaco_db.targets ADD COLUMN IF NOT EXISTS uniprot_enriched_at TIMESTAMPTZ;

-- Add useful indexes for common queries
CREATE INDEX IF NOT EXISTS idx_targets_ec ON pharmaco_db.targets(ec_number);
CREATE INDEX IF NOT EXISTS idx_targets_existence ON pharmaco_db.targets(protein_existence);
CREATE INDEX IF NOT EXISTS idx_targets_seq_length ON pharmaco_db.targets(sequence_length);
CREATE INDEX IF NOT EXISTS idx_targets_keywords ON pharmaco_db.targets USING gin(keywords);
CREATE INDEX IF NOT EXISTS idx_targets_diseases ON pharmaco_db.targets USING gin(disease_associations);
CREATE INDEX IF NOT EXISTS idx_targets_domains ON pharmaco_db.targets USING gin(domains);
CREATE INDEX IF NOT EXISTS idx_targets_tissue ON pharmaco_db.targets USING gin(tissue_expression);
