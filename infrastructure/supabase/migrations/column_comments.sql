-- ============================================================
-- PharmacoDB — Column descriptions for all tables
-- Visible in Supabase dashboard + any SQL client
-- ============================================================

-- ─── TARGETS ─────────────────────────────────────────────────
COMMENT ON TABLE pharmaco_db.targets IS 'Protein targets from ChEMBL + UniProt enrichment. 56 columns covering function, structure, pathways, disease, and expression.';

COMMENT ON COLUMN pharmaco_db.targets.id IS 'Internal primary key';
COMMENT ON COLUMN pharmaco_db.targets.chembl_id IS 'ChEMBL target identifier (e.g. CHEMBL203)';
COMMENT ON COLUMN pharmaco_db.targets.uniprot_id IS 'UniProt accession (e.g. P00533)';
COMMENT ON COLUMN pharmaco_db.targets.pref_name IS 'Preferred target name from ChEMBL';
COMMENT ON COLUMN pharmaco_db.targets.gene_name IS 'HGNC gene symbol (e.g. EGFR, BRAF)';
COMMENT ON COLUMN pharmaco_db.targets.organism IS 'Source organism (e.g. Homo sapiens)';
COMMENT ON COLUMN pharmaco_db.targets.tax_id IS 'NCBI taxonomy ID (9606 = human)';
COMMENT ON COLUMN pharmaco_db.targets.target_type IS 'ChEMBL target type (SINGLE PROTEIN, PROTEIN COMPLEX, etc.)';
COMMENT ON COLUMN pharmaco_db.targets.protein_family IS 'Broad protein family (Kinase, GPCR, Ion channel, Protease, etc.)';
COMMENT ON COLUMN pharmaco_db.targets.protein_class_l1 IS 'ChEMBL protein classification level 1';
COMMENT ON COLUMN pharmaco_db.targets.protein_class_l2 IS 'ChEMBL protein classification level 2';
COMMENT ON COLUMN pharmaco_db.targets.protein_class_l3 IS 'ChEMBL protein classification level 3';
COMMENT ON COLUMN pharmaco_db.targets.uniprot_function IS 'Protein function description from UniProt';
COMMENT ON COLUMN pharmaco_db.targets.uniprot_subcellular IS 'Subcellular localization from UniProt';
COMMENT ON COLUMN pharmaco_db.targets.go_molecular_function IS 'Gene Ontology molecular function terms (array)';
COMMENT ON COLUMN pharmaco_db.targets.go_biological_process IS 'Gene Ontology biological process terms (array)';
COMMENT ON COLUMN pharmaco_db.targets.go_cellular_component IS 'Gene Ontology cellular component terms (array)';
COMMENT ON COLUMN pharmaco_db.targets.pathway_kegg IS 'KEGG pathway memberships (array)';
COMMENT ON COLUMN pharmaco_db.targets.pathway_reactome IS 'Reactome pathway memberships (array)';
COMMENT ON COLUMN pharmaco_db.targets.pdb_ids IS 'PDB structure IDs (array). Count = structural tractability proxy';
COMMENT ON COLUMN pharmaco_db.targets.ncbi_gene_id IS 'NCBI Gene/Entrez ID';
COMMENT ON COLUMN pharmaco_db.targets.ncbi_gene_summary IS 'Gene summary from NCBI';
COMMENT ON COLUMN pharmaco_db.targets.is_druggable IS 'Has at least 1 known drug mechanism in ChEMBL';
COMMENT ON COLUMN pharmaco_db.targets.num_approved_drugs IS 'Count of distinct compounds with drug mechanisms for this target';
COMMENT ON COLUMN pharmaco_db.targets.sequence_length IS 'Amino acid sequence length';
COMMENT ON COLUMN pharmaco_db.targets.mass IS 'Protein molecular mass in Da';
COMMENT ON COLUMN pharmaco_db.targets.ec_number IS 'Enzyme Commission number (e.g. 2.7.10.1 for tyrosine kinases)';
COMMENT ON COLUMN pharmaco_db.targets.signal_peptide IS 'Signal peptide annotation from UniProt';
COMMENT ON COLUMN pharmaco_db.targets.transmembrane_regions IS 'Transmembrane region annotations (array)';
COMMENT ON COLUMN pharmaco_db.targets.active_site IS 'Active site residue annotation from UniProt';
COMMENT ON COLUMN pharmaco_db.targets.binding_sites IS 'Binding site annotations (array)';
COMMENT ON COLUMN pharmaco_db.targets.disulfide_bonds IS 'Disulfide bond annotations';
COMMENT ON COLUMN pharmaco_db.targets.glycosylation_sites IS 'Glycosylation site annotations';
COMMENT ON COLUMN pharmaco_db.targets.phosphorylation_sites IS 'Phosphorylation site annotations';
COMMENT ON COLUMN pharmaco_db.targets.domains IS 'Protein domain annotations (array)';
COMMENT ON COLUMN pharmaco_db.targets.interpro_ids IS 'InterPro domain IDs (array)';
COMMENT ON COLUMN pharmaco_db.targets.pfam_ids IS 'Pfam domain IDs (array)';
COMMENT ON COLUMN pharmaco_db.targets.tissue_expression IS 'Tissue expression from UniProt (array of tissue names)';
COMMENT ON COLUMN pharmaco_db.targets.tissue_specificity IS 'Tissue specificity description from UniProt';
COMMENT ON COLUMN pharmaco_db.targets.disease_associations IS 'Disease associations from UniProt (array)';
COMMENT ON COLUMN pharmaco_db.targets.involvement_in_disease IS 'Disease involvement narrative from UniProt';
COMMENT ON COLUMN pharmaco_db.targets.keywords IS 'UniProt keyword annotations (array)';
COMMENT ON COLUMN pharmaco_db.targets.protein_existence IS 'UniProt protein existence evidence level (1=protein, 5=uncertain)';
COMMENT ON COLUMN pharmaco_db.targets.cross_refs_omim IS 'OMIM cross-reference IDs (array)';
COMMENT ON COLUMN pharmaco_db.targets.cross_refs_orphanet IS 'Orphanet cross-reference IDs (array)';
COMMENT ON COLUMN pharmaco_db.targets.cross_refs_pharmgkb IS 'PharmGKB cross-reference IDs (array)';
COMMENT ON COLUMN pharmaco_db.targets.cross_refs_reactome IS 'Reactome cross-reference IDs (array)';
COMMENT ON COLUMN pharmaco_db.targets.cross_refs_string IS 'STRING PPI network IDs (array)';
COMMENT ON COLUMN pharmaco_db.targets.cross_refs_intact IS 'IntAct interaction IDs (array)';
COMMENT ON COLUMN pharmaco_db.targets.isoforms IS 'Protein isoform annotations (array)';
COMMENT ON COLUMN pharmaco_db.targets.alternative_names IS 'Alternative protein names from UniProt (array)';
COMMENT ON COLUMN pharmaco_db.targets.chromosome IS 'Chromosomal location';
COMMENT ON COLUMN pharmaco_db.targets.gene_location IS 'Genomic location details';
COMMENT ON COLUMN pharmaco_db.targets.uniprot_enriched_at IS 'Timestamp of UniProt enrichment completion';

-- ─── COMPOUNDS ───────────────────────────────────────────────
COMMENT ON TABLE pharmaco_db.compounds IS 'Small molecule compounds from ChEMBL. MW 100-900 Da with SMILES. Physicochemical properties pre-computed.';

COMMENT ON COLUMN pharmaco_db.compounds.id IS 'Internal primary key';
COMMENT ON COLUMN pharmaco_db.compounds.canonical_smiles IS 'Canonical SMILES string (RDKit-parseable)';
COMMENT ON COLUMN pharmaco_db.compounds.inchi IS 'IUPAC InChI structural representation';
COMMENT ON COLUMN pharmaco_db.compounds.inchi_key IS 'InChI key (27-char hash for fast lookup)';
COMMENT ON COLUMN pharmaco_db.compounds.chembl_id IS 'ChEMBL compound identifier (e.g. CHEMBL25)';
COMMENT ON COLUMN pharmaco_db.compounds.pubchem_cid IS 'PubChem Compound ID';
COMMENT ON COLUMN pharmaco_db.compounds.drugbank_id IS 'DrugBank identifier';
COMMENT ON COLUMN pharmaco_db.compounds.bindingdb_id IS 'BindingDB identifier';
COMMENT ON COLUMN pharmaco_db.compounds.pref_name IS 'Preferred compound name (INN/USAN where available)';
COMMENT ON COLUMN pharmaco_db.compounds.molecular_weight IS 'Molecular weight of freebase (Da)';
COMMENT ON COLUMN pharmaco_db.compounds.alogp IS 'Calculated partition coefficient (ALogP). Drug-like: -0.4 to 5.6';
COMMENT ON COLUMN pharmaco_db.compounds.hba IS 'Hydrogen bond acceptors. Lipinski: <= 10';
COMMENT ON COLUMN pharmaco_db.compounds.hbd IS 'Hydrogen bond donors. Lipinski: <= 5';
COMMENT ON COLUMN pharmaco_db.compounds.psa IS 'Polar surface area (A^2). Oral absorption: < 140';
COMMENT ON COLUMN pharmaco_db.compounds.rtb IS 'Rotatable bonds. Oral bioavailability: <= 10';
COMMENT ON COLUMN pharmaco_db.compounds.num_ro5_violations IS 'Number of Lipinski Rule of Five violations (0-4). Drug-like: <= 1';
COMMENT ON COLUMN pharmaco_db.compounds.aromatic_rings IS 'Number of aromatic rings';
COMMENT ON COLUMN pharmaco_db.compounds.heavy_atoms IS 'Number of heavy (non-hydrogen) atoms';
COMMENT ON COLUMN pharmaco_db.compounds.qed_weighted IS 'Quantitative Estimate of Drug-likeness (0-1). Higher = more drug-like';
COMMENT ON COLUMN pharmaco_db.compounds.molecular_formula IS 'Molecular formula (e.g. C20H25N3O)';
COMMENT ON COLUMN pharmaco_db.compounds.is_drug IS 'Has reached clinical Phase IV (approved drug)';
COMMENT ON COLUMN pharmaco_db.compounds.max_phase IS 'Maximum clinical trial phase reached (0-4). 4 = approved';
COMMENT ON COLUMN pharmaco_db.compounds.is_natural_product IS 'Classified as natural product in ChEMBL';

-- ─── ASSAYS ──────────────────────────────────────────────────
COMMENT ON TABLE pharmaco_db.assays IS 'Biological assays from ChEMBL. Human assays with confidence >= 4. Links compounds to targets via bioactivities.';

COMMENT ON COLUMN pharmaco_db.assays.id IS 'Internal primary key';
COMMENT ON COLUMN pharmaco_db.assays.chembl_id IS 'ChEMBL assay identifier (e.g. CHEMBL674637)';
COMMENT ON COLUMN pharmaco_db.assays.description IS 'Assay description (max 500 chars)';
COMMENT ON COLUMN pharmaco_db.assays.assay_type IS 'B=Binding, F=Functional, A=ADMET, T=Toxicity, P=Physicochemical';
COMMENT ON COLUMN pharmaco_db.assays.assay_category IS 'Assay category classification';
COMMENT ON COLUMN pharmaco_db.assays.src_id IS 'ChEMBL data source ID';
COMMENT ON COLUMN pharmaco_db.assays.src_description IS 'Data source description';
COMMENT ON COLUMN pharmaco_db.assays.journal IS 'Publication journal';
COMMENT ON COLUMN pharmaco_db.assays.year IS 'Publication year';
COMMENT ON COLUMN pharmaco_db.assays.doi IS 'Publication DOI';
COMMENT ON COLUMN pharmaco_db.assays.target_id IS 'FK to targets table';
COMMENT ON COLUMN pharmaco_db.assays.target_chembl_id IS 'ChEMBL target ID for reference';
COMMENT ON COLUMN pharmaco_db.assays.assay_organism IS 'Organism used in the assay';
COMMENT ON COLUMN pharmaco_db.assays.assay_cell_type IS 'Cell type used in the assay';
COMMENT ON COLUMN pharmaco_db.assays.confidence_score IS 'ChEMBL confidence score (0-9). 9 = direct single protein target';

-- ─── BIOACTIVITIES ───────────────────────────────────────────
COMMENT ON TABLE pharmaco_db.bioactivities IS 'Quantitative bioactivity measurements. Core SAR data linking compounds to targets with pChEMBL values. 4.8M rows with pChEMBL.';

COMMENT ON COLUMN pharmaco_db.bioactivities.id IS 'Internal primary key';
COMMENT ON COLUMN pharmaco_db.bioactivities.compound_id IS 'FK to compounds table';
COMMENT ON COLUMN pharmaco_db.bioactivities.target_id IS 'FK to targets table';
COMMENT ON COLUMN pharmaco_db.bioactivities.assay_id IS 'FK to assays table';
COMMENT ON COLUMN pharmaco_db.bioactivities.activity_type IS 'Measurement type: Ki, IC50, EC50, Kd, etc.';
COMMENT ON COLUMN pharmaco_db.bioactivities.relation IS 'Value relation: =, <, >, <=, >=, ~';
COMMENT ON COLUMN pharmaco_db.bioactivities.value IS 'Standardized activity value (typically nM)';
COMMENT ON COLUMN pharmaco_db.bioactivities.units IS 'Units of the value (nM, uM, etc.)';
COMMENT ON COLUMN pharmaco_db.bioactivities.pchembl_value IS 'Standardized -log10(molar IC50/Ki/Kd/EC50). Higher = more potent. >=6 is active, >=7 is potent';
COMMENT ON COLUMN pharmaco_db.bioactivities.activity_class IS 'Derived: active (pChEMBL>=7), intermediate (5-7), inactive (<5)';
COMMENT ON COLUMN pharmaco_db.bioactivities.source IS 'Data source (ChEMBL)';
COMMENT ON COLUMN pharmaco_db.bioactivities.source_id IS 'Original activity ID in source database';
COMMENT ON COLUMN pharmaco_db.bioactivities.data_validity IS 'ChEMBL data validity comment (flags potential issues)';
COMMENT ON COLUMN pharmaco_db.bioactivities.potential_duplicate IS 'Flagged as potential duplicate by ChEMBL';

-- ─── DRUG MECHANISMS ─────────────────────────────────────────
COMMENT ON TABLE pharmaco_db.drug_mechanisms IS 'Known drug mechanisms of action from ChEMBL. Links approved/clinical drugs to their target and MoA.';

COMMENT ON COLUMN pharmaco_db.drug_mechanisms.id IS 'Internal primary key';
COMMENT ON COLUMN pharmaco_db.drug_mechanisms.compound_id IS 'FK to compounds table';
COMMENT ON COLUMN pharmaco_db.drug_mechanisms.target_id IS 'FK to targets table';
COMMENT ON COLUMN pharmaco_db.drug_mechanisms.mechanism_of_action IS 'Free text MoA description (e.g. Epidermal growth factor receptor inhibitor)';
COMMENT ON COLUMN pharmaco_db.drug_mechanisms.action_type IS 'Standardized action: INHIBITOR, AGONIST, ANTAGONIST, MODULATOR, etc.';
COMMENT ON COLUMN pharmaco_db.drug_mechanisms.direct_interaction IS 'True if drug directly interacts with target protein';
COMMENT ON COLUMN pharmaco_db.drug_mechanisms.molecular_mechanism IS 'Molecular-level mechanism detail';
COMMENT ON COLUMN pharmaco_db.drug_mechanisms.selectivity_comment IS 'Selectivity notes (e.g. Selective, Non-selective)';

-- ─── DRUG INDICATIONS ────────────────────────────────────────
COMMENT ON TABLE pharmaco_db.drug_indications IS 'Therapeutic indications for drugs. MeSH and EFO disease mappings with clinical phase.';

COMMENT ON COLUMN pharmaco_db.drug_indications.id IS 'Internal primary key';
COMMENT ON COLUMN pharmaco_db.drug_indications.compound_id IS 'FK to compounds table';
COMMENT ON COLUMN pharmaco_db.drug_indications.mesh_id IS 'MeSH disease identifier';
COMMENT ON COLUMN pharmaco_db.drug_indications.mesh_heading IS 'MeSH disease name (e.g. Breast Neoplasms)';
COMMENT ON COLUMN pharmaco_db.drug_indications.efo_id IS 'Experimental Factor Ontology disease ID';
COMMENT ON COLUMN pharmaco_db.drug_indications.efo_term IS 'EFO disease term';
COMMENT ON COLUMN pharmaco_db.drug_indications.max_phase IS 'Maximum clinical phase for this indication (1-4)';
COMMENT ON COLUMN pharmaco_db.drug_indications.indication_refs IS 'Reference sources for this indication';

-- ─── CROSS REFERENCES ────────────────────────────────────────
COMMENT ON TABLE pharmaco_db.cross_references IS 'Cross-database links for targets and compounds. PDB, STRING, OMIM, IntAct, etc.';

COMMENT ON COLUMN pharmaco_db.cross_references.id IS 'Internal primary key';
COMMENT ON COLUMN pharmaco_db.cross_references.compound_id IS 'FK to compounds (nullable)';
COMMENT ON COLUMN pharmaco_db.cross_references.target_id IS 'FK to targets (nullable)';
COMMENT ON COLUMN pharmaco_db.cross_references.db_name IS 'Source database name (PDB, STRING, OMIM, IntAct, Orphanet, NCBI_Gene)';
COMMENT ON COLUMN pharmaco_db.cross_references.db_id IS 'Identifier in the source database';
COMMENT ON COLUMN pharmaco_db.cross_references.url IS 'Direct URL to the entry in the source database';

-- ─── INGESTION LOG ───────────────────────────────────────────
COMMENT ON TABLE pharmaco_db.ingestion_log IS 'Pipeline execution log tracking ingestion steps, status, and row counts.';
