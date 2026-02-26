# PharmacoDB Research Mission
## Computational Drug Discovery Pattern Mining & Novel Rule Extraction

> **Classification**: Internal Research Protocol
> **Target Publication**: Nature Methods / Nature Computational Science
> **PI**: BindX Research Division
> **Dataset**: PharmacoDB (~2M compounds, ~11K targets, bioactivities, mechanisms, indications, cross-refs)
> **Compute**: RunPod GPU cluster (unlimited), local HPC
> **Timeline**: Unlimited (iterative until exhaustion of all leads)

---

## 1. EXECUTIVE SUMMARY

This document defines an exhaustive, iterative research program to mine PharmacoDB for **validated known patterns** and **novel, unpublished rules** for in silico drug discovery. Every finding must be:

1. **Statistically rigorous** (corrected p-values, effect sizes, confidence intervals)
2. **Validated on blind holdout data** (20% never seen during exploration)
3. **Translatable** into scoring functions or filters for BindX pipeline
4. **Publication-ready** with figures, tables, and methodology sections

The program is organized into 12 Work Packages (WPs), each producing deliverables that feed into the next. Execution is **iterative**: results from early WPs reshape the scope of later ones. Dead ends are documented; unexpected signals are pursued to exhaustion.

---

## 2. DATA FOUNDATION

### 2.1 Available Data (Current State)

| Table | Rows | Key Fields |
|-------|------|------------|
| targets | 10,962 | UniProt function, GO terms, KEGG/Reactome pathways, PDB IDs, domains, InterPro/Pfam, disease associations, tissue expression, EC numbers, 56 total columns |
| compounds | ~2M | SMILES, InChI, MW, ALogP, HBA, HBD, PSA, RTB, Ro5 violations, aromatic rings, heavy atoms, QED, molecular formula |
| cross_references | 87,930 | PDB, STRING, NCBI Gene, OMIM, IntAct, Orphanet linkages |
| bioactivities | *pending* | compound_id, target_id, activity_type (Ki/IC50/EC50), pChEMBL, activity_class |
| assays | *pending* | assay type/category, organism, cell type, confidence score, journal/DOI |
| drug_mechanisms | *pending* | mechanism_of_action, action_type, direct_interaction |
| drug_indications | *pending* | MeSH/EFO disease terms, max_phase |

### 2.2 Enrichment Pipeline (In Progress)

- **Phase 2**: PubChem CIDs, ATC classification, GtoPdb interactions, compound/target tags, compound_structures (fingerprints, Murcko scaffolds)
- **Phase 3**: DGIdb, SIDER side effects, OpenFDA adverse events, ClinicalTrials.gov, STRING PPI, Reactome pathways, DisGeNET, KEGG

### 2.3 Blind Holdout Protocol

Before ANY analysis begins:

```
HOLDOUT SPLIT (by target, stratified by protein_family):
- Training/Exploration set: 80% of targets + all their compounds/bioactivities
- Blind Validation set: 20% of targets (held in sealed envelope)
  - Stratified so every protein_family is represented proportionally
  - NEVER used for hypothesis generation, only for final validation
  - Split recorded with random seed for reproducibility

TEMPORAL SPLIT (by compound registration date where available):
- Temporal training: compounds registered before 2020
- Temporal validation: compounds registered 2020+
- Tests whether discovered patterns generalize to novel chemistry
```

---

## 3. WORK PACKAGES

### WP1: Chemical Space Cartography & Dark Matter Analysis
**Priority: HIGH | Data needed: compounds (AVAILABLE NOW)**

#### 1.1 Objectives
- Map the complete chemical space of PharmacoDB using multiple dimensionality reduction methods
- Identify systematic voids (underexplored regions) and dense clusters
- Characterize "dark chemical matter" (compounds with no recorded activity despite being tested)
- Generate coverage maps per physicochemical property space

#### 1.2 Methods
- **UMAP/t-SNE** on Morgan FP (2048-bit), MACCS keys, and RDKit descriptors
- **HDBSCAN clustering** to identify natural compound groupings
- **Property space tiling**: MW x ALogP x PSA x QED 4D grid, count compounds per cell
- **Murcko scaffold decomposition**: frequency distribution, scaffold diversity index
- **Fragment analysis**: BRICS decomposition, frequency of fragments, fragment-property correlations
- **Chemical space density estimation**: KDE on fingerprint-based 2D projections
- **Natural product analysis**: Compare is_natural_product=true vs synthetic compound space coverage
- **Named compound enrichment**: Analyze the 37,464 named compounds -- do they cluster differently?

#### 1.3 Deliverables
- Chemical space atlas (interactive, exportable)
- List of void regions with suggested generative targets
- Scaffold frequency catalog with family annotations
- Fragment importance ranking
- **BindX Rule**: Chemical space density-based novelty score for new compounds

#### 1.4 Iterative Follow-up
- If voids correlate with specific property ranges -> investigate if those ranges are intrinsically difficult or simply unexplored
- If scaffold clusters align with target families (once bioactivities arrive) -> this feeds WP4

---

### WP2: Physicochemical Rule Mining Beyond Lipinski
**Priority: HIGH | Data needed: compounds + bioactivities**

#### 2.1 Objectives
- Validate Lipinski/Veber/Ghose/Egan rules against real bioactivity data at scale
- Discover NEW property rules that predict active vs. inactive compounds
- Mine target-family-specific property windows
- Identify emergent multi-property trade-offs (AMODO-EO style)

#### 2.2 Methods
- **Rule validation matrix**: For each known rule (Lipinski Ro5, Veber, Ghose, Egan, Muegge, lead-like, fragment-like, PPI-like, CNS MPO, BBB):
  - Sensitivity/specificity/PPV/NPV against pChEMBL >= 6 (active), < 5 (inactive)
  - Stratified by target family, by activity type (Ki vs IC50 vs EC50)
  - ROC-AUC of rule vs. activity class

- **Decision tree mining**: Train shallow decision trees (depth 3-5) on compound properties -> activity class
  - Extract interpretable rules
  - Compare discovered rules vs. known rules
  - Identify novel property combinations with high discriminative power

- **Emergent trade-off discovery**:
  - Compute all pairwise property ratios (MW/PSA, HBA/RTB, ALogP*aromatic_rings, etc.)
  - Mutual information between ratios and activity class
  - Surprise metric: ratios with high MI but not in established rules

- **Target-family windows**:
  - For each protein_family (kinase, GPCR, ion channel, protease, etc.):
    - Property distributions of actives vs. inactives
    - Optimal property windows (Bayesian optimization of thresholds)
    - Compare windows across families -> identify family-specific vs. universal rules

- **QED deconstruction**:
  - Which QED sub-components matter most for which target families?
  - Can we build a better "target-aware QED"?

#### 2.3 Deliverables
- Comprehensive rule validation table (publishable supplementary)
- Novel rules with effect sizes and validation on blind set
- Target-family-specific property windows
- Improved QED formula per target class
- **BindX Rule**: Target-aware compound scoring function

---

### WP3: Target Druggability Prediction from Protein Features
**Priority: HIGH | Data needed: targets (AVAILABLE NOW)**

#### 3.1 Objectives
- Build a druggability prediction model using the 56 target feature columns
- Identify which protein features most strongly predict successful drug discovery
- Discover novel feature combinations that predict "hidden druggable" targets
- Validate against known approved drug targets

#### 3.2 Methods
- **Feature engineering from existing columns**:
  - GO term embeddings (Word2Vec on GO term sets)
  - KEGG pathway graph features (centrality, connectivity)
  - Domain composition vectors (Pfam/InterPro one-hot)
  - Disease association count and diversity
  - PDB structure count (proxy for research interest and tractability)
  - Tissue specificity entropy (Shannon entropy of tissue expression)
  - Sequence features: length, mass, transmembrane count, PTM density

- **Druggability label construction**:
  - Positive: targets with >= 1 compound at pChEMBL >= 7 AND max_phase >= 3
  - Negative: targets with >= 10 tested compounds but none at pChEMBL >= 6
  - Unknown: insufficient data

- **Models**:
  - Random Forest + SHAP explanations
  - XGBoost with hyperparameter tuning
  - Logistic regression (for interpretability)
  - ESM-2 protein embeddings (via HuggingFace) as additional features

- **Novel angles**:
  - **Cross-species druggability conservation**: Is a target druggable in human iff druggable in mouse?
  - **Pathway context**: Targets in well-drugged pathways are more likely druggable?
  - **Binding site geometry** (from PDB): pocket volume, hydrophobicity, flexibility
  - **Disease association strength**: Is Open Targets genetic evidence predictive of druggability?

#### 3.3 Deliverables
- Druggability prediction model (AUC > 0.85 on blind set)
- Feature importance ranking (SHAP-based)
- List of "hidden druggable" targets (predicted druggable, not yet drugged)
- **BindX Rule**: Target prioritization score for campaign creation

---

### WP4: Privileged Scaffold Mining by Target Family
**Priority: CRITICAL | Data needed: compounds + bioactivities + targets**

#### 4.1 Objectives
- Systematically identify privileged scaffolds for each target family
- Quantify scaffold-family affinity enrichment
- Discover cross-family scaffolds (true privileged structures) vs. family-specific scaffolds
- Build scaffold recommendation engine

#### 4.2 Methods
- **Murcko scaffold extraction** for all compounds with pChEMBL >= 6
- **Enrichment analysis**: For each scaffold S and family F:
  - P(active | scaffold=S, family=F) / P(active | family=F) = Enrichment Factor
  - Fisher exact test with BH correction
  - Effect size (odds ratio + 95% CI)

- **Scaffold-target family heatmap**: Scaffolds (rows) x Families (columns), colored by enrichment
  - Cluster to find privileged groups (high enrichment across multiple families)
  - vs. specialist scaffolds (high enrichment in exactly one family)

- **Substructure mining** (more granular than Murcko):
  - BRICS fragments enrichment per family
  - Molecular matched pair analysis within families
  - RECAP decomposition
  - R-group enrichment at specific positions

- **Anti-privileged scaffolds**: Scaffolds that are systematically INACTIVE across all families
  - These are PAINS-like but data-driven rather than rule-based
  - Compare to known PAINS alerts -> identify novel problematic substructures

- **Temporal privileged scaffold evolution**: Do newly discovered scaffolds (post-2020) differ from historical ones?

#### 4.3 Deliverables
- Privileged scaffold catalog per target family (with enrichment, frequency, confidence)
- Cross-family privileged structure atlas
- Anti-privileged scaffold blacklist (data-driven PAINS replacement)
- Substructure-family affinity rules
- **BindX Rule**: Scaffold-based compound prioritization score per target family

---

### WP5: Activity Cliff Cartography & SAR Transfer Rules
**Priority: CRITICAL | Data needed: compounds + bioactivities**

#### 5.1 Objectives
- Map all activity cliffs in PharmacoDB (structurally similar pairs with large activity differences)
- Identify structural modifications responsible for cliffs
- Discover cross-target cliff conservation (same modification causes cliff on related targets)
- Extract "SAR transfer rules" -- modifications that reliably improve potency

#### 5.2 Methods
- **Pairwise similarity**: Tanimoto on Morgan FP for all compounds within each target
  - Threshold: Tanimoto >= 0.8, |delta pChEMBL| >= 2.0
  - Result: activity cliff compound pairs

- **SALI (Structure-Activity Landscape Index)**: |delta_activity| / (1 - similarity)
  - Top 1% SALI pairs = most extreme cliffs

- **Cliff decomposition**:
  - For each cliff pair: compute molecular matched pair (MMP) transformation
  - Catalog all transformations causing cliffs
  - Frequency x magnitude = "transformation impact score"

- **Cross-target cliff analysis** (NOVEL):
  - For each MMP transformation found in cliffs on target A:
    - Does the same transformation cause cliffs on targets B, C, D within the same family?
    - Does it cause cliffs across families?
  - "Universal cliff transformations" vs "target-specific cliff transformations"

- **SAR transfer rules** (NOVEL):
  - MMP transformations that CONSISTENTLY improve potency (not just create cliffs)
  - Regression: delta_pChEMBL ~ transformation + target_family + property_changes
  - Goal: "Transformation X improves potency by Y units on kinases with 90% reliability"

- **Cliff prediction model**:
  - Given two compounds, predict if they form a cliff
  - Features: structural diff, property diff, target family, scaffold
  - Can we predict cliffs before measuring activity?

#### 5.3 Deliverables
- Complete activity cliff atlas (interactive, filterable)
- MMP transformation catalog with impact scores
- Cross-target cliff conservation matrix
- SAR transfer rules per target family
- Cliff prediction model
- **BindX Rule**: SAR-transfer-aware compound scoring; penalty for cliff-prone modifications

---

### WP6: Polypharmacology & Selectivity Pattern Mining
**Priority: CRITICAL | Data needed: bioactivities + targets**

#### 6.1 Objectives
- Build compound polypharmacology profiles (activity across multiple targets)
- Identify selectivity determinants (what makes a compound selective?)
- Discover "selectivity switches" -- small modifications that toggle target selectivity
- Mine therapeutic polypharmacology opportunities

#### 6.2 Methods
- **Activity matrix construction**: Compounds (rows) x Targets (columns), values = pChEMBL
  - Sparsity analysis: what % of matrix is filled?
  - Sub-matrices for well-tested compound-target blocks

- **Selectivity metrics**:
  - Selectivity entropy: -sum(p_i * log(p_i)) over target activity profile
  - Gini coefficient of activity distribution
  - Max selectivity window: best_pChEMBL - second_best_pChEMBL
  - Target count at pChEMBL >= 6

- **Selectivity determinants** (NOVEL):
  - For each target pair (A, B) with shared compounds:
    - What structural features predict A-selective vs B-selective?
    - Decision tree on molecular descriptors -> selectivity class
    - Extract interpretable selectivity rules

- **Polypharmacology profiles**:
  - Cluster compounds by their activity vectors (correlation-based distance)
  - "Polypharmacology fingerprint" -- binary vector of targets hit at pChEMBL >= 6
  - Association rules mining: {hits target A, hits target B} -> {hits target C}

- **Selectivity switches** (NOVEL):
  - MMP analysis where one compound is selective, the other promiscuous
  - Catalog "selectivity-inducing" and "selectivity-breaking" transformations

- **Therapeutic polypharmacology**:
  - Cross-reference polypharmacology profiles with disease-target associations
  - Score: sum of disease relevance scores for all targets hit
  - Identify compounds with therapeutically synergistic multi-target profiles

#### 6.3 Deliverables
- Selectivity rule catalog per target pair
- Polypharmacology opportunity ranking
- Selectivity switch transformation library
- **BindX Rule**: Selectivity-aware scoring; polypharmacology opportunity flag

---

### WP7: Bioactivity Matrix Imputation & Virtual Screening Enhancement
**Priority: HIGH | Data needed: bioactivities**

#### 7.1 Objectives
- Impute missing values in the compound-target bioactivity matrix
- Discover "hidden interactions" -- untested compound-target pairs predicted as strong binders
- Build confidence-aware imputation for BindX virtual screening prioritization

#### 7.2 Methods
- **Matrix factorization**:
  - NMF (Non-negative Matrix Factorization) on pChEMBL matrix
  - Bayesian PMF (Probabilistic Matrix Factorization) with uncertainty
  - Deep matrix factorization (neural collaborative filtering)

- **Side information integration**:
  - Compound features (properties, fingerprints) -> compound embeddings
  - Target features (GO terms, domains, family) -> target embeddings
  - Graph neural network on compound-target bipartite graph
  - Knowledge graph embeddings (cross-references, pathways, diseases)

- **Confidence calibration**:
  - For each predicted interaction, compute:
    - Prediction variance (from Bayesian methods)
    - Neighborhood support (how many similar compounds/targets have measured data?)
    - Model agreement (ensemble disagreement)

- **Validation**:
  - Leave-one-out on known interactions
  - Leave-one-target-out (can we predict a target's profile from related targets?)
  - Temporal validation (predict 2020+ interactions from pre-2020 data)
  - Blind holdout validation

- **Active learning simulation**:
  - If we could measure N new datapoints, which ones would most improve the model?
  - Quantify information gain per measurement

#### 7.3 Deliverables
- Imputed bioactivity matrix with confidence scores
- Top 1000 highest-confidence novel interaction predictions
- Active learning prioritization for experimental validation
- **BindX Rule**: Imputation-based virtual screening pre-filter

---

### WP8: Drug Repurposing via Target-Disease-Compound Triangulation
**Priority: HIGH | Data needed: bioactivities + drug_mechanisms + drug_indications + disease associations**

#### 8.1 Objectives
- Build a knowledge graph from PharmacoDB
- Identify drug repurposing candidates by triangulating compound-target-disease relationships
- Score repurposing opportunities by strength of evidence
- Discover non-obvious repurposing paths (multi-hop reasoning)

#### 8.2 Methods
- **Knowledge graph construction**:
  - Nodes: compounds, targets, diseases (from indications + disease_associations), pathways, protein families
  - Edges: bioactivity (weighted by pChEMBL), mechanism, indication, disease-target association, PPI (STRING), pathway membership
  - Edge weights: confidence x relevance

- **Path-based reasoning**:
  - Drug -> Target -> Disease: direct evidence
  - Drug -> Target -> Pathway -> Target2 -> Disease: pathway-mediated
  - Drug -> Target -> PPI_partner -> Disease: network-mediated
  - Score each path by product of edge confidences

- **Embedding-based prediction**:
  - TransE/RotatE/ComplEx embeddings of the KG
  - Predict missing Drug->Disease edges
  - Compare with path-based scores

- **Consilience scoring**: Agreement across multiple methods (path-based, embedding-based, literature co-occurrence) -> strongest candidates

- **Negative validation**: Cross-reference with known failed clinical trials -> do our methods correctly avoid predicting failed drugs?

#### 8.3 Deliverables
- Knowledge graph (exportable, queryable)
- Ranked repurposing candidates with multi-evidence scores
- Novel disease-target-compound triangles
- **BindX Rule**: Repurposing opportunity score per compound-disease pair

---

### WP9: Mechanism-of-Action Prediction from Bioactivity Fingerprints
**Priority: MEDIUM | Data needed: bioactivities + drug_mechanisms**

#### 9.1 Objectives
- Predict mechanism of action from a compound's bioactivity profile
- Identify "MoA signatures" -- characteristic activity patterns per mechanism
- Discover compounds with unexpected MoA predictions (potential novel mechanisms)

#### 9.2 Methods
- **Bioactivity fingerprint construction**:
  - For each compound: vector of pChEMBL values across all tested targets
  - Sparse representation: only targets with measured data
  - Dense representation: imputed values from WP7

- **MoA classification**:
  - Labels: from drug_mechanisms.action_type (inhibitor, agonist, antagonist, etc.)
  - Models: kNN on bioactivity fingerprints, Random Forest, neural network
  - Feature importance: which targets' activities are most informative for each MoA?

- **MoA signature extraction**:
  - For each MoA class: characteristic activity pattern (mean + variance per target)
  - Discriminative targets: targets where activity differs most between MoA classes
  - "MoA panel": minimal set of targets to assay for MoA determination

- **Novel MoA detection**:
  - Compounds whose bioactivity profiles don't match any known MoA signature
  - Outlier detection in MoA space -> candidates for novel mechanisms
  - Cross-reference with structural novelty (from WP1)

#### 9.3 Deliverables
- MoA prediction model with confidence
- MoA signature catalog
- Minimal MoA assay panel recommendation
- Novel MoA candidate list
- **BindX Rule**: MoA-informed compound classification

---

### WP10: Safety & Toxicity Signal Mining
**Priority: HIGH | Data needed: bioactivities + side_effects + adverse_events + toxicology_data**

#### 10.1 Objectives
- Correlate structural features with safety signals
- Identify "toxicophores" -- substructures associated with adverse events
- Build predictive safety scoring from compound features
- Discover target-mediated safety risks (on-target toxicity patterns)

#### 10.2 Methods
- **Structure-toxicity correlation**:
  - Substructure mining: enrichment of fragments in compounds with reported adverse events
  - MMP analysis: transformations that introduce/remove safety liabilities
  - Property windows for safety: MW, ALogP, PSA ranges associated with fewer AEs

- **Target-mediated toxicity**:
  - For each target: proportion of active compounds with safety signals
  - "Safety risk targets" -- targets where activity correlates with adverse events
  - Pathway-level safety analysis: pathways enriched in safety-risk targets

- **Predictive safety model**:
  - Features: structural (fingerprints, descriptors), target profile (from WP6), properties
  - Labels: binary (any serious AE) and multi-label (specific AE categories)
  - Calibrated probabilities for risk assessment

- **Selectivity-safety relationship** (NOVEL):
  - Are selective compounds safer than promiscuous ones?
  - Is there a "safety selectivity window" -- minimum selectivity needed for acceptable safety?
  - Target-pair-specific: hitting target A + target B = safe, but target A + target C = toxic

#### 10.3 Deliverables
- Data-driven toxicophore catalog (beyond PAINS)
- Target-mediated safety risk map
- Predictive safety model
- Selectivity-safety rules
- **BindX Rule**: Safety scoring function for compound prioritization

---

### WP11: Multi-Objective Pareto Frontier Analysis
**Priority: HIGH | Data needed: compounds + bioactivities + ADMET data**

#### 11.1 Objectives
- Map Pareto frontiers for all relevant objective pairs/triples
- Identify target-family-specific trade-off shapes
- Discover "Pareto sweet spots" -- regions where multiple objectives can be simultaneously optimized
- Find emergent trade-offs not in established frameworks

#### 11.2 Methods
- **Pairwise Pareto frontiers**:
  - For each target family x {potency vs. MW, potency vs. ALogP, potency vs. PSA, potency vs. QED, potency vs. selectivity, potency vs. Ro5 violations}
  - Compute empirical Pareto front
  - Fit smooth frontier curve (kernel regression)
  - Quantify "room for improvement" -- distance of average compound from frontier

- **Multi-objective Pareto analysis (3D+)**:
  - Potency x QED x selectivity: 3D Pareto surface per family
  - NSGA-II style non-dominated sorting of real compounds
  - Pareto rank as compound quality metric

- **Emergent trade-offs** (AMODO-EO inspired):
  - Compute all pairwise property ratios and products
  - For each pair: correlation with bioactivity that is NOT explained by individual properties
  - Information gain of combined feature over individual features
  - Statistical test: interaction term significance in regression

- **Historical frontier movement**:
  - Compounds binned by year of publication (from assay metadata)
  - Pareto frontier per decade -> has the frontier advanced? Where? Where not?
  - Identifies "stalled" optimization campaigns

#### 11.3 Deliverables
- Pareto frontier atlas per target family
- Emergent trade-off catalog
- Pareto-based compound quality scores
- Historical frontier evolution analysis
- **BindX Rule**: Multi-objective Pareto ranking function

---

### WP12: Integrated Meta-Analysis & BindX Rule Synthesis
**Priority: FINAL | Data needed: all WP results**

#### 12.1 Objectives
- Integrate all discovered rules into a unified scoring framework
- Validate the integrated framework on blind holdout data
- Quantify the improvement over standard methods (Lipinski alone, single-objective scoring, etc.)
- Package as BindX algorithm module

#### 12.2 Methods
- **Rule integration**:
  - Each WP produces scoring rules; combine into a unified multi-criteria decision function
  - Weight optimization: learn optimal weights on training data, validate on blind set
  - Ablation study: contribution of each WP's rules to final performance

- **Benchmark against baselines**:
  - Baseline 1: Lipinski Ro5 + PAINS filter
  - Baseline 2: QED-based ranking
  - Baseline 3: Single-objective pChEMBL prediction
  - Our framework vs. baselines on multiple metrics (enrichment factor, BEDROC, AUC)

- **Sensitivity analysis**:
  - Which rules are most impactful? (SHAP on final model)
  - Which rules are most robust? (stability across bootstrap samples)
  - Which rules are most novel? (not captured by existing scoring functions)

- **Edge case analysis**:
  - Where does the framework fail? What types of compounds/targets?
  - Are failures systematic or random?
  - Suggest future data needs to address gaps

#### 12.3 Deliverables
- Unified BindX scoring framework with learned weights
- Comprehensive benchmark results
- Publication manuscript (Nature Methods format)
- Code package for BindX integration
- **BindX Algorithm**: Complete scoring pipeline deployable in production

---

## 4. EXECUTION PROTOCOL

### 4.1 Iterative Methodology

```
FOR each WP:
  1. EXPLORE: Run initial analysis, examine distributions, identify signals
  2. HYPOTHESIZE: Form specific testable hypotheses from observed patterns
  3. TEST: Statistical tests on exploration set with proper corrections
  4. PURSUE: For confirmed signals, dig deeper:
     - What mechanism explains this pattern?
     - Does it hold across subgroups (families, organisms, assay types)?
     - What are the edge cases and exceptions?
     - Can we find this signal in external databases?
  5. VALIDATE: Test final rules on blind holdout
  6. DOCUMENT: Full methodology + results + figures
  7. ITERATE: New questions arising from results -> back to EXPLORE
```

### 4.2 Statistical Rigor

- All hypothesis tests: Benjamini-Hochberg FDR correction at alpha=0.05
- Effect sizes mandatory (Cohen's d, odds ratio, enrichment factor)
- Confidence intervals on all point estimates (95%)
- Bootstrapped standard errors for complex statistics
- Multiple comparison correction within and across WPs
- Pre-registration of hypotheses before blind set validation

### 4.3 Reproducibility

- All analyses as Python scripts or Jupyter notebooks
- Random seeds fixed and recorded
- Complete dependency specification (requirements.txt)
- Data versioning with checksums
- Results cached with provenance tracking

### 4.4 Computational Stack

| Tool | Use |
|------|-----|
| **RDKit** | Molecular descriptors, fingerprints, MMP, scaffolds, substructure search |
| **scikit-learn** | ML models, clustering, dimensionality reduction |
| **XGBoost/LightGBM** | Gradient boosting models |
| **PyTorch** | Deep learning (graph neural networks, matrix factorization) |
| **PyTorch Geometric** | GNNs for molecular graphs and knowledge graphs |
| **UMAP/openTSNE** | Dimensionality reduction for visualization |
| **HDBSCAN** | Density-based clustering |
| **statsmodels** | Statistical tests, regression |
| **networkx/igraph** | Graph analysis for KG and PPI networks |
| **mols2grid** | Interactive molecule visualization |
| **Plotly/Seaborn/Matplotlib** | Publication-quality figures |
| **DuckDB/Polars** | Fast analytical queries on large datasets |
| **ESM-2** (HuggingFace) | Protein language model embeddings |
| **RunPod GPU** | Heavy compute (GNN training, large MF, ESM-2 inference) |

---

## 5. SUCCESS CRITERIA

### Publication-Ready Results
- Minimum 3 novel, statistically validated findings not in existing literature
- At least 1 finding that challenges or refines an established rule
- Blind holdout validation AUC improvement >= 5% over standard methods
- All figures at Nature-quality resolution and clarity

### BindX Integration
- Scoring function implementable in < 100ms per compound
- Configurable per target family
- Graceful degradation when data is missing (some features unavailable)
- API-compatible with existing BindX pipeline

### Data Coverage
- Analyses must cover >= 80% of protein families with >= 50 compounds
- Minimum 3 independent validation approaches per major finding
- Edge case documentation for all rules

---

## 6. RISK MANAGEMENT

| Risk | Mitigation |
|------|-----------|
| Bioactivities not yet loaded | Start with WP1 (compounds only) and WP3 (targets only) immediately. Queue WP2/4/5/6 for when bioactivities arrive. |
| Sparse bioactivity matrix (< 1% filled) | Use collaborative filtering (WP7), restrict analyses to well-tested compound-target blocks |
| Confounding by assay type | Always stratify by activity_type (Ki vs IC50 vs EC50); use pChEMBL for normalization |
| Multiple comparison inflation | Pre-register primary analyses; BH-FDR correction; replicate in holdout |
| Computational cost | GPU-accelerated fingerprint computation; DuckDB for fast analytics; incremental computation |
| Overfitting to ChEMBL biases | Temporal validation; diversity of validation sets; cross-database checks where possible |

---

## 7. TIMELINE (ADAPTIVE)

| Phase | WPs | Trigger |
|-------|-----|---------|
| **Immediate** | WP1 (chemical space), WP3 (druggability) | Compounds + targets available NOW |
| **After bioactivities load** | WP2 (property rules), WP4 (scaffolds), WP5 (cliffs), WP6 (selectivity) | bioactivities table > 0 rows |
| **After Phase 2/3 enrichment** | WP8 (repurposing), WP9 (MoA), WP10 (safety) | drug_mechanisms + indications + side_effects populated |
| **After imputation** | WP7 (matrix imputation) | WP6 completed (need activity matrix) |
| **After all WPs** | WP11 (Pareto), WP12 (integration) | All previous WPs produce preliminary results |
| **Iterative loops** | Any WP | Unexpected findings trigger re-exploration |

---

## 8. APPENDIX: LITERATURE REFERENCES

### Foundational
1. Lipinski CA et al. "Rule of Five" -- Adv Drug Deliv Rev, 2001
2. Leeson PD, Springthorpe B. "Influence of drug-like properties" -- Nature Rev Drug Discov, 2007
3. Bajorath J. "Activity landscape representations" -- J Med Chem, 2017

### Cutting Edge (2024-2026)
4. Activity cliffs: ACARL framework, J Cheminformatics 2025; MTPNet 2025
5. Privileged scaffolds: ChemBounce, Bioinformatics 2025; ACS target class repurposing, 2025
6. Polypharmacology: PCM rigorous evaluation, JCIM 2025; polypharmacology new drugs, Pharmacol Reports 2025
7. Matrix imputation: QComp, JCIM 2025; hyperbolic MF, Scientific Reports 2023
8. Drug repurposing: K-Paths 2025; OREGANO KG, 2024
9. Chemical space: BioReCS, Frontiers 2025; ML-guided docking, Nature Comp Sci 2025
10. Pareto: PMMG, Advanced Science 2025; AMODO-EO, ChemRxiv 2025
11. Druggability: DrugProtAI, Briefings Bioinf 2025; DrugTar, bioRxiv 2024
12. MoA prediction: Cell Painting, Nature Comms 2024; XGDP, Scientific Reports 2025
13. Foundation models: Boltz-2, 2025; MolE, 2025; VideoMol, Nature Comms 2024

---

*This document is a living specification. It will be updated as results from each WP reshape priorities and reveal new directions.*
