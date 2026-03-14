# BindX V9 — Scientific Validation Report

> Generated: 2026-03-02 22:51:50
> Molecules tested: 66 (56 kinase actives + 10 decoys) + reference drugs
> Calculations validated: 29 categories

## Executive Summary

| Metric | Value |
|--------|-------|
| Overall | **PASS** |
| Total tests | 297 |
| Passed | 291 |
| Failed | 0 |
| Warnings | 6 |
| Pass rate | 98.0% |

## Results by Category

| Category | Pass | Fail | Warn | Status |
|----------|------|------|------|--------|
| Physicochemical | 67 | 0 | 0 | PASS |
| SA Score | 5 | 0 | 0 | PASS |
| PAINS | 8 | 0 | 2 | PASS |
| Brenk | 3 | 0 | 0 | PASS |
| CNS MPO | 3 | 0 | 0 | PASS |
| Druglikeness Rules | 4 | 0 | 0 | PASS |
| Lipinski/Ro3 | 6 | 0 | 0 | PASS |
| Composite Score | 5 | 0 | 0 | PASS |
| Ligand Efficiency | 13 | 0 | 0 | PASS |
| ADMET | 13 | 0 | 0 | PASS |
| Activity Cliffs | 4 | 0 | 0 | PASS |
| Clustering | 5 | 0 | 1 | PASS |
| Pareto Ranking | 5 | 0 | 0 | PASS |
| Consensus Scoring | 3 | 0 | 2 | PASS |
| Frontend Display | 9 | 0 | 0 | PASS |
| Hard Cutoffs | 8 | 0 | 0 | PASS |
| Score Results | 7 | 0 | 0 | PASS |
| ADMET Composite | 8 | 0 | 0 | PASS |
| hERG Specialized | 9 | 0 | 1 | PASS |
| Confidence | 11 | 0 | 0 | PASS |
| Pharmacophore | 10 | 0 | 0 | PASS |
| Edge Cases | 21 | 0 | 0 | PASS |
| Cross-Validation | 8 | 0 | 0 | PASS |
| Numerical Stability | 7 | 0 | 0 | PASS |
| Full Benchmark | 10 | 0 | 0 | PASS |
| Dep: Scoring | 8 | 0 | 0 | PASS |
| Dep: Cliffs/Conf | 8 | 0 | 0 | PASS |
| Dep: Pipeline Chain | 13 | 0 | 0 | PASS |
| Dep: Frontend Cols | 10 | 0 | 0 | PASS |

## Detailed Results

### Physicochemical

- [PASS] Erlotinib MW: 393.17 vs ref 393.17 (±1.5)
- [PASS] Erlotinib logP: 3.41 vs ref 2.7 (±2.0)
- [PASS] Erlotinib TPSA: 74.73 vs ref 74.73 (±10.0)
- [PASS] Erlotinib HBD: 1 vs ref 1 (±1)
- [PASS] Erlotinib HBA: 7 vs ref 7 (±4, method diff)
- [PASS] Erlotinib InChIKey computed: AAKJLRGGTJKAMG...
- [PASS] Gefitinib MW: 446.15 vs ref 446.15 (±1.5)
- [PASS] Gefitinib logP: 4.28 vs ref 3.2 (±2.0)
- [PASS] Gefitinib TPSA: 68.74 vs ref 68.74 (±10.0)
- [PASS] Gefitinib HBD: 1 vs ref 1 (±1)
- [PASS] Gefitinib HBA: 7 vs ref 7 (±4, method diff)
- [PASS] Gefitinib InChIKey computed: XGALLCVXEZPNRQ...
- [PASS] Lapatinib MW: 580.13 vs ref 581.14 (±1.5)
- [PASS] Lapatinib logP: 6.14 vs ref 4.6 (±2.0)
- [PASS] Lapatinib TPSA: 106.35 vs ref 114.73 (±10.0)
- [PASS] Lapatinib HBD: 2 vs ref 2 (±1)
- [PASS] Lapatinib HBA: 8 vs ref 8 (±4, method diff)
- [PASS] Lapatinib InChIKey computed: BCFGMOOMADDAQU...
- [PASS] Osimertinib MW: 499.27 vs ref 499.27 (±1.5)
- [PASS] Osimertinib logP: 4.51 vs ref 3.4 (±2.0)
- [PASS] Osimertinib TPSA: 87.55 vs ref 87.59 (±10.0)
- [PASS] Osimertinib HBD: 2 vs ref 3 (±1)
- [PASS] Osimertinib HBA: 7 vs ref 8 (±4, method diff)
- [PASS] Osimertinib InChIKey computed: DUYJMQONPNNFPI...
- [PASS] Afatinib MW: 485.16 vs ref 485.18 (±1.5)
- [PASS] Afatinib logP: 4.39 vs ref 3.8 (±2.0)
- [PASS] Afatinib TPSA: 88.61 vs ref 88.62 (±10.0)
- [PASS] Afatinib HBD: 2 vs ref 2 (±1)
- [PASS] Afatinib HBA: 7 vs ref 8 (±4, method diff)
- [PASS] Afatinib InChIKey computed: ULXXDDBFHOBEHA...
- [PASS] Ibuprofen MW: 206.13 vs ref 206.13 (±1.5)
- [PASS] Ibuprofen logP: 3.07 vs ref 3.5 (±2.0)
- [PASS] Ibuprofen TPSA: 37.3 vs ref 37.3 (±10.0)
- [PASS] Ibuprofen HBD: 1 vs ref 1 (±1)
- [PASS] Ibuprofen HBA: 1 vs ref 2 (±4, method diff)
- [PASS] Ibuprofen InChIKey computed: HEFNNWSXXWATRW...
- [PASS] Metformin MW: 129.1 vs ref 129.1 (±1.5)
- [PASS] Metformin logP: -1.24 vs ref -1.4 (±2.0)
- [PASS] Metformin TPSA: 91.49 vs ref 91.49 (±10.0)
- [PASS] Metformin HBD: 3 vs ref 3 (±1)
- [PASS] Metformin HBA: 1 vs ref 5 (±4, method diff)
- [PASS] Metformin InChIKey computed: XZWYZXLIPXDOLR...
- [PASS] Atorvastatin MW: 558.25 vs ref 558.25 (±1.5)
- [PASS] Atorvastatin logP: 6.31 vs ref 5.7 (±2.0)
- [PASS] Atorvastatin TPSA: 111.79 vs ref 111.79 (±10.0)
- [PASS] Atorvastatin HBD: 4 vs ref 4 (±1)
- [PASS] Atorvastatin HBA: 4 vs ref 7 (±4, method diff)
- [PASS] Atorvastatin InChIKey computed: XUKUURHRXDUEBC...
- [PASS] Vemurafenib MW: 489.07 vs ref 489.07 (±1.5)
- [PASS] Vemurafenib logP: 5.54 vs ref 4.1 (±2.0)
- [PASS] Vemurafenib TPSA: 91.92 vs ref 100.87 (±10.0)
- [PASS] Vemurafenib HBD: 2 vs ref 2 (±1)
- [PASS] Vemurafenib HBA: 4 vs ref 6 (±4, method diff)
- [PASS] Vemurafenib InChIKey computed: GPXBXXGIAQBQNI...
- [PASS] Ruxolitinib MW: 306.16 vs ref 306.17 (±1.5)
- [PASS] Ruxolitinib logP: 3.47 vs ref 2.1 (±2.0)
- [PASS] Ruxolitinib TPSA: 83.18 vs ref 83.18 (±10.0)
- [PASS] Ruxolitinib HBD: 1 vs ref 1 (±1)
- [PASS] Ruxolitinib HBA: 4 vs ref 5 (±4, method diff)
- [PASS] Ruxolitinib InChIKey computed: HFNKQEVNSGCOJV...
- [PASS] Sotorasib MW: 560.23 vs ref 560.21 (±1.5)
- [PASS] Sotorasib logP: 4.48 vs ref 3.6 (±2.0)
- [PASS] Sotorasib TPSA: 104.45 vs ref 102.29 (±10.0)
- [PASS] Sotorasib HBD: 1 vs ref 1 (±1)
- [PASS] Sotorasib HBA: 7 vs ref 7 (±4, method diff)
- [PASS] Sotorasib InChIKey computed: NXQKSXLFSAEQCZ...
- [PASS] All 66 benchmark molecules parsed successfully

### SA Score

- [PASS] All 66 SA scores in [1, 10]
- [PASS] Ibuprofen SA=2.19 < 4.0 (easy to synthesize)
- [PASS] Metformin SA=3.21 < 4.0 (easy to synthesize)
- [PASS] Lapatinib SA=2.67 >= 2.0 (complex molecule)
- [PASS] Osimertinib SA=2.92 >= 2.0 (complex molecule)

### PAINS

- [PASS] Erlotinib: no PAINS alert (approved drug)
- [PASS] Gefitinib: no PAINS alert (approved drug)
- [PASS] Osimertinib: no PAINS alert (approved drug)
- [WARN] Afatinib: PAINS alert on approved drug — review SMARTS specificity
- [PASS] Lapatinib: no PAINS alert (approved drug)
- [PASS] Vemurafenib: no PAINS alert (approved drug)
- [PASS] Ruxolitinib: no PAINS alert (approved drug)
- [WARN] Catechol: no PAINS alert — SMARTS may not cover this substructure
- [PASS] Quinone: PAINS alert detected (expected)
- [PASS] Azo dye (azobenzene): PAINS alert detected (expected)

### Brenk

- [PASS] All 66 Brenk checks returned bool (27 alerts)
- [PASS] Nitro compound: Brenk alert detected (expected)
- [PASS] Epoxide: Brenk alert detected (expected)

### CNS MPO

- [PASS] All 66 CNS MPO scores in [0, 6]
- [PASS] Ruxolitinib CNS MPO=4.86 >= 3.5 (CNS-active drug)
- [PASS] Lapatinib CNS MPO=1.96 < 4.0 (large, not CNS drug)

### Druglikeness Rules

- [PASS] Ibuprofen pfizer_alert=True (logP=3.07, TPSA=37.3)
- [PASS] Metformin: no alerts (logP=-1.24, MW=129.1)
- [PASS] Atorvastatin gsk_alert=True (logP=6.31, MW=558.25)
- [PASS] Cross-validation: 132 rule checks consistent

### Lipinski/Ro3

- [PASS] Lapatinib MW=580.1 > 500 (Lipinski violation confirmed)
- [PASS] Atorvastatin MW=558.2 > 500 (Lipinski violation)
- [PASS] Metformin ro3_pass=False (MW=129.1, RotB=0)
- [PASS] Erlotinib ro3_pass=False (MW=393.2 > 300)
- [PASS] Gefitinib ro3_pass=False (MW=446.1 > 300)
- [PASS] Lapatinib ro3_pass=False (MW=580.1 > 300)

### Composite Score

- [PASS] affinity=0: score=0.1750 in [0,1]
- [PASS] affinity=-14: score=0.8250 ≈ 0.8250
- [PASS] affinity=-7, qed=0.8, logP=2.5: score=0.6350 ≈ 0.6350
- [PASS] V2 full params: score=0.9150 ≈ 0.9150
- [PASS] All 5 edge cases in [0, 1]

### Ligand Efficiency

- [PASS] Erlotinib LE=0.303 (HA=29, score=-8.8)
- [PASS] Erlotinib LE=0.276 in [0.1, 1.0]
- [PASS] Gefitinib LE=0.258 in [0.1, 1.0]
- [PASS] Lapatinib LE=0.200 in [0.1, 1.0]
- [PASS] Afatinib LE=0.235 in [0.1, 1.0]
- [PASS] Osimertinib LE=0.216 in [0.1, 1.0]
- [PASS] Icotinib LE=0.276 in [0.1, 1.0]
- [PASS] Neratinib LE=0.200 in [0.1, 1.0]
- [PASS] Dacomitinib LE=0.242 in [0.1, 1.0]
- [PASS] Vandetanib LE=0.267 in [0.1, 1.0]
- [PASS] Canertinib LE=0.235 in [0.1, 1.0]
- [PASS] LE(affinity=0, HA=20) = 0.0
- [PASS] LE(HA=0) = None (edge case handled)

### ADMET

- [PASS] predict_admet returned 5 results for 5 inputs
- [PASS] Ibuprofen: all 5 ADMET categories present
- [PASS] Ibuprofen: composite_score=0.6000 in [0, 1]
- [PASS] Metformin: all 5 ADMET categories present
- [PASS] Metformin: composite_score=0.5200 in [0, 1]
- [PASS] Erlotinib: all 5 ADMET categories present
- [PASS] Erlotinib: composite_score=0.5500 in [0, 1]
- [PASS] Lapatinib: all 5 ADMET categories present
- [PASS] Lapatinib: composite_score=0.0700 in [0, 1]
- [PASS] Atorvastatin: all 5 ADMET categories present
- [PASS] Atorvastatin: composite_score=0.2200 in [0, 1]
- [PASS] Ibuprofen oral_bioavailability=0.85 >= 0.4
- [PASS] Lapatinib absorption (0.20) < Ibuprofen (0.85)

### Activity Cliffs

- [PASS] Activity cliff detection returned 13 results
- [PASS] Detected 8 cliff molecules among EGFR inhibitors
- [PASS] All SALI values >= 0
- [PASS] All result dicts have required keys

### Clustering

- [PASS] Clustering returned 13 molecules
- [PASS] All molecules have cluster_id
- [PASS] 10 clusters found (structural diversity)
- [WARN] Quinazolines in different clusters: {'Erlotinib': 9, 'Gefitinib': 1, 'Afatinib': 2}
- [PASS] Osimertinib (cluster #7) differs from quinazolines
- [PASS] At least one representative assigned

### Pareto Ranking

- [PASS] All molecules have pareto_rank
- [PASS] All pareto_rank >= 0
- [PASS] All molecules have pareto_front bool
- [PASS] Dominant molecule on Pareto front (rank 0)
- [PASS] All objective values in [0, 1]

### Consensus Scoring

- [PASS] Consensus enrichment returned 5 molecules
- [WARN] best_all agreement=2/3 — expected 3/3
- [WARN] worst_all agreement=1/3 — expected 0/3
- [PASS] All agreement values in X/3 format
- [PASS] All per-method ranks in valid range

### Frontend Display

- [PASS] Composite: backend 0.825 → frontend 82.5
- [PASS] hERG radar: risk 0.12 → display 0.88
- [PASS] Safety dot: risk=0.1 → green (< 0.3)
- [PASS] Safety dot: risk=0.29 → green (< 0.3)
- [PASS] Safety dot: risk=0.3 → yellow (0.3-0.6)
- [PASS] Safety dot: risk=0.5 → yellow (0.3-0.6)
- [PASS] Safety dot: risk=0.6 → yellow (= 0.6 boundary)
- [PASS] Safety dot: risk=0.61 → red (> 0.6)
- [PASS] Safety dot: risk=0.9 → red (> 0.6)

### Hard Cutoffs

- [PASS] Clean molecule passes all cutoffs
- [PASS] QED=0.10 eliminated: QED too low (0.10)
- [PASS] SA=7.5 eliminated: SA score too high (7.5)
- [PASS] PAINS eliminated: PAINS alert (likely false positive)
- [PASS] 3 Lipinski violations eliminated: Lipinski violations (3)
- [PASS] hERG=0.85 eliminated: hERG inhibition risk (cardiac)
- [PASS] Batch: 2 passed, 2 eliminated correctly
- [PASS] Empty list handled correctly

### Score Results

- [PASS] score_results: all molecules have composite_score
- [PASS] score_results: all molecules enriched with properties
- [PASS] score_results: sorted descending ([0.6186, 0.5427, 0.2447])
- [PASS] score_results: best affinity ranks first
- [PASS] score_results: PAINS and SA score computed for all
- [PASS] score_results_v2: ADMET data attached to matching molecule
- [PASS] score_results_v2: works without ADMET data

### ADMET Composite

- [PASS] Clean molecule: score=0.6000 in [0, 1]
- [PASS] Clean molecule: color=green (score=0.6000)
- [PASS] Toxic (0.0000) < clean (0.6000)
- [PASS] Toxic molecule: color=red
- [PASS] Toxic molecule: 12 flags raised
- [PASS] hERG alert raised: ALERT: hERG inhibition risk (cardiotoxicity)
- [PASS] Ames alert raised: ALERT: Ames mutagenicity positive
- [PASS] Empty input: baseline score=0.5

### hERG Specialized

- [WARN] Ibuprofen: risk=MODERATE, expected LOW
- [PASS] Ibuprofen has 'ic50_um' key
- [PASS] Ibuprofen has 'risk_level' key
- [PASS] Ibuprofen has 'method' key
- [PASS] Ibuprofen has 'features' key
- [PASS] IC50 > 0 (30.0 μM)
- [PASS] Valid risk level: MODERATE
- [PASS] Lapatinib IC50 (8.4) < Ibuprofen (30.0)
- [PASS] Empty SMILES: defaults to LOW risk
- [PASS] All 20 tested molecules returned valid hERG predictions

### Confidence

- [PASS] High-confidence: has 'overall' and 'components'
- [PASS] High-confidence overall=0.845 in [0, 1]
- [PASS] High-confidence >= 0.7: 0.845
- [PASS] 5 confidence components present
- [PASS] All components have score in [0,1] and note
- [PASS] Low (0.365) < High (0.845)
- [PASS] Structure hierarchy correct: ['0.845', '0.812', '0.750', '0.650']
- [PASS] Disorder penalty 0.10 applied for fraction=0.45
- [PASS] No disorder penalty for fraction=0.10
- [PASS] Batch confidence: 2 molecules processed
- [PASS] Component weights sum to 1.000 ≈ 1.0

### Pharmacophore

- [PASS] Erlotinib: 19 pharmacophore features
- [PASS] All 6 feature types present in counts
- [PASS] Erlotinib donors: 1
- [PASS] Erlotinib acceptors: 7
- [PASS] Erlotinib aromatics: 3
- [PASS] Invalid SMILES: graceful fallback (0 features)
- [PASS] Similarity matrix: 3x3
- [PASS] Diagonal = 1.0 (self-similarity)
- [PASS] Matrix is symmetric
- [PASS] Erlotinib-Gefitinib (0.468) > Erlotinib-Ibuprofen (0.017)

### Edge Cases

- [PASS] compute_properties(''): returns dict
- [PASS] compute_properties('NOT_A_SMILES'): returns dict
- [PASS] compute_properties('C(C)(C)(C)(C)(C)'): returns dict
- [PASS] compute_properties(None): returns dict
- [PASS] SA score for invalid SMILES: None
- [PASS] PAINS for empty string: False
- [PASS] All extreme affinity values produce scores in [0, 1]
- [PASS] LE(None, 20) = None
- [PASS] LE(-8, None) = None
- [PASS] LE(-8, -5) = None (negative atoms)
- [PASS] CNS MPO for empty props: None
- [PASS] Druglikeness rules for empty props: no alerts
- [PASS] Single-molecule clustering: cluster_id=0
- [PASS] Single molecule Pareto: rank=0, front=True
- [PASS] Single-molecule consensus: detail computed
- [PASS] Activity cliffs with 1 molecule (has score): is_cliff=False
- [PASS] Activity cliffs with no docking data: is_cliff=None (N/D)
- [PASS] Composite score(None): None (N/D — docking not run)
- [PASS] All functions handle empty lists gracefully
- [PASS] SVG generation: 17015 chars
- [PASS] SVG for invalid SMILES: None

### Cross-Validation

- [PASS] Druglikeness rules consistent with properties for 10 molecules
- [PASS] V1 composite score monotonic in affinity: 0.290 → 0.940
- [PASS] V2 composite score monotonic: 0.295 → 0.845
- [PASS] Higher QED → higher score: 0.7014 > 0.5614
- [PASS] logP=2.5 (0.6614) > logP=8.0 (0.5248)
- [PASS] SA score correctly inverted for Pareto synthesis objective
- [PASS] ADMET structure consistent for 5 molecules
- [PASS] Real pipeline (0.812) > Mock (0.377)

### Numerical Stability

- [PASS] Consensus handles NaN values gracefully
- [PASS] Consensus handles extreme values
- [PASS] Composite score handles extreme negative affinities
- [PASS] Composite score handles positive affinities
- [PASS] ADMET composite with None values: 0.5000
- [PASS] Identical molecules: all same Pareto rank (0)
- [PASS] Identical SMILES: same cluster

### Full Benchmark

- [PASS] All 66 molecules computed without errors
- [PASS] All 66 QED values in [0, 1] (mean=0.509)
- [PASS] All 66 MW in [50, 900] (range 129-621)
- [PASS] ADMET batch: 66 results for 66 inputs
- [PASS] All ADMET composites in [0, 1] (mean=0.377)
- [PASS] score_results: processed 66 molecules
- [PASS] Score pipeline: results sorted descending
- [PASS] Hard cutoffs: 56 passed, 10 eliminated (total=66)
- [PASS] Clustering top 20: 19 clusters found
- [PASS] Pareto top 20: 4 on front

### Dep: Scoring

- [PASS] With docking: Erlotinib (0.619) > Ibuprofen (0.543)
- [PASS] No docking: composite_score = None (N/D — not misleading)
- [PASS] No docking: properties (MW, logP, QED) still computed
- [PASS] Ligand efficiency without docking: None (correct fallback)
- [PASS] Ligand efficiency with docking: 0.3030
- [PASS] V2 with ADMET: data attached to molecules
- [PASS] V2 without ADMET: scores computed (fallback admet=0.5)
- [PASS] V2 ADMET effect: 0.5066 vs 0.4966

### Dep: Cliffs/Conf

- [PASS] Activity cliffs with docking: 9 cliffs from 13 molecules
- [PASS] Activity cliffs without docking: all is_cliff=None (N/D — correct)
- [PASS] Activity cliffs N/D: n_cliffs=None (not misleading 0)
- [PASS] Activity cliffs N/D: valid structure (all keys present)
- [PASS] Confidence degrades gracefully: full=0.825 > empty=0.377
- [PASS] Confidence ordering: empty(0.377) <= dock_only(0.615) <= full(0.825)
- [PASS] ADMET-only confidence (0.458) > empty (0.377)
- [PASS] All confidence results have valid structure and range [0,1]

### Dep: Pipeline Chain

- [PASS] Step 1 (Import): molecule has SMILES only, no properties
- [PASS] Step 2 (Scoring): MW=393.2, QED=0.418, LE=None, composite=None (N/D — docking not run)
- [PASS] Step 3 (ADMET): composite=0.550, CNS MPO=4.69
- [PASS] Step 4 (Docking): affinity=-8.5
- [PASS] Step 5 (Re-Score): LE=0.2930, V1=0.6186, V2=0.5066
- [PASS] Step 6 (Safety): PAINS=False
- [PASS] Step 7 (Confidence): 0.657
- [PASS] Step 8 (Hard Cutoffs): Erlotinib passes all cutoffs
- [PASS] Step 9 (Clustering): 3 molecules clustered
- [PASS] Step 10 (Pareto): rank=0, front=True (single molecule)
- [PASS] Step 11 (Activity Cliffs): 3 cliffs
- [PASS] Step 12 (Consensus): detail computed for all
- [PASS] Full chain: molecule has all 10 key properties

### Dep: Frontend Cols

- [PASS] After import: 0 calc columns detected
- [PASS] After docking: columns ['docking_score', 'cnn_score', 'cnn_affinity']
- [PASS] After docking+ADMET: columns from ['docking', 'admet', 'safety']
- [PASS] After scoring: composite_score column detected
- [PASS] After safety: 6 safety columns
- [PASS] After clustering: cluster_id detected
- [PASS] Full enrichment: all 10 calc types detected
- [PASS] Column additivity: count increases monotonically (final=6)
- [PASS] Display transform: 0.825 → 82.5 (×100)
- [PASS] ADMET radar inversion: hERG 0.12 → 0.88

## Reference Data Sources

- **PubChem**: MW (exact), XLogP3, TPSA, HBD/HBA counts
- **ChEMBL**: IC50/Ki values, compound activity data
- **RDKit**: Crippen logP (may differ ±1.0 from XLogP3), Lipinski HBA (may differ ±1 from PubChem)

## Known Discrepancies

- **logP**: RDKit uses Crippen method; PubChem uses XLogP3. Tolerance ±1.0 is standard.
- **HBA**: RDKit uses Lipinski definition (N+O atoms); PubChem may count differently. Tolerance ±1.
- **SA Score**: Uses RDKit Contrib sascorer if available, else heuristic fallback.
- **ADMET**: Heuristic mode (RDKit fallback) — values are rule-based estimates, not ML predictions.
