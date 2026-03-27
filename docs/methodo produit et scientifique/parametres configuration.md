# BindX definitive parameterization guide for in silico drug discovery

**Every computational method in the BindX stack — GNINA, AutoDock Vina, RDKit, REINVENT 4, AiZynthFinder, ADMET-AI, P2Rank, and OpenFE — has parameters whose optimal values are well-established in the literature but rarely compiled in one place.** This guide distills systematic benchmarks, sensitivity analyses, and best-practices papers into a single actionable configuration manual. Each recommendation is backed by specific benchmark numbers and citations. Where the literature disagrees, both positions are presented. The guide follows a consistent format: safe defaults, context-specific variations, performance benchmarks, pitfalls, and citations.

---

## 1. AutoDock Vina 1.2+ produces best results at exhaustiveness 32 for pose prediction

**Exhaustiveness** is Vina's most impactful parameter. Agarwal & Smith (2023, *Mol. Inform.* 42:e2200188) systematically tested values of 1, 8, 25, 50, 75, and 100 on PDBbind v2017 and found that **exhaustiveness = 8 is sufficient for virtual screening**, while gains plateau above 25. Eberhardt et al. (2021, *J. Chem. Inf. Model.* 61:3891) recommend **exhaustiveness = 32** for challenging ligands such as imatinib. Computation scales linearly — doubling exhaustiveness doubles wall-clock time (Trott & Olson, 2010, *J. Comput. Chem.* 31:455).

| Use case | Exhaustiveness | num_modes | energy_range |
|---|---|---|---|
| Virtual screening (HTS) | 8 | 9 | 3 kcal/mol |
| Pose prediction / re-docking | 32 | 9–20 | 3 kcal/mol |
| Flexible/large ligands (>8 RotB) | 32–64 | 20 | 4 kcal/mol |
| Peptide docking | 128–256 | 20–50 | 5 kcal/mol |

**Box size** should follow the Feinstein & Brylinski (2015, *J. Cheminform.* 7:18) formula: **box dimension = 2.857 × ligand radius of gyration** in each axis, corresponding to an Rg/box ratio of 0.35. This improved EF1% from 7.67 to **8.20** and improved ranking for ~2/3 of targets. The Vina manual warns at volumes exceeding 27,000 Å³ (~30×30×30 Å).

**Scoring function selection** depends on the task. On CASF-2016, Vina achieves **90.2% top-1 docking success** (RMSD < 2 Å) but only Pearson r ≈ 0.60 for affinity prediction (Su et al., 2019, *J. Chem. Inf. Model.* 59:895). Vinardo (Quiroga & Villarreal, 2016, *PLOS ONE* 11:e0155183) outperforms Vina in ranking power (**70.8% high-success** vs 58.5% on CASF-2013). The AD4 scoring function excels at metal-containing sites due to explicit electrostatics but is **3× slower** (Eberhardt 2021). On DUD-E, Vina AUC = 68 ± 4.7 and AD4 AUC = 66.4 ± 10.2.

**Ligand flexibility** degrades accuracy above **8 rotatable bonds** across all docking programs (Erickson et al., 2004, *J. Med. Chem.* 47:45). Below 8 rotatable bonds, success rates reach 95%; above 8, only CDOCKER maintained 71% accuracy. **Receptor flexibility** should be limited to **1–4 side chains** — an Amaro lab study on VEGFR-2 found 1 flexible residue (Glu885) improved Pearson r to 0.568 with 51% increased processing time, while adding more residues did not consistently help.

**Reproducibility** requires fixing the random seed, but Jaghoori et al. (2016, *J. Comput. Aided Mol. Des.* 30:237) showed that even identical seeds on different operating systems produce energy differences of **0.7 kcal/mol** — always record seed, OS, compiler, and Vina version. For protonation, use **Dimorphite-DL** (Ropp et al., 2019, *J. Cheminform.* 11:14) to enumerate ionization states at pH 6.8–7.4, then convert to PDBQT using Meeko.

---

## 2. GNINA's CNN rescoring doubles enrichment over Vina at negligible speed cost

GNINA's default **`--cnn_scoring rescore`** mode runs Vina's Monte Carlo sampling, then rescores and re-ranks final poses with a CNN ensemble. This achieves **73% top-1 success** in redocking versus 58% for Vina (McNutt et al., 2021, *J. Cheminform.* 13:43) — at essentially the same speed as Vina alone when a GPU is available. The full CNN-driven mode (`--cnn_scoring all`) is **~1000× slower** and not recommended.

**CNN model selection** matters for different tasks. On DUD-E, the Dense ensemble achieves median AUC = **0.80** versus 0.79 for the default ensemble (Sunseri & Koes, 2021, *Molecules* 26:7369). GNINA outperforms Vina on **89 of 117 targets** across DUD-E and LIT-PCBA, with **median EF1% more than twice that of Vina**. The CNN affinity score provides better virtual screening enrichment than the pose score; the combined CNN_VS score (pose × affinity) provides a modest additional boost.

**GNINA 1.3** (McNutt et al., 2025, *J. Cheminform.* 17:28) migrated from Caffe to PyTorch, cutting CPU runtime from 129s to 23s (5.6× speedup) while improving cross-docking from 37% to **40% top-1**. The new `--cnn=fast` single-model option runs in **16s on CPU** — only 1.3s slower than Vina. Redocking showed a slight regression (69% → 67%), indicating the 1.3 ensemble traded marginal redocking accuracy for better generalization.

| Configuration | Speed (CPU, 4 cores) | Speed (GPU, RTX 3080) | Redocking top-1 | Cross-docking top-1 |
|---|---|---|---|---|
| Vina (no CNN) | 14.7s | N/A | 58% | 27% |
| GNINA 1.3 default ensemble | 23s | ~8s | 67% | 40% |
| GNINA 1.3 `--cnn=fast` | 16s | ~15s | ~65% | ~38% |
| GNINA 1.0 (Caffe) | 129s | ~8s | 69% | 37% |

**`autobox_add` of 4 Å** is optimal when the binding site is known — increasing it actually *decreases* pose accuracy (McNutt 2021). For uncertain sites, increase to 8–12 Å. The **`min_rmsd_filter`** of 1.0 Å (default) does not significantly affect docking performance. For exhaustiveness, **GNINA at 8 outperforms Vina at 16** because the CNN rescore corrects for suboptimal sampling.

**Uni-Dock** (Yu et al., 2023, *J. Chem. Theory Comput.* 19:3336) delivers **>1000× speedup** over Vina on a single V100 GPU, docking at ~0.10 s/ligand in fast mode with **no measurable accuracy loss** on CASF-2016 and DUD-E. A hierarchical strategy (Fast → top 10% → Balanced → top 10% → Detailed) screened **38.2M molecules in ~12 hours on 100 V100 GPUs**. The optimal pipeline for ultra-large campaigns: Uni-Dock fast mode → GNINA CNN rescore of top 1–5%.

---

## 3. Consensus scoring improves EF1% by 30–50% when combining complementary functions

Rank-based consensus methods are preferred over score-based methods because scores have variable units, scales, and offsets (Palacio-Rodriguez et al., 2019, *Sci. Rep.* 9:5142). **Exponential Consensus Ranking (ECR)** assigns P_j(i) = exp(−rank_j(i)/σ), where **σ ≈ 5% of dataset size**, and sums across programs. ECR results are essentially independent of σ, a major advantage over vote-based methods.

The optimal number of programs is **3–5**, with saturation beyond 5–7 (Towards Effective Consensus Scoring, *Interdiscip. Sci.* 2023). Ericksen et al. (2017, *J. Chem. Inf. Model.* 57:1579) showed mean consensus EF1 = **26** versus best individual program EF1 = 18 — a **44% improvement**. Only complementary/orthogonal scoring functions should be combined; correlated functions provide redundant information (Yang et al., 2005). The recommended open-source combination for BindX is **AutoDock Vina + Vinardo + GNINA CNN** — these span empirical, modified-empirical, and ML-based scoring with low correlation.

For **ensemble docking**, **3–5 receptor conformations** is the practical sweet spot. Larger ensembles increase false positives (Sci. Rep. 2021, s41598-021-04448-5). AlphaFold structures perform **consistently worse** than experimental PDB structures for virtual screening (Scardino et al., *iScience* 2023), though state-specific AF2 modeling (e.g., DFGin/DFGout for kinases) can improve results. Score aggregation should use **best-score (minimum)** for initial enrichment, coupled with ECR for consensus across conformations. Including energy penalties from MD simulations is essential — otherwise high-energy states dominate docking (Fischer et al., *PNAS* 2021).

---

## 4. ML scoring functions dominate empirical ones on CASF-2016 but require careful validation

The CASF-2016 benchmark (Su et al., 2019) on 285 complexes established the definitive comparison. Classical scoring functions peak at Pearson r ≈ 0.63 (X-Score HM) for affinity prediction, while ML-based functions reach r = **0.86** (OnionNet-2, TopBP). For docking power, ChemPLP@GOLD leads at ~72–74% success rate among empirical functions. GNINA CNN rescore achieves **73% redocking success** and median DUD-E AUC = 0.79.

| Scoring Function | Scoring Power (r) | Docking Power | Screening EF1% | Type |
|---|---|---|---|---|
| OnionNet-2 | 0.864 | — | — | ML (CNN) |
| TopBP | 0.861 | — | — | ML (topology) |
| GNINA Default | ~0.6 (estimated) | 73% (redock) | 2× Vina | ML (CNN) |
| RTMScore | High | Excellent | High | ML (GNN+MDN) |
| X-Score(HM) | 0.631 | — | — | Empirical |
| ChemPLP@GOLD | ~0.58 | 72–74% | ~11.9 | Empirical |
| AutoDock Vina | ~0.48 | ~58% | ~7.7 | Empirical |

**Task-specific recommendations**: For **pose prediction**, use GNINA CNN or ChemPLP. For **affinity ranking**, use OnionNet-2 or ΔvinaXGB. For **virtual screening**, use GNINA CNN_VS with OnionNet-SFCT+Vina as rescorer (EF1% = **15.5** on CASF-2016). The critical caveat: ML scoring functions show inflated performance when training/test similarity is high (Li & Yang, 2017), though recent work shows they still outperform classical counterparts on truly blind benchmarks.

---

## 5. ECFP4 at 2048 bits with Tanimoto similarity remains the optimal default fingerprint

Riniker & Landrum (2013, *J. Cheminform.* 5:26) benchmarked 14 fingerprints across 88 targets and found **ECFP4 (radius 2) had the best mean rank**, though the difference from ECFP6 was not statistically significant (r² = 0.999 correlation). Probst & Reymond (2018, *J. Cheminform.* 10:49) confirmed this across the same benchmark set.

For bit length, Landrum's 2023 collision analysis on 4M PubChem compounds showed that **2048-bit Morgan r=2** maintains Spearman r = 0.995 with the unfolded fingerprint (mean similarity deviation = 0.0045). At 1024 bits, r drops to 0.990 — still excellent for most applications. At 512 bits (r = 0.979), distortion becomes meaningful. Morgan fingerprints have inherently low collision rates (~2.4% bit density at 2048 bits) compared to path-based fingerprints.

**Count vs binary**: Use **binary for tree-based models** (RF, XGBoost, LightGBM) where splits only need thresholds. Use **count for Bayesian and linear models** where frequency information adds discriminative power. FCFP (functional-class) fingerprints are preferred over ECFP for **scaffold hopping** due to pharmacophoric abstraction (Rogers & Hahn, 2010, *J. Chem. Inf. Model.* 50:742).

**Tanimoto similarity** is the gold standard metric — Bajusz et al. (2015, *J. Cheminform.* 7:20) compared 8 metrics and identified Tanimoto, Dice, Cosine, and Soergel distance as the best choices. Tanimoto and Dice produce identical rankings. For substructure-like searching, use asymmetric Tversky (α=0.9, β=0.1).

---

## 6. Random Forest QSAR works best with 500 trees, unlimited depth, and min_samples_leaf = 5

Svetnik et al. (2003, *J. Chem. Inf. Comput. Sci.* 43:1947) established RF as the QSAR workhorse. Performance plateaus at **100–300 trees** (Kensert et al., 2018, *J. Cheminform.* 10:49), with 500 providing a robust default. For `max_features`, use **√n** for classification and **n/3** for regression — both well-validated defaults for ECFP features (Sheridan et al., 2021, arXiv:2105.08626).

**Do not cap `max_depth`** in Random Forests — unlike boosting, RF controls overfitting through ensemble averaging, and fully grown trees reduce bias (Probst et al., 2018, arXiv:1804.03515). Instead, set **`min_samples_leaf = 5`** to prevent terminal nodes from memorizing individual compounds (Sheridan 2021 used this across all 30 pharmaceutical datasets). For class imbalance, **`class_weight='balanced'`** is the simplest effective approach. SMOTEENN (SMOTE + Edited Nearest Neighbors) outperformed all other resampling methods on 12 imbalanced Tox21 bioassays (Sun et al., 2020, *J. Cheminform.* 12:54), but adds complexity. Use **permutation importance** over Gini importance for interpreting ECFP bits (Strobl et al., 2007), and always verify important bits by examining their substructural meaning via RDKit's `bitInfo`.

**XGBoost/LightGBM** requires more careful tuning. Sheridan et al. (2021) found **learning_rate = 0.05** with n_estimators = 2000 (with early stopping) works across 30 pharmaceutical datasets. Boldini et al. (2023, *J. Cheminform.* 15:73) benchmarked 157,590 gradient boosting models and found XGBoost outperforms LightGBM by ~5% but LightGBM is **~100× faster** for large datasets. For ECFP features, **max_depth = 3–6** is optimal; for continuous descriptors, 6–8. The recommended Optuna search space: learning_rate [0.01, 0.3] log-uniform, max_depth [3, 8], subsample [0.5, 1.0], colsample_bytree [0.3, 1.0], reg_alpha and reg_lambda [1e-8, 10.0] log-uniform. Use **≥50–100 Optuna trials** with early stopping (patience = 30–50).

---

## 7. UMAP-based splits now surpass scaffold splits as the gold standard for realistic benchmarking

Scaffold splitting using **generic (CSK) Bemis-Murcko scaffolds** (Bemis & Murcko, 1996, *J. Med. Chem.* 39:2887) is the minimum standard — specific scaffolds fragment datasets excessively (average ~1.93 molecules per scaffold). However, Guo et al. (2025, *J. Cheminform.* 17:79) demonstrated that UMAP-based clustering splits (ECFP4, 2048 bits → UMAP with Jaccard distance, n_neighbors=80, min_dist=0.95 → agglomerative clustering with K=7) produce the **most challenging and realistic evaluations**. GEM median hit rate with Butina split approximately *doubles* that with UMAP split. Tetko & Clevert (*JCIM* 2025) independently confirmed UMAP as most challenging, with the largest ΔROC-AUC (median 0.088).

**Temporal splitting is preferred when timestamps are available** — Sheridan (2013, *J. Chem. Inf. Model.* 53:783) showed time-split R² better approximates prospective prediction than random splits (optimistic) or leave-class-out (pessimistic). The SIMPD algorithm (Landrum et al., 2023, *J. Cheminform.* 15:119) can simulate temporal-like splits for public datasets. The hierarchy is: **temporal > UMAP > scaffold > random**. Default split ratio: **80/10/10**; for small datasets (<1000) or conformal prediction, use **70/15/15** or **60/20/20**.

For **applicability domain**, use Tanimoto similarity **≥ 0.4 on ECFP4** as the primary threshold (Sheridan et al., 2004, *J. Chem. Inf. Comput. Sci.* 44:1912). Below 0.3 is background noise for structurally unrelated molecules (Godden et al., 2000). Implement Sheridan's **3D applicability domain** (similarity × RF tree variance × predicted value) as a calibrated lookup table — this was shown to be more discriminative than any single metric across diverse and narrow training sets (Sheridan, 2012, *JCIM* 52:814; 2013, *JCIM* 53:2837; 2015, *JCIM* 55:1098).

---

## 8. Conformal prediction should always use Mondrian CP for drug discovery classification

Standard (non-Mondrian) conformal prediction systematically under-covers minority classes in imbalanced bioactivity datasets. **Mondrian CP** generates separate nonconformity score distributions per class, ensuring validity holds per-class (Sun, Carlsson, Ahlberg et al., 2017, *JCIM* 57:1591). For classification, always use Mondrian CCP or ACP with 10 iterations.

**Significance level α** should be context-dependent: **α = 0.20** for hit triage/virtual screening (tolerate 20% errors for throughput), **α = 0.10** for lead optimization (standard), **α = 0.05** for safety/toxicity assessment (Alvarsson et al., 2021, *J. Pharm. Sci.* 110:42). Always present calibration plots across multiple significance levels. The minimum calibration set is **≥200 total examples** (≥100 per class for Mondrian CP). For regression, use **normalized nonconformity scores** with a difficulty estimator to produce adaptive-width intervals (Svensson et al., 2018, *JCIM* 58:1132).

For iterative **DMTA cycles**, standard CP's exchangeability assumption breaks under distribution shift. **Adaptive Conformal Inference (ACI)** (Gibbs & Candès, 2021, NeurIPS) continuously adjusts the working significance level via online gradient descent, achieving desired coverage regardless of distribution shifts. The simplest practical approach is to **update calibration sets with recent experimental data at each cycle** — Morger et al. (2021, *J. Cheminform.* 13:35) showed this restores validity when temporal distribution shift breaks CP calibration.

---

## 9. REINVENT 4's default σ = 128 with diversity filters provides balanced generation

REINVENT 4 (Loeffler et al., 2024, *J. Cheminform.* 16:20) uses the augmented likelihood formulation: log P_aug(T) = log P_prior(T) + σ·S(T). The **default σ = 128** balances exploration and exploitation. Higher σ (200–256) increases exploitation at the risk of mode collapse; lower σ (64) keeps the agent closer to the prior. The Adam optimizer learning rate defaults to 0.0001.

The **diversity filter** with IdenticalMurckoScaffold type, **bucket_size = 25**, and minscore = 0.4 prevents mode collapse while allowing adequate scaffold exploration. Each scaffold bucket holds up to 25 molecules; once full, further instances score zero. Additionally, a global SMILES memory of size 1 prevents exact duplicate generation. For Mol2Mol optimization, use PenalizeSameSmiles with penalty_multiplier = 0.5 instead.

**Scoring function design** is critical. Always include drug-likeness constraints (QED or explicit Lipinski/Veber rules) when optimizing docking scores — docking rewards large, lipophobic molecules (Guo et al., 2021). Use **geometric mean** aggregation (recommended) over arithmetic mean. Penalty components (PAINS, structural alerts) are applied multiplicatively. Transform functions must match the property: reverse_sigmoid for "lower is better" (docking), double_sigmoid for range constraints (MW 200–500).

**Augmented Memory** (Guo & Schwaller, 2024, *JACS Au* 4:2160) dramatically prevents mode collapse by combining experience replay with SMILES augmentation and selective memory purge. Key parameters: replay buffer top-k = 100, SMILES augmentation repetitions = 5–10, multiple gradient updates per batch = 8. This achieves state-of-the-art on the PMO benchmark with only 9,600 oracle calls for docking-based MPO.

| Parameter | Default | Recommended | Citation |
|---|---|---|---|
| `batch_size` | 128 | 128 (256 for complex scoring) | Loeffler 2024 |
| `sigma` | 128 | 128 | Loeffler 2024 |
| `rate` (Adam lr) | 0.0001 | 0.0001 | Source code |
| Diversity bucket_size | 25 | 25–50 | Source code |
| `min_steps` | 25 | 25 | Source code |
| `max_steps` | 100 | 100–1000 (task-dependent) | Loeffler 2024 |

For **staged learning**, define multiple scoring stages of increasing complexity: Stage 1 (drug-likeness + structural alerts, 100 steps) → Stage 2 (add docking/QSAR, 500 steps) → Stage 3 (add ADMET, 300 steps). Transfer learning before RL (e.g., 35–200 epochs on known actives) can match experience replay performance.

---

## 10. AiZynthFinder finds most drug routes in 3–5 steps with default MCTS parameters

AiZynthFinder (Genheden et al., 2020, *J. Cheminform.* 12:70) uses MCTS with **C = 1.4** (UCB exploration constant), **max_transforms = 6**, **iteration_limit = 100**, **time_limit = 120s**, and **cutoff_number = 50** (expansion width). These defaults find solutions for the majority of ChEMBL drug-like molecules, typically in **3–5 retrosynthetic steps**. The search finds solutions in <10 seconds and completes in <1 minute for typical targets.

The **RAscore** (Thakkar et al., 2021, *Chem. Sci.* 12:3339) provides a binary synthesizability classifier at **4,500× the speed** of full retrosynthetic analysis. Use a threshold of **0.5** for binary classification. Reference values: imatinib RAscore = 0.995 (easily synthesizable), morphine RAscore = 8.3×10⁻⁷ (complex ring systems). Stock configuration supports ZINC (default), Enamine REAL, and custom catalogs via InChI key HDF5 files or MolBloom.

---

## 11. ADMET thresholds follow well-established pharmacological cutoffs

**ADMET-AI** (Swanson et al., 2024, *Bioinformatics* 40:btae416) uses Chemprop-RDKit with an ensemble of **5 models** across 41 TDC ADMET datasets. It achieved the best average rank on the TDC ADMET Benchmark Group with 45% speed reduction versus competitors.

Key endpoint thresholds for BindX configuration:

| Endpoint | Threshold | Interpretation | Citation |
|---|---|---|---|
| hERG IC50 | >10 μM | Low cardiac risk; safety margin ≥30× free Cmax | Redfern 2003; Leishman 2020 |
| logBB (BBB) | ≥ −1 (screening), >0 (CNS drugs) | BBB-permeable | MoleculeNet standard |
| Caco-2 Papp | >10 × 10⁻⁶ cm/s | High permeability | Pharmacological consensus |
| CYP inhibition | IC50 > 10 μM (all isoforms) | Low DDI risk | FDA guidance 2020 |
| HLM CLint | <8.6 μL/min/mg | Low clearance | In vitro-in vivo scaling |
| PPB (fu) | >0.1 | Adequate free fraction | Pharmacological consensus |

**CNS MPO** (Wager et al., 2010, *ACS Chem. Neurosci.* 1:435) scores 6 parameters with linear desirability functions summing to 0–6: ClogP (≤3 = 1.0, ≥5 = 0.0), ClogD (≤2 = 1.0, ≥4 = 0.0), MW (≤360 = 1.0, ≥500 = 0.0), TPSA (40–90 optimal), HBD (≤0.5 = 1.0, ≥3.5 = 0.0), pKa (≤8 = 1.0, ≥10 = 0.0). The threshold **CNS MPO ≥ 4.0** captures 74% of marketed CNS drugs and 94% of post-implementation Pfizer CNS candidates.

**Lipophilic Efficiency** (LipE = pIC₅₀ − cLogP) targets **≥5 for lead compounds** and ≥6 for clinical candidates (Leeson & Springthorpe, 2007, *Nat. Rev. Drug Discov.* 6:881). QED (Bickerton et al., 2012, *Nat. Chem.* 4:90) with **QEDw,mo ≥ 0.40** offers 48% greater specificity than Ro5. For Pareto optimization, NSGA-II with population size = 100 and 150–200 generations is the standard (Deb et al., 2002; Fromer & Coley, 2022, *Chem. Sci.*).

---

## 12. OpenFE achieves MUE ~1.2 kcal/mol with Sage and 5 ns per window

**Force field**: OpenFF 2.0 Sage with AM1-BCC charges paired with AMBER ff14SB for protein and TIP3P water is the recommended open-source stack. Hahn et al. (2024, *JCIM*) found Sage achieves ΔΔG RMSE = **1.7 kcal/mol** and MUE = **1.2 kcal/mol** across 22 targets. A consensus averaging Sage + GAFF2.1x + CGenFF approaches the accuracy of proprietary OPLS3e. Lee et al. (2022, *Front. Mol. Biosci.*) specifically found ff14SB + GAFF2.2 + TIP3P with AM1-BCC charges gave MUE = **0.89 kcal/mol** — and RESP charges did *not* improve over AM1-BCC.

**Lambda windows**: Use **11 windows** (OpenFE default with HREX) for standard perturbations (≤4 heavy atom changes). Increase to **16–20 windows** for larger perturbations — Clarke et al. (2022, *JCIM* 62:2068) found 20 windows optimal for BRD4 inhibitors. Monitor overlap matrix superdiagonal elements: values should be **> 0.03** (Midgley et al., 2025, *JCIM*). **Simulation length** of 5 ns/window is the standard minimum (Wang et al., 2015, *JACS* 137:2695). Running **3–5 independent replicas** at 5 ns is preferable to a single 20 ns run.

**Perturbation networks** should use LOMAP with cycle closure (Liu et al., 2013, *J. Comput.-Aided Mol. Des.* 27:755). Set LOMAP score threshold > **0.4** and limit perturbations to **≤ 4–6 heavy atom changes** per edge. Avoid ring-breaking transformations and charge-changing perturbations without Rocklin-style PB corrections (Rocklin et al., 2013, *J. Chem. Phys.* 139:184103). FEP should **not** be used for scaffold hopping, large R-group changes (>6–8 heavy atoms), metal binding sites, covalent ligands, or cases with different binding modes between ligands.

The maximal achievable accuracy is limited by **experimental reproducibility** (~0.4–0.64 kcal/mol MUE; Hahn et al., 2023, *Commun. Chem.* 6:202). Current best FEP methods achieve pairwise RMSE ~1.1 kcal/mol, approaching this limit.

---

## 13. MM-GBSA works best as a ranking tool with igb=5, single trajectory, and no entropy

**GB model**: Use **igb=5 (GBOBC2)** with mbondi radii as the robust default (Wang et al., *Chem. Rev.* 2019, 119:9478) or **igb=8 (GBneck2)** with mbondi3 radii as the modern alternative (Nguyen et al., 2013, *J. Chem. Theory Comput.* 9:2020). For docking pose rescoring, igb=2 (GBOBC1) with internal dielectric εin = 2 gave the highest success rate of **69.4%** in a 98-complex benchmark (Hou et al., 2011, *JCIM* 51:69).

**Internal dielectric** is context-dependent: **εin = 1** for hydrophobic binding sites, **εin = 2** as a general default for docking rescoring, **εin = 4** for highly charged interfaces. Chen et al. (2020, *JCIM* 60:5353) showed that variable dielectric (εin=4 for charged residues, 1 for nonpolar) outperformed uniform dielectric.

**Simulation length** of **1–5 ns MD** with extraction of 100–200 frames is the production standard. Multiple independent short runs (5 × 1 ns) outperform a single long trajectory (Genheden & Ryde, 2015, *Expert Opin. Drug Discov.* 10:449). Use the **single-trajectory** approach as default — ΔEinternal cancels between complex, receptor, and ligand terms, significantly reducing noise. The 3-trajectory approach is 3× more expensive and typically *worsens* ranking precision.

**Entropy should be omitted** for ranking congeneric series of similar-sized ligands — it largely cancels (Sun et al., 2018, *Phys. Chem. Chem. Phys.* 20:14450). For diverse ligands, the **interaction entropy** method (Duan et al., 2016, *JACS* 138:5722) requires zero additional computation, but Ryde & Genheden (2021, *JCTC* 17:5187) showed it is poorly conditioned when σ_IE > 15 kJ/mol. Typical ranking accuracy: Pearson r = **0.4–0.7** for well-curated congeneric series, vastly inferior to FEP (r = 0.75; Wang 2015).

---

## 14. Protein preparation requires careful histidine assignment and restrained minimization

Use **PropKa 3.1 at pH 7.4** for standard drug targets. Adjust for compartment-specific targets: pH 4.5–5.0 for lysosomal enzymes (cathepsins), pH 1–2 for gastric targets (pepsin), pH 5.5–6.5 for endosomal targets (Olsson et al., 2011, *JCTC* 7:525). Histidine protonation state significantly affects virtual screening enrichment — AUC differences of 0.87 to 0.97 were reported for a single target (Kim et al., 2013, PMC3639364). Decision criteria: if PropKa pKa > pH → HIP (+1); if pKa < pH → neutral (HID or HIE based on local hydrogen bond network). Always inspect the local environment of catalytic and binding-site histidines.

For **missing loops**, generate **≥100–300 MODELLER models** for research quality (up to 10,000 with clustering for highly flexible loops) and select the model with the **lowest DOPE score** — zDOPE < −1.0 indicates high confidence. Validate with Ramachandran analysis (>90% favored).

**Minimization protocol**: 500–5,000 steps steepest descent (removes worst clashes) followed by 500–10,000 steps conjugate gradient (refines geometry), with heavy-atom restraints of 25–50 kcal/mol/Å² progressively released. The Schrödinger standard is restrained minimization to RMSD convergence of **0.30 Å**. Use **AMBER ff14SB** as the force field for docking preparation (Maier et al., 2015, *JCTC* 11:3696).

For **water retention**, keep waters making **≥2 hydrogen bonds to non-water atoms** and within **5.0 Å of any ligand atom**. Remove waters with B-factor > 40–60 Å². WaterMap (Abel et al., 2008, *JACS* 130:2817) or open-source GIST provide thermodynamic analysis of hydration sites. Conserved waters (positions within 1.0–1.5 Å across multiple crystal structures) should always be retained.

---

## 15. P2Rank with the AlphaFold model achieves >80% success on predicted structures

P2Rank (Krivák & Hoksza, 2018, *J. Cheminform.* 10:39) uses the **`-c alphafold`** configuration for AlphaFold/cryo-EM/NMR structures, as the default model uses B-factor features that are repurposed for pLDDT in AF2. On the HOLO4K dataset, the default model achieves ~73% top-n success rate and ~85% top-(n+2). On AlphaFold structures: **81% success rate** (206/256), with >90% when considering top-3 predicted pockets.

For docking box construction from P2Rank output, center the box on the predicted pocket center with **10 Å padding** per side. Run **fpocket + P2Rank rescoring** (`prank rescore fpocket.ds`) for consensus — this achieves 72% recall and 47% precision for real pockets (PocketVec, *Nature Commun.* 2024). Druggable pockets should have fpocket drug score > **0.5** and contain **≥14 residues** with volume **300–1500 Å³**.

---

## 16. Apply PAINS + Brenk as minimum; relax Ro5 for macrocycles and PROTACs

The PAINS filter set (Baell & Holloway, 2010, *J. Med. Chem.* 53:2719) contains **480 substructure patterns** across three families (A: ~16 highest-confidence, B: ~55, C: ~409). Apply all three families as initial screen but **flag rather than auto-reject** — the Eli Lilly study (PMC6088356) found only 3 alerts show true pan-assay promiscuity. The **Brenk filter** (Brenk et al., 2008, *ChemMedChem* 3:435) adds 105 structural alerts for reactivity, toxicity, and metabolic instability. In RDKit, combine PAINS + BRENK catalogs as the minimum filter stack; add NIH and CHEMBL_MLSMR for HTS triage.

**Ro5 relaxation** for beyond-Rule-of-5 (bRo5) space (Doak et al., 2014, *Chem. Biol.* 21:1115): MW up to 3000, cLogP −2 to +10, but **HBD ≤ 6 remains the most restrictive hard limit** for oral absorption. Macrocycles achieve oral bioavailability through chameleonicity (intramolecular H-bonding in non-polar environments). PROTACs (MW 700–1200) can achieve oral absorption with compact E3 ligase ligands. For fragment screening, apply the Rule of Three (MW < 300, cLogP ≤ 3, HBD ≤ 3, HBA ≤ 3).

---

## 17. Activity cliff detection uses SALI with ΔpKi ≥ 2 and fingerprint-appropriate Tc thresholds

The **Structure-Activity Landscape Index** (SALI(i,j) = |A_i − A_j| / (1 − sim(i,j))) quantifies activity cliffs, with higher values indicating more dramatic potency differences between structurally similar compounds (Guha & Van Drie, 2008, *JCIM* 48:646). The significance threshold is the **90th or 95th percentile** of the SALI distribution within the dataset.

For activity cliff detection, Stumpfe & Bajorath (2012, *J. Med. Chem.* 55:2932; 2014, *J. Med. Chem.* 57:18) established **ΔpKi ≥ 2** (100-fold potency difference) as the consensus threshold. The Tanimoto similarity threshold depends on the fingerprint: **Tc ≥ 0.55 for MACCS keys** or **Tc ≥ 0.40 for ECFP4** — these yield comparable chemical neighborhoods. For MMP-based cliffs, use single-cut Hussain-Rea fragmentation (Hussain & Rea, 2010, *JCIM* 50:339) via the **mmpdb** tool, which handles stereochemistry and provides environment fingerprints (Dalke et al., 2018, *JCIM*).

Second-generation thresholds are target-set-dependent: compute the **mean + 2σ of the potency difference distribution** within the target's compound pairs. This adapts thresholds to datasets where activity ranges vary — for example, neurokinin 1 receptor requires only ΔpKi = 1.4 (Stumpfe & Bajorath, 2019, *ACS Omega* 4:14360).

---

## Conclusion: an integrated BindX configuration philosophy

The overarching pattern across all 22 methods is that **literature-validated defaults work well for >80% of cases**, but systematic parameter studies reveal 10–50% performance gains from context-specific optimization. Three principles unify the optimal BindX configuration: First, **consensus beats any single method** — combining Vina + GNINA + Vinardo scoring, multiple ADMET predictors, and ensemble receptor conformations consistently improves reliability even when individual components are imperfect. Second, **ML-based methods require honest validation** — scaffold splits are a minimum standard, but UMAP or temporal splits with conformal prediction provide the most realistic performance estimates. Third, **computational budget allocation matters** — Uni-Dock for initial screening (0.1 s/ligand), GNINA CNN rescoring for hit confirmation (8 s/ligand on GPU), and OpenFE for lead optimization (days per compound) represent the correct escalation of computational rigor matching the decision stage.

The single most common pitfall across methods is **overconfidence from inadequate splitting**: random splits overestimate QSAR model performance, rescoring on the same poses used for training inflates ML scoring function accuracy, and ADMET predictions outside their applicability domain should be flagged rather than trusted. Every prediction in the BindX pipeline should carry an uncertainty estimate — whether from conformal prediction, ensemble variance, or Sheridan's multi-dimensional applicability domain — and compounds outside these bounds should be triaged accordingly.