// ---------------------------------------------------------------------------
// BindX V9 — Mock data for frontend development
// All mock data follows CDC_V9.md schema: Project > Campaign > Phase > Run
// ---------------------------------------------------------------------------

// ---- Molecules with progressive enrichment ----
// Molecules start with basic info (import), then gain columns from subsequent runs

function mol(id, smiles, name, extra = {}) {
  return {
    id: `mol-${id}`,
    smiles,
    canonical_smiles: smiles,
    name,
    bookmarked: false,
    generation_level: 0,
    parent_molecule_id: null,
    source_run_id: 'run-a1',
    ...extra,
  }
}

// 20 realistic molecules for EGFR Hit Discovery (Phase A)
const PHASE_A_MOLECULES = [
  mol(1, 'COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1', 'Gefitinib', {
    bookmarked: true, source_run_id: 'run-a1',
    docking_score: -9.2, cnn_score: 0.87, cnn_affinity: 7.1, cnn_vs: 6.18,
    logP: 3.2, MW: 446.9, HBD: 1, HBA: 7, TPSA: 68.7, QED: 0.72,
    solubility: 0.65, BBB: 0.3, hERG: 0.12, metabolic_stability: 0.7, lipinski_pass: true,
    composite_score: 82.5, cluster_id: 1, scaffold: 'quinazoline', interactions_count: 6,
    admet: { absorption: 0.82, distribution: 0.65, metabolism: 0.72, excretion: 0.55, toxicity: 0.88, pharmacodynamics: 0.7 },
    svg_2d: null,
  }),
  mol(2, 'C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1', 'Erlotinib', {
    bookmarked: true, source_run_id: 'run-a1',
    docking_score: -8.8, cnn_score: 0.82, cnn_affinity: 6.9, cnn_vs: 5.66,
    logP: 2.9, MW: 393.4, HBD: 1, HBA: 7, TPSA: 74.7, QED: 0.68,
    solubility: 0.55, BBB: 0.25, hERG: 0.08, metabolic_stability: 0.65, lipinski_pass: true,
    composite_score: 79.1, cluster_id: 1, scaffold: 'quinazoline', interactions_count: 5,
    admet: { absorption: 0.78, distribution: 0.6, metabolism: 0.68, excretion: 0.5, toxicity: 0.92, pharmacodynamics: 0.65 },
  }),
  mol(3, 'CS(=O)(=O)CCNCc1ccc(-c2ccc3ncnc(Nc4ccc(OCc5cccc(F)c5)c(Cl)c4)c3c2)o1', 'Lapatinib', {
    bookmarked: true, source_run_id: 'run-a1',
    docking_score: -9.5, cnn_score: 0.91, cnn_affinity: 7.4, cnn_vs: 6.73,
    logP: 4.6, MW: 581.1, HBD: 2, HBA: 7, TPSA: 114.7, QED: 0.45,
    solubility: 0.35, BBB: 0.15, hERG: 0.22, metabolic_stability: 0.55, lipinski_pass: false,
    composite_score: 74.8, cluster_id: 2, scaffold: 'quinazoline-furan', interactions_count: 8,
    admet: { absorption: 0.55, distribution: 0.45, metabolism: 0.6, excretion: 0.4, toxicity: 0.72, pharmacodynamics: 0.8 },
  }),
  mol(4, 'COc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCC1CCN(C)CC1', 'Vandetanib', {
    source_run_id: 'run-a1',
    docking_score: -8.1, cnn_score: 0.75, cnn_affinity: 6.5, cnn_vs: 4.88,
    logP: 3.8, MW: 475.4, HBD: 1, HBA: 5, TPSA: 59.5, QED: 0.58,
    solubility: 0.42, BBB: 0.35, hERG: 0.28, metabolic_stability: 0.6, lipinski_pass: true,
    composite_score: 68.3, cluster_id: 1, scaffold: 'quinazoline', interactions_count: 4,
    admet: { absorption: 0.7, distribution: 0.55, metabolism: 0.62, excretion: 0.48, toxicity: 0.68, pharmacodynamics: 0.6 },
  }),
  mol(5, 'C=CC(=O)Nc1cc(Nc2nccc(-c3cn(C)c4ccccc34)n2)c(OC)cc1N(C)CCN(C)C', 'Osimertinib', {
    bookmarked: true, source_run_id: 'run-a1',
    docking_score: -9.8, cnn_score: 0.93, cnn_affinity: 7.8, cnn_vs: 7.25,
    logP: 3.4, MW: 499.6, HBD: 2, HBA: 6, TPSA: 87.6, QED: 0.52,
    solubility: 0.48, BBB: 0.42, hERG: 0.18, metabolic_stability: 0.72, lipinski_pass: true,
    composite_score: 88.2, cluster_id: 3, scaffold: 'pyrimidine-indole', interactions_count: 7,
    admet: { absorption: 0.75, distribution: 0.7, metabolism: 0.74, excretion: 0.52, toxicity: 0.82, pharmacodynamics: 0.85 },
  }),
  mol(6, 'CCOc1cc2ncc(C#N)c(Nc3ccc(OCc4ccccn4)c(Cl)c3)c2cc1NC(=O)/C=C/CN(C)C', 'Neratinib', {
    source_run_id: 'run-a1',
    docking_score: -9.1, cnn_score: 0.84, cnn_affinity: 7.0, cnn_vs: 5.88,
    logP: 4.2, MW: 557.0, HBD: 2, HBA: 8, TPSA: 120.7, QED: 0.38,
    solubility: 0.3, BBB: 0.1, hERG: 0.25, metabolic_stability: 0.5, lipinski_pass: false,
    composite_score: 71.5, cluster_id: 2, scaffold: 'quinoline', interactions_count: 6,
    admet: { absorption: 0.52, distribution: 0.4, metabolism: 0.55, excretion: 0.38, toxicity: 0.7, pharmacodynamics: 0.72 },
  }),
  mol(7, 'CN(C)C/C=C/C(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1OC1CCOC1', 'Afatinib', {
    bookmarked: true, source_run_id: 'run-a1',
    docking_score: -9.4, cnn_score: 0.89, cnn_affinity: 7.3, cnn_vs: 6.5,
    logP: 3.1, MW: 485.9, HBD: 2, HBA: 8, TPSA: 88.6, QED: 0.55,
    solubility: 0.52, BBB: 0.28, hERG: 0.14, metabolic_stability: 0.68, lipinski_pass: true,
    composite_score: 85.1, cluster_id: 1, scaffold: 'quinazoline', interactions_count: 7,
    admet: { absorption: 0.8, distribution: 0.62, metabolism: 0.7, excretion: 0.53, toxicity: 0.86, pharmacodynamics: 0.78 },
  }),
  mol(8, 'Fc1ccc(Nc2ncnc3cc4c(cc23)OCO4)cc1Cl', 'CHEMBL_501234', {
    source_run_id: 'run-a1',
    docking_score: -7.5, cnn_score: 0.62, cnn_affinity: 5.8, cnn_vs: 3.6,
    logP: 2.8, MW: 318.7, HBD: 2, HBA: 4, TPSA: 55.3, QED: 0.82,
    solubility: 0.78, BBB: 0.55, hERG: 0.06, metabolic_stability: 0.82, lipinski_pass: true,
    composite_score: 62.4, cluster_id: 1, scaffold: 'quinazoline', interactions_count: 3,
    admet: { absorption: 0.88, distribution: 0.72, metabolism: 0.8, excretion: 0.65, toxicity: 0.94, pharmacodynamics: 0.55 },
  }),
  mol(9, 'CCc1c(N)ncnc1Nc1ccc2c(cnn2C)c1', 'CHEMBL_602891', {
    source_run_id: 'run-a1',
    docking_score: -6.8, cnn_score: 0.55, cnn_affinity: 5.2, cnn_vs: 2.86,
    logP: 1.9, MW: 255.3, HBD: 3, HBA: 5, TPSA: 92.1, QED: 0.88,
    solubility: 0.85, BBB: 0.62, hERG: 0.04, metabolic_stability: 0.88, lipinski_pass: true,
    composite_score: 55.8, cluster_id: 4, scaffold: 'pyrimidine-pyrazole', interactions_count: 4,
    admet: { absorption: 0.9, distribution: 0.78, metabolism: 0.85, excretion: 0.72, toxicity: 0.96, pharmacodynamics: 0.48 },
  }),
  mol(10, 'O=C(Nc1cccc(C(F)(F)F)c1)c1cnc2ccc(-c3ccc(F)cc3)cc2n1', 'CHEMBL_789012', {
    source_run_id: 'run-a1',
    docking_score: -8.3, cnn_score: 0.78, cnn_affinity: 6.6, cnn_vs: 5.15,
    logP: 4.1, MW: 411.4, HBD: 1, HBA: 4, TPSA: 62.8, QED: 0.65,
    solubility: 0.4, BBB: 0.38, hERG: 0.19, metabolic_stability: 0.62, lipinski_pass: true,
    composite_score: 70.2, cluster_id: 5, scaffold: 'quinoxaline', interactions_count: 5,
    admet: { absorption: 0.72, distribution: 0.58, metabolism: 0.64, excretion: 0.45, toxicity: 0.78, pharmacodynamics: 0.68 },
  }),
  mol(11, 'Nc1ncnc2c1cnn2-c1ccc(O)c(O)c1', 'CHEMBL_345678', {
    source_run_id: 'run-a1',
    docking_score: -7.2, cnn_score: 0.58, cnn_affinity: 5.5, cnn_vs: 3.19,
    logP: 0.8, MW: 268.3, HBD: 4, HBA: 6, TPSA: 118.5, QED: 0.76,
    solubility: 0.92, BBB: 0.12, hERG: 0.03, metabolic_stability: 0.9, lipinski_pass: true,
    composite_score: 58.1, cluster_id: 4, scaffold: 'purine', interactions_count: 4,
    admet: { absorption: 0.65, distribution: 0.42, metabolism: 0.82, excretion: 0.78, toxicity: 0.97, pharmacodynamics: 0.5 },
  }),
  mol(12, 'COc1ccc(-c2cc3c(N)ncnc3[nH]2)cc1OC', 'CHEMBL_456789', {
    source_run_id: 'run-a1',
    docking_score: -7.8, cnn_score: 0.68, cnn_affinity: 6.1, cnn_vs: 4.15,
    logP: 2.1, MW: 296.3, HBD: 3, HBA: 5, TPSA: 85.4, QED: 0.81,
    solubility: 0.72, BBB: 0.48, hERG: 0.07, metabolic_stability: 0.78, lipinski_pass: true,
    composite_score: 66.7, cluster_id: 3, scaffold: 'benzimidazole', interactions_count: 5,
    admet: { absorption: 0.84, distribution: 0.65, metabolism: 0.76, excretion: 0.6, toxicity: 0.92, pharmacodynamics: 0.62 },
  }),
  // Molecules with ONLY import + docking (no ADMET yet) — simulates partial enrichment
  mol(13, 'CC(C)n1c(/C=C/c2cc(F)c(O)c(F)c2)nc2cnc(N)nc21', 'ZINC_001234', {
    source_run_id: 'run-a1',
    docking_score: -8.6, cnn_score: 0.8, cnn_affinity: 6.7, cnn_vs: 5.36,
    logP: null, MW: 358.4, HBD: null, HBA: null, TPSA: null, QED: null,
    composite_score: null, cluster_id: null, scaffold: null,
  }),
  mol(14, 'O=C(Nc1cccnc1)Nc1ccc(-c2ccnc3cc(Cl)ccc23)cc1', 'ZINC_005678', {
    source_run_id: 'run-a1',
    docking_score: -7.9, cnn_score: 0.71, cnn_affinity: 6.3, cnn_vs: 4.47,
    logP: null, MW: 388.9, HBD: null, HBA: null, TPSA: null, QED: null,
    composite_score: null, cluster_id: null, scaffold: null,
  }),
  // Molecules with ONLY import (no docking, no ADMET) — just imported
  mol(15, 'CC(=O)Nc1ccc(O)c(-c2nc3ccccc3[nH]2)c1', 'ZINC_009012', {
    source_run_id: 'run-a1',
    docking_score: null, cnn_score: null, cnn_affinity: null, cnn_vs: null,
  }),
  mol(16, 'Oc1ccc(-c2nc(-c3ccccc3)c(-c3ccccc3)[nH]2)cc1', 'ZINC_003456', {
    source_run_id: 'run-a1',
    docking_score: null, cnn_score: null, cnn_affinity: null, cnn_vs: null,
  }),
  mol(17, 'COc1cc(Nc2ncc3cc(-c4ccccc4F)c(=O)n(C)c3n2)cc(OC)c1OC', 'PubChem_12345', {
    source_run_id: 'run-a1',
    docking_score: -8.0, cnn_score: 0.73, cnn_affinity: 6.4, cnn_vs: 4.67,
    logP: 3.5, MW: 437.5, HBD: 1, HBA: 6, TPSA: 75.3, QED: 0.6,
    solubility: 0.45, BBB: 0.32, hERG: 0.16, metabolic_stability: 0.64, lipinski_pass: true,
    composite_score: 69.8, cluster_id: 5, scaffold: 'pyrimidine', interactions_count: 5,
    admet: { absorption: 0.74, distribution: 0.58, metabolism: 0.66, excretion: 0.48, toxicity: 0.82, pharmacodynamics: 0.64 },
  }),
  mol(18, 'CC(C)(O)c1ccc(-c2ccc3ncnc(Nc4ccc(F)c(Cl)c4)c3c2)cn1', 'PubChem_67890', {
    source_run_id: 'run-a1',
    docking_score: -8.7, cnn_score: 0.83, cnn_affinity: 6.8, cnn_vs: 5.64,
    logP: 3.0, MW: 412.9, HBD: 2, HBA: 5, TPSA: 72.1, QED: 0.67,
    solubility: 0.58, BBB: 0.4, hERG: 0.11, metabolic_stability: 0.71, lipinski_pass: true,
    composite_score: 76.4, cluster_id: 1, scaffold: 'quinazoline', interactions_count: 6,
    admet: { absorption: 0.8, distribution: 0.64, metabolism: 0.72, excretion: 0.55, toxicity: 0.88, pharmacodynamics: 0.7 },
  }),
  mol(19, 'O=c1[nH]c2cc(F)ccc2n1-c1ccc(-c2ccnc(N)n2)cc1', 'Enamine_T5612', {
    source_run_id: 'run-a1',
    docking_score: -7.4, cnn_score: 0.6, cnn_affinity: 5.6, cnn_vs: 3.36,
    logP: 2.4, MW: 336.3, HBD: 2, HBA: 5, TPSA: 88.2, QED: 0.78,
    solubility: 0.7, BBB: 0.35, hERG: 0.09, metabolic_stability: 0.76, lipinski_pass: true,
    composite_score: 61.2, cluster_id: 4, scaffold: 'benzimidazolone', interactions_count: 4,
    admet: { absorption: 0.82, distribution: 0.6, metabolism: 0.75, excretion: 0.62, toxicity: 0.9, pharmacodynamics: 0.55 },
  }),
  mol(20, 'CC1(C)CNc2cc3c(Nc4cccc(Cl)c4F)ncnc3cc21', 'Enamine_Z8901', {
    source_run_id: 'run-a1',
    docking_score: -8.9, cnn_score: 0.85, cnn_affinity: 7.0, cnn_vs: 5.95,
    logP: 3.6, MW: 368.8, HBD: 2, HBA: 4, TPSA: 56.8, QED: 0.74,
    solubility: 0.5, BBB: 0.45, hERG: 0.13, metabolic_stability: 0.69, lipinski_pass: true,
    composite_score: 77.9, cluster_id: 1, scaffold: 'quinazoline', interactions_count: 6,
    admet: { absorption: 0.78, distribution: 0.68, metabolism: 0.7, excretion: 0.52, toxicity: 0.86, pharmacodynamics: 0.72 },
  }),
]

// Phase B molecules — fewer, bookmarked hits from Phase A + generated variants
const PHASE_B_MOLECULES = [
  // Promoted from Phase A (bookmarked)
  { ...PHASE_A_MOLECULES[0], id: 'mol-b1', generation_level: 0, bookmarked: false },
  { ...PHASE_A_MOLECULES[4], id: 'mol-b2', generation_level: 0, bookmarked: false },
  { ...PHASE_A_MOLECULES[6], id: 'mol-b3', generation_level: 0, bookmarked: true },
  // Generated variants
  mol('b4', 'COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCNCC1', 'Gefitinib-v1', {
    generation_level: 1, parent_molecule_id: 'mol-b1', source_run_id: 'run-b2',
    docking_score: -9.6, cnn_score: 0.9, cnn_affinity: 7.5, cnn_vs: 6.75,
    logP: 2.8, MW: 431.9, HBD: 2, HBA: 7, TPSA: 78.2, QED: 0.74,
    composite_score: 86.4, cluster_id: 1, scaffold: 'quinazoline', bookmarked: true,
    admet: { absorption: 0.85, distribution: 0.68, metabolism: 0.74, excretion: 0.58, toxicity: 0.9, pharmacodynamics: 0.78 },
  }),
  mol('b5', 'C=CC(=O)Nc1cc(Nc2nccc(-c3cn(C)c4ccccc34)n2)c(OC)cc1NC', 'Osimertinib-v1', {
    generation_level: 1, parent_molecule_id: 'mol-b2', source_run_id: 'run-b2',
    docking_score: -10.1, cnn_score: 0.95, cnn_affinity: 8.0, cnn_vs: 7.6,
    logP: 3.0, MW: 470.6, HBD: 3, HBA: 6, TPSA: 95.2, QED: 0.56,
    composite_score: 91.3, cluster_id: 3, scaffold: 'pyrimidine-indole', bookmarked: true,
    admet: { absorption: 0.78, distribution: 0.72, metabolism: 0.76, excretion: 0.55, toxicity: 0.84, pharmacodynamics: 0.88 },
  }),
  mol('b6', 'CN(C)C/C=C/C(=O)Nc1cc2c(Nc3ccc(F)cc3)ncnc2cc1OC1CCOC1', 'Afatinib-v1', {
    generation_level: 1, parent_molecule_id: 'mol-b3', source_run_id: 'run-b2',
    docking_score: -9.0, cnn_score: 0.86, cnn_affinity: 7.1, cnn_vs: 6.11,
    logP: 2.9, MW: 476.5, HBD: 2, HBA: 8, TPSA: 88.6, QED: 0.58,
    composite_score: 83.7, cluster_id: 1, scaffold: 'quinazoline', bookmarked: false,
    admet: { absorption: 0.82, distribution: 0.64, metabolism: 0.71, excretion: 0.54, toxicity: 0.87, pharmacodynamics: 0.76 },
  }),
]

// ---- Runs ----

const RUNS_PHASE_A = [
  {
    id: 'run-a1', phase_id: 'phase-a', type: 'import', status: 'completed',
    config: { source: 'chembl', query: 'EGFR kinase inhibitors', max_molecules: 500 },
    started_at: '2026-02-20T10:15:00Z', completed_at: '2026-02-20T10:16:30Z',
    error_message: null, archived: false, created_at: '2026-02-20T10:15:00Z',
  },
  {
    id: 'run-a2', phase_id: 'phase-a', type: 'docking', status: 'completed',
    config: { engine: 'gnina_gpu', exhaustiveness: 32, num_modes: 9, seed: 0 },
    started_at: '2026-02-20T10:20:00Z', completed_at: '2026-02-20T10:23:45Z',
    error_message: null, archived: false, created_at: '2026-02-20T10:19:00Z',
  },
  {
    id: 'run-a3', phase_id: 'phase-a', type: 'admet', status: 'completed',
    config: { properties: ['logP', 'solubility', 'BBB', 'hERG', 'metabolic_stability'] },
    started_at: '2026-02-20T10:25:00Z', completed_at: '2026-02-20T10:25:30Z',
    error_message: null, archived: false, created_at: '2026-02-20T10:24:00Z',
  },
  {
    id: 'run-a4', phase_id: 'phase-a', type: 'scoring', status: 'completed',
    config: { weights: { docking_score: 0.3, CNNscore: 0.2, logP: 0.15, solubility: 0.1, selectivity: 0.15, novelty: 0.1 } },
    started_at: '2026-02-20T10:26:00Z', completed_at: '2026-02-20T10:26:15Z',
    error_message: null, archived: false, created_at: '2026-02-20T10:25:45Z',
  },
  {
    id: 'run-a5', phase_id: 'phase-a', type: 'enrichment', status: 'running',
    config: { analyses: ['prolif', 'clustering'] },
    started_at: '2026-02-20T10:30:00Z', completed_at: null,
    error_message: null, archived: false, created_at: '2026-02-20T10:29:00Z',
    progress: 65,
  },
]

const RUNS_PHASE_B = [
  {
    id: 'run-b1', phase_id: 'phase-b', type: 'import', status: 'completed',
    config: { source: 'phase_selection', source_phase_id: 'phase-a', selection: 'bookmarked' },
    started_at: '2026-02-21T09:00:00Z', completed_at: '2026-02-21T09:00:05Z',
    error_message: null, archived: false, created_at: '2026-02-21T09:00:00Z',
  },
  {
    id: 'run-b2', phase_id: 'phase-b', type: 'generation', status: 'completed',
    config: { method: 'scaffold_hopping', iterations: 3, variants_per_iteration: 5 },
    started_at: '2026-02-21T09:05:00Z', completed_at: '2026-02-21T09:12:30Z',
    error_message: null, archived: false, created_at: '2026-02-21T09:04:00Z',
  },
]

// ---- Phases ----

const PHASES = [
  {
    id: 'phase-a',
    campaign_id: 'camp-1',
    type: 'hit_discovery',
    label: 'Phase A',
    status: 'frozen',
    frozen_at: '2026-02-20T11:00:00Z',
    column_presets: ['name', 'docking_score', 'cnn_score', 'logP', 'MW', 'TPSA', 'lipinski_pass', 'composite_score'],
    created_at: '2026-02-20T10:10:00Z',
    runs: RUNS_PHASE_A,
    molecules: PHASE_A_MOLECULES,
    stats: { total_molecules: 20, bookmarked: 5, runs_completed: 4, runs_running: 1 },
  },
  {
    id: 'phase-b',
    campaign_id: 'camp-1',
    type: 'hit_to_lead',
    label: 'Phase B',
    status: 'active',
    frozen_at: null,
    column_presets: ['name', 'docking_score', 'cnn_score', 'composite_score', 'generation_level', 'cluster_id', 'scaffold'],
    created_at: '2026-02-21T09:00:00Z',
    runs: RUNS_PHASE_B,
    molecules: PHASE_B_MOLECULES,
    stats: { total_molecules: 6, bookmarked: 2, runs_completed: 2, runs_running: 0 },
  },
  {
    id: 'phase-c',
    campaign_id: 'camp-1',
    type: 'lead_optimization',
    label: 'Phase C',
    status: 'active',
    frozen_at: null,
    column_presets: ['name', 'composite_score', 'logP', 'solubility', 'BBB', 'hERG', 'metabolic_stability', 'interactions_count'],
    created_at: null, // Not created yet
    runs: [],
    molecules: [],
    stats: { total_molecules: 0, bookmarked: 0, runs_completed: 0, runs_running: 0 },
  },
]

// ---- Campaigns ----

const CAMPAIGNS = [
  {
    id: 'camp-1',
    project_id: 'proj-1',
    name: 'ATP Pocket — ChEMBL Screening',
    pocket_config: {
      center: [10.5, 20.3, 15.1],
      size: [22, 22, 22],
      residues: ['THR790', 'MET793', 'LYS745', 'ASP855', 'LEU718'],
      druggability: 0.92,
    },
    scoring_weights: { docking_score: 0.3, CNNscore: 0.2, logP: 0.15, solubility: 0.1, selectivity: 0.15, novelty: 0.1 },
    docking_defaults: { engine: 'gnina_gpu', exhaustiveness: 32, num_modes: 9, seed: 0 },
    rules: { lipinski: true, pains: true, max_MW: 600 },
    created_at: '2026-02-20T10:05:00Z',
    phases: PHASES,
  },
]

// ---- Projects ----

export const MOCK_PROJECTS = [
  {
    id: 'proj-1',
    name: 'EGFR Inhibitors',
    description: 'Targeting EGFR kinase domain (T790M mutant) for non-small cell lung cancer. Focus on third-generation covalent inhibitors.',
    target_name: 'EGFR',
    target_pdb_id: '1M17',
    uniprot_id: 'P00533',
    status: 'active',
    created_at: '2026-02-20T10:00:00Z',
    updated_at: '2026-02-21T09:15:00Z',
    campaigns: CAMPAIGNS,
  },
  {
    id: 'proj-2',
    name: 'BRAF V600E',
    description: 'Selective BRAF V600E inhibitors for melanoma. Avoiding paradoxical ERK activation.',
    target_name: 'BRAF',
    target_pdb_id: '4MNE',
    uniprot_id: 'P15056',
    status: 'active',
    created_at: '2026-02-18T14:00:00Z',
    updated_at: '2026-02-19T16:30:00Z',
    campaigns: [{
      id: 'camp-2', project_id: 'proj-2', name: 'DFG-out Pocket',
      pocket_config: { center: [5.2, 12.1, 8.7], size: [20, 20, 20], residues: ['V600E', 'C532', 'D594'] },
      scoring_weights: { docking_score: 0.3, CNNscore: 0.2, logP: 0.15, solubility: 0.1, selectivity: 0.15, novelty: 0.1 },
      created_at: '2026-02-18T14:05:00Z',
      phases: [{
        id: 'phase-braf-a', campaign_id: 'camp-2', type: 'hit_discovery', label: 'Phase A',
        status: 'active', frozen_at: null,
        column_presets: ['name', 'docking_score', 'cnn_score', 'logP', 'MW', 'composite_score'],
        created_at: '2026-02-18T14:10:00Z',
        runs: [
          { id: 'run-braf-1', phase_id: 'phase-braf-a', type: 'import', status: 'completed',
            config: { source: 'chembl', max_molecules: 200 },
            started_at: '2026-02-18T14:15:00Z', completed_at: '2026-02-18T14:15:20Z',
            created_at: '2026-02-18T14:15:00Z' },
        ],
        molecules: [
          mol('braf-1', 'Clc1ccc(-c2cc(C(F)(F)F)nn2-c2ccncc2)c(NS(=O)(=O)c2ccc(F)cc2)c1', 'Vemurafenib', {
            docking_score: -10.2, cnn_score: 0.92, cnn_affinity: 7.8, cnn_vs: 7.18,
            logP: 4.5, MW: 489.9, composite_score: 85.1,
          }),
        ],
        stats: { total_molecules: 1, bookmarked: 0, runs_completed: 1, runs_running: 0 },
      }],
    }],
  },
  {
    id: 'proj-3',
    name: 'JAK2 Selectivity',
    description: 'JAK2-selective inhibitors for myeloproliferative neoplasms. Must spare JAK1/JAK3.',
    target_name: 'JAK2',
    target_pdb_id: '4AQC',
    uniprot_id: 'O60674',
    status: 'active',
    created_at: '2026-02-15T09:00:00Z',
    updated_at: '2026-02-15T09:00:00Z',
    campaigns: [],
  },
]

// ---- Column definitions ----
// All possible columns a molecule can have, with metadata

export const ALL_COLUMNS = [
  { key: 'name', label: 'Name', type: 'text', width: 140, sortable: true },
  { key: 'smiles', label: 'SMILES', type: 'smiles', width: 200, sortable: false },
  { key: 'source_run_id', label: 'Source', type: 'text', width: 80, sortable: true },
  { key: 'docking_score', label: 'Docking', type: 'number', unit: 'kcal/mol', width: 90, sortable: true, colorScale: 'lower-better' },
  { key: 'cnn_score', label: 'CNN Score', type: 'number', width: 90, sortable: true, colorScale: 'higher-better' },
  { key: 'cnn_affinity', label: 'CNN Aff.', type: 'number', width: 80, sortable: true, colorScale: 'higher-better' },
  { key: 'cnn_vs', label: 'CNN VS', type: 'number', width: 80, sortable: true, colorScale: 'higher-better' },
  { key: 'composite_score', label: 'Composite', type: 'number', width: 90, sortable: true, colorScale: 'higher-better' },
  { key: 'logP', label: 'LogP', type: 'number', width: 70, sortable: true },
  { key: 'MW', label: 'MW', type: 'number', unit: 'Da', width: 70, sortable: true },
  { key: 'HBD', label: 'HBD', type: 'number', width: 55, sortable: true },
  { key: 'HBA', label: 'HBA', type: 'number', width: 55, sortable: true },
  { key: 'TPSA', label: 'TPSA', type: 'number', unit: 'A2', width: 70, sortable: true },
  { key: 'QED', label: 'QED', type: 'number', width: 60, sortable: true, colorScale: 'higher-better' },
  { key: 'lipinski_pass', label: 'Lipinski', type: 'boolean', width: 70, sortable: true },
  { key: 'solubility', label: 'Solubility', type: 'number', width: 80, sortable: true, colorScale: 'higher-better' },
  { key: 'BBB', label: 'BBB', type: 'number', width: 60, sortable: true },
  { key: 'hERG', label: 'hERG', type: 'number', width: 60, sortable: true, colorScale: 'lower-better' },
  { key: 'metabolic_stability', label: 'Met. Stab.', type: 'number', width: 80, sortable: true, colorScale: 'higher-better' },
  { key: 'cluster_id', label: 'Cluster', type: 'number', width: 65, sortable: true },
  { key: 'scaffold', label: 'Scaffold', type: 'text', width: 110, sortable: true },
  { key: 'interactions_count', label: 'Contacts', type: 'number', width: 75, sortable: true, colorScale: 'higher-better' },
  { key: 'generation_level', label: 'Gen. Level', type: 'number', width: 80, sortable: true },
  { key: 'parent_molecule_id', label: 'Parent', type: 'text', width: 80, sortable: false },
]

// Column presets per phase type (CDC section 4.4)
export const COLUMN_PRESETS = {
  hit_discovery: ['name', 'docking_score', 'cnn_score', 'logP', 'MW', 'HBD', 'HBA', 'TPSA', 'lipinski_pass', 'composite_score'],
  hit_to_lead: ['name', 'docking_score', 'cnn_score', 'composite_score', 'generation_level', 'cluster_id', 'scaffold'],
  lead_optimization: ['name', 'composite_score', 'logP', 'solubility', 'BBB', 'hERG', 'metabolic_stability', 'interactions_count'],
}

// ---- Run type definitions ----

export const RUN_TYPES = [
  { type: 'import', label: 'Import Molecules', icon: 'upload', description: 'Import from SDF, SMILES file, or internal selection', phases: ['hit_discovery', 'hit_to_lead', 'lead_optimization'] },
  { type: 'docking', label: 'Molecular Docking', icon: 'target', description: 'Dock molecules against the target protein', phases: ['hit_discovery', 'hit_to_lead', 'lead_optimization'] },
  { type: 'admet', label: 'ADMET Properties', icon: 'shield', description: 'Predict absorption, distribution, metabolism, excretion, toxicity', phases: ['hit_discovery', 'hit_to_lead', 'lead_optimization'] },
  { type: 'scoring', label: 'Composite Scoring', icon: 'star', description: 'Calculate weighted composite score', phases: ['hit_discovery', 'hit_to_lead', 'lead_optimization'] },
  { type: 'enrichment', label: 'Enrichment Analysis', icon: 'layers', description: 'ProLIF interactions, clustering, scaffold analysis', phases: ['hit_discovery', 'hit_to_lead', 'lead_optimization'] },
  { type: 'generation', label: 'De Novo Generation', icon: 'sparkles', description: 'Generate new molecules from selected hits', phases: ['hit_to_lead', 'lead_optimization'] },
  { type: 'clustering', label: 'Diversity Clustering', icon: 'grid', description: 'Cluster by scaffold, compute Tanimoto similarity', phases: ['hit_discovery', 'hit_to_lead', 'lead_optimization'] },
]

// ---- Phase type definitions ----

export const PHASE_TYPES = {
  hit_discovery: { label: 'Hit Discovery', short: 'A', color: 'blue', description: 'Find hits in a compound library' },
  hit_to_lead: { label: 'Hit-to-Lead', short: 'B', color: 'purple', description: 'Optimize hits, generate analogues' },
  lead_optimization: { label: 'Lead Optimization', short: 'C', color: 'green', description: 'Refine leads with full ADMET profiling' },
}

// ---- Rich molecule detail data (for drawer / detail panels) ----
// Keyed by molecule id — lookup at render time

export const MOLECULE_DETAILS = {
  'mol-1': {
    interactions: {
      residues: [
        { name: 'MET793', type: 'HBond', distance: 2.1, chain: 'A' },
        { name: 'THR790', type: 'HBond', distance: 2.8, chain: 'A' },
        { name: 'LYS745', type: 'Salt Bridge', distance: 3.2, chain: 'A' },
        { name: 'ASP855', type: 'HBond', distance: 2.5, chain: 'A' },
        { name: 'LEU718', type: 'Hydrophobic', distance: 3.8, chain: 'A' },
        { name: 'VAL726', type: 'Hydrophobic', distance: 4.0, chain: 'A' },
      ],
      total_contacts: 12, hbonds: 3, hydrophobic: 5, ionic: 1, pi_stacking: 2, water_bridges: 1,
    },
    synthesis: {
      feasibility: 0.82, total_cost: '$180', num_steps: 3,
      steps: [
        { product: 'Quinazoline core', reagent: '2-amino-4-chlorobenzoic acid', conditions: 'DMF, 120C, 4h', yield: 0.85 },
        { product: 'N-aryl intermediate', reagent: '3-Chloro-4-fluoroaniline', conditions: 'Pd(PPh3)4, K2CO3, dioxane', yield: 0.72 },
        { product: 'Gefitinib', reagent: 'Morpholine side chain', conditions: 'DIPEA, MeCN, RT, 12h', yield: 0.88 },
      ],
    },
    safety: {
      pains_alerts: [], pains_pass: true,
      brenk_alerts: [],
      off_target: [
        { target: 'HER2 (ERBB2)', similarity: 0.82, risk: 'medium', family: 'Kinase' },
        { target: 'HER3 (ERBB3)', similarity: 0.65, risk: 'low', family: 'Kinase' },
      ],
      herg_risk: 0.12, herg_pass: true,
      ames_risk: 0.05, hepatotox_risk: 0.08,
      confidence: { overall: 0.85, binding: 0.92, admet: 0.78, selectivity: 0.82, safety: 0.88 },
    },
  },
  'mol-2': {
    interactions: {
      residues: [
        { name: 'MET793', type: 'HBond', distance: 2.3, chain: 'A' },
        { name: 'LYS745', type: 'Salt Bridge', distance: 3.1, chain: 'A' },
        { name: 'LEU788', type: 'Hydrophobic', distance: 3.6, chain: 'A' },
        { name: 'ALA743', type: 'Hydrophobic', distance: 3.9, chain: 'A' },
        { name: 'GLY796', type: 'HBond', distance: 2.7, chain: 'A' },
      ],
      total_contacts: 10, hbonds: 2, hydrophobic: 4, ionic: 1, pi_stacking: 2, water_bridges: 1,
    },
    synthesis: {
      feasibility: 0.91, total_cost: '$120', num_steps: 2,
      steps: [
        { product: 'Quinazoline core', reagent: '2-aminobenzaldehyde', conditions: 'AcOH, reflux, 2h', yield: 0.9 },
        { product: 'Erlotinib', reagent: 'Alkoxy substitution', conditions: 'K2CO3, DMF, 80C', yield: 0.85 },
      ],
    },
    safety: {
      pains_alerts: [], pains_pass: true, brenk_alerts: [],
      off_target: [
        { target: 'HER2 (ERBB2)', similarity: 0.71, risk: 'low', family: 'Kinase' },
      ],
      herg_risk: 0.08, herg_pass: true, ames_risk: 0.03, hepatotox_risk: 0.05,
      confidence: { overall: 0.88, binding: 0.89, admet: 0.85, selectivity: 0.90, safety: 0.92 },
    },
  },
  'mol-5': {
    interactions: {
      residues: [
        { name: 'CYS797', type: 'Covalent', distance: 1.8, chain: 'A' },
        { name: 'MET793', type: 'HBond', distance: 2.0, chain: 'A' },
        { name: 'THR790', type: 'HBond', distance: 2.6, chain: 'A' },
        { name: 'LYS745', type: 'Salt Bridge', distance: 3.0, chain: 'A' },
        { name: 'ASP855', type: 'HBond', distance: 2.4, chain: 'A' },
        { name: 'LEU718', type: 'Hydrophobic', distance: 3.5, chain: 'A' },
        { name: 'VAL726', type: 'Hydrophobic', distance: 3.7, chain: 'A' },
      ],
      total_contacts: 15, hbonds: 3, hydrophobic: 4, ionic: 1, pi_stacking: 3, water_bridges: 2, covalent: 1,
    },
    synthesis: {
      feasibility: 0.68, total_cost: '$450', num_steps: 5,
      steps: [
        { product: 'Pyrimidine core', reagent: 'Malononitrile + guanidine', conditions: 'EtOH, reflux', yield: 0.78 },
        { product: 'Indole coupling', reagent: '1-methylindole-3-boronic acid', conditions: 'Suzuki, Pd, 80C', yield: 0.7 },
        { product: 'Aniline coupling', reagent: '4-methoxy-3-aminoaniline', conditions: 'Buchwald-Hartwig', yield: 0.65 },
        { product: 'Side chain', reagent: 'Dimethylaminoethyl', conditions: 'Reductive amination', yield: 0.82 },
        { product: 'Osimertinib', reagent: 'Acryloyl chloride', conditions: 'Et3N, DCM, 0C', yield: 0.75 },
      ],
    },
    safety: {
      pains_alerts: [], pains_pass: true, brenk_alerts: ['Michael acceptor (acrylamide warhead — intended covalent)'],
      off_target: [
        { target: 'HER2 (ERBB2)', similarity: 0.55, risk: 'low', family: 'Kinase' },
        { target: 'EGFR WT', similarity: 0.3, risk: 'low', family: 'Kinase' },
        { target: 'BTK', similarity: 0.2, risk: 'low', family: 'Kinase' },
      ],
      herg_risk: 0.18, herg_pass: true, ames_risk: 0.07, hepatotox_risk: 0.12,
      confidence: { overall: 0.91, binding: 0.96, admet: 0.82, selectivity: 0.88, safety: 0.85 },
    },
  },
  'mol-7': {
    interactions: {
      residues: [
        { name: 'CYS797', type: 'Covalent', distance: 1.9, chain: 'A' },
        { name: 'MET793', type: 'HBond', distance: 2.2, chain: 'A' },
        { name: 'LYS745', type: 'Salt Bridge', distance: 3.3, chain: 'A' },
        { name: 'THR854', type: 'HBond', distance: 2.9, chain: 'A' },
        { name: 'PHE856', type: 'Pi-Stacking', distance: 3.8, chain: 'A' },
        { name: 'LEU792', type: 'Hydrophobic', distance: 4.1, chain: 'A' },
      ],
      total_contacts: 13, hbonds: 2, hydrophobic: 4, ionic: 1, pi_stacking: 2, water_bridges: 1, covalent: 1,
    },
    synthesis: {
      feasibility: 0.75, total_cost: '$280', num_steps: 4,
      steps: [
        { product: 'Quinazoline core', reagent: '2-amino-4,5-dimethoxybenzoic acid', conditions: 'Formamide, 180C', yield: 0.8 },
        { product: 'Aniline coupling', reagent: '3-chloro-4-fluoroaniline', conditions: 'DIPEA, iPrOH, reflux', yield: 0.75 },
        { product: 'Ether linkage', reagent: 'Tetrahydrofuryl alcohol', conditions: 'NaH, DMF', yield: 0.7 },
        { product: 'Afatinib', reagent: 'Dimethylaminocrotonoyl chloride', conditions: 'Et3N, DCM', yield: 0.72 },
      ],
    },
    safety: {
      pains_alerts: [], pains_pass: true, brenk_alerts: ['Michael acceptor (intended covalent)'],
      off_target: [
        { target: 'HER2 (ERBB2)', similarity: 0.78, risk: 'medium', family: 'Kinase' },
        { target: 'HER4 (ERBB4)', similarity: 0.62, risk: 'low', family: 'Kinase' },
      ],
      herg_risk: 0.14, herg_pass: true, ames_risk: 0.06, hepatotox_risk: 0.09,
      confidence: { overall: 0.87, binding: 0.93, admet: 0.8, selectivity: 0.78, safety: 0.86 },
    },
  },
}

// ---- Agent mock recommendations (campaign-level) ----

export const MOCK_AGENT_INSIGHTS = [
  {
    id: 'insight-1', type: 'recommendation', priority: 'high', timestamp: '2026-02-21T09:30:00Z',
    title: 'Run ADMET on Phase B molecules',
    body: 'Phase B has 6 molecules from scaffold hopping but no ADMET data yet. The generated variants Osimertinib-v1 and Gefitinib-v1 show excellent docking scores (CNN VS > 6.5) and should be prioritized for ADMET profiling before proceeding to Lead Optimization.',
    action: { label: 'Launch ADMET Run', phaseId: 'phase-b', runType: 'admet' },
  },
  {
    id: 'insight-2', type: 'observation', priority: 'medium', timestamp: '2026-02-20T11:05:00Z',
    title: 'Quinazoline scaffold dominates Phase A hits',
    body: 'Among the 5 bookmarked hits in Phase A, 4 share a quinazoline core scaffold. Consider diversifying with pyrimidine-indole or benzimidazole scaffolds from cluster 3 and 4 to reduce attrition risk.',
  },
  {
    id: 'insight-3', type: 'alert', priority: 'high', timestamp: '2026-02-21T09:15:00Z',
    title: 'Lapatinib fails Lipinski Rule of 5',
    body: 'Lapatinib (MW=581, logP=4.6) violates Lipinski rules. While bookmarked as a hit, consider deprioritizing or generating lower-MW analogues in Phase B.',
  },
  {
    id: 'insight-4', type: 'observation', priority: 'low', timestamp: '2026-02-20T10:45:00Z',
    title: 'Enrichment factor exceeds random baseline',
    body: 'Phase A screening achieved EF=2.1x at top-10% cutoff (CNN VS ranking), confirming the docking protocol correctly prioritizes known EGFR actives over decoys.',
  },
]

// ---- Target assessment data ----

export const MOCK_TARGET_ASSESSMENT = {
  'proj-1': {
    druggability_score: 0.92,
    evidence_score: 0.95,
    confidence: 0.88,
    disease_context: 'Non-small cell lung cancer (NSCLC)',
    known_drugs: ['Gefitinib (Iressa)', 'Erlotinib (Tarceva)', 'Osimertinib (Tagrisso)', 'Afatinib (Gilotrif)'],
    protein_family: 'Receptor Tyrosine Kinase (ErbB family)',
    gene_name: 'EGFR',
    organism: 'Homo sapiens',
    sequence_length: 1210,
    pdb_structures: 48,
    chembl_compounds: 15234,
    binding_site: {
      residues: ['THR790', 'MET793', 'LYS745', 'ASP855', 'LEU718', 'VAL726', 'ALA743', 'GLY796', 'CYS797'],
      volume: 482.5,
      druggability_probability: 0.94,
    },
  },
}

// ---- Activity log ----

export const MOCK_ACTIVITY = [
  { id: 'act-1', type: 'run_completed', message: 'Enrichment analysis running', project_id: 'proj-1', timestamp: '2026-02-20T10:30:00Z' },
  { id: 'act-2', type: 'phase_frozen', message: 'Phase A frozen with 5 bookmarked hits', project_id: 'proj-1', timestamp: '2026-02-20T11:00:00Z' },
  { id: 'act-3', type: 'run_completed', message: 'Generation completed — 3 new variants', project_id: 'proj-1', timestamp: '2026-02-21T09:12:30Z' },
  { id: 'act-4', type: 'agent_insight', message: 'Agent recommended ADMET run for Phase B', project_id: 'proj-1', timestamp: '2026-02-21T09:30:00Z' },
]
