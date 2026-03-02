// ---------------------------------------------------------------------------
// BindX V9 — Column definitions & phase metadata
// Defines the molecule table column structure for PhaseDashboard.
// CDC §4.1: "All columns with results are shown by default"
// ---------------------------------------------------------------------------

export const PHASE_TYPES = {
  hit_discovery: { label: 'Hit Discovery', short: 'A', color: 'blue', description: 'Find hits in a compound library' },
  hit_to_lead: { label: 'Hit-to-Lead', short: 'B', color: 'purple', description: 'Optimize hits, generate analogues' },
  lead_optimization: { label: 'Lead Optimization', short: 'C', color: 'green', description: 'Refine leads with full ADMET profiling' },
}

// All possible columns in the molecule table.
// key = flat field name on the molecule row (after property flattening)
export const ALL_COLUMNS = [
  // Identity (always available from import)
  { key: 'name', label: 'Name', type: 'text', group: 'identity', width: 140, sortable: true },
  { key: 'smiles', label: 'SMILES', type: 'smiles', group: 'identity', width: 200, sortable: false },
  { key: 'source_run_id', label: 'Source', type: 'text', group: 'identity', width: 80, sortable: true },
  // Docking scores
  { key: 'docking_score', label: 'Docking', type: 'number', group: 'docking', unit: 'kcal/mol', width: 90, sortable: true, colorScale: 'lower-better' },
  { key: 'cnn_score', label: 'CNN Score', type: 'number', group: 'docking', width: 90, sortable: true, colorScale: 'higher-better' },
  { key: 'cnn_affinity', label: 'CNN Aff.', type: 'number', group: 'docking', width: 80, sortable: true, colorScale: 'higher-better' },
  { key: 'cnn_vs', label: 'CNN VS', type: 'number', group: 'docking', width: 80, sortable: true, colorScale: 'higher-better' },
  // Drug properties / ADMET
  { key: 'logP', label: 'LogP', type: 'number', group: 'admet', width: 70, sortable: true },
  { key: 'MW', label: 'MW', type: 'number', group: 'admet', unit: 'Da', width: 70, sortable: true },
  { key: 'HBD', label: 'HBD', type: 'number', group: 'admet', width: 55, sortable: true },
  { key: 'HBA', label: 'HBA', type: 'number', group: 'admet', width: 55, sortable: true },
  { key: 'TPSA', label: 'TPSA', type: 'number', group: 'admet', unit: 'A2', width: 70, sortable: true },
  { key: 'QED', label: 'QED', type: 'number', group: 'admet', width: 60, sortable: true, colorScale: 'higher-better' },
  { key: 'lipinski_pass', label: 'Lipinski', type: 'boolean', group: 'admet', width: 70, sortable: true },
  { key: 'solubility', label: 'Solubility', type: 'number', group: 'admet', width: 80, sortable: true, colorScale: 'higher-better' },
  { key: 'BBB', label: 'BBB', type: 'number', group: 'admet', width: 60, sortable: true },
  { key: 'hERG', label: 'hERG', type: 'number', group: 'admet', width: 60, sortable: true, colorScale: 'lower-better' },
  { key: 'metabolic_stability', label: 'Met. Stab.', type: 'number', group: 'admet', width: 80, sortable: true, colorScale: 'higher-better' },
  { key: 'oral_bioavailability', label: 'Oral Bioavail.', type: 'number', group: 'admet', width: 95, sortable: true, colorScale: 'higher-better' },
  { key: 'plasma_protein_binding', label: 'PPB', type: 'number', group: 'admet', width: 60, sortable: true },
  // Scoring
  { key: 'composite_score', label: 'Composite', type: 'number', group: 'scoring', width: 90, sortable: true, colorScale: 'higher-better' },
  // Enrichment
  { key: 'interactions_count', label: 'Contacts', type: 'number', group: 'enrichment', width: 75, sortable: true, colorScale: 'higher-better' },
  { key: 'scaffold', label: 'Scaffold', type: 'text', group: 'enrichment', width: 110, sortable: true },
  // Clustering
  { key: 'cluster_id', label: 'Cluster', type: 'number', group: 'clustering', width: 65, sortable: true },
  { key: 'scaffold_smiles', label: 'Scaffold SMILES', type: 'smiles', group: 'clustering', width: 150, sortable: false },
  { key: 'tanimoto_to_centroid', label: 'Tanimoto', type: 'number', group: 'clustering', width: 80, sortable: true, colorScale: 'higher-better' },
  { key: 'is_representative', label: 'Representative', type: 'boolean', group: 'clustering', width: 90, sortable: true },
  // Off-target
  { key: 'selectivity_score', label: 'Selectivity', type: 'number', group: 'off_target', width: 85, sortable: true, colorScale: 'higher-better' },
  { key: 'off_target_hits', label: 'Off-targets', type: 'number', group: 'off_target', width: 80, sortable: true, colorScale: 'lower-better' },
  { key: 'selectivity_ratio', label: 'Select. Ratio', type: 'number', group: 'off_target', width: 90, sortable: true, colorScale: 'higher-better' },
  // Confidence
  { key: 'confidence_score', label: 'Confidence', type: 'number', group: 'confidence', width: 85, sortable: true, colorScale: 'higher-better', popup: 'confidence' },
  { key: 'pains_alert', label: 'PAINS', type: 'boolean', group: 'confidence', width: 60, sortable: true },
  { key: 'applicability_domain', label: 'Appl. Domain', type: 'boolean', group: 'confidence', width: 90, sortable: true },
  // Retrosynthesis
  { key: 'n_synth_steps', label: 'Synth Steps', type: 'number', group: 'retrosynthesis', width: 85, sortable: true, colorScale: 'lower-better' },
  { key: 'synth_confidence', label: 'Synth Conf.', type: 'number', group: 'retrosynthesis', width: 85, sortable: true, colorScale: 'higher-better', popup: 'retrosynthesis' },
  { key: 'synth_cost_estimate', label: 'Synth Cost', type: 'text', group: 'retrosynthesis', width: 80, sortable: true },
  { key: 'reagents_available', label: 'Reagents', type: 'boolean', group: 'retrosynthesis', width: 70, sortable: true },
  // Safety
  { key: 'herg_risk', label: 'hERG Risk', type: 'number', group: 'safety', width: 75, sortable: true, colorScale: 'lower-better' },
  { key: 'ames_mutagenicity', label: 'AMES', type: 'boolean', group: 'safety', width: 60, sortable: true },
  { key: 'hepatotoxicity', label: 'Hepatotox', type: 'number', group: 'safety', width: 80, sortable: true, colorScale: 'lower-better' },
  { key: 'skin_sensitization', label: 'Skin Sens.', type: 'boolean', group: 'safety', width: 80, sortable: true },
  { key: 'carcinogenicity', label: 'Carcino.', type: 'number', group: 'safety', width: 70, sortable: true, colorScale: 'lower-better' },
  { key: 'safety_color_code', label: 'Safety', type: 'text', group: 'safety', width: 65, sortable: true, popup: 'safety' },
  // Generation
  { key: 'generation_level', label: 'Gen. Level', type: 'number', group: 'generation', width: 80, sortable: true },
  { key: 'parent_molecule_id', label: 'Parent', type: 'text', group: 'generation', width: 80, sortable: false },
  // User annotations
  { key: '_detail', label: '', type: 'action', group: 'annotation', width: 36, sortable: false },
  { key: 'tags', label: 'Tags', type: 'tags', group: 'annotation', width: 140, sortable: false },
  { key: 'invalidated', label: 'Invalid', type: 'invalidation', group: 'annotation', width: 60, sortable: true },
  { key: 'user_comment', label: 'Notes', type: 'editable_text', group: 'annotation', width: 160, sortable: false },
  { key: 'ai_comment', label: 'AI Notes', type: 'text', group: 'annotation', width: 180, sortable: false },
]

// Column key lookup for fast access
const COLUMN_MAP = Object.fromEntries(ALL_COLUMNS.map(c => [c.key, c]))

// Column presets per phase type
export const COLUMN_PRESETS = {
  hit_discovery: ['name', 'docking_score', 'cnn_score', 'logP', 'MW', 'HBD', 'HBA', 'TPSA', 'lipinski_pass', 'composite_score', 'confidence_score', 'safety_color_code'],
  hit_to_lead: ['name', 'docking_score', 'cnn_score', 'composite_score', 'generation_level', 'cluster_id', 'scaffold', 'synth_confidence', 'safety_color_code'],
  lead_optimization: ['name', 'composite_score', 'logP', 'solubility', 'BBB', 'hERG', 'metabolic_stability', 'interactions_count', 'selectivity_score', 'synth_confidence', 'safety_color_code'],
}

// Backend → frontend key mapping (backend names that differ from column keys)
const PROP_ALIASES = {
  hbd: 'HBD', hba: 'HBA', qed: 'QED', tpsa: 'TPSA',
  bbb_permeability: 'BBB', herg_inhibition: 'hERG', herg_risk: 'herg_risk',
  color_code: 'safety_color_code',
}

/**
 * Recursively flatten a nested object into a single-level dict.
 * Skips arrays and special keys (flags, note, status, smiles, confidence_modifier).
 */
function deepFlatten(obj, out = {}) {
  const SKIP_KEYS = new Set(['flags', 'note', 'status', 'smiles', 'confidence_modifier', 'nearest_tanimoto', 'docking_status', 'tree', 'children', 'reaction', 'reactants', 'reactant_names', 'conditions'])
  for (const [k, v] of Object.entries(obj)) {
    if (SKIP_KEYS.has(k)) continue
    if (v && typeof v === 'object' && !Array.isArray(v)) {
      deepFlatten(v, out)
    } else {
      const alias = PROP_ALIASES[k] || k
      out[alias] = v
    }
  }
  return out
}

/**
 * Flatten a molecule's `properties` dict from the API into a flat row.
 * Handles deeply nested formats like:
 *   { properties: { admet: { toxicity: { hepatotoxicity: 0.07 } }, physicochemical: { MW: 84 } } }
 * Output shape: { hepatotoxicity: 0.07, MW: 84, ... }
 */
export function flattenMoleculeProperties(mol) {
  if (!mol) return mol
  const flat = { ...mol }
  if (mol.properties && typeof mol.properties === 'object') {
    Object.assign(flat, deepFlatten(mol.properties))
  }
  return flat
}

/**
 * Detect which columns actually have data across a set of molecules.
 * Returns the keys that appear in at least one molecule with a non-null value.
 * Identity columns (name, smiles) are always included.
 */
export function detectAvailableColumns(molecules) {
  const alwaysShow = new Set(['name'])
  const available = new Set(alwaysShow)

  for (const mol of molecules) {
    for (const col of ALL_COLUMNS) {
      if (available.has(col.key)) continue
      const val = mol[col.key]
      if (val !== undefined && val !== null && val !== '') {
        available.add(col.key)
      }
    }
  }
  return ALL_COLUMNS.filter(c => available.has(c.key))
}

/**
 * Get column definition by key. Returns undefined for unknown keys.
 */
export function getColumnDef(key) {
  return COLUMN_MAP[key]
}

// ---------------------------------------------------------------------------
// Run type definitions (CDC §3.4)
// ---------------------------------------------------------------------------

export const RUN_TYPES = [
  { type: 'import', label: 'Import Molecules', icon: 'upload', description: 'Import from SDF, SMILES file, or internal selection', phases: ['hit_discovery', 'hit_to_lead', 'lead_optimization'] },
  { type: 'calculation', label: 'Run Calculations', icon: 'calculator', description: 'Select one or more calculation types to run on selected molecules', phases: ['hit_discovery', 'hit_to_lead', 'lead_optimization'] },
  { type: 'generation', label: 'De Novo Generation', icon: 'sparkles', description: 'Generate new molecules from selected hits', phases: ['hit_to_lead', 'lead_optimization'] },
]

export const CALCULATION_SUBTYPES = [
  { key: 'docking', label: 'Molecular Docking', icon: 'target', description: 'Dock molecules against the target protein', columns: ['docking_score', 'cnn_score', 'cnn_affinity', 'poses'] },
  { key: 'admet', label: 'ADMET Properties', icon: 'shield', description: 'Predict absorption, distribution, metabolism, excretion, toxicity', columns: ['logP', 'solubility', 'BBB', 'hERG', 'metabolic_stability', 'oral_bioavailability', 'plasma_protein_binding'] },
  { key: 'scoring', label: 'Composite Scoring', icon: 'star', description: 'Calculate weighted composite score', columns: ['composite_score'] },
  { key: 'enrichment', label: 'Enrichment Analysis', icon: 'layers', description: 'ProLIF interactions, scaffold analysis', columns: ['interactions_count', 'scaffold'] },
  { key: 'clustering', label: 'Diversity Clustering', icon: 'grid', description: 'Cluster by scaffold, compute Tanimoto similarity', columns: ['cluster_id', 'scaffold_smiles', 'tanimoto_to_centroid'] },
  { key: 'off_target', label: 'Off-target Selectivity', icon: 'crosshair', description: 'Assess selectivity against off-target proteins', columns: ['selectivity_score', 'off_target_hits', 'selectivity_ratio'] },
  { key: 'confidence', label: 'Confidence Analysis', icon: 'check-circle', description: 'PAINS filters, applicability domain, convergence', columns: ['confidence_score', 'pains_alert', 'applicability_domain', 'confidence_flags'] },
  { key: 'retrosynthesis', label: 'Retrosynthesis', icon: 'git-branch', description: 'Synthesis feasibility, cost estimation, reagent availability', columns: ['n_synth_steps', 'synth_confidence', 'synth_cost_estimate', 'reagents_available'] },
  { key: 'safety', label: 'Safety Profile', icon: 'alert-triangle', description: 'Full safety: hERG, AMES, hepatotoxicity, carcinogenicity', columns: ['herg_risk', 'ames_mutagenicity', 'hepatotoxicity', 'skin_sensitization', 'carcinogenicity', 'safety_color_code'] },
]
