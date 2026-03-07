// ---------------------------------------------------------------------------
// BindX V9 — Column definitions & phase metadata
// Defines the molecule table column structure for PhaseDashboard.
// CDC §4.1: "All columns with results are shown by default"
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Group metadata — labels & colors for the grouped header row
// ---------------------------------------------------------------------------
export const GROUP_META = {
  molecule:          { label: 'Molecule',             bg: 'bg-slate-50',    text: 'text-slate-500',    border: 'border-slate-300' },
  annotation:        { label: 'Annotations',          bg: 'bg-amber-50',   text: 'text-amber-600',    border: 'border-amber-300' },
  docking:           { label: 'Molecular Docking',    bg: 'bg-blue-50',    text: 'text-blue-600',     border: 'border-blue-400' },
  admet:             { label: 'ADMET Properties',     bg: 'bg-emerald-50', text: 'text-emerald-600',  border: 'border-emerald-400' },
  scoring:           { label: 'Molecular Properties',  bg: 'bg-violet-50',  text: 'text-violet-600',   border: 'border-violet-400' },
  enrichment:        { label: 'Enrichment Analysis',  bg: 'bg-cyan-50',    text: 'text-cyan-600',     border: 'border-cyan-400' },
  clustering:        { label: 'Diversity Clustering', bg: 'bg-teal-50',    text: 'text-teal-600',     border: 'border-teal-400' },
  off_target:        { label: 'Off-target Selectivity', bg: 'bg-orange-50', text: 'text-orange-600',  border: 'border-orange-400' },
  confidence:        { label: 'Confidence Analysis',  bg: 'bg-indigo-50',  text: 'text-indigo-600',   border: 'border-indigo-400' },
  retrosynthesis:    { label: 'Retrosynthesis',       bg: 'bg-pink-50',    text: 'text-pink-600',     border: 'border-pink-400' },
  activity_cliffs:   { label: 'Activity Cliffs',      bg: 'bg-yellow-50',  text: 'text-yellow-600',   border: 'border-yellow-400' },
  pharmacophore:     { label: 'Pharmacophore Mapping', bg: 'bg-fuchsia-50', text: 'text-fuchsia-600', border: 'border-fuchsia-400' },
  generation:        { label: 'De Novo Generation',   bg: 'bg-sky-50',     text: 'text-sky-600',      border: 'border-sky-400' },
}

export const PHASE_TYPES = {
  hit_discovery: { label: 'Hit Discovery', short: 'A', color: 'blue', description: 'Find hits in a compound library' },
  hit_to_lead: { label: 'Hit-to-Lead', short: 'B', color: 'purple', description: 'Optimize hits, generate analogues' },
  lead_optimization: { label: 'Lead Optimization', short: 'C', color: 'green', description: 'Refine leads with full ADMET profiling' },
}

// All possible columns in the molecule table.
// key = flat field name on the molecule row (after property flattening)
export const ALL_COLUMNS = [
  // Molecule identity (always available from import)
  { key: 'name', label: 'Name', type: 'text', group: 'molecule', width: 140, sortable: true },
  { key: 'smiles', label: 'SMILES', type: 'smiles', group: 'molecule', width: 200, sortable: false },
  { key: 'source_run_id', label: 'Source', type: 'source', group: 'molecule', width: 110, sortable: true },
  // User annotations
  { key: 'tags', label: 'Tags', type: 'tags', group: 'annotation', width: 160, sortable: false, hiddenByDefault: true },
  { key: 'user_comment', label: 'Notes', type: 'editable_text', group: 'annotation', width: 180, sortable: false, hiddenByDefault: true },
  { key: 'invalidated', label: 'Invalid', type: 'invalidation', group: 'annotation', width: 65, sortable: true, hiddenByDefault: true },
  { key: 'ai_comment', label: 'AI Notes', type: 'text', group: 'annotation', width: 180, sortable: false, hiddenByDefault: true },
  // Docking scores
  { key: 'docking_score', label: 'Docking', type: 'number', group: 'docking', unit: 'kcal/mol', width: 90, sortable: true, colorScale: 'lower-better' },
  { key: 'cnn_score', label: 'CNN Score', type: 'number', group: 'docking', width: 90, sortable: true, colorScale: 'higher-better' },
  { key: 'cnn_affinity', label: 'CNN Aff.', type: 'number', group: 'docking', width: 80, sortable: true, colorScale: 'higher-better' },
  { key: 'cnn_vs', label: 'CNN VS', type: 'number', group: 'docking', width: 80, sortable: true, colorScale: 'higher-better' },
  { key: 'consensus_ecr', label: 'ECR', type: 'number', group: 'docking', width: 65, sortable: true, colorScale: 'higher-better' },
  // Molecular Properties (from scoring run: physicochemical descriptors, drug-likeness, SA score)
  { key: 'logP', label: 'LogP', type: 'number', group: 'scoring', width: 70, sortable: true },
  { key: 'MW', label: 'MW', type: 'number', group: 'scoring', unit: 'Da', width: 70, sortable: true },
  { key: 'HBD', label: 'HBD', type: 'number', group: 'scoring', width: 55, sortable: true },
  { key: 'HBA', label: 'HBA', type: 'number', group: 'scoring', width: 55, sortable: true },
  { key: 'TPSA', label: 'TPSA', type: 'number', group: 'scoring', unit: 'A2', width: 70, sortable: true },
  { key: 'QED', label: 'QED', type: 'number', group: 'scoring', width: 60, sortable: true, colorScale: 'higher-better' },
  { key: 'lipinski_pass', label: 'Lipinski', type: 'boolean', group: 'scoring', width: 70, sortable: true },
  { key: 'heavy_atom_count', label: 'HA', type: 'number', group: 'scoring', width: 45, sortable: true },
  { key: 'sa_score', label: 'SA Score', type: 'number', group: 'scoring', width: 70, sortable: true, colorScale: 'lower-better' },
  { key: 'ro3_pass', label: 'Ro3', type: 'boolean', group: 'scoring', width: 50, sortable: true },
  { key: 'inchikey', label: 'InChIKey', type: 'text', group: 'scoring', width: 180, sortable: true },
  { key: 'ligand_efficiency', label: 'LE', type: 'number', group: 'docking', width: 55, sortable: true, colorScale: 'higher-better' },
  // ADMET Properties (from admet run: ADME + toxicity + druglikeness rules + ADMET score)
  { key: 'solubility', label: 'Solubility', type: 'number', group: 'admet', width: 80, sortable: true, colorScale: 'higher-better' },
  { key: 'BBB', label: 'BBB', type: 'number', group: 'admet', width: 60, sortable: true },
  { key: 'hERG', label: 'hERG', type: 'number', group: 'admet', width: 60, sortable: true, colorScale: 'lower-better' },
  { key: 'metabolic_stability', label: 'Met. Stab.', type: 'number', group: 'admet', width: 80, sortable: true, colorScale: 'higher-better' },
  { key: 'oral_bioavailability', label: 'Oral Bioavail.', type: 'number', group: 'admet', width: 95, sortable: true, colorScale: 'higher-better' },
  { key: 'plasma_protein_binding', label: 'PPB', type: 'number', group: 'admet', width: 60, sortable: true },
  { key: 'cns_mpo', label: 'CNS MPO', type: 'number', group: 'admet', width: 70, sortable: true, colorScale: 'higher-better' },
  { key: 'pfizer_alert', label: 'Pfizer 3/75', type: 'boolean', group: 'admet', width: 80, sortable: true },
  { key: 'gsk_alert', label: 'GSK 4/400', type: 'boolean', group: 'admet', width: 75, sortable: true },
  { key: 'brenk_alert', label: 'Brenk', type: 'boolean', group: 'admet', width: 60, sortable: true },
  { key: 'ames_mutagenicity', label: 'AMES', type: 'boolean', group: 'admet', width: 60, sortable: true },
  { key: 'hepatotoxicity', label: 'Hepatotox', type: 'number', group: 'admet', width: 80, sortable: true, colorScale: 'lower-better' },
  { key: 'skin_sensitization', label: 'Skin Sens.', type: 'boolean', group: 'admet', width: 80, sortable: true },
  { key: 'carcinogenicity', label: 'Carcino.', type: 'number', group: 'admet', width: 70, sortable: true, colorScale: 'lower-better' },
  { key: 'safety_color_code', label: 'Safety', type: 'text', group: 'admet', width: 65, sortable: true, popup: 'safety' },
  { key: 'composite_score', label: 'ADMET Score', type: 'number', group: 'admet', width: 95, sortable: true, colorScale: 'higher-better' },
  // Enrichment Analysis
  { key: 'interactions_count', label: 'Contacts', type: 'number', group: 'enrichment', width: 75, sortable: true, colorScale: 'higher-better' },
  { key: 'scaffold', label: 'Scaffold', type: 'text', group: 'enrichment', width: 110, sortable: true },
  // Diversity Clustering
  { key: 'cluster_id', label: 'Cluster', type: 'number', group: 'clustering', width: 65, sortable: true },
  { key: 'scaffold_smiles', label: 'Scaffold SMILES', type: 'smiles', group: 'clustering', width: 150, sortable: false },
  { key: 'tanimoto_to_centroid', label: 'Tanimoto', type: 'number', group: 'clustering', width: 80, sortable: true, colorScale: 'higher-better' },
  { key: 'is_representative', label: 'Representative', type: 'boolean', group: 'clustering', width: 90, sortable: true },
  // Off-target Selectivity
  { key: 'selectivity_score', label: 'Selectivity', type: 'number', group: 'off_target', width: 85, sortable: true, colorScale: 'higher-better' },
  { key: 'off_target_hits', label: 'Off-targets', type: 'number', group: 'off_target', width: 80, sortable: true, colorScale: 'lower-better' },
  { key: 'selectivity_ratio', label: 'Select. Ratio', type: 'number', group: 'off_target', width: 90, sortable: true, colorScale: 'higher-better' },
  // Confidence Analysis
  { key: 'confidence_score', label: 'Confidence', type: 'number', group: 'confidence', width: 85, sortable: true, colorScale: 'higher-better', popup: 'confidence' },
  { key: 'pains_alert', label: 'PAINS', type: 'boolean', group: 'confidence', width: 60, sortable: true },
  { key: 'applicability_domain', label: 'Appl. Domain', type: 'boolean', group: 'confidence', width: 90, sortable: true },
  // Retrosynthesis
  { key: 'n_synth_steps', label: 'Synth Steps', type: 'number', group: 'retrosynthesis', width: 85, sortable: true, colorScale: 'lower-better' },
  { key: 'synth_confidence', label: 'Synth Conf.', type: 'number', group: 'retrosynthesis', width: 85, sortable: true, colorScale: 'higher-better', popup: 'retrosynthesis' },
  { key: 'synth_cost_estimate', label: 'Synth Cost', type: 'text', group: 'retrosynthesis', width: 80, sortable: true },
  { key: 'reagents_available', label: 'Reagents', type: 'boolean', group: 'retrosynthesis', width: 70, sortable: true },
  // Activity Cliffs
  { key: 'is_cliff', label: 'Cliff', type: 'boolean', group: 'activity_cliffs', width: 55, sortable: true },
  { key: 'sali_max', label: 'SALI', type: 'number', group: 'activity_cliffs', width: 65, sortable: true, colorScale: 'higher-better' },
  { key: 'n_cliffs', label: '# Cliffs', type: 'number', group: 'activity_cliffs', width: 65, sortable: true },
  // Pharmacophore
  { key: 'pharmacophore_features', label: 'Pharma. Features', type: 'number', group: 'pharmacophore', width: 95, sortable: true, colorScale: 'higher-better' },
  { key: 'pharmacophore_similarity', label: 'Pharma. Sim.', type: 'number', group: 'pharmacophore', width: 85, sortable: true, colorScale: 'higher-better' },
  // Generation
  { key: 'generation_level', label: 'Gen. Level', type: 'number', group: 'generation', width: 80, sortable: true },
  { key: 'parent_molecule_id', label: 'Parent', type: 'text', group: 'generation', width: 80, sortable: false },
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
  bbb_permeability: 'BBB', herg_inhibition: 'hERG',
  color_code: 'safety_color_code',
  affinity: 'docking_score', vina_score: 'docking_score',
}

/**
 * Recursively flatten a nested object into a single-level dict.
 * Skips arrays and special keys (flags, note, status, smiles, confidence_modifier).
 */
function deepFlatten(obj, out = {}) {
  const SKIP_KEYS = new Set(['flags', 'note', 'status', 'smiles', 'confidence_modifier', 'nearest_tanimoto', 'docking_status', 'tree', 'children', 'reaction', 'reactants', 'reactant_names', 'conditions', 'results', 'warnings', 'pose_molblock'])
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
      if (available.has(col.key) || col.hiddenByDefault) continue
      const val = mol[col.key]
      if (val !== undefined && val !== null && val !== '') {
        available.add(col.key)
      }
    }
  }
  return ALL_COLUMNS.filter(c => available.has(c.key))
}

/** Annotation columns — always available in ColumnSelector, hidden by default */
export const ANNOTATION_COLUMNS = ALL_COLUMNS.filter(c => c.hiddenByDefault)

/**
 * Get column definition by key. Returns undefined for unknown keys.
 */
export function getColumnDef(key) {
  return COLUMN_MAP[key]
}

/**
 * Get effective color scale config for a column, considering user overrides.
 * @param {string} key — column key
 * @param {Object} overrides — columnColorOverrides from settingsStore
 * @returns {{ mode: string, intervals?: Array }} — resolved color config
 *   mode: 'none' | 'higher-better' | 'lower-better' | 'custom'
 *   intervals (custom only): [{ min: number|null, max: number|null, color: string }]
 */
export function getEffectiveColorScale(key, overrides = {}) {
  if (overrides[key]) return overrides[key]
  const col = COLUMN_MAP[key]
  if (!col || !col.colorScale) return { mode: 'none' }
  return { mode: col.colorScale }
}

// ---------------------------------------------------------------------------
// Run type definitions (CDC §3.4)
// ---------------------------------------------------------------------------

export const RUN_TYPES = [
  { type: 'import', label: 'Import Molecules', icon: 'upload', description: 'Import from SDF, SMILES file, or internal selection', phases: ['hit_discovery', 'hit_to_lead', 'lead_optimization'] },
  { type: 'calculation', label: 'Run Calculation', icon: 'calculator', description: 'Select a calculation type to run on selected molecules', phases: ['hit_discovery', 'hit_to_lead', 'lead_optimization'] },
  { type: 'generation', label: 'De Novo Generation', icon: 'sparkles', description: 'Generate new molecules from selected hits', phases: ['hit_to_lead', 'lead_optimization'] },
]

// Labels come from GROUP_META — single source of truth for naming
export const CALCULATION_SUBTYPES = [
  { key: 'docking', icon: 'target', description: 'Dock molecules against the target protein', columns: ['docking_score', 'cnn_score', 'cnn_affinity', 'ligand_efficiency', 'cnn_vs', 'consensus_ecr'] },
  { key: 'admet', icon: 'shield', description: 'Predict absorption, distribution, metabolism, excretion, toxicity', columns: ['solubility', 'BBB', 'hERG', 'metabolic_stability', 'oral_bioavailability', 'plasma_protein_binding', 'cns_mpo', 'pfizer_alert', 'gsk_alert', 'brenk_alert', 'ames_mutagenicity', 'hepatotoxicity', 'skin_sensitization', 'carcinogenicity', 'safety_color_code', 'composite_score'] },
  { key: 'scoring', icon: 'star', description: 'Physicochemical descriptors (MW, LogP, TPSA), drug-likeness rules (Lipinski, QED, Ro3), SA score', columns: ['logP', 'MW', 'HBD', 'HBA', 'TPSA', 'QED', 'lipinski_pass', 'heavy_atom_count', 'sa_score', 'ro3_pass', 'inchikey'] },
  { key: 'enrichment', icon: 'layers', description: 'ProLIF interactions, scaffold analysis', columns: ['interactions_count', 'scaffold'] },
  { key: 'clustering', icon: 'grid', description: 'Cluster by scaffold, compute Tanimoto similarity', columns: ['cluster_id', 'scaffold_smiles', 'tanimoto_to_centroid'] },
  { key: 'off_target', icon: 'crosshair', description: 'Assess selectivity against off-target proteins', columns: ['selectivity_score', 'off_target_hits', 'selectivity_ratio'] },
  { key: 'confidence', icon: 'check-circle', description: 'PAINS filters, applicability domain, convergence', columns: ['confidence_score', 'pains_alert', 'applicability_domain', 'confidence_flags'] },
  { key: 'retrosynthesis', icon: 'git-branch', description: 'Synthesis feasibility, cost estimation, reagent availability', columns: ['n_synth_steps', 'synth_confidence', 'synth_cost_estimate', 'reagents_available'] },
  { key: 'pharmacophore', icon: 'hexagon', description: 'Map 3D pharmacophoric features and compute pairwise similarity', columns: [] },
  { key: 'activity_cliffs', icon: 'trending-up', description: 'Detect structure-activity cliffs — similar molecules with large activity differences', columns: ['is_cliff', 'sali_max', 'n_cliffs'] },
].map(st => ({ ...st, label: GROUP_META[st.key].label }))

// Estimated run times (displayed in RunCreator confirmation)
export const ESTIMATED_TIMES = {
  import: '~5 seconds', calculation: '~1-5 minutes', generation: '~5-10 minutes',
  docking: '~3-5 min', admet: '~30 sec', scoring: '~15 sec',
  enrichment: '~1-2 min', clustering: '~30 sec', off_target: '~1 min',
  confidence: '~15 sec', retrosynthesis: '~1 min',
  pharmacophore: '~30 sec', activity_cliffs: '~15 sec',
}

// ---------------------------------------------------------------------------
// Column defaults for RunCreator checklist (Step 3)
// Columns NOT listed here: checked by default, no dependency.
// ---------------------------------------------------------------------------
export const COLUMN_DEFAULTS = {
  // Scoring — optional
  heavy_atom_count:  { checked: false },
  ro3_pass:          { checked: false },
  inchikey:          { checked: false },
  ligand_efficiency: { checked: true },
  // ADMET — piggyback optionals
  cns_mpo:      { checked: false },
  pfizer_alert: { checked: false },
  gsk_alert:    { checked: false },
  brenk_alert:  { checked: false },
  // Docking — key metrics (GNINA only)
  cnn_score:     { checked: true, requiresEngine: 'gnina' },
  cnn_affinity:  { checked: true, requiresEngine: 'gnina' },
  cnn_vs:        { checked: true, requiresEngine: 'gnina' },
  consensus_ecr: { checked: true, requiresEngine: 'gnina' },
}

/**
 * Check if a column is available given the current run config.
 * Columns with requiresEngine are only available when the engine matches.
 */
export function isColumnAvailable(colKey, config = {}) {
  const def = COLUMN_DEFAULTS[colKey]
  if (!def?.requiresEngine) return true
  const engine = config?.docking?.engine || config?.engine || ''
  // 'gnina' matches gnina_gpu and gnina_cpu
  if (def.requiresEngine === 'gnina') {
    return engine.startsWith('gnina')
  }
  return engine === def.requiresEngine
}

/**
 * Get the list of available columns for a calculation subtype given the config.
 */
export function getAvailableColumns(calcKey, config = {}) {
  const sub = CALCULATION_SUBTYPES.find(s => s.key === calcKey)
  if (!sub) return []
  return sub.columns.filter(k => isColumnAvailable(k, config))
}

/**
 * Get default checked column keys for a given calculation subtype.
 */
export function getDefaultCheckedColumns(calcKey, config = {}) {
  const available = getAvailableColumns(calcKey, config)
  return available.filter(k => COLUMN_DEFAULTS[k]?.checked ?? true)
}
