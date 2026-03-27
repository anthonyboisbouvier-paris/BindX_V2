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
  adme:              { label: 'ADME',                  bg: 'bg-emerald-50', text: 'text-emerald-600',  border: 'border-emerald-400' },
  toxicity:          { label: 'Toxicity',              bg: 'bg-red-50',    text: 'text-red-600',      border: 'border-red-400' },
  scoring:           { label: 'Molecular Properties',  bg: 'bg-violet-50',  text: 'text-violet-600',   border: 'border-violet-400' },
  interactions:      { label: 'Interactions',          bg: 'bg-lime-50',    text: 'text-lime-600',     border: 'border-lime-400' },
  scaffold:          { label: 'Structural Analysis',  bg: 'bg-emerald-50', text: 'text-emerald-600',  border: 'border-emerald-400' },
  clustering:        { label: 'Diversity Clustering', bg: 'bg-teal-50',    text: 'text-teal-600',     border: 'border-teal-400' },
  off_target:        { label: 'Off-target Selectivity', bg: 'bg-orange-50', text: 'text-orange-600',  border: 'border-orange-400' },
  confidence:        { label: 'Confidence Analysis',  bg: 'bg-indigo-50',  text: 'text-indigo-600',   border: 'border-indigo-400' },
  retrosynthesis:    { label: 'Retrosynthesis',       bg: 'bg-pink-50',    text: 'text-pink-600',     border: 'border-pink-400' },
  pharmacophore:     { label: 'Reference Similarity',  bg: 'bg-fuchsia-50', text: 'text-fuchsia-600', border: 'border-fuchsia-400' },
  composite:         { label: 'Composite Score',      bg: 'bg-amber-50',   text: 'text-amber-600',    border: 'border-amber-400' },
  afvs:              { label: 'AFVS Screening',       bg: 'bg-cyan-50',    text: 'text-cyan-600',     border: 'border-cyan-400' },
  generation:        { label: 'Gen.',                  bg: 'bg-sky-50',     text: 'text-sky-600',      border: 'border-sky-400' },
  import_data:       { label: 'Import Data',          bg: 'bg-rose-50',   text: 'text-rose-600',     border: 'border-rose-400' },
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
  { key: 'name', label: 'Name', type: 'text', group: 'molecule', width: 120, sortable: true },
  { key: 'smiles', label: 'SMILES', type: 'smiles', group: 'molecule', width: 160, sortable: false },
  { key: 'source_run_id', label: 'Source', type: 'source', group: 'molecule', width: 90, sortable: true, sortKey: 'source' },
  // User annotations
  { key: 'tags', label: 'Tags', type: 'tags', group: 'annotation', width: 160, sortable: false, hiddenByDefault: true },
  { key: 'user_comment', label: 'Notes', type: 'editable_text', group: 'annotation', width: 180, sortable: false, hiddenByDefault: true },
  { key: 'invalidated', label: 'Invalid', type: 'invalidation', group: 'annotation', width: 65, sortable: true, hiddenByDefault: true },
  { key: 'ai_comment', label: 'AI Notes', type: 'text', group: 'annotation', width: 180, sortable: false, hiddenByDefault: true },
  // Docking scores
  { key: 'docking_score', label: 'Docking', type: 'number', group: 'docking', unit: 'kcal/mol', width: 75, sortable: true, colorScale: 'lower-better' },
  { key: 'cnn_score', label: 'CNN Score', type: 'number', group: 'docking', width: 70, sortable: true, colorScale: 'higher-better' },
  { key: 'cnn_affinity', label: 'CNN Aff.', type: 'number', group: 'docking', width: 65, sortable: true, colorScale: 'higher-better' },
  { key: 'cnn_vs', label: 'CNN VS', type: 'number', group: 'docking', width: 65, sortable: true, colorScale: 'higher-better' },
  { key: 'consensus_ecr', label: 'ECR', type: 'number', group: 'docking', width: 55, sortable: true, colorScale: 'higher-better' },
  { key: 'pocket_distance', label: 'Pocket Dist.', type: 'number', group: 'docking', unit: 'Å', width: 72, sortable: true, colorScale: 'lower-better' },
  // Molecular Properties (from scoring run: physicochemical descriptors, drug-likeness, SA score)
  { key: 'logP', label: 'LogP', type: 'number', group: 'scoring', width: 55, sortable: true },
  { key: 'MW', label: 'MW', type: 'number', group: 'scoring', unit: 'Da', width: 60, sortable: true },
  { key: 'HBD', label: 'HBD', type: 'number', group: 'scoring', width: 45, sortable: true },
  { key: 'HBA', label: 'HBA', type: 'number', group: 'scoring', width: 45, sortable: true },
  { key: 'TPSA', label: 'TPSA', type: 'number', group: 'scoring', unit: 'A2', width: 55, sortable: true },
  { key: 'QED', label: 'QED', type: 'number', group: 'scoring', width: 50, sortable: true, colorScale: 'higher-better' },
  { key: 'lipinski_pass', label: 'Lipinski', type: 'boolean', group: 'scoring', width: 60, sortable: true },
  { key: 'heavy_atom_count', label: 'HA', type: 'number', group: 'scoring', width: 40, sortable: true },
  { key: 'sa_score', label: 'SA Score', type: 'number', group: 'scoring', width: 62, sortable: true, colorScale: 'lower-better' },
  { key: 'ro3_pass', label: 'Ro3', type: 'boolean', group: 'scoring', width: 45, sortable: true },
  { key: 'inchikey', label: 'InChIKey', type: 'text', group: 'scoring', width: 160, sortable: true },
  { key: 'n_hydrophobic', label: 'Hydroph.', type: 'number', group: 'scoring', width: 55, sortable: true },
  { key: 'n_aromatic', label: 'Arom.', type: 'number', group: 'scoring', width: 48, sortable: true },
  { key: 'n_pos_ionizable', label: 'Pos. Ion.', type: 'number', group: 'scoring', width: 58, sortable: true },
  { key: 'n_neg_ionizable', label: 'Neg. Ion.', type: 'number', group: 'scoring', width: 58, sortable: true },
  { key: 'ligand_efficiency', label: 'LE', type: 'number', group: 'docking', width: 45, sortable: true, colorScale: 'higher-better' },
  // ADME Properties (from adme run: absorption, distribution, metabolism, excretion + druglikeness rules)
  { key: 'solubility', label: 'logS', type: 'number', group: 'adme', width: 50, sortable: true, colorScale: 'higher-better' },
  { key: 'BBB', label: 'BBB', type: 'number', group: 'adme', width: 50, sortable: true },
  { key: 'oral_bioavailability', label: 'Oral Bio.', type: 'number', group: 'adme', width: 65, sortable: true, colorScale: 'higher-better' },
  { key: 'plasma_protein_binding', label: 'PPB', type: 'number', group: 'adme', width: 50, sortable: true },
  { key: 'half_life', label: 'Half-life', type: 'number', group: 'adme', unit: 'h', width: 62, sortable: true },
  { key: 'cyp_inhibitions', label: 'CYP Inh.', type: 'number', group: 'adme', width: 60, sortable: true, colorScale: 'lower-better', computed: true },
  { key: 'cns_mpo', label: 'CNS MPO', type: 'number', group: 'adme', width: 62, sortable: true, colorScale: 'higher-better' },
  { key: 'pfizer_alert', label: 'Pfizer', type: 'boolean', group: 'adme', width: 55, sortable: true, invertBoolean: true },
  { key: 'gsk_alert', label: 'GSK', type: 'boolean', group: 'adme', width: 50, sortable: true, invertBoolean: true },
  { key: 'brenk_alert', label: 'Brenk', type: 'boolean', group: 'adme', width: 50, sortable: true, invertBoolean: true },
  { key: 'pains_alert', label: 'PAINS', type: 'boolean', group: 'adme', width: 50, sortable: true, invertBoolean: true },
  { key: 'applicability_domain', label: 'Appl. Dom.', type: 'boolean', group: 'adme', width: 65, sortable: true },
  // Toxicity Properties (from toxicity run: safety endpoints + composite score)
  { key: 'hERG', label: 'hERG', type: 'number', group: 'toxicity', width: 50, sortable: true, colorScale: 'lower-better' },
  { key: 'herg_risk_level', label: 'hERG Risk', type: 'text', group: 'toxicity', width: 70, sortable: true },
  { key: 'herg_ic50', label: 'hERG IC50', type: 'number', group: 'toxicity', width: 70, sortable: true, colorScale: 'higher-better', unit: 'µM' },
  { key: 'ames_mutagenicity', label: 'AMES', type: 'number', group: 'toxicity', width: 50, sortable: true, colorScale: 'lower-better' },
  { key: 'hepatotoxicity', label: 'Hepato.', type: 'number', group: 'toxicity', width: 58, sortable: true, colorScale: 'lower-better' },
  { key: 'skin_sensitization', label: 'Skin', type: 'number', group: 'toxicity', width: 50, sortable: true, colorScale: 'lower-better' },
  { key: 'carcinogenicity', label: 'Carcino.', type: 'number', group: 'toxicity', width: 60, sortable: true, colorScale: 'lower-better' },
  { key: 'safety_color_code', label: 'Safety', type: 'text', group: 'toxicity', width: 55, sortable: true, popup: 'safety' },
  { key: 'composite_score', label: 'Score', type: 'number', group: 'toxicity', width: 55, sortable: true, colorScale: 'higher-better' },
  // Interactions (protein-ligand contacts)
  { key: 'n_interactions', label: 'Contacts', type: 'number', group: 'interactions', width: 62, sortable: true, colorScale: 'higher-better' },
  { key: 'functional_contacts', label: 'Func.', type: 'number', group: 'interactions', width: 50, sortable: true, colorScale: 'higher-better' },
  { key: 'interaction_quality', label: 'Quality', type: 'number', group: 'interactions', width: 58, sortable: true, colorScale: 'higher-better' },
  { key: 'total_functional', label: 'Tot. Func.', type: 'number', group: 'interactions', width: 65, sortable: true },
  { key: 'key_hbonds', label: 'Key HB', type: 'number', group: 'interactions', width: 55, sortable: true, colorScale: 'higher-better' },
  // Scaffold (Murcko decomposition + BRICS)
  { key: 'scaffold_smiles', label: 'Scaffold', type: 'smiles', group: 'scaffold', width: 160, sortable: true },
  { key: 'n_modifiable_positions', label: 'R-Pos.', type: 'number', group: 'scaffold', width: 55, sortable: true },
  { key: 'brics_bond_count', label: 'BRICS', type: 'number', group: 'scaffold', width: 55, sortable: true },
  { key: 'scaffold_n_rings', label: 'Rings', type: 'number', group: 'scaffold', width: 48, sortable: true },
  { key: 'scaffold_group', label: 'Series', type: 'number', group: 'scaffold', width: 55, sortable: true },
  { key: 'scaffold_group_size', label: 'Series Size', type: 'number', group: 'scaffold', width: 70, sortable: true, colorScale: 'higher-better' },
  // Diversity Clustering (merged into scaffold / Structural Analysis)
  { key: 'cluster_id', label: 'Cluster', type: 'number', group: 'scaffold', width: 55, sortable: true },
  { key: 'cluster_size', label: 'Size', type: 'number', group: 'scaffold', width: 45, sortable: true },
  { key: 'tanimoto_to_centroid', label: 'Tanimoto', type: 'number', group: 'scaffold', width: 65, sortable: true, colorScale: 'higher-better' },
  { key: 'is_representative', label: 'Repr.', type: 'boolean', group: 'scaffold', width: 50, sortable: true },
  // Off-target Selectivity
  { key: 'selectivity_score', label: 'Select.', type: 'number', group: 'off_target', width: 55, sortable: true, colorScale: 'higher-better' },
  { key: 'off_target_hits', label: 'Off-tgt', type: 'number', group: 'off_target', width: 55, sortable: true, colorScale: 'lower-better' },
  { key: 'selectivity_ratio', label: 'Sel. Ratio', type: 'number', group: 'off_target', width: 60, sortable: true, colorScale: 'higher-better' },
  // Confidence Analysis (legacy — kept for backward compat with old runs)
  { key: 'confidence_score', label: 'Confid.', type: 'number', group: 'confidence', width: 58, sortable: true, colorScale: 'higher-better', popup: 'confidence' },
  // Retrosynthesis
  { key: 'n_synth_steps', label: 'Steps', type: 'number', group: 'retrosynthesis', width: 50, sortable: true, colorScale: 'lower-better' },
  { key: 'synth_confidence', label: 'Conf.', type: 'number', group: 'retrosynthesis', width: 50, sortable: true, colorScale: 'higher-better', popup: 'retrosynthesis' },
  { key: 'synth_cost_estimate', label: 'Cost', type: 'text', group: 'retrosynthesis', width: 55, sortable: true },
  { key: 'reagents_available', label: 'Reag.', type: 'boolean', group: 'retrosynthesis', width: 50, sortable: true },
  // Composite Score (weighted multi-criteria)
  { key: 'weighted_score', label: 'Score', type: 'number', group: 'composite', width: 55, sortable: true, colorScale: 'higher-better' },
  // Reference Similarity (pharmacophoric fingerprint scoring — up to 5 references)
  { key: 'pharmacophore_fit_1', label: 'Ref. 1', type: 'number', group: 'pharmacophore', width: 60, sortable: true, colorScale: 'higher-better' },
  { key: 'pharmacophore_fit_2', label: 'Ref. 2', type: 'number', group: 'pharmacophore', width: 60, sortable: true, colorScale: 'higher-better' },
  { key: 'pharmacophore_fit_3', label: 'Ref. 3', type: 'number', group: 'pharmacophore', width: 60, sortable: true, colorScale: 'higher-better' },
  { key: 'pharmacophore_fit_4', label: 'Ref. 4', type: 'number', group: 'pharmacophore', width: 60, sortable: true, colorScale: 'higher-better' },
  { key: 'pharmacophore_fit_5', label: 'Ref. 5', type: 'number', group: 'pharmacophore', width: 60, sortable: true, colorScale: 'higher-better' },
  { key: 'pharmacophore_fit_best', label: 'Best Score', type: 'number', group: 'pharmacophore', width: 70, sortable: true, colorScale: 'higher-better' },
  { key: 'pharmacophore_fit_avg', label: 'Avg Score', type: 'number', group: 'pharmacophore', width: 70, sortable: true, colorScale: 'higher-better' },
  // Generation
  { key: 'generation_level', label: 'Gen.', type: 'number', group: 'generation', width: 50, sortable: true },
  { key: 'parent_molecule_id', label: 'Parent', type: 'text', group: 'generation', width: 70, sortable: false },
  // Import Data (metadata from source databases)
  { key: 'pchembl_value', label: 'pChEMBL', type: 'number', group: 'import_data', width: 65, sortable: true, colorScale: 'higher-better' },
  { key: 'activity_value_nM', label: 'Act. (nM)', type: 'number', group: 'import_data', width: 70, sortable: true, colorScale: 'lower-better' },
  { key: 'activity_type', label: 'Act. Type', type: 'text', group: 'import_data', width: 70, sortable: true },
  { key: 'assay_name', label: 'Assay', type: 'text', group: 'import_data', width: 140, sortable: true },
  { key: 'import_mwt', label: 'MW (src)', type: 'number', group: 'import_data', width: 65, sortable: true },
  { key: 'import_logp', label: 'LogP (src)', type: 'number', group: 'import_data', width: 65, sortable: true },
  { key: 'import_hbd', label: 'HBD (src)', type: 'number', group: 'import_data', width: 60, sortable: true },
  { key: 'import_hba', label: 'HBA (src)', type: 'number', group: 'import_data', width: 60, sortable: true },
  { key: 'import_tpsa', label: 'TPSA (src)', type: 'number', group: 'import_data', width: 65, sortable: true },
]

// Column key lookup for fast access
const COLUMN_MAP = Object.fromEntries(ALL_COLUMNS.map(c => [c.key, c]))

// Column presets per phase type
export const COLUMN_PRESETS = {
  hit_discovery: ['name', 'docking_score', 'cnn_score', 'logP', 'MW', 'HBD', 'HBA', 'TPSA', 'lipinski_pass', 'composite_score', 'safety_color_code'],
  hit_to_lead: ['name', 'docking_score', 'cnn_score', 'composite_score', 'generation_level', 'cluster_id', 'scaffold_smiles', 'synth_confidence', 'safety_color_code'],
  lead_optimization: ['name', 'composite_score', 'logP', 'solubility', 'BBB', 'hERG', 'oral_bioavailability', 'n_interactions', 'selectivity_score', 'synth_confidence', 'safety_color_code'],
}

// Backend → frontend key mapping (backend names that differ from column keys)
const PROP_ALIASES = {
  hbd: 'HBD', hba: 'HBA', qed: 'QED', tpsa: 'TPSA',
  bbb_permeability: 'BBB', herg_inhibition: 'hERG',
  color_code: 'safety_color_code',
  affinity: 'docking_score', vina_score: 'docking_score',
  mwt: 'import_mwt', mw: 'import_mwt',
  logp: 'import_logp',
}

/**
 * Recursively flatten a nested object into a single-level dict.
 * Skips arrays and special keys (flags, note, status, smiles, confidence_modifier).
 */
function deepFlatten(obj, out = {}) {
  const SKIP_KEYS = new Set(['flags', 'note', 'status', 'smiles', 'confidence_modifier', 'nearest_tanimoto', 'docking_status', 'tree', 'children', 'reaction', 'reactants', 'reactant_names', 'conditions', 'results', 'warnings', 'pose_molblock', 'preparation', 'breakdown', 'weights_used', 'scaffold_positions', 'scaffold_svg', 'interactions_detail', 'interactions_method', 'mapped_functional', 'herg_specialized'])
  for (const [k, v] of Object.entries(obj)) {
    if (SKIP_KEYS.has(k)) continue
    if (Array.isArray(v)) continue // skip arrays — they can't render as table cells
    if (v && typeof v === 'object') {
      deepFlatten(v, out)
    } else {
      const alias = PROP_ALIASES[k] || k
      // Parse string numbers from EAV storage (e.g. "6.5" → 6.5)
      if (typeof v === 'string' && v !== '' && !isNaN(v)) {
        out[alias] = parseFloat(v)
      } else {
        out[alias] = v
      }
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
    // Preserve preparation metadata as object (skipped by deepFlatten to avoid column pollution)
    if (mol.properties.preparation && typeof mol.properties.preparation === 'object') {
      flat.preparation = mol.properties.preparation
    }
    // Preserve composite breakdown/weights for detail panel
    const compositeData = mol.properties.composite
    if (compositeData && typeof compositeData === 'object') {
      if (compositeData.breakdown) flat.breakdown = compositeData.breakdown
      if (compositeData.weights_used) flat.weights_used = compositeData.weights_used
    }
    // Preserve interactions detail for detail panel
    const interData = mol.properties.interactions
    if (interData && typeof interData === 'object') {
      if (interData.interactions_detail) flat.interactions_detail = interData.interactions_detail
      if (interData.interactions_method) flat.interactions_method = interData.interactions_method
      if (interData.mapped_functional) flat.mapped_functional = interData.mapped_functional
    }
    // Preserve scaffold metadata for detail panel
    const scaffoldData = mol.properties.scaffold
    if (scaffoldData && typeof scaffoldData === 'object') {
      if (scaffoldData.scaffold_positions) flat.scaffold_positions = scaffoldData.scaffold_positions
      if (scaffoldData.scaffold_svg) flat.scaffold_svg = scaffoldData.scaffold_svg
    }
    // Preserve pharmacophore metadata for detail panel + extract indexed fits
    const pharmaData = mol.properties.pharmacophore
    if (pharmaData && typeof pharmaData === 'object') {
      if (pharmaData.matched_types) flat.matched_types = pharmaData.matched_types
      if (pharmaData.missing_types) flat.missing_types = pharmaData.missing_types
      if (pharmaData.feature_counts) flat.feature_counts = pharmaData.feature_counts
      if (pharmaData.ref_feature_counts) flat.ref_feature_counts = pharmaData.ref_feature_counts
      // Multi-reference indexed fits (pharmacophore_fit_1..5 + pharmacophore_fit_best)
      for (let i = 1; i <= 5; i++) {
        const k = `pharmacophore_fit_${i}`
        if (pharmaData[k] != null) flat[k] = pharmaData[k]
      }
      if (pharmaData.pharmacophore_fit_best != null) flat.pharmacophore_fit_best = pharmaData.pharmacophore_fit_best
      if (pharmaData.pharmacophore_fit_avg != null) flat.pharmacophore_fit_avg = pharmaData.pharmacophore_fit_avg
    }
    // Extract hERG specialized data (IC50 + risk level) into flat keys
    const hergSpec = mol.properties.toxicity?.herg_specialized || mol.properties.herg_specialized
    if (hergSpec && typeof hergSpec === 'object') {
      flat.herg_specialized = hergSpec
      flat.herg_risk_level = hergSpec.risk_level || null
      flat.herg_ic50 = hergSpec.ic50_um != null ? Math.round(hergSpec.ic50_um * 10) / 10 : null
    }
  }
  // Computed: CYP inhibitions count (number of CYPs with probability > 0.5)
  const _CYP_KEYS = ['cyp1a2_inhibitor', 'cyp2c9_inhibitor', 'cyp2c19_inhibitor', 'cyp2d6_inhibitor', 'cyp3a4_inhibitor']
  const cypVals = _CYP_KEYS.map(k => flat[k]).filter(v => v != null && typeof v === 'number')
  if (cypVals.length > 0) {
    flat.cyp_inhibitions = cypVals.filter(v => v > 0.5).length
  }
  // Safety: purge stray object/array values that would crash React rendering
  const _ALLOWED_OBJECTS = new Set(['properties', 'preparation', 'breakdown', 'weights_used', 'interactions_detail', 'mapped_functional', 'scaffold_positions', 'herg_specialized', 'matched_types', 'missing_types', 'feature_counts', 'ref_feature_counts'])
  for (const k of Object.keys(flat)) {
    const v = flat[k]
    if (v != null && typeof v === 'object' && !_ALLOWED_OBJECTS.has(k)) {
      if (process.env.NODE_ENV !== 'production') {
        console.warn(`[flattenMoleculeProperties] Removing unexpected object at key "${k}":`, v)
      }
      delete flat[k]
    }
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
// Pharmacophore feature constants (shared between ProjectHome & MoleculeDetailPanel)
// ---------------------------------------------------------------------------

export const FEATURE_TYPES_ORDER = ['donor', 'acceptor', 'hydrophobic', 'aromatic', 'positive', 'negative']

export const FEATURE_LABELS = {
  donor: 'H-Bond Donor',
  acceptor: 'H-Bond Acceptor',
  hydrophobic: 'Hydrophobic',
  aromatic: 'Aromatic',
  positive: 'Positive',
  negative: 'Negative',
}

// Colors aligned with backend _FEATURE_COLORS (pharmacophore.py SVG rendering)
export const FEATURE_COLORS_HEX = {
  donor: '#3B82F5',       // (0.23, 0.51, 0.96) — blue
  acceptor: '#F04545',    // (0.94, 0.27, 0.27) — red
  hydrophobic: '#FABF26', // (0.98, 0.75, 0.15) — yellow
  aromatic: '#A854F7',    // (0.66, 0.33, 0.97) — violet
  positive: '#0FBA82',    // (0.06, 0.73, 0.51) — green
  negative: '#FA7838',    // (0.98, 0.47, 0.22) — orange
}

// ---------------------------------------------------------------------------
// Run type definitions (CDC §3.4)
// ---------------------------------------------------------------------------

export const RUN_TYPES = [
  { type: 'import', label: 'Import Molecules', icon: 'upload', description: 'Import from SDF, SMILES file, or internal selection', phases: ['hit_discovery', 'hit_to_lead', 'lead_optimization'] },
  { type: 'calculation', label: 'Run Calculation', icon: 'calculator', description: 'Select a calculation type to run on selected molecules', phases: ['hit_discovery', 'hit_to_lead', 'lead_optimization'] },
  { type: 'generation', label: 'De Novo Generation', icon: 'sparkles', description: 'Generate new molecules from selected hits', phases: ['hit_to_lead', 'lead_optimization'] },
  { type: 'afvs', label: 'AFVS — Ultra-Large Screening', icon: 'globe', description: '69B molecules (Enamine REAL Space) via AdaptiveFlow + AWS Batch', phases: ['hit_discovery'] },
]

// Labels come from GROUP_META — single source of truth for naming
export const CALCULATION_SUBTYPES = [
  { key: 'docking', icon: 'target', description: 'Dock molecules against the target protein', columns: ['docking_score', 'cnn_score', 'cnn_affinity', 'ligand_efficiency', 'cnn_vs', 'consensus_ecr', 'pocket_distance'] },
  { key: 'adme', icon: 'activity', description: 'Pharmacokinetics: how the drug circulates in the body (absorption, distribution, metabolism, excretion)', columns: ['solubility', 'BBB', 'oral_bioavailability', 'plasma_protein_binding', 'half_life', 'cyp_inhibitions', 'cns_mpo', 'pfizer_alert', 'gsk_alert', 'brenk_alert', 'pains_alert'] },
  { key: 'toxicity', icon: 'shield-alert', description: 'Toxicological risk: is the drug dangerous for the patient? (hERG, Ames, hepatotoxicity)', columns: ['hERG', 'herg_risk_level', 'herg_ic50', 'ames_mutagenicity', 'hepatotoxicity', 'skin_sensitization', 'carcinogenicity', 'safety_color_code', 'composite_score'] },
  { key: 'scoring', icon: 'star', description: 'Physicochemical descriptors (MW, LogP, TPSA), drug-likeness rules (Lipinski, QED, Ro3), SA score, pharmacophore features', columns: ['logP', 'MW', 'HBD', 'HBA', 'TPSA', 'QED', 'lipinski_pass', 'heavy_atom_count', 'sa_score', 'ro3_pass', 'inchikey', 'n_hydrophobic', 'n_aromatic', 'n_pos_ionizable', 'n_neg_ionizable'] },
  { key: 'interactions', icon: 'share-2', description: 'Protein-ligand interaction analysis (H-bonds, hydrophobic, pi-stacking)', columns: ['n_interactions', 'functional_contacts', 'total_functional', 'interaction_quality', 'key_hbonds'] },
  { key: 'scaffold', icon: 'hexagon', description: 'Murcko scaffold decomposition, BRICS bonds, R-group positions, chemical series grouping, and Butina diversity clustering', columns: ['scaffold_smiles', 'n_modifiable_positions', 'brics_bond_count', 'scaffold_n_rings', 'scaffold_group', 'scaffold_group_size', 'cluster_id', 'cluster_size', 'tanimoto_to_centroid', 'is_representative'] },
  { key: 'off_target', icon: 'crosshair', description: 'Assess selectivity against off-target proteins', columns: ['selectivity_score', 'off_target_hits', 'selectivity_ratio'] },
  // confidence retired — PAINS moved to ADME (structural alert), applicability domain to ADME
  { key: 'retrosynthesis', icon: 'git-branch', description: 'Retrosynthetic planning via AiZynthFinder (USPTO) or RDKit heuristic fallback. Steps, confidence, cost category. Reagent availability only with AiZynthFinder.', columns: ['n_synth_steps', 'synth_confidence', 'synth_cost_estimate', 'reagents_available'] },
  { key: 'pharmacophore', icon: 'hexagon', description: 'Score molecules against reference ligands using pharmacophoric fingerprint similarity (2D Gobbi)', columns: ['pharmacophore_fit_1', 'pharmacophore_fit_2', 'pharmacophore_fit_3', 'pharmacophore_fit_4', 'pharmacophore_fit_5', 'pharmacophore_fit_best', 'pharmacophore_fit_avg'] },
  { key: 'composite', icon: 'trophy', description: 'Weighted multi-criteria score combining docking, ADME, selectivity, and drug-likeness results from previous runs', columns: ['weighted_score'] },
].map(st => ({ ...st, label: GROUP_META[st.key].label }))

// Estimated run times (displayed in RunCreator confirmation)
export const ESTIMATED_TIMES = {
  import: '~5 seconds', calculation: '~1-5 minutes', generation: '~5-10 minutes',
  docking: '~3-5 min', adme: '~30 sec', toxicity: '~15 sec', scoring: '~15 sec',
  interactions: '~1-2 min', scaffold: '~30 sec', off_target: '~1 min',
  retrosynthesis: '~1 min',
  pharmacophore: '~30 sec',
  composite: '~10 sec',
}

// ---------------------------------------------------------------------------
// Column defaults for RunCreator checklist (Step 3)
// Columns NOT listed here: checked by default, no dependency.
// ---------------------------------------------------------------------------
export const COLUMN_DEFAULTS = {
  // Scoring — optional
  heavy_atom_count:    { checked: false },
  ro3_pass:            { checked: false },
  inchikey:            { checked: false },
  n_hydrophobic:       { checked: false },
  n_aromatic:          { checked: false },
  n_pos_ionizable:     { checked: false },
  n_neg_ionizable:     { checked: false },
  ligand_efficiency: { checked: true },
  // ADME — piggyback optionals
  cns_mpo:      { checked: false },
  pfizer_alert: { checked: false },
  gsk_alert:    { checked: false },
  brenk_alert:  { checked: false },
  pains_alert:  { checked: false },
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
