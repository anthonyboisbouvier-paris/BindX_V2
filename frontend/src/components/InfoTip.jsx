import React, { useState, useRef, useEffect } from 'react'

/**
 * InfoTip — Pedagogical tooltip component
 *
 * Shows a (i) icon that reveals an explanation popup on hover/click.
 * Used across the dashboard to explain scientific concepts, column meanings,
 * calculation types, and features.
 *
 * Props:
 *   text     — tooltip content (string or JSX)
 *   variant  — 'dark' (default, for light backgrounds) | 'light' (for dark backgrounds)
 *   size     — 'sm' (default) | 'xs'
 */
export default function InfoTip({ text, variant = 'dark', size = 'sm', className = '' }) {
  const [open, setOpen] = useState(false)
  const [position, setPosition] = useState('bottom')
  const [hAlign, setHAlign] = useState('center') // 'center' | 'left' | 'right'
  const tipRef = useRef(null)
  const popupRef = useRef(null)

  // Close on click outside
  useEffect(() => {
    if (!open) return
    const handler = (e) => {
      if (tipRef.current && !tipRef.current.contains(e.target)) {
        setOpen(false)
      }
    }
    document.addEventListener('mousedown', handler)
    return () => document.removeEventListener('mousedown', handler)
  }, [open])

  // Compute fixed position for popup (escapes overflow-hidden containers)
  const [popupStyle, setPopupStyle] = useState({})
  useEffect(() => {
    if (!open || !popupRef.current || !tipRef.current) return
    const rect = tipRef.current.getBoundingClientRect()
    const popupRect = popupRef.current.getBoundingClientRect()

    // Vertical: prefer below, flip to above if no space
    const goUp = rect.bottom + popupRect.height + 8 > window.innerHeight
    setPosition(goUp ? 'top' : 'bottom')

    // Horizontal: center on icon, push left/right if overflowing
    let left = rect.left + rect.width / 2 - popupRect.width / 2
    let hA = 'center'
    if (left < 8) { left = rect.left; hA = 'left' }
    else if (left + popupRect.width > window.innerWidth - 8) { left = rect.right - popupRect.width; hA = 'right' }
    setHAlign(hA)

    const top = goUp ? rect.top - popupRect.height - 8 : rect.bottom + 8
    setPopupStyle({ left, top })
  }, [open])

  if (!text) return null

  const iconSize = size === 'xs' ? 'w-3.5 h-3.5 text-[8px]' : 'w-4 h-4 text-[10px]'

  const isDark = variant === 'dark'
  const popupClasses = isDark
    ? 'bg-gray-800 border-gray-700 text-gray-200 shadow-xl shadow-black/20'
    : 'bg-bx-s2 border-bx-s3/60 text-slate-300 shadow-xl shadow-black/30'

  return (
    <span
      ref={tipRef}
      className={`inline-flex items-center relative ${className}`}
      onMouseEnter={() => setOpen(true)}
      onMouseLeave={() => setOpen(false)}
    >
      <button
        type="button"
        onClick={(e) => { e.stopPropagation(); setOpen(v => !v) }}
        className={`ml-1 ${iconSize} rounded-full
                   bg-gray-200 hover:bg-blue-100
                   text-gray-400 hover:text-blue-500
                   font-bold flex items-center justify-center transition-colors cursor-help
                   focus:outline-none focus:ring-1 focus:ring-blue-300`}
        aria-label="More information"
      >
        i
      </button>

      {open && (
        <div
          ref={popupRef}
          className={`fixed z-[9999] max-w-[220px] px-3 py-2.5 rounded-lg border
                     text-xs leading-relaxed break-words whitespace-normal
                     ${popupClasses}`}
          style={{ pointerEvents: 'none', ...popupStyle }}
          role="tooltip"
        >
          {text}
        </div>
      )}
    </span>
  )
}

// ---------------------------------------------------------------------------
// Tooltip dictionary — centralized explanations for all concepts
// English, accessible, scientifically accurate. 2-4 sentences per entry.
// Explains HOW each value is computed (which columns/properties feed into it).
// ---------------------------------------------------------------------------
export const TIPS = {
  // --- Columns: Docking ---
  docking_score: 'Predicted binding energy between the molecule and the protein pocket (kcal/mol). Computed by GNINA using CNN scoring on the docked pose. More negative = stronger binding. Good candidates typically score below -7 kcal/mol.',
  cnn_score: 'Neural network pose confidence score (0-1). GNINA\'s CNN evaluates whether the molecule is correctly oriented in the binding pocket. Above 0.7, the pose is considered reliable.',
  cnn_affinity: 'CNN-predicted binding affinity expressed as pKd. Computed by GNINA\'s deep learning model from the 3D docked pose. Higher values mean stronger binding; above 6 is promising.',
  cnn_vs: 'Virtual screening score combining pose quality (cnn_score) and predicted affinity (cnn_affinity). Used to rank molecules in large-scale screening campaigns.',
  consensus_ecr: 'Exponential Consensus Ranking: aggregates multiple docking scores (docking_score, cnn_score, cnn_affinity) into a single ranking. Higher values indicate consistently good ranking across all scoring methods.',

  // --- Columns: ADMET ---
  logP: 'Lipophilicity measure — the molecule\'s ability to cross lipid membranes. Computed from the molecular structure (atom contributions). Ideal range: -0.4 to 5. Above 5, poor absorption and increased toxicity risk.',
  MW: 'Molecular weight in Daltons. Computed directly from the molecular formula. Oral drugs are typically 150-500 Da. Above 500, intestinal absorption drops significantly (Lipinski Rule of 5).',
  HBD: 'Number of hydrogen bond donors (OH and NH groups). Counted from the molecular structure. Lipinski\'s Rule of 5 recommends a maximum of 5 for good oral absorption.',
  HBA: 'Number of hydrogen bond acceptors (O and N atoms). Counted from the molecular structure. Lipinski recommends a maximum of 10. Excess acceptors reduce membrane permeability.',
  TPSA: 'Topological Polar Surface Area in Angstroms squared. Calculated from the sum of polar atom surfaces (O, N, attached H). Below 140 for oral drugs; below 90 for blood-brain barrier penetration.',
  QED: 'Quantitative Estimate of Drug-likeness (0-1). Combines MW, logP, HBD, HBA, TPSA, rotatable bonds, aromatic rings, and alerts into a single score. Above 0.5 indicates a favorable drug-like profile.',
  solubility: 'Predicted aqueous solubility. Estimated from molecular descriptors (logP, MW, aromatic ring count). Essential for drug dissolution and absorption. Higher values are better.',
  BBB: 'Blood-Brain Barrier penetration prediction. Computed from logP, TPSA, MW, and HBD. High score = molecule can reach the brain. Desirable for neurological targets, undesirable otherwise.',
  hERG: 'Predicted hERG potassium channel inhibition risk. Computed from molecular descriptors. High values signal cardiac toxicity risk (QT prolongation). One of the most critical safety endpoints in drug development.',
  half_life: 'Predicted elimination half-life in hours. Estimated from molecular descriptors and ADMET-AI models. Indicates how long the drug stays active in the body. Longer half-life allows for less frequent dosing.',
  cyp_inhibitions: 'Number of CYP450 enzymes inhibited (probability > 50%). Counts CYP1A2, CYP2C9, CYP2C19, CYP2D6, CYP3A4. 0 = clean metabolic profile, ≥2 = high drug-drug interaction risk. Detail per CYP available in the ADME panel.',
  oral_bioavailability: 'Fraction of the drug reaching systemic circulation after oral administration (0-100%). Predicted from absorption and first-pass metabolism models. Above 30% is generally acceptable.',
  lipinski_pass: 'Lipinski\'s Rule of 5: checks if MW <= 500, logP <= 5, HBD <= 5, HBA <= 10. All four properties are read from the corresponding dashboard columns. Pass = the molecule is likely orally bioavailable.',
  plasma_protein_binding: 'Fraction bound to blood plasma proteins. Predicted from molecular descriptors. Highly bound molecules (>95%) have less free drug available to reach the target, potentially reducing efficacy.',
  cns_mpo: 'CNS Multi-Parameter Optimization score (0-6). Computed from logP, TPSA, MW, HBD, pKa, and CLint. Each property contributes 0-1 points. Above 4 indicates good CNS drug potential.',
  pfizer_alert: 'Pfizer 3/75 alert: flags molecules with logP > 3 AND TPSA < 75 (both from dashboard columns). This lipophilic + low-polarity profile is associated with increased toxicity risk.',
  gsk_alert: 'GSK 4/400 alert: flags molecules with logP > 4 AND MW > 400 (both from dashboard columns). This profile is associated with higher clinical failure rates.',
  heavy_atom_count: 'Number of non-hydrogen atoms. Counted directly from the molecular structure. Typical drugs contain 15-35 heavy atoms. Used to compute ligand efficiency (docking_score / heavy_atom_count).',
  sa_score: 'Synthetic Accessibility score (1-10). Computed from molecular complexity, ring systems, and stereocenters. Lower = easier to synthesize. Below 4 is considered readily accessible.',
  ro3_pass: 'Rule of 3 for fragments: checks MW <= 300, logP <= 3, HBD <= 3, HBA <= 3. Used to evaluate whether a molecule is a good starting fragment for optimization campaigns.',

  // --- Columns: Scoring ---
  composite_score: 'ADMET safety score combining toxicity endpoints into a single metric (0-100). Aggregates hERG, Ames, hepatotoxicity, skin sensitization, and carcinogenicity. Higher = safer compound. Not to be confused with the Weighted Score.',
  weighted_score: 'Weighted multi-criteria score combining your selected metrics from previous runs (docking, ADME, selectivity, drug-likeness). Configure weights when creating the Composite Score run. Higher = better overall candidate.',
  ligand_efficiency: 'Ligand Efficiency = -docking_score / heavy_atom_count. Normalizes binding affinity by molecular size. Allows fair comparison between molecules of different sizes. Higher is better.',

  // --- Columns: Interactions ---
  n_interactions: 'Number of specific protein-ligand contacts (H-bonds, hydrophobic, pi-stacking, salt bridges). Computed by ProLIF or RDKit distance-based analysis from the docked pose. More contacts generally indicate a more stable binding mode.',
  functional_contacts: 'Number of contacts with functionally important residues in the binding site. Identified by cross-referencing detected interactions with known catalytic/binding residues from UniProt.',
  interaction_quality: 'Quality score (0-1) based on the proportion of functionally important residues contacted. Higher quality means the molecule interacts with key residues for target activity.',
  key_hbonds: 'Number of hydrogen bonds with key binding site residues. H-bonds with catalytic or conserved residues are the strongest indicator of specific, high-affinity binding.',

  // --- Columns: Scaffold ---
  scaffold_smiles: 'Murcko scaffold — the core ring system after removing side chains. Computed using RDKit\'s MurckoScaffold decomposition. Molecules sharing the same scaffold belong to the same chemical series.',
  n_modifiable_positions: 'Number of R-group positions where chemical modifications are possible. Identified by BRICS decomposition and substituent detection. More positions = more room for structure-activity optimization.',
  brics_bond_count: 'Number of BRICS retrosynthetically accessible bonds. These are bonds that can be cleaved to generate building blocks for combinatorial chemistry. Higher count indicates a more modular molecule.',
  scaffold_n_rings: 'Number of rings in the Murcko scaffold. Ring count affects molecular complexity, rigidity, and drug-likeness. Most oral drugs have 2-4 rings.',

  // --- Columns: Clustering ---
  cluster_id: 'Chemical cluster number assigned by Butina clustering algorithm. Molecules are grouped by Tanimoto similarity of Morgan fingerprints (radius 2). Same cluster = similar chemical scaffold.',
  scaffold_smiles: 'SMILES notation of the scaffold shared by molecules in this cluster. Represents the common chemical motif of the group.',
  tanimoto_to_centroid: 'Tanimoto similarity (0-1) between this molecule\'s Morgan fingerprint and the cluster centroid. Closer to 1 = more representative of the chemical group.',
  is_representative: 'Whether this molecule is the centroid (most representative) of its cluster. Computed as the molecule with highest average similarity to all other cluster members. Useful for selecting one exemplar per group.',

  // --- Columns: Off-target ---
  selectivity_score: 'Selectivity score (0-1). Computed by comparing the docking score on the intended target vs. scores on anti-targets (hERG, CYP450, etc.). Higher = more selective for the intended target.',
  off_target_hits: 'Number of off-target proteins with significant predicted binding. Computed from docking against a panel of anti-targets. Fewer hits = more specific and safer molecule.',
  selectivity_ratio: 'Ratio of target affinity to the best off-target affinity (docking_score_target / docking_score_best_off_target). A ratio above 10 indicates excellent selectivity.',

  // --- Columns: Confidence ---
  confidence_score: 'Overall confidence score (0-1). Aggregates reliability indicators: PAINS alerts, Brenk alerts, applicability domain check, and prediction convergence across pipeline steps. Low score = interpret results with caution.',
  pains_alert: 'PAINS (Pan-Assay Interference Compounds) alert. Detected by matching the molecule\'s structure against known problematic substructure patterns. These molecules are prone to false positives in biological assays.',
  applicability_domain: 'Whether the molecule falls within the training domain of the prediction models. Computed by comparing molecular descriptors to the models\' training set. Outside the domain = less reliable predictions.',
  brenk_alert: 'Brenk structural alert: detects reactive or unstable chemical groups known to cause problems in drug development. Computed by substructure matching against Brenk\'s published filter list. An alert does not disqualify but warrants attention.',

  // --- Columns: Retrosynthesis ---
  n_synth_steps: 'Number of synthetic steps needed to make this molecule from commercially available reagents. Computed by AI retrosynthesis (AiZynthFinder). Fewer steps = faster and cheaper synthesis.',
  synth_confidence: 'Confidence in the proposed synthesis route feasibility (0-1). Computed by the retrosynthesis AI based on route plausibility and reaction precedent. Above 0.7 = realistic route; below 0.3 = uncertain.',
  synth_cost_estimate: 'Estimated synthesis cost: low, medium, or high. Derived from number of steps (n_synth_steps), reagent availability (reagents_available), and reaction complexity.',
  reagents_available: 'Whether all required starting materials are commercially available. Checked against supplier catalogs by the retrosynthesis pipeline. If not, additional preparation steps are needed.',

  // --- Columns: Safety ---
  herg_risk: 'hERG channel inhibition risk (0-1). Predicted from molecular descriptors related to cation-pi interactions and lipophilicity. This potassium channel regulates heartbeat; risk >0.5 can cause fatal arrhythmias. Elimination criterion in pharma.',
  ames_mutagenicity: 'Ames mutagenicity test prediction. Computed by matching structural alerts and ML models trained on Ames test data. A positive result is a serious safety red flag indicating potential DNA damage.',
  hepatotoxicity: 'Liver toxicity risk (0-1). Predicted from molecular descriptors using models trained on DILI (Drug-Induced Liver Injury) datasets. High score signals danger of hepatic damage.',
  skin_sensitization: 'Skin sensitization risk. Predicted from reactive group detection and molecular descriptors. A positive result means the molecule could trigger allergic reactions on skin contact.',
  carcinogenicity: 'Carcinogenicity risk (0-1). Predicted from structural alerts and ML models trained on rodent carcinogenicity data. Any significant score should be investigated as a priority.',
  safety_color_code: 'Overall safety summary. Computed from herg_risk, ames_mutagenicity, hepatotoxicity, skin_sensitization, and carcinogenicity columns. Green = safe profile, Yellow = moderate signals, Red = at least one critical risk.',

  // --- Columns: Activity Cliffs ---
  is_cliff: 'Activity cliff flag: two molecules with high structural similarity (Tanimoto > 0.85) but very different activity. Computed by comparing docking_score differences vs. fingerprint similarity across all pairs.',
  sali_max: 'Maximum SALI (Structure-Activity Landscape Index). Computed as |activity_diff| / (1 - similarity) for the most contrasting neighbor pair. Higher values indicate a sharper activity cliff.',
  n_cliffs: 'Number of cliff pairs involving this molecule. Computed from pairwise SALI analysis. High count means the molecule sits in a sensitive region of the structure-activity landscape.',

  // --- Columns: Pharmacophore ---
  pharmacophore_features: 'Number of pharmacophoric features detected (H-bond donors/acceptors, hydrophobic zones, charges, aromatic rings). Computed from the 3D molecular structure. More features = richer target interaction potential.',
  pharmacophore_similarity: 'Pharmacophore similarity to the reference model (0-1). Computed by aligning 3D pharmacophore features with the best-known active compounds. High score = same 3D interaction pattern as proven hits.',

  // --- Columns: Generation ---
  generation_level: 'Generation number: how many rounds of generative AI design since the original molecule. Generation 0 = imported molecule, 1 = first optimization cycle, 2 = second cycle, etc.',
  parent_molecule_id: 'ID of the parent molecule from which this one was derived by generative design. Allows tracing the chemical genealogy across optimization cycles.',

  // --- Columns: Identity ---
  smiles: 'SMILES: universal text notation for chemical structures. Each molecule has a unique SMILES encoding its atoms, bonds, and stereochemistry. Used as the primary molecular identifier.',
  inchikey: 'InChIKey: a unique 27-character hash identifier for each molecule. Enables lookup across global databases (PubChem, ChEMBL, ZINC). Computed deterministically from the molecular structure.',

  // --- Run types ---
  run_import: 'Imports molecules from databases (ChEMBL, ZINC), files (SDF/CSV), or bookmarked hits from a previous phase. This is the first step of every phase to populate the molecule table.',
  run_docking: 'Molecular docking with GNINA (GPU-accelerated). Simulates how each molecule fits into the protein binding pocket and predicts binding energy. Produces docking_score, cnn_score, cnn_affinity columns.',
  run_admet: 'Predicts ADMET properties: Absorption, Distribution, Metabolism, Excretion, Toxicity. Computes logP, MW, HBD, HBA, TPSA, QED, solubility, BBB, hERG, and Lipinski columns from the molecular structure.',
  run_scoring: 'Computes a weighted composite score combining docking_score, QED, and ADMET properties into a single 0-100 ranking. Also computes ligand_efficiency from docking_score and heavy_atom_count.',
  run_interactions: 'Protein-ligand interaction analysis: computes H-bonds, hydrophobic contacts, pi-stacking, and salt bridges from docked poses using ProLIF or RDKit. Requires a prior docking run.',
  run_scaffold: 'Scaffold decomposition: extracts the Murcko core ring system, identifies BRICS bonds and R-group modification positions. Only requires SMILES — no docking needed.',
  run_generation: 'AI de novo generation. Creates new molecules by modifying top hits: fragment replacement, R-group exploration, scaffold hopping. Produces new molecules with generation_level and parent_molecule_id.',
  run_clustering: 'Chemical clustering using Butina algorithm on Morgan fingerprints (Tanimoto similarity). Groups similar molecules to identify chemical series and ensure structural diversity in the selection.',
  run_off_target: 'Off-target selectivity: docks molecules against anti-target proteins (hERG, CYP450, etc.) to verify they don\'t bind unintended targets. Produces selectivity_score, off_target_hits, selectivity_ratio.',
  run_confidence: 'Confidence analysis: evaluates prediction reliability via PAINS alerts, Brenk filters, applicability domain checks, and cross-method convergence. Produces confidence_score and alert columns.',
  run_retrosynthesis: 'AI retrosynthesis: plans the chemical synthesis route, estimates step count (n_synth_steps), cost (synth_cost_estimate), confidence, and checks reagent availability.',
  run_composite: 'Weighted composite score: aggregates results from your previous runs (docking, ADME, scoring, off-target) into a single multi-criteria ranking. You configure the weight of each metric. Produces weighted_score column.',
  run_safety: 'Full safety profile: hERG risk (cardiac), Ames mutagenicity, hepatotoxicity, skin sensitization, and carcinogenicity. Aggregated into a color-coded safety summary (safety_color_code).',

  // --- Concepts ---
  pareto: 'Pareto front analysis: identifies optimal molecules when two objectives compete (e.g., potency vs. selectivity). Points on the front cannot be improved on one axis without degrading the other. Select any two numeric columns as axes.',
  bookmark: 'Bookmarked molecules are flagged as hits of interest. They can be promoted to the next pipeline phase via the "Next Phase" action. Bookmarks are preserved when the phase is frozen.',
  freeze: 'Freezing a phase locks its data. Bookmarked molecules are preserved for the next phase. Unfreezing is possible but requires confirmation, especially if downstream phases have runs.',
  phase_a: 'Hit Discovery: initial screening to find molecules that bind the target. Uses large compound libraries and rapid filtering (docking + ADMET) to identify first hits from thousands of candidates.',
  phase_b: 'Hit-to-Lead: refines initial hits with detailed analysis. Evaluates ADMET, off-target selectivity, and structure-activity relationships to select the most promising lead compounds.',
  phase_c: 'Lead Optimization: optimizes leads for potency, selectivity, and drug-like properties. Uses generative AI, retrosynthesis, and safety profiling. Best candidates advance to preclinical studies.',
  campaign: 'A campaign groups all phases for one target + pocket combination. It defines the screening strategy and the composite score weights used to rank molecules.',
  pocket: 'Binding pocket: a cavity on the protein surface where a drug molecule can bind. Automatically detected by P2Rank or manually selected from the 3D structure during target setup.',

  // --- Annotation columns ---
  tags: 'Custom tags to organize your molecules. Add labels like "promising", "backup", "review needed" to structure your analysis workflow.',
  invalidated: 'Marks a molecule as invalid (e.g., PAINS alert, synthesis issue, reviewer decision). Invalidated rows are visually dimmed in the table but not deleted.',
  user_comment: 'Free-text notes for your analysis. Add observations, reminders, or decisions directly in the table. Visible to all project members.',
  ai_comment: 'Notes generated by the AI agent during automated analysis. Populated when the Campaign Agent reviews molecules and identifies patterns or concerns.',
}
