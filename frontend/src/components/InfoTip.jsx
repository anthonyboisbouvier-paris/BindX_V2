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

  // Auto-position: check if popup would overflow viewport
  useEffect(() => {
    if (!open || !popupRef.current || !tipRef.current) return
    const rect = tipRef.current.getBoundingClientRect()
    const popupHeight = popupRef.current.offsetHeight
    if (rect.bottom + popupHeight + 8 > window.innerHeight) {
      setPosition('top')
    } else {
      setPosition('bottom')
    }
  }, [open])

  if (!text) return null

  const iconSize = size === 'xs' ? 'w-3.5 h-3.5 text-[8px]' : 'w-4 h-4 text-[10px]'

  const isDark = variant === 'dark'
  const popupClasses = isDark
    ? 'bg-gray-800 border-gray-700 text-gray-200 shadow-xl shadow-black/20'
    : 'bg-bx-s2 border-bx-s3/60 text-slate-300 shadow-xl shadow-black/30'
  const arrowBg = isDark ? 'bg-gray-800 border-gray-700' : 'bg-bx-s2 border-bx-s3/60'

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
          className={`absolute z-[60] w-64 px-3 py-2.5 rounded-lg border
                     text-xs leading-relaxed
                     ${popupClasses}
                     ${position === 'top' ? 'bottom-full mb-2' : 'top-full mt-2'}
                     left-1/2 -translate-x-1/2`}
          style={{ pointerEvents: 'none' }}
          role="tooltip"
        >
          <div
            className={`absolute w-2 h-2 rotate-45 left-1/2 -translate-x-1/2 ${arrowBg} ${
              position === 'top'
                ? 'bottom-[-5px] border-r border-b'
                : 'top-[-5px] border-l border-t'
            }`}
          />
          {text}
        </div>
      )}
    </span>
  )
}

// ---------------------------------------------------------------------------
// Tooltip dictionary — centralized explanations for all concepts
// ---------------------------------------------------------------------------
export const TIPS = {
  // Column tooltips
  docking_score: 'Docking score from GNINA. Lower (more negative) = stronger predicted binding affinity. Typically ranges from -12 to 0 kcal/mol.',
  cnn_score: 'Convolutional Neural Network pose score (0-1). Higher = more likely to be the correct binding pose. Above 0.7 is good.',
  cnn_affinity: 'CNN-predicted binding affinity in pKd. Higher = stronger binding. Above 6 is considered good.',
  cnn_vs: 'CNN Virtual Screening score. A composite ranking metric combining pose quality and affinity.',
  composite_score: 'Weighted composite score combining docking, ADMET, and drug-likeness. Higher = better overall drug candidate.',
  smiles: 'SMILES (Simplified Molecular Input Line Entry System). A text notation for chemical structure.',
  MW: 'Molecular Weight in Daltons. Drug-like molecules typically fall between 150-500 Da (Lipinski Ro5).',
  logP: 'Octanol-water partition coefficient. Measures lipophilicity. Ideal range: -0.4 to 5.6 (Lipinski: ≤5).',
  HBD: 'Hydrogen Bond Donors. Number of OH and NH groups. Lipinski limit: ≤5.',
  HBA: 'Hydrogen Bond Acceptors. Number of O and N atoms. Lipinski limit: ≤10.',
  TPSA: 'Topological Polar Surface Area (Å²). Predicts oral absorption and BBB penetration. Below 140 Å² for oral drugs.',
  QED: 'Quantitative Estimate of Drug-likeness (0-1). Higher = more drug-like. Above 0.5 is generally favorable.',
  solubility: 'Predicted aqueous solubility. Higher is better for oral bioavailability.',
  BBB: 'Blood-Brain Barrier permeability prediction. Higher score = more likely to cross the BBB.',
  hERG: 'hERG channel inhibition risk. High values indicate potential cardiac toxicity (QT prolongation).',
  metabolic_stability: 'Predicted metabolic stability. Higher = longer half-life, less rapid clearance.',
  oral_bioavailability: 'Predicted fraction of drug reaching systemic circulation after oral administration.',
  lipinski_pass: 'Lipinski Rule of 5: MW ≤ 500, logP ≤ 5, HBD ≤ 5, HBA ≤ 10. Predicts oral absorption.',
  safety_color_code: 'Overall safety assessment. Green = low risk, Yellow = moderate caution, Red = high concern.',
  cluster_id: 'Chemical cluster assignment from structural similarity analysis (Butina clustering).',
  generation_level: 'Number of generative design iterations from the original scaffold.',
  interactions_count: 'Total number of protein-ligand interactions detected by ProLIF analysis.',

  // Run types
  run_import: 'Import molecules from external databases (ChEMBL, ZINC), file upload (SDF/CSV), or from a previous phase\'s bookmarked hits.',
  run_docking: 'Molecular docking using GNINA (GPU-accelerated). Predicts binding pose and affinity of molecules in the target\'s binding pocket.',
  run_admet: 'ADMET property prediction (Absorption, Distribution, Metabolism, Excretion, Toxicity). Evaluates drug-likeness and safety profiles.',
  run_scoring: 'Composite scoring: combines docking scores, ADMET properties, and drug-likeness into a single weighted ranking.',
  run_enrichment: 'Enrichment analysis: ProLIF interaction fingerprints, clustering, and structural analysis.',
  run_generation: 'AI-driven molecule generation. Creates novel analogs via scaffold hopping, fragment linking, or R-group exploration.',
  run_clustering: 'Chemical clustering using Butina algorithm on molecular fingerprints. Groups similar molecules together.',

  // Concepts
  pareto: 'Pareto front analysis: identifies molecules that are optimal across two competing objectives (e.g., potency vs. selectivity). Points on the Pareto front cannot be improved on one axis without worsening the other.',
  bookmark: 'Bookmarked molecules are flagged as hits of interest. They can be promoted to the next phase of the drug discovery pipeline.',
  freeze: 'Freezing a phase locks its data. Bookmarked molecules are preserved for the next phase. Unfreezing is possible but requires confirmation.',
  phase_a: 'Hit Discovery: Initial screening to find molecules that bind to the target. Uses large libraries and fast filtering.',
  phase_b: 'Hit-to-Lead: Refine initial hits with detailed analysis. Evaluate ADMET, selectivity, and SAR relationships.',
  phase_c: 'Lead Optimization: Optimize lead compounds for potency, selectivity, and drug-like properties. Final candidates for preclinical.',
  campaign: 'A campaign groups phases for one target-pocket combination. It defines the screening strategy and scoring weights.',
  pocket: 'A binding pocket is a cavity on the protein surface where a drug molecule can bind. Identified by P2Rank or manual selection.',

  // Annotation columns
  tags: 'Custom tags for organizing and categorizing molecules. Add labels like "promising", "backup", "review" to track your analysis.',
  invalidated: 'Mark a molecule as invalid (e.g., PAINS alert, synthesis issue, reviewer decision). Invalidated rows are dimmed in the table.',
  user_comment: 'Free-text notes for your analysis. Add observations, reminders, or decisions directly in the table.',
  ai_comment: 'AI-generated notes from the agent analysis. Automatically filled when the AI Agent reviews molecules.',
}
