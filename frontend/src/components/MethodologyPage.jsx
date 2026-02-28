import React, { useState } from 'react'
import { useNavigate } from 'react-router-dom'

// --------------------------------------------------
// Section icons
// --------------------------------------------------
const ICONS = {
  structure: (
    <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
        d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
    </svg>
  ),
  pocket: (
    <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
        d="M9 20l-5.447-2.724A1 1 0 013 16.382V5.618a1 1 0 011.447-.894L9 7m0 13l6-3m-6 3V7m6 10l4.553 2.276A1 1 0 0021 18.382V7.618a1 1 0 00-.553-.894L15 4m0 13V4m0 0L9 7" />
    </svg>
  ),
  docking: (
    <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
        d="M13.828 10.172a4 4 0 00-5.656 0l-4 4a4 4 0 105.656 5.656l1.102-1.101m-.758-4.899a4 4 0 005.656 0l4-4a4 4 0 00-5.656-5.656l-1.1 1.1" />
    </svg>
  ),
  generation: (
    <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
        d="M9.663 17h4.673M12 3v1m6.364 1.636l-.707.707M21 12h-1M4 12H3m3.343-5.657l-.707-.707m2.828 9.9a5 5 0 117.072 0l-.548.547A3.374 3.374 0 0014 18.469V19a2 2 0 11-4 0v-.531c0-.895-.356-1.754-.988-2.386l-.548-.547z" />
    </svg>
  ),
  admet: (
    <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
        d="M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 01-2 2h-2a2 2 0 01-2-2z" />
    </svg>
  ),
  retrosynthesis: (
    <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
        d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
    </svg>
  ),
  safety: (
    <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
        d="M9 12l2 2 4-4m5.618-4.016A11.955 11.955 0 0112 2.944a11.955 11.955 0 01-8.618 3.04A12.02 12.02 0 003 9c0 5.591 3.824 10.29 9 11.622 5.176-1.332 9-6.03 9-11.622 0-1.042-.133-2.052-.382-3.016z" />
    </svg>
  ),
  optimization: (
    <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
        d="M13 7h8m0 0v8m0-8l-8 8-4-4-6 6" />
    </svg>
  ),
}

// --------------------------------------------------
// Section data
// --------------------------------------------------
const SECTIONS = [
  {
    id: 'structure',
    icon: 'structure',
    title: 'Structure Prediction',
    tool: 'AlphaFold2 (DeepMind) or ESMFold (Meta)',
    description:
      'Predicts the 3D atomic structure of a protein from its amino acid sequence using deep learning. AlphaFold2 is used first (via the EBI AlphaFold DB); ESMFold provides a fallback when no precomputed structure is available. The predicted structure is provided as a PDB file, which is used in all downstream steps.',
    accuracy: 'Median GDT-TS 92.4 on CASP14 for AlphaFold2. ESMFold achieves similar performance in single-chain prediction.',
    limitation:
      'Less reliable for intrinsically disordered regions, large multi-chain complexes, and proteins with few homologs in the training set (orphan proteins).',
    reference: 'Jumper et al., Nature 2021 (AlphaFold2); Lin et al., Science 2023 (ESMFold)',
    badge: 'Step 1',
    badgeColor: 'bg-blue-100 text-blue-700',
  },
  {
    id: 'pocket',
    icon: 'pocket',
    title: 'Pocket Detection',
    tool: 'fpocket 4.0',
    description:
      'Identifies druggable binding pockets on the protein surface using Voronoi tessellation and alpha sphere methods. The largest or highest-scoring pocket is selected automatically as the docking box. Pocket druggability scores are reported alongside volume and hydrophobicity.',
    accuracy: 'Detects known binding sites in >80% of benchmark structures when a co-crystal ligand is present.',
    limitation:
      'Performance decreases for allosteric sites, cryptic pockets, and protein-protein interaction interfaces that require conformational change to become accessible.',
    reference: 'Le Guilloux et al., BMC Bioinformatics 2009; Schmidtke et al., J Chem Inf Model 2010',
    badge: 'Step 2',
    badgeColor: 'bg-indigo-100 text-indigo-700',
  },
  {
    id: 'docking',
    icon: 'docking',
    title: 'Molecular Docking',
    tool: 'GNINA 1.1 (primary) / AutoDock Vina 1.2 (fallback)',
    description:
      'Performs structure-based molecular docking using GNINA with CNN rescoring (primary) or AutoDock Vina (fallback). ' +
      'Protein structures are prepared by stripping all HETATM records (waters, co-crystallized ligands, crystallization artifacts) to ensure a clean binding pocket. ' +
      'Binding pockets are defined prior to docking via P2Rank/fpocket, and ligands are docked within a fixed search box centered on the pocket. ' +
      'GNINA extends Vina with convolutional neural network scoring, producing three complementary scores: Vina affinity (kcal/mol), CNN pose score (0-1), and CNN predicted affinity (pK units). ' +
      'Docking results provide atom-level 3D poses directly from the docking engine, without any post-processing or coordinate manipulation. ' +
      'Ligand positions are expected to lie within or near the binding cavity rather than at the geometric center of the pocket. ' +
      'When no valid docking pose is available, only 2D representations are shown to avoid misleading interpretations. ' +
      'Docking scores are used for relative ranking among compounds, not as absolute affinity predictions.',
    accuracy:
      'GNINA reproduces co-crystal poses within 2A RMSD in ~70-75% of benchmark cases (vs ~60-70% for Vina alone). ' +
      'CNN rescoring significantly improves pose selection accuracy. ' +
      'Relative ranking of known actives vs. decoys shows enrichment factors (EF1%) of 5-15x across standard benchmarks (DUD-E, CASF-2016). ' +
      'A Spearman rank correlation of 0.3-0.5 between docking scores and experimental IC50 values is typical and expected for structure-based methods.',
    limitation:
      'Rigid receptor approximation does not capture induced-fit binding. ' +
      'Scoring accuracy decreases for metalloproteins, covalent binders, and highly flexible ligands (>15 rotatable bonds). ' +
      'Absolute docking scores should NOT be interpreted as binding free energies — use only for relative ranking within the same target.',
    reference: 'McNutt et al., J Cheminform 2021 (GNINA); Eberhardt et al., J Chem Inf Model 2021 (Vina 1.2); Trott & Olson, J Comput Chem 2010',
    badge: 'Step 3',
    badgeColor: 'bg-purple-100 text-purple-700',
  },
  {
    id: 'generation',
    icon: 'generation',
    title: 'AI Molecule Generation',
    tool: 'REINVENT4 (generative RL-based SMILES model)',
    description:
      'Generates novel drug-like molecules optimized toward the target pocket using reinforcement learning. A transformer-based SMILES generator is fine-tuned with a reward signal derived from the docking score and ADMET properties. Generated molecules are then re-docked to confirm their predicted affinity.',
    accuracy: 'Generated molecules typically show 10-30% improvement in docking score versus ChEMBL screening hits. Novelty (Tanimoto < 0.4 to training set) is typically >80%.',
    limitation:
      'Generated molecules may be chemically unusual or difficult to synthesize. SMILES validity can drop for very long sequences. Experimental confirmation is required for all AI-generated candidates.',
    reference: 'Loeffler et al., J Chem Inf Model 2024 (REINVENT4)',
    badge: 'Step 4',
    badgeColor: 'bg-pink-100 text-pink-700',
  },
  {
    id: 'admet',
    icon: 'admet',
    title: 'ADMET Prediction',
    tool: 'ADMET-AI (graph neural network ensemble)',
    description:
      'Predicts absorption, distribution, metabolism, excretion, and toxicity (ADMET) properties from molecular structure. The model ensemble covers >40 endpoints including Caco-2 permeability, hERG inhibition, CYP450 metabolism, aqueous solubility, plasma protein binding, and LD50 toxicity. A composite ADMET score (0-1) is computed and color-coded green/yellow/red.',
    accuracy: 'Typical AUROC 0.80-0.90 on internal test sets. Performance varies by endpoint; simpler endpoints (e.g. LogD) are more accurate than complex ones (e.g. CNS penetration).',
    limitation:
      'Predictions are based on molecular fingerprints/graphs and may not capture metabolic soft spots, reactive intermediates, or species-specific differences. In vitro confirmation is required before progressing candidates.',
    reference: 'Swanson et al., Bioinformatics 2023 (ADMET-AI)',
    badge: 'Step 5',
    badgeColor: 'bg-green-100 text-green-700',
  },
  {
    id: 'retrosynthesis',
    icon: 'retrosynthesis',
    title: 'Retrosynthesis Planning',
    tool: 'AiZynthFinder (Monte Carlo Tree Search + template library)',
    description:
      'Plans a synthetic route for each candidate molecule by recursively applying retrosynthetic disconnection rules until commercial starting materials are reached. Each step uses a trained reaction template selector. Routes are ranked by number of steps, availability of reagents, and estimated yield.',
    accuracy: 'Finds a valid route in 80-90% of drug-like molecules in 60 seconds (Monte Carlo Tree Search, up to 5 steps). Route quality may be lower for novel scaffolds not in the template library.',
    limitation:
      'Does not consider stereoselectivity, regioselectivity, or actual reaction conditions (solvent, temperature, catalyst). Route feasibility should be validated by a medicinal chemist.',
    reference: 'Genheden et al., J Cheminform 2020 (AiZynthFinder)',
    badge: 'Step 6',
    badgeColor: 'bg-yellow-100 text-yellow-700',
  },
  {
    id: 'safety',
    icon: 'safety',
    title: 'Off-Target Safety Screening',
    tool: 'DockIt V5 — multi-target docking panel',
    description:
      'Screens each candidate against a panel of 10 key anti-targets: hERG (cardiac arrhythmia), CYP3A4, CYP2D6, CYP2C9 (drug metabolism), hNAV1.5 (cardiac), PXR (drug induction), COX-1/2 (GI toxicity), carbonic anhydrase II, and acetylcholinesterase. Uses AutoDock Vina with pre-prepared anti-target structures from the Protein Data Bank. A selectivity score (0-1) summarizes overall off-target safety.',
    accuracy: 'Off-target binding predictions have lower accuracy than primary target docking due to less experimental data for validation. Use as a triage tool only.',
    limitation:
      'Computational selectivity predictions must be confirmed with experimental assays (e.g., Cerep selectivity panel, hERG patch-clamp, CYP inhibition IC50). False negatives and positives are common.',
    reference: 'DockIt V5 internal pipeline. Anti-target structures from PDB (updated quarterly).',
    badge: 'V5 New',
    badgeColor: 'bg-orange-100 text-orange-700',
  },
  {
    id: 'optimization',
    icon: 'optimization',
    title: 'Lead Optimization',
    tool: 'DockIt V5 — multi-objective RL optimization',
    description:
      'Iteratively improves a selected lead molecule using multi-objective reinforcement learning. A generative model proposes structural variants each iteration; each variant is scored against a weighted combination of binding affinity, predicted toxicity, oral bioavailability, and synthesis complexity. The Pareto front is tracked across iterations. Users can adjust objective weights to reflect their project priorities.',
    accuracy: 'Improvement of 10-30% in composite score is typical over 10-20 iterations. Quality of results depends strongly on the quality of scoring functions used.',
    limitation:
      'Optimization is constrained to the scoring functions\' ability to capture desired properties. Unexpected failure modes (e.g., novel toxicophores introduced during optimization) may not be detected. Wet lab follow-up is mandatory.',
    reference: 'DockIt V5 internal pipeline. Based on REINVENT4 and multi-objective optimization literature.',
    badge: 'V5 New',
    badgeColor: 'bg-orange-100 text-orange-700',
  },
]

// --------------------------------------------------
// Docking integrity rules
// --------------------------------------------------
const INTEGRITY_RULES = [
  {
    icon: (
      <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12l2 2 4-4m6 2a9 9 0 11-18 0 9 9 0 0118 0z" />
      </svg>
    ),
    rule: 'Atom-level 3D poses come directly from GNINA/Vina docking output',
    detail: 'No post-processing, no artificial translation, no coordinate manipulation. What the docking engine produces is what you see.',
  },
  {
    icon: (
      <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M18.364 18.364A9 9 0 005.636 5.636m12.728 12.728A9 9 0 015.636 5.636m12.728 12.728L5.636 5.636" />
      </svg>
    ),
    rule: 'No 3D pose displayed without real docking',
    detail: 'If docking fails or is not performed, the interface shows a 2D structural diagram only. No misleading 3D visualization is ever generated.',
  },
  {
    icon: (
      <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 7h8m0 0v8m0-8l-8 8-4-4-6 6" />
      </svg>
    ),
    rule: 'Scores are for relative ranking, not absolute prediction',
    detail: 'Docking scores (Vina affinity, CNN score, CNN affinity) rank compounds within the same target. They are NOT predictions of experimental binding affinity (IC50, Kd).',
  },
  {
    icon: (
      <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M17.657 16.657L13.414 20.9a1.998 1.998 0 01-2.827 0l-4.244-4.243a8 8 0 1111.314 0z" />
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 11a3 3 0 11-6 0 3 3 0 016 0z" />
      </svg>
    ),
    rule: 'Pocket center is a reference point, not the expected ligand position',
    detail: 'Docked ligands are expected to lie 5-10 Angstroms from the geometric pocket center, within or at the surface of the binding cavity. This is physically normal.',
  },
  {
    icon: (
      <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2" />
      </svg>
    ),
    rule: 'Receptor preparation strips all non-protein atoms',
    detail: 'Waters (HOH), co-crystallized ligands, and crystallization artifacts are automatically removed before docking. Only protein ATOM records are kept to ensure a clean binding pocket.',
  },
  {
    icon: (
      <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
      </svg>
    ),
    rule: 'Pipeline is stable and reproducible',
    detail: 'GNINA uses deterministic random seeds. Repeated runs on the same input produce consistent rankings (top-5 overlap > 80%). Score variation between runs is typically < 0.5 kcal/mol.',
  },
]

// --------------------------------------------------
// Section card
// --------------------------------------------------
function SectionCard({ section, isOpen, onToggle }) {
  return (
    <div className="card overflow-hidden">
      {/* Header — always visible */}
      <button
        type="button"
        onClick={onToggle}
        className="w-full flex items-center gap-4 px-5 py-4 text-left hover:bg-gray-50 transition-colors"
      >
        <div className="w-10 h-10 rounded-lg bg-bx-surface flex items-center justify-center text-bx-mint flex-shrink-0">
          {ICONS[section.icon]}
        </div>
        <div className="flex-1 min-w-0">
          <div className="flex items-center gap-2 flex-wrap">
            <h3 className="font-bold text-bx-light-text text-base">{section.title}</h3>
            <span className={`inline-flex items-center px-2 py-0.5 rounded-full text-xs font-semibold ${section.badgeColor}`}>
              {section.badge}
            </span>
          </div>
          <p className="text-xs text-gray-500 mt-0.5 truncate">{section.tool}</p>
        </div>
        <svg
          className={`w-5 h-5 text-gray-400 flex-shrink-0 transition-transform ${isOpen ? 'rotate-180' : ''}`}
          fill="none" stroke="currentColor" viewBox="0 0 24 24"
        >
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
        </svg>
      </button>

      {/* Expanded content */}
      {isOpen && (
        <div className="px-5 pb-5 border-t border-gray-50">
          {/* Tool name */}
          <div className="mt-4 mb-3">
            <span className="text-xs font-semibold text-gray-400 uppercase tracking-wide">Tool</span>
            <p className="text-sm font-medium text-bx-light-text mt-0.5">{section.tool}</p>
          </div>

          {/* Description */}
          <div className="mb-4">
            <span className="text-xs font-semibold text-gray-400 uppercase tracking-wide">What it does</span>
            <p className="text-sm text-gray-600 mt-1 leading-relaxed">{section.description}</p>
          </div>

          {/* Accuracy / Limitation / Reference */}
          <div className="grid grid-cols-1 sm:grid-cols-3 gap-4">
            <div className="bg-green-50 rounded-lg p-3">
              <div className="flex items-center gap-1.5 mb-1">
                <svg className="w-3.5 h-3.5 text-green-600 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12l2 2 4-4m6 2a9 9 0 11-18 0 9 9 0 0118 0z" />
                </svg>
                <span className="text-xs font-semibold text-green-700 uppercase tracking-wide">Accuracy</span>
              </div>
              <p className="text-xs text-green-800 leading-relaxed">{section.accuracy}</p>
            </div>

            <div className="bg-amber-50 rounded-lg p-3">
              <div className="flex items-center gap-1.5 mb-1">
                <svg className="w-3.5 h-3.5 text-amber-600 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 9v2m0 4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
                </svg>
                <span className="text-xs font-semibold text-amber-700 uppercase tracking-wide">Limitation</span>
              </div>
              <p className="text-xs text-amber-800 leading-relaxed">{section.limitation}</p>
            </div>

            <div className="bg-gray-50 rounded-lg p-3">
              <div className="flex items-center gap-1.5 mb-1">
                <svg className="w-3.5 h-3.5 text-gray-500 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 6.253v13m0-13C10.832 5.477 9.246 5 7.5 5S4.168 5.477 3 6.253v13C4.168 18.477 5.754 18 7.5 18s3.332.477 4.5 1.253m0-13C13.168 5.477 14.754 5 16.5 5c1.746 0 3.332.477 4.5 1.253v13C19.832 18.477 18.246 18 16.5 18c-1.746 0-3.332.477-4.5 1.253" />
                </svg>
                <span className="text-xs font-semibold text-gray-500 uppercase tracking-wide">Reference</span>
              </div>
              <p className="text-xs text-gray-600 leading-relaxed italic">{section.reference}</p>
            </div>
          </div>
        </div>
      )}
    </div>
  )
}

// --------------------------------------------------
// MethodologyPage
// --------------------------------------------------
export default function MethodologyPage({ onBack }) {
  const [openSection, setOpenSection] = useState(null)
  const navigate = useNavigate()
  const navigateBack = onBack || (() => navigate(-1))

  const toggleSection = (id) => {
    setOpenSection((prev) => (prev === id ? null : id))
  }

  return (
    <div className="space-y-6 max-w-4xl mx-auto">
      {/* Back + header */}
      <div className="flex items-center gap-4">
        <button
          onClick={navigateBack}
          className="flex items-center gap-1.5 px-3 py-2 bg-gray-100 hover:bg-gray-200 text-gray-700 text-sm rounded-lg transition-colors flex-shrink-0"
        >
          <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10 19l-7-7m0 0l7-7m-7 7h18" />
          </svg>
          Back
        </button>
        <div>
          <h1 className="text-2xl font-extrabold text-bx-light-text">Methodology</h1>
          <p className="text-sm text-gray-500 mt-0.5">
            Scientific description of each tool and algorithm used in the DockIt pipeline
          </p>
        </div>
      </div>

      {/* Expand all / Collapse all */}
      <div className="flex justify-end gap-2">
        <button
          onClick={() => setOpenSection('all')}
          className="text-xs text-bx-light-text hover:underline"
        >
          Expand all
        </button>
        <span className="text-gray-300">|</span>
        <button
          onClick={() => setOpenSection(null)}
          className="text-xs text-bx-light-text hover:underline"
        >
          Collapse all
        </button>
      </div>

      {/* Pipeline overview */}
      <div className="bg-bx-surface rounded-xl p-5 text-white">
        <h2 className="font-bold text-lg mb-2">DockIt Computational Pipeline</h2>
        <p className="text-white/70 text-sm leading-relaxed mb-4">
          DockIt automates a complete structure-based virtual screening and hit-to-lead workflow.
          Each step is performed by a best-in-class computational tool with fallback mechanisms to ensure robustness.
        </p>
        <div className="flex flex-wrap gap-2">
          {SECTIONS.map((s, i) => (
            <div key={s.id} className="flex items-center gap-1.5 text-xs text-white/80">
              {i > 0 && <svg className="w-3 h-3 text-white/30" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
              </svg>}
              <span className={`px-2 py-0.5 rounded-full text-xs font-semibold ${s.badgeColor}`}>{i + 1}</span>
              <span>{s.title}</span>
            </div>
          ))}
        </div>
      </div>

      {/* Section cards */}
      <div className="space-y-3">
        {SECTIONS.map((section) => (
          <SectionCard
            key={section.id}
            section={section}
            isOpen={openSection === 'all' || openSection === section.id}
            onToggle={() => toggleSection(section.id)}
          />
        ))}
      </div>

      {/* Docking Integrity & Validation */}
      <div className="card overflow-hidden">
        <div className="px-5 py-4 bg-bx-surface">
          <h2 className="font-bold text-lg text-white">Docking Integrity & Validation Rules</h2>
          <p className="text-white/70 text-sm mt-1">
            Principles ensuring scientific rigor and transparency in all docking results
          </p>
        </div>
        <div className="p-5 space-y-4">
          {INTEGRITY_RULES.map((item, idx) => (
            <div key={idx} className="flex items-start gap-3">
              <div className="w-8 h-8 rounded-lg bg-blue-50 flex items-center justify-center text-bx-light-text flex-shrink-0 mt-0.5">
                {item.icon}
              </div>
              <div>
                <p className="text-sm font-semibold text-bx-light-text">{item.rule}</p>
                <p className="text-xs text-gray-500 mt-0.5 leading-relaxed">{item.detail}</p>
              </div>
            </div>
          ))}
        </div>

        {/* Summary statement */}
        <div className="px-5 pb-5">
          <div className="bg-blue-50 rounded-lg p-4">
            <p className="text-sm text-bx-light-text font-medium italic leading-relaxed">
              "DockIt performs structure-based virtual screening using GNINA/Vina.
              Protein structures are prepared from experimental PDB or predicted AlphaFold data,
              binding pockets are defined prior to docking, and ligands are docked within a fixed search box.
              Docking results provide atom-level 3D poses directly from the docking engine, without post-processing or coordinate manipulation.
              Ligand positions are expected to lie within or near the binding cavity rather than at the geometric center of the pocket.
              When no docking pose is available, only 2D representations are shown to avoid misleading interpretations.
              Docking scores are used for relative ranking, not absolute affinity prediction."
            </p>
          </div>
        </div>
      </div>

      {/* Benchmark Validation Results */}
      <div className="card overflow-hidden">
        <div className="px-5 py-4 bg-green-700">
          <h2 className="font-bold text-lg text-white">Benchmark Validation Results</h2>
          <p className="text-white/70 text-sm mt-1">
            Five-target benchmark with known actives vs. decoys (February 2026)
          </p>
        </div>
        <div className="p-5 space-y-6">

          {/* Protocol summary */}
          <div className="bg-gray-50 rounded-lg p-4">
            <h3 className="text-sm font-bold text-gray-700 mb-2">Benchmark Protocol</h3>
            <p className="text-xs text-gray-600 leading-relaxed">
              Five kinase targets (EGFR, CDK2, BRAF V600E, JAK2, KRAS G12C) were screened using GNINA 1.1
              with experimental PDB structures. Each set included 9-13 known inhibitors (actives) and
              10 non-target drugs (decoys) to evaluate enrichment. Enrichment Factor (EF) is the ratio of
              active fraction in top-ranked compounds vs. the random expectation. EF10% counts actives
              in the top 10% of the ranked list. All receptor structures were prepared with automatic
              HETATM stripping. Ligands were provided as SMILES and converted to SDF for native GNINA input.
            </p>
          </div>

          {/* Enrichment factors table */}
          <div>
            <h3 className="text-sm font-bold text-gray-700 mb-2">Enrichment Factors — Actives vs. Decoys</h3>
            <div className="overflow-x-auto">
              <table className="w-full text-xs border-collapse">
                <thead>
                  <tr className="bg-bx-surface text-white">
                    <th className="px-3 py-2 text-left font-semibold rounded-tl-lg">Target</th>
                    <th className="px-3 py-2 text-center font-semibold">Actives</th>
                    <th className="px-3 py-2 text-center font-semibold">Decoys</th>
                    <th className="px-3 py-2 text-center font-semibold">EF Vina</th>
                    <th className="px-3 py-2 text-center font-semibold">EF CNN</th>
                    <th className="px-3 py-2 text-center font-semibold">EF10% Vina</th>
                    <th className="px-3 py-2 text-center font-semibold">EF10% CNN</th>
                    <th className="px-3 py-2 text-center font-semibold rounded-tr-lg">GNINA success</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-gray-100">
                  <tr className="hover:bg-gray-50">
                    <td className="px-3 py-2 font-semibold text-bx-light-text">EGFR (P00533)</td>
                    <td className="px-3 py-2 text-center text-gray-600">13</td>
                    <td className="px-3 py-2 text-center text-gray-600">10</td>
                    <td className="px-3 py-2 text-center"><span className="text-green-700 font-bold">1.05x</span></td>
                    <td className="px-3 py-2 text-center"><span className="text-green-700 font-bold">1.34x</span></td>
                    <td className="px-3 py-2 text-center text-gray-700">1.5</td>
                    <td className="px-3 py-2 text-center text-gray-700">1.5</td>
                    <td className="px-3 py-2 text-center"><span className="text-green-600 font-bold">100%</span></td>
                  </tr>
                  <tr className="bg-gray-50/40 hover:bg-gray-50">
                    <td className="px-3 py-2 font-semibold text-bx-light-text">CDK2 (P24941)</td>
                    <td className="px-3 py-2 text-center text-gray-600">11</td>
                    <td className="px-3 py-2 text-center text-gray-600">10</td>
                    <td className="px-3 py-2 text-center"><span className="text-green-700 font-bold">1.30x</span></td>
                    <td className="px-3 py-2 text-center"><span className="text-green-700 font-bold">1.57x</span></td>
                    <td className="px-3 py-2 text-center text-gray-700">0.9</td>
                    <td className="px-3 py-2 text-center text-gray-700">2.7</td>
                    <td className="px-3 py-2 text-center"><span className="text-green-600 font-bold">100%</span></td>
                  </tr>
                  <tr className="hover:bg-gray-50">
                    <td className="px-3 py-2 font-semibold text-bx-light-text">BRAF V600E (P15056)</td>
                    <td className="px-3 py-2 text-center text-gray-600">11</td>
                    <td className="px-3 py-2 text-center text-gray-600">10</td>
                    <td className="px-3 py-2 text-center"><span className="text-green-700 font-bold">1.61x</span></td>
                    <td className="px-3 py-2 text-center"><span className="text-green-700 font-bold">2.10x</span></td>
                    <td className="px-3 py-2 text-center text-gray-700">1.8</td>
                    <td className="px-3 py-2 text-center text-gray-700">2.7</td>
                    <td className="px-3 py-2 text-center"><span className="text-green-600 font-bold">100%</span></td>
                  </tr>
                  <tr className="bg-gray-50/40 hover:bg-gray-50">
                    <td className="px-3 py-2 font-semibold text-bx-light-text">JAK2 (O60674)</td>
                    <td className="px-3 py-2 text-center text-gray-600">12</td>
                    <td className="px-3 py-2 text-center text-gray-600">10</td>
                    <td className="px-3 py-2 text-center"><span className="text-green-700 font-bold">1.19x</span></td>
                    <td className="px-3 py-2 text-center"><span className="text-green-700 font-bold">1.36x</span></td>
                    <td className="px-3 py-2 text-center text-gray-700">0.8</td>
                    <td className="px-3 py-2 text-center text-gray-700">0.8</td>
                    <td className="px-3 py-2 text-center"><span className="text-green-600 font-bold">100%</span></td>
                  </tr>
                  <tr className="hover:bg-gray-50">
                    <td className="px-3 py-2 font-semibold text-bx-light-text">KRAS G12C (P01116)</td>
                    <td className="px-3 py-2 text-center text-gray-600">9</td>
                    <td className="px-3 py-2 text-center text-gray-600">10</td>
                    <td className="px-3 py-2 text-center"><span className="text-green-700 font-bold">1.59x</span></td>
                    <td className="px-3 py-2 text-center"><span className="text-green-700 font-bold">1.92x</span></td>
                    <td className="px-3 py-2 text-center text-gray-700">1.3</td>
                    <td className="px-3 py-2 text-center text-amber-600 font-semibold">0.0</td>
                    <td className="px-3 py-2 text-center"><span className="text-amber-600 font-bold">77%</span></td>
                  </tr>
                  {/* Summary row */}
                  <tr className="bg-green-50 border-t-2 border-green-200">
                    <td className="px-3 py-2 font-bold text-green-800">Mean (5 targets)</td>
                    <td className="px-3 py-2 text-center font-semibold text-green-800">11.2</td>
                    <td className="px-3 py-2 text-center font-semibold text-green-800">10</td>
                    <td className="px-3 py-2 text-center"><span className="text-green-800 font-bold">1.35x</span></td>
                    <td className="px-3 py-2 text-center"><span className="text-green-800 font-bold">1.66x</span></td>
                    <td className="px-3 py-2 text-center font-semibold text-green-800">1.26</td>
                    <td className="px-3 py-2 text-center font-semibold text-green-800">1.54</td>
                    <td className="px-3 py-2 text-center"><span className="text-green-800 font-bold">95%</span></td>
                  </tr>
                </tbody>
              </table>
            </div>
            <p className="text-xs text-gray-400 mt-1.5 italic">
              EF = mean active score / mean decoy score (overall enrichment). EF10% = actives recovered in top 10% of ranked list (expected = 1.0 at random).
            </p>
          </div>

          {/* Spearman correlation table */}
          <div>
            <h3 className="text-sm font-bold text-gray-700 mb-2">Rank Correlation with Experimental IC50</h3>
            <div className="overflow-x-auto">
              <table className="w-full text-xs border-collapse">
                <thead>
                  <tr className="bg-gray-100 text-gray-600 uppercase tracking-wide">
                    <th className="px-3 py-2 text-left font-semibold">Target</th>
                    <th className="px-3 py-2 text-center font-semibold">n actives</th>
                    <th className="px-3 py-2 text-center font-semibold">Vina rho</th>
                    <th className="px-3 py-2 text-center font-semibold">Vina p-value</th>
                    <th className="px-3 py-2 text-center font-semibold">CNN rho</th>
                    <th className="px-3 py-2 text-center font-semibold">CNN p-value</th>
                    <th className="px-3 py-2 text-center font-semibold">Significant?</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-gray-100">
                  <tr className="hover:bg-gray-50">
                    <td className="px-3 py-2 font-medium text-gray-700">EGFR</td>
                    <td className="px-3 py-2 text-center text-gray-600">13</td>
                    <td className="px-3 py-2 text-center text-gray-600">-0.265</td>
                    <td className="px-3 py-2 text-center text-gray-500">0.381</td>
                    <td className="px-3 py-2 text-center text-gray-600">0.221</td>
                    <td className="px-3 py-2 text-center text-gray-500">0.468</td>
                    <td className="px-3 py-2 text-center text-gray-400">No</td>
                  </tr>
                  <tr className="bg-gray-50/40 hover:bg-gray-50">
                    <td className="px-3 py-2 font-medium text-gray-700">CDK2</td>
                    <td className="px-3 py-2 text-center text-gray-600">11</td>
                    <td className="px-3 py-2 text-center text-gray-600">0.282</td>
                    <td className="px-3 py-2 text-center text-gray-500">0.401</td>
                    <td className="px-3 py-2 text-center text-gray-600">-0.082</td>
                    <td className="px-3 py-2 text-center text-gray-500">0.811</td>
                    <td className="px-3 py-2 text-center text-gray-400">No</td>
                  </tr>
                  <tr className="hover:bg-gray-50">
                    <td className="px-3 py-2 font-medium text-gray-700">BRAF V600E</td>
                    <td className="px-3 py-2 text-center text-gray-600">11</td>
                    <td className="px-3 py-2 text-center text-gray-600">0.077</td>
                    <td className="px-3 py-2 text-center text-gray-500">0.821</td>
                    <td className="px-3 py-2 text-center text-gray-600">0.296</td>
                    <td className="px-3 py-2 text-center text-gray-500">0.377</td>
                    <td className="px-3 py-2 text-center text-gray-400">No</td>
                  </tr>
                  <tr className="bg-gray-50/40 hover:bg-gray-50">
                    <td className="px-3 py-2 font-medium text-gray-700">JAK2</td>
                    <td className="px-3 py-2 text-center text-gray-600">12</td>
                    <td className="px-3 py-2 text-center text-gray-600">-0.141</td>
                    <td className="px-3 py-2 text-center text-gray-500">0.662</td>
                    <td className="px-3 py-2 text-center text-gray-600">-0.063</td>
                    <td className="px-3 py-2 text-center text-gray-500">0.845</td>
                    <td className="px-3 py-2 text-center text-gray-400">No</td>
                  </tr>
                  <tr className="hover:bg-gray-50">
                    <td className="px-3 py-2 font-medium text-gray-700">KRAS G12C</td>
                    <td className="px-3 py-2 text-center text-gray-600">9</td>
                    <td className="px-3 py-2 text-center font-semibold text-green-700">-0.817</td>
                    <td className="px-3 py-2 text-center">
                      <span className="bg-green-100 text-green-800 font-bold px-1.5 py-0.5 rounded">0.007*</span>
                    </td>
                    <td className="px-3 py-2 text-center text-gray-600">-0.400</td>
                    <td className="px-3 py-2 text-center text-gray-500">0.286</td>
                    <td className="px-3 py-2 text-center text-green-600 font-semibold">Vina only</td>
                  </tr>
                </tbody>
              </table>
            </div>
            <p className="text-xs text-gray-400 mt-1.5 italic">
              Spearman rank correlation between docking score and experimental IC50 (actives only, with IC50 data). * p &lt; 0.05.
              Negative rho for Vina is expected (more negative score = better rank, lower IC50 = more potent).
            </p>
          </div>

          {/* Reproducibility + GPU/CPU */}
          <div className="grid grid-cols-1 sm:grid-cols-2 gap-4">
            <div className="bg-green-50 rounded-lg p-4">
              <h4 className="text-xs font-bold text-green-800 uppercase tracking-wide mb-2">Reproducibility (CPU runs)</h4>
              <ul className="text-xs text-green-700 space-y-1.5">
                <li className="flex justify-between">
                  <span>Top-5 overlap between runs</span>
                  <span className="font-bold">100%</span>
                </li>
                <li className="flex justify-between">
                  <span>Top-10 overlap between runs</span>
                  <span className="font-bold">100%</span>
                </li>
                <li className="flex justify-between">
                  <span>Mean score difference</span>
                  <span className="font-bold">0.69 kcal/mol</span>
                </li>
              </ul>
              <p className="text-xs text-green-600 mt-2 leading-relaxed">
                GNINA uses deterministic seeds. Rankings are highly stable across repeated CPU runs.
              </p>
            </div>
            <div className="bg-amber-50 rounded-lg p-4">
              <h4 className="text-xs font-bold text-amber-800 uppercase tracking-wide mb-2">GPU vs. CPU Consistency</h4>
              <ul className="text-xs text-amber-700 space-y-1.5">
                <li className="flex justify-between">
                  <span>Rank correlation (rho)</span>
                  <span className="font-bold text-amber-600">-0.046 (poor)</span>
                </li>
                <li className="flex justify-between">
                  <span>Top-5 overlap</span>
                  <span className="font-bold text-amber-600">0%</span>
                </li>
                <li className="flex justify-between">
                  <span>Mean score difference</span>
                  <span className="font-bold text-amber-600">6.27 kcal/mol</span>
                </li>
              </ul>
              <p className="text-xs text-amber-600 mt-2 leading-relaxed">
                GPU and CPU GNINA inference modes differ substantially. DockIt runs CPU mode consistently
                to ensure reproducibility. GPU results should not be mixed with CPU results for ranking.
              </p>
            </div>
          </div>

          {/* Key findings */}
          <div className="grid grid-cols-1 sm:grid-cols-2 gap-4">
            <div className="bg-blue-50 rounded-lg p-3">
              <h4 className="text-xs font-bold text-blue-800 uppercase tracking-wide mb-1">Interpretation</h4>
              <ul className="text-xs text-blue-700 space-y-1 list-disc list-inside">
                <li>CNN affinity consistently outperforms Vina (mean EF 1.66x vs 1.35x across 5 targets)</li>
                <li>Rank correlations with IC50 are weak across most targets — expected for docking-based scoring</li>
                <li>KRAS significant correlation (rho=-0.817, p=0.007) likely reflects wide IC50 range in the active set</li>
                <li>KRAS EF10% CNN = 0.0 is consistent with lower GNINA success rate (77%) for this target</li>
                <li>All 4 targets with 100% GNINA success show positive enrichment in both Vina and CNN</li>
              </ul>
            </div>
            <div className="bg-amber-50 rounded-lg p-3">
              <h4 className="text-xs font-bold text-amber-800 uppercase tracking-wide mb-1">Limitations noted</h4>
              <ul className="text-xs text-amber-700 space-y-1 list-disc list-inside">
                <li>EF values (1.05-2.10x) are modest vs. DUD-E benchmarks (5-15x) — decoys here are real drugs, not random compounds</li>
                <li>Small active sets (9-13 per target) limit statistical power for rank correlation</li>
                <li>GPU/CPU inconsistency indicates CNN scores are hardware-dependent; use CPU for comparability</li>
                <li>Composite score includes druglikeness terms; use raw Vina/CNN for enrichment evaluation</li>
              </ul>
            </div>
          </div>

          {/* Scoring formula */}
          <div className="bg-blue-50 rounded-lg p-4">
            <h4 className="text-xs font-bold text-blue-800 uppercase tracking-wide mb-2">Composite Scoring Formula</h4>
            <p className="text-xs text-blue-700 leading-relaxed font-mono">
              Score = 0.65 x norm_affinity + 0.20 x QED + 0.15 x logP_penalty
            </p>
            <p className="text-xs text-blue-600 mt-2 leading-relaxed">
              Binding affinity is the dominant term (65% weight). QED and logP provide secondary
              druglikeness filtering. Affinity is normalized over [-14, 0] kcal/mol. For V2 mode with ADMET:
              Score = 0.55 x affinity + 0.20 x ADMET + 0.15 x QED + 0.10 x novelty.
              Hard cutoffs: Lipinski &gt;2 violations, QED &lt;0.25, SA &gt;6.0, PAINS alerts, hERG risk.
            </p>
          </div>
        </div>
      </div>

      {/* General disclaimer */}
      <div className="bg-amber-50 border border-amber-100 rounded-xl p-5">
        <div className="flex items-start gap-3">
          <svg className="w-5 h-5 text-amber-500 flex-shrink-0 mt-0.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
          </svg>
          <div>
            <p className="text-sm font-semibold text-amber-800 mb-1">General Disclaimer</p>
            <p className="text-xs text-amber-700 leading-relaxed">
              All predictions produced by DockIt are computational and exploratory in nature.
              They do not constitute medical advice and should not be used to make clinical decisions.
              All candidates must be synthesized and validated experimentally (biochemical assays, cell-based assays, in vivo studies)
              before any consideration for further development.
              Computational scores are indicative only and may not correlate with experimental activity.
              DockIt is a research tool for early-stage drug discovery support.
            </p>
          </div>
        </div>
      </div>
    </div>
  )
}
