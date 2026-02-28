import React from 'react'
import InfoTip from './InfoTip.jsx'

// --------------------------------------------------
// Component icons (SVG paths)
// --------------------------------------------------
const COMPONENT_ICONS = {
  structure: (
    <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
        d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
    </svg>
  ),
  pocket: (
    <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
        d="M9 20l-5.447-2.724A1 1 0 013 16.382V5.618a1 1 0 011.447-.894L9 7m0 13l6-3m-6 3V7m6 10l4.553 2.276A1 1 0 0021 18.382V7.618a1 1 0 00-.553-.894L15 4m0 13V4m0 0L9 7" />
    </svg>
  ),
  docking: (
    <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
        d="M13.828 10.172a4 4 0 00-5.656 0l-4 4a4 4 0 105.656 5.656l1.102-1.101m-.758-4.899a4 4 0 005.656 0l4-4a4 4 0 00-5.656-5.656l-1.1 1.1" />
    </svg>
  ),
  admet: (
    <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
        d="M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 01-2 2h-2a2 2 0 01-2-2z" />
    </svg>
  ),
  synthesis: (
    <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
        d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
    </svg>
  ),
}

const COMPONENT_LABELS = {
  structure: 'Structure Prediction',
  pocket: 'Pocket Detection',
  docking: 'Molecular Docking',
  admet: 'ADMET Prediction',
  synthesis: 'Retrosynthesis',
}

const COMPONENT_TIPS = {
  structure: 'Confidence in the predicted 3D protein structure. High for AlphaFold structures with pLDDT > 90.',
  pocket: 'Confidence that the binding pocket detected by fpocket is druggable and well-defined.',
  docking: 'Confidence in the docking pose. Lower for AI-generated molecules without experimental validation data.',
  admet: 'Confidence in ADMET predictions (absorption, distribution, metabolism, excretion, toxicity). Based on ADMET-AI model performance (typical AUROC 0.80-0.90).',
  synthesis: 'Confidence that a feasible synthetic route exists. Based on AiZynthFinder retrosynthesis analysis.',
}

// --------------------------------------------------
// Structure source hierarchy definition
// Priority order: pdb_experimental > alphafold > esmfold > homology
// --------------------------------------------------
const SOURCE_HIERARCHY = [
  {
    key: 'pdb_experimental',
    label: 'PDB Experimental',
    description: 'Crystal or cryo-EM structure from the Protein Data Bank',
    confidenceRange: '95-99%',
  },
  {
    key: 'alphafold',
    label: 'AlphaFold 2',
    description: 'AI prediction by DeepMind — pLDDT confidence per residue',
    confidenceRange: '70-90%',
  },
  {
    key: 'esmfold',
    label: 'ESMFold',
    description: 'Language-model prediction by Meta AI — fallback if AlphaFold unavailable',
    confidenceRange: '60-80%',
  },
  {
    key: 'homology',
    label: 'Homology Model',
    description: 'Template-based model — lowest confidence, used as last resort',
    confidenceRange: '40-65%',
  },
]

// Determine status of each source given the active source
function getSourceStatus(sourceKey, activeSource) {
  const activeIndex = SOURCE_HIERARCHY.findIndex((s) => s.key === activeSource)
  const thisIndex = SOURCE_HIERARCHY.findIndex((s) => s.key === sourceKey)

  if (activeIndex === -1) {
    // Unknown active source — mark all as unknown
    return 'unknown'
  }
  if (thisIndex === activeIndex) return 'used'
  if (thisIndex < activeIndex) return 'not_tried' // better source — not needed
  return 'failed'                                  // worse source — tried before this was found
}

// --------------------------------------------------
// Score color helper
// --------------------------------------------------
function scoreColor(val) {
  if (val >= 0.7) return '#00e6a0'   // green
  if (val >= 0.5) return '#eab308'   // yellow
  return '#ef4444'                    // red
}

function scoreTextClass(val) {
  if (val >= 0.7) return 'text-green-600'
  if (val >= 0.5) return 'text-yellow-500'
  return 'text-red-500'
}

function overallTextClass(val) {
  if (val >= 0.7) return 'text-dockit-green'
  if (val >= 0.5) return 'text-yellow-500'
  return 'text-red-500'
}

// --------------------------------------------------
// Status icon for source hierarchy
// --------------------------------------------------
function SourceStatusIcon({ status }) {
  if (status === 'used') {
    return (
      <span className="flex-shrink-0 w-5 h-5 rounded-full bg-green-100 flex items-center justify-center">
        <svg className="w-3 h-3 text-green-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M5 13l4 4L19 7" />
        </svg>
      </span>
    )
  }
  if (status === 'failed') {
    return (
      <span className="flex-shrink-0 w-5 h-5 rounded-full bg-red-100 flex items-center justify-center">
        <svg className="w-3 h-3 text-red-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M6 18L18 6M6 6l12 12" />
        </svg>
      </span>
    )
  }
  if (status === 'not_tried') {
    return (
      <span className="flex-shrink-0 w-5 h-5 rounded-full bg-gray-100 flex items-center justify-center">
        <span className="text-gray-400 text-sm font-bold leading-none">—</span>
      </span>
    )
  }
  // unknown
  return (
    <span className="flex-shrink-0 w-5 h-5 rounded-full bg-gray-100 flex items-center justify-center">
      <span className="text-gray-300 text-xs font-bold leading-none">?</span>
    </span>
  )
}

// --------------------------------------------------
// Structure Source Hierarchy section
// --------------------------------------------------
function StructureSourceHierarchy({ pipeline_summary }) {
  if (!pipeline_summary) return null

  const structureSource = pipeline_summary?.effective_structure_source || pipeline_summary?.structure_source || null
  const pdbInfo = pipeline_summary.pdb_info || {}
  const pdbId = pdbInfo.pdb_id || pipeline_summary.pdb_id || null
  const pdbResolution = pdbInfo.resolution || pipeline_summary.pdb_resolution || null

  // Resolve the normalised active source key
  let activeSource = structureSource
  if (structureSource === 'pdb' || structureSource === 'pdb_experimental') activeSource = 'pdb_experimental'
  else if (structureSource === 'alphafold' || structureSource === 'alphafold2') activeSource = 'alphafold'
  else if (structureSource === 'esmfold') activeSource = 'esmfold'
  else if (structureSource === 'homology' || structureSource === 'homology_model') activeSource = 'homology'

  // Disorder info
  const disorderInfo = pipeline_summary.disorder_info ?? pipeline_summary.disordered_regions ?? null
  const disorderPct = (() => {
    if (!disorderInfo) return null
    if (typeof disorderInfo === 'number') return Math.round(disorderInfo * 100)
    if (typeof disorderInfo.fraction === 'number') return Math.round(disorderInfo.fraction * 100)
    if (typeof disorderInfo.percent === 'number') return Math.round(disorderInfo.percent)
    return null
  })()

  return (
    <div className="px-5 pt-4 pb-2">
      <p className="text-xs font-semibold text-gray-400 uppercase tracking-wide mb-3">
        Structure Source Hierarchy
        <InfoTip text="Shows which protein structure source was used. Sources are ordered by reliability: experimental data is most accurate, computational predictions are fallbacks." />
      </p>

      <div className="space-y-1.5">
        {SOURCE_HIERARCHY.map((src) => {
          const status = getSourceStatus(src.key, activeSource)
          const isUsed = status === 'used'
          const isFailed = status === 'failed'

          return (
            <div
              key={src.key}
              className={`flex items-center gap-3 px-3 py-2 rounded-lg ${
                isUsed
                  ? 'bg-green-50 border border-green-200'
                  : isFailed
                  ? 'bg-red-50/50 border border-red-100'
                  : 'bg-gray-50 border border-gray-100'
              }`}
            >
              <SourceStatusIcon status={status} />

              <div className="flex-1 min-w-0">
                <div className="flex items-center gap-2">
                  <span
                    className={`text-xs font-semibold ${
                      isUsed ? 'text-green-800' : isFailed ? 'text-red-600' : 'text-gray-400'
                    }`}
                  >
                    {src.label}
                  </span>
                  {isUsed && pdbId && src.key === 'pdb_experimental' && (
                    <span className="text-xs font-mono text-green-700 bg-green-100 px-1.5 py-0.5 rounded">
                      {pdbId}
                    </span>
                  )}
                  {isUsed && pdbResolution && src.key === 'pdb_experimental' && (
                    <span className="text-xs text-green-600">{pdbResolution} A</span>
                  )}
                </div>
                <p className={`text-xs leading-tight mt-0.5 ${
                  isUsed ? 'text-green-700' : isFailed ? 'text-red-500' : 'text-gray-400'
                }`}>
                  {isFailed ? 'Not available for this protein' : src.description}
                </p>
              </div>

              <span
                className={`flex-shrink-0 text-xs font-mono ${
                  isUsed ? 'text-green-600' : 'text-gray-300'
                }`}
              >
                {src.confidenceRange}
              </span>
            </div>
          )
        })}
      </div>

      {/* Disorder warning */}
      {disorderPct !== null && disorderPct > 0 && (
        <div className="mt-3 flex items-start gap-2 p-2.5 bg-yellow-50 border border-yellow-200 rounded-lg">
          <svg className="w-4 h-4 text-yellow-500 flex-shrink-0 mt-0.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
          </svg>
          <p className="text-xs text-yellow-700 leading-relaxed">
            {disorderPct}% of residues are in disordered regions. Predictions in these areas may be less reliable.
          </p>
        </div>
      )}
    </div>
  )
}

// --------------------------------------------------
// Progress bar for one component
// --------------------------------------------------
function ComponentRow({ name, data }) {
  const pct = Math.round((data.score ?? 0) * 100)
  const color = scoreColor(data.score ?? 0)
  const icon = COMPONENT_ICONS[name]
  const label = COMPONENT_LABELS[name] || name
  const tip = COMPONENT_TIPS[name] || ''

  return (
    <div className="flex items-start gap-3 py-3 border-b border-gray-50 last:border-0">
      {/* Icon */}
      <div className="flex-shrink-0 w-9 h-9 rounded-lg bg-gray-50 flex items-center justify-center text-gray-400">
        {icon}
      </div>

      {/* Content */}
      <div className="flex-1 min-w-0">
        <div className="flex items-center justify-between mb-1">
          <span className="text-sm font-medium text-gray-700 flex items-center">
            {label}
            {tip && <InfoTip text={tip} />}
          </span>
          <span className={`text-sm font-bold ${scoreTextClass(data.score ?? 0)}`}>
            {pct}%
          </span>
        </div>
        {/* Progress bar */}
        <div className="w-full h-2 bg-gray-100 rounded-full overflow-hidden">
          <div
            className="h-2 rounded-full transition-all duration-500"
            style={{ width: `${pct}%`, backgroundColor: color }}
          />
        </div>
        {/* Note */}
        {data.note && (
          <p className="text-xs text-gray-400 mt-1 leading-relaxed">{data.note}</p>
        )}
      </div>
    </div>
  )
}

// --------------------------------------------------
// ConfidenceBreakdown
// --------------------------------------------------
export default function ConfidenceBreakdown({ confidence, moleculeName, pipeline_summary }) {
  if (!confidence) {
    return (
      <div className="bg-white rounded-xl border border-gray-100 shadow-sm p-6 text-center text-gray-400 text-sm">
        No confidence data available.
      </div>
    )
  }

  const overall = confidence.overall ?? 0
  const overallPct = Math.round(overall * 100)
  const components = confidence.components || {}

  // Find lowest component
  const componentEntries = Object.entries(components)
  const lowestEntry = componentEntries.length > 0
    ? componentEntries.reduce((min, entry) => (entry[1].score ?? 1) < (min[1].score ?? 1) ? entry : min)
    : null

  return (
    <div className="bg-white rounded-xl border border-gray-100 shadow-sm overflow-hidden">
      {/* Header */}
      <div className="bg-dockit-blue px-5 py-4 text-white">
        <h3 className="font-bold text-base">
          Confidence Breakdown
          {moleculeName && (
            <span className="text-dockit-green ml-2 font-mono">{moleculeName}</span>
          )}
        </h3>
        <p className="text-white/60 text-xs mt-0.5">
          Multi-component reliability assessment
          <InfoTip text="The overall confidence score aggregates uncertainty from each pipeline step. A high score means the prediction is well-supported across all computational tools." />
        </p>
      </div>

      {/* Structure Source Hierarchy — shown first */}
      <StructureSourceHierarchy pipeline_summary={pipeline_summary} />

      {/* Divider */}
      {pipeline_summary && <div className="mx-5 border-t border-gray-100" />}

      {/* Overall score */}
      <div className="flex items-center gap-6 px-5 py-5 border-b border-gray-100">
        {/* Circular score display */}
        <div className="flex-shrink-0 relative w-20 h-20">
          <svg viewBox="0 0 80 80" className="w-20 h-20 -rotate-90">
            <circle
              cx="40" cy="40" r="32"
              fill="none"
              stroke="#f3f4f6"
              strokeWidth="10"
            />
            <circle
              cx="40" cy="40" r="32"
              fill="none"
              stroke={scoreColor(overall)}
              strokeWidth="10"
              strokeDasharray={`${2 * Math.PI * 32}`}
              strokeDashoffset={`${2 * Math.PI * 32 * (1 - overall)}`}
              strokeLinecap="round"
              style={{ transition: 'stroke-dashoffset 0.6s ease' }}
            />
          </svg>
          <div className="absolute inset-0 flex flex-col items-center justify-center">
            <span className={`text-xl font-extrabold leading-none ${overallTextClass(overall)}`}>
              {overallPct}%
            </span>
          </div>
        </div>

        <div>
          <p className="text-sm font-semibold text-gray-700">Overall Confidence</p>
          <p className={`text-xs mt-1 font-medium ${overallTextClass(overall)}`}>
            {overall >= 0.7 ? 'High confidence' : overall >= 0.5 ? 'Moderate confidence' : 'Low confidence — review carefully'}
          </p>
          <p className="text-xs text-gray-400 mt-1 leading-relaxed max-w-xs">
            {overall >= 0.7
              ? 'Prediction is well-supported across all pipeline components.'
              : overall >= 0.5
              ? 'Some components introduce uncertainty. Experimental validation recommended.'
              : 'Multiple pipeline components have low confidence. Treat results as exploratory.'}
          </p>
        </div>
      </div>

      {/* Component breakdown */}
      <div className="px-5 py-2">
        <p className="text-xs font-semibold text-gray-400 uppercase tracking-wide py-2">Component Scores</p>
        {componentEntries.length === 0 ? (
          <p className="text-xs text-gray-400 italic py-4 text-center">No component data available.</p>
        ) : (
          componentEntries.map(([name, data]) => (
            <ComponentRow key={name} name={name} data={data} />
          ))
        )}
      </div>

      {/* Key limitation callout */}
      {lowestEntry && (lowestEntry[1].score ?? 1) < 0.7 && (
        <div className="mx-5 mb-4 p-3 bg-amber-50 border border-amber-100 rounded-lg">
          <div className="flex items-start gap-2">
            <svg className="w-4 h-4 text-amber-500 flex-shrink-0 mt-0.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
            </svg>
            <div>
              <p className="text-xs font-semibold text-amber-700">
                Key limitation: {COMPONENT_LABELS[lowestEntry[0]] || lowestEntry[0]}
              </p>
              <p className="text-xs text-amber-600 mt-0.5 leading-relaxed">
                {lowestEntry[1].note || `The ${COMPONENT_LABELS[lowestEntry[0]] || lowestEntry[0]} step has the lowest confidence (${Math.round((lowestEntry[1].score ?? 0) * 100)}%). This is the primary source of uncertainty in this prediction.`}
              </p>
            </div>
          </div>
        </div>
      )}
    </div>
  )
}
