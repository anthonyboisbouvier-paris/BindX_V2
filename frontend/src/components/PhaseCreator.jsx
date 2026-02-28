import React, { useEffect, useRef } from 'react'
import { PHASE_TYPES } from '../mock/data.js'

// ---------------------------------------------------------------------------
// Phase type details for the modal
// ---------------------------------------------------------------------------

const PHASE_DETAILS = {
  hit_discovery: {
    focuses: [
      'High-throughput virtual screening',
      'Docking-based prioritization',
      'Initial ADMET filtering',
      'Hit list generation',
    ],
    defaultColumns: ['name', 'docking_score', 'cnn_score', 'logP', 'MW', 'TPSA', 'lipinski_pass', 'composite_score'],
    colorClass: 'bg-blue-600',
    borderClass: 'border-blue-200',
    bgClass: 'bg-blue-50',
    textClass: 'text-blue-700',
  },
  hit_to_lead: {
    focuses: [
      'Full ADMET profiling',
      'Selectivity screening',
      'Interaction analysis',
      'Final lead selection',
    ],
    defaultColumns: ['name', 'docking_score', 'cnn_score', 'composite_score', 'generation_level', 'cluster_id', 'scaffold'],
    colorClass: 'bg-purple-600',
    borderClass: 'border-purple-200',
    bgClass: 'bg-purple-50',
    textClass: 'text-purple-700',
  },
  lead_optimization: {
    focuses: [
      'Full ADMET profiling',
      'Selectivity screening',
      'Interaction analysis',
      'Final lead selection',
    ],
    defaultColumns: ['name', 'composite_score', 'logP', 'solubility', 'BBB', 'hERG', 'metabolic_stability', 'interactions_count'],
    colorClass: 'bg-green-600',
    borderClass: 'border-green-200',
    bgClass: 'bg-green-50',
    textClass: 'text-green-700',
  },
}

// Phase sequence
const PHASE_SEQUENCE = ['hit_discovery', 'hit_to_lead', 'lead_optimization']

// ---------------------------------------------------------------------------
// PhaseCreator
// ---------------------------------------------------------------------------
// Props:
//   campaignId     — string
//   existingPhases — array of existing phase objects
//   onCreatePhase  — (phaseConfig) => void
//   onClose        — () => void

export default function PhaseCreator({ campaignId, existingPhases = [], onCreatePhase, onClose }) {
  const overlayRef = useRef(null)

  // Determine next phase type
  const existingTypes = new Set(existingPhases.filter(p => p.created_at).map(p => p.type))
  const nextType = PHASE_SEQUENCE.find(t => !existingTypes.has(t)) || 'lead_optimization'
  const phaseTypeMeta = PHASE_TYPES[nextType] || { label: 'Unknown Phase', short: '?', description: '' }
  const phaseDetails = PHASE_DETAILS[nextType] || PHASE_DETAILS.lead_optimization

  // Previous phase bookmarked count
  const prevPhaseIdx = PHASE_SEQUENCE.indexOf(nextType) - 1
  const prevType = prevPhaseIdx >= 0 ? PHASE_SEQUENCE[prevPhaseIdx] : null
  const prevPhase = prevType ? existingPhases.find(p => p.type === prevType && p.created_at) : null
  const availableMols = prevPhase?.stats?.bookmarked ?? 0

  // Close on escape
  useEffect(() => {
    const handler = (e) => { if (e.key === 'Escape') onClose?.() }
    document.addEventListener('keydown', handler)
    return () => document.removeEventListener('keydown', handler)
  }, [onClose])

  const handleCreate = () => {
    if (onCreatePhase) {
      onCreatePhase({
        campaign_id: campaignId,
        type: nextType,
        label: `Phase ${phaseTypeMeta.short}`,
        column_presets: phaseDetails.defaultColumns,
        created_at: new Date().toISOString(),
        status: 'active',
        runs: [],
        molecules: [],
        stats: { total_molecules: 0, bookmarked: 0, runs_completed: 0, runs_running: 0 },
      })
    }
    onClose?.()
  }

  return (
    /* Backdrop */
    <div
      ref={overlayRef}
      className="fixed inset-0 bg-black/40 backdrop-blur-sm flex items-center justify-center z-50 px-4"
      onClick={(e) => { if (e.target === overlayRef.current) onClose?.() }}
    >
      <div className="bg-white rounded-2xl shadow-2xl w-full max-w-md overflow-hidden">
        {/* Header */}
        <div className={`${phaseDetails.bgClass} ${phaseDetails.borderClass} border-b px-6 py-4 flex items-center justify-between`}>
          <div className="flex items-center gap-3">
            <span className={`w-8 h-8 rounded-lg ${phaseDetails.colorClass} text-white flex items-center justify-center font-bold text-sm`}>
              {phaseTypeMeta.short}
            </span>
            <div>
              <h3 className="font-semibold text-gray-800 text-sm">Create New Phase</h3>
              <p className={`text-xs font-medium ${phaseDetails.textClass}`}>
                Phase {phaseTypeMeta.short} — {phaseTypeMeta.label}
              </p>
            </div>
          </div>
          <button
            onClick={onClose}
            className="w-7 h-7 rounded-full bg-white/50 hover:bg-white flex items-center justify-center transition-colors text-gray-500"
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
            </svg>
          </button>
        </div>

        {/* Body */}
        <div className="px-6 py-5 space-y-5">
          {/* Description */}
          <div>
            <p className="text-xs font-semibold text-gray-500 uppercase tracking-wider mb-2">
              This phase focuses on
            </p>
            <ul className="space-y-1.5">
              {phaseDetails.focuses.map((item, i) => (
                <li key={i} className="flex items-start gap-2 text-sm text-gray-600">
                  <svg className="w-4 h-4 text-bx-mint shrink-0 mt-0.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 13l4 4L19 7" />
                  </svg>
                  {item}
                </li>
              ))}
            </ul>
          </div>

          {/* Input source */}
          {prevPhase && (
            <div className={`${phaseDetails.bgClass} ${phaseDetails.borderClass} border rounded-lg px-4 py-3`}>
              <p className="text-xs font-semibold text-gray-500 uppercase tracking-wider mb-1">Input</p>
              <p className="text-sm text-gray-700">
                Bookmarked molecules from Phase {PHASE_TYPES[prevType]?.short}
              </p>
              <p className={`text-xs ${phaseDetails.textClass} font-semibold mt-0.5`}>
                {availableMols} molecule{availableMols !== 1 ? 's' : ''} available
              </p>
            </div>
          )}

          {/* Default columns */}
          <div>
            <p className="text-xs font-semibold text-gray-500 uppercase tracking-wider mb-2">Default columns</p>
            <div className="flex flex-wrap gap-1.5">
              {phaseDetails.defaultColumns.map(col => (
                <span
                  key={col}
                  className="text-xs font-mono px-2 py-0.5 bg-gray-100 text-gray-600 rounded"
                >
                  {col}
                </span>
              ))}
            </div>
          </div>
        </div>

        {/* Footer */}
        <div className="px-6 py-4 bg-gray-50 border-t border-gray-100 flex items-center justify-end gap-3">
          <button
            onClick={onClose}
            className="px-4 py-2 text-sm text-gray-500 hover:text-gray-700 rounded-lg hover:bg-gray-100 transition-colors"
          >
            Cancel
          </button>
          <button
            onClick={handleCreate}
            className={`px-5 py-2 ${phaseDetails.colorClass} hover:opacity-90 text-white text-sm font-semibold rounded-lg transition-opacity flex items-center gap-2`}
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 4v16m8-8H4" />
            </svg>
            Create Phase {phaseTypeMeta.short}
          </button>
        </div>
      </div>
    </div>
  )
}
