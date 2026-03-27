import React, { useState, useCallback } from 'react'
import { v9GenerateReport } from '../api.js'

const SECTIONS = [
  { id: 'summary', label: 'Summary', desc: 'Project info, molecule counts, run stats' },
  { id: 'top_molecules', label: 'Top Molecules', desc: 'Best molecules ranked by composite score' },
  { id: 'property_distributions', label: 'Property Distributions', desc: 'Histograms of key properties (MW, logP, QED...)' },
  { id: 'admet_profiles', label: 'ADMET Profiles', desc: 'Absorption, distribution, metabolism, excretion data' },
  { id: 'safety_overview', label: 'Safety Overview', desc: 'Safety color distribution, PAINS/Brenk alerts' },
  { id: 'pareto_analysis', label: 'Pareto Analysis', desc: 'Multi-objective optimization summary' },
  { id: 'retrosynthesis', label: 'Retrosynthesis', desc: 'Synthesis routes, steps, costs for top molecules' },
  { id: 'scaffold', label: 'Structural Analysis', desc: 'Scaffold decomposition + diversity clustering' },
]

const TOP_N_OPTIONS = [5, 10, 20, 50]
const MAX_MOLECULES = 200

export default function ReportBuilderModal({ isOpen, onClose, phaseId, selectedMoleculeIds, onToast }) {
  const [selectedSections, setSelectedSections] = useState(new Set(SECTIONS.map(s => s.id)))
  const [topN, setTopN] = useState(10)
  const [loading, setLoading] = useState(false)

  const selCount = selectedMoleculeIds?.size || 0
  const useSelection = selCount > 0
  const tooMany = selCount > MAX_MOLECULES

  const toggleSection = useCallback((id) => {
    setSelectedSections(prev => {
      const next = new Set(prev)
      if (next.has(id)) next.delete(id)
      else next.add(id)
      return next
    })
  }, [])

  const toggleAll = useCallback(() => {
    setSelectedSections(prev => {
      if (prev.size === SECTIONS.length) return new Set()
      return new Set(SECTIONS.map(s => s.id))
    })
  }, [])

  const handleGenerate = useCallback(async () => {
    if (!phaseId || selectedSections.size === 0) return
    if (tooMany) return
    setLoading(true)
    try {
      const opts = {
        sections: [...selectedSections],
        top_n: topN,
        format: 'pdf',
      }
      if (useSelection) {
        opts.molecule_ids = [...selectedMoleculeIds]
      }
      const blob = await v9GenerateReport(phaseId, opts)
      const url = URL.createObjectURL(blob)
      const a = document.createElement('a')
      a.href = url
      a.download = `BindX_Report_${new Date().toISOString().slice(0, 10)}.pdf`
      document.body.appendChild(a)
      a.click()
      document.body.removeChild(a)
      URL.revokeObjectURL(url)
      onToast?.('Report downloaded successfully', 'success')
      onClose()
    } catch (err) {
      onToast?.(err.userMessage || 'Failed to generate report', 'error')
    } finally {
      setLoading(false)
    }
  }, [phaseId, selectedSections, topN, useSelection, selectedMoleculeIds, tooMany, onClose, onToast])

  if (!isOpen) return null

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center">
      <div className="absolute inset-0 bg-black/40 backdrop-blur-sm" onClick={onClose} />

      <div className="relative bg-white rounded-2xl shadow-2xl w-full max-w-lg mx-4 overflow-hidden">
        {/* Header */}
        <div className="px-6 py-4 border-b border-gray-100 flex items-center justify-between">
          <div>
            <h2 className="text-lg font-bold text-gray-800">Generate Report</h2>
            <p className="text-xs text-gray-400 mt-0.5">Professional PDF report for this phase</p>
          </div>
          <button onClick={onClose} className="text-gray-400 hover:text-gray-600 transition-colors">
            <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
            </svg>
          </button>
        </div>

        {/* Body */}
        <div className="px-6 py-4 space-y-4 max-h-[60vh] overflow-y-auto">
          {/* Scope indicator */}
          <div className={`rounded-lg p-3 text-xs flex items-center gap-2 ${
            tooMany ? 'bg-red-50 border border-red-200' : useSelection ? 'bg-blue-50 border border-blue-200' : 'bg-gray-50 border border-gray-200'
          }`}>
            <svg className="w-4 h-4 shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                d={useSelection
                  ? "M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2m-6 9l2 2 4-4"
                  : "M19 11H5m14 0a2 2 0 012 2v6a2 2 0 01-2 2H5a2 2 0 01-2-2v-6a2 2 0 012-2m14 0V9a2 2 0 00-2-2M5 11V9a2 2 0 012-2m0 0V5a2 2 0 012-2h6a2 2 0 012 2v2M7 7h10"
                } />
            </svg>
            {tooMany ? (
              <span className="text-red-600 font-medium">
                {selCount} molecules selected — max {MAX_MOLECULES}. Reduce your selection.
              </span>
            ) : useSelection ? (
              <span className="text-blue-700 font-medium">
                Report on <span className="font-bold">{selCount}</span> selected molecule{selCount > 1 ? 's' : ''}
              </span>
            ) : (
              <span className="text-gray-500">Report on <span className="font-medium">all molecules</span> in this phase</span>
            )}
          </div>

          {/* Section Checklist */}
          <div>
            <div className="flex items-center justify-between mb-2">
              <label className="text-sm font-semibold text-gray-700">Sections</label>
              <button onClick={toggleAll} className="text-xs text-blue-500 hover:text-blue-700">
                {selectedSections.size === SECTIONS.length ? 'Deselect all' : 'Select all'}
              </button>
            </div>
            <div className="grid grid-cols-2 gap-2">
              {SECTIONS.map(section => (
                <label
                  key={section.id}
                  className={`flex items-start gap-2 p-2.5 rounded-lg border cursor-pointer transition-colors ${
                    selectedSections.has(section.id)
                      ? 'border-blue-300 bg-blue-50/50'
                      : 'border-gray-200 hover:border-gray-300'
                  }`}
                >
                  <input
                    type="checkbox"
                    checked={selectedSections.has(section.id)}
                    onChange={() => toggleSection(section.id)}
                    className="mt-0.5 rounded border-gray-300 text-blue-500 focus:ring-blue-500"
                  />
                  <div>
                    <p className="text-xs font-medium text-gray-700">{section.label}</p>
                    <p className="text-[10px] text-gray-400 mt-0.5 leading-tight">{section.desc}</p>
                  </div>
                </label>
              ))}
            </div>
          </div>

          {/* Top N Slider */}
          <div>
            <label className="text-sm font-semibold text-gray-700 block mb-2">
              Top N molecules: <span className="text-blue-500">{topN}</span>
            </label>
            <div className="flex items-center gap-3">
              <input
                type="range"
                min={5}
                max={50}
                step={5}
                value={topN}
                onChange={e => setTopN(Number(e.target.value))}
                className="flex-1 accent-blue-500"
              />
              <div className="flex gap-1">
                {TOP_N_OPTIONS.map(n => (
                  <button
                    key={n}
                    onClick={() => setTopN(n)}
                    className={`px-2 py-0.5 text-xs rounded ${
                      topN === n ? 'bg-blue-500 text-white' : 'bg-gray-100 text-gray-600 hover:bg-gray-200'
                    }`}
                  >
                    {n}
                  </button>
                ))}
              </div>
            </div>
          </div>

          {/* Preview */}
          <div className="bg-gray-50 rounded-lg p-3">
            <p className="text-xs text-gray-500">
              <span className="font-medium">{selectedSections.size}</span> sections selected —
              Top <span className="font-medium">{topN}</span> molecules by composite score
            </p>
          </div>
        </div>

        {/* Footer */}
        <div className="px-6 py-4 border-t border-gray-100 flex items-center justify-end gap-3">
          <button
            onClick={onClose}
            className="px-4 py-2 text-sm font-medium text-gray-600 hover:text-gray-800 transition-colors"
          >
            Cancel
          </button>
          <button
            onClick={handleGenerate}
            disabled={loading || selectedSections.size === 0 || tooMany}
            className="flex items-center gap-2 px-5 py-2 bg-bx-surface text-white rounded-lg text-sm font-semibold hover:bg-bx-elevated transition-colors disabled:opacity-50 disabled:cursor-not-allowed"
          >
            {loading ? (
              <>
                <div className="w-4 h-4 border-2 border-white/30 border-t-white rounded-full animate-spin" />
                Generating...
              </>
            ) : (
              <>
                <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                    d="M12 10v6m0 0l-3-3m3 3l3-3m2 8H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
                </svg>
                Generate PDF
              </>
            )}
          </button>
        </div>
      </div>
    </div>
  )
}
