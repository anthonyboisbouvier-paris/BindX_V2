import React, { useState, useRef, useEffect, useCallback } from 'react'

function SmilesPreview({ smiles, width = 80, height = 55 }) {
  const canvasRef = useRef(null)
  const [failed, setFailed] = useState(false)
  useEffect(() => {
    if (!smiles || !canvasRef.current) return
    if (smiles.length < 2) { setFailed(true); return }
    let cancelled = false
    setFailed(false)
    import('smiles-drawer').then(mod => {
      if (cancelled) return
      const SmilesDrawer = mod.default || mod
      const drawer = new SmilesDrawer.Drawer({ width, height, bondThickness: 1 })
      SmilesDrawer.parse(smiles, (tree) => {
        if (!cancelled) drawer.draw(tree, canvasRef.current, 'light')
      }, () => { if (!cancelled) setFailed(true) })
    }).catch(() => { if (!cancelled) setFailed(true) })
    return () => { cancelled = true }
  }, [smiles, width, height])
  if (failed) return <span className="text-[10px] font-mono text-gray-400 px-1">{smiles?.slice(0, 20)}</span>
  return <canvas ref={canvasRef} width={width} height={height} className="bg-white rounded border border-gray-100" />
}

function formatVal(v, d = 2) {
  if (v == null) return '\u2014'
  return typeof v === 'number' ? v.toFixed(d) : String(v)
}

/**
 * Modal to review and import MMP-generated suggestions.
 */
export default function MMPSuggestionsModal({
  open,
  onClose,
  suggestions,
  transform,
  nInputMolecules,
  loading,
  onImport,
  importing,
}) {
  const [selected, setSelected] = useState(() => new Set())
  const allIds = suggestions?.map((_, i) => i) || []

  const toggleAll = useCallback(() => {
    if (selected.size === allIds.length) setSelected(new Set())
    else setSelected(new Set(allIds))
  }, [selected.size, allIds.length])

  const toggleOne = useCallback((idx) => {
    setSelected(prev => {
      const next = new Set(prev)
      if (next.has(idx)) next.delete(idx)
      else next.add(idx)
      return next
    })
  }, [])

  // Select all by default when suggestions arrive
  useEffect(() => {
    if (suggestions?.length) setSelected(new Set(suggestions.map((_, i) => i)))
  }, [suggestions])

  if (!open) return null

  const selectedSuggestions = suggestions?.filter((_, i) => selected.has(i)) || []

  const handleImport = () => {
    if (!selectedSuggestions.length) return
    onImport(selectedSuggestions)
  }

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/40 backdrop-blur-sm" onClick={onClose}>
      <div
        className="bg-white rounded-xl shadow-2xl w-[900px] max-w-[95vw] max-h-[85vh] flex flex-col"
        onClick={e => e.stopPropagation()}
      >
        {/* Header */}
        <div className="px-6 py-4 border-b border-gray-100 flex items-center justify-between">
          <div>
            <h2 className="text-lg font-bold text-bx-light-text flex items-center gap-2">
              <svg className="w-5 h-5 text-violet-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M12 4v16m8-8H4" />
              </svg>
              MMP Suggestions
            </h2>
            <p className="text-xs text-gray-400 mt-0.5">
              Transform: <span className="font-mono text-gray-500">{transform}</span>
              {nInputMolecules != null && <span className="ml-2">({nInputMolecules} input molecules)</span>}
            </p>
          </div>
          <button onClick={onClose} className="p-1.5 rounded-lg hover:bg-gray-100 transition text-gray-400 hover:text-gray-600">
            <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
            </svg>
          </button>
        </div>

        {/* Body */}
        <div className="flex-1 overflow-auto px-6 py-3">
          {loading ? (
            <div className="flex flex-col items-center justify-center py-16 gap-3">
              <div className="w-8 h-8 border-2 border-violet-400 border-t-transparent rounded-full animate-spin" />
              <p className="text-sm text-gray-400">Generating suggestions...</p>
            </div>
          ) : !suggestions?.length ? (
            <div className="flex flex-col items-center justify-center py-16 gap-3 text-gray-400">
              <svg className="w-10 h-10" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M9.172 16.172a4 4 0 015.656 0M9 10h.01M15 10h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
              </svg>
              <p className="text-sm">No valid suggestions generated for this transformation.</p>
            </div>
          ) : (
            <table className="w-full text-xs border-collapse">
              <thead className="sticky top-0 bg-white z-10">
                <tr className="border-b border-gray-200">
                  <th className="p-2 w-8">
                    <input
                      type="checkbox"
                      checked={selected.size === allIds.length && allIds.length > 0}
                      onChange={toggleAll}
                      className="rounded border-gray-300"
                    />
                  </th>
                  <th className="p-2 text-left font-medium text-gray-500">Structure</th>
                  <th className="p-2 text-left font-medium text-gray-500">Source</th>
                  <th className="p-2 text-left font-medium text-gray-500 font-mono">SMILES</th>
                  <th className="p-2 text-right font-medium text-gray-500">QED</th>
                  <th className="p-2 text-right font-medium text-gray-500">MW</th>
                  <th className="p-2 text-right font-medium text-gray-500">logP</th>
                  <th className="p-2 text-right font-medium text-gray-500">HBD</th>
                  <th className="p-2 text-right font-medium text-gray-500">HBA</th>
                </tr>
              </thead>
              <tbody>
                {suggestions.map((s, i) => (
                  <tr
                    key={i}
                    className={`border-b border-gray-100 hover:bg-violet-50/30 cursor-pointer ${selected.has(i) ? 'bg-violet-50/20' : ''}`}
                    onClick={() => toggleOne(i)}
                  >
                    <td className="p-2 text-center" onClick={e => e.stopPropagation()}>
                      <input
                        type="checkbox"
                        checked={selected.has(i)}
                        onChange={() => toggleOne(i)}
                        className="rounded border-gray-300"
                      />
                    </td>
                    <td className="p-2">
                      <SmilesPreview smiles={s.smiles} />
                    </td>
                    <td className="p-2 text-gray-500 truncate max-w-[100px]" title={s.source_mol_name}>
                      {s.source_mol_name || '\u2014'}
                    </td>
                    <td className="p-2 font-mono text-[10px] text-gray-500 truncate max-w-[180px]" title={s.smiles}>
                      {s.smiles}
                    </td>
                    <td className="p-2 text-right font-mono">{formatVal(s.qed, 3)}</td>
                    <td className="p-2 text-right font-mono">{formatVal(s.mw, 1)}</td>
                    <td className="p-2 text-right font-mono">{formatVal(s.logp, 2)}</td>
                    <td className="p-2 text-right font-mono">{formatVal(s.hbd, 0)}</td>
                    <td className="p-2 text-right font-mono">{formatVal(s.hba, 0)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          )}
        </div>

        {/* Footer */}
        {suggestions?.length > 0 && (
          <div className="px-6 py-3 border-t border-gray-100 flex items-center justify-between">
            <span className="text-xs text-gray-400">
              {selected.size} of {suggestions.length} selected
            </span>
            <div className="flex items-center gap-3">
              <button
                onClick={onClose}
                className="px-4 py-2 rounded-lg text-sm font-medium border border-gray-200 text-gray-600 hover:bg-gray-50 transition"
              >
                Cancel
              </button>
              <button
                onClick={handleImport}
                disabled={!selected.size || importing}
                className="px-4 py-2 rounded-lg text-sm font-medium bg-violet-600 text-white hover:bg-violet-700 transition disabled:opacity-50 disabled:cursor-not-allowed flex items-center gap-2"
              >
                {importing ? (
                  <>
                    <div className="w-4 h-4 border-2 border-white border-t-transparent rounded-full animate-spin" />
                    Importing...
                  </>
                ) : (
                  <>
                    <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 4v16m8-8H4" />
                    </svg>
                    Import {selected.size} molecule{selected.size !== 1 ? 's' : ''}
                  </>
                )}
              </button>
            </div>
          </div>
        )}
      </div>
    </div>
  )
}
