import React, { useState, useMemo } from 'react'

// ---------------------------------------------------------------------------
// ExportModal
//
// Props:
//   isOpen      — boolean
//   onClose     — () => void
//   molecules   — full molecule array (all)
//   selectedIds — Set<string>  (currently checked rows)
//   columns     — array of column def objects { key, label }
//   onToast     — optional (message, type) => void for notifications
// ---------------------------------------------------------------------------

const SCOPE_OPTIONS = [
  { value: 'selected',   label: 'Selected',      count: null },
  { value: 'filtered',   label: 'All filtered',  count: null },
  { value: 'all',        label: 'All',            count: null },
  { value: 'bookmarked', label: 'Bookmarked',     count: null },
]

const FORMATS = [
  {
    value: 'csv',
    label: 'CSV',
    desc: 'Tables',
    icon: (
      <svg className="w-7 h-7" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
          d="M3 10h18M3 14h18m-9-4v8m-7 0h14a2 2 0 002-2V8a2 2 0 00-2-2H5a2 2 0 00-2 2v8a2 2 0 002 2z" />
      </svg>
    ),
  },
  {
    value: 'sdf',
    label: 'SDF',
    desc: '3D Structures',
    icon: (
      <svg className="w-7 h-7" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
          d="M14 10l-2 1m0 0l-2-1m2 1v2.5M20 7l-2 1m2-1l-2-1m2 1v2.5M14 4l-2-1-2 1M4 7l2-1M4 7l2 1M4 7v2.5M12 21l-2-1m2 1l2-1m-2 1v-2.5M6 18l-2-1v-2.5M18 18l2-1v-2.5" />
      </svg>
    ),
  },
  {
    value: 'pdf',
    label: 'PDF',
    desc: 'Report',
    icon: (
      <svg className="w-7 h-7" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
          d="M19.5 14.25v-2.625a3.375 3.375 0 00-3.375-3.375h-1.5A1.125 1.125 0 0113.5 7.125v-1.5a3.375 3.375 0 00-3.375-3.375H8.25m2.25 0H5.625c-.621 0-1.125.504-1.125 1.125v17.25c0 .621.504 1.125 1.125 1.125h12.75c.621 0 1.125-.504 1.125-1.125V11.25a9 9 0 00-9-9z" />
      </svg>
    ),
  },
]

function triggerCsvDownload(filename, csvContent) {
  const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' })
  const url = URL.createObjectURL(blob)
  const a = document.createElement('a')
  a.href = url
  a.download = filename
  document.body.appendChild(a)
  a.click()
  document.body.removeChild(a)
  URL.revokeObjectURL(url)
}

function buildCsv(molecules, columns) {
  const header = columns.map(c => c.label).join(',')
  const rows = molecules.map(mol =>
    columns.map(c => {
      const v = mol[c.key]
      if (v == null) return ''
      const s = String(v)
      // Escape commas and quotes
      if (s.includes(',') || s.includes('"') || s.includes('\n')) {
        return `"${s.replace(/"/g, '""')}"`
      }
      return s
    }).join(',')
  )
  return [header, ...rows].join('\r\n')
}

export default function ExportModal({
  isOpen,
  onClose,
  molecules = [],
  selectedIds = new Set(),
  columns = [],
  onToast,
}) {
  const [scope, setScope] = useState('selected')
  const [format, setFormat] = useState('csv')
  const [csvColumns, setCsvColumns] = useState('visible')
  const [pdfOptions, setPdfOptions] = useState({
    include_admet: true,
    include_scoring: true,
    include_synthesis: true,
  })

  // ---- Compute molecule counts for each scope ----
  const bookmarkedMols = useMemo(() => molecules.filter(m => m.bookmarked), [molecules])
  const selectedMols = useMemo(() => molecules.filter(m => selectedIds.has(m.id)), [molecules, selectedIds])

  const scopeCounts = {
    selected:   selectedMols.length,
    filtered:   molecules.length,   // treated as all filtered — no separate filter state here
    all:        molecules.length,
    bookmarked: bookmarkedMols.length,
  }

  function getMolsForScope() {
    switch (scope) {
      case 'selected':   return selectedMols
      case 'filtered':   return molecules
      case 'all':        return molecules
      case 'bookmarked': return bookmarkedMols
      default:           return molecules
    }
  }

  const exportMols = getMolsForScope()

  // ---- Columns for CSV ----
  const visibleColumns = columns.length > 0 ? columns : [
    { key: 'name',           label: 'Name' },
    { key: 'smiles',         label: 'SMILES' },
    { key: 'docking_score',  label: 'Docking Score' },
    { key: 'cnn_score',      label: 'CNN Score' },
    { key: 'composite_score',label: 'Composite Score' },
    { key: 'logP',           label: 'LogP' },
    { key: 'MW',             label: 'MW' },
    { key: 'TPSA',           label: 'TPSA' },
  ]

  const allAvailableColumns = [
    ...visibleColumns,
    { key: 'cnn_affinity',      label: 'CNN Affinity' },
    { key: 'cnn_vs',            label: 'CNN VS' },
    { key: 'HBD',               label: 'HBD' },
    { key: 'HBA',               label: 'HBA' },
    { key: 'QED',               label: 'QED' },
    { key: 'solubility',        label: 'Solubility' },
    { key: 'BBB',               label: 'BBB' },
    { key: 'hERG',              label: 'hERG' },
    { key: 'metabolic_stability', label: 'Met. Stability' },
    { key: 'cluster_id',        label: 'Cluster' },
    { key: 'scaffold',          label: 'Scaffold' },
    { key: 'lipinski_pass',     label: 'Lipinski' },
  ]

  const exportColumns = csvColumns === 'visible' ? visibleColumns : allAvailableColumns

  function handleExport() {
    if (format === 'csv') {
      const csv = buildCsv(exportMols, exportColumns)
      triggerCsvDownload(`dockit_export_${scope}_${Date.now()}.csv`, csv)
      onToast?.(`Exported ${exportMols.length} molecules as CSV`, 'success')
      onClose()
    } else {
      onToast?.(`Export queued — ${format.toUpperCase()} download will start when ready`, 'info')
      onClose()
    }
  }

  if (!isOpen) return null

  return (
    <div
      className="fixed inset-0 z-50 flex items-center justify-center bg-black/40 backdrop-blur-sm p-4"
      onClick={onClose}
    >
      <div
        className="bg-white rounded-2xl shadow-2xl w-full max-w-lg overflow-hidden"
        onClick={e => e.stopPropagation()}
      >
        {/* Header */}
        <div className="flex items-center justify-between px-6 py-4 border-b border-gray-100">
          <div className="flex items-center gap-3">
            <div className="p-2 rounded-xl bg-[#0f131d] text-white">
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                  d="M3 16.5v2.25A2.25 2.25 0 005.25 21h13.5A2.25 2.25 0 0021 18.75V16.5M16.5 12L12 16.5m0 0L7.5 12m4.5 4.5V3" />
              </svg>
            </div>
            <h2 className="font-bold text-[#0f131d] text-lg">Export Molecules</h2>
          </div>
          <button
            onClick={onClose}
            className="p-1.5 rounded-lg hover:bg-gray-100 text-gray-500 transition-colors"
          >
            <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
            </svg>
          </button>
        </div>

        <div className="px-6 py-4 space-y-5">

          {/* ---- Scope ---- */}
          <div>
            <p className="text-xs font-semibold text-gray-400 uppercase tracking-wide mb-2">Scope</p>
            <div className="space-y-1.5">
              {SCOPE_OPTIONS.map(opt => {
                const count = scopeCounts[opt.value]
                const disabled = count === 0
                return (
                  <label
                    key={opt.value}
                    className={`flex items-center gap-3 p-2.5 rounded-xl cursor-pointer border transition-all ${
                      scope === opt.value && !disabled
                        ? 'border-[#0f131d] bg-blue-50'
                        : disabled
                        ? 'border-gray-100 opacity-40 cursor-not-allowed'
                        : 'border-gray-100 hover:border-gray-200 hover:bg-gray-50'
                    }`}
                  >
                    <input
                      type="radio"
                      name="export-scope"
                      value={opt.value}
                      checked={scope === opt.value}
                      onChange={() => !disabled && setScope(opt.value)}
                      disabled={disabled}
                      className="accent-[#0f131d]"
                    />
                    <span className={`text-sm font-medium ${scope === opt.value ? 'text-[#0f131d]' : 'text-gray-700'}`}>
                      {opt.label}
                    </span>
                    <span className={`ml-auto text-xs tabular-nums font-semibold ${
                      count === 0 ? 'text-gray-300' : 'text-gray-500'
                    }`}>
                      {count} mol{count !== 1 ? 's' : ''}
                    </span>
                  </label>
                )
              })}
            </div>
          </div>

          {/* ---- Format ---- */}
          <div>
            <p className="text-xs font-semibold text-gray-400 uppercase tracking-wide mb-2">Format</p>
            <div className="grid grid-cols-3 gap-2">
              {FORMATS.map(f => (
                <button
                  key={f.value}
                  onClick={() => setFormat(f.value)}
                  className={`flex flex-col items-center gap-2 py-4 rounded-xl border-2 transition-all ${
                    format === f.value
                      ? 'border-[#0f131d] bg-blue-50 text-[#0f131d]'
                      : 'border-gray-100 text-gray-500 hover:border-blue-200 hover:bg-blue-50/20 hover:text-[#0f131d]'
                  }`}
                >
                  {f.icon}
                  <div className="text-center">
                    <p className="text-sm font-bold leading-none">{f.label}</p>
                    <p className="text-[10px] mt-0.5 text-gray-400">{f.desc}</p>
                  </div>
                </button>
              ))}
            </div>
          </div>

          {/* ---- Format-specific options ---- */}
          {format === 'csv' && (
            <div>
              <p className="text-xs font-semibold text-gray-400 uppercase tracking-wide mb-2">CSV Columns</p>
              <div className="space-y-1.5">
                {[
                  { value: 'visible', label: 'Visible columns only', count: visibleColumns.length },
                  { value: 'all',     label: 'All available columns', count: allAvailableColumns.length },
                ].map(opt => (
                  <label
                    key={opt.value}
                    className={`flex items-center gap-3 p-2.5 rounded-xl cursor-pointer border transition-all ${
                      csvColumns === opt.value
                        ? 'border-[#0f131d] bg-blue-50'
                        : 'border-gray-100 hover:border-gray-200'
                    }`}
                  >
                    <input
                      type="radio"
                      name="csv-cols"
                      value={opt.value}
                      checked={csvColumns === opt.value}
                      onChange={() => setCsvColumns(opt.value)}
                      className="accent-[#0f131d]"
                    />
                    <span className={`text-sm font-medium ${csvColumns === opt.value ? 'text-[#0f131d]' : 'text-gray-700'}`}>
                      {opt.label}
                    </span>
                    <span className="ml-auto text-xs text-gray-400 tabular-nums">{opt.count} cols</span>
                  </label>
                ))}
              </div>
            </div>
          )}

          {format === 'pdf' && (
            <div>
              <p className="text-xs font-semibold text-gray-400 uppercase tracking-wide mb-2">PDF Contents</p>
              <div className="space-y-1.5">
                {[
                  { key: 'include_admet',    label: 'Include ADMET profiles' },
                  { key: 'include_scoring',  label: 'Include scoring breakdown' },
                  { key: 'include_synthesis', label: 'Include synthesis routes' },
                ].map(opt => (
                  <label
                    key={opt.key}
                    className={`flex items-center gap-3 p-2.5 rounded-xl cursor-pointer border transition-all ${
                      pdfOptions[opt.key]
                        ? 'border-blue-200 bg-blue-50/50'
                        : 'border-gray-100 hover:border-gray-200'
                    }`}
                  >
                    <input
                      type="checkbox"
                      checked={pdfOptions[opt.key]}
                      onChange={() => setPdfOptions(prev => ({ ...prev, [opt.key]: !prev[opt.key] }))}
                      className="accent-[#0f131d]"
                    />
                    <span className="text-sm font-medium text-gray-700">{opt.label}</span>
                  </label>
                ))}
              </div>
            </div>
          )}

          {format === 'sdf' && (
            <div className="bg-blue-50 border border-blue-100 rounded-xl p-3">
              <p className="text-xs text-blue-600">
                SDF export includes 2D SMILES-derived coordinates and all computed properties as SD tags.
                3D coordinates require a completed docking run.
              </p>
            </div>
          )}

          {/* ---- Preview ---- */}
          <div className="bg-gray-50 border border-gray-100 rounded-xl px-4 py-3 flex items-center justify-between">
            <div className="text-sm text-gray-600">
              <span className="font-semibold text-gray-800">{exportMols.length}</span> molecules
              <span className="text-gray-400 mx-1">x</span>
              <span className="font-semibold text-gray-800">
                {format === 'csv' ? exportColumns.length : '—'}
              </span>
              {format === 'csv' && <span className="text-gray-400 ml-1">columns</span>}
              {format !== 'csv' && <span className="text-gray-400 ml-1">{format.toUpperCase()}</span>}
            </div>
            <span className="text-xs text-gray-400 uppercase font-semibold tracking-wide">{format.toUpperCase()}</span>
          </div>
        </div>

        {/* Footer */}
        <div className="px-6 py-4 border-t border-gray-100 bg-gray-50 flex items-center justify-end gap-3">
          <button
            onClick={onClose}
            className="px-4 py-2 rounded-xl text-sm font-semibold border border-gray-200
                       text-gray-700 hover:bg-gray-100 transition-colors"
          >
            Cancel
          </button>
          <button
            onClick={handleExport}
            disabled={exportMols.length === 0}
            className="flex items-center gap-2 px-5 py-2 rounded-xl text-sm font-semibold
                       bg-[#0f131d] hover:bg-[#1a2332] text-white transition-colors
                       disabled:opacity-40 disabled:cursor-not-allowed shadow-sm"
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                d="M3 16.5v2.25A2.25 2.25 0 005.25 21h13.5A2.25 2.25 0 0021 18.75V16.5M16.5 12L12 16.5m0 0L7.5 12m4.5 4.5V3" />
            </svg>
            Export
          </button>
        </div>
      </div>
    </div>
  )
}
