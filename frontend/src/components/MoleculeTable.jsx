import React, { useState, useMemo, useRef, useCallback } from 'react'

// ---------------------------------------------------------------------------
// Number formatter per column key
// ---------------------------------------------------------------------------
function formatNumber(key, value) {
  if (value == null || typeof value !== 'number') return null
  // Key-specific precision
  const decimals = {
    docking_score: 1,
    cnn_score: 2,
    cnn_affinity: 1,
    cnn_vs: 2,
    composite_score: 1,
    logP: 1,
    MW: 0,
    HBD: 0,
    HBA: 0,
    TPSA: 1,
    QED: 2,
    solubility: 2,
    BBB: 2,
    hERG: 2,
    metabolic_stability: 2,
    cluster_id: 0,
    generation_level: 0,
    interactions_count: 0,
  }
  const d = decimals[key] !== undefined ? decimals[key] : 2
  return d === 0 ? Math.round(value).toString() : value.toFixed(d)
}

// ---------------------------------------------------------------------------
// Color scale: returns Tailwind classes for a value in a column
// ---------------------------------------------------------------------------
function getValueColorClasses(value, allValues, scale) {
  if (value == null || !allValues || !allValues.length) return { cell: '', text: '' }
  const valid = allValues.filter(v => v != null)
  if (valid.length < 4) return { cell: '', text: '' }
  const sorted = [...valid].sort((a, b) => a - b)
  const n = sorted.length
  // Rank from 0 (lowest) to 1 (highest)
  const lteCount = sorted.filter(v => v <= value).length
  const rank = lteCount / n

  if (scale === 'higher-better') {
    if (rank >= 0.75) return { cell: 'bg-green-50', text: 'text-green-700 font-semibold' }
    if (rank <= 0.25) return { cell: 'bg-red-50', text: 'text-red-600' }
  } else if (scale === 'lower-better') {
    if (rank <= 0.25) return { cell: 'bg-green-50', text: 'text-green-700 font-semibold' }
    if (rank >= 0.75) return { cell: 'bg-red-50', text: 'text-red-600' }
  }
  return { cell: '', text: '' }
}

// ---------------------------------------------------------------------------
// Sort icon
// ---------------------------------------------------------------------------
function SortIcon({ direction, rank }) {
  if (!direction) {
    return (
      <svg className="w-3 h-3 text-gray-300 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
          d="M7 16V4m0 0L3 8m4-4l4 4m6 0v12m0 0l4-4m-4 4l-4-4" />
      </svg>
    )
  }
  return (
    <div className="flex items-center gap-0.5 flex-shrink-0">
      {rank && <span className="text-[8px] text-bx-light-text font-bold">{rank}</span>}
      <svg className="w-3 h-3 text-bx-light-text" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        {direction === 'asc' ? (
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 15l7-7 7 7" />
        ) : (
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
        )}
      </svg>
    </div>
  )
}

// ---------------------------------------------------------------------------
// Cell value renderer
// ---------------------------------------------------------------------------
function CellValue({ col, value, colorClasses }) {
  if (value == null || value === undefined) {
    return <span className="text-gray-300 select-none">—</span>
  }

  if (col.type === 'boolean') {
    return value ? (
      <span className="inline-flex items-center justify-center w-5 h-5 rounded-full bg-green-100">
        <svg className="w-3 h-3 text-green-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M5 13l4 4L19 7" />
        </svg>
      </span>
    ) : (
      <span className="inline-flex items-center justify-center w-5 h-5 rounded-full bg-red-50">
        <svg className="w-3 h-3 text-red-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M6 18L18 6M6 6l12 12" />
        </svg>
      </span>
    )
  }

  if (col.type === 'smiles') {
    const truncated = value.length > 20 ? value.slice(0, 20) + '…' : value
    return (
      <span className="font-mono text-[11px] text-gray-500" title={value}>
        {truncated}
      </span>
    )
  }

  if (col.type === 'number' && typeof value === 'number') {
    const formatted = formatNumber(col.key, value)
    return (
      <span className={`tabular-nums ${colorClasses?.text || 'text-gray-700'}`}>
        {formatted}
      </span>
    )
  }

  // Text
  const str = String(value)
  return (
    <span className="text-gray-700 block truncate" title={str} style={{ maxWidth: col.width ? `${col.width - 16}px` : 120 }}>
      {str}
    </span>
  )
}

// ---------------------------------------------------------------------------
// MoleculeTable — feature-complete data table
// Props:
//   molecules        — array of molecule objects
//   columns          — column definitions (from ALL_COLUMNS, filtered to visible)
//   selectedIds      — Set of selected molecule ids
//   onToggleSelect   — (molId) => void
//   onSelectAll      — () => void
//   onRowClick       — (molecule) => void
//   onToggleBookmark — (molId) => void — undefined means disabled
//   activeRowId      — id of molecule in detail panel (highlighted differently)
// ---------------------------------------------------------------------------
export default function MoleculeTable({
  molecules = [],
  columns = [],
  selectedIds = new Set(),
  onToggleSelect,
  onSelectAll,
  onRowClick,
  onToggleBookmark,
  activeRowId = null,
}) {
  const [sorts, setSorts] = useState([]) // [{key, dir}] max 2
  const tableRef = useRef(null)

  // Pre-compute column value arrays for color scaling
  const colValues = useMemo(() => {
    const map = {}
    columns.forEach(col => {
      if (col.type === 'number' && col.colorScale) {
        map[col.key] = molecules
          .map(m => m[col.key])
          .filter(v => v != null && typeof v === 'number')
      }
    })
    return map
  }, [molecules, columns])

  // Sorted molecules
  const sorted = useMemo(() => {
    if (!sorts.length) return molecules
    return [...molecules].sort((a, b) => {
      for (const { key, dir } of sorts) {
        const av = a[key], bv = b[key]
        if (av == null && bv == null) continue
        if (av == null) return 1
        if (bv == null) return -1
        const cmp = av < bv ? -1 : av > bv ? 1 : 0
        if (cmp !== 0) return dir === 'asc' ? cmp : -cmp
      }
      return 0
    })
  }, [molecules, sorts])

  // Header click: shift for multi-sort (max 2 levels)
  const handleHeaderClick = useCallback((col, e) => {
    if (!col.sortable) return
    setSorts(prev => {
      const existing = prev.find(s => s.key === col.key)
      if (e.shiftKey) {
        if (existing) {
          if (existing.dir === 'asc') {
            return prev.map(s => s.key === col.key ? { ...s, dir: 'desc' } : s)
          }
          return prev.filter(s => s.key !== col.key)
        }
        const next = [...prev, { key: col.key, dir: 'asc' }]
        return next.slice(-2)
      } else {
        if (existing) {
          if (existing.dir === 'asc') return [{ key: col.key, dir: 'desc' }]
          return []
        }
        return [{ key: col.key, dir: 'asc' }]
      }
    })
  }, [])

  function getSortDir(key) {
    const s = sorts.find(s => s.key === key)
    return s ? s.dir : null
  }

  function getSortRank(key) {
    if (sorts.length < 2) return null
    const idx = sorts.findIndex(s => s.key === key)
    return idx >= 0 ? (idx + 1).toString() : null
  }

  const allVisible = molecules.length > 0 && molecules.every(m => selectedIds.has(m.id))
  const someVisible = !allVisible && molecules.some(m => selectedIds.has(m.id))

  if (!molecules.length) {
    return (
      <div className="card p-10 text-center">
        <div className="flex flex-col items-center gap-3">
          <div className="w-12 h-12 rounded-full bg-gray-100 flex items-center justify-center">
            <svg className="w-6 h-6 text-gray-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
            </svg>
          </div>
          <p className="text-sm font-medium text-gray-600">No molecules in this view</p>
          <p className="text-sm text-gray-400">Try adjusting your filters or run an Import to add molecules</p>
        </div>
      </div>
    )
  }

  return (
    <div className="flex flex-col card overflow-hidden">
      {/* Scrollable table */}
      <div
        ref={tableRef}
        className="overflow-auto"
        style={{
          maxHeight: 'calc(100vh - 420px)',
          minHeight: 200,
          scrollbarWidth: 'thin',
          scrollbarColor: '#e5e7eb transparent',
        }}
      >
        <table
          className="w-full text-sm border-collapse"
          style={{ minWidth: Math.max(600, columns.reduce((sum, c) => sum + (c.width || 80), 0) + 90) }}
        >
          <thead className="sticky top-0 z-10">
            <tr className="bg-gray-50 border-b border-gray-200">
              {/* Checkbox column */}
              <th className="w-9 px-2 py-2.5 text-left border-b border-gray-200" style={{ width: 40 }}>
                <input
                  type="checkbox"
                  checked={allVisible}
                  ref={el => { if (el) el.indeterminate = someVisible }}
                  onChange={onSelectAll}
                  className="accent-bx-mint cursor-pointer w-3.5 h-3.5"
                  aria-label="Select all"
                />
              </th>

              {/* Bookmark column */}
              <th className="px-1 py-2.5 border-b border-gray-200" style={{ width: 32 }}>
                <svg className="w-3.5 h-3.5 text-gray-300 mx-auto" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                    d="M5 5a2 2 0 012-2h10a2 2 0 012 2v16l-7-3.5L5 21V5z" />
                </svg>
              </th>

              {/* Data columns */}
              {columns.map(col => {
                const dir = getSortDir(col.key)
                const rank = getSortRank(col.key)
                const isActive = dir !== null
                return (
                  <th
                    key={col.key}
                    className={`px-3 py-2.5 border-b border-gray-200 font-semibold uppercase tracking-wide text-[10px] whitespace-nowrap select-none ${
                      col.type === 'number' ? 'text-right' : 'text-left'
                    } ${col.sortable ? 'cursor-pointer transition-colors' : ''} ${
                      isActive ? 'text-bx-light-text bg-blue-50/50' : 'text-gray-500 hover:text-bx-light-text hover:bg-gray-100'
                    }`}
                    style={{ minWidth: col.width || 80 }}
                    onClick={e => handleHeaderClick(col, e)}
                    title={col.sortable ? 'Click to sort · Shift+click to add sort level' : col.label}
                  >
                    <span className={`inline-flex items-center gap-1 ${col.type === 'number' ? 'flex-row-reverse' : ''}`}>
                      {col.label}
                      {col.unit && (
                        <span className="font-normal text-gray-300 text-[9px]">({col.unit})</span>
                      )}
                      {col.sortable && <SortIcon direction={dir} rank={rank} />}
                    </span>
                  </th>
                )
              })}
            </tr>
          </thead>

          <tbody className="divide-y divide-gray-50">
            {sorted.map((mol, idx) => {
              const isSelected = selectedIds.has(mol.id)
              const isActive = activeRowId === mol.id
              const isBookmarked = mol.bookmarked

              let rowBg = idx % 2 === 0 ? 'bg-white' : 'bg-slate-50/40'
              if (isSelected) rowBg = 'bg-blue-50'
              if (isActive) rowBg = 'bg-blue-100/80'

              return (
                <tr
                  key={mol.id}
                  className={`cursor-pointer transition-colors duration-100 group ${rowBg} hover:bg-blue-50/60 ${
                    isSelected ? 'border-l-[3px] border-bx-surface' :
                    isBookmarked ? 'border-l-2 border-yellow-300' :
                    'border-l-2 border-transparent'
                  }`}
                  onClick={() => onRowClick && onRowClick(mol)}
                >
                  {/* Checkbox */}
                  <td
                    className="px-2 py-2 w-9"
                    onClick={e => e.stopPropagation()}
                  >
                    <input
                      type="checkbox"
                      checked={isSelected}
                      onChange={() => onToggleSelect && onToggleSelect(mol.id)}
                      className="accent-bx-mint cursor-pointer w-3.5 h-3.5"
                      aria-label={`Select ${mol.name || mol.id}`}
                    />
                  </td>

                  {/* Bookmark star */}
                  <td
                    className="px-1 py-2 text-center w-8"
                    onClick={e => e.stopPropagation()}
                  >
                    <button
                      onClick={() => onToggleBookmark && onToggleBookmark(mol.id)}
                      disabled={!onToggleBookmark}
                      className={`p-0.5 rounded transition-colors ${
                        onToggleBookmark ? 'hover:bg-yellow-50 cursor-pointer' : 'cursor-default opacity-50'
                      }`}
                      aria-label={isBookmarked ? 'Remove bookmark' : 'Bookmark'}
                    >
                      <svg
                        className={`w-3.5 h-3.5 transition-all ${
                          isBookmarked ? 'text-yellow-400 fill-yellow-400 drop-shadow-sm' :
                          'text-gray-200 fill-none group-hover:text-gray-300'
                        }`}
                        stroke="currentColor" viewBox="0 0 24 24"
                      >
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                          d="M5 5a2 2 0 012-2h10a2 2 0 012 2v16l-7-3.5L5 21V5z" />
                      </svg>
                    </button>
                  </td>

                  {/* Data cells */}
                  {columns.map(col => {
                    const value = mol[col.key]
                    const colorClasses = col.colorScale && col.type === 'number'
                      ? getValueColorClasses(value, colValues[col.key] || [], col.colorScale)
                      : { cell: '', text: '' }

                    return (
                      <td
                        key={col.key}
                        className={`px-3 py-2 ${col.type === 'number' ? 'text-right tabular-nums' : 'text-left'} ${colorClasses.cell}`}
                      >
                        <CellValue col={col} value={value} colorClasses={colorClasses} />
                      </td>
                    )
                  })}
                </tr>
              )
            })}
          </tbody>
        </table>
      </div>

      {/* Footer */}
      <div className="px-4 py-2 bg-gray-50 border-t border-gray-100 flex items-center justify-between">
        <span className="text-sm text-gray-400">
          Showing{' '}
          <strong className="text-gray-600 tabular-nums">{sorted.length}</strong>
          {' '}of{' '}
          <strong className="text-gray-600 tabular-nums">{molecules.length}</strong>
          {' '}molecules
        </span>
        <div className="flex items-center gap-3">
          {sorts.length > 0 && (
            <div className="flex items-center gap-1.5 text-sm text-gray-500">
              <span>Sorted by</span>
              {sorts.map((s, i) => (
                <span key={s.key} className="flex items-center gap-0.5">
                  {i > 0 && <span className="text-gray-300">·</span>}
                  <span className="font-medium text-bx-light-text">
                    {columns.find(c => c.key === s.key)?.label || s.key}
                  </span>
                  <span className="text-[10px] text-gray-400">{s.dir === 'asc' ? '↑' : '↓'}</span>
                </span>
              ))}
              <button
                onClick={() => setSorts([])}
                className="ml-1 text-gray-400 hover:text-red-500 transition-colors text-[10px] underline underline-offset-2"
              >
                Clear
              </button>
            </div>
          )}
          {selectedIds.size > 0 && (
            <span className="text-sm text-bx-light-text font-medium tabular-nums">
              {selectedIds.size} selected
            </span>
          )}
        </div>
      </div>
    </div>
  )
}
