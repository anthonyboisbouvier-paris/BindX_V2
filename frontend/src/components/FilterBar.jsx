import React, { useState, useMemo, useEffect, useRef, useCallback } from 'react'

// ---------------------------------------------------------------------------
// FilterBar — compact inline filtering: search + quick toggles + column filters
// Props:
//   molecules      — array of molecule objects
//   columns        — array of column definitions (from ALL_COLUMNS)
//   onFilteredChange — callback(filteredMolecules)
// ---------------------------------------------------------------------------

const QUICK_FILTERS = [
  {
    id: 'bookmarked', label: 'Bookmarked',
    icon: <svg className="w-3 h-3" fill="currentColor" viewBox="0 0 24 24"><path d="M5 5a2 2 0 012-2h10a2 2 0 012 2v16l-7-3.5L5 21V5z" /></svg>,
    test: m => m.bookmarked === true,
  },
  {
    id: 'has_docking', label: 'Docking',
    icon: <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24"><circle cx="12" cy="12" r="3" strokeWidth={1.5}/><path strokeLinecap="round" strokeWidth={1.5} d="M12 2v4m0 12v4M2 12h4m12 0h4"/></svg>,
    test: m => m.docking_score != null,
  },
  {
    id: 'has_admet', label: 'ADMET',
    icon: <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M9 12l2 2 4-4m5.618-4.016A11.955 11.955 0 0112 2.944a11.955 11.955 0 01-8.618 3.04A12.02 12.02 0 003 9c0 5.591 3.824 10.29 9 11.622 5.176-1.332 9-6.03 9-11.622 0-1.042-.133-2.052-.382-3.016z"/></svg>,
    test: m => m.admet != null || m.QED != null,
  },
  {
    id: 'lipinski', label: 'Lipinski',
    icon: <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M5 13l4 4L19 7"/></svg>,
    test: m => m.lipinski_pass === true,
  },
]

function XIcon({ className = 'w-3 h-3' }) {
  return (
    <svg className={className} fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M6 18L18 6M6 6l12 12" />
    </svg>
  )
}

// ---------------------------------------------------------------------------
// Inline filter popover — appears when user picks a column from the select
// ---------------------------------------------------------------------------
function InlineFilterConfig({ col, molecules, onApply, onCancel }) {
  const [numMin, setNumMin] = useState('')
  const [numMax, setNumMax] = useState('')
  const [textVal, setTextVal] = useState('')
  const [boolVal, setBoolVal] = useState(true)
  const ref = useRef(null)

  const colRange = useMemo(() => {
    if (col.type !== 'number') return null
    const vals = molecules.map(m => m[col.key]).filter(v => v != null && typeof v === 'number')
    if (!vals.length) return null
    return { min: Math.min(...vals), max: Math.max(...vals) }
  }, [col, molecules])

  useEffect(() => {
    function handler(e) { if (ref.current && !ref.current.contains(e.target)) onCancel() }
    document.addEventListener('mousedown', handler)
    return () => document.removeEventListener('mousedown', handler)
  }, [onCancel])

  function handleApply() {
    let label, test
    if (col.type === 'number') {
      const mn = numMin !== '' ? parseFloat(numMin) : null
      const mx = numMax !== '' ? parseFloat(numMax) : null
      if (mn == null && mx == null) return
      const parts = []
      if (mn != null) parts.push(`≥ ${mn}`)
      if (mx != null) parts.push(`≤ ${mx}`)
      label = `${col.label} ${parts.join(', ')}`
      test = m => {
        const v = m[col.key]
        if (v == null) return false
        if (mn != null && v < mn) return false
        if (mx != null && v > mx) return false
        return true
      }
    } else if (col.type === 'boolean') {
      label = `${col.label} = ${boolVal ? 'Yes' : 'No'}`
      test = m => m[col.key] === boolVal
    } else {
      if (!textVal.trim()) return
      const tv = textVal.trim().toLowerCase()
      label = `${col.label} ~ "${textVal.trim()}"`
      test = m => {
        const v = m[col.key]
        return v != null && String(v).toLowerCase().includes(tv)
      }
    }
    onApply({ id: `col_${col.key}_${Date.now()}`, label, test, colKey: col.key })
  }

  function handleKeyDown(e) {
    if (e.key === 'Enter') handleApply()
    if (e.key === 'Escape') onCancel()
  }

  return (
    <div ref={ref} className="absolute left-0 top-full mt-1 z-50 bg-white border border-gray-200 rounded-xl shadow-xl w-64 overflow-hidden">
      <div className="px-3 py-2 border-b border-gray-100 flex items-center justify-between">
        <span className="text-sm font-semibold text-gray-700">{col.label}</span>
        {colRange && (
          <span className="text-[10px] text-gray-400 tabular-nums">
            {colRange.min.toFixed(1)} — {colRange.max.toFixed(1)}
          </span>
        )}
      </div>
      <div className="px-3 py-2.5">
        {col.type === 'number' && (
          <div className="flex gap-2">
            <input
              type="number" value={numMin} onChange={e => setNumMin(e.target.value)}
              onKeyDown={handleKeyDown}
              placeholder={colRange ? `Min (${colRange.min.toFixed(1)})` : 'Min'}
              className="flex-1 text-sm border border-gray-200 rounded-lg px-2 py-1.5 focus:outline-none focus:ring-1 focus:ring-bx-mint w-0"
              step="0.1" autoFocus
            />
            <input
              type="number" value={numMax} onChange={e => setNumMax(e.target.value)}
              onKeyDown={handleKeyDown}
              placeholder={colRange ? `Max (${colRange.max.toFixed(1)})` : 'Max'}
              className="flex-1 text-sm border border-gray-200 rounded-lg px-2 py-1.5 focus:outline-none focus:ring-1 focus:ring-bx-mint w-0"
              step="0.1"
            />
          </div>
        )}
        {col.type === 'boolean' && (
          <div className="flex gap-2">
            {[true, false].map(v => (
              <button key={String(v)} onClick={() => setBoolVal(v)}
                className={`flex-1 py-1.5 rounded-lg text-sm font-medium border transition-colors ${
                  boolVal === v
                    ? 'bg-bx-surface text-white border-bx-surface'
                    : 'border-gray-200 text-gray-600 hover:bg-gray-50'
                }`}
              >{v ? 'Yes' : 'No'}</button>
            ))}
          </div>
        )}
        {(col.type === 'text' || col.type === 'smiles') && (
          <input
            type="text" value={textVal} onChange={e => setTextVal(e.target.value)}
            onKeyDown={handleKeyDown}
            placeholder={`Contains...`}
            className="w-full text-sm border border-gray-200 rounded-lg px-2 py-1.5 focus:outline-none focus:ring-1 focus:ring-bx-mint"
            autoFocus
          />
        )}
      </div>
      <div className="px-3 py-2 border-t border-gray-100 flex justify-end gap-2">
        <button onClick={onCancel} className="text-xs text-gray-400 hover:text-gray-600 px-2 py-1">Cancel</button>
        <button onClick={handleApply}
          className="px-3 py-1 bg-bx-surface text-white text-xs rounded-lg hover:bg-bx-elevated transition-colors font-medium"
        >Apply</button>
      </div>
    </div>
  )
}

// ---------------------------------------------------------------------------
// Main FilterBar
// ---------------------------------------------------------------------------
export default function FilterBar({ molecules = [], columns = [], onFilteredChange }) {
  const [activeQuickFilters, setActiveQuickFilters] = useState(new Set())
  const [customFilters, setCustomFilters] = useState([])
  const [searchText, setSearchText] = useState('')
  const [configCol, setConfigCol] = useState(null) // column being configured
  const addBtnRef = useRef(null)

  // Filterable columns for the dropdown
  const filterableCols = useMemo(
    () => columns.filter(c => c.sortable && c.type !== 'action' && c.type !== 'invalidation' && c.type !== 'tags' && c.type !== 'editable_text'),
    [columns]
  )

  // Compute filtered molecules
  const filteredMolecules = useMemo(() => {
    let result = molecules

    // Global text search
    if (searchText.trim()) {
      const q = searchText.trim().toLowerCase()
      result = result.filter(m =>
        (m.name && m.name.toLowerCase().includes(q)) ||
        (m.smiles && m.smiles.toLowerCase().includes(q)) ||
        (m.id && String(m.id).toLowerCase().includes(q))
      )
    }

    // Quick filters (AND)
    for (const qfId of activeQuickFilters) {
      const qf = QUICK_FILTERS.find(f => f.id === qfId)
      if (qf) result = result.filter(qf.test)
    }

    // Custom column filters (AND)
    for (const cf of customFilters) {
      result = result.filter(cf.test)
    }

    return result
  }, [molecules, activeQuickFilters, customFilters, searchText])

  // Notify parent
  useEffect(() => {
    onFilteredChange && onFilteredChange(filteredMolecules)
  }, [filteredMolecules, onFilteredChange])

  const hasAnyFilter = activeQuickFilters.size > 0 || customFilters.length > 0 || searchText.trim() !== ''
  const filterCount = activeQuickFilters.size + customFilters.length + (searchText.trim() ? 1 : 0)

  const toggleQuickFilter = useCallback((id) => {
    setActiveQuickFilters(prev => {
      const next = new Set(prev)
      if (next.has(id)) next.delete(id)
      else next.add(id)
      return next
    })
  }, [])

  const removeCustomFilter = useCallback((id) => {
    setCustomFilters(prev => prev.filter(f => f.id !== id))
  }, [])

  const clearAll = useCallback(() => {
    setActiveQuickFilters(new Set())
    setCustomFilters([])
    setSearchText('')
  }, [])

  const handleColSelect = useCallback((e) => {
    const key = e.target.value
    if (!key) return
    const col = filterableCols.find(c => c.key === key)
    if (col) setConfigCol(col)
    e.target.value = ''
  }, [filterableCols])

  return (
    <div className="card px-4 py-2.5">
      {/* Single compact row */}
      <div className="flex items-center gap-2 flex-wrap">
        {/* Search input */}
        <div className="relative">
          <svg className="absolute left-2 top-1/2 -translate-y-1/2 w-3.5 h-3.5 text-gray-400 pointer-events-none" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
          </svg>
          <input
            value={searchText}
            onChange={e => setSearchText(e.target.value)}
            placeholder="Search..."
            className="text-sm border border-gray-200 rounded-lg pl-7 pr-2 py-1.5 w-40 focus:outline-none focus:ring-1 focus:ring-bx-mint focus:w-52 transition-all"
          />
          {searchText && (
            <button onClick={() => setSearchText('')} className="absolute right-1.5 top-1/2 -translate-y-1/2 text-gray-400 hover:text-gray-600">
              <XIcon className="w-2.5 h-2.5" />
            </button>
          )}
        </div>

        <span className="text-gray-200 text-sm select-none">|</span>

        {/* Quick filter toggle pills */}
        {QUICK_FILTERS.map(qf => {
          const isActive = activeQuickFilters.has(qf.id)
          const count = molecules.filter(qf.test).length
          if (count === 0) return null
          return (
            <button
              key={qf.id}
              onClick={() => toggleQuickFilter(qf.id)}
              className={`flex items-center gap-1 px-2 py-1 rounded-full text-xs font-medium border transition-all duration-150 ${
                isActive
                  ? 'bg-bx-surface text-white border-bx-surface shadow-sm'
                  : 'bg-white text-gray-500 border-gray-200 hover:border-gray-300 hover:text-gray-700'
              }`}
              title={`${qf.label} (${count})`}
            >
              {qf.icon}
              {qf.label}
              <span className={`text-[9px] tabular-nums ${isActive ? 'text-white/60' : 'text-gray-400'}`}>{count}</span>
            </button>
          )
        })}

        <span className="text-gray-200 text-sm select-none">|</span>

        {/* Column filter selector */}
        <div className="relative" ref={addBtnRef}>
          <select
            onChange={handleColSelect}
            defaultValue=""
            className="text-xs border border-dashed border-gray-300 rounded-lg pl-2 pr-6 py-1.5 text-gray-500 bg-white hover:border-gray-400 focus:outline-none focus:ring-1 focus:ring-bx-mint cursor-pointer appearance-none"
            style={{ backgroundImage: `url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 20 20' fill='%239ca3af'%3E%3Cpath fill-rule='evenodd' d='M5.23 7.21a.75.75 0 011.06.02L10 11.168l3.71-3.938a.75.75 0 111.08 1.04l-4.25 4.5a.75.75 0 01-1.08 0l-4.25-4.5a.75.75 0 01.02-1.06z'/%3E%3C/svg%3E")`, backgroundRepeat: 'no-repeat', backgroundPosition: 'right 4px center', backgroundSize: '14px' }}
          >
            <option value="">+ Column filter</option>
            {filterableCols.map(c => (
              <option key={c.key} value={c.key}>{c.label} ({c.type})</option>
            ))}
          </select>

          {configCol && (
            <InlineFilterConfig
              col={configCol}
              molecules={molecules}
              onApply={(filter) => {
                setCustomFilters(prev => [...prev, filter])
                setConfigCol(null)
              }}
              onCancel={() => setConfigCol(null)}
            />
          )}
        </div>

        {/* Spacer */}
        <div className="flex-1" />

        {/* Status + clear */}
        <span className={`text-xs tabular-nums ${hasAnyFilter ? 'text-bx-light-text font-medium' : 'text-gray-400'}`}>
          {filteredMolecules.length}/{molecules.length}
        </span>
        {hasAnyFilter && (
          <button onClick={clearAll} className="text-xs text-gray-400 hover:text-red-500 transition-colors">
            Clear
          </button>
        )}
      </div>

      {/* Active custom filter chips (only shown when present) */}
      {customFilters.length > 0 && (
        <div className="flex flex-wrap items-center gap-1.5 mt-2 pt-2 border-t border-gray-100">
          {customFilters.map(cf => (
            <span key={cf.id}
              className="inline-flex items-center gap-1 px-2 py-0.5 bg-blue-50 border border-blue-200 text-blue-700 rounded-full text-[11px] font-medium"
            >
              {cf.label}
              <button onClick={() => removeCustomFilter(cf.id)}
                className="text-blue-400 hover:text-blue-700 transition-colors"
              ><XIcon className="w-2.5 h-2.5" /></button>
            </span>
          ))}
        </div>
      )}
    </div>
  )
}
