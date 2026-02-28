import React, { useState, useMemo, useEffect, useRef, useCallback } from 'react'

// ---------------------------------------------------------------------------
// FilterBar — smart filtering with quick toggles and custom filter chips
// Props:
//   molecules      — array of molecule objects
//   columns        — array of column definitions (from ALL_COLUMNS)
//   onFilteredChange — callback(filteredMolecules)
// ---------------------------------------------------------------------------

const QUICK_FILTERS = [
  { id: 'bookmarked',   label: 'Bookmarked',    test: m => m.bookmarked === true },
  { id: 'has_docking',  label: 'Has docking',   test: m => m.docking_score != null },
  { id: 'has_admet',    label: 'Has ADMET',      test: m => m.admet != null },
  { id: 'lipinski',     label: 'Lipinski pass',  test: m => m.lipinski_pass === true },
]

function ChevronDown({ className = 'w-3.5 h-3.5' }) {
  return (
    <svg className={className} fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
    </svg>
  )
}

function XIcon({ className = 'w-3 h-3' }) {
  return (
    <svg className={className} fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M6 18L18 6M6 6l12 12" />
    </svg>
  )
}

// Add-filter dropdown — lets the user pick a column then configure a filter
function AddFilterDropdown({ columns, molecules, onAdd, onClose }) {
  const [step, setStep] = useState('pick_column') // 'pick_column' | 'configure'
  const [selectedCol, setSelectedCol] = useState(null)
  const [numMin, setNumMin] = useState('')
  const [numMax, setNumMax] = useState('')
  const [textVal, setTextVal] = useState('')
  const [boolVal, setBoolVal] = useState(true)
  const ref = useRef(null)

  // Close on outside click
  useEffect(() => {
    function handler(e) {
      if (ref.current && !ref.current.contains(e.target)) onClose()
    }
    document.addEventListener('mousedown', handler)
    return () => document.removeEventListener('mousedown', handler)
  }, [onClose])

  // Filterable columns — exclude boolean columns handled by quick filters and non-filterable
  const filterableCols = columns.filter(c => c.sortable && c.type !== 'smiles')

  // Compute data range for selected numeric column
  const colRange = useMemo(() => {
    if (!selectedCol || selectedCol.type !== 'number') return null
    const vals = molecules.map(m => m[selectedCol.key]).filter(v => v != null && typeof v === 'number')
    if (!vals.length) return null
    return { min: Math.min(...vals), max: Math.max(...vals) }
  }, [selectedCol, molecules])

  function handleApply() {
    if (!selectedCol) return
    let label, test
    if (selectedCol.type === 'number') {
      const mn = numMin !== '' ? parseFloat(numMin) : null
      const mx = numMax !== '' ? parseFloat(numMax) : null
      if (mn == null && mx == null) return
      const parts = []
      if (mn != null) parts.push(`>= ${mn}`)
      if (mx != null) parts.push(`<= ${mx}`)
      label = `${selectedCol.label} ${parts.join(' and ')}`
      test = m => {
        const v = m[selectedCol.key]
        if (v == null) return false
        if (mn != null && v < mn) return false
        if (mx != null && v > mx) return false
        return true
      }
    } else if (selectedCol.type === 'boolean') {
      label = `${selectedCol.label} = ${boolVal ? 'Yes' : 'No'}`
      test = m => m[selectedCol.key] === boolVal
    } else {
      if (!textVal.trim()) return
      label = `${selectedCol.label} contains "${textVal.trim()}"`
      const tv = textVal.trim().toLowerCase()
      test = m => {
        const v = m[selectedCol.key]
        return v != null && String(v).toLowerCase().includes(tv)
      }
    }
    onAdd({ id: `custom_${Date.now()}`, label, test })
    onClose()
  }

  return (
    <div
      ref={ref}
      className="absolute left-0 top-full mt-1 z-40 bg-white border border-gray-200 rounded-xl shadow-xl w-72 overflow-hidden"
    >
      {step === 'pick_column' ? (
        <>
          <div className="px-3 pt-3 pb-2 border-b border-gray-100">
            <p className="text-xs font-semibold text-gray-600 uppercase tracking-wide">Add Filter — Select Column</p>
          </div>
          <div className="overflow-y-auto max-h-64 py-1">
            {filterableCols.map(col => (
              <button
                key={col.key}
                onClick={() => { setSelectedCol(col); setStep('configure') }}
                className="w-full flex items-center justify-between px-3 py-2 text-sm text-gray-700 hover:bg-blue-50 hover:text-[#0f131d] transition-colors text-left"
              >
                <span>{col.label}</span>
                <span className="text-[10px] text-gray-400 uppercase">{col.type}</span>
              </button>
            ))}
          </div>
          <div className="px-3 py-2 border-t border-gray-100">
            <button onClick={onClose} className="text-xs text-gray-400 hover:text-gray-600">Cancel</button>
          </div>
        </>
      ) : (
        <>
          <div className="px-3 pt-3 pb-2 border-b border-gray-100 flex items-center gap-2">
            <button onClick={() => setStep('pick_column')} className="text-[#0f131d] hover:text-[#1a2332]">
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 19l-7-7 7-7" />
              </svg>
            </button>
            <p className="text-xs font-semibold text-gray-700">Filter: {selectedCol?.label}</p>
            {colRange && (
              <span className="ml-auto text-[10px] text-gray-400">
                {colRange.min.toFixed(1)} – {colRange.max.toFixed(1)}
              </span>
            )}
          </div>
          <div className="px-3 py-3 space-y-3">
            {selectedCol?.type === 'number' && (
              <>
                <div className="flex gap-2">
                  <div className="flex-1">
                    <label className="block text-[10px] text-gray-400 mb-1">Min</label>
                    <input
                      type="number"
                      value={numMin}
                      onChange={e => setNumMin(e.target.value)}
                      placeholder={colRange ? colRange.min.toFixed(1) : ''}
                      className="w-full text-xs border border-gray-200 rounded-lg px-2.5 py-1.5 focus:outline-none focus:ring-1 focus:ring-[#0f131d]"
                      step="0.1"
                    />
                  </div>
                  <div className="flex-1">
                    <label className="block text-[10px] text-gray-400 mb-1">Max</label>
                    <input
                      type="number"
                      value={numMax}
                      onChange={e => setNumMax(e.target.value)}
                      placeholder={colRange ? colRange.max.toFixed(1) : ''}
                      className="w-full text-xs border border-gray-200 rounded-lg px-2.5 py-1.5 focus:outline-none focus:ring-1 focus:ring-[#0f131d]"
                      step="0.1"
                    />
                  </div>
                </div>
              </>
            )}
            {selectedCol?.type === 'boolean' && (
              <div className="flex gap-2">
                {[true, false].map(v => (
                  <button
                    key={String(v)}
                    onClick={() => setBoolVal(v)}
                    className={`flex-1 py-1.5 rounded-lg text-xs font-medium border transition-colors ${
                      boolVal === v
                        ? 'bg-[#0f131d] text-white border-[#0f131d]'
                        : 'border-gray-200 text-gray-600 hover:bg-gray-50'
                    }`}
                  >
                    {v ? 'Yes' : 'No'}
                  </button>
                ))}
              </div>
            )}
            {selectedCol?.type === 'text' && (
              <input
                type="text"
                value={textVal}
                onChange={e => setTextVal(e.target.value)}
                placeholder={`Search ${selectedCol.label}...`}
                className="w-full text-xs border border-gray-200 rounded-lg px-2.5 py-1.5 focus:outline-none focus:ring-1 focus:ring-[#0f131d]"
                autoFocus
              />
            )}
          </div>
          <div className="px-3 pb-3 flex items-center justify-between">
            <button onClick={onClose} className="text-xs text-gray-400 hover:text-gray-600">Cancel</button>
            <button
              onClick={handleApply}
              className="px-3 py-1.5 bg-[#0f131d] text-white text-xs rounded-lg hover:bg-[#1a2332] transition-colors font-medium"
            >
              Apply filter
            </button>
          </div>
        </>
      )}
    </div>
  )
}

export default function FilterBar({ molecules = [], columns = [], onFilteredChange }) {
  const [activeQuickFilters, setActiveQuickFilters] = useState(new Set())
  const [customFilters, setCustomFilters] = useState([]) // [{id, label, test}]
  const [showAddFilter, setShowAddFilter] = useState(false)
  const [isExpanded, setIsExpanded] = useState(true)
  const addRef = useRef(null)

  // Compute filtered molecules whenever filters change
  const filteredMolecules = useMemo(() => {
    let result = molecules
    // Apply quick filters
    for (const qfId of activeQuickFilters) {
      const qf = QUICK_FILTERS.find(f => f.id === qfId)
      if (qf) result = result.filter(qf.test)
    }
    // Apply custom filters
    for (const cf of customFilters) {
      result = result.filter(cf.test)
    }
    return result
  }, [molecules, activeQuickFilters, customFilters])

  // Notify parent
  useEffect(() => {
    onFilteredChange && onFilteredChange(filteredMolecules)
  }, [filteredMolecules, onFilteredChange])

  const hasAnyFilter = activeQuickFilters.size > 0 || customFilters.length > 0
  const filterCount = activeQuickFilters.size + customFilters.length

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

  const clearAllFilters = useCallback(() => {
    setActiveQuickFilters(new Set())
    setCustomFilters([])
  }, [])

  const addCustomFilter = useCallback((filter) => {
    setCustomFilters(prev => [...prev, filter])
  }, [])

  return (
    <div className="bg-white rounded-xl border border-gray-100 shadow-sm overflow-hidden">
      {/* Header row */}
      <div className="flex items-center justify-between px-4 py-2.5 bg-gray-50 border-b border-gray-100">
        <button
          onClick={() => setIsExpanded(v => !v)}
          className="flex items-center gap-2 text-xs font-semibold text-gray-500 uppercase tracking-wide hover:text-[#0f131d] transition-colors"
        >
          <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
              d="M3 4a1 1 0 011-1h16a1 1 0 011 1v2a1 1 0 01-.293.707L13 13.414V19a1 1 0 01-.553.894l-4 2A1 1 0 017 21v-7.586L3.293 6.707A1 1 0 013 6V4z" />
          </svg>
          Filters
          {filterCount > 0 && (
            <span className="inline-flex items-center justify-center w-4 h-4 rounded-full bg-[#0f131d] text-white text-[9px] font-bold">
              {filterCount}
            </span>
          )}
          <ChevronDown className={`w-3 h-3 transition-transform ${isExpanded ? 'rotate-180' : ''}`} />
        </button>
        <div className="flex items-center gap-3 text-xs">
          <span className={hasAnyFilter ? 'text-[#0f131d] font-medium' : 'text-gray-400'}>
            Showing <strong className="tabular-nums">{filteredMolecules.length}</strong> of{' '}
            <strong className="tabular-nums">{molecules.length}</strong>
          </span>
          {hasAnyFilter && (
            <button
              onClick={clearAllFilters}
              className="text-gray-400 hover:text-red-500 transition-colors flex items-center gap-1"
            >
              Clear all
            </button>
          )}
        </div>
      </div>

      {isExpanded && (
        <div className="px-4 py-3 space-y-2.5">
          {/* Quick filter pills */}
          <div className="flex flex-wrap items-center gap-2">
            <span className="text-[10px] font-semibold text-gray-400 uppercase tracking-wider mr-1">Quick:</span>
            {QUICK_FILTERS.map(qf => {
              const isActive = activeQuickFilters.has(qf.id)
              // Count how many molecules pass this filter alone
              const count = molecules.filter(qf.test).length
              return (
                <button
                  key={qf.id}
                  onClick={() => toggleQuickFilter(qf.id)}
                  className={`flex items-center gap-1.5 px-2.5 py-1 rounded-full text-xs font-medium border transition-all duration-150 ${
                    isActive
                      ? 'bg-[#0f131d] text-white border-[#0f131d] shadow-sm'
                      : 'bg-white text-gray-600 border-gray-200 hover:border-[#0f131d] hover:text-[#0f131d]'
                  }`}
                >
                  {qf.label}
                  <span className={`text-[9px] tabular-nums ${isActive ? 'text-white/70' : 'text-gray-400'}`}>
                    {count}
                  </span>
                </button>
              )
            })}
          </div>

          {/* Active custom filter chips */}
          {customFilters.length > 0 && (
            <div className="flex flex-wrap items-center gap-2">
              <span className="text-[10px] font-semibold text-gray-400 uppercase tracking-wider mr-1">Active:</span>
              {customFilters.map(cf => (
                <span
                  key={cf.id}
                  className="inline-flex items-center gap-1.5 px-2.5 py-1 bg-blue-50 border border-blue-200 text-blue-700 rounded-full text-xs font-medium"
                >
                  {cf.label}
                  <button
                    onClick={() => removeCustomFilter(cf.id)}
                    className="ml-0.5 text-blue-400 hover:text-blue-700 transition-colors"
                    aria-label="Remove filter"
                  >
                    <XIcon className="w-2.5 h-2.5" />
                  </button>
                </span>
              ))}
            </div>
          )}

          {/* Add filter button */}
          <div className="relative" ref={addRef}>
            <button
              onClick={() => setShowAddFilter(v => !v)}
              className={`flex items-center gap-1.5 px-2.5 py-1 rounded-full text-xs border transition-all duration-150 ${
                showAddFilter
                  ? 'bg-[#0f131d] text-white border-[#0f131d]'
                  : 'text-[#0f131d] border-[#0f131d]/40 hover:bg-blue-50'
              }`}
            >
              <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M12 4v16m8-8H4" />
              </svg>
              Add filter
            </button>
            {showAddFilter && (
              <AddFilterDropdown
                columns={columns}
                molecules={molecules}
                onAdd={addCustomFilter}
                onClose={() => setShowAddFilter(false)}
              />
            )}
          </div>
        </div>
      )}
    </div>
  )
}
