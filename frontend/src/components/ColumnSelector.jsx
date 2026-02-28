import React, { useState, useRef, useEffect } from 'react'

// ---------------------------------------------------------------------------
// Column groups for organized display
// ---------------------------------------------------------------------------
const GROUPS = [
  { label: 'Identity',       keys: ['name', 'smiles', 'source_run_id'] },
  { label: 'Docking Scores', keys: ['docking_score', 'cnn_score', 'cnn_affinity', 'cnn_vs'] },
  { label: 'Drug Properties',keys: ['logP', 'MW', 'HBD', 'HBA', 'TPSA', 'QED', 'lipinski_pass'] },
  { label: 'ADMET',          keys: ['solubility', 'BBB', 'hERG', 'metabolic_stability'] },
  { label: 'Scoring',        keys: ['composite_score'] },
  { label: 'Analysis',       keys: ['cluster_id', 'scaffold', 'interactions_count'] },
  { label: 'Generation',     keys: ['generation_level', 'parent_molecule_id'] },
]

const TYPE_COLORS = {
  number:  'text-blue-500 bg-blue-50',
  text:    'text-gray-500 bg-gray-100',
  boolean: 'text-green-600 bg-green-50',
  smiles:  'text-purple-500 bg-purple-50',
}

// ---------------------------------------------------------------------------
// ColumnSelector — polished dropdown to toggle visible columns
// Props:
//   allColumns   — full column definitions array
//   visibleKeys  — currently visible column key array
//   onChange     — callback(newVisibleKeysArray)
// ---------------------------------------------------------------------------
export default function ColumnSelector({ allColumns = [], visibleKeys = [], onChange }) {
  const [open, setOpen] = useState(false)
  const [search, setSearch] = useState('')
  const [collapsedGroups, setCollapsedGroups] = useState(new Set())
  const panelRef = useRef(null)
  const buttonRef = useRef(null)

  const visibleSet = new Set(visibleKeys)

  // Close panel on outside click
  useEffect(() => {
    if (!open) return
    function handle(e) {
      if (
        panelRef.current && !panelRef.current.contains(e.target) &&
        buttonRef.current && !buttonRef.current.contains(e.target)
      ) {
        setOpen(false)
        setSearch('')
      }
    }
    document.addEventListener('mousedown', handle)
    return () => document.removeEventListener('mousedown', handle)
  }, [open])

  function toggle(key) {
    const next = new Set(visibleSet)
    if (next.has(key)) next.delete(key)
    else next.add(key)
    onChange && onChange([...next])
  }

  function toggleGroup(groupKeys, show) {
    const next = new Set(visibleSet)
    groupKeys.forEach(k => show ? next.add(k) : next.delete(k))
    onChange && onChange([...next])
  }

  function toggleGroupCollapse(label) {
    setCollapsedGroups(prev => {
      const next = new Set(prev)
      if (next.has(label)) next.delete(label)
      else next.add(label)
      return next
    })
  }

  function showAll() {
    onChange && onChange(allColumns.map(c => c.key))
  }

  function hideAll() {
    // Always keep 'name' visible
    onChange && onChange(['name'])
  }

  // Filter by search
  const searchLower = search.toLowerCase().trim()
  const getGroupKeys = (group) =>
    group.keys.filter(k => {
      if (!searchLower) return true
      const col = allColumns.find(c => c.key === k)
      return col && (
        col.label.toLowerCase().includes(searchLower) ||
        col.key.toLowerCase().includes(searchLower)
      )
    })

  return (
    <div className="relative">
      <button
        ref={buttonRef}
        onClick={() => { setOpen(v => !v); if (!open) setSearch('') }}
        className={`flex items-center gap-1.5 px-3 py-1.5 rounded-lg text-sm border transition-all duration-150 ${
          open
            ? 'border-[#0f131d] text-[#0f131d] bg-blue-50 shadow-sm'
            : 'border-gray-200 text-gray-600 hover:bg-gray-50 hover:border-gray-300'
        }`}
      >
        <svg className="w-3.5 h-3.5 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
            d="M9 17V7m0 10a2 2 0 01-2 2H5a2 2 0 01-2-2V7a2 2 0 012-2h2a2 2 0 012 2m0 10a2 2 0 002 2h2a2 2 0 002-2M9 7a2 2 0 012-2h2a2 2 0 012 2m0 10V7m0 10a2 2 0 002 2h2a2 2 0 002-2V7a2 2 0 00-2-2h-2a2 2 0 00-2 2" />
        </svg>
        <span className="font-medium">Columns</span>
        <span className="text-[10px] font-normal tabular-nums px-1.5 py-0.5 bg-gray-100 rounded-full text-gray-500">
          {visibleSet.size}/{allColumns.length}
        </span>
        <svg
          className={`w-3 h-3 transition-transform duration-150 ${open ? 'rotate-180' : ''}`}
          fill="none" stroke="currentColor" viewBox="0 0 24 24"
        >
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
        </svg>
      </button>

      {open && (
        <div
          ref={panelRef}
          className="absolute left-0 top-full mt-1.5 z-50 bg-white border border-gray-200 rounded-xl shadow-xl w-76"
          style={{ width: '300px' }}
        >
          {/* Panel header */}
          <div className="flex items-center justify-between px-3 pt-3 pb-2 border-b border-gray-100">
            <span className="text-xs font-bold text-gray-700">Column Visibility</span>
            <div className="flex items-center gap-2">
              <button
                onClick={hideAll}
                className="text-[10px] text-gray-400 hover:text-gray-600 transition-colors"
              >
                Hide all
              </button>
              <span className="text-gray-200 text-xs">|</span>
              <button
                onClick={showAll}
                className="text-[10px] text-[#0f131d] hover:text-[#1a2332] transition-colors font-medium"
              >
                Show all
              </button>
            </div>
          </div>

          {/* Search input */}
          <div className="px-3 py-2 border-b border-gray-100">
            <div className="relative">
              <svg className="absolute left-2.5 top-1/2 -translate-y-1/2 w-3 h-3 text-gray-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
              </svg>
              <input
                type="text"
                value={search}
                onChange={e => setSearch(e.target.value)}
                placeholder="Search columns..."
                className="w-full text-xs border border-gray-200 rounded-lg pl-7 pr-2.5 py-1.5 focus:outline-none focus:ring-1 focus:ring-[#0f131d] text-gray-700 placeholder-gray-300"
                autoFocus
              />
              {search && (
                <button
                  onClick={() => setSearch('')}
                  className="absolute right-2 top-1/2 -translate-y-1/2 text-gray-300 hover:text-gray-500"
                >
                  <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M6 18L18 6M6 6l12 12" />
                  </svg>
                </button>
              )}
            </div>
          </div>

          {/* Column groups */}
          <div className="overflow-y-auto max-h-72 py-1">
            {GROUPS.map(group => {
              const keys = getGroupKeys(group)
              if (!keys.length) return null
              const isCollapsed = collapsedGroups.has(group.label) && !searchLower
              const groupVisible = keys.filter(k => visibleSet.has(k)).length
              const allGroupVisible = groupVisible === keys.length
              const someGroupVisible = groupVisible > 0 && groupVisible < keys.length

              return (
                <div key={group.label} className="border-b border-gray-50 last:border-b-0">
                  {/* Group header */}
                  <div className="flex items-center gap-1.5 px-3 py-1.5 bg-gray-50/80">
                    <button
                      onClick={() => toggleGroupCollapse(group.label)}
                      className="flex-1 flex items-center gap-1.5 text-left"
                      disabled={!!searchLower}
                    >
                      <svg
                        className={`w-3 h-3 text-gray-400 transition-transform ${isCollapsed ? '-rotate-90' : ''}`}
                        fill="none" stroke="currentColor" viewBox="0 0 24 24"
                      >
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
                      </svg>
                      <span className="text-[10px] font-semibold text-gray-500 uppercase tracking-wider">
                        {group.label}
                      </span>
                      <span className="text-[9px] text-gray-400 tabular-nums">
                        ({groupVisible}/{keys.length})
                      </span>
                    </button>
                    <div className="flex items-center gap-1">
                      <button
                        onClick={() => toggleGroup(keys, true)}
                        className="text-[9px] text-gray-400 hover:text-[#0f131d] transition-colors px-1"
                        title="Show all in group"
                      >
                        All
                      </button>
                      <button
                        onClick={() => toggleGroup(keys, false)}
                        className="text-[9px] text-gray-400 hover:text-red-500 transition-colors px-1"
                        title="Hide all in group"
                      >
                        None
                      </button>
                    </div>
                  </div>

                  {/* Column items */}
                  {!isCollapsed && (
                    <div className="py-0.5">
                      {keys.map(k => {
                        const col = allColumns.find(c => c.key === k)
                        if (!col) return null
                        const checked = visibleSet.has(k)
                        const typeColor = TYPE_COLORS[col.type] || TYPE_COLORS.text
                        return (
                          <label
                            key={k}
                            className={`flex items-center gap-2 px-3 py-1.5 cursor-pointer transition-colors hover:bg-gray-50 ${
                              checked ? 'opacity-100' : 'opacity-60'
                            }`}
                          >
                            <input
                              type="checkbox"
                              checked={checked}
                              onChange={() => toggle(k)}
                              className="accent-[#0f131d] w-3.5 h-3.5 flex-shrink-0"
                            />
                            <span className={`text-xs flex-1 ${checked ? 'text-gray-700 font-medium' : 'text-gray-500'}`}>
                              {col.label}
                            </span>
                            <div className="flex items-center gap-1.5">
                              {col.unit && (
                                <span className="text-[9px] text-gray-300">({col.unit})</span>
                              )}
                              <span className={`text-[8px] font-semibold px-1 py-0.5 rounded uppercase ${typeColor}`}>
                                {col.type === 'boolean' ? 'bool' : col.type}
                              </span>
                            </div>
                          </label>
                        )
                      })}
                    </div>
                  )}
                </div>
              )
            })}
            {searchLower && GROUPS.every(g => getGroupKeys(g).length === 0) && (
              <div className="py-6 text-center text-xs text-gray-400">
                No columns matching "{search}"
              </div>
            )}
          </div>

          {/* Footer */}
          <div className="px-3 py-2 border-t border-gray-100 flex items-center justify-between bg-gray-50/50 rounded-b-xl">
            <span className="text-[10px] text-gray-400">
              <strong className="text-gray-600 tabular-nums">{visibleSet.size}</strong> of {allColumns.length} visible
            </span>
            <button
              onClick={() => { setOpen(false); setSearch('') }}
              className="text-[10px] text-[#0f131d] hover:underline font-medium"
            >
              Done
            </button>
          </div>
        </div>
      )}
    </div>
  )
}
