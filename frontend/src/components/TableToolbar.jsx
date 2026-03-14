import React, { useState, useMemo, useEffect, useRef, useCallback } from 'react'
import { QUICK_FILTERS, InlineFilterConfig } from './FilterBar.jsx'
import ColumnSelector from './ColumnSelector.jsx'
import ColumnColorSettings from './ColumnColorSettings.jsx'
import BindXLogo from './BindXLogo.jsx'
import { ALL_COLUMNS } from '../lib/columns.js'

function XIcon({ className = 'w-3 h-3' }) {
  return (
    <svg className={className} fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M6 18L18 6M6 6l12 12" />
    </svg>
  )
}

// ---------------------------------------------------------------------------
// TableToolbar — unified toolbar replacing SelectionToolbar + FilterBar + Column controls
// 1 permanent row (filter + columns + count) + 1 conditional row (selection bar)
// ---------------------------------------------------------------------------
export default function TableToolbar({
  // Filter
  molecules = [],
  columns = [],
  onFilteredChange,
  onServerFilterChange,
  // Advanced filter (controlled by parent)
  showAdvancedFilter,
  onToggleAdvancedFilter,
  advancedFilterActive,
  onClearAdvancedFilter,
  // Column display
  visibleKeys,
  onVisibleKeysChange,
  // Count display (final, after all filters including advanced)
  filteredMoleculeCount,
  moleculeTotal,
  moleculesLoading,
  // Selection
  selectedCount = 0,
  totalCount = 0,
  bookmarkedCount = 0,
  filteredCount = null,
  onSelectAll,
  onSelectNone,
  onSelectBookmarked,
  onSelectFiltered,
  // Bulk actions
  onExport,
  onBookmarkSelected,
  onSendToNextPhase,
  isFrozen = false,
}) {
  const [showColorSettings, setShowColorSettings] = useState(false)

  // --- Filter state (absorbed from FilterBar) ---
  const [activeQuickFilters, setActiveQuickFilters] = useState(new Set())
  const [customFilters, setCustomFilters] = useState([])
  const [searchText, setSearchText] = useState('')
  const [configCol, setConfigCol] = useState(null)
  const addBtnRef = useRef(null)
  const debounceRef = useRef(null)

  // Filterable columns for the dropdown
  const filterableCols = useMemo(
    () => columns.filter(c => c.sortable && c.type !== 'action' && c.type !== 'invalidation' && c.type !== 'tags' && c.type !== 'editable_text'),
    [columns]
  )

  // Compute filtered molecules (basic filters only — parent applies advanced on top)
  // Local filtering: quick filters + custom column filters only.
  // Search text is handled server-side (sent via onServerFilterChange debounce).
  const filteredMolecules = useMemo(() => {
    let result = molecules
    for (const qfId of activeQuickFilters) {
      const qf = QUICK_FILTERS.find(f => f.id === qfId)
      if (qf) result = result.filter(qf.test)
    }
    for (const cf of customFilters) {
      result = result.filter(cf.test)
    }
    return result
  }, [molecules, activeQuickFilters, customFilters])

  // Notify parent of basic filter result
  useEffect(() => {
    onFilteredChange && onFilteredChange(filteredMolecules)
  }, [filteredMolecules, onFilteredChange])

  // Debounce search/bookmarked changes to the server (always — pagination is server-side)
  useEffect(() => {
    if (!onServerFilterChange) return
    if (debounceRef.current) clearTimeout(debounceRef.current)
    debounceRef.current = setTimeout(() => {
      onServerFilterChange({
        search: searchText.trim(),
        bookmarked_only: activeQuickFilters.has('bookmarked'),
      })
    }, 300)
    return () => { if (debounceRef.current) clearTimeout(debounceRef.current) }
  }, [searchText, activeQuickFilters, onServerFilterChange])

  const hasAnyFilter = activeQuickFilters.size > 0 || customFilters.length > 0 || searchText.trim() !== ''
  const hasAnyFilterOrAdvanced = hasAnyFilter || advancedFilterActive
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
    onClearAdvancedFilter && onClearAdvancedFilter()
  }, [onClearAdvancedFilter])

  const handleColSelect = useCallback((e) => {
    const key = e.target.value
    if (!key) return
    const col = filterableCols.find(c => c.key === key)
    if (col) setConfigCol(col)
    e.target.value = ''
  }, [filterableCols])

  return (
    <div className="space-y-1.5">
      {/* ── Row 1: Unified toolbar (always visible) ── */}
      <div className="card px-3 py-2">
        <div className="flex items-center gap-2 flex-wrap">
          {/* Search input */}
          <div className="relative flex-shrink-0">
            <svg className="absolute left-2 top-1/2 -translate-y-1/2 w-3.5 h-3.5 text-gray-400 pointer-events-none" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
            </svg>
            <input
              value={searchText}
              onChange={e => setSearchText(e.target.value)}
              placeholder="Search..."
              className="text-sm border border-gray-200 rounded-lg pl-7 pr-2 py-1.5 w-32 focus:outline-none focus:ring-1 focus:ring-bx-mint focus:w-44 transition-all"
            />
            {searchText && (
              <button onClick={() => setSearchText('')} className="absolute right-1.5 top-1/2 -translate-y-1/2 text-gray-400 hover:text-gray-600">
                <XIcon className="w-2.5 h-2.5" />
              </button>
            )}
          </div>

          {/* Quick filter toggle pills */}
          {QUICK_FILTERS.map(qf => {
            const isActive = activeQuickFilters.has(qf.id)
            const count = molecules.filter(qf.test).length
            if (count === 0) return null
            return (
              <button
                key={qf.id}
                onClick={() => toggleQuickFilter(qf.id)}
                className={`flex items-center gap-1 px-2 py-1 rounded-full text-xs font-medium border transition-all duration-150 flex-shrink-0 ${
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
          <div className="relative flex-shrink-0" ref={addBtnRef}>
            <select
              onChange={handleColSelect}
              defaultValue=""
              className="text-xs border border-dashed border-gray-300 rounded-lg pl-2 pr-6 py-1.5 text-gray-500 bg-white hover:border-gray-400 focus:outline-none focus:ring-1 focus:ring-bx-mint cursor-pointer appearance-none"
              style={{ backgroundImage: `url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 20 20' fill='%239ca3af'%3E%3Cpath fill-rule='evenodd' d='M5.23 7.21a.75.75 0 011.06.02L10 11.168l3.71-3.938a.75.75 0 111.08 1.04l-4.25 4.5a.75.75 0 01-1.08 0l-4.25-4.5a.75.75 0 01.02-1.06z'/%3E%3C/svg%3E")`, backgroundRepeat: 'no-repeat', backgroundPosition: 'right 4px center', backgroundSize: '14px' }}
            >
              <option value="">+ Filter</option>
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

          {/* Advanced filter button */}
          <button
            onClick={onToggleAdvancedFilter}
            className={`flex items-center gap-1.5 px-2.5 py-1.5 rounded-lg text-xs font-medium border transition-colors flex-shrink-0 ${
              showAdvancedFilter
                ? 'bg-bx-surface text-white border-bx-surface'
                : 'bg-white text-gray-500 border-gray-200 hover:border-gray-300 hover:text-gray-700'
            }`}
          >
            <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                d="M3 4a1 1 0 011-1h16a1 1 0 011 1v2.586a1 1 0 01-.293.707l-6.414 6.414a1 1 0 00-.293.707V17l-4 4v-6.586a1 1 0 00-.293-.707L3.293 7.293A1 1 0 013 6.586V4z" />
            </svg>
            Advanced
            {advancedFilterActive && (
              <span className="w-1.5 h-1.5 rounded-full bg-bx-cyan" />
            )}
          </button>

          <span className="text-gray-200 text-sm select-none">|</span>

          {/* Column selector + color settings */}
          <div className="flex items-center gap-1.5 flex-shrink-0">
            <ColumnSelector
              allColumns={ALL_COLUMNS}
              visibleKeys={visibleKeys}
              onChange={onVisibleKeysChange}
            />
            <div className="relative">
              <button
                onClick={() => setShowColorSettings(v => !v)}
                className={`flex items-center gap-1.5 px-2 py-1.5 rounded-lg text-xs font-medium border transition-all duration-150 ${
                  showColorSettings
                    ? 'border-bx-surface text-bx-light-text bg-blue-50 shadow-sm'
                    : 'border-gray-200 text-gray-500 hover:bg-gray-50 hover:border-gray-300'
                }`}
                title="Column color settings"
              >
                <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                    d="M7 21a4 4 0 01-4-4V5a2 2 0 012-2h4a2 2 0 012 2v12a4 4 0 01-4 4zm0 0h12a2 2 0 002-2v-4a2 2 0 00-2-2h-2.343M11 7.343l1.657-1.657a2 2 0 012.828 0l2.829 2.829a2 2 0 010 2.828l-8.486 8.485M7 17h.01" />
                </svg>
                Colors
              </button>
              <ColumnColorSettings
                isOpen={showColorSettings}
                onClose={() => setShowColorSettings(false)}
                availableColumns={columns}
              />
            </div>
          </div>

          {/* Spacer */}
          <div className="flex-1" />

          {/* Molecule count + clear */}
          <div className="flex items-center gap-2 flex-shrink-0">
            <span className="text-xs text-gray-400 tabular-nums flex items-center gap-1.5 whitespace-nowrap">
              {moleculesLoading && <BindXLogo variant="loading" size={16} />}
              {hasAnyFilterOrAdvanced
                ? <><span className="text-bx-light-text font-medium">{filteredMoleculeCount}</span> / {moleculeTotal.toLocaleString()} mol.</>
                : <>{moleculeTotal.toLocaleString()} mol.</>
              }
            </span>
            {hasAnyFilterOrAdvanced && (
              <button onClick={clearAll} className="text-gray-400 hover:text-red-500 transition-colors" title="Clear all filters">
                <XIcon className="w-3 h-3" />
              </button>
            )}
          </div>
        </div>

        {/* Custom filter chips (only shown when present) */}
        {customFilters.length > 0 && (
          <div className="flex flex-wrap items-center gap-1.5 mt-2 pt-2 border-t border-gray-100">
            {customFilters.map(cf => (
              <span key={cf.id}
                className="inline-flex items-center gap-1 px-2 py-0.5 bg-blue-50 border border-blue-200 text-blue-700 rounded-full text-xs font-medium"
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

      {/* ── Row 2: Selection action bar (conditional — only when selectedCount > 0) ── */}
      {selectedCount > 0 && (
        <div className="card px-3 py-2 bg-blue-50/40 border-l-4 border-l-bx-accent">
          <div className="flex items-center gap-2 flex-wrap">
            {/* Selection state */}
            <div className="flex items-center gap-1.5 text-sm font-semibold text-bx-light-text">
              <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M5 13l4 4L19 7" />
              </svg>
              {selectedCount} selected
            </div>

            <span className="text-gray-300 text-sm select-none">|</span>

            {/* Quick select buttons */}
            <div className="flex items-center gap-2 text-xs">
              <button onClick={onSelectAll}
                className="text-gray-500 hover:text-bx-light-text hover:underline underline-offset-2 transition-colors">
                All ({totalCount})
              </button>

              {bookmarkedCount > 0 && (
                <button onClick={onSelectBookmarked}
                  className="text-gray-500 hover:text-yellow-600 hover:underline underline-offset-2 transition-colors flex items-center gap-0.5">
                  <svg className="w-2.5 h-2.5 text-yellow-400" fill="currentColor" viewBox="0 0 24 24">
                    <path d="M5 5a2 2 0 012-2h10a2 2 0 012 2v16l-7-3.5L5 21V5z" />
                  </svg>
                  Bookmarked ({bookmarkedCount})
                </button>
              )}

              {(filterCount > 0 || advancedFilterActive) && filteredCount != null && (
                <button onClick={onSelectFiltered}
                  className="text-blue-600 hover:text-bx-light-text hover:underline underline-offset-2 transition-colors">
                  Filtered ({filteredCount})
                </button>
              )}

              <button onClick={onSelectNone}
                className="text-gray-400 hover:text-red-500 transition-colors flex items-center gap-0.5">
                <svg className="w-2.5 h-2.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M6 18L18 6M6 6l12 12" />
                </svg>
                Clear
              </button>
            </div>

            {/* Spacer */}
            <div className="flex-1" />

            {/* Bulk actions */}
            <div className="flex items-center gap-2">
              <button
                onClick={onExport}
                className="flex items-center gap-1 px-2.5 py-1 rounded-lg text-xs font-medium border border-gray-200 text-gray-600 hover:bg-gray-50 hover:border-gray-300 transition-all"
              >
                <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                    d="M12 10v6m0 0l-3-3m3 3l3-3m2 8H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
                </svg>
                Export
              </button>

              <div title={isFrozen ? 'Phase frozen' : 'Bookmark selected'}>
                <button
                  onClick={!isFrozen ? onBookmarkSelected : undefined}
                  disabled={isFrozen}
                  className={`flex items-center gap-1 px-2.5 py-1 rounded-lg text-xs font-medium border transition-all ${
                    isFrozen
                      ? 'border-gray-200 text-gray-300 cursor-not-allowed bg-gray-50'
                      : 'border-yellow-300 text-yellow-700 hover:bg-yellow-50 hover:border-yellow-400'
                  }`}
                >
                  <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                      d="M5 5a2 2 0 012-2h10a2 2 0 012 2v16l-7-3.5L5 21V5z" />
                  </svg>
                  Bookmark ({selectedCount})
                </button>
              </div>


            </div>
          </div>
        </div>
      )}
    </div>
  )
}
