import React, { useState, useRef, useEffect, useMemo } from 'react'
import { ALL_COLUMNS, GROUP_META, getEffectiveColorScale } from '../lib/columns.js'
import useSettingsStore from '../stores/settingsStore.js'

const MODES = [
  { value: 'none', label: 'None' },
  { value: 'higher-better', label: 'Higher ▲' },
  { value: 'lower-better', label: 'Lower ▼' },
  { value: 'custom', label: 'Custom' },
]

const BAND_COLORS = [
  { value: 'green', label: 'Green', tw: 'bg-green-400' },
  { value: 'yellow-green', label: 'Lime', tw: 'bg-lime-400' },
  { value: 'yellow', label: 'Yellow', tw: 'bg-yellow-400' },
  { value: 'orange', label: 'Orange', tw: 'bg-orange-400' },
  { value: 'red', label: 'Red', tw: 'bg-red-400' },
  { value: 'gray', label: 'Gray', tw: 'bg-gray-400' },
]

const DEFAULT_INTERVALS = [
  { min: null, max: 0, color: 'green' },
  { min: 0, max: 1, color: 'yellow' },
  { min: 1, max: null, color: 'red' },
]

function IntervalEditor({ intervals, onChange }) {
  const update = (idx, field, val) => {
    const next = intervals.map((band, i) =>
      i === idx ? { ...band, [field]: val } : band
    )
    onChange(next)
  }

  const addBand = () => {
    const last = intervals[intervals.length - 1]
    onChange([...intervals, { min: last?.max ?? 0, max: null, color: 'gray' }])
  }

  const removeBand = (idx) => {
    if (intervals.length <= 1) return
    onChange(intervals.filter((_, i) => i !== idx))
  }

  return (
    <div className="space-y-1.5 pl-2 border-l-2 border-gray-200 ml-1">
      {intervals.map((band, idx) => (
        <div key={idx} className="flex items-center gap-1.5 text-[11px]">
          <input
            type="number"
            value={band.min ?? ''}
            onChange={e => update(idx, 'min', e.target.value === '' ? null : Number(e.target.value))}
            placeholder="-∞"
            className="w-14 px-1 py-0.5 border border-gray-200 rounded text-center text-gray-600 focus:outline-none focus:ring-1 focus:ring-bx-mint"
          />
          <span className="text-gray-400">to</span>
          <input
            type="number"
            value={band.max ?? ''}
            onChange={e => update(idx, 'max', e.target.value === '' ? null : Number(e.target.value))}
            placeholder="+∞"
            className="w-14 px-1 py-0.5 border border-gray-200 rounded text-center text-gray-600 focus:outline-none focus:ring-1 focus:ring-bx-mint"
          />
          <select
            value={band.color}
            onChange={e => update(idx, 'color', e.target.value)}
            className="text-[10px] border border-gray-200 rounded px-1 py-0.5 bg-white"
          >
            {BAND_COLORS.map(c => (
              <option key={c.value} value={c.value}>{c.label}</option>
            ))}
          </select>
          <span className={`w-3 h-3 rounded-full flex-shrink-0 ${BAND_COLORS.find(c => c.value === band.color)?.tw || 'bg-gray-300'}`} />
          {intervals.length > 1 && (
            <button
              onClick={() => removeBand(idx)}
              className="text-gray-300 hover:text-red-400 transition-colors"
              title="Remove band"
            >
              <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
              </svg>
            </button>
          )}
        </div>
      ))}
      <button
        onClick={addBand}
        className="text-[10px] text-bx-light-text hover:underline flex items-center gap-0.5"
      >
        <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 4v16m8-8H4" />
        </svg>
        Add band
      </button>
    </div>
  )
}

export default function ColumnColorSettings({ isOpen, onClose, availableColumns }) {
  const panelRef = useRef(null)
  const [search, setSearch] = useState('')
  const [expandedKey, setExpandedKey] = useState(null)

  const columnColorOverrides = useSettingsStore(s => s.columnColorOverrides)
  const setColumnColor = useSettingsStore(s => s.setColumnColor)
  const resetColumnColor = useSettingsStore(s => s.resetColumnColor)

  // Close on click outside
  useEffect(() => {
    if (!isOpen) return
    const handler = (e) => {
      if (panelRef.current && !panelRef.current.contains(e.target)) {
        onClose()
      }
    }
    document.addEventListener('mousedown', handler)
    return () => document.removeEventListener('mousedown', handler)
  }, [isOpen, onClose])

  // Numeric columns grouped by category — filtered to available columns if provided
  const groups = useMemo(() => {
    const availableKeys = availableColumns ? new Set(availableColumns.map(c => c.key)) : null
    const numCols = ALL_COLUMNS.filter(c => c.type === 'number' && (!availableKeys || availableKeys.has(c.key)))
    const groupOrder = Object.keys(GROUP_META)
    const map = new Map()
    for (const col of numCols) {
      const g = col.group || 'molecule'
      if (!map.has(g)) map.set(g, [])
      map.get(g).push(col)
    }
    return groupOrder
      .filter(g => map.has(g))
      .map(g => ({ group: g, label: GROUP_META[g]?.label || g, meta: GROUP_META[g], cols: map.get(g) }))
  }, [availableColumns])

  const searchLower = search.toLowerCase().trim()

  const overrideCount = Object.keys(columnColorOverrides).length

  if (!isOpen) return null

  return (
    <div
      ref={panelRef}
      className="absolute left-0 top-full mt-1.5 z-50 bg-white border border-gray-200 rounded-xl shadow-xl"
      style={{ width: '380px' }}
    >
      {/* Header */}
      <div className="flex items-center justify-between px-3 pt-3 pb-2 border-b border-gray-100">
        <div className="flex items-center gap-2">
          <svg className="w-3.5 h-3.5 text-gray-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
              d="M7 21a4 4 0 01-4-4V5a2 2 0 012-2h4a2 2 0 012 2v12a4 4 0 01-4 4zm0 0h12a2 2 0 002-2v-4a2 2 0 00-2-2h-2.343M11 7.343l1.657-1.657a2 2 0 012.828 0l2.829 2.829a2 2 0 010 2.828l-8.486 8.485M7 17h.01" />
          </svg>
          <span className="text-xs font-bold text-gray-700">Column Colors</span>
          {overrideCount > 0 && (
            <span className="text-[10px] px-1.5 py-0.5 bg-bx-surface/10 text-bx-light-text rounded-full font-medium">
              {overrideCount} custom
            </span>
          )}
        </div>
        <button onClick={onClose} className="text-gray-400 hover:text-gray-600">
          <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
          </svg>
        </button>
      </div>

      {/* Search */}
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
            className="w-full text-xs border border-gray-200 rounded-lg pl-7 pr-2.5 py-1.5 focus:outline-none focus:ring-1 focus:ring-bx-mint text-gray-700 placeholder-gray-300"
            autoFocus
          />
        </div>
      </div>

      {/* Column list */}
      <div className="overflow-y-auto max-h-80 py-1">
        {groups.map(group => {
          const cols = group.cols.filter(col =>
            !searchLower ||
            col.label.toLowerCase().includes(searchLower) ||
            col.key.toLowerCase().includes(searchLower) ||
            group.label.toLowerCase().includes(searchLower)
          )
          if (!cols.length) return null

          return (
            <div key={group.group} className="border-b border-gray-50 last:border-b-0">
              <div className={`px-3 py-1.5 ${group.meta?.bg || 'bg-gray-50'}`}>
                <span className={`text-[10px] font-semibold uppercase tracking-wider ${group.meta?.text || 'text-gray-500'}`}>
                  {group.label}
                </span>
              </div>
              <div className="py-0.5">
                {cols.map(col => {
                  const effective = getEffectiveColorScale(col.key, columnColorOverrides)
                  const isOverridden = !!columnColorOverrides[col.key]
                  const isCustomExpanded = expandedKey === col.key && effective.mode === 'custom'

                  return (
                    <div key={col.key} className="px-3 py-1.5 hover:bg-gray-50 transition-colors">
                      <div className="flex items-center gap-2">
                        <span className={`text-xs flex-1 ${isOverridden ? 'text-gray-800 font-medium' : 'text-gray-600'}`}>
                          {col.label}
                          {col.unit && <span className="text-[9px] text-gray-300 ml-1">({col.unit})</span>}
                        </span>
                        {isOverridden && (
                          <button
                            onClick={() => resetColumnColor(col.key)}
                            className="text-[9px] text-gray-400 hover:text-red-500 transition-colors"
                            title="Reset to default"
                          >
                            Reset
                          </button>
                        )}
                        <select
                          value={effective.mode}
                          onChange={e => {
                            const mode = e.target.value
                            if (mode === 'custom') {
                              setColumnColor(col.key, { mode, intervals: effective.intervals || [...DEFAULT_INTERVALS] })
                              setExpandedKey(col.key)
                            } else if (mode === 'none' && !col.colorScale) {
                              resetColumnColor(col.key)
                            } else {
                              setColumnColor(col.key, { mode })
                              if (expandedKey === col.key) setExpandedKey(null)
                            }
                          }}
                          className="text-[11px] border border-gray-200 rounded px-1.5 py-0.5 bg-white text-gray-600 w-24"
                        >
                          {MODES.map(m => (
                            <option key={m.value} value={m.value}>{m.label}</option>
                          ))}
                        </select>
                        {effective.mode === 'custom' && (
                          <button
                            onClick={() => setExpandedKey(expandedKey === col.key ? null : col.key)}
                            className="text-gray-400 hover:text-gray-600"
                          >
                            <svg className={`w-3 h-3 transition-transform ${isCustomExpanded ? 'rotate-180' : ''}`}
                              fill="none" stroke="currentColor" viewBox="0 0 24 24">
                              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
                            </svg>
                          </button>
                        )}
                      </div>
                      {isCustomExpanded && (
                        <div className="mt-2">
                          <IntervalEditor
                            intervals={effective.intervals || [...DEFAULT_INTERVALS]}
                            onChange={intervals => setColumnColor(col.key, { mode: 'custom', intervals })}
                          />
                        </div>
                      )}
                    </div>
                  )
                })}
              </div>
            </div>
          )
        })}
      </div>

      {/* Footer */}
      <div className="px-3 py-2 border-t border-gray-100 flex items-center justify-between bg-gray-50/50 rounded-b-xl">
        <span className="text-[10px] text-gray-400">
          {overrideCount > 0
            ? <><strong className="text-gray-600">{overrideCount}</strong> overridden</>
            : 'All defaults'
          }
        </span>
        <button
          onClick={onClose}
          className="text-[10px] text-bx-light-text hover:underline font-medium"
        >
          Done
        </button>
      </div>
    </div>
  )
}
