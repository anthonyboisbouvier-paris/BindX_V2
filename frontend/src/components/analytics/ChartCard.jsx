import React, { useState, useRef, useEffect, useMemo } from 'react'
import CHART_REGISTRY, { CHART_TYPES } from './chartRegistry.js'
import { ALL_COLUMNS } from '../../lib/columns.js'

function useContainerWidth() {
  const ref = useRef(null)
  const [width, setWidth] = useState(0)
  useEffect(() => {
    if (!ref.current) return
    const observer = new ResizeObserver(entries => {
      const w = entries[0]?.contentRect?.width || 0
      if (w > 0) setWidth(Math.floor(w))
    })
    observer.observe(ref.current)
    setWidth(Math.floor(ref.current.getBoundingClientRect().width))
    return () => observer.disconnect()
  }, [])
  return [ref, width]
}

/** Build a human-readable title from the chart config */
function autoTitle(config) {
  const colLabel = (key) => ALL_COLUMNS.find(c => c.key === key)?.label || key
  switch (config.type) {
    case 'histogram':   return config.groupKey ? `${colLabel(config.key)} by ${colLabel(config.groupKey)}` : `${colLabel(config.key)} Distribution`
    case 'bar':         return `${colLabel(config.key)} Breakdown`
    case 'pie':         return `${colLabel(config.key)} Proportions`
    case 'scatter':     return config.colorKey ? `${colLabel(config.xKey)} vs ${colLabel(config.yKey)} (color: ${colLabel(config.colorKey)})` : `${colLabel(config.xKey)} vs ${colLabel(config.yKey)}`
    case 'box':         return `${colLabel(config.metricKey)} by ${colLabel(config.groupKey)}`
    case 'correlation': return 'Correlation Matrix'
    case 'bubble':      return config.colorKey ? `${colLabel(config.xKey)} vs ${colLabel(config.yKey)} (size: ${colLabel(config.sizeKey)}, color: ${colLabel(config.colorKey)})` : `${colLabel(config.xKey)} vs ${colLabel(config.yKey)} (size: ${colLabel(config.sizeKey)})`
    case 'topn':        return `Top ${config.n || 15} by ${colLabel(config.key)}`
    default:            return config.type
  }
}

export default function ChartCard({ config, molecules, onUpdate, onRemove, selectableColumns }) {
  const [editing, setEditing] = useState(false)
  const [fullscreen, setFullscreen] = useState(false)
  const [containerRef, containerWidth] = useContainerWidth()
  const [fsWidth, setFsWidth] = useState(0)

  // Close fullscreen on Escape + measure fullscreen width
  useEffect(() => {
    if (!fullscreen) return
    const onKey = (e) => { if (e.key === 'Escape') setFullscreen(false) }
    window.addEventListener('keydown', onKey)
    // Compute fullscreen content width (95vw - padding - scrollbar)
    setFsWidth(Math.floor(window.innerWidth * 0.95 - 48))
    const onResize = () => setFsWidth(Math.floor(window.innerWidth * 0.95 - 48))
    window.addEventListener('resize', onResize)
    return () => { window.removeEventListener('keydown', onKey); window.removeEventListener('resize', onResize) }
  }, [fullscreen])

  const entry = CHART_REGISTRY[config.type]
  if (!entry) return null

  const ChartComponent = entry.component
  const title = autoTitle(config)

  const cardContent = (isFs) => {
    const w = isFs ? fsWidth : containerWidth
    const h = isFs ? Math.floor(window.innerHeight * 0.75) : 220
    return (
      <>
        {/* Header */}
        <div className="flex items-center justify-between px-3 py-2 border-b border-gray-100">
          <h4 className={`${isFs ? 'text-sm' : 'text-xs'} font-semibold text-gray-600 truncate`}>{title}</h4>
          <div className="flex items-center gap-1">
            {!isFs && (
              <button
                onClick={() => setEditing(v => !v)}
                className="p-1 rounded hover:bg-gray-100 text-gray-400 hover:text-gray-600 transition-colors"
                title="Configure"
              >
                <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                    d="M10.325 4.317c.426-1.756 2.924-1.756 3.35 0a1.724 1.724 0 002.573 1.066c1.543-.94 3.31.826 2.37 2.37a1.724 1.724 0 001.066 2.573c1.756.426 1.756 2.924 0 3.35a1.724 1.724 0 00-1.066 2.573c.94 1.543-.826 3.31-2.37 2.37a1.724 1.724 0 00-2.573 1.066c-.426 1.756-2.924 1.756-3.35 0a1.724 1.724 0 00-2.573-1.066c-1.543.94-3.31-.826-2.37-2.37a1.724 1.724 0 00-1.066-2.573c-1.756-.426-1.756-2.924 0-3.35a1.724 1.724 0 001.066-2.573c-.94-1.543.826-3.31 2.37-2.37.996.608 2.296.07 2.572-1.065z" />
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M15 12a3 3 0 11-6 0 3 3 0 016 0z" />
                </svg>
              </button>
            )}
            <button
              onClick={() => setFullscreen(v => !v)}
              className="p-1 rounded hover:bg-gray-100 text-gray-400 hover:text-gray-600 transition-colors"
              title={isFs ? 'Exit fullscreen' : 'Fullscreen'}
            >
              {isFs ? (
                <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                    d="M9 9L4 4m0 0v4m0-4h4m6 6l5 5m0 0v-4m0 4h-4M9 15l-5 5m0 0v-4m0 4h4m6-6l5-5m0 0v4m0-4h-4" />
                </svg>
              ) : (
                <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                    d="M4 8V4m0 0h4M4 4l5 5m11-5h-4m4 0v4m0 0l-5-5m-7 14H4m0 0v-4m0 4l5-5m11 5h-4m4 0v-4m0 0l-5 5" />
                </svg>
              )}
            </button>
            {!isFs && (
              <button
                onClick={() => onRemove(config.id)}
                className="p-1 rounded hover:bg-red-50 text-gray-400 hover:text-red-500 transition-colors"
                title="Remove chart"
              >
                <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
                </svg>
              </button>
            )}
          </div>
        </div>

        {/* Config editor (inline card only) */}
        {!isFs && editing && (
          <ChartConfigEditor config={config} onUpdate={onUpdate} onClose={() => setEditing(false)}
            selectableColumns={selectableColumns} />
        )}

        {/* Chart */}
        <div ref={isFs ? undefined : containerRef} className="p-2" style={{ width: '100%', minHeight: isFs ? h + 16 : (config.type === 'correlation' ? 'auto' : 236) }}>
          {w > 0 ? (
            <ChartComponent config={config} molecules={molecules} width={w - 16} height={h} />
          ) : (
            <div style={{ height: h }} />
          )}
        </div>
      </>
    )
  }

  return (
    <>
      <div className="card overflow-hidden" style={{ minWidth: 0 }}>
        {cardContent(false)}
      </div>

      {/* Fullscreen overlay */}
      {fullscreen && (
        <div className="fixed inset-0 z-50 bg-black/50 flex items-center justify-center p-6"
          onClick={(e) => { if (e.target === e.currentTarget) setFullscreen(false) }}>
          <div className="bg-white rounded-xl shadow-2xl w-full max-w-[95vw] max-h-[95vh] overflow-auto">
            {cardContent(true)}
          </div>
        </div>
      )}
    </>
  )
}

function ChartConfigEditor({ config, onUpdate, onClose, selectableColumns = [] }) {
  // Continuous numeric columns (for axes, metrics, histograms)
  const numericCols = useMemo(
    () => selectableColumns.filter(c => c.type === 'number'),
    [selectableColumns]
  )
  // Categorical columns (for grouping, bar, pie — text, boolean, source)
  const categoricalCols = useMemo(
    () => selectableColumns.filter(c => c.type === 'text' || c.type === 'boolean' || c.type === 'source'),
    [selectableColumns]
  )

  return (
    <div className="px-3 py-2 bg-gray-50 border-b border-gray-100 space-y-2">
      {/* Type */}
      <div className="flex items-center gap-2">
        <label className="text-xs text-gray-500 w-12">Type</label>
        <select
          className="flex-1 text-xs border border-gray-200 rounded px-2 py-1 bg-white"
          value={config.type}
          onChange={e => onUpdate(config.id, { type: e.target.value })}
        >
          {CHART_TYPES.map(t => <option key={t.type} value={t.type}>{t.label}</option>)}
        </select>
      </div>

      {/* Type-specific config */}
      {(config.type === 'histogram') && (
        <>
          <ColumnPicker label="Metric" value={config.key} columns={numericCols}
            onChange={key => onUpdate(config.id, { key })} />
          <ColumnPicker label="Group" value={config.groupKey} columns={categoricalCols}
            onChange={groupKey => onUpdate(config.id, { groupKey })} optional />
          <div className="flex items-center gap-2">
            <label className="text-xs text-gray-500 w-12">Bins</label>
            <input type="number" min={5} max={50} className="w-16 text-xs border border-gray-200 rounded px-2 py-1 bg-white"
              value={config.bins || 20} onChange={e => onUpdate(config.id, { bins: Number(e.target.value) })} />
          </div>
        </>
      )}

      {(config.type === 'bar' || config.type === 'pie') && (
        <ColumnPicker label="Key" value={config.key} columns={categoricalCols}
          onChange={key => onUpdate(config.id, { key })} />
      )}

      {config.type === 'scatter' && (
        <>
          <ColumnPicker label="X" value={config.xKey} columns={numericCols}
            onChange={xKey => onUpdate(config.id, { xKey })} />
          <ColumnPicker label="Y" value={config.yKey} columns={numericCols}
            onChange={yKey => onUpdate(config.id, { yKey })} />
          <ColumnPicker label="Color" value={config.colorKey} columns={numericCols}
            onChange={colorKey => onUpdate(config.id, { colorKey })} optional />
          <ScaleConfig config={config} onUpdate={onUpdate} />
        </>
      )}

      {config.type === 'box' && (
        <>
          <ColumnPicker label="Metric" value={config.metricKey} columns={numericCols}
            onChange={metricKey => onUpdate(config.id, { metricKey })} />
          <ColumnPicker label="Group" value={config.groupKey} columns={categoricalCols}
            onChange={groupKey => onUpdate(config.id, { groupKey })} />
        </>
      )}

      {config.type === 'correlation' && (
        <CorrelationColumnPicker
          selectableColumns={selectableColumns}
          selectedKeys={config.corrColumns}
          onChange={corrColumns => onUpdate(config.id, { corrColumns })}
        />
      )}

      {config.type === 'bubble' && (
        <>
          <ColumnPicker label="X" value={config.xKey} columns={numericCols}
            onChange={xKey => onUpdate(config.id, { xKey })} />
          <ColumnPicker label="Y" value={config.yKey} columns={numericCols}
            onChange={yKey => onUpdate(config.id, { yKey })} />
          <ColumnPicker label="Size" value={config.sizeKey} columns={numericCols}
            onChange={sizeKey => onUpdate(config.id, { sizeKey })} />
          <ColumnPicker label="Color" value={config.colorKey} columns={numericCols}
            onChange={colorKey => onUpdate(config.id, { colorKey })} optional />
        </>
      )}

      {config.type === 'topn' && (
        <>
          <ColumnPicker label="Metric" value={config.key} columns={numericCols}
            onChange={key => onUpdate(config.id, { key })} />
          <div className="flex items-center gap-2">
            <label className="text-xs text-gray-500 w-12">Top N</label>
            <input type="number" min={5} max={50} className="w-16 text-xs border border-gray-200 rounded px-2 py-1 bg-white"
              value={config.n || 15} onChange={e => onUpdate(config.id, { n: Number(e.target.value) })} />
          </div>
        </>
      )}

      <button onClick={onClose} className="text-xs text-bx-light-text hover:underline">Done</button>
    </div>
  )
}

function CorrelationColumnPicker({ selectableColumns, selectedKeys, onChange }) {
  const numCols = useMemo(
    () => selectableColumns.filter(c => c.type === 'number'),
    [selectableColumns]
  )
  // undefined/null = all selected
  const activeKeys = selectedKeys || numCols.map(c => c.key)
  const activeSet = new Set(activeKeys)

  const toggle = (key) => {
    const next = activeSet.has(key)
      ? activeKeys.filter(k => k !== key)
      : [...activeKeys, key]
    // If all selected, store undefined (= all)
    onChange(next.length === numCols.length ? undefined : next)
  }

  const allChecked = activeKeys.length === numCols.length
  const toggleAll = () => onChange(allChecked ? [] : undefined)

  return (
    <div className="space-y-1">
      <div className="flex items-center gap-2">
        <label className="text-xs text-gray-500">Columns</label>
        <button onClick={toggleAll} className="text-xs text-blue-500 hover:underline">
          {allChecked ? 'Deselect all' : 'Select all'}
        </button>
      </div>
      <div className="flex flex-wrap gap-x-3 gap-y-0.5 max-h-24 overflow-y-auto">
        {numCols.map(c => (
          <label key={c.key} className="flex items-center gap-1 text-xs text-gray-600 cursor-pointer">
            <input type="checkbox" checked={activeSet.has(c.key)} onChange={() => toggle(c.key)}
              className="w-3 h-3 rounded border-gray-300" />
            {c.label}
          </label>
        ))}
      </div>
    </div>
  )
}

function ScaleConfig({ config, onUpdate }) {
  const isManual = config.scaleMode === 'manual'
  const inputCls = "w-20 text-xs border border-gray-200 rounded px-2 py-1 bg-white"

  return (
    <div className="space-y-1">
      <div className="flex items-center gap-2">
        <label className="text-xs text-gray-500 w-12">Scale</label>
        <select
          className="flex-1 text-xs border border-gray-200 rounded px-2 py-1 bg-white"
          value={config.scaleMode || 'auto'}
          onChange={e => {
            const mode = e.target.value
            onUpdate(config.id, { scaleMode: mode === 'auto' ? undefined : mode })
          }}
        >
          <option value="auto">Auto</option>
          <option value="manual">Manual</option>
        </select>
      </div>
      {isManual && (
        <div className="flex items-center gap-2 flex-wrap">
          <label className="text-xs text-gray-400 w-12">X</label>
          <input type="number" step="any" placeholder="min" className={inputCls}
            value={config.xDomainMin ?? ''} onChange={e => onUpdate(config.id, { xDomainMin: e.target.value === '' ? undefined : Number(e.target.value) })} />
          <span className="text-xs text-gray-400">→</span>
          <input type="number" step="any" placeholder="max" className={inputCls}
            value={config.xDomainMax ?? ''} onChange={e => onUpdate(config.id, { xDomainMax: e.target.value === '' ? undefined : Number(e.target.value) })} />
          <label className="text-xs text-gray-400 w-12 ml-2">Y</label>
          <input type="number" step="any" placeholder="min" className={inputCls}
            value={config.yDomainMin ?? ''} onChange={e => onUpdate(config.id, { yDomainMin: e.target.value === '' ? undefined : Number(e.target.value) })} />
          <span className="text-xs text-gray-400">→</span>
          <input type="number" step="any" placeholder="max" className={inputCls}
            value={config.yDomainMax ?? ''} onChange={e => onUpdate(config.id, { yDomainMax: e.target.value === '' ? undefined : Number(e.target.value) })} />
        </div>
      )}
    </div>
  )
}

function ColumnPicker({ label, value, columns, onChange, optional }) {
  return (
    <div className="flex items-center gap-2">
      <label className="text-xs text-gray-500 w-12">{label}</label>
      <select
        className="flex-1 text-xs border border-gray-200 rounded px-2 py-1 bg-white"
        value={value || ''}
        onChange={e => onChange(e.target.value || undefined)}
      >
        <option value="">{optional ? '— None —' : '--'}</option>
        {columns.map(c => <option key={c.key} value={c.key}>{c.label}</option>)}
      </select>
    </div>
  )
}
