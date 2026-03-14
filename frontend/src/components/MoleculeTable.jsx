import React, { useState, useMemo, useRef, useCallback, useEffect } from 'react'
import { createPortal } from 'react-dom'
import { useVirtualizer } from '@tanstack/react-virtual'
import InfoTip, { TIPS } from './InfoTip.jsx'
import { GROUP_META, getEffectiveColorScale } from '../lib/columns.js'
import useSettingsStore from '../stores/settingsStore.js'

// ---------------------------------------------------------------------------
// Number formatter per column key
// ---------------------------------------------------------------------------
function formatNumber(key, value) {
  if (value == null || typeof value !== 'number') return null
  // Composite score: backend stores 0-1, display as 0-100
  if (key === 'composite_score') return Math.round(value * 100).toString()
  // Half-life: clamp negative predictions to 0
  if (key === 'half_life' && value < 0) value = 0
  // Key-specific precision
  const decimals = {
    docking_score: 1,
    cnn_score: 2,
    cnn_affinity: 1,
    cnn_vs: 2,
    logP: 1,
    MW: 0,
    HBD: 0,
    HBA: 0,
    TPSA: 1,
    QED: 2,
    solubility: 2,
    BBB: 2,
    hERG: 2,
    half_life: 1,
    cyp_inhibitions: 0,
    cluster_id: 0,
    generation_level: 0,
    interactions_count: 0,
  }
  const d = decimals[key] !== undefined ? decimals[key] : 2
  return d === 0 ? Math.round(value).toString() : value.toFixed(d)
}

// ---------------------------------------------------------------------------
// Color scale: returns Tailwind classes for a value in a column
// colorConfig: { mode, intervals? } from getEffectiveColorScale
// ---------------------------------------------------------------------------
const CUSTOM_COLOR_MAP = {
  green:        { cell: 'bg-green-50',  text: 'text-green-700 font-semibold' },
  'yellow-green': { cell: 'bg-lime-50', text: 'text-lime-700' },
  yellow:       { cell: 'bg-yellow-50', text: 'text-yellow-700' },
  orange:       { cell: 'bg-orange-50', text: 'text-orange-600' },
  red:          { cell: 'bg-red-50',    text: 'text-red-600' },
  gray:         { cell: 'bg-gray-50',   text: 'text-gray-500' },
}

function getValueColorClasses(value, allValues, colorConfig, globalRange) {
  if (value == null) return { cell: '', text: '' }

  // Accept string for backward compat (scale name) or object { mode, intervals }
  const config = typeof colorConfig === 'string' ? { mode: colorConfig } : colorConfig
  if (!config || config.mode === 'none') return { cell: '', text: '' }

  // Custom interval mode: match first interval where value fits
  if (config.mode === 'custom' && config.intervals?.length) {
    for (const band of config.intervals) {
      const aboveMin = band.min == null || value >= band.min
      const belowMax = band.max == null || value < band.max
      if (aboveMin && belowMax) {
        return CUSTOM_COLOR_MAP[band.color] || { cell: '', text: '' }
      }
    }
    return { cell: '', text: '' }
  }

  // Compute rank: use global min/max if available, else fall back to page-local percentile
  let rank
  if (globalRange && globalRange.min != null && globalRange.max != null && globalRange.max !== globalRange.min) {
    rank = (value - globalRange.min) / (globalRange.max - globalRange.min)
  } else {
    if (!allValues || !allValues.length) return { cell: '', text: '' }
    const valid = allValues.filter(v => v != null)
    if (valid.length < 4) return { cell: '', text: '' }
    const sorted = [...valid].sort((a, b) => a - b)
    const n = sorted.length
    const lteCount = sorted.filter(v => v <= value).length
    rank = lteCount / n
  }

  if (config.mode === 'higher-better') {
    if (rank >= 0.75) return { cell: 'bg-green-50', text: 'text-green-700 font-semibold' }
    if (rank <= 0.25) return { cell: 'bg-red-50', text: 'text-red-600' }
  } else if (config.mode === 'lower-better') {
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
// SMILES 2D structure renderer (SvgDrawer for React compatibility)
// ---------------------------------------------------------------------------
let _smilesDrawerPromise = null
function getSmilesDrawer() {
  if (!_smilesDrawerPromise) {
    _smilesDrawerPromise = import('smiles-drawer').then(mod => mod.default || mod)
  }
  return _smilesDrawerPromise
}

// SVG cache to avoid re-parsing the same SMILES
const _svgCache = new Map()

function SmilesImage({ smiles, width = 120, height = 60 }) {
  const thumbRef = useRef(null)
  const [svgHtml, setSvgHtml] = useState(() => _svgCache.get(smiles) || null)
  const [failed, setFailed] = useState(false)
  const [popupPos, setPopupPos] = useState(null) // { x, y } or null

  useEffect(() => {
    if (!smiles) return
    if (_svgCache.has(smiles)) { setSvgHtml(_svgCache.get(smiles)); return }
    let cancelled = false
    getSmilesDrawer().then(SD => {
      if (cancelled) return
      try {
        SD.parse(smiles, (tree) => {
          if (cancelled) return
          const svgEl = document.createElementNS('http://www.w3.org/2000/svg', 'svg')
          svgEl.setAttribute('width', String(width))
          svgEl.setAttribute('height', String(height))
          document.body.appendChild(svgEl)
          const drawer = new SD.SvgDrawer({ width, height })
          drawer.draw(tree, svgEl, 'light')
          const html = svgEl.outerHTML
          document.body.removeChild(svgEl)
          _svgCache.set(smiles, html)
          if (!cancelled) setSvgHtml(html)
        }, () => {
          if (!cancelled) setFailed(true)
        })
      } catch {
        if (!cancelled) setFailed(true)
      }
    }).catch(() => { if (!cancelled) setFailed(true) })
    return () => { cancelled = true }
  }, [smiles, width, height])

  const handleEnter = useCallback(() => {
    if (!thumbRef.current) return
    const rect = thumbRef.current.getBoundingClientRect()
    const popW = 360, popH = 240
    let x = rect.left + rect.width / 2 - popW / 2
    let y = rect.top - popH - 8
    // Clamp inside viewport
    if (x < 8) x = 8
    if (x + popW > window.innerWidth - 8) x = window.innerWidth - 8 - popW
    if (y < 8) y = rect.bottom + 8 // flip below if no room above
    setPopupPos({ x, y })
  }, [])

  if (failed || !svgHtml) {
    if (failed) {
      const truncated = smiles.length > 20 ? smiles.slice(0, 20) + '…' : smiles
      return <span className="font-mono text-xs text-gray-500" title={smiles}>{truncated}</span>
    }
    return <span className="text-gray-300 text-[10px]">loading...</span>
  }

  return (
    <span
      ref={thumbRef}
      className="inline-block cursor-pointer"
      onMouseEnter={handleEnter}
      onMouseLeave={() => setPopupPos(null)}
    >
      <span
        dangerouslySetInnerHTML={{ __html: svgHtml }}
        style={{ display: 'block', width, height }}
      />
      {/* Hover zoom — fixed portal to escape overflow clipping */}
      {popupPos && createPortal(
        <div
          className="fixed z-[9999] bg-white border border-gray-200 rounded-xl shadow-2xl pointer-events-none"
          style={{ left: popupPos.x, top: popupPos.y, width: 360, padding: 12 }}
        >
          <div
            dangerouslySetInnerHTML={{ __html: svgHtml }}
            style={{ width: 336, height: 190, display: 'flex', alignItems: 'center', justifyContent: 'center' }}
            className="[&>svg]:!w-full [&>svg]:!h-full"
          />
          <div className="text-[9px] text-gray-400 font-mono mt-1 truncate" style={{ maxWidth: 336 }}>{smiles}</div>
        </div>,
        document.body
      )}
    </span>
  )
}

// ---------------------------------------------------------------------------
// Inline-editable text cell
// ---------------------------------------------------------------------------
function EditableTextCell({ value, onSave, placeholder = 'Add note…' }) {
  const [editing, setEditing] = useState(false)
  const [draft, setDraft] = useState(value || '')
  const inputRef = useRef(null)

  useEffect(() => { setDraft(value || '') }, [value])
  useEffect(() => { if (editing && inputRef.current) inputRef.current.focus() }, [editing])

  const commit = () => {
    setEditing(false)
    const trimmed = draft.trim()
    if (trimmed !== (value || '')) onSave(trimmed)
  }

  if (editing) {
    return (
      <input
        ref={inputRef}
        type="text"
        value={draft}
        onChange={e => setDraft(e.target.value)}
        onBlur={commit}
        onKeyDown={e => {
          if (e.key === 'Enter') commit()
          if (e.key === 'Escape') { setEditing(false); setDraft(value || '') }
        }}
        className="w-full text-xs px-2 py-1 border border-indigo-400 rounded-md bg-white outline-none ring-2 ring-indigo-200/60 shadow-sm"
        style={{ maxWidth: 170 }}
        onClick={e => e.stopPropagation()}
      />
    )
  }

  return (
    <span
      className={`text-xs cursor-pointer block truncate rounded px-1.5 py-0.5 transition-colors ${
        value ? 'text-gray-700 bg-slate-50 hover:bg-indigo-50' : 'text-gray-300 hover:text-indigo-400 hover:bg-indigo-50/50 italic'
      }`}
      style={{ maxWidth: 170 }}
      title={value || 'Click to add a note'}
      onClick={e => { e.stopPropagation(); setEditing(true) }}
    >
      {value || placeholder}
    </span>
  )
}

// ---------------------------------------------------------------------------
// Tags cell — renders chips with add/remove
// ---------------------------------------------------------------------------
// Stable color per tag name (hash-based, same tag = same color everywhere)
const TAG_PALETTES = [
  'bg-blue-100 text-blue-700 border-blue-200',
  'bg-purple-100 text-purple-700 border-purple-200',
  'bg-emerald-100 text-emerald-700 border-emerald-200',
  'bg-amber-100 text-amber-700 border-amber-200',
  'bg-rose-100 text-rose-700 border-rose-200',
  'bg-cyan-100 text-cyan-700 border-cyan-200',
  'bg-indigo-100 text-indigo-700 border-indigo-200',
  'bg-teal-100 text-teal-700 border-teal-200',
  'bg-orange-100 text-orange-700 border-orange-200',
  'bg-pink-100 text-pink-700 border-pink-200',
  'bg-lime-100 text-lime-700 border-lime-200',
  'bg-fuchsia-100 text-fuchsia-700 border-fuchsia-200',
]
function tagColor(tag) {
  let h = 0
  for (let i = 0; i < tag.length; i++) h = ((h << 5) - h + tag.charCodeAt(i)) | 0
  return TAG_PALETTES[Math.abs(h) % TAG_PALETTES.length]
}

function TagsCell({ tags = [], onUpdate, allKnownTags = [] }) {
  const [adding, setAdding] = useState(false)
  const [draft, setDraft] = useState('')
  const inputRef = useRef(null)
  const dropdownRef = useRef(null)

  useEffect(() => { if (adding && inputRef.current) inputRef.current.focus() }, [adding])

  // Suggestions: existing tags not already on this molecule, filtered by draft
  const suggestions = useMemo(() => {
    const tagsSet = new Set(tags)
    const available = allKnownTags.filter(t => !tagsSet.has(t))
    if (!draft.trim()) return available
    const lower = draft.toLowerCase()
    return available.filter(t => t.toLowerCase().includes(lower))
  }, [allKnownTags, tags, draft])

  const handleAdd = useCallback(() => {
    const t = draft.trim()
    if (t && !tags.includes(t)) {
      onUpdate([...tags, t])
    }
    setDraft('')
    setAdding(false)
  }, [draft, tags, onUpdate])

  const handlePickSuggestion = (tag) => {
    if (!tags.includes(tag)) {
      onUpdate([...tags, tag])
    }
    setDraft('')
    setAdding(false)
  }

  const handleRemove = (tag) => {
    onUpdate(tags.filter(t => t !== tag))
  }

  // Compute dropdown position from input ref
  const [dropdownPos, setDropdownPos] = useState(null)
  useEffect(() => {
    if (adding && inputRef.current && suggestions.length > 0) {
      const rect = inputRef.current.getBoundingClientRect()
      setDropdownPos({ top: rect.bottom + 4, left: rect.left })
    } else {
      setDropdownPos(null)
    }
  }, [adding, suggestions.length, draft])

  // Close portal dropdown on outside click
  useEffect(() => {
    if (!adding) return
    function handle(e) {
      if (dropdownRef.current?.contains(e.target)) return
      if (inputRef.current?.contains(e.target)) return
      handleAdd()
    }
    document.addEventListener('mousedown', handle)
    return () => document.removeEventListener('mousedown', handle)
  }, [adding, handleAdd])

  return (
    <span className="inline-flex flex-wrap items-center gap-1" onClick={e => e.stopPropagation()}>
      {tags.map(tag => (
        <span key={tag} className={`inline-flex items-center gap-0.5 px-2 py-0.5 rounded-full text-[10px] font-semibold border ${tagColor(tag)} shadow-sm`}>
          {tag}
          <button
            onClick={e => { e.stopPropagation(); handleRemove(tag) }}
            className="ml-0.5 opacity-40 hover:opacity-100 transition-opacity text-xs leading-none"
            title={`Remove "${tag}"`}
          >×</button>
        </span>
      ))}
      {adding ? (
        <>
          <input
            ref={inputRef}
            type="text"
            value={draft}
            onChange={e => setDraft(e.target.value)}
            onKeyDown={e => {
              if (e.key === 'Enter') handleAdd()
              if (e.key === 'Escape') { setAdding(false); setDraft('') }
            }}
            className="text-[10px] px-2 py-0.5 border border-indigo-400 rounded-full bg-white outline-none ring-1 ring-indigo-200 w-20 shadow-sm"
            placeholder="new tag…"
            onClick={e => e.stopPropagation()}
          />
          {/* Suggestions dropdown — rendered as portal to escape overflow:hidden */}
          {dropdownPos && suggestions.length > 0 && createPortal(
            <div
              ref={dropdownRef}
              className="fixed z-[9999] bg-white border border-gray-200 rounded-lg shadow-xl py-1 min-w-[130px] max-h-40 overflow-y-auto"
              style={{ top: dropdownPos.top, left: dropdownPos.left }}
              onClick={e => e.stopPropagation()}
            >
              {suggestions.slice(0, 10).map(s => (
                <button
                  key={s}
                  className="flex items-center gap-1.5 w-full text-left text-[11px] px-2.5 py-1.5 hover:bg-indigo-50 text-gray-700 truncate"
                  onMouseDown={e => { e.preventDefault(); e.stopPropagation(); handlePickSuggestion(s) }}
                >
                  <span className={`w-2.5 h-2.5 rounded-full flex-shrink-0 border ${tagColor(s)}`} />
                  {s}
                </button>
              ))}
            </div>,
            document.body
          )}
        </>
      ) : (
        <button
          onClick={e => { e.stopPropagation(); setAdding(true) }}
          className="w-5 h-5 rounded-full bg-indigo-50 hover:bg-indigo-100 border border-indigo-200 flex items-center justify-center text-indigo-400 hover:text-indigo-600 text-xs transition-colors shadow-sm"
          title="Add tag"
        >+</button>
      )}
    </span>
  )
}

// ---------------------------------------------------------------------------
// Cell value renderer
// ---------------------------------------------------------------------------
function CellValue({ col, value, colorClasses, mol, onAnnotation, allKnownTags, runNameMap }) {
  // Invalidation toggle
  if (col.type === 'invalidation') {
    const checked = !!value
    return (
      <span className="flex items-center justify-center" onClick={e => e.stopPropagation()}>
        <button
          onClick={e => { e.stopPropagation(); onAnnotation && onAnnotation(mol.id, 'invalidated', !checked) }}
          className={`w-3.5 h-3.5 rounded border flex items-center justify-center transition-all cursor-pointer ${
            checked
              ? 'bg-red-500 border-red-500 text-white'
              : 'bg-white border-gray-300 hover:border-red-300 hover:bg-red-50 text-transparent'
          }`}
          title={checked ? 'Click to mark as valid' : 'Click to mark as invalid'}
        >
          <svg className="w-2.5 h-2.5" fill="none" stroke="currentColor" viewBox="0 0 24 24" strokeWidth={3}>
            <path strokeLinecap="round" strokeLinejoin="round" d="M6 18L18 6M6 6l12 12" />
          </svg>
        </button>
      </span>
    )
  }

  // Tags
  if (col.type === 'tags') {
    return (
      <TagsCell
        tags={Array.isArray(value) ? value : []}
        onUpdate={(tags) => onAnnotation && onAnnotation(mol.id, 'tags', tags)}
        allKnownTags={allKnownTags || []}
      />
    )
  }

  // Editable text (user comments)
  if (col.type === 'editable_text') {
    return (
      <EditableTextCell
        value={value || ''}
        onSave={(text) => onAnnotation && onAnnotation(mol.id, col.key, text)}
        placeholder="Add note…"
      />
    )
  }

  // Source column — resolve run ID to label + detail tooltip
  if (col.type === 'source') {
    if (value == null) return <span className="text-gray-300 select-none">—</span>
    const entry = runNameMap?.[value]
    const label = entry?.label || String(value).slice(0, 8)
    const detail = entry?.detail || label
    return (
      <span className="text-gray-600 text-xs block truncate" title={detail} style={{ maxWidth: col.width ? `${col.width - 16}px` : 90 }}>
        {label}
      </span>
    )
  }

  if (value == null || value === undefined) {
    return <span className="text-gray-300 select-none">—</span>
  }

  if (col.type === 'boolean') {
    // For alert columns (invertBoolean), true = bad (red), false = good (green)
    const isGood = col.invertBoolean ? !value : !!value
    return isGood ? (
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
    return <SmilesImage smiles={value} width={120} height={60} />
  }

  if (col.type === 'number' && typeof value === 'number') {
    const formatted = formatNumber(col.key, value)
    return (
      <span className={`tabular-nums ${colorClasses?.text || 'text-gray-700'}`}>
        {formatted}
      </span>
    )
  }

  // Safety color code — render colored dot only (no label)
  if (col.key === 'safety_color_code') {
    const lc = String(value).toLowerCase()
    const dotColor = lc === 'green' ? 'bg-green-500' : lc === 'yellow' || lc === 'orange' ? 'bg-yellow-400' : lc === 'red' ? 'bg-red-500' : 'bg-gray-400'
    return (
      <span className="inline-flex items-center justify-center w-full">
        <span className={`w-3 h-3 rounded-full ${dotColor} flex-shrink-0`} title={lc} />
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
const ROW_HEIGHT = 44

export default function MoleculeTable({
  molecules = [],
  columns = [],
  selectedIds = new Set(),
  onToggleSelect,
  onSelectAll,
  onRowClick,
  onToggleBookmark,
  onCellPopup,
  onAnnotation,
  activeRowId = null,
  highlightedIds = null,
  totalCount,
  runs = [],
  // Pagination
  currentPage = 1,
  totalPages = 1,
  pageSize = 100,
  onPageChange,
  onPageSizeChange,
  // Server sort
  onServerSort,
  serverSortKey = null,
  serverSortDir = null,
  globalColumnRanges = {},
}) {
  const columnColorOverrides = useSettingsStore(s => s.columnColorOverrides)
  const [sorts, setSorts] = useState([]) // local secondary sort (shift-click)
  const [columnFilters, setColumnFilters] = useState({}) // {key: {type:'text',value} | {type:'number',min,max} | {type:'boolean',value}}
  const [showFilterCol, setShowFilterCol] = useState(null) // key of column with open filter
  const tableRef = useRef(null)
  const filterDropdownRef = useRef(null)

  // Collect all known tags across all molecules for tag suggestions
  const allKnownTags = useMemo(() => {
    const tagSet = new Set()
    for (const mol of molecules) {
      if (Array.isArray(mol.tags)) mol.tags.forEach(t => tagSet.add(t))
    }
    return [...tagSet].sort()
  }, [molecules])

  // Run ID → { label, detail } (e.g. label="ChEMBL", detail="Import from ChEMBL (max 50)")
  const runNameMap = useMemo(() => {
    const map = {}
    for (const run of runs) {
      let label, detail
      if (run.type === 'import') {
        const cfg = run.config || {}
        const src = cfg.source || cfg.database || ''
        const dbNames = cfg.databases || (src ? [src] : [])
        if (dbNames.length > 0) {
          label = dbNames.map(d => d.charAt(0).toUpperCase() + d.slice(1)).join('+')
          const parts = [label]
          if (cfg.max_compounds || cfg.max_per_source) parts.push(`max ${cfg.max_compounds || cfg.max_per_source}`)
          if (cfg.filters) {
            const f = cfg.filters
            if (f.mwt_min || f.mwt_max) parts.push(`MW ${f.mwt_min || ''}–${f.mwt_max || ''}`)
            if (f.logp_min || f.logp_max) parts.push(`LogP ${f.logp_min || ''}–${f.logp_max || ''}`)
          }
          detail = parts.join(' · ')
        } else if (cfg.source === 'phase_selection') {
          label = 'Phase import'
          detail = 'Bookmarked from previous phase'
        } else {
          label = 'Import'
          detail = cfg.original_filename || 'Import'
        }
      } else if (run.type === 'calculation' && run.calculation_types?.length) {
        label = run.calculation_types.map(t => t.charAt(0).toUpperCase() + t.slice(1)).join(' + ')
        detail = label
      } else {
        label = run.type ? run.type.charAt(0).toUpperCase() + run.type.slice(1) : run.id
        detail = label
      }
      map[run.id] = { label, detail }
    }
    return map
  }, [runs])

  // --- Resizable columns ---
  const STORAGE_KEY = 'bx-col-widths'
  const [colWidthOverrides, setColWidthOverrides] = useState(() => {
    try { return JSON.parse(localStorage.getItem(STORAGE_KEY)) || {} } catch { return {} }
  })
  const resizeRef = useRef(null) // { colKey, startX, startWidth }

  const handleResizeStart = useCallback((e, colKey, currentWidth) => {
    e.preventDefault()
    e.stopPropagation()
    resizeRef.current = { colKey, startX: e.clientX, startWidth: currentWidth }
    const onMove = (ev) => {
      if (!resizeRef.current) return
      const delta = ev.clientX - resizeRef.current.startX
      const newWidth = Math.max(40, resizeRef.current.startWidth + delta)
      setColWidthOverrides(prev => {
        const next = { ...prev, [resizeRef.current.colKey]: newWidth }
        localStorage.setItem(STORAGE_KEY, JSON.stringify(next))
        return next
      })
    }
    const onUp = () => {
      resizeRef.current = null
      document.removeEventListener('mousemove', onMove)
      document.removeEventListener('mouseup', onUp)
      document.body.style.cursor = ''
      document.body.style.userSelect = ''
    }
    document.addEventListener('mousemove', onMove)
    document.addEventListener('mouseup', onUp)
    document.body.style.cursor = 'col-resize'
    document.body.style.userSelect = 'none'
  }, [])

  const handleResizeReset = useCallback((colKey) => {
    setColWidthOverrides(prev => {
      const next = { ...prev }
      delete next[colKey]
      localStorage.setItem(STORAGE_KEY, JSON.stringify(next))
      return next
    })
  }, [])

  // Close column filter dropdown on click outside
  useEffect(() => {
    if (!showFilterCol) return
    const handler = (e) => {
      if (filterDropdownRef.current && !filterDropdownRef.current.contains(e.target)) {
        setShowFilterCol(null)
      }
    }
    document.addEventListener('mousedown', handler)
    return () => document.removeEventListener('mousedown', handler)
  }, [showFilterCol])

  // Effective color config per column (merges user overrides with defaults)
  const colorConfigs = useMemo(() => {
    const map = {}
    columns.forEach(col => {
      if (col.type === 'number') {
        map[col.key] = getEffectiveColorScale(col.key, columnColorOverrides)
      }
    })
    return map
  }, [columns, columnColorOverrides])

  // Pre-compute column value arrays for color scaling (all numeric columns that have a color mode)
  const colValues = useMemo(() => {
    const map = {}
    columns.forEach(col => {
      const cfg = colorConfigs[col.key]
      if (col.type === 'number' && cfg && cfg.mode !== 'none') {
        map[col.key] = molecules
          .map(m => m[col.key])
          .filter(v => v != null && typeof v === 'number')
      }
    })
    return map
  }, [molecules, columns, colorConfigs])

  // Column-filtered molecules
  const colFiltered = useMemo(() => {
    const filterKeys = Object.keys(columnFilters).filter(k => {
      const f = columnFilters[k]
      if (f == null) return false
      // Backward compat: plain string (old format)
      if (typeof f === 'string') return f !== ''
      // Typed filter objects
      if (f.type === 'text') return f.value && f.value !== ''
      if (f.type === 'number') return f.min != null || f.max != null
      if (f.type === 'boolean') return f.value !== null && f.value !== undefined
      if (f.type === 'enum') return f.value != null
      return false
    })
    if (!filterKeys.length) return molecules
    return molecules.filter(mol => {
      return filterKeys.every(key => {
        const val = mol[key]
        const f = columnFilters[key]
        // Backward compat: plain string (old format) → text search
        if (typeof f === 'string') {
          if (val == null) return false
          const filterVal = f.toLowerCase()
          if (typeof val === 'number') return String(val).includes(filterVal)
          if (typeof val === 'boolean') return (val ? 'yes' : 'no').includes(filterVal)
          return String(val).toLowerCase().includes(filterVal)
        }
        // Typed filters
        if (f.type === 'enum') {
          if (val == null) return false
          return String(val).toLowerCase() === f.value.toLowerCase()
        }
        if (f.type === 'text') {
          if (val == null) return false
          // Tags: search within array items
          if (Array.isArray(val)) {
            return val.some(t => String(t).toLowerCase().includes(f.value.toLowerCase()))
          }
          return String(val).toLowerCase().includes(f.value.toLowerCase())
        }
        if (f.type === 'number') {
          if (val == null || typeof val !== 'number') return false
          if (f.min != null && val < f.min) return false
          if (f.max != null && val > f.max) return false
          return true
        }
        if (f.type === 'boolean') {
          if (f.value === null) return true // "All" selected
          return val === f.value
        }
        return true
      })
    })
  }, [molecules, columnFilters])

  // Sorted molecules
  const sorted = useMemo(() => {
    if (!sorts.length) return colFiltered
    return [...colFiltered].sort((a, b) => {
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
  }, [colFiltered, sorts])

  // Virtual scrolling
  const rowVirtualizer = useVirtualizer({
    count: sorted.length,
    getScrollElement: () => tableRef.current,
    estimateSize: () => ROW_HEIGHT,
    overscan: 15,
  })

  // Scroll to active row when activeRowId changes (e.g. prev/next in detail panel)
  useEffect(() => {
    if (!activeRowId) return
    const idx = sorted.findIndex(m => m.id === activeRowId)
    if (idx >= 0) rowVirtualizer.scrollToIndex(idx, { align: 'auto' })
  }, [activeRowId]) // eslint-disable-line react-hooks/exhaustive-deps

  // Compute effective column width: override > max(declared width, label-based minimum)
  const effectiveWidths = useMemo(() => columns.map(col => {
    if (colWidthOverrides[col.key]) return Math.max(40, colWidthOverrides[col.key])
    const labelChars = (col.label || '').length
    const unitChars = col.unit ? col.unit.length + 2 : 0
    const textWidth = (labelChars + unitChars) * 7.5
    const extras = 24
    const padding = 20
    const minFromLabel = Math.ceil(textWidth + extras + padding)
    return Math.max(col.width || 80, minFromLabel)
  }), [columns, colWidthOverrides])

  // Grid template: checkbox(40px) + bookmark(32px) + data columns
  const gridTemplate = useMemo(() => {
    const colWidths = effectiveWidths.map(w => `${w}px`)
    return `40px 32px ${colWidths.join(' ')}`
  }, [effectiveWidths])

  const minTableWidth = useMemo(
    () => Math.max(600, effectiveWidths.reduce((sum, w) => sum + w, 0) + 90),
    [effectiveWidths]
  )

  // No custom label mapping — GROUP_META is the single source of truth

  // Group runs: merge consecutive columns with the same group for the group header row
  const groupRuns = useMemo(() => {
    const result = []
    for (const col of columns) {
      const g = col.group || 'molecule'
      if (result.length > 0 && result[result.length - 1].group === g) {
        result[result.length - 1].count++
      } else {
        result.push({ group: g, count: 1 })
      }
    }
    return result
  }, [columns])

  // Header click: primary sort goes to server, shift-click adds local secondary sort
  const handleHeaderClick = useCallback((col, e) => {
    if (!col.sortable) return
    if (e.shiftKey) {
      // Shift-click: local secondary sort on top of server sort
      setSorts(prev => {
        const existing = prev.find(s => s.key === col.key)
        if (existing) {
          if (existing.dir === 'asc') return prev.map(s => s.key === col.key ? { ...s, dir: 'desc' } : s)
          return prev.filter(s => s.key !== col.key)
        }
        return [...prev, { key: col.key, dir: 'asc' }].slice(-2)
      })
    } else {
      // Normal click: server-side sort
      setSorts([]) // clear local secondary sorts
      if (onServerSort) {
        const effectiveKey = col.sortKey || col.key
        const currentDir = serverSortKey === effectiveKey ? serverSortDir : null
        if (!currentDir) onServerSort(effectiveKey, 'asc')
        else if (currentDir === 'asc') onServerSort(effectiveKey, 'desc')
        else onServerSort('created_at', 'desc') // third click resets to default
      }
    }
  }, [onServerSort, serverSortKey, serverSortDir])

  function getSortDir(key, sortKey) {
    // Check server sort first (use sortKey if available), then local secondary sorts
    const effectiveKey = sortKey || key
    if (serverSortKey === effectiveKey) return serverSortDir
    const s = sorts.find(s => s.key === key)
    return s ? s.dir : null
  }

  function getSortRank(key, sortKey) {
    // Only show rank when there are multiple sort levels (server + local)
    const effectiveKey = sortKey || key
    const totalSorts = (serverSortKey ? 1 : 0) + sorts.length
    if (totalSorts < 2) return null
    if (serverSortKey === effectiveKey) return '1'
    const idx = sorts.findIndex(s => s.key === key)
    return idx >= 0 ? (idx + (serverSortKey ? 2 : 1)).toString() : null
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

  // Effective total for footer display
  const displayTotal = totalCount != null ? totalCount : molecules.length

  return (
    <div className="flex flex-col card overflow-hidden">
      {/* Scrollable virtual table */}
      <div
        ref={tableRef}
        className="overflow-auto"
        style={{
          maxHeight: 'calc(100vh - 170px)',
          minHeight: 200,
          scrollbarWidth: 'thin',
          scrollbarColor: '#e5e7eb transparent',
        }}
      >
        <div style={{ minWidth: minTableWidth }}>
          {/* Group header row — sticky */}
          <div className="sticky top-0 z-10 bg-gray-50">
          <div
            style={{ display: 'grid', gridTemplateColumns: gridTemplate }}
          >
            {/* Empty cells for checkbox + bookmark */}
            <div className="col-span-2" />
            {/* Group spans */}
            {groupRuns.map((run, i) => {
              const meta = GROUP_META[run.group] || GROUP_META.molecule
              const label = meta.label
              return (
                <div
                  key={`${run.group}-${i}`}
                  className={`${meta.bg} ${meta.text} ${meta.border} border-b-2 text-[10px] uppercase tracking-wider font-semibold text-center py-1 whitespace-nowrap`}
                  style={{ gridColumn: `span ${run.count}` }}
                >
                  {label}
                </div>
              )
            })}
          </div>
          {/* Column header row */}
          <div
            className="bg-gray-50 border-b border-gray-200"
            style={{ display: 'grid', gridTemplateColumns: gridTemplate }}
          >
            {/* Checkbox column */}
            <div className="px-2 py-2.5 flex items-center">
              <input
                type="checkbox"
                checked={allVisible}
                ref={el => { if (el) el.indeterminate = someVisible }}
                onChange={onSelectAll}
                className="accent-bx-mint cursor-pointer w-3.5 h-3.5"
                aria-label="Select all"
              />
            </div>

            {/* Bookmark column header */}
            <div className="px-1 py-2.5 flex items-center justify-center">
              <button
                onClick={() => {
                  if (!onToggleBookmark) return
                  const targets = selectedIds.size > 0
                    ? molecules.filter(m => selectedIds.has(m.id))
                    : molecules
                  const allBookmarked = targets.every(m => m.bookmarked)
                  targets.forEach(m => {
                    if (allBookmarked || !m.bookmarked) onToggleBookmark(m.id)
                  })
                }}
                disabled={!onToggleBookmark}
                className={`p-0.5 rounded transition-colors ${onToggleBookmark ? 'hover:bg-yellow-50 cursor-pointer' : 'cursor-default opacity-50'}`}
                title={selectedIds.size > 0 ? `Toggle bookmark on ${selectedIds.size} selected` : 'Toggle bookmark on all'}
              >
                <svg className="w-3.5 h-3.5 text-gray-300" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                    d="M5 5a2 2 0 012-2h10a2 2 0 012 2v16l-7-3.5L5 21V5z" />
                </svg>
              </button>
            </div>

            {/* Data column headers */}
            {columns.map(col => {
              const dir = getSortDir(col.key, col.sortKey)
              const rank = getSortRank(col.key, col.sortKey)
              const isActive = dir !== null
              const hasFilter = (() => {
                const f = columnFilters[col.key]
                if (f == null) return false
                if (typeof f === 'string') return f !== ''
                if (f.type === 'text') return f.value && f.value !== ''
                if (f.type === 'number') return f.min != null || f.max != null
                if (f.type === 'boolean') return f.value !== null && f.value !== undefined
                if (f.type === 'enum') return f.value != null
                return false
              })()
              const isFilterable = col.type === 'number' || col.type === 'text' || col.type === 'boolean' || col.type === 'smiles' || col.type === 'tags' || col.type === 'editable_text' || col.type === 'invalidation' || col.type === 'source'

              return (
                <div
                  key={col.key}
                  className={`px-2 py-1.5 font-semibold uppercase tracking-wide text-xs select-none relative flex flex-col gap-0.5 items-center ${
                    col.sortable ? 'cursor-pointer transition-colors' : ''
                  } ${
                    isActive ? 'text-bx-light-text bg-blue-50/50' : 'text-gray-500 hover:text-bx-light-text hover:bg-gray-100'
                  }`}
                  onClick={e => handleHeaderClick(col, e)}
                >
                  {/* Line 1: label + unit */}
                  <span className="flex items-center gap-1 whitespace-nowrap leading-tight">
                    {col.label}
                    {col.unit && (
                      <span className="font-normal text-gray-300 text-[9px]">({col.unit})</span>
                    )}
                  </span>
                  {/* Line 2: sort + info tip + filter — all centered */}
                  {(col.sortable || TIPS[col.key] || isFilterable) && (
                    <span className="flex items-center gap-1">
                      {col.sortable && <SortIcon direction={dir} rank={rank} />}
                      {TIPS[col.key] && <InfoTip text={TIPS[col.key]} size="xs" />}
                      {isFilterable && (
                        <button
                          onClick={e => { e.stopPropagation(); setShowFilterCol(prev => prev === col.key ? null : col.key) }}
                          className={`p-0.5 rounded transition-colors ${hasFilter ? 'text-blue-500' : 'text-gray-300 hover:text-gray-500'}`}
                          title="Filter column"
                        >
                          <svg className="w-2.5 h-2.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M3 4a1 1 0 011-1h16a1 1 0 011 1v2.586a1 1 0 01-.293.707l-6.414 6.414a1 1 0 00-.293.707V17l-4 4v-6.586a1 1 0 00-.293-.707L3.293 7.293A1 1 0 013 6.586V4z" />
                          </svg>
                        </button>
                      )}
                    </span>
                  )}
                  {/* Column filter dropdown */}
                  {showFilterCol === col.key && (
                    <div
                      ref={filterDropdownRef}
                      className="absolute top-full left-0 z-20 mt-1 bg-white border border-gray-200 rounded-lg shadow-lg p-1.5"
                      onClick={e => e.stopPropagation()}
                    >
                      {/* Number filter: min/max inputs */}
                      {col.type === 'number' && (() => {
                        const f = columnFilters[col.key]
                        const minVal = (f && typeof f === 'object' && f.type === 'number') ? (f.min ?? '') : ''
                        const maxVal = (f && typeof f === 'object' && f.type === 'number') ? (f.max ?? '') : ''
                        return (
                          <div className="flex items-center gap-1">
                            <input
                              type="number"
                              autoFocus
                              value={minVal}
                              onChange={e => {
                                const v = e.target.value === '' ? null : Number(e.target.value)
                                setColumnFilters(prev => ({
                                  ...prev,
                                  [col.key]: { type: 'number', min: v, max: (prev[col.key]?.type === 'number' ? prev[col.key]?.max : null) ?? null }
                                }))
                              }}
                              onKeyDown={e => { if (e.key === 'Escape') setShowFilterCol(null) }}
                              placeholder="Min"
                              className="text-xs px-1.5 py-1 border border-gray-200 rounded w-16 outline-none focus:border-blue-300 tabular-nums"
                            />
                            <span className="text-gray-300 text-xs">–</span>
                            <input
                              type="number"
                              value={maxVal}
                              onChange={e => {
                                const v = e.target.value === '' ? null : Number(e.target.value)
                                setColumnFilters(prev => ({
                                  ...prev,
                                  [col.key]: { type: 'number', min: (prev[col.key]?.type === 'number' ? prev[col.key]?.min : null) ?? null, max: v }
                                }))
                              }}
                              onKeyDown={e => { if (e.key === 'Escape') setShowFilterCol(null) }}
                              placeholder="Max"
                              className="text-xs px-1.5 py-1 border border-gray-200 rounded w-16 outline-none focus:border-blue-300 tabular-nums"
                            />
                          </div>
                        )
                      })()}

                      {/* Boolean filter: All / Yes / No buttons */}
                      {col.type === 'boolean' && (() => {
                        const f = columnFilters[col.key]
                        const currentVal = (f && typeof f === 'object' && f.type === 'boolean') ? f.value : null
                        const btnBase = 'text-xs px-2 py-1 rounded transition-colors'
                        const btnActive = 'bg-blue-100 text-blue-700 font-medium'
                        const btnInactive = 'text-gray-500 hover:bg-gray-100'
                        return (
                          <div className="flex items-center gap-0.5">
                            <button
                              className={`${btnBase} ${currentVal === null || currentVal === undefined ? btnActive : btnInactive}`}
                              onClick={() => setColumnFilters(prev => { const n = { ...prev }; delete n[col.key]; return n })}
                            >All</button>
                            <button
                              className={`${btnBase} ${currentVal === true ? btnActive : btnInactive}`}
                              onClick={() => setColumnFilters(prev => ({ ...prev, [col.key]: { type: 'boolean', value: true } }))}
                            >Yes</button>
                            <button
                              className={`${btnBase} ${currentVal === false ? btnActive : btnInactive}`}
                              onClick={() => setColumnFilters(prev => ({ ...prev, [col.key]: { type: 'boolean', value: false } }))}
                            >No</button>
                          </div>
                        )
                      })()}

                      {/* Safety color code filter: colored buttons */}
                      {col.key === 'safety_color_code' && (() => {
                        const f = columnFilters[col.key]
                        const currentVal = (f && typeof f === 'object' && f.type === 'enum') ? f.value : null
                        const btnBase = 'text-xs px-2 py-1 rounded transition-colors inline-flex items-center gap-1'
                        const btnActive = 'ring-2 ring-blue-400 font-medium'
                        const btnInactive = 'opacity-60 hover:opacity-100'
                        const colors = [
                          { value: 'green', dot: 'bg-green-500', label: '' },
                          { value: 'yellow', dot: 'bg-yellow-400', label: '' },
                          { value: 'red', dot: 'bg-red-500', label: '' },
                        ]
                        return (
                          <div className="flex items-center gap-1">
                            <button
                              className={`${btnBase} ${currentVal === null ? 'bg-blue-100 text-blue-700 font-medium' : 'text-gray-500 hover:bg-gray-100'}`}
                              onClick={() => setColumnFilters(prev => { const n = { ...prev }; delete n[col.key]; return n })}
                            >All</button>
                            {colors.map(c => (
                              <button
                                key={c.value}
                                className={`${btnBase} bg-gray-50 ${currentVal === c.value ? btnActive : btnInactive}`}
                                onClick={() => setColumnFilters(prev => ({ ...prev, [col.key]: { type: 'enum', value: c.value } }))}
                                title={c.value}
                              >
                                <span className={`w-3 h-3 rounded-full ${c.dot}`} />
                              </button>
                            ))}
                          </div>
                        )
                      })()}

                      {/* Tags filter: text search on tag values */}
                      {col.type === 'tags' && (() => {
                        const f = columnFilters[col.key]
                        const textVal = (f && f.type === 'text') ? f.value : ''
                        return (
                          <input
                            type="text"
                            autoFocus
                            value={textVal}
                            onChange={e => setColumnFilters(prev => ({ ...prev, [col.key]: { type: 'text', value: e.target.value } }))}
                            onKeyDown={e => { if (e.key === 'Escape') setShowFilterCol(null) }}
                            placeholder="Filter by tag…"
                            className="text-xs px-2 py-1 border border-gray-200 rounded w-28 outline-none focus:border-blue-300"
                          />
                        )
                      })()}

                      {/* Editable text filter (notes) */}
                      {col.type === 'editable_text' && (() => {
                        const f = columnFilters[col.key]
                        const textVal = (f && f.type === 'text') ? f.value : ''
                        return (
                          <input
                            type="text"
                            autoFocus
                            value={textVal}
                            onChange={e => setColumnFilters(prev => ({ ...prev, [col.key]: { type: 'text', value: e.target.value } }))}
                            onKeyDown={e => { if (e.key === 'Escape') setShowFilterCol(null) }}
                            placeholder={`Filter ${col.label}…`}
                            className="text-xs px-2 py-1 border border-gray-200 rounded w-28 outline-none focus:border-blue-300"
                          />
                        )
                      })()}

                      {/* Invalidation filter: All / Yes / No */}
                      {col.type === 'invalidation' && (() => {
                        const f = columnFilters[col.key]
                        const currentVal = (f && typeof f === 'object' && f.type === 'boolean') ? f.value : null
                        const btnBase = 'text-xs px-2 py-1 rounded transition-colors'
                        const btnActive = 'bg-blue-100 text-blue-700 font-medium'
                        const btnInactive = 'text-gray-500 hover:bg-gray-100'
                        return (
                          <div className="flex items-center gap-0.5">
                            <button
                              className={`${btnBase} ${currentVal === null || currentVal === undefined ? btnActive : btnInactive}`}
                              onClick={() => setColumnFilters(prev => { const n = { ...prev }; delete n[col.key]; return n })}
                            >All</button>
                            <button
                              className={`${btnBase} ${currentVal === true ? btnActive : btnInactive}`}
                              onClick={() => setColumnFilters(prev => ({ ...prev, [col.key]: { type: 'boolean', value: true } }))}
                            >Yes</button>
                            <button
                              className={`${btnBase} ${currentVal === false ? btnActive : btnInactive}`}
                              onClick={() => setColumnFilters(prev => ({ ...prev, [col.key]: { type: 'boolean', value: false } }))}
                            >No</button>
                          </div>
                        )
                      })()}

                      {/* Text / SMILES filter: text input (skip safety_color_code) */}
                      {(col.type === 'text' || col.type === 'smiles' || col.type === 'source') && col.key !== 'safety_color_code' && (() => {
                        const f = columnFilters[col.key]
                        // Backward compat: plain string → extract value
                        const textVal = typeof f === 'string' ? f : (f && f.type === 'text' ? f.value : '')
                        return (
                          <input
                            type="text"
                            autoFocus
                            value={textVal}
                            onChange={e => setColumnFilters(prev => ({ ...prev, [col.key]: { type: 'text', value: e.target.value } }))}
                            onKeyDown={e => { if (e.key === 'Escape') setShowFilterCol(null) }}
                            placeholder={`Filter ${col.label}…`}
                            className="text-xs px-2 py-1 border border-gray-200 rounded w-28 outline-none focus:border-blue-300"
                          />
                        )
                      })()}

                      {hasFilter && (
                        <button
                          onClick={() => setColumnFilters(prev => { const n = { ...prev }; delete n[col.key]; return n })}
                          className="text-xs text-red-400 hover:text-red-600 mt-0.5 block"
                        >Clear</button>
                      )}
                    </div>
                  )}
                  {/* Resize handle */}
                  <div
                    className="absolute top-0 right-0 w-1 h-full cursor-col-resize hover:bg-blue-400/50 active:bg-blue-500/70 transition-colors z-10"
                    onMouseDown={e => handleResizeStart(e, col.key, effectiveWidths[columns.indexOf(col)])}
                    onDoubleClick={e => { e.stopPropagation(); handleResizeReset(col.key) }}
                    title="Drag to resize · Double-click to reset"
                  />
                </div>
              )
            })}
          </div>
          </div>{/* end sticky group+header wrapper */}

          {/* Virtualized body */}
          <div
            style={{
              height: `${rowVirtualizer.getTotalSize()}px`,
              width: '100%',
              position: 'relative',
            }}
          >
            {rowVirtualizer.getVirtualItems().map(virtualRow => {
              const mol = sorted[virtualRow.index]
              const idx = virtualRow.index
              const isSelected = selectedIds.has(mol.id)
              const isActive = activeRowId === mol.id
              const isHighlighted = highlightedIds ? highlightedIds.has(mol.id) : false
              const isBookmarked = mol.bookmarked

              let rowBg = idx % 2 === 0 ? 'bg-white' : 'bg-slate-50/40'
              if (isSelected) rowBg = 'bg-blue-50'
              if (isHighlighted) rowBg = 'bg-emerald-50/70'
              if (isActive) rowBg = 'bg-blue-100/80'

              const isInvalidated = !!mol.invalidated

              return (
                <div
                  key={mol.id}
                  className={`cursor-pointer transition-colors duration-100 group ${rowBg} hover:bg-blue-50/60 border-b border-gray-50 ${
                    isInvalidated ? 'border-l-[3px] border-l-red-300' :
                    isSelected ? 'border-l-[3px] border-l-bx-surface' :
                    isBookmarked ? 'border-l-2 border-l-yellow-300' :
                    'border-l-2 border-l-transparent'
                  }`}
                  style={{
                    display: 'grid',
                    gridTemplateColumns: gridTemplate,
                    position: 'absolute',
                    top: 0,
                    left: 0,
                    width: '100%',
                    height: `${virtualRow.size}px`,
                    transform: `translateY(${virtualRow.start}px)`,
                  }}
                  onClick={(e) => onRowClick && onRowClick(mol, e)}
                >
                  {/* Checkbox */}
                  <div
                    className="px-2 flex items-center"
                    onClick={e => e.stopPropagation()}
                  >
                    <input
                      type="checkbox"
                      checked={isSelected}
                      onChange={() => onToggleSelect && onToggleSelect(mol.id)}
                      className="accent-bx-mint cursor-pointer w-3.5 h-3.5"
                      aria-label={`Select ${mol.name || mol.id}`}
                    />
                  </div>

                  {/* Bookmark star */}
                  <div
                    className="px-1 flex items-center justify-center"
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
                  </div>

                  {/* Data cells */}
                  {columns.map(col => {
                    const value = mol[col.key]
                    const colorClasses = col.type === 'number' && colorConfigs[col.key]
                      ? getValueColorClasses(value, colValues[col.key] || [], colorConfigs[col.key], globalColumnRanges[col.key])
                      : { cell: '', text: '' }
                    const hasPopup = col.popup && onCellPopup && value != null

                    return (
                      <div
                        key={col.key}
                        className={`px-3 flex items-center max-h-[${ROW_HEIGHT}px] ${col.type === 'tags' ? 'overflow-visible' : 'overflow-hidden'} ${col.type === 'number' ? 'justify-end tabular-nums' : 'justify-start'} ${colorClasses.cell} ${hasPopup ? 'cursor-pointer hover:bg-blue-50/80 transition-colors' : ''} ${mol.invalidated && col.key !== 'invalidated' ? 'opacity-40' : ''}`}
                        onClick={hasPopup ? (e) => { e.stopPropagation(); onCellPopup(col.popup, mol) } : undefined}
                        title={hasPopup ? `Click for ${col.label} details` : undefined}
                      >
                        <span className={`text-sm ${hasPopup ? 'underline decoration-dotted underline-offset-2 decoration-gray-300' : ''}`}>
                          <CellValue col={col} value={value} colorClasses={colorClasses} mol={mol} onAnnotation={onAnnotation} allKnownTags={allKnownTags} runNameMap={runNameMap} />
                        </span>
                      </div>
                    )
                  })}
                </div>
              )
            })}
          </div>

        </div>
      </div>

      {/* Footer — pagination controls */}
      <div className="px-4 py-2 bg-gray-50 border-t border-gray-100 flex items-center justify-between gap-4">
        {/* Left: status text */}
        <span className="text-xs text-gray-400 tabular-nums whitespace-nowrap">
          Page <strong className="text-gray-600">{currentPage}</strong> of <strong className="text-gray-600">{totalPages}</strong>
          {' '}&mdash; {(displayTotal || 0).toLocaleString()} molecules
          {Object.keys(columnFilters).filter(k => {
            const f = columnFilters[k]
            if (f == null) return false
            if (typeof f === 'string') return f !== ''
            if (f.type === 'text') return f.value && f.value !== ''
            if (f.type === 'number') return f.min != null || f.max != null
            if (f.type === 'boolean') return f.value !== null && f.value !== undefined
            return false
          }).length > 0 && (
            <button
              onClick={() => setColumnFilters({})}
              className="ml-2 text-xs text-blue-500 hover:text-blue-700 underline"
            >Clear column filters</button>
          )}
          {(serverSortKey || sorts.length > 0) && (
            <span className="ml-2 text-gray-400">
              Sorted by{' '}
              {serverSortKey && serverSortKey !== 'created_at' && (
                <span>
                  <span className="font-medium text-bx-light-text">{columns.find(c => c.key === serverSortKey)?.label || serverSortKey}</span>
                  <span className="text-[10px]">{serverSortDir === 'asc' ? ' ↑' : ' ↓'}</span>
                </span>
              )}
              {sorts.map((s, i) => (
                <span key={s.key}>
                  {(i > 0 || (serverSortKey && serverSortKey !== 'created_at')) && ' · '}
                  <span className="font-medium text-bx-light-text">{columns.find(c => c.key === s.key)?.label || s.key}</span>
                  <span className="text-[10px]">{s.dir === 'asc' ? ' ↑' : ' ↓'}</span>
                </span>
              ))}
              <button onClick={() => { setSorts([]); onServerSort?.('created_at', 'desc') }} className="ml-1 text-gray-400 hover:text-red-500 text-xs underline underline-offset-2">Clear</button>
            </span>
          )}
        </span>

        {/* Center: page navigation */}
        {totalPages > 1 && onPageChange && (
          <div className="flex items-center gap-1">
            <button
              onClick={() => onPageChange(currentPage - 1)}
              disabled={currentPage <= 1}
              className="px-2 py-1 text-xs font-medium rounded-md border border-gray-200 text-gray-600 hover:bg-gray-100 disabled:opacity-30 disabled:cursor-not-allowed transition-colors"
            >
              Prev
            </button>
            {(() => {
              const pages = []
              const addPage = (p) => pages.push({ type: 'page', value: p })
              const addEllipsis = (key) => pages.push({ type: 'ellipsis', key })

              addPage(1)
              if (currentPage > 3) addEllipsis('start')
              for (let p = Math.max(2, currentPage - 1); p <= Math.min(totalPages - 1, currentPage + 1); p++) {
                addPage(p)
              }
              if (currentPage < totalPages - 2) addEllipsis('end')
              if (totalPages > 1) addPage(totalPages)

              // Deduplicate
              const seen = new Set()
              return pages.filter(item => {
                const key = item.type === 'page' ? `p${item.value}` : item.key
                if (seen.has(key)) return false
                seen.add(key)
                return true
              }).map(item => {
                if (item.type === 'ellipsis') {
                  return <span key={item.key} className="px-1 text-xs text-gray-400 select-none">&hellip;</span>
                }
                const isActive = item.value === currentPage
                return (
                  <button
                    key={item.value}
                    onClick={() => onPageChange(item.value)}
                    className={`min-w-[28px] px-1.5 py-1 text-xs font-medium rounded-md border transition-colors ${
                      isActive
                        ? 'bg-bx-surface text-white border-bx-surface shadow-sm'
                        : 'border-gray-200 text-gray-600 hover:bg-gray-100'
                    }`}
                  >
                    {item.value}
                  </button>
                )
              })
            })()}
            <button
              onClick={() => onPageChange(currentPage + 1)}
              disabled={currentPage >= totalPages}
              className="px-2 py-1 text-xs font-medium rounded-md border border-gray-200 text-gray-600 hover:bg-gray-100 disabled:opacity-30 disabled:cursor-not-allowed transition-colors"
            >
              Next
            </button>
          </div>
        )}

        {/* Right: page size selector + selected count */}
        <div className="flex items-center gap-3 flex-shrink-0">
          {onPageSizeChange && (
            <div className="flex items-center gap-1.5">
              <span className="text-xs text-gray-400">Show</span>
              {[50, 100, 200].map(size => (
                <button
                  key={size}
                  onClick={() => onPageSizeChange(size)}
                  className={`px-2 py-0.5 text-xs font-medium rounded-md border transition-colors ${
                    pageSize === size
                      ? 'bg-bx-surface text-white border-bx-surface'
                      : 'border-gray-200 text-gray-500 hover:bg-gray-100'
                  }`}
                >
                  {size}
                </button>
              ))}
            </div>
          )}
          {selectedIds.size > 0 && (
            <span className="text-xs text-bx-light-text font-medium tabular-nums">
              {selectedIds.size} selected
            </span>
          )}
        </div>
      </div>
    </div>
  )
}
