import React, { useState, useMemo, useRef, useCallback, useEffect } from 'react'
import InfoTip, { TIPS } from './InfoTip.jsx'

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
  const containerRef = useRef(null)
  const [svgHtml, setSvgHtml] = useState(() => _svgCache.get(smiles) || null)
  const [failed, setFailed] = useState(false)
  const [hovered, setHovered] = useState(false)

  useEffect(() => {
    if (!smiles) return
    if (_svgCache.has(smiles)) { setSvgHtml(_svgCache.get(smiles)); return }
    let cancelled = false
    getSmilesDrawer().then(SD => {
      if (cancelled) return
      try {
        SD.parse(smiles, (tree) => {
          if (cancelled) return
          // Create a temporary SVG element
          const svgEl = document.createElementNS('http://www.w3.org/2000/svg', 'svg')
          svgEl.setAttribute('width', String(width))
          svgEl.setAttribute('height', String(height))
          document.body.appendChild(svgEl) // must be in DOM briefly for SvgDrawer
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

  if (failed || !svgHtml) {
    if (failed) {
      const truncated = smiles.length > 20 ? smiles.slice(0, 20) + '…' : smiles
      return <span className="font-mono text-[11px] text-gray-500" title={smiles}>{truncated}</span>
    }
    return <span className="text-gray-300 text-[10px]">loading...</span>
  }

  return (
    <span
      className="relative inline-block cursor-pointer"
      onMouseEnter={() => setHovered(true)}
      onMouseLeave={() => setHovered(false)}
    >
      <span
        ref={containerRef}
        dangerouslySetInnerHTML={{ __html: svgHtml }}
        style={{ display: 'block', width, height }}
        title={smiles}
      />
      {/* Hover zoom */}
      {hovered && (
        <span
          className="absolute z-50 bg-white border border-gray-200 rounded-xl shadow-2xl p-2 pointer-events-none"
          style={{ bottom: '110%', left: '50%', transform: 'translateX(-50%)', width: 240, height: 140 }}
          dangerouslySetInnerHTML={{
            __html: svgHtml.replace(`width="${width}"`, 'width="224"').replace(`height="${height}"`, 'height="124"')
          }}
        />
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

  useEffect(() => { if (editing && inputRef.current) inputRef.current.focus() }, [editing])

  if (editing) {
    return (
      <input
        ref={inputRef}
        type="text"
        value={draft}
        onChange={e => setDraft(e.target.value)}
        onBlur={() => { setEditing(false); if (draft !== (value || '')) onSave(draft) }}
        onKeyDown={e => {
          if (e.key === 'Enter') { setEditing(false); if (draft !== (value || '')) onSave(draft) }
          if (e.key === 'Escape') { setEditing(false); setDraft(value || '') }
        }}
        className="w-full text-xs px-1.5 py-0.5 border border-blue-300 rounded bg-white outline-none ring-1 ring-blue-200"
        style={{ maxWidth: 140 }}
      />
    )
  }

  return (
    <span
      className="text-xs text-gray-500 cursor-pointer hover:text-gray-700 block truncate"
      style={{ maxWidth: 140 }}
      title={value || 'Click to add'}
      onClick={e => { e.stopPropagation(); setEditing(true) }}
    >
      {value || <span className="text-gray-300 italic">{placeholder}</span>}
    </span>
  )
}

// ---------------------------------------------------------------------------
// Tags cell — renders chips with add/remove
// ---------------------------------------------------------------------------
function TagsCell({ tags = [], onUpdate }) {
  const [adding, setAdding] = useState(false)
  const [draft, setDraft] = useState('')
  const inputRef = useRef(null)

  useEffect(() => { if (adding && inputRef.current) inputRef.current.focus() }, [adding])

  const colors = ['bg-blue-100 text-blue-700', 'bg-purple-100 text-purple-700', 'bg-green-100 text-green-700', 'bg-amber-100 text-amber-700', 'bg-pink-100 text-pink-700']

  const handleAdd = () => {
    const t = draft.trim()
    if (t && !tags.includes(t)) {
      onUpdate([...tags, t])
    }
    setDraft('')
    setAdding(false)
  }

  const handleRemove = (tag) => {
    onUpdate(tags.filter(t => t !== tag))
  }

  return (
    <span className="inline-flex flex-wrap items-center gap-1" onClick={e => e.stopPropagation()}>
      {tags.map((tag, i) => (
        <span key={tag} className={`inline-flex items-center gap-0.5 px-1.5 py-0.5 rounded-full text-[10px] font-medium ${colors[i % colors.length]}`}>
          {tag}
          <button onClick={() => handleRemove(tag)} className="ml-0.5 opacity-50 hover:opacity-100">×</button>
        </span>
      ))}
      {adding ? (
        <input
          ref={inputRef}
          type="text"
          value={draft}
          onChange={e => setDraft(e.target.value)}
          onBlur={handleAdd}
          onKeyDown={e => { if (e.key === 'Enter') handleAdd(); if (e.key === 'Escape') { setAdding(false); setDraft('') } }}
          className="text-[10px] px-1 py-0.5 border border-blue-300 rounded bg-white outline-none w-14"
          placeholder="tag"
        />
      ) : (
        <button
          onClick={() => setAdding(true)}
          className="w-4 h-4 rounded-full bg-gray-100 hover:bg-gray-200 flex items-center justify-center text-gray-400 hover:text-gray-600 text-xs"
          title="Add tag"
        >+</button>
      )}
    </span>
  )
}

// ---------------------------------------------------------------------------
// Cell value renderer
// ---------------------------------------------------------------------------
function CellValue({ col, value, colorClasses, mol, onAnnotation }) {
  // Action column (detail link) — rendered by parent <td>, returns null here
  if (col.type === 'action') return null

  // Invalidation toggle
  if (col.type === 'invalidation') {
    const checked = !!value
    return (
      <span className="flex items-center justify-center" onClick={e => e.stopPropagation()}>
        <input
          type="checkbox"
          checked={checked}
          onChange={() => onAnnotation && onAnnotation(mol.id, 'invalidated', !checked)}
          className="accent-red-500 cursor-pointer w-3.5 h-3.5"
          title={checked ? 'Mark as valid' : 'Mark as invalid'}
        />
      </span>
    )
  }

  // Tags
  if (col.type === 'tags') {
    return (
      <TagsCell
        tags={Array.isArray(value) ? value : []}
        onUpdate={(tags) => onAnnotation && onAnnotation(mol.id, 'tags', tags)}
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

  // Safety color code — render colored dot + label
  if (col.key === 'safety_color_code') {
    const lc = String(value).toLowerCase()
    const dotColor = lc === 'green' ? 'bg-green-500' : lc === 'yellow' || lc === 'orange' ? 'bg-yellow-400' : lc === 'red' ? 'bg-red-500' : 'bg-gray-400'
    return (
      <span className="inline-flex items-center gap-1.5">
        <span className={`w-2.5 h-2.5 rounded-full ${dotColor} flex-shrink-0`} />
        <span className="text-xs text-gray-600 capitalize">{lc}</span>
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
  onCellPopup,
  onAnnotation,
  activeRowId = null,
}) {
  const [sorts, setSorts] = useState([]) // [{key, dir}] max 2
  const [columnFilters, setColumnFilters] = useState({}) // {key: filterValue}
  const [showFilterCol, setShowFilterCol] = useState(null) // key of column with open filter
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

  // Column-filtered molecules
  const colFiltered = useMemo(() => {
    const filterKeys = Object.keys(columnFilters).filter(k => columnFilters[k] !== '')
    if (!filterKeys.length) return molecules
    return molecules.filter(mol => {
      return filterKeys.every(key => {
        const val = mol[key]
        const filterVal = columnFilters[key].toLowerCase()
        if (val == null) return false
        if (typeof val === 'number') return String(val).includes(filterVal)
        if (typeof val === 'boolean') return (val ? 'yes' : 'no').includes(filterVal)
        return String(val).toLowerCase().includes(filterVal)
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

              {/* Bookmark column — click to toggle bookmark on all selected */}
              <th className="px-1 py-2.5 border-b border-gray-200" style={{ width: 32 }}>
                <button
                  onClick={() => {
                    if (!onToggleBookmark) return
                    // Bookmark all selected, or all visible if none selected
                    const targets = selectedIds.size > 0
                      ? molecules.filter(m => selectedIds.has(m.id))
                      : molecules
                    const allBookmarked = targets.every(m => m.bookmarked)
                    targets.forEach(m => {
                      if (allBookmarked || !m.bookmarked) onToggleBookmark(m.id)
                    })
                  }}
                  disabled={!onToggleBookmark}
                  className={`p-0.5 rounded transition-colors mx-auto block ${onToggleBookmark ? 'hover:bg-yellow-50 cursor-pointer' : 'cursor-default opacity-50'}`}
                  title={selectedIds.size > 0 ? `Toggle bookmark on ${selectedIds.size} selected` : 'Toggle bookmark on all'}
                >
                  <svg className="w-3.5 h-3.5 text-gray-300" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                      d="M5 5a2 2 0 012-2h10a2 2 0 012 2v16l-7-3.5L5 21V5z" />
                  </svg>
                </button>
              </th>

              {/* Data columns */}
              {columns.map(col => {
                const dir = getSortDir(col.key)
                const rank = getSortRank(col.key)
                const isActive = dir !== null
                const hasFilter = columnFilters[col.key] && columnFilters[col.key] !== ''
                const isFilterable = col.type === 'number' || col.type === 'text' || col.type === 'boolean' || col.type === 'smiles'

                // Action column has no header label
                if (col.type === 'action') {
                  return <th key={col.key} className="w-9 px-1 py-2.5 border-b border-gray-200" style={{ width: 36 }} />
                }

                return (
                  <th
                    key={col.key}
                    className={`px-3 py-2.5 border-b border-gray-200 font-semibold uppercase tracking-wide text-[10px] whitespace-nowrap select-none relative ${
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
                      {TIPS[col.key] && <InfoTip text={TIPS[col.key]} size="xs" />}
                      {col.sortable && <SortIcon direction={dir} rank={rank} />}
                      {isFilterable && (
                        <button
                          onClick={e => { e.stopPropagation(); setShowFilterCol(prev => prev === col.key ? null : col.key) }}
                          className={`ml-0.5 p-0.5 rounded transition-colors ${hasFilter ? 'text-blue-500' : 'text-gray-300 hover:text-gray-500'}`}
                          title="Filter column"
                        >
                          <svg className="w-2.5 h-2.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M3 4a1 1 0 011-1h16a1 1 0 011 1v2.586a1 1 0 01-.293.707l-6.414 6.414a1 1 0 00-.293.707V17l-4 4v-6.586a1 1 0 00-.293-.707L3.293 7.293A1 1 0 013 6.586V4z" />
                          </svg>
                        </button>
                      )}
                    </span>
                    {/* Column filter dropdown */}
                    {showFilterCol === col.key && (
                      <div
                        className="absolute top-full left-0 z-20 mt-1 bg-white border border-gray-200 rounded-lg shadow-lg p-1.5"
                        onClick={e => e.stopPropagation()}
                      >
                        <input
                          type="text"
                          autoFocus
                          value={columnFilters[col.key] || ''}
                          onChange={e => setColumnFilters(prev => ({ ...prev, [col.key]: e.target.value }))}
                          onKeyDown={e => { if (e.key === 'Escape') setShowFilterCol(null) }}
                          placeholder={`Filter ${col.label}…`}
                          className="text-xs px-2 py-1 border border-gray-200 rounded w-28 outline-none focus:border-blue-300"
                        />
                        {hasFilter && (
                          <button
                            onClick={() => setColumnFilters(prev => { const n = { ...prev }; delete n[col.key]; return n })}
                            className="text-[10px] text-red-400 hover:text-red-600 mt-0.5 block"
                          >Clear</button>
                        )}
                      </div>
                    )}
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

              const isInvalidated = !!mol.invalidated

              return (
                <tr
                  key={mol.id}
                  className={`cursor-pointer transition-colors duration-100 group ${rowBg} hover:bg-blue-50/60 ${
                    isInvalidated ? 'border-l-[3px] border-red-300' :
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
                    const hasPopup = col.popup && onCellPopup && value != null

                    // Action column — detail link icon
                    if (col.type === 'action') {
                      return (
                        <td key={col.key} className="px-1 py-2 text-center w-9">
                          <button
                            onClick={e => { e.stopPropagation(); onRowClick && onRowClick(mol) }}
                            className="p-0.5 rounded hover:bg-blue-50 text-gray-300 hover:text-blue-500 transition-colors"
                            title="View details"
                          >
                            <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M15 12a3 3 0 11-6 0 3 3 0 016 0z" />
                              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M2.458 12C3.732 7.943 7.523 5 12 5c4.478 0 8.268 2.943 9.542 7-1.274 4.057-5.064 7-9.542 7-4.477 0-8.268-2.943-9.542-7z" />
                            </svg>
                          </button>
                        </td>
                      )
                    }

                    return (
                      <td
                        key={col.key}
                        className={`px-3 py-2 ${col.type === 'number' ? 'text-right tabular-nums' : 'text-left'} ${colorClasses.cell} ${hasPopup ? 'cursor-pointer hover:bg-blue-50/80 transition-colors' : ''} ${mol.invalidated && col.key !== 'invalidated' ? 'opacity-40' : ''}`}
                        onClick={hasPopup ? (e) => { e.stopPropagation(); onCellPopup(col.popup, mol) } : undefined}
                        title={hasPopup ? `Click for ${col.label} details` : undefined}
                      >
                        <span className={hasPopup ? 'underline decoration-dotted underline-offset-2 decoration-gray-300' : ''}>
                          <CellValue col={col} value={value} colorClasses={colorClasses} mol={mol} onAnnotation={onAnnotation} />
                        </span>
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
          {sorted.length !== molecules.length && (
            <>{' '}of{' '}<strong className="text-gray-600 tabular-nums">{molecules.length}</strong></>
          )}
          {' '}molecules
          {Object.keys(columnFilters).filter(k => columnFilters[k]).length > 0 && (
            <button
              onClick={() => setColumnFilters({})}
              className="ml-2 text-[10px] text-blue-500 hover:text-blue-700 underline"
            >Clear column filters</button>
          )}
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
