import React, { useState, useRef, useCallback, useMemo, useEffect } from 'react'
import useSettingsStore from '../stores/settingsStore.js'
import ParetoSettings from './ParetoSettings.jsx'

// --------------------------------------------------
// Extended objectives config
// --------------------------------------------------
export const ALL_OBJECTIVES = [
  { key: 'docking_score',   label: 'Docking Score',    extract: m => m.docking_score ?? null,    higher_is_better: false },
  { key: 'cnn_score',       label: 'CNN Score',        extract: m => m.cnn_score ?? null,        higher_is_better: true },
  { key: 'cnn_affinity',    label: 'CNN Affinity',     extract: m => m.cnn_affinity ?? null,     higher_is_better: true },
  { key: 'composite_score', label: 'Composite',        extract: m => m.composite_score ?? null,  higher_is_better: true },
  { key: 'QED',             label: 'QED',              extract: m => m.QED ?? null,              higher_is_better: true },
  { key: 'logP',            label: 'LogP',             extract: m => m.logP ?? null,             higher_is_better: false },
  { key: 'MW',              label: 'MW (Da)',          extract: m => m.MW ?? null,               higher_is_better: false },
  { key: 'TPSA',            label: 'TPSA',             extract: m => m.TPSA ?? null,             higher_is_better: false },
  { key: 'solubility',      label: 'Solubility',       extract: m => m.solubility ?? null,       higher_is_better: true },
  { key: 'synth_confidence',label: 'Synth. Confidence',extract: m => m.synth_confidence ?? null, higher_is_better: true },
  { key: 'confidence_score',label: 'Confidence',       extract: m => m.confidence_score ?? null, higher_is_better: true },
]

// Compute Pareto front (2-objective dominance)
function computeParetoFront(mols, xKey, yKey, objectives = ALL_OBJECTIVES) {
  const xObj = objectives.find(o => o.key === xKey)
  const yObj = objectives.find(o => o.key === yKey)
  if (!xObj || !yObj) return new Set()

  // Determine direction: we normalize so "better" = higher
  const xSign = xObj.higher_is_better === false ? -1 : 1
  const ySign = yObj.higher_is_better === false ? -1 : 1

  const points = mols.map((m, i) => ({
    idx: i,
    x: (xObj.extract(m) ?? -Infinity) * xSign,
    y: (yObj.extract(m) ?? -Infinity) * ySign,
  })).filter(p => p.x > -Infinity && p.y > -Infinity)

  const paretoIds = new Set()
  for (const p of points) {
    let dominated = false
    for (const q of points) {
      if (q === p) continue
      if (q.x >= p.x && q.y >= p.y && (q.x > p.x || q.y > p.y)) {
        dominated = true
        break
      }
    }
    if (!dominated) paretoIds.add(mols[p.idx].id)
  }
  return paretoIds
}

const COLOR_SCHEMES = {
  pareto_rank: { label: 'Pareto Rank', fn: (m, paretoIds) => paretoIds?.has(m.id) ? '#00e6a0' : '#9ca3af' },
  bookmarked:  { label: 'Bookmarked',  fn: m => m.bookmarked ? '#f59e0b' : '#d1d5db' },
  source:      { label: 'Source',      fn: m => {
    const s = (m.source || '').toLowerCase()
    if (s.includes('zinc')) return '#3b82f6'
    if (s.includes('chembl')) return '#00e6a0'
    if (s.includes('ai') || s.includes('generated')) return '#a855f7'
    return '#9ca3af'
  }},
  cluster:     { label: 'Cluster',     fn: m => {
    const colors = ['#ef4444','#3b82f6','#00e6a0','#f59e0b','#a855f7','#06b6d4','#ec4899','#84cc16']
    return colors[(m.cluster_id || 0) % colors.length]
  }},
  safety:      { label: 'Safety',      fn: m => {
    const c = (m.safety_color_code || '').toLowerCase()
    return c === 'green' ? '#00e6a0' : c === 'red' ? '#ef4444' : c === 'yellow' || c === 'orange' ? '#eab308' : '#9ca3af'
  }},
  invalidated: { label: 'Invalidated', fn: m => m.invalidated ? '#ef4444' : '#00e6a0' },
}

// --------------------------------------------------
// Inline SMILES preview for Pareto tooltip (lightweight, uses smiles-drawer)
// --------------------------------------------------
let _sdPromise = null
const _paretoSvgCache = new Map()

function ParetoSmilesPreview({ smiles }) {
  const [html, setHtml] = React.useState(() => _paretoSvgCache.get(smiles) || null)
  React.useEffect(() => {
    if (!smiles || _paretoSvgCache.has(smiles)) return
    let cancelled = false
    if (!_sdPromise) _sdPromise = import('smiles-drawer').then(m => m.default || m)
    _sdPromise.then(SD => {
      if (cancelled) return
      try {
        SD.parse(smiles, (tree) => {
          if (cancelled) return
          const svgEl = document.createElementNS('http://www.w3.org/2000/svg', 'svg')
          svgEl.setAttribute('width', '160')
          svgEl.setAttribute('height', '60')
          document.body.appendChild(svgEl)
          const drawer = new SD.SvgDrawer({ width: 160, height: 60 })
          drawer.draw(tree, svgEl, 'dark')
          const out = svgEl.outerHTML
          document.body.removeChild(svgEl)
          _paretoSvgCache.set(smiles, out)
          if (!cancelled) setHtml(out)
        }, () => {})
      } catch {}
    }).catch(() => {})
    return () => { cancelled = true }
  }, [smiles])

  if (!html) return <span style={{ color: '#6b7280', fontSize: '9px' }}>{smiles.length > 30 ? smiles.slice(0, 30) + '…' : smiles}</span>
  return <span dangerouslySetInnerHTML={{ __html: html }} />
}

// --------------------------------------------------
// Layout constants
// --------------------------------------------------
const VW = 560, VH = 340
const ML = 58, MR = 20, MT = 20, MB = 52
const PW = VW - ML - MR, PH = VH - MT - MB

// --------------------------------------------------
// Helpers
// --------------------------------------------------
function getVal(mol, key) {
  const obj = ALL_OBJECTIVES.find(o => o.key === key)
  return obj ? obj.extract(mol) : null
}

function linspace(a, b, n) {
  return Array.from({ length: n }, (_, i) => a + (i / (n - 1)) * (b - a))
}

function pearsonCorrelation(xs, ys) {
  const n = xs.length
  if (n < 3) return null
  const mx = xs.reduce((a, b) => a + b, 0) / n
  const my = ys.reduce((a, b) => a + b, 0) / n
  let num = 0, dx2 = 0, dy2 = 0
  for (let i = 0; i < n; i++) {
    const dx = xs[i] - mx, dy = ys[i] - my
    num += dx * dy; dx2 += dx * dx; dy2 += dy * dy
  }
  const den = Math.sqrt(dx2 * dy2)
  return den > 0 ? num / den : 0
}

function median(arr) {
  if (!arr.length) return 0
  const s = [...arr].sort((a, b) => a - b)
  const m = Math.floor(s.length / 2)
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2
}

function std(arr) {
  if (arr.length < 2) return 0
  const m = arr.reduce((a, b) => a + b, 0) / arr.length
  return Math.sqrt(arr.reduce((s, v) => s + (v - m) ** 2, 0) / (arr.length - 1))
}


// --------------------------------------------------
// Main ParetoFront component
// --------------------------------------------------
export default function ParetoFront({ molecules, onSelect, onToggleBookmark, visibleKeys }) {
  const paretoOverrides = useSettingsStore(s => s.paretoOverrides)
  const visibleSet = useMemo(() => visibleKeys ? new Set(visibleKeys) : null, [visibleKeys])

  // Effective objectives: filter by enabled + visible columns, override direction
  const effectiveObjectives = useMemo(() =>
    ALL_OBJECTIVES
      .filter(o => paretoOverrides[o.key]?.enabled !== false)
      .filter(o => !visibleSet || visibleSet.has(o.key))
      .map(o => {
        const ov = paretoOverrides[o.key]
        if (ov?.higher_is_better !== undefined) {
          return { ...o, higher_is_better: ov.higher_is_better }
        }
        return o
      }),
    [paretoOverrides, visibleSet]
  )

  const [xKey, setXKey] = useState('docking_score')
  const [yKey, setYKey] = useState('composite_score')
  const [colorBy, setColorBy] = useState('pareto_rank')
  const [sizeBy, setSizeBy] = useState('none')
  const [showDominated, setShowDominated] = useState(true)
  const [showStats, setShowStats] = useState(false)
  const [hovered, setHovered] = useState(null)
  const [hoveredPos, setHoveredPos] = useState({ x: 0, y: 0 })
  const [selected, setSelected] = useState([])
  const [zoom, setZoom] = useState(1)
  const [brushing, setBrushing] = useState(false)
  const [brushRect, setBrushRect] = useState(null)
  const [dragOrigin, setDragOrigin] = useState(null) // left-click drag tracking
  const [contextMenu, setContextMenu] = useState(null) // { x, y } screen coords
  const svgRef = useRef(null)
  const contextMenuRef = useRef(null)
  const hitMolRef = useRef(null) // molecule under cursor at mousedown

  // Auto-switch keys if objective disabled
  useEffect(() => {
    const keys = effectiveObjectives.map(o => o.key)
    if (keys.length === 0) return
    if (!keys.includes(xKey)) setXKey(keys[0])
    if (!keys.includes(yKey)) setYKey(keys[Math.min(1, keys.length - 1)])
  }, [effectiveObjectives]) // eslint-disable-line react-hooks/exhaustive-deps

  // Filter molecules that have at least one numeric property
  const mols = useMemo(() =>
    (molecules || []).filter(m => {
      return ALL_OBJECTIVES.some(o => o.extract(m) != null)
    }),
    [molecules]
  )

  // Compute pareto front locally (2-objective) using effective directions
  const paretoIds = useMemo(() => computeParetoFront(mols, xKey, yKey, effectiveObjectives), [mols, xKey, yKey, effectiveObjectives])

  const filteredMols = useMemo(() =>
    showDominated ? mols : mols.filter(m => paretoIds.has(m.id)),
    [mols, showDominated, paretoIds]
  )

  const paretoMols = mols.filter(m => paretoIds.has(m.id))
  const nonParetoMols = filteredMols.filter(m => !paretoIds.has(m.id))
  const paretoCount = paretoMols.length
  const selectedIdSet = useMemo(() => new Set(selected.map(m => m.id)), [selected])

  // Dynamic ranges
  const xValues = filteredMols.map(m => getVal(m, xKey)).filter(v => v != null)
  const yValues = filteredMols.map(m => getVal(m, yKey)).filter(v => v != null)
  const xMin = xValues.length ? Math.min(...xValues) : 0
  const xMax = xValues.length ? Math.max(...xValues) : 1
  const yMin = yValues.length ? Math.min(...yValues) : 0
  const yMax = yValues.length ? Math.max(...yValues) : 1
  const xPad = (xMax - xMin) * 0.05 || 0.1
  const yPad = (yMax - yMin) * 0.05 || 0.1
  const xRange = [xMin - xPad, xMax + xPad]
  const yRange = [yMin - yPad, yMax + yPad]

  const toX = useCallback(v => ML + ((v - xRange[0]) / (xRange[1] - xRange[0])) * PW, [xRange[0], xRange[1]]) // eslint-disable-line
  const toY = useCallback(v => MT + PH - ((v - yRange[0]) / (yRange[1] - yRange[0])) * PH, [yRange[0], yRange[1]]) // eslint-disable-line
  const xTicks = linspace(xRange[0], xRange[1], 5)
  const yTicks = linspace(yRange[0], yRange[1], 5)

  const colorFnRaw = COLOR_SCHEMES[colorBy]?.fn || COLOR_SCHEMES.pareto_rank.fn
  const colorFn = useCallback(mol => colorFnRaw(mol, paretoIds), [colorFnRaw, paretoIds])

  const getRadius = useCallback((mol) => {
    if (sizeBy === 'none') return 5
    const val = getVal(mol, sizeBy)
    if (val == null) return 4
    const allVals = filteredMols.map(m => getVal(m, sizeBy)).filter(v => v != null)
    if (!allVals.length) return 5
    const mn = Math.min(...allVals), mx = Math.max(...allVals)
    const range = mx - mn || 1
    return 3 + ((val - mn) / range) * 8
  }, [sizeBy, filteredMols])

  // Zoom handler
  const handleWheel = useCallback((e) => {
    e.preventDefault()
    const delta = e.deltaY > 0 ? 0.9 : 1.1
    setZoom(z => Math.max(0.5, Math.min(5, z * delta)))
  }, [])

  // Attach wheel as non-passive so preventDefault works
  useEffect(() => {
    const el = svgRef.current
    if (!el) return
    el.addEventListener('wheel', handleWheel, { passive: false })
    return () => el.removeEventListener('wheel', handleWheel)
  }, [handleWheel])

  // Close context menu on click outside
  useEffect(() => {
    if (!contextMenu) return
    const handler = (e) => {
      if (contextMenuRef.current && !contextMenuRef.current.contains(e.target)) {
        setContextMenu(null)
      }
    }
    document.addEventListener('mousedown', handler)
    return () => document.removeEventListener('mousedown', handler)
  }, [contextMenu])

  // Selection: click=toggle, drag=brush area, right-click=context menu
  const handleMouseDown = useCallback((e) => {
    if (e.button !== 0) return // only left button
    setContextMenu(null)
    const rect = svgRef.current.getBoundingClientRect()
    const svgX = (e.clientX - rect.left) / rect.width * VW
    const svgY = (e.clientY - rect.top) / rect.height * VH
    setDragOrigin({ clientX: e.clientX, clientY: e.clientY, svgX, svgY })
  }, [])

  const handleMouseMove = useCallback((e) => {
    if (!dragOrigin) return
    const dx = e.clientX - dragOrigin.clientX
    const dy = e.clientY - dragOrigin.clientY
    if (!brushing && (Math.abs(dx) > 5 || Math.abs(dy) > 5)) {
      setBrushing(true)
      hitMolRef.current = null // dragging — not a click
    }
    if (brushing || (Math.abs(dx) > 5 || Math.abs(dy) > 5)) {
      const rect = svgRef.current.getBoundingClientRect()
      const cx = (e.clientX - rect.left) / rect.width * VW
      const cy = (e.clientY - rect.top) / rect.height * VH
      setBrushRect({
        x: dragOrigin.svgX, y: dragOrigin.svgY,
        w: cx - dragOrigin.svgX, h: cy - dragOrigin.svgY,
      })
    }
  }, [dragOrigin, brushing])

  const handleMouseUp = useCallback((e) => {
    if (brushing && brushRect) {
      // Finish brush — add enclosed points to selection (cumulative)
      const bx = Math.min(brushRect.x, brushRect.x + brushRect.w)
      const by = Math.min(brushRect.y, brushRect.y + brushRect.h)
      const bw = Math.abs(brushRect.w)
      const bh = Math.abs(brushRect.h)
      if (bw > 5 && bh > 5) {
        // Account for zoom when comparing coords
        const z = zoom || 1
        const gbx = bx / z, gby = by / z, gbw = bw / z, gbh = bh / z
        const brushed = filteredMols.filter(mol => {
          const xv = getVal(mol, xKey), yv = getVal(mol, yKey)
          if (xv == null || yv == null) return false
          const px = toX(xv), py = toY(yv)
          return px >= gbx && px <= gbx + gbw && py >= gby && py <= gby + gbh
        })
        if (brushed.length > 0) {
          setSelected(prev => {
            const ids = new Set(prev.map(m => m.id))
            const toAdd = brushed.filter(m => !ids.has(m.id))
            return [...prev, ...toAdd]
          })
          setTimeout(() => setContextMenu({ x: e.clientX, y: e.clientY }), 50)
        }
      }
    } else if (hitMolRef.current) {
      // Click on a point — toggle it in/out of selection
      const mol = hitMolRef.current
      setSelected(prev =>
        prev.some(m => m.id === mol.id)
          ? prev.filter(m => m.id !== mol.id)
          : [...prev, mol]
      )
    } else if (dragOrigin && !brushing) {
      // Click on empty area — clear selection
      setSelected([])
    }
    hitMolRef.current = null
    setBrushing(false)
    setBrushRect(null)
    setDragOrigin(null)
  }, [brushing, brushRect, dragOrigin, filteredMols, xKey, yKey, toX, toY, zoom])

  // Double-click on a point: open detail panel
  const handlePointDblClick = useCallback((mol) => {
    if (onSelect) onSelect(mol)
  }, [onSelect])

  const handleBookmarkSelected = useCallback(() => {
    if (!onToggleBookmark) return
    // Only bookmark non-bookmarked molecules
    selected.filter(mol => !mol.bookmarked).forEach(mol => onToggleBookmark(mol.id))
    setContextMenu(null)
  }, [selected, onToggleBookmark])

  const handleUnbookmarkSelected = useCallback(() => {
    if (!onToggleBookmark) return
    // Only unbookmark bookmarked molecules
    selected.filter(mol => mol.bookmarked).forEach(mol => onToggleBookmark(mol.id))
    setContextMenu(null)
  }, [selected, onToggleBookmark])

  const handleSvgContextMenu = useCallback((e) => {
    e.preventDefault()
    if (selected.length > 0) setContextMenu({ x: e.clientX, y: e.clientY })
  }, [selected])

  const clearSelection = useCallback(() => setSelected([]), [])
  const handleResetZoom = useCallback(() => { setZoom(1) }, [])

  // Export CSV
  const exportCSV = useCallback(() => {
    const data = selected.length > 0 ? selected : filteredMols
    const keys = ALL_OBJECTIVES.map(o => o.key)
    const header = ['name', 'smiles', ...keys, 'pareto_optimal'].join(',')
    const rows = data.map(m => [
      `"${(m.name || '').replace(/"/g, '""')}"`,
      `"${(m.smiles || '').replace(/"/g, '""')}"`,
      ...keys.map(k => getVal(m, k) ?? ''),
      paretoIds.has(m.id) ? 'yes' : 'no'
    ].join(','))
    const csv = [header, ...rows].join('\n')
    const blob = new Blob([csv], { type: 'text/csv' })
    const url = URL.createObjectURL(blob)
    const a = document.createElement('a'); a.href = url; a.download = 'pareto_data.csv'; a.click()
    URL.revokeObjectURL(url)
  }, [selected, filteredMols])

  // Stats
  const stats = useMemo(() => {
    const data = selected.length > 0 ? selected : filteredMols
    const xs = data.map(m => getVal(m, xKey)).filter(v => v != null)
    const ys = data.map(m => getVal(m, yKey)).filter(v => v != null)
    return {
      n: data.length,
      xMean: xs.length ? (xs.reduce((a, b) => a + b, 0) / xs.length).toFixed(3) : 'N/A',
      xMedian: xs.length ? median(xs).toFixed(3) : 'N/A',
      xStd: xs.length > 1 ? std(xs).toFixed(3) : 'N/A',
      yMean: ys.length ? (ys.reduce((a, b) => a + b, 0) / ys.length).toFixed(3) : 'N/A',
      yMedian: ys.length ? median(ys).toFixed(3) : 'N/A',
      yStd: ys.length > 1 ? std(ys).toFixed(3) : 'N/A',
      correlation: xs.length > 2 && ys.length > 2 ? pearsonCorrelation(xs, ys)?.toFixed(3) : 'N/A',
      paretoPercent: mols.length > 0 ? ((paretoCount / mols.length) * 100).toFixed(1) : '0',
      paretoCount,
    }
  }, [selected, filteredMols, xKey, yKey, mols, paretoCount])

  if (!mols.length) {
    return (
      <div className="flex flex-col items-center justify-center py-10 text-gray-300">
        <svg className="w-10 h-10 mb-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
            d="M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 01-2 2h-2a2 2 0 01-2-2z" />
        </svg>
        <p className="text-sm">No Pareto data available</p>
      </div>
    )
  }

  const xLabel = ALL_OBJECTIVES.find(o => o.key === xKey)?.label || xKey
  const yLabel = ALL_OBJECTIVES.find(o => o.key === yKey)?.label || yKey

  // Scatter view
  const svgTransform = zoom !== 1 ? `scale(${zoom})` : undefined

  return (
    <div className="space-y-3">
      <Toolbar {...{ xKey, setXKey, yKey, setYKey, colorBy, setColorBy,
        sizeBy, setSizeBy, showDominated, setShowDominated, showStats, setShowStats,
        handleResetZoom, exportCSV, paretoCount, zoom, selected, effectiveObjectives, visibleSet }} />

      <div className="card overflow-hidden"
        style={{ cursor: 'crosshair' }}>
        <svg ref={svgRef} viewBox={`0 0 ${VW} ${VH}`} className="w-full" style={{ display: 'block' }}
          onMouseDown={handleMouseDown} onMouseMove={handleMouseMove}
          onMouseUp={handleMouseUp} onMouseLeave={handleMouseUp}
          onContextMenu={handleSvgContextMenu}>
          <rect x={0} y={0} width={VW} height={VH} fill="#fff" />

          <g transform={svgTransform}>
            <rect x={ML} y={MT} width={PW} height={PH} fill="#f9fafb" rx={4} />

            {/* Grid */}
            {xTicks.map((t, i) => {
              const x = toX(t)
              return (
                <g key={`xg-${i}`}>
                  <line x1={x} y1={MT} x2={x} y2={MT + PH} stroke="#e5e7eb"
                    strokeWidth={i === 0 || i === 4 ? 1.5 : 0.75} strokeDasharray={i === 0 || i === 4 ? '0' : '3,3'} />
                  <text x={x} y={MT + PH + 16} textAnchor="middle" fontSize={9} fill="#9ca3af">
                    {Math.abs(t) < 10 ? t.toFixed(2) : t.toFixed(1)}
                  </text>
                </g>
              )
            })}
            {yTicks.map((t, i) => {
              const y = toY(t)
              return (
                <g key={`yg-${i}`}>
                  <line x1={ML} y1={y} x2={ML + PW} y2={y} stroke="#e5e7eb"
                    strokeWidth={i === 0 || i === 4 ? 1.5 : 0.75} strokeDasharray={i === 0 || i === 4 ? '0' : '3,3'} />
                  <text x={ML - 8} y={y + 3.5} textAnchor="end" fontSize={9} fill="#9ca3af">
                    {Math.abs(t) < 10 ? t.toFixed(2) : t.toFixed(1)}
                  </text>
                </g>
              )
            })}

            {/* Pareto front line + shaded area */}
            {(() => {
              if (paretoMols.length < 2) return null
              const sorted = [...paretoMols].sort((a, b) => (getVal(a, xKey) ?? 0) - (getVal(b, xKey) ?? 0))
              const pts = sorted.map(m => ({ x: toX(getVal(m, xKey) ?? 0), y: toY(getVal(m, yKey) ?? 0) }))
              const linePoints = pts.map(p => `${p.x},${p.y}`).join(' ')
              // Shaded area under Pareto front
              const areaPath = pts.map((p, i) => `${i === 0 ? 'M' : 'L'} ${p.x} ${p.y}`).join(' ')
                + ` L ${pts[pts.length-1].x} ${MT + PH} L ${pts[0].x} ${MT + PH} Z`
              return (
                <>
                  <path d={areaPath} fill="#00e6a0" fillOpacity={0.06} />
                  <polyline points={linePoints} fill="none" stroke="#00e6a0" strokeWidth={1.5}
                    strokeDasharray="5,4" strokeLinecap="round" opacity={0.7} />
                </>
              )
            })()}

            {/* Non-pareto dots */}
            {nonParetoMols.map((mol, idx) => {
              const xv = getVal(mol, xKey), yv = getVal(mol, yKey)
              if (xv == null || yv == null) return null
              const cx = toX(xv), cy = toY(yv)
              const r = getRadius(mol)
              const fill = colorFn(mol)
              const isHov = hovered === mol
              const isSel = selectedIdSet.has(mol.id)
              return (
                <g key={`np-${idx}`}>
                  {isSel && <circle cx={cx} cy={cy} r={r + 5} fill="#3b82f6" fillOpacity={0.2} stroke="#3b82f6" strokeWidth={2} strokeDasharray="3,2" />}
                  <circle cx={cx} cy={cy} r={Math.max(r + 6, 10)} fill="transparent"
                    style={{ cursor: 'pointer' }}
                    onMouseEnter={() => { setHovered(mol); setHoveredPos({ x: cx, y: cy }) }}
                    onMouseLeave={() => setHovered(null)}
                    onMouseDown={(e) => { if (e.button === 0) hitMolRef.current = mol }}
                    onDoubleClick={() => handlePointDblClick(mol)} />
                  <circle cx={cx} cy={cy} r={isHov ? r + 2 : isSel ? r + 1 : r}
                    fill={isSel ? '#3b82f6' : fill} fillOpacity={isHov ? 0.9 : isSel ? 0.9 : 0.5}
                    stroke={isSel ? '#1d4ed8' : isHov ? '#4b5563' : 'none'} strokeWidth={isSel ? 2 : 1.5}
                    style={{ cursor: 'pointer', transition: 'r 0.1s', pointerEvents: 'none' }} />
                </g>
              )
            })}

            {/* Pareto dots */}
            {paretoMols.map((mol, idx) => {
              const xv = getVal(mol, xKey), yv = getVal(mol, yKey)
              if (xv == null || yv == null) return null
              const cx = toX(xv), cy = toY(yv)
              const r = getRadius(mol)
              const isHov = hovered === mol
              const isSel = selectedIdSet.has(mol.id)
              return (
                <g key={`p-${idx}`}>
                  {isSel && <circle cx={cx} cy={cy} r={r + 8} fill="#3b82f6" fillOpacity={0.2} stroke="#3b82f6" strokeWidth={2} strokeDasharray="3,2" />}
                  <circle cx={cx} cy={cy} r={Math.max(r + 8, 12)} fill="transparent"
                    style={{ cursor: 'pointer' }}
                    onMouseEnter={() => { setHovered(mol); setHoveredPos({ x: cx, y: cy }) }}
                    onMouseLeave={() => setHovered(null)}
                    onMouseDown={(e) => { if (e.button === 0) hitMolRef.current = mol }}
                    onDoubleClick={() => handlePointDblClick(mol)} />
                  <circle cx={cx} cy={cy} r={isHov ? r + 6 : r + 3} fill={isSel ? '#3b82f6' : '#00e6a0'} fillOpacity={isSel ? 0.25 : 0.15}
                    style={{ pointerEvents: 'none' }} />
                  <circle cx={cx} cy={cy} r={isHov ? r + 3 : isSel ? r + 2 : r + 1}
                    fill={isSel ? '#3b82f6' : '#00e6a0'} stroke={isSel ? '#1d4ed8' : '#fff'} strokeWidth={isSel ? 2 : 1.5}
                    style={{ cursor: 'pointer', transition: 'r 0.1s', pointerEvents: 'none' }} />
                  {mol.pareto_rank != null && (
                    <text x={cx} y={cy + 3.5} textAnchor="middle" fontSize={7} fontWeight={700} fill="#fff"
                      style={{ pointerEvents: 'none', userSelect: 'none' }}>{mol.pareto_rank}</text>
                  )}
                </g>
              )
            })}

            {/* Axis labels */}
            <text x={ML + PW / 2} y={VH - 6} textAnchor="middle" fontSize={11} fontWeight={600} fill="#374151">{xLabel}</text>
            <text x={12} y={MT + PH / 2} textAnchor="middle" fontSize={11} fontWeight={600} fill="#374151"
              transform={`rotate(-90, 12, ${MT + PH / 2})`}>{yLabel}</text>
          </g>

          {/* Brush rectangle */}
          {brushing && brushRect && (
            <rect x={Math.min(brushRect.x, brushRect.x + brushRect.w)}
              y={Math.min(brushRect.y, brushRect.y + brushRect.h)}
              width={Math.abs(brushRect.w)} height={Math.abs(brushRect.h)}
              fill="#3b82f6" fillOpacity={0.1} stroke="#3b82f6" strokeWidth={1} strokeDasharray="4,2" />
          )}

          {/* Selection count badge */}
          {selected.length > 0 && (
            <g>
              <rect x={ML + 6} y={MT + 6} width={84} height={24} rx={12} fill="#3b82f6" fillOpacity={0.9} />
              <text x={ML + 48} y={MT + 22} textAnchor="middle" fontSize={11} fontWeight={700} fill="#fff">
                {selected.length} selected
              </text>
            </g>
          )}

          {/* R² correlation badge */}
          {stats.correlation !== 'N/A' && (
            <g>
              <rect x={ML + PW - 80} y={MT + 6} width={74} height={22} rx={4} fill="#f1f5f9" stroke="#e2e8f0" strokeWidth={0.5} />
              <text x={ML + PW - 43} y={MT + 21} textAnchor="middle" fontSize={10} fill="#64748b" fontWeight={600}>
                R² = {(Number(stats.correlation) ** 2).toFixed(3)}
              </text>
            </g>
          )}

          {/* Tooltip */}
          {hovered && !brushing && (() => {
            const W = 200, H = 170
            let tx = hoveredPos.x + 14, ty = hoveredPos.y - H / 2
            if (tx + W > VW - 4) tx = hoveredPos.x - W - 14
            if (ty < MT) ty = MT
            if (ty + H > VH - 4) ty = VH - H - 4
            return (
              <foreignObject x={tx} y={ty} width={W} height={H} style={{ overflow: 'visible', pointerEvents: 'none' }}>
                <div xmlns="http://www.w3.org/1999/xhtml" style={{
                  background: '#0b1120', color: '#fff', borderRadius: '8px', padding: '8px 10px',
                  fontSize: '11px', lineHeight: '1.5', boxShadow: '0 4px 16px rgba(0,0,0,0.25)', width: `${W}px`,
                }}>
                  <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', marginBottom: '4px' }}>
                    <span style={{ fontWeight: 700, color: '#00e6a0', fontSize: '12px' }}>
                      {hovered.name || 'Molecule'}
                    </span>
                    {hovered.bookmarked && (
                      <span style={{ color: '#f59e0b', fontSize: '10px' }}>★</span>
                    )}
                  </div>
                  {paretoIds.has(hovered.id) && (
                    <div style={{ color: '#86efac', fontSize: '10px', marginBottom: '3px' }}>
                      Pareto-optimal
                    </div>
                  )}
                  {/* 2D structure preview */}
                  {hovered.smiles && (
                    <div style={{ marginBottom: '4px', background: '#1a2332', borderRadius: '4px', padding: '3px', textAlign: 'center' }}>
                      <ParetoSmilesPreview smiles={hovered.smiles} />
                    </div>
                  )}
                  <div style={{ fontSize: '10px', marginBottom: '4px' }}>
                    {ALL_OBJECTIVES.slice(0, 6).map(({ key, label, extract }) => {
                      const val = extract(hovered)
                      return (
                        <div key={key} style={{ display: 'flex', justifyContent: 'space-between', gap: '4px' }}>
                          <span style={{ opacity: 0.7 }}>{label}</span>
                          <span style={{ fontWeight: 600 }}>{val != null ? Number(val).toFixed(3) : 'N/A'}</span>
                        </div>
                      )
                    })}
                  </div>
                  <div style={{ fontSize: '10px', opacity: 0.6 }}>
                    {hovered.source && <span>Source: {hovered.source} </span>}
                    {hovered.cluster_id != null && <span>Family {hovered.cluster_id + 1}</span>}
                  </div>
                </div>
              </foreignObject>
            )
          })()}
        </svg>
      </div>

      {/* Context menu */}
      {contextMenu && selected.length > 0 && (() => {
        const unbookmarkedCount = selected.filter(m => !m.bookmarked).length
        const bookmarkedCount = selected.filter(m => m.bookmarked).length
        return (
          <div
            ref={contextMenuRef}
            className="fixed z-50 bg-white border border-gray-200 rounded-lg shadow-xl py-1 min-w-[180px]"
            style={{ left: contextMenu.x, top: contextMenu.y }}
          >
            <div className="px-3 py-1.5 text-xs text-gray-400 border-b border-gray-100">
              {selected.length} molecule{selected.length > 1 ? 's' : ''} selected
            </div>
            {unbookmarkedCount > 0 && (
              <button
                onClick={handleBookmarkSelected}
                className="w-full flex items-center gap-2 px-3 py-2 text-sm text-gray-700 hover:bg-yellow-50 hover:text-yellow-700 transition-colors"
              >
                <svg className="w-4 h-4 text-yellow-400" fill="currentColor" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                    d="M5 5a2 2 0 012-2h10a2 2 0 012 2v16l-7-3.5L5 21V5z" />
                </svg>
                Bookmark ({unbookmarkedCount})
              </button>
            )}
            {bookmarkedCount > 0 && (
              <button
                onClick={handleUnbookmarkSelected}
                className="w-full flex items-center gap-2 px-3 py-2 text-sm text-gray-700 hover:bg-red-50 hover:text-red-600 transition-colors"
              >
                <svg className="w-4 h-4 text-gray-300" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                    d="M5 5a2 2 0 012-2h10a2 2 0 012 2v16l-7-3.5L5 21V5z" />
                </svg>
                Remove bookmark ({bookmarkedCount})
              </button>
            )}
            <button
              onClick={() => { setSelected([]); setContextMenu(null) }}
              className="w-full flex items-center gap-2 px-3 py-2 text-sm text-gray-500 hover:bg-gray-50 transition-colors border-t border-gray-100"
            >
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
              </svg>
              Clear selection
            </button>
          </div>
        )
      })()}

      {/* Legend */}
      <div className="flex flex-wrap items-center gap-5 text-sm text-gray-500 px-1">
        <div className="flex items-center gap-2">
          <svg width={24} height={12}>
            <circle cx={6} cy={6} r={5} fill="#00e6a0" stroke="#fff" strokeWidth={1.5} />
            <line x1={10} y1={6} x2={22} y2={6} stroke="#00e6a0" strokeWidth={1.5} strokeDasharray="3,2" />
          </svg>
          <span>Pareto front</span>
        </div>
        <div className="flex items-center gap-2">
          <svg width={14} height={12}><circle cx={6} cy={6} r={4} fill="#9ca3af" fillOpacity={0.5} /></svg>
          <span>Dominated</span>
        </div>
        <span className="text-gray-300">Click: select/deselect | Drag: area select | Right-click: actions</span>
      </div>

      {showStats && <StatsPanel stats={stats} xLabel={xLabel} yLabel={yLabel} />}
    </div>
  )
}

// --------------------------------------------------
// Toolbar component
// --------------------------------------------------
function Toolbar({ xKey, setXKey, yKey, setYKey, colorBy, setColorBy,
  sizeBy, setSizeBy, showDominated, setShowDominated, showStats, setShowStats,
  handleResetZoom, exportCSV, paretoCount, zoom, selected, clearSelection, effectiveObjectives = ALL_OBJECTIVES, visibleSet }) {
  const [showParetoSettings, setShowParetoSettings] = React.useState(false)
  const paretoOverrides = useSettingsStore(s => s.paretoOverrides)
  const hasOverrides = Object.keys(paretoOverrides).length > 0

  return (
    <div className="flex flex-wrap items-center gap-3">
      <div className="flex items-center gap-1.5">
        <label className="text-sm font-medium text-gray-500">X:</label>
        <select value={xKey} onChange={e => setXKey(e.target.value)}
          className="text-sm border border-gray-200 rounded-md px-1.5 py-1 bg-white text-gray-700 focus:outline-none focus:ring-1 focus:ring-bx-mint">
          {effectiveObjectives.map(o => <option key={o.key} value={o.key}>{o.label}</option>)}
        </select>
      </div>
      <div className="flex items-center gap-1.5">
        <label className="text-sm font-medium text-gray-500">Y:</label>
        <select value={yKey} onChange={e => setYKey(e.target.value)}
          className="text-sm border border-gray-200 rounded-md px-1.5 py-1 bg-white text-gray-700 focus:outline-none focus:ring-1 focus:ring-bx-mint">
          {effectiveObjectives.map(o => <option key={o.key} value={o.key}>{o.label}</option>)}
        </select>
      </div>
      <select value={colorBy} onChange={e => setColorBy(e.target.value)}
        className="text-sm border border-gray-200 rounded-md px-1.5 py-1 bg-white text-gray-700">
        {Object.entries(COLOR_SCHEMES).map(([k, v]) => <option key={k} value={k}>{v.label}</option>)}
      </select>
      <select value={sizeBy} onChange={e => setSizeBy(e.target.value)}
        className="text-sm border border-gray-200 rounded-md px-1.5 py-1 bg-white text-gray-700">
        <option value="none">Size: off</option>
        {effectiveObjectives.map(o => <option key={o.key} value={o.key}>Size: {o.label}</option>)}
      </select>
      <label className="flex items-center gap-1 text-sm text-gray-500 cursor-pointer">
        <input type="checkbox" checked={showDominated} onChange={e => setShowDominated(e.target.checked)}
          className="accent-bx-mint" />
        Dominated
      </label>
      <div className="ml-auto flex items-center gap-2">
        {zoom > 1 && (
          <button onClick={handleResetZoom} className="text-sm text-blue-600 hover:underline">Reset zoom</button>
        )}
        <button onClick={exportCSV}
          className="flex items-center gap-1 px-2 py-1 text-sm font-medium text-gray-600 bg-gray-100 hover:bg-gray-200 rounded-lg transition-colors">
          <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 10v6m0 0l-3-3m3 3l3-3m2 8H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
          </svg>
          CSV
        </button>
        <button onClick={() => setShowStats(v => !v)}
          className={`px-2 py-1 text-sm font-medium rounded-lg transition-colors ${showStats ? 'bg-bx-surface text-white' : 'bg-gray-100 text-gray-600 hover:bg-gray-200'}`}>
          Stats
        </button>
        <div className="relative">
          <button
            onClick={() => setShowParetoSettings(v => !v)}
            className={`flex items-center gap-1 px-2 py-1 text-sm font-medium rounded-lg transition-colors ${
              showParetoSettings ? 'bg-bx-surface text-white' : 'bg-gray-100 text-gray-600 hover:bg-gray-200'
            }`}
            title="Configure Pareto objectives"
          >
            <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                d="M10.325 4.317c.426-1.756 2.924-1.756 3.35 0a1.724 1.724 0 002.573 1.066c1.543-.94 3.31.826 2.37 2.37a1.724 1.724 0 001.066 2.573c1.756.426 1.756 2.924 0 3.35a1.724 1.724 0 00-1.066 2.573c.94 1.543-.826 3.31-2.37 2.37a1.724 1.724 0 00-2.573 1.066c-.426 1.756-2.924 1.756-3.35 0a1.724 1.724 0 00-2.573-1.066c-1.543.94-3.31-.826-2.37-2.37a1.724 1.724 0 00-1.066-2.573c-1.756-.426-1.756-2.924 0-3.35a1.724 1.724 0 001.066-2.573c-.94-1.543.826-3.31 2.37-2.37.996.608 2.296.07 2.572-1.065z" />
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M15 12a3 3 0 11-6 0 3 3 0 016 0z" />
            </svg>
            {hasOverrides && <span className="w-1.5 h-1.5 rounded-full bg-bx-cyan" />}
          </button>
          <ParetoSettings
            isOpen={showParetoSettings}
            onClose={() => setShowParetoSettings(false)}
            allObjectives={visibleSet ? ALL_OBJECTIVES.filter(o => visibleSet.has(o.key)) : ALL_OBJECTIVES}
          />
        </div>
        {paretoCount > 0 && (
          <span className="text-sm text-gray-400">
            <span className="font-semibold text-bx-mint">{paretoCount}</span> on front
          </span>
        )}
        {selected.length > 0 && (
          <span className="flex items-center gap-1.5 text-sm text-blue-600 font-medium">
            {selected.length} selected
            <button onClick={clearSelection} className="text-blue-400 hover:text-blue-600" title="Clear selection">
              <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
              </svg>
            </button>
          </span>
        )}
      </div>
    </div>
  )
}

// --------------------------------------------------
// Stats panel
// --------------------------------------------------
function StatsPanel({ stats, xLabel, yLabel }) {
  return (
    <div className="card p-4">
      <h4 className="text-sm font-semibold text-gray-400 uppercase tracking-wider mb-3">Statistics</h4>
      <div className="grid grid-cols-2 sm:grid-cols-4 gap-4 text-sm">
        <div>
          <p className="text-gray-400 mb-1">Points</p>
          <p className="font-bold text-bx-light-text">{stats.n}</p>
        </div>
        <div>
          <p className="text-gray-400 mb-1">Correlation (r)</p>
          <p className="font-bold text-bx-light-text">{stats.correlation}</p>
        </div>
        <div>
          <p className="text-gray-400 mb-1">Pareto front</p>
          <p className="font-bold text-bx-mint">{stats.paretoCount} ({stats.paretoPercent}%)</p>
        </div>
        <div>
          <p className="text-gray-400 mb-1">{xLabel} mean</p>
          <p className="font-bold text-gray-700">{stats.xMean}</p>
        </div>
        <div>
          <p className="text-gray-400 mb-1">{xLabel} median</p>
          <p className="font-bold text-gray-700">{stats.xMedian}</p>
        </div>
        <div>
          <p className="text-gray-400 mb-1">{xLabel} std</p>
          <p className="font-bold text-gray-700">{stats.xStd}</p>
        </div>
        <div>
          <p className="text-gray-400 mb-1">{yLabel} mean</p>
          <p className="font-bold text-gray-700">{stats.yMean}</p>
        </div>
        <div>
          <p className="text-gray-400 mb-1">{yLabel} std</p>
          <p className="font-bold text-gray-700">{stats.yStd}</p>
        </div>
      </div>
    </div>
  )
}
