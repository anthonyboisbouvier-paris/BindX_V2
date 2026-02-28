import React, { useState, useRef, useCallback, useMemo, useEffect } from 'react'

// --------------------------------------------------
// Extended objectives config
// --------------------------------------------------
const ALL_OBJECTIVES = [
  { key: 'affinity',        label: 'Affinity',         extract: m => m.pareto_objectives?.affinity ?? null },
  { key: 'safety',          label: 'Safety',           extract: m => m.pareto_objectives?.safety ?? null },
  { key: 'bioavailability', label: 'Bioavailability',  extract: m => m.pareto_objectives?.bioavailability ?? null },
  { key: 'synthesis',       label: 'Synthesizability', extract: m => m.pareto_objectives?.synthesis ?? null },
  { key: 'qed',             label: 'QED',              extract: m => m.qed ?? null },
  { key: 'logp',            label: 'LogP',             extract: m => m.logp ?? null },
  { key: 'mw',              label: 'MW (Da)',          extract: m => m.mw ?? null },
  { key: 'sa_score',        label: 'SA Score',         extract: m => m.sa_score ?? null },
  { key: 'composite_score', label: 'Composite',        extract: m => m.composite_score ?? null },
  { key: 'cnn_score',       label: 'CNN Score',        extract: m => m.cnn_score ?? null },
  { key: 'selectivity',     label: 'Selectivity',      extract: m => m.off_target?.selectivity_score ?? null },
]

const COLOR_SCHEMES = {
  pareto_rank: { label: 'Pareto Rank', fn: m => m.pareto_front ? '#00e6a0' : '#9ca3af' },
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
  admet:       { label: 'ADMET',       fn: m => {
    const c = m.admet?.color_code
    return c === 'green' ? '#00e6a0' : c === 'red' ? '#ef4444' : '#eab308'
  }},
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
// Mini radar for tooltip
// --------------------------------------------------
function MiniRadar({ mol, size = 60 }) {
  const axes = ALL_OBJECTIVES.slice(0, 6)
  const cx = size / 2, cy = size / 2, r = size / 2 - 6
  const points = axes.map((a, i) => {
    const angle = (Math.PI * 2 * i) / axes.length - Math.PI / 2
    const val = Math.max(0, Math.min(1, a.extract(mol) ?? 0))
    return { x: cx + Math.cos(angle) * r * val, y: cy + Math.sin(angle) * r * val }
  })
  const path = points.map((p, i) => `${i === 0 ? 'M' : 'L'} ${p.x} ${p.y}`).join(' ') + ' Z'
  return (
    <svg width={size} height={size}>
      {/* Grid */}
      {[0.25, 0.5, 0.75, 1].map(s => (
        <polygon key={s} points={axes.map((_, i) => {
          const a = (Math.PI * 2 * i) / axes.length - Math.PI / 2
          return `${cx + Math.cos(a) * r * s},${cy + Math.sin(a) * r * s}`
        }).join(' ')} fill="none" stroke="#e5e7eb" strokeWidth={0.5} />
      ))}
      {/* Data polygon */}
      <path d={path} fill="#00e6a0" fillOpacity={0.25} stroke="#00e6a0" strokeWidth={1.5} />
      {points.map((p, i) => (
        <circle key={i} cx={p.x} cy={p.y} r={2} fill="#00e6a0" />
      ))}
    </svg>
  )
}

// --------------------------------------------------
// Distribution histogram
// --------------------------------------------------
function Histogram({ values, width, height, color = '#3b82f6' }) {
  if (!values.length) return null
  const nBins = Math.min(15, Math.max(5, Math.ceil(Math.sqrt(values.length))))
  const min = Math.min(...values), max = Math.max(...values)
  const range = max - min || 1
  const binW = range / nBins
  const bins = Array(nBins).fill(0)
  values.forEach(v => {
    const idx = Math.min(nBins - 1, Math.floor((v - min) / binW))
    bins[idx]++
  })
  const maxCount = Math.max(...bins)
  const barW = width / nBins - 1
  return (
    <g>
      {bins.map((count, i) => {
        const h = maxCount > 0 ? (count / maxCount) * height : 0
        return (
          <rect key={i} x={i * (barW + 1)} y={height - h} width={barW} height={h}
            fill={color} fillOpacity={0.6} rx={1} />
        )
      })}
    </g>
  )
}

// --------------------------------------------------
// Radar overlay view
// --------------------------------------------------
function RadarOverlay({ molecules, width, height }) {
  const axes = ALL_OBJECTIVES.slice(0, 8)
  const cx = width / 2, cy = height / 2, r = Math.min(width, height) / 2 - 30
  const colors = ['#00e6a0', '#3b82f6', '#a855f7']

  return (
    <svg width={width} height={height} className="mx-auto">
      {/* Grid circles */}
      {[0.25, 0.5, 0.75, 1].map(s => (
        <polygon key={s} points={axes.map((_, i) => {
          const a = (Math.PI * 2 * i) / axes.length - Math.PI / 2
          return `${cx + Math.cos(a) * r * s},${cy + Math.sin(a) * r * s}`
        }).join(' ')} fill="none" stroke="#e5e7eb" strokeWidth={0.5} />
      ))}
      {/* Axis lines + labels */}
      {axes.map((ax, i) => {
        const a = (Math.PI * 2 * i) / axes.length - Math.PI / 2
        const ex = cx + Math.cos(a) * r, ey = cy + Math.sin(a) * r
        const lx = cx + Math.cos(a) * (r + 18), ly = cy + Math.sin(a) * (r + 18)
        return (
          <g key={ax.key}>
            <line x1={cx} y1={cy} x2={ex} y2={ey} stroke="#d1d5db" strokeWidth={0.5} />
            <text x={lx} y={ly} textAnchor="middle" fontSize={8} fill="#6b7280">{ax.label}</text>
          </g>
        )
      })}
      {/* Molecule polygons */}
      {molecules.map((mol, mi) => {
        const pts = axes.map((ax, i) => {
          const a = (Math.PI * 2 * i) / axes.length - Math.PI / 2
          const val = Math.max(0, Math.min(1, ax.extract(mol) ?? 0))
          return { x: cx + Math.cos(a) * r * val, y: cy + Math.sin(a) * r * val }
        })
        const path = pts.map((p, i) => `${i === 0 ? 'M' : 'L'} ${p.x} ${p.y}`).join(' ') + ' Z'
        return (
          <g key={mi}>
            <path d={path} fill={colors[mi]} fillOpacity={0.15} stroke={colors[mi]} strokeWidth={2} />
            {pts.map((p, i) => (
              <circle key={i} cx={p.x} cy={p.y} r={3} fill={colors[mi]} />
            ))}
          </g>
        )
      })}
      {/* Legend */}
      {molecules.map((mol, mi) => (
        <text key={mi} x={10} y={20 + mi * 14} fontSize={10} fill={colors[mi]} fontWeight={600}>
          {mol.name || `Molecule ${mi + 1}`}
        </text>
      ))}
    </svg>
  )
}

// --------------------------------------------------
// Main ParetoFront component
// --------------------------------------------------
export default function ParetoFront({ molecules, onSelect }) {
  const [xKey, setXKey] = useState('affinity')
  const [yKey, setYKey] = useState('safety')
  const [viewMode, setViewMode] = useState('scatter')
  const [colorBy, setColorBy] = useState('pareto_rank')
  const [sizeBy, setSizeBy] = useState('none')
  const [showDominated, setShowDominated] = useState(true)
  const [showStats, setShowStats] = useState(false)
  const [hovered, setHovered] = useState(null)
  const [hoveredPos, setHoveredPos] = useState({ x: 0, y: 0 })
  const [selected, setSelected] = useState([])
  const [zoom, setZoom] = useState(1)
  const [pan, setPan] = useState({ x: 0, y: 0 })
  const [isPanning, setIsPanning] = useState(false)
  const [panStart, setPanStart] = useState(null)
  const [brushing, setBrushing] = useState(false)
  const [brushRect, setBrushRect] = useState(null)
  const svgRef = useRef(null)

  const mols = useMemo(() =>
    (molecules || []).filter(m => m.pareto_objectives && typeof m.pareto_objectives === 'object'),
    [molecules]
  )

  const filteredMols = useMemo(() =>
    showDominated ? mols : mols.filter(m => m.pareto_front === true),
    [mols, showDominated]
  )

  const paretoMols = mols.filter(m => m.pareto_front === true)
  const nonParetoMols = filteredMols.filter(m => m.pareto_front !== true)
  const paretoCount = paretoMols.length

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

  const colorFn = COLOR_SCHEMES[colorBy]?.fn || COLOR_SCHEMES.pareto_rank.fn

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

  // Pan handlers
  const handleMouseDown = useCallback((e) => {
    if (e.shiftKey) {
      // Brush mode
      const rect = svgRef.current.getBoundingClientRect()
      const sx = (e.clientX - rect.left) / rect.width * VW
      const sy = (e.clientY - rect.top) / rect.height * VH
      setBrushing(true)
      setBrushRect({ x: sx, y: sy, w: 0, h: 0 })
      return
    }
    if (zoom > 1) {
      setIsPanning(true)
      setPanStart({ x: e.clientX - pan.x, y: e.clientY - pan.y })
    }
  }, [zoom, pan])

  const handleMouseMove = useCallback((e) => {
    if (brushing && brushRect) {
      const rect = svgRef.current.getBoundingClientRect()
      const cx = (e.clientX - rect.left) / rect.width * VW
      const cy = (e.clientY - rect.top) / rect.height * VH
      setBrushRect(prev => ({ ...prev, w: cx - prev.x, h: cy - prev.y }))
      return
    }
    if (isPanning && panStart) {
      setPan({ x: e.clientX - panStart.x, y: e.clientY - panStart.y })
    }
  }, [isPanning, panStart, brushing, brushRect])

  const handleMouseUp = useCallback(() => {
    if (brushing && brushRect) {
      // Select points in brush rect
      const bx = Math.min(brushRect.x, brushRect.x + brushRect.w)
      const by = Math.min(brushRect.y, brushRect.y + brushRect.h)
      const bw = Math.abs(brushRect.w)
      const bh = Math.abs(brushRect.h)
      if (bw > 5 && bh > 5) {
        const brushed = filteredMols.filter(mol => {
          const xv = getVal(mol, xKey), yv = getVal(mol, yKey)
          if (xv == null || yv == null) return false
          const px = toX(xv), py = toY(yv)
          return px >= bx && px <= bx + bw && py >= by && py <= by + bh
        })
        setSelected(brushed.slice(0, 3))
      }
      setBrushing(false)
      setBrushRect(null)
      return
    }
    setIsPanning(false)
    setPanStart(null)
  }, [brushing, brushRect, filteredMols, xKey, yKey, toX, toY])

  const handlePointClick = useCallback((mol, e) => {
    if (e.ctrlKey || e.metaKey) {
      setSelected(prev => {
        if (prev.includes(mol)) return prev.filter(m => m !== mol)
        if (prev.length >= 3) return [...prev.slice(1), mol]
        return [...prev, mol]
      })
    } else {
      setSelected([mol])
      if (onSelect) onSelect(mol)
    }
  }, [onSelect])

  const handleResetZoom = useCallback(() => { setZoom(1); setPan({ x: 0, y: 0 }) }, [])

  // Export CSV
  const exportCSV = useCallback(() => {
    const data = selected.length > 0 ? selected : filteredMols
    const keys = ALL_OBJECTIVES.map(o => o.key)
    const header = ['name', 'smiles', ...keys, 'pareto_rank', 'pareto_front'].join(',')
    const rows = data.map(m => [
      `"${(m.name || '').replace(/"/g, '""')}"`,
      `"${(m.smiles || '').replace(/"/g, '""')}"`,
      ...keys.map(k => getVal(m, k) ?? ''),
      m.pareto_rank ?? '', m.pareto_front ?? ''
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

  // Radar view
  if (viewMode === 'radar') {
    const radarMols = selected.length > 0 ? selected : paretoMols.slice(0, 3)
    return (
      <div className="space-y-3">
        <Toolbar {...{ xKey, setXKey, yKey, setYKey, viewMode, setViewMode, colorBy, setColorBy,
          sizeBy, setSizeBy, showDominated, setShowDominated, showStats, setShowStats,
          handleResetZoom, exportCSV, paretoCount, zoom, selected }} />
        <div className="bg-white rounded-xl border border-gray-100 overflow-hidden shadow-sm p-4">
          <RadarOverlay molecules={radarMols} width={400} height={350} />
          <p className="text-xs text-gray-400 text-center mt-2">
            {radarMols.length === 0 ? 'Select molecules to compare' : `Comparing ${radarMols.length} molecule(s)`}
            {' — Ctrl+click points in scatter to select up to 3'}
          </p>
        </div>
        {showStats && <StatsPanel stats={stats} xLabel={xLabel} yLabel={yLabel} />}
      </div>
    )
  }

  // Distribution view
  if (viewMode === 'distribution') {
    return (
      <div className="space-y-3">
        <Toolbar {...{ xKey, setXKey, yKey, setYKey, viewMode, setViewMode, colorBy, setColorBy,
          sizeBy, setSizeBy, showDominated, setShowDominated, showStats, setShowStats,
          handleResetZoom, exportCSV, paretoCount, zoom, selected }} />
        <div className="bg-white rounded-xl border border-gray-100 overflow-hidden shadow-sm">
          <div className="grid grid-cols-2 gap-4 p-4">
            <div>
              <p className="text-xs font-semibold text-gray-500 mb-2">{xLabel} Distribution</p>
              <svg width="100%" viewBox={`0 0 ${PW} 80`}>
                <Histogram values={xValues} width={PW} height={70} color="#3b82f6" />
              </svg>
            </div>
            <div>
              <p className="text-xs font-semibold text-gray-500 mb-2">{yLabel} Distribution</p>
              <svg width="100%" viewBox={`0 0 ${PW} 80`}>
                <Histogram values={yValues} width={PW} height={70} color="#00e6a0" />
              </svg>
            </div>
          </div>
        </div>
        {showStats && <StatsPanel stats={stats} xLabel={xLabel} yLabel={yLabel} />}
      </div>
    )
  }

  // Scatter view (default)
  const svgTransform = zoom !== 1 || pan.x !== 0 || pan.y !== 0
    ? `translate(${pan.x},${pan.y}) scale(${zoom})` : undefined

  return (
    <div className="space-y-3">
      <Toolbar {...{ xKey, setXKey, yKey, setYKey, viewMode, setViewMode, colorBy, setColorBy,
        sizeBy, setSizeBy, showDominated, setShowDominated, showStats, setShowStats,
        handleResetZoom, exportCSV, paretoCount, zoom, selected }} />

      <div className="bg-white rounded-xl border border-gray-100 overflow-hidden shadow-sm"
        style={{ cursor: zoom > 1 ? (isPanning ? 'grabbing' : 'grab') : 'crosshair' }}>
        <svg ref={svgRef} viewBox={`0 0 ${VW} ${VH}`} className="w-full" style={{ display: 'block' }}
          onMouseDown={handleMouseDown} onMouseMove={handleMouseMove}
          onMouseUp={handleMouseUp} onMouseLeave={handleMouseUp}>
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
              const isSel = selected.includes(mol)
              return (
                <circle key={`np-${idx}`} cx={cx} cy={cy} r={isHov ? r + 2 : isSel ? r + 1 : r}
                  fill={fill} fillOpacity={isHov ? 0.9 : isSel ? 0.8 : 0.5}
                  stroke={isSel ? '#0f131d' : isHov ? '#4b5563' : 'none'} strokeWidth={isSel ? 2 : 1.5}
                  style={{ cursor: 'pointer', transition: 'r 0.1s' }}
                  onMouseEnter={() => { setHovered(mol); setHoveredPos({ x: cx, y: cy }) }}
                  onMouseLeave={() => setHovered(null)}
                  onClick={(e) => handlePointClick(mol, e)} />
              )
            })}

            {/* Pareto dots */}
            {paretoMols.map((mol, idx) => {
              const xv = getVal(mol, xKey), yv = getVal(mol, yKey)
              if (xv == null || yv == null) return null
              const cx = toX(xv), cy = toY(yv)
              const r = getRadius(mol)
              const isHov = hovered === mol
              const isSel = selected.includes(mol)
              return (
                <g key={`p-${idx}`}>
                  <circle cx={cx} cy={cy} r={isHov ? r + 6 : r + 3} fill="#00e6a0" fillOpacity={0.15} />
                  <circle cx={cx} cy={cy} r={isHov ? r + 3 : isSel ? r + 1 : r + 1}
                    fill="#00e6a0" stroke={isSel ? '#0f131d' : '#fff'} strokeWidth={isSel ? 2 : 1.5}
                    style={{ cursor: 'pointer', transition: 'r 0.1s' }}
                    onMouseEnter={() => { setHovered(mol); setHoveredPos({ x: cx, y: cy }) }}
                    onMouseLeave={() => setHovered(null)}
                    onClick={(e) => handlePointClick(mol, e)} />
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

          {/* Tooltip */}
          {hovered && !brushing && (() => {
            const W = 190, H = 160
            let tx = hoveredPos.x + 14, ty = hoveredPos.y - H / 2
            if (tx + W > VW - 4) tx = hoveredPos.x - W - 14
            if (ty < MT) ty = MT
            if (ty + H > VH - 4) ty = VH - H - 4
            return (
              <foreignObject x={tx} y={ty} width={W} height={H} style={{ overflow: 'visible', pointerEvents: 'none' }}>
                <div xmlns="http://www.w3.org/1999/xhtml" style={{
                  background: '#0f131d', color: '#fff', borderRadius: '8px', padding: '8px 10px',
                  fontSize: '11px', lineHeight: '1.5', boxShadow: '0 4px 16px rgba(0,0,0,0.25)', width: `${W}px`,
                }}>
                  <div style={{ fontWeight: 700, marginBottom: '4px', color: '#00e6a0', fontSize: '12px' }}>
                    {hovered.name || 'Molecule'}
                  </div>
                  {hovered.pareto_front && (
                    <div style={{ color: '#86efac', fontSize: '10px', marginBottom: '3px' }}>
                      Pareto front — rank {hovered.pareto_rank ?? '?'}
                    </div>
                  )}
                  <div style={{ display: 'flex', gap: '8px', marginBottom: '4px' }}>
                    <MiniRadar mol={hovered} size={50} />
                    <div style={{ flex: 1, fontSize: '10px' }}>
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

      {/* Legend */}
      <div className="flex flex-wrap items-center gap-5 text-xs text-gray-500 px-1">
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
        <span className="text-gray-300">Shift+drag: brush select | Ctrl+click: multi-select | Scroll: zoom</span>
      </div>

      {showStats && <StatsPanel stats={stats} xLabel={xLabel} yLabel={yLabel} />}
    </div>
  )
}

// --------------------------------------------------
// Toolbar component
// --------------------------------------------------
function Toolbar({ xKey, setXKey, yKey, setYKey, viewMode, setViewMode, colorBy, setColorBy,
  sizeBy, setSizeBy, showDominated, setShowDominated, showStats, setShowStats,
  handleResetZoom, exportCSV, paretoCount, zoom, selected }) {
  return (
    <div className="flex flex-wrap items-center gap-3">
      <div className="flex items-center gap-1.5">
        <label className="text-xs font-medium text-gray-500">X:</label>
        <select value={xKey} onChange={e => setXKey(e.target.value)}
          className="text-xs border border-gray-200 rounded-md px-1.5 py-1 bg-white text-gray-700 focus:outline-none focus:ring-1 focus:ring-[#0f131d]">
          {ALL_OBJECTIVES.map(o => <option key={o.key} value={o.key}>{o.label}</option>)}
        </select>
      </div>
      <div className="flex items-center gap-1.5">
        <label className="text-xs font-medium text-gray-500">Y:</label>
        <select value={yKey} onChange={e => setYKey(e.target.value)}
          className="text-xs border border-gray-200 rounded-md px-1.5 py-1 bg-white text-gray-700 focus:outline-none focus:ring-1 focus:ring-[#0f131d]">
          {ALL_OBJECTIVES.map(o => <option key={o.key} value={o.key}>{o.label}</option>)}
        </select>
      </div>
      <div className="flex gap-0.5 p-0.5 bg-gray-100 rounded-lg">
        {[{ v: 'scatter', l: 'Scatter' }, { v: 'radar', l: 'Radar' }, { v: 'distribution', l: 'Dist.' }].map(({ v, l }) => (
          <button key={v} onClick={() => setViewMode(v)}
            className={`px-2 py-0.5 rounded text-xs font-medium transition-colors ${viewMode === v ? 'bg-white text-[#0f131d] shadow-sm' : 'text-gray-500 hover:text-gray-700'}`}>
            {l}
          </button>
        ))}
      </div>
      <select value={colorBy} onChange={e => setColorBy(e.target.value)}
        className="text-xs border border-gray-200 rounded-md px-1.5 py-1 bg-white text-gray-700">
        {Object.entries(COLOR_SCHEMES).map(([k, v]) => <option key={k} value={k}>{v.label}</option>)}
      </select>
      <select value={sizeBy} onChange={e => setSizeBy(e.target.value)}
        className="text-xs border border-gray-200 rounded-md px-1.5 py-1 bg-white text-gray-700">
        <option value="none">Size: off</option>
        {ALL_OBJECTIVES.map(o => <option key={o.key} value={o.key}>Size: {o.label}</option>)}
      </select>
      <label className="flex items-center gap-1 text-xs text-gray-500 cursor-pointer">
        <input type="checkbox" checked={showDominated} onChange={e => setShowDominated(e.target.checked)}
          className="accent-[#0f131d]" />
        Dominated
      </label>
      <div className="ml-auto flex items-center gap-2">
        {zoom > 1 && (
          <button onClick={handleResetZoom} className="text-xs text-blue-600 hover:underline">Reset zoom</button>
        )}
        <button onClick={exportCSV}
          className="flex items-center gap-1 px-2 py-1 text-xs font-medium text-gray-600 bg-gray-100 hover:bg-gray-200 rounded-lg transition-colors">
          <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 10v6m0 0l-3-3m3 3l3-3m2 8H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
          </svg>
          CSV
        </button>
        <button onClick={() => setShowStats(v => !v)}
          className={`px-2 py-1 text-xs font-medium rounded-lg transition-colors ${showStats ? 'bg-[#0f131d] text-white' : 'bg-gray-100 text-gray-600 hover:bg-gray-200'}`}>
          Stats
        </button>
        {paretoCount > 0 && (
          <span className="text-xs text-gray-400">
            <span className="font-semibold text-[#00e6a0]">{paretoCount}</span> on front
          </span>
        )}
        {selected.length > 0 && (
          <span className="text-xs text-blue-600 font-medium">{selected.length} selected</span>
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
    <div className="bg-white rounded-xl border border-gray-100 shadow-sm p-4">
      <h4 className="text-xs font-semibold text-gray-400 uppercase tracking-wider mb-3">Statistics</h4>
      <div className="grid grid-cols-2 sm:grid-cols-4 gap-4 text-xs">
        <div>
          <p className="text-gray-400 mb-1">Points</p>
          <p className="font-bold text-[#0f131d]">{stats.n}</p>
        </div>
        <div>
          <p className="text-gray-400 mb-1">Correlation (r)</p>
          <p className="font-bold text-[#0f131d]">{stats.correlation}</p>
        </div>
        <div>
          <p className="text-gray-400 mb-1">Pareto front</p>
          <p className="font-bold text-[#00e6a0]">{stats.paretoCount} ({stats.paretoPercent}%)</p>
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
