import React, { useEffect, useMemo, useState, useRef, useCallback } from 'react'
import { useParams, Link, useNavigate } from 'react-router-dom'
import { useWorkspace } from '../contexts/WorkspaceContext.jsx'
import { ALL_COLUMNS, GROUP_META, PHASE_TYPES, flattenMoleculeProperties, detectAvailableColumns } from '../lib/columns.js'
import { v9GetTsne } from '../api.js'
import {
  ResponsiveContainer, BarChart, Bar, XAxis, YAxis, Tooltip, CartesianGrid,
  ScatterChart, Scatter, ZAxis, Cell,
  RadarChart, Radar, PolarGrid, PolarAngleAxis, PolarRadiusAxis, Legend,
} from 'recharts'

// ──────────────────────────────────────────────────────────────────────────
// Constants
// ──────────────────────────────────────────────────────────────────────────

const PALETTE = ['#3b82f6', '#ef4444', '#10b981', '#f59e0b', '#8b5cf6', '#ec4899', '#06b6d4', '#f97316']
const HISTOGRAM_KEYS = ['docking_score', 'composite_score', 'MW', 'logP', 'QED', 'TPSA']
const RADAR_AXES = [
  { key: 'docking_score', label: 'Docking', invert: true },
  { key: 'QED', label: 'QED', invert: false },
  { key: 'logP', label: 'LogP', invert: true },
  { key: 'TPSA', label: 'TPSA', invert: true },
  { key: 'solubility', label: 'Solubility', invert: false },
  { key: 'composite_score', label: 'Composite', invert: false },
]

const COL_MAP = Object.fromEntries(ALL_COLUMNS.map(c => [c.key, c]))

// ──────────────────────────────────────────────────────────────────────────
// SMILES 2D drawing helper
// ──────────────────────────────────────────────────────────────────────────

function SmilesCanvas({ smiles, width = 100, height = 80 }) {
  const canvasRef = useRef(null)
  useEffect(() => {
    if (!smiles || !canvasRef.current) return
    let cancelled = false
    import('smiles-drawer').then(mod => {
      if (cancelled) return
      const SmilesDrawer = mod.default || mod
      const drawer = new SmilesDrawer.Drawer({ width, height, bondThickness: 1 })
      SmilesDrawer.parse(smiles, (tree) => {
        if (!cancelled) drawer.draw(tree, canvasRef.current, 'light')
      }, () => {})
    }).catch(() => {})
    return () => { cancelled = true }
  }, [smiles, width, height])
  return <canvas ref={canvasRef} width={width} height={height} className="bg-white rounded" />
}

// ──────────────────────────────────────────────────────────────────────────
// Histogram helper
// ──────────────────────────────────────────────────────────────────────────

function buildHistogramData(molecules, key, bins = 20) {
  const values = molecules.map(m => m[key]).filter(v => v != null && typeof v === 'number')
  if (values.length < 2) return null
  const min = Math.min(...values)
  const max = Math.max(...values)
  if (min === max) return null
  const binWidth = (max - min) / bins
  const counts = Array(bins).fill(0)
  for (const v of values) {
    const idx = Math.min(Math.floor((v - min) / binWidth), bins - 1)
    counts[idx]++
  }
  return counts.map((count, i) => ({
    range: `${(min + i * binWidth).toFixed(1)}`,
    count,
  }))
}

// ──────────────────────────────────────────────────────────────────────────
// Correlation helper
// ──────────────────────────────────────────────────────────────────────────

function computeCorrelationMatrix(molecules, keys) {
  const matrix = []
  for (const k1 of keys) {
    const row = {}
    row._key = k1
    for (const k2 of keys) {
      const pairs = molecules
        .map(m => [m[k1], m[k2]])
        .filter(([a, b]) => a != null && b != null && typeof a === 'number' && typeof b === 'number')
      if (pairs.length < 3) { row[k2] = null; continue }
      const n = pairs.length
      const meanX = pairs.reduce((s, p) => s + p[0], 0) / n
      const meanY = pairs.reduce((s, p) => s + p[1], 0) / n
      let num = 0, denX = 0, denY = 0
      for (const [x, y] of pairs) {
        const dx = x - meanX, dy = y - meanY
        num += dx * dy; denX += dx * dx; denY += dy * dy
      }
      const den = Math.sqrt(denX * denY)
      row[k2] = den === 0 ? 0 : num / den
    }
    matrix.push(row)
  }
  return matrix
}

// ──────────────────────────────────────────────────────────────────────────
// Section Card wrapper
// ──────────────────────────────────────────────────────────────────────────

function SectionCard({ title, subtitle, fullWidth, children }) {
  return (
    <div className={`card overflow-hidden ${fullWidth ? 'col-span-2' : ''}`}>
      <div className="px-5 py-3 border-b border-gray-100">
        <h3 className="text-sm font-semibold text-gray-700">{title}</h3>
        {subtitle && <p className="text-xs text-gray-400 mt-0.5">{subtitle}</p>}
      </div>
      <div className="p-5">{children}</div>
    </div>
  )
}

// ──────────────────────────────────────────────────────────────────────────
// Main Component
// ──────────────────────────────────────────────────────────────────────────

export default function AnalyticsPage() {
  const { projectId, phaseId } = useParams()
  const navigate = useNavigate()
  const {
    currentProject,
    currentCampaign,
    currentPhase,
    phaseMolecules,
    selectProject,
    selectPhase,
  } = useWorkspace()

  useEffect(() => { if (projectId) selectProject(projectId) }, [projectId]) // eslint-disable-line react-hooks/exhaustive-deps
  useEffect(() => { if (phaseId) selectPhase(phaseId) }, [phaseId]) // eslint-disable-line react-hooks/exhaustive-deps

  const flatMolecules = useMemo(() => phaseMolecules.map(flattenMoleculeProperties), [phaseMolecules])

  const numericKeys = useMemo(() => {
    const keys = new Set()
    for (const col of ALL_COLUMNS) {
      if (col.type !== 'number') continue
      if (flatMolecules.some(m => m[col.key] != null && typeof m[col.key] === 'number')) {
        keys.add(col.key)
      }
    }
    return [...keys]
  }, [flatMolecules])

  const phaseTypeMeta = currentPhase ? (PHASE_TYPES[currentPhase.type] || {}) : {}

  // Scatter plot state
  const [scatterX, setScatterX] = useState('docking_score')
  const [scatterY, setScatterY] = useState('composite_score')

  // Radar molecule selection
  const [radarIds, setRadarIds] = useState([])
  const bookmarkedMols = useMemo(() => flatMolecules.filter(m => m.bookmarked), [flatMolecules])

  // t-SNE state
  const [tsneData, setTsneData] = useState(null)
  const [tsneLoading, setTsneLoading] = useState(false)

  // Load t-SNE
  useEffect(() => {
    if (!phaseId) return
    setTsneLoading(true)
    v9GetTsne(phaseId)
      .then(data => setTsneData(data))
      .catch(() => setTsneData(null))
      .finally(() => setTsneLoading(false))
  }, [phaseId])

  // Auto-select radar molecules from bookmarked
  useEffect(() => {
    if (radarIds.length === 0 && bookmarkedMols.length >= 2) {
      setRadarIds(bookmarkedMols.slice(0, 6).map(m => m.id))
    }
  }, [bookmarkedMols]) // eslint-disable-line react-hooks/exhaustive-deps

  // ─── KPI Computations ───
  const kpis = useMemo(() => {
    const total = flatMolecules.length
    const bm = flatMolecules.filter(m => m.bookmarked).length
    const scores = flatMolecules.map(m => m.composite_score ?? m.weighted_score).filter(v => v != null)
    const avgScore = scores.length ? scores.reduce((s, v) => s + v, 0) / scores.length : null
    const lipinski = flatMolecules.filter(m => m.lipinski_pass === true).length
    const lipTotal = flatMolecules.filter(m => m.lipinski_pass != null).length
    const safetyGreen = flatMolecules.filter(m => m.safety_color_code === 'green').length
    const safetyTotal = flatMolecules.filter(m => m.safety_color_code != null).length
    return {
      total, bookmarked: bm,
      avgScore: avgScore?.toFixed(3),
      lipinskiPct: lipTotal ? (lipinski / lipTotal * 100).toFixed(1) : null,
      safetyGreenPct: safetyTotal ? (safetyGreen / safetyTotal * 100).toFixed(1) : null,
    }
  }, [flatMolecules])

  // ─── Scatter data ───
  const scatterData = useMemo(() => {
    return flatMolecules
      .filter(m => m[scatterX] != null && m[scatterY] != null)
      .map(m => ({
        x: m[scatterX],
        y: m[scatterY],
        name: m.name || m.id,
        cluster: m.cluster_id,
        score: m.composite_score,
        smiles: m.smiles,
        id: m.id,
      }))
  }, [flatMolecules, scatterX, scatterY])

  // ─── Radar data ───
  const radarData = useMemo(() => {
    const mols = flatMolecules.filter(m => radarIds.includes(m.id))
    if (mols.length < 2) return null
    const axes = RADAR_AXES.filter(a => mols.some(m => m[a.key] != null))
    if (axes.length < 3) return null

    const ranges = {}
    for (const a of axes) {
      const vals = mols.map(m => m[a.key]).filter(v => v != null)
      ranges[a.key] = { min: Math.min(...vals), max: Math.max(...vals) }
    }

    return axes.map(axis => {
      const r = ranges[axis.key]
      const span = r.max - r.min || 1
      const pt = { axis: axis.label }
      mols.forEach((m, i) => {
        const raw = m[axis.key]
        let norm = raw != null ? (raw - r.min) / span : 0
        if (axis.invert) norm = 1 - norm
        pt[`mol_${i}`] = Math.round(norm * 100) / 100
      })
      return pt
    })
  }, [flatMolecules, radarIds])

  const radarMols = useMemo(() => flatMolecules.filter(m => radarIds.includes(m.id)), [flatMolecules, radarIds])

  // ─── Correlation matrix ───
  const corrKeys = useMemo(() => numericKeys.slice(0, 12), [numericKeys])
  const corrMatrix = useMemo(() => {
    if (corrKeys.length < 3) return null
    return computeCorrelationMatrix(flatMolecules, corrKeys)
  }, [flatMolecules, corrKeys])

  // ─── Activity cliffs ───
  const cliffMols = useMemo(() => {
    return flatMolecules
      .filter(m => m.is_cliff === true && m.sali_max != null)
      .sort((a, b) => (b.sali_max || 0) - (a.sali_max || 0))
      .slice(0, 20)
  }, [flatMolecules])

  // ─── t-SNE points with molecule data ───
  const tsnePoints = useMemo(() => {
    if (!tsneData?.points?.length) return null
    const molMap = new Map(flatMolecules.map(m => [String(m.id), m]))
    return tsneData.points.map(p => {
      const mol = molMap.get(p.molecule_id)
      return {
        x: p.x,
        y: p.y,
        name: mol?.name || p.molecule_id.slice(0, 8),
        score: mol?.composite_score,
        cluster: mol?.cluster_id,
      }
    })
  }, [tsneData, flatMolecules])

  if (!currentPhase) {
    return (
      <div className="flex items-center justify-center h-64">
        <div className="w-8 h-8 border-2 border-bx-accent border-t-transparent rounded-full animate-spin" />
      </div>
    )
  }

  return (
    <div className="space-y-4 pb-8">
      {/* Breadcrumb */}
      <nav className="flex items-center gap-1.5 text-sm text-gray-400">
        <Link to="/" className="hover:text-bx-mint transition-colors">Projects</Link>
        <ChevronIcon />
        <Link to={`/project/${projectId}`} className="hover:text-bx-light-text transition-colors truncate max-w-[120px]">
          {currentProject?.name || 'Project'}
        </Link>
        <ChevronIcon />
        <Link to={`/project/${projectId}/phase/${phaseId}`} className="hover:text-bx-light-text transition-colors">
          Phase {phaseTypeMeta.short || '?'}
        </Link>
        <ChevronIcon />
        <span className="text-gray-600 font-medium">Analytics</span>
      </nav>

      {/* Header */}
      <div className="card px-5 py-4 flex items-center justify-between">
        <div>
          <h1 className="text-xl font-bold text-bx-light-text">Analytics</h1>
          <p className="text-sm text-gray-400 mt-0.5">
            Phase {phaseTypeMeta.short} — {phaseTypeMeta.label} — {kpis.total} molecules
          </p>
        </div>
        <button
          onClick={() => navigate(`/project/${projectId}/phase/${phaseId}`)}
          className="flex items-center gap-1.5 px-3 py-2 rounded-lg text-sm font-medium border border-gray-200 text-gray-600 hover:bg-gray-50 transition-colors"
        >
          <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10 19l-7-7m0 0l7-7m-7 7h18" />
          </svg>
          Dashboard
        </button>
      </div>

      {/* (a) KPI Cards */}
      <div className="grid grid-cols-5 gap-3">
        <KpiCard label="Total" value={kpis.total} color="text-bx-mint" />
        <KpiCard label="Bookmarked" value={kpis.bookmarked} color="text-amber-500" />
        <KpiCard label="Avg Score" value={kpis.avgScore || '—'} color="text-blue-500" />
        <KpiCard label="Lipinski Pass" value={kpis.lipinskiPct ? `${kpis.lipinskiPct}%` : '—'} color="text-emerald-500" />
        <KpiCard label="Safety Green" value={kpis.safetyGreenPct ? `${kpis.safetyGreenPct}%` : '—'} color="text-green-500" />
      </div>

      {/* Chart Grid */}
      <div className="grid grid-cols-1 md:grid-cols-2 gap-4">

        {/* (b) Histograms */}
        {HISTOGRAM_KEYS.map(key => {
          const data = buildHistogramData(flatMolecules, key)
          if (!data) return null
          const colDef = COL_MAP[key]
          const group = colDef?.group
          const color = group === 'docking' ? '#3b82f6' : group === 'scoring' ? '#8b5cf6' : group === 'adme' ? '#10b981' : group === 'toxicity' ? '#ef4444' : '#64748b'
          return (
            <SectionCard key={key} title={`${colDef?.label || key} Distribution`}>
              <ResponsiveContainer width="100%" height={180}>
                <BarChart data={data}>
                  <CartesianGrid strokeDasharray="3 3" stroke="#f1f5f9" />
                  <XAxis dataKey="range" tick={{ fontSize: 9 }} interval={Math.floor(data.length / 5)} />
                  <YAxis tick={{ fontSize: 9 }} />
                  <Tooltip formatter={(v) => [v, 'Count']} labelFormatter={(l) => `Range: ${l}`} />
                  <Bar dataKey="count" fill={color} radius={[2, 2, 0, 0]} />
                </BarChart>
              </ResponsiveContainer>
            </SectionCard>
          )
        })}

        {/* (c) Scatter Plot Explorer */}
        <SectionCard title="Scatter Plot Explorer" subtitle="Select X/Y axes" fullWidth>
          <div className="flex items-center gap-3 mb-3">
            <label className="text-xs text-gray-500">X:</label>
            <select value={scatterX} onChange={e => setScatterX(e.target.value)} className="text-xs border border-gray-200 rounded px-2 py-1">
              {numericKeys.map(k => <option key={k} value={k}>{COL_MAP[k]?.label || k}</option>)}
            </select>
            <label className="text-xs text-gray-500">Y:</label>
            <select value={scatterY} onChange={e => setScatterY(e.target.value)} className="text-xs border border-gray-200 rounded px-2 py-1">
              {numericKeys.map(k => <option key={k} value={k}>{COL_MAP[k]?.label || k}</option>)}
            </select>
          </div>
          <ResponsiveContainer width="100%" height={320}>
            <ScatterChart>
              <CartesianGrid strokeDasharray="3 3" stroke="#f1f5f9" />
              <XAxis type="number" dataKey="x" name={COL_MAP[scatterX]?.label || scatterX} tick={{ fontSize: 10 }} />
              <YAxis type="number" dataKey="y" name={COL_MAP[scatterY]?.label || scatterY} tick={{ fontSize: 10 }} />
              <ZAxis range={[30, 30]} />
              <Tooltip
                content={({ payload }) => {
                  if (!payload?.[0]) return null
                  const d = payload[0].payload
                  return (
                    <div className="bg-white shadow-lg rounded-lg p-2 border text-xs">
                      <p className="font-semibold">{d.name}</p>
                      <p>{COL_MAP[scatterX]?.label}: {d.x?.toFixed(2)}</p>
                      <p>{COL_MAP[scatterY]?.label}: {d.y?.toFixed(2)}</p>
                      {d.cluster != null && <p>Cluster: {d.cluster}</p>}
                    </div>
                  )
                }}
              />
              <Scatter data={scatterData} fill="#3b82f6" fillOpacity={0.6}>
                {scatterData.map((entry, i) => (
                  <Cell
                    key={i}
                    fill={entry.cluster != null ? PALETTE[entry.cluster % PALETTE.length] : '#3b82f6'}
                  />
                ))}
              </Scatter>
            </ScatterChart>
          </ResponsiveContainer>
        </SectionCard>

        {/* (d) Radar Chart Comparison */}
        {bookmarkedMols.length >= 2 && (
          <SectionCard title="Radar Comparison" subtitle="Bookmarked molecules">
            <div className="flex flex-wrap gap-1 mb-3">
              {bookmarkedMols.slice(0, 10).map(m => {
                const selected = radarIds.includes(m.id)
                return (
                  <button
                    key={m.id}
                    onClick={() => {
                      setRadarIds(prev => {
                        if (selected) return prev.filter(id => id !== m.id)
                        if (prev.length >= 6) return prev
                        return [...prev, m.id]
                      })
                    }}
                    className={`px-2 py-0.5 text-[10px] rounded-full border transition-colors ${
                      selected ? 'bg-blue-50 border-blue-300 text-blue-700' : 'border-gray-200 text-gray-500 hover:border-gray-300'
                    }`}
                  >
                    {m.name || String(m.id).slice(0, 8)}
                  </button>
                )
              })}
            </div>
            {radarData && radarData.length > 2 ? (
              <ResponsiveContainer width="100%" height={280}>
                <RadarChart data={radarData}>
                  <PolarGrid stroke="#e2e8f0" />
                  <PolarAngleAxis dataKey="axis" tick={{ fontSize: 10, fill: '#64748b' }} />
                  <PolarRadiusAxis angle={90} domain={[0, 1]} tick={{ fontSize: 8 }} />
                  {radarMols.map((mol, i) => (
                    <Radar
                      key={mol.id}
                      name={mol.name || `Mol ${i + 1}`}
                      dataKey={`mol_${i}`}
                      stroke={PALETTE[i % PALETTE.length]}
                      fill={PALETTE[i % PALETTE.length]}
                      fillOpacity={0.08}
                      strokeWidth={2}
                    />
                  ))}
                  <Legend wrapperStyle={{ fontSize: 10 }} />
                  <Tooltip />
                </RadarChart>
              </ResponsiveContainer>
            ) : (
              <p className="text-xs text-gray-400 text-center py-8">Select 2+ bookmarked molecules</p>
            )}
          </SectionCard>
        )}

        {/* (e) Correlation Heatmap */}
        {corrMatrix && (
          <SectionCard title="Correlation Heatmap" subtitle="Pearson coefficients">
            <div className="overflow-x-auto">
              <table className="text-[9px]">
                <thead>
                  <tr>
                    <th className="px-1 py-0.5" />
                    {corrKeys.map(k => (
                      <th key={k} className="px-1 py-0.5 font-medium text-gray-500 whitespace-nowrap" style={{ writingMode: 'vertical-rl', transform: 'rotate(180deg)' }}>
                        {COL_MAP[k]?.label || k}
                      </th>
                    ))}
                  </tr>
                </thead>
                <tbody>
                  {corrMatrix.map((row) => (
                    <tr key={row._key}>
                      <td className="px-1 py-0.5 font-medium text-gray-500 whitespace-nowrap">{COL_MAP[row._key]?.label || row._key}</td>
                      {corrKeys.map(k => {
                        const val = row[k]
                        if (val == null) return <td key={k} className="px-1 py-0.5 bg-gray-100 text-center">—</td>
                        const intensity = Math.abs(val)
                        const bg = val > 0
                          ? `rgba(239, 68, 68, ${intensity * 0.5})`
                          : `rgba(59, 130, 246, ${intensity * 0.5})`
                        return (
                          <td key={k} className="px-1 py-0.5 text-center w-8 h-8" style={{ backgroundColor: bg }} title={`${row._key} × ${k}: ${val.toFixed(2)}`}>
                            {val.toFixed(1)}
                          </td>
                        )
                      })}
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </SectionCard>
        )}

        {/* (f) Activity Cliffs */}
        {cliffMols.length > 0 && (
          <SectionCard title="Activity Cliffs" subtitle={`${cliffMols.length} molecules with activity cliffs`} fullWidth>
            <div className="overflow-x-auto">
              <table className="w-full text-xs">
                <thead>
                  <tr className="bg-gray-50">
                    <th className="px-3 py-2 text-left">Name</th>
                    <th className="px-3 py-2 text-left">SMILES</th>
                    <th className="px-3 py-2 text-center">SALI</th>
                    <th className="px-3 py-2 text-center">Score</th>
                    <th className="px-3 py-2 text-center">Cliffs</th>
                  </tr>
                </thead>
                <tbody>
                  {cliffMols.map(m => (
                    <tr key={m.id} className="border-t border-gray-50 hover:bg-gray-50/50">
                      <td className="px-3 py-1.5 font-medium text-gray-700">{m.name || '—'}</td>
                      <td className="px-3 py-1.5">
                        <SmilesCanvas smiles={m.smiles} width={100} height={50} />
                      </td>
                      <td className="px-3 py-1.5 text-center font-semibold text-amber-600">{m.sali_max?.toFixed(2)}</td>
                      <td className="px-3 py-1.5 text-center tabular-nums">{(m.composite_score ?? m.docking_score)?.toFixed(2) || '—'}</td>
                      <td className="px-3 py-1.5 text-center">{m.n_cliffs ?? '—'}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </SectionCard>
        )}

        {/* (g) Chemical Space (t-SNE) */}
        {(tsneLoading || tsnePoints) && (
          <SectionCard title="Chemical Space (t-SNE)" subtitle="Morgan fingerprint 2D projection" fullWidth>
            {tsneLoading ? (
              <div className="flex items-center justify-center h-40 gap-2 text-gray-400">
                <div className="w-5 h-5 border-2 border-blue-400 border-t-transparent rounded-full animate-spin" />
                <span className="text-sm">Computing t-SNE...</span>
              </div>
            ) : tsnePoints && tsnePoints.length > 0 ? (
              <ResponsiveContainer width="100%" height={320}>
                <ScatterChart>
                  <CartesianGrid strokeDasharray="3 3" stroke="#f1f5f9" />
                  <XAxis type="number" dataKey="x" tick={{ fontSize: 9 }} name="t-SNE 1" />
                  <YAxis type="number" dataKey="y" tick={{ fontSize: 9 }} name="t-SNE 2" />
                  <ZAxis range={[25, 25]} />
                  <Tooltip
                    content={({ payload }) => {
                      if (!payload?.[0]) return null
                      const d = payload[0].payload
                      return (
                        <div className="bg-white shadow-lg rounded-lg p-2 border text-xs">
                          <p className="font-semibold">{d.name}</p>
                          {d.score != null && <p>Score: {d.score.toFixed(3)}</p>}
                          {d.cluster != null && <p>Cluster: {d.cluster}</p>}
                        </div>
                      )
                    }}
                  />
                  <Scatter data={tsnePoints} fillOpacity={0.7}>
                    {tsnePoints.map((p, i) => (
                      <Cell key={i} fill={p.cluster != null ? PALETTE[p.cluster % PALETTE.length] : '#3b82f6'} />
                    ))}
                  </Scatter>
                </ScatterChart>
              </ResponsiveContainer>
            ) : (
              <p className="text-xs text-gray-400 text-center py-8">t-SNE not available for this phase</p>
            )}
          </SectionCard>
        )}
      </div>
    </div>
  )
}

// ──────────────────────────────────────────────────────────────────────────
// Small helper components
// ──────────────────────────────────────────────────────────────────────────

function KpiCard({ label, value, color }) {
  return (
    <div className="card px-4 py-3 text-center">
      <p className={`text-2xl font-bold ${color}`}>{value}</p>
      <p className="text-[10px] uppercase tracking-wider text-gray-400 mt-0.5">{label}</p>
    </div>
  )
}

function ChevronIcon() {
  return (
    <svg className="w-3 h-3 text-gray-300" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
    </svg>
  )
}
