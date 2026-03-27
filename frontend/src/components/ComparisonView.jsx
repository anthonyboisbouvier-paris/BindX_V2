import React, { useMemo, useRef, useEffect, useCallback } from 'react'
import { RadarChart, Radar, PolarGrid, PolarAngleAxis, PolarRadiusAxis, Legend, ResponsiveContainer, Tooltip } from 'recharts'
import { ALL_COLUMNS, GROUP_META } from '../lib/columns.js'

// ──────────────────────────────────────────────────────────────────────────
// Constants
// ──────────────────────────────────────────────────────────────────────────

const PALETTE = ['#3b82f6', '#ef4444', '#10b981', '#f59e0b', '#8b5cf6', '#ec4899']

const RADAR_AXES = [
  { key: 'docking_score', label: 'Docking', invert: true },
  { key: 'cnn_score', label: 'CNN Score', invert: false },
  { key: 'QED', label: 'QED', invert: false },
  { key: 'logP', label: 'LogP', invert: true },
  { key: 'TPSA', label: 'TPSA', invert: true },
  { key: 'solubility', label: 'Solubility', invert: false },
  { key: 'sa_score', label: 'SA Score', invert: true },
  { key: 'composite_score', label: 'Composite', invert: false },
]

const COL_MAP = Object.fromEntries(ALL_COLUMNS.map(c => [c.key, c]))

// ──────────────────────────────────────────────────────────────────────────
// SMILES 2D drawing helper (smiles-drawer)
// ──────────────────────────────────────────────────────────────────────────

function SmilesImage({ smiles, width = 160, height = 120 }) {
  const canvasRef = useRef(null)

  useEffect(() => {
    if (!smiles || !canvasRef.current) return
    let cancelled = false

    import('smiles-drawer').then(mod => {
      if (cancelled) return
      const SmilesDrawer = mod.default || mod
      const drawer = new SmilesDrawer.Drawer({ width, height, bondThickness: 1.2 })
      SmilesDrawer.parse(smiles, (tree) => {
        if (!cancelled) drawer.draw(tree, canvasRef.current, 'light')
      }, () => {})
    }).catch(() => {})

    return () => { cancelled = true }
  }, [smiles, width, height])

  return <canvas ref={canvasRef} width={width} height={height} className="bg-white rounded" />
}

// ──────────────────────────────────────────────────────────────────────────
// Normalize for radar
// ──────────────────────────────────────────────────────────────────────────

function normalizeRadar(molecules) {
  // Find min/max for each axis
  const ranges = {}
  for (const axis of RADAR_AXES) {
    const vals = molecules.map(m => m[axis.key]).filter(v => v != null && typeof v === 'number')
    if (vals.length === 0) { ranges[axis.key] = null; continue }
    ranges[axis.key] = { min: Math.min(...vals), max: Math.max(...vals) }
  }

  return RADAR_AXES.filter(a => ranges[a.key] !== null).map(axis => {
    const r = ranges[axis.key]
    const span = r.max - r.min || 1
    const point = { axis: axis.label }
    molecules.forEach((mol, i) => {
      const raw = mol[axis.key]
      if (raw == null) { point[`mol_${i}`] = 0; return }
      let norm = (raw - r.min) / span
      if (axis.invert) norm = 1 - norm
      point[`mol_${i}`] = Math.round(norm * 100) / 100
    })
    return point
  })
}

// ──────────────────────────────────────────────────────────────────────────
// Main component
// ──────────────────────────────────────────────────────────────────────────

export default function ComparisonView({ molecules, onBack, onBookmark, onExportCSV }) {
  if (!molecules || molecules.length < 2) {
    return (
      <div className="text-center py-12 text-gray-400">
        Select 2-6 molecules to compare.
      </div>
    )
  }

  const radarData = useMemo(() => normalizeRadar(molecules), [molecules])

  // Determine which numeric columns have data across any molecule
  const numericCols = useMemo(() => {
    return ALL_COLUMNS.filter(c => {
      if (c.type !== 'number') return false
      return molecules.some(m => m[c.key] != null && typeof m[c.key] === 'number')
    })
  }, [molecules])

  // Group numeric columns by group
  const groupedCols = useMemo(() => {
    const groups = {}
    for (const col of numericCols) {
      const g = col.group || 'other'
      if (!groups[g]) groups[g] = []
      groups[g].push(col)
    }
    return groups
  }, [numericCols])

  // Find best/worst per column
  const rankings = useMemo(() => {
    const result = {}
    for (const col of numericCols) {
      const vals = molecules.map((m, i) => ({ i, v: m[col.key] })).filter(x => x.v != null && typeof x.v === 'number')
      if (vals.length < 2) continue
      const sorted = [...vals].sort((a, b) => a.v - b.v)
      const isLowerBetter = col.colorScale === 'lower-better'
      result[col.key] = {
        best: isLowerBetter ? sorted[0].i : sorted[sorted.length - 1].i,
        worst: isLowerBetter ? sorted[sorted.length - 1].i : sorted[0].i,
      }
    }
    return result
  }, [numericCols, molecules])

  // Best composite molecule
  const bestIdx = useMemo(() => {
    let best = 0
    let bestScore = -Infinity
    molecules.forEach((m, i) => {
      const s = m.composite_score ?? m.weighted_score ?? -Infinity
      if (s > bestScore) { bestScore = s; best = i }
    })
    return best
  }, [molecules])

  // Export comparison as CSV
  const handleExportCSV = useCallback(() => {
    const rows = [['Property', ...molecules.map(m => m.name || m.id)]]
    for (const col of numericCols) {
      rows.push([col.label, ...molecules.map(m => m[col.key] ?? '')])
    }
    const csv = rows.map(r => r.join(',')).join('\n')
    const blob = new Blob([csv], { type: 'text/csv' })
    const url = URL.createObjectURL(blob)
    const a = document.createElement('a')
    a.href = url
    a.download = `bindx_comparison_${Date.now()}.csv`
    a.click()
    URL.revokeObjectURL(url)
  }, [molecules, numericCols])

  return (
    <div className="space-y-6">
      {/* (a) Header Cards */}
      <div className="flex gap-4 overflow-x-auto pb-2">
        {molecules.map((mol, i) => (
          <div key={mol.id} className="card p-4 min-w-[200px] flex-shrink-0 border-t-4" style={{ borderTopColor: PALETTE[i % PALETTE.length] }}>
            <SmilesImage smiles={mol.smiles} width={160} height={120} />
            <h3 className="font-semibold text-sm text-gray-700 mt-2 truncate">{mol.name || `Mol ${i + 1}`}</h3>
            <div className="flex items-center gap-2 mt-1">
              <span className="text-2xl font-bold" style={{ color: PALETTE[i % PALETTE.length] }}>
                {(mol.composite_score ?? mol.weighted_score)?.toFixed(3) || '—'}
              </span>
              {mol.bookmarked && (
                <svg className="w-4 h-4 text-yellow-400" fill="currentColor" viewBox="0 0 24 24">
                  <path d="M5 5a2 2 0 012-2h10a2 2 0 012 2v16l-7-3.5L5 21V5z" />
                </svg>
              )}
            </div>
            <p className="text-[10px] text-gray-400 mt-0.5">Composite Score</p>
          </div>
        ))}
      </div>

      {/* (b) Radar Chart */}
      {radarData.length > 2 && (
        <div className="card p-5">
          <h3 className="text-sm font-semibold text-gray-700 mb-3">Radar Comparison</h3>
          <ResponsiveContainer width="100%" height={340}>
            <RadarChart data={radarData}>
              <PolarGrid stroke="#e2e8f0" />
              <PolarAngleAxis dataKey="axis" tick={{ fontSize: 11, fill: '#64748b' }} />
              <PolarRadiusAxis angle={90} domain={[0, 1]} tick={{ fontSize: 9 }} />
              {molecules.map((mol, i) => (
                <Radar
                  key={mol.id}
                  name={mol.name || `Mol ${i + 1}`}
                  dataKey={`mol_${i}`}
                  stroke={PALETTE[i % PALETTE.length]}
                  fill={PALETTE[i % PALETTE.length]}
                  fillOpacity={0.1}
                  strokeWidth={2}
                />
              ))}
              <Legend wrapperStyle={{ fontSize: 11 }} />
              <Tooltip />
            </RadarChart>
          </ResponsiveContainer>
        </div>
      )}

      {/* (c) Property Comparison Table */}
      <div className="card overflow-hidden">
        <h3 className="text-sm font-semibold text-gray-700 px-5 py-3 border-b border-gray-100">
          Property Comparison
        </h3>
        <div className="overflow-x-auto">
          <table className="w-full text-xs">
            <thead>
              <tr className="bg-gray-50">
                <th className="sticky left-0 bg-gray-50 z-10 text-left px-4 py-2 font-semibold text-gray-500">Property</th>
                {molecules.map((mol, i) => (
                  <th key={mol.id} className="px-3 py-2 text-center font-semibold" style={{ color: PALETTE[i % PALETTE.length] }}>
                    {mol.name || `Mol ${i + 1}`}
                  </th>
                ))}
              </tr>
            </thead>
            <tbody>
              {Object.entries(groupedCols).map(([group, cols]) => (
                <React.Fragment key={group}>
                  <tr>
                    <td colSpan={molecules.length + 1} className={`px-4 py-1.5 text-[10px] font-bold uppercase tracking-wider ${GROUP_META[group]?.text || 'text-gray-400'} ${GROUP_META[group]?.bg || 'bg-gray-50'}`}>
                      {GROUP_META[group]?.label || group}
                    </td>
                  </tr>
                  {cols.map(col => (
                    <tr key={col.key} className="border-t border-gray-50 hover:bg-gray-50/50">
                      <td className="sticky left-0 bg-white z-10 px-4 py-1.5 text-gray-600 font-medium whitespace-nowrap">
                        {col.label} {col.unit ? <span className="text-gray-300">({col.unit})</span> : null}
                      </td>
                      {molecules.map((mol, i) => {
                        const val = mol[col.key]
                        const rank = rankings[col.key]
                        const isBest = rank?.best === i
                        const isWorst = rank?.worst === i
                        return (
                          <td key={mol.id} className={`px-3 py-1.5 text-center tabular-nums ${isBest ? 'bg-emerald-50 text-emerald-700 font-semibold' : isWorst ? 'bg-red-50 text-red-600' : 'text-gray-700'}`}>
                            {val == null ? '—' : typeof val === 'number' ? val.toFixed(2) : typeof val === 'boolean' ? (val ? 'Yes' : 'No') : String(val)}
                          </td>
                        )
                      })}
                    </tr>
                  ))}
                </React.Fragment>
              ))}
            </tbody>
          </table>
        </div>
      </div>

      {/* (d) ADMET Profile Comparison */}
      {molecules.some(m => m.safety_color_code != null) && (
        <div className="card p-5">
          <h3 className="text-sm font-semibold text-gray-700 mb-3">ADMET Profiles</h3>
          <div className="grid grid-cols-2 md:grid-cols-3 lg:grid-cols-6 gap-3">
            {molecules.map((mol, i) => (
              <div key={mol.id} className="border rounded-lg p-3" style={{ borderColor: PALETTE[i % PALETTE.length] + '40' }}>
                <p className="text-xs font-semibold truncate" style={{ color: PALETTE[i % PALETTE.length] }}>{mol.name || `Mol ${i + 1}`}</p>
                <div className="flex items-center gap-1.5 mt-2">
                  <span className={`w-3 h-3 rounded-full ${
                    mol.safety_color_code === 'green' ? 'bg-emerald-400' :
                    mol.safety_color_code === 'yellow' ? 'bg-amber-400' :
                    mol.safety_color_code === 'red' ? 'bg-red-400' : 'bg-gray-300'
                  }`} />
                  <span className="text-[10px] text-gray-500">{mol.safety_color_code || '—'}</span>
                </div>
                <div className="mt-2 space-y-1 text-[10px] text-gray-500">
                  {mol.pains_alert != null && <p>PAINS: {mol.pains_alert ? 'Alert' : 'Clear'}</p>}
                  {mol.brenk_alert != null && <p>Brenk: {mol.brenk_alert ? 'Alert' : 'Clear'}</p>}
                  {mol.cyp_inhibitions != null && <p>CYP inh: {mol.cyp_inhibitions}/5</p>}
                  {mol.hERG != null && <p>hERG: {mol.hERG.toFixed(2)}</p>}
                </div>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* (e) Synthesis Feasibility */}
      {molecules.some(m => m.n_synth_steps != null) && (
        <div className="card p-5">
          <h3 className="text-sm font-semibold text-gray-700 mb-3">Synthesis Feasibility</h3>
          <div className="space-y-3">
            {['n_synth_steps', 'synth_confidence', 'synth_cost_estimate'].map(key => {
              const colDef = COL_MAP[key]
              if (!colDef) return null
              const vals = molecules.map(m => m[key]).filter(v => v != null)
              if (vals.length === 0) return null
              const maxVal = typeof vals[0] === 'number' ? Math.max(...vals.filter(v => typeof v === 'number'), 1) : null
              return (
                <div key={key}>
                  <p className="text-xs text-gray-500 mb-1">{colDef.label}</p>
                  <div className="space-y-1">
                    {molecules.map((mol, i) => {
                      const v = mol[key]
                      const pct = maxVal && typeof v === 'number' ? (v / maxVal) * 100 : 0
                      return (
                        <div key={mol.id} className="flex items-center gap-2">
                          <span className="text-[10px] w-16 truncate" style={{ color: PALETTE[i % PALETTE.length] }}>
                            {mol.name || `Mol ${i + 1}`}
                          </span>
                          {typeof v === 'number' ? (
                            <>
                              <div className="flex-1 h-3 bg-gray-100 rounded-full overflow-hidden">
                                <div className="h-full rounded-full" style={{ width: `${pct}%`, backgroundColor: PALETTE[i % PALETTE.length] }} />
                              </div>
                              <span className="text-xs text-gray-600 w-12 text-right tabular-nums">{v.toFixed(key === 'n_synth_steps' ? 0 : 2)}</span>
                            </>
                          ) : (
                            <span className="text-xs text-gray-500">{v ?? '—'}</span>
                          )}
                        </div>
                      )
                    })}
                  </div>
                </div>
              )
            })}
          </div>
        </div>
      )}

      {/* Action buttons */}
      <div className="flex items-center gap-3">
        {onBack && (
          <button onClick={onBack} className="flex items-center gap-1.5 px-4 py-2 rounded-lg text-sm font-medium border border-gray-200 text-gray-600 hover:bg-gray-50 transition-colors">
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10 19l-7-7m0 0l7-7m-7 7h18" />
            </svg>
            Back to Dashboard
          </button>
        )}
        {onBookmark && (
          <button onClick={() => onBookmark(molecules[bestIdx]?.id)} className="flex items-center gap-1.5 px-4 py-2 rounded-lg text-sm font-medium border border-yellow-300 text-yellow-700 hover:bg-yellow-50 transition-colors">
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M5 5a2 2 0 012-2h10a2 2 0 012 2v16l-7-3.5L5 21V5z" />
            </svg>
            Bookmark Best
          </button>
        )}
        <button onClick={handleExportCSV} className="flex items-center gap-1.5 px-4 py-2 rounded-lg text-sm font-medium border border-gray-200 text-gray-600 hover:bg-gray-50 transition-colors">
          <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 10v6m0 0l-3-3m3 3l3-3m2 8H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
          </svg>
          Export CSV
        </button>
      </div>
    </div>
  )
}
