import React, { useEffect, useState, useMemo, useCallback, useRef } from 'react'
import { useParams, Link } from 'react-router-dom'
import { useWorkspace } from '../contexts/WorkspaceContext.jsx'
import { PHASE_TYPES } from '../lib/columns.js'
import { v9GetSARAnalysis, v9SmilesToSvgDiff } from '../api.js'
import BindXLogo from '../components/BindXLogo.jsx'

// ──────────────────────────────────────────────────────────────────────────
// Helpers
// ──────────────────────────────────────────────────────────────────────────

function fmt(v, d = 2) { return v == null ? '—' : v.toFixed(d) }

function ChevronRight() {
  return <svg className="w-3 h-3 text-gray-300" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" /></svg>
}

// ──────────────────────────────────────────────────────────────────────────
// SALI badge
// ──────────────────────────────────────────────────────────────────────────

function SaliBadge({ value }) {
  if (value == null) return <span className="text-gray-300">—</span>
  let bg
  if (value >= 0.8) { bg = 'bg-orange-100 text-orange-700' }
  else if (value >= 0.5) { bg = 'bg-amber-50 text-amber-600' }
  else { bg = 'bg-gray-100 text-gray-500' }
  return (
    <span className={`inline-block px-2 py-0.5 rounded-full text-[11px] font-semibold ${bg}`}>
      {value.toFixed(1)}
    </span>
  )
}

// ──────────────────────────────────────────────────────────────────────────
// CliffPairRow — loads diff SVGs lazily via IntersectionObserver
// ──────────────────────────────────────────────────────────────────────────

const SVG_W = 150
const SVG_H = 100

function CliffPairRow({ pair, selectedProps, moleculeProperties }) {
  const [svgA, setSvgA] = useState(null)
  const [svgB, setSvgB] = useState(null)
  const [visible, setVisible] = useState(false)
  const ref = useRef(null)

  useEffect(() => {
    if (!ref.current) return
    const observer = new IntersectionObserver(
      ([entry]) => { if (entry.isIntersecting) { setVisible(true); observer.disconnect() } },
      { rootMargin: '200px' }
    )
    observer.observe(ref.current)
    return () => observer.disconnect()
  }, [])

  useEffect(() => {
    if (!visible || !pair.mol_a_smiles || !pair.mol_b_smiles) return
    v9SmilesToSvgDiff(pair.mol_a_smiles, pair.mol_b_smiles, SVG_W, SVG_H)
      .then(({ svg_a, svg_b }) => { setSvgA(svg_a); setSvgB(svg_b) })
      .catch(() => { setSvgA(null); setSvgB(null) })
  }, [visible, pair.mol_a_smiles, pair.mol_b_smiles])

  const labelA = pair.mol_a_name || pair.mol_a_id?.slice(0, 8)
  const labelB = pair.mol_b_name || pair.mol_b_id?.slice(0, 8)

  // Look up molecule properties for multi-score deltas
  const propsA = moleculeProperties?.[pair.mol_a_id] || {}
  const propsB = moleculeProperties?.[pair.mol_b_id] || {}

  const svgPlaceholder = (isVisible) => (
    <div className="flex items-center justify-center rounded bg-gray-50 border border-gray-100"
      style={{ width: SVG_W, height: SVG_H }}>
      {isVisible ? (
        <div className="w-4 h-4 border-2 border-gray-300 border-t-transparent rounded-full animate-spin" />
      ) : (
        <span className="text-[9px] text-gray-300">...</span>
      )}
    </div>
  )

  return (
    <tr ref={ref} className="border-b border-gray-100 hover:bg-gray-50/50">
      <td className="px-4 py-3">
        <div className="flex flex-col items-center gap-0.5">
          {svgA ? (
            <div className="rounded bg-white border border-gray-100" style={{ width: SVG_W, height: SVG_H, overflow: 'hidden' }}
              dangerouslySetInnerHTML={{ __html: svgA }} />
          ) : svgPlaceholder(visible)}
          {labelA && <div className="text-[11px] text-gray-400 truncate max-w-[140px]" title={labelA}>{labelA}</div>}
        </div>
      </td>
      <td className="px-4 py-3">
        <div className="flex flex-col items-center gap-0.5">
          {svgB ? (
            <div className="rounded bg-white border border-gray-100" style={{ width: SVG_W, height: SVG_H, overflow: 'hidden' }}
              dangerouslySetInnerHTML={{ __html: svgB }} />
          ) : svgPlaceholder(visible)}
          {labelB && <div className="text-[11px] text-gray-400 truncate max-w-[140px]" title={labelB}>{labelB}</div>}
        </div>
      </td>
      <td className="px-4 py-3 text-right font-mono text-xs text-gray-600">
        {fmt(pair.similarity, 2)}
      </td>
      {selectedProps.map(sp => {
        const valA = propsA[sp.key] ?? pair.activity_a
        const valB = propsB[sp.key] ?? pair.activity_b
        const delta = (valA != null && valB != null) ? valA - valB : null
        const color = delta == null ? 'text-gray-300' : delta >= 0 ? 'text-green-600' : 'text-red-600'
        return (
          <td key={sp.key} className={`px-3 py-3 text-right font-semibold text-xs ${color}`}>
            {delta == null ? '—' : `${delta >= 0 ? '+' : ''}${fmt(delta, 3)}`}
          </td>
        )
      })}
      <td className="px-4 py-3 text-right">
        <SaliBadge value={pair.sali} />
      </td>
    </tr>
  )
}

// ──────────────────────────────────────────────────────────────────────────
// Multi-select score chips
// ──────────────────────────────────────────────────────────────────────────

function ScoreSelector({ properties, selected, onToggle }) {
  return (
    <div className="flex flex-wrap items-center gap-1.5">
      {properties.map(p => {
        const active = selected.includes(p.key)
        return (
          <button
            key={p.key}
            onClick={() => onToggle(p.key)}
            className={`px-2.5 py-1 rounded-full text-[11px] font-medium border transition-all ${
              active
                ? 'bg-bx-surface text-white border-bx-surface'
                : 'bg-white text-gray-500 border-gray-200 hover:border-gray-400'
            }`}
          >
            {p.label}
          </button>
        )
      })}
    </div>
  )
}

// ──────────────────────────────────────────────────────────────────────────
// Main page
// ──────────────────────────────────────────────────────────────────────────

export default function ActivityExplorerPage() {
  const { projectId, phaseId } = useParams()
  const { currentProject, currentPhase, selectProject, selectPhase } = useWorkspace()

  useEffect(() => { if (projectId) selectProject(projectId) }, [projectId]) // eslint-disable-line
  useEffect(() => { if (phaseId) selectPhase(phaseId) }, [phaseId]) // eslint-disable-line

  const [sarData, setSarData] = useState(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState(null)

  // Controls
  const [propertyKey, setPropertyKey] = useState('auto')
  const [selectedKeys, setSelectedKeys] = useState([])
  const [saliMin, setSaliMin] = useState(0.5)

  const fetchData = useCallback((propKey = 'auto') => {
    if (!phaseId) return
    setLoading(true)
    setError(null)
    v9GetSARAnalysis(phaseId, propKey)
      .then(result => {
        setSarData(result)
        if (propKey === 'auto' && result.property_key) {
          setPropertyKey(result.property_key)
          // Auto-select the primary property if nothing selected yet
          setSelectedKeys(prev => prev.length ? prev : [result.property_key])
        }
      })
      .catch(err => setError(err.userMessage || err.message || 'Analysis failed'))
      .finally(() => setLoading(false))
  }, [phaseId])

  useEffect(() => { fetchData() }, [phaseId]) // eslint-disable-line

  // Available properties
  const numericProps = useMemo(() =>
    (sarData?.available_properties || []).filter(p => p.coverage >= 0.2),
  [sarData])

  // Toggle a property on/off
  const handleToggle = useCallback((key) => {
    setSelectedKeys(prev => {
      if (prev.includes(key)) {
        // Don't allow deselecting the last one
        if (prev.length <= 1) return prev
        return prev.filter(k => k !== key)
      }
      return [...prev, key]
    })
  }, [])

  // The selected property objects in order
  const selectedProps = useMemo(() =>
    selectedKeys.map(k => numericProps.find(p => p.key === k)).filter(Boolean),
  [selectedKeys, numericProps])

  // Filter cliffs by SALI threshold
  const cliffPairs = useMemo(() => {
    const pairs = sarData?.cliffs?.pairs || []
    return pairs.filter(p => p.sali >= saliMin)
  }, [sarData, saliMin])

  // Molecule properties map for multi-score deltas
  const moleculeProperties = sarData?.molecule_properties || {}

  const phaseTypeMeta = currentPhase ? (PHASE_TYPES[currentPhase.type] || {}) : {}

  // Loading
  if (loading && !sarData) {
    return (
      <div className="min-h-screen bg-bx-bg flex flex-col items-center justify-center gap-4">
        <BindXLogo variant="loading" size={48} />
        <p className="text-sm text-gray-400">Analyzing activity cliffs...</p>
      </div>
    )
  }

  // Error
  if (error && !sarData) {
    return (
      <div className="min-h-screen bg-bx-bg flex flex-col items-center justify-center gap-4">
        <p className="text-sm text-red-400">{error}</p>
        <button onClick={() => fetchData()} className="px-4 py-2 rounded-lg text-sm font-medium bg-bx-surface text-white hover:bg-bx-elevated transition">Retry</button>
      </div>
    )
  }

  if (!sarData) return null

  return (
    <div className="min-h-screen bg-bx-bg">
      {/* Header */}
      <div className="px-6 py-4 border-b border-gray-100">
        <nav className="flex items-center gap-1.5 text-sm text-gray-400 mb-2">
          <Link to="/" className="hover:text-bx-mint transition-colors">Projects</Link>
          <ChevronRight />
          <Link to={`/project/${projectId}`} className="hover:text-bx-light-text transition-colors truncate max-w-[120px]">{currentProject?.name || 'Project'}</Link>
          <ChevronRight />
          <Link to={`/project/${projectId}/phase/${phaseId}`} className="hover:text-bx-light-text transition-colors">Phase {phaseTypeMeta.short || '?'}</Link>
          <ChevronRight />
          <span className="text-gray-600 font-medium">Activity Explorer</span>
        </nav>

        <div className="flex items-center justify-between">
          <div>
            <h1 className="text-xl font-bold text-bx-light-text">Activity Explorer</h1>
            <p className="text-xs text-gray-400 mt-0.5 max-w-xl">
              Pairs of structurally similar molecules with large score differences.
              A small structural change with a big effect is worth investigating.
            </p>
          </div>
          <Link to={`/project/${projectId}/phase/${phaseId}`} className="flex items-center gap-1.5 px-3 py-2 rounded-lg text-sm font-medium border border-gray-200 text-gray-600 hover:bg-gray-50 transition">
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10 19l-7-7m0 0l7-7m-7 7h18" /></svg>
            Back to Dashboard
          </Link>
        </div>
      </div>

      {/* Controls */}
      <div className="px-6 pt-4">
        <div className="card px-5 py-3 space-y-3">
          <div className="flex flex-wrap items-center gap-6">
            {/* Score chips */}
            <div className="flex items-center gap-2.5">
              <label className="text-xs text-gray-500 font-medium whitespace-nowrap">Scores</label>
              <ScoreSelector properties={numericProps} selected={selectedKeys} onToggle={handleToggle} />
            </div>
          </div>

          {/* SALI slider */}
          <div className="flex items-center gap-3">
            <label className="text-xs text-gray-500 font-medium whitespace-nowrap">Min. SALI</label>
            <input
              type="range"
              min={0.3}
              max={1.5}
              step={0.05}
              value={saliMin}
              onChange={e => setSaliMin(parseFloat(e.target.value))}
              className="flex-1 h-1.5 rounded-full accent-amber-500 max-w-xs"
            />
            <span className="text-xs text-gray-600 font-mono whitespace-nowrap">
              SALI ≥ {saliMin.toFixed(2)} · {cliffPairs.length} pair{cliffPairs.length !== 1 ? 's' : ''}
            </span>
          </div>
        </div>
      </div>

      {/* Table */}
      <div className="px-6 py-4">
        <div className="card overflow-hidden">
          {cliffPairs.length === 0 ? (
            <div className="px-6 py-12 text-center">
              <p className="text-sm text-gray-400">No activity cliffs above SALI ≥ {saliMin.toFixed(2)}.</p>
              <p className="text-xs text-gray-300 mt-1">Try lowering the similarity threshold.</p>
            </div>
          ) : (
            <div className="overflow-x-auto">
              <table className="w-full text-sm border-collapse">
                <thead>
                  <tr className="border-b-2 border-gray-200 bg-gray-50/50">
                    <th className="px-4 py-3 text-left text-xs text-gray-500 font-medium">Molecule A</th>
                    <th className="px-4 py-3 text-left text-xs text-gray-500 font-medium">Molecule B</th>
                    <th className="px-4 py-3 text-right text-xs text-gray-500 font-medium">Similarity</th>
                    {selectedProps.map(sp => (
                      <th key={sp.key} className="px-3 py-3 text-right text-xs text-gray-500 font-medium whitespace-nowrap">
                        Δ {sp.label}
                      </th>
                    ))}
                    <th className="px-4 py-3 text-right text-xs text-gray-500 font-medium">SALI</th>
                  </tr>
                </thead>
                <tbody>
                  {cliffPairs.map((p, i) => (
                    <CliffPairRow
                      key={i}
                      pair={p}
                      selectedProps={selectedProps}
                      moleculeProperties={moleculeProperties}
                    />
                  ))}
                </tbody>
              </table>
            </div>
          )}
        </div>
      </div>

      {/* Footer */}
      <div className="px-6 py-3 border-t border-gray-100">
        <p className="text-[10px] text-gray-300">
          Activity cliffs: SALI (Guha & Van Drie, JCIM 2008) | ECFP4 fingerprints (Morgan radius 2, 2048 bits)
        </p>
      </div>
    </div>
  )
}
