import React, { useEffect, useState, useMemo, useRef, useCallback } from 'react'
import { useParams, Link, useNavigate } from 'react-router-dom'
import { useWorkspace } from '../contexts/WorkspaceContext.jsx'
import { PHASE_TYPES } from '../lib/columns.js'
import { v9GetSARAnalysis, v9ApplyMMPTransform, v9ImportMMPSuggestions } from '../api.js'
import BindXLogo from '../components/BindXLogo.jsx'
import SARMatrixView from '../components/SARMatrixView.jsx'
import MMPSuggestionsModal from '../components/MMPSuggestionsModal.jsx'

// ──────────────────────────────────────────────────────────────────────────
// Helpers
// ──────────────────────────────────────────────────────────────────────────

function fmt(v, d = 2) { return v == null ? '—' : v.toFixed(d) }
function fmtDelta(v) { return v == null ? '—' : `${v >= 0 ? '+' : ''}${v.toFixed(3)}` }

function SmilesCanvas({ smiles, width = 80, height = 60 }) {
  const canvasRef = useRef(null)
  const [failed, setFailed] = useState(false)
  useEffect(() => {
    if (!smiles || !canvasRef.current) return
    if (smiles === 'H' || smiles === '[H]' || smiles.length < 2) { setFailed(true); return }
    let cancelled = false
    setFailed(false)
    import('smiles-drawer').then(mod => {
      if (cancelled) return
      const SD = mod.default || mod
      const d = new SD.Drawer({ width, height, bondThickness: 1 })
      SD.parse(smiles, t => { if (!cancelled) d.draw(t, canvasRef.current, 'light') }, () => { if (!cancelled) setFailed(true) })
    }).catch(() => { if (!cancelled) setFailed(true) })
    return () => { cancelled = true }
  }, [smiles, width, height])
  if (failed || !smiles || smiles === 'H') return <span className="text-[10px] font-mono text-gray-400 px-1">{smiles || 'H'}</span>
  return <canvas ref={canvasRef} width={width} height={height} className="bg-white rounded" />
}

function ChevronRight() {
  return <svg className="w-3 h-3 text-gray-300" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" /></svg>
}

// ──────────────────────────────────────────────────────────────────────────
// Series Selector
// ──────────────────────────────────────────────────────────────────────────

function SeriesSelector({ series, selectedScaffold, onSelect }) {
  if (!series || series.length <= 1) return null
  return (
    <div className="card px-5 py-4">
      <h2 className="text-sm font-semibold text-bx-light-text mb-1">Chemical Series</h2>
      <p className="text-xs text-gray-400 mb-3">
        Your molecules have been grouped by their core structure (scaffold). Select a series to focus the analysis.
        Larger series give more reliable results.
      </p>
      <div className="flex gap-3 overflow-x-auto pb-1">
        {series.map(s => {
          const sel = s.scaffold_smiles === selectedScaffold
          return (
            <button
              key={s.scaffold_smiles ?? '__other__'}
              onClick={() => onSelect(s.scaffold_smiles)}
              className={`relative flex-shrink-0 flex flex-col items-center gap-1 p-2.5 rounded-lg border-2 transition-all cursor-pointer hover:shadow-md ${
                sel ? 'border-amber-400 ring-2 ring-amber-200 bg-amber-50/50' : 'border-gray-200 bg-white hover:border-gray-300'
              }`}
              style={{ minWidth: 120 }}
            >
              {s.scaffold_svg ? (
                <div className="rounded bg-white" dangerouslySetInnerHTML={{ __html: s.scaffold_svg }} style={{ width: 100, height: 70, overflow: 'hidden' }} />
              ) : (
                <div className="flex items-center justify-center rounded bg-gray-50 text-gray-400 text-xs" style={{ width: 100, height: 70 }}>No scaffold</div>
              )}
              <div className={`text-sm font-bold ${sel ? 'text-amber-600' : 'text-bx-light-text'}`}>{s.n_molecules} mol</div>
              {sel && (
                <div className="absolute -top-1 -right-1 w-4 h-4 bg-amber-400 rounded-full flex items-center justify-center">
                  <svg className="w-2.5 h-2.5 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={3} d="M5 13l4 4L19 7" /></svg>
                </div>
              )}
            </button>
          )
        })}
      </div>
    </div>
  )
}

// ──────────────────────────────────────────────────────────────────────────
// Overview (stats only, no graph)
// ──────────────────────────────────────────────────────────────────────────

function OverviewSection({ data }) {
  const rg = data.rgroup
  return (
    <div className="card px-5 py-4">
      <h2 className="text-base font-semibold text-bx-light-text mb-3">Analysis Summary</h2>

      <div className="grid grid-cols-2 sm:grid-cols-4 gap-3 mb-4">
        <Stat label="Molecules analyzed" value={data.n_molecules} sub={data.n_molecules_total !== data.n_molecules ? `of ${data.n_molecules_total} total` : null} />
        <Stat label="With activity data" value={data.n_with_activity} />
        <Stat label="R-group match" value={rg ? `${Math.round(rg.match_rate * 100)}%` : '—'} sub={rg ? `${rg.n_matched} of ${rg.n_total}` : null} />
        <Stat label="Variable positions" value={rg?.variable_r_groups?.length ?? 0} sub={rg?.variable_r_groups?.join(', ')} />
      </div>

      {data.scaffold_svg && (
        <div className="flex items-start gap-4">
          <div className="flex-shrink-0 border border-gray-100 rounded-lg p-2 bg-white" dangerouslySetInnerHTML={{ __html: data.scaffold_svg }} style={{ maxWidth: 250 }} />
          <div className="text-xs text-gray-500 space-y-1.5 pt-1">
            <p><strong>Core scaffold</strong> = the common structural backbone shared by all molecules in this series.</p>
            <p>The analysis identifies which <em>substituent positions</em> (R-groups) vary between molecules, and how those variations affect the measured properties.</p>
          </div>
        </div>
      )}

      {data.error && (
        <div className="mt-3 flex items-center gap-2 px-3 py-2 rounded-lg bg-amber-50 border border-amber-200 text-amber-700 text-xs">
          Partial analysis — {data.error}
        </div>
      )}
    </div>
  )
}

function Stat({ label, value, sub }) {
  return (
    <div className="dark-inset text-center py-3 px-2 rounded-lg">
      <p className="text-lg font-bold text-bx-light-text">{value ?? '—'}</p>
      <p className="text-[10px] text-gray-400 uppercase tracking-wide">{label}</p>
      {sub && <p className="text-[10px] text-gray-400 mt-0.5">{sub}</p>}
    </div>
  )
}

// ──────────────────────────────────────────────────────────────────────────
// Activity Cliffs — TABLE ONLY, no scatter chart
// ──────────────────────────────────────────────────────────────────────────

function ActivityCliffsSection({ data, propLabel }) {
  const cliffs = data?.cliffs
  if (!cliffs?.pairs?.length) return null

  return (
    <div className="card px-5 py-4">
      <h2 className="text-base font-semibold text-bx-light-text mb-1">Activity Cliffs</h2>

      <div className="text-xs text-gray-500 bg-blue-50 border border-blue-100 rounded-lg px-4 py-3 mb-4 space-y-1.5">
        <p className="font-semibold text-blue-700">What are activity cliffs?</p>
        <p>Two molecules that look <strong>almost identical</strong> structurally but have <strong>very different activity</strong>. This is important because it means a tiny structural change has a big effect — useful for optimization, but also a risk.</p>
        <p><strong>Similarity</strong> = how structurally close the two molecules are (0 = different, 1 = identical). <strong>|Δ {propLabel}|</strong> = how much their activity differs. <strong>SALI</strong> = cliff severity score (higher = more surprising).</p>
        <p>{cliffs.n_cliffs} cliff pair{cliffs.n_cliffs > 1 ? 's' : ''} detected. Threshold: SALI ≥ {cliffs.sali_threshold} (top 5% most extreme).</p>
      </div>

      <div className="overflow-x-auto max-h-[400px] overflow-y-auto">
        <table className="w-full text-xs border-collapse">
          <thead className="sticky top-0 bg-white z-10">
            <tr className="border-b-2 border-gray-200">
              <th className="p-2 text-left text-gray-500 font-medium">Molecule A</th>
              <th className="p-2 text-left text-gray-500 font-medium">Molecule B</th>
              <th className="p-2 text-right text-gray-500 font-medium" title="Structural similarity (0–1)">Similarity</th>
              <th className="p-2 text-right text-gray-500 font-medium" title="Absolute difference in activity">|Δ {propLabel}|</th>
              <th className="p-2 text-right text-gray-500 font-medium" title="Structure-Activity Landscape Index — higher = bigger cliff">SALI</th>
            </tr>
          </thead>
          <tbody>
            {cliffs.pairs.slice(0, 20).map((p, i) => (
              <tr key={i} className="border-b border-gray-100 hover:bg-red-50/20">
                <td className="p-2">
                  <div className="font-medium truncate max-w-[140px]" title={p.mol_a_name}>{p.mol_a_name || '—'}</div>
                </td>
                <td className="p-2">
                  <div className="font-medium truncate max-w-[140px]" title={p.mol_b_name}>{p.mol_b_name || '—'}</div>
                </td>
                <td className="p-2 text-right font-mono">{fmt(p.similarity, 3)}</td>
                <td className="p-2 text-right font-semibold text-amber-600">{fmt(Math.abs(p.delta_activity), 3)}</td>
                <td className="p-2 text-right font-bold text-red-600">{fmt(p.sali, 1)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      {cliffs.pairs.length > 20 && (
        <p className="text-[10px] text-gray-400 mt-2">Showing top 20 of {cliffs.pairs.length} cliff pairs (sorted by severity).</p>
      )}
    </div>
  )
}

// ──────────────────────────────────────────────────────────────────────────
// MMP Rules — TABLE ONLY, no bar chart, with Apply button
// ──────────────────────────────────────────────────────────────────────────

function MMPSection({ data, propLabel, onApply }) {
  const mmps = data?.mmps
  if (!mmps?.rules?.length) return null

  const CONF = { high: { bg: '#dcfce7', text: '#15803d', label: 'High' }, medium: { bg: '#fef3c7', text: '#b45309', label: 'Medium' }, low: { bg: '#f3f4f6', text: '#6b7280', label: 'Low' } }

  return (
    <div className="card px-5 py-4">
      <h2 className="text-base font-semibold text-bx-light-text mb-1">Matched Molecular Pair Rules</h2>

      <div className="text-xs text-gray-500 bg-blue-50 border border-blue-100 rounded-lg px-4 py-3 mb-4 space-y-1.5">
        <p className="font-semibold text-blue-700">What are MMP rules?</p>
        <p>A Matched Molecular Pair = two molecules that differ by exactly one structural change (e.g. replacing -H by -F). When we observe this same change across multiple molecule pairs, we get a <strong>rule</strong>.</p>
        <p><strong>Transformation</strong> = the structural change (A → B). <strong>Mean Δ</strong> = average effect on {propLabel}. <strong>Pairs</strong> = how many times this change was observed (more = more reliable).</p>
        <p><strong>Confidence</strong>: <span className="text-green-700">High</span> = ≥5 pairs with consistent direction. <span className="text-amber-700">Medium</span> = ≥3 pairs. <span className="text-gray-500">Low</span> = fewer pairs.</p>
        <p>Click <strong>Apply</strong> on a rule to generate new molecule suggestions using that transformation.</p>
      </div>

      <div className="overflow-x-auto max-h-[400px] overflow-y-auto">
        <table className="w-full text-xs border-collapse">
          <thead className="sticky top-0 bg-white z-10">
            <tr className="border-b-2 border-gray-200">
              <th className="p-2 text-left text-gray-500 font-medium">Transformation</th>
              <th className="p-2 text-center text-gray-500 font-medium">Confidence</th>
              <th className="p-2 text-right text-gray-500 font-medium">Pairs</th>
              <th className="p-2 text-right text-gray-500 font-medium" title="Average change in property value">Mean Δ {propLabel}</th>
              <th className="p-2 text-center text-gray-500 font-medium">Action</th>
            </tr>
          </thead>
          <tbody>
            {mmps.rules.map((r, i) => {
              const c = CONF[r.confidence] || CONF.low
              return (
                <tr key={i} className="border-b border-gray-100 hover:bg-gray-50">
                  <td className="p-2 font-mono text-[11px] max-w-[220px] truncate" title={r.transform}>{r.transform}</td>
                  <td className="p-2 text-center">
                    <span className="inline-block px-1.5 py-0.5 rounded text-[10px] font-medium" style={{ backgroundColor: c.bg, color: c.text }}>
                      {c.label}
                    </span>
                  </td>
                  <td className="p-2 text-right">{r.n_pairs}</td>
                  <td className={`p-2 text-right font-semibold ${r.mean_delta != null && r.mean_delta >= 0 ? 'text-green-600' : 'text-red-600'}`}>
                    {fmtDelta(r.mean_delta)}
                  </td>
                  <td className="p-2 text-center">
                    {r.from_smiles && r.to_smiles && (
                      <button
                        onClick={() => onApply(r)}
                        className="inline-flex items-center gap-1 px-2.5 py-1 rounded text-[10px] font-medium bg-violet-50 text-violet-600 hover:bg-violet-100 transition"
                      >
                        Apply
                      </button>
                    )}
                  </td>
                </tr>
              )
            })}
          </tbody>
        </table>
      </div>
      <p className="text-[10px] text-gray-400 mt-2">{mmps.n_rules} rules from {mmps.n_pairs_total} molecule pairs.</p>
    </div>
  )
}

// ──────────────────────────────────────────────────────────────────────────
// Main page
// ──────────────────────────────────────────────────────────────────────────

export default function SARExplorerPage() {
  const { projectId, phaseId } = useParams()
  const navigate = useNavigate()
  const { currentProject, currentPhase, selectProject, selectPhase } = useWorkspace()

  useEffect(() => { if (projectId) selectProject(projectId) }, [projectId]) // eslint-disable-line react-hooks/exhaustive-deps
  useEffect(() => { if (phaseId) selectPhase(phaseId) }, [phaseId]) // eslint-disable-line react-hooks/exhaustive-deps

  const [sarData, setSarData] = useState(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState(null)
  const [selectedScaffold, setSelectedScaffold] = useState(null)

  const [mmpModal, setMmpModal] = useState({ open: false, loading: false, suggestions: [], transform: '', nInput: 0 })
  const [importing, setImporting] = useState(false)

  const propLabel = useMemo(() => {
    if (!sarData) return 'Score'
    const p = sarData.available_properties?.find(p => p.key === sarData.property_key)
    return p?.label || sarData.property_key || 'Score'
  }, [sarData])

  const fetchSAR = useCallback((scaffold = null) => {
    if (!phaseId) return
    setLoading(true)
    setError(null)
    v9GetSARAnalysis(phaseId, 'auto', scaffold)
      .then(result => {
        setSarData(result)
        if (result.selected_scaffold !== undefined) setSelectedScaffold(result.selected_scaffold)
      })
      .catch(err => setError(err.userMessage || err.message || 'SAR analysis failed'))
      .finally(() => setLoading(false))
  }, [phaseId])

  useEffect(() => { fetchSAR() }, [phaseId]) // eslint-disable-line react-hooks/exhaustive-deps

  const handleSeriesChange = scaffoldSmiles => {
    setSelectedScaffold(scaffoldSmiles)
    setSarData(null)
    fetchSAR(scaffoldSmiles)
  }

  const handleApplyMMP = useCallback(rule => {
    if (!phaseId) return
    setMmpModal({ open: true, loading: true, suggestions: [], transform: rule.transform, nInput: 0 })
    v9ApplyMMPTransform(phaseId, rule.from_smiles, rule.to_smiles)
      .then(r => setMmpModal(prev => ({ ...prev, loading: false, suggestions: r.suggestions || [], nInput: r.n_input_molecules || 0 })))
      .catch(() => setMmpModal(prev => ({ ...prev, loading: false, suggestions: [] })))
  }, [phaseId])

  const handleImport = useCallback(selected => {
    if (!phaseId || !selected.length) return
    setImporting(true)
    v9ImportMMPSuggestions(phaseId, selected.map(s => s.smiles), selected.map((_, i) => `MMP-${i + 1}`), mmpModal.transform)
      .then(r => { setMmpModal(prev => ({ ...prev, open: false })); if (r.imported > 0) navigate(`/project/${projectId}/phase/${phaseId}`) })
      .catch(err => console.error('MMP import failed:', err))
      .finally(() => setImporting(false))
  }, [phaseId, projectId, mmpModal.transform, navigate])

  const phaseTypeMeta = currentPhase ? (PHASE_TYPES[currentPhase.type] || {}) : {}

  // Loading
  if (loading && !sarData) {
    return (
      <div className="min-h-screen bg-bx-bg flex flex-col items-center justify-center gap-4">
        <BindXLogo variant="loading" size={48} />
        <p className="text-sm text-gray-400">Analyzing structure-activity relationships...</p>
      </div>
    )
  }

  // Error
  if (error && !sarData) {
    return (
      <div className="min-h-screen bg-bx-bg flex flex-col items-center justify-center gap-4">
        <p className="text-sm text-red-400">{error}</p>
        <button onClick={() => fetchSAR()} className="px-4 py-2 rounded-lg text-sm font-medium bg-bx-surface text-white hover:bg-bx-elevated transition">Retry</button>
      </div>
    )
  }

  if (!sarData) return null

  const rgroup = sarData.rgroup
  const variableRGroups = rgroup?.variable_r_groups || []

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
          <span className="text-gray-600 font-medium">SAR Explorer</span>
        </nav>
        <div className="flex items-center justify-between">
          <div>
            <h1 className="text-xl font-bold text-bx-light-text">SAR Explorer</h1>
            <p className="text-xs text-gray-400 mt-0.5">
              Understand how structural changes affect molecule properties. Identify optimization opportunities.
            </p>
          </div>
          <Link to={`/project/${projectId}/phase/${phaseId}`} className="flex items-center gap-1.5 px-3 py-2 rounded-lg text-sm font-medium border border-gray-200 text-gray-600 hover:bg-gray-50 transition">
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10 19l-7-7m0 0l7-7m-7 7h18" /></svg>
            Back to Dashboard
          </Link>
        </div>
      </div>

      {/* What is this page — top-level explanation */}
      <div className="px-6 pt-4">
        <div className="text-xs text-gray-500 bg-gray-50 border border-gray-200 rounded-lg px-4 py-3 space-y-1">
          <p className="font-semibold text-gray-600">What does this page do?</p>
          <p>This page analyzes your molecules to find patterns between their <strong>structure</strong> (shape) and their <strong>activity</strong> (how well they work). It answers three key questions:</p>
          <ol className="list-decimal list-inside space-y-0.5 ml-2">
            <li><strong>SAR Matrix</strong> — What happens when I change one substituent while keeping the rest fixed?</li>
            <li><strong>Activity Cliffs</strong> — Are there cases where a tiny structural change causes a huge activity jump?</li>
            <li><strong>MMP Rules</strong> — What structural changes consistently improve (or worsen) the activity across many molecules?</li>
          </ol>
        </div>
      </div>

      {/* Series selector */}
      {sarData.series?.length > 1 && (
        <div className="px-6 pt-4">
          <SeriesSelector series={sarData.series} selectedScaffold={selectedScaffold} onSelect={handleSeriesChange} />
        </div>
      )}

      {/* Sections */}
      <div className="px-6 py-4 space-y-4">
        <OverviewSection data={sarData} />

        {rgroup && variableRGroups.length >= 2 && (
          <SARMatrixView
            rgroup={rgroup}
            moleculeProperties={sarData.molecule_properties || {}}
            availableProperties={sarData.available_properties || []}
            variableRGroups={variableRGroups}
          />
        )}

        <ActivityCliffsSection data={sarData} propLabel={propLabel} />
        <MMPSection data={sarData} propLabel={propLabel} onApply={handleApplyMMP} />
      </div>

      {/* Footer */}
      <div className="px-6 py-4 border-t border-gray-100">
        <div className="text-[10px] text-gray-300 space-y-0.5">
          <p>Scaffolds: Bemis-Murcko generic (J. Med. Chem. 1996) | R-groups: RDKit RGroupDecompose | Activity cliffs: SALI (Guha & Van Drie, JCIM 2008) | ECFP4 fingerprints</p>
          <p>MMP: Hussain-Rea algorithm (rdMMPA) | SAR Matrix: Wawer & Bajorath 2007</p>
        </div>
      </div>

      <MMPSuggestionsModal
        open={mmpModal.open}
        onClose={() => setMmpModal(prev => ({ ...prev, open: false }))}
        suggestions={mmpModal.suggestions}
        transform={mmpModal.transform}
        nInputMolecules={mmpModal.nInput}
        loading={mmpModal.loading}
        onImport={handleImport}
        importing={importing}
      />
    </div>
  )
}
