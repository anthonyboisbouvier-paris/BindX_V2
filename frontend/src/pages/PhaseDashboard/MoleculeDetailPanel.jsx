import React, { useState, useRef, useEffect, useMemo, useCallback } from 'react'
import Badge from '../../components/Badge.jsx'
import ProteinViewer from '../../components/ProteinViewer.jsx'
import LigandViewer3D from '../../components/LigandViewer3D.jsx'
import SafetyReport from '../../components/SafetyReport.jsx'
import SynthesisTree from '../../components/SynthesisTree.jsx'
import ConfidenceBreakdown from '../../components/ConfidenceBreakdown.jsx'
import InfoTip, { TIPS } from '../../components/InfoTip.jsx'
import InteractionDiagram from '../../components/InteractionDiagram.jsx'
import InteractionView from '../../components/InteractionView.jsx'

// ---------------------------------------------------------------------------
// Error boundary — catches render errors in detail panel content
// ---------------------------------------------------------------------------
class PanelErrorBoundary extends React.Component {
  constructor(props) {
    super(props)
    this.state = { error: null }
  }
  static getDerivedStateFromError(error) {
    return { error }
  }
  componentDidCatch(error, info) {
    console.error('[PanelErrorBoundary] Render error:', error?.message, '\nComponent stack:', info?.componentStack)
  }
  componentDidUpdate(prevProps) {
    // Auto-reset when molecule changes
    if (prevProps.resetKey !== this.props.resetKey && this.state.error) {
      this.setState({ error: null })
    }
  }
  render() {
    if (this.state.error) {
      return (
        <div className="p-4 bg-red-50 rounded-lg border border-red-200 text-center space-y-2">
          <p className="text-sm font-semibold text-red-700">Render error in detail panel</p>
          <p className="text-xs text-red-500 font-mono break-all">{String(this.state.error?.message || this.state.error)}</p>
          <button
            onClick={() => this.setState({ error: null })}
            className="px-3 py-1 text-xs bg-red-100 hover:bg-red-200 text-red-700 rounded-lg transition-colors"
          >
            Retry
          </button>
        </div>
      )
    }
    return this.props.children
  }
}

// ---------------------------------------------------------------------------
// Detail panel tab content helpers
// ---------------------------------------------------------------------------
function ScoresTab({ mol }) {
  const scores = [
    { label: 'Docking', value: mol.docking_score, unit: 'kcal/mol', color: 'bg-bx-surface', textColor: 'text-bx-light-text', format: v => Number(v).toFixed(1) },
    { label: 'CNN Score', value: mol.cnn_score, unit: '', color: 'bg-green-500', textColor: 'text-green-600', format: v => Number(v).toFixed(2) },
    { label: 'CNN Aff.', value: mol.cnn_affinity, unit: 'pKi', color: 'bg-teal-500', textColor: 'text-teal-600', format: v => Number(v).toFixed(1) },
  ]
  return (
    <div className="space-y-3">
      <div className="grid grid-cols-2 gap-2">
        {scores.filter(s => s.value != null).map(s => (
          <div key={s.label} className="bg-gray-50 rounded-xl p-3 border border-gray-100">
            <p className="text-[9px] text-gray-400 uppercase tracking-wider mb-1">{s.label}</p>
            <p className={`text-xl font-bold tabular-nums ${s.textColor}`}>
              {s.format(s.value)}
            </p>
            {s.unit && <p className="text-[9px] text-gray-400">{s.unit}</p>}
            <div className="mt-2 h-1 bg-gray-200 rounded-full overflow-hidden">
              <div
                className={`h-1 ${s.color} rounded-full`}
                style={{ width: (() => {
                  if (s.label === 'CNN Score') return `${Math.min(100, Number(s.value) * 100)}%`
                  if (s.label === 'CNN Affinity') return `${Math.min(100, Number(s.value) / 10 * 100)}%`
                  // Docking score: range -12 to 0, more negative = better
                  return `${Math.min(100, Math.abs(Number(s.value)) / 12 * 100)}%`
                })() }}
              />
            </div>
          </div>
        ))}
      </div>
      {scores.every(s => s.value == null) && (
        <p className="text-sm text-gray-400 text-center py-4">No docking scores yet — run a Docking analysis.</p>
      )}

      {/* Ligand Preparation metadata */}
      {mol.preparation && typeof mol.preparation === 'object' && (
        <div className="mt-4 pt-3 border-t border-gray-100">
          <p className="text-[9px] text-gray-400 uppercase tracking-wider mb-2">Ligand Preparation</p>
          <div className="flex flex-wrap gap-1.5 mb-2">
            {mol.preparation.tautomer_changed && (
              <span className="inline-flex items-center px-2 py-0.5 rounded-full text-[10px] font-medium bg-amber-100 text-amber-700">
                Tautomer changed
              </span>
            )}
            {mol.preparation.protonated && (
              <span className="inline-flex items-center px-2 py-0.5 rounded-full text-[10px] font-medium bg-green-100 text-green-700">
                Protonated pH 7.4
              </span>
            )}
            {mol.preparation.standardized && (
              <span className="inline-flex items-center px-2 py-0.5 rounded-full text-[10px] font-medium bg-blue-100 text-blue-700">
                Standardized
              </span>
            )}
            {mol.preparation.minimization && mol.preparation.minimization !== 'none' && (
              <span className="inline-flex items-center px-2 py-0.5 rounded-full text-[10px] font-medium bg-gray-100 text-gray-600">
                {mol.preparation.minimization}
              </span>
            )}
            {mol.preparation.rotatable_bonds != null && (
              <span className="inline-flex items-center px-2 py-0.5 rounded-full text-[10px] font-medium bg-gray-100 text-gray-600">
                {mol.preparation.rotatable_bonds} rot. bonds
              </span>
            )}
          </div>
          {mol.preparation.prepared_smiles && mol.preparation.input_smiles &&
           mol.preparation.prepared_smiles !== mol.preparation.input_smiles && (
            <div className="text-[10px] text-gray-500 space-y-0.5">
              <p className="truncate"><span className="text-gray-400">Prepared:</span>{' '}
                <code className="font-mono text-[9px] bg-gray-50 px-1 rounded">{mol.preparation.prepared_smiles}</code>
              </p>
            </div>
          )}
        </div>
      )}
    </div>
  )
}

function CompositeScoreTab({ mol }) {
  const score = mol.weighted_score
  const breakdown = mol.breakdown || {}
  const weightsUsed = mol.weights_used || {}

  const METRIC_LABELS = {
    docking_score: 'Docking',
    cnn_score: 'CNN Score',
    logP: 'LogP',
    solubility: 'Solubility',
    selectivity: 'Selectivity',
    qed: 'QED',
    safety: 'Safety',
    novelty: 'Novelty',
  }

  const METRIC_COLORS = {
    docking_score: 'bg-indigo-500',
    cnn_score: 'bg-green-500',
    logP: 'bg-blue-500',
    solubility: 'bg-cyan-500',
    selectivity: 'bg-purple-500',
    qed: 'bg-amber-500',
    safety: 'bg-red-400',
    novelty: 'bg-teal-500',
  }

  if (score == null) {
    return <p className="text-sm text-gray-400 text-center py-4">No composite score yet — run a Composite Score analysis.</p>
  }

  const breakdownEntries = Object.entries(breakdown).filter(([, v]) => v != null)

  return (
    <div className="space-y-3">
      {/* Main score */}
      <div className="bg-amber-50 rounded-xl p-4 border border-amber-200 text-center">
        <p className="text-[9px] text-amber-500 uppercase tracking-wider mb-1">Weighted Score</p>
        <p className="text-3xl font-bold tabular-nums text-amber-600">
          {Math.round(score * 100)}
        </p>
        <p className="text-[10px] text-amber-400">/100</p>
        <div className="mt-2 h-1.5 bg-amber-100 rounded-full overflow-hidden">
          <div className="h-1.5 bg-amber-500 rounded-full transition-all" style={{ width: `${score * 100}%` }} />
        </div>
      </div>

      {/* Breakdown bars */}
      {breakdownEntries.length > 0 && (
        <div className="space-y-2">
          <p className="text-[9px] text-gray-400 uppercase tracking-wider">Breakdown</p>
          {breakdownEntries.map(([key, val]) => {
            const weight = weightsUsed[key]
            return (
              <div key={key} className="flex items-center gap-2">
                <span className="text-[10px] text-gray-500 w-20 truncate">{METRIC_LABELS[key] || key}</span>
                <div className="flex-1 h-2 bg-gray-100 rounded-full overflow-hidden">
                  <div
                    className={`h-2 rounded-full ${METRIC_COLORS[key] || 'bg-gray-400'}`}
                    style={{ width: `${val * 100}%` }}
                  />
                </div>
                <span className="text-[10px] text-gray-500 tabular-nums w-10 text-right">{Math.round(val * 100)}%</span>
                {weight != null && (
                  <span className="text-[9px] text-gray-300 w-8 text-right">×{weight}</span>
                )}
              </div>
            )
          })}
        </div>
      )}
    </div>
  )
}

function PropertiesTab({ mol }) {
  const props = [
    { key: 'MW',   label: 'MW',    unit: 'Da',    threshold: { max: 500 }, decimals: 0 },
    { key: 'logP', label: 'LogP',  unit: '',      threshold: { min: -2, max: 5 }, decimals: 1 },
    { key: 'HBD',  label: 'HBD',   unit: '',      threshold: { max: 5 }, decimals: 0 },
    { key: 'HBA',  label: 'HBA',   unit: '',      threshold: { max: 10 }, decimals: 0 },
    { key: 'TPSA', label: 'TPSA',  unit: 'A²',    threshold: { max: 140 }, decimals: 1 },
    { key: 'QED',  label: 'QED',   unit: '',      threshold: { min: 0.5 }, decimals: 2 },
  ]

  const available = props.filter(p => mol[p.key] != null)
  if (!available.length) {
    return <p className="text-sm text-gray-400 text-center py-4">No physicochemical properties available yet.</p>
  }

  return (
    <div className="space-y-3">
      {mol.lipinski_pass != null && (
        <div className={`flex items-center gap-2 px-3 py-2 rounded-lg text-sm font-medium ${
          mol.lipinski_pass ? 'bg-green-50 text-green-700 border border-green-200' : 'bg-red-50 text-red-600 border border-red-200'
        }`}>
          <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            {mol.lipinski_pass
              ? <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M5 13l4 4L19 7" />
              : <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M6 18L18 6M6 6l12 12" />}
          </svg>
          Lipinski Rule of 5 — {mol.lipinski_pass ? 'Pass' : 'Fail'}
        </div>
      )}
      {available.length > 0 && (
        <div className="grid grid-cols-2 gap-2">
          {available.map(p => {
            const val = mol[p.key]
            const formatted = p.decimals === 0 ? Math.round(val) : val.toFixed(p.decimals)
            const inRange = (() => {
              if (!p.threshold) return null
              if (p.threshold.min != null && val < p.threshold.min) return false
              if (p.threshold.max != null && val > p.threshold.max) return false
              return true
            })()
            return (
              <div key={p.key} className={`bg-gray-50 rounded-lg px-3 py-2 border ${
                inRange === false ? 'border-red-200 bg-red-50' : 'border-gray-100'
              }`}>
                <p className="text-[9px] text-gray-400 uppercase tracking-wider">{p.label}</p>
                <p className={`text-base font-bold tabular-nums ${inRange === false ? 'text-red-600' : 'text-gray-800'}`}>
                  {formatted}
                  {p.unit && <span className="text-[10px] font-normal text-gray-400 ml-0.5">{p.unit}</span>}
                </p>
                {p.threshold && (
                  <p className="text-[9px] text-gray-400">
                    {p.threshold.min != null && p.threshold.max != null ? `${p.threshold.min}–${p.threshold.max}` :
                     p.threshold.max != null ? `≤ ${p.threshold.max}` : `≥ ${p.threshold.min}`}
                  </p>
                )}
              </div>
            )
          })}
        </div>
      )}
    </div>
  )
}

// Helper: format ADMET value for display (raw, no unnecessary transformation)
function formatAdmet(key, val) {
  if (val == null) return null
  // Solubility: logS (mol/L) — show as-is with unit
  if (key === 'solubility') {
    const label = val > -2 ? 'High' : val > -4 ? 'Moderate' : val > -6 ? 'Low' : 'Insoluble'
    return { display: val.toFixed(1), unit: 'logS', label, color: val > -2 ? 'text-green-600' : val > -4 ? 'text-amber-600' : 'text-red-600' }
  }
  // Probabilities (0-1): show as percentage
  if (typeof val === 'number' && val >= 0 && val <= 1) {
    const pct = Math.round(val * 100)
    return { display: `${pct}%`, unit: '', label: null, color: null }
  }
  // Fallback: raw number
  return { display: typeof val === 'number' ? val.toFixed(2) : String(val), unit: '', label: null, color: null }
}

// Helper: risk bar for a probability value (0-1)
function RiskBar({ label, value, invert }) {
  if (value == null || typeof value !== 'number') return null
  const pct = Math.round(value * 100)
  const isRisk = invert ? value > 0.3 : value < 0.3
  const isMid = invert ? value > 0.15 : value < 0.5
  const barColor = isRisk ? 'bg-red-400' : isMid ? 'bg-amber-400' : 'bg-green-400'
  const textColor = isRisk ? 'text-red-600' : isMid ? 'text-amber-600' : 'text-green-600'
  return (
    <div>
      <div className="flex justify-between text-[10px] mb-0.5">
        <span className="text-gray-500">{label}</span>
        <span className={`font-medium tabular-nums ${textColor}`}>{pct}%</span>
      </div>
      <div className="w-full bg-gray-100 rounded-full h-1.5">
        <div className={`h-1.5 rounded-full transition-all ${barColor}`}
          style={{ width: `${Math.min(100, pct)}%` }} />
      </div>
    </div>
  )
}

// Helper: pass/fail pill
function PassFailPill({ label, alert }) {
  if (alert == null) return null
  return (
    <div className={`rounded-lg px-2.5 py-2 border text-center ${alert ? 'bg-red-50 border-red-200' : 'bg-green-50 border-green-200'}`}>
      <p className="text-[9px] text-gray-400 uppercase">{label}</p>
      <p className={`font-bold ${alert ? 'text-red-600' : 'text-green-600'}`}>
        {alert ? 'Alert' : 'Pass'}
      </p>
    </div>
  )
}

function AdmeTab({ mol, details }) {
  const safety = details?.safety

  // Extract CYP data from multiple sources (flat keys OR nested properties)
  const metabSource = mol.properties?.admet?.metabolism || mol.properties?.adme?.metabolism || {}
  const getCyp = (key) => mol[key] ?? metabSource[key] ?? null

  const ADME_KEYS = ['oral_bioavailability', 'intestinal_permeability', 'solubility', 'BBB',
    'plasma_protein_binding', 'half_life', 'pgp_substrate', 'vd']
  const CYP_KEYS = ['cyp1a2_inhibitor', 'cyp2c9_inhibitor', 'cyp2c19_inhibitor', 'cyp2d6_inhibitor', 'cyp3a4_inhibitor']
  const hasAdme = ADME_KEYS.some(k => mol[k] != null) || CYP_KEYS.some(k => getCyp(k) != null)
  const hasDruglikeness = safety?.pains_pass != null || safety?.brenk_alert != null ||
    safety?.pfizer_alert != null || safety?.gsk_alert != null || safety?.cns_mpo != null

  if (!hasAdme && !hasDruglikeness) {
    return <p className="text-sm text-gray-400 text-center py-4">No ADME data — run an ADME analysis first.</p>
  }

  return (
    <div className="space-y-3">
      {/* ── DRUGLIKENESS ALERTS ── */}
      {hasDruglikeness && (
        <div>
          <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wider mb-1.5">Druglikeness Rules</p>
          {safety.pains_pass != null && (
            <div className={`flex items-center justify-between px-3 py-1.5 rounded-lg text-sm border mb-1.5 ${
              safety.pains_pass ? 'bg-green-50 border-green-200 text-green-700' : 'bg-red-50 border-red-200 text-red-600'
            }`}>
              <span className="font-semibold">PAINS</span>
              <span>{safety.pains_pass ? 'Pass' : `Fail (${safety.pains_alerts?.length || 0})`}</span>
            </div>
          )}
          <div className="grid grid-cols-3 gap-1.5">
            <PassFailPill label="Brenk" alert={safety.brenk_alert} />
            <PassFailPill label="Pfizer 3/75" alert={safety.pfizer_alert} />
            <PassFailPill label="GSK 4/400" alert={safety.gsk_alert} />
          </div>
          {safety.cns_mpo != null && (
            <div className="mt-1.5 bg-gray-50 rounded-lg px-2.5 py-2 border border-gray-100">
              <div className="flex justify-between text-[10px] mb-0.5">
                <span className="text-gray-500 font-semibold">CNS MPO</span>
                <span className={`font-bold ${safety.cns_mpo >= 4 ? 'text-green-600' : safety.cns_mpo >= 3 ? 'text-amber-600' : 'text-red-600'}`}>
                  {safety.cns_mpo.toFixed(1)} / 6
                </span>
              </div>
              <div className="w-full bg-gray-100 rounded-full h-1.5">
                <div className={`h-1.5 rounded-full ${safety.cns_mpo >= 4 ? 'bg-green-400' : safety.cns_mpo >= 3 ? 'bg-amber-400' : 'bg-red-400'}`}
                  style={{ width: `${Math.min(100, (safety.cns_mpo / 6) * 100)}%` }} />
              </div>
            </div>
          )}
        </div>
      )}

      {/* ── ABSORPTION ── */}
      {(mol.oral_bioavailability != null || mol.intestinal_permeability != null ||
        mol.solubility != null || mol.pgp_substrate != null) && (
        <div className="space-y-1.5">
          <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wider">Absorption</p>
          <RiskBar label="Oral bioavailability" value={mol.oral_bioavailability} />
          <RiskBar label="Intestinal permeability" value={mol.intestinal_permeability} />
          {mol.solubility != null && (() => {
            const f = formatAdmet('solubility', mol.solubility)
            return (
              <div className="flex items-center justify-between bg-gray-50 rounded-lg px-2.5 py-1.5 border border-gray-100">
                <span className="text-[10px] text-gray-500">Solubility</span>
                <div className="flex items-center gap-1.5">
                  <span className={`text-[10px] font-bold tabular-nums ${f.color}`}>{f.display}</span>
                  <span className="text-[9px] text-gray-400">{f.unit}</span>
                  <span className={`text-[9px] font-semibold px-1.5 py-0.5 rounded-full border ${
                    f.label === 'High' ? 'bg-green-50 border-green-200 text-green-600' :
                    f.label === 'Moderate' ? 'bg-amber-50 border-amber-200 text-amber-600' :
                    'bg-red-50 border-red-200 text-red-600'
                  }`}>{f.label}</span>
                </div>
              </div>
            )
          })()}
          <RiskBar label="P-gp substrate" value={mol.pgp_substrate} invert />
        </div>
      )}

      {/* ── DISTRIBUTION ── */}
      {(mol.plasma_protein_binding != null || mol.BBB != null || mol.vd != null) && (
        <div className="space-y-1.5">
          <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wider">Distribution</p>
          {mol.plasma_protein_binding != null && (
            <div className="flex items-center justify-between bg-gray-50 rounded-lg px-2.5 py-1.5 border border-gray-100">
              <span className="text-[10px] text-gray-500">Plasma protein binding</span>
              <span className={`text-[10px] font-bold tabular-nums ${
                mol.plasma_protein_binding > 0.9 ? 'text-red-600' : mol.plasma_protein_binding > 0.7 ? 'text-amber-600' : 'text-green-600'
              }`}>{Math.round(mol.plasma_protein_binding * 100)}%</span>
            </div>
          )}
          <RiskBar label="BBB permeability" value={mol.BBB} />
          {mol.vd != null && (
            <div className="flex items-center justify-between bg-gray-50 rounded-lg px-2.5 py-1.5 border border-gray-100">
              <span className="text-[10px] text-gray-500">Volume of distribution</span>
              <span className="text-[10px] font-bold tabular-nums text-gray-700">{typeof mol.vd === 'number' ? mol.vd.toFixed(2) : mol.vd}</span>
            </div>
          )}
        </div>
      )}

      {/* ── METABOLISM (CYP) ── */}
      {(() => {
        const cyps = [
          { key: 'cyp1a2_inhibitor', label: 'CYP1A2' },
          { key: 'cyp2c9_inhibitor', label: 'CYP2C9' },
          { key: 'cyp2c19_inhibitor', label: 'CYP2C19' },
          { key: 'cyp2d6_inhibitor', label: 'CYP2D6' },
          { key: 'cyp3a4_inhibitor', label: 'CYP3A4' },
        ].map(c => ({ ...c, value: getCyp(c.key) })).filter(c => c.value != null)
        if (cyps.length === 0) return null
        return (
          <div>
            <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wider mb-1.5">Metabolism — CYP inhibition</p>
            <div className="grid grid-cols-5 gap-1">
              {cyps.map(c => {
                const isInhibitor = c.value > 0.5
                return (
                  <div key={c.key} className={`rounded-lg px-1 py-1.5 border text-center ${isInhibitor ? 'bg-red-50 border-red-200' : 'bg-green-50 border-green-200'}`}>
                    <p className="text-[8px] text-gray-400">{c.label}</p>
                    <p className={`text-[10px] font-bold ${isInhibitor ? 'text-red-600' : 'text-green-600'}`}>
                      {Math.round(c.value * 100)}%
                    </p>
                  </div>
                )
              })}
            </div>
          </div>
        )
      })()}

      {/* ── EXCRETION ── */}
      {mol.half_life != null && (
        <div className="space-y-1.5">
          <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wider">Excretion</p>
          <div className="flex items-center justify-between bg-gray-50 rounded-lg px-2.5 py-1.5 border border-gray-100">
            <span className="text-[10px] text-gray-500">Half-life</span>
            {(() => {
              const hl = typeof mol.half_life === 'number' ? Math.max(0, mol.half_life) : 0
              return (
                <span className={`text-[10px] font-bold tabular-nums ${
                  hl > 12 ? 'text-green-600' : hl > 4 ? 'text-amber-600' : 'text-red-600'
                }`}>{hl.toFixed(1)}<span className="text-[9px] font-normal text-gray-400 ml-0.5">h</span></span>
              )
            })()}
          </div>
        </div>
      )}
    </div>
  )
}

function ToxicityDetailTab({ mol, details }) {
  const safety = details?.safety || (() => {
    const hasAny = mol.herg_risk != null || mol.ames_mutagenicity != null || mol.hepatotoxicity != null ||
                   mol.skin_sensitization != null || mol.carcinogenicity != null || mol.safety_color_code != null
    if (!hasAny) return null
    return {
      herg_risk: mol.herg_risk ?? mol.hERG ?? null,
      ames_risk: mol.ames_mutagenicity ?? null,
      hepatotox_risk: mol.hepatotoxicity ?? null,
    }
  })()

  const hasTox = safety?.herg_risk != null || safety?.ames_risk != null || safety?.hepatotox_risk != null ||
    mol.skin_sensitization != null || mol.carcinogenicity != null
  const hasScore = mol.composite_score != null || mol.safety_color_code != null

  if (!hasTox && !hasScore) {
    return <p className="text-sm text-gray-400 text-center py-4">No toxicity data — run a Toxicity analysis first.</p>
  }

  return (
    <div className="space-y-3">
      {/* ── SAFETY SCORE SUMMARY ── */}
      {hasScore && (
        <div className="flex items-center gap-3 bg-gray-50 rounded-lg px-3 py-2.5 border border-gray-100">
          {mol.safety_color_code && (
            <div className={`w-3 h-3 rounded-full flex-shrink-0 ${
              mol.safety_color_code === 'green' ? 'bg-green-500' :
              mol.safety_color_code === 'yellow' ? 'bg-yellow-500' : 'bg-red-500'
            }`} />
          )}
          <div className="flex-1">
            <p className="text-[10px] font-semibold text-gray-400 uppercase">Safety Score</p>
            {mol.composite_score != null && (
              <p className={`text-lg font-bold tabular-nums ${
                mol.composite_score >= 0.7 ? 'text-green-600' : mol.composite_score >= 0.4 ? 'text-amber-600' : 'text-red-600'
              }`}>{Math.round(mol.composite_score * 100)}<span className="text-sm text-gray-400">/100</span></p>
            )}
          </div>
        </div>
      )}

      {/* ── TOXICITY ENDPOINTS ── */}
      {hasTox && (
        <div className="space-y-1.5">
          <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wider">Toxicity Endpoints</p>
          <RiskBar label="hERG inhibition" value={safety?.herg_risk} invert />
          {/* hERG Specialized: IC50 + risk level when available */}
          {mol.herg_specialized && (
            <div className="flex items-center gap-2 pl-1 -mt-0.5">
              <span className={`inline-flex items-center px-1.5 py-0.5 rounded text-[9px] font-semibold uppercase ${
                mol.herg_specialized.risk_level === 'LOW' ? 'bg-green-100 text-green-700' :
                mol.herg_specialized.risk_level === 'MODERATE' ? 'bg-amber-100 text-amber-700' :
                mol.herg_specialized.risk_level === 'HIGH' ? 'bg-red-100 text-red-700' : 'bg-gray-100 text-gray-500'
              }`}>{mol.herg_specialized.risk_level}</span>
              <span className="text-[10px] text-gray-500">
                IC50 ≈ {Math.round(mol.herg_specialized.ic50_um * 10) / 10} µM
              </span>
            </div>
          )}
          <RiskBar label="Ames mutagenicity" value={safety?.ames_risk} invert />
          <RiskBar label="Hepatotoxicity" value={safety?.hepatotox_risk} invert />
          <RiskBar label="Skin sensitization" value={mol.skin_sensitization} invert />
          <RiskBar label="Carcinogenicity" value={mol.carcinogenicity} invert />
        </div>
      )}
    </div>
  )
}

// Off-Target Tab — separate run type
function OffTargetTab({ mol, details }) {
  const safety = details?.safety
  const hasDetailedResults = safety?.off_target?.length > 0

  // Fallback: show summary from flat molecule data even when per-target results are missing
  const selectivity = mol.selectivity_score
  const hits = mol.off_target_hits
  const ratio = mol.selectivity_ratio

  if (!hasDetailedResults && selectivity == null && hits == null) {
    return <p className="text-sm text-gray-400 text-center py-4">No off-target data — run the Off-Target analysis first.</p>
  }

  const riskColor = r => r === 'high' ? 'text-red-600 bg-red-50 border-red-200' :
                         r === 'medium' ? 'text-amber-600 bg-amber-50 border-amber-200' :
                         'text-green-600 bg-green-50 border-green-200'

  // Summary bar (works with both detailed and flat data)
  const summaryScore = safety?.off_target_summary?.selectivity_score ?? selectivity
  const nSafe = safety?.off_target_summary?.n_safe
  const nTotal = safety?.off_target_summary?.n_total

  return (
    <div className="space-y-3">
      {/* Summary */}
      {(summaryScore != null || hits != null) && (
        <div className="flex items-center justify-between bg-gray-50 rounded-lg px-3 py-2 border border-gray-100">
          <div className="flex items-center gap-2">
            <span className="text-[10px] font-semibold text-gray-400 uppercase">Selectivity</span>
            {summaryScore != null && (
              <span className={`text-sm font-bold tabular-nums ${
                summaryScore >= 0.7 ? 'text-green-600' :
                summaryScore >= 0.4 ? 'text-amber-600' : 'text-red-600'
              }`}>{(summaryScore * 100).toFixed(0)}%</span>
            )}
          </div>
          <div className="flex items-center gap-3">
            {hits != null && (
              <span className="text-[10px] text-gray-500">
                {hits} off-target hit{hits !== 1 ? 's' : ''}
              </span>
            )}
            {nSafe != null && nTotal != null && (
              <span className="text-[10px] text-gray-500">
                {nSafe}/{nTotal} safe
              </span>
            )}
          </div>
        </div>
      )}

      {/* Selectivity ratio */}
      {ratio != null && !hasDetailedResults && (
        <div className="bg-gray-50 rounded-lg px-3 py-2 border border-gray-100">
          <div className="flex items-center justify-between mb-1">
            <span className="text-[10px] text-gray-400 uppercase">Selectivity Ratio</span>
            <span className="text-sm font-bold tabular-nums text-gray-700">{ratio.toFixed(2)}</span>
          </div>
          <div className="h-1.5 bg-gray-200 rounded-full overflow-hidden">
            <div className={`h-1.5 rounded-full ${ratio >= 0.7 ? 'bg-green-500' : ratio >= 0.4 ? 'bg-amber-500' : 'bg-red-500'}`}
              style={{ width: `${Math.min(ratio * 100, 100)}%` }} />
          </div>
        </div>
      )}

      {/* Per-target details (when available) */}
      {hasDetailedResults && (
        <div className="space-y-1.5">
          {safety.off_target.map((ot, i) => (
            <div key={i} className="flex items-center justify-between bg-gray-50 rounded-lg px-2.5 py-1.5">
              <div className="min-w-0 flex-1">
                <p className="text-sm font-medium text-gray-700">{ot.target}</p>
                <p className="text-[10px] text-gray-400 truncate">
                  {ot.risk_description && <span>{ot.risk_description} · </span>}
                  {ot.score != null && <span>{ot.score.toFixed(1)} kcal/mol (threshold {ot.threshold})</span>}
                  {ot.family && <span>{ot.family} · sim. {((ot.similarity || 0) * 100).toFixed(0)}%</span>}
                </p>
              </div>
              <span className={`text-[10px] font-semibold px-2 py-0.5 rounded-full border capitalize ml-2 flex-shrink-0 ${riskColor(ot.risk)}`}>
                {ot.status || ot.risk}
              </span>
            </div>
          ))}
        </div>
      )}

      {/* Warnings */}
      {safety?.off_target_summary?.warnings?.length > 0 && (
        <div className="bg-amber-50 border border-amber-200 rounded-lg px-2.5 py-1.5">
          <p className="text-[10px] font-semibold text-amber-700 mb-0.5">Warnings</p>
          {safety.off_target_summary.warnings.map((w, i) => (
            <p key={i} className="text-[10px] text-amber-600">{w}</p>
          ))}
        </div>
      )}
    </div>
  )
}

function SynthesisTab({ details, onOpenPopup }) {
  const synth = details?.synthesis
  if (!synth) {
    return <p className="text-sm text-gray-400 text-center py-4">No retrosynthesis data for this molecule.</p>
  }

  const isRdkit = synth.method === 'rdkit_heuristic'

  return (
    <div className="space-y-3">
      {/* Method badge + Open full tree button */}
      <div className="flex items-center gap-2">
        {synth.method && (
          <span className={`inline-flex items-center px-2 py-0.5 rounded-full text-[10px] font-semibold ${
            isRdkit
              ? 'bg-amber-100 text-amber-700 border border-amber-200'
              : 'bg-indigo-100 text-indigo-700 border border-indigo-200'
          }`}>
            {isRdkit ? 'RDKit Heuristic' : 'AiZynthFinder'}
          </span>
        )}
        {onOpenPopup && (
          <button
            onClick={onOpenPopup}
            className="flex-1 flex items-center justify-center gap-2 px-4 py-2.5 bg-pink-50 hover:bg-pink-100 border border-pink-200 rounded-xl text-pink-700 font-semibold text-sm transition-colors"
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2" />
            </svg>
            View full retrosynthesis route
          </button>
        )}
      </div>

      {/* Summary */}
      <div className="grid grid-cols-2 gap-2 text-sm">
        {synth.num_steps != null && (
          <div className="bg-gray-50 rounded-lg px-2.5 py-2 text-center border border-gray-100">
            <p className="text-[9px] text-gray-400 uppercase">Steps</p>
            <p className="font-bold text-gray-800 text-base">{synth.num_steps}</p>
          </div>
        )}
        {synth.total_cost != null && (
          <div className="bg-gray-50 rounded-lg px-2.5 py-2 text-center border border-gray-100">
            <p className="text-[9px] text-gray-400 uppercase">Cost</p>
            <p className="font-bold text-gray-800 text-base">${typeof synth.total_cost === 'number' ? synth.total_cost.toFixed(0) : synth.total_cost}</p>
          </div>
        )}
        {synth.feasibility != null && (
          <div className={`rounded-lg px-2.5 py-2 text-center border ${synth.feasibility >= 0.7 ? 'bg-green-50 border-green-200' : synth.feasibility >= 0.5 ? 'bg-amber-50 border-amber-200' : 'bg-red-50 border-red-200'}`}>
            <p className="text-[9px] text-gray-400 uppercase">Feasibility</p>
            <p className={`font-bold text-base ${synth.feasibility >= 0.7 ? 'text-green-600' : synth.feasibility >= 0.5 ? 'text-amber-600' : 'text-red-600'}`}>
              {(synth.feasibility * 100).toFixed(0)}%
          </p>
        </div>
        )}
        {synth.reagents_available != null && (
          <div className={`rounded-lg px-2.5 py-2 text-center border ${synth.reagents_available ? 'bg-green-50 border-green-200' : 'bg-amber-50 border-amber-200'}`}>
            <p className="text-[9px] text-gray-400 uppercase">Reagents</p>
            <p className={`font-bold text-base ${synth.reagents_available ? 'text-green-600' : 'text-amber-600'}`}>
              {synth.reagents_available ? 'Available' : 'Partial'}
            </p>
          </div>
        )}
      </div>

      {/* Reagent availability detail — hidden for RDKit (no ZINC data) */}
      {!isRdkit && synth.reagent_availability?.length > 0 && (
        <div className="space-y-1">
          <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wider">Reagent availability</p>
          {synth.reagent_availability.map((r, i) => (
            <div key={i} className="flex items-center justify-between bg-gray-50 rounded-lg px-2.5 py-1.5 text-[11px]">
              <span className="text-gray-700 font-mono break-all">{r.name || r.smiles || `Reagent ${i + 1}`}</span>
              <span className={`font-semibold flex-shrink-0 ml-2 ${r.available ? 'text-green-600' : 'text-red-500'}`}>
                {r.available ? 'Yes' : 'No'}
              </span>
            </div>
          ))}
        </div>
      )}

      {/* Steps */}
      <div className="space-y-2">
        <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wider">Retrosynthesis route</p>
        {synth.steps.map((step, i) => {
          const conf = step.confidence ?? step.yield
          const reaction = step.reaction || step.product || `Step ${i + 1}`
          const reagents = step.reactant_names?.join(', ') || step.reactants?.join(', ') || step.reagent || ''
          return (
            <div key={i} className="flex items-start gap-2.5">
              <div className="flex-shrink-0 w-5 h-5 rounded-full bg-bx-surface text-white text-[9px] font-bold flex items-center justify-center mt-0.5">
                {i + 1}
              </div>
              <div className="flex-1 bg-gray-50 rounded-lg p-2.5 border border-gray-100">
                <div className="flex items-start justify-between gap-2">
                  <div className="min-w-0 flex-1">
                    <p className="text-sm font-semibold text-gray-700">{reaction}</p>
                    {reagents && <p className="text-[10px] text-gray-500 mt-0.5 break-all">{reagents}</p>}
                    {step.conditions && <p className="text-[9px] text-gray-400 italic mt-0.5">{step.conditions}</p>}
                  </div>
                  {conf != null && (
                    <div className="text-right flex-shrink-0">
                      <p className="text-sm font-bold text-bx-mint tabular-nums">{(conf * 100).toFixed(0)}%</p>
                      <p className="text-[9px] text-gray-400">confidence</p>
                    </div>
                  )}
                </div>
                {conf != null && (
                  <div className="mt-2 w-full bg-gray-200 rounded-full h-1">
                    <div className="h-1 bg-bx-mint rounded-full" style={{ width: `${conf * 100}%` }} />
                  </div>
                )}
              </div>
            </div>
          )
        })}
      </div>
    </div>
  )
}

// Interaction type badge colors
const INTER_TYPE_STYLES = {
  Hbond:         { bg: 'bg-blue-100', text: 'text-blue-700' },
  HBDonor:       { bg: 'bg-blue-100', text: 'text-blue-700' },
  HBAcceptor:    { bg: 'bg-sky-100', text: 'text-sky-700' },
  Hydrophobic:   { bg: 'bg-gray-100', text: 'text-gray-600' },
  PiStacking:    { bg: 'bg-purple-100', text: 'text-purple-700' },
  CationPi:      { bg: 'bg-orange-100', text: 'text-orange-700' },
  PiCation:      { bg: 'bg-orange-100', text: 'text-orange-700' },
  SaltBridge:    { bg: 'bg-teal-100', text: 'text-teal-700' },
  Ionic:         { bg: 'bg-teal-100', text: 'text-teal-700' },
  Cationic:      { bg: 'bg-teal-100', text: 'text-teal-700' },
  Anionic:       { bg: 'bg-teal-100', text: 'text-teal-700' },
  VDW:           { bg: 'bg-gray-50', text: 'text-gray-400' },
  VdWContact:    { bg: 'bg-gray-50', text: 'text-gray-400' },
  MetalAcceptor: { bg: 'bg-violet-100', text: 'text-violet-700' },
}

const ORIGIN_ROLES = new Set(['active_site', 'binding', 'metal_binding', 'custom', 'key'])

function InteractionsTab({ molecule, details, onHighlightResidue, onShowInteractions, functionalResidues, computedDistances }) {
  const [viewMode, setViewMode] = useState('diagram')

  const nInteractions = molecule?.n_interactions
  const functionalContacts = molecule?.functional_contacts
  const totalFunctional = molecule?.total_functional || 0
  const interactionQuality = molecule?.interaction_quality
  const keyHbonds = molecule?.key_hbonds
  const method = molecule?.interactions_method
  const rawInterDetail = molecule?.interactions_detail || []

  // Enrich with 3D-computed distances when backend data lacks them
  const interDetail = useMemo(() => {
    if (!computedDistances || Object.keys(computedDistances).length === 0) return rawInterDetail
    return rawInterDetail.map(row => {
      if (row.distance != null) return row
      const cd = computedDistances[row.residue_number]
      return cd != null ? { ...row, distance: cd } : row
    })
  }, [rawInterDetail, computedDistances])

  if (nInteractions == null && interDetail.length === 0) {
    return <p className="text-sm text-gray-400 text-center py-4">No interaction data — run an Interactions analysis first.</p>
  }

  // Build the interaction data object for InteractionDiagram / InteractionView
  const interactionData = {
    interactions: interDetail,
    functional_contacts: functionalContacts,
    total_functional: totalFunctional,
    interaction_quality: interactionQuality,
    key_hbonds: keyHbonds,
    method: method,
  }

  // Group interactions by type for badges
  const typeGroups = {}
  interDetail.forEach(row => {
    const t = row.type || 'VDW'
    if (!typeGroups[t]) typeGroups[t] = { total: 0, functional: 0 }
    typeGroups[t].total++
    if (row.is_functional) typeGroups[t].functional++
  })

  // Missing functional contacts: functional residues NOT contacted
  const contactedFuncNums = new Set(
    interDetail.filter(r => r.is_functional && r.residue_number != null).map(r => r.residue_number)
  )

  // Use mapped_functional from backend (PDB-numbered, offset-corrected) when available
  const mappedFunc = molecule?.mapped_functional || []

  const missingFunc = useMemo(() => {
    if (mappedFunc.length > 0) {
      return mappedFunc.filter(mf => mf.number != null && !contactedFuncNums.has(mf.number))
    }
    // Fallback for old runs without mapped_functional
    return (functionalResidues || []).filter(
      fr => fr.enabled !== false && fr.number != null
        && ORIGIN_ROLES.has(fr.role || 'key')
        && !contactedFuncNums.has(fr.number)
    )
  }, [mappedFunc, functionalResidues, contactedFuncNums])

  // Build residue origins map for the Origin column
  // Uses PDB-numbered data from backend to match interaction table residue numbers
  const residueOrigins = useMemo(() => {
    const map = {}
    if (mappedFunc.length > 0) {
      mappedFunc.forEach(mf => {
        if (mf.number != null) {
          map[mf.number] = { role: mf.role || 'key', label: mf.description || mf.role || 'Key' }
        }
      })
    } else {
      // Fallback for old runs: use raw functionalResidues (assumes offset=0)
      ;(functionalResidues || []).forEach(fr => {
        if (fr.number != null && ORIGIN_ROLES.has(fr.role || 'key')) {
          map[fr.number] = { role: fr.role || 'key', label: fr.description || fr.role || 'Key' }
        }
      })
    }
    return map
  }, [mappedFunc, functionalResidues])

  const handleResidueClick = (resNum, resName) => {
    if (onHighlightResidue && resNum != null) {
      onHighlightResidue({ number: resNum, aa: resName, isMissing: false })
    }
  }

  return (
    <div className="space-y-3">
      {/* Summary metrics */}
      <div className="flex flex-wrap gap-2">
        {[
          { label: 'Contacts', value: nInteractions, color: 'text-gray-700' },
          totalFunctional > 0
            ? { label: 'Functional', value: `${functionalContacts ?? 0}/${totalFunctional}`, color: 'text-lime-600' }
            : { label: 'Functional', value: functionalContacts, color: 'text-lime-600' },
          interactionQuality != null ? { label: 'Quality', value: (interactionQuality * 100).toFixed(0) + '%', color: 'text-emerald-600' } : null,
          { label: 'Key HBonds', value: keyHbonds, color: 'text-blue-600' },
        ].filter(s => s && s.value != null).map(s => (
          <div key={s.label} className="bg-gray-50 rounded-lg px-2.5 py-1.5 text-center border border-gray-100 min-w-[52px]">
            <p className={`text-base font-bold tabular-nums ${s.color}`}>{s.value}</p>
            <p className="text-[9px] text-gray-400">{s.label}</p>
          </div>
        ))}
      </div>

      {/* Type badges */}
      {Object.keys(typeGroups).length > 0 && (
        <div className="flex flex-wrap gap-1.5">
          {Object.entries(typeGroups).map(([type, counts]) => {
            const st = INTER_TYPE_STYLES[type] || { bg: 'bg-gray-100', text: 'text-gray-600' }
            return (
              <span key={type} className={`inline-flex items-center gap-1 px-2 py-0.5 rounded-full text-[10px] font-semibold ${st.bg} ${st.text}`}>
                {({ Hbond:'H-Bond', HBDonor:'H-Bond Donor', HBAcceptor:'H-Bond Acceptor', PiStacking:'Pi Stacking', CationPi:'Cation-Pi', PiCation:'Cation-Pi', SaltBridge:'Salt Bridge', Ionic:'Ionic', Cationic:'Ionic (+)', Anionic:'Ionic (-)', VDW:'VdW', VdWContact:'VdW', MetalAcceptor:'Metal Coord.' })[type] || type}
                <span className="font-bold">{counts.total}</span>
                {counts.functional > 0 && (
                  <span className="text-green-600 ml-0.5">({counts.functional} func.)</span>
                )}
              </span>
            )
          })}
        </div>
      )}

      {/* Show on 3D button */}
      {interDetail.length > 0 && onShowInteractions && (
        <button
          onClick={() => onShowInteractions(interDetail)}
          className="w-full flex items-center justify-center gap-2 px-3 py-2 bg-lime-50 hover:bg-lime-100 border border-lime-200 rounded-xl text-lime-700 font-semibold text-sm transition-colors"
        >
          <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
              d="M14 10l-2 1m0 0l-2-1m2 1v2.5M20 7l-2 1m2-1l-2-1m2 1v2.5M14 4l-2-1-2 1M4 7l2-1M4 7l2 1M4 7v2.5M12 21l-2-1m2 1l2-1m-2 1v-2.5M6 18l-2-1v-2.5M18 18l2-1v-2.5" />
          </svg>
          Show on 3D
        </button>
      )}

      {/* Missing functional contacts */}
      {missingFunc.length > 0 && (
        <div className="p-2.5 bg-red-50 rounded-lg border border-red-200">
          <p className="text-[10px] font-semibold text-red-600 uppercase tracking-wider mb-1.5">Missing functional contacts</p>
          <div className="flex flex-wrap gap-1">
            {missingFunc.map(fr => (
              <button
                key={fr.number}
                onClick={() => onHighlightResidue && onHighlightResidue({ number: fr.number, aa: fr.aa, isMissing: true })}
                className="inline-flex items-center gap-1 px-2 py-0.5 rounded-full text-[10px] font-semibold bg-red-100 text-red-700 hover:bg-red-200 transition-colors cursor-pointer"
              >
                {fr.aa || '?'}{fr.number}
              </button>
            ))}
          </div>
        </div>
      )}

      {/* View mode toggle */}
      <div className="flex items-center gap-1 bg-gray-100 rounded-lg p-0.5 w-fit">
        <button
          onClick={() => setViewMode('diagram')}
          className={`px-3 py-1 text-xs font-medium rounded-md transition-colors ${
            viewMode === 'diagram' ? 'bg-white text-bx-light-text shadow-sm' : 'text-gray-500 hover:text-gray-700'
          }`}
        >
          Diagram
        </button>
        <button
          onClick={() => setViewMode('table')}
          className={`px-3 py-1 text-xs font-medium rounded-md transition-colors ${
            viewMode === 'table' ? 'bg-white text-bx-light-text shadow-sm' : 'text-gray-500 hover:text-gray-700'
          }`}
        >
          Table
        </button>
      </div>

      {/* Interaction visualization */}
      {viewMode === 'diagram' ? (
        <InteractionDiagram interactions={interactionData} />
      ) : (
        <InteractionView
          interactions={interactionData}
          onResidueClick={onHighlightResidue ? handleResidueClick : undefined}
          residueOrigins={residueOrigins}
        />
      )}
    </div>
  )
}

function ScaffoldTab({ molecule }) {
  const scaffoldSmiles = molecule?.scaffold_smiles
  const nPositions = molecule?.n_modifiable_positions
  const bricsBonds = molecule?.brics_bond_count
  const nRings = molecule?.scaffold_n_rings
  const scaffoldMw = molecule?.scaffold_mw
  const positions = molecule?.scaffold_positions || []
  const svgContent = molecule?.scaffold_svg

  if (!scaffoldSmiles && nPositions == null) {
    return <p className="text-sm text-gray-400 text-center py-4">No structural data — run a Structural Analysis first.</p>
  }

  return (
    <div className="space-y-3">
      {/* Summary stats */}
      <div className="flex flex-wrap gap-2">
        {[
          nPositions != null && { label: 'R-Positions', value: nPositions, color: 'text-emerald-600' },
          bricsBonds != null && { label: 'BRICS Bonds', value: bricsBonds, color: 'text-teal-600' },
          nRings != null && { label: 'Rings', value: nRings, color: 'text-blue-600' },
          scaffoldMw != null && { label: 'MW', value: Number(scaffoldMw).toFixed(0), color: 'text-gray-600' },
        ].filter(Boolean).map(s => (
          <div key={s.label} className="bg-gray-50 rounded-lg px-2.5 py-1.5 text-center border border-gray-100 min-w-[52px]">
            <p className={`text-base font-bold tabular-nums ${s.color}`}>{s.value}</p>
            <p className="text-[9px] text-gray-400">{s.label}</p>
          </div>
        ))}
      </div>

      {/* Scaffold SMILES */}
      {scaffoldSmiles && (
        <div className="bg-gray-50 rounded-xl border border-gray-100 p-3">
          <p className="text-[9px] font-semibold text-gray-400 uppercase tracking-wider mb-1">Scaffold SMILES</p>
          <p className="text-xs font-mono text-gray-600 break-all">{scaffoldSmiles}</p>
        </div>
      )}

      {/* SVG visualization */}
      {svgContent && (
        <div className="bg-white rounded-xl border border-gray-100 p-2">
          <p className="text-[9px] font-semibold text-gray-400 uppercase tracking-wider mb-2">Annotated Structure</p>
          <div className="flex justify-center" dangerouslySetInnerHTML={{ __html: svgContent }} />
        </div>
      )}

      {/* Modifiable positions */}
      {positions.length > 0 && (
        <div className="space-y-1.5">
          <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wider">Modifiable Positions</p>
          {positions.map((pos, i) => (
            <div key={i} className="flex items-center gap-2.5 py-1.5 border-b border-gray-50 last:border-0">
              <span className={`text-[9px] font-semibold px-1.5 py-0.5 rounded-full flex-shrink-0 ${
                pos.type === 'R-group' ? 'bg-green-100 text-green-700'
                : pos.type === 'BRICS' ? 'bg-red-100 text-red-700'
                : 'bg-gray-100 text-gray-600'
              }`}>
                {pos.type || 'R-group'}
              </span>
              <span className="text-sm text-gray-700">{pos.label || pos.smiles || `Position ${i + 1}`}</span>
              {pos.atom_idx != null && (
                <span className="text-[9px] text-gray-400 ml-auto">atom {pos.atom_idx}</span>
              )}
            </div>
          ))}
        </div>
      )}

      {/* Diversity Cluster */}
      {molecule?.cluster_id != null && (
        <div className="bg-gray-50 rounded-xl border border-gray-100 p-3">
          <p className="text-[9px] font-semibold text-gray-400 uppercase tracking-wider mb-2">Diversity Cluster</p>
          <div className="flex flex-wrap gap-2">
            {[
              { label: 'Cluster', value: molecule.cluster_id },
              molecule.cluster_size != null && { label: 'Size', value: molecule.cluster_size },
              molecule.tanimoto_to_centroid != null && { label: 'Tanimoto', value: Number(molecule.tanimoto_to_centroid).toFixed(3) },
              molecule.is_representative && { label: 'Repr.', value: 'Yes' },
            ].filter(Boolean).map(s => (
              <div key={s.label} className="bg-white rounded-lg px-2.5 py-1.5 text-center border border-gray-100 min-w-[52px]">
                <p className="text-base font-bold tabular-nums text-teal-600">{s.value}</p>
                <p className="text-[9px] text-gray-400">{s.label}</p>
              </div>
            ))}
          </div>
        </div>
      )}
    </div>
  )
}

// ---------------------------------------------------------------------------
// Detail Popup Modal — Safety / Synthesis / Confidence
// ---------------------------------------------------------------------------
function DetailPopupModal({ type, molecule, onClose }) {
  if (!molecule) return null

  const props = molecule.properties || {}
  const flat = molecule // already flattened

  const titleMap = {
    safety: 'Safety Report',
    retrosynthesis: 'Retrosynthesis Analysis',
    confidence: 'Confidence Breakdown',
  }

  const content = (() => {
    switch (type) {
      case 'safety':
        return (
          <SafetyReport
            offTargetResults={props.safety?.off_target || props.off_target || null}
            moleculeName={molecule.name || molecule.id}
          />
        )
      case 'retrosynthesis':
        return (
          <SynthesisTree
            synthesisRoute={props.retrosynthesis || null}
          />
        )
      case 'confidence': {
        const confidence = props.confidence || null
        return (
          <ConfidenceBreakdown
            confidence={confidence}
            moleculeName={molecule.name || molecule.id}
            pipeline_summary={props.pipeline_summary || null}
          />
        )
      }
      default:
        return <p className="text-gray-400 text-sm p-6">Unknown popup type: {type}</p>
    }
  })()

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center p-4" onClick={onClose}>
      <div className="absolute inset-0 bg-black/40 backdrop-blur-sm" />
      <div
        className="relative bg-white rounded-2xl shadow-2xl border border-gray-200 max-w-2xl w-full max-h-[85vh] overflow-hidden flex flex-col animate-in"
        onClick={e => e.stopPropagation()}
      >
        {/* Header */}
        <div className="flex items-center justify-between px-6 py-4 border-b border-gray-100 bg-gradient-to-r from-gray-50 to-white flex-shrink-0">
          <div className="flex items-center gap-3">
            <h3 className="text-base font-bold text-gray-800">{titleMap[type] || 'Details'}</h3>
            <span className="text-sm text-gray-400">{molecule.name || molecule.id}</span>
          </div>
          <button
            onClick={onClose}
            className="p-1.5 rounded-lg hover:bg-gray-100 text-gray-400 hover:text-gray-700 transition-colors"
          >
            <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
            </svg>
          </button>
        </div>

        {/* Content */}
        <div className="overflow-y-auto flex-1 p-6" style={{ scrollbarWidth: 'thin' }}>
          {content}
        </div>
      </div>
    </div>
  )
}


// ---------------------------------------------------------------------------
// Molecule Detail Panel (inline drawer)
// ---------------------------------------------------------------------------
const TABS = [
  { id: 'scores',       label: 'Docking',      tip: 'Docking scores from computational screening' },
  { id: 'composite',    label: 'Score',        tip: 'Weighted multi-criteria score combining docking, ADME, selectivity, and drug-likeness' },
  { id: 'properties',   label: 'Props',       tip: 'Physicochemical properties: MW, logP, TPSA, Lipinski Ro5' },
  { id: 'adme',         label: 'ADME',        tip: 'Pharmacokinetics: absorption, distribution, metabolism, excretion' },
  { id: 'toxicity',     label: 'Toxicity',    tip: 'Safety endpoints: hERG, Ames, hepatotoxicity, carcinogenicity' },
  { id: 'offtarget',    label: 'Off-Target',  tip: 'Selectivity profiling against off-target proteins' },
  { id: 'synthesis',    label: 'Synth.',      tip: 'Retrosynthesis feasibility, estimated cost, and synthetic route' },
  { id: 'scaffold',     label: 'Structure',   tip: 'Murcko scaffold decomposition, BRICS bonds, R-group positions, and diversity clustering' },
]

function MoleculeDetailPanel({ molecule, molecules, onClose, onToggleBookmark, onRowClick, isFrozen, project, campaign, selectedMolecules, onCellPopup, panelPosition, onTogglePosition }) {
  const [activeTab, setActiveTab] = useState('scores')
  const [copied, setCopied] = useState(false)
  const [viewerMode, setViewerMode] = useState(null) // null = auto-detect: 'protein' | 'ligand' | 'results'
  const [highlightedInteraction, setHighlightedInteraction] = useState(null)
  const [interactionContacts, setInteractionContacts] = useState(null)
  const [computedDistances, setComputedDistances] = useState({}) // { residue_number: distance }
  const [panelFullscreen, setPanelFullscreen] = useState(false)

  // Escape key exits panel fullscreen
  useEffect(() => {
    if (!panelFullscreen) return
    const handleKey = (e) => { if (e.key === 'Escape') setPanelFullscreen(false) }
    window.addEventListener('keydown', handleKey)
    return () => window.removeEventListener('keydown', handleKey)
  }, [panelFullscreen])

  // Auto-update interaction contacts when molecule changes (if interactions are shown)
  const molId = molecule?.id
  useEffect(() => {
    if (interactionContacts) {
      const newDetail = molecule?.interactions_detail || []
      if (newDetail.length > 0) {
        setInteractionContacts(newDetail)
      } else {
        setInteractionContacts(null)
      }
      setHighlightedInteraction(null)
      setComputedDistances({})
    }
  }, [molId]) // eslint-disable-line react-hooks/exhaustive-deps

  // Extract target structure info from project
  const targetPreview = project?.target_preview || {}
  const funcResidues = targetPreview.functional_residues?.residues || []
  const pdbUrl = targetPreview.structure?.download_url || null
  const selectedPocket = useMemo(() => {
    const pockets = targetPreview.pockets || []
    const idx = targetPreview.selected_pocket_index
    return idx != null && pockets[idx] ? pockets[idx] : null
  }, [targetPreview])

  // Build docked molblocks: Ctrl+click highlights → overlay multiple poses
  const dockedMolblocks = useMemo(() => {
    if (!selectedMolecules || selectedMolecules.length === 0) return null
    // Merge highlighted molecules + current active molecule (avoid duplicates)
    const allMols = [...selectedMolecules]
    if (molecule && !allMols.find(m => m.id === molecule.id)) {
      allMols.unshift(molecule)
    }
    if (allMols.length <= 1) return null
    const blocks = allMols
      .map(m => m.properties?.docking?.pose_molblock)
      .filter(Boolean)
    return blocks.length > 1 ? blocks.map(mb => ({ molblock: mb })) : null
  }, [selectedMolecules, molecule])

  const currentIdx = molecules.findIndex(m => m.id === molecule.id)
  const prevMol = currentIdx > 0 ? molecules[currentIdx - 1] : null
  const nextMol = currentIdx < molecules.length - 1 ? molecules[currentIdx + 1] : null

  // Build details from molecule properties (API data, not mock)
  const details = useMemo(() => {
    const props = molecule.properties || {}
    const retro = props.retrosynthesis
    const ot = props.off_target
    // Build safety section: merge safety run + off_target run + druglikeness_rules data
    const safetyRun = props.safety || {}
    const dlRules = props.druglikeness_rules || {}
    const confRun = props.confidence || {}
    const safetyObj = {
      herg_risk: safetyRun.herg_risk ?? molecule.herg_risk ?? null,
      ames_risk: safetyRun.ames_mutagenicity ?? molecule.ames_mutagenicity ?? null,
      hepatotox_risk: safetyRun.hepatotoxicity ?? molecule.hepatotoxicity ?? null,
      pains_pass: confRun.pains_alert != null ? !confRun.pains_alert
                : safetyRun.pains_alert != null ? !safetyRun.pains_alert
                : molecule.pains_alert != null ? !molecule.pains_alert : null,
      pains_alerts: safetyRun.pains_alerts || [],
      brenk_alert: dlRules.brenk_alert ?? molecule.brenk_alert ?? null,
      pfizer_alert: dlRules.pfizer_alert ?? molecule.pfizer_alert ?? null,
      gsk_alert: dlRules.gsk_alert ?? molecule.gsk_alert ?? null,
      cns_mpo: dlRules.cns_mpo ?? molecule.cns_mpo ?? null,
    }
    // Inject off-target panel results (from off_target run)
    if (ot?.results) {
      safetyObj.off_target = Object.entries(ot.results).map(([name, data]) => ({
        target: name,
        score: data.score,
        threshold: data.threshold,
        status: data.status,
        risk: data.status === 'risk' ? 'high' : 'low',
        risk_description: data.risk_description,
      }))
      safetyObj.off_target_summary = {
        selectivity_score: ot.selectivity_score,
        n_safe: ot.n_safe,
        n_total: ot.n_total,
        warnings: ot.warnings || [],
      }
    }
    const hasSafety = safetyObj.herg_risk != null || safetyObj.ames_risk != null ||
                      safetyObj.hepatotox_risk != null || safetyObj.pains_pass != null ||
                      safetyObj.brenk_alert != null || safetyObj.pfizer_alert != null ||
                      safetyObj.gsk_alert != null || safetyObj.off_target?.length > 0
    return {
      interactions: props.interactions || props.enrichment || null,
      synthesis: retro ? {
        feasibility: retro.confidence ?? retro.synth_confidence ?? null,
        total_cost: retro.estimated_cost ?? retro.synth_cost_estimate ?? retro.cost_estimate ?? null,
        num_steps: retro.n_steps ?? retro.n_synth_steps ?? null,
        steps: retro.steps || [],
        reagents_available: retro.all_reagents_available ?? retro.reagents_available ?? null,
        reagent_availability: retro.reagent_availability || [],
        method: retro.method || null,
      } : null,
      safety: hasSafety ? safetyObj : null,
    }
  }, [molecule])

  function copySmiles() {
    if (!molecule.smiles) return
    navigator.clipboard.writeText(molecule.smiles).then(() => {
      setCopied(true)
      setTimeout(() => setCopied(false), 1500)
    })
  }

  const panelContent = (
    <div className={
      panelFullscreen
        ? "fixed inset-0 z-50 bg-white flex flex-col"
        : panelPosition === 'bottom'
          ? "bg-white rounded-xl border border-gray-200 shadow-xl overflow-hidden flex flex-col"
          : "bg-white rounded-xl border border-gray-200 shadow-xl overflow-hidden flex flex-col sticky top-4"
    }
      style={panelFullscreen ? undefined : panelPosition === 'bottom' ? { maxHeight: '80vh' } : { maxHeight: 'calc(100vh - 120px)' }}>

      {/* Panel header */}
      <div className="flex items-center justify-between px-4 py-3 border-b border-gray-100 bg-gradient-to-r from-gray-50 to-white flex-shrink-0">
        <div className="flex items-center gap-2 min-w-0 flex-1">
          <h3 className="font-bold text-bx-light-text text-sm truncate">
            {molecule.name || molecule.id}
          </h3>
          {!isFrozen && (
            <button
              onClick={() => onToggleBookmark && onToggleBookmark(molecule.id)}
              className="flex-shrink-0 p-0.5 rounded transition-colors hover:bg-yellow-50"
              aria-label={molecule.bookmarked ? 'Remove bookmark' : 'Bookmark'}
            >
              <svg
                className={`w-4 h-4 transition-colors ${molecule.bookmarked ? 'text-yellow-400 fill-yellow-400' : 'text-gray-300 fill-none hover:text-yellow-300'}`}
                stroke="currentColor" viewBox="0 0 24 24"
              >
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                  d="M5 5a2 2 0 012-2h10a2 2 0 012 2v16l-7-3.5L5 21V5z" />
              </svg>
            </button>
          )}
          {molecule.generation_level > 0 && (
            <Badge variant="purple" size="sm">Gen {molecule.generation_level}</Badge>
          )}
          {molecule.parent_molecule_id && molecules && (() => {
            const parent = molecules.find(m => m.id === molecule.parent_molecule_id)
            if (!parent) return null
            return (
              <button
                onClick={() => onRowClick(parent)}
                className="text-[10px] text-purple-500 hover:text-purple-700 transition-colors"
                title={`Parent: ${parent.name || parent.id}`}
              >
                ← {parent.name || 'Parent'}
              </button>
            )
          })()}
          {molecule.bookmarked && (
            <Badge variant="yellow" size="sm">Bookmarked</Badge>
          )}
        </div>
        <div className="flex items-center gap-1 flex-shrink-0 ml-2">
          <button
            onClick={() => setPanelFullscreen(f => !f)}
            className="p-1.5 rounded-lg hover:bg-gray-100 text-gray-400 hover:text-gray-700 transition-colors"
            aria-label={panelFullscreen ? 'Exit fullscreen' : 'Fullscreen panel'}
            title={panelFullscreen ? 'Exit fullscreen (Esc)' : 'Fullscreen panel'}
          >
            {panelFullscreen ? (
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                  d="M9 9L4 4m0 0v4m0-4h4m6 6l5 5m0 0v-4m0 4h-4M9 15l-5 5m0 0h4m-4 0v-4m11-6l5-5m0 0h-4m4 0v4" />
              </svg>
            ) : (
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                  d="M4 8V4m0 0h4M4 4l5 5m11-1V4m0 0h-4m4 0l-5 5M4 16v4m0 0h4m-4 0l5-5m11 5v-4m0 4h-4m4 0l-5-5" />
              </svg>
            )}
          </button>
          {onTogglePosition && !panelFullscreen && (
            <button
              onClick={onTogglePosition}
              className="p-1.5 rounded-lg hover:bg-gray-100 text-gray-400 hover:text-gray-700 transition-colors"
              aria-label={panelPosition === 'bottom' ? 'Move panel to right' : 'Move panel below'}
              title={panelPosition === 'bottom' ? 'Panel to right' : 'Panel below'}
            >
              {panelPosition === 'bottom' ? (
                <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                    d="M9 3h6m-6 0v18m0-18H5a2 2 0 00-2 2v14a2 2 0 002 2h4m6-18v18m0-18h4a2 2 0 012 2v14a2 2 0 01-2 2h-4" />
                </svg>
              ) : (
                <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                    d="M3 9h18M3 9V5a2 2 0 012-2h14a2 2 0 012 2v4M3 9v10a2 2 0 002 2h14a2 2 0 002-2V9" />
                </svg>
              )}
            </button>
          )}
          <button
            onClick={() => { setPanelFullscreen(false); onClose() }}
            className="p-1.5 rounded-lg hover:bg-gray-100 text-gray-400 hover:text-gray-700 transition-colors"
            aria-label="Close detail panel"
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
            </svg>
          </button>
        </div>
      </div>

      {/* Scrollable content */}
      <div className="overflow-y-auto flex-1 p-4"
        style={{ scrollbarWidth: 'thin', scrollbarColor: '#e5e7eb transparent' }}>
        <PanelErrorBoundary resetKey={molId}>

        {/* View toggle: Protein / Ligand / Results */}
        {(() => {
          const hasProtein = !!pdbUrl
          const hasLigand = !!molecule.smiles
          const effectiveMode = viewerMode || (hasProtein ? 'protein' : hasLigand ? 'ligand' : 'results')
          const isBottomLayout = panelPosition === 'bottom' && !panelFullscreen

          const viewButtons = [
            hasProtein && { id: 'protein', label: 'Protein', icon: (
              <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                  d="M14 10l-2 1m0 0l-2-1m2 1v2.5M20 7l-2 1m2-1l-2-1m2 1v2.5M14 4l-2-1-2 1M4 7l2-1M4 7l2 1M4 7v2.5M12 21l-2-1m2 1l2-1m-2 1v-2.5M6 18l-2-1v-2.5M18 18l2-1v-2.5" />
              </svg>
            )},
            hasLigand && { id: 'ligand', label: 'Ligand', icon: (
              <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                  d="M9.75 3.104v5.714a2.25 2.25 0 01-.659 1.591L5 14.5M9.75 3.104c-.251.023-.501.05-.75.082m.75-.082a24.301 24.301 0 014.5 0m0 0v5.714a2.25 2.25 0 00.659 1.591L19 14.5M14.25 3.104c.251.023.501.05.75.082M19 14.5l-2.47 2.47a3.749 3.749 0 01-5.06 0L9 14.5m10 0a4.5 4.5 0 01-9 0" />
              </svg>
            )},
            { id: 'results', label: 'Results', icon: (
              <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                  d="M3 10h18M3 14h18m-9-4v8m-7 0h14a2 2 0 002-2V8a2 2 0 00-2-2H5a2 2 0 00-2 2v8a2 2 0 002 2z" />
              </svg>
            )},
          ].filter(Boolean)

          const toggleBar = (
            <div className="flex gap-1 mb-2">
              {viewButtons.map(btn => (
                <button
                  key={btn.id}
                  onClick={() => setViewerMode(btn.id)}
                  className={`flex items-center gap-1.5 px-3 py-1.5 rounded-lg text-[11px] font-semibold transition-all ${
                    effectiveMode === btn.id
                      ? 'bg-bx-surface text-white shadow-sm'
                      : 'bg-gray-100 text-gray-500 hover:bg-gray-200'
                  }`}
                >
                  {btn.icon}
                  {btn.label}
                </button>
              ))}
            </div>
          )

          // --- Viewer block (3D) ---
          const viewerBlock = (
            <>
              {effectiveMode === 'protein' && hasProtein && (
                <ProteinViewer
                  pdbUrl={pdbUrl}
                  selectedPocket={selectedPocket}
                  height={panelFullscreen ? 500 : isBottomLayout ? 380 : 280}
                  dockedMolblock={!dockedMolblocks ? (molecule.properties?.docking?.pose_molblock || null) : null}
                  dockedMolblocks={dockedMolblocks}
                  highlightedResidue={highlightedInteraction}
                  interactionContacts={interactionContacts}
                  onEnrichedContacts={setComputedDistances}
                />
              )}
              {effectiveMode === 'ligand' && hasLigand && (
                <LigandViewer3D smiles={molecule.smiles} height={panelFullscreen ? 500 : isBottomLayout ? 380 : 280} />
              )}
            </>
          )

          // --- Side content block (interactions, tabs, SMILES) ---
          const sideContent = (
            <div className="space-y-3">
              {effectiveMode === 'protein' && hasProtein && (
                <InteractionsTab
                  molecule={molecule}
                  details={details}
                  onHighlightResidue={(val) => {
                    setHighlightedInteraction(prev =>
                      prev && prev.number === val?.number ? null : val
                    )
                  }}
                  onShowInteractions={(contacts) => {
                    setInteractionContacts(contacts)
                  }}
                  functionalResidues={funcResidues}
                  computedDistances={computedDistances}
                />
              )}
              {effectiveMode === 'results' && (
                <div>
                  <div className="flex border-b border-gray-100 gap-0.5 flex-wrap">
                    {TABS.map(tab => (
                      <button
                        key={tab.id}
                        onClick={() => setActiveTab(tab.id)}
                        className={`px-2 py-1.5 text-[11px] font-semibold border-b-2 transition-all duration-150 whitespace-nowrap inline-flex flex-col items-center gap-0 ${
                          activeTab === tab.id
                            ? 'border-bx-surface text-bx-light-text'
                            : 'border-transparent text-gray-400 hover:text-gray-600 hover:border-gray-200'
                        }`}
                      >
                        <span>{tab.label}</span>
                        {tab.tip && <InfoTip text={tab.tip} size="xs" />}
                      </button>
                    ))}
                  </div>
                  <div className="pt-3">
                    {activeTab === 'scores' && <ScoresTab mol={molecule} />}
                    {activeTab === 'composite' && <CompositeScoreTab mol={molecule} />}
                    {activeTab === 'properties' && <PropertiesTab mol={molecule} />}
                    {activeTab === 'adme' && <AdmeTab mol={molecule} details={details} />}
                    {activeTab === 'toxicity' && <ToxicityDetailTab mol={molecule} details={details} />}
                    {activeTab === 'offtarget' && <OffTargetTab mol={molecule} details={details} />}
                    {activeTab === 'synthesis' && <SynthesisTab details={details} onOpenPopup={onCellPopup ? () => onCellPopup('retrosynthesis', molecule) : null} />}
                    {activeTab === 'scaffold' && <ScaffoldTab molecule={molecule} />}
                  </div>
                </div>
              )}
              {/* SMILES block */}
              {molecule.smiles && (
                <div className="bg-gray-50 rounded-xl border border-gray-100 p-3">
                  <div className="flex items-center justify-between mb-1.5">
                    <p className="text-[9px] font-semibold text-gray-400 uppercase tracking-wider">SMILES</p>
                    <button
                      onClick={copySmiles}
                      className="flex items-center gap-1 text-[9px] text-gray-400 hover:text-bx-light-text transition-colors"
                    >
                      {copied ? (
                        <>
                          <svg className="w-3 h-3 text-green-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M5 13l4 4L19 7" />
                          </svg>
                          Copied
                        </>
                      ) : (
                        <>
                          <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                              d="M8 16H6a2 2 0 01-2-2V6a2 2 0 012-2h8a2 2 0 012 2v2m-6 12h8a2 2 0 002-2v-8a2 2 0 00-2-2h-8a2 2 0 00-2 2v8a2 2 0 002 2z" />
                          </svg>
                          Copy
                        </>
                      )}
                    </button>
                  </div>
                  <p className="font-mono text-[9px] text-gray-500 break-all leading-relaxed">{molecule.smiles}</p>
                </div>
              )}
            </div>
          )

          // --- Layout: 2-col (bottom) or 1-col (right/fullscreen) ---
          if (isBottomLayout) {
            return (
              <div>
                {toggleBar}
                <div className="flex gap-4">
                  <div className="w-1/2 min-w-0">{viewerBlock}</div>
                  <div className="w-1/2 min-w-0 overflow-y-auto" style={{ maxHeight: '70vh', scrollbarWidth: 'thin', scrollbarColor: '#e5e7eb transparent' }}>
                    {sideContent}
                  </div>
                </div>
              </div>
            )
          }

          return (
            <div className="space-y-4">
              {toggleBar}
              {viewerBlock}
              {effectiveMode === 'protein' && hasProtein && <div className="mt-3">{sideContent}</div>}
              {effectiveMode === 'ligand' && hasLigand && sideContent}
              {effectiveMode === 'results' && sideContent}
            </div>
          )
        })()}
        </PanelErrorBoundary>
      </div>

      {/* Prev / Next navigation */}
      <div className="flex items-center justify-between px-4 py-2.5 border-t border-gray-100 bg-gray-50/50 flex-shrink-0">
        <button
          onClick={() => prevMol && onRowClick && onRowClick(prevMol)}
          disabled={!prevMol}
          className={`flex items-center gap-1 text-sm font-medium transition-colors ${
            prevMol ? 'text-bx-light-text hover:text-[#1a2332]' : 'text-gray-300 cursor-not-allowed'
          }`}
        >
          <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 19l-7-7 7-7" />
          </svg>
          {prevMol ? prevMol.name || 'Prev' : 'Prev'}
        </button>
        <span className="text-[10px] text-gray-400 tabular-nums">
          {currentIdx + 1} / {molecules.length}
        </span>
        <button
          onClick={() => nextMol && onRowClick && onRowClick(nextMol)}
          disabled={!nextMol}
          className={`flex items-center gap-1 text-sm font-medium transition-colors ${
            nextMol ? 'text-bx-light-text hover:text-[#1a2332]' : 'text-gray-300 cursor-not-allowed'
          }`}
        >
          {nextMol ? nextMol.name || 'Next' : 'Next'}
          <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
          </svg>
        </button>
      </div>
    </div>
  )

  // In fullscreen: render with backdrop overlay
  if (panelFullscreen) {
    return (
      <>
        <div className="fixed inset-0 z-40 bg-black/30" onClick={() => setPanelFullscreen(false)} />
        {panelContent}
      </>
    )
  }

  return panelContent
}

export default MoleculeDetailPanel
export { DetailPopupModal }
