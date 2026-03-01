import React, { useEffect, useState, useMemo, useCallback, useRef } from 'react'
import { useParams, Link, useNavigate } from 'react-router-dom'

import { useWorkspace } from '../contexts/WorkspaceContext.jsx'
import { useToast } from '../contexts/ToastContext.jsx'
import { ALL_COLUMNS, COLUMN_PRESETS, PHASE_TYPES, flattenMoleculeProperties, detectAvailableColumns } from '../lib/columns.js'

import MoleculeTable from '../components/MoleculeTable.jsx'
import SelectionToolbar from '../components/SelectionToolbar.jsx'
import ColumnSelector from '../components/ColumnSelector.jsx'
import FilterBar from '../components/FilterBar.jsx'
import RunHistory from '../components/RunHistory.jsx'
import RunCreator from '../components/RunCreator.jsx'
import RunProgress from '../components/RunProgress.jsx'
import FreezeDialog from '../components/FreezeDialog.jsx'
import ParetoFront from '../components/ParetoFront.jsx'
import Badge from '../components/Badge.jsx'
import BindXLogo from '../components/BindXLogo.jsx'
import ProteinViewer from '../components/ProteinViewer.jsx'
import SafetyReport from '../components/SafetyReport.jsx'
import SynthesisTree from '../components/SynthesisTree.jsx'
import ConfidenceBreakdown from '../components/ConfidenceBreakdown.jsx'
import ExportModal from '../components/ExportModal.jsx'
import CampaignAgentPanel from '../components/CampaignAgentPanel.jsx'
import InfoTip, { TIPS } from '../components/InfoTip.jsx'

// ---------------------------------------------------------------------------
// Freeze Banner
// ---------------------------------------------------------------------------
function FreezeBanner() {
  return (
    <div className="flex items-center gap-2.5 rounded-[9px] px-4 py-2.5 text-sm" style={{ background: 'var(--blue-soft)', color: 'var(--blue)', border: '1px solid rgba(59,130,246,.12)' }}>
      <svg className="w-4 h-4 text-blue-500 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
          d="M12 15v2m-6 4h12a2 2 0 002-2v-6a2 2 0 00-2-2H6a2 2 0 00-2 2v6a2 2 0 002 2zm10-10V7a4 4 0 00-8 0v4h8z" />
      </svg>
      <div>
        <span className="font-semibold">Phase frozen</span>
        <span className="ml-1 font-normal">— molecules locked for downstream analysis. No new runs or edits allowed.</span>
      </div>
    </div>
  )
}

// ---------------------------------------------------------------------------
// Stats KPI cards
// ---------------------------------------------------------------------------
function StatsBar({ stats }) {
  const cards = [
    { value: stats.total_molecules ?? 0, label: 'Molecules', colorClass: 'text-bx-mint' },
    { value: stats.bookmarked ?? 0, label: 'Bookmarked', colorClass: 'text-bx-amber' },
    { value: stats.runs_completed ?? 0, label: 'Runs done', colorClass: 'text-bx-cyan' },
    {
      value: stats.runs_running ?? 0,
      label: 'Running',
      colorClass: stats.runs_running > 0 ? 'text-bx-blue' : 'text-bx-dim',
      animated: stats.runs_running > 0,
    },
  ]

  return (
    <div className="dark-inset">
      <div className="grid grid-cols-4 gap-4">
        {cards.map(card => (
          <div key={card.label} className="text-center py-2">
            <div className="flex items-center justify-center gap-1.5 mb-0.5">
              {card.animated && <BindXLogo variant="loading" size={14} />}
              <p className={`stat-value ${card.colorClass}`}>{card.value}</p>
            </div>
            <p className="stat-label">{card.label}</p>
          </div>
        ))}
      </div>
    </div>
  )
}

// ---------------------------------------------------------------------------
// ADMET Radar SVG (standalone, used in detail panel)
// ---------------------------------------------------------------------------
function AdmetRadar({ admet }) {
  const keys = ['absorption', 'distribution', 'metabolism', 'excretion', 'toxicity', 'pharmacodynamics']
  const labels = ['Absorption', 'Distribution', 'Metabolism', 'Excretion', 'Toxicity', 'Pharm.']
  const values = keys.map(k => (typeof admet[k] === 'number' ? admet[k] : 0.5))
  const N = values.length, CX = 70, CY = 70, R = 52

  function pt(v, i) {
    const a = (Math.PI * 2 * i) / N - Math.PI / 2
    return { x: CX + v * R * Math.cos(a), y: CY + v * R * Math.sin(a) }
  }

  const dataPoints = values.map((v, i) => pt(v, i))
  const polygon = dataPoints.map(p => `${p.x.toFixed(1)},${p.y.toFixed(1)}`).join(' ')

  return (
    <div className="flex items-center gap-4">
      <svg width={140} height={140} viewBox="0 0 140 140" className="flex-shrink-0">
        {/* Grid rings */}
        {[0.25, 0.5, 0.75, 1.0].map((r, ri) => {
          const pts = values.map((_, i) => {
            const p = pt(r, i)
            return `${p.x.toFixed(1)},${p.y.toFixed(1)}`
          }).join(' ')
          return (
            <polygon key={ri} points={pts} fill="none"
              stroke={r === 1.0 ? '#d1d5db' : '#e5e7eb'}
              strokeWidth={r === 1.0 ? 1 : 0.5} />
          )
        })}
        {/* Axes */}
        {values.map((_, i) => {
          const tip = pt(1, i)
          return <line key={i} x1={CX} y1={CY} x2={tip.x.toFixed(1)} y2={tip.y.toFixed(1)}
            stroke="#d1d5db" strokeWidth={0.75} />
        })}
        {/* Data polygon */}
        <polygon points={polygon} fill="#00e6a0" fillOpacity={0.15} stroke="#00e6a0" strokeWidth={1.5} />
        {/* Data points */}
        {dataPoints.map((p, i) => (
          <circle key={i} cx={p.x.toFixed(1)} cy={p.y.toFixed(1)} r={3}
            fill="#00e6a0" stroke="white" strokeWidth={1.5} />
        ))}
        {/* Labels */}
        {values.map((_, i) => {
          const lp = pt(1.3, i)
          const anchor = lp.x < CX - 5 ? 'end' : lp.x > CX + 5 ? 'start' : 'middle'
          return (
            <text key={i} x={lp.x.toFixed(1)} y={lp.y.toFixed(1)}
              fontSize={7} fill="#6b7280" textAnchor={anchor} dominantBaseline="central">
              {labels[i]}
            </text>
          )
        })}
      </svg>
      <div className="flex-1 space-y-1.5">
        {keys.map((k, i) => {
          const v = values[i]
          const pct = Math.round(v * 100)
          const color = v >= 0.7 ? 'bg-green-400' : v >= 0.5 ? 'bg-yellow-400' : 'bg-red-400'
          return (
            <div key={k}>
              <div className="flex justify-between text-[10px] mb-0.5">
                <span className="text-gray-500 capitalize">{k}</span>
                <span className="font-medium text-gray-700 tabular-nums">{pct}%</span>
              </div>
              <div className="w-full bg-gray-100 rounded-full h-1">
                <div className={`h-1 rounded-full ${color} transition-all`} style={{ width: `${pct}%` }} />
              </div>
            </div>
          )
        })}
      </div>
    </div>
  )
}

// ---------------------------------------------------------------------------
// Detail panel tab content helpers
// ---------------------------------------------------------------------------
function ScoresTab({ mol }) {
  const scores = [
    { label: 'Docking', value: mol.docking_score, unit: 'kcal/mol', color: 'bg-bx-surface', textColor: 'text-bx-light-text', format: v => v.toFixed(1) },
    { label: 'CNN Score', value: mol.cnn_score, unit: '', color: 'bg-green-500', textColor: 'text-green-600', format: v => v.toFixed(2) },
    { label: 'CNN Aff.', value: mol.cnn_affinity, unit: 'pKi', color: 'bg-teal-500', textColor: 'text-teal-600', format: v => v.toFixed(1) },
    { label: 'Composite', value: mol.composite_score, unit: '/100', color: 'bg-amber-500', textColor: 'text-amber-700', format: v => v.toFixed(1) },
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
                style={{ width: s.label === 'Composite' ? `${s.value}%` : s.label === 'CNN Score' ? `${s.value * 100}%` : '60%' }}
              />
            </div>
          </div>
        ))}
      </div>
      {scores.every(s => s.value == null) && (
        <p className="text-sm text-gray-400 text-center py-4">No scores available yet — run Docking and Scoring analyses.</p>
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
    </div>
  )
}

function AdmetTab({ mol }) {
  // Try mol.admet first (nested), then reconstruct from flat properties
  const admet = mol.admet || (() => {
    const keys = ['oral_bioavailability', 'solubility', 'BBB', 'metabolic_stability', 'hERG', 'QED']
    const found = keys.filter(k => mol[k] != null)
    if (!found.length) return null
    const obj = {}
    found.forEach(k => { obj[k] = mol[k] })
    return obj
  })()
  if (!admet) {
    return <p className="text-sm text-gray-400 text-center py-4">No ADMET data — run the ADMET analysis first.</p>
  }
  return <AdmetRadar admet={admet} />
}

function SafetyTab({ mol, details }) {
  // Try details.safety first, then reconstruct from flat molecule props
  const safety = details?.safety || (() => {
    const hasAny = mol.herg_risk != null || mol.ames_mutagenicity != null || mol.hepatotoxicity != null ||
                   mol.pains_alert != null || mol.safety_color_code != null
    if (!hasAny) return null
    return {
      herg_risk: mol.herg_risk ?? mol.hERG ?? null,
      ames_risk: mol.ames_mutagenicity ?? null,
      hepatotox_risk: mol.hepatotoxicity ?? null,
      pains_pass: mol.pains_alert != null ? !mol.pains_alert : null,
      pains_alerts: mol.pains_alert ? [mol.pains_alert] : [],
    }
  })()
  if (!safety) {
    return <p className="text-sm text-gray-400 text-center py-4">No safety data available for this molecule.</p>
  }

  const riskColor = r => r === 'high' ? 'text-red-600 bg-red-50 border-red-200' :
                         r === 'medium' ? 'text-amber-600 bg-amber-50 border-amber-200' :
                         'text-green-600 bg-green-50 border-green-200'

  return (
    <div className="space-y-3">
      {/* PAINS */}
      <div className={`flex items-center justify-between px-3 py-2 rounded-lg text-sm border ${
        safety.pains_pass ? 'bg-green-50 border-green-200 text-green-700' : 'bg-red-50 border-red-200 text-red-600'
      }`}>
        <span className="font-semibold">PAINS Filter</span>
        <span>{safety.pains_pass ? 'Pass' : `Fail (${safety.pains_alerts?.length || 0} alerts)`}</span>
      </div>

      {/* Brenk alerts */}
      {safety.brenk_alerts?.length > 0 && (
        <div className="bg-amber-50 border border-amber-200 rounded-lg px-3 py-2">
          <p className="text-[10px] font-semibold text-amber-700 uppercase tracking-wider mb-1">Brenk Alerts</p>
          {safety.brenk_alerts.map((a, i) => (
            <p key={i} className="text-sm text-amber-700">{a}</p>
          ))}
        </div>
      )}

      {/* Off-target */}
      {safety.off_target?.length > 0 && (
        <div>
          <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wider mb-1.5">Off-Target Risks</p>
          <div className="space-y-1.5">
            {safety.off_target.map((ot, i) => (
              <div key={i} className="flex items-center justify-between bg-gray-50 rounded-lg px-2.5 py-1.5">
                <div>
                  <p className="text-sm font-medium text-gray-700">{ot.target}</p>
                  <p className="text-[10px] text-gray-400">{ot.family} · similarity {(ot.similarity * 100).toFixed(0)}%</p>
                </div>
                <span className={`text-[10px] font-semibold px-2 py-0.5 rounded-full border capitalize ${riskColor(ot.risk)}`}>
                  {ot.risk}
                </span>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* hERG */}
      {safety.herg_risk != null && (
        <div>
          <div className="flex justify-between text-[10px] mb-1">
            <span className="text-gray-500 font-semibold">hERG Risk</span>
            <span className={`font-bold ${safety.herg_risk > 0.3 ? 'text-red-600' : 'text-green-600'}`}>
              {(safety.herg_risk * 100).toFixed(0)}%
            </span>
          </div>
          <div className="w-full bg-gray-100 rounded-full h-1.5">
            <div
              className={`h-1.5 rounded-full ${safety.herg_risk > 0.3 ? 'bg-red-400' : safety.herg_risk > 0.15 ? 'bg-amber-400' : 'bg-green-400'}`}
              style={{ width: `${safety.herg_risk * 100}%` }}
            />
          </div>
        </div>
      )}

      {/* Ames + Hepatotox */}
      <div className="grid grid-cols-2 gap-2 text-sm">
        {safety.ames_risk != null && (
          <div className="bg-gray-50 rounded-lg px-2.5 py-2 border border-gray-100">
            <p className="text-[9px] text-gray-400 uppercase tracking-wider">Ames risk</p>
            <p className={`font-bold tabular-nums ${safety.ames_risk > 0.3 ? 'text-red-600' : 'text-green-600'}`}>
              {(safety.ames_risk * 100).toFixed(0)}%
            </p>
          </div>
        )}
        {safety.hepatotox_risk != null && (
          <div className="bg-gray-50 rounded-lg px-2.5 py-2 border border-gray-100">
            <p className="text-[9px] text-gray-400 uppercase tracking-wider">Hepatotox.</p>
            <p className={`font-bold tabular-nums ${safety.hepatotox_risk > 0.3 ? 'text-red-600' : 'text-green-600'}`}>
              {(safety.hepatotox_risk * 100).toFixed(0)}%
            </p>
          </div>
        )}
      </div>
    </div>
  )
}

function SynthesisTab({ details }) {
  const synth = details?.synthesis
  if (!synth) {
    return <p className="text-sm text-gray-400 text-center py-4">No retrosynthesis data for this molecule.</p>
  }

  return (
    <div className="space-y-3">
      {/* Summary */}
      <div className="grid grid-cols-3 gap-2 text-sm">
        <div className="bg-gray-50 rounded-lg px-2.5 py-2 text-center border border-gray-100">
          <p className="text-[9px] text-gray-400 uppercase">Steps</p>
          <p className="font-bold text-gray-800 text-base">{synth.num_steps}</p>
        </div>
        <div className="bg-gray-50 rounded-lg px-2.5 py-2 text-center border border-gray-100">
          <p className="text-[9px] text-gray-400 uppercase">Cost</p>
          <p className="font-bold text-gray-800 text-base">{synth.total_cost}</p>
        </div>
        <div className={`rounded-lg px-2.5 py-2 text-center border ${synth.feasibility >= 0.7 ? 'bg-green-50 border-green-200' : synth.feasibility >= 0.5 ? 'bg-amber-50 border-amber-200' : 'bg-red-50 border-red-200'}`}>
          <p className="text-[9px] text-gray-400 uppercase">Feasibility</p>
          <p className={`font-bold text-base ${synth.feasibility >= 0.7 ? 'text-green-600' : synth.feasibility >= 0.5 ? 'text-amber-600' : 'text-red-600'}`}>
            {(synth.feasibility * 100).toFixed(0)}%
          </p>
        </div>
      </div>

      {/* Steps */}
      <div className="space-y-2">
        <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wider">Retrosynthesis route</p>
        {synth.steps.map((step, i) => (
          <div key={i} className="flex items-start gap-2.5">
            <div className="flex-shrink-0 w-5 h-5 rounded-full bg-bx-surface text-white text-[9px] font-bold flex items-center justify-center mt-0.5">
              {i + 1}
            </div>
            <div className="flex-1 bg-gray-50 rounded-lg p-2.5 border border-gray-100">
              <div className="flex items-start justify-between gap-2">
                <div>
                  <p className="text-sm font-semibold text-gray-700">{step.product}</p>
                  <p className="text-[10px] text-gray-500 mt-0.5">{step.reagent}</p>
                  <p className="text-[9px] text-gray-400 italic mt-0.5">{step.conditions}</p>
                </div>
                <div className="text-right flex-shrink-0">
                  <p className="text-sm font-bold text-bx-mint tabular-nums">{(step.yield * 100).toFixed(0)}%</p>
                  <p className="text-[9px] text-gray-400">yield</p>
                </div>
              </div>
              <div className="mt-2 w-full bg-gray-200 rounded-full h-1">
                <div className="h-1 bg-bx-mint rounded-full" style={{ width: `${step.yield * 100}%` }} />
              </div>
            </div>
          </div>
        ))}
      </div>
    </div>
  )
}

function InteractionsTab({ details }) {
  const inter = details?.interactions
  if (!inter) {
    return <p className="text-sm text-gray-400 text-center py-4">No interaction data — run Enrichment to compute ProLIF interactions.</p>
  }

  const typeColors = {
    'HBond': 'bg-blue-100 text-blue-700',
    'Hydrophobic': 'bg-orange-100 text-orange-700',
    'Salt Bridge': 'bg-purple-100 text-purple-700',
    'Pi-Stacking': 'bg-indigo-100 text-indigo-700',
    'Covalent': 'bg-red-100 text-red-700',
    'Water Bridge': 'bg-cyan-100 text-cyan-700',
  }

  const maxDist = inter.residues ? Math.max(...inter.residues.map(r => r.distance)) : 5

  return (
    <div className="space-y-3">
      {/* Summary counts */}
      <div className="flex flex-wrap gap-2">
        {[
          { label: 'Total', value: inter.total_contacts, color: 'text-gray-700' },
          { label: 'H-bonds', value: inter.hbonds, color: 'text-blue-600' },
          { label: 'Hydrophobic', value: inter.hydrophobic, color: 'text-orange-600' },
          { label: 'Ionic', value: inter.ionic, color: 'text-purple-600' },
          { label: 'Pi-stack', value: inter.pi_stacking, color: 'text-indigo-600' },
          inter.covalent ? { label: 'Covalent', value: inter.covalent, color: 'text-red-600' } : null,
        ].filter(Boolean).map(s => (
          <div key={s.label} className="bg-gray-50 rounded-lg px-2.5 py-1.5 text-center border border-gray-100 min-w-[52px]">
            <p className={`text-base font-bold tabular-nums ${s.color}`}>{s.value}</p>
            <p className="text-[9px] text-gray-400">{s.label}</p>
          </div>
        ))}
      </div>

      {/* Residue list */}
      <div className="space-y-1.5">
        <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wider">Interacting residues</p>
        {inter.residues?.map((r, i) => (
          <div key={i} className="flex items-center gap-2.5 py-1.5 border-b border-gray-50 last:border-0">
            <span className="text-sm font-mono font-semibold text-gray-700 w-16 flex-shrink-0">{r.name}</span>
            <span className={`text-[9px] font-semibold px-1.5 py-0.5 rounded-full ${typeColors[r.type] || 'bg-gray-100 text-gray-600'} flex-shrink-0`}>
              {r.type}
            </span>
            <div className="flex-1 flex items-center gap-1.5">
              <div className="flex-1 bg-gray-100 rounded-full h-1">
                <div
                  className="h-1 bg-bx-surface rounded-full"
                  style={{ width: `${((maxDist - r.distance) / maxDist) * 100}%`, opacity: 0.6 }}
                />
              </div>
              <span className="text-[9px] text-gray-400 tabular-nums flex-shrink-0">{r.distance.toFixed(1)} Å</span>
            </div>
            <span className="text-[9px] text-gray-400 flex-shrink-0">{r.chain}</span>
          </div>
        ))}
      </div>
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
      case 'confidence':
        return (
          <ConfidenceBreakdown
            confidence={props.confidence || (flat.confidence_score != null ? {
              overall: flat.confidence_score,
              components: {
                structure: { score: flat.structure_confidence, label: 'Structure' },
                docking: { score: flat.docking_confidence, label: 'Docking' },
                admet: { score: flat.admet_confidence, label: 'ADMET' },
              },
              flags: flat.confidence_flags || [],
            } : null)}
            moleculeName={molecule.name || molecule.id}
            pipeline_summary={props.pipeline_summary || null}
          />
        )
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
  { id: 'scores',       label: 'Scores',     tip: 'Docking and composite scores from computational screening' },
  { id: 'properties',   label: 'Props',      tip: 'Physicochemical properties: MW, logP, TPSA, Lipinski Ro5' },
  { id: 'admet',        label: 'ADMET',      tip: 'Absorption, Distribution, Metabolism, Excretion, Toxicity predictions' },
  { id: 'safety',       label: 'Safety',     tip: 'Safety profile: hERG, AMES mutagenicity, hepatotoxicity, PAINS filters' },
  { id: 'synthesis',    label: 'Synth.',     tip: 'Retrosynthesis feasibility, estimated cost, and synthetic route' },
  { id: 'interactions', label: 'Interact.',  tip: 'Protein-ligand interaction fingerprints (ProLIF): H-bonds, hydrophobic, ionic contacts' },
]

function MoleculeDetailPanel({ molecule, molecules, onClose, onToggleBookmark, onRowClick, isFrozen, project }) {
  const [activeTab, setActiveTab] = useState('scores')
  const [copied, setCopied] = useState(false)

  // Extract target structure info from project
  const targetPreview = project?.target_preview || {}
  const pdbUrl = targetPreview.structure?.download_url || null
  const selectedPocket = useMemo(() => {
    const pockets = targetPreview.pockets || []
    const idx = targetPreview.selected_pocket_index
    return idx != null && pockets[idx] ? pockets[idx] : null
  }, [targetPreview])

  const currentIdx = molecules.findIndex(m => m.id === molecule.id)
  const prevMol = currentIdx > 0 ? molecules[currentIdx - 1] : null
  const nextMol = currentIdx < molecules.length - 1 ? molecules[currentIdx + 1] : null

  // Build details from molecule properties (API data, not mock)
  const details = useMemo(() => {
    const props = molecule.properties || {}
    return {
      interactions: props.enrichment || null,
      synthesis: props.retrosynthesis ? {
        feasibility: props.retrosynthesis.synth_confidence,
        total_cost: props.retrosynthesis.synth_cost_estimate,
        num_steps: props.retrosynthesis.n_synth_steps,
        steps: props.retrosynthesis.steps || [],
      } : null,
    }
  }, [molecule])

  function copySmiles() {
    if (!molecule.smiles) return
    navigator.clipboard.writeText(molecule.smiles).then(() => {
      setCopied(true)
      setTimeout(() => setCopied(false), 1500)
    })
  }

  return (
    <div className="bg-white rounded-xl border border-gray-200 shadow-xl overflow-hidden flex flex-col sticky top-4"
      style={{ maxHeight: 'calc(100vh - 120px)' }}>

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
          {molecule.bookmarked && (
            <Badge variant="yellow" size="sm">Bookmarked</Badge>
          )}
        </div>
        <button
          onClick={onClose}
          className="flex-shrink-0 ml-2 p-1.5 rounded-lg hover:bg-gray-100 text-gray-400 hover:text-gray-700 transition-colors"
          aria-label="Close detail panel"
        >
          <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
          </svg>
        </button>
      </div>

      {/* Scrollable content */}
      <div className="overflow-y-auto flex-1 p-4 space-y-4"
        style={{ scrollbarWidth: 'thin', scrollbarColor: '#e5e7eb transparent' }}>

        {/* 3D Viewer */}
        {pdbUrl ? (
          <ProteinViewer
            pdbUrl={pdbUrl}
            selectedPocket={selectedPocket}
            height={280}
          />
        ) : (
          <div className="rounded-xl bg-[#0a0d15] border border-bx-surface/30 overflow-hidden">
            <div className="px-3 py-2 border-b border-white/5 flex items-center justify-between">
              <div className="flex items-center gap-2">
                <div className="w-2 h-2 rounded-full bg-red-500" />
                <div className="w-2 h-2 rounded-full bg-amber-400" />
                <div className="w-2 h-2 rounded-full bg-green-400" />
                <span className="ml-1 text-[10px] font-semibold text-white/50 uppercase tracking-wide">3D Structure</span>
              </div>
            </div>
            <div className="h-28 flex flex-col items-center justify-center gap-2">
              <svg className="w-10 h-10 text-bx-light-text/60" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={0.75}
                  d="M14 10l-2 1m0 0l-2-1m2 1v2.5M20 7l-2 1m2-1l-2-1m2 1v2.5M14 4l-2-1-2 1M4 7l2-1M4 7l2 1M4 7v2.5M12 21l-2-1m2 1l2-1m-2 1v-2.5M6 18l-2-1v-2.5M18 18l2-1v-2.5" />
              </svg>
              <p className="text-[10px] text-white/25">Configure target in Target Setup to view structure</p>
            </div>
          </div>
        )}

        {/* Tabs */}
        <div>
          <div className="flex border-b border-gray-100 -mx-4 px-4 gap-0.5 overflow-x-auto"
            style={{ scrollbarWidth: 'none' }}>
            {TABS.map(tab => (
              <button
                key={tab.id}
                onClick={() => setActiveTab(tab.id)}
                className={`flex-shrink-0 px-2.5 py-2 text-sm font-semibold border-b-2 transition-all duration-150 whitespace-nowrap inline-flex items-center ${
                  activeTab === tab.id
                    ? 'border-bx-surface text-bx-light-text'
                    : 'border-transparent text-gray-400 hover:text-gray-600 hover:border-gray-200'
                }`}
              >
                {tab.label}
                {tab.tip && <InfoTip text={tab.tip} size="xs" />}
              </button>
            ))}
          </div>

          <div className="pt-3">
            {activeTab === 'scores' && <ScoresTab mol={molecule} />}
            {activeTab === 'properties' && <PropertiesTab mol={molecule} />}
            {activeTab === 'admet' && <AdmetTab mol={molecule} />}
            {activeTab === 'safety' && <SafetyTab mol={molecule} details={details} />}
            {activeTab === 'synthesis' && <SynthesisTab details={details} />}
            {activeTab === 'interactions' && <InteractionsTab details={details} />}
          </div>
        </div>

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
}

// ---------------------------------------------------------------------------
// PhaseDashboard — main view
// ---------------------------------------------------------------------------
export default function PhaseDashboard() {
  const { projectId, phaseId } = useParams()

  const {
    currentProject,
    currentCampaign,
    currentPhase,
    phaseMolecules,
    currentPhaseRuns,
    runsLoading,
    moleculesLoading,
    selectProject,
    selectPhase,
    selectedMoleculeIds,
    toggleSelection,
    selectAll,
    selectNone,
    selectBookmarked,
    toggleBookmark,
    bookmarkSelected,
    freezePhase,
    unfreezePhase,
    getPhaseStatus,
    createRun,
    cancelRun,
    archiveRun,
    importFile,
    importDatabase,
  } = useWorkspace()

  const { addToast } = useToast()

  // Sync URL params to workspace state
  useEffect(() => { if (projectId) selectProject(projectId) }, [projectId, selectProject])
  useEffect(() => { if (phaseId) selectPhase(phaseId) }, [phaseId, selectPhase])

  // --- Flatten molecule properties from API (EAV → flat columns) ---
  const flatMolecules = useMemo(
    () => phaseMolecules.map(flattenMoleculeProperties),
    [phaseMolecules]
  )

  // --- Detect which columns have data ---
  const availableColumns = useMemo(
    () => detectAvailableColumns(flatMolecules),
    [flatMolecules]
  )

  // --- User annotations (tags, invalidation, comments) — local state keyed by mol id ---
  const [annotations, setAnnotations] = useState({}) // { molId: { tags: [...], invalidated: bool, user_comment: '...' } }

  const handleAnnotation = useCallback((molId, key, value) => {
    setAnnotations(prev => ({
      ...prev,
      [molId]: { ...(prev[molId] || {}), [key]: value },
    }))
  }, [])

  // Merge annotations into flat molecules
  const annotatedMolecules = useMemo(
    () => flatMolecules.map(m => ({ ...m, ...(annotations[m.id] || {}) })),
    [flatMolecules, annotations]
  )

  // --- Local UI state ---
  const [selectedMolecule, setSelectedMolecule] = useState(null)
  const [showRunCreator, setShowRunCreator] = useState(false)
  const [showPareto, setShowPareto] = useState(true)
  const [showFreezeDialog, setShowFreezeDialog] = useState(false)
  const [filteredMolecules, setFilteredMolecules] = useState(annotatedMolecules)
  const [popupState, setPopupState] = useState(null) // { type: 'safety'|'retrosynthesis'|'confidence', molecule }
  const [showExportModal, setShowExportModal] = useState(false)
  const [showAgentPanel, setShowAgentPanel] = useState(false)
  const [showRunLogs, setShowRunLogs] = useState(false)

  // Update filtered whenever phase molecules change
  useEffect(() => { setFilteredMolecules(annotatedMolecules) }, [annotatedMolecules])

  // Auto-select first molecule when molecules load and none is selected
  useEffect(() => {
    if (flatMolecules.length > 0 && !selectedMolecule) {
      setSelectedMolecule(flatMolecules[0])
    }
  }, [flatMolecules.length])

  // Active run (first created/running run) + queued count
  const activeRun = useMemo(
    () => currentPhaseRuns.find(r => r.status === 'created' || r.status === 'running') || null,
    [currentPhaseRuns]
  )
  const queuedRunCount = useMemo(
    () => Math.max(0, currentPhaseRuns.filter(r => r.status === 'created' || r.status === 'running').length - 1),
    [currentPhaseRuns]
  )

  // Track run completions/failures for toast notifications
  const prevRunsRef = useRef(currentPhaseRuns)
  useEffect(() => {
    const prev = prevRunsRef.current
    prevRunsRef.current = currentPhaseRuns
    if (!prev.length) return

    for (const run of currentPhaseRuns) {
      const prevRun = prev.find(r => r.id === run.id)
      if (!prevRun) continue
      if (prevRun.status !== 'completed' && run.status === 'completed') {
        addToast('Run completed successfully', 'success')
      }
      if (prevRun.status !== 'failed' && run.status === 'failed') {
        addToast(run.error_message || 'Run failed', 'error')
      }
    }
  }, [currentPhaseRuns, addToast])

  // Column visibility — show all available columns by default (CDC §4.1)
  const [visibleKeys, setVisibleKeys] = useState(() => COLUMN_PRESETS.hit_discovery)
  const availableKeysStr = useMemo(
    () => availableColumns.map(c => c.key).join(','),
    [availableColumns]
  )
  useEffect(() => {
    if (currentPhase) {
      // CDC: "all columns with results are shown by default"
      // Use available columns if we have data, otherwise fall back to phase preset
      if (availableColumns.length > 1) {
        setVisibleKeys(availableColumns.map(c => c.key))
      } else {
        const preset = COLUMN_PRESETS[currentPhase.type] || COLUMN_PRESETS.hit_discovery
        setVisibleKeys(preset)
      }
      setSelectedMolecule(null)
    }
  }, [currentPhase?.id, availableKeysStr])

  // Freeze status
  const freezeOverride = currentPhase ? getPhaseStatus(currentPhase.id) : null
  const effectiveStatus = freezeOverride || currentPhase?.status || 'active'
  const isFrozen = effectiveStatus === 'frozen'

  // Phase type meta
  const phaseTypeMeta = currentPhase ? (PHASE_TYPES[currentPhase.type] || {}) : {}

  // Visible columns
  const visibleColumns = useMemo(
    () => ALL_COLUMNS.filter(c => visibleKeys.includes(c.key)),
    [visibleKeys]
  )

  // Bookmarked count (live)
  const bookmarkedCount = useMemo(
    () => flatMolecules.filter(m => m.bookmarked).length,
    [flatMolecules]
  )

  // Live stats
  const liveStats = useMemo(() => {
    if (!currentPhase) return { total_molecules: 0, bookmarked: 0, runs_completed: 0, runs_running: 0 }
    const completedRuns = currentPhaseRuns.filter(r => r.status === 'completed').length
    const runningRuns = currentPhaseRuns.filter(r => r.status === 'running' || r.status === 'created').length
    return {
      total_molecules: flatMolecules.length,
      bookmarked: bookmarkedCount,
      runs_completed: completedRuns,
      runs_running: runningRuns,
    }
  }, [currentPhase, currentPhaseRuns, flatMolecules, bookmarkedCount])

  // Number of active filters
  const [activeFilterCount, setActiveFilterCount] = useState(0)

  const handleFilteredChange = useCallback((result) => {
    setFilteredMolecules(result)
  }, [])

  // Row click — opens detail panel
  const handleRowClick = useCallback((mol) => {
    setSelectedMolecule(prev => prev?.id === mol.id ? null : mol)
  }, [])

  const handleCloseDetail = useCallback(() => setSelectedMolecule(null), [])

  const handleCellPopup = useCallback((type, mol) => {
    // Find the original molecule with full properties (not flattened)
    const original = phaseMolecules.find(m => m.id === mol.id)
    setPopupState({ type, molecule: { ...mol, properties: original?.properties || {} } })
  }, [phaseMolecules])

  const [runSubmitting, setRunSubmitting] = useState(false)

  const handleNewRun = useCallback(() => {
    if (activeRun) {
      addToast('A run is already in progress. Wait for it to finish before launching a new one.', 'warning')
      return
    }
    setShowRunCreator(true)
  }, [activeRun, addToast])

  const handleRunSubmit = useCallback(async (runConfig) => {
    if (!phaseId) return
    try {
      setRunSubmitting(true)

      if (runConfig.type === 'import' && runConfig.config?.sourceMode === 'external' && runConfig.config?._file) {
        // File upload import
        await importFile(phaseId, runConfig.config._file)
        addToast('File uploaded and import started', 'success')
      } else if (runConfig.type === 'import' && runConfig.config?.sourceMode === 'internal') {
        // Internal import — bookmarked molecules from another phase
        await createRun(phaseId, {
          type: 'import',
          config: {
            source: 'phase_selection',
            selection: 'bookmarked',
            source_phase_id: runConfig.config.source_phase_id,
          },
        })
        addToast('Importing bookmarked molecules...', 'success')
      } else if (runConfig.type === 'import' && runConfig.config?.sourceMode === 'database') {
        // Database import — get uniprot_id from project target
        const uniprotId = currentProject?.target_input_value || currentProject?.target_preview?.uniprot?.id || ''
        await importDatabase(phaseId, {
          databases: runConfig.config.databases,
          uniprot_id: uniprotId,
          max_per_source: runConfig.config.maxPerSource || 50,
          filters: runConfig.config.filters || {},
        })
        addToast(`Fetching compounds from ${runConfig.config.databases.join(', ')}...`, 'success')
      } else {
        // Standard run (calculation, generation, SMILES import)
        await createRun(phaseId, runConfig)
        addToast('Run created and dispatched', 'success')
      }

      setShowRunCreator(false)
    } catch (err) {
      addToast(err.userMessage || err.message || 'Failed to create run', 'error')
    } finally {
      setRunSubmitting(false)
    }
  }, [phaseId, createRun, importFile, importDatabase, currentProject, addToast])

  const handleCancelRun = useCallback(async (runId) => {
    try {
      await cancelRun(runId)
      addToast('Run cancelled', 'info')
    } catch (err) {
      addToast(err.userMessage || 'Failed to cancel run', 'error')
    }
  }, [cancelRun, addToast])

  const handleArchiveRun = useCallback(async (runId) => {
    try {
      await archiveRun(runId)
      addToast('Run archived', 'info')
    } catch (err) {
      addToast(err.userMessage || 'Failed to archive run', 'error')
    }
  }, [archiveRun, addToast])

  const handleExport = useCallback(() => {
    setShowExportModal(true)
  }, [])

  const handleFreezeToggle = useCallback(() => setShowFreezeDialog(true), [])
  const handleFreezeConfirm = useCallback(() => {
    if (!currentPhase) return
    if (isFrozen) {
      unfreezePhase(currentPhase.id)
      addToast('Phase unfrozen — you can now modify molecules and run analyses', 'info')
    } else {
      const bCount = flatMolecules.filter(m => m.bookmarked).length
      freezePhase(currentPhase.id)
      addToast(`Phase frozen — ${bCount} bookmarked molecule${bCount !== 1 ? 's' : ''} locked for next phase`, 'success')
    }
  }, [isFrozen, currentPhase, freezePhase, unfreezePhase, addToast, flatMolecules])

  const handleSelectAll = useCallback(() => {
    if (flatMolecules.every(m => selectedMoleculeIds.has(m.id))) selectNone()
    else selectAll()
  }, [flatMolecules, selectedMoleculeIds, selectAll, selectNone])

  const handleSelectFiltered = useCallback(() => {
    filteredMolecules.forEach(m => {
      if (!selectedMoleculeIds.has(m.id)) toggleSelection(m.id)
    })
  }, [filteredMolecules, selectedMoleculeIds, toggleSelection])

  const hasDownstreamRuns = useMemo(() => {
    if (!currentCampaign || !currentPhase) return false
    const phases = currentCampaign.phases || []
    const myIdx = phases.findIndex(p => p.id === currentPhase.id)
    if (myIdx < 0) return false
    return phases.slice(myIdx + 1).some(p => (p.runs || []).length > 0)
  }, [currentCampaign, currentPhase])

  // Send bookmarks to next phase
  const navigate = useNavigate()
  const handleSendToNextPhase = useCallback(async () => {
    if (!currentCampaign || !currentPhase) return
    const phases = currentCampaign.phases || []
    const myIdx = phases.findIndex(p => p.id === currentPhase.id)
    if (myIdx < 0) return

    const nextPhase = phases[myIdx + 1]
    if (!nextPhase) {
      addToast('No next phase in this campaign. Create a new phase from the project page first.', 'warning')
      return
    }

    try {
      await createRun(nextPhase.id, {
        type: 'import',
        config: {
          source: 'phase_selection',
          selection: 'bookmarked',
          source_phase_id: currentPhase.id,
        },
      })
      addToast(`Sending ${bookmarkedCount} bookmarked molecules to ${nextPhase.label || 'next phase'}...`, 'success')
      navigate(`/project/${projectId}/phase/${nextPhase.id}`)
    } catch (err) {
      addToast(err.userMessage || 'Failed to send bookmarks to next phase', 'error')
    }
  }, [currentCampaign, currentPhase, createRun, bookmarkedCount, addToast, navigate, projectId])

  const showDetail = !!selectedMolecule

  // --- Loading / not found ---
  if (!currentPhase) {
    return (
      <div className="flex flex-col items-center justify-center h-64 gap-3">
        <div className="w-12 h-12 rounded-full bg-gray-100 flex items-center justify-center">
          <svg className="w-6 h-6 text-gray-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
              d="M9.172 16.172a4 4 0 015.656 0M9 10h.01M15 10h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
          </svg>
        </div>
        <p className="text-gray-400 text-sm">
          {currentProject ? 'Phase not found.' : 'Loading workspace...'}
        </p>
      </div>
    )
  }

  // --- Empty phase (no molecules) ---
  if (flatMolecules.length === 0) {
    return (
      <div className="space-y-4 pb-8">
        {/* Breadcrumb */}
        <Breadcrumb projectId={projectId} projectName={currentProject?.name} phase={currentPhase} phaseTypeMeta={phaseTypeMeta} />

        {/* Phase header */}
        <PhaseHeader
          phase={currentPhase}
          phaseTypeMeta={phaseTypeMeta}
          isFrozen={isFrozen}
          stats={liveStats}
          onFreezeToggle={handleFreezeToggle}
          onNewRun={handleNewRun}
          hasActiveRun={!!activeRun}
          onOpenAgent={() => setShowAgentPanel(true)}
        />

        {/* Active run progress */}
        <RunProgress run={activeRun} onCancel={handleCancelRun} queuedCount={queuedRunCount} />

        {/* Empty state */}
        <div className="card p-12 text-center">
          <div className="flex flex-col items-center gap-4 max-w-sm mx-auto">
            <div className="w-16 h-16 rounded-full bg-gray-100 flex items-center justify-center">
              <svg className="w-8 h-8 text-gray-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1}
                  d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
              </svg>
            </div>
            <div>
              <h3 className="text-base font-bold text-gray-700 mb-1">No molecules yet</h3>
              <p className="text-sm text-gray-500 leading-relaxed">
                Run an Import to add molecules to this phase. You can import from ChEMBL, ZINC, a file, or promote hits from a previous phase.
              </p>
            </div>
            {!isFrozen && (
              <button
                onClick={handleNewRun}
                className="flex items-center gap-2 px-5 py-2.5 bg-bx-surface text-white rounded-xl text-sm font-semibold hover:bg-bx-elevated transition-colors shadow-sm"
              >
                <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 4v16m8-8H4" />
                </svg>
                Import Molecules
              </button>
            )}
          </div>
        </div>

        <RunHistory runs={currentPhaseRuns} onCancel={handleCancelRun} onArchive={handleArchiveRun} />

        <RunCreator
          phaseId={phaseId}
          phaseType={currentPhase.type}
          isOpen={showRunCreator}
          onClose={() => setShowRunCreator(false)}
          onSubmit={handleRunSubmit}
          selectedMoleculeIds={selectedMoleculeIds}
          submitting={runSubmitting}
        />
        <FreezeDialog
          isOpen={showFreezeDialog}
          onClose={() => setShowFreezeDialog(false)}
          onConfirm={handleFreezeConfirm}
          action={isFrozen ? 'unfreeze' : 'freeze'}
          phaseName={currentPhase.label}
          bookmarkedMols={flatMolecules.filter(m => m.bookmarked)}
          hasDownstreamRuns={hasDownstreamRuns}
        />
      </div>
    )
  }

  const phase = currentPhase

  return (
    <div className="space-y-4 pb-8">
      {/* Breadcrumb */}
      <Breadcrumb projectId={projectId} projectName={currentProject?.name} phase={phase} phaseTypeMeta={phaseTypeMeta} />

      {/* Phase header card */}
      <PhaseHeader
        phase={phase}
        phaseTypeMeta={phaseTypeMeta}
        isFrozen={isFrozen}
        stats={liveStats}
        onFreezeToggle={handleFreezeToggle}
        onNewRun={handleNewRun}
        hasActiveRun={!!activeRun}
        onOpenAgent={() => setShowAgentPanel(true)}
      />

      {/* Active run progress */}
      <RunProgress run={activeRun} onCancel={handleCancelRun} queuedCount={queuedRunCount} />

      {/* Selection Toolbar */}
      <SelectionToolbar
        selectedCount={selectedMoleculeIds.size}
        totalCount={flatMolecules.length}
        bookmarkedCount={bookmarkedCount}
        filteredCount={filteredMolecules.length}
        activeFilterCount={activeFilterCount}
        onSelectAll={handleSelectAll}
        onSelectNone={selectNone}
        onSelectBookmarked={selectBookmarked}
        onSelectFiltered={handleSelectFiltered}
        onExport={handleExport}
        onBookmarkSelected={bookmarkSelected}
        onSendToNextPhase={!isFrozen && bookmarkedCount > 0 ? handleSendToNextPhase : undefined}
        isFrozen={isFrozen}
      />

      {/* Filter bar */}
      <FilterBarWithCount
        molecules={flatMolecules}
        columns={availableColumns}
        onFilteredChange={handleFilteredChange}
        onFilterCountChange={setActiveFilterCount}
      />

      {/* Table + Detail panel split layout */}
      <div className="flex gap-4 items-start">
        {/* Table (60% when drawer open, 100% when closed) */}
        <div className={`min-w-0 space-y-2 transition-all duration-200 ${showDetail ? 'flex-1' : 'w-full'}`}>
          {/* Column controls row */}
          <div className="flex items-center justify-between gap-3">
            <div className="flex items-center gap-2">
              <ColumnSelector
                allColumns={ALL_COLUMNS}
                visibleKeys={visibleKeys}
                onChange={setVisibleKeys}
              />
              <span className="text-sm text-gray-400 tabular-nums flex items-center gap-1.5">
                {moleculesLoading && <BindXLogo variant="loading" size={20} />}
                {filteredMolecules.length !== flatMolecules.length
                  ? `${filteredMolecules.length} of ${flatMolecules.length} molecules`
                  : `${flatMolecules.length} molecules`}
              </span>
            </div>
          </div>

          <MoleculeTable
            molecules={filteredMolecules}
            columns={visibleColumns}
            selectedIds={selectedMoleculeIds}
            onToggleSelect={toggleSelection}
            onSelectAll={handleSelectAll}
            onRowClick={handleRowClick}
            onToggleBookmark={isFrozen ? undefined : toggleBookmark}
            onCellPopup={handleCellPopup}
            onAnnotation={isFrozen ? undefined : handleAnnotation}
            activeRowId={selectedMolecule?.id || null}
          />
        </div>

        {/* Detail panel (40% fixed width when open) */}
        {showDetail && selectedMolecule && (
          <div className="flex-shrink-0 w-[40%] min-w-[300px] max-w-[520px]">
            <MoleculeDetailPanel
              molecule={flatMolecules.find(m => m.id === selectedMolecule.id) || selectedMolecule}
              molecules={filteredMolecules}
              onClose={handleCloseDetail}
              onToggleBookmark={isFrozen ? undefined : toggleBookmark}
              onRowClick={handleRowClick}
              isFrozen={isFrozen}
              project={currentProject}
            />
          </div>
        )}
      </div>

      {/* Run History timeline */}
      <RunHistory runs={currentPhaseRuns} onCancel={handleCancelRun} onArchive={handleArchiveRun} />

      {/* Pareto Analysis (collapsible) */}
      <div className="card overflow-hidden">
        <button
          onClick={() => setShowPareto(v => !v)}
          className="w-full flex items-center justify-between px-5 py-3.5 hover:bg-gray-50 transition-colors text-left"
        >
          <div className="flex items-center gap-2.5">
            <svg className="w-4 h-4 text-bx-light-text" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                d="M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 01-2 2h-2a2 2 0 01-2-2z" />
            </svg>
            <span className="font-semibold text-gray-700 text-sm">Pareto Analysis</span>
            <InfoTip text={TIPS.pareto} size="xs" />
            <span className="text-sm text-gray-400">2D objective scatter plot</span>
          </div>
          <svg
            className={`w-4 h-4 text-gray-400 transition-transform duration-150 ${showPareto ? 'rotate-180' : ''}`}
            fill="none" stroke="currentColor" viewBox="0 0 24 24"
          >
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
          </svg>
        </button>
        {showPareto && (
          <div className="px-5 pb-5 border-t border-gray-100 pt-4">
            <ParetoFront
              molecules={filteredMolecules}
              onSelect={mol => mol && handleRowClick(mol)}
              onToggleBookmark={isFrozen ? undefined : toggleBookmark}
            />
          </div>
        )}
      </div>

      {/* Run Logs (collapsible terminal) */}
      <div className="card overflow-hidden">
        <button
          onClick={() => setShowRunLogs(v => !v)}
          className="w-full flex items-center justify-between px-5 py-3 hover:bg-gray-50 transition-colors text-left"
        >
          <div className="flex items-center gap-2.5">
            <svg className="w-4 h-4 text-gray-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                d="M8 9l3 3-3 3m5 0h3M5 20h14a2 2 0 002-2V6a2 2 0 00-2-2H5a2 2 0 00-2 2v12a2 2 0 002 2z" />
            </svg>
            <span className="font-semibold text-gray-700 text-sm">Run Logs</span>
            {currentPhaseRuns.some(r => r.logs?.length > 0) && (
              <span className="text-[10px] text-gray-400 tabular-nums">
                {currentPhaseRuns.reduce((s, r) => s + (r.logs?.length || 0), 0)} entries
              </span>
            )}
          </div>
          <svg
            className={`w-4 h-4 text-gray-400 transition-transform duration-150 ${showRunLogs ? 'rotate-180' : ''}`}
            fill="none" stroke="currentColor" viewBox="0 0 24 24"
          >
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
          </svg>
        </button>
        {showRunLogs && (
          <div className="border-t border-gray-100">
            <div className="bg-gray-900 rounded-b-lg max-h-64 overflow-y-auto px-4 py-3 font-mono text-xs space-y-0.5 scrollbar-thin">
              {(() => {
                const allLogs = currentPhaseRuns
                  .filter(r => r.logs?.length > 0)
                  .flatMap(r => (r.logs || []).map(log => ({
                    ...log,
                    runType: r.type,
                    runId: r.id,
                    ts: log.timestamp ? new Date(log.timestamp).getTime() : 0,
                  })))
                  .sort((a, b) => a.ts - b.ts)

                if (allLogs.length === 0) {
                  return (
                    <div className="text-gray-500 text-center py-6">
                      No run logs available yet. Logs will appear here when runs produce output.
                    </div>
                  )
                }

                return allLogs.map((log, i) => {
                  const levelColor = log.level === 'error' ? 'text-red-400'
                    : log.level === 'warn' || log.level === 'warning' ? 'text-amber-400'
                    : log.level === 'success' ? 'text-emerald-400'
                    : 'text-gray-400'
                  return (
                    <div key={i} className={`flex gap-2 ${levelColor}`}>
                      {log.timestamp && (
                        <span className="text-gray-600 flex-shrink-0 w-[70px]">
                          {new Date(log.timestamp).toLocaleTimeString()}
                        </span>
                      )}
                      <span className="text-gray-500 flex-shrink-0 w-[65px] truncate">[{log.runType || 'run'}]</span>
                      <span className="break-all">{log.message || String(log)}</span>
                    </div>
                  )
                })
              })()}
            </div>
          </div>
        )}
      </div>

      {/* Modals */}
      <RunCreator
        phaseId={phaseId}
        phaseType={phase.type}
        isOpen={showRunCreator}
        onClose={() => setShowRunCreator(false)}
        onSubmit={handleRunSubmit}
        selectedMoleculeIds={selectedMoleculeIds}
        submitting={runSubmitting}
      />

      <FreezeDialog
        isOpen={showFreezeDialog}
        onClose={() => setShowFreezeDialog(false)}
        onConfirm={handleFreezeConfirm}
        action={isFrozen ? 'unfreeze' : 'freeze'}
        phaseName={phase.label}
        bookmarkedMols={flatMolecules.filter(m => m.bookmarked)}
        hasDownstreamRuns={hasDownstreamRuns}
      />

      {/* AI Agent panel */}
      <CampaignAgentPanel
        isOpen={showAgentPanel}
        onClose={() => setShowAgentPanel(false)}
        project={currentProject}
        campaign={currentCampaign}
        phase={currentPhase}
        molecules={flatMolecules}
      />

      {/* Export modal */}
      <ExportModal
        isOpen={showExportModal}
        onClose={() => setShowExportModal(false)}
        molecules={filteredMolecules}
        selectedIds={selectedMoleculeIds}
        columns={visibleColumns}
        onToast={addToast}
      />

      {/* Detail popups (Safety / Synthesis / Confidence) */}
      {popupState && (
        <DetailPopupModal
          type={popupState.type}
          molecule={popupState.molecule}
          onClose={() => setPopupState(null)}
        />
      )}
    </div>
  )
}

// ---------------------------------------------------------------------------
// FilterBar wrapper that also tracks active filter count
// ---------------------------------------------------------------------------
function FilterBarWithCount({ molecules, columns, onFilteredChange, onFilterCountChange }) {
  const [activeQuickFilters] = useState(new Set())

  const handleFilteredChange = useCallback((filtered) => {
    onFilteredChange(filtered)
    // Estimate filter count via difference (rough heuristic)
    onFilterCountChange(filtered.length < molecules.length ? 1 : 0)
  }, [onFilteredChange, onFilterCountChange, molecules.length])

  return (
    <FilterBar
      molecules={molecules}
      columns={columns}
      onFilteredChange={handleFilteredChange}
    />
  )
}

// ---------------------------------------------------------------------------
// Breadcrumb
// ---------------------------------------------------------------------------
function Breadcrumb({ projectId, projectName, phase, phaseTypeMeta }) {
  return (
    <nav className="flex items-center gap-1.5 text-sm text-gray-400" aria-label="Breadcrumb">
      <Link to="/" className="hover:text-bx-mint transition-colors">Projects</Link>
      <svg className="w-3 h-3 text-gray-300" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
      </svg>
      <Link to={`/project/${projectId}`} className="hover:text-bx-light-text transition-colors truncate max-w-[120px]">
        {projectName || 'Project'}
      </Link>
      <svg className="w-3 h-3 text-gray-300" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
      </svg>
      <span className="text-gray-600 font-medium truncate">
        {phase.label} — {phaseTypeMeta.label || phase.type}
      </span>
    </nav>
  )
}

// ---------------------------------------------------------------------------
// Phase Header card
// ---------------------------------------------------------------------------
function PhaseHeader({ phase, phaseTypeMeta, isFrozen, stats, onFreezeToggle, onNewRun, hasActiveRun, onOpenAgent }) {
  return (
    <div className="card overflow-hidden">
      <div className={`h-1.5 ${isFrozen ? 'bg-blue-400' : 'bg-bx-mint'}`} />
      <div className="px-5 py-4">
        <div className="flex flex-wrap items-start justify-between gap-4">
          {/* Title + stats */}
          <div className="space-y-3">
            <div className="flex items-center gap-2.5">
              <h1 className="text-xl font-bold text-bx-light-text">
                {phase.label}
                <span className="ml-2 text-base font-normal text-gray-400">
                  — {phaseTypeMeta.label || phase.type}
                </span>
              </h1>
              {isFrozen && <Badge variant="blue" size="sm">Frozen</Badge>}
            </div>
            <StatsBar stats={stats} />
          </div>

          {/* Action buttons */}
          <div className="flex items-center gap-2">
            <button
              onClick={onOpenAgent}
              className="flex items-center gap-1.5 px-3 py-2 rounded-lg text-sm font-medium border border-emerald-200 text-emerald-700 bg-emerald-50 hover:bg-emerald-100 transition-all duration-150"
            >
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                  d="M9.75 17L9 20l-1 1h8l-1-1-.75-3M3 13h18M5 17H3a2 2 0 01-2-2V5a2 2 0 012-2h14a2 2 0 012 2v10a2 2 0 01-2 2h-2" />
              </svg>
              AI Agent
            </button>
            <button
              onClick={onFreezeToggle}
              className={`flex items-center gap-1.5 px-3 py-2 rounded-lg text-sm font-medium border transition-all duration-150 ${
                isFrozen
                  ? 'border-blue-300 text-blue-700 bg-blue-50 hover:bg-blue-100'
                  : 'border-gray-200 text-gray-600 hover:bg-gray-50 hover:border-gray-300'
              }`}
            >
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                  d={isFrozen
                    ? "M8 11V7a4 4 0 118 0m-4 8v2m-6 4h12a2 2 0 002-2v-6a2 2 0 00-2-2H6a2 2 0 00-2 2v6a2 2 0 002 2z"
                    : "M12 15v2m-6 4h12a2 2 0 002-2v-6a2 2 0 00-2-2H6a2 2 0 00-2 2v6a2 2 0 002 2zm10-10V7a4 4 0 00-8 0v4h8z"}
                />
              </svg>
              {isFrozen ? 'Unfreeze' : 'Freeze Phase'}
            </button>
            <button
              onClick={onNewRun}
              disabled={isFrozen}
              className={`flex items-center gap-1.5 px-4 py-2 rounded-lg text-sm font-semibold transition-all duration-150 ${
                isFrozen
                  ? 'bg-gray-100 text-gray-400 cursor-not-allowed'
                  : hasActiveRun
                    ? 'bg-gray-100 text-gray-400 cursor-not-allowed'
                    : 'bg-bx-surface hover:bg-bx-elevated text-white shadow-sm hover:shadow-md'
              }`}
            >
              {hasActiveRun ? (
                <>
                  <BindXLogo variant="loading" size={16} />
                  Run in progress
                </>
              ) : (
                <>
                  <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 4v16m8-8H4" />
                  </svg>
                  New Run
                </>
              )}
            </button>
          </div>
        </div>

        {isFrozen && (
          <div className="mt-4">
            <FreezeBanner />
          </div>
        )}
      </div>
    </div>
  )
}
