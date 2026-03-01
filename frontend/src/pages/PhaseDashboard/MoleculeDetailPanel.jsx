import React, { useState, useRef, useEffect, useMemo, useCallback } from 'react'
import Badge from '../../components/Badge.jsx'
import ProteinViewer from '../../components/ProteinViewer.jsx'
import LigandViewer3D from '../../components/LigandViewer3D.jsx'
import SafetyReport from '../../components/SafetyReport.jsx'
import SynthesisTree from '../../components/SynthesisTree.jsx'
import ConfidenceBreakdown from '../../components/ConfidenceBreakdown.jsx'
import InfoTip, { TIPS } from '../../components/InfoTip.jsx'

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
  const [viewerMode, setViewerMode] = useState(null) // null = auto-detect on render

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

        {/* 3D Viewer with Protein / Ligand toggle */}
        {(() => {
          const hasProtein = !!pdbUrl
          const hasLigand = !!molecule.smiles
          // Auto-detect default mode
          const effectiveMode = viewerMode || (hasProtein ? 'protein' : hasLigand ? 'ligand' : null)

          return (
            <div>
              {/* Toggle buttons — show only if both views are available */}
              {hasProtein && hasLigand && (
                <div className="flex gap-1 mb-2">
                  <button
                    onClick={() => setViewerMode('protein')}
                    className={`flex items-center gap-1.5 px-3 py-1.5 rounded-lg text-[11px] font-semibold transition-all ${
                      effectiveMode === 'protein'
                        ? 'bg-bx-surface text-white shadow-sm'
                        : 'bg-gray-100 text-gray-500 hover:bg-gray-200'
                    }`}
                  >
                    <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                        d="M14 10l-2 1m0 0l-2-1m2 1v2.5M20 7l-2 1m2-1l-2-1m2 1v2.5M14 4l-2-1-2 1M4 7l2-1M4 7l2 1M4 7v2.5M12 21l-2-1m2 1l2-1m-2 1v-2.5M6 18l-2-1v-2.5M18 18l2-1v-2.5" />
                    </svg>
                    Protein
                  </button>
                  <button
                    onClick={() => setViewerMode('ligand')}
                    className={`flex items-center gap-1.5 px-3 py-1.5 rounded-lg text-[11px] font-semibold transition-all ${
                      effectiveMode === 'ligand'
                        ? 'bg-bx-surface text-white shadow-sm'
                        : 'bg-gray-100 text-gray-500 hover:bg-gray-200'
                    }`}
                  >
                    <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                        d="M9.75 3.104v5.714a2.25 2.25 0 01-.659 1.591L5 14.5M9.75 3.104c-.251.023-.501.05-.75.082m.75-.082a24.301 24.301 0 014.5 0m0 0v5.714a2.25 2.25 0 00.659 1.591L19 14.5M14.25 3.104c.251.023.501.05.75.082M19 14.5l-2.47 2.47a3.749 3.749 0 01-5.06 0L9 14.5m10 0a4.5 4.5 0 01-9 0" />
                    </svg>
                    Ligand
                  </button>
                </div>
              )}

              {/* Viewer content */}
              {effectiveMode === 'protein' && hasProtein && (
                <ProteinViewer
                  pdbUrl={pdbUrl}
                  selectedPocket={selectedPocket}
                  height={280}
                />
              )}

              {effectiveMode === 'ligand' && hasLigand && (
                <LigandViewer3D smiles={molecule.smiles} height={280} />
              )}

              {/* Fallback: no protein AND no ligand */}
              {!effectiveMode && (
                <div className="rounded-xl bg-[#0e1628] border border-bx-surface/30 overflow-hidden">
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
            </div>
          )
        })()}

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

export default MoleculeDetailPanel
export { DetailPopupModal }
