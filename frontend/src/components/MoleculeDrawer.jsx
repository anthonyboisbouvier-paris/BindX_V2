import React, { useState, useEffect, useCallback } from 'react'
import Badge from './Badge.jsx'
import { MOLECULE_DETAILS } from '../mock/data.js'

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
function fmt(v, decimals = 2) {
  if (v == null) return '—'
  if (typeof v === 'boolean') return v ? 'Yes' : 'No'
  if (typeof v === 'number') return Number(v).toFixed(decimals)
  return String(v)
}

function pct(v) {
  if (v == null) return null
  return Math.round(Number(v) * 100)
}

/** Return a percentile rank (0-100) of `value` within `arr` (ascending or descending). */
function percentileRank(value, arr, lowerIsBetter = false) {
  if (value == null || !arr || arr.length === 0) return null
  const valid = arr.filter(v => v != null)
  if (valid.length === 0) return null
  const sorted = [...valid].sort((a, b) => a - b)
  const rank = sorted.filter(v => (lowerIsBetter ? v <= value : v >= value)).length
  return Math.round((rank / sorted.length) * 100)
}

// Drug-like range checks per property
const PROP_RANGES = {
  MW:   { min: 0,   max: 500, warn: 500, color: v => v <= 500 ? 'green' : v <= 600 ? 'yellow' : 'red' },
  logP: { min: -2,  max: 5,   warn: 5,   color: v => v <= 5 && v >= -2 ? 'green' : v <= 6 ? 'yellow' : 'red' },
  HBD:  { min: 0,   max: 5,              color: v => v <= 5 ? 'green' : v <= 7 ? 'yellow' : 'red' },
  HBA:  { min: 0,   max: 10,             color: v => v <= 10 ? 'green' : v <= 12 ? 'yellow' : 'red' },
  TPSA: { min: 0,   max: 140,            color: v => v <= 140 ? 'green' : v <= 160 ? 'yellow' : 'red' },
  QED:  { min: 0.6, max: 1,              color: v => v >= 0.6 ? 'green' : v >= 0.4 ? 'yellow' : 'red' },
}

const PROP_COLORS = {
  green:  { bg: 'bg-green-50',  text: 'text-green-700',  border: 'border-green-100' },
  yellow: { bg: 'bg-yellow-50', text: 'text-yellow-700', border: 'border-yellow-100' },
  red:    { bg: 'bg-red-50',    text: 'text-red-700',    border: 'border-red-100' },
  gray:   { bg: 'bg-gray-50',   text: 'text-gray-600',   border: 'border-gray-100' },
}

function propColorClass(key, value) {
  const def = PROP_RANGES[key]
  if (!def || value == null) return 'gray'
  return def.color(value)
}

// ---------------------------------------------------------------------------
// Interaction type config
// ---------------------------------------------------------------------------
const INTERACTION_COLORS = {
  'HBond':       { dot: 'bg-green-500',  bar: 'bg-green-400',  text: 'text-green-700',  label: 'H-Bond' },
  'Covalent':    { dot: 'bg-red-500',    bar: 'bg-red-400',    text: 'text-red-700',    label: 'Covalent' },
  'Salt Bridge': { dot: 'bg-blue-500',   bar: 'bg-blue-400',   text: 'text-blue-700',   label: 'Salt Bridge' },
  'Hydrophobic': { dot: 'bg-orange-400', bar: 'bg-orange-400', text: 'text-orange-700', label: 'Hydrophobic' },
  'Pi-Stacking': { dot: 'bg-purple-500', bar: 'bg-purple-400', text: 'text-purple-700', label: 'Pi-Stack' },
  'Water Bridge':{ dot: 'bg-gray-400',   bar: 'bg-gray-400',   text: 'text-gray-600',   label: 'Water' },
}

function interactionConfig(type) {
  return INTERACTION_COLORS[type] || { dot: 'bg-gray-400', bar: 'bg-gray-400', text: 'text-gray-600', label: type }
}

// ---------------------------------------------------------------------------
// Sub-components
// ---------------------------------------------------------------------------

function TabButton({ active, onClick, children }) {
  return (
    <button
      onClick={onClick}
      className={`px-3 py-2 text-sm font-semibold whitespace-nowrap transition-colors border-b-2 ${
        active
          ? 'border-bx-surface text-bx-light-text'
          : 'border-transparent text-gray-500 hover:text-gray-700 hover:border-gray-300'
      }`}
    >
      {children}
    </button>
  )
}

function NoDataCard({ runType, message }) {
  return (
    <div className="flex flex-col items-center justify-center py-10 text-gray-300">
      <svg className="w-10 h-10 mb-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.2}
          d="M9 17v-2m3 2v-4m3 4v-6m2 10H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
      </svg>
      <p className="text-sm text-gray-400 font-medium">{message || `No data available`}</p>
      {runType && (
        <p className="text-sm text-gray-300 mt-1">Run a <span className="font-semibold">{runType}</span> analysis to generate this data</p>
      )}
    </div>
  )
}

function MiniBar({ value, max, colorClass = 'bg-blue-500' }) {
  const pct = max > 0 ? Math.min(100, (value / max) * 100) : 0
  return (
    <div className="w-full bg-gray-100 rounded-full h-1.5 overflow-hidden">
      <div className={`h-full rounded-full transition-all ${colorClass}`} style={{ width: `${pct}%` }} />
    </div>
  )
}

function ProgressBar({ value, label, description, colorClass = 'bg-blue-500' }) {
  return (
    <div className="space-y-1">
      <div className="flex items-center justify-between text-sm">
        <span className="font-medium text-gray-700">{label}</span>
        <span className="text-gray-500">{value}%</span>
      </div>
      <div className="w-full bg-gray-100 rounded-full h-2 overflow-hidden">
        <div className={`h-full rounded-full ${colorClass}`} style={{ width: `${value}%` }} />
      </div>
      {description && <p className="text-[10px] text-gray-400">{description}</p>}
    </div>
  )
}

// ---------------------------------------------------------------------------
// Tab 1 — Scores
// ---------------------------------------------------------------------------
function TabScores({ molecule, allMolecules }) {
  const scores = [
    {
      key: 'docking_score',
      label: 'Docking Score',
      unit: 'kcal/mol',
      color: { bg: 'bg-blue-50', border: 'border-blue-100', value: 'text-bx-light-text', accent: 'bg-blue-500' },
      lowerIsBetter: true,
      decimals: 1,
    },
    {
      key: 'cnn_score',
      label: 'CNN Score',
      unit: '0–1',
      color: { bg: 'bg-green-50', border: 'border-green-100', value: 'text-green-700', accent: 'bg-green-500' },
      lowerIsBetter: false,
      decimals: 2,
    },
    {
      key: 'cnn_vs',
      label: 'CNN VS',
      unit: 'pKd',
      color: { bg: 'bg-purple-50', border: 'border-purple-100', value: 'text-purple-700', accent: 'bg-purple-500' },
      lowerIsBetter: false,
      decimals: 2,
    },
    {
      key: 'composite_score',
      label: 'Composite Score',
      unit: '0–100',
      color: { bg: 'bg-amber-50', border: 'border-amber-100', value: 'text-amber-700', accent: 'bg-amber-500' },
      lowerIsBetter: false,
      decimals: 1,
    },
  ]

  const allValues = {}
  if (allMolecules) {
    scores.forEach(s => {
      allValues[s.key] = allMolecules.map(m => m[s.key]).filter(v => v != null)
    })
  }

  const phaseAvg = {}
  scores.forEach(s => {
    const vals = allValues[s.key]
    if (vals && vals.length > 0) {
      phaseAvg[s.key] = vals.reduce((a, b) => a + b, 0) / vals.length
    }
  })

  const hasAny = scores.some(s => molecule[s.key] != null)
  if (!hasAny) {
    return <NoDataCard runType="docking or scoring" message="No scores computed yet" />
  }

  return (
    <div className="space-y-4">
      <div className="grid grid-cols-2 gap-3">
        {scores.map(s => {
          const val = molecule[s.key]
          const rank = val != null ? percentileRank(val, allValues[s.key], s.lowerIsBetter) : null
          const c = s.color
          return (
            <div key={s.key} className={`rounded-xl border ${c.border} ${c.bg} p-3`}>
              {val == null ? (
                <div className="text-center text-gray-300">
                  <p className="text-sm text-gray-400 font-medium mb-1">{s.label}</p>
                  <p className="text-sm text-gray-300">Not computed</p>
                </div>
              ) : (
                <>
                  <p className="text-sm text-gray-500 font-medium mb-1">{s.label}</p>
                  <p className={`text-xl font-bold tabular-nums ${c.value}`}>{fmt(val, s.decimals)}</p>
                  <p className="text-[10px] text-gray-400 mb-2">{s.unit}</p>
                  {rank != null && (
                    <div className="space-y-1">
                      <MiniBar value={rank} max={100} colorClass={c.accent} />
                      <p className="text-[10px] text-gray-500">Top {rank}% in phase</p>
                    </div>
                  )}
                </>
              )}
            </div>
          )
        })}
      </div>

      {/* Bar chart vs phase average */}
      {Object.keys(phaseAvg).length > 0 && (
        <div className="bg-gray-50 rounded-xl p-3 space-y-2">
          <p className="text-sm font-semibold text-gray-500 uppercase tracking-wide mb-3">
            vs Phase Average
          </p>
          {scores.map(s => {
            const val = molecule[s.key]
            const avg = phaseAvg[s.key]
            if (val == null || avg == null) return null
            const ref = s.lowerIsBetter ? Math.abs(val) : val
            const refAvg = s.lowerIsBetter ? Math.abs(avg) : avg
            const maxVal = Math.max(ref, refAvg) * 1.2 || 1
            return (
              <div key={s.key} className="space-y-1">
                <div className="flex items-center justify-between text-[10px] text-gray-500">
                  <span className="font-medium">{s.label}</span>
                  <span className="text-gray-400">avg {fmt(avg, s.decimals)}</span>
                </div>
                <div className="flex items-center gap-2">
                  <div className="flex-1 bg-gray-200 rounded-full h-2 overflow-hidden">
                    <div
                      className={s.color.accent + ' h-full rounded-full'}
                      style={{ width: `${(ref / maxVal) * 100}%` }}
                    />
                  </div>
                  <span className={`text-sm font-bold tabular-nums ${s.color.value}`}>
                    {fmt(val, s.decimals)}
                  </span>
                </div>
                {/* avg marker */}
                <div className="flex items-center gap-2">
                  <div className="flex-1 bg-gray-200 rounded-full h-1 overflow-hidden relative">
                    <div
                      className="bg-gray-400 h-full rounded-full"
                      style={{ width: `${(refAvg / maxVal) * 100}%` }}
                    />
                  </div>
                  <span className="text-[10px] text-gray-400 w-8">avg</span>
                </div>
              </div>
            )
          }).filter(Boolean)}
        </div>
      )}
    </div>
  )
}

// ---------------------------------------------------------------------------
// Tab 2 — Properties
// ---------------------------------------------------------------------------
function TabProperties({ molecule }) {
  const props = [
    { key: 'MW',   label: 'MW',   unit: 'Da',   decimals: 1 },
    { key: 'logP', label: 'LogP', unit: '',      decimals: 2 },
    { key: 'HBD',  label: 'HBD',  unit: '',      decimals: 0 },
    { key: 'HBA',  label: 'HBA',  unit: '',      decimals: 0 },
    { key: 'TPSA', label: 'TPSA', unit: 'A2',    decimals: 1 },
    { key: 'QED',  label: 'QED',  unit: '',      decimals: 2 },
  ]

  const hasProps = props.some(p => molecule[p.key] != null)
  if (!hasProps) {
    return <NoDataCard runType="ADMET" message="No physicochemical data available" />
  }

  // Lipinski
  const lipinskiChecks = [
    { label: 'MW < 500', pass: molecule.MW != null && molecule.MW < 500, value: molecule.MW != null ? `${fmt(molecule.MW, 1)} Da` : '—' },
    { label: 'LogP < 5',  pass: molecule.logP != null && molecule.logP < 5, value: molecule.logP != null ? fmt(molecule.logP) : '—' },
    { label: 'HBD <= 5',  pass: molecule.HBD != null && molecule.HBD <= 5, value: molecule.HBD != null ? molecule.HBD : '—' },
    { label: 'HBA <= 10', pass: molecule.HBA != null && molecule.HBA <= 10, value: molecule.HBA != null ? molecule.HBA : '—' },
  ]
  const lipinskiPossible = lipinskiChecks.every(c => c.pass !== false || molecule[c.label.split(' ')[0]] == null)
  const lipinskiPass = molecule.lipinski_pass != null
    ? molecule.lipinski_pass
    : lipinskiChecks.filter(c => c.pass !== null).every(c => c.pass)

  return (
    <div className="space-y-4">
      <div className="grid grid-cols-3 gap-2">
        {props.map(p => {
          const val = molecule[p.key]
          const colorKey = propColorClass(p.key, val)
          const c = PROP_COLORS[colorKey]
          return (
            <div key={p.key} className={`rounded-lg border ${c.border} ${c.bg} px-3 py-2`}>
              <p className="text-[10px] text-gray-400 font-medium uppercase tracking-wide">{p.label}</p>
              <p className={`text-sm font-bold tabular-nums mt-0.5 ${c.text}`}>
                {val != null ? fmt(val, p.decimals) : '—'}
                {val != null && p.unit ? <span className="text-[10px] font-normal ml-0.5">{p.unit}</span> : null}
              </p>
            </div>
          )
        })}
      </div>

      {/* Lipinski card */}
      {molecule.lipinski_pass != null && (
        <div className={`rounded-xl border px-4 py-3 ${
          lipinskiPass
            ? 'bg-green-50 border-green-100'
            : 'bg-red-50 border-red-100'
        }`}>
          <div className="flex items-center gap-2 mb-2">
            <span className={`text-sm font-bold ${lipinskiPass ? 'text-green-700' : 'text-red-700'}`}>
              Lipinski Rule of 5:
            </span>
            <span className={`px-2 py-0.5 rounded-full text-sm font-bold ${
              lipinskiPass ? 'bg-green-100 text-green-700' : 'bg-red-100 text-red-700'
            }`}>
              {lipinskiPass ? 'PASS' : 'FAIL'}
            </span>
          </div>
          <div className="space-y-1 pl-2 border-l-2 border-gray-200">
            {lipinskiChecks.map((c, i) => (
              <div key={i} className="flex items-center gap-2 text-sm">
                <span className={c.pass ? 'text-green-500' : 'text-red-500'}>
                  {c.pass ? (
                    <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M5 13l4 4L19 7" />
                    </svg>
                  ) : (
                    <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M6 18L18 6M6 6l12 12" />
                    </svg>
                  )}
                </span>
                <span className="text-gray-600">{c.label}:</span>
                <span className="font-semibold text-gray-800">{c.value}</span>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* Additional */}
      {(molecule.scaffold || molecule.cluster_id != null || molecule.interactions_count != null) && (
        <div className="bg-gray-50 rounded-xl p-3 space-y-2">
          <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wide">Additional</p>
          <div className="grid grid-cols-2 gap-2 text-sm">
            {molecule.scaffold && (
              <div>
                <span className="text-gray-400">Scaffold</span>
                <p className="font-semibold text-gray-700 capitalize">{molecule.scaffold}</p>
              </div>
            )}
            {molecule.cluster_id != null && (
              <div>
                <span className="text-gray-400">Cluster</span>
                <p className="font-semibold text-gray-700">#{molecule.cluster_id}</p>
              </div>
            )}
            {molecule.interactions_count != null && (
              <div>
                <span className="text-gray-400">Contacts</span>
                <p className="font-semibold text-gray-700">{molecule.interactions_count}</p>
              </div>
            )}
            {molecule.solubility != null && (
              <div>
                <span className="text-gray-400">Solubility</span>
                <p className="font-semibold text-gray-700">{pct(molecule.solubility)}%</p>
              </div>
            )}
          </div>
        </div>
      )}
    </div>
  )
}

// ---------------------------------------------------------------------------
// Tab 3 — ADMET (custom radar + breakdown bars)
// ---------------------------------------------------------------------------
function RadarSVG({ values, color }) {
  const N = values.length
  const CX = 90, CY = 90, R = 68
  const LABELS = ['Abs', 'Dist', 'Met', 'Exc', 'Tox', 'PD']
  const RINGS = [0.25, 0.5, 0.75, 1.0]

  function pt(v, i) {
    const a = (Math.PI * 2 * i) / N - Math.PI / 2
    return { x: CX + v * R * Math.cos(a), y: CY + v * R * Math.sin(a) }
  }

  const pts = values.map((v, i) => pt(v, i))
  const polygon = pts.map(p => `${p.x.toFixed(1)},${p.y.toFixed(1)}`).join(' ')

  return (
    <svg width={180} height={180} viewBox="0 0 180 180" className="overflow-visible">
      {RINGS.map((r, ri) => {
        const rpts = values.map((_, i) => {
          const p = pt(r, i)
          return `${p.x.toFixed(1)},${p.y.toFixed(1)}`
        }).join(' ')
        return (
          <polygon key={ri} points={rpts} fill="none"
            stroke={ri === RINGS.length - 1 ? '#d1d5db' : '#e5e7eb'}
            strokeWidth={ri === RINGS.length - 1 ? 1.5 : 1}
            strokeDasharray={ri < RINGS.length - 1 ? '3,3' : undefined}
          />
        )
      })}
      {values.map((_, i) => {
        const tip = pt(1, i)
        return <line key={i} x1={CX} y1={CY} x2={tip.x} y2={tip.y} stroke="#e5e7eb" strokeWidth={1} />
      })}
      <polygon points={polygon} fill={color} fillOpacity={0.18} stroke={color} strokeWidth={2} strokeLinejoin="round" />
      {pts.map((p, i) => (
        <circle key={i} cx={p.x} cy={p.y} r={4} fill={color} stroke="white" strokeWidth={1.5} />
      ))}
      {values.map((_, i) => {
        const lp = pt(1.32, i)
        const anchor = lp.x < CX - 8 ? 'end' : lp.x > CX + 8 ? 'start' : 'middle'
        return (
          <text key={i} x={lp.x} y={lp.y} fontSize={8} fill="#6b7280" textAnchor={anchor} dominantBaseline="central">
            {LABELS[i]}
          </text>
        )
      })}
    </svg>
  )
}

function TabADMET({ molecule }) {
  const admet = molecule.admet
  if (!admet) {
    return <NoDataCard runType="ADMET" message="No ADMET data available" />
  }

  const axes = [
    { key: 'absorption',       label: 'Absorption',       desc: 'Good oral bioavailability',       color: 'bg-blue-500' },
    { key: 'distribution',     label: 'Distribution',     desc: 'Moderate BBB permeation',         color: 'bg-purple-500' },
    { key: 'metabolism',       label: 'Metabolism',       desc: 'Low CYP inhibition',              color: 'bg-amber-500' },
    { key: 'excretion',        label: 'Excretion',        desc: 'Moderate renal clearance',        color: 'bg-teal-500' },
    { key: 'toxicity',         label: 'Toxicity',         desc: 'Low toxicity risk',               color: 'bg-green-500' },
    { key: 'pharmacodynamics', label: 'Pharmacodynamics', desc: 'Target engagement predicted',     color: 'bg-indigo-500' },
  ]

  const values = axes.map(a => typeof admet[a.key] === 'number' ? admet[a.key] : 0.5)
  const avg = values.reduce((a, b) => a + b, 0) / values.length
  const overallColor = avg >= 0.72 ? '#00e6a0' : avg >= 0.55 ? '#eab308' : '#ef4444'
  const overallLabel = avg >= 0.72 ? 'Good' : avg >= 0.55 ? 'Mixed' : 'Poor'
  const overallVariant = avg >= 0.72 ? 'green' : avg >= 0.55 ? 'yellow' : 'red'

  return (
    <div className="space-y-4">
      <div className="flex flex-col items-center gap-2">
        <RadarSVG values={values} color={overallColor} />
        <div className="text-center">
          <span className={`text-sm font-semibold px-2 py-1 rounded-full ${
            overallVariant === 'green' ? 'bg-green-100 text-green-700' :
            overallVariant === 'yellow' ? 'bg-yellow-100 text-yellow-700' :
            'bg-red-100 text-red-700'
          }`}>
            Overall ADMET: {overallLabel} ({Math.round(avg * 100)}%)
          </span>
        </div>
      </div>

      <div className="space-y-2">
        {axes.map((a, i) => (
          <ProgressBar
            key={a.key}
            value={Math.round(values[i] * 100)}
            label={a.label}
            description={a.desc}
            colorClass={a.color}
          />
        ))}
      </div>

      {(molecule.hERG != null || molecule.metabolic_stability != null || molecule.BBB != null) && (
        <div className="bg-gray-50 rounded-xl p-3 space-y-1">
          <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wide mb-2">Key ADMET Metrics</p>
          {molecule.hERG != null && (
            <div className="flex items-center justify-between text-sm">
              <span className="text-gray-500">hERG Risk</span>
              <span className={`font-bold tabular-nums ${molecule.hERG < 0.2 ? 'text-green-600' : molecule.hERG < 0.4 ? 'text-yellow-600' : 'text-red-600'}`}>
                {pct(molecule.hERG)}% {molecule.hERG < 0.2 ? '(Safe)' : molecule.hERG < 0.4 ? '(Moderate)' : '(High)'}
              </span>
            </div>
          )}
          {molecule.metabolic_stability != null && (
            <div className="flex items-center justify-between text-sm">
              <span className="text-gray-500">Metabolic Stability</span>
              <span className={`font-bold tabular-nums ${molecule.metabolic_stability > 0.6 ? 'text-green-600' : 'text-yellow-600'}`}>
                {pct(molecule.metabolic_stability)}%
              </span>
            </div>
          )}
          {molecule.BBB != null && (
            <div className="flex items-center justify-between text-sm">
              <span className="text-gray-500">BBB Penetration</span>
              <span className="font-bold tabular-nums text-gray-700">{pct(molecule.BBB)}%</span>
            </div>
          )}
          {molecule.solubility != null && (
            <div className="flex items-center justify-between text-sm">
              <span className="text-gray-500">Solubility</span>
              <span className={`font-bold tabular-nums ${molecule.solubility > 0.5 ? 'text-green-600' : 'text-yellow-600'}`}>
                {pct(molecule.solubility)}%
              </span>
            </div>
          )}
        </div>
      )}
    </div>
  )
}

// ---------------------------------------------------------------------------
// Tab 4 — Safety
// ---------------------------------------------------------------------------
function TabSafety({ molecule, details }) {
  const safety = details?.safety
  if (!safety) {
    return <NoDataCard runType="enrichment" message="No safety data available" />
  }

  const riskVariant = r => r === 'high' ? 'red' : r === 'medium' ? 'yellow' : 'green'

  return (
    <div className="space-y-4">
      {/* PAINS + Brenk */}
      <div className="space-y-2">
        <div className={`flex items-center gap-2 p-2 rounded-lg text-sm ${
          safety.pains_pass ? 'bg-green-50' : 'bg-red-50'
        }`}>
          <span className={safety.pains_pass ? 'text-green-500' : 'text-red-500'}>
            {safety.pains_pass ? (
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M5 13l4 4L19 7" />
              </svg>
            ) : (
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M6 18L18 6M6 6l12 12" />
              </svg>
            )}
          </span>
          <span className="font-semibold text-gray-700">PAINS Check:</span>
          <span className={safety.pains_pass ? 'text-green-700 font-semibold' : 'text-red-700 font-semibold'}>
            {safety.pains_pass ? 'PASS (no alerts)' : `FAIL (${safety.pains_alerts?.length || 0} alerts)`}
          </span>
        </div>

        {safety.brenk_alerts && safety.brenk_alerts.length > 0 && (
          <div className="bg-amber-50 border border-amber-100 rounded-lg p-2">
            <p className="text-sm font-semibold text-amber-700 mb-1">
              Brenk Alerts ({safety.brenk_alerts.length})
            </p>
            {safety.brenk_alerts.map((a, i) => (
              <p key={i} className="text-sm text-amber-600 flex items-start gap-1">
                <span className="mt-0.5 flex-shrink-0">--</span>
                {a}
              </p>
            ))}
          </div>
        )}
        {(!safety.brenk_alerts || safety.brenk_alerts.length === 0) && (
          <div className="flex items-center gap-2 text-sm text-green-600">
            <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M5 13l4 4L19 7" />
            </svg>
            Brenk: No structural alerts
          </div>
        )}
      </div>

      {/* Off-target */}
      {safety.off_target && safety.off_target.length > 0 && (
        <div>
          <p className="text-sm font-semibold text-gray-500 uppercase tracking-wide mb-2">Off-Target Screening</p>
          <div className="rounded-xl border border-gray-100 overflow-hidden">
            <table className="w-full text-sm">
              <thead>
                <tr className="bg-gray-50 text-left">
                  <th className="px-3 py-2 font-semibold text-gray-500">Target</th>
                  <th className="px-3 py-2 font-semibold text-gray-500 text-right">Similarity</th>
                  <th className="px-3 py-2 font-semibold text-gray-500 text-center">Risk</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-gray-50">
                {safety.off_target.map((ot, i) => (
                  <tr key={i} className="hover:bg-gray-50">
                    <td className="px-3 py-2 text-gray-700 font-medium">{ot.target}</td>
                    <td className="px-3 py-2 text-right font-mono text-gray-600">
                      {Math.round(ot.similarity * 100)}%
                    </td>
                    <td className="px-3 py-2 text-center">
                      <Badge variant={riskVariant(ot.risk)} size="sm">
                        {ot.risk?.toUpperCase()}
                      </Badge>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}

      {/* Risk bars */}
      <div className="space-y-2">
        {[
          { label: 'hERG Risk',      value: pct(safety.herg_risk),     pass: safety.herg_pass,     threshold: 20 },
          { label: 'Ames Mutagenicity', value: pct(safety.ames_risk),   pass: safety.ames_risk < 0.2,  threshold: 20 },
          { label: 'Hepatotoxicity', value: pct(safety.hepatotox_risk), pass: safety.hepatotox_risk < 0.2, threshold: 20 },
        ].filter(r => r.value != null).map(r => (
          <div key={r.label} className="space-y-1">
            <div className="flex items-center justify-between text-sm">
              <span className="font-medium text-gray-700">{r.label}</span>
              <span className={`font-bold ${r.pass ? 'text-green-600' : 'text-red-600'}`}>
                {r.value}% — {r.pass ? 'Safe' : 'Risk'} (threshold: {r.threshold}%)
              </span>
            </div>
            <div className="w-full bg-gray-100 rounded-full h-2 overflow-hidden">
              <div
                className={`h-full rounded-full ${r.pass ? 'bg-green-400' : 'bg-red-400'}`}
                style={{ width: `${Math.min(100, r.value)}%` }}
              />
            </div>
          </div>
        ))}
      </div>

      {/* Confidence breakdown */}
      {safety.confidence && (
        <div className="bg-bx-surface/5 rounded-xl p-3 space-y-2">
          <div className="flex items-center justify-between mb-1">
            <p className="text-sm font-semibold text-bx-light-text">Confidence</p>
            <span className="text-sm font-bold text-bx-light-text">{pct(safety.confidence.overall)}% overall</span>
          </div>
          {[
            { key: 'binding',     label: 'Binding' },
            { key: 'admet',       label: 'ADMET' },
            { key: 'selectivity', label: 'Selectivity' },
            { key: 'safety',      label: 'Safety' },
          ].map(c => (
            <div key={c.key} className="flex items-center gap-2 text-sm">
              <span className="w-20 text-gray-500 flex-shrink-0">{c.label}</span>
              <div className="flex-1 bg-white/60 rounded-full h-1.5 overflow-hidden">
                <div
                  className="bg-bx-surface h-full rounded-full"
                  style={{ width: `${pct(safety.confidence[c.key]) || 0}%` }}
                />
              </div>
              <span className="w-8 text-right font-semibold text-bx-light-text">
                {pct(safety.confidence[c.key]) || '—'}%
              </span>
            </div>
          ))}
        </div>
      )}
    </div>
  )
}

// ---------------------------------------------------------------------------
// Tab 5 — Synthesis
// ---------------------------------------------------------------------------
function TabSynthesis({ details }) {
  const synthesis = details?.synthesis
  if (!synthesis) {
    return <NoDataCard runType="enrichment" message="No synthesis data available" />
  }

  return (
    <div className="space-y-4">
      {synthesis.steps && synthesis.steps.length > 0 ? (
        <>
          <div className="space-y-3">
            {synthesis.steps.map((step, i) => (
              <div key={i}>
                <div className="bg-white rounded-xl border border-gray-100 p-3 shadow-sm">
                  <div className="flex items-start gap-2 mb-2">
                    <span className="flex-shrink-0 w-5 h-5 rounded-full bg-bx-surface text-white text-[10px] font-bold flex items-center justify-center">
                      {i + 1}
                    </span>
                    <p className="font-semibold text-sm text-gray-800">{step.product}</p>
                  </div>
                  <div className="pl-7 space-y-1 text-sm text-gray-600">
                    <div className="flex items-start gap-1">
                      <span className="text-gray-400 flex-shrink-0">Reagent:</span>
                      <span className="font-medium">{step.reagent}</span>
                    </div>
                    <div className="flex items-start gap-1">
                      <span className="text-gray-400 flex-shrink-0">Conditions:</span>
                      <span className="font-mono text-sm">{step.conditions}</span>
                    </div>
                    <div className="flex items-center gap-2 mt-2">
                      <span className="text-gray-400">Yield:</span>
                      <div className="flex-1 bg-gray-100 rounded-full h-1.5 overflow-hidden">
                        <div
                          className={`h-full rounded-full ${step.yield >= 0.8 ? 'bg-green-500' : step.yield >= 0.6 ? 'bg-amber-400' : 'bg-red-400'}`}
                          style={{ width: `${Math.round(step.yield * 100)}%` }}
                        />
                      </div>
                      <span className={`font-bold tabular-nums ${step.yield >= 0.8 ? 'text-green-600' : step.yield >= 0.6 ? 'text-amber-600' : 'text-red-600'}`}>
                        {Math.round(step.yield * 100)}%
                      </span>
                    </div>
                  </div>
                </div>
                {i < synthesis.steps.length - 1 && (
                  <div className="flex justify-center my-1">
                    <svg className="w-4 h-4 text-gray-300" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
                    </svg>
                  </div>
                )}
              </div>
            ))}
          </div>

          <div className="bg-gray-50 rounded-xl p-3 border border-gray-100">
            <div className="grid grid-cols-3 gap-3 text-center text-sm">
              <div>
                <p className="text-gray-400">Total Steps</p>
                <p className="font-bold text-gray-800 text-sm">{synthesis.num_steps}</p>
              </div>
              <div>
                <p className="text-gray-400">Est. Cost</p>
                <p className="font-bold text-gray-800 text-sm">{synthesis.total_cost}</p>
              </div>
              <div>
                <p className="text-gray-400">Feasibility</p>
                <p className={`font-bold text-sm ${
                  synthesis.feasibility >= 0.8 ? 'text-green-600' :
                  synthesis.feasibility >= 0.6 ? 'text-amber-600' : 'text-red-600'
                }`}>
                  {Math.round(synthesis.feasibility * 100)}%
                </p>
              </div>
            </div>
          </div>
        </>
      ) : (
        <NoDataCard runType="enrichment" message="Synthesis route not available" />
      )}
    </div>
  )
}

// ---------------------------------------------------------------------------
// Tab 6 — Interactions
// ---------------------------------------------------------------------------
function TabInteractions({ details }) {
  const interactions = details?.interactions
  if (!interactions) {
    return <NoDataCard runType="enrichment" message="No interaction data available" />
  }

  const { residues = [], total_contacts, hbonds, hydrophobic, ionic, pi_stacking, water_bridges, covalent } = interactions

  return (
    <div className="space-y-4">
      {/* Summary */}
      <div className="bg-bx-surface/5 rounded-xl p-3">
        <div className="grid grid-cols-3 gap-2 text-center text-sm">
          {[
            { label: 'Total',       value: total_contacts, color: 'text-bx-light-text' },
            { label: 'H-Bonds',     value: hbonds,          color: 'text-green-600' },
            { label: 'Hydrophobic', value: hydrophobic,     color: 'text-orange-500' },
            { label: 'Ionic',       value: ionic,           color: 'text-blue-600' },
            { label: 'Pi-Stack',    value: pi_stacking,     color: 'text-purple-600' },
            { label: 'Water',       value: water_bridges,   color: 'text-gray-500' },
          ].filter(x => x.value != null).map(x => (
            <div key={x.label}>
              <p className={`text-base font-bold tabular-nums ${x.color}`}>{x.value}</p>
              <p className="text-gray-400">{x.label}</p>
            </div>
          ))}
        </div>
        {covalent > 0 && (
          <div className="mt-2 text-center">
            <Badge variant="red" size="sm">{covalent} Covalent</Badge>
          </div>
        )}
      </div>

      {/* Residues table */}
      {residues.length > 0 && (
        <div className="rounded-xl border border-gray-100 overflow-hidden">
          <table className="w-full text-sm">
            <thead>
              <tr className="bg-gray-50 text-left">
                <th className="px-3 py-2 font-semibold text-gray-500">Residue</th>
                <th className="px-3 py-2 font-semibold text-gray-500">Type</th>
                <th className="px-3 py-2 font-semibold text-gray-500 text-right">Dist.</th>
                <th className="px-3 py-2 font-semibold text-gray-500">Chain</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-gray-50">
              {residues.map((r, i) => {
                const cfg = interactionConfig(r.type)
                const distPct = Math.max(0, Math.min(100, ((5 - r.distance) / 3.5) * 100))
                return (
                  <tr key={i} className="hover:bg-gray-50">
                    <td className="px-3 py-2 font-mono font-semibold text-gray-800">{r.name}</td>
                    <td className="px-3 py-2">
                      <div className="flex items-center gap-1.5">
                        <span className={`w-2 h-2 rounded-full flex-shrink-0 ${cfg.dot}`} />
                        <span className={`font-medium ${cfg.text}`}>{cfg.label}</span>
                      </div>
                    </td>
                    <td className="px-3 py-2 text-right">
                      <div className="flex items-center gap-1 justify-end">
                        <div className="w-12 bg-gray-100 rounded-full h-1.5 overflow-hidden">
                          <div className={`h-full rounded-full ${cfg.bar}`} style={{ width: `${distPct}%` }} />
                        </div>
                        <span className="font-mono text-gray-700 w-10 text-right">{r.distance.toFixed(1)} A</span>
                      </div>
                    </td>
                    <td className="px-3 py-2 text-gray-500">{r.chain}</td>
                  </tr>
                )
              })}
            </tbody>
          </table>
        </div>
      )}
    </div>
  )
}

// ---------------------------------------------------------------------------
// MoleculeDrawer — main component
// ---------------------------------------------------------------------------
const TABS = ['Scores', 'Properties', 'ADMET', 'Safety', 'Synthesis', 'Interactions']

export default function MoleculeDrawer({
  molecule,
  isOpen,
  onClose,
  onNext,
  onPrev,
  onToggleBookmark,
  allMolecules,
}) {
  const [tab, setTab] = useState('Scores')
  const [copied, setCopied] = useState(false)

  const details = molecule ? MOLECULE_DETAILS[molecule.id] : null

  // Index label
  const idx = allMolecules && molecule
    ? allMolecules.findIndex(m => m.id === molecule.id) + 1
    : null
  const total = allMolecules ? allMolecules.length : null

  // Keyboard navigation
  useEffect(() => {
    if (!isOpen) return
    function handler(e) {
      if (e.key === 'ArrowRight' && onNext) onNext()
      if (e.key === 'ArrowLeft' && onPrev) onPrev()
      if (e.key === 'Escape') onClose()
    }
    window.addEventListener('keydown', handler)
    return () => window.removeEventListener('keydown', handler)
  }, [isOpen, onNext, onPrev, onClose])

  const handleCopy = useCallback(() => {
    if (!molecule?.smiles) return
    navigator.clipboard.writeText(molecule.smiles).then(() => {
      setCopied(true)
      setTimeout(() => setCopied(false), 1800)
    }).catch(() => {})
  }, [molecule])

  if (!isOpen || !molecule) return null

  return (
    <>
      <div
        className="fixed inset-0 z-30 bg-black/10"
        onClick={onClose}
        aria-hidden="true"
      />

      <div
        className="fixed right-0 top-0 bottom-0 z-40 bg-white flex flex-col overflow-hidden
                   border-l border-gray-200 shadow-2xl transition-transform duration-200"
        style={{ width: '40%', minWidth: '380px', maxWidth: '660px' }}
        onClick={e => e.stopPropagation()}
      >
        {/* ----------------------------------------------------------------- */}
        {/* Header                                                              */}
        {/* ----------------------------------------------------------------- */}
        <div className="flex items-center justify-between px-4 py-3 border-b border-gray-100 bg-bx-surface flex-shrink-0">
          <div className="flex items-center gap-2 min-w-0">
            <h3 className="font-bold text-white text-sm truncate">
              {molecule.name || molecule.id}
            </h3>
            {molecule.bookmarked && (
              <svg className="w-4 h-4 text-yellow-400 fill-yellow-400 flex-shrink-0" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                  d="M5 5a2 2 0 012-2h10a2 2 0 012 2v16l-7-3.5L5 21V5z" />
              </svg>
            )}
            {molecule.generation_level > 0 && (
              <Badge variant="purple" size="sm">Gen {molecule.generation_level}</Badge>
            )}
          </div>
          <div className="flex items-center gap-1 flex-shrink-0">
            {onToggleBookmark && (
              <button
                onClick={() => onToggleBookmark(molecule.id)}
                className="p-1.5 rounded hover:bg-white/10 text-white/70 hover:text-yellow-400 transition-colors"
                title="Toggle bookmark"
              >
                <svg className={`w-4 h-4 ${molecule.bookmarked ? 'fill-yellow-400 text-yellow-400' : ''}`}
                  fill={molecule.bookmarked ? 'currentColor' : 'none'} stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                    d="M5 5a2 2 0 012-2h10a2 2 0 012 2v16l-7-3.5L5 21V5z" />
                </svg>
              </button>
            )}
            <button
              onClick={onClose}
              className="p-1.5 rounded hover:bg-white/10 text-white/70 transition-colors"
              title="Close"
            >
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
              </svg>
            </button>
          </div>
        </div>

        {/* Quick props strip */}
        <div className="flex items-center gap-3 px-4 py-2 bg-gray-50 border-b border-gray-100 flex-shrink-0">
          {molecule.scaffold && <span className="text-sm font-medium text-gray-600 capitalize">{molecule.scaffold}</span>}
          {molecule.scaffold && <span className="text-gray-200">|</span>}
          {molecule.MW != null && <span className="text-sm text-gray-500">MW {fmt(molecule.MW, 1)} Da</span>}
          {molecule.logP != null && <span className="text-sm text-gray-500">LogP {fmt(molecule.logP)}</span>}
          {molecule.lipinski_pass != null && (
            <Badge variant={molecule.lipinski_pass ? 'green' : 'red'} size="sm">
              {molecule.lipinski_pass ? 'Lipinski OK' : 'Lipinski Fail'}
            </Badge>
          )}
        </div>

        {/* ----------------------------------------------------------------- */}
        {/* 3D Viewer placeholder                                               */}
        {/* ----------------------------------------------------------------- */}
        <div className="flex-shrink-0 bg-[#0d1b2a] border-b border-bx-surface/30">
          <div className="h-36 flex flex-col items-center justify-center gap-2 relative">
            <svg className="w-10 h-10 text-bx-light-text" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={0.8}
                d="M14 10l-2 1m0 0l-2-1m2 1v2.5M20 7l-2 1m2-1l-2-1m2 1v2.5M14 4l-2-1-2 1M4 7l2-1M4 7l2 1M4 7v2.5M12 21l-2-1m2 1l2-1m-2 1v-2.5M6 18l-2-1v-2.5M18 18l2-1v-2.5" />
            </svg>
            <p className="text-[11px] text-white/30 text-center">3D viewer — connect backend for structure</p>
            {molecule.smiles && (
              <p className="text-[9px] font-mono text-white/20 text-center max-w-full truncate px-4"
                title={molecule.smiles}>
                {molecule.smiles.slice(0, 70)}{molecule.smiles.length > 70 ? '...' : ''}
              </p>
            )}
            <span className="absolute top-2 right-3 text-[10px] text-white/20">Backend not connected</span>
          </div>
        </div>

        {/* ----------------------------------------------------------------- */}
        {/* Tabs                                                                */}
        {/* ----------------------------------------------------------------- */}
        <div className="flex border-b border-gray-100 overflow-x-auto flex-shrink-0 bg-white">
          {TABS.map(t => (
            <TabButton key={t} active={tab === t} onClick={() => setTab(t)}>
              {t}
            </TabButton>
          ))}
        </div>

        {/* ----------------------------------------------------------------- */}
        {/* Tab content                                                         */}
        {/* ----------------------------------------------------------------- */}
        <div className="flex-1 overflow-y-auto p-4">
          {tab === 'Scores'       && <TabScores molecule={molecule} allMolecules={allMolecules} />}
          {tab === 'Properties'   && <TabProperties molecule={molecule} />}
          {tab === 'ADMET'        && <TabADMET molecule={molecule} />}
          {tab === 'Safety'       && <TabSafety molecule={molecule} details={details} />}
          {tab === 'Synthesis'    && <TabSynthesis details={details} />}
          {tab === 'Interactions' && <TabInteractions details={details} />}
        </div>

        {/* ----------------------------------------------------------------- */}
        {/* SMILES + Navigation footer                                          */}
        {/* ----------------------------------------------------------------- */}
        <div className="flex-shrink-0 border-t border-gray-100 bg-gray-50 px-4 py-3 space-y-3">
          {molecule.smiles && (
            <div className="flex items-start gap-2">
              <div className="flex-1 min-w-0">
                <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wide mb-0.5">SMILES</p>
                <p className="font-mono text-[10px] text-gray-600 break-all leading-relaxed line-clamp-2">
                  {molecule.smiles}
                </p>
              </div>
              <button
                onClick={handleCopy}
                className={`flex-shrink-0 px-2 py-1.5 text-[10px] font-semibold rounded-lg transition-colors ${
                  copied
                    ? 'bg-green-100 text-green-700'
                    : 'bg-white border border-gray-200 text-gray-600 hover:bg-gray-100'
                }`}
              >
                {copied ? 'Copied!' : 'Copy'}
              </button>
            </div>
          )}

          <div className="flex items-center justify-between">
            <button
              onClick={onPrev}
              disabled={!onPrev}
              className="flex items-center gap-1.5 px-3 py-1.5 rounded-lg text-sm font-medium
                         border border-gray-200 text-gray-700 hover:bg-gray-100 disabled:opacity-30
                         disabled:cursor-not-allowed transition-colors"
            >
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 19l-7-7 7-7" />
              </svg>
              Prev
            </button>

            {idx != null && total != null && (
              <span className="text-sm text-gray-400 tabular-nums">{idx} of {total}</span>
            )}

            <button
              onClick={onNext}
              disabled={!onNext}
              className="flex items-center gap-1.5 px-3 py-1.5 rounded-lg text-sm font-medium
                         border border-gray-200 text-gray-700 hover:bg-gray-100 disabled:opacity-30
                         disabled:cursor-not-allowed transition-colors"
            >
              Next
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
              </svg>
            </button>
          </div>
        </div>
      </div>
    </>
  )
}
