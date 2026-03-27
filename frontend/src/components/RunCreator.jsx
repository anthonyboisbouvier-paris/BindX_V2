import React, { useState, useRef, useEffect, lazy, Suspense } from 'react'
import { RUN_TYPES, CALCULATION_SUBTYPES, GROUP_META, ESTIMATED_TIMES, COLUMN_DEFAULTS, ALL_COLUMNS, getDefaultCheckedColumns, isColumnAvailable, getAvailableColumns } from '../lib/columns.js'
import ScoringWeightsEditor from './ScoringWeightsEditor.jsx'
import AFVSConfigForm from './AFVSConfigForm.jsx'
import Badge from './Badge.jsx'
import BindXLogo from './BindXLogo.jsx'
import InfoTip, { TIPS } from './InfoTip.jsx'
import { v9GpuHealth } from '../api.js'

const ScaffoldAnalyzer = lazy(() => import('./ScaffoldAnalyzer.jsx'))

// ---------------------------------------------------------------------------
// Run type icon map
// ---------------------------------------------------------------------------
function RunTypeIcon({ icon, className = 'w-5 h-5' }) {
  const p = { className, fill: 'none', stroke: 'currentColor', viewBox: '0 0 24 24' }
  switch (icon) {
    case 'upload':
      return <svg {...p}><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
        d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-8l-4-4m0 0L8 8m4-4v12" /></svg>
    case 'target':
      return <svg {...p}><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
        d="M9 12.75L11.25 15 15 9.75M21 12c0 4.97-4.03 9-9 9s-9-4.03-9-9 4.03-9 9-9 9 4.03 9 9z" /></svg>
    case 'shield':
      return <svg {...p}><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
        d="M9 12l2 2 4-4m5.618-4.016A11.955 11.955 0 0112 2.944a11.955 11.955 0 01-8.618 3.04A12.02 12.02 0 003 9c0 5.591 3.824 10.29 9 11.622 5.176-1.332 9-6.03 9-11.622 0-1.042-.133-2.052-.382-3.016z" /></svg>
    case 'star':
      return <svg {...p}><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
        d="M11.049 2.927c.3-.921 1.603-.921 1.902 0l1.519 4.674a1 1 0 00.95.69h4.915c.969 0 1.371 1.24.588 1.81l-3.976 2.888a1 1 0 00-.363 1.118l1.518 4.674c.3.922-.755 1.688-1.538 1.118l-3.976-2.888a1 1 0 00-1.176 0l-3.976 2.888c-.783.57-1.838-.197-1.538-1.118l1.518-4.674a1 1 0 00-.363-1.118l-3.976-2.888c-.784-.57-.38-1.81.588-1.81h4.914a1 1 0 00.951-.69l1.519-4.674z" /></svg>
    case 'layers':
      return <svg {...p}><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
        d="M12 6.042A8.967 8.967 0 006 3.75c-1.052 0-2.062.18-3 .512v14.25A8.987 8.987 0 016 18c2.305 0 4.408.867 6 2.292m0-14.25a8.966 8.966 0 016-2.292c1.052 0 2.062.18 3 .512v14.25A8.987 8.987 0 0018 18a8.967 8.967 0 00-6 2.292m0-14.25v14.25" /></svg>
    case 'sparkles':
      return <svg {...p}><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
        d="M9.813 15.904L9 18.75l-.813-2.846a4.5 4.5 0 00-3.09-3.09L2.25 12l2.846-.813a4.5 4.5 0 003.09-3.09L9 5.25l.813 2.846a4.5 4.5 0 003.09 3.09L15.75 12l-2.846.813a4.5 4.5 0 00-3.09 3.09zM18.259 8.715L18 9.75l-.259-1.035a3.375 3.375 0 00-2.455-2.456L14.25 6l1.036-.259a3.375 3.375 0 002.455-2.456L18 2.25l.259 1.035a3.375 3.375 0 002.456 2.456L21.75 6l-1.035.259a3.375 3.375 0 00-2.456 2.456z" /></svg>
    case 'trophy':
      return <svg {...p}><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
        d="M16.5 18.75h-9m9 0a3 3 0 013 3h-15a3 3 0 013-3m9 0v-4.5A3.375 3.375 0 0013.125 10.875h-2.25A3.375 3.375 0 007.5 14.25v4.5m9-13.5V3.375c0-.621-.504-1.125-1.125-1.125h-6.75C8.004 2.25 7.5 2.754 7.5 3.375V5.25m9 0h1.875c.621 0 1.125.504 1.125 1.125v2.25c0 1.035-.84 1.875-1.875 1.875H16.5m-9-5.25H5.625C5.004 5.25 4.5 5.754 4.5 6.375v2.25c0 1.035.84 1.875 1.875 1.875H7.5" /></svg>
    case 'grid':
      return <svg {...p}><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
        d="M3.75 6A2.25 2.25 0 016 3.75h2.25A2.25 2.25 0 0110.5 6v2.25a2.25 2.25 0 01-2.25 2.25H6a2.25 2.25 0 01-2.25-2.25V6zM3.75 15.75A2.25 2.25 0 016 13.5h2.25a2.25 2.25 0 012.25 2.25V18a2.25 2.25 0 01-2.25 2.25H6A2.25 2.25 0 013.75 18v-2.25zM13.5 6a2.25 2.25 0 012.25-2.25H18A2.25 2.25 0 0120.25 6v2.25A2.25 2.25 0 0118 10.5h-2.25a2.25 2.25 0 01-2.25-2.25V6zM13.5 15.75a2.25 2.25 0 012.25-2.25H18a2.25 2.25 0 012.25 2.25V18A2.25 2.25 0 0118 20.25h-2.25A2.25 2.25 0 0113.5 18v-2.25z" /></svg>
    case 'calculator':
      return <svg {...p}><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
        d="M15.75 15.75V18m-7.5-6.75h.008v.008H8.25v-.008zm0 2.25h.008v.008H8.25V13.5zm0 2.25h.008v.008H8.25v-.008zm0 2.25h.008v.008H8.25V18zm2.498-6.75h.008v.008h-.008v-.008zm0 2.25h.008v.008h-.008V13.5zm0 2.25h.008v.008h-.008v-.008zm0 2.25h.008v.008h-.008V18zm2.504-6.75h.008v.008h-.008v-.008zm0 2.25h.008v.008h-.008V13.5zm0 2.25h.008v.008h-.008v-.008zm2.498-6.75h.008v.008H15.75v-.008zm0 2.25h.008v.008H15.75V13.5zM8.25 6h7.5v2.25h-7.5V6zM6 20.25h12A1.5 1.5 0 0019.5 18.75V5.25A1.5 1.5 0 0018 3.75H6A1.5 1.5 0 004.5 5.25v13.5A1.5 1.5 0 006 20.25z" /></svg>
    case 'globe':
      return <svg {...p}><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
        d="M12 21a9.004 9.004 0 008.716-6.747M12 21a9.004 9.004 0 01-8.716-6.747M12 21c2.485 0 4.5-4.03 4.5-9S14.485 3 12 3m0 18c-2.485 0-4.5-4.03-4.5-9S9.515 3 12 3m0 0a8.997 8.997 0 017.843 4.582M12 3a8.997 8.997 0 00-7.843 4.582m15.686 0A11.953 11.953 0 0112 10.5c-2.998 0-5.74-1.1-7.843-2.918m15.686 0A8.959 8.959 0 0121 12c0 .778-.099 1.533-.284 2.253m0 0A17.919 17.919 0 0112 16.5a17.92 17.92 0 01-8.716-2.247m0 0A8.966 8.966 0 013 12c0-1.777.514-3.431 1.4-4.834" /></svg>
    default:
      return <svg {...p}><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
        d="M5.25 5.653c0-.856.917-1.398 1.667-.986l11.54 6.347a1.125 1.125 0 010 1.972l-11.54 6.347a1.125 1.125 0 01-1.667-.986V5.653z" /></svg>
  }
}

// ---------------------------------------------------------------------------
// Color per run type
// ---------------------------------------------------------------------------
// Run type colors — calc subtypes derived from GROUP_META (single source of truth)
const _baseRunColors = {
  import:      { text: 'text-gray-500',  bg: 'bg-gray-100',  ring: 'ring-gray-200',  fill: 'bg-gray-500'  },
  calculation: { text: 'text-blue-700',  bg: 'bg-blue-100',  ring: 'ring-blue-200',  fill: 'bg-blue-600'  },
  generation:  { text: 'text-pink-700',  bg: 'bg-pink-100',  ring: 'ring-pink-200',  fill: 'bg-pink-600'  },
}
// Auto-generate subtype colors from GROUP_META border color (e.g. 'border-blue-400' → base 'blue')
const RUN_TYPE_COLORS = { ..._baseRunColors }
for (const [key, meta] of Object.entries(GROUP_META)) {
  if (_baseRunColors[key]) continue // skip import/calculation/generation
  const color = meta.border.replace('border-', '').replace(/-\d+$/, '') // e.g. 'blue', 'emerald'
  RUN_TYPE_COLORS[key] = {
    text: `text-${color}-700`, bg: `bg-${color}-100`,
    ring: `ring-${color}-200`, fill: `bg-${color}-600`,
  }
}

function runColor(type) {
  return RUN_TYPE_COLORS[type] || RUN_TYPE_COLORS.import
}

// ---------------------------------------------------------------------------
// Step indicator
// ---------------------------------------------------------------------------
function StepIndicator({ current, labels }) {
  const total = labels.length
  return (
    <div className="flex items-center justify-center mb-6">
      {labels.map((label, i) => {
        const step = i + 1
        const done = step < current
        const active = step === current
        return (
          <React.Fragment key={step}>
            <div className="flex flex-col items-center gap-1">
              <div className={`w-8 h-8 rounded-full flex items-center justify-center text-sm font-bold transition-all ${
                done   ? 'bg-green-500 text-white ring-2 ring-green-200'
                : active ? 'bg-bx-surface text-white ring-2 ring-bx-mint/30'
                : 'bg-gray-100 text-gray-400'
              }`}>
                {done ? (
                  <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M5 13l4 4L19 7" />
                  </svg>
                ) : step}
              </div>
              <span className={`text-[10px] font-semibold ${active ? 'text-bx-light-text' : done ? 'text-green-600' : 'text-gray-400'}`}>
                {label}
              </span>
            </div>
            {i < total - 1 && (
              <div className={`w-20 h-0.5 mx-2 mb-4 rounded-full transition-colors ${done ? 'bg-green-400' : 'bg-gray-200'}`} />
            )}
          </React.Fragment>
        )
      })}
    </div>
  )
}

// ---------------------------------------------------------------------------
// Exhaustiveness slider — conformational search depth
// ---------------------------------------------------------------------------
const EXHAUST_STEPS = [
  { value: 8,  label: 'Fast',       desc: 'Large-scale screening' },
  { value: 16, label: 'Standard',   desc: 'Balanced exploration' },
  { value: 32, label: 'Thorough',   desc: 'Best accuracy/speed trade-off' },
  { value: 64, label: 'Intensive',  desc: 'Final validation, max conformations' },
]

function ExhaustivenessSlider({ value, onChange }) {
  const idx = EXHAUST_STEPS.findIndex(s => s.value >= value)
  const activeIdx = idx === -1 ? EXHAUST_STEPS.length - 1 : idx
  const active = EXHAUST_STEPS[activeIdx]

  return (
    <div>
      <div className="flex items-center justify-between mb-1">
        <label className="text-sm font-medium text-gray-600">Search Depth</label>
        <span className="text-xs font-semibold text-bx-accent bg-bx-accent/10 px-2 py-0.5 rounded-full">
          {active.label} ({active.value})
        </span>
      </div>
      <p className="text-[10px] text-gray-400 mb-2">Conformations explored per ligand</p>
      <input
        type="range" min={0} max={EXHAUST_STEPS.length - 1} step={1}
        value={activeIdx}
        onChange={e => onChange(EXHAUST_STEPS[Number(e.target.value)].value)}
        className="w-full accent-bx-mint"
      />
      <div className="flex justify-between mt-1">
        {EXHAUST_STEPS.map((s, i) => (
          <button key={s.value} type="button"
            onClick={() => onChange(s.value)}
            className={`text-[10px] transition-colors ${i === activeIdx ? 'text-bx-accent font-semibold' : 'text-gray-400 hover:text-gray-600'}`}>
            {s.label}
          </button>
        ))}
      </div>
      <p className="text-[10px] text-gray-400 mt-1 text-center italic">{active.desc}</p>
    </div>
  )
}

// ---------------------------------------------------------------------------
// Stepper control (integer increment/decrement)
// ---------------------------------------------------------------------------
function Stepper({ value, min, max, onChange }) {
  return (
    <div className="flex items-center gap-2">
      <button
        onClick={() => onChange(Math.max(min, value - 1))}
        disabled={value <= min}
        className="w-7 h-7 rounded-lg border border-gray-200 flex items-center justify-center
                   text-gray-600 hover:bg-gray-100 disabled:opacity-30 disabled:cursor-not-allowed transition-colors"
      >
        <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M20 12H4" />
        </svg>
      </button>
      <span className="w-8 text-center font-bold text-bx-light-text text-sm tabular-nums">{value}</span>
      <button
        onClick={() => onChange(Math.min(max, value + 1))}
        disabled={value >= max}
        className="w-7 h-7 rounded-lg border border-gray-200 flex items-center justify-center
                   text-gray-600 hover:bg-gray-100 disabled:opacity-30 disabled:cursor-not-allowed transition-colors"
      >
        <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M12 4v16m8-8H4" />
        </svg>
      </button>
    </div>
  )
}

// ---------------------------------------------------------------------------
// Default configs
// ---------------------------------------------------------------------------
const DEFAULT_CONFIGS = {
  import:      {
    sourceMode: 'database', database: null, databases: [], maxPerSource: 100,
    zinc_subset: 'in-stock',
    filters: {},
  },
  calculation: {
    calculation_types: [],
    // Per-subtype configs
    docking: { engine: 'gnina_gpu', exhaustiveness: 32, num_modes: 9, seed: 0, boxSizeX: '', boxSizeY: '', boxSizeZ: '' },
    adme: {},
    toxicity: {},
    scoring: {},
    composite: { weights: { docking_score: 0.30, cnn_score: 0.20, logP: 0.15, solubility: 0.10, selectivity: 0.15, novelty: 0.10 } },
    interactions: {},
    scaffold: { cutoff: 0.65 },
    off_target: {},
    confidence: {},
    retrosynthesis: { method: 'aizynthfinder' },
  },
  afvs: {
    strategy: 'balanced', docking_scenario: 'smina_vinardo',
    exhaustiveness_prescreen: 1, exhaustiveness_primary: 8,
    reps_per_tranche: 1, tranche_pct: 0.2,
    mw_max: 500, logp_max: 5, hbd_max: 5, hba_max: 10,
    top_n_results: 10000, budget_cap_usd: 700,
  },
  generation:  { mode: 'batch', method: 'scaffold_hopping', iterations: 3, variants_per_iteration: 5, qed_min: 0.4, lipinski: true, pains_filter: true, include_docking: true, include_admet: true, include_scoring: true,
    // Optimization defaults
    opt_iterations: 5, opt_variants_per_iteration: 50,
    weights: { binding_affinity: 0.35, toxicity: 0.25, bioavailability: 0.20, synthesis_ease: 0.20 },
  },
}


// ---------------------------------------------------------------------------
// Step 2: Configuration forms
// ---------------------------------------------------------------------------
function FormSection({ title, children }) {
  return (
    <div className="space-y-2">
      <p className="text-sm font-semibold text-gray-400 uppercase tracking-wide">{title}</p>
      {children}
    </div>
  )
}

function EngineCard({ id, label, badges = [], selected, onClick, disabled = false, disabledReason = '' }) {
  return (
    <button
      onClick={() => !disabled && onClick(id)}
      disabled={disabled}
      title={disabled ? disabledReason : ''}
      className={`w-full flex items-start gap-3 p-3 rounded-xl border-2 text-left transition-all ${
        disabled
          ? 'border-gray-100 bg-gray-50 opacity-50 cursor-not-allowed'
          : selected
            ? 'border-bx-mint bg-blue-50 ring-2 ring-bx-mint/20'
            : 'border-gray-100 hover:border-blue-200 hover:shadow-sm'
      }`}
    >
      <div className={`w-2.5 h-2.5 rounded-full mt-1.5 flex-shrink-0 border-2 transition-colors ${
        disabled ? 'border-gray-200' : selected ? 'border-bx-mint bg-bx-surface' : 'border-gray-300'
      }`} />
      <div className="flex-1 min-w-0">
        <p className={`text-sm font-semibold ${disabled ? 'text-gray-400' : selected ? 'text-bx-light-text' : 'text-gray-800'}`}>{label}</p>
        <div className="flex flex-wrap gap-1 mt-1">
          {badges.map(b => (
            <span key={b.text} className={`text-[10px] px-1.5 py-0.5 rounded-full font-semibold ${b.color}`}>
              {b.text}
            </span>
          ))}
          {disabled && disabledReason && (
            <span className="text-[10px] px-1.5 py-0.5 rounded-full font-semibold bg-gray-100 text-gray-400">
              {disabledReason}
            </span>
          )}
        </div>
      </div>
    </button>
  )
}

// ---------------------------------------------------------------------------
// Range slider helper
// ---------------------------------------------------------------------------
function RangeInput({ label, min, max, step = 1, value, onChange, unit = '' }) {
  return (
    <div>
      <div className="flex justify-between text-[10px] mb-1">
        <span className="text-gray-500 font-medium">{label}</span>
        <span className="font-semibold text-gray-700 tabular-nums">{value}{unit}</span>
      </div>
      <input
        type="range" min={min} max={max} step={step} value={value}
        onChange={e => onChange(parseFloat(e.target.value))}
        className="w-full h-1.5 bg-gray-200 rounded-full appearance-none cursor-pointer accent-bx-mint"
      />
    </div>
  )
}

function RangeInputDual({ label, min, max, step = 1, valueMin, valueMax, onChangeMin, onChangeMax, unit = '' }) {
  return (
    <div>
      <div className="flex justify-between text-[10px] mb-1">
        <span className="text-gray-500 font-medium">{label}</span>
        <span className="font-semibold text-gray-700 tabular-nums">{valueMin}–{valueMax}{unit}</span>
      </div>
      <div className="flex gap-2 items-center">
        <input
          type="number" min={min} max={valueMax} step={step} value={valueMin}
          onChange={e => onChangeMin(parseFloat(e.target.value))}
          className="w-20 text-xs border border-gray-200 rounded-lg px-2 py-1.5 bg-white tabular-nums
                     focus:outline-none focus:ring-2 focus:ring-bx-mint/30 focus:border-bx-mint"
        />
        <span className="text-gray-300 text-xs">to</span>
        <input
          type="number" min={valueMin} max={max} step={step} value={valueMax}
          onChange={e => onChangeMax(parseFloat(e.target.value))}
          className="w-20 text-xs border border-gray-200 rounded-lg px-2 py-1.5 bg-white tabular-nums
                     focus:outline-none focus:ring-2 focus:ring-bx-mint/30 focus:border-bx-mint"
        />
      </div>
    </div>
  )
}

function FilterToggle({ label, checked, onChange, description }) {
  return (
    <label className={`flex items-center gap-2.5 p-2 rounded-lg cursor-pointer transition-colors ${
      checked ? 'bg-emerald-50' : 'hover:bg-gray-50'
    }`}>
      <input type="checkbox" checked={checked} onChange={e => onChange(e.target.checked)} className="accent-bx-mint" />
      <div>
        <span className="text-sm font-medium text-gray-700">{label}</span>
        {description && <p className="text-[10px] text-gray-400">{description}</p>}
      </div>
    </label>
  )
}

// ---------------------------------------------------------------------------
// Import Filters — Accordion sections
// ---------------------------------------------------------------------------
function ImportFilters({ filters, onChange, selectedDBs }) {
  const [showCommon, setShowCommon] = useState(true)
  const [showSpecific, setShowSpecific] = useState(false)
  const f = filters
  const set = (key, val) => onChange({ ...f, [key]: val })

  const hasChembl = selectedDBs.includes('chembl')
  const hasPubchem = selectedDBs.includes('pubchem')
  const hasEnamine = selectedDBs.includes('enamine')
  const hasFragments = selectedDBs.includes('fragments')
  const hasSpecificFilters = hasChembl || hasPubchem || hasEnamine || hasFragments

  return (
    <div className="space-y-2">
      {/* Common Filters Accordion */}
      <button
        onClick={() => setShowCommon(v => !v)}
        className="w-full flex items-center justify-between p-3 rounded-xl border border-gray-200 bg-gray-50 hover:bg-gray-100 transition-colors text-left"
      >
        <div className="flex items-center gap-2">
          <svg className="w-4 h-4 text-gray-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
              d="M12 3c2.755 0 5.455.232 8.083.678.533.09.917.556.917 1.096v1.044a2.25 2.25 0 01-.659 1.591l-5.432 5.432a2.25 2.25 0 00-.659 1.591v2.927a2.25 2.25 0 01-1.244 2.013L9.75 21v-6.568a2.25 2.25 0 00-.659-1.591L3.659 7.409A2.25 2.25 0 013 5.818V4.774c0-.54.384-1.006.917-1.096A48.32 48.32 0 0112 3z" />
          </svg>
          <span className="text-sm font-semibold text-gray-700">Drug-likeness Filters</span>
          <span className="text-[10px] px-1.5 py-0.5 bg-bx-mint/10 text-bx-mint rounded font-semibold">
            {(f.lipinski ? 1 : 0) + (f.pains ? 1 : 0) + 6} active
          </span>
        </div>
        <svg className={`w-4 h-4 text-gray-400 transition-transform ${showCommon ? 'rotate-180' : ''}`}
          fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
        </svg>
      </button>

      {showCommon && (
        <div className="border border-gray-200 rounded-xl p-4 space-y-3 bg-white">
          <div className="grid grid-cols-2 gap-3">
            <RangeInputDual label="MW (Da)" min={0} max={1000} step={10}
              valueMin={f.mw_min ?? 200} valueMax={f.mw_max ?? 500}
              onChangeMin={v => set('mw_min', v)} onChangeMax={v => set('mw_max', v)} unit=" Da" />
            <RangeInputDual label="LogP" min={-5} max={10} step={0.5}
              valueMin={f.logp_min ?? -1} valueMax={f.logp_max ?? 5}
              onChangeMin={v => set('logp_min', v)} onChangeMax={v => set('logp_max', v)} />
          </div>
          <div className="grid grid-cols-3 gap-3">
            <RangeInput label="HBD max" min={0} max={15} value={f.hbd_max ?? 5} onChange={v => set('hbd_max', v)} />
            <RangeInput label="HBA max" min={0} max={20} value={f.hba_max ?? 10} onChange={v => set('hba_max', v)} />
            <RangeInput label="RotBonds max" min={0} max={20} value={f.rotatable_max ?? 10} onChange={v => set('rotatable_max', v)} />
          </div>
          <div className="grid grid-cols-2 gap-3">
            <RangeInputDual label="TPSA (A²)" min={0} max={300} step={5}
              valueMin={f.tpsa_min ?? 0} valueMax={f.tpsa_max ?? 140}
              onChangeMin={v => set('tpsa_min', v)} onChangeMax={v => set('tpsa_max', v)} unit=" A²" />
            <RangeInput label="QED min" min={0} max={1} step={0.05}
              value={f.qed_min ?? 0.3} onChange={v => set('qed_min', v)} />
          </div>
          <div className="grid grid-cols-2 gap-2 pt-1">
            <FilterToggle label="Lipinski Ro5" checked={f.lipinski ?? true} onChange={v => set('lipinski', v)}
              description="Max 1 violation allowed" />
            <FilterToggle label="PAINS Filter" checked={f.pains ?? true} onChange={v => set('pains', v)}
              description="Reject reactive/false positive motifs" />
          </div>
        </div>
      )}

      {/* Database-specific Filters */}
      {hasSpecificFilters && (
        <>
          <button
            onClick={() => setShowSpecific(v => !v)}
            className="w-full flex items-center justify-between p-3 rounded-xl border border-gray-200 bg-gray-50 hover:bg-gray-100 transition-colors text-left"
          >
            <div className="flex items-center gap-2">
              <svg className="w-4 h-4 text-gray-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                  d="M10.5 6h9.75M10.5 6a1.5 1.5 0 11-3 0m3 0a1.5 1.5 0 10-3 0M3.75 6H7.5m3 12h9.75m-9.75 0a1.5 1.5 0 01-3 0m3 0a1.5 1.5 0 00-3 0m-3.75 0H7.5m9-6h3.75m-3.75 0a1.5 1.5 0 01-3 0m3 0a1.5 1.5 0 00-3 0m-9.75 0h9.75" />
              </svg>
              <span className="text-sm font-semibold text-gray-700">Database-Specific Filters</span>
            </div>
            <svg className={`w-4 h-4 text-gray-400 transition-transform ${showSpecific ? 'rotate-180' : ''}`}
              fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
            </svg>
          </button>

          {showSpecific && (
            <div className="border border-gray-200 rounded-xl p-4 space-y-4 bg-white">
              {/* ChEMBL */}
              {hasChembl && (
                <div className="space-y-2">
                  <p className="text-[10px] font-semibold text-blue-600 uppercase tracking-wider flex items-center gap-1">
                    <span>🧬</span> ChEMBL Filters
                  </p>
                  <div className="grid grid-cols-2 gap-3">
                    <div>
                      <p className="text-[10px] text-gray-500 font-medium mb-1">Activity types</p>
                      <div className="flex flex-wrap gap-1.5">
                        {['IC50', 'Ki', 'EC50', 'Kd'].map(at => {
                          const active = (f.activity_types || []).includes(at)
                          return (
                            <button key={at} onClick={() => {
                              const arr = f.activity_types || []
                              set('activity_types', active ? arr.filter(x => x !== at) : [...arr, at])
                            }}
                              className={`px-2 py-0.5 rounded text-[10px] font-semibold border transition-colors ${
                                active ? 'bg-blue-100 border-blue-300 text-blue-700' : 'bg-gray-50 border-gray-200 text-gray-500 hover:border-blue-200'
                              }`}
                            >{at}</button>
                          )
                        })}
                      </div>
                    </div>
                    <div>
                      <p className="text-[10px] text-gray-500 font-medium mb-1">Activity cutoff (nM)</p>
                      <input type="number" min={1} max={100000} step={100}
                        value={f.activity_cutoff ?? 10000}
                        onChange={e => set('activity_cutoff', parseInt(e.target.value) || 10000)}
                        className="w-28 text-xs border border-gray-200 rounded-lg px-2 py-1.5 bg-white tabular-nums
                                   focus:outline-none focus:ring-2 focus:ring-bx-mint/30 focus:border-bx-mint" />
                    </div>
                  </div>
                  <RangeInput label="pChEMBL min" min={3} max={10} step={0.1}
                    value={f.pchembl_min ?? 5.0} onChange={v => set('pchembl_min', v)} />
                </div>
              )}

              {/* PubChem */}
              {hasPubchem && (
                <div className="space-y-2">
                  <p className="text-[10px] font-semibold text-purple-600 uppercase tracking-wider flex items-center gap-1">
                    <span>🔬</span> PubChem Filters
                  </p>
                  <div className="grid grid-cols-2 gap-3">
                    <div>
                      <p className="text-[10px] text-gray-500 font-medium mb-1">Activity type</p>
                      <select
                        value={f.pubchem_activity_type || 'IC50'}
                        onChange={e => set('pubchem_activity_type', e.target.value)}
                        className="text-xs border border-gray-200 rounded-lg px-2 py-1.5 bg-white
                                   focus:outline-none focus:ring-2 focus:ring-bx-mint/30 focus:border-bx-mint"
                      >
                        {['IC50', 'Ki', 'EC50', 'Potency'].map(t => <option key={t} value={t}>{t}</option>)}
                      </select>
                    </div>
                    <FilterToggle label="Active only" checked={f.pubchem_active_only ?? true}
                      onChange={v => set('pubchem_active_only', v)} description="Only fetch 'active' compounds" />
                  </div>
                </div>
              )}

              {/* Enamine */}
              {hasEnamine && (
                <div className="space-y-2">
                  <p className="text-[10px] font-semibold text-orange-600 uppercase tracking-wider flex items-center gap-1">
                    <span>🏭</span> Enamine REAL Filters
                  </p>
                  <div>
                    <p className="text-[10px] text-gray-500 font-medium mb-1">Scaffold complexity</p>
                    <div className="flex gap-2">
                      {['2-fragment', '3-fragment', 'mixed'].map(sc => (
                        <button key={sc} onClick={() => set('scaffold_complexity', sc)}
                          className={`px-2.5 py-1 rounded-lg text-xs font-semibold border transition-colors ${
                            f.scaffold_complexity === sc ? 'bg-orange-100 border-orange-300 text-orange-700' : 'bg-gray-50 border-gray-200 text-gray-500 hover:border-orange-200'
                          }`}
                        >{sc}</button>
                      ))}
                    </div>
                  </div>
                </div>
              )}

              {/* Fragment Library */}
              {hasFragments && (
                <div className="space-y-2">
                  <p className="text-[10px] font-semibold text-teal-600 uppercase tracking-wider flex items-center gap-1">
                    <span>🧩</span> Fragment Library Filters
                  </p>
                  <FilterToggle label="Rule-of-3 strict" checked={f.ro3_strict ?? true}
                    onChange={v => set('ro3_strict', v)}
                    description="MW≤300, logP≤3, HBD≤3, HBA≤3, RotBonds≤3" />
                </div>
              )}
            </div>
          )}
        </>
      )}
    </div>
  )
}

// ---------------------------------------------------------------------------
// Step 2: Configuration forms
// ---------------------------------------------------------------------------
function ConfigForm({ runType, config, onChange, phase, selectedCount = 0, gpuEngines, project, campaign, hasDockingRun = false }) {
  const set = (key, value) => onChange({ ...config, [key]: value })

  // Auto-switch engine if GPU is unavailable and current selection is gnina_gpu
  useEffect(() => {
    if (!gpuEngines || gpuEngines.gnina_gpu) return
    // Standalone docking run
    if (runType === 'docking' && config.engine === 'gnina_gpu') {
      onChange({ ...config, engine: 'gnina_cpu' })
    }
    // Calculation sub-type docking
    if (runType === 'calculation' && config.docking?.engine === 'gnina_gpu') {
      onChange({ ...config, docking: { ...config.docking, engine: 'gnina_cpu' } })
    }
  }, [gpuEngines]) // eslint-disable-line react-hooks/exhaustive-deps

  const toggleArr = (key, item) => {
    const arr = config[key] || []
    onChange({ ...config, [key]: arr.includes(item) ? arr.filter(v => v !== item) : [...arr, item] })
  }

  const fileInputRef = useRef(null)

  switch (runType) {
    case 'import': {
      const IMPORT_SOURCES = [
        { key: 'zinc', label: 'ZINC20', desc: 'Drug-like commercially available molecules (~389M)', icon: '💊', sourceMode: 'database' },
        { key: 'chembl', label: 'ChEMBL', desc: 'Known bioactive compounds for your target', icon: '🧬', sourceMode: 'database', needsTarget: true },
        { key: 'pubchem', label: 'PubChem', desc: 'NCBI compound database with bioactivity data', icon: '🔬', sourceMode: 'database', needsTarget: true },
        { key: 'enamine', label: 'Enamine REAL', desc: 'Virtual combinatorial library (37B+ compounds)', icon: '🏭', sourceMode: 'database' },
        { key: 'fragments', label: 'Fragment Library', desc: 'Rule-of-3 fragments for FBDD screening', icon: '🧩', sourceMode: 'database' },
        { key: 'external', label: 'Upload File', desc: 'Import from SDF, SMILES, or CSV file', icon: '📁', sourceMode: 'external' },
        { key: 'internal', label: 'Phase Bookmarks', desc: 'Import bookmarked molecules from previous phase', icon: '🔖', sourceMode: 'internal' },
      ]
      const selectedSource = config.database || null

      return (
        <div className="space-y-4">
          <FormSection title="Select a source">
            <div className="space-y-2">
              {IMPORT_SOURCES.map(src => (
                <button key={src.key}
                  onClick={() => onChange({ ...config, database: src.key, sourceMode: src.sourceMode, databases: [src.key] })}
                  className={`w-full flex items-start gap-3 p-3 rounded-xl border-2 text-left transition-all ${
                    selectedSource === src.key
                      ? 'border-bx-mint bg-emerald-50'
                      : 'border-gray-100 hover:border-gray-200 hover:shadow-sm'
                  }`}
                >
                  <div className={`w-3 h-3 rounded-full mt-1 flex-shrink-0 border-2 transition-colors ${
                    selectedSource === src.key ? 'border-bx-mint bg-bx-surface' : 'border-gray-300'
                  }`} />
                  <div className="flex-1 min-w-0">
                    <div className="flex items-center gap-2">
                      <span>{src.icon}</span>
                      <span className="text-sm font-semibold text-gray-800">{src.label}</span>
                      {src.needsTarget && (
                        <span className="text-[10px] px-1.5 py-0.5 bg-blue-100 text-blue-600 rounded font-medium">target-aware</span>
                      )}
                    </div>
                    <p className="text-xs text-gray-500 mt-0.5">{src.desc}</p>
                  </div>
                </button>
              ))}
            </div>
          </FormSection>
        </div>
      )
    }

    case 'calculation': {
      const selectedCalc = (config.calculation_types || [])[0] || null
      const selectCalcType = (key) => {
        onChange({ ...config, calculation_types: [key] })
      }
      return (
        <div className="space-y-4">
          <FormSection title="Select a calculation to run">
            <div className="grid grid-cols-1 gap-2">
              {CALCULATION_SUBTYPES.map(sub => {
                const checked = selectedCalc === sub.key
                const subColor = RUN_TYPE_COLORS[sub.key] || RUN_TYPE_COLORS.calculation
                const disabled = (sub.key === 'interactions' && !hasDockingRun) || sub.key === 'off_target'
                return (
                  <button key={sub.key}
                    onClick={() => !disabled && selectCalcType(sub.key)}
                    disabled={disabled}
                    title={disabled ? 'Run docking first to enable interactions analysis' : undefined}
                    className={`flex items-start gap-3 p-3 rounded-xl border-2 text-left transition-all ${
                      disabled ? 'opacity-50 cursor-not-allowed border-gray-100 bg-gray-50'
                      : checked ? `border-blue-400 ${subColor.bg}` : 'border-gray-100 hover:border-gray-200'
                    }`}
                  >
                    <div className={`w-3 h-3 rounded-full mt-1 flex-shrink-0 border-2 transition-colors ${
                      checked && !disabled ? 'border-blue-500 bg-blue-500' : 'border-gray-300'
                    }`} />
                    <div className="flex-1 min-w-0">
                      <div className="flex items-center gap-2">
                        <p className={`text-sm font-semibold ${disabled ? 'text-gray-400' : checked ? 'text-bx-light-text' : 'text-gray-700'}`}>
                          {sub.label}
                          {TIPS[`run_${sub.key}`] && <InfoTip text={TIPS[`run_${sub.key}`]} size="xs" />}
                        </p>
                        <span className="text-[9px] text-gray-400">{ESTIMATED_TIMES[sub.key]}</span>
                        {disabled && (
                          <span className="text-[9px] px-1.5 py-0.5 rounded-full bg-red-100 text-red-600 font-semibold">
                            {sub.key === 'off_target' ? 'Not available' : 'Requires docking'}
                          </span>
                        )}
                      </div>
                      <p className="text-sm text-gray-500 mt-0.5">{sub.description}</p>
                      <div className="flex flex-wrap gap-1 mt-1">
                        {sub.columns.slice(0, 4).map(c => (
                          <span key={c} className="text-[9px] px-1.5 py-0.5 rounded bg-white/70 text-gray-500 border border-gray-200">{c}</span>
                        ))}
                        {sub.columns.length > 4 && <span className="text-[9px] text-gray-400">+{sub.columns.length - 4}</span>}
                      </div>
                    </div>
                  </button>
                )
              })}
            </div>
          </FormSection>

          {/* Per-subtype config panels */}
          {/* Docking + scoring config panels moved to Step 3 (Review) */}

          {selectedCalc === 'docking' && !project?.receptor_prep_report && (
            <div className="bg-red-50 border border-red-200 rounded-xl px-4 py-3 text-sm text-red-700 flex items-start gap-2">
              <svg className="w-5 h-5 flex-shrink-0 mt-0.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
              </svg>
              <span>
                <strong>Receptor not prepared.</strong> Go to Target Setup and click "Prepare Receptor" before running docking.
              </span>
            </div>
          )}

          {selectedCalc === 'pharmacophore' && !campaign?.reference_ligands?.length && (
            <div className="bg-red-50 border border-red-200 rounded-xl px-4 py-3 text-sm text-red-700 flex items-start gap-2">
              <svg className="w-5 h-5 flex-shrink-0 mt-0.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
              </svg>
              <span>
                <strong>No reference ligands set.</strong> Go to Project Home → Campaign to add at least one reference ligand.
              </span>
            </div>
          )}

          {/* Interactions info cards */}
          {selectedCalc === 'interactions' && (
            <div className="space-y-3">
              <div className="p-3 rounded-xl border border-lime-200 bg-lime-50">
                <p className="text-sm font-semibold text-lime-700">Protein-Ligand Interaction Analysis</p>
                <p className="text-sm text-lime-600 mt-1">Computes H-bonds, hydrophobic contacts, pi-stacking and salt bridges using ProLIF/RDKit.</p>
              </div>
            </div>
          )}

          {/* Scaffold / Structural Analysis info */}
          {selectedCalc === 'scaffold' && (
            <div className="p-3 rounded-xl border border-emerald-200 bg-emerald-50">
              <p className="text-sm font-semibold text-emerald-700">Structural Analysis</p>
              <p className="text-sm text-emerald-600 mt-1">Murcko scaffold extraction, BRICS bond analysis, R-group positions, and Butina diversity clustering on Morgan fingerprints. Clustering cutoff is configurable in the next step.</p>
            </div>
          )}

          {!selectedCalc && (
            <div className="bg-amber-50 border border-amber-100 rounded-xl px-4 py-3 text-sm text-amber-700">
              Select a calculation type to proceed.
            </div>
          )}
        </div>
      )
    }

    case 'docking':
      return (
        <div className="space-y-5">
          <FormSection title="Docking Engine">
            <div className="space-y-2">
              <EngineCard id="gnina_gpu" label="GNINA (GPU)"
                badges={[{ text: 'Fastest', color: 'bg-green-100 text-green-700' }, { text: 'CNN rescoring', color: 'bg-blue-100 text-blue-700' }]}
                selected={config.engine === 'gnina_gpu'} onClick={v => set('engine', v)}
                disabled={gpuEngines && !gpuEngines.gnina_gpu}
                disabledReason="RunPod not configured" />
              <EngineCard id="gnina_cpu" label="GNINA (CPU)"
                badges={[{ text: 'CNN rescoring', color: 'bg-blue-100 text-blue-700' }]}
                selected={config.engine === 'gnina_cpu'} onClick={v => set('engine', v)} />
              <EngineCard id="vina" label="AutoDock Vina"
                badges={[{ text: 'Classic', color: 'bg-gray-100 text-gray-600' }]}
                selected={config.engine === 'vina'} onClick={v => set('engine', v)} />
              <EngineCard id="diffdock" label="DiffDock"
                badges={[{ text: 'Deep Learning', color: 'bg-purple-100 text-purple-700' }, { text: 'Blind docking', color: 'bg-amber-100 text-amber-700' }]}
                selected={config.engine === 'diffdock'} onClick={v => set('engine', v)}
                disabled={true}
                disabledReason="Coming soon" />
            </div>
          </FormSection>

          <FormSection title="Search Parameters">
            <div className="space-y-3">
              <ExhaustivenessSlider
                value={config.exhaustiveness ?? 32}
                onChange={v => set('exhaustiveness', v)}
              />

              <div className="flex items-center justify-between">
                <div>
                  <label className="text-sm font-medium text-gray-600">Num poses</label>
                  <p className="text-[10px] text-gray-400">per molecule</p>
                </div>
                <Stepper
                  value={config.num_modes ?? 9}
                  min={1} max={20}
                  onChange={v => set('num_modes', v)}
                />
              </div>

              <div className="flex items-center justify-between">
                <div>
                  <label className="text-sm font-medium text-gray-600">Seed</label>
                  <p className="text-[10px] text-gray-400">0 = random</p>
                </div>
                <input
                  type="number" min={0}
                  value={config.seed ?? 0}
                  onChange={e => set('seed', Number(e.target.value))}
                  className="w-20 text-sm text-right border border-gray-200 rounded-lg px-2 py-1.5
                             focus:outline-none focus:ring-2 focus:ring-bx-mint/30 focus:border-bx-mint tabular-nums"
                />
              </div>

              <div>
                <div className="flex items-center justify-between text-sm mb-1">
                  <label className="font-medium text-gray-600">Energy range</label>
                  <span className="font-bold text-bx-light-text tabular-nums">{config.energy_range ?? 3} kcal/mol</span>
                </div>
                <input
                  type="range" min={1} max={5} step={1}
                  value={config.energy_range ?? 3}
                  onChange={e => set('energy_range', Number(e.target.value))}
                  className="w-full accent-bx-mint"
                />
                <div className="flex justify-between text-[10px] text-gray-400 mt-0.5">
                  {[1, 2, 3, 4, 5].map(v => <span key={v}>{v}</span>)}
                </div>
              </div>

              <div>
                <div className="flex items-center justify-between text-sm mb-1">
                  <label className="font-medium text-gray-600">Box padding</label>
                  <span className="font-bold text-bx-light-text tabular-nums">{config.autobox_add ?? 4} Å</span>
                </div>
                <input
                  type="range" min={2} max={12} step={2}
                  value={config.autobox_add ?? 4}
                  onChange={e => set('autobox_add', Number(e.target.value))}
                  className="w-full accent-bx-mint"
                />
                <div className="flex justify-between text-[10px] text-gray-400 mt-0.5">
                  {[2, 4, 6, 8, 10, 12].map(v => <span key={v}>{v}</span>)}
                </div>
                <p className="text-[10px] text-gray-400 mt-0.5">Padding around pocket. 4Å for known sites, 8-12Å if uncertain.</p>
              </div>

              {(config.engine === 'gnina_gpu' || config.engine === 'gnina_cpu' || !config.engine) && (
                <div>
                  <label className="text-sm font-medium text-gray-600 mb-1 block">CNN scoring mode</label>
                  <div className="flex gap-2">
                    {[
                      { id: 'rescore', label: 'Rescore', desc: 'Fast — same speed as Vina + CNN rescore' },
                      { id: 'refinement', label: 'Refinement', desc: 'Better poses — CNN guides optimization' },
                    ].map(mode => (
                      <button key={mode.id}
                        onClick={() => set('cnn_scoring', mode.id)}
                        className={`flex-1 text-left p-2 rounded-lg border text-xs transition-all ${
                          (config.cnn_scoring || 'rescore') === mode.id
                            ? 'border-blue-300 bg-blue-50 text-blue-700'
                            : 'border-gray-200 text-gray-500 hover:border-gray-300'
                        }`}>
                        <div className="font-medium">{mode.label}</div>
                        <div className="text-[10px] opacity-70">{mode.desc}</div>
                      </button>
                    ))}
                  </div>
                  {config.cnn_scoring === 'refinement' && (
                    <div className="flex items-center gap-1.5 mt-2 text-[10px] text-amber-600 bg-amber-50 rounded-lg px-2.5 py-1.5 border border-amber-200">
                      <svg className="w-3.5 h-3.5 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
                      </svg>
                      ~3x slower, may reject ~17% of molecules (positive Vina scores)
                    </div>
                  )}
                </div>
              )}
            </div>
          </FormSection>

          <FormSection title="Box Size (Angstroms)">
            <div className="grid grid-cols-3 gap-2">
              {['X', 'Y', 'Z'].map(axis => (
                <div key={axis}>
                  <label className="block text-[10px] text-gray-400 mb-0.5">{axis}</label>
                  <input
                    type="number" min={1}
                    value={config[`boxSize${axis}`] || ''}
                    onChange={e => set(`boxSize${axis}`, e.target.value)}
                    placeholder="Auto"
                    className="w-full text-sm text-center border border-gray-200 rounded-lg px-2 py-1.5
                               focus:outline-none focus:ring-2 focus:ring-bx-mint/30 focus:border-bx-mint"
                  />
                </div>
              ))}
            </div>
            <p className="text-[10px] text-gray-400">Leave blank to use pocket-derived box</p>
          </FormSection>
        </div>
      )

    case 'adme':
    case 'toxicity':
      return null  // No extra config — column selection handled in Step 3 checklist

    case 'scoring':
      return null  // No extra config — physicochemical descriptors only

    case 'composite': {
      const weights = config.weights || DEFAULT_CONFIGS.composite.weights
      const total = Object.values(weights).reduce((a, b) => a + b, 0)
      const totalOk = Math.abs(total - 1.0) < 0.01
      return (
        <div className="space-y-4">
          <FormSection title="Score Weights">
            <p className="text-xs text-gray-400 mb-2">
              Adjust the importance of each metric. The composite score will aggregate results from your previous runs.
            </p>
            {Object.entries(weights).map(([key, val]) => (
              <div key={key}>
                <div className="flex items-center justify-between text-sm mb-1">
                  <span className="font-medium text-gray-600 capitalize">{key.replace(/_/g, ' ')}</span>
                  <span className="font-bold text-bx-light-text tabular-nums">{Number(val).toFixed(2)}</span>
                </div>
                <input
                  type="range" min={0} max={1} step={0.05}
                  value={val}
                  onChange={e => onChange({
                    ...config,
                    weights: { ...weights, [key]: parseFloat(e.target.value) }
                  })}
                  className="w-full accent-bx-mint"
                />
              </div>
            ))}
          </FormSection>

          <div className="flex items-center justify-between">
            <div className={`flex items-center gap-2 text-sm font-semibold ${totalOk ? 'text-green-600' : 'text-red-600'}`}>
              <span className={`w-2 h-2 rounded-full flex-shrink-0 ${totalOk ? 'bg-green-500' : 'bg-red-500'}`} />
              Total weight: {total.toFixed(2)}
              {totalOk ? ' — OK' : ' — must sum to 1.0'}
            </div>
            <button
              onClick={() => onChange({ ...config, weights: DEFAULT_CONFIGS.composite.weights })}
              className="text-sm px-3 py-1.5 rounded-lg bg-gray-100 text-gray-600 hover:bg-gray-200 transition-colors font-medium"
            >
              Reset defaults
            </button>
          </div>
        </div>
      )
    }

    case 'scaffold': {
      const cutoff = config.cutoff ?? 0.65
      const similarity = Math.round((1 - cutoff) * 100)
      return (
        <div className="space-y-5">
          <div className="p-3 rounded-xl border border-emerald-200 bg-emerald-50">
            <p className="text-sm font-semibold text-emerald-700">Structural Analysis</p>
            <p className="text-sm text-emerald-600 mt-1">Murcko scaffold extraction, BRICS bond analysis, R-group position identification, and Butina diversity clustering on Morgan fingerprints. Only requires SMILES — no docking needed.</p>
          </div>

          <FormSection title="Clustering Cutoff">
            <div>
              <div className="flex items-center justify-between text-sm mb-1">
                <span className="font-medium text-gray-600">Cutoff</span>
                <span className="font-bold text-bx-light-text tabular-nums">{Number(cutoff).toFixed(2)}</span>
              </div>
              <input
                type="range" min={0.3} max={0.8} step={0.05}
                value={cutoff}
                onChange={e => set('cutoff', parseFloat(e.target.value))}
                className="w-full accent-bx-mint"
              />
              <div className="flex justify-between text-[10px] text-gray-400 mt-0.5">
                <span>0.30 (strict)</span>
                <span>0.80 (loose)</span>
              </div>
              <div className="mt-2 bg-emerald-50 border border-emerald-100 rounded-lg px-3 py-2 text-sm text-emerald-700">
                Molecules with <span className="font-bold">&ge;{similarity}%</span> Tanimoto similarity will be grouped together
              </div>
            </div>
          </FormSection>
        </div>
      )
    }

    case 'generation': {
      const methods = [
        { value: 'scaffold_hopping',  label: 'Scaffold Hopping',  desc: 'Replace core scaffold while keeping substituents', implemented: true },
        { value: 'fragment_growing',  label: 'Fragment Growing',  desc: 'Extend molecule at reactive positions (BRICS fragments)', implemented: true },
        { value: 'de_novo',           label: 'De Novo SMILES',    desc: 'Generate novel structures from scratch (unconstrained)', implemented: true },
        { value: 'bioisosteric',      label: 'Bioisosteric',      desc: 'Replace functional groups with bioisosteric equivalents', implemented: true },
        { value: 'fragment_linking',  label: 'Fragment Linking',   desc: 'Recombine BRICS fragments from multiple parents', implemented: true },
        { value: 'matched_pairs',     label: 'Matched Pairs',     desc: 'Apply curated medicinal chemistry transformations', implemented: true },
      ]
      const mode = config.mode || 'batch'
      const isOptimization = mode === 'optimization'
      const iters = isOptimization ? (config.opt_iterations ?? 5) : (config.iterations ?? 3)
      const vars = isOptimization ? (config.opt_variants_per_iteration ?? 50) : (config.variants_per_iteration ?? 5)
      const estTotal = isOptimization ? `~${vars * iters} tested → top 10` : `~${iters * vars * (selectedCount || 1)} molecules`

      // Optimization weights
      const optWeights = config.weights || { binding_affinity: 0.35, toxicity: 0.25, bioavailability: 0.20, synthesis_ease: 0.20 }
      const optWeightTotal = Object.values(optWeights).reduce((a, b) => a + b, 0)
      const optWeightOk = Math.abs(optWeightTotal - 1.0) < 0.02

      const WEIGHT_META = [
        { key: 'binding_affinity', label: 'Binding Affinity', color: 'accent-blue-500', desc: 'Docking score priority' },
        { key: 'toxicity', label: 'Low Toxicity', color: 'accent-green-500', desc: 'hERG, hepatotoxicity, AMES' },
        { key: 'bioavailability', label: 'Bioavailability', color: 'accent-purple-500', desc: 'Oral absorption, BBB' },
        { key: 'synthesis_ease', label: 'Synthesis Ease', color: 'accent-orange-500', desc: 'Synthetic accessibility' },
      ]

      return (
        <div className="space-y-5">
          {/* Mode toggle: Batch vs Optimization vs Molecule */}
          <FormSection title="Generation mode">
            <div className="flex gap-2">
              {[
                { value: 'batch', label: 'Batch', desc: 'Apply to all selected molecules', implemented: true },
                { value: 'optimization', label: 'Optimization', desc: 'Multi-iteration lead optimization', implemented: true },
                { value: 'molecule', label: 'Per-Molecule', desc: 'R-group control per molecule', implemented: true },
              ].map(m => (
                <button
                  key={m.value}
                  onClick={() => m.implemented ? set('mode', m.value) : null}
                  disabled={!m.implemented}
                  className={`flex-1 p-3 rounded-xl border-2 text-left transition-all ${
                    !m.implemented
                      ? 'border-gray-100 bg-gray-50/50 opacity-60 cursor-not-allowed'
                      : mode === m.value
                        ? 'border-pink-400 bg-pink-50'
                        : 'border-gray-100 hover:border-pink-200'
                  }`}
                >
                  <p className={`text-sm font-semibold flex items-center gap-2 ${mode === m.value && m.implemented ? 'text-pink-700' : 'text-gray-700'}`}>
                    {m.label}
                    {!m.implemented && (
                      <span className="text-[9px] font-bold uppercase px-1.5 py-0.5 rounded bg-amber-100 text-amber-600 border border-amber-200">
                        Coming Soon
                      </span>
                    )}
                  </p>
                  <p className="text-[10px] text-gray-400 mt-0.5">{m.desc}</p>
                </button>
              ))}
            </div>
          </FormSection>

          {/* Optimization mode — weights + parameters */}
          {isOptimization && (
            <>
              <FormSection title="Optimization Objectives">
                <div className="space-y-3 bg-gradient-to-b from-pink-50/50 to-white rounded-xl border border-pink-100 p-4">
                  <p className="text-[10px] text-gray-400 mb-2">
                    Set the relative importance of each objective. The optimizer will balance these goals across iterations.
                  </p>
                  {WEIGHT_META.map(wm => (
                    <div key={wm.key}>
                      <div className="flex items-center justify-between text-sm mb-1">
                        <div>
                          <span className="font-medium text-gray-700">{wm.label}</span>
                          <span className="text-[10px] text-gray-400 ml-2">{wm.desc}</span>
                        </div>
                        <span className="font-bold text-pink-700 tabular-nums w-10 text-right">
                          {(optWeights[wm.key] ?? 0).toFixed(2)}
                        </span>
                      </div>
                      <input
                        type="range" min={0} max={1} step={0.05}
                        value={optWeights[wm.key] ?? 0}
                        onChange={e => set('weights', { ...optWeights, [wm.key]: parseFloat(e.target.value) })}
                        className={`w-full h-1.5 rounded-full appearance-none cursor-pointer ${wm.color}`}
                      />
                    </div>
                  ))}
                  <div className="flex items-center justify-between pt-2 border-t border-pink-100">
                    <div className={`flex items-center gap-2 text-sm font-semibold ${optWeightOk ? 'text-green-600' : 'text-red-500'}`}>
                      <span className={`w-2 h-2 rounded-full ${optWeightOk ? 'bg-green-500' : 'bg-red-500'}`} />
                      Total: {optWeightTotal.toFixed(2)}
                      {optWeightOk ? '' : ' — adjust to sum to 1.0'}
                    </div>
                    <button
                      onClick={() => set('weights', { binding_affinity: 0.35, toxicity: 0.25, bioavailability: 0.20, synthesis_ease: 0.20 })}
                      className="text-xs px-2.5 py-1 rounded-lg bg-gray-100 text-gray-500 hover:bg-gray-200 transition-colors font-medium"
                    >
                      Reset
                    </button>
                  </div>
                </div>
              </FormSection>

              <FormSection title="Optimization Parameters">
                <div className="space-y-4">
                  <div className="flex items-center justify-between">
                    <div>
                      <p className="text-sm font-medium text-gray-600">Iterations</p>
                      <p className="text-[10px] text-gray-400">optimization rounds (more = better but slower)</p>
                    </div>
                    <Stepper value={config.opt_iterations ?? 5} min={2} max={20}
                      onChange={v => set('opt_iterations', v)} />
                  </div>
                  <div className="flex items-center justify-between">
                    <div>
                      <p className="text-sm font-medium text-gray-600">Variants per iteration</p>
                      <p className="text-[10px] text-gray-400">structural modifications per round</p>
                    </div>
                    <Stepper value={config.opt_variants_per_iteration ?? 50} min={10} max={100}
                      onChange={v => set('opt_variants_per_iteration', v)} />
                  </div>
                </div>
              </FormSection>

              <div className="bg-pink-50 border border-pink-100 rounded-xl p-3 space-y-1">
                <div className="flex items-center justify-between">
                  <span className="text-sm text-pink-600">Estimated scope</span>
                  <span className="text-sm font-bold text-pink-700 tabular-nums">{estTotal}</span>
                </div>
                <p className="text-[10px] text-gray-400">
                  Generates variants, docks, scores with multi-objective function, keeps top 5, iterates. Best molecules stored.
                </p>
              </div>
            </>
          )}

          {/* Batch mode: Method selection */}
          {!isOptimization && (
            <FormSection title="Method">
              <div className="space-y-2">
                {methods.map(m => (
                  <button
                    key={m.value}
                    onClick={() => m.implemented ? set('method', m.value) : null}
                    disabled={!m.implemented}
                    className={`w-full flex items-start gap-3 p-3 rounded-xl border-2 text-left transition-all ${
                      !m.implemented
                        ? 'border-gray-100 bg-gray-50/50 opacity-60 cursor-not-allowed'
                        : config.method === m.value
                          ? 'border-pink-400 bg-pink-50'
                          : 'border-gray-100 hover:border-pink-200 hover:bg-pink-50/20'
                    }`}
                  >
                    <div className={`w-2.5 h-2.5 rounded-full mt-1.5 border-2 flex-shrink-0 ${
                      config.method === m.value && m.implemented ? 'border-pink-500 bg-pink-500' : 'border-gray-300'
                    }`} />
                    <div className="flex-1">
                      <p className={`text-sm font-semibold flex items-center gap-2 ${config.method === m.value && m.implemented ? 'text-pink-700' : 'text-gray-800'}`}>
                        {m.label}
                        {!m.implemented && (
                          <span className="text-[9px] font-bold uppercase px-1.5 py-0.5 rounded bg-amber-100 text-amber-600 border border-amber-200">
                            Coming Soon
                          </span>
                        )}
                      </p>
                      <p className="text-sm text-gray-400 mt-0.5">{m.desc}</p>
                    </div>
                  </button>
                ))}
              </div>
            </FormSection>
          )}

          {/* Molecule mode: ScaffoldAnalyzer */}
          {mode === 'molecule' && config.molecule_smiles && (
            <FormSection title="R-group control">
              <Suspense fallback={<div className="h-32 flex items-center justify-center"><BindXLogo variant="loading" size={24} label="Loading analyzer..." /></div>}>
                <ScaffoldAnalyzer
                  smiles={config.molecule_smiles}
                  onRulesChange={(rules) => set('scaffold_rules', rules)}
                />
              </Suspense>
            </FormSection>
          )}

          {mode === 'molecule' && !config.molecule_smiles && (
            <div className="bg-pink-50 border border-pink-100 rounded-xl p-4 text-sm text-pink-600">
              Select a single molecule in the dashboard to enable per-molecule R-group control. The molecule's SMILES will be analyzed for modifiable positions.
            </div>
          )}

          {/* Batch mode: Parameters */}
          {!isOptimization && (
            <>
              <FormSection title="Generation Parameters">
                <div className="space-y-4">
                  <div className="flex items-center justify-between">
                    <div>
                      <p className="text-sm font-medium text-gray-600">Iterations</p>
                      <p className="text-[10px] text-gray-400">rounds of generation</p>
                    </div>
                    <Stepper value={iters} min={1} max={5}
                      onChange={v => set('iterations', v)} />
                  </div>
                  <div className="flex items-center justify-between">
                    <div>
                      <p className="text-sm font-medium text-gray-600">Variants per iteration</p>
                      <p className="text-[10px] text-gray-400">new molecules per round</p>
                    </div>
                    <Stepper value={vars} min={3} max={20}
                      onChange={v => set('variants_per_iteration', v)} />
                  </div>
                </div>
              </FormSection>

              {/* Filters */}
              <FormSection title="Quality filters">
                <div className="space-y-3">
                  <div className="flex items-center justify-between">
                    <div>
                      <p className="text-sm font-medium text-gray-600">Min QED</p>
                      <p className="text-[10px] text-gray-400">drug-likeness threshold (0-1)</p>
                    </div>
                    <div className="flex items-center gap-2">
                      <input
                        type="range" min="0" max="0.9" step="0.05"
                        value={config.qed_min ?? 0.4}
                        onChange={e => set('qed_min', parseFloat(e.target.value))}
                        className="w-24 h-1.5 accent-pink-500"
                      />
                      <span className="text-sm font-mono font-semibold text-gray-700 w-8 text-right">
                        {(config.qed_min ?? 0.4).toFixed(2)}
                      </span>
                    </div>
                  </div>
                  <label className="flex items-center justify-between cursor-pointer">
                    <div>
                      <p className="text-sm font-medium text-gray-600">Lipinski Rule of 5</p>
                      <p className="text-[10px] text-gray-400">enforce drug-like properties</p>
                    </div>
                    <input type="checkbox" checked={config.lipinski !== false}
                      onChange={() => set('lipinski', !config.lipinski)}
                      className="accent-pink-500 w-4 h-4" />
                  </label>
                  <label className="flex items-center justify-between cursor-pointer">
                    <div>
                      <p className="text-sm font-medium text-gray-600">PAINS filter</p>
                      <p className="text-[10px] text-gray-400">reject pan-assay interference compounds</p>
                    </div>
                    <input type="checkbox" checked={config.pains_filter !== false}
                      onChange={() => set('pains_filter', !config.pains_filter)}
                      className="accent-pink-500 w-4 h-4" />
                  </label>
                </div>
              </FormSection>

              {/* Estimate preview */}
              <div className="bg-pink-50 border border-pink-100 rounded-xl p-3 flex items-center justify-between">
                <span className="text-sm text-pink-600">Estimated output</span>
                <span className="text-sm font-bold text-pink-700 tabular-nums">{estTotal}</span>
              </div>

              <FormSection title="Post-generation analyses">
                <div className="flex flex-wrap gap-2">
                  {[
                    { key: 'include_docking', label: 'Docking' },
                    { key: 'include_admet',   label: 'ADMET' },
                    { key: 'include_scoring', label: 'Scoring' },
                  ].map(item => {
                    const checked = config[item.key] !== false
                    return (
                      <label key={item.key}
                        className={`flex items-center gap-1.5 px-3 py-1.5 rounded-full border text-sm font-semibold cursor-pointer transition-all ${
                          checked
                            ? 'border-pink-300 bg-pink-50 text-pink-700'
                            : 'border-gray-200 text-gray-400 hover:border-pink-200'
                        }`}
                      >
                        <input
                          type="checkbox"
                          checked={checked}
                          onChange={() => set(item.key, !checked)}
                          className="accent-pink-500"
                        />
                        {item.label}
                      </label>
                    )
                  })}
                </div>
                <p className="text-[10px] text-gray-400 mt-1">Auto-run docking/ADMET/scoring on generated molecules</p>
              </FormSection>
            </>
          )}
        </div>
      )
    }

    case 'afvs':
      return (
        <AFVSConfigForm
          config={config}
          onChange={onChange}
          project={project}
          campaign={campaign}
        />
      )

    default:
      return <p className="text-sm text-gray-400">No configuration needed for this run type.</p>
  }
}

// ---------------------------------------------------------------------------
// Column checklist for calculation runs (Step 3)
// ---------------------------------------------------------------------------
function ColumnChecklist({ calcGroupKey, includedColumns, onChange, config = {} }) {
  const sub = CALCULATION_SUBTYPES.find(s => s.key === calcGroupKey)
  if (!sub || !sub.columns.length) return null

  // Filter columns based on current config (e.g. engine selection)
  const availableKeys = getAvailableColumns(calcGroupKey, config)
  const colDefs = availableKeys.map(k => ALL_COLUMNS.find(c => c.key === k)).filter(Boolean)

  // Auto-clean includedColumns when config changes (e.g. engine switch removes CNN cols)
  const availableKeyStr = availableKeys.join(',')
  useEffect(() => {
    const avail = new Set(availableKeys)
    const cleaned = includedColumns.filter(k => avail.has(k) || !sub.columns.includes(k))
    if (cleaned.length !== includedColumns.length) {
      onChange(cleaned)
    }
  }, [availableKeyStr]) // eslint-disable-line react-hooks/exhaustive-deps

  const total = colDefs.length
  const checkedCount = colDefs.filter(c => includedColumns.includes(c.key)).length
  const allChecked = checkedCount === total

  function toggle(key) {
    const def = COLUMN_DEFAULTS[key]
    if (def?.requires) return // disabled
    onChange(
      includedColumns.includes(key)
        ? includedColumns.filter(k => k !== key)
        : [...includedColumns, key]
    )
  }

  function toggleAll() {
    const enabledKeys = colDefs.filter(c => !COLUMN_DEFAULTS[c.key]?.requires).map(c => c.key)
    if (allChecked) {
      onChange(includedColumns.filter(k => !enabledKeys.includes(k)))
    } else {
      const newSet = new Set([...includedColumns, ...enabledKeys])
      onChange([...newSet])
    }
  }

  return (
    <div className="bg-gray-50 rounded-xl border border-gray-100 p-4 space-y-3">
      <div className="flex items-center justify-between">
        <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wide">Columns to compute</p>
        <button
          onClick={toggleAll}
          className="text-[11px] font-medium text-bx-accent hover:underline"
        >
          {allChecked ? 'Deselect all' : 'Select all'}
        </button>
      </div>

      <div className="grid grid-cols-2 gap-x-4 gap-y-1.5">
        {colDefs.map(col => {
          const def = COLUMN_DEFAULTS[col.key]
          const disabled = !!def?.requires
          const checked = includedColumns.includes(col.key)
          return (
            <label
              key={col.key}
              className={`flex items-center gap-2 py-1 px-1.5 rounded-lg text-sm cursor-pointer select-none transition-colors
                ${disabled ? 'opacity-50 cursor-not-allowed' : 'hover:bg-gray-100'}`}
              title={disabled ? def.requires : ''}
            >
              <input
                type="checkbox"
                checked={checked}
                disabled={disabled}
                onChange={() => toggle(col.key)}
                className="rounded border-gray-300 text-bx-accent focus:ring-bx-accent/30 w-3.5 h-3.5"
              />
              <span className={`${disabled ? 'text-gray-400' : 'text-gray-700'}`}>{col.label}</span>
              {disabled && (
                <span className="text-[10px] text-amber-500 flex items-center gap-0.5">
                  <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
                  </svg>
                  {def.requires}
                </span>
              )}
            </label>
          )
        })}
      </div>

      <p className="text-[11px] text-gray-400 text-right">{checkedCount} / {total} columns selected</p>
    </div>
  )
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ZINC20 subsets (synced with backend ZINC20_SUBSETS)
// ---------------------------------------------------------------------------
const ZINC20_SUBSETS = [
  { key: 'in-stock', label: 'In Stock', desc: '~12M', count: 12000000 },
  { key: 'for-sale', label: 'For Sale', desc: '~389M', count: 389000000 },
  { key: 'fda', label: 'FDA Approved', desc: '~1.4K', count: 1400 },
  { key: 'in-trials', label: 'In Trials', desc: '~5.8K', count: 5800 },
  { key: 'in-vivo', label: 'In Vivo', desc: '~115K', count: 115000 },
  { key: 'in-vitro', label: 'In Vitro Active', desc: '~276K', count: 276000 },
  { key: 'natural-products', label: 'Natural Products', desc: '~80K', count: 80000 },
  { key: 'biogenic', label: 'Biogenic', desc: '~135K', count: 135000 },
  { key: 'clean', label: 'PAINS-Filtered', desc: '~125M', count: 125000000 },
  { key: 'endogenous', label: 'Metabolites', desc: '~50K', count: 50000 },
  { key: 'world', label: 'World Approved', desc: '~3.4K', count: 3400 },
]

// ---------------------------------------------------------------------------
// Import source-specific configuration (displayed in Step 3 / Review)
// ---------------------------------------------------------------------------
function ImportSourceConfig({ config, onConfigChange }) {
  const db = config.database
  const f = config.filters || {}
  const setFilter = (key, val) => onConfigChange({ ...config, filters: { ...f, [key]: val } })
  const set = (key, val) => onConfigChange({ ...config, [key]: val })

  if (db === 'zinc') {
    // Multi-subset: zinc_subset is a "+" joined string (e.g. "in-stock+fda")
    const selectedSubsets = (config.zinc_subset || 'in-stock').split('+').filter(Boolean)
    const toggleSubset = (key) => {
      let next
      if (selectedSubsets.includes(key)) {
        next = selectedSubsets.filter(k => k !== key)
      } else {
        next = [...selectedSubsets, key]
      }
      if (next.length === 0) next = ['in-stock'] // at least one required
      set('zinc_subset', next.join('+'))
    }
    // Approximate total count (union — may overlap, shown as "up to")
    const approxCount = selectedSubsets.reduce((sum, k) => {
      const s = ZINC20_SUBSETS.find(x => x.key === k)
      return sum + (s?.count || 0)
    }, 0)
    const formatCount = (n) => {
      if (n >= 1e9) return `${(n / 1e9).toFixed(1)}B`
      if (n >= 1e6) return `${(n / 1e6).toFixed(1)}M`
      if (n >= 1e3) return `${(n / 1e3).toFixed(1)}K`
      return String(n)
    }

    return (
      <div className="bg-gray-50 rounded-xl border border-gray-100 p-4 space-y-4">
        <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wide">ZINC20 Configuration</p>

        {/* Subsets — multi-select checkboxes */}
        <div>
          <div className="flex items-center justify-between mb-2">
            <p className="text-xs font-medium text-gray-600">Subsets</p>
            <p className="text-[10px] tabular-nums text-gray-500">
              ~<span className="font-semibold text-gray-700">{formatCount(approxCount)}</span> compounds available
            </p>
          </div>
          <div className="grid grid-cols-2 gap-1.5">
            {ZINC20_SUBSETS.map(s => {
              const checked = selectedSubsets.includes(s.key)
              return (
                <label key={s.key}
                  className={`flex items-center gap-2 px-2.5 py-1.5 rounded-lg border cursor-pointer transition-colors text-xs ${
                    checked
                      ? 'bg-blue-50 border-blue-200 text-blue-800'
                      : 'bg-white border-gray-150 text-gray-600 hover:border-gray-300'
                  }`}
                >
                  <input type="checkbox" checked={checked} onChange={() => toggleSubset(s.key)}
                    className="w-3.5 h-3.5 rounded border-gray-300 text-blue-600 focus:ring-blue-500/30" />
                  <span className="font-medium">{s.label}</span>
                  <span className="text-[10px] text-gray-400 ml-auto">{s.desc}</span>
                </label>
              )
            })}
          </div>
        </div>

        {/* Property filters */}
        <div>
          <p className="text-xs font-medium text-gray-600 mb-2">Property Filters</p>
          <div className="grid grid-cols-2 gap-3">
            <RangeInputDual label="MW (Da)" min={0} max={1000} step={10}
              valueMin={f.mwt_min ?? 200} valueMax={f.mwt_max ?? 500}
              onChangeMin={v => setFilter('mwt_min', v)} onChangeMax={v => setFilter('mwt_max', v)} unit=" Da" />
            <RangeInputDual label="LogP" min={-5} max={10} step={0.5}
              valueMin={f.logp_min ?? -1} valueMax={f.logp_max ?? 5}
              onChangeMin={v => setFilter('logp_min', v)} onChangeMax={v => setFilter('logp_max', v)} />
          </div>
          <div className="grid grid-cols-2 gap-3 mt-2">
            <RangeInput label="HBD max" min={0} max={15} value={f.hbd_max ?? 5} onChange={v => setFilter('hbd_max', v)} />
            <RangeInput label="HBA max" min={0} max={20} value={f.hba_max ?? 10} onChange={v => setFilter('hba_max', v)} />
          </div>
          <div className="grid grid-cols-2 gap-3 mt-2">
            <RangeInputDual label="TPSA (A²)" min={0} max={300} step={5}
              valueMin={f.tpsa_min ?? 0} valueMax={f.tpsa_max ?? 140}
              onChangeMin={v => setFilter('tpsa_min', v)} onChangeMax={v => setFilter('tpsa_max', v)} unit=" A²" />
            <RangeInput label="RotBonds max" min={0} max={20} value={f.rotatable_max ?? 10} onChange={v => setFilter('rotatable_max', v)} />
          </div>
          <div className="grid grid-cols-2 gap-3 mt-2">
            <RangeInputDual label="Rings" min={0} max={10} step={1}
              valueMin={f.num_rings_min ?? 0} valueMax={f.num_rings_max ?? 6}
              onChangeMin={v => setFilter('num_rings_min', v)} onChangeMax={v => setFilter('num_rings_max', v)} />
            <RangeInputDual label="Heavy atoms" min={0} max={80} step={1}
              valueMin={f.num_heavy_atoms_min ?? 10} valueMax={f.num_heavy_atoms_max ?? 50}
              onChangeMin={v => setFilter('num_heavy_atoms_min', v)} onChangeMax={v => setFilter('num_heavy_atoms_max', v)} />
          </div>
        </div>

        {/* Max compounds — after filters */}
        <div>
          <p className="text-xs font-medium text-gray-600 mb-1">Max compounds to import</p>
          <input type="number" min={10} max={1000000} step={10}
            value={config.maxPerSource || 100}
            onChange={e => set('maxPerSource', Math.max(10, parseInt(e.target.value) || 100))}
            className="w-40 text-sm border border-gray-200 rounded-lg px-3 py-2 bg-white tabular-nums
                       focus:outline-none focus:ring-2 focus:ring-bx-mint/30 focus:border-bx-mint" />
          <p className="text-[10px] text-gray-400 mt-1">Up to 1,000,000 — large imports may take several minutes</p>
        </div>
      </div>
    )
  }

  if (db === 'chembl') {
    return (
      <div className="bg-gray-50 rounded-xl border border-gray-100 p-4 space-y-4">
        <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wide">ChEMBL Configuration</p>
        <div>
          <p className="text-xs font-medium text-gray-600 mb-1">Max compounds</p>
          <input type="number" min={10} max={10000} step={10}
            value={config.maxPerSource || 50}
            onChange={e => set('maxPerSource', Math.max(10, parseInt(e.target.value) || 50))}
            className="w-32 text-sm border border-gray-200 rounded-lg px-3 py-2 bg-white tabular-nums
                       focus:outline-none focus:ring-2 focus:ring-bx-mint/30 focus:border-bx-mint" />
        </div>
        <div className="grid grid-cols-2 gap-3">
          <div>
            <p className="text-[10px] text-gray-500 font-medium mb-1">Activity types</p>
            <div className="flex flex-wrap gap-1.5">
              {['IC50', 'Ki', 'EC50', 'Kd'].map(at => {
                const active = (f.activity_types || ['IC50', 'Ki']).includes(at)
                return (
                  <button key={at} onClick={() => {
                    const arr = f.activity_types || ['IC50', 'Ki']
                    setFilter('activity_types', active ? arr.filter(x => x !== at) : [...arr, at])
                  }}
                    className={`px-2 py-0.5 rounded text-[10px] font-semibold border transition-colors ${
                      active ? 'bg-blue-100 border-blue-300 text-blue-700' : 'bg-gray-50 border-gray-200 text-gray-500 hover:border-blue-200'
                    }`}
                  >{at}</button>
                )
              })}
            </div>
          </div>
          <div>
            <p className="text-[10px] text-gray-500 font-medium mb-1">Activity cutoff (nM)</p>
            <input type="number" min={1} max={100000} step={100}
              value={f.activity_cutoff ?? 10000}
              onChange={e => setFilter('activity_cutoff', parseInt(e.target.value) || 10000)}
              className="w-28 text-xs border border-gray-200 rounded-lg px-2 py-1.5 bg-white tabular-nums
                         focus:outline-none focus:ring-2 focus:ring-bx-mint/30 focus:border-bx-mint" />
          </div>
        </div>
        <RangeInput label="pChEMBL min" min={3} max={10} step={0.1}
          value={f.pchembl_min ?? 5.0} onChange={v => setFilter('pchembl_min', v)} />
      </div>
    )
  }

  if (db === 'pubchem') {
    return (
      <div className="bg-gray-50 rounded-xl border border-gray-100 p-4 space-y-4">
        <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wide">PubChem Configuration</p>
        <div>
          <p className="text-xs font-medium text-gray-600 mb-1">Max compounds</p>
          <input type="number" min={10} max={10000} step={10}
            value={config.maxPerSource || 50}
            onChange={e => set('maxPerSource', Math.max(10, parseInt(e.target.value) || 50))}
            className="w-32 text-sm border border-gray-200 rounded-lg px-3 py-2 bg-white tabular-nums
                       focus:outline-none focus:ring-2 focus:ring-bx-mint/30 focus:border-bx-mint" />
        </div>
        <div className="grid grid-cols-2 gap-3">
          <div>
            <p className="text-[10px] text-gray-500 font-medium mb-1">Activity type</p>
            <select
              value={f.pubchem_activity_type || 'IC50'}
              onChange={e => setFilter('pubchem_activity_type', e.target.value)}
              className="text-xs border border-gray-200 rounded-lg px-2 py-1.5 bg-white
                         focus:outline-none focus:ring-2 focus:ring-bx-mint/30 focus:border-bx-mint"
            >
              {['IC50', 'Ki', 'EC50', 'Potency'].map(t => <option key={t} value={t}>{t}</option>)}
            </select>
          </div>
          <FilterToggle label="Active only" checked={f.pubchem_active_only ?? true}
            onChange={v => setFilter('pubchem_active_only', v)} description="Only fetch 'active' compounds" />
        </div>
      </div>
    )
  }

  if (db === 'enamine') {
    return (
      <div className="bg-gray-50 rounded-xl border border-gray-100 p-4 space-y-4">
        <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wide">Enamine REAL Configuration</p>
        <div>
          <p className="text-xs font-medium text-gray-600 mb-1">Max compounds</p>
          <input type="number" min={10} max={10000} step={10}
            value={config.maxPerSource || 500}
            onChange={e => set('maxPerSource', Math.max(10, parseInt(e.target.value) || 500))}
            className="w-32 text-sm border border-gray-200 rounded-lg px-3 py-2 bg-white tabular-nums
                       focus:outline-none focus:ring-2 focus:ring-bx-mint/30 focus:border-bx-mint" />
        </div>
        <div>
          <p className="text-[10px] text-gray-500 font-medium mb-1">Scaffold complexity</p>
          <div className="flex gap-2">
            {['2-fragment', '3-fragment', 'mixed'].map(sc => (
              <button key={sc} onClick={() => setFilter('scaffold_complexity', sc)}
                className={`px-2.5 py-1 rounded-lg text-xs font-semibold border transition-colors ${
                  (f.scaffold_complexity || 'mixed') === sc ? 'bg-orange-100 border-orange-300 text-orange-700' : 'bg-gray-50 border-gray-200 text-gray-500 hover:border-orange-200'
                }`}
              >{sc}</button>
            ))}
          </div>
        </div>
      </div>
    )
  }

  if (db === 'fragments') {
    return (
      <div className="bg-gray-50 rounded-xl border border-gray-100 p-4 space-y-4">
        <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wide">Fragment Library Configuration</p>
        <FilterToggle label="Rule-of-3 strict" checked={f.ro3_strict ?? true}
          onChange={v => setFilter('ro3_strict', v)}
          description="MW<=300, logP<=3, HBD<=3, HBA<=3, RotBonds<=3" />
      </div>
    )
  }

  return null
}

// ---------------------------------------------------------------------------
// Docking time/cost estimate (benchmark-based)
// ---------------------------------------------------------------------------
function DockingEstimate({ nMols, exhaustiveness = 8, cnnScoring = 'rescore' }) {
  if (!nMols || nMols <= 0) return null

  // Throughput model (from GNINA GPU benchmark 2026-03-10)
  const baseThroughput = Math.min(30, 10 + 0.068 * nMols) // mol/min, sub-linear scaling
  const exhMultiplier = 1 + (exhaustiveness - 8) * 0.0125 // exh=32 → 1.3x
  const refMultiplier = cnnScoring === 'refinement' ? 2.7 : 1.0

  const effectiveThroughput = baseThroughput / (exhMultiplier * refMultiplier)
  const timeMin = nMols / effectiveThroughput
  const timeSec = timeMin * 60
  const cost = timeSec * 0.00016

  const fmtTime = timeMin < 1 ? `${Math.round(timeSec)}s` : timeMin < 60 ? `~${Math.round(timeMin)} min` : `~${(timeMin / 60).toFixed(1)}h`
  const fmtCost = cost < 0.01 ? `<$0.01` : `~$${cost.toFixed(2)}`
  const fmtThroughput = `${effectiveThroughput.toFixed(0)} mol/min`

  return (
    <div className="bg-indigo-50 border border-indigo-100 rounded-xl px-4 py-3 space-y-1">
      <p className="text-[10px] font-semibold text-indigo-400 uppercase tracking-wide">GPU Estimate</p>
      <div className="flex items-center gap-4 text-sm">
        <div>
          <span className="text-indigo-500">Time: </span>
          <span className="font-bold text-indigo-800">{fmtTime}</span>
        </div>
        <div>
          <span className="text-indigo-500">Cost: </span>
          <span className="font-bold text-indigo-800">{fmtCost}</span>
        </div>
        <div className="ml-auto text-xs text-indigo-400">{fmtThroughput}</div>
      </div>
    </div>
  )
}

// Step 3: Confirmation
// ---------------------------------------------------------------------------
function ConfirmationView({ runType, config, selectedCount, includedColumns, onIncludedColumnsChange, runAllMolecules, onRunAllMoleculesChange, onConfigChange, gpuEngines, project }) {
  const rt = RUN_TYPES.find(t => t.type === runType)
  if (!rt) return null
  const c = runColor(runType)

  const configLines = []
  if (runType === 'import') {
    const db = config.database
    const sourceLabels = { zinc: 'ZINC20', chembl: 'ChEMBL', pubchem: 'PubChem', enamine: 'Enamine REAL', fragments: 'Fragment Library', external: 'File Upload', internal: 'Phase Bookmarks' }
    if (config.sourceMode === 'database' && db) {
      configLines.push(`Source: ${sourceLabels[db] || db}`)
    } else if (config.sourceMode === 'internal') {
      configLines.push('Source: Phase bookmarks')
    } else if (config.sourceMode === 'external') {
      configLines.push(`Source: File upload${config._file ? ` (${config._file.name})` : ''}`)
    }
    if (config.sourceName) configLines.push(`Label: ${config.sourceName}`)
  }
  if (runType === 'calculation') {
    const calcTypes = config.calculation_types || []
    const calcLabel = calcTypes[0] ? (CALCULATION_SUBTYPES.find(s => s.key === calcTypes[0])?.label || calcTypes[0]) : 'None'
    configLines.push(`Calculation: ${calcLabel}`)
    // Docking/scoring details shown as interactive panels below — no text summary needed
  }
  if (runType === 'docking') {
    configLines.push(`Engine: ${config.engine === 'gnina_gpu' ? 'GNINA (GPU)' : config.engine === 'gnina_cpu' ? 'GNINA (CPU)' : 'AutoDock Vina'}`)
    configLines.push(`Exhaustiveness: ${config.exhaustiveness ?? 32}  |  Poses: ${config.num_modes ?? 9}`)
    const boxX = config.boxSizeX || 'auto', boxY = config.boxSizeY || 'auto', boxZ = config.boxSizeZ || 'auto'
    configLines.push(`Box: ${boxX} x ${boxY} x ${boxZ} A`)
  }
  if (runType === 'adme') {
    configLines.push('Full pharmacokinetics profiling (absorption, distribution, metabolism, excretion)')
  }
  if (runType === 'toxicity') {
    configLines.push('Full toxicological risk assessment (hERG, Ames, hepatotoxicity, carcinogenicity)')
  }
  if (runType === 'scoring') {
    configLines.push('Physicochemical descriptors (MW, LogP, TPSA, QED, Lipinski, SA Score)')
  }
  if (runType === 'composite') {
    const ws = config.weights || {}
    const total = Object.values(ws).reduce((a, b) => a + b, 0)
    configLines.push(`${Object.keys(ws).length} score weights, total: ${total.toFixed(2)}`)
    configLines.push('Aggregates results from previous runs (docking, ADME, scoring, etc.)')
  }
  if (runType === 'interactions') {
    configLines.push('Protein-ligand interaction analysis (H-bonds, hydrophobic, pi-stacking)')
  }
  if (runType === 'scaffold' || (runType === 'calculation' && config.calculation_types?.[0] === 'scaffold')) {
    configLines.push('Murcko scaffold decomposition + BRICS bonds + R-group positions')
    configLines.push(`Butina diversity clustering (cutoff: ${Number(config.scaffold?.cutoff ?? 0.65).toFixed(2)})`)
  }
  if (runType === 'generation') {
    if (config.mode === 'optimization') {
      configLines.push('Mode: Lead Optimization (multi-iteration)')
      const w = config.weights || {}
      const wLabels = Object.entries(w).map(([k, v]) => `${k.replace(/_/g, ' ')}: ${(v * 100).toFixed(0)}%`).join(', ')
      configLines.push(`Weights: ${wLabels}`)
      configLines.push(`${config.opt_iterations ?? 5} iterations x ${config.opt_variants_per_iteration ?? 50} variants/iter`)
    } else {
      configLines.push(`Method: ${config.method}`)
      configLines.push(`${config.iterations ?? 3} iterations x ${config.variants_per_iteration ?? 5} variants = ~${(config.iterations ?? 3) * (config.variants_per_iteration ?? 5)} molecules`)
    }
  }
  if (runType === 'retrosynthesis' || (runType === 'calculation' && config.calculation_types?.[0] === 'retrosynthesis')) {
    const m = config.retrosynthesis?.method || 'aizynthfinder'
    configLines.push(`Method: ${m === 'aizynthfinder' ? 'AiZynthFinder (AI)' : 'RDKit Heuristic (fast)'}`)
  }

  return (
    <div className="space-y-4">
      <div className={`rounded-xl border p-4 flex items-start gap-4 ${c.bg} border-current/10`}>
        <div className={`flex-shrink-0 p-2.5 rounded-xl ${c.fill} text-white`}>
          <RunTypeIcon icon={rt.icon} className="w-6 h-6" />
        </div>
        <div className="min-w-0">
          <p className="font-bold text-gray-800 text-base">{rt.label}</p>
          <p className="text-sm text-gray-500 mt-0.5">{rt.description}</p>
        </div>
      </div>

      <div className="bg-gray-50 rounded-xl border border-gray-100 p-4 space-y-2">
        <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wide">Configuration</p>
        {configLines.map((l, i) => (
          <div key={i} className="flex items-start gap-2 text-sm text-gray-700">
            <span className="text-gray-300 flex-shrink-0 mt-0.5">--</span>
            <span>{l}</span>
          </div>
        ))}
      </div>

      {/* Import source-specific config (editable in review/Step 3) */}
      {runType === 'import' && config.sourceMode === 'database' && onConfigChange && (
        <ImportSourceConfig config={config} onConfigChange={onConfigChange} />
      )}

      {/* Import file upload (in Step 3 for external source) */}
      {runType === 'import' && config.sourceMode === 'external' && onConfigChange && (
        <div className="bg-gray-50 rounded-xl border border-gray-100 p-4 space-y-3">
          <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wide">File Upload</p>
          <div
            onClick={() => document.getElementById('import-file-input')?.click()}
            className="border-2 border-dashed border-gray-200 rounded-xl p-6 text-center cursor-pointer
                       hover:border-bx-mint/40 hover:bg-blue-50/30 transition-all group"
          >
            {config._file ? (
              <div className="flex items-center justify-center gap-2">
                <svg className="w-5 h-5 text-bx-mint" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 13l4 4L19 7" />
                </svg>
                <span className="text-sm font-semibold text-gray-700">{config._file.name}</span>
                <span className="text-xs text-gray-400">({(config._file.size / 1024).toFixed(1)} KB)</span>
              </div>
            ) : (
              <>
                <svg className="w-9 h-9 text-gray-300 mx-auto mb-2 group-hover:text-bx-light-text/40 transition-colors"
                  fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                    d="M7 16a4 4 0 01-.88-7.903A5 5 0 1115.9 6L16 6a5 5 0 011 9.9M15 13l-3-3m0 0l-3 3m3-3v12" />
                </svg>
                <p className="text-sm font-medium text-gray-500 mb-1">Drop file here or click to browse</p>
                <div className="flex justify-center gap-2 mt-2">
                  {['SDF', 'SMILES', 'CSV'].map(f => (
                    <span key={f} className="px-2 py-0.5 rounded text-[10px] font-bold bg-gray-100 text-gray-500">.{f.toLowerCase()}</span>
                  ))}
                </div>
              </>
            )}
            <input
              id="import-file-input"
              type="file"
              accept=".sdf,.smi,.smiles,.csv"
              className="hidden"
              onChange={e => {
                const f = e.target.files?.[0]
                if (f) onConfigChange({ ...config, _file: f, sourceName: config.sourceName || f.name })
              }}
            />
          </div>
        </div>
      )}

      {/* Per-subtype config panels (editable in review) */}
      {runType === 'calculation' && config.calculation_types?.[0] === 'docking' && onConfigChange && (
        <div className="bg-gray-50 rounded-xl border border-gray-100 p-4 space-y-3">
          <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wide">Docking Configuration</p>
          <div className="space-y-2">
            <EngineCard id="gnina_gpu" label="GNINA (GPU)"
              badges={[{ text: 'Fastest', color: 'bg-green-100 text-green-700' }]}
              selected={(config.docking?.engine || 'gnina_gpu') === 'gnina_gpu'}
              onClick={v => onConfigChange({ ...config, docking: { ...config.docking, engine: v } })}
              disabled={gpuEngines && !gpuEngines.gnina_gpu}
              disabledReason="RunPod not configured" />
            <EngineCard id="gnina_cpu" label="GNINA (CPU)"
              badges={[{ text: 'CNN rescoring', color: 'bg-blue-100 text-blue-700' }]}
              selected={(config.docking?.engine) === 'gnina_cpu'}
              onClick={v => onConfigChange({ ...config, docking: { ...config.docking, engine: v } })} />
            <EngineCard id="vina" label="AutoDock Vina"
              badges={[{ text: 'Classic', color: 'bg-gray-100 text-gray-600' }]}
              selected={(config.docking?.engine) === 'vina'}
              onClick={v => onConfigChange({ ...config, docking: { ...config.docking, engine: v } })} />
          </div>
          <div className="mt-2">
            <ExhaustivenessSlider
              value={config.docking?.exhaustiveness ?? 32}
              onChange={v => onConfigChange({ ...config, docking: { ...config.docking, exhaustiveness: v } })}
            />
          </div>
          {((config.docking?.engine || 'gnina_gpu') !== 'vina') && (
            <div className="mt-3">
              <label className="text-sm font-medium text-gray-600 mb-1 block">CNN scoring mode</label>
              <div className="flex gap-2">
                {[
                  { id: 'rescore', label: 'Rescore', desc: 'Fast — CNN rescores Vina poses' },
                  { id: 'refinement', label: 'Refinement', desc: 'Better quality — CNN optimizes poses (+50% time)' },
                ].map(mode => (
                  <button key={mode.id}
                    onClick={() => onConfigChange({ ...config, docking: { ...config.docking, cnn_scoring: mode.id } })}
                    className={`flex-1 text-left p-2 rounded-lg border text-xs transition-all ${
                      (config.docking?.cnn_scoring || 'rescore') === mode.id
                        ? 'border-blue-300 bg-blue-50 text-blue-700'
                        : 'border-gray-200 text-gray-500 hover:border-gray-300'
                    }`}>
                    <div className="font-medium">{mode.label}</div>
                    <div className="text-[10px] opacity-70">{mode.desc}</div>
                  </button>
                ))}
              </div>
              {(config.docking?.cnn_scoring) === 'refinement' && (
                <div className="flex items-center gap-1.5 mt-2 text-[10px] text-amber-600 bg-amber-50 rounded-lg px-2.5 py-1.5 border border-amber-200">
                  <svg className="w-3.5 h-3.5 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
                  </svg>
                  ~3x slower, may reject ~17% of molecules (positive Vina scores)
                </div>
              )}
            </div>
          )}
        </div>
      )}

      {/* Docking estimate — GPU runs only */}
      {runType === 'calculation' && config.calculation_types?.[0] === 'docking' && (config.docking?.engine || 'gnina_gpu') !== 'vina' && (
        <DockingEstimate
          nMols={runAllMolecules ? null : selectedCount}
          exhaustiveness={config.docking?.exhaustiveness ?? 32}
          cnnScoring={config.docking?.cnn_scoring || 'rescore'}
        />
      )}

      {runType === 'calculation' && config.calculation_types?.[0] === 'composite' && onConfigChange && (
        <div className="bg-gray-50 rounded-xl border border-gray-100 p-4 space-y-2">
          <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wide">Score Weights</p>
          <p className="text-xs text-gray-400 mb-2">
            Adjust the importance of each metric. The composite score aggregates results from your previous runs.
          </p>
          <ScoringWeightsEditor
            weights={config.composite?.weights || DEFAULT_CONFIGS.calculation.composite.weights}
            onChange={(w) => onConfigChange({ ...config, composite: { ...config.composite, weights: w } })}
          />
        </div>
      )}

      {/* Interactions: functional residues checklist */}
      {runType === 'calculation' && config.calculation_types?.[0] === 'interactions' && onConfigChange && (() => {
        const funcRes = project?.target_preview?.functional_residues?.residues || []
        const keyHbondNums = new Set(project?.target_preview?.functional_residues?.key_hbond_residues || [])
        const excluded = config.interactions?.excluded_residues || []
        const activeCount = funcRes.filter(r => r.enabled !== false && !excluded.includes(r.number)).length

        const toggleExclude = (num) => {
          const prev = config.interactions?.excluded_residues || []
          const next = prev.includes(num) ? prev.filter(n => n !== num) : [...prev, num]
          onConfigChange({ ...config, interactions: { ...config.interactions, excluded_residues: next } })
        }

        return (
          <div className="bg-gray-50 rounded-xl border border-gray-100 p-4 space-y-2">
            <div className="flex items-center justify-between">
              <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wide">Functional Residues</p>
              {funcRes.length > 0 && (
                <span className="text-[10px] font-medium px-2 py-0.5 rounded-full bg-lime-100 text-lime-700">
                  {activeCount}/{funcRes.length} active
                </span>
              )}
            </div>
            {funcRes.length === 0 ? (
              <p className="text-xs text-gray-400">No functional residues defined — all contacts will be reported without quality scoring.</p>
            ) : (
              <>
                <p className="text-[10px] text-gray-400 mb-1">Uncheck residues to exclude from quality scoring</p>
                <div className="grid grid-cols-2 gap-1.5 max-h-48 overflow-y-auto pr-1">
                  {funcRes.filter(r => r.enabled !== false).map(r => {
                    const isExcluded = excluded.includes(r.number)
                    const isKeyHB = keyHbondNums.has(r.number)
                    return (
                      <label
                        key={r.number}
                        className={`flex items-center gap-2 px-2 py-1.5 rounded-lg border cursor-pointer transition-colors ${
                          isExcluded
                            ? 'border-gray-200 bg-white opacity-50'
                            : 'border-lime-200 bg-lime-50/50 hover:bg-lime-50'
                        }`}
                      >
                        <input
                          type="checkbox"
                          checked={!isExcluded}
                          onChange={() => toggleExclude(r.number)}
                          className="w-3.5 h-3.5 rounded border-gray-300 text-lime-600 focus:ring-lime-500"
                        />
                        <span className="text-xs font-mono font-medium text-gray-700">
                          {r.aa || '?'}{r.number}
                        </span>
                        {isKeyHB && (
                          <span className="text-[9px] font-semibold px-1 py-0.5 rounded bg-blue-100 text-blue-600 ml-auto">HB</span>
                        )}
                      </label>
                    )
                  })}
                </div>
              </>
            )}
          </div>
        )
      })()}

      {/* Retrosynthesis method selector */}
      {(runType === 'retrosynthesis' || (runType === 'calculation' && config.calculation_types?.[0] === 'retrosynthesis')) && onConfigChange && (
        <div className="bg-gray-50 rounded-xl border border-gray-100 p-4 space-y-3">
          <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wide">Retrosynthesis Method</p>
          <div className="space-y-2">
            <EngineCard id="aizynthfinder" label="AiZynthFinder (AI)"
              badges={[
                { text: 'Recommended', color: 'bg-green-100 text-green-700' },
                { text: '~1-2 min/mol', color: 'bg-amber-100 text-amber-700' },
              ]}
              selected={(config.retrosynthesis?.method || 'aizynthfinder') === 'aizynthfinder'}
              onClick={v => onConfigChange({ ...config, retrosynthesis: { ...config.retrosynthesis, method: v } })}
            />
            {(config.retrosynthesis?.method || 'aizynthfinder') === 'aizynthfinder' && (
              <p className="text-[11px] text-gray-500 pl-6 -mt-1">
                Neural network trained on ~1M USPTO reactions. Plans multi-step routes with stock verification (ZINC).
              </p>
            )}
            <EngineCard id="rdkit" label="RDKit Heuristic (fast)"
              badges={[
                { text: 'Instant', color: 'bg-blue-100 text-blue-700' },
                { text: '~1 sec/mol', color: 'bg-gray-100 text-gray-600' },
              ]}
              selected={config.retrosynthesis?.method === 'rdkit'}
              onClick={v => onConfigChange({ ...config, retrosynthesis: { ...config.retrosynthesis, method: v } })}
            />
            {config.retrosynthesis?.method === 'rdkit' && (
              <p className="text-[11px] text-gray-500 pl-6 -mt-1">
                SMARTS pattern matching with 7 reaction types (Suzuki, amidation, etc.). Fast but limited to common disconnections.
              </p>
            )}
          </div>
        </div>
      )}

      {/* Scaffold: clustering cutoff slider */}
      {runType === 'calculation' && config.calculation_types?.[0] === 'scaffold' && onConfigChange && (() => {
        const cutoff = config.scaffold?.cutoff ?? 0.65
        const similarity = Math.round((1 - cutoff) * 100)
        return (
          <div className="bg-gray-50 rounded-xl border border-gray-100 p-4 space-y-3">
            <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wide">Clustering Cutoff</p>
            <div>
              <div className="flex items-center justify-between text-sm mb-1">
                <span className="font-medium text-gray-600">Cutoff</span>
                <span className="font-bold text-bx-light-text tabular-nums">{Number(cutoff).toFixed(2)}</span>
              </div>
              <input
                type="range" min={0.3} max={0.8} step={0.05}
                value={cutoff}
                onChange={e => onConfigChange({ ...config, scaffold: { ...config.scaffold, cutoff: parseFloat(e.target.value) } })}
                className="w-full accent-bx-mint"
              />
              <div className="flex justify-between text-[10px] text-gray-400 mt-0.5">
                <span>0.30 (strict)</span>
                <span>0.80 (loose)</span>
              </div>
              <div className="mt-2 bg-emerald-50 border border-emerald-100 rounded-lg px-3 py-2 text-sm text-emerald-700">
                Molecules with <span className="font-bold">&ge;{similarity}%</span> Tanimoto similarity will be grouped together
              </div>
            </div>
          </div>
        )
      })()}

      {/* Column checklist for calculation runs */}
      {runType === 'calculation' && config.calculation_types?.length === 1 && includedColumns && onIncludedColumnsChange && (
        <ColumnChecklist
          calcGroupKey={config.calculation_types[0]}
          includedColumns={includedColumns}
          onChange={onIncludedColumnsChange}
          config={config}
        />
      )}

      {/* Validation warnings */}
      {runType === 'calculation' && (!config.calculation_types || config.calculation_types.length === 0) && (
        <div className="flex items-center gap-2 bg-red-50 border border-red-200 rounded-xl px-4 py-2.5 text-sm text-red-600">
          <svg className="w-4 h-4 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
          </svg>
          No calculation type selected. Go back and select one.
        </div>
      )}
      {/* Molecule input: status + run-all toggle */}
      {(runType === 'calculation' || runType === 'generation') && (
        <div className="space-y-2">
          {runAllMolecules ? (
            <div className="flex items-center gap-2 bg-green-50 border border-green-200 rounded-xl px-4 py-2.5 text-sm text-green-700">
              <svg className="w-4 h-4 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M5 13l4 4L19 7" />
              </svg>
              All molecules in phase will be processed
            </div>
          ) : selectedCount > 0 ? (
            <div className="flex items-center gap-2 bg-green-50 border border-green-200 rounded-xl px-4 py-2.5 text-sm text-green-700">
              <svg className="w-4 h-4 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M5 13l4 4L19 7" />
              </svg>
              {selectedCount} molecule{selectedCount > 1 ? 's' : ''} selected as input
            </div>
          ) : (
            <div className="flex items-center gap-2 bg-amber-50 border border-amber-200 rounded-xl px-4 py-2.5 text-sm text-amber-700">
              <svg className="w-4 h-4 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
              </svg>
              No molecules selected. Use the toggle below or select molecules in the dashboard.
            </div>
          )}
          {runType === 'calculation' && onRunAllMoleculesChange && (
            <label className="flex items-center gap-2 px-1 py-1 text-sm text-gray-600 cursor-pointer select-none hover:text-gray-800">
              <input
                type="checkbox"
                checked={runAllMolecules || false}
                onChange={e => onRunAllMoleculesChange(e.target.checked)}
                className="rounded border-gray-300 text-bx-accent focus:ring-bx-accent/30 w-3.5 h-3.5"
              />
              Run on all molecules in phase (ignore selection)
            </label>
          )}
        </div>
      )}

      <div className="flex items-center gap-3 bg-blue-50 border border-blue-100 rounded-xl px-4 py-3">
        <svg className="w-4 h-4 text-blue-400 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
            d="M12 8v4l3 3m6-3a9 9 0 11-18 0 9 9 0 0118 0z" />
        </svg>
        <div className="text-sm">
          <span className="text-blue-600">Estimated time: </span>
          <span className="font-bold text-blue-800">{ESTIMATED_TIMES[runType] || '~1 minute'}</span>
        </div>
        <div className="ml-auto text-[10px] text-blue-400">
          {(runType === 'docking' && config.engine === 'gnina_gpu') ? 'GPU (RunPod)' : 'CPU (local)'}
        </div>
      </div>
    </div>
  )
}

// ---------------------------------------------------------------------------
// RunCreator — main component
// ---------------------------------------------------------------------------
export default function RunCreator({ phaseId, phaseType, isOpen, onClose, onSubmit, selectedMoleculeIds, submitting, project, campaign, runs = [] }) {
  const [step, setStep] = useState(1)
  const [selectedType, setSelectedType] = useState(null)
  const [config, setConfig] = useState({})
  const [includedColumns, setIncludedColumns] = useState(null)
  const [runAllMolecules, setRunAllMolecules] = useState(false)
  const [gpuEngines, setGpuEngines] = useState(null)

  useEffect(() => {
    if (!isOpen) return
    v9GpuHealth()
      .then(data => setGpuEngines(data.engines || {}))
      .catch(() => setGpuEngines({ gnina_gpu: false, gnina_cpu: true, vina: true, diffdock: false }))
  }, [isOpen])

  const hasDockingRun = runs.some(
    r => r.status === 'completed' && r.calculation_types?.includes('docking')
  )

  if (!isOpen) return null

  const allRunTypes = RUN_TYPES
  const STEP_LABELS = ['Type', 'Configure', 'Launch']

  function handleSelectType(type) {
    setSelectedType(type)
    setConfig(DEFAULT_CONFIGS[type] || {})
    setStep(2)
  }

  const needsReceptorPrep = selectedType === 'calculation' && config.calculation_types?.[0] === 'docking' && !project?.receptor_prep_report
  const needsPharmaModel = selectedType === 'calculation' && config.calculation_types?.[0] === 'pharmacophore' && !campaign?.reference_ligands?.length

  const canSubmit = (() => {
    if (submitting) return false
    if (needsReceptorPrep) return false
    if (needsPharmaModel) return false
    if (selectedType === 'import') {
      if (config.sourceMode === 'external' && !config._file) return false
      if (config.sourceMode === 'database' && !config.database) return false
    }
    if (selectedType === 'calculation') {
      if (!config.calculation_types?.length) return false
      if (!runAllMolecules && !selectedMoleculeIds?.size) return false
    }
    if (selectedType === 'generation' && !selectedMoleculeIds?.size) return false
    if (selectedType === 'afvs' && !(config.budget_cap_usd > 0)) return false
    return true
  })()

  function handleSubmit() {
    if (!canSubmit) return
    const payload = { type: selectedType, config, phase_id: phaseId }
    // For calculation runs, include calculation_types and selected molecules
    if (selectedType === 'calculation') {
      payload.calculation_types = config.calculation_types || []
      if (!runAllMolecules && selectedMoleculeIds?.size > 0) {
        payload.input_molecule_ids = [...selectedMoleculeIds]
      }
    }
    if (selectedType === 'afvs') {
      payload.config = { ...payload.config, _afvs: true }
    }
    if (selectedType === 'generation') {
      if (selectedMoleculeIds?.size > 0) {
        payload.input_molecule_ids = [...selectedMoleculeIds]
      }
      // Map optimization-specific config keys for backend
      if (config.mode === 'optimization') {
        payload.config = {
          ...config,
          iterations: config.opt_iterations ?? 5,
          variants_per_iteration: config.opt_variants_per_iteration ?? 50,
        }
      }
    }
    // Flatten composite weights to top-level config for backend
    if (selectedType === 'calculation' && (config.calculation_types || [])[0] === 'composite') {
      payload.config = { ...payload.config, weights: config.composite?.weights || DEFAULT_CONFIGS.calculation.composite.weights }
    }
    if (includedColumns) {
      payload.config = { ...payload.config, included_columns: includedColumns }
    }
    if (runAllMolecules) {
      payload.config = { ...payload.config, run_all_molecules: true }
    }
    const excludedRes = config.interactions?.excluded_residues
    if (excludedRes?.length > 0) {
      payload.config = { ...payload.config, excluded_residues: excludedRes }
    }
    if (onSubmit) onSubmit(payload)
    handleClose()
  }

  function handleClose() {
    setStep(1)
    setSelectedType(null)
    setConfig({})
    setIncludedColumns(null)
    setRunAllMolecules(false)
    onClose()
  }

  const selectedRt = selectedType ? RUN_TYPES.find(t => t.type === selectedType) : null

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/50 backdrop-blur-sm p-4">
      <div
        className="bg-white rounded-2xl shadow-2xl w-full max-w-2xl max-h-[90vh] flex flex-col overflow-hidden"
        onClick={e => e.stopPropagation()}
      >
        {/* Header */}
        <div className="flex items-center justify-between px-6 py-4 border-b border-gray-100 flex-shrink-0">
          <div className="flex items-center gap-3">
            {selectedRt && step > 1 && (
              <div className={`p-1.5 rounded-lg ${runColor(selectedType).fill} text-white`}>
                <RunTypeIcon icon={selectedRt.icon} className="w-4 h-4" />
              </div>
            )}
            <div>
              <h2 className="font-bold text-bx-light-text text-lg">
                {step === 1 ? 'New Run' : step === 2 ? `Configure: ${selectedRt?.label}` : 'Review & Launch'}
              </h2>
              {phaseType && (
                <p className="text-sm text-gray-400">
                  {phaseType.replace(/_/g, ' ')}
                </p>
              )}
            </div>
          </div>
          <button
            onClick={handleClose}
            className="p-1.5 rounded-lg hover:bg-gray-100 text-gray-500 transition-colors"
          >
            <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
            </svg>
          </button>
        </div>

        {/* Body */}
        <div className="flex-1 overflow-y-auto px-6 py-5">
          <StepIndicator current={step} labels={STEP_LABELS} />

          {/* Step 1: Choose type */}
          {step === 1 && (
            <div>
              <p className="text-sm text-gray-500 mb-4 text-center">
                Choose the type of analysis to run on your molecules.
              </p>
              <div className="grid grid-cols-2 gap-3">
                {allRunTypes.map(rt => {
                  const available = !phaseType || rt.phases.includes(phaseType)
                  const c = runColor(rt.type)
                  return (
                    <button
                      key={rt.type}
                      onClick={() => available && handleSelectType(rt.type)}
                      className={`flex items-start gap-3 p-4 rounded-xl border-2 text-left transition-all ${
                        available
                          ? `border-gray-100 hover:border-blue-300 hover:shadow-sm hover:bg-blue-50/30 group`
                          : 'border-gray-100 opacity-40 cursor-not-allowed'
                      }`}
                    >
                      <span className={`flex-shrink-0 p-2 rounded-xl ${c.bg} ${c.text} transition-colors`}>
                        <RunTypeIcon icon={rt.icon} className="w-5 h-5" />
                      </span>
                      <div className="min-w-0 flex-1">
                        <p className="font-semibold text-gray-800 text-sm">{rt.label}</p>
                        <p className="text-sm text-gray-500 mt-0.5 leading-relaxed">{rt.description}</p>
                        {!available && (
                          <p className="text-[10px] text-gray-400 mt-1 italic">
                            Available in {rt.phases.filter(p => p !== phaseType).map(p => p.replace(/_/g, ' ')).join(', ')}
                          </p>
                        )}
                      </div>
                    </button>
                  )
                })}
              </div>
            </div>
          )}

          {/* Step 2: Configure */}
          {step === 2 && selectedType && (
            <ConfigForm
              runType={selectedType}
              config={config}
              onChange={setConfig}
              phase={phaseType}
              selectedCount={selectedMoleculeIds?.size || 0}
              gpuEngines={gpuEngines}
              project={project}
              campaign={campaign}
              hasDockingRun={hasDockingRun}
            />
          )}

          {/* Step 3: Confirm */}
          {step === 3 && selectedType && (
            <ConfirmationView runType={selectedType} config={config} selectedCount={selectedMoleculeIds?.size || 0} includedColumns={includedColumns} onIncludedColumnsChange={setIncludedColumns} runAllMolecules={runAllMolecules} onRunAllMoleculesChange={setRunAllMolecules} onConfigChange={setConfig} gpuEngines={gpuEngines} project={project} />
          )}
        </div>

        {/* Footer */}
        <div className="flex items-center justify-between px-6 py-4 border-t border-gray-100 bg-gray-50 flex-shrink-0">
          <button
            onClick={() => step > 1 ? setStep(s => s - 1) : handleClose()}
            className="px-4 py-2 rounded-xl text-sm font-medium border border-gray-200 text-gray-700 hover:bg-gray-100 transition-colors"
          >
            {step === 1 ? 'Cancel' : 'Back'}
          </button>

          {step === 2 && (
            <button
              onClick={() => {
                // Initialize column checklist for calculation runs (filtered by config options like engine)
                if (selectedType === 'calculation' && config.calculation_types?.length === 1) {
                  setIncludedColumns(getDefaultCheckedColumns(config.calculation_types[0], config))
                }
                setStep(3)
              }}
              className="flex items-center gap-2 px-5 py-2 rounded-xl text-sm font-semibold bg-bx-surface hover:bg-bx-elevated text-white transition-colors"
            >
              Review
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
              </svg>
            </button>
          )}

          {step === 3 && (
            <button
              onClick={handleSubmit}
              disabled={!canSubmit}
              className={`flex items-center gap-2 px-5 py-2 rounded-xl text-sm font-semibold transition-colors shadow-sm shadow-green-200 ${
                !canSubmit
                  ? 'bg-gray-300 text-gray-500 cursor-not-allowed'
                  : 'bg-bx-mint hover:bg-green-600 text-white'
              }`}
            >
              {submitting ? (
                <>
                  <BindXLogo variant="loading" size={16} />
                  Launching...
                </>
              ) : (
                <>
                  <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                      d="M5.25 5.653c0-.856.917-1.398 1.667-.986l11.54 6.347a1.125 1.125 0 010 1.972l-11.54 6.347a1.125 1.125 0 01-1.667-.986V5.653z" />
                  </svg>
                  Launch Run
                </>
              )}
            </button>
          )}
        </div>
      </div>
    </div>
  )
}
