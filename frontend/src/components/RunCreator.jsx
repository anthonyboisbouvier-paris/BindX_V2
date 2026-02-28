import React, { useState, useRef } from 'react'
import { RUN_TYPES, CALCULATION_SUBTYPES } from '../mock/data.js'
import Badge from './Badge.jsx'
import BindXLogo from './BindXLogo.jsx'

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
    case 'grid':
      return <svg {...p}><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
        d="M3.75 6A2.25 2.25 0 016 3.75h2.25A2.25 2.25 0 0110.5 6v2.25a2.25 2.25 0 01-2.25 2.25H6a2.25 2.25 0 01-2.25-2.25V6zM3.75 15.75A2.25 2.25 0 016 13.5h2.25a2.25 2.25 0 012.25 2.25V18a2.25 2.25 0 01-2.25 2.25H6A2.25 2.25 0 013.75 18v-2.25zM13.5 6a2.25 2.25 0 012.25-2.25H18A2.25 2.25 0 0120.25 6v2.25A2.25 2.25 0 0118 10.5h-2.25a2.25 2.25 0 01-2.25-2.25V6zM13.5 15.75a2.25 2.25 0 012.25-2.25H18a2.25 2.25 0 012.25 2.25V18A2.25 2.25 0 0118 20.25h-2.25A2.25 2.25 0 0113.5 18v-2.25z" /></svg>
    case 'calculator':
      return <svg {...p}><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
        d="M15.75 15.75V18m-7.5-6.75h.008v.008H8.25v-.008zm0 2.25h.008v.008H8.25V13.5zm0 2.25h.008v.008H8.25v-.008zm0 2.25h.008v.008H8.25V18zm2.498-6.75h.008v.008h-.008v-.008zm0 2.25h.008v.008h-.008V13.5zm0 2.25h.008v.008h-.008v-.008zm0 2.25h.008v.008h-.008V18zm2.504-6.75h.008v.008h-.008v-.008zm0 2.25h.008v.008h-.008V13.5zm0 2.25h.008v.008h-.008v-.008zm2.498-6.75h.008v.008H15.75v-.008zm0 2.25h.008v.008H15.75V13.5zM8.25 6h7.5v2.25h-7.5V6zM6 20.25h12A1.5 1.5 0 0019.5 18.75V5.25A1.5 1.5 0 0018 3.75H6A1.5 1.5 0 004.5 5.25v13.5A1.5 1.5 0 006 20.25z" /></svg>
    default:
      return <svg {...p}><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
        d="M5.25 5.653c0-.856.917-1.398 1.667-.986l11.54 6.347a1.125 1.125 0 010 1.972l-11.54 6.347a1.125 1.125 0 01-1.667-.986V5.653z" /></svg>
  }
}

// ---------------------------------------------------------------------------
// Color per run type
// ---------------------------------------------------------------------------
const RUN_TYPE_COLORS = {
  import:      { text: 'text-gray-500',   bg: 'bg-gray-100',    ring: 'ring-gray-200',   fill: 'bg-gray-500'   },
  calculation: { text: 'text-blue-700',   bg: 'bg-blue-100',    ring: 'ring-blue-200',   fill: 'bg-blue-600'   },
  generation:  { text: 'text-pink-700',   bg: 'bg-pink-100',    ring: 'ring-pink-200',   fill: 'bg-pink-600'   },
  // Subtypes (used in RunHistory display)
  docking:     { text: 'text-blue-700',   bg: 'bg-blue-100',    ring: 'ring-blue-200',   fill: 'bg-blue-600'   },
  admet:       { text: 'text-green-700',  bg: 'bg-green-100',   ring: 'ring-green-200',  fill: 'bg-green-600'  },
  scoring:     { text: 'text-amber-700',  bg: 'bg-amber-100',   ring: 'ring-amber-200',  fill: 'bg-amber-500'  },
  enrichment:  { text: 'text-purple-700', bg: 'bg-purple-100',  ring: 'ring-purple-200', fill: 'bg-purple-600' },
  clustering:  { text: 'text-teal-700',   bg: 'bg-teal-100',    ring: 'ring-teal-200',   fill: 'bg-teal-600'   },
  off_target:  { text: 'text-red-700',    bg: 'bg-red-100',     ring: 'ring-red-200',    fill: 'bg-red-600'    },
  confidence:  { text: 'text-indigo-700', bg: 'bg-indigo-100',  ring: 'ring-indigo-200', fill: 'bg-indigo-600' },
  retrosynthesis: { text: 'text-orange-700', bg: 'bg-orange-100', ring: 'ring-orange-200', fill: 'bg-orange-500' },
  safety:      { text: 'text-rose-700',   bg: 'bg-rose-100',    ring: 'ring-rose-200',   fill: 'bg-rose-600'   },
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
  import:      { sourceMode: 'external', format: 'sdf', sourceName: '' },
  calculation: {
    calculation_types: [],
    // Per-subtype configs
    docking: { engine: 'gnina_gpu', exhaustiveness: 32, num_modes: 9, seed: 0, boxSizeX: '', boxSizeY: '', boxSizeZ: '' },
    admet: { properties: ['logP', 'solubility', 'BBB', 'hERG', 'metabolic_stability', 'CYP', 'oral_bioavailability', 'ames'] },
    scoring: { weights: { docking_score: 0.30, cnn_score: 0.20, logP: 0.15, solubility: 0.10, selectivity: 0.15, novelty: 0.10 } },
    enrichment: { analyses: ['prolif', 'clustering', 'scaffold', 'pharmacophore'] },
    clustering: { method: 'butina', cutoff: 0.5 },
    off_target: {},
    confidence: {},
    retrosynthesis: {},
    safety: {},
  },
  generation:  { method: 'scaffold_hopping', iterations: 3, variants_per_iteration: 5, include_docking: true, include_admet: true, include_scoring: true },
}

// Keep legacy keys for ConfigForm compatibility
const LEGACY_DEFAULT_CONFIGS = {
  docking:    { engine: 'gnina_gpu', exhaustiveness: 32, num_modes: 9, seed: 0, boxSizeX: '', boxSizeY: '', boxSizeZ: '' },
  admet:      { properties: ['logP', 'solubility', 'BBB', 'hERG', 'metabolic_stability', 'CYP', 'oral_bioavailability', 'ames'] },
  scoring:    { weights: { docking_score: 0.30, cnn_score: 0.20, logP: 0.15, solubility: 0.10, selectivity: 0.15, novelty: 0.10 } },
  enrichment: { analyses: ['prolif', 'clustering', 'scaffold', 'pharmacophore'] },
  clustering: { method: 'butina', cutoff: 0.5 },
}

const ESTIMATED_TIMES = {
  import: '~5 seconds', calculation: '~1-5 minutes', generation: '~5-10 minutes',
  docking: '~3-5 min', admet: '~30 sec', scoring: '~15 sec',
  enrichment: '~1-2 min', clustering: '~30 sec', off_target: '~1 min',
  confidence: '~15 sec', retrosynthesis: '~1 min', safety: '~30 sec',
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

function EngineCard({ id, label, badges = [], selected, onClick }) {
  return (
    <button
      onClick={() => onClick(id)}
      className={`w-full flex items-start gap-3 p-3 rounded-xl border-2 text-left transition-all ${
        selected
          ? 'border-bx-mint bg-blue-50 ring-2 ring-bx-mint/20'
          : 'border-gray-100 hover:border-blue-200 hover:shadow-sm'
      }`}
    >
      <div className={`w-2.5 h-2.5 rounded-full mt-1.5 flex-shrink-0 border-2 transition-colors ${
        selected ? 'border-bx-mint bg-bx-surface' : 'border-gray-300'
      }`} />
      <div className="flex-1 min-w-0">
        <p className={`text-sm font-semibold ${selected ? 'text-bx-light-text' : 'text-gray-800'}`}>{label}</p>
        <div className="flex flex-wrap gap-1 mt-1">
          {badges.map(b => (
            <span key={b.text} className={`text-[10px] px-1.5 py-0.5 rounded-full font-semibold ${b.color}`}>
              {b.text}
            </span>
          ))}
        </div>
      </div>
    </button>
  )
}

function ConfigForm({ runType, config, onChange, phase }) {
  const set = (key, value) => onChange({ ...config, [key]: value })

  const toggleArr = (key, item) => {
    const arr = config[key] || []
    onChange({ ...config, [key]: arr.includes(item) ? arr.filter(v => v !== item) : [...arr, item] })
  }

  const fileInputRef = useRef(null)

  switch (runType) {
    case 'import':
      return (
        <div className="space-y-4">
          <FormSection title="Source">
            <div className="flex gap-2">
              {[
                { value: 'external', label: 'External file' },
                { value: 'internal', label: 'Internal selection' },
              ].map(opt => (
                <label key={opt.value}
                  className={`flex-1 flex items-center gap-2 p-3 rounded-xl border-2 cursor-pointer transition-all ${
                    config.sourceMode === opt.value
                      ? 'border-bx-mint bg-blue-50'
                      : 'border-gray-100 hover:border-gray-200'
                  }`}
                >
                  <input
                    type="radio"
                    checked={config.sourceMode === opt.value}
                    onChange={() => set('sourceMode', opt.value)}
                    className="accent-bx-mint"
                  />
                  <span className="text-sm font-medium text-gray-700">{opt.label}</span>
                </label>
              ))}
            </div>
          </FormSection>

          {config.sourceMode === 'external' ? (
            <>
              <div
                onClick={() => fileInputRef.current?.click()}
                className="border-2 border-dashed border-gray-200 rounded-xl p-6 text-center cursor-pointer
                           hover:border-bx-mint/40 hover:bg-blue-50/30 transition-all group"
              >
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
                <input ref={fileInputRef} type="file" accept=".sdf,.smi,.csv" className="hidden" />
              </div>
              <div>
                <label className="block text-sm font-medium text-gray-500 mb-1">Source name (optional)</label>
                <input
                  type="text"
                  value={config.sourceName || ''}
                  onChange={e => set('sourceName', e.target.value)}
                  placeholder="e.g. ChEMBL EGFR batch 1"
                  className="w-full text-sm border border-gray-200 rounded-lg px-3 py-2 bg-white
                             focus:outline-none focus:ring-2 focus:ring-bx-mint/30 focus:border-bx-mint"
                />
              </div>
            </>
          ) : (
            <div className="bg-blue-50 border border-blue-100 rounded-xl p-4 space-y-2">
              <div className="flex items-center gap-2">
                <svg className="w-4 h-4 text-blue-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                    d="M5 5a2 2 0 012-2h10a2 2 0 012 2v16l-7-3.5L5 21V5z" />
                </svg>
                <p className="text-sm font-semibold text-blue-700">From previous phase (bookmarked)</p>
              </div>
              <p className="text-sm text-blue-500">
                Imports all bookmarked molecules from the previous phase into this one.
              </p>
            </div>
          )}
        </div>
      )

    case 'calculation': {
      const calcTypes = config.calculation_types || []
      const toggleCalcType = (key) => {
        const updated = calcTypes.includes(key)
          ? calcTypes.filter(k => k !== key)
          : [...calcTypes, key]
        onChange({ ...config, calculation_types: updated })
      }
      return (
        <div className="space-y-4">
          <FormSection title="Select calculations to run">
            <div className="grid grid-cols-1 gap-2">
              {CALCULATION_SUBTYPES.map(sub => {
                const checked = calcTypes.includes(sub.key)
                const subColor = RUN_TYPE_COLORS[sub.key] || RUN_TYPE_COLORS.calculation
                return (
                  <label key={sub.key}
                    className={`flex items-start gap-3 p-3 rounded-xl border-2 cursor-pointer transition-all ${
                      checked ? `border-blue-400 ${subColor.bg}` : 'border-gray-100 hover:border-gray-200'
                    }`}
                  >
                    <input type="checkbox" checked={checked} onChange={() => toggleCalcType(sub.key)}
                      className="accent-bx-mint mt-0.5" />
                    <div className="flex-1 min-w-0">
                      <div className="flex items-center gap-2">
                        <p className={`text-sm font-semibold ${checked ? 'text-bx-light-text' : 'text-gray-700'}`}>{sub.label}</p>
                        <span className="text-[9px] text-gray-400">{ESTIMATED_TIMES[sub.key]}</span>
                      </div>
                      <p className="text-sm text-gray-500 mt-0.5">{sub.description}</p>
                      <div className="flex flex-wrap gap-1 mt-1">
                        {sub.columns.slice(0, 4).map(c => (
                          <span key={c} className="text-[9px] px-1.5 py-0.5 rounded bg-white/70 text-gray-500 border border-gray-200">{c}</span>
                        ))}
                        {sub.columns.length > 4 && <span className="text-[9px] text-gray-400">+{sub.columns.length - 4}</span>}
                      </div>
                    </div>
                  </label>
                )
              })}
            </div>
          </FormSection>

          {/* Per-subtype config panels */}
          {calcTypes.includes('docking') && (
            <FormSection title="Docking configuration">
              <div className="space-y-2 bg-blue-50/50 rounded-xl p-3">
                <EngineCard id="gnina_gpu" label="GNINA (GPU)"
                  badges={[{ text: 'Fastest', color: 'bg-green-100 text-green-700' }]}
                  selected={(config.docking?.engine || 'gnina_gpu') === 'gnina_gpu'}
                  onClick={v => onChange({ ...config, docking: { ...config.docking, engine: v } })} />
                <EngineCard id="gnina_cpu" label="GNINA (CPU)"
                  badges={[{ text: 'CNN rescoring', color: 'bg-blue-100 text-blue-700' }]}
                  selected={(config.docking?.engine) === 'gnina_cpu'}
                  onClick={v => onChange({ ...config, docking: { ...config.docking, engine: v } })} />
                <EngineCard id="vina" label="AutoDock Vina"
                  badges={[{ text: 'Classic', color: 'bg-gray-100 text-gray-600' }]}
                  selected={(config.docking?.engine) === 'vina'}
                  onClick={v => onChange({ ...config, docking: { ...config.docking, engine: v } })} />
                <div className="flex items-center gap-4 mt-2">
                  <div>
                    <label className="text-[10px] font-medium text-gray-500">Exhaustiveness</label>
                    <Stepper value={config.docking?.exhaustiveness ?? 32} min={8} max={64}
                      onChange={v => onChange({ ...config, docking: { ...config.docking, exhaustiveness: v } })} />
                  </div>
                </div>
              </div>
            </FormSection>
          )}

          {calcTypes.includes('scoring') && (
            <FormSection title="Scoring weights">
              <p className="text-sm text-gray-400">Configure in the scoring weights editor when run starts.</p>
            </FormSection>
          )}

          {calcTypes.length === 0 && (
            <div className="bg-amber-50 border border-amber-100 rounded-xl px-4 py-3 text-sm text-amber-700">
              Select at least one calculation type to proceed.
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
                selected={config.engine === 'gnina_gpu'} onClick={v => set('engine', v)} />
              <EngineCard id="gnina_cpu" label="GNINA (CPU)"
                badges={[{ text: 'CNN rescoring', color: 'bg-blue-100 text-blue-700' }]}
                selected={config.engine === 'gnina_cpu'} onClick={v => set('engine', v)} />
              <EngineCard id="vina" label="AutoDock Vina"
                badges={[{ text: 'Classic', color: 'bg-gray-100 text-gray-600' }]}
                selected={config.engine === 'vina'} onClick={v => set('engine', v)} />
            </div>
          </FormSection>

          <FormSection title="Search Parameters">
            <div className="space-y-3">
              <div>
                <div className="flex items-center justify-between text-sm mb-1">
                  <label className="font-medium text-gray-600">Exhaustiveness</label>
                  <span className="font-bold text-bx-light-text tabular-nums">{config.exhaustiveness ?? 32}</span>
                </div>
                <input
                  type="range" min={8} max={64} step={8}
                  value={config.exhaustiveness ?? 32}
                  onChange={e => set('exhaustiveness', Number(e.target.value))}
                  className="w-full accent-bx-mint"
                />
                <div className="flex justify-between text-[10px] text-gray-400 mt-0.5">
                  {[8, 16, 32, 64].map(v => <span key={v}>{v}</span>)}
                </div>
              </div>

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

    case 'admet': {
      const allProps = [
        { key: 'logP',                 label: 'LogP',                  desc: 'Lipophilicity' },
        { key: 'solubility',           label: 'Solubility',            desc: 'Aqueous solubility' },
        { key: 'BBB',                  label: 'BBB Penetration',       desc: 'Blood-brain barrier' },
        { key: 'hERG',                 label: 'hERG Risk',             desc: 'Cardiotoxicity risk' },
        { key: 'metabolic_stability',  label: 'Metabolic Stability',   desc: 'CYP stability' },
        { key: 'CYP',                  label: 'CYP Inhibition',        desc: '5 isoforms' },
        { key: 'oral_bioavailability', label: 'Oral Bioavailability',  desc: 'F% prediction' },
        { key: 'ames',                 label: 'Ames Mutagenicity',     desc: 'Genotoxicity flag' },
      ]
      const all = config.properties || []
      const allSelected = allProps.every(p => all.includes(p.key))
      return (
        <div className="space-y-3">
          <div className="flex items-center justify-between">
            <label className="text-sm font-semibold text-gray-500 uppercase tracking-wide">Properties</label>
            <div className="flex gap-2">
              <button
                onClick={() => onChange({ ...config, properties: allProps.map(p => p.key) })}
                className={`text-sm px-2 py-1 rounded font-medium transition-colors ${
                  allSelected ? 'bg-bx-surface text-white' : 'bg-gray-100 text-gray-600 hover:bg-gray-200'
                }`}
              >
                All
              </button>
              <button
                onClick={() => onChange({ ...config, properties: [] })}
                className="text-sm px-2 py-1 rounded font-medium bg-gray-100 text-gray-600 hover:bg-gray-200 transition-colors"
              >
                None
              </button>
            </div>
          </div>
          <div className="grid grid-cols-2 gap-1.5">
            {allProps.map(prop => {
              const checked = all.includes(prop.key)
              return (
                <label key={prop.key}
                  className={`flex items-start gap-2 p-2.5 rounded-xl cursor-pointer border transition-all ${
                    checked
                      ? 'border-green-200 bg-green-50'
                      : 'border-gray-100 hover:border-gray-200 hover:bg-gray-50'
                  }`}
                >
                  <input
                    type="checkbox"
                    checked={checked}
                    onChange={() => toggleArr('properties', prop.key)}
                    className="accent-green-600 mt-0.5 flex-shrink-0"
                  />
                  <div className="min-w-0">
                    <p className="text-sm font-semibold text-gray-700">{prop.label}</p>
                    <p className="text-[10px] text-gray-400">{prop.desc}</p>
                  </div>
                </label>
              )
            })}
          </div>
        </div>
      )
    }

    case 'scoring': {
      const weights = config.weights || DEFAULT_CONFIGS.scoring.weights
      const total = Object.values(weights).reduce((a, b) => a + b, 0)
      const totalOk = Math.abs(total - 1.0) < 0.01
      return (
        <div className="space-y-4">
          <FormSection title="Score Weights">
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
              onClick={() => onChange({ ...config, weights: DEFAULT_CONFIGS.scoring.weights })}
              className="text-sm px-3 py-1.5 rounded-lg bg-gray-100 text-gray-600 hover:bg-gray-200 transition-colors font-medium"
            >
              Reset defaults
            </button>
          </div>
        </div>
      )
    }

    case 'enrichment': {
      const allAnalyses = [
        { value: 'prolif',        label: 'ProLIF Interactions',  desc: 'Protein-ligand interaction fingerprints' },
        { value: 'clustering',    label: 'Tanimoto Clustering',  desc: 'Group molecules by structural similarity' },
        { value: 'scaffold',      label: 'Scaffold Analysis',    desc: 'Identify core scaffolds across hits' },
        { value: 'pharmacophore', label: 'Pharmacophore',        desc: 'Map 3D pharmacophoric features' },
      ]
      return (
        <div className="space-y-2">
          <label className="block text-sm font-semibold text-gray-500 uppercase tracking-wide mb-3">
            Analyses to run
          </label>
          {allAnalyses.map(a => {
            const checked = (config.analyses || []).includes(a.value)
            return (
              <label key={a.value}
                className={`flex items-start gap-3 p-3 rounded-xl cursor-pointer border transition-all ${
                  checked
                    ? 'border-purple-200 bg-purple-50'
                    : 'border-gray-100 hover:border-gray-200 hover:bg-gray-50'
                }`}
              >
                <input
                  type="checkbox"
                  checked={checked}
                  onChange={() => toggleArr('analyses', a.value)}
                  className="accent-purple-600 mt-0.5 flex-shrink-0"
                />
                <div>
                  <p className="text-sm font-semibold text-gray-700">{a.label}</p>
                  <p className="text-sm text-gray-400">{a.desc}</p>
                </div>
              </label>
            )
          })}
        </div>
      )
    }

    case 'generation': {
      const methods = [
        { value: 'scaffold_hopping', label: 'Scaffold Hopping',  desc: 'Replace core scaffold while keeping substituents' },
        { value: 'fragment_growing', label: 'Fragment Growing',  desc: 'Extend molecule at reactive positions' },
        { value: 'de_novo',          label: 'De Novo SMILES',    desc: 'Generate novel structures from scratch' },
      ]
      return (
        <div className="space-y-5">
          <FormSection title="Method">
            <div className="space-y-2">
              {methods.map(m => (
                <button
                  key={m.value}
                  onClick={() => set('method', m.value)}
                  className={`w-full flex items-start gap-3 p-3 rounded-xl border-2 text-left transition-all ${
                    config.method === m.value
                      ? 'border-pink-400 bg-pink-50'
                      : 'border-gray-100 hover:border-pink-200 hover:bg-pink-50/20'
                  }`}
                >
                  <div className={`w-2.5 h-2.5 rounded-full mt-1.5 border-2 flex-shrink-0 ${
                    config.method === m.value ? 'border-pink-500 bg-pink-500' : 'border-gray-300'
                  }`} />
                  <div>
                    <p className={`text-sm font-semibold ${config.method === m.value ? 'text-pink-700' : 'text-gray-800'}`}>
                      {m.label}
                    </p>
                    <p className="text-sm text-gray-400 mt-0.5">{m.desc}</p>
                  </div>
                </button>
              ))}
            </div>
          </FormSection>

          <FormSection title="Generation Parameters">
            <div className="space-y-4">
              <div className="flex items-center justify-between">
                <div>
                  <p className="text-sm font-medium text-gray-600">Iterations</p>
                  <p className="text-[10px] text-gray-400">rounds of generation</p>
                </div>
                <Stepper value={config.iterations ?? 3} min={1} max={5}
                  onChange={v => set('iterations', v)} />
              </div>
              <div className="flex items-center justify-between">
                <div>
                  <p className="text-sm font-medium text-gray-600">Variants per iteration</p>
                  <p className="text-[10px] text-gray-400">new molecules per round</p>
                </div>
                <Stepper value={config.variants_per_iteration ?? 5} min={3} max={20}
                  onChange={v => set('variants_per_iteration', v)} />
              </div>
            </div>
          </FormSection>

          <FormSection title="Include analyses">
            <div className="flex flex-wrap gap-2">
              {[
                { key: 'include_docking', label: 'Docking' },
                { key: 'include_admet',   label: 'ADMET' },
                { key: 'include_scoring', label: 'Scoring' },
              ].map(item => (
                <label key={item.key}
                  className={`flex items-center gap-1.5 px-3 py-1.5 rounded-full border cursor-pointer text-sm font-semibold transition-all ${
                    config[item.key]
                      ? 'border-pink-300 bg-pink-50 text-pink-700'
                      : 'border-gray-200 text-gray-500 hover:border-gray-300'
                  }`}
                >
                  <input
                    type="checkbox"
                    checked={!!config[item.key]}
                    onChange={() => set(item.key, !config[item.key])}
                    className="accent-pink-500"
                  />
                  {item.label}
                </label>
              ))}
            </div>
          </FormSection>
        </div>
      )
    }

    case 'clustering': {
      const methods = [
        { value: 'butina',   label: 'Butina',   desc: 'Sphere exclusion — fast, deterministic' },
        { value: 'tanimoto', label: 'Tanimoto', desc: 'Pairwise similarity matrix' },
        { value: 'ecfp4',    label: 'ECFP4',    desc: 'Morgan fingerprint based' },
      ]
      const cutoff = config.cutoff ?? 0.5
      const estClusters = Math.max(1, Math.round(20 * (1 - cutoff) * 2))
      return (
        <div className="space-y-5">
          <FormSection title="Clustering Method">
            <div className="space-y-2">
              {methods.map(m => (
                <button
                  key={m.value}
                  onClick={() => set('method', m.value)}
                  className={`w-full flex items-start gap-3 p-3 rounded-xl border-2 text-left transition-all ${
                    config.method === m.value
                      ? 'border-teal-400 bg-teal-50'
                      : 'border-gray-100 hover:border-teal-200'
                  }`}
                >
                  <div className={`w-2.5 h-2.5 rounded-full mt-1.5 border-2 flex-shrink-0 ${
                    config.method === m.value ? 'border-teal-500 bg-teal-500' : 'border-gray-300'
                  }`} />
                  <div>
                    <p className={`text-sm font-semibold ${config.method === m.value ? 'text-teal-700' : 'text-gray-800'}`}>
                      {m.label}
                    </p>
                    <p className="text-sm text-gray-400">{m.desc}</p>
                  </div>
                </button>
              ))}
            </div>
          </FormSection>

          <FormSection title="Similarity Cutoff">
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
                <span>0.30 (loose)</span>
                <span>0.80 (strict)</span>
              </div>
              <div className="mt-2 bg-teal-50 border border-teal-100 rounded-lg px-3 py-2 text-sm text-teal-700">
                Estimated clusters: <span className="font-bold">~{estClusters}</span>
              </div>
            </div>
          </FormSection>
        </div>
      )
    }

    default:
      return <p className="text-sm text-gray-400">No configuration needed for this run type.</p>
  }
}

// ---------------------------------------------------------------------------
// Step 3: Confirmation
// ---------------------------------------------------------------------------
function ConfirmationView({ runType, config, selectedCount }) {
  const rt = RUN_TYPES.find(t => t.type === runType)
  if (!rt) return null
  const c = runColor(runType)

  const configLines = []
  if (runType === 'import') {
    configLines.push(`Source: ${config.sourceMode === 'internal' ? 'Phase bookmarks' : 'External file'}`)
    if (config.sourceName) configLines.push(`Label: ${config.sourceName}`)
  }
  if (runType === 'calculation') {
    const calcTypes = config.calculation_types || []
    const labels = calcTypes.map(k => CALCULATION_SUBTYPES.find(s => s.key === k)?.label || k)
    configLines.push(`Calculations: ${labels.join(', ')}`)
    if (calcTypes.includes('docking')) {
      const dk = config.docking || {}
      configLines.push(`Docking engine: ${dk.engine === 'gnina_gpu' ? 'GNINA (GPU)' : dk.engine === 'gnina_cpu' ? 'GNINA (CPU)' : dk.engine || 'GNINA (GPU)'}`)
    }
  }
  if (runType === 'docking') {
    configLines.push(`Engine: ${config.engine === 'gnina_gpu' ? 'GNINA (GPU)' : config.engine === 'gnina_cpu' ? 'GNINA (CPU)' : 'AutoDock Vina'}`)
    configLines.push(`Exhaustiveness: ${config.exhaustiveness ?? 32}  |  Poses: ${config.num_modes ?? 9}`)
    const boxX = config.boxSizeX || 'auto', boxY = config.boxSizeY || 'auto', boxZ = config.boxSizeZ || 'auto'
    configLines.push(`Box: ${boxX} x ${boxY} x ${boxZ} A`)
  }
  if (runType === 'admet') {
    const props = config.properties || []
    configLines.push(`${props.length} properties: ${props.slice(0, 3).join(', ')}${props.length > 3 ? ` +${props.length - 3} more` : ''}`)
  }
  if (runType === 'scoring') {
    const ws = config.weights || {}
    const total = Object.values(ws).reduce((a, b) => a + b, 0)
    configLines.push(`${Object.keys(ws).length} score weights, total: ${total.toFixed(2)}`)
  }
  if (runType === 'enrichment') {
    const a = config.analyses || []
    configLines.push(`Analyses: ${a.join(', ') || 'none selected'}`)
  }
  if (runType === 'generation') {
    configLines.push(`Method: ${config.method}`)
    configLines.push(`${config.iterations ?? 3} iterations x ${config.variants_per_iteration ?? 5} variants = ~${(config.iterations ?? 3) * (config.variants_per_iteration ?? 5)} molecules`)
  }
  if (runType === 'clustering') {
    configLines.push(`Method: ${config.method}, cutoff: ${Number(config.cutoff ?? 0.5).toFixed(2)}`)
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

      {/* Validation warnings */}
      {runType === 'calculation' && (!config.calculation_types || config.calculation_types.length === 0) && (
        <div className="flex items-center gap-2 bg-red-50 border border-red-200 rounded-xl px-4 py-2.5 text-sm text-red-600">
          <svg className="w-4 h-4 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
          </svg>
          No calculation types selected. Go back and select at least one.
        </div>
      )}
      {(runType === 'calculation' || runType === 'generation') && selectedCount === 0 && (
        <div className="flex items-center gap-2 bg-amber-50 border border-amber-200 rounded-xl px-4 py-2.5 text-sm text-amber-700">
          <svg className="w-4 h-4 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
          </svg>
          No molecules selected in the dashboard. Select molecules before launching.
        </div>
      )}

      {/* Molecule selection count */}
      {(runType === 'calculation' || runType === 'generation') && selectedCount > 0 && (
        <div className="flex items-center gap-2 bg-green-50 border border-green-200 rounded-xl px-4 py-2.5 text-sm text-green-700">
          <svg className="w-4 h-4 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M5 13l4 4L19 7" />
          </svg>
          {selectedCount} molecule{selectedCount > 1 ? 's' : ''} selected as input
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
export default function RunCreator({ phaseId, phaseType, isOpen, onClose, onSubmit, selectedMoleculeIds, submitting }) {
  const [step, setStep] = useState(1)
  const [selectedType, setSelectedType] = useState(null)
  const [config, setConfig] = useState({})

  if (!isOpen) return null

  const allRunTypes = RUN_TYPES
  const STEP_LABELS = ['Type', 'Configure', 'Launch']

  function handleSelectType(type) {
    setSelectedType(type)
    setConfig(DEFAULT_CONFIGS[type] || {})
    setStep(2)
  }

  const canSubmit = (() => {
    if (submitting) return false
    if (selectedType === 'calculation') {
      if (!config.calculation_types?.length) return false
      if (!selectedMoleculeIds?.size) return false
    }
    if (selectedType === 'generation' && !selectedMoleculeIds?.size) return false
    return true
  })()

  function handleSubmit() {
    if (!canSubmit) return
    const payload = { type: selectedType, config, phase_id: phaseId }
    // For calculation runs, include calculation_types and selected molecules
    if (selectedType === 'calculation') {
      payload.calculation_types = config.calculation_types || []
      if (selectedMoleculeIds?.size > 0) {
        payload.input_molecule_ids = [...selectedMoleculeIds]
      }
    }
    if (selectedType === 'generation' && selectedMoleculeIds?.size > 0) {
      payload.input_molecule_ids = [...selectedMoleculeIds]
    }
    if (onSubmit) onSubmit(payload)
    handleClose()
  }

  function handleClose() {
    setStep(1)
    setSelectedType(null)
    setConfig({})
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
            />
          )}

          {/* Step 3: Confirm */}
          {step === 3 && selectedType && (
            <ConfirmationView runType={selectedType} config={config} selectedCount={selectedMoleculeIds?.size || 0} />
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
              onClick={() => setStep(3)}
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
