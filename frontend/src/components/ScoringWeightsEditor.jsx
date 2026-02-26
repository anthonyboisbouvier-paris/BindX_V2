import React, { useState, useCallback } from 'react'

// ---------------------------------------------------------------------------
// Default weight labels / display names
// ---------------------------------------------------------------------------

const WEIGHT_META = {
  docking_score: { label: 'Docking Score', color: '#3b82f6' },
  CNNscore: { label: 'CNN Score', color: '#8b5cf6' },
  logP: { label: 'LogP', color: '#06b6d4' },
  solubility: { label: 'Solubility', color: '#22c55e' },
  selectivity: { label: 'Selectivity', color: '#f59e0b' },
  novelty: { label: 'Novelty', color: '#ec4899' },
}

// ---------------------------------------------------------------------------
// Mini radar / pie chart using SVG
// ---------------------------------------------------------------------------

function WeightPieChart({ weights }) {
  const entries = Object.entries(weights).filter(([, v]) => v > 0)
  const total = entries.reduce((s, [, v]) => s + v, 0)
  if (total === 0 || entries.length === 0) return null

  const cx = 40
  const cy = 40
  const r = 32
  let startAngle = -Math.PI / 2

  const slices = entries.map(([key, val]) => {
    const frac = val / total
    const angle = frac * 2 * Math.PI
    const endAngle = startAngle + angle

    const x1 = cx + r * Math.cos(startAngle)
    const y1 = cy + r * Math.sin(startAngle)
    const x2 = cx + r * Math.cos(endAngle)
    const y2 = cy + r * Math.sin(endAngle)
    const largeArc = angle > Math.PI ? 1 : 0

    const slice = { key, x1, y1, x2, y2, cx, cy, r, largeArc, startAngle, endAngle }
    startAngle = endAngle
    return slice
  })

  return (
    <svg width="80" height="80" viewBox="0 0 80 80" className="shrink-0">
      {slices.map(({ key, x1, y1, x2, y2, largeArc }) => (
        <path
          key={key}
          d={`M ${cx} ${cy} L ${x1} ${y1} A ${r} ${r} 0 ${largeArc} 1 ${x2} ${y2} Z`}
          fill={(WEIGHT_META[key] || { color: '#94a3b8' }).color}
          stroke="white"
          strokeWidth="1.5"
        />
      ))}
    </svg>
  )
}

// ---------------------------------------------------------------------------
// Single weight row
// ---------------------------------------------------------------------------

function WeightRow({ name, value, onChange }) {
  const meta = WEIGHT_META[name] || { label: name.replace(/_/g, ' '), color: '#94a3b8' }

  return (
    <div className="flex items-center gap-3 py-1.5">
      <div className="flex items-center gap-2 w-32 shrink-0">
        <span
          className="w-2.5 h-2.5 rounded-full shrink-0"
          style={{ backgroundColor: meta.color }}
        />
        <span className="text-xs text-gray-600 truncate">{meta.label}</span>
      </div>
      <div className="flex-1 relative">
        <input
          type="range"
          min="0"
          max="1"
          step="0.01"
          value={value}
          onChange={e => onChange(name, parseFloat(e.target.value))}
          className="w-full h-1.5 rounded-full appearance-none cursor-pointer"
          style={{
            background: `linear-gradient(to right, ${meta.color} 0%, ${meta.color} ${value * 100}%, #e5e7eb ${value * 100}%, #e5e7eb 100%)`,
            accentColor: meta.color,
          }}
        />
      </div>
      <span className="text-xs font-mono font-semibold text-gray-700 w-10 text-right shrink-0">
        {value.toFixed(2)}
      </span>
    </div>
  )
}

// ---------------------------------------------------------------------------
// ScoringWeightsEditor
// ---------------------------------------------------------------------------
// Props:
//   weights  — { docking_score: 0.3, CNNscore: 0.2, logP: 0.15, ... }
//   onChange — (newWeights) => void
//   onReset  — () => void

const DEFAULT_WEIGHTS = {
  docking_score: 0.30,
  CNNscore: 0.20,
  logP: 0.15,
  solubility: 0.10,
  selectivity: 0.15,
  novelty: 0.10,
}

export default function ScoringWeightsEditor({ weights: externalWeights, onChange, onReset }) {
  const [localWeights, setLocalWeights] = useState(externalWeights || DEFAULT_WEIGHTS)
  const weights = externalWeights || localWeights

  const total = Object.values(weights).reduce((s, v) => s + v, 0)
  const isValid = Math.abs(total - 1.0) < 0.05

  const handleChange = useCallback((name, value) => {
    const updated = { ...weights, [name]: value }
    setLocalWeights(updated)
    if (onChange) onChange(updated)
  }, [weights, onChange])

  const handleReset = useCallback(() => {
    setLocalWeights(DEFAULT_WEIGHTS)
    if (onChange) onChange(DEFAULT_WEIGHTS)
    if (onReset) onReset()
  }, [onChange, onReset])

  return (
    <div className="bg-white rounded-xl border border-gray-100 shadow-sm overflow-hidden">
      <div className="px-5 py-4 border-b border-gray-50">
        <div className="flex items-center justify-between">
          <h4 className="text-sm font-semibold text-gray-700">Scoring Weights</h4>
          <button
            onClick={handleReset}
            className="text-xs text-gray-400 hover:text-gray-600 transition-colors flex items-center gap-1"
          >
            <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
            </svg>
            Reset defaults
          </button>
        </div>
      </div>

      <div className="px-5 py-4">
        <div className="flex gap-5">
          {/* Sliders */}
          <div className="flex-1 space-y-0.5">
            {Object.keys(weights).map(name => (
              <WeightRow
                key={name}
                name={name}
                value={weights[name]}
                onChange={handleChange}
              />
            ))}
          </div>

          {/* Pie chart */}
          <div className="flex flex-col items-center justify-center gap-1 shrink-0">
            <WeightPieChart weights={weights} />
            <span className="text-xs text-gray-400">Distribution</span>
          </div>
        </div>

        {/* Total row */}
        <div className={`mt-4 flex items-center justify-between px-3 py-2 rounded-lg border ${
          isValid
            ? 'bg-green-50 border-green-100'
            : 'bg-red-50 border-red-100'
        }`}>
          <span className="text-xs font-medium text-gray-600">Total</span>
          <div className="flex items-center gap-1.5">
            <span className={`text-sm font-bold font-mono ${isValid ? 'text-green-600' : 'text-red-600'}`}>
              {total.toFixed(2)}
            </span>
            {isValid ? (
              <svg className="w-4 h-4 text-green-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 13l4 4L19 7" />
              </svg>
            ) : (
              <svg className="w-4 h-4 text-red-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                  d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
              </svg>
            )}
          </div>
        </div>
      </div>
    </div>
  )
}
