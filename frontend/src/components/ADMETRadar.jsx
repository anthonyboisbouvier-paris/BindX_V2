import React, { useMemo, useState } from 'react'
import InfoTip from './InfoTip.jsx'

// --------------------------------------------------
// Geometry helpers
// --------------------------------------------------

/**
 * Convert a normalised radius value + axis index to SVG x,y coordinates.
 * Axis 0 starts at the top (−π/2 offset).
 */
function polarToXY(value, axisIndex, totalAxes, cx, cy, radius) {
  const angle = (Math.PI * 2 * axisIndex) / totalAxes - Math.PI / 2
  return {
    x: cx + value * radius * Math.cos(angle),
    y: cy + value * radius * Math.sin(angle),
  }
}

/**
 * Build an SVG polygon `points` string from an array of normalised values (0-1).
 */
function buildPoints(values, cx, cy, radius, totalAxes) {
  return values
    .map((v, i) => {
      const { x, y } = polarToXY(v, i, totalAxes, cx, cy, radius)
      return `${x.toFixed(2)},${y.toFixed(2)}`
    })
    .join(' ')
}

// --------------------------------------------------
// Color helpers
// --------------------------------------------------
const FILL_COLORS = {
  green:  '#00e6a0',
  yellow: '#eab308',
  red:    '#ef4444',
}

function resolveColor(colorCode) {
  if (!colorCode) return FILL_COLORS.yellow
  return FILL_COLORS[colorCode.toLowerCase()] || FILL_COLORS.yellow
}

// --------------------------------------------------
// Extract 6 normalised ADMET axis values (0-1).
// Handles both flat and nested payload shapes.
// --------------------------------------------------
function extractAxisValues(admet) {
  if (!admet) return Array(6).fill(0.5)

  const getFlat   = (key)          => admet[key]
  const getNested = (sec, key)     => admet[sec]?.[key]
  const get       = (fk, sec, nk, fallback = 0.5) => {
    const v = getFlat(fk) ?? getNested(sec, nk)
    return v !== null && v !== undefined ? Number(v) : fallback
  }

  // 1. Absorption — oral bioavailability (0-1, high = good)
  const absorption = get('oral_bioavailability', 'absorption', 'oral_bioavailability', 0.5)

  // 2. Distribution — BBB permeability (0-1)
  const distribution = get('bbb_permeability', 'distribution', 'bbb_permeability', 0.5)

  // 3. Metabolism — average CYP inhibition then INVERTED (low inhibition = high score)
  const cypPairs = [
    ['cyp1a2_inhibition',  'metabolism', 'cyp1a2_inhibition'],
    ['cyp2c9_inhibition',  'metabolism', 'cyp2c9_inhibition'],
    ['cyp2c19_inhibition', 'metabolism', 'cyp2c19_inhibition'],
    ['cyp2d6_inhibition',  'metabolism', 'cyp2d6_inhibition'],
    ['cyp3a4_inhibition',  'metabolism', 'cyp3a4_inhibition'],
  ]
  const cypVals = cypPairs
    .map(([fk, sec, nk]) => getFlat(fk) ?? getNested(sec, nk))
    .filter((v) => v !== null && v !== undefined)
    .map(Number)
  const avgCyp = cypVals.length
    ? cypVals.reduce((a, b) => a + b, 0) / cypVals.length
    : 0.3
  const metabolism = Math.max(0, Math.min(1, 1 - avgCyp))

  // 4. Excretion — renal clearance mapped to favour moderate values
  const rawClear = getFlat('clearance')
    ?? getFlat('renal_clearance')
    ?? getNested('excretion', 'clearance')
    ?? getNested('excretion', 'renal_clearance')
  let excretion
  if (rawClear !== null && rawClear !== undefined) {
    const norm = Math.min(1, Number(rawClear) / 100)
    // Best score around 0.2-0.3 normalised (moderate clearance)
    excretion = Math.max(0, Math.min(1, 1 - Math.abs(norm - 0.25) * 2))
  } else {
    excretion = 0.5
  }

  // 5. Toxicity — invert average of hERG + ames + hepatotox (low risk = high score)
  const toxRaw = [
    getFlat('herg_inhibition')  ?? getNested('toxicity', 'herg_inhibition'),
    getFlat('ames_mutagenicity') ?? getNested('toxicity', 'ames_mutagenicity'),
    getFlat('hepatotoxicity')   ?? getNested('toxicity', 'hepatotoxicity'),
  ].filter((v) => v !== null && v !== undefined).map(Number)
  const avgTox = toxRaw.length
    ? toxRaw.reduce((a, b) => a + b, 0) / toxRaw.length
    : 0.2
  const toxicity = Math.max(0, Math.min(1, 1 - avgTox))

  // 6. Drug-likeness — composite score (0-1)
  const drugLikeness = get('composite_score', null, null, 0.5)

  return [absorption, distribution, metabolism, excretion, toxicity, drugLikeness]
}

// --------------------------------------------------
// Flag badge
// --------------------------------------------------
function FlagBadge({ label, warn }) {
  return (
    <span
      className={`inline-flex items-center gap-1 px-2 py-0.5 rounded-full text-sm font-medium ${
        warn
          ? 'bg-red-100 text-red-700 border border-red-200'
          : 'bg-green-100 text-green-700 border border-green-200'
      }`}
    >
      <span className={`w-1.5 h-1.5 rounded-full flex-shrink-0 ${warn ? 'bg-red-500' : 'bg-green-500'}`} />
      {label}
    </span>
  )
}

// ADMET axis tooltips
const AXIS_TIPS = {
  Absorption: 'Oral bioavailability — how well the molecule is absorbed through the digestive tract (0-100%). Higher is better.',
  Distribution: 'Blood-Brain Barrier permeability — ability to cross into the brain. High values may be undesirable for non-CNS drugs.',
  Metabolism: 'Metabolic stability — low CYP enzyme inhibition means fewer drug interactions. Higher score = safer metabolism.',
  Excretion: 'Renal clearance — moderate excretion rate is optimal; too fast or too slow can be problematic.',
  Toxicity: 'Safety score based on hERG inhibition, Ames mutagenicity, and hepatotoxicity. Higher = lower toxicity risk.',
  'Drug-likeness': 'Overall drug-likeness composite score. Combines multiple pharmacological properties into a single metric.',
}

// --------------------------------------------------
// ADMETRadar — SVG radar chart for 6 ADMET categories
// --------------------------------------------------
export default function ADMETRadar({ admet }) {
  const AXES = [
    { label: 'Absorption',    short: 'Abs' },
    { label: 'Distribution',  short: 'Dist' },
    { label: 'Metabolism',    short: 'Met' },
    { label: 'Excretion',     short: 'Exc' },
    { label: 'Toxicity',      short: 'Tox' },
    { label: 'Drug-likeness', short: 'DL' },
  ]

  const N   = AXES.length
  const CX  = 115
  const CY  = 115
  const R   = 85
  const RINGS = [0.25, 0.5, 0.75, 1.0]

  const values     = useMemo(() => extractAxisValues(admet), [admet])
  const colorCode  = admet?.color_code || 'yellow'
  const fillColor  = resolveColor(colorCode)
  const scoreValue = admet?.composite_score
  const flags      = admet?.flags || []
  const [showDetailedDMPK, setShowDetailedDMPK] = useState(false)

  // Axis tip positions (value = 1)
  const axisTips = AXES.map((_, i) => polarToXY(1, i, N, CX, CY, R))

  // Label positions slightly outside the radar
  const labelPos = AXES.map((axis, i) => {
    const { x, y } = polarToXY(1.26, i, N, CX, CY, R)
    return { x, y, ...axis }
  })

  if (!admet) {
    return (
      <div className="card p-6 flex flex-col items-center justify-center h-64 text-gray-300">
        <svg className="w-10 h-10 mb-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
            d="M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 01-2 2h-2a2 2 0 01-2-2z" />
        </svg>
        <p className="text-sm">ADMET data not available</p>
      </div>
    )
  }

  return (
    <div className="card overflow-hidden">
      {/* Header */}
      <div className="px-4 py-3 bg-bx-surface flex items-center justify-between">
        <h3 className="text-white font-semibold text-sm flex items-center gap-2">
          <svg className="w-4 h-4 text-bx-mint" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
              d="M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 01-2 2h-2a2 2 0 01-2-2z" />
          </svg>
          ADMET Profile
          <InfoTip text="Absorption, Distribution, Metabolism, Excretion, Toxicity — key pharmacokinetic properties for drug candidates." />
        </h3>
        {scoreValue !== undefined && scoreValue !== null && (
          <div className="flex items-center gap-2">
            <span className="w-2.5 h-2.5 rounded-full" style={{ backgroundColor: fillColor }} />
            <span className="text-white text-sm font-bold">{Number(scoreValue).toFixed(3)}</span>
            <span className="text-white/50 text-sm">composite score</span>
          </div>
        )}
      </div>

      {/* Body */}
      <div className="p-4 flex flex-col items-center gap-4">

        {/* SVG radar */}
        <svg
          width="230"
          height="230"
          viewBox="0 0 230 230"
          aria-label="ADMET Radar"
          className="overflow-visible"
        >
          {/* Grid rings */}
          {RINGS.map((ring, ri) => (
            <polygon
              key={ri}
              points={buildPoints(Array(N).fill(ring), CX, CY, R, N)}
              fill="none"
              stroke="#e5e7eb"
              strokeWidth={ri === RINGS.length - 1 ? 1.5 : 1}
              strokeDasharray={ri < RINGS.length - 1 ? '4,4' : undefined}
            />
          ))}

          {/* Ring percentage labels (on axis 0, going up) */}
          {[0.25, 0.5, 0.75].map((ring) => {
            const { x, y } = polarToXY(ring, 0, N, CX, CY, R)
            return (
              <text key={ring} x={x + 4} y={y - 2} fontSize="7" fill="#9ca3af" textAnchor="start">
                {Math.round(ring * 100)}%
              </text>
            )
          })}

          {/* Axis spokes */}
          {axisTips.map((tip, i) => (
            <line
              key={i}
              x1={CX} y1={CY}
              x2={tip.x} y2={tip.y}
              stroke="#d1d5db"
              strokeWidth={1}
            />
          ))}

          {/* Data polygon */}
          <polygon
            points={buildPoints(values, CX, CY, R, N)}
            fill={fillColor}
            fillOpacity={0.18}
            stroke={fillColor}
            strokeWidth={2}
            strokeLinejoin="round"
          />

          {/* Data dots */}
          {values.map((v, i) => {
            const { x, y } = polarToXY(v, i, N, CX, CY, R)
            return (
              <circle
                key={i}
                cx={x} cy={y} r={4.5}
                fill={fillColor}
                stroke="white"
                strokeWidth={1.5}
              />
            )
          })}

          {/* Axis labels */}
          {labelPos.map((lp, i) => {
            const anchor = lp.x < CX - 6 ? 'end' : lp.x > CX + 6 ? 'start' : 'middle'
            return (
              <g key={i}>
                <text
                  x={lp.x} y={lp.y - 5}
                  fontSize="9" fontWeight="600"
                  fill="#0f131d"
                  textAnchor={anchor}
                  dominantBaseline="central"
                >
                  {lp.short}
                </text>
                <text
                  x={lp.x} y={lp.y + 7}
                  fontSize="7.5"
                  fill="#6b7280"
                  textAnchor={anchor}
                  dominantBaseline="central"
                >
                  {Math.round(values[i] * 100)}%
                </text>
              </g>
            )
          })}
        </svg>

        {/* Axes legend grid */}
        <div className="grid grid-cols-3 gap-x-6 gap-y-1.5 w-full px-2">
          {AXES.map((axis, i) => (
            <div key={i} className="flex items-center gap-1.5">
              <span
                className="w-2 h-2 rounded-full flex-shrink-0"
                style={{ backgroundColor: fillColor }}
              />
              <span className="text-sm text-gray-500 truncate flex items-center">
                {axis.label}
                <InfoTip text={AXIS_TIPS[axis.label] || axis.label} />
              </span>
              <span className="text-sm font-semibold text-gray-700 ml-auto">
                {Math.round(values[i] * 100)}%
              </span>
            </div>
          ))}
        </div>

        {/* ADMET flags */}
        {flags.length > 0 && (
          <div className="w-full border-t border-gray-100 pt-3">
            <p className="text-sm font-medium text-gray-400 uppercase tracking-wide mb-2">
              ADMET Alerts
            </p>
            <div className="flex flex-wrap gap-1.5">
              {flags.map((flag, i) => {
                const rawLabel = typeof flag === 'string' ? flag : (flag?.label || flag?.message || String(flag))
                const isWarn = typeof flag === 'string'
                  ? /alert|warn|risk|fail|herg|toxic|mutagenic/i.test(flag)
                  : flag?.type === 'warning' || flag?.type === 'alert'
                return <FlagBadge key={i} label={rawLabel} warn={isWarn} />
              })}
            </div>
          </div>
        )}

        {/* Detailed DMPK table (collapsible) */}
        <div className="w-full border-t border-gray-100 pt-3">
          <button
            onClick={() => setShowDetailedDMPK(v => !v)}
            className="flex items-center gap-1.5 text-sm font-medium text-gray-500 hover:text-bx-light-text transition-colors"
          >
            <svg className={`w-3.5 h-3.5 transition-transform ${showDetailedDMPK ? 'rotate-90' : ''}`} fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
            </svg>
            {showDetailedDMPK ? 'Hide' : 'Show'} detailed DMPK values
          </button>
          {showDetailedDMPK && (
            <div className="mt-2 overflow-hidden rounded-lg border border-gray-200">
              <table className="w-full text-sm">
                <thead>
                  <tr className="bg-gray-50 text-left">
                    <th className="px-3 py-1.5 font-semibold text-gray-600">Property</th>
                    <th className="px-3 py-1.5 font-semibold text-gray-600 text-right">Value</th>
                    <th className="px-3 py-1.5 font-semibold text-gray-600 text-center">Status</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-gray-100">
                  {[
                    { key: 'oral_bioavailability', label: 'Oral Bioavailability', section: 'absorption', format: v => `${(v * 100).toFixed(0)}%`, thresholds: [0.3, 0.7] },
                    { key: 'intestinal_permeability', label: 'Intestinal Permeability', section: 'absorption', format: v => Number(v).toFixed(2), thresholds: [0.3, 0.7] },
                    { key: 'ppb', label: 'Plasma Protein Binding', section: 'distribution', format: v => `${(v * 100).toFixed(0)}%`, thresholds: [0.5, 0.95], invert: true },
                    { key: 'bbb_permeability', label: 'BBB Permeability', section: 'distribution', format: v => Number(v).toFixed(2), thresholds: [0.3, 0.7] },
                    { key: 'clearance', label: 'Clearance', section: 'excretion', format: v => `${Number(v).toFixed(1)} mL/min/kg`, thresholds: [5, 30], invert: true },
                    { key: 'half_life', label: 'Half-life', section: 'excretion', format: v => `${Number(v).toFixed(1)} h`, thresholds: [2, 8] },
                    { key: 'cyp1a2_inhibition', label: 'CYP1A2 Inhibitor', section: 'metabolism', format: v => Number(v).toFixed(2), thresholds: [0.5, 0.3], invert: true },
                    { key: 'cyp2c9_inhibition', label: 'CYP2C9 Inhibitor', section: 'metabolism', format: v => Number(v).toFixed(2), thresholds: [0.5, 0.3], invert: true },
                    { key: 'cyp2c19_inhibition', label: 'CYP2C19 Inhibitor', section: 'metabolism', format: v => Number(v).toFixed(2), thresholds: [0.5, 0.3], invert: true },
                    { key: 'cyp2d6_inhibition', label: 'CYP2D6 Inhibitor', section: 'metabolism', format: v => Number(v).toFixed(2), thresholds: [0.5, 0.3], invert: true },
                    { key: 'cyp3a4_inhibition', label: 'CYP3A4 Inhibitor', section: 'metabolism', format: v => Number(v).toFixed(2), thresholds: [0.5, 0.3], invert: true },
                    { key: 'herg_inhibition', label: 'hERG Inhibition', section: 'toxicity', format: v => Number(v).toFixed(2), thresholds: [0.5, 0.3], invert: true },
                    { key: 'ames_mutagenicity', label: 'Ames Mutagenicity', section: 'toxicity', format: v => Number(v).toFixed(2), thresholds: [0.5, 0.3], invert: true },
                    { key: 'hepatotoxicity', label: 'Hepatotoxicity', section: 'toxicity', format: v => Number(v).toFixed(2), thresholds: [0.5, 0.3], invert: true },
                  ].map(({ key, label, section, format, thresholds, invert }) => {
                    const val = admet[key] ?? admet[section]?.[key]
                    if (val === null || val === undefined) return null
                    const numVal = Number(val)
                    let status, statusColor, statusBg
                    if (invert) {
                      // Lower is better (e.g. CYP inhibition, toxicity)
                      if (numVal <= thresholds[1]) { status = 'Low'; statusColor = 'text-green-700'; statusBg = 'bg-green-100' }
                      else if (numVal <= thresholds[0]) { status = 'Moderate'; statusColor = 'text-yellow-700'; statusBg = 'bg-yellow-100' }
                      else { status = 'High'; statusColor = 'text-red-700'; statusBg = 'bg-red-100' }
                    } else {
                      // Higher is better (e.g. bioavailability)
                      if (numVal >= thresholds[1]) { status = 'Good'; statusColor = 'text-green-700'; statusBg = 'bg-green-100' }
                      else if (numVal >= thresholds[0]) { status = 'Moderate'; statusColor = 'text-yellow-700'; statusBg = 'bg-yellow-100' }
                      else { status = 'Low'; statusColor = 'text-red-700'; statusBg = 'bg-red-100' }
                    }
                    return (
                      <tr key={key} className="hover:bg-gray-50">
                        <td className="px-3 py-1.5 text-gray-600">{label}</td>
                        <td className="px-3 py-1.5 text-right font-mono text-gray-800">{format(numVal)}</td>
                        <td className="px-3 py-1.5 text-center">
                          <span className={`inline-block px-1.5 py-0.5 rounded text-[10px] font-semibold ${statusColor} ${statusBg}`}>
                            {status}
                          </span>
                        </td>
                      </tr>
                    )
                  }).filter(Boolean)}
                </tbody>
              </table>
              {/* Show message if no detailed data */}
              {[
                'oral_bioavailability', 'intestinal_permeability', 'ppb', 'bbb_permeability',
                'clearance', 'half_life', 'cyp1a2_inhibition', 'cyp2c9_inhibition',
                'cyp2c19_inhibition', 'cyp2d6_inhibition', 'cyp3a4_inhibition',
                'herg_inhibition', 'ames_mutagenicity', 'hepatotoxicity',
              ].every(k => (admet[k] ?? admet?.absorption?.[k] ?? admet?.distribution?.[k] ?? admet?.metabolism?.[k] ?? admet?.excretion?.[k] ?? admet?.toxicity?.[k]) == null) && (
                <div className="px-3 py-4 text-center text-sm text-gray-400">
                  No detailed DMPK values available for this molecule
                </div>
              )}
            </div>
          )}
        </div>

        {/* Color-coded summary band */}
        <div
          className="w-full rounded-lg px-3 py-2 flex items-center justify-between"
          style={{ backgroundColor: fillColor + '1a', borderLeft: `3px solid ${fillColor}` }}
        >
          <div className="flex items-center gap-2">
            <span className="w-3 h-3 rounded-full" style={{ backgroundColor: fillColor }} />
            <span className="text-sm font-semibold" style={{ color: fillColor }}>
              {colorCode.charAt(0).toUpperCase() + colorCode.slice(1).toLowerCase()}
            </span>
          </div>
          <span className="text-sm text-gray-500">
            {colorCode.toLowerCase() === 'green'  && 'Favorable ADMET profile'}
            {colorCode.toLowerCase() === 'yellow' && 'Acceptable ADMET profile'}
            {colorCode.toLowerCase() === 'red'    && 'ADMET issues detected'}
          </span>
        </div>

        {/* Legend */}
        <div className="flex items-center gap-4 text-sm text-gray-400">
          <span className="flex items-center gap-1">
            <span className="w-2 h-2 rounded-full bg-green-500 inline-block" />
            Favorable
          </span>
          <span className="flex items-center gap-1">
            <span className="w-2 h-2 rounded-full bg-yellow-400 inline-block" />
            Acceptable
          </span>
          <span className="flex items-center gap-1">
            <span className="w-2 h-2 rounded-full bg-red-500 inline-block" />
            Unfavorable
          </span>
        </div>
      </div>
    </div>
  )
}
