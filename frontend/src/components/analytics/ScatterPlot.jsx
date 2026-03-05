import React, { useMemo } from 'react'
import { ScatterChart, Scatter, XAxis, YAxis, Tooltip, CartesianGrid, ZAxis } from 'recharts'
import { ALL_COLUMNS } from '../../lib/columns.js'

/** Format tick values: max 3 decimal places, strip trailing zeros */
const fmtTick = v => typeof v !== 'number' ? v : Number.isInteger(v) ? v : parseFloat(v.toFixed(3))

/** Red (#ef4444) → Yellow (#eab308) → Green (#22c55e) gradient based on normalized 0→1 value */
function colorScale(norm) {
  if (norm <= 0.5) {
    const t = norm * 2
    const r = Math.round(239 + (234 - 239) * t)
    const g = Math.round(68 + (179 - 68) * t)
    const b = Math.round(68 + (8 - 68) * t)
    return `rgb(${r},${g},${b})`
  } else {
    const t = (norm - 0.5) * 2
    const r = Math.round(234 + (34 - 234) * t)
    const g = Math.round(179 + (197 - 179) * t)
    const b = Math.round(8 + (94 - 8) * t)
    return `rgb(${r},${g},${b})`
  }
}

export default function ScatterPlot({ config, molecules, width = 400, height = 220 }) {
  const { xKey, yKey, colorKey, scaleMode, xDomainMin, xDomainMax, yDomainMin, yDomainMax } = config
  const xDef = ALL_COLUMNS.find(c => c.key === xKey)
  const yDef = ALL_COLUMNS.find(c => c.key === yKey)
  const colorDef = colorKey ? ALL_COLUMNS.find(c => c.key === colorKey) : null
  const xLabel = xDef?.label || xKey
  const yLabel = yDef?.label || yKey
  const colorLabel = colorDef?.label || colorKey
  const lowerBetter = colorDef?.colorScale === 'lower-better'

  const { data, cMin, cRange } = useMemo(() => {
    const mapped = molecules
      .map(m => ({
        x: m[xKey], y: m[yKey], name: m.name || '',
        c: colorKey ? m[colorKey] : undefined,
      }))
      .filter(d => d.x != null && d.y != null && typeof d.x === 'number' && typeof d.y === 'number')

    if (!colorKey) return { data: mapped, cMin: 0, cRange: 1 }
    const cVals = mapped.filter(d => d.c != null && typeof d.c === 'number').map(d => d.c)
    const cMin = cVals.length ? Math.min(...cVals) : 0
    const cMax = cVals.length ? Math.max(...cVals) : 1
    return { data: mapped, cMin, cRange: cMax - cMin || 1 }
  }, [molecules, xKey, yKey, colorKey])

  // Compute domains — must be before any early return (hooks rule)
  const xDomain = useMemo(() => {
    if (!data.length) return ['auto', 'auto']
    if (scaleMode === 'manual' && xDomainMin != null && xDomainMax != null) return [xDomainMin, xDomainMax]
    const xs = data.map(d => d.x)
    const min = Math.min(...xs), max = Math.max(...xs)
    const pad = (max - min) * 0.05 || 0.5
    return [
      scaleMode === 'manual' && xDomainMin != null ? xDomainMin : min - pad,
      scaleMode === 'manual' && xDomainMax != null ? xDomainMax : max + pad,
    ]
  }, [data, scaleMode, xDomainMin, xDomainMax])

  const yDomain = useMemo(() => {
    if (!data.length) return ['auto', 'auto']
    if (scaleMode === 'manual' && yDomainMin != null && yDomainMax != null) return [yDomainMin, yDomainMax]
    const ys = data.map(d => d.y)
    const min = Math.min(...ys), max = Math.max(...ys)
    const pad = (max - min) * 0.05 || 0.5
    return [
      scaleMode === 'manual' && yDomainMin != null ? yDomainMin : min - pad,
      scaleMode === 'manual' && yDomainMax != null ? yDomainMax : max + pad,
    ]
  }, [data, scaleMode, yDomainMin, yDomainMax])

  if (!data.length) return <div className="text-gray-400 text-xs text-center py-6">No data for {xLabel} vs {yLabel}</div>

  const hasColor = colorKey && data.some(d => d.c != null)
  const legendW = hasColor ? 50 : 0

  return (
    <div style={{ width, position: 'relative' }}>
      <ScatterChart width={width - legendW} height={height} margin={{ top: 10, right: 10, left: 0, bottom: 20 }}>
        <CartesianGrid strokeDasharray="3 3" stroke="#f0f0f0" />
        <XAxis dataKey="x" type="number" name={xLabel} tick={{ fontSize: 10 }}
          domain={xDomain} allowDataOverflow tickFormatter={fmtTick}
          label={{ value: xLabel, position: 'insideBottom', offset: -10, fontSize: 11, fill: '#6b7280' }} />
        <YAxis dataKey="y" type="number" name={yLabel} tick={{ fontSize: 10 }}
          domain={yDomain} allowDataOverflow tickFormatter={fmtTick}
          label={{ value: yLabel, angle: -90, position: 'insideLeft', offset: 10, fontSize: 11, fill: '#6b7280' }} />
        <ZAxis range={[40, 40]} />
        <Tooltip
          contentStyle={{ fontSize: 11, borderRadius: 8 }}
          formatter={(val, name) => [typeof val === 'number' ? val.toFixed(2) : val, name]}
          labelFormatter={(_, payload) => {
            const p = payload?.[0]?.payload
            if (!p) return ''
            const parts = [p.name]
            if (hasColor && p.c != null) parts.push(`${colorLabel}: ${typeof p.c === 'number' ? p.c.toFixed(2) : p.c}`)
            return parts.join(' | ')
          }}
        />
        <Scatter data={data} isAnimationActive={false}
          shape={(props) => {
            const { cx, cy, payload } = props
            let fill = '#3b82f6'
            if (hasColor && payload.c != null && typeof payload.c === 'number') {
              let norm = (payload.c - cMin) / cRange
              if (lowerBetter) norm = 1 - norm
              fill = colorScale(Math.max(0, Math.min(1, norm)))
            }
            return <circle cx={cx} cy={cy} r={4} fill={fill} stroke="#374151" strokeWidth={0.5} fillOpacity={0.8} />
          }}
        />
      </ScatterChart>
      {/* Color legend */}
      {hasColor && (
        <svg width={40} height={height - 40} style={{ position: 'absolute', right: 0, top: 10 }}>
          <defs>
            <linearGradient id="scatterColorGrad" x1="0" y1="0" x2="0" y2="1">
              <stop offset="0%" stopColor={lowerBetter ? colorScale(0) : colorScale(1)} />
              <stop offset="50%" stopColor={colorScale(0.5)} />
              <stop offset="100%" stopColor={lowerBetter ? colorScale(1) : colorScale(0)} />
            </linearGradient>
          </defs>
          <rect x={4} y={0} width={12} height={height - 70} fill="url(#scatterColorGrad)" rx={3} stroke="#d1d5db" strokeWidth={0.5} />
          <text x={20} y={10} fontSize={8} fill="#6b7280">{(cMin + cRange).toFixed(1)}</text>
          <text x={20} y={(height - 70) / 2 + 3} fontSize={8} fill="#6b7280">{(cMin + cRange / 2).toFixed(1)}</text>
          <text x={20} y={height - 72} fontSize={8} fill="#6b7280">{cMin.toFixed(1)}</text>
          <text x={0} y={height - 50} fontSize={8} fill="#6b7280" fontWeight={500}>{colorLabel}</text>
        </svg>
      )}
    </div>
  )
}
