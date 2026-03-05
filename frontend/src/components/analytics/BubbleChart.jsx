import React, { useMemo } from 'react'
import { ScatterChart, Scatter, XAxis, YAxis, ZAxis, Tooltip, CartesianGrid } from 'recharts'
import { ALL_COLUMNS } from '../../lib/columns.js'

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

export default function BubbleChart({ config, molecules, width = 400, height = 220 }) {
  const { xKey, yKey, sizeKey, colorKey } = config
  const xDef = ALL_COLUMNS.find(c => c.key === xKey)
  const yDef = ALL_COLUMNS.find(c => c.key === yKey)
  const sizeDef = ALL_COLUMNS.find(c => c.key === sizeKey)
  const colorDef = colorKey ? ALL_COLUMNS.find(c => c.key === colorKey) : null
  const xLabel = xDef?.label || xKey
  const yLabel = yDef?.label || yKey
  const sizeLabel = sizeDef?.label || sizeKey
  const colorLabel = colorDef?.label || colorKey
  const lowerBetter = colorDef?.colorScale === 'lower-better'

  const { data, cMin, cRange } = useMemo(() => {
    const mapped = molecules
      .map(m => ({
        x: m[xKey],
        y: m[yKey],
        z: m[sizeKey],
        c: colorKey ? m[colorKey] : undefined,
        name: m.name || '',
        bookmarked: m.bookmarked,
      }))
      .filter(d => d.x != null && d.y != null && d.z != null
        && typeof d.x === 'number' && typeof d.y === 'number' && typeof d.z === 'number')

    if (!colorKey) return { data: mapped, cMin: 0, cRange: 1 }
    const cVals = mapped.filter(d => d.c != null && typeof d.c === 'number').map(d => d.c)
    const cMin = cVals.length ? Math.min(...cVals) : 0
    const cMax = cVals.length ? Math.max(...cVals) : 1
    return { data: mapped, cMin, cRange: cMax - cMin || 1 }
  }, [molecules, xKey, yKey, sizeKey, colorKey])

  if (!data.length) return <div className="text-gray-400 text-xs text-center py-6">No data for {xLabel} / {yLabel} / {sizeLabel}</div>

  const zMin = Math.min(...data.map(d => d.z))
  const zMax = Math.max(...data.map(d => d.z))
  const hasColor = colorKey && data.some(d => d.c != null)
  const fmtZ = v => typeof v === 'number' ? parseFloat(v.toFixed(2)) : v
  const legendW = hasColor ? 50 : 0

  return (
    <div style={{ width, position: 'relative' }}>
      <ScatterChart width={width - legendW} height={height - 30} margin={{ top: 10, right: 10, left: 0, bottom: 20 }}>
        <CartesianGrid strokeDasharray="3 3" stroke="#f0f0f0" />
        <XAxis dataKey="x" type="number" name={xLabel} tick={{ fontSize: 10 }} tickFormatter={fmtTick}
          label={{ value: xLabel, position: 'insideBottom', offset: -10, fontSize: 11, fill: '#6b7280' }} />
        <YAxis dataKey="y" type="number" name={yLabel} tick={{ fontSize: 10 }} tickFormatter={fmtTick}
          label={{ value: yLabel, angle: -90, position: 'insideLeft', offset: 10, fontSize: 11, fill: '#6b7280' }} />
        <ZAxis dataKey="z" type="number" name={sizeLabel} range={[20, 400]} domain={[zMin, zMax]} />
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
            const zRange = zMax - zMin || 1
            const norm = (payload.z - zMin) / zRange
            const r = 3 + norm * 12
            let fill = '#3b82f6'
            let opacity = 0.5
            if (hasColor && payload.c != null && typeof payload.c === 'number') {
              let cNorm = (payload.c - cMin) / cRange
              if (lowerBetter) cNorm = 1 - cNorm
              fill = colorScale(Math.max(0, Math.min(1, cNorm)))
              opacity = 0.75
            }
            return (
              <g>
                <circle cx={cx} cy={cy} r={r} fill={fill} fillOpacity={opacity} stroke="#374151" strokeWidth={0.5} />
                {payload.bookmarked && (
                  <text x={cx} y={cy + 1} textAnchor="middle" dominantBaseline="middle"
                    fontSize={Math.max(8, r)} fill="#f59e0b" stroke="#fff" strokeWidth={0.3}>★</text>
                )}
              </g>
            )
          }}
        />
      </ScatterChart>
      {/* Color gradient legend (right) */}
      {hasColor && (
        <svg width={40} height={height - 60} style={{ position: 'absolute', right: 0, top: 10 }}>
          <defs>
            <linearGradient id="bubbleColorGrad" x1="0" y1="0" x2="0" y2="1">
              <stop offset="0%" stopColor={lowerBetter ? colorScale(0) : colorScale(1)} />
              <stop offset="50%" stopColor={colorScale(0.5)} />
              <stop offset="100%" stopColor={lowerBetter ? colorScale(1) : colorScale(0)} />
            </linearGradient>
          </defs>
          <rect x={4} y={0} width={12} height={height - 90} fill="url(#bubbleColorGrad)" rx={3} stroke="#d1d5db" strokeWidth={0.5} />
          <text x={20} y={10} fontSize={8} fill="#6b7280">{fmtZ(cMin + cRange)}</text>
          <text x={20} y={(height - 90) / 2 + 3} fontSize={8} fill="#6b7280">{fmtZ(cMin + cRange / 2)}</text>
          <text x={20} y={height - 92} fontSize={8} fill="#6b7280">{fmtZ(cMin)}</text>
          <text x={0} y={height - 70} fontSize={8} fill="#6b7280" fontWeight={500}>{colorLabel}</text>
        </svg>
      )}
      {/* Size legend (bottom) */}
      <div className="flex items-center justify-center gap-1" style={{ fontSize: 10, color: '#6b7280' }}>
        <span style={{ color: '#9ca3af' }}>{sizeLabel}:</span>
        <svg width={8} height={8}><circle cx={4} cy={4} r={3} fill="#94a3b8" fillOpacity={0.4}/></svg>
        <span>{fmtZ(zMin)}</span>
        <span style={{ color: '#d1d5db' }}>→</span>
        <svg width={16} height={16}><circle cx={8} cy={8} r={7} fill="#94a3b8" fillOpacity={0.4}/></svg>
        <span>{fmtZ(zMax)}</span>
      </div>
    </div>
  )
}
