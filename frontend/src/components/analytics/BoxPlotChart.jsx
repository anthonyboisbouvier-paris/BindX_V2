import React, { useMemo } from 'react'
import { BarChart, Bar, XAxis, YAxis, Tooltip, CartesianGrid } from 'recharts'
import { ALL_COLUMNS } from '../../lib/columns.js'

const PALETTE = ['#3b82f6', '#00e6a0', '#a855f7', '#f59e0b', '#ef4444', '#06b6d4', '#ec4899', '#84cc16']

function computeBoxStats(values) {
  if (!values.length) return null
  const sorted = [...values].sort((a, b) => a - b)
  const n = sorted.length
  const q1 = sorted[Math.floor(n * 0.25)]
  const median = sorted[Math.floor(n * 0.5)]
  const q3 = sorted[Math.floor(n * 0.75)]
  const min = sorted[0]
  const max = sorted[n - 1]
  const iqr = q3 - q1
  const whiskerLo = Math.max(min, q1 - 1.5 * iqr)
  const whiskerHi = Math.min(max, q3 + 1.5 * iqr)
  return { q1, median, q3, min, max, whiskerLo, whiskerHi, count: n }
}

// Truncate long group names
function truncate(str, max = 10) {
  return str.length > max ? str.slice(0, max) + '…' : str
}

export default function BoxPlotChart({ config, molecules, width = 400, height = 220 }) {
  const { metricKey, groupKey } = config
  const metricDef = ALL_COLUMNS.find(c => c.key === metricKey)
  const groupDef = ALL_COLUMNS.find(c => c.key === groupKey)
  const metricLabel = metricDef?.label || metricKey
  const groupLabel = groupDef?.label || groupKey

  const data = useMemo(() => {
    const groups = {}
    molecules.forEach(m => {
      const g = m[groupKey]
      const v = m[metricKey]
      if (g == null || v == null || typeof v !== 'number') return
      const gKey = String(g)
      if (!groups[gKey]) groups[gKey] = []
      groups[gKey].push(v)
    })

    return Object.entries(groups)
      .map(([name, values]) => {
        const stats = computeBoxStats(values)
        if (!stats) return null
        return { name, displayName: truncate(name), ...stats }
      })
      .filter(Boolean)
      .sort((a, b) => {
        const na = Number(a.name), nb = Number(b.name)
        if (!isNaN(na) && !isNaN(nb)) return na - nb
        return a.name.localeCompare(b.name)
      })
      .slice(0, 20)
  }, [molecules, metricKey, groupKey])

  if (!data.length) return <div className="text-gray-400 text-xs text-center py-6">No data for {metricLabel} by {groupLabel}</div>

  const yMin = Math.min(...data.map(d => d.whiskerLo))
  const yMax = Math.max(...data.map(d => d.whiskerHi))
  const yPad = (yMax - yMin) * 0.1 || 1

  // Determine if labels need rotation
  const needsRotation = data.length > 6 || data.some(d => d.name.length > 6)

  return (
    <BarChart width={width} height={height} data={data}
      margin={{ top: 10, right: 10, left: 0, bottom: needsRotation ? 30 : 10 }}>
      <CartesianGrid strokeDasharray="3 3" stroke="#f0f0f0" />
      <XAxis dataKey="displayName" tick={{ fontSize: 9 }}
        angle={needsRotation ? -40 : 0}
        textAnchor={needsRotation ? 'end' : 'middle'}
        height={needsRotation ? 50 : 30}
        interval={0}
        label={{ value: groupLabel, position: 'insideBottom', offset: needsRotation ? -20 : -2, fontSize: 10, fill: '#6b7280' }} />
      <YAxis tick={{ fontSize: 10 }} domain={[yMin - yPad, yMax + yPad]}
        tickFormatter={v => Math.abs(v) >= 1000 ? `${(v / 1000).toFixed(1)}k` : v.toFixed(1)}
        label={{ value: metricLabel, angle: -90, position: 'insideLeft', offset: 10, fontSize: 11, fill: '#6b7280' }} />
      <Tooltip
        contentStyle={{ fontSize: 11, borderRadius: 8 }}
        formatter={(_, __, props) => {
          const d = props.payload
          return [`Q1: ${d.q1?.toFixed(2)} | Med: ${d.median?.toFixed(2)} | Q3: ${d.q3?.toFixed(2)} | n=${d.count}`, metricLabel]
        }}
        labelFormatter={(_, payload) => payload?.[0]?.payload?.name || ''}
      />
      <Bar dataKey="q3" isAnimationActive={false}
        shape={(props) => {
          const { x, y, width: bw, payload, index } = props
          const color = PALETTE[index % PALETTE.length]
          const chartArea = props.background || {}
          const areaH = chartArea.height || (height - 40)
          const areaY = chartArea.y || 10
          const domain = [yMin - yPad, yMax + yPad]
          const scale = areaH / (domain[1] - domain[0])
          const toPixelY = (val) => areaY + areaH - (val - domain[0]) * scale

          const q1Y = toPixelY(payload.q1)
          const q3Y = toPixelY(payload.q3)
          const medY = toPixelY(payload.median)
          const loY = toPixelY(payload.whiskerLo)
          const hiY = toPixelY(payload.whiskerHi)
          const cx = x + bw / 2

          return (
            <g>
              {/* IQR box */}
              <rect x={x + 4} y={q3Y} width={bw - 8} height={q1Y - q3Y}
                fill={color} stroke="#374151" strokeWidth={1} rx={3} />
              {/* Median line */}
              <line x1={x + 4} y1={medY} x2={x + bw - 4} y2={medY}
                stroke="#fff" strokeWidth={2} />
              {/* Whisker lines */}
              <line x1={cx} y1={q3Y} x2={cx} y2={hiY} stroke="#6b7280" strokeWidth={1.5} />
              <line x1={cx} y1={q1Y} x2={cx} y2={loY} stroke="#6b7280" strokeWidth={1.5} />
              {/* Whisker caps */}
              <line x1={cx - 6} y1={hiY} x2={cx + 6} y2={hiY} stroke="#6b7280" strokeWidth={1.5} />
              <line x1={cx - 6} y1={loY} x2={cx + 6} y2={loY} stroke="#6b7280" strokeWidth={1.5} />
            </g>
          )
        }}
      />
    </BarChart>
  )
}
