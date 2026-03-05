import React, { useMemo } from 'react'
import { BarChart, Bar, XAxis, YAxis, Tooltip, CartesianGrid } from 'recharts'
import { ALL_COLUMNS } from '../../lib/columns.js'

const CATEGORY_COLORS = {
  green: '#22c55e', red: '#ef4444', yellow: '#eab308', orange: '#f97316',
  true: '#22c55e', false: '#ef4444',
}
const PALETTE = ['#3b82f6', '#00e6a0', '#a855f7', '#f59e0b', '#ef4444', '#06b6d4', '#ec4899', '#84cc16']

export default function CategoryBarChart({ config, molecules, width = 400, height = 220 }) {
  const { key } = config
  const colDef = ALL_COLUMNS.find(c => c.key === key)
  const label = colDef?.label || key

  const data = useMemo(() => {
    const counts = {}
    molecules.forEach(m => {
      const val = m[key]
      const cat = val == null ? 'N/A' : String(val)
      counts[cat] = (counts[cat] || 0) + 1
    })
    return Object.entries(counts)
      .sort((a, b) => b[1] - a[1])
      .map(([name, count], i) => ({
        name,
        count,
        barFill: CATEGORY_COLORS[name.toLowerCase()] || PALETTE[i % PALETTE.length],
      }))
  }, [molecules, key])

  if (!data.length) return <div className="text-gray-400 text-xs text-center py-6">No data for {label}</div>

  return (
    <BarChart width={width} height={height} data={data} margin={{ top: 5, right: 10, left: 0, bottom: 5 }}>
      <CartesianGrid strokeDasharray="3 3" stroke="#f0f0f0" />
      <XAxis dataKey="name" tick={{ fontSize: 10 }} />
      <YAxis tick={{ fontSize: 10 }} allowDecimals={false} />
      <Tooltip contentStyle={{ fontSize: 11, borderRadius: 8 }} />
      <Bar dataKey="count" isAnimationActive={false}
        shape={(props) => {
          const { x, y, width: w, height: h, payload } = props
          return <rect x={x} y={y} width={w} height={h} fill={payload.barFill} stroke="#374151" strokeWidth={0.5} rx={4} />
        }}
      />
    </BarChart>
  )
}
