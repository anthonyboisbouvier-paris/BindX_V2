import React, { useMemo } from 'react'
import { BarChart, Bar, XAxis, YAxis, Tooltip, CartesianGrid } from 'recharts'
import { ALL_COLUMNS } from '../../lib/columns.js'

const PALETTE = ['#3b82f6', '#00e6a0', '#a855f7', '#f59e0b', '#ef4444', '#06b6d4', '#ec4899', '#84cc16',
  '#14b8a6', '#8b5cf6', '#f97316', '#0ea5e9', '#e11d48', '#65a30d', '#7c3aed', '#0891b2']

const fmtTick = v => typeof v !== 'number' ? v : Number.isInteger(v) ? v : parseFloat(v.toFixed(3))

export default function TopNRanked({ config, molecules, width = 400, height = 220 }) {
  const { key, n = 15 } = config
  const colDef = ALL_COLUMNS.find(c => c.key === key)
  const label = colDef?.label || key
  const unit = colDef?.unit || ''
  const lowerBetter = colDef?.colorScale === 'lower-better'

  const data = useMemo(() => {
    const withVal = molecules
      .filter(m => m[key] != null && typeof m[key] === 'number')
      .map(m => ({
        name: m.name || m.id?.slice(0, 8) || '?',
        value: m[key],
        bookmarked: m.bookmarked,
      }))

    // Sort: lower-better = ascending, otherwise descending
    withVal.sort((a, b) => lowerBetter ? a.value - b.value : b.value - a.value)
    return withVal.slice(0, n)
  }, [molecules, key, n, lowerBetter])

  if (!data.length) return <div className="text-gray-400 text-xs text-center py-6">No data for {label}</div>

  // Horizontal bar chart — swap axes
  const barH = Math.max(180, data.length * 18)
  const effectiveH = Math.min(barH, height)

  return (
    <BarChart width={width} height={effectiveH} data={data} layout="vertical"
      margin={{ top: 5, right: 15, left: 5, bottom: 5 }}>
      <CartesianGrid strokeDasharray="3 3" stroke="#f0f0f0" horizontal={false} />
      <XAxis type="number" tick={{ fontSize: 9 }} tickFormatter={fmtTick}
        label={{ value: `${label}${unit ? ` (${unit})` : ''}`, position: 'insideBottom', offset: -2, fontSize: 10, fill: '#6b7280' }} />
      <YAxis type="category" dataKey="name" tick={{ fontSize: 8 }} width={70} />
      <Tooltip
        contentStyle={{ fontSize: 11, borderRadius: 8 }}
        formatter={(val) => [typeof val === 'number' ? val.toFixed(2) : val, label]}
      />
      <Bar dataKey="value" isAnimationActive={false}
        shape={(props) => {
          const { x, y, width: w, height: h, index, payload } = props
          const fill = payload.bookmarked ? '#f59e0b' : PALETTE[index % PALETTE.length]
          return (
            <g>
              <rect x={x} y={y} width={Math.max(0, w)} height={h} fill={fill} rx={3} />
              {payload.bookmarked && (
                <text x={x + Math.max(0, w) + 3} y={y + h / 2 + 3} fontSize={8} fill="#f59e0b">★</text>
              )}
            </g>
          )
        }}
      />
    </BarChart>
  )
}
