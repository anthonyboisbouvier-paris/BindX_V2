import React, { useMemo } from 'react'
import { PieChart, Pie, Cell, Tooltip, Legend } from 'recharts'
import { ALL_COLUMNS } from '../../lib/columns.js'

const COLORS = ['#22c55e', '#ef4444', '#3b82f6', '#f59e0b', '#a855f7', '#06b6d4', '#ec4899']

export default function PieBreakdown({ config, molecules, width = 400, height = 220 }) {
  const { key } = config
  const colDef = ALL_COLUMNS.find(c => c.key === key)
  const label = colDef?.label || key

  const data = useMemo(() => {
    const counts = {}
    molecules.forEach(m => {
      const raw = m[key]
      let cat
      if (raw === true || raw === 1) cat = 'Pass'
      else if (raw === false || raw === 0) cat = 'Fail'
      else if (raw == null) cat = 'N/A'
      else cat = String(raw)
      counts[cat] = (counts[cat] || 0) + 1
    })
    return Object.entries(counts)
      .sort((a, b) => b[1] - a[1])
      .slice(0, 10) // Limit to top 10 categories for readability
      .map(([name, value]) => ({ name, value }))
  }, [molecules, key])

  if (!data.length) return <div className="text-gray-400 text-xs text-center py-6">No data for {label}</div>

  return (
    <PieChart width={width} height={height}>
      <Pie
        data={data}
        cx="50%"
        cy="45%"
        innerRadius={40}
        outerRadius={70}
        paddingAngle={2}
        dataKey="value"
        isAnimationActive={false}
        label={({ name, percent }) => `${name} ${(percent * 100).toFixed(0)}%`}
        labelLine={{ strokeWidth: 1 }}
        style={{ fontSize: 10 }}
      >
        {data.map((_, i) => (
          <Cell key={i} fill={COLORS[i % COLORS.length]} stroke="#fff" strokeWidth={2} />
        ))}
      </Pie>
      <Tooltip contentStyle={{ fontSize: 11, borderRadius: 8 }} />
      <Legend iconSize={8} wrapperStyle={{ fontSize: 10 }} />
    </PieChart>
  )
}
