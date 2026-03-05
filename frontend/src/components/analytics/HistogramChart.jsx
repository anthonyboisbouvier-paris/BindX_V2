import React, { useMemo } from 'react'
import { BarChart, Bar, XAxis, YAxis, Tooltip, CartesianGrid, Legend } from 'recharts'
import { ALL_COLUMNS } from '../../lib/columns.js'

const PALETTE = ['#3b82f6', '#00e6a0', '#a855f7', '#f59e0b', '#ef4444', '#06b6d4', '#ec4899', '#84cc16',
  '#14b8a6', '#8b5cf6']

export default function HistogramChart({ config, molecules, width = 400, height = 220 }) {
  const { key, bins = 20, groupKey } = config
  const colDef = ALL_COLUMNS.find(c => c.key === key)
  const label = colDef?.label || key
  const unit = colDef?.unit || ''
  const groupDef = groupKey ? ALL_COLUMNS.find(c => c.key === groupKey) : null

  // Simple histogram (no grouping)
  const simpleData = useMemo(() => {
    if (groupKey) return null
    const values = molecules.map(m => m[key]).filter(v => v != null && typeof v === 'number')
    if (!values.length) return []

    const min = Math.min(...values)
    const max = Math.max(...values)
    const range = max - min || 1
    const nBins = Math.min(bins, Math.max(5, values.length))
    const binWidth = range / nBins

    const buckets = Array.from({ length: nBins }, (_, i) => ({
      binStart: min + i * binWidth,
      binEnd: min + (i + 1) * binWidth,
      count: 0,
    }))

    values.forEach(v => {
      const idx = Math.min(nBins - 1, Math.floor((v - min) / binWidth))
      buckets[idx].count++
    })

    return buckets.map(b => ({
      name: b.binStart.toFixed(1),
      count: b.count,
      range: `${b.binStart.toFixed(2)} – ${b.binEnd.toFixed(2)}`,
    }))
  }, [molecules, key, bins, groupKey])

  // Grouped histogram
  const { groupedData, groups } = useMemo(() => {
    if (!groupKey) return { groupedData: null, groups: [] }
    const values = molecules
      .map(m => ({ v: m[key], g: m[groupKey] }))
      .filter(d => d.v != null && typeof d.v === 'number' && d.g != null)
    if (!values.length) return { groupedData: [], groups: [] }

    // Discover groups (limit to top 8)
    const groupCounts = {}
    values.forEach(d => { groupCounts[String(d.g)] = (groupCounts[String(d.g)] || 0) + 1 })
    const topGroups = Object.entries(groupCounts)
      .sort((a, b) => b[1] - a[1])
      .slice(0, 8)
      .map(([g]) => g)
    const groupSet = new Set(topGroups)

    // Compute bins
    const allVals = values.filter(d => groupSet.has(String(d.g))).map(d => d.v)
    const min = Math.min(...allVals)
    const max = Math.max(...allVals)
    const range = max - min || 1
    const nBins = Math.min(bins, Math.max(5, allVals.length))
    const binWidth = range / nBins

    const buckets = Array.from({ length: nBins }, (_, i) => {
      const row = {
        name: (min + i * binWidth).toFixed(1),
        range: `${(min + i * binWidth).toFixed(2)} – ${(min + (i + 1) * binWidth).toFixed(2)}`,
      }
      topGroups.forEach(g => { row[g] = 0 })
      return row
    })

    values.forEach(d => {
      const gStr = String(d.g)
      if (!groupSet.has(gStr)) return
      const idx = Math.min(nBins - 1, Math.floor((d.v - min) / binWidth))
      buckets[idx][gStr]++
    })

    return { groupedData: buckets, groups: topGroups }
  }, [molecules, key, bins, groupKey])

  const data = groupKey ? groupedData : simpleData
  if (!data?.length) return <div className="text-gray-400 text-xs text-center py-6">No numeric data for {label}</div>

  // Tick formatting: skip every other tick if too many
  const showEveryN = data.length > 15 ? 3 : data.length > 8 ? 2 : 1

  return (
    <BarChart width={width} height={height} data={data} margin={{ top: 5, right: 10, left: 0, bottom: 25 }}>
      <CartesianGrid strokeDasharray="3 3" stroke="#f0f0f0" />
      <XAxis dataKey="name" tick={{ fontSize: 9 }} angle={-35} textAnchor="end" height={45}
        interval={showEveryN - 1}
        label={{ value: `${label}${unit ? ` (${unit})` : ''}`, position: 'insideBottom', offset: -15, fontSize: 10, fill: '#6b7280' }} />
      <YAxis tick={{ fontSize: 10 }} allowDecimals={false}
        label={{ value: 'Count', angle: -90, position: 'insideLeft', offset: 10, fontSize: 11, fill: '#6b7280' }} />
      <Tooltip
        formatter={(val, name) => [val, groupKey ? name : 'Count']}
        labelFormatter={(_, payload) => payload?.[0]?.payload?.range || ''}
        contentStyle={{ fontSize: 11, borderRadius: 8 }}
      />
      {groupKey && groups.length > 0 ? (
        <>
          {groups.map((g, i) => (
            <Bar key={g} dataKey={g} stackId="stack" isAnimationActive={false}
              shape={(props) => {
                const { x, y, width: w, height: h } = props
                return <rect x={x} y={y} width={w} height={Math.max(0, h)} fill={PALETTE[i % PALETTE.length]} stroke="#fff" strokeWidth={0.5} rx={1} />
              }}
            />
          ))}
          <Legend iconSize={8} wrapperStyle={{ fontSize: 9, paddingTop: 2 }} />
        </>
      ) : (
        <Bar dataKey="count" isAnimationActive={false}
          shape={(props) => {
            const { x, y, width: w, height: h } = props
            return <rect x={x} y={y} width={w} height={h} fill="#3b82f6" stroke="#2563eb" strokeWidth={1} rx={2} />
          }}
        />
      )}
    </BarChart>
  )
}
