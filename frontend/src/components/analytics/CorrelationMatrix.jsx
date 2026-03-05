import React, { useMemo, useState } from 'react'
import { ALL_COLUMNS } from '../../lib/columns.js'

function pearson(xs, ys) {
  const n = xs.length
  if (n < 3) return null
  const mx = xs.reduce((a, b) => a + b, 0) / n
  const my = ys.reduce((a, b) => a + b, 0) / n
  let num = 0, dx2 = 0, dy2 = 0
  for (let i = 0; i < n; i++) {
    const dx = xs[i] - mx, dy = ys[i] - my
    num += dx * dy; dx2 += dx * dx; dy2 += dy * dy
  }
  const den = Math.sqrt(dx2 * dy2)
  return den > 0 ? num / den : 0
}

function corrColor(r) {
  if (r == null) return '#f3f4f6'
  const abs = Math.min(1, Math.abs(r))
  if (r >= 0) {
    const g = Math.round(255 - abs * 120)
    const rb = Math.round(255 - abs * 200)
    return `rgb(${rb}, ${g}, 255)`
  } else {
    const g = Math.round(255 - abs * 120)
    const b = Math.round(255 - abs * 200)
    return `rgb(255, ${g}, ${b})`
  }
}

export default function CorrelationMatrix({ config, molecules, width = 400 }) {
  const [hovered, setHovered] = useState(null)

  // Filter to selected columns (corrColumns = undefined means all)
  const numCols = useMemo(() => {
    const candidates = ALL_COLUMNS.filter(c => c.type === 'number')
    const withData = candidates.filter(c =>
      molecules.some(m => m[c.key] != null && typeof m[c.key] === 'number')
    )
    if (config.corrColumns?.length) {
      const allowed = new Set(config.corrColumns)
      return withData.filter(c => allowed.has(c.key))
    }
    return withData.slice(0, 12)
  }, [molecules, config.corrColumns])

  const matrix = useMemo(() => {
    const n = numCols.length
    const result = Array.from({ length: n }, () => Array(n).fill(null))
    const vectors = numCols.map(col =>
      molecules.map(m => m[col.key]).map(v => (v != null && typeof v === 'number') ? v : null)
    )
    for (let i = 0; i < n; i++) {
      for (let j = i; j < n; j++) {
        const xs = [], ys = []
        for (let k = 0; k < molecules.length; k++) {
          if (vectors[i][k] != null && vectors[j][k] != null) {
            xs.push(vectors[i][k])
            ys.push(vectors[j][k])
          }
        }
        const r = pearson(xs, ys)
        result[i][j] = r
        result[j][i] = r
      }
    }
    return result
  }, [numCols, molecules])

  if (numCols.length < 2) return <div className="text-gray-400 text-xs text-center py-6">Need at least 2 numeric columns</div>

  const n = numCols.length
  // Layout: row labels left, grid, column labels bottom, legend right
  const labelW = 80
  const legendW = 50
  const bottomLabelH = 80
  const padding = 4
  const availW = width - labelW - legendW - padding * 2
  // Smaller cells: cap at 36px
  const cellSize = Math.max(16, Math.min(Math.floor(availW / n), 36))
  const gridW = cellSize * n
  const gridH = cellSize * n
  const gridTop = padding
  const svgW = labelW + gridW + legendW + padding * 2
  const svgH = gridTop + gridH + bottomLabelH

  return (
    <div style={{ width, display: 'flex', flexDirection: 'column', alignItems: 'center' }}>
      <svg width={svgW} height={svgH} style={{ display: 'block' }}>
        {/* Row labels (left) */}
        {numCols.map((col, i) => (
          <text key={`rl-${i}`}
            x={labelW - 6}
            y={gridTop + i * cellSize + cellSize / 2 + 4}
            textAnchor="end" fontSize={Math.min(10, cellSize - 2)} fill="#374151" fontWeight={500}>
            {col.label}
          </text>
        ))}
        {/* Cells */}
        {matrix.map((row, i) =>
          row.map((r, j) => {
            const x = labelW + j * cellSize
            const y = gridTop + i * cellSize
            const isHov = hovered?.i === i && hovered?.j === j
            return (
              <g key={`${i}-${j}`}
                onMouseEnter={() => setHovered({ i, j, r })}
                onMouseLeave={() => setHovered(null)}>
                <rect x={x} y={y} width={cellSize - 1} height={cellSize - 1}
                  fill={corrColor(r)} stroke={isHov ? '#1e293b' : '#e5e7eb'} strokeWidth={isHov ? 2 : 0.5} rx={2}
                  style={{ cursor: 'default' }} />
                {cellSize >= 28 && r != null && (
                  <text x={x + cellSize / 2 - 0.5} y={y + cellSize / 2 + 3}
                    textAnchor="middle" fontSize={cellSize >= 34 ? 9 : 7}
                    fontWeight={600} fill={Math.abs(r) > 0.5 ? '#fff' : '#374151'}>
                    {r.toFixed(cellSize >= 34 ? 2 : 1)}
                  </text>
                )}
              </g>
            )
          })
        )}
        {/* Column labels (bottom, rotated) */}
        {numCols.map((col, j) => {
          const cx = labelW + j * cellSize + cellSize / 2
          const cy = gridTop + gridH + 6
          return (
            <text key={`cl-${j}`}
              x={cx} y={cy}
              textAnchor="start" fontSize={Math.min(10, cellSize - 2)} fill="#374151" fontWeight={500}
              transform={`rotate(55, ${cx}, ${cy})`}>
              {col.label}
            </text>
          )
        })}
        {/* Color legend bar (right side) */}
        <defs>
          <linearGradient id="corrGrad" x1="0" y1="0" x2="0" y2="1">
            <stop offset="0%" stopColor={corrColor(1)} />
            <stop offset="50%" stopColor={corrColor(0)} />
            <stop offset="100%" stopColor={corrColor(-1)} />
          </linearGradient>
        </defs>
        <rect x={labelW + gridW + 14} y={gridTop} width={14} height={gridH}
          fill="url(#corrGrad)" rx={3} stroke="#d1d5db" strokeWidth={0.5} />
        <text x={labelW + gridW + 34} y={gridTop + 10} fontSize={9} fill="#6b7280" fontWeight={600}>+1</text>
        <text x={labelW + gridW + 34} y={gridTop + gridH / 2 + 3} fontSize={9} fill="#6b7280" fontWeight={600}>0</text>
        <text x={labelW + gridW + 34} y={gridTop + gridH - 3} fontSize={9} fill="#6b7280" fontWeight={600}>-1</text>
      </svg>
      {/* Tooltip info bar */}
      <div className="text-xs text-gray-600 px-1 py-1 flex items-center justify-center" style={{ minHeight: 22 }}>
        {hovered && hovered.r != null ? (
          <>
            <strong>{numCols[hovered.i].label}</strong>
            {' vs '}
            <strong>{numCols[hovered.j].label}</strong>
            {`: r = ${hovered.r.toFixed(3)} (R² = ${(hovered.r ** 2).toFixed(3)})`}
          </>
        ) : (
          <span className="text-gray-400">Hover a cell to see details</span>
        )}
      </div>
    </div>
  )
}
