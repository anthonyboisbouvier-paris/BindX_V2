import React, { useState, useRef, useCallback } from 'react'

/**
 * OptimizationChart — SVG line chart for score evolution during lead optimization.
 * Props:
 *   iterations: Array<{ iteration: number, best_score: number }>
 */
export default function OptimizationChart({ iterations = [] }) {
  const [tooltip, setTooltip] = useState(null)
  const svgRef = useRef(null)

  const WIDTH = 480
  const HEIGHT = 200
  const PAD = { top: 16, right: 20, bottom: 36, left: 52 }
  const INNER_W = WIDTH - PAD.left - PAD.right
  const INNER_H = HEIGHT - PAD.top - PAD.bottom

  if (!iterations || iterations.length === 0) {
    return (
      <div
        className="flex items-center justify-center bg-gray-50 rounded-xl border border-gray-100 text-gray-300 text-sm"
        style={{ height: HEIGHT }}
      >
        No data yet
      </div>
    )
  }

  // Compute scales
  const scores = iterations.map((d) => d.best_score)
  const minScore = Math.min(...scores)
  const maxScore = Math.max(...scores)
  const scoreRange = maxScore - minScore || 1
  const minIter = iterations[0].iteration
  const maxIter = iterations[iterations.length - 1].iteration
  const iterRange = maxIter - minIter || 1

  const xScale = (iter) => ((iter - minIter) / iterRange) * INNER_W
  const yScale = (score) => INNER_H - ((score - minScore) / scoreRange) * INNER_H

  // Build path points
  const points = iterations.map((d) => ({
    x: xScale(d.iteration),
    y: yScale(d.best_score),
    ...d,
  }))

  const linePath = points
    .map((p, i) => `${i === 0 ? 'M' : 'L'} ${p.x} ${p.y}`)
    .join(' ')

  const areaPath =
    `M ${points[0].x} ${INNER_H}` +
    points.map((p) => ` L ${p.x} ${p.y}`).join('') +
    ` L ${points[points.length - 1].x} ${INNER_H} Z`

  // Y-axis ticks (4 ticks)
  const yTicks = [0, 0.33, 0.67, 1].map((t) => ({
    y: t * INNER_H,
    label: (minScore + (1 - t) * scoreRange).toFixed(2),
  }))

  // X-axis ticks (up to 6 ticks)
  const xTickCount = Math.min(6, iterations.length)
  const xTickStep = Math.max(1, Math.floor(iterations.length / xTickCount))
  const xTicks = iterations
    .filter((_, i) => i % xTickStep === 0 || i === iterations.length - 1)
    .map((d) => ({ x: xScale(d.iteration), label: `#${d.iteration}` }))

  const handleMouseMove = useCallback(
    (e) => {
      if (!svgRef.current) return
      const rect = svgRef.current.getBoundingClientRect()
      const mouseX = ((e.clientX - rect.left) / rect.width) * WIDTH - PAD.left
      // Find closest point
      let closest = points[0]
      let minDist = Infinity
      for (const p of points) {
        const dist = Math.abs(p.x - mouseX)
        if (dist < minDist) {
          minDist = dist
          closest = p
        }
      }
      if (minDist < 30) {
        setTooltip(closest)
      } else {
        setTooltip(null)
      }
    },
    [points]
  )

  const gradientId = 'opt-chart-gradient'

  return (
    <div className="relative">
      <svg
        ref={svgRef}
        viewBox={`0 0 ${WIDTH} ${HEIGHT}`}
        className="w-full rounded-xl bg-gray-50"
        style={{ maxHeight: HEIGHT + 20 }}
        onMouseMove={handleMouseMove}
        onMouseLeave={() => setTooltip(null)}
      >
        <defs>
          <linearGradient id={gradientId} x1="0" y1="0" x2="0" y2="1">
            <stop offset="0%" stopColor="#22c55e" stopOpacity="0.25" />
            <stop offset="100%" stopColor="#22c55e" stopOpacity="0.02" />
          </linearGradient>
        </defs>

        <g transform={`translate(${PAD.left},${PAD.top})`}>
          {/* Y grid lines */}
          {yTicks.map((t, i) => (
            <g key={i}>
              <line
                x1={0} y1={t.y}
                x2={INNER_W} y2={t.y}
                stroke="#e5e7eb"
                strokeWidth={0.8}
                strokeDasharray={i === 0 ? undefined : '3,3'}
              />
              <text
                x={-8} y={t.y + 4}
                textAnchor="end"
                fontSize={9}
                fill="#9ca3af"
                fontFamily="monospace"
              >
                {t.label}
              </text>
            </g>
          ))}

          {/* X axis ticks */}
          {xTicks.map((t, i) => (
            <g key={i}>
              <line
                x1={t.x} y1={INNER_H}
                x2={t.x} y2={INNER_H + 4}
                stroke="#d1d5db"
                strokeWidth={0.8}
              />
              <text
                x={t.x} y={INNER_H + 14}
                textAnchor="middle"
                fontSize={9}
                fill="#9ca3af"
                fontFamily="monospace"
              >
                {t.label}
              </text>
            </g>
          ))}

          {/* Area fill */}
          <path d={areaPath} fill={`url(#${gradientId})`} />

          {/* Line */}
          <path
            d={linePath}
            fill="none"
            stroke="#1e3a5f"
            strokeWidth={2}
            strokeLinejoin="round"
            strokeLinecap="round"
          />

          {/* Dots */}
          {points.map((p, i) => (
            <circle
              key={i}
              cx={p.x}
              cy={p.y}
              r={tooltip && tooltip.iteration === p.iteration ? 5 : 3}
              fill={tooltip && tooltip.iteration === p.iteration ? '#22c55e' : '#1e3a5f'}
              stroke="white"
              strokeWidth={1.5}
              style={{ transition: 'r 0.1s, fill 0.1s' }}
            />
          ))}

          {/* Axis labels */}
          <text
            x={INNER_W / 2} y={INNER_H + 30}
            textAnchor="middle"
            fontSize={10}
            fill="#6b7280"
          >
            Iteration
          </text>
          <text
            x={-INNER_H / 2} y={-40}
            textAnchor="middle"
            fontSize={10}
            fill="#6b7280"
            transform="rotate(-90)"
          >
            Best Score
          </text>

          {/* Tooltip vertical line */}
          {tooltip && (
            <line
              x1={tooltip.x} y1={0}
              x2={tooltip.x} y2={INNER_H}
              stroke="#22c55e"
              strokeWidth={1}
              strokeDasharray="4,3"
              opacity={0.6}
            />
          )}
        </g>
      </svg>

      {/* Tooltip box */}
      {tooltip && (
        <div
          className="absolute pointer-events-none bg-gray-900 text-white text-xs px-2.5 py-1.5 rounded-lg shadow-lg"
          style={{
            left: `calc(${((tooltip.x + PAD.left) / WIDTH) * 100}% + 8px)`,
            top: `calc(${((tooltip.y + PAD.top) / HEIGHT) * 100}%)`,
            transform: 'translateY(-50%)',
          }}
        >
          <div className="font-semibold">Iteration {tooltip.iteration}</div>
          <div className="text-green-400 font-mono">{Number(tooltip.best_score).toFixed(4)}</div>
        </div>
      )}
    </div>
  )
}
