import React, { useState, useMemo, useRef, useEffect } from 'react'

// ──────────────────────────────────────────────────────────────────────────
// Helpers
// ──────────────────────────────────────────────────────────────────────────

function SmilesCanvas({ smiles, width = 60, height = 40 }) {
  const canvasRef = useRef(null)
  const [failed, setFailed] = useState(false)
  useEffect(() => {
    if (!smiles || !canvasRef.current) return
    if (smiles === 'H' || smiles === '[H]' || smiles.length < 2) { setFailed(true); return }
    let cancelled = false
    setFailed(false)
    import('smiles-drawer').then(mod => {
      if (cancelled) return
      const SD = mod.default || mod
      const d = new SD.Drawer({ width, height, bondThickness: 1 })
      SD.parse(smiles, t => { if (!cancelled) d.draw(t, canvasRef.current, 'light') }, () => { if (!cancelled) setFailed(true) })
    }).catch(() => { if (!cancelled) setFailed(true) })
    return () => { cancelled = true }
  }, [smiles, width, height])
  if (failed || !smiles || smiles === 'H') return <span className="text-[10px] font-mono text-gray-400">{smiles || 'H'}</span>
  return <canvas ref={canvasRef} width={width} height={height} className="bg-white rounded" />
}

function fmt(v, d = 2) { return v == null ? '—' : v.toFixed(d) }

const LOWER_IS_BETTER = new Set(['docking_score', 'logP', 'MW', 'TPSA'])

/** Simple background: green if good, red if bad, white if null */
function cellBg(val, min, max, propKey) {
  if (val == null || min === max) return '#ffffff'
  let t = (val - min) / (max - min)            // 0 = min, 1 = max
  if (LOWER_IS_BETTER.has(propKey)) t = 1 - t  // invert: lower = good
  // t: 0 = bad, 1 = good → green
  const r = Math.round(255 - t * 80)
  const g = Math.round(220 + t * 35)
  const b = Math.round(220 - t * 60)
  return `rgb(${r},${g},${b})`
}

// ──────────────────────────────────────────────────────────────────────────
// Component
// ──────────────────────────────────────────────────────────────────────────

export default function SARMatrixView({
  rgroup,
  moleculeProperties = {},
  availableProperties = [],
  variableRGroups = [],
}) {
  const molecules = rgroup?.molecules || []
  const matrix = rgroup?.sar_matrix

  // Only one control: which property to display
  const numericProps = useMemo(() =>
    (availableProperties || []).filter(p => p.coverage >= 0.2),
  [availableProperties])
  const [displayProp, setDisplayProp] = useState(null)

  useEffect(() => {
    if (!numericProps.length) return
    if (!displayProp || !numericProps.find(p => p.key === displayProp)) {
      setDisplayProp(numericProps[0].key)
    }
  }, [numericProps]) // eslint-disable-line react-hooks/exhaustive-deps

  const propLabel = numericProps.find(p => p.key === displayProp)?.label || displayProp || 'Score'

  // Compute min/max for coloring
  const range = useMemo(() => {
    let min = Infinity, max = -Infinity
    for (const m of molecules) {
      const v = (moleculeProperties[m.id] || m.properties || {})[displayProp]
      if (typeof v === 'number') { if (v < min) min = v; if (v > max) max = v }
    }
    return { min: min === Infinity ? 0 : min, max: max === -Infinity ? 1 : max }
  }, [molecules, moleculeProperties, displayProp])

  // Guard
  if (!variableRGroups.length || molecules.length < 3) return null
  if (!matrix || !matrix.rows?.length || !matrix.cols?.length) {
    // Flat fallback
    return (
      <div className="card px-5 py-4">
        <h2 className="text-base font-semibold text-bx-light-text mb-2">SAR Table</h2>
        <HelpBlock />
        <p className="text-xs text-amber-600 bg-amber-50 border border-amber-200 rounded px-3 py-2 mb-3">
          Not enough structural variation to build a 2D matrix. Showing flat table instead.
        </p>
        <FlatView molecules={molecules} variableRGroups={variableRGroups} moleculeProperties={moleculeProperties} displayProp={displayProp} propLabel={propLabel} />
      </div>
    )
  }

  const { rows, cols, cells, r1_key, r2_key } = matrix

  return (
    <div className="card px-5 py-4">
      <div className="flex items-center justify-between mb-3">
        <h2 className="text-base font-semibold text-bx-light-text">
          SAR Matrix — {r1_key} vs {r2_key}
        </h2>
        <div className="flex items-center gap-2">
          <label className="text-xs text-gray-500">Property displayed:</label>
          <select
            value={displayProp || ''}
            onChange={e => setDisplayProp(e.target.value)}
            className="text-xs border border-gray-200 rounded px-2 py-1 bg-white"
          >
            {numericProps.map(p => <option key={p.key} value={p.key}>{p.label}</option>)}
          </select>
        </div>
      </div>

      <HelpBlock r1={r1_key} r2={r2_key} propLabel={propLabel} />

      {/* Matrix table */}
      <div className="overflow-x-auto mt-3">
        <table className="text-xs border-collapse w-auto">
          <thead>
            <tr>
              <th className="p-2 border border-gray-200 bg-gray-50 text-[10px] text-gray-500 font-medium">
                {r1_key} ↓ \ {r2_key} →
              </th>
              {cols.map((c, i) => (
                <th key={i} className="p-2 border border-gray-200 bg-gray-50 text-center" style={{ minWidth: 90 }}>
                  <SmilesCanvas smiles={c} width={55} height={38} />
                  <div className="text-[9px] text-gray-400 mt-0.5 truncate max-w-[80px]" title={c}>{c}</div>
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {rows.map((r, ri) => (
              <tr key={ri}>
                <td className="p-2 border border-gray-200 bg-gray-50 text-center" style={{ minWidth: 90 }}>
                  <SmilesCanvas smiles={r} width={55} height={38} />
                  <div className="text-[9px] text-gray-400 mt-0.5 truncate max-w-[80px]" title={r}>{r}</div>
                </td>
                {(cells[ri] || []).map((cell, ci) => {
                  if (!cell) {
                    return (
                      <td key={ci} className="p-2 border border-dashed border-gray-200 text-center bg-gray-50/30" title="Combination not tested — potential design opportunity">
                        <span className="text-gray-300 text-[10px]">not tested</span>
                      </td>
                    )
                  }
                  const props = cell.properties || {}
                  const val = props[displayProp]
                  const bg = cellBg(val, range.min, range.max, displayProp)
                  return (
                    <td key={ci} className="p-2 border border-gray-200 text-center relative" style={{ backgroundColor: bg }}
                      title={`${cell.name || cell.smiles}\n${propLabel}: ${fmt(val, 3)}${cell.n_molecules > 1 ? `\n(${cell.n_molecules} molecules)` : ''}`}
                    >
                      {cell.n_molecules > 1 && (
                        <span className="absolute top-0.5 right-0.5 text-[8px] bg-blue-500 text-white rounded-full w-4 h-4 flex items-center justify-center font-bold">
                          {cell.n_molecules}
                        </span>
                      )}
                      <div className="text-[11px] font-semibold">{fmt(val, 2)}</div>
                      <div className="text-[9px] text-gray-500 truncate max-w-[70px]" title={cell.name}>{cell.name || ''}</div>
                    </td>
                  )
                })}
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      {/* Legend */}
      <div className="mt-3 flex items-center gap-4 text-[10px] text-gray-400">
        <span>
          <span className="inline-block w-3 h-3 rounded-sm mr-1" style={{ backgroundColor: cellBg(range.min, range.min, range.max, displayProp) }} />
          {LOWER_IS_BETTER.has(displayProp) ? 'Best' : 'Worst'} ({fmt(range.min)})
        </span>
        <span>→</span>
        <span>
          <span className="inline-block w-3 h-3 rounded-sm mr-1" style={{ backgroundColor: cellBg(range.max, range.min, range.max, displayProp) }} />
          {LOWER_IS_BETTER.has(displayProp) ? 'Worst' : 'Best'} ({fmt(range.max)})
        </span>
        <span className="ml-4 border-l pl-4 border-gray-200">
          <span className="inline-block w-4 h-3 border border-dashed border-gray-300 rounded-sm mr-1 bg-gray-50/30" />
          = combination not yet tested (design opportunity)
        </span>
      </div>
    </div>
  )
}

// ──────────────────────────────────────────────────────────────────────────
// Explanatory block
// ──────────────────────────────────────────────────────────────────────────

function HelpBlock({ r1, r2, propLabel }) {
  return (
    <div className="text-xs text-gray-500 bg-blue-50 border border-blue-100 rounded-lg px-4 py-3 space-y-1.5">
      <p className="font-semibold text-blue-700">How to read this matrix</p>
      {r1 && r2 ? (
        <>
          <p>Each <strong>row</strong> = one substituent at position <strong>{r1}</strong>. Each <strong>column</strong> = one substituent at position <strong>{r2}</strong>.</p>
          <p>Each cell shows a molecule that has that specific combination of {r1} + {r2}. The number displayed is <strong>{propLabel || 'the selected property'}</strong>. The greener the cell, the better.</p>
          <p><strong>Read a row</strong> to see the effect of changing {r2} while keeping {r1} constant (no confounding).</p>
          <p><strong>Read a column</strong> to see the effect of changing {r1} while keeping {r2} constant.</p>
          <p>Cells marked <em>"not tested"</em> = a combination that doesn't exist yet in your dataset. These are opportunities for new molecule design.</p>
        </>
      ) : (
        <>
          <p>This table shows each molecule with its variable R-group substituents and the selected property value.</p>
          <p>Look for patterns: which substituents are associated with better values?</p>
        </>
      )}
    </div>
  )
}

// ──────────────────────────────────────────────────────────────────────────
// Flat fallback for <5 molecules or 1 variable R-group
// ──────────────────────────────────────────────────────────────────────────

function FlatView({ molecules, variableRGroups, moleculeProperties, displayProp, propLabel }) {
  return (
    <div className="overflow-x-auto">
      <table className="w-full text-xs border-collapse">
        <thead>
          <tr>
            <th className="p-2 border-b border-gray-200 text-left font-medium text-gray-500">Molecule</th>
            {variableRGroups.map(rk => (
              <th key={rk} className="p-2 border-b border-gray-200 text-center font-medium text-gray-500">{rk}</th>
            ))}
            <th className="p-2 border-b border-gray-200 text-right font-medium text-gray-500">{propLabel}</th>
          </tr>
        </thead>
        <tbody>
          {molecules.map((m, i) => {
            const val = (moleculeProperties?.[m.id] || m.properties || {})[displayProp]
            return (
              <tr key={m.id || i} className="hover:bg-gray-50">
                <td className="p-2 border-b border-gray-100 font-medium">{m.name || '—'}</td>
                {variableRGroups.map(rk => (
                  <td key={rk} className="p-2 border-b border-gray-100 text-center">
                    {m.r_groups_display?.[rk] ? <SmilesCanvas smiles={m.r_groups_display[rk]} width={50} height={35} /> : '—'}
                  </td>
                ))}
                <td className="p-2 border-b border-gray-100 text-right font-semibold">{fmt(val, 3)}</td>
              </tr>
            )
          })}
        </tbody>
      </table>
    </div>
  )
}
