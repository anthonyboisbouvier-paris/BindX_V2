import React, { useEffect, useRef, useMemo } from 'react'
import useSettingsStore from '../stores/settingsStore.js'
import { generateDefaultCharts } from './analytics/defaultCharts.js'
import { CHART_TYPES } from './analytics/chartRegistry.js'
import { ALL_COLUMNS } from '../lib/columns.js'
import ChartCard from './analytics/ChartCard.jsx'

export default function AnalyticsPanel({ molecules, availableColumns, visibleKeys }) {
  const charts = useSettingsStore(s => s.analyticsCharts)
  const addChart = useSettingsStore(s => s.addAnalyticsChart)
  const updateChart = useSettingsStore(s => s.updateAnalyticsChart)
  const removeChart = useSettingsStore(s => s.removeAnalyticsChart)
  const resetCharts = useSettingsStore(s => s.resetAnalyticsCharts)

  // Columns selectable = only those visible in the dashboard
  const selectableColumns = useMemo(() => {
    if (!visibleKeys?.length) return availableColumns || []
    const vSet = new Set(visibleKeys)
    return ALL_COLUMNS.filter(c => vSet.has(c.key))
  }, [visibleKeys, availableColumns])

  // Auto-populate default charts on first render when store is empty
  const initialized = useRef(false)
  useEffect(() => {
    if (initialized.current || charts.length > 0 || !selectableColumns?.length) return
    initialized.current = true
    const defaults = generateDefaultCharts(selectableColumns)
    defaults.forEach(c => addChart(c))
  }, [selectableColumns]) // eslint-disable-line react-hooks/exhaustive-deps

  const handleAdd = (type) => {
    // Pick sensible defaults from visible columns
    const numCols = selectableColumns.filter(c => c.type === 'number')
    const firstNum = numCols[0]?.key || 'docking_score'
    const secondNum = numCols[1]?.key || numCols[0]?.key || 'QED'
    const firstAny = selectableColumns[0]?.key || 'safety_color_code'

    const thirdNum = numCols[2]?.key || secondNum
    const defaults = {
      histogram: { type: 'histogram', key: firstNum, bins: 20 },
      bar: { type: 'bar', key: firstAny },
      pie: { type: 'pie', key: firstAny },
      scatter: { type: 'scatter', xKey: firstNum, yKey: secondNum },
      box: { type: 'box', metricKey: firstNum, groupKey: firstAny },
      correlation: { type: 'correlation' },
      bubble: { type: 'bubble', xKey: firstNum, yKey: secondNum, sizeKey: thirdNum },
      topn: { type: 'topn', key: firstNum, n: 15 },
    }
    addChart(defaults[type] || defaults.histogram)
  }

  return (
    <div className="space-y-3">
      {/* Toolbar */}
      <div className="flex items-center justify-between gap-2">
        <div className="flex items-center gap-2 flex-wrap">
          <span className="text-xs text-gray-400">{charts.length} chart{charts.length !== 1 ? 's' : ''}</span>
          <div className="flex gap-1">
            {CHART_TYPES.map(t => (
              <button
                key={t.type}
                onClick={() => handleAdd(t.type)}
                className="px-2 py-0.5 text-xs bg-gray-100 text-gray-600 hover:bg-gray-200 rounded-md transition-colors"
                title={`Add ${t.label}`}
              >
                + {t.label}
              </button>
            ))}
          </div>
        </div>
        {charts.length > 0 && (
          <button
            onClick={resetCharts}
            className="text-xs text-gray-400 hover:text-red-500 transition-colors"
            title="Remove all charts"
          >
            Clear all
          </button>
        )}
      </div>

      {/* Chart grid */}
      {charts.length === 0 ? (
        <div className="text-center py-8 text-gray-400 text-sm">
          No charts configured. Click a button above to add one, or they will auto-populate when data is available.
        </div>
      ) : (
        <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
          {charts.map(chart => (
            <div key={chart.id} className={chart.type === 'correlation' ? 'md:col-span-2' : ''}>
              <ChartCard
                config={chart}
                molecules={molecules}
                onUpdate={updateChart}
                onRemove={removeChart}
                selectableColumns={selectableColumns}
              />
            </div>
          ))}
        </div>
      )}
    </div>
  )
}
