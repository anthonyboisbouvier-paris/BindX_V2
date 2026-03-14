import React, { useEffect, useRef, useMemo, useState, useCallback } from 'react'
import useSettingsStore from '../stores/settingsStore.js'
import { generateDefaultCharts } from './analytics/defaultCharts.js'
import { CHART_TYPES } from './analytics/chartRegistry.js'
import { ALL_COLUMNS } from '../lib/columns.js'
import ChartCard from './analytics/ChartCard.jsx'
import { queryAgent } from '../api.js'

export default function AnalyticsPanel({ molecules, availableColumns, visibleKeys }) {
  const charts = useSettingsStore(s => s.analyticsCharts)
  const addChart = useSettingsStore(s => s.addAnalyticsChart)
  const updateChart = useSettingsStore(s => s.updateAnalyticsChart)
  const removeChart = useSettingsStore(s => s.removeAnalyticsChart)
  const resetCharts = useSettingsStore(s => s.resetAnalyticsCharts)

  // AI advisor state
  const [aiPrompt, setAiPrompt] = useState('')
  const [aiLoading, setAiLoading] = useState(false)
  const [aiError, setAiError] = useState(null)
  const [aiExplanation, setAiExplanation] = useState(null)

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

  // AI Chart Advisor — submit prompt
  const handleAiSubmit = useCallback(async () => {
    const prompt = aiPrompt.trim()
    if (!prompt || aiLoading) return

    setAiLoading(true)
    setAiError(null)
    setAiExplanation(null)

    try {
      const columns = selectableColumns.map(c => ({
        key: c.key,
        label: c.label || c.key,
        type: c.type || 'string',
        group: c.group || null,
      }))

      const result = await queryAgent('chart_advisor', {
        prompt,
        columns,
        molecule_count: molecules?.length || 0,
      })

      if (!result.available) {
        setAiError(result.error || 'AI advisor unavailable. Check that the OpenAI API key is configured.')
        return
      }

      const aiCharts = result.charts || []
      if (aiCharts.length === 0) {
        setAiError('No charts could be generated for this prompt. Try being more specific.')
        return
      }

      // Keep correlation matrix, replace everything else
      resetCharts()
      const numericCount = selectableColumns.filter(c => c.type === 'number').length
      if (numericCount >= 3) {
        addChart({ type: 'correlation' })
      }
      aiCharts.forEach(c => addChart(c))

      setAiExplanation(result.explanation || null)
      setAiPrompt('')
    } catch (err) {
      setAiError(err?.response?.data?.detail || err.message || 'Failed to generate charts')
    } finally {
      setAiLoading(false)
    }
  }, [aiPrompt, aiLoading, selectableColumns, molecules, resetCharts, addChart])

  // Handle Enter key in prompt input
  const handleKeyDown = (e) => {
    if (e.key === 'Enter' && !e.shiftKey) {
      e.preventDefault()
      handleAiSubmit()
    }
  }

  // Reset to defaults (correlation matrix only)
  const handleResetToDefaults = () => {
    resetCharts()
    const defaults = generateDefaultCharts(selectableColumns)
    defaults.forEach(c => addChart(c))
    setAiExplanation(null)
    setAiError(null)
  }

  return (
    <div className="space-y-3">
      {/* AI Chart Advisor */}
      <div className="bg-bx-s1/30 border border-bx-surface/30 rounded-lg p-3 space-y-2">
        <div className="flex items-center gap-2">
          <div className="flex-1 relative">
            <input
              type="text"
              value={aiPrompt}
              onChange={e => setAiPrompt(e.target.value)}
              onKeyDown={handleKeyDown}
              placeholder="Ask AI to generate charts... (e.g. 'show docking vs toxicity relationship')"
              disabled={aiLoading}
              className="w-full px-3 py-1.5 text-sm bg-bx-bg border border-bx-surface/30 rounded-md
                         text-bx-text placeholder:text-bx-dim
                         focus:outline-none focus:border-bx-mint focus:ring-1 focus:ring-bx-mint/30
                         disabled:opacity-50 transition-colors"
            />
            {aiLoading && (
              <div className="absolute right-3 top-1/2 -translate-y-1/2">
                <div className="w-4 h-4 border-2 border-bx-mint/30 border-t-bx-mint rounded-full animate-spin" />
              </div>
            )}
          </div>
          <button
            onClick={handleAiSubmit}
            disabled={!aiPrompt.trim() || aiLoading}
            className="px-3 py-1.5 text-sm font-medium bg-bx-mint text-white rounded-md
                       hover:bg-bx-mint-dim disabled:opacity-40 disabled:cursor-not-allowed
                       transition-colors whitespace-nowrap"
          >
            Generate
          </button>
          <button
            onClick={handleResetToDefaults}
            className="px-2 py-1.5 text-xs text-bx-dim hover:text-bx-sub
                       transition-colors whitespace-nowrap"
            title="Reset to default (correlation matrix only)"
          >
            Reset
          </button>
        </div>

        {aiError && (
          <p className="text-xs text-red-400">{aiError}</p>
        )}
        {aiExplanation && (
          <p className="text-xs text-bx-dim italic">{aiExplanation}</p>
        )}
      </div>

      {/* Toolbar — manual add buttons */}
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
          No charts configured. Use the AI advisor above or click a button to add one manually.
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
