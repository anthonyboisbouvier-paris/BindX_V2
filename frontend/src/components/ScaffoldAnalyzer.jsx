import React, { useState, useEffect, useCallback, useMemo } from 'react'
import { analyzeScaffold } from '../api.js'

// ---------------------------------------------------------------------------
// Strategy labels
// ---------------------------------------------------------------------------

const STRATEGY_OPTIONS = [
  { value: 'any', label: 'Any' },
  { value: 'add_fg', label: 'Add functional group' },
  { value: 'swap_halogen', label: 'Swap halogen' },
  { value: 'swap_atom', label: 'Swap atom' },
  { value: 'modify_chain', label: 'Modify chain' },
]

const ALL_STRATEGIES = ['add_fg', 'swap_halogen', 'swap_atom', 'modify_chain']

// ---------------------------------------------------------------------------
// PositionCard
// ---------------------------------------------------------------------------

function PositionCard({ position, state, onToggleFreeze, onStrategyChange, onToggleGroup }) {
  const isFrozen = state.frozen
  const selectedGroups = state.selectedGroups || []

  return (
    <div className={`rounded-lg border p-3 transition-all ${
      isFrozen
        ? 'bg-gray-50 border-gray-200 opacity-75'
        : 'bg-white border-green-200'
    }`}>
      <div className="flex items-center justify-between mb-2">
        <div className="flex items-center gap-2">
          <span className={`inline-flex items-center justify-center w-8 h-8 rounded-full text-xs font-bold ${
            isFrozen ? 'bg-gray-200 text-gray-500' : 'bg-green-100 text-green-700'
          }`}>
            {position.label}
          </span>
          <div>
            <span className="text-sm font-semibold text-gray-700">{position.atom_symbol}</span>
            <span className="text-xs text-gray-400 ml-1.5 font-mono">{position.current_group}</span>
          </div>
          {position.is_brics_site && (
            <span className="text-[10px] font-semibold px-1.5 py-0.5 rounded bg-red-50 text-red-600 uppercase">BRICS</span>
          )}
        </div>
        <button
          onClick={() => onToggleFreeze(position.position_idx)}
          className={`flex items-center gap-1 px-2.5 py-1 rounded-md text-xs font-semibold transition-colors ${
            isFrozen
              ? 'bg-gray-200 text-gray-600 hover:bg-gray-300'
              : 'bg-green-100 text-green-700 hover:bg-green-200'
          }`}
        >
          {isFrozen ? (
            <>
              <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                  d="M12 15v2m-6 4h12a2 2 0 002-2v-6a2 2 0 00-2-2H6a2 2 0 00-2 2v6a2 2 0 002 2zm10-10V7a4 4 0 00-8 0v4h8z" />
              </svg>
              Frozen
            </>
          ) : (
            <>
              <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                  d="M11 5H6a2 2 0 00-2 2v11a2 2 0 002 2h11a2 2 0 002-2v-5m-1.414-9.414a2 2 0 112.828 2.828L11.828 15H9v-2.828l8.586-8.586z" />
              </svg>
              Modify
            </>
          )}
        </button>
      </div>

      {!isFrozen && (
        <div className="space-y-2 mt-2">
          <div>
            <label className="text-[10px] font-semibold text-gray-400 uppercase">Strategy</label>
            <select
              value={state.strategy}
              onChange={(e) => onStrategyChange(position.position_idx, e.target.value)}
              className="mt-0.5 w-full text-xs border border-gray-200 rounded-md px-2 py-1.5 bg-white text-gray-700 focus:ring-1 focus:ring-green-400 focus:border-green-400"
            >
              {STRATEGY_OPTIONS
                .filter(opt => opt.value === 'any' || position.applicable_strategies.includes(opt.value))
                .map(opt => (
                  <option key={opt.value} value={opt.value}>{opt.label}</option>
                ))}
            </select>
          </div>

          {position.suggested_replacements.length > 0 && (
            <div>
              <label className="text-[10px] font-semibold text-gray-400 uppercase">Groups</label>
              <div className="flex flex-wrap gap-1 mt-1">
                {position.suggested_replacements.map(group => {
                  const isActive = selectedGroups.includes(group)
                  return (
                    <button
                      key={group}
                      onClick={() => onToggleGroup(position.position_idx, group)}
                      className={`px-2 py-0.5 rounded text-[11px] font-mono font-semibold transition-colors ${
                        isActive
                          ? 'bg-[#0f131d] text-white'
                          : 'bg-gray-100 text-gray-500 hover:bg-gray-200'
                      }`}
                    >
                      {group}
                    </button>
                  )
                })}
              </div>
              <p className="text-[10px] text-gray-400 mt-0.5">
                {selectedGroups.length === 0 ? 'All groups allowed' : `${selectedGroups.length} selected`}
              </p>
            </div>
          )}
        </div>
      )}
    </div>
  )
}

// ---------------------------------------------------------------------------
// GlobalControls
// ---------------------------------------------------------------------------

function GlobalControls({ preserveScaffold, setPreserveScaffold, similarity, setSimilarity, maxMwChange, setMaxMwChange, globalStrategies, setGlobalStrategies }) {
  return (
    <div className="space-y-3">
      <div className="flex items-center gap-3">
        <label className="flex items-center gap-2 cursor-pointer">
          <input type="checkbox" checked={preserveScaffold} onChange={(e) => setPreserveScaffold(e.target.checked)}
            className="w-4 h-4 rounded border-gray-300 text-[#00e6a0] focus:ring-[#00e6a0]" />
          <span className="text-xs font-semibold text-gray-600">Preserve scaffold</span>
        </label>
      </div>

      <div className="grid grid-cols-2 gap-4">
        <div>
          <div className="flex items-center justify-between mb-1">
            <label className="text-[10px] font-semibold text-gray-400 uppercase">Min similarity</label>
            <span className="text-xs font-bold text-[#0f131d] tabular-nums">{similarity.toFixed(2)}</span>
          </div>
          <input type="range" min={0.1} max={0.9} step={0.05} value={similarity}
            onChange={(e) => setSimilarity(parseFloat(e.target.value))}
            className="w-full h-1.5 bg-gray-200 rounded-lg appearance-none cursor-pointer accent-[#00e6a0]" />
        </div>
        <div>
          <div className="flex items-center justify-between mb-1">
            <label className="text-[10px] font-semibold text-gray-400 uppercase">Max MW change</label>
            <span className="text-xs font-bold text-[#0f131d] tabular-nums">{maxMwChange} Da</span>
          </div>
          <input type="range" min={10} max={500} step={10} value={maxMwChange}
            onChange={(e) => setMaxMwChange(parseInt(e.target.value))}
            className="w-full h-1.5 bg-gray-200 rounded-lg appearance-none cursor-pointer accent-[#00e6a0]" />
        </div>
      </div>

      <div>
        <label className="text-[10px] font-semibold text-gray-400 uppercase mb-1.5 block">Allowed strategies</label>
        <div className="flex flex-wrap gap-2">
          {ALL_STRATEGIES.map(strat => {
            const active = globalStrategies.length === 0 || globalStrategies.includes(strat)
            const label = STRATEGY_OPTIONS.find(o => o.value === strat)?.label || strat
            return (
              <label key={strat} className="flex items-center gap-1.5 cursor-pointer">
                <input
                  type="checkbox"
                  checked={active}
                  onChange={(e) => {
                    if (e.target.checked) {
                      // If currently filtering, add this one
                      if (globalStrategies.length > 0) {
                        setGlobalStrategies([...globalStrategies, strat])
                      }
                      // If all were checked (empty = all), do nothing
                    } else {
                      // Remove this one
                      if (globalStrategies.length === 0) {
                        // Was "all", now filter out this one
                        setGlobalStrategies(ALL_STRATEGIES.filter(s => s !== strat))
                      } else {
                        const next = globalStrategies.filter(s => s !== strat)
                        setGlobalStrategies(next.length === ALL_STRATEGIES.length ? [] : next)
                      }
                    }
                  }}
                  className="w-3.5 h-3.5 rounded border-gray-300 text-[#00e6a0] focus:ring-[#00e6a0]"
                />
                <span className="text-xs text-gray-600">{label}</span>
              </label>
            )
          })}
        </div>
      </div>
    </div>
  )
}

// ---------------------------------------------------------------------------
// ScaffoldAnalyzer (main)
// ---------------------------------------------------------------------------

export default function ScaffoldAnalyzer({ smiles, onRulesChange }) {
  const [analysis, setAnalysis] = useState(null)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState(null)

  // Per-position states: { [position_idx]: { frozen, strategy, selectedGroups } }
  const [positionStates, setPositionStates] = useState({})

  // Global controls
  const [preserveScaffold, setPreserveScaffold] = useState(true)
  const [similarity, setSimilarity] = useState(0.3)
  const [maxMwChange, setMaxMwChange] = useState(100)
  const [globalStrategies, setGlobalStrategies] = useState([]) // empty = all

  // Fetch scaffold analysis
  useEffect(() => {
    if (!smiles) return
    let cancelled = false
    setLoading(true)
    setError(null)
    analyzeScaffold(smiles)
      .then(data => {
        if (cancelled) return
        setAnalysis(data)
        // Initialize position states
        const states = {}
        for (const pos of (data.positions || [])) {
          states[pos.position_idx] = {
            frozen: false,
            strategy: 'any',
            selectedGroups: [],
          }
        }
        setPositionStates(states)
      })
      .catch(err => {
        if (cancelled) return
        setError(err.userMessage || err.message || 'Failed to analyze scaffold')
        setAnalysis(null)
      })
      .finally(() => {
        if (!cancelled) setLoading(false)
      })
    return () => { cancelled = true }
  }, [smiles])

  // Build rules object whenever state changes
  const rules = useMemo(() => {
    if (!analysis) return null
    const posRules = []
    const frozenPositions = []

    for (const pos of (analysis.positions || [])) {
      const st = positionStates[pos.position_idx]
      if (!st) continue
      if (st.frozen) {
        frozenPositions.push(pos.position_idx)
      } else {
        posRules.push({
          position_idx: pos.position_idx,
          strategy: st.strategy || 'any',
          allowed_groups: st.selectedGroups || [],
          frozen: false,
        })
      }
    }

    return {
      rules: posRules,
      frozen_positions: frozenPositions,
      allowed_strategies: globalStrategies,
      preserve_scaffold: preserveScaffold,
      min_similarity: similarity,
      max_mw_change: maxMwChange,
      core_atom_indices: analysis.core_atom_indices || [],
    }
  }, [analysis, positionStates, preserveScaffold, similarity, maxMwChange, globalStrategies])

  // Notify parent
  useEffect(() => {
    if (onRulesChange) onRulesChange(rules)
  }, [rules, onRulesChange])

  const handleToggleFreeze = useCallback((idx) => {
    setPositionStates(prev => ({
      ...prev,
      [idx]: { ...prev[idx], frozen: !prev[idx]?.frozen },
    }))
  }, [])

  const handleStrategyChange = useCallback((idx, strategy) => {
    setPositionStates(prev => ({
      ...prev,
      [idx]: { ...prev[idx], strategy },
    }))
  }, [])

  const handleToggleGroup = useCallback((idx, group) => {
    setPositionStates(prev => {
      const current = prev[idx]?.selectedGroups || []
      const next = current.includes(group)
        ? current.filter(g => g !== group)
        : [...current, group]
      return { ...prev, [idx]: { ...prev[idx], selectedGroups: next } }
    })
  }, [])

  if (loading) {
    return (
      <div className="bg-white rounded-xl border border-gray-100 shadow-sm p-5">
        <div className="flex items-center gap-3 text-gray-400">
          <svg className="w-5 h-5 animate-spin" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
              d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
          </svg>
          <span className="text-sm">Analyzing scaffold...</span>
        </div>
      </div>
    )
  }

  if (error) {
    return (
      <div className="bg-white rounded-xl border border-red-100 shadow-sm p-5">
        <p className="text-xs text-red-600">{error}</p>
      </div>
    )
  }

  if (!analysis || !analysis.positions?.length) return null

  const nPositions = analysis.positions.length
  const nBrics = analysis.brics_bond_count || 0

  return (
    <div className="bg-white rounded-xl border border-gray-100 shadow-sm overflow-hidden">
      {/* Header */}
      <div className="flex items-center justify-between px-5 py-3 border-b border-gray-100">
        <h2 className="text-xs font-semibold text-gray-400 uppercase tracking-wider">Structural Analysis</h2>
        <div className="flex items-center gap-2">
          <span className="text-xs font-semibold text-[#0f131d]">{nPositions} position{nPositions !== 1 ? 's' : ''}</span>
          <span className="text-xs text-gray-300">|</span>
          <span className="text-xs font-semibold text-red-500">{nBrics} BRICS</span>
        </div>
      </div>

      {/* Main body: SVG + Positions */}
      <div className="flex flex-col md:flex-row">
        {/* Left: SVG */}
        <div className="md:w-1/2 p-4 flex items-center justify-center border-b md:border-b-0 md:border-r border-gray-100 bg-gray-50/50">
          {analysis.annotated_svg ? (
            <div
              className="max-w-full"
              dangerouslySetInnerHTML={{ __html: analysis.annotated_svg }}
            />
          ) : (
            <div className="text-gray-300 text-sm">No SVG available</div>
          )}
        </div>

        {/* Right: Position cards */}
        <div className="md:w-1/2 p-4 space-y-2 max-h-[400px] overflow-y-auto">
          {analysis.positions.map(pos => (
            <PositionCard
              key={pos.position_idx}
              position={pos}
              state={positionStates[pos.position_idx] || { frozen: false, strategy: 'any', selectedGroups: [] }}
              onToggleFreeze={handleToggleFreeze}
              onStrategyChange={handleStrategyChange}
              onToggleGroup={handleToggleGroup}
            />
          ))}
        </div>
      </div>

      {/* Global controls */}
      <div className="px-5 py-4 border-t border-gray-100 bg-gray-50/30">
        <GlobalControls
          preserveScaffold={preserveScaffold}
          setPreserveScaffold={setPreserveScaffold}
          similarity={similarity}
          setSimilarity={setSimilarity}
          maxMwChange={maxMwChange}
          setMaxMwChange={setMaxMwChange}
          globalStrategies={globalStrategies}
          setGlobalStrategies={setGlobalStrategies}
        />
      </div>
    </div>
  )
}
