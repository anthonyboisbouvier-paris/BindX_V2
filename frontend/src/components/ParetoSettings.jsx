import React, { useRef, useEffect } from 'react'
import useSettingsStore from '../stores/settingsStore.js'

export default function ParetoSettings({ isOpen, onClose, allObjectives = [] }) {
  const panelRef = useRef(null)
  const paretoOverrides = useSettingsStore(s => s.paretoOverrides)
  const setParetoObjective = useSettingsStore(s => s.setParetoObjective)
  const resetParetoObjective = useSettingsStore(s => s.resetParetoObjective)

  // Close on click outside
  useEffect(() => {
    if (!isOpen) return
    const handler = (e) => {
      if (panelRef.current && !panelRef.current.contains(e.target)) {
        onClose()
      }
    }
    document.addEventListener('mousedown', handler)
    return () => document.removeEventListener('mousedown', handler)
  }, [isOpen, onClose])

  if (!isOpen) return null

  const enabledCount = allObjectives.filter(o => {
    const ov = paretoOverrides[o.key]
    return ov?.enabled !== false
  }).length

  const overrideCount = Object.keys(paretoOverrides).length

  return (
    <div
      ref={panelRef}
      className="absolute right-0 top-full mt-1.5 z-50 bg-white border border-gray-200 rounded-xl shadow-xl"
      style={{ width: '340px' }}
    >
      {/* Header */}
      <div className="flex items-center justify-between px-3 pt-3 pb-2 border-b border-gray-100">
        <div className="flex items-center gap-2">
          <svg className="w-3.5 h-3.5 text-gray-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
              d="M10.325 4.317c.426-1.756 2.924-1.756 3.35 0a1.724 1.724 0 002.573 1.066c1.543-.94 3.31.826 2.37 2.37a1.724 1.724 0 001.066 2.573c1.756.426 1.756 2.924 0 3.35a1.724 1.724 0 00-1.066 2.573c.94 1.543-.826 3.31-2.37 2.37a1.724 1.724 0 00-2.573 1.066c-.426 1.756-2.924 1.756-3.35 0a1.724 1.724 0 00-2.573-1.066c-1.543.94-3.31-.826-2.37-2.37a1.724 1.724 0 00-1.066-2.573c-1.756-.426-1.756-2.924 0-3.35a1.724 1.724 0 001.066-2.573c-.94-1.543.826-3.31 2.37-2.37.996.608 2.296.07 2.572-1.065z" />
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M15 12a3 3 0 11-6 0 3 3 0 016 0z" />
          </svg>
          <span className="text-xs font-bold text-gray-700">Pareto Objectives</span>
          <span className="text-[10px] px-1.5 py-0.5 bg-gray-100 text-gray-500 rounded-full font-medium tabular-nums">
            {enabledCount}/{allObjectives.length}
          </span>
        </div>
        <button onClick={onClose} className="text-gray-400 hover:text-gray-600">
          <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
          </svg>
        </button>
      </div>

      {/* Objectives list */}
      <div className="overflow-y-auto max-h-80 py-1">
        {allObjectives.map(obj => {
          const ov = paretoOverrides[obj.key] || {}
          const enabled = ov.enabled !== false
          const direction = ov.higher_is_better ?? obj.higher_is_better
          const isOverridden = ov.enabled !== undefined || ov.higher_is_better !== undefined

          return (
            <div
              key={obj.key}
              className={`flex items-center gap-2.5 px-3 py-2 hover:bg-gray-50 transition-colors ${
                !enabled ? 'opacity-50' : ''
              }`}
            >
              {/* Enable toggle */}
              <input
                type="checkbox"
                checked={enabled}
                onChange={e => setParetoObjective(obj.key, { enabled: e.target.checked })}
                className="accent-bx-mint w-3.5 h-3.5 flex-shrink-0"
              />

              {/* Label */}
              <span className={`text-xs flex-1 ${enabled ? 'text-gray-700 font-medium' : 'text-gray-400'}`}>
                {obj.label}
              </span>

              {/* Override indicator */}
              {isOverridden && (
                <span className="w-1.5 h-1.5 rounded-full bg-bx-cyan flex-shrink-0" title="Modified from default" />
              )}

              {/* Direction toggle: Min / Max */}
              <div className="flex items-center bg-gray-100 rounded-md overflow-hidden">
                <button
                  onClick={() => setParetoObjective(obj.key, { higher_is_better: false })}
                  disabled={!enabled}
                  className={`px-2 py-0.5 text-[10px] font-semibold transition-colors ${
                    !direction
                      ? 'bg-bx-surface text-white'
                      : 'text-gray-400 hover:text-gray-600'
                  }`}
                >
                  Min
                </button>
                <button
                  onClick={() => setParetoObjective(obj.key, { higher_is_better: true })}
                  disabled={!enabled}
                  className={`px-2 py-0.5 text-[10px] font-semibold transition-colors ${
                    direction
                      ? 'bg-bx-surface text-white'
                      : 'text-gray-400 hover:text-gray-600'
                  }`}
                >
                  Max
                </button>
              </div>

              {/* Reset button */}
              {isOverridden && (
                <button
                  onClick={() => resetParetoObjective(obj.key)}
                  className="text-[9px] text-gray-400 hover:text-red-500 transition-colors flex-shrink-0"
                  title="Reset to default"
                >
                  <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                      d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
                  </svg>
                </button>
              )}
            </div>
          )
        })}
      </div>

      {/* Footer */}
      <div className="px-3 py-2 border-t border-gray-100 flex items-center justify-between bg-gray-50/50 rounded-b-xl">
        <span className="text-[10px] text-gray-400">
          {overrideCount > 0
            ? <><strong className="text-gray-600">{overrideCount}</strong> modified</>
            : 'All defaults'
          }
        </span>
        <button
          onClick={onClose}
          className="text-[10px] text-bx-light-text hover:underline font-medium"
        >
          Done
        </button>
      </div>
    </div>
  )
}
