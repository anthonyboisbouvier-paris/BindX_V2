import React, { useState } from 'react'

// ---------------------------------------------------------------------------
// Run type icons
// ---------------------------------------------------------------------------
function RunTypeIcon({ type, className = 'w-4 h-4' }) {
  const props = { className, fill: 'none', stroke: 'currentColor', viewBox: '0 0 24 24' }
  switch (type) {
    case 'import':
      return (
        <svg {...props}>
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
            d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-4l-4 4m0 0l-4-4m4 4V4" />
        </svg>
      )
    case 'docking':
      return (
        <svg {...props}>
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
            d="M7 16a4 4 0 01-.88-7.903A5 5 0 1115.9 6L16 6a5 5 0 011 9.9M15 13l-3-3m0 0l-3 3m3-3v12" />
        </svg>
      )
    case 'admet':
      return (
        <svg {...props}>
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
            d="M9 12l2 2 4-4m5.618-4.016A11.955 11.955 0 0112 2.944a11.955 11.955 0 01-8.618 3.04A12.02 12.02 0 003 9c0 5.591 3.824 10.29 9 11.622 5.176-1.332 9-6.03 9-11.622 0-1.042-.133-2.052-.382-3.016z" />
        </svg>
      )
    case 'scoring':
      return (
        <svg {...props}>
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
            d="M11.049 2.927c.3-.921 1.603-.921 1.902 0l1.519 4.674a1 1 0 00.95.69h4.915c.969 0 1.371 1.24.588 1.81l-3.976 2.888a1 1 0 00-.363 1.118l1.518 4.674c.3.922-.755 1.688-1.538 1.118l-3.976-2.888a1 1 0 00-1.176 0l-3.976 2.888c-.783.57-1.838-.197-1.538-1.118l1.518-4.674a1 1 0 00-.363-1.118l-3.976-2.888c-.784-.57-.38-1.81.588-1.81h4.914a1 1 0 00.951-.69l1.519-4.674z" />
        </svg>
      )
    case 'enrichment':
      return (
        <svg {...props}>
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
            d="M5 3v4M3 5h4M6 17v4m-2-2h4m5-16l2.286 6.857L21 12l-5.714 2.143L13 21l-2.286-6.857L5 12l5.714-2.143L13 3z" />
        </svg>
      )
    case 'generation':
      return (
        <svg {...props}>
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
            d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
        </svg>
      )
    case 'clustering':
      return (
        <svg {...props}>
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
            d="M4 6a2 2 0 012-2h2a2 2 0 012 2v2a2 2 0 01-2 2H6a2 2 0 01-2-2V6zM14 6a2 2 0 012-2h2a2 2 0 012 2v2a2 2 0 01-2 2h-2a2 2 0 01-2-2V6zM4 16a2 2 0 012-2h2a2 2 0 012 2v2a2 2 0 01-2 2H6a2 2 0 01-2-2v-2zM14 16a2 2 0 012-2h2a2 2 0 012 2v2a2 2 0 01-2 2h-2a2 2 0 01-2-2v-2z" />
        </svg>
      )
    default:
      return (
        <svg {...props}>
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
            d="M14.752 11.168l-3.197-2.132A1 1 0 0010 9.87v4.263a1 1 0 001.555.832l3.197-2.132a1 1 0 000-1.664z" />
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
            d="M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
        </svg>
      )
  }
}

// ---------------------------------------------------------------------------
// Duration formatter
// ---------------------------------------------------------------------------
function formatDuration(started, completed) {
  if (!started) return null
  if (!completed) return 'running...'
  const ms = new Date(completed) - new Date(started)
  if (isNaN(ms) || ms < 0) return null
  const secs = Math.floor(ms / 1000)
  if (secs < 60) return `${secs}s`
  const mins = Math.floor(secs / 60)
  const rem = secs % 60
  return rem > 0 ? `${mins}m ${rem}s` : `${mins}m`
}

// ---------------------------------------------------------------------------
// Config summary one-liner
// ---------------------------------------------------------------------------
function configSummary(type, config) {
  if (!config) return null
  switch (type) {
    case 'import':
      if (config.source === 'phase_selection')
        return `From ${config.source_phase_id || 'phase'}, ${config.selection || 'all'}`
      return `${config.source || 'file'} · ${config.max_molecules ?? '?'} mols`
    case 'docking':
      return `${(config.engine || 'gnina').toUpperCase()} · exh=${config.exhaustiveness ?? 32}`
    case 'admet':
      if (config.properties) {
        const shown = config.properties.slice(0, 3).join(', ')
        return shown + (config.properties.length > 3 ? ' +' + (config.properties.length - 3) : '')
      }
      return 'All properties'
    case 'scoring':
      return 'Weighted composite'
    case 'enrichment':
      return config.analyses ? config.analyses.join(', ') : 'ProLIF · Clustering'
    case 'generation':
      return `${config.method || 'Scaffold hopping'} · ${config.iterations ?? 3} iter`
    case 'clustering':
      return `Butina · cutoff=${config.cutoff ?? 0.5}`
    default:
      return null
  }
}

// ---------------------------------------------------------------------------
// Status node colors and icons
// ---------------------------------------------------------------------------
const STATUS_CONFIG = {
  completed: {
    ring: 'ring-2 ring-green-400',
    bg: 'bg-green-500',
    icon: (
      <svg className="w-4 h-4 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M5 13l4 4L19 7" />
      </svg>
    ),
    label: 'Completed',
    labelColor: 'text-green-600',
    connectorColor: 'bg-green-200',
  },
  running: {
    ring: 'ring-2 ring-blue-400 ring-offset-1',
    bg: 'bg-blue-500',
    icon: (
      <svg className="w-4 h-4 text-white animate-spin" fill="none" viewBox="0 0 24 24">
        <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
        <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
      </svg>
    ),
    label: 'Running',
    labelColor: 'text-blue-600',
    connectorColor: 'bg-blue-200',
  },
  failed: {
    ring: 'ring-2 ring-red-400',
    bg: 'bg-red-500',
    icon: (
      <svg className="w-4 h-4 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M6 18L18 6M6 6l12 12" />
      </svg>
    ),
    label: 'Failed',
    labelColor: 'text-red-500',
    connectorColor: 'bg-gray-200',
  },
  queued: {
    ring: 'ring-2 ring-gray-300',
    bg: 'bg-gray-300',
    icon: (
      <svg className="w-4 h-4 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
          d="M12 8v4l3 3m6-3a9 9 0 11-18 0 9 9 0 0118 0z" />
      </svg>
    ),
    label: 'Queued',
    labelColor: 'text-gray-400',
    connectorColor: 'bg-gray-100',
  },
}

// ---------------------------------------------------------------------------
// RunHistory — horizontal stepper timeline
// ---------------------------------------------------------------------------
export default function RunHistory({ runs = [] }) {
  const [activeRunId, setActiveRunId] = useState(null)

  if (!runs.length) {
    return (
      <div className="bg-white rounded-xl border border-gray-100 shadow-sm p-8 text-center">
        <div className="flex flex-col items-center gap-3">
          <div className="w-12 h-12 rounded-full bg-gray-100 flex items-center justify-center">
            <svg className="w-6 h-6 text-gray-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                d="M14.752 11.168l-3.197-2.132A1 1 0 0010 9.87v4.263a1 1 0 001.555.832l3.197-2.132a1 1 0 000-1.664z" />
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                d="M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
            </svg>
          </div>
          <p className="text-sm font-medium text-gray-600">No runs yet</p>
          <p className="text-xs text-gray-400">Create your first run to start analyzing molecules</p>
        </div>
      </div>
    )
  }

  const completedCount = runs.filter(r => r.status === 'completed').length
  const runningCount = runs.filter(r => r.status === 'running').length
  const failedCount = runs.filter(r => r.status === 'failed').length

  return (
    <div className="bg-white rounded-xl border border-gray-100 shadow-sm overflow-hidden">
      {/* Header */}
      <div className="px-5 py-3 border-b border-gray-100 bg-gray-50 flex items-center justify-between">
        <div className="flex items-center gap-2">
          <svg className="w-3.5 h-3.5 text-gray-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
              d="M13 10V3L4 14h7v7l9-11h-7z" />
          </svg>
          <h3 className="text-xs font-semibold text-gray-500 uppercase tracking-wide">Run History</h3>
        </div>
        <div className="flex items-center gap-3 text-xs">
          {completedCount > 0 && (
            <span className="flex items-center gap-1 text-green-600 font-medium">
              <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M5 13l4 4L19 7" />
              </svg>
              {completedCount} done
            </span>
          )}
          {runningCount > 0 && (
            <span className="flex items-center gap-1 text-blue-600 font-medium">
              <svg className="w-3 h-3 animate-spin" fill="none" viewBox="0 0 24 24">
                <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
                <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
              </svg>
              {runningCount} running
            </span>
          )}
          {failedCount > 0 && (
            <span className="flex items-center gap-1 text-red-500 font-medium">
              <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M6 18L18 6M6 6l12 12" />
              </svg>
              {failedCount} failed
            </span>
          )}
          <span className="text-gray-300">|</span>
          <span className="text-gray-400">{runs.length} total</span>
        </div>
      </div>

      {/* Timeline */}
      <div className="px-5 py-5 overflow-x-auto" style={{ scrollbarWidth: 'thin' }}>
        <div className="flex items-start gap-0 min-w-max">
          {runs.map((run, i) => {
            const sc = STATUS_CONFIG[run.status] || STATUS_CONFIG.queued
            const duration = formatDuration(run.started_at, run.completed_at)
            const summary = configSummary(run.type, run.config)
            const typeLabel = run.type.charAt(0).toUpperCase() + run.type.slice(1)
            const isActive = activeRunId === run.id
            const isLast = i === runs.length - 1

            return (
              <div key={run.id} className="flex items-start">
                {/* Node */}
                <div className="flex flex-col items-center" style={{ width: 140 }}>
                  {/* Step indicator + circle + type icon */}
                  <button
                    onClick={() => setActiveRunId(isActive ? null : run.id)}
                    className={`flex flex-col items-center gap-2 w-full px-2 py-2 rounded-xl transition-all duration-150 ${
                      isActive ? 'bg-blue-50' : 'hover:bg-gray-50'
                    }`}
                  >
                    {/* Circle */}
                    <div className={`relative w-10 h-10 rounded-full flex items-center justify-center ${sc.bg} ${sc.ring} transition-all duration-200`}>
                      <RunTypeIcon type={run.type} className="w-4.5 h-4.5 text-white" />
                      {run.status === 'running' && (
                        <span className="absolute -top-0.5 -right-0.5 w-3 h-3 bg-blue-400 rounded-full border-2 border-white animate-pulse" />
                      )}
                    </div>

                    {/* Run label */}
                    <div className="text-center">
                      <p className="text-[10px] text-gray-400 tabular-nums">#{i + 1}</p>
                      <p className="text-xs font-semibold text-gray-700 leading-tight">{typeLabel}</p>
                    </div>

                    {/* Status */}
                    <span className={`text-[10px] font-semibold ${sc.labelColor}`}>
                      {run.status === 'running' && run.progress != null
                        ? `${run.progress}%`
                        : sc.label}
                    </span>
                  </button>

                  {/* Expanded detail below node */}
                  {isActive && (
                    <div className="w-full mt-1 px-2 py-2 bg-blue-50 rounded-xl border border-blue-100 space-y-1.5">
                      {summary && (
                        <p className="text-[10px] text-gray-500 text-center">{summary}</p>
                      )}
                      {duration && (
                        <p className="text-[10px] text-gray-400 text-center flex items-center justify-center gap-1">
                          <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                              d="M12 8v4l3 3m6-3a9 9 0 11-18 0 9 9 0 0118 0z" />
                          </svg>
                          {duration}
                        </p>
                      )}
                      {run.status === 'running' && run.progress != null && (
                        <div className="w-full bg-blue-100 rounded-full h-1.5 overflow-hidden">
                          <div
                            className="h-full bg-blue-500 rounded-full transition-all duration-700"
                            style={{ width: `${run.progress}%` }}
                          />
                        </div>
                      )}
                      {run.error_message && (
                        <p className="text-[10px] text-red-500 text-center truncate" title={run.error_message}>
                          {run.error_message}
                        </p>
                      )}
                    </div>
                  )}

                  {/* Summary always shown below */}
                  {!isActive && summary && (
                    <p className="text-[9px] text-gray-400 text-center px-1 mt-1 leading-snug">
                      {summary}
                    </p>
                  )}
                  {!isActive && duration && (
                    <p className="text-[9px] text-gray-400 text-center tabular-nums">{duration}</p>
                  )}

                  {/* Progress bar always visible for running */}
                  {run.status === 'running' && run.progress != null && !isActive && (
                    <div className="w-full mt-1.5 px-2">
                      <div className="w-full bg-blue-100 rounded-full h-1 overflow-hidden">
                        <div
                          className="h-full bg-blue-500 rounded-full transition-all duration-700"
                          style={{ width: `${run.progress}%` }}
                        />
                      </div>
                    </div>
                  )}
                </div>

                {/* Connector line between nodes */}
                {!isLast && (
                  <div className="flex-shrink-0 flex items-center self-center mt-0" style={{ height: 2, width: 24 }}>
                    <div
                      className={`h-0.5 w-full ${
                        run.status === 'completed' ? 'bg-green-300' :
                        run.status === 'running' ? 'bg-blue-200' :
                        'bg-gray-200'
                      } ${runs[i + 1]?.status === 'queued' ? 'border-dashed' : ''}`}
                      style={runs[i + 1]?.status === 'queued' ? {
                        borderTop: '2px dashed #e5e7eb',
                        background: 'none',
                        height: 0,
                      } : {}}
                    />
                    {/* Arrow head */}
                    <svg className={`w-2 h-2 -ml-1 flex-shrink-0 ${
                      run.status === 'completed' ? 'text-green-400' :
                      run.status === 'running' ? 'text-blue-300' :
                      'text-gray-300'
                    }`} fill="currentColor" viewBox="0 0 8 8">
                      <polygon points="0,0 8,4 0,8" />
                    </svg>
                  </div>
                )}
              </div>
            )
          })}
        </div>
      </div>
    </div>
  )
}
