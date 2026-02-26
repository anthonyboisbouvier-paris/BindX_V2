import React from 'react'
import { MOCK_PROJECTS } from '../mock/data.js'

// ---------------------------------------------------------------------------
// Relative time formatter
// ---------------------------------------------------------------------------

function relativeTime(iso) {
  if (!iso) return ''
  const diff = Date.now() - new Date(iso).getTime()
  const mins = Math.floor(diff / 60000)
  if (mins < 60) return `${mins}m ago`
  const hrs = Math.floor(mins / 60)
  if (hrs < 24) return `${hrs}h ago`
  const days = Math.floor(hrs / 24)
  return `${days}d ago`
}

// ---------------------------------------------------------------------------
// Color by activity type
// ---------------------------------------------------------------------------

function dotColor(type) {
  switch (type) {
    case 'run_completed': return 'bg-green-400'
    case 'phase_frozen': return 'bg-blue-400'
    case 'agent_insight': return 'bg-amber-400'
    case 'run_started': return 'bg-blue-300'
    case 'import': return 'bg-gray-400'
    default: return 'bg-gray-300'
  }
}

function typeLabel(type) {
  switch (type) {
    case 'run_completed': return 'Run'
    case 'phase_frozen': return 'Phase'
    case 'agent_insight': return 'Agent'
    case 'run_started': return 'Run'
    default: return ''
  }
}

function typeLabelColor(type) {
  switch (type) {
    case 'run_completed': return 'text-green-600 bg-green-50'
    case 'phase_frozen': return 'text-blue-600 bg-blue-50'
    case 'agent_insight': return 'text-amber-600 bg-amber-50'
    default: return 'text-gray-500 bg-gray-50'
  }
}

// ---------------------------------------------------------------------------
// ActivityTimeline
// ---------------------------------------------------------------------------
// Props:
//   activities — array of { id, type, message, project_id, timestamp }
//   limit      — max items to show (default: all)
//   showProject — show project name (default: true)

export default function ActivityTimeline({ activities = [], limit, showProject = true }) {
  const items = limit ? activities.slice(0, limit) : activities

  if (items.length === 0) {
    return (
      <div className="text-center py-8 text-gray-400 text-sm">
        No recent activity.
      </div>
    )
  }

  return (
    <div className="space-y-0">
      {items.map((act, idx) => {
        const isLast = idx === items.length - 1
        const project = MOCK_PROJECTS.find(p => p.id === act.project_id)
        const label = typeLabel(act.type)
        const labelColor = typeLabelColor(act.type)

        return (
          <div key={act.id} className="flex gap-3">
            {/* Left column: dot + line */}
            <div className="flex flex-col items-center">
              <div className={`w-2.5 h-2.5 rounded-full mt-1 shrink-0 ${dotColor(act.type)}`} />
              {!isLast && (
                <div className="w-px flex-1 bg-gray-100 my-1" style={{ minHeight: 20 }} />
              )}
            </div>

            {/* Right column: content */}
            <div className={`pb-4 ${isLast ? '' : ''}`}>
              <div className="flex items-center gap-2 flex-wrap">
                {label && (
                  <span className={`text-xs font-medium px-1.5 py-0.5 rounded ${labelColor}`}>
                    {label}
                  </span>
                )}
                <span className="text-xs text-gray-400">{relativeTime(act.timestamp)}</span>
              </div>
              <p className="text-sm text-gray-700 mt-0.5 leading-snug">{act.message}</p>
              {showProject && project && (
                <p className="text-xs text-gray-400 mt-0.5 font-medium">{project.name}</p>
              )}
            </div>
          </div>
        )
      })}
    </div>
  )
}
