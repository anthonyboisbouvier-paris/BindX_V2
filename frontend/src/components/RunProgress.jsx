import React from 'react'

// ---------------------------------------------------------------------------
// RunProgress — persistent banner showing active run status
// ---------------------------------------------------------------------------
export default function RunProgress({ run, onCancel }) {
  if (!run) return null

  const isCreated = run.status === 'created'
  const progress = run.progress ?? 0

  const typeLabel = run.type === 'calculation' && run.calculation_types?.length
    ? run.calculation_types.map(t => t.charAt(0).toUpperCase() + t.slice(1)).join(' + ')
    : run.type.charAt(0).toUpperCase() + run.type.slice(1)

  return (
    <div className="bg-gradient-to-r from-blue-50 to-indigo-50 border border-blue-200 rounded-xl px-4 py-3 shadow-sm animate-in">
      <div className="flex items-center justify-between gap-4">
        {/* Left: status + type */}
        <div className="flex items-center gap-3 min-w-0">
          {/* Spinner */}
          <div className="flex-shrink-0">
            <svg className="w-5 h-5 text-blue-500 animate-spin" fill="none" viewBox="0 0 24 24">
              <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
              <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
            </svg>
          </div>

          <div className="min-w-0">
            <div className="flex items-center gap-2">
              <span className="text-sm font-semibold text-blue-800">
                {isCreated ? 'Starting' : 'Running'}: {typeLabel}
              </span>
              {run.current_step && (
                <span className="text-xs text-blue-500 truncate">
                  — {run.current_step}
                </span>
              )}
            </div>

            {/* Progress bar */}
            <div className="flex items-center gap-2 mt-1.5">
              <div className="flex-1 bg-blue-100 rounded-full h-1.5 overflow-hidden min-w-[120px]">
                <div
                  className="h-full bg-blue-500 rounded-full transition-all duration-700 ease-out"
                  style={{ width: `${isCreated ? 0 : progress}%` }}
                />
              </div>
              <span className="text-xs font-medium text-blue-600 tabular-nums flex-shrink-0">
                {isCreated ? 'Queued' : `${progress}%`}
              </span>
            </div>
          </div>
        </div>

        {/* Right: cancel button */}
        {onCancel && (
          <button
            onClick={() => onCancel(run.id)}
            className="flex-shrink-0 flex items-center gap-1.5 px-3 py-1.5 text-xs font-medium text-red-600 bg-red-50 hover:bg-red-100 border border-red-200 rounded-lg transition-colors"
          >
            <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
            </svg>
            Cancel
          </button>
        )}
      </div>
    </div>
  )
}
