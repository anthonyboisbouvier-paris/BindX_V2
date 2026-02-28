import React from 'react'

// ---------------------------------------------------------------------------
// FreezeDialog — confirmation modal for freeze / unfreeze phase actions
//
// Props:
//   isOpen            — boolean
//   onClose           — () => void
//   onConfirm         — () => void  (called before onClose)
//   action            — 'freeze' | 'unfreeze'
//   phaseName         — string  e.g. "Phase A"
//   bookmarkedMols    — array of { id, name } for preview list
//   hasDownstreamRuns — boolean  (true = risky unfreeze)
//   downstreamCount   — number   (# downstream runs/phases)
// ---------------------------------------------------------------------------
export default function FreezeDialog({
  isOpen,
  onClose,
  onConfirm,
  action = 'freeze',
  phaseName = 'Phase',
  bookmarkedMols = [],
  hasDownstreamRuns = false,
  downstreamCount = 0,
}) {
  if (!isOpen) return null

  const isFreeze = action === 'freeze'
  const isRiskyUnfreeze = !isFreeze && hasDownstreamRuns

  function handleConfirm() {
    onConfirm?.()
    onClose?.()
  }

  // -------------------------------------------------------------------------
  // Header
  // -------------------------------------------------------------------------
  const headerBg = isRiskyUnfreeze
    ? 'bg-red-600'
    : isFreeze
    ? 'bg-bx-surface'
    : 'bg-gray-700'

  // -------------------------------------------------------------------------
  // Confirm button
  // -------------------------------------------------------------------------
  const confirmClass = isRiskyUnfreeze
    ? 'bg-red-600 hover:bg-red-700 text-white'
    : isFreeze
    ? 'bg-bx-surface hover:bg-bx-elevated text-white'
    : 'bg-gray-600 hover:bg-gray-700 text-white'

  const confirmLabel = isFreeze
    ? 'Freeze Phase'
    : isRiskyUnfreeze
    ? 'Unfreeze Anyway'
    : 'Unfreeze'

  return (
    <div
      className="fixed inset-0 z-50 flex items-center justify-center bg-black/40 backdrop-blur-sm p-4"
      onClick={onClose}
    >
      <div
        className="bg-white rounded-2xl shadow-2xl max-w-md w-full overflow-hidden"
        onClick={e => e.stopPropagation()}
      >
        {/* ----------------------------------------------------------------- */}
        {/* Colored header bar                                                  */}
        {/* ----------------------------------------------------------------- */}
        <div className={`${headerBg} px-6 py-5`}>
          <div className="flex items-center gap-4">
            {/* Icon */}
            <div className="flex-shrink-0 w-12 h-12 rounded-xl bg-white/10 flex items-center justify-center">
              {isFreeze ? (
                <svg className="w-7 h-7 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                    d="M16.5 10.5V6.75a4.5 4.5 0 10-9 0v3.75m-.75 11.25h10.5a2.25 2.25 0 002.25-2.25v-6.75a2.25 2.25 0 00-2.25-2.25H6.75a2.25 2.25 0 00-2.25 2.25v6.75a2.25 2.25 0 002.25 2.25z" />
                </svg>
              ) : isRiskyUnfreeze ? (
                <svg className="w-7 h-7 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                    d="M12 9v3.75m-9.303 3.376c-.866 1.5.217 3.374 1.948 3.374h14.71c1.73 0 2.813-1.874 1.948-3.374L13.949 3.378c-.866-1.5-3.032-1.5-3.898 0L2.697 16.126zM12 15.75h.007v.008H12v-.008z" />
                </svg>
              ) : (
                <svg className="w-7 h-7 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                    d="M13.5 10.5V6.75a4.5 4.5 0 119 0v3.75M3.75 21.75h10.5a2.25 2.25 0 002.25-2.25v-6.75a2.25 2.25 0 00-2.25-2.25H3.75a2.25 2.25 0 00-2.25 2.25v6.75a2.25 2.25 0 002.25 2.25z" />
                </svg>
              )}
            </div>

            <div>
              <h2 className="text-white font-bold text-lg leading-tight">
                {isFreeze ? `Freeze ${phaseName}` : `Unfreeze ${phaseName}`}
                {isFreeze ? ' — Hit Discovery?' : '?'}
              </h2>
              <p className="text-white/60 text-sm mt-0.5">
                {isFreeze
                  ? 'Lock the bookmarked selection for the next phase'
                  : isRiskyUnfreeze
                  ? 'This action may impact downstream data'
                  : 'Unlock this phase for editing'
                }
              </p>
            </div>
          </div>
        </div>

        {/* ----------------------------------------------------------------- */}
        {/* Body                                                                */}
        {/* ----------------------------------------------------------------- */}
        <div className="px-6 py-5 space-y-4">

          {/* Freeze: description + molecule preview */}
          {isFreeze && (
            <>
              <p className="text-sm text-gray-600 leading-relaxed">
                Freezing locks the current bookmarked selection
                <strong> ({bookmarkedMols.length} {bookmarkedMols.length === 1 ? 'molecule' : 'molecules'})</strong>.
                These will be available as input for the next phase. No further changes
                can be made to this phase while frozen.
              </p>

              {bookmarkedMols.length > 0 && (
                <div className="bg-blue-50 border border-blue-100 rounded-xl p-3">
                  <p className="text-[10px] font-semibold text-blue-500 uppercase tracking-wide mb-2">
                    Bookmarked molecules to lock
                  </p>
                  <div className="space-y-1">
                    {bookmarkedMols.slice(0, 6).map((mol, i) => (
                      <div key={mol.id || i} className="flex items-center gap-2 text-xs text-blue-700">
                        <svg className="w-3.5 h-3.5 text-blue-400 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                            d="M5 5a2 2 0 012-2h10a2 2 0 012 2v16l-7-3.5L5 21V5z" />
                        </svg>
                        <span className="font-medium">{mol.name || mol.id}</span>
                      </div>
                    ))}
                    {bookmarkedMols.length > 6 && (
                      <p className="text-xs text-blue-400 pl-5">
                        +{bookmarkedMols.length - 6} more
                      </p>
                    )}
                  </div>
                </div>
              )}

              {bookmarkedMols.length === 0 && (
                <div className="bg-amber-50 border border-amber-100 rounded-xl p-3 flex items-start gap-2">
                  <svg className="w-4 h-4 text-amber-500 flex-shrink-0 mt-0.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                      d="M12 9v3.75m-9.303 3.376c-.866 1.5.217 3.374 1.948 3.374h14.71c1.73 0 2.813-1.874 1.948-3.374L13.949 3.378c-.866-1.5-3.032-1.5-3.898 0L2.697 16.126zM12 15.75h.007v.008H12v-.008z" />
                  </svg>
                  <p className="text-xs text-amber-700">
                    No bookmarked molecules. Freeze anyway? The phase will lock with an empty selection.
                  </p>
                </div>
              )}
            </>
          )}

          {/* Normal unfreeze */}
          {!isFreeze && !hasDownstreamRuns && (
            <p className="text-sm text-gray-600 leading-relaxed">
              <strong>{phaseName}</strong> will be unlocked. You can add or modify molecules
              and re-run analyses. Previous results will be preserved.
            </p>
          )}

          {/* Risky unfreeze: warning */}
          {isRiskyUnfreeze && (
            <>
              <p className="text-sm text-gray-600 leading-relaxed">
                <strong>{phaseName}</strong> has downstream phases that already use its frozen
                selection.{downstreamCount > 0 && ` ${downstreamCount} run${downstreamCount > 1 ? 's' : ''} depend on this data.`}
              </p>

              <div className="bg-red-50 border border-red-200 rounded-xl p-4 space-y-2">
                <div className="flex items-center gap-2">
                  <svg className="w-4 h-4 text-red-500 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                      d="M12 9v3.75m-9.303 3.376c-.866 1.5.217 3.374 1.948 3.374h14.71c1.73 0 2.813-1.874 1.948-3.374L13.949 3.378c-.866-1.5-3.032-1.5-3.898 0L2.697 16.126zM12 15.75h.007v.008H12v-.008z" />
                  </svg>
                  <p className="text-xs font-bold text-red-700 uppercase tracking-wide">
                    Warning: Destructive action
                  </p>
                </div>
                <ul className="text-xs text-red-600 space-y-1 pl-1">
                  <li className="flex items-start gap-1.5">
                    <span className="mt-0.5 flex-shrink-0">--</span>
                    Downstream phase molecules may lose their source reference
                  </li>
                  <li className="flex items-start gap-1.5">
                    <span className="mt-0.5 flex-shrink-0">--</span>
                    Scoring and ADMET data in subsequent phases may be invalidated
                  </li>
                  <li className="flex items-start gap-1.5">
                    <span className="mt-0.5 flex-shrink-0">--</span>
                    Consider creating a new phase branch instead of unfreezing
                  </li>
                </ul>
              </div>
            </>
          )}
        </div>

        {/* ----------------------------------------------------------------- */}
        {/* Footer                                                              */}
        {/* ----------------------------------------------------------------- */}
        <div className="px-6 py-4 border-t border-gray-100 bg-gray-50 flex justify-end gap-3">
          <button
            onClick={onClose}
            className="px-4 py-2 rounded-xl text-sm font-semibold border border-gray-200
                       text-gray-700 hover:bg-gray-100 transition-colors"
          >
            Cancel
          </button>
          <button
            onClick={handleConfirm}
            className={`px-5 py-2 rounded-xl text-sm font-semibold transition-colors ${confirmClass}`}
          >
            {confirmLabel}
          </button>
        </div>
      </div>
    </div>
  )
}
