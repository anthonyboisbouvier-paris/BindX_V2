import React from 'react'

// ---------------------------------------------------------------------------
// SelectionToolbar — smart action bar above the molecule table
// Props:
//   selectedCount        — number of currently selected molecules
//   totalCount           — total molecules in view
//   bookmarkedCount      — bookmarked molecules count
//   filteredCount        — molecules passing active filters (for "Select filtered")
//   activeFilterCount    — number of active filters
//   onSelectAll          — select all molecules
//   onSelectNone         — clear selection
//   onSelectBookmarked   — select bookmarked
//   onSelectFiltered     — select filtered (if filters active)
//   onExport             — export callback(type: 'csv_visible'|'csv_all'|'sdf'|'pdf')
//   onBookmarkSelected   — bookmark currently selected
//   onSendToNextPhase    — send bookmarked molecules to next phase
//   isFrozen             — phase is frozen
// ---------------------------------------------------------------------------
export default function SelectionToolbar({
  selectedCount = 0,
  totalCount = 0,
  bookmarkedCount = 0,
  filteredCount = null,
  activeFilterCount = 0,
  onSelectAll,
  onSelectNone,
  onSelectBookmarked,
  onSelectFiltered,
  onExport,
  onBookmarkSelected,
  onSendToNextPhase,
  isFrozen = false,
}) {
  return (
    <div className="flex flex-wrap items-center gap-3 card px-4 py-2.5">
      {/* Left — selection state + quick select */}
      <div className="flex items-center gap-2 flex-wrap">
        {/* Selection pill */}
        <div
          className={`flex items-center gap-1.5 px-2.5 py-1 rounded-full text-sm font-semibold border transition-colors ${
            selectedCount > 0
              ? 'bg-bx-surface text-white border-bx-surface'
              : 'bg-gray-100 text-gray-500 border-gray-200'
          }`}
        >
          {selectedCount > 0 ? (
            <>
              <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M5 13l4 4L19 7" />
              </svg>
              {selectedCount} selected
            </>
          ) : (
            'No selection'
          )}
        </div>

        <span className="text-gray-200 text-sm">|</span>

        {/* Quick select buttons */}
        <button
          onClick={onSelectAll}
          className="text-sm text-gray-500 hover:text-bx-light-text transition-colors hover:underline underline-offset-2"
        >
          All ({totalCount})
        </button>

        {bookmarkedCount > 0 && (
          <button
            onClick={onSelectBookmarked}
            className="text-sm text-gray-500 hover:text-yellow-600 transition-colors hover:underline underline-offset-2 flex items-center gap-0.5"
          >
            <svg className="w-3 h-3 text-yellow-400" fill="currentColor" viewBox="0 0 24 24">
              <path d="M5 5a2 2 0 012-2h10a2 2 0 012 2v16l-7-3.5L5 21V5z" />
            </svg>
            Bookmarked ({bookmarkedCount})
          </button>
        )}

        {activeFilterCount > 0 && filteredCount != null && (
          <button
            onClick={onSelectFiltered}
            className="text-sm text-blue-600 hover:text-bx-light-text transition-colors hover:underline underline-offset-2 flex items-center gap-0.5"
          >
            Filtered ({filteredCount})
          </button>
        )}

        {selectedCount > 0 && (
          <button
            onClick={onSelectNone}
            className="text-sm text-gray-400 hover:text-red-500 transition-colors flex items-center gap-0.5"
          >
            <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M6 18L18 6M6 6l12 12" />
            </svg>
            Clear
          </button>
        )}

        {activeFilterCount > 0 && (
          <span className="flex items-center gap-1 px-2 py-0.5 bg-blue-50 text-blue-700 text-[10px] font-semibold rounded-full border border-blue-200">
            <svg className="w-2.5 h-2.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                d="M3 4a1 1 0 011-1h16a1 1 0 011 1v2a1 1 0 01-.293.707L13 13.414V19a1 1 0 01-.553.894l-4 2A1 1 0 017 21v-7.586L3.293 6.707A1 1 0 013 6V4z" />
            </svg>
            {activeFilterCount} filter{activeFilterCount !== 1 ? 's' : ''} active
          </span>
        )}
      </div>

      {/* Spacer */}
      <div className="flex-1" />

      {/* Right — action buttons */}
      <div className="flex items-center gap-2">
        {isFrozen && (
          <span className="flex items-center gap-1.5 px-2.5 py-1 bg-blue-50 text-blue-600 text-sm font-medium rounded-full border border-blue-200">
            <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                d="M12 15v2m-6 4h12a2 2 0 002-2v-6a2 2 0 00-2-2H6a2 2 0 00-2 2v6a2 2 0 002 2zm10-10V7a4 4 0 00-8 0v4h8z" />
            </svg>
            Phase frozen
          </span>
        )}

        {/* Export */}
        <button
          onClick={onExport}
          className="flex items-center gap-1.5 px-3 py-1.5 rounded-lg text-sm font-medium border border-gray-200 text-gray-700 hover:bg-gray-50 hover:border-gray-300 transition-all duration-150"
        >
          <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
              d="M12 10v6m0 0l-3-3m3 3l3-3m2 8H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
          </svg>
          Export
        </button>

        {/* Bookmark selected */}
        <div title={selectedCount === 0 ? 'Select molecules first' : isFrozen ? 'Phase frozen' : 'Bookmark selected'}>
          <button
            onClick={selectedCount > 0 && !isFrozen ? onBookmarkSelected : undefined}
            disabled={selectedCount === 0 || isFrozen}
            className={`flex items-center gap-1.5 px-3 py-1.5 rounded-lg text-sm font-medium border transition-all duration-150 ${
              selectedCount === 0 || isFrozen
                ? 'border-gray-200 text-gray-300 cursor-not-allowed bg-gray-50'
                : 'border-yellow-300 text-yellow-700 hover:bg-yellow-50 hover:border-yellow-400'
            }`}
          >
            <svg className="w-3.5 h-3.5" fill={selectedCount > 0 && !isFrozen ? 'none' : 'none'} stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                d="M5 5a2 2 0 012-2h10a2 2 0 012 2v16l-7-3.5L5 21V5z" />
            </svg>
            Bookmark
            {selectedCount > 0 && (
              <span className="tabular-nums text-sm opacity-70">({selectedCount})</span>
            )}
          </button>
        </div>

        {/* Send bookmarks to next phase */}
        {onSendToNextPhase && bookmarkedCount > 0 && (
          <button
            onClick={onSendToNextPhase}
            className="flex items-center gap-1.5 px-3 py-1.5 rounded-lg text-sm font-medium border border-emerald-300 text-emerald-700 bg-emerald-50 hover:bg-emerald-100 hover:border-emerald-400 transition-all duration-150"
            title={`Send ${bookmarkedCount} bookmarked molecules to the next phase`}
          >
            <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                d="M13 7l5 5m0 0l-5 5m5-5H6" />
            </svg>
            Next Phase
            <span className="tabular-nums text-sm opacity-70">({bookmarkedCount})</span>
          </button>
        )}
      </div>
    </div>
  )
}
