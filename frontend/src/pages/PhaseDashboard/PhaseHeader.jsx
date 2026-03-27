import React, { useCallback, useState, useEffect, useRef } from 'react'
import { Link } from 'react-router-dom'
import FilterBar from '../../components/FilterBar.jsx'
import Badge from '../../components/Badge.jsx'
import BindXLogo from '../../components/BindXLogo.jsx'
import InfoTip, { TIPS } from '../../components/InfoTip.jsx'

function FreezeBanner() {
  return (
    <div className="flex items-center gap-2.5 rounded-[9px] px-4 py-2.5 text-sm" style={{ background: 'var(--blue-soft)', color: 'var(--blue)', border: '1px solid rgba(59,130,246,.12)' }}>
      <svg className="w-4 h-4 text-blue-500 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
          d="M12 15v2m-6 4h12a2 2 0 002-2v-6a2 2 0 00-2-2H6a2 2 0 00-2 2v6a2 2 0 002 2zm10-10V7a4 4 0 00-8 0v4h8z" />
      </svg>
      <div>
        <span className="font-semibold">Phase frozen</span>
        <span className="ml-1 font-normal">— molecules locked for downstream analysis. No new runs or edits allowed.</span>
      </div>
    </div>
  )
}

// ---------------------------------------------------------------------------

function StatsBar({ stats }) {
  const cards = [
    { value: stats.total_molecules ?? 0, label: 'Molecules', colorClass: 'text-bx-mint' },
    { value: stats.bookmarked ?? 0, label: 'Bookmarked', colorClass: 'text-bx-amber' },
    { value: stats.runs_completed ?? 0, label: 'Runs done', colorClass: 'text-bx-cyan' },
    {
      value: stats.runs_running ?? 0,
      label: 'Running',
      colorClass: stats.runs_running > 0 ? 'text-bx-blue' : 'text-bx-dim',
      animated: stats.runs_running > 0,
    },
  ]

  return (
    <div className="dark-inset">
      <div className="grid grid-cols-4 gap-4">
        {cards.map(card => (
          <div key={card.label} className="text-center py-2">
            <div className="flex items-center justify-center gap-1.5 mb-0.5">
              {card.animated && <BindXLogo variant="loading" size={14} />}
              <p className={`stat-value ${card.colorClass}`}>{card.value}</p>
            </div>
            <p className="stat-label">{card.label}</p>
          </div>
        ))}
      </div>
    </div>
  )
}


function FilterBarWithCount({ molecules, columns, onFilteredChange, onFilterCountChange, serverMode, onServerFilterChange, totalFromServer }) {
  const handleFilteredChange = useCallback((filtered) => {
    onFilteredChange(filtered)
    // Estimate filter count via difference (rough heuristic)
    onFilterCountChange(filtered.length < molecules.length ? 1 : 0)
  }, [onFilteredChange, onFilterCountChange, molecules.length])

  return (
    <FilterBar
      molecules={molecules}
      columns={columns}
      onFilteredChange={handleFilteredChange}
      serverMode={serverMode}
      onServerFilterChange={onServerFilterChange}
      totalFromServer={totalFromServer}
    />
  )
}

function Breadcrumb({ projectId, projectName, phase, phaseTypeMeta }) {
  return (
    <nav className="flex items-center gap-1.5 text-sm text-gray-400" aria-label="Breadcrumb">
      <Link to="/" className="hover:text-bx-mint transition-colors">Projects</Link>
      <svg className="w-3 h-3 text-gray-300" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
      </svg>
      <Link to={`/project/${projectId}`} className="hover:text-bx-light-text transition-colors truncate max-w-[120px]">
        {projectName || 'Project'}
      </Link>
      <svg className="w-3 h-3 text-gray-300" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
      </svg>
      <span className="text-gray-600 font-medium truncate">
        Phase {phaseTypeMeta.short || '?'} — {phaseTypeMeta.label || phase.type}
      </span>
    </nav>
  )
}

function PhaseDescription({ phaseId }) {
  const storageKey = `bindx_phase_desc_${phaseId}`
  const [value, setValue] = useState(() => localStorage.getItem(storageKey) || '')

  useEffect(() => {
    setValue(localStorage.getItem(`bindx_phase_desc_${phaseId}`) || '')
  }, [phaseId])

  const handleChange = (e) => {
    setValue(e.target.value)
    localStorage.setItem(storageKey, e.target.value)
  }

  return (
    <div className="flex flex-col h-full">
      <label className="text-[10px] uppercase tracking-wider text-gray-400 font-semibold mb-1 shrink-0">Phase Description</label>
      <textarea
        value={value}
        onChange={handleChange}
        placeholder="Goals, strategy, notes..."
        className="flex-1 text-xs text-gray-600 border border-gray-200 rounded-lg px-3 py-2 resize-none overflow-y-auto focus:outline-none focus:ring-1 focus:ring-bx-mint focus:border-bx-mint placeholder-gray-300"
      />
    </div>
  )
}

function PhaseHeader({ phase, phaseTypeMeta, isFrozen, stats, onFreezeToggle, onNewRun, hasActiveRun, onOpenAgent, onDeletePhase, onSendToNextPhase, bookmarkedCount, onOpenReport, onOpenSAR }) {
  return (
    <div className="card overflow-hidden">
      <div className={`h-1.5 ${isFrozen ? 'bg-blue-400' : 'bg-bx-mint'}`} />
      <div className="px-5 py-4">
        <div className="flex flex-wrap items-start justify-between gap-4">
          {/* Title */}
          <div className="flex items-center gap-2.5">
            <h1 className="text-xl font-bold text-bx-light-text">
              Phase {phaseTypeMeta.short || '?'}
              <span className="ml-2 text-base font-normal text-gray-400">
                — {phaseTypeMeta.label || phase.type}
                <InfoTip text={TIPS[`phase_${phase.type === 'hit_discovery' ? 'a' : phase.type === 'hit_to_lead' ? 'b' : 'c'}`]} size="xs" />
              </span>
            </h1>
            {isFrozen && <Badge variant="blue" size="sm">Frozen</Badge>}
          </div>

          {/* Action buttons */}
          <div className="flex items-center gap-2">
            <button
              onClick={onDeletePhase}
              className="flex items-center gap-1.5 px-3 py-2 rounded-lg text-sm font-medium border border-red-200 text-red-600 hover:bg-red-50 hover:border-red-300 transition-all duration-150"
              title="Delete this phase"
            >
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                  d="M19 7l-.867 12.142A2 2 0 0116.138 21H7.862a2 2 0 01-1.995-1.858L5 7m5 4v6m4-6v6m1-10V4a1 1 0 00-1-1h-4a1 1 0 00-1 1v3M4 7h16" />
              </svg>
            </button>
            {onOpenSAR && (
              <button
                onClick={onOpenSAR}
                className="flex items-center gap-1.5 px-3 py-2 rounded-lg text-sm font-medium border border-amber-200 text-amber-700 bg-amber-50 hover:bg-amber-100 transition-all duration-150"
                title="Activity Explorer — find activity cliffs"
              >
                <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                    d="M13 10V3L4 14h7v7l9-11h-7z" />
                </svg>
                Activity
              </button>
            )}
            {onOpenReport && (
              <button
                onClick={onOpenReport}
                className="flex items-center gap-1.5 px-3 py-2 rounded-lg text-sm font-medium border border-sky-200 text-sky-700 bg-sky-50 hover:bg-sky-100 transition-all duration-150"
                title="Generate PDF report"
              >
                <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                    d="M9 12h6m-6 4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
                </svg>
                Report
              </button>
            )}
            <button
              onClick={onOpenAgent}
              className="flex items-center gap-1.5 px-3 py-2 rounded-lg text-sm font-medium border border-emerald-200 text-emerald-700 bg-emerald-50 hover:bg-emerald-100 transition-all duration-150"
            >
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                  d="M9.75 17L9 20l-1 1h8l-1-1-.75-3M3 13h18M5 17H3a2 2 0 01-2-2V5a2 2 0 012-2h14a2 2 0 012 2v10a2 2 0 01-2 2h-2" />
              </svg>
              AI Agent
              <InfoTip text="AI agent analyzes your molecules and provides insights on chemical series, property trends, safety alerts, and optimization suggestions." size="xs" />
            </button>
            <button
              onClick={onFreezeToggle}
              className={`flex items-center gap-1.5 px-3 py-2 rounded-lg text-sm font-medium border transition-all duration-150 ${
                isFrozen
                  ? 'border-blue-300 text-blue-700 bg-blue-50 hover:bg-blue-100'
                  : 'border-gray-200 text-gray-600 hover:bg-gray-50 hover:border-gray-300'
              }`}
            >
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                  d={isFrozen
                    ? "M8 11V7a4 4 0 118 0m-4 8v2m-6 4h12a2 2 0 002-2v-6a2 2 0 00-2-2H6a2 2 0 00-2 2v6a2 2 0 002 2z"
                    : "M12 15v2m-6 4h12a2 2 0 002-2v-6a2 2 0 00-2-2H6a2 2 0 00-2 2v6a2 2 0 002 2zm10-10V7a4 4 0 00-8 0v4h8z"}
                />
              </svg>
              {isFrozen ? 'Unfreeze' : 'Freeze Phase'}
              <InfoTip text={TIPS.freeze} size="xs" />
            </button>
            {onSendToNextPhase && bookmarkedCount > 0 && !isFrozen && (
              <button
                onClick={onSendToNextPhase}
                className="flex items-center gap-1.5 px-3 py-2 rounded-lg text-sm font-medium border border-bx-mint/30 text-bx-mint-dim bg-bx-mint/5 hover:bg-bx-mint/10 transition-all duration-150"
              >
                <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 7l5 5m0 0l-5 5m5-5H6" />
                </svg>
                Next Phase ({bookmarkedCount})
              </button>
            )}
            <button
              onClick={onNewRun}
              disabled={isFrozen}
              className={`flex items-center gap-1.5 px-4 py-2 rounded-lg text-sm font-semibold transition-all duration-150 ${
                isFrozen
                  ? 'bg-gray-100 text-gray-400 cursor-not-allowed'
                  : hasActiveRun
                    ? 'bg-gray-100 text-gray-400 cursor-not-allowed'
                    : 'bg-bx-surface hover:bg-bx-elevated text-white shadow-sm hover:shadow-md'
              }`}
            >
              {hasActiveRun ? (
                <>
                  <BindXLogo variant="loading" size={16} />
                  Run in progress
                </>
              ) : (
                <>
                  <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 4v16m8-8H4" />
                  </svg>
                  New Run
                </>
              )}
            </button>
          </div>
        </div>

        {/* Stats + Description row */}
        <div className="flex items-stretch gap-4 mt-3">
          <StatsBar stats={stats} />
          <div className="flex-1 min-w-[200px]" style={{ height: 80 }}>
            <PhaseDescription phaseId={phase.id} />
          </div>
        </div>

        {isFrozen && (
          <div className="mt-4">
            <FreezeBanner />
          </div>
        )}
      </div>
    </div>
  )
}


export { FreezeBanner, StatsBar, FilterBarWithCount, Breadcrumb }
export default PhaseHeader
