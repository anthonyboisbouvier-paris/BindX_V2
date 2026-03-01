import React, { useEffect, useState, useMemo, useCallback, useRef } from 'react'
import { useParams, Link, useNavigate } from 'react-router-dom'

import { useWorkspace } from '../../contexts/WorkspaceContext.jsx'
import { useToast } from '../../contexts/ToastContext.jsx'
import { ALL_COLUMNS, COLUMN_PRESETS, PHASE_TYPES, flattenMoleculeProperties, detectAvailableColumns } from '../../lib/columns.js'

import MoleculeTable from '../../components/MoleculeTable.jsx'
import SelectionToolbar from '../../components/SelectionToolbar.jsx'
import ColumnSelector from '../../components/ColumnSelector.jsx'
import RunHistory from '../../components/RunHistory.jsx'
import RunCreator from '../../components/RunCreator.jsx'
import RunProgress from '../../components/RunProgress.jsx'
import FreezeDialog from '../../components/FreezeDialog.jsx'
import ParetoFront from '../../components/ParetoFront.jsx'
import ExportModal from '../../components/ExportModal.jsx'
import CampaignAgentPanel from '../../components/CampaignAgentPanel.jsx'
import InfoTip, { TIPS } from '../../components/InfoTip.jsx'
import BindXLogo from '../../components/BindXLogo.jsx'

import AdvancedFilterBuilder from '../../components/AdvancedFilterBuilder.jsx'
import MoleculeDetailPanel, { DetailPopupModal } from './MoleculeDetailPanel.jsx'
import PhaseHeader, { Breadcrumb, FilterBarWithCount, StatsBar } from './PhaseHeader.jsx'

export default function PhaseDashboard() {
  const { projectId, phaseId } = useParams()

  const {
    currentProject,
    currentCampaign,
    currentPhase,
    phaseMolecules,
    currentPhaseRuns,
    runsLoading,
    moleculesLoading,
    selectProject,
    selectPhase,
    selectedMoleculeIds,
    toggleSelection,
    selectAll,
    selectNone,
    selectBookmarked,
    toggleBookmark,
    bookmarkSelected,
    freezePhase,
    unfreezePhase,
    getPhaseStatus,
    createRun,
    cancelRun,
    archiveRun,
    importFile,
    importDatabase,
    moleculeTotal,
    moleculeHasMore,
    moleculeStats,
    loadMoreMolecules,
    updateServerSort,
    updateServerFilters,
  } = useWorkspace()

  const { addToast } = useToast()

  // Sync URL params to workspace state
  useEffect(() => { if (projectId) selectProject(projectId) }, [projectId, selectProject])
  useEffect(() => { if (phaseId) selectPhase(phaseId) }, [phaseId, selectPhase])

  // --- Flatten molecule properties from API (EAV → flat columns) ---
  const flatMolecules = useMemo(
    () => phaseMolecules.map(flattenMoleculeProperties),
    [phaseMolecules]
  )

  // --- Detect which columns have data ---
  const availableColumns = useMemo(
    () => detectAvailableColumns(flatMolecules),
    [flatMolecules]
  )

  // --- User annotations (tags, invalidation, comments) — local state keyed by mol id ---
  const [annotations, setAnnotations] = useState({}) // { molId: { tags: [...], invalidated: bool, user_comment: '...' } }

  const handleAnnotation = useCallback((molId, key, value) => {
    setAnnotations(prev => ({
      ...prev,
      [molId]: { ...(prev[molId] || {}), [key]: value },
    }))
  }, [])

  // Merge annotations into flat molecules
  const annotatedMolecules = useMemo(
    () => flatMolecules.map(m => ({ ...m, ...(annotations[m.id] || {}) })),
    [flatMolecules, annotations]
  )

  // --- Local UI state ---
  const [selectedMolecule, setSelectedMolecule] = useState(null)
  const [showRunCreator, setShowRunCreator] = useState(false)
  const [showPareto, setShowPareto] = useState(true)
  const [showFreezeDialog, setShowFreezeDialog] = useState(false)
  const [filteredMolecules, setFilteredMolecules] = useState(annotatedMolecules)
  const [popupState, setPopupState] = useState(null) // { type: 'safety'|'retrosynthesis'|'confidence', molecule }
  const [showExportModal, setShowExportModal] = useState(false)
  const [showAgentPanel, setShowAgentPanel] = useState(false)
  const [showRunLogs, setShowRunLogs] = useState(false)
  const [advancedFilterFn, setAdvancedFilterFn] = useState(null)
  const [showAdvancedFilter, setShowAdvancedFilter] = useState(false)

  // Update filtered whenever phase molecules or advanced filter change
  useEffect(() => {
    if (advancedFilterFn) {
      setFilteredMolecules(annotatedMolecules.filter(advancedFilterFn))
    } else {
      setFilteredMolecules(annotatedMolecules)
    }
  }, [annotatedMolecules, advancedFilterFn])

  // Auto-select first molecule when molecules load and none is selected
  useEffect(() => {
    if (flatMolecules.length > 0 && !selectedMolecule) {
      setSelectedMolecule(flatMolecules[0])
    }
  }, [flatMolecules.length])

  // Active run (first created/running run) + queued count
  const activeRun = useMemo(
    () => currentPhaseRuns.find(r => r.status === 'created' || r.status === 'running') || null,
    [currentPhaseRuns]
  )
  const queuedRunCount = useMemo(
    () => Math.max(0, currentPhaseRuns.filter(r => r.status === 'created' || r.status === 'running').length - 1),
    [currentPhaseRuns]
  )

  // Track run completions/failures for toast notifications
  const prevRunsRef = useRef(currentPhaseRuns)
  useEffect(() => {
    const prev = prevRunsRef.current
    prevRunsRef.current = currentPhaseRuns
    if (!prev.length) return

    for (const run of currentPhaseRuns) {
      const prevRun = prev.find(r => r.id === run.id)
      if (!prevRun) continue
      if (prevRun.status !== 'completed' && run.status === 'completed') {
        addToast('Run completed successfully', 'success')
      }
      if (prevRun.status !== 'failed' && run.status === 'failed') {
        addToast(run.error_message || 'Run failed', 'error')
      }
    }
  }, [currentPhaseRuns, addToast])

  // Column visibility — show all available columns by default (CDC §4.1)
  const [visibleKeys, setVisibleKeys] = useState(() => COLUMN_PRESETS.hit_discovery)
  const availableKeysStr = useMemo(
    () => availableColumns.map(c => c.key).join(','),
    [availableColumns]
  )
  useEffect(() => {
    if (currentPhase) {
      // CDC: "all columns with results are shown by default"
      // Use available columns if we have data, otherwise fall back to phase preset
      if (availableColumns.length > 1) {
        setVisibleKeys(availableColumns.map(c => c.key))
      } else {
        const preset = COLUMN_PRESETS[currentPhase.type] || COLUMN_PRESETS.hit_discovery
        setVisibleKeys(preset)
      }
      setSelectedMolecule(null)
    }
  }, [currentPhase?.id, availableKeysStr])

  // Freeze status
  const freezeOverride = currentPhase ? getPhaseStatus(currentPhase.id) : null
  const effectiveStatus = freezeOverride || currentPhase?.status || 'active'
  const isFrozen = effectiveStatus === 'frozen'

  // Phase type meta
  const phaseTypeMeta = currentPhase ? (PHASE_TYPES[currentPhase.type] || {}) : {}

  // Visible columns
  const visibleColumns = useMemo(
    () => ALL_COLUMNS.filter(c => visibleKeys.includes(c.key)),
    [visibleKeys]
  )

  // Server mode threshold: use server-side filtering when total > 5000
  const useServerFilters = moleculeTotal > 5000

  // Bookmarked count — use server stats if in server mode, else compute locally
  const bookmarkedCount = useMemo(
    () => useServerFilters ? (moleculeStats.bookmarked || 0) : flatMolecules.filter(m => m.bookmarked).length,
    [useServerFilters, moleculeStats.bookmarked, flatMolecules]
  )

  // Live stats — use server stats for totals when in server mode
  const liveStats = useMemo(() => {
    if (!currentPhase) return { total_molecules: 0, bookmarked: 0, runs_completed: 0, runs_running: 0 }
    const completedRuns = currentPhaseRuns.filter(r => r.status === 'completed').length
    const runningRuns = currentPhaseRuns.filter(r => r.status === 'running' || r.status === 'created').length
    return {
      total_molecules: useServerFilters ? moleculeTotal : flatMolecules.length,
      bookmarked: bookmarkedCount,
      runs_completed: completedRuns,
      runs_running: runningRuns,
    }
  }, [currentPhase, currentPhaseRuns, flatMolecules, bookmarkedCount, useServerFilters, moleculeTotal])

  // Number of active filters
  const [activeFilterCount, setActiveFilterCount] = useState(0)

  const handleFilteredChange = useCallback((result) => {
    setFilteredMolecules(result)
  }, [])

  // Row click — opens detail panel
  const handleRowClick = useCallback((mol) => {
    setSelectedMolecule(prev => prev?.id === mol.id ? null : mol)
  }, [])

  const handleCloseDetail = useCallback(() => setSelectedMolecule(null), [])

  const handleCellPopup = useCallback((type, mol) => {
    // Find the original molecule with full properties (not flattened)
    const original = phaseMolecules.find(m => m.id === mol.id)
    setPopupState({ type, molecule: { ...mol, properties: original?.properties || {} } })
  }, [phaseMolecules])

  const [runSubmitting, setRunSubmitting] = useState(false)

  const handleNewRun = useCallback(() => {
    if (activeRun) {
      addToast('A run is already in progress. Wait for it to finish before launching a new one.', 'warning')
      return
    }
    setShowRunCreator(true)
  }, [activeRun, addToast])

  const handleRunSubmit = useCallback(async (runConfig) => {
    if (!phaseId) return
    try {
      setRunSubmitting(true)

      if (runConfig.type === 'import' && runConfig.config?.sourceMode === 'external' && runConfig.config?._file) {
        // File upload import
        await importFile(phaseId, runConfig.config._file)
        addToast('File uploaded and import started', 'success')
      } else if (runConfig.type === 'import' && runConfig.config?.sourceMode === 'internal') {
        // Internal import — bookmarked molecules from another phase
        await createRun(phaseId, {
          type: 'import',
          config: {
            source: 'phase_selection',
            selection: 'bookmarked',
            source_phase_id: runConfig.config.source_phase_id,
          },
        })
        addToast('Importing bookmarked molecules...', 'success')
      } else if (runConfig.type === 'import' && runConfig.config?.sourceMode === 'database') {
        // Database import — get uniprot_id from project target
        const uniprotId = currentProject?.target_input_value || currentProject?.target_preview?.uniprot?.id || ''
        await importDatabase(phaseId, {
          databases: runConfig.config.databases,
          uniprot_id: uniprotId,
          max_per_source: runConfig.config.maxPerSource || 50,
          filters: runConfig.config.filters || {},
        })
        addToast(`Fetching compounds from ${runConfig.config.databases.join(', ')}...`, 'success')
      } else {
        // Standard run (calculation, generation, SMILES import)
        await createRun(phaseId, runConfig)
        addToast('Run created and dispatched', 'success')
      }

      setShowRunCreator(false)
    } catch (err) {
      addToast(err.userMessage || err.message || 'Failed to create run', 'error')
    } finally {
      setRunSubmitting(false)
    }
  }, [phaseId, createRun, importFile, importDatabase, currentProject, addToast])

  const handleCancelRun = useCallback(async (runId) => {
    try {
      await cancelRun(runId)
      addToast('Run cancelled', 'info')
    } catch (err) {
      addToast(err.userMessage || 'Failed to cancel run', 'error')
    }
  }, [cancelRun, addToast])

  const handleArchiveRun = useCallback(async (runId) => {
    try {
      await archiveRun(runId)
      addToast('Run archived', 'info')
    } catch (err) {
      addToast(err.userMessage || 'Failed to archive run', 'error')
    }
  }, [archiveRun, addToast])

  const handleExport = useCallback(() => {
    setShowExportModal(true)
  }, [])

  const handleFreezeToggle = useCallback(() => setShowFreezeDialog(true), [])
  const handleFreezeConfirm = useCallback(() => {
    if (!currentPhase) return
    if (isFrozen) {
      unfreezePhase(currentPhase.id)
      addToast('Phase unfrozen — you can now modify molecules and run analyses', 'info')
    } else {
      const bCount = flatMolecules.filter(m => m.bookmarked).length
      freezePhase(currentPhase.id)
      addToast(`Phase frozen — ${bCount} bookmarked molecule${bCount !== 1 ? 's' : ''} locked for next phase`, 'success')
    }
  }, [isFrozen, currentPhase, freezePhase, unfreezePhase, addToast, flatMolecules])

  const handleSelectAll = useCallback(() => {
    if (flatMolecules.every(m => selectedMoleculeIds.has(m.id))) selectNone()
    else selectAll()
  }, [flatMolecules, selectedMoleculeIds, selectAll, selectNone])

  const handleSelectFiltered = useCallback(() => {
    filteredMolecules.forEach(m => {
      if (!selectedMoleculeIds.has(m.id)) toggleSelection(m.id)
    })
  }, [filteredMolecules, selectedMoleculeIds, toggleSelection])

  const hasDownstreamRuns = useMemo(() => {
    if (!currentCampaign || !currentPhase) return false
    const phases = currentCampaign.phases || []
    const myIdx = phases.findIndex(p => p.id === currentPhase.id)
    if (myIdx < 0) return false
    return phases.slice(myIdx + 1).some(p => (p.runs || []).length > 0)
  }, [currentCampaign, currentPhase])

  // Send bookmarks to next phase
  const navigate = useNavigate()
  const handleSendToNextPhase = useCallback(async () => {
    if (!currentCampaign || !currentPhase) return
    const phases = currentCampaign.phases || []
    const myIdx = phases.findIndex(p => p.id === currentPhase.id)
    if (myIdx < 0) return

    const nextPhase = phases[myIdx + 1]
    if (!nextPhase) {
      addToast('No next phase in this campaign. Create a new phase from the project page first.', 'warning')
      return
    }

    try {
      await createRun(nextPhase.id, {
        type: 'import',
        config: {
          source: 'phase_selection',
          selection: 'bookmarked',
          source_phase_id: currentPhase.id,
        },
      })
      addToast(`Sending ${bookmarkedCount} bookmarked molecules to ${nextPhase.label || 'next phase'}...`, 'success')
      navigate(`/project/${projectId}/phase/${nextPhase.id}`)
    } catch (err) {
      addToast(err.userMessage || 'Failed to send bookmarks to next phase', 'error')
    }
  }, [currentCampaign, currentPhase, createRun, bookmarkedCount, addToast, navigate, projectId])

  const showDetail = !!selectedMolecule

  // --- Loading / not found ---
  if (!currentPhase) {
    return (
      <div className="flex flex-col items-center justify-center h-64 gap-3">
        <div className="w-12 h-12 rounded-full bg-gray-100 flex items-center justify-center">
          <svg className="w-6 h-6 text-gray-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
              d="M9.172 16.172a4 4 0 015.656 0M9 10h.01M15 10h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
          </svg>
        </div>
        <p className="text-gray-400 text-sm">
          {currentProject ? 'Phase not found.' : 'Loading workspace...'}
        </p>
      </div>
    )
  }

  // --- Empty phase (no molecules) ---
  if (flatMolecules.length === 0) {
    return (
      <div className="space-y-4 pb-8">
        {/* Breadcrumb */}
        <Breadcrumb projectId={projectId} projectName={currentProject?.name} phase={currentPhase} phaseTypeMeta={phaseTypeMeta} />

        {/* Phase header */}
        <PhaseHeader
          phase={currentPhase}
          phaseTypeMeta={phaseTypeMeta}
          isFrozen={isFrozen}
          stats={liveStats}
          onFreezeToggle={handleFreezeToggle}
          onNewRun={handleNewRun}
          hasActiveRun={!!activeRun}
          onOpenAgent={() => setShowAgentPanel(true)}
        />

        {/* Active run progress */}
        <RunProgress run={activeRun} onCancel={handleCancelRun} queuedCount={queuedRunCount} />

        {/* Empty state */}
        <div className="card p-12 text-center">
          <div className="flex flex-col items-center gap-4 max-w-sm mx-auto">
            <div className="w-16 h-16 rounded-full bg-gray-100 flex items-center justify-center">
              <svg className="w-8 h-8 text-gray-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1}
                  d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
              </svg>
            </div>
            <div>
              <h3 className="text-base font-bold text-gray-700 mb-1">No molecules yet</h3>
              <p className="text-sm text-gray-500 leading-relaxed">
                Run an Import to add molecules to this phase. You can import from ChEMBL, ZINC, a file, or promote hits from a previous phase.
              </p>
            </div>
            {!isFrozen && (
              <button
                onClick={handleNewRun}
                className="flex items-center gap-2 px-5 py-2.5 bg-bx-surface text-white rounded-xl text-sm font-semibold hover:bg-bx-elevated transition-colors shadow-sm"
              >
                <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 4v16m8-8H4" />
                </svg>
                Import Molecules
              </button>
            )}
          </div>
        </div>

        <RunHistory runs={currentPhaseRuns} onCancel={handleCancelRun} onArchive={handleArchiveRun} />

        <RunCreator
          phaseId={phaseId}
          phaseType={currentPhase.type}
          isOpen={showRunCreator}
          onClose={() => setShowRunCreator(false)}
          onSubmit={handleRunSubmit}
          selectedMoleculeIds={selectedMoleculeIds}
          submitting={runSubmitting}
        />
        <FreezeDialog
          isOpen={showFreezeDialog}
          onClose={() => setShowFreezeDialog(false)}
          onConfirm={handleFreezeConfirm}
          action={isFrozen ? 'unfreeze' : 'freeze'}
          phaseName={currentPhase.label}
          bookmarkedMols={flatMolecules.filter(m => m.bookmarked)}
          hasDownstreamRuns={hasDownstreamRuns}
        />
      </div>
    )
  }

  const phase = currentPhase

  return (
    <div className="space-y-4 pb-8">
      {/* Breadcrumb */}
      <Breadcrumb projectId={projectId} projectName={currentProject?.name} phase={phase} phaseTypeMeta={phaseTypeMeta} />

      {/* Phase header card */}
      <PhaseHeader
        phase={phase}
        phaseTypeMeta={phaseTypeMeta}
        isFrozen={isFrozen}
        stats={liveStats}
        onFreezeToggle={handleFreezeToggle}
        onNewRun={handleNewRun}
        hasActiveRun={!!activeRun}
        onOpenAgent={() => setShowAgentPanel(true)}
      />

      {/* Active run progress */}
      <RunProgress run={activeRun} onCancel={handleCancelRun} queuedCount={queuedRunCount} />

      {/* Selection Toolbar */}
      <SelectionToolbar
        selectedCount={selectedMoleculeIds.size}
        totalCount={flatMolecules.length}
        bookmarkedCount={bookmarkedCount}
        filteredCount={filteredMolecules.length}
        activeFilterCount={activeFilterCount}
        onSelectAll={handleSelectAll}
        onSelectNone={selectNone}
        onSelectBookmarked={selectBookmarked}
        onSelectFiltered={handleSelectFiltered}
        onExport={handleExport}
        onBookmarkSelected={bookmarkSelected}
        onSendToNextPhase={!isFrozen && bookmarkedCount > 0 ? handleSendToNextPhase : undefined}
        isFrozen={isFrozen}
      />

      {/* Filter bar */}
      <div className="space-y-2">
        <div className="flex items-center gap-2">
          <div className="flex-1">
            <FilterBarWithCount
              molecules={flatMolecules}
              columns={availableColumns}
              onFilteredChange={handleFilteredChange}
              onFilterCountChange={setActiveFilterCount}
              serverMode={useServerFilters}
              onServerFilterChange={updateServerFilters}
              totalFromServer={useServerFilters ? moleculeTotal : undefined}
            />
          </div>
          <button
            onClick={() => setShowAdvancedFilter(v => !v)}
            className={`flex items-center gap-1.5 px-3 py-1.5 rounded-lg text-xs font-medium border transition-colors flex-shrink-0 ${
              showAdvancedFilter
                ? 'bg-bx-surface text-white border-bx-surface'
                : 'bg-white text-gray-500 border-gray-200 hover:border-gray-300 hover:text-gray-700'
            }`}
          >
            <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                d="M3 4a1 1 0 011-1h16a1 1 0 011 1v2.586a1 1 0 01-.293.707l-6.414 6.414a1 1 0 00-.293.707V17l-4 4v-6.586a1 1 0 00-.293-.707L3.293 7.293A1 1 0 013 6.586V4z" />
            </svg>
            Advanced
            {advancedFilterFn && (
              <span className="w-1.5 h-1.5 rounded-full bg-bx-cyan" />
            )}
          </button>
        </div>

        <AdvancedFilterBuilder
          columns={availableColumns}
          onFilterChange={setAdvancedFilterFn}
          isOpen={showAdvancedFilter}
          onToggle={() => setShowAdvancedFilter(v => !v)}
        />
      </div>

      {/* Table + Detail panel split layout */}
      <div className="flex gap-4 items-start">
        {/* Table (60% when drawer open, 100% when closed) */}
        <div className={`min-w-0 space-y-2 transition-all duration-200 ${showDetail ? 'flex-1' : 'w-full'}`}>
          {/* Column controls row */}
          <div className="flex items-center justify-between gap-3">
            <div className="flex items-center gap-2">
              <ColumnSelector
                allColumns={ALL_COLUMNS}
                visibleKeys={visibleKeys}
                onChange={setVisibleKeys}
              />
              <span className="text-sm text-gray-400 tabular-nums flex items-center gap-1.5">
                {moleculesLoading && <BindXLogo variant="loading" size={20} />}
                {filteredMolecules.length !== moleculeTotal
                  ? `${filteredMolecules.length} of ${moleculeTotal.toLocaleString()} molecules`
                  : `${moleculeTotal.toLocaleString()} molecules`}
                {moleculeHasMore && <span className="text-[10px] text-gray-300">(scroll for more)</span>}
              </span>
            </div>
          </div>

          <MoleculeTable
            molecules={filteredMolecules}
            columns={visibleColumns}
            selectedIds={selectedMoleculeIds}
            onToggleSelect={toggleSelection}
            onSelectAll={handleSelectAll}
            onRowClick={handleRowClick}
            onToggleBookmark={isFrozen ? undefined : toggleBookmark}
            onCellPopup={handleCellPopup}
            onAnnotation={isFrozen ? undefined : handleAnnotation}
            activeRowId={selectedMolecule?.id || null}
            onLoadMore={moleculeHasMore ? loadMoreMolecules : undefined}
            hasMore={moleculeHasMore}
            totalCount={moleculeTotal}
          />
        </div>

        {/* Detail panel (40% fixed width when open) */}
        {showDetail && selectedMolecule && (
          <div className="flex-shrink-0 w-[40%] min-w-[300px] max-w-[520px]">
            <MoleculeDetailPanel
              molecule={flatMolecules.find(m => m.id === selectedMolecule.id) || selectedMolecule}
              molecules={filteredMolecules}
              onClose={handleCloseDetail}
              onToggleBookmark={isFrozen ? undefined : toggleBookmark}
              onRowClick={handleRowClick}
              isFrozen={isFrozen}
              project={currentProject}
            />
          </div>
        )}
      </div>

      {/* Run History timeline */}
      <RunHistory runs={currentPhaseRuns} onCancel={handleCancelRun} onArchive={handleArchiveRun} />

      {/* Pareto Analysis (collapsible) */}
      <div className="card overflow-hidden">
        <button
          onClick={() => setShowPareto(v => !v)}
          className="w-full flex items-center justify-between px-5 py-3.5 hover:bg-gray-50 transition-colors text-left"
        >
          <div className="flex items-center gap-2.5">
            <svg className="w-4 h-4 text-bx-light-text" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                d="M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 01-2 2h-2a2 2 0 01-2-2z" />
            </svg>
            <span className="font-semibold text-gray-700 text-sm">Pareto Analysis</span>
            <InfoTip text={TIPS.pareto} size="xs" />
            <span className="text-sm text-gray-400">2D objective scatter plot</span>
          </div>
          <svg
            className={`w-4 h-4 text-gray-400 transition-transform duration-150 ${showPareto ? 'rotate-180' : ''}`}
            fill="none" stroke="currentColor" viewBox="0 0 24 24"
          >
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
          </svg>
        </button>
        {showPareto && (
          <div className="px-5 pb-5 border-t border-gray-100 pt-4">
            <ParetoFront
              molecules={filteredMolecules}
              onSelect={mol => mol && handleRowClick(mol)}
              onToggleBookmark={isFrozen ? undefined : toggleBookmark}
            />
          </div>
        )}
      </div>

      {/* Run Logs (collapsible terminal) */}
      <div className="card overflow-hidden">
        <button
          onClick={() => setShowRunLogs(v => !v)}
          className="w-full flex items-center justify-between px-5 py-3 hover:bg-gray-50 transition-colors text-left"
        >
          <div className="flex items-center gap-2.5">
            <svg className="w-4 h-4 text-gray-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                d="M8 9l3 3-3 3m5 0h3M5 20h14a2 2 0 002-2V6a2 2 0 00-2-2H5a2 2 0 00-2 2v12a2 2 0 002 2z" />
            </svg>
            <span className="font-semibold text-gray-700 text-sm">Run Logs</span>
            {currentPhaseRuns.some(r => r.logs?.length > 0) && (
              <span className="text-[10px] text-gray-400 tabular-nums">
                {currentPhaseRuns.reduce((s, r) => s + (r.logs?.length || 0), 0)} entries
              </span>
            )}
          </div>
          <svg
            className={`w-4 h-4 text-gray-400 transition-transform duration-150 ${showRunLogs ? 'rotate-180' : ''}`}
            fill="none" stroke="currentColor" viewBox="0 0 24 24"
          >
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
          </svg>
        </button>
        {showRunLogs && (
          <div className="border-t border-gray-100">
            <div className="bg-gray-900 rounded-b-lg max-h-64 overflow-y-auto px-4 py-3 font-mono text-xs space-y-0.5 scrollbar-thin">
              {(() => {
                const allLogs = currentPhaseRuns
                  .filter(r => r.logs?.length > 0)
                  .flatMap(r => (r.logs || []).map(log => ({
                    ...log,
                    runType: r.type,
                    runId: r.id,
                    ts: log.timestamp ? new Date(log.timestamp).getTime() : 0,
                  })))
                  .sort((a, b) => a.ts - b.ts)

                if (allLogs.length === 0) {
                  return (
                    <div className="text-gray-500 text-center py-6">
                      No run logs available yet. Logs will appear here when runs produce output.
                    </div>
                  )
                }

                return allLogs.map((log, i) => {
                  const levelColor = log.level === 'error' ? 'text-red-400'
                    : log.level === 'warn' || log.level === 'warning' ? 'text-amber-400'
                    : log.level === 'success' ? 'text-emerald-400'
                    : 'text-gray-400'
                  return (
                    <div key={i} className={`flex gap-2 ${levelColor}`}>
                      {log.timestamp && (
                        <span className="text-gray-600 flex-shrink-0 w-[70px]">
                          {new Date(log.timestamp).toLocaleTimeString()}
                        </span>
                      )}
                      <span className="text-gray-500 flex-shrink-0 w-[65px] truncate">[{log.runType || 'run'}]</span>
                      <span className="break-all">{log.message || String(log)}</span>
                    </div>
                  )
                })
              })()}
            </div>
          </div>
        )}
      </div>

      {/* Modals */}
      <RunCreator
        phaseId={phaseId}
        phaseType={phase.type}
        isOpen={showRunCreator}
        onClose={() => setShowRunCreator(false)}
        onSubmit={handleRunSubmit}
        selectedMoleculeIds={selectedMoleculeIds}
        submitting={runSubmitting}
      />

      <FreezeDialog
        isOpen={showFreezeDialog}
        onClose={() => setShowFreezeDialog(false)}
        onConfirm={handleFreezeConfirm}
        action={isFrozen ? 'unfreeze' : 'freeze'}
        phaseName={phase.label}
        bookmarkedMols={flatMolecules.filter(m => m.bookmarked)}
        hasDownstreamRuns={hasDownstreamRuns}
      />

      {/* AI Agent panel */}
      <CampaignAgentPanel
        isOpen={showAgentPanel}
        onClose={() => setShowAgentPanel(false)}
        project={currentProject}
        campaign={currentCampaign}
        phase={currentPhase}
        molecules={flatMolecules}
      />

      {/* Export modal */}
      <ExportModal
        isOpen={showExportModal}
        onClose={() => setShowExportModal(false)}
        molecules={filteredMolecules}
        selectedIds={selectedMoleculeIds}
        columns={visibleColumns}
        onToast={addToast}
      />

      {/* Detail popups (Safety / Synthesis / Confidence) */}
      {popupState && (
        <DetailPopupModal
          type={popupState.type}
          molecule={popupState.molecule}
          onClose={() => setPopupState(null)}
        />
      )}
    </div>
  )
}
