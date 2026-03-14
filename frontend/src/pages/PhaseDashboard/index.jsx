import React, { useEffect, useState, useMemo, useCallback, useRef } from 'react'
import { useParams, Link, useNavigate } from 'react-router-dom'

import { useWorkspace } from '../../contexts/WorkspaceContext.jsx'
import { useToast } from '../../contexts/ToastContext.jsx'
import { v9UpdateAnnotations } from '../../api.js'
import { ALL_COLUMNS, COLUMN_PRESETS, PHASE_TYPES, flattenMoleculeProperties, detectAvailableColumns } from '../../lib/columns.js'

import MoleculeTable from '../../components/MoleculeTable.jsx'
import TableToolbar from '../../components/TableToolbar.jsx'
import RunHistory from '../../components/RunHistory.jsx'
import RunCreator from '../../components/RunCreator.jsx'
import RunProgress from '../../components/RunProgress.jsx'
import FreezeDialog from '../../components/FreezeDialog.jsx'
import ParetoFront from '../../components/ParetoFront.jsx'
import AnalyticsPanel from '../../components/AnalyticsPanel.jsx'
import ExportModal from '../../components/ExportModal.jsx'
import CampaignAgentPanel from '../../components/CampaignAgentPanel.jsx'
import InfoTip, { TIPS } from '../../components/InfoTip.jsx'

import AdvancedFilterBuilder from '../../components/AdvancedFilterBuilder.jsx'
import MoleculeDetailPanel, { DetailPopupModal } from './MoleculeDetailPanel.jsx'
import PhaseHeader, { Breadcrumb } from './PhaseHeader.jsx'

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
    deletePhase,
    getPhaseStatus,
    createRun,
    cancelRun,
    archiveRun,
    importFile,
    importDatabase,
    moleculeTotal,
    moleculePage,
    moleculePageSize,
    moleculeStats,
    goToPage,
    setPageSize,
    serverSort,
    updateServerSort,
    updateServerFilters,
    createPhase,
  } = useWorkspace()

  const { addToast } = useToast()
  const navigate = useNavigate()

  // Sync URL params to workspace state (only when URL params change, not on function identity change)
  // eslint-disable-next-line react-hooks/exhaustive-deps
  useEffect(() => { if (projectId) selectProject(projectId) }, [projectId])
  // eslint-disable-next-line react-hooks/exhaustive-deps
  useEffect(() => { if (phaseId) selectPhase(phaseId) }, [phaseId])

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

  // --- User annotations (tags, invalidation, comments) — optimistic local + API persist ---
  const [annotations, setAnnotations] = useState({})
  const debounceTimers = useRef({})

  const handleAnnotation = useCallback((molId, key, value) => {
    // Optimistic local update
    setAnnotations(prev => ({
      ...prev,
      [molId]: { ...(prev[molId] || {}), [key]: value },
    }))
    // Debounced API persist (300ms for text, immediate for tags/boolean)
    const delay = key === 'user_comment' ? 600 : 100
    const timerKey = `${molId}:${key}`
    if (debounceTimers.current[timerKey]) clearTimeout(debounceTimers.current[timerKey])
    debounceTimers.current[timerKey] = setTimeout(() => {
      v9UpdateAnnotations(molId, { [key]: value }).catch(err => {
        console.error('Failed to save annotation:', err)
      })
      delete debounceTimers.current[timerKey]
    }, delay)
  }, [])

  // Merge annotations into flat molecules (API fields + local overrides)
  const annotatedMolecules = useMemo(
    () => flatMolecules.map(m => ({ ...m, ...(annotations[m.id] || {}) })),
    [flatMolecules, annotations]
  )

  // --- Local UI state ---
  const [selectedMolecule, setSelectedMolecule] = useState(null)
  const [highlightedIds, setHighlightedIds] = useState(new Set())
  const [showRunCreator, setShowRunCreator] = useState(false)
  const [showPareto, setShowPareto] = useState(true)
  const [showFreezeDialog, setShowFreezeDialog] = useState(false)
  const [filterBarResult, setFilterBarResult] = useState(null)
  const [filteredMolecules, setFilteredMolecules] = useState(annotatedMolecules)
  const [popupState, setPopupState] = useState(null) // { type: 'safety'|'retrosynthesis'|'confidence', molecule }
  const [showExportModal, setShowExportModal] = useState(false)
  const [showAgentPanel, setShowAgentPanel] = useState(false)
  const [showRunLogs, setShowRunLogs] = useState(false)
  const [showAnalytics, setShowAnalytics] = useState(true)
  const [advancedFilterFn, setAdvancedFilterFn] = useState(null)
  const [showAdvancedFilter, setShowAdvancedFilter] = useState(false)

  // Reset stale filterBarResult when new molecules arrive from API
  const flatMolsRef = useRef(flatMolecules)
  if (flatMolsRef.current !== flatMolecules) {
    flatMolsRef.current = flatMolecules
    // filterBarResult holds old page data — null it so we use fresh flatMolecules
    if (filterBarResult) setFilterBarResult(null)
  }

  // Combine FilterBar result + advanced filter + annotations overlay
  useEffect(() => {
    const base = filterBarResult || flatMolecules
    // Re-apply annotations on top (filterBarResult has stale snapshots without annotations)
    const withAnnotations = base.map(m => ({ ...m, ...(annotations[m.id] || {}) }))
    if (advancedFilterFn) {
      setFilteredMolecules(withAnnotations.filter(advancedFilterFn))
    } else {
      setFilteredMolecules(withAnnotations)
    }
  }, [filterBarResult, flatMolecules, annotations, advancedFilterFn])

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

  // Bookmarked count — always from server stats (pagination = only 1 page loaded)
  const bookmarkedCount = moleculeStats.bookmarked || 0

  // Live stats — always use server totals (client only has current page)
  const liveStats = useMemo(() => {
    if (!currentPhase) return { total_molecules: 0, bookmarked: 0, runs_completed: 0, runs_running: 0 }
    const completedRuns = currentPhaseRuns.filter(r => r.status === 'completed').length
    const runningRuns = currentPhaseRuns.filter(r => r.status === 'running' || r.status === 'created').length
    return {
      total_molecules: moleculeTotal,
      bookmarked: bookmarkedCount,
      runs_completed: completedRuns,
      runs_running: runningRuns,
    }
  }, [currentPhase, currentPhaseRuns, bookmarkedCount, moleculeTotal])

  const handleFilteredChange = useCallback((result) => {
    setFilterBarResult(result)
  }, [])

  // Row click — opens detail panel. Ctrl/Cmd+click accumulates highlights for 3D overlay.
  const handleRowClick = useCallback((mol, e) => {
    if (e && (e.ctrlKey || e.metaKey)) {
      // Ctrl+click: toggle this molecule in highlighted set
      setHighlightedIds(prev => {
        const next = new Set(prev)
        if (next.has(mol.id)) next.delete(mol.id)
        else next.add(mol.id)
        return next
      })
      // Keep detail panel open on the clicked molecule
      setSelectedMolecule(mol)
    } else {
      // Normal click: single select, clear highlights
      setSelectedMolecule(prev => prev?.id === mol.id ? null : mol)
      setHighlightedIds(new Set())
    }
  }, [])

  const handleCloseDetail = useCallback(() => {
    setSelectedMolecule(null)
    setHighlightedIds(new Set())
  }, [])

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
        // Resolve source_phase_id by phase type order (hit_discovery → hit_to_lead → lead_optimization)
        let sourcePhaseId = runConfig.config.source_phase_id
        if (!sourcePhaseId) {
          const PHASE_ORDER = ['hit_discovery', 'hit_to_lead', 'lead_optimization']
          const phases = currentCampaign?.phases || []
          const myType = currentPhase?.type || phases.find(p => String(p.id) === String(phaseId))?.type
          const myOrderIdx = PHASE_ORDER.indexOf(myType)
          if (myOrderIdx > 0) {
            // Find the phase of the previous type
            const prevType = PHASE_ORDER[myOrderIdx - 1]
            const prevPhase = phases.find(p => p.type === prevType)
            if (prevPhase) sourcePhaseId = prevPhase.id
          }
          // Fallback: try index-based lookup
          if (!sourcePhaseId && phases.length > 1) {
            const myIdx = phases.findIndex(p => String(p.id) === String(phaseId))
            if (myIdx > 0) sourcePhaseId = phases[myIdx - 1].id
          }
        }
        if (!sourcePhaseId) {
          addToast('No previous phase found to import from. This is the first phase.', 'error')
          return
        }
        await createRun(phaseId, {
          type: 'import',
          config: {
            source: 'phase_selection',
            selection: 'bookmarked',
            source_phase_id: sourcePhaseId,
          },
        })
        addToast('Importing bookmarked molecules...', 'success')
      } else if (runConfig.type === 'import' && runConfig.config?.sourceMode === 'database') {
        // Database import — single source with specific config
        const uniprotId = currentProject?.target_input_value || currentProject?.target_preview?.uniprot?.id || ''
        const db = runConfig.config.database
        const payload = {
          database: db,
          uniprot_id: uniprotId,
          max_compounds: runConfig.config.maxPerSource || 100,
          filters: runConfig.config.filters || {},
        }
        // Source-specific config
        if (db === 'zinc') payload.zinc_subset = runConfig.config.zinc_subset || 'in-stock'
        if (db === 'enamine') payload.scaffold_complexity = runConfig.config.filters?.scaffold_complexity || 'mixed'

        await importDatabase(phaseId, payload)
        const sourceLabels = { zinc: 'ZINC20', chembl: 'ChEMBL', pubchem: 'PubChem', enamine: 'Enamine REAL', fragments: 'Fragment Library' }
        addToast(`Fetching compounds from ${sourceLabels[db] || db}...`, 'success')
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
  }, [phaseId, createRun, importFile, importDatabase, currentProject, currentCampaign, currentPhase, addToast])

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

  const handleDeletePhase = useCallback(async () => {
    if (!currentPhase) return
    const phaseName = PHASE_TYPES[currentPhase.type]?.label || currentPhase.type
    const molCount = moleculeTotal || flatMolecules.length
    const runCount = (currentPhaseRuns || []).length
    const msg = `Delete phase "${phaseName}"? This will permanently remove ${molCount} molecule${molCount !== 1 ? 's' : ''} and ${runCount} run${runCount !== 1 ? 's' : ''}. This cannot be undone.`
    if (!window.confirm(msg)) return
    try {
      await deletePhase(currentPhase.id)
      addToast('Phase deleted', 'success')
      navigate(`/project/${projectId}`)
    } catch (err) {
      addToast(err.userMessage || 'Failed to delete phase', 'error')
    }
  }, [currentPhase, flatMolecules.length, currentPhaseRuns, deletePhase, addToast, navigate, projectId])

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
  const handleSendToNextPhase = useCallback(async () => {
    if (!currentCampaign || !currentPhase) return
    const phases = currentCampaign.phases || []
    const typeOrder = ['hit_discovery', 'hit_to_lead', 'lead_optimization']
    const currentTypeIdx = typeOrder.indexOf(currentPhase.type)

    if (currentTypeIdx < 0 || currentTypeIdx >= typeOrder.length - 1) {
      addToast('This is the final phase (Lead Optimization). No further phase can be created.', 'warning')
      return
    }

    const nextType = typeOrder[currentTypeIdx + 1]

    // Find existing next phase by type (more reliable than index)
    let nextPhase = phases.find(p => p.type === nextType)

    // If next phase already exists, just navigate to it
    if (nextPhase) {
      navigate(`/project/${projectId}/phase/${nextPhase.id}`)
      return
    }

    // Auto-create next phase if it doesn't exist
    const nextLabel = PHASE_TYPES[nextType]?.label || nextType
    try {
      nextPhase = await createPhase(currentCampaign.id, { type: nextType })
      addToast(`Created ${nextLabel} phase`, 'success')
    } catch (err) {
      // If 409 (already exists), find it and navigate
      if (err?.status === 409) {
        const existing = phases.find(p => p.type === nextType)
        if (existing) {
          navigate(`/project/${projectId}/phase/${existing.id}`)
          return
        }
      }
      addToast(err.userMessage || 'Failed to create next phase', 'error')
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
      addToast(`Sending ${bookmarkedCount} bookmarked molecules to ${PHASE_TYPES[nextType]?.label || 'next phase'}...`, 'success')
      navigate(`/project/${projectId}/phase/${nextPhase.id}`)
    } catch (err) {
      addToast(err.userMessage || 'Failed to send bookmarks to next phase', 'error')
    }
  }, [currentCampaign, currentPhase, createRun, createPhase, bookmarkedCount, addToast, navigate, projectId])

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
          onDeletePhase={handleDeletePhase}
          onSendToNextPhase={handleSendToNextPhase}
          bookmarkedCount={bookmarkedCount}
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
          project={currentProject}
        />
        <FreezeDialog
          isOpen={showFreezeDialog}
          onClose={() => setShowFreezeDialog(false)}
          onConfirm={handleFreezeConfirm}
          action={isFrozen ? 'unfreeze' : 'freeze'}
          phaseName={PHASE_TYPES[currentPhase.type]?.label || currentPhase.type}
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
        onDeletePhase={handleDeletePhase}
        onSendToNextPhase={handleSendToNextPhase}
        bookmarkedCount={bookmarkedCount}
      />

      {/* Active run progress */}
      <RunProgress run={activeRun} onCancel={handleCancelRun} queuedCount={queuedRunCount} />

      {/* Unified Toolbar (filters + columns + selection) */}
      <TableToolbar
        molecules={annotatedMolecules}
        columns={availableColumns}
        onFilteredChange={handleFilteredChange}
        onServerFilterChange={updateServerFilters}
        showAdvancedFilter={showAdvancedFilter}
        onToggleAdvancedFilter={() => setShowAdvancedFilter(v => !v)}
        advancedFilterActive={!!advancedFilterFn}
        onClearAdvancedFilter={() => setAdvancedFilterFn(null)}
        visibleKeys={visibleKeys}
        onVisibleKeysChange={setVisibleKeys}
        filteredMoleculeCount={filteredMolecules.length}
        moleculeTotal={moleculeTotal}
        moleculesLoading={moleculesLoading}
        selectedCount={selectedMoleculeIds.size}
        totalCount={flatMolecules.length}
        bookmarkedCount={bookmarkedCount}
        filteredCount={filteredMolecules.length}
        onSelectAll={handleSelectAll}
        onSelectNone={selectNone}
        onSelectBookmarked={selectBookmarked}
        onSelectFiltered={handleSelectFiltered}
        onExport={handleExport}
        onBookmarkSelected={bookmarkSelected}
        onSendToNextPhase={!isFrozen && bookmarkedCount > 0 ? handleSendToNextPhase : undefined}
        isFrozen={isFrozen}
      />

      <AdvancedFilterBuilder
        columns={availableColumns}
        onFilterChange={(fn) => setAdvancedFilterFn(() => fn)}
        isOpen={showAdvancedFilter}
        onToggle={() => setShowAdvancedFilter(v => !v)}
      />

      {/* Table + Detail panel split layout */}
      <div className="flex gap-4 items-start">
        {/* Table (60% when drawer open, 100% when closed) */}
        <div className={`min-w-0 space-y-2 transition-all duration-200 ${showDetail ? 'flex-1' : 'w-full'}`}>

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
            highlightedIds={highlightedIds.size > 0 ? highlightedIds : null}
            totalCount={moleculeTotal}
            runs={currentPhaseRuns}
            currentPage={moleculePage}
            totalPages={Math.max(1, Math.ceil(moleculeTotal / moleculePageSize))}
            pageSize={moleculePageSize}
            onPageChange={goToPage}
            onPageSizeChange={setPageSize}
            onServerSort={updateServerSort}
            serverSortKey={serverSort.sort_by}
            serverSortDir={serverSort.sort_dir}
            globalColumnRanges={moleculeStats.column_ranges}
          />
        </div>

        {/* Detail panel (40% fixed width when open) */}
        {showDetail && selectedMolecule && (
          <div className="flex-shrink-0 w-[40%] min-w-[300px] max-w-[520px]">
            <MoleculeDetailPanel
              molecule={(() => {
                const flat = flatMolecules.find(m => m.id === selectedMolecule.id) || selectedMolecule
                const original = phaseMolecules.find(m => m.id === selectedMolecule.id)
                return original ? { ...flat, properties: original.properties } : flat
              })()}
              molecules={filteredMolecules}
              onClose={handleCloseDetail}
              onToggleBookmark={isFrozen ? undefined : toggleBookmark}
              onRowClick={handleRowClick}
              isFrozen={isFrozen}
              project={currentProject}
              selectedMolecules={highlightedIds.size > 0 ? phaseMolecules.filter(m => highlightedIds.has(m.id)) : null}
              onCellPopup={handleCellPopup}
            />
          </div>
        )}
      </div>

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
              visibleKeys={visibleKeys}
            />
          </div>
        )}
      </div>

      {/* Analytics (collapsible) */}
      <div className="card overflow-hidden">
        <button
          onClick={() => setShowAnalytics(v => !v)}
          className="w-full flex items-center justify-between px-5 py-3.5 hover:bg-gray-50 transition-colors text-left"
        >
          <div className="flex items-center gap-2.5">
            <svg className="w-4 h-4 text-bx-light-text" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                d="M16 8v8m-4-5v5m-4-2v2m-2 4h12a2 2 0 002-2V6a2 2 0 00-2-2H6a2 2 0 00-2 2v12a2 2 0 002 2z" />
            </svg>
            <span className="font-semibold text-gray-700 text-sm">Analytics</span>
            <span className="text-sm text-gray-400">Distributions, correlations, breakdowns</span>
          </div>
          <svg
            className={`w-4 h-4 text-gray-400 transition-transform duration-150 ${showAnalytics ? 'rotate-180' : ''}`}
            fill="none" stroke="currentColor" viewBox="0 0 24 24"
          >
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
          </svg>
        </button>
        {showAnalytics && (
          <div className="px-5 pb-5 border-t border-gray-100 pt-4">
            <AnalyticsPanel
              molecules={filteredMolecules}
              availableColumns={availableColumns}
              visibleKeys={visibleKeys}
            />
          </div>
        )}
      </div>

      {/* Run History timeline */}
      <RunHistory runs={currentPhaseRuns} onCancel={handleCancelRun} onArchive={handleArchiveRun} />

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
        project={currentProject}
      />

      <FreezeDialog
        isOpen={showFreezeDialog}
        onClose={() => setShowFreezeDialog(false)}
        onConfirm={handleFreezeConfirm}
        action={isFrozen ? 'unfreeze' : 'freeze'}
        phaseName={PHASE_TYPES[phase.type]?.label || phase.type}
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
