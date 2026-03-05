/**
 * BindX V9 — Workspace bridge and compat layer.
 *
 * Bridges React-only state (auth, toast) into the Zustand store and
 * provides a backward-compatible useWorkspace() hook so existing
 * consumers don't need to change.
 *
 * For new components or performance-critical paths, import
 * useWorkspaceStore directly from '../stores/workspaceStore'.
 */

import React, { useEffect, useRef } from 'react'
import { useAuth } from './AuthContext.jsx'
import { useToast } from './ToastContext.jsx'
import useWorkspaceStore from '../stores/workspaceStore'

/**
 * WorkspaceBridge — renders nothing, syncs React context into Zustand.
 * Must be placed inside AuthProvider and ToastProvider.
 */
export function WorkspaceBridge() {
  const { isAuthenticated, loading: authLoading } = useAuth()
  const { addToast } = useToast()
  const prevPhaseId = useRef(null)

  // Sync auth state into Zustand store
  useEffect(() => {
    useWorkspaceStore.getState()._syncAuth(isAuthenticated)
  }, [isAuthenticated])

  // Sync toast function
  useEffect(() => {
    useWorkspaceStore.getState()._setToast(addToast)
  }, [addToast])

  // Trigger initial project load once auth settles
  useEffect(() => {
    if (authLoading) return
    useWorkspaceStore.getState().refreshProjects()
  }, [authLoading, isAuthenticated])

  // Watch phase changes → trigger data refresh
  useEffect(() => {
    return useWorkspaceStore.subscribe(
      (state) => {
        const phaseId = state.currentPhaseId
        if (phaseId !== prevPhaseId.current) {
          prevPhaseId.current = phaseId
          state._onPhaseChange(phaseId)
        }
      }
    )
  }, [])

  // Watch server sort/filter changes → re-fetch page 1
  const prevSortRef = useRef(null)
  const prevFiltersRef = useRef(null)
  useEffect(() => {
    return useWorkspaceStore.subscribe(
      (state) => {
        const { serverSort, serverFilters, currentPhaseId, _isAuthenticated } = state
        const sortKey = `${serverSort.sort_by}:${serverSort.sort_dir}`
        const filterKey = `${serverFilters.search}:${serverFilters.bookmarked_only}`

        if (prevSortRef.current === null) {
          prevSortRef.current = sortKey
          prevFiltersRef.current = filterKey
          return
        }

        if (sortKey !== prevSortRef.current || filterKey !== prevFiltersRef.current) {
          prevSortRef.current = sortKey
          prevFiltersRef.current = filterKey
          if (currentPhaseId && _isAuthenticated) {
            state.refreshMolecules(currentPhaseId)
          }
        }
      }
    )
  }, [])

  // Run polling: auto-refresh while runs are active
  useEffect(() => {
    let interval = null

    const check = () => {
      const { phaseRuns, currentPhaseId, _isAuthenticated } = useWorkspaceStore.getState()
      return _isAuthenticated && currentPhaseId && phaseRuns.some(r => r.status === 'created' || r.status === 'running')
    }

    const startPolling = () => {
      if (interval) return
      interval = setInterval(() => {
        if (!check()) {
          clearInterval(interval)
          interval = null
          return
        }
        const s = useWorkspaceStore.getState()
        s.refreshRuns()
        s.refreshStats()
        s.refreshMolecules(s.currentPhaseId)
      }, 5000)
    }

    // Subscribe to any state change, check if polling should start/stop
    const unsub = useWorkspaceStore.subscribe(() => {
      if (check() && !interval) startPolling()
    })

    if (check()) startPolling()

    return () => {
      unsub()
      if (interval) clearInterval(interval)
    }
  }, [])

  return null
}

/**
 * useWorkspace — backward-compatible hook.
 *
 * Returns the same shape as the old WorkspaceContext value.
 * Components re-render on any store change (same as Context did).
 *
 * For performance, use useWorkspaceStore with selectors directly.
 */
export function useWorkspace() {
  const store = useWorkspaceStore()

  // Compute derived values (these were useMemo'd in the old context)
  const currentProject = store.getCurrentProject()
  const currentCampaign = store.getCurrentCampaign()
  const currentPhase = store.getCurrentPhase()
  const phaseMolecules = store.getPhaseMolecules()
  const currentPhaseRuns = store.getCurrentPhaseRuns()

  return {
    // Data
    projects: store.projects,
    currentProject,
    currentCampaign,
    currentPhase,
    phaseMolecules,
    currentPhaseRuns,
    loading: store.loading,
    error: store.error,
    runsLoading: store.runsLoading,
    moleculesLoading: store.moleculesLoading,

    // Paginated molecule metadata
    moleculeTotal: store.moleculeTotal,
    moleculeHasMore: store.moleculeHasMore,
    moleculeStats: store.moleculeStats,
    serverSort: store.serverSort,
    serverFilters: store.serverFilters,

    // CRUD
    createProject: store.createProject,
    updateProject: store.updateProject,
    deleteProject: store.deleteProject,
    createCampaign: store.createCampaign,
    updateCampaign: store.updateCampaign,
    createPhase: store.createPhase,
    refreshProjects: store.refreshProjects,

    // Navigation
    selectProject: store.selectProject,
    selectCampaign: store.selectCampaign,
    selectPhase: store.selectPhase,

    // Selection
    selectedMoleculeIds: store.selectedMoleculeIds,
    toggleSelection: store.toggleSelection,
    selectAll: store.selectAll,
    selectNone: store.selectNone,
    selectBookmarked: store.selectBookmarked,

    // Bookmarks
    toggleBookmark: store.toggleBookmark,
    bookmarkSelected: store.bookmarkSelected,

    // Phase actions
    freezePhase: store.freezePhase,
    unfreezePhase: store.unfreezePhase,
    deletePhase: store.deletePhase,
    getPhaseStatus: store.getPhaseStatus,

    // Runs & Molecules actions
    createRun: store.createRun,
    cancelRun: store.cancelRun,
    archiveRun: store.archiveRun,
    importFile: store.importFile,
    importDatabase: store.importDatabase,
    refreshRuns: store.refreshRuns,
    refreshMolecules: store.refreshMolecules,
    loadMoreMolecules: store.loadMoreMolecules,
    updateServerSort: store.updateServerSort,
    updateServerFilters: store.updateServerFilters,
    refreshStats: store.refreshStats,
  }
}

// Legacy export for backward compat — WorkspaceProvider is now just the bridge
export function WorkspaceProvider({ children }) {
  return (
    <>
      <WorkspaceBridge />
      {children}
    </>
  )
}
