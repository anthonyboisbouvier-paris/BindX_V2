import React, { createContext, useContext, useState, useCallback, useMemo, useEffect, useRef } from 'react'
import { useAuth } from './AuthContext.jsx'
import { useToast } from './ToastContext.jsx'
import {
  v9ListProjects,
  v9CreateProject,
  v9UpdateProject,
  v9DeleteProject,
  v9CreateCampaign,
  v9UpdateCampaign,
  v9CreatePhase,
  v9FreezePhase,
  v9UnfreezePhase,
  v9ListRuns,
  v9CreateRun,
  v9CancelRun,
  v9ArchiveRun,
  v9ImportFile,
  v9ImportDatabase,
  v9ListMolecules,
  v9MoleculeStats,
  v9BookmarkMolecule,
  v9BookmarkBatch,
} from '../api'

const WorkspaceContext = createContext(null)

export function WorkspaceProvider({ children }) {
  const { isAuthenticated, loading: authLoading } = useAuth()
  const { addToast } = useToast()

  // --- Core data state ---
  const [projects, setProjects] = useState([])
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState(null)

  // --- Core navigation state ---
  const [currentProjectId, setCurrentProjectId] = useState(null)
  const [currentCampaignId, setCurrentCampaignId] = useState(null)
  const [currentPhaseId, setCurrentPhaseId] = useState(null)

  // --- Molecule selection (checkboxes in dashboard) ---
  const [selectedMoleculeIds, setSelectedMoleculeIds] = useState(new Set())

  // --- Bookmarks (persisted per phase in a real app) ---
  const [bookmarkOverrides, setBookmarkOverrides] = useState({})

  // --- Phase status overrides (local until API-backed) ---
  const [phaseStatusOverrides, setPhaseStatusOverrides] = useState({})

  // --- Runs from API ---
  const [phaseRuns, setPhaseRuns] = useState([])
  const [runsLoading, setRunsLoading] = useState(false)

  // --- Paginated molecule store ---
  const [moleculePages, setMoleculePages] = useState([])   // [[...], [...]]
  const [moleculeTotal, setMoleculeTotal] = useState(0)
  const [moleculeNextCursor, setMoleculeNextCursor] = useState(null)
  const [moleculeHasMore, setMoleculeHasMore] = useState(false)
  const [moleculeStats, setMoleculeStats] = useState({ total: 0, bookmarked: 0, ai_generated: 0 })
  const [moleculesLoading, setMoleculesLoading] = useState(false)
  const [serverSort, setServerSort] = useState({ sort_by: 'created_at', sort_dir: 'desc' })
  const [serverFilters, setServerFilters] = useState({ search: '', bookmarked_only: false })
  const loadMoreLockRef = useRef(false)

  // --- Fetch projects from API ---
  const refreshProjects = useCallback(async () => {
    if (!isAuthenticated) {
      setProjects([])
      setLoading(false)
      return
    }
    try {
      setLoading(true)
      setError(null)
      const data = await v9ListProjects()
      setProjects(data)
    } catch (err) {
      console.warn('[WorkspaceContext] API error:', err.message)
      setError(err.message)
      setProjects([])
    } finally {
      setLoading(false)
    }
  }, [isAuthenticated])

  // Wait for auth to settle before fetching projects (avoids race condition)
  useEffect(() => {
    if (authLoading) return // Auth still initializing — wait
    refreshProjects()
  }, [authLoading, refreshProjects])

  // --- CRUD actions ---
  const createProject = useCallback(async (data) => {
    if (!isAuthenticated) {
      addToast('Please sign in to create a project', 'error')
      return null
    }
    try {
      const created = await v9CreateProject(data)
      setProjects(prev => [created, ...prev])
      return created
    } catch (err) {
      addToast(err.userMessage || err.message || 'Failed to create project', 'error')
      return null
    }
  }, [isAuthenticated, addToast])

  const updateProject = useCallback(async (projectId, data) => {
    if (!isAuthenticated) return null
    try {
      const updated = await v9UpdateProject(projectId, data)
      setProjects(prev => prev.map(p => p.id === projectId ? updated : p))
      return updated
    } catch (err) {
      addToast(err.userMessage || err.message || 'Failed to update project', 'error')
      throw err
    }
  }, [isAuthenticated, addToast])

  const deleteProject = useCallback(async (projectId) => {
    if (!isAuthenticated) return
    await v9DeleteProject(projectId)
    setProjects(prev => prev.filter(p => p.id !== projectId))
  }, [isAuthenticated])

  const createCampaign = useCallback(async (projectId, data) => {
    if (!isAuthenticated) return null
    const created = await v9CreateCampaign(projectId, data)
    await refreshProjects()
    return created
  }, [isAuthenticated, refreshProjects])

  const updateCampaign = useCallback(async (campaignId, data) => {
    if (!isAuthenticated) return null
    const updated = await v9UpdateCampaign(campaignId, data)
    await refreshProjects()
    return updated
  }, [isAuthenticated, refreshProjects])

  const createPhase = useCallback(async (campaignId, data) => {
    if (!isAuthenticated) return null
    const created = await v9CreatePhase(campaignId, data)
    await refreshProjects()
    return created
  }, [isAuthenticated, refreshProjects])

  const freezePhase = useCallback(async (phaseId) => {
    if (!isAuthenticated) return
    await v9FreezePhase(phaseId)
    await refreshProjects()
  }, [isAuthenticated, refreshProjects])

  const unfreezePhase = useCallback(async (phaseId) => {
    if (!isAuthenticated) return
    const result = await v9UnfreezePhase(phaseId)
    await refreshProjects()
    return result
  }, [isAuthenticated, refreshProjects])

  // --- Runs API actions ---
  const refreshRuns = useCallback(async (pId) => {
    const id = pId || currentPhaseId
    if (!id || !isAuthenticated) return
    try {
      setRunsLoading(true)
      const data = await v9ListRuns(id)
      setPhaseRuns(data)
    } catch (err) {
      console.warn('[WorkspaceContext] Failed to fetch runs:', err.message)
    } finally {
      setRunsLoading(false)
    }
  }, [currentPhaseId, isAuthenticated])

  // --- Molecules: paginated fetch ---
  const refreshMolecules = useCallback(async (pId, { append = false, cursor: cursorArg } = {}) => {
    const id = pId || currentPhaseId
    if (!id || !isAuthenticated) return
    try {
      if (!append) setMoleculesLoading(true)

      const params = {
        sort_by: serverSort.sort_by,
        sort_dir: serverSort.sort_dir,
        limit: 200,
      }
      if (serverFilters.search) params.search = serverFilters.search
      if (serverFilters.bookmarked_only) params.bookmarked_only = true
      if (append && cursorArg) {
        params.cursor = cursorArg
      }

      const data = await v9ListMolecules(id, params)
      const mols = Array.isArray(data) ? data : data.molecules || []

      if (append) {
        setMoleculePages(prev => [...prev, mols])
      } else {
        setMoleculePages([mols])
      }

      setMoleculeTotal(data.total || mols.length)
      setMoleculeNextCursor(data.next_cursor || null)
      setMoleculeHasMore(data.has_more || false)
    } catch (err) {
      console.warn('[WorkspaceContext] Failed to fetch molecules:', err.message)
    } finally {
      setMoleculesLoading(false)
      loadMoreLockRef.current = false
    }
  }, [currentPhaseId, isAuthenticated, serverSort, serverFilters])

  // --- Stats fetch (lightweight) ---
  const refreshStats = useCallback(async (pId) => {
    const id = pId || currentPhaseId
    if (!id || !isAuthenticated) return
    try {
      const stats = await v9MoleculeStats(id)
      setMoleculeStats(stats)
    } catch (err) {
      console.warn('[WorkspaceContext] Failed to fetch stats:', err.message)
    }
  }, [currentPhaseId, isAuthenticated])

  // --- Load more (infinite scroll trigger) ---
  const loadMoreMolecules = useCallback(() => {
    if (loadMoreLockRef.current || !moleculeHasMore || !moleculeNextCursor) return
    loadMoreLockRef.current = true
    refreshMolecules(currentPhaseId, { append: true, cursor: moleculeNextCursor })
  }, [moleculeHasMore, moleculeNextCursor, currentPhaseId, refreshMolecules])

  // --- Server sort/filter updates (trigger fresh page 1) ---
  const updateServerSort = useCallback((sort_by, sort_dir) => {
    setServerSort({ sort_by, sort_dir })
  }, [])

  const updateServerFilters = useCallback((filters) => {
    setServerFilters(prev => ({ ...prev, ...filters }))
  }, [])

  // Re-fetch page 1 when sort or filters change
  useEffect(() => {
    if (currentPhaseId && isAuthenticated) {
      refreshMolecules(currentPhaseId)
    }
  }, [serverSort, serverFilters]) // eslint-disable-line react-hooks/exhaustive-deps

  const createRun = useCallback(async (phaseId, data) => {
    if (!isAuthenticated) return null
    const created = await v9CreateRun(phaseId, data)
    await refreshRuns(phaseId)
    return created
  }, [isAuthenticated, refreshRuns])

  const cancelRun = useCallback(async (runId) => {
    if (!isAuthenticated) return
    await v9CancelRun(runId)
    await refreshRuns()
  }, [isAuthenticated, refreshRuns])

  const archiveRun = useCallback(async (runId) => {
    if (!isAuthenticated) return
    await v9ArchiveRun(runId)
    await refreshRuns()
  }, [isAuthenticated, refreshRuns])

  const importFile = useCallback(async (phaseId, file) => {
    if (!isAuthenticated) return null
    const created = await v9ImportFile(phaseId, file)
    await refreshRuns(phaseId)
    return created
  }, [isAuthenticated, refreshRuns])

  const importDatabase = useCallback(async (phaseId, config) => {
    if (!isAuthenticated) return null
    const created = await v9ImportDatabase(phaseId, config)
    await refreshRuns(phaseId)
    return created
  }, [isAuthenticated, refreshRuns])

  // --- Run polling: auto-refresh while runs are active ---
  useEffect(() => {
    if (!isAuthenticated || !currentPhaseId) return
    const hasActive = phaseRuns.some(r => r.status === 'created' || r.status === 'running')
    if (!hasActive) return

    const interval = setInterval(() => {
      refreshRuns()
      // During active runs: refresh stats (lightweight) + page 1 only
      refreshStats()
      refreshMolecules(currentPhaseId)
    }, 5000)
    return () => clearInterval(interval)
  }, [isAuthenticated, currentPhaseId, phaseRuns, refreshRuns, refreshMolecules, refreshStats])

  // --- Fetch runs & molecules when phase changes ---
  useEffect(() => {
    if (currentPhaseId && isAuthenticated) {
      refreshRuns(currentPhaseId)
      refreshMolecules(currentPhaseId)
      refreshStats(currentPhaseId)
    } else {
      setPhaseRuns([])
      setMoleculePages([])
      setMoleculeTotal(0)
      setMoleculeNextCursor(null)
      setMoleculeHasMore(false)
      setMoleculeStats({ total: 0, bookmarked: 0, ai_generated: 0 })
    }
  }, [currentPhaseId, isAuthenticated]) // eslint-disable-line react-hooks/exhaustive-deps

  // --- Derived data ---
  const currentProject = useMemo(
    () => projects.find(p => p.id === currentProjectId) || null,
    [projects, currentProjectId]
  )

  const currentCampaign = useMemo(() => {
    if (!currentProject) return null
    if (currentCampaignId) return currentProject.campaigns?.find(c => c.id === currentCampaignId) || null
    return currentProject.campaigns?.[0] || null
  }, [currentProject, currentCampaignId])

  const currentPhase = useMemo(() => {
    if (!currentCampaign) return null
    return currentCampaign.phases?.find(p => p.id === currentPhaseId) || null
  }, [currentCampaign, currentPhaseId])

  // Flatten all pages into a single array, applying bookmark overrides
  const phaseMolecules = useMemo(() => {
    if (!currentPhase) return []
    const allMols = moleculePages.flat()
    return allMols.map(m => ({
      ...m,
      bookmarked: bookmarkOverrides[m.id] !== undefined ? bookmarkOverrides[m.id] : m.bookmarked,
    }))
  }, [currentPhase, moleculePages, bookmarkOverrides])

  // Runs from API
  const currentPhaseRuns = useMemo(() => {
    if (!currentPhase) return []
    return phaseRuns
  }, [currentPhase, phaseRuns])

  // --- Navigation actions ---
  const selectProject = useCallback((projectId) => {
    setCurrentProjectId(projectId)
    setCurrentCampaignId(null)
    setCurrentPhaseId(null)
    setSelectedMoleculeIds(new Set())
  }, [])

  const selectCampaign = useCallback((campaignId) => {
    setCurrentCampaignId(campaignId)
    setCurrentPhaseId(null)
    setSelectedMoleculeIds(new Set())
  }, [])

  const selectPhase = useCallback((phaseId) => {
    setCurrentPhaseId(phaseId)
    setSelectedMoleculeIds(new Set())
    setBookmarkOverrides({})
    // Reset server state for new phase
    setServerSort({ sort_by: 'created_at', sort_dir: 'desc' })
    setServerFilters({ search: '', bookmarked_only: false })
  }, [])

  // --- Selection actions ---
  const toggleSelection = useCallback((molId) => {
    setSelectedMoleculeIds(prev => {
      const next = new Set(prev)
      if (next.has(molId)) next.delete(molId)
      else next.add(molId)
      return next
    })
  }, [])

  const selectAll = useCallback(() => {
    setSelectedMoleculeIds(new Set(phaseMolecules.map(m => m.id)))
  }, [phaseMolecules])

  const selectNone = useCallback(() => {
    setSelectedMoleculeIds(new Set())
  }, [])

  const selectBookmarked = useCallback(() => {
    setSelectedMoleculeIds(new Set(
      phaseMolecules.filter(m => m.bookmarked).map(m => m.id)
    ))
  }, [phaseMolecules])

  // --- Bookmark actions ---
  const toggleBookmark = useCallback(async (molId) => {
    const mol = phaseMolecules.find(m => m.id === molId)
    const currentVal = bookmarkOverrides[molId]
    const originalVal = mol?.bookmarked || false
    const effectiveVal = currentVal !== undefined ? currentVal : originalVal
    const newVal = !effectiveVal

    // Optimistic update
    setBookmarkOverrides(prev => ({ ...prev, [molId]: newVal }))

    if (isAuthenticated) {
      try {
        await v9BookmarkMolecule(molId, newVal)
      } catch (err) {
        console.warn('[WorkspaceContext] Bookmark failed, reverting:', err.message)
        setBookmarkOverrides(prev => ({ ...prev, [molId]: effectiveVal }))
      }
    }
  }, [phaseMolecules, bookmarkOverrides, isAuthenticated])

  const bookmarkSelected = useCallback(async () => {
    const ids = [...selectedMoleculeIds]
    // Optimistic update
    setBookmarkOverrides(prev => {
      const next = { ...prev }
      ids.forEach(id => { next[id] = true })
      return next
    })

    if (isAuthenticated && currentPhaseId && ids.length > 0) {
      try {
        await v9BookmarkBatch(currentPhaseId, ids, true)
      } catch (err) {
        console.warn('[WorkspaceContext] Batch bookmark failed:', err.message)
      }
    }
  }, [selectedMoleculeIds, isAuthenticated, currentPhaseId])

  // --- Phase status getter (for local overrides) ---
  const getPhaseStatus = useCallback((phaseId) => {
    return phaseStatusOverrides[phaseId] || null
  }, [phaseStatusOverrides])

  const value = useMemo(() => ({
    // Data
    projects,
    currentProject,
    currentCampaign,
    currentPhase,
    phaseMolecules,
    currentPhaseRuns,
    loading,
    error,
    runsLoading,
    moleculesLoading,

    // Paginated molecule metadata
    moleculeTotal,
    moleculeHasMore,
    moleculeStats,
    serverSort,
    serverFilters,

    // CRUD
    createProject,
    updateProject,
    deleteProject,
    createCampaign,
    updateCampaign,
    createPhase,
    refreshProjects,

    // Navigation
    selectProject,
    selectCampaign,
    selectPhase,

    // Selection
    selectedMoleculeIds,
    toggleSelection,
    selectAll,
    selectNone,
    selectBookmarked,

    // Bookmarks
    toggleBookmark,
    bookmarkSelected,

    // Phase actions
    freezePhase,
    unfreezePhase,
    getPhaseStatus,

    // Runs & Molecules actions
    createRun,
    cancelRun,
    archiveRun,
    importFile,
    importDatabase,
    refreshRuns,
    refreshMolecules,
    loadMoreMolecules,
    updateServerSort,
    updateServerFilters,
    refreshStats,
  }), [
    projects, currentProject, currentCampaign, currentPhase, phaseMolecules, currentPhaseRuns,
    loading, error, runsLoading, moleculesLoading,
    moleculeTotal, moleculeHasMore, moleculeStats, serverSort, serverFilters,
    createProject, updateProject, deleteProject, createCampaign, updateCampaign, createPhase, refreshProjects,
    selectProject, selectCampaign, selectPhase,
    selectedMoleculeIds, toggleSelection, selectAll, selectNone, selectBookmarked,
    toggleBookmark, bookmarkSelected,
    freezePhase, unfreezePhase, getPhaseStatus,
    createRun, cancelRun, archiveRun, importFile, importDatabase, refreshRuns, refreshMolecules,
    loadMoreMolecules, updateServerSort, updateServerFilters, refreshStats,
  ])

  return (
    <WorkspaceContext.Provider value={value}>
      {children}
    </WorkspaceContext.Provider>
  )
}

export function useWorkspace() {
  const ctx = useContext(WorkspaceContext)
  if (!ctx) throw new Error('useWorkspace must be used within WorkspaceProvider')
  return ctx
}
