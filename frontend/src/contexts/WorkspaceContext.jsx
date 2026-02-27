import React, { createContext, useContext, useState, useCallback, useMemo, useEffect } from 'react'
import { useAuth } from './AuthContext.jsx'
import { MOCK_PROJECTS } from '../mock/data'
import {
  v9ListProjects,
  v9CreateProject,
  v9UpdateProject,
  v9DeleteProject,
  v9CreateCampaign,
  v9CreatePhase,
  v9FreezePhase,
  v9UnfreezePhase,
  v9ListRuns,
  v9CreateRun,
  v9CancelRun,
  v9ArchiveRun,
  v9ImportFile,
  v9ListMolecules,
  v9BookmarkMolecule,
  v9BookmarkBatch,
} from '../api'

const WorkspaceContext = createContext(null)

export function WorkspaceProvider({ children }) {
  const { isAuthenticated } = useAuth()

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

  // --- Runs & Molecules from API ---
  const [phaseRuns, setPhaseRuns] = useState([])
  const [apiMolecules, setApiMolecules] = useState([])
  const [runsLoading, setRunsLoading] = useState(false)
  const [moleculesLoading, setMoleculesLoading] = useState(false)

  // --- Fetch projects from API or fallback to mock ---
  const refreshProjects = useCallback(async () => {
    if (!isAuthenticated) {
      setProjects(MOCK_PROJECTS)
      setLoading(false)
      return
    }
    try {
      setLoading(true)
      setError(null)
      const data = await v9ListProjects()
      setProjects(data)
    } catch (err) {
      console.warn('[WorkspaceContext] API unavailable, using mock data:', err.message)
      setError(err.message)
      setProjects(MOCK_PROJECTS)
    } finally {
      setLoading(false)
    }
  }, [isAuthenticated])

  useEffect(() => {
    refreshProjects()
  }, [refreshProjects])

  // --- CRUD actions ---
  const createProject = useCallback(async (data) => {
    if (!isAuthenticated) {
      // Mock fallback
      const fake = {
        id: `proj-${Date.now()}`,
        ...data,
        status: 'active',
        created_at: new Date().toISOString(),
        updated_at: new Date().toISOString(),
        campaigns: [],
      }
      setProjects(prev => [fake, ...prev])
      return fake
    }
    const created = await v9CreateProject(data)
    setProjects(prev => [created, ...prev])
    return created
  }, [isAuthenticated])

  const updateProject = useCallback(async (projectId, data) => {
    if (!isAuthenticated) {
      setProjects(prev => prev.map(p => p.id === projectId ? { ...p, ...data } : p))
      return
    }
    const updated = await v9UpdateProject(projectId, data)
    setProjects(prev => prev.map(p => p.id === projectId ? updated : p))
    return updated
  }, [isAuthenticated])

  const deleteProject = useCallback(async (projectId) => {
    if (!isAuthenticated) {
      setProjects(prev => prev.filter(p => p.id !== projectId))
      return
    }
    await v9DeleteProject(projectId)
    setProjects(prev => prev.filter(p => p.id !== projectId))
  }, [isAuthenticated])

  const createCampaign = useCallback(async (projectId, data) => {
    if (!isAuthenticated) return null
    const created = await v9CreateCampaign(projectId, data)
    await refreshProjects()
    return created
  }, [isAuthenticated, refreshProjects])

  const createPhase = useCallback(async (campaignId, data) => {
    if (!isAuthenticated) return null
    const created = await v9CreatePhase(campaignId, data)
    await refreshProjects()
    return created
  }, [isAuthenticated, refreshProjects])

  const freezePhase = useCallback(async (phaseId) => {
    if (!isAuthenticated) {
      setPhaseStatusOverrides(prev => ({ ...prev, [phaseId]: 'frozen' }))
      return
    }
    await v9FreezePhase(phaseId)
    await refreshProjects()
  }, [isAuthenticated, refreshProjects])

  const unfreezePhase = useCallback(async (phaseId) => {
    if (!isAuthenticated) {
      setPhaseStatusOverrides(prev => ({ ...prev, [phaseId]: 'active' }))
      return
    }
    const result = await v9UnfreezePhase(phaseId)
    await refreshProjects()
    return result
  }, [isAuthenticated, refreshProjects])

  // --- Runs & Molecules API actions ---
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

  const refreshMolecules = useCallback(async (pId) => {
    const id = pId || currentPhaseId
    if (!id || !isAuthenticated) return
    try {
      setMoleculesLoading(true)
      const data = await v9ListMolecules(id)
      setApiMolecules(Array.isArray(data) ? data : data.molecules || [])
    } catch (err) {
      console.warn('[WorkspaceContext] Failed to fetch molecules:', err.message)
    } finally {
      setMoleculesLoading(false)
    }
  }, [currentPhaseId, isAuthenticated])

  const createRun = useCallback(async (phaseId, data) => {
    if (!isAuthenticated) {
      console.log('[mock] createRun:', data)
      return null
    }
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

  // --- Run polling: auto-refresh while runs are active ---
  useEffect(() => {
    if (!isAuthenticated || !currentPhaseId) return
    const hasActive = phaseRuns.some(r => r.status === 'created' || r.status === 'running')
    if (!hasActive) return

    const interval = setInterval(() => {
      refreshRuns()
      refreshMolecules()
    }, 5000)
    return () => clearInterval(interval)
  }, [isAuthenticated, currentPhaseId, phaseRuns, refreshRuns, refreshMolecules])

  // --- Fetch runs & molecules when phase changes ---
  useEffect(() => {
    if (currentPhaseId && isAuthenticated) {
      refreshRuns(currentPhaseId)
      refreshMolecules(currentPhaseId)
    } else {
      setPhaseRuns([])
      setApiMolecules([])
    }
  }, [currentPhaseId, isAuthenticated])

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

  // Molecules: API when authenticated, mock otherwise
  const phaseMolecules = useMemo(() => {
    if (!currentPhase) return []
    const mols = isAuthenticated && apiMolecules.length > 0
      ? apiMolecules
      : currentPhase.molecules || []
    return mols.map(m => ({
      ...m,
      bookmarked: bookmarkOverrides[m.id] !== undefined ? bookmarkOverrides[m.id] : m.bookmarked,
    }))
  }, [currentPhase, isAuthenticated, apiMolecules, bookmarkOverrides])

  // Runs: API when authenticated, mock otherwise
  const currentPhaseRuns = useMemo(() => {
    if (!currentPhase) return []
    return isAuthenticated ? phaseRuns : (currentPhase.runs || [])
  }, [currentPhase, isAuthenticated, phaseRuns])

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

    // CRUD
    createProject,
    updateProject,
    deleteProject,
    createCampaign,
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
    refreshRuns,
    refreshMolecules,
  }), [
    projects, currentProject, currentCampaign, currentPhase, phaseMolecules, currentPhaseRuns,
    loading, error, runsLoading, moleculesLoading,
    createProject, updateProject, deleteProject, createCampaign, createPhase, refreshProjects,
    selectProject, selectCampaign, selectPhase,
    selectedMoleculeIds, toggleSelection, selectAll, selectNone, selectBookmarked,
    toggleBookmark, bookmarkSelected,
    freezePhase, unfreezePhase, getPhaseStatus,
    createRun, cancelRun, archiveRun, importFile, refreshRuns, refreshMolecules,
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
