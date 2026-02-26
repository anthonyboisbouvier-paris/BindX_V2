import React, { createContext, useContext, useState, useCallback, useMemo } from 'react'
import { MOCK_PROJECTS } from '../mock/data'

const WorkspaceContext = createContext(null)

export function WorkspaceProvider({ children }) {
  // --- Core navigation state ---
  const [projects] = useState(MOCK_PROJECTS)
  const [currentProjectId, setCurrentProjectId] = useState(null)
  const [currentCampaignId, setCurrentCampaignId] = useState(null)
  const [currentPhaseId, setCurrentPhaseId] = useState(null)

  // --- Molecule selection (checkboxes in dashboard) ---
  const [selectedMoleculeIds, setSelectedMoleculeIds] = useState(new Set())

  // --- Bookmarks (persisted per phase in a real app) ---
  const [bookmarkOverrides, setBookmarkOverrides] = useState({})

  // --- Derived data ---
  const currentProject = useMemo(
    () => projects.find(p => p.id === currentProjectId) || null,
    [projects, currentProjectId]
  )

  const currentCampaign = useMemo(() => {
    if (!currentProject) return null
    if (currentCampaignId) return currentProject.campaigns.find(c => c.id === currentCampaignId) || null
    return currentProject.campaigns[0] || null
  }, [currentProject, currentCampaignId])

  const currentPhase = useMemo(() => {
    if (!currentCampaign) return null
    return currentCampaign.phases.find(p => p.id === currentPhaseId) || null
  }, [currentCampaign, currentPhaseId])

  // Molecules with bookmark overrides applied
  const phaseMolecules = useMemo(() => {
    if (!currentPhase) return []
    return currentPhase.molecules.map(m => ({
      ...m,
      bookmarked: bookmarkOverrides[m.id] !== undefined ? bookmarkOverrides[m.id] : m.bookmarked,
    }))
  }, [currentPhase, bookmarkOverrides])

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
    if (!currentPhase) return
    setSelectedMoleculeIds(new Set(currentPhase.molecules.map(m => m.id)))
  }, [currentPhase])

  const selectNone = useCallback(() => {
    setSelectedMoleculeIds(new Set())
  }, [])

  const selectBookmarked = useCallback(() => {
    setSelectedMoleculeIds(new Set(
      phaseMolecules.filter(m => m.bookmarked).map(m => m.id)
    ))
  }, [phaseMolecules])

  // --- Bookmark actions ---
  const toggleBookmark = useCallback((molId) => {
    setBookmarkOverrides(prev => {
      const currentVal = prev[molId]
      const mol = currentPhase?.molecules.find(m => m.id === molId)
      const originalVal = mol?.bookmarked || false
      const effectiveVal = currentVal !== undefined ? currentVal : originalVal
      return { ...prev, [molId]: !effectiveVal }
    })
  }, [currentPhase])

  const bookmarkSelected = useCallback(() => {
    setBookmarkOverrides(prev => {
      const next = { ...prev }
      selectedMoleculeIds.forEach(id => { next[id] = true })
      return next
    })
  }, [selectedMoleculeIds])

  // --- Phase actions ---
  // These are mocked — in prod they'd call Supabase
  const [phaseStatusOverrides, setPhaseStatusOverrides] = useState({})

  const freezePhase = useCallback((phaseId) => {
    setPhaseStatusOverrides(prev => ({ ...prev, [phaseId]: 'frozen' }))
  }, [])

  const unfreezePhase = useCallback((phaseId) => {
    setPhaseStatusOverrides(prev => ({ ...prev, [phaseId]: 'active' }))
  }, [])

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
  }), [
    projects, currentProject, currentCampaign, currentPhase, phaseMolecules,
    selectProject, selectCampaign, selectPhase,
    selectedMoleculeIds, toggleSelection, selectAll, selectNone, selectBookmarked,
    toggleBookmark, bookmarkSelected,
    freezePhase, unfreezePhase, getPhaseStatus,
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
