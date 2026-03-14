/**
 * BindX V9 — Workspace Zustand store.
 *
 * Single store managing projects, navigation, molecules, runs,
 * selection, and bookmarks. Replaces the old WorkspaceContext.
 *
 * Auth integration: call _syncAuth(isAuthenticated) from WorkspaceBridge.
 * Toast integration: call _setToast(addToast) from WorkspaceBridge.
 */

import { create } from 'zustand'
import {
  v9ListProjects,
  v9CreateProject,
  v9UpdateProject,
  v9DeleteProject,
  v9CreateCampaign,
  v9UpdateCampaign,
  v9CreatePhase,
  v9DeletePhase,
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

// ---------------------------------------------------------------------------
// Store
// ---------------------------------------------------------------------------

const useWorkspaceStore = create((set, get) => ({
  // ===== Bridge state (synced from React) =====
  _isAuthenticated: false,
  _addToast: null,
  _syncAuth: (isAuth) => set({ _isAuthenticated: isAuth }),
  _setToast: (fn) => set({ _addToast: fn }),

  // ===== Projects =====
  projects: [],
  loading: true,
  error: null,

  refreshProjects: async () => {
    const { _isAuthenticated } = get()
    if (!_isAuthenticated) {
      set({ projects: [], loading: false })
      return
    }
    try {
      set({ loading: true, error: null })
      const data = await v9ListProjects()
      set({ projects: data, loading: false })
    } catch (err) {
      console.warn('[Workspace] API error:', err.message)
      set({ projects: [], error: err.message, loading: false })
    }
  },

  createProject: async (data) => {
    const { _isAuthenticated, _addToast } = get()
    if (!_isAuthenticated) {
      _addToast?.('Please sign in to create a project', 'error')
      return null
    }
    try {
      const created = await v9CreateProject(data)
      set(s => ({ projects: [created, ...s.projects] }))
      return created
    } catch (err) {
      _addToast?.(err.userMessage || err.message || 'Failed to create project', 'error')
      return null
    }
  },

  updateProject: async (projectId, data) => {
    const { _isAuthenticated, _addToast } = get()
    if (!_isAuthenticated) return null
    try {
      const updated = await v9UpdateProject(projectId, data)
      set(s => ({ projects: s.projects.map(p => p.id === projectId ? updated : p) }))
      return updated
    } catch (err) {
      _addToast?.(err.userMessage || err.message || 'Failed to update project', 'error')
      throw err
    }
  },

  deleteProject: async (projectId) => {
    if (!get()._isAuthenticated) return
    try {
      await v9DeleteProject(projectId)
      set(s => ({ projects: s.projects.filter(p => p.id !== projectId) }))
    } catch (err) {
      get()._addToast?.(err.userMessage || err.message || 'Failed to delete project', 'error')
      throw err
    }
  },

  createCampaign: async (projectId, data) => {
    if (!get()._isAuthenticated) return null
    try {
      const created = await v9CreateCampaign(projectId, data)
      await get().refreshProjects()
      return created
    } catch (err) {
      get()._addToast?.(err.userMessage || err.message || 'Failed to create campaign', 'error')
      return null
    }
  },

  updateCampaign: async (campaignId, data) => {
    if (!get()._isAuthenticated) return null
    try {
      const updated = await v9UpdateCampaign(campaignId, data)
      await get().refreshProjects()
      return updated
    } catch (err) {
      get()._addToast?.(err.userMessage || err.message || 'Failed to update campaign', 'error')
      return null
    }
  },

  createPhase: async (campaignId, data) => {
    if (!get()._isAuthenticated) return null
    try {
      const created = await v9CreatePhase(campaignId, data)
      await get().refreshProjects()
      return created
    } catch (err) {
      get()._addToast?.(err.userMessage || err.message || 'Failed to create phase', 'error')
      return null
    }
  },

  freezePhase: async (phaseId) => {
    if (!get()._isAuthenticated) return
    try {
      await v9FreezePhase(phaseId)
      await get().refreshProjects()
    } catch (err) {
      get()._addToast?.(err.userMessage || err.message || 'Failed to freeze phase', 'error')
    }
  },

  unfreezePhase: async (phaseId) => {
    if (!get()._isAuthenticated) return
    try {
      const result = await v9UnfreezePhase(phaseId)
      await get().refreshProjects()
      return result
    } catch (err) {
      get()._addToast?.(err.userMessage || err.message || 'Failed to unfreeze phase', 'error')
      return null
    }
  },

  deletePhase: async (phaseId) => {
    if (!get()._isAuthenticated) return
    await v9DeletePhase(phaseId)
    // Clear current phase immediately to stop polling on deleted phase
    set({ currentPhaseId: null, molecules: [], currentPhaseRuns: [] })
    await get().refreshProjects()
  },

  // ===== Navigation =====
  currentProjectId: null,
  currentCampaignId: null,
  currentPhaseId: null,

  selectProject: (projectId) => {
    if (get().currentProjectId === projectId) return
    set({
      currentProjectId: projectId,
      currentCampaignId: null,
      currentPhaseId: null,
      selectedMoleculeIds: new Set(),
    })
  },

  selectCampaign: (campaignId) => {
    if (get().currentCampaignId === campaignId) return
    set({
      currentCampaignId: campaignId,
      currentPhaseId: null,
      selectedMoleculeIds: new Set(),
    })
  },

  selectPhase: (phaseId) => set({ currentPhaseId: phaseId }),

  // ===== Derived (computed from state) =====
  // NOTE: Zustand doesn't have native computed — we use getter functions
  // that consumers call, or compute inline in selectors.

  getCurrentProject: () => {
    const { projects, currentProjectId } = get()
    return projects.find(p => p.id === currentProjectId) || null
  },

  getCurrentCampaign: () => {
    const project = get().getCurrentProject()
    if (!project) return null
    const { currentCampaignId } = get()
    if (currentCampaignId) return project.campaigns?.find(c => c.id === currentCampaignId) || null
    return project.campaigns?.[0] || null
  },

  getCurrentPhase: () => {
    const campaign = get().getCurrentCampaign()
    if (!campaign) return null
    return campaign.phases?.find(p => p.id === get().currentPhaseId) || null
  },

  // ===== Molecule selection =====
  selectedMoleculeIds: new Set(),

  toggleSelection: (molId) => set(s => {
    const next = new Set(s.selectedMoleculeIds)
    if (next.has(molId)) next.delete(molId)
    else next.add(molId)
    return { selectedMoleculeIds: next }
  }),

  selectAll: () => {
    const mols = get().getPhaseMolecules()
    set({ selectedMoleculeIds: new Set(mols.map(m => m.id)) })
  },

  selectNone: () => set({ selectedMoleculeIds: new Set() }),

  selectBookmarked: () => {
    const mols = get().getPhaseMolecules()
    set({ selectedMoleculeIds: new Set(mols.filter(m => m.bookmarked).map(m => m.id)) })
  },

  // ===== Bookmarks =====
  bookmarkOverrides: {},

  toggleBookmark: async (molId) => {
    const { bookmarkOverrides, _isAuthenticated } = get()
    const mols = get().getPhaseMolecules()
    const mol = mols.find(m => m.id === molId)
    const currentVal = bookmarkOverrides[molId]
    const originalVal = mol?.bookmarked || false
    const effectiveVal = currentVal !== undefined ? currentVal : originalVal
    const newVal = !effectiveVal

    // Optimistic update
    set(s => ({ bookmarkOverrides: { ...s.bookmarkOverrides, [molId]: newVal } }))

    if (_isAuthenticated) {
      try {
        await v9BookmarkMolecule(molId, newVal)
      } catch (err) {
        console.warn('[Workspace] Bookmark failed, reverting:', err.message)
        set(s => ({ bookmarkOverrides: { ...s.bookmarkOverrides, [molId]: effectiveVal } }))
      }
    }
  },

  bookmarkSelected: async () => {
    const { selectedMoleculeIds, _isAuthenticated, currentPhaseId, bookmarkOverrides: prev } = get()
    const ids = [...selectedMoleculeIds]

    // Optimistic update
    set(s => {
      const next = { ...s.bookmarkOverrides }
      ids.forEach(id => { next[id] = true })
      return { bookmarkOverrides: next }
    })

    if (_isAuthenticated && currentPhaseId && ids.length > 0) {
      try {
        await v9BookmarkBatch(currentPhaseId, ids, true)
      } catch (err) {
        console.warn('[Workspace] Batch bookmark failed, reverting:', err.message)
        set({ bookmarkOverrides: prev })
        get()._addToast?.('Batch bookmark failed', 'error')
      }
    }
  },

  // ===== Phase status overrides =====
  phaseStatusOverrides: {},

  getPhaseStatus: (phaseId) => {
    return get().phaseStatusOverrides[phaseId] || null
  },

  // ===== Runs =====
  phaseRuns: [],
  runsLoading: false,
  _runsRequestId: 0,

  refreshRuns: async (pId) => {
    const id = pId || get().currentPhaseId
    if (!id || !get()._isAuthenticated) return
    const requestId = (get()._runsRequestId || 0) + 1
    set({ _runsRequestId: requestId, runsLoading: true })
    try {
      const data = await v9ListRuns(id)
      // Discard stale response if a newer request was initiated
      if (get()._runsRequestId !== requestId) return
      set({ phaseRuns: data, runsLoading: false })
    } catch (err) {
      if (get()._runsRequestId !== requestId) return
      console.warn('[Workspace] Failed to fetch runs:', err.message)
      set({ runsLoading: false })
    }
  },

  createRun: async (phaseId, data) => {
    if (!get()._isAuthenticated) return null
    const created = await v9CreateRun(phaseId, data)
    await get().refreshRuns(phaseId)
    return created
  },

  cancelRun: async (runId) => {
    if (!get()._isAuthenticated) return
    await v9CancelRun(runId)
    await get().refreshRuns()
  },

  archiveRun: async (runId) => {
    if (!get()._isAuthenticated) return
    await v9ArchiveRun(runId)
    await get().refreshRuns()
  },

  importFile: async (phaseId, file) => {
    if (!get()._isAuthenticated) return null
    const created = await v9ImportFile(phaseId, file)
    await get().refreshRuns(phaseId)
    return created
  },

  importDatabase: async (phaseId, config) => {
    if (!get()._isAuthenticated) return null
    const created = await v9ImportDatabase(phaseId, config)
    await get().refreshRuns(phaseId)
    return created
  },

  // ===== Molecules (page-based pagination) =====
  molecules: [],
  moleculeTotal: 0,
  moleculePage: 1,
  moleculePageSize: 100,
  moleculeStats: { total: 0, bookmarked: 0, ai_generated: 0, column_ranges: {} },
  moleculesLoading: false,
  serverSort: { sort_by: 'created_at', sort_dir: 'desc' },
  serverFilters: { search: '', bookmarked_only: false },
  _molRequestId: 0,

  refreshMolecules: async (pId, { page } = {}) => {
    const id = pId || get().currentPhaseId
    if (!id || !get()._isAuthenticated) return
    // Increment request counter — only the latest request's response is applied
    const requestId = (get()._molRequestId || 0) + 1
    set({ _molRequestId: requestId, moleculesLoading: true })
    try {
      const { serverSort, serverFilters, moleculePage, moleculePageSize } = get()
      const effectivePage = page || moleculePage
      const params = {
        sort_by: serverSort.sort_by,
        sort_dir: serverSort.sort_dir,
        limit: moleculePageSize,
        offset: (effectivePage - 1) * moleculePageSize,
      }
      if (serverFilters.search) params.search = serverFilters.search
      if (serverFilters.bookmarked_only) params.bookmarked_only = true

      const data = await v9ListMolecules(id, params)

      // Discard stale response if a newer request was initiated
      if (get()._molRequestId !== requestId) return

      const mols = Array.isArray(data) ? data : data.molecules || []

      set({
        molecules: mols,
        moleculeTotal: data.total || mols.length,
        moleculePage: effectivePage,
        moleculesLoading: false,
        bookmarkOverrides: {},  // Clear overrides — server data is source of truth
      })
    } catch (err) {
      if (get()._molRequestId !== requestId) return
      console.warn('[Workspace] Failed to fetch molecules:', err.message)
      set({ moleculesLoading: false })
    }
  },

  refreshStats: async (pId) => {
    const id = pId || get().currentPhaseId
    if (!id || !get()._isAuthenticated) return
    try {
      const stats = await v9MoleculeStats(id)
      set({ moleculeStats: stats })
    } catch (err) {
      console.warn('[Workspace] Failed to fetch stats:', err.message)
    }
  },

  goToPage: (page) => {
    const { moleculeTotal, moleculePageSize, currentPhaseId } = get()
    const maxPage = Math.max(1, Math.ceil(moleculeTotal / moleculePageSize))
    const clamped = Math.max(1, Math.min(page, maxPage))
    set({ moleculePage: clamped, selectedMoleculeIds: new Set() })
    get().refreshMolecules(currentPhaseId, { page: clamped })
  },

  setPageSize: (size) => {
    set({ moleculePageSize: size, moleculePage: 1, selectedMoleculeIds: new Set() })
    get().refreshMolecules(get().currentPhaseId, { page: 1 })
  },

  updateServerSort: (sort_by, sort_dir) => {
    set({ serverSort: { sort_by, sort_dir }, moleculePage: 1, selectedMoleculeIds: new Set() })
    get().refreshMolecules(get().currentPhaseId, { page: 1 })
  },

  updateServerFilters: (filters) => {
    set(s => ({
      serverFilters: { ...s.serverFilters, ...filters },
      moleculePage: 1,
      selectedMoleculeIds: new Set(),
    }))
    get().refreshMolecules(get().currentPhaseId, { page: 1 })
  },

  // ===== Derived getters =====
  getPhaseMolecules: () => {
    const { molecules, bookmarkOverrides } = get()
    const phase = get().getCurrentPhase()
    if (!phase) return []
    return molecules.map(m => ({
      ...m,
      bookmarked: bookmarkOverrides[m.id] !== undefined ? bookmarkOverrides[m.id] : m.bookmarked,
    }))
  },

  getCurrentPhaseRuns: () => {
    const phase = get().getCurrentPhase()
    if (!phase) return []
    return get().phaseRuns
  },

  // ===== Phase change handler =====
  _onPhaseChange: async (phaseId) => {
    if (phaseId && get()._isAuthenticated) {
      set({
        moleculePage: 1,
        selectedMoleculeIds: new Set(),
        bookmarkOverrides: {},
        serverSort: { sort_by: 'created_at', sort_dir: 'desc' },
        serverFilters: { search: '', bookmarked_only: false },
      })
      get().refreshRuns(phaseId)
      get().refreshMolecules(phaseId, { page: 1 })
      get().refreshStats(phaseId)
    } else {
      set({
        phaseRuns: [],
        molecules: [],
        moleculeTotal: 0,
        moleculePage: 1,
        moleculeStats: { total: 0, bookmarked: 0, ai_generated: 0, column_ranges: {} },
      })
    }
  },
}))

export default useWorkspaceStore
