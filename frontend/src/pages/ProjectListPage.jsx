import React, { useState, useMemo, useCallback, useEffect } from 'react'
import { useNavigate } from 'react-router-dom'
import { useWorkspace } from '../contexts/WorkspaceContext.jsx'
import { MOCK_ACTIVITY, PHASE_TYPES } from '../mock/data.js'
import ActivityTimeline from '../components/ActivityTimeline.jsx'

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function relativeTime(iso) {
  if (!iso) return ''
  const diff = Date.now() - new Date(iso).getTime()
  const mins = Math.floor(diff / 60000)
  if (mins < 60) return `${mins}m ago`
  const hrs = Math.floor(mins / 60)
  if (hrs < 24) return `${hrs}h ago`
  const days = Math.floor(hrs / 24)
  return `${days}d ago`
}

function highlight(text, query) {
  if (!query || !text) return text
  const idx = text.toLowerCase().indexOf(query.toLowerCase())
  if (idx === -1) return text
  return (
    <>
      {text.slice(0, idx)}
      <mark className="bg-yellow-100 text-yellow-800 rounded px-0.5">{text.slice(idx, idx + query.length)}</mark>
      {text.slice(idx + query.length)}
    </>
  )
}

// ---------------------------------------------------------------------------
// Phase progress indicator (3 circles + lines)
// ---------------------------------------------------------------------------

function PhaseProgress({ campaigns }) {
  const allPhases = (campaigns || []).flatMap(c => c.phases || [])

  const phaseStatus = (type) => {
    const phase = allPhases.find(p => p.type === type)
    if (!phase || !phase.created_at) return 'empty'
    const hasMols = (phase.stats?.total_molecules || 0) > 0
    const hasRuns = (phase.stats?.runs_completed || 0) > 0
    if (phase.status === 'frozen' && hasMols) return 'done'
    if (hasMols && hasRuns) return 'active'
    if (hasMols) return 'partial'
    return 'empty'
  }

  const types = ['hit_discovery', 'hit_to_lead', 'lead_optimization']
  const labels = ['A', 'B', 'C']

  return (
    <div className="flex items-center gap-1">
      {types.map((type, i) => {
        const status = phaseStatus(type)
        const isLast = i === types.length - 1
        const meta = PHASE_TYPES[type]

        const circleClass = {
          done: 'bg-[#22c55e] border-[#22c55e] text-white',
          active: 'bg-[#22c55e]/20 border-[#22c55e] text-[#22c55e]',
          partial: 'bg-blue-100 border-blue-400 text-blue-600',
          empty: 'bg-white border-gray-200 text-gray-300',
        }[status]

        const lineClass = status === 'done' || status === 'active'
          ? 'bg-[#22c55e]/40'
          : 'bg-gray-100'

        return (
          <React.Fragment key={type}>
            <div className="flex flex-col items-center gap-0.5">
              <div className={`w-6 h-6 rounded-full border-2 flex items-center justify-center text-xs font-bold transition-colors ${circleClass}`}>
                {status === 'done' ? (
                  <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={3} d="M5 13l4 4L19 7" />
                  </svg>
                ) : labels[i]}
              </div>
              <span className="text-xs text-gray-400" style={{ fontSize: 10 }}>{labels[i]}</span>
            </div>
            {!isLast && (
              <div className={`h-0.5 w-5 mb-3.5 rounded ${lineClass}`} />
            )}
          </React.Fragment>
        )
      })}
    </div>
  )
}

// ---------------------------------------------------------------------------
// Project status badge
// ---------------------------------------------------------------------------

function StatusBadge({ project }) {
  const allPhases = (project.campaigns || []).flatMap(c => c.phases || [])
  const totalMols = allPhases.reduce((s, p) => s + (p.stats?.total_molecules || 0), 0)

  if (totalMols === 0) {
    return (
      <span className="inline-flex items-center gap-1.5 px-2 py-0.5 rounded-full text-xs bg-gray-100 text-gray-400">
        <span className="w-1.5 h-1.5 rounded-full bg-gray-300" />
        Empty
      </span>
    )
  }
  if (project.status === 'active') {
    return (
      <span className="inline-flex items-center gap-1.5 px-2 py-0.5 rounded-full text-xs bg-green-50 text-green-700 font-medium">
        <span className="w-1.5 h-1.5 rounded-full bg-[#22c55e] animate-pulse" />
        Active
      </span>
    )
  }
  return (
    <span className="inline-flex items-center gap-1.5 px-2 py-0.5 rounded-full text-xs bg-gray-100 text-gray-500">
      <span className="w-1.5 h-1.5 rounded-full bg-gray-400" />
      {project.status || 'Unknown'}
    </span>
  )
}

// ---------------------------------------------------------------------------
// Project card
// ---------------------------------------------------------------------------

function ProjectCard({ project, onClick, searchQuery }) {
  const allPhases = (project.campaigns || []).flatMap(c => c.phases || [])
  const totalMolecules = allPhases.reduce((s, p) => s + (p.stats?.total_molecules || 0), 0)
  const totalRuns = allPhases.reduce((s, p) => s + (p.stats?.runs_completed || 0) + (p.stats?.runs_running || 0), 0)

  return (
    <button
      onClick={() => onClick(project.id)}
      className="group bg-white rounded-xl border border-gray-100 shadow-sm hover:shadow-lg hover:-translate-y-0.5 transition-all duration-200 text-left w-full overflow-hidden focus:outline-none focus:ring-2 focus:ring-[#1e3a5f]/30 focus:ring-offset-2"
    >
      {/* Top accent bar: blue=active, gray=empty */}
      <div className={`h-1 w-full ${
        totalMolecules > 0 && project.status === 'active'
          ? 'bg-gradient-to-r from-[#1e3a5f] to-[#22c55e]'
          : 'bg-gray-100'
      }`} />

      <div className="p-5">
        {/* Header */}
        <div className="flex items-start justify-between gap-2 mb-3">
          <h3 className="font-bold text-[#1e3a5f] text-base leading-tight group-hover:text-[#22c55e] transition-colors">
            {highlight(project.name, searchQuery)}
          </h3>
          <StatusBadge project={project} />
        </div>

        {/* Target row */}
        <div className="flex items-center gap-2 mb-3 flex-wrap">
          {project.target_name && (
            <span className="inline-flex items-center gap-1 text-xs text-[#1e3a5f] font-semibold">
              <svg className="w-3.5 h-3.5 text-gray-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <circle cx="12" cy="12" r="9" strokeWidth={2} />
                <circle cx="12" cy="12" r="5" strokeWidth={2} />
                <circle cx="12" cy="12" r="1.5" fill="currentColor" strokeWidth={0} />
              </svg>
              {highlight(project.target_name, searchQuery)}
            </span>
          )}
          {project.uniprot_id && (
            <span className="text-xs font-mono text-gray-400 bg-gray-50 px-1.5 py-0.5 rounded border border-gray-100">
              {project.uniprot_id}
            </span>
          )}
          {project.target_pdb_id && (
            <span className="text-xs font-mono text-gray-400">
              PDB: {project.target_pdb_id}
            </span>
          )}
        </div>

        {/* Phase progress */}
        <div className="mb-4">
          <PhaseProgress campaigns={project.campaigns} />
        </div>

        {/* Stats row */}
        <div className="grid grid-cols-3 gap-2 mb-3 pt-3 border-t border-gray-50">
          <div className="text-center">
            <p className="text-base font-bold text-[#1e3a5f]">{totalMolecules}</p>
            <p className="text-xs text-gray-400">molecules</p>
          </div>
          <div className="text-center">
            <p className="text-base font-bold text-gray-700">{totalRuns}</p>
            <p className="text-xs text-gray-400">runs</p>
          </div>
          <div className="text-center">
            <p className="text-xs font-semibold text-gray-500 truncate">
              {project.updated_at ? relativeTime(project.updated_at) : '—'}
            </p>
            <p className="text-xs text-gray-400">last update</p>
          </div>
        </div>

        {/* Description (if searched) */}
        {searchQuery && project.description?.toLowerCase().includes(searchQuery.toLowerCase()) && (
          <p className="text-xs text-gray-400 line-clamp-1 italic">
            {highlight(project.description, searchQuery)}
          </p>
        )}

        {/* Footer: campaign count + chevron */}
        <div className="flex items-center justify-between mt-1">
          <span className="text-xs text-gray-400">
            {project.campaigns.length} campaign{project.campaigns.length !== 1 ? 's' : ''}
          </span>
          <svg className="w-4 h-4 text-gray-300 group-hover:text-[#22c55e] transition-colors" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
          </svg>
        </div>
      </div>
    </button>
  )
}

// ---------------------------------------------------------------------------
// New Project modal
// ---------------------------------------------------------------------------

function NewProjectModal({ onClose, onCreated }) {
  const [name, setName] = useState('')
  const [targetName, setTargetName] = useState('')
  const [pdbId, setPdbId] = useState('')
  const [uniprotId, setUniprotId] = useState('')
  const [description, setDescription] = useState('')
  const [saving, setSaving] = useState(false)

  useEffect(() => {
    const handler = (e) => { if (e.key === 'Escape') onClose() }
    document.addEventListener('keydown', handler)
    return () => document.removeEventListener('keydown', handler)
  }, [onClose])

  const handleSubmit = async (e) => {
    e.preventDefault()
    if (!name.trim() || !targetName.trim()) return
    setSaving(true)
    try {
      const project = await onCreated({
        name: name.trim(),
        description: description.trim() || undefined,
        target_name: targetName.trim(),
        target_pdb_id: pdbId.trim() || undefined,
        target_input_type: uniprotId.trim() ? 'uniprot_id' : 'pdb_id',
        target_input_value: uniprotId.trim() || pdbId.trim() || targetName.trim(),
      })
      return project
    } catch (err) {
      console.error('[NewProjectModal] Create failed:', err)
      setSaving(false)
    }
  }

  return (
    <div
      className="fixed inset-0 bg-black/40 backdrop-blur-sm flex items-center justify-center z-50 px-4"
      onClick={(e) => { if (e.target === e.currentTarget) onClose() }}
    >
      <div className="bg-white rounded-2xl shadow-2xl w-full max-w-lg overflow-hidden">
        {/* Header */}
        <div className="bg-[#1e3a5f] px-6 py-4 flex items-center justify-between">
          <div>
            <h3 className="text-white font-bold text-base">New Project</h3>
            <p className="text-blue-300 text-xs mt-0.5">Set up a new drug discovery project</p>
          </div>
          <button
            onClick={onClose}
            className="w-8 h-8 rounded-full bg-white/10 hover:bg-white/20 flex items-center justify-center text-white transition-colors"
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
            </svg>
          </button>
        </div>

        {/* Form */}
        <form onSubmit={handleSubmit} className="px-6 py-5 space-y-4">
          {/* Name */}
          <div>
            <label className="block text-xs font-semibold text-gray-500 uppercase tracking-wider mb-1.5">
              Project Name <span className="text-red-400">*</span>
            </label>
            <input
              value={name}
              onChange={e => setName(e.target.value)}
              required
              placeholder="e.g. EGFR Inhibitors Phase III"
              className="w-full text-sm px-3 py-2.5 rounded-lg border border-gray-200 focus:outline-none focus:ring-2 focus:ring-[#1e3a5f]/20 focus:border-[#1e3a5f] transition-colors"
            />
          </div>

          {/* Target + PDB/UniProt row */}
          <div className="grid grid-cols-3 gap-3">
            <div className="col-span-1">
              <label className="block text-xs font-semibold text-gray-500 uppercase tracking-wider mb-1.5">
                Target <span className="text-red-400">*</span>
              </label>
              <input
                value={targetName}
                onChange={e => setTargetName(e.target.value)}
                required
                placeholder="e.g. EGFR"
                className="w-full text-sm px-3 py-2.5 rounded-lg border border-gray-200 focus:outline-none focus:ring-2 focus:ring-[#1e3a5f]/20 focus:border-[#1e3a5f] transition-colors"
              />
            </div>
            <div>
              <label className="block text-xs font-semibold text-gray-500 uppercase tracking-wider mb-1.5">
                PDB ID
              </label>
              <input
                value={pdbId}
                onChange={e => setPdbId(e.target.value.toUpperCase())}
                placeholder="e.g. 1M17"
                className="w-full text-sm px-3 py-2.5 rounded-lg border border-gray-200 font-mono focus:outline-none focus:ring-2 focus:ring-[#1e3a5f]/20 focus:border-[#1e3a5f] transition-colors"
              />
            </div>
            <div>
              <label className="block text-xs font-semibold text-gray-500 uppercase tracking-wider mb-1.5">
                UniProt ID
              </label>
              <input
                value={uniprotId}
                onChange={e => setUniprotId(e.target.value.toUpperCase())}
                placeholder="e.g. P00533"
                className="w-full text-sm px-3 py-2.5 rounded-lg border border-gray-200 font-mono focus:outline-none focus:ring-2 focus:ring-[#1e3a5f]/20 focus:border-[#1e3a5f] transition-colors"
              />
            </div>
          </div>

          {/* Description */}
          <div>
            <label className="block text-xs font-semibold text-gray-500 uppercase tracking-wider mb-1.5">
              Description
            </label>
            <textarea
              value={description}
              onChange={e => setDescription(e.target.value)}
              rows={3}
              placeholder="Brief description of the drug discovery goals and context..."
              className="w-full text-sm px-3 py-2.5 rounded-lg border border-gray-200 resize-none focus:outline-none focus:ring-2 focus:ring-[#1e3a5f]/20 focus:border-[#1e3a5f] transition-colors"
            />
          </div>

          {/* Actions */}
          <div className="flex items-center justify-end gap-3 pt-2">
            <button
              type="button"
              onClick={onClose}
              className="px-4 py-2 text-sm text-gray-500 hover:text-gray-700 rounded-lg hover:bg-gray-100 transition-colors"
            >
              Cancel
            </button>
            <button
              type="submit"
              disabled={saving || !name.trim() || !targetName.trim()}
              className="px-5 py-2.5 bg-[#22c55e] hover:bg-[#16a34a] text-white text-sm font-semibold rounded-lg transition-colors disabled:opacity-50 flex items-center gap-2"
            >
              {saving ? (
                <svg className="w-4 h-4 animate-spin" fill="none" viewBox="0 0 24 24">
                  <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
                  <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
                </svg>
              ) : (
                <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 4v16m8-8H4" />
                </svg>
              )}
              Create Project
            </button>
          </div>
        </form>
      </div>
    </div>
  )
}

// ---------------------------------------------------------------------------
// ProjectListPage
// ---------------------------------------------------------------------------

export default function ProjectListPage() {
  const { projects, createProject, loading } = useWorkspace()
  const navigate = useNavigate()

  const [searchRaw, setSearchRaw] = useState('')
  const [searchQuery, setSearchQuery] = useState('')
  const [showNewProject, setShowNewProject] = useState(false)

  // Debounce search
  useEffect(() => {
    const t = setTimeout(() => setSearchQuery(searchRaw), 300)
    return () => clearTimeout(t)
  }, [searchRaw])

  // Filtered projects
  const filtered = useMemo(() => {
    if (!searchQuery) return projects
    const q = searchQuery.toLowerCase()
    return projects.filter(p =>
      p.name?.toLowerCase().includes(q) ||
      p.description?.toLowerCase().includes(q) ||
      p.target_name?.toLowerCase().includes(q) ||
      p.target_input_value?.toLowerCase().includes(q)
    )
  }, [projects, searchQuery])

  // Stats
  const activeCount = projects.filter(p => {
    const allMols = (p.campaigns || []).flatMap(c => c.phases || []).reduce((s, ph) => s + (ph.stats?.total_molecules || 0), 0)
    return p.status === 'active' && allMols > 0
  }).length
  const emptyCount = projects.filter(p => {
    const allMols = (p.campaigns || []).flatMap(c => c.phases || []).reduce((s, ph) => s + (ph.stats?.total_molecules || 0), 0)
    return allMols === 0
  }).length

  const handleCreated = useCallback(async (data) => {
    const project = await createProject(data)
    setShowNewProject(false)
    navigate(`/project/${project.id}`)
    return project
  }, [createProject, navigate])

  return (
    <div className="space-y-6">
      {/* Page header */}
      <div className="flex items-start justify-between gap-4 flex-wrap">
        <div>
          <h1 className="text-2xl font-bold text-[#1e3a5f]">My Projects</h1>
          <p className="text-sm text-gray-400 mt-0.5">
            <span className="font-medium text-gray-600">{projects.length}</span> projects
            {activeCount > 0 && <> &middot; <span className="text-[#22c55e] font-medium">{activeCount} active</span></>}
            {emptyCount > 0 && <> &middot; <span className="text-gray-400">{emptyCount} empty</span></>}
          </p>
        </div>
        <button
          onClick={() => setShowNewProject(true)}
          className="flex items-center gap-2 px-4 py-2.5 bg-[#22c55e] hover:bg-[#16a34a] text-white text-sm font-semibold rounded-xl transition-colors shadow-sm"
        >
          <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 4v16m8-8H4" />
          </svg>
          New Project
        </button>
      </div>

      {/* Search */}
      <div className="relative">
        <svg className="absolute left-3.5 top-1/2 -translate-y-1/2 w-4 h-4 text-gray-400 pointer-events-none" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
        </svg>
        <input
          value={searchRaw}
          onChange={e => setSearchRaw(e.target.value)}
          placeholder="Search projects by name, target, or description..."
          className="w-full pl-10 pr-10 py-2.5 text-sm rounded-xl border border-gray-200 bg-white focus:outline-none focus:ring-2 focus:ring-[#1e3a5f]/20 focus:border-[#1e3a5f] transition-colors placeholder-gray-400"
        />
        {searchRaw && (
          <button
            onClick={() => { setSearchRaw(''); setSearchQuery('') }}
            className="absolute right-3.5 top-1/2 -translate-y-1/2 text-gray-400 hover:text-gray-600"
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
            </svg>
          </button>
        )}
      </div>

      {/* Project grid */}
      {loading ? (
        <div className="flex items-center justify-center py-20">
          <svg className="w-8 h-8 animate-spin text-[#1e3a5f]" fill="none" viewBox="0 0 24 24">
            <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
            <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
          </svg>
          <span className="ml-3 text-sm text-gray-500">Loading projects...</span>
        </div>
      ) : filtered.length === 0 ? (
        <div className="text-center py-20 bg-white rounded-xl border border-dashed border-gray-200">
          {searchQuery ? (
            <>
              <p className="text-gray-500 font-medium">No projects matching &ldquo;{searchQuery}&rdquo;</p>
              <button onClick={() => { setSearchRaw(''); setSearchQuery('') }} className="text-sm text-[#1e3a5f] hover:underline mt-2">
                Clear search
              </button>
            </>
          ) : (
            <>
              <svg className="w-16 h-16 mx-auto mb-4 text-gray-200" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.2}
                  d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
              </svg>
              <p className="text-lg font-semibold text-gray-500 mb-1">No projects yet</p>
              <p className="text-sm text-gray-400 mb-6 max-w-sm mx-auto">
                Create your first project to start your drug discovery campaign.
              </p>
              <button
                onClick={() => setShowNewProject(true)}
                className="inline-flex items-center gap-2 px-4 py-2.5 bg-[#22c55e] hover:bg-[#16a34a] text-white text-sm font-semibold rounded-xl transition-colors"
              >
                <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 4v16m8-8H4" />
                </svg>
                Create First Project
              </button>
            </>
          )}
        </div>
      ) : (
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
          {filtered.map(project => (
            <ProjectCard
              key={project.id}
              project={project}
              onClick={(id) => navigate(`/project/${id}`)}
              searchQuery={searchQuery}
            />
          ))}
        </div>
      )}

      {/* Recent Activity section */}
      {MOCK_ACTIVITY.length > 0 && (
        <div className="bg-white rounded-xl border border-gray-100 shadow-sm px-5 py-4">
          <h3 className="text-sm font-semibold text-gray-700 mb-4 flex items-center gap-2">
            <svg className="w-4 h-4 text-gray-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
                d="M12 8v4l3 3m6-3a9 9 0 11-18 0 9 9 0 0118 0z" />
            </svg>
            Recent Activity
          </h3>
          <ActivityTimeline activities={MOCK_ACTIVITY} limit={5} showProject />
        </div>
      )}

      {/* New Project modal */}
      {showNewProject && (
        <NewProjectModal
          onClose={() => setShowNewProject(false)}
          onCreated={handleCreated}
        />
      )}
    </div>
  )
}
