import React, { useState, useMemo, useCallback, useEffect } from 'react'
import { useNavigate } from 'react-router-dom'
import { useWorkspace } from '../contexts/WorkspaceContext.jsx'
import { PHASE_TYPES } from '../mock/data.js'
import BindXLogo from '../components/BindXLogo.jsx'

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

  // Phase brand colors per design system
  const phaseColors = [
    { border: 'border-bx-mint', text: 'text-bx-mint', bg: 'bg-bx-mint/10', fill: 'bg-bx-mint' },
    { border: 'border-bx-cyan', text: 'text-bx-cyan', bg: 'bg-bx-cyan/10', fill: 'bg-bx-cyan' },
    { border: 'border-bx-blue', text: 'text-bx-blue', bg: 'bg-bx-blue/10', fill: 'bg-bx-blue' },
  ]

  return (
    <div className="dark-inset inline-flex items-center px-3 py-2 gap-0">
      {types.map((type, i) => {
        const status = phaseStatus(type)
        const isLast = i === types.length - 1
        const pc = phaseColors[i]

        const isDone = status === 'done' || status === 'active'
        const nodeClass = isDone
          ? `${pc.border} ${pc.text} ${pc.bg}`
          : 'border-white/[.07] text-bx-dim bg-bx-bg'

        return (
          <React.Fragment key={type}>
            <div className="flex flex-col items-center gap-0.5">
              <div className={`w-6 h-6 rounded-full border-2 flex items-center justify-center text-[.55rem] font-bold ${nodeClass}`}>
                {status === 'done' ? (
                  <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={3} d="M5 13l4 4L19 7" />
                  </svg>
                ) : labels[i]}
              </div>
              <span className="text-bx-sub" style={{ fontSize: '.55rem', fontWeight: 500 }}>{labels[i]}</span>
            </div>
            {!isLast && (
              <div className={`h-0.5 w-5 mb-3.5 mx-0.5 rounded ${
                isDone ? 'bg-gradient-to-r from-bx-mint to-bx-cyan' : 'bg-white/[.07]'
              }`} />
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
    return <span className="badge badge-empty">Empty</span>
  }
  if (project.status === 'active') {
    return <span className="badge badge-active">Active</span>
  }
  return <span className="badge badge-empty">{project.status || 'Unknown'}</span>
}

// ---------------------------------------------------------------------------
// Project card
// ---------------------------------------------------------------------------

function ProjectCard({ project, onClick, onDelete, searchQuery }) {
  const allPhases = (project.campaigns || []).flatMap(c => c.phases || [])
  const totalMolecules = allPhases.reduce((s, p) => s + (p.stats?.total_molecules || 0), 0)
  const totalRuns = allPhases.reduce((s, p) => s + (p.stats?.runs_completed || 0) + (p.stats?.runs_running || 0), 0)

  return (
    <button
      onClick={() => onClick(project.id)}
      className="group bg-white rounded-card border border-bx-light-border shadow-card hover:shadow-card-hover hover:-translate-y-0.5 transition-all duration-200 text-left w-full overflow-hidden focus:outline-none focus:ring-2 focus:ring-bx-mint/20 focus:ring-offset-2"
    >
      {/* Top accent bar — visible only on active/hover */}
      <div className={`h-[3px] w-full transition-all ${
        totalMolecules > 0 && project.status === 'active'
          ? 'bg-gradient-to-r from-bx-mint to-bx-cyan'
          : 'bg-transparent group-hover:bg-bx-mint/30'
      }`} />

      <div className="p-5">
        {/* Header */}
        <div className="flex items-start justify-between gap-2 mb-3">
          <h3 className="font-bold text-bx-light-text text-base leading-tight group-hover:text-bx-mint-dim transition-colors font-serif">
            {highlight(project.name, searchQuery)}
          </h3>
          <StatusBadge project={project} />
        </div>

        {/* Target row — design system tags */}
        <div className="flex items-center gap-1.5 mb-3 flex-wrap">
          {project.target_name && (
            <span className="tag tag-gene-l">{highlight(project.target_name, searchQuery)}</span>
          )}
          {project.uniprot_id && (
            <span className="tag tag-uni-l">{project.uniprot_id}</span>
          )}
          {project.target_pdb_id && (
            <span className="tag tag-pdb-l">PDB: {project.target_pdb_id}</span>
          )}
        </div>

        {/* Phase progress */}
        <div className="mb-4">
          <PhaseProgress campaigns={project.campaigns} />
        </div>

        {/* Stats row */}
        <div className="grid grid-cols-3 gap-2 mb-3 pt-3 border-t border-bx-light-border-s">
          <div className="text-center">
            <p className="text-base font-extrabold text-bx-light-text">{totalMolecules}</p>
            <p className="text-[.55rem] text-bx-light-muted">molecules</p>
          </div>
          <div className="text-center">
            <p className="text-base font-extrabold text-bx-light-text2">{totalRuns}</p>
            <p className="text-[.55rem] text-bx-light-muted">runs</p>
          </div>
          <div className="text-center">
            <p className="text-sm font-mono font-medium text-bx-light-muted truncate">
              {project.updated_at ? relativeTime(project.updated_at) : '—'}
            </p>
            <p className="text-[.55rem] text-bx-light-muted">last update</p>
          </div>
        </div>

        {/* Description (if searched) */}
        {searchQuery && project.description?.toLowerCase().includes(searchQuery.toLowerCase()) && (
          <p className="text-sm text-gray-400 line-clamp-1 italic">
            {highlight(project.description, searchQuery)}
          </p>
        )}

        {/* Footer: campaign count + delete + chevron */}
        <div className="flex items-center justify-between mt-1">
          <span className="text-[.55rem] text-bx-light-muted font-mono">
            {project.campaigns.length} campaign{project.campaigns.length !== 1 ? 's' : ''}
          </span>
          <div className="flex items-center gap-2">
            <span
              role="button"
              tabIndex={0}
              onClick={(e) => {
                e.stopPropagation()
                if (window.confirm(`Delete "${project.name}"?`)) onDelete(project.id)
              }}
              onKeyDown={(e) => { if (e.key === 'Enter') e.currentTarget.click() }}
              className="w-6 h-6 rounded-md flex items-center justify-center text-bx-light-muted hover:text-bx-red hover:bg-bx-red/[.06] opacity-0 group-hover:opacity-100 transition-all cursor-pointer"
            >
              <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 7l-.867 12.142A2 2 0 0116.138 21H7.862a2 2 0 01-1.995-1.858L5 7m5 4v6m4-6v6m1-10V4a1 1 0 00-1-1h-4a1 1 0 00-1 1v3M4 7h16" />
              </svg>
            </span>
            <svg className="w-4 h-4 text-bx-light-muted group-hover:text-bx-mint transition-colors" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
            </svg>
          </div>
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
  const [description, setDescription] = useState('')
  const [saving, setSaving] = useState(false)

  useEffect(() => {
    const handler = (e) => { if (e.key === 'Escape') onClose() }
    document.addEventListener('keydown', handler)
    return () => document.removeEventListener('keydown', handler)
  }, [onClose])

  const handleSubmit = async (e) => {
    e.preventDefault()
    if (!name.trim()) return
    setSaving(true)
    try {
      await onCreated({
        name: name.trim(),
        description: description.trim() || undefined,
      })
    } catch (err) {
      console.error('[NewProjectModal] Create failed:', err)
    } finally {
      setSaving(false)
    }
  }

  return (
    <div
      className="fixed inset-0 bg-black/40 backdrop-blur-sm flex items-center justify-center z-50 px-4"
      onClick={(e) => { if (e.target === e.currentTarget) onClose() }}
    >
      <div className="bg-white rounded-[14px] shadow-lg w-full max-w-lg overflow-hidden">
        {/* Header */}
        <div className="bg-bx-surface px-6 py-4 flex items-center justify-between">
          <div>
            <h3 className="text-bx-text font-bold text-base">New Project</h3>
            <p className="text-bx-sub text-sm mt-0.5">Set up a new drug discovery project</p>
          </div>
          <button
            onClick={onClose}
            className="w-8 h-8 rounded-full bg-white/[.07] hover:bg-white/[.12] flex items-center justify-center text-white transition-colors"
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
            <label className="label">
              Project Name <span className="text-bx-red">*</span>
            </label>
            <input
              value={name}
              onChange={e => setName(e.target.value)}
              required
              autoFocus
              placeholder="e.g. EGFR Inhibitor Discovery"
              className="input-field text-sm py-2.5"
            />
          </div>

          {/* Description */}
          <div>
            <label className="label">
              Description
            </label>
            <textarea
              value={description}
              onChange={e => setDescription(e.target.value)}
              rows={3}
              placeholder="Brief description of the drug discovery goals and context..."
              className="input-field text-sm py-2.5 resize-none"
            />
          </div>

          {/* Actions */}
          <div className="flex items-center justify-end gap-3 pt-2">
            <button
              type="button"
              onClick={onClose}
              className="btn-ghost"
            >
              Cancel
            </button>
            <button
              type="submit"
              disabled={saving || !name.trim()}
              className="btn-primary text-sm py-2.5 px-5"
            >
              {saving ? (
                <BindXLogo variant="loading" size={16} />
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
  const { projects, createProject, deleteProject, loading } = useWorkspace()
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
    if (!project) return // auth error, toast already shown
    setShowNewProject(false)
    navigate(`/project/${project.id}`)
    return project
  }, [createProject, navigate])

  return (
    <div className="space-y-6">
      {/* Page header */}
      <div className="flex items-start justify-between gap-4 flex-wrap">
        <div>
          <h1 className="text-[1.6rem] font-bold text-bx-light-text font-serif">My Projects</h1>
          <p className="text-[.78rem] text-bx-light-muted mt-0.5">
            <span className="font-medium text-bx-light-text2">{projects.length}</span> projects
            {activeCount > 0 && <> &middot; <span className="text-bx-mint-dim font-medium">{activeCount} active</span></>}
            {emptyCount > 0 && <> &middot; <span className="text-bx-light-muted">{emptyCount} empty</span></>}
          </p>
        </div>
        <button
          onClick={() => setShowNewProject(true)}
          className="btn-primary text-sm py-2.5 px-4"
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
          className="input-field pl-10 pr-10 py-2.5"
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
          <BindXLogo variant="loading" size={64} label="Loading projects..." />
        </div>
      ) : filtered.length === 0 ? (
        <div className="text-center py-20 card border-dashed">
          {searchQuery ? (
            <>
              <p className="text-bx-light-text2 font-medium">No projects matching &ldquo;{searchQuery}&rdquo;</p>
              <button onClick={() => { setSearchRaw(''); setSearchQuery('') }} className="text-sm text-bx-mint-dim hover:underline mt-2">
                Clear search
              </button>
            </>
          ) : (
            <>
              <div className="mb-4">
                <BindXLogo variant="idle" size={64} />
              </div>
              <p className="text-lg font-semibold text-bx-light-text2 mb-1">No projects yet</p>
              <p className="text-[.78rem] text-bx-light-muted mb-6 max-w-sm mx-auto">
                Create your first project to start your drug discovery campaign.
              </p>
              <button
                onClick={() => setShowNewProject(true)}
                className="btn-primary text-sm py-2.5 px-4"
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
              onDelete={(id) => deleteProject(id)}
              searchQuery={searchQuery}
            />
          ))}
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
