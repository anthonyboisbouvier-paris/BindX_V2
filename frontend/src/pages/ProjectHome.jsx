import React, { useEffect, useState } from 'react'
import { useParams, useNavigate, Link } from 'react-router-dom'
import { useWorkspace } from '../contexts/WorkspaceContext.jsx'
import { useToast } from '../contexts/ToastContext.jsx'
import {
  PHASE_TYPES,
  MOCK_AGENT_INSIGHTS,
  MOCK_TARGET_ASSESSMENT,
  MOCK_ACTIVITY,
} from '../mock/data.js'

import TargetInfoCard from '../components/TargetInfoCard.jsx'
import CampaignAgentPanel from '../components/CampaignAgentPanel.jsx'
import ScoringWeightsEditor from '../components/ScoringWeightsEditor.jsx'
import PhaseCreator from '../components/PhaseCreator.jsx'
import ActivityTimeline from '../components/ActivityTimeline.jsx'

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function fmtDate(iso) {
  if (!iso) return null
  return new Date(iso).toLocaleDateString('en-US', { month: 'short', day: 'numeric', year: 'numeric' })
}

// ---------------------------------------------------------------------------
// Phase color palette
// ---------------------------------------------------------------------------

const PHASE_COLORS = {
  hit_discovery: {
    border: 'border-blue-200',
    leftBar: 'bg-blue-500',
    header: 'bg-blue-50',
    accent: 'text-blue-700',
    badge: 'bg-blue-100 text-blue-700',
    btn: 'bg-blue-600 hover:bg-blue-700 text-white',
    dotFilled: 'bg-blue-500',
    frozenBg: 'bg-blue-50/50',
  },
  hit_to_lead: {
    border: 'border-purple-200',
    leftBar: 'bg-purple-500',
    header: 'bg-purple-50',
    accent: 'text-purple-700',
    badge: 'bg-purple-100 text-purple-700',
    btn: 'bg-purple-600 hover:bg-purple-700 text-white',
    dotFilled: 'bg-purple-500',
    frozenBg: 'bg-purple-50/50',
  },
  lead_optimization: {
    border: 'border-green-200',
    leftBar: 'bg-[#22c55e]',
    header: 'bg-green-50',
    accent: 'text-green-700',
    badge: 'bg-green-100 text-green-700',
    btn: 'bg-[#22c55e] hover:bg-[#16a34a] text-white',
    dotFilled: 'bg-[#22c55e]',
    frozenBg: 'bg-green-50/50',
  },
}

// Funnel widths per phase (visually narrowing)
const FUNNEL_WIDTHS = ['w-full', 'max-w-[90%]', 'max-w-[78%]']

// ---------------------------------------------------------------------------
// Phase status badge (inline)
// ---------------------------------------------------------------------------

function PhaseBadge({ status, notCreated }) {
  if (notCreated) {
    return (
      <span className="inline-flex items-center gap-1 px-2 py-0.5 rounded-full text-xs bg-gray-100 text-gray-400">
        <span className="w-1.5 h-1.5 rounded-full bg-gray-300" />
        Not started
      </span>
    )
  }
  if (status === 'frozen') {
    return (
      <span className="inline-flex items-center gap-1.5 px-2 py-0.5 rounded-full text-xs bg-blue-50 text-blue-600 font-medium">
        <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
            d="M12 15v2m-6 4h12a2 2 0 002-2v-6a2 2 0 00-2-2H6a2 2 0 00-2 2v6a2 2 0 002 2zm10-10V7a4 4 0 00-8 0v4h8z" />
        </svg>
        Frozen
      </span>
    )
  }
  return (
    <span className="inline-flex items-center gap-1.5 px-2 py-0.5 rounded-full text-xs bg-green-50 text-green-700 font-medium">
      <span className="w-1.5 h-1.5 rounded-full bg-[#22c55e]" />
      Active
    </span>
  )
}

// ---------------------------------------------------------------------------
// Funnel Phase Card
// ---------------------------------------------------------------------------

function FunnelPhaseCard({ phase, phaseIndex, projectId, navigate, onCreatePhase }) {
  const notCreated = !phase.created_at
  const colors = PHASE_COLORS[phase.type] || PHASE_COLORS.hit_discovery
  const phaseType = PHASE_TYPES[phase.type] || { label: phase.label, short: '?', color: 'blue' }
  const stats = phase.stats || {}
  const isFrozen = phase.status === 'frozen'
  const width = FUNNEL_WIDTHS[phaseIndex] || 'w-full'

  if (notCreated) {
    return (
      <div className={`mx-auto ${width}`}>
        <div
          className={`
            border-2 border-dashed border-gray-200 rounded-xl p-5
            flex items-center justify-between gap-4
            bg-gray-50/60
          `}
        >
          <div className="flex items-center gap-3">
            <div className="w-8 h-8 rounded-lg bg-gray-100 flex items-center justify-center">
              <span className="text-sm font-bold text-gray-400">{phaseType.short}</span>
            </div>
            <div>
              <p className="text-sm font-semibold text-gray-400">{phase.label}</p>
              <p className="text-xs text-gray-400">{phaseType.label}</p>
            </div>
          </div>
          <button
            className={`flex items-center gap-1.5 px-3 py-1.5 text-xs font-semibold rounded-lg transition-colors ${colors.btn}`}
            onClick={() => onCreatePhase && onCreatePhase(phase)}
          >
            <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 4v16m8-8H4" />
            </svg>
            Create Phase
          </button>
        </div>
      </div>
    )
  }

  return (
    <div className={`mx-auto ${width}`}>
      <div
        className={`
          rounded-xl border shadow-sm overflow-hidden
          transition-all duration-200 hover:shadow-md
          ${colors.border}
          ${isFrozen ? colors.frozenBg : 'bg-white'}
        `}
      >
        {/* Left color bar */}
        <div className="flex">
          <div className={`w-1 shrink-0 ${colors.leftBar}`} />
          <div className="flex-1">
            {/* Card header */}
            <div className={`px-4 py-3 ${colors.header} flex items-center justify-between`}>
              <div className="flex items-center gap-2.5">
                <span className={`w-7 h-7 rounded-lg flex items-center justify-center text-sm font-bold ${colors.badge}`}>
                  {phaseType.short}
                </span>
                <div>
                  <p className={`text-sm font-semibold ${colors.accent}`}>{phase.label}</p>
                  <p className="text-xs text-gray-500">{phaseType.label}</p>
                </div>
                {isFrozen && (
                  <svg className="w-3.5 h-3.5 text-blue-400 ml-1" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                      d="M12 15v2m-6 4h12a2 2 0 002-2v-6a2 2 0 00-2-2H6a2 2 0 00-2 2v6a2 2 0 002 2zm10-10V7a4 4 0 00-8 0v4h8z" />
                  </svg>
                )}
              </div>
              <div className="flex items-center gap-2">
                <PhaseBadge status={phase.status} notCreated={false} />
                <button
                  className={`flex items-center gap-1 px-2.5 py-1 text-xs font-medium rounded-lg transition-colors ${colors.btn}`}
                  onClick={() => navigate(`/project/${projectId}/phase/${phase.id}`)}
                >
                  Open
                  <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
                  </svg>
                </button>
              </div>
            </div>

            {/* Card body */}
            <div className="px-4 py-3 flex items-center gap-6 flex-wrap">
              {/* Main stats */}
              <div className="flex items-center gap-5">
                <div className="text-center">
                  <p className="text-xl font-bold text-[#1e3a5f]">{stats.total_molecules ?? 0}</p>
                  <p className="text-xs text-gray-400">molecules</p>
                </div>
                <div className="text-gray-200 text-lg font-light">|</div>
                <div className="text-center">
                  <p className="text-xl font-bold text-gray-700">
                    {stats.bookmarked ?? 0}
                    <span className="text-sm text-gray-400 font-normal ml-0.5">bookmarked</span>
                  </p>
                </div>
                <div className="text-gray-200 text-lg font-light">|</div>
                <div className="text-center">
                  <p className="text-sm font-semibold text-gray-600">
                    {stats.runs_completed ?? 0}
                    <span className="text-xs text-gray-400 font-normal ml-0.5">
                      {(stats.runs_running ?? 0) > 0 ? `+ ${stats.runs_running} running` : 'runs'}
                    </span>
                  </p>
                </div>
              </div>

              {/* Running indicator */}
              {(stats.runs_running ?? 0) > 0 && (
                <div className="flex items-center gap-1.5 text-xs text-amber-600 bg-amber-50 px-2.5 py-1 rounded-lg ml-auto">
                  <svg className="w-3 h-3 animate-spin" fill="none" viewBox="0 0 24 24">
                    <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
                    <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
                  </svg>
                  {stats.runs_running} run{stats.runs_running !== 1 ? 's' : ''} in progress
                </div>
              )}
            </div>
          </div>
        </div>
      </div>
    </div>
  )
}

// ---------------------------------------------------------------------------
// Funnel connector arrows
// ---------------------------------------------------------------------------

function FunnelConnector({ prevPhase }) {
  if (!prevPhase?.created_at) return null
  const bookmarked = prevPhase.stats?.bookmarked ?? 0
  return (
    <div className="flex items-center justify-center gap-2 py-1">
      <div className="flex flex-col items-center gap-0.5">
        <div className="w-px h-3 bg-gray-200" />
        <div className="flex items-center gap-2">
          <div className="w-2 h-px bg-gray-200" />
          <span className="text-xs text-gray-400 bg-white px-2 py-0.5 rounded-full border border-gray-100">
            {bookmarked} promoted
          </span>
          <div className="w-2 h-px bg-gray-200" />
        </div>
        <svg className="w-3 h-3 text-gray-300" fill="currentColor" viewBox="0 0 12 12">
          <path d="M6 9L1 3h10z" />
        </svg>
      </div>
    </div>
  )
}

// ---------------------------------------------------------------------------
// Campaign Settings (collapsible, wraps ScoringWeightsEditor)
// ---------------------------------------------------------------------------

function CampaignSettings({ campaign }) {
  const [open, setOpen] = useState(false)

  if (!campaign) return null
  const { scoring_weights, docking_defaults, rules } = campaign

  return (
    <div className="bg-white rounded-xl border border-gray-100 shadow-sm overflow-hidden">
      <button
        onClick={() => setOpen(o => !o)}
        className="w-full flex items-center justify-between px-5 py-4 text-left hover:bg-gray-50 transition-colors"
      >
        <div className="flex items-center gap-2">
          <svg className="w-4 h-4 text-gray-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
              d="M10.325 4.317c.426-1.756 2.924-1.756 3.35 0a1.724 1.724 0 002.573 1.066c1.543-.94 3.31.826 2.37 2.37a1.724 1.724 0 001.065 2.572c1.756.426 1.756 2.924 0 3.35a1.724 1.724 0 00-1.066 2.573c.94 1.543-.826 3.31-2.37 2.37a1.724 1.724 0 00-2.572 1.065c-.426 1.756-2.924 1.756-3.35 0a1.724 1.724 0 00-2.573-1.066c-1.543.94-3.31-.826-2.37-2.37a1.724 1.724 0 00-1.065-2.572c-1.756-.426-1.756-2.924 0-3.35a1.724 1.724 0 001.066-2.573c-.94-1.543.826-3.31 2.37-2.37.996.608 2.296.07 2.572-1.065z" />
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 12a3 3 0 11-6 0 3 3 0 016 0z" />
          </svg>
          <span className="text-sm font-semibold text-gray-700">Campaign Settings</span>
          <span className="text-xs text-gray-400">Scoring Weights, Docking Defaults, Filter Rules</span>
        </div>
        <svg
          className={`w-4 h-4 text-gray-400 transition-transform duration-200 ${open ? 'rotate-180' : ''}`}
          fill="none" stroke="currentColor" viewBox="0 0 24 24"
        >
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
        </svg>
      </button>

      {open && (
        <div className="border-t border-gray-50 px-5 py-4 space-y-5">
          {/* Scoring weights interactive editor */}
          {scoring_weights && (
            <div>
              <ScoringWeightsEditor
                weights={scoring_weights}
                onChange={() => {}}
              />
            </div>
          )}

          {/* Docking defaults */}
          {docking_defaults && (
            <div>
              <p className="text-xs font-semibold text-gray-500 uppercase tracking-wider mb-2">Docking Defaults</p>
              <div className="flex flex-wrap gap-2">
                {Object.entries(docking_defaults).map(([k, v]) => (
                  <span key={k} className="px-2.5 py-1 bg-blue-50 text-blue-700 text-xs rounded-md font-mono">
                    {k}: {String(v)}
                  </span>
                ))}
              </div>
            </div>
          )}

          {/* Filter rules */}
          {rules && (
            <div>
              <p className="text-xs font-semibold text-gray-500 uppercase tracking-wider mb-2">Filter Rules</p>
              <div className="flex flex-wrap gap-2">
                {rules.lipinski && (
                  <span className="px-2.5 py-1 bg-green-50 text-green-700 text-xs rounded-md flex items-center gap-1">
                    <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 13l4 4L19 7" />
                    </svg>
                    Lipinski Ro5
                  </span>
                )}
                {rules.pains && (
                  <span className="px-2.5 py-1 bg-green-50 text-green-700 text-xs rounded-md flex items-center gap-1">
                    <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 13l4 4L19 7" />
                    </svg>
                    PAINS Filter
                  </span>
                )}
                {rules.max_MW && (
                  <span className="px-2.5 py-1 bg-blue-50 text-blue-700 text-xs rounded-md font-mono">
                    MW &lt; {rules.max_MW} Da
                  </span>
                )}
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  )
}

// ---------------------------------------------------------------------------
// Edit project inline form
// ---------------------------------------------------------------------------

function EditProjectPanel({ project, onSave, onCancel }) {
  const [name, setName] = useState(project.name || '')
  const [description, setDescription] = useState(project.description || '')

  return (
    <div className="bg-blue-50 border border-blue-100 rounded-xl px-5 py-4 space-y-3">
      <div className="grid grid-cols-1 sm:grid-cols-2 gap-3">
        <div>
          <label className="block text-xs font-semibold text-gray-500 mb-1">Project Name</label>
          <input
            value={name}
            onChange={e => setName(e.target.value)}
            className="w-full text-sm px-3 py-2 rounded-lg border border-gray-200 focus:outline-none focus:ring-2 focus:ring-[#1e3a5f]/20 focus:border-[#1e3a5f]"
          />
        </div>
        <div>
          <label className="block text-xs font-semibold text-gray-500 mb-1">Description</label>
          <textarea
            value={description}
            onChange={e => setDescription(e.target.value)}
            rows={2}
            className="w-full text-sm px-3 py-2 rounded-lg border border-gray-200 resize-none focus:outline-none focus:ring-2 focus:ring-[#1e3a5f]/20 focus:border-[#1e3a5f]"
          />
        </div>
      </div>
      <div className="flex gap-2">
        <button
          onClick={() => onSave({ name, description })}
          className="px-4 py-1.5 bg-[#1e3a5f] text-white text-xs font-semibold rounded-lg hover:bg-[#1e4a7f] transition-colors"
        >
          Save
        </button>
        <button
          onClick={onCancel}
          className="px-4 py-1.5 bg-gray-100 text-gray-600 text-xs font-medium rounded-lg hover:bg-gray-200 transition-colors"
        >
          Cancel
        </button>
      </div>
    </div>
  )
}

// ---------------------------------------------------------------------------
// ProjectHome
// ---------------------------------------------------------------------------

export default function ProjectHome() {
  const { projectId } = useParams()
  const navigate = useNavigate()
  const { projects, selectProject, updateProject, deleteProject, createCampaign, createPhase } = useWorkspace()
  const { addToast } = useToast()

  const [editing, setEditing] = useState(false)
  const [showPhaseCreator, setShowPhaseCreator] = useState(false)
  const [creatingForPhase, setCreatingForPhase] = useState(null)
  const [campaignCreating, setCampaignCreating] = useState(false)

  useEffect(() => {
    if (projectId) selectProject(projectId)
  }, [projectId, selectProject])

  const project = projects.find(p => p.id === projectId)

  if (!project) {
    return (
      <div className="text-center py-20">
        <p className="text-gray-400 text-sm">Project not found.</p>
        <Link to="/" className="text-[#1e3a5f] text-sm hover:underline mt-2 inline-block">
          Back to projects
        </Link>
      </div>
    )
  }

  const campaign = project.campaigns?.[0] || null
  const phases = campaign?.phases || []
  const assessment = MOCK_TARGET_ASSESSMENT[project.id] || null

  // Filter activity to this project
  const projectActivity = MOCK_ACTIVITY.filter(a => a.project_id === project.id)

  return (
    <div className="space-y-5">
      {/* Back navigation */}
      <div>
        <Link
          to="/"
          className="inline-flex items-center gap-1.5 text-xs text-gray-400 hover:text-[#1e3a5f] transition-colors"
        >
          <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 19l-7-7 7-7" />
          </svg>
          All Projects
        </Link>
      </div>

      {/* ---- Project header ---- */}
      <div className="bg-white rounded-xl border border-gray-100 shadow-sm px-6 py-5">
        <div className="flex items-start justify-between gap-4">
          <div className="flex-1 min-w-0">
            <div className="flex items-center gap-3 flex-wrap mb-1.5">
              <h1 className="text-2xl font-bold text-[#1e3a5f] leading-tight">{project.name}</h1>
              <span className={`inline-flex items-center gap-1.5 px-2.5 py-0.5 rounded-full text-xs font-medium ${
                project.status === 'active'
                  ? 'bg-green-50 text-green-700'
                  : 'bg-gray-100 text-gray-500'
              }`}>
                <span className={`w-1.5 h-1.5 rounded-full ${project.status === 'active' ? 'bg-[#22c55e]' : 'bg-gray-400'}`} />
                {project.status === 'active' ? 'Active' : project.status || 'Active'}
              </span>
            </div>

            {project.description && (
              <p className="text-gray-500 text-sm leading-relaxed mb-3 max-w-2xl">
                {project.description}
              </p>
            )}

            {/* Target metadata row */}
            <div className="flex items-center flex-wrap gap-3 text-sm">
              {project.target_name && (
                <span className="inline-flex items-center gap-1.5 text-gray-700 font-medium">
                  <svg className="w-4 h-4 text-gray-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <circle cx="12" cy="12" r="9" strokeWidth={1.8} />
                    <circle cx="12" cy="12" r="5" strokeWidth={1.8} />
                    <circle cx="12" cy="12" r="1.5" fill="currentColor" strokeWidth={0} />
                  </svg>
                  {project.target_name}
                </span>
              )}
              {project.uniprot_id && (
                <span className="px-2 py-0.5 bg-blue-50 text-blue-700 rounded font-mono text-xs font-semibold">
                  {project.uniprot_id}
                </span>
              )}
              {project.target_pdb_id && (
                <span className="text-gray-400 text-xs">
                  PDB: <span className="font-mono text-gray-600 font-semibold">{project.target_pdb_id}</span>
                </span>
              )}
              {project.updated_at && (
                <span className="text-gray-400 text-xs ml-auto">
                  Updated {fmtDate(project.updated_at)}
                </span>
              )}
            </div>
          </div>

          {/* Action buttons */}
          <div className="flex items-center gap-2 shrink-0">
            <button
              className="px-3 py-1.5 text-xs text-gray-500 border border-gray-200 rounded-lg hover:bg-gray-50 transition-colors flex items-center gap-1.5"
              onClick={() => setEditing(e => !e)}
            >
              <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                  d="M11 5H6a2 2 0 00-2 2v11a2 2 0 002 2h11a2 2 0 002-2v-5m-1.414-9.414a2 2 0 112.828 2.828L11.828 15H9v-2.828l8.586-8.586z" />
              </svg>
              Edit
            </button>
            <button
              className="px-3 py-1.5 text-xs text-gray-400 border border-gray-200 rounded-lg hover:bg-gray-50 transition-colors"
              onClick={async () => {
                if (!confirm('Archive this project? It will be hidden from the project list.')) return
                try {
                  await updateProject(projectId, { status: 'archived' })
                  addToast('Project archived', 'info')
                  navigate('/')
                } catch (err) {
                  addToast(err.userMessage || 'Failed to archive project', 'error')
                }
              }}
            >
              Archive
            </button>
          </div>
        </div>

        {/* Inline edit panel */}
        {editing && (
          <div className="mt-4">
            <EditProjectPanel
              project={project}
              onSave={async (data) => {
                try {
                  await updateProject(projectId, data)
                } catch (err) {
                  console.error('[ProjectHome] Update failed:', err)
                }
                setEditing(false)
              }}
              onCancel={() => setEditing(false)}
            />
          </div>
        )}
      </div>

      {/* ---- Target Info Card ---- */}
      {assessment && (
        <TargetInfoCard assessment={assessment} project={project} />
      )}

      {/* ---- Campaign section ---- */}
      {campaign ? (
        <>
          {/* Campaign header */}
          <div className="flex items-center justify-between">
            <div>
              <p className="text-xs text-gray-400 uppercase tracking-wider font-semibold mb-0.5">Campaign</p>
              <h2 className="text-lg font-bold text-gray-800">{campaign.name}</h2>
            </div>
            <span className="text-xs text-gray-400">
              {phases.filter(p => p.created_at).length}/{phases.length} phases created
            </span>
          </div>

          {/* Drug Discovery Funnel */}
          <div className="bg-white rounded-xl border border-gray-100 shadow-sm px-6 py-6">
            <div className="flex items-center gap-2 mb-5">
              <svg className="w-4 h-4 text-gray-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
                  d="M3 4a1 1 0 011-1h16a1 1 0 011 1v2.586a1 1 0 01-.293.707l-6.414 6.414a1 1 0 00-.293.707V17l-4 4v-6.586a1 1 0 00-.293-.707L3.293 7.293A1 1 0 013 6.586V4z" />
              </svg>
              <h3 className="text-sm font-semibold text-gray-700">Drug Discovery Funnel</h3>
            </div>

            <div className="space-y-2">
              {phases.map((phase, i) => (
                <React.Fragment key={phase.id}>
                  {i > 0 && (
                    <FunnelConnector prevPhase={phases[i - 1]} />
                  )}
                  <FunnelPhaseCard
                    phase={phase}
                    phaseIndex={i}
                    projectId={projectId}
                    navigate={navigate}
                    onCreatePhase={(ph) => {
                      setCreatingForPhase(ph)
                      setShowPhaseCreator(true)
                    }}
                  />
                </React.Fragment>
              ))}
            </div>
          </div>

          {/* AI Agent Insights */}
          <CampaignAgentPanel
            insights={MOCK_AGENT_INSIGHTS}
            onAskAgent={(q) => addToast('AI Agent is not yet connected to the backend', 'info')}
          />

          {/* Campaign settings */}
          <CampaignSettings campaign={campaign} />
        </>
      ) : (
        /* No campaign yet */
        <div className="bg-white rounded-xl border border-dashed border-gray-200 px-6 py-16 text-center">
          <svg className="w-12 h-12 text-gray-200 mx-auto mb-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
              d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
          </svg>
          <p className="text-gray-500 font-medium mb-1">No campaigns yet</p>
          <p className="text-sm text-gray-400 mb-5 max-w-sm mx-auto">
            Create a campaign to define a binding pocket and start screening molecules.
          </p>
          <button
            className={`inline-flex items-center gap-2 px-4 py-2.5 text-sm font-semibold rounded-xl transition-colors ${
              campaignCreating
                ? 'bg-gray-300 text-gray-500 cursor-not-allowed'
                : 'bg-[#1e3a5f] text-white hover:bg-[#1e4a7f]'
            }`}
            disabled={campaignCreating}
            onClick={async () => {
              try {
                setCampaignCreating(true)
                await createCampaign(projectId, {
                  name: `Campaign — ${project.target_name || project.name}`,
                })
                addToast('Campaign created', 'success')
              } catch (err) {
                addToast(err.userMessage || 'Failed to create campaign', 'error')
              } finally {
                setCampaignCreating(false)
              }
            }}
          >
            {campaignCreating ? (
              <>
                <svg className="w-4 h-4 animate-spin" fill="none" viewBox="0 0 24 24">
                  <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
                  <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
                </svg>
                Creating...
              </>
            ) : (
              <>
                <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 4v16m8-8H4" />
                </svg>
                New Campaign
              </>
            )}
          </button>
        </div>
      )}

      {/* Recent Activity */}
      {projectActivity.length > 0 && (
        <div className="bg-white rounded-xl border border-gray-100 shadow-sm px-5 py-4">
          <h3 className="text-sm font-semibold text-gray-700 mb-4 flex items-center gap-2">
            <svg className="w-4 h-4 text-gray-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
                d="M12 8v4l3 3m6-3a9 9 0 11-18 0 9 9 0 0118 0z" />
            </svg>
            Recent Activity
          </h3>
          <ActivityTimeline activities={projectActivity} showProject={false} />
        </div>
      )}

      {/* Phase creator modal */}
      {showPhaseCreator && (
        <PhaseCreator
          campaignId={campaign?.id}
          existingPhases={phases}
          onCreatePhase={async (newPhaseConfig) => {
            try {
              await createPhase(newPhaseConfig.campaign_id, { type: newPhaseConfig.type })
              addToast('Phase created', 'success')
            } catch (err) {
              addToast(err.userMessage || 'Failed to create phase', 'error')
            }
            setShowPhaseCreator(false)
            setCreatingForPhase(null)
          }}
          onClose={() => {
            setShowPhaseCreator(false)
            setCreatingForPhase(null)
          }}
        />
      )}
    </div>
  )
}
