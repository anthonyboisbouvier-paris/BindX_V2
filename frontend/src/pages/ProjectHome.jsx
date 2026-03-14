import React, { useEffect, useState, useMemo } from 'react'
import { useParams, useNavigate, Link } from 'react-router-dom'
import { useWorkspace } from '../contexts/WorkspaceContext.jsx'
import { useToast } from '../contexts/ToastContext.jsx'
import { PHASE_TYPES } from '../lib/columns.js'

import PhaseCreator from '../components/PhaseCreator.jsx'
import ProteinViewer from '../components/ProteinViewer.jsx'
import BindXLogo from '../components/BindXLogo.jsx'
import CampaignAgentPanel from '../components/CampaignAgentPanel.jsx'

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
    border: 'border-emerald-200',
    leftBar: 'bg-bx-mint',
    header: 'bg-emerald-50',
    accent: 'text-emerald-700',
    badge: 'badge-phase-a',
    btn: 'btn-primary',
    dotFilled: 'bg-bx-mint',
    frozenBg: 'bg-emerald-50/50',
  },
  hit_to_lead: {
    border: 'border-cyan-200',
    leftBar: 'bg-bx-cyan',
    header: 'bg-cyan-50',
    accent: 'text-cyan-700',
    badge: 'badge-phase-b',
    btn: 'bg-bx-cyan hover:brightness-110 text-white text-sm font-bold rounded-btn px-3 py-1.5',
    dotFilled: 'bg-bx-cyan',
    frozenBg: 'bg-cyan-50/50',
  },
  lead_optimization: {
    border: 'border-blue-200',
    leftBar: 'bg-bx-blue',
    header: 'bg-blue-50',
    accent: 'text-blue-700',
    badge: 'badge-phase-c',
    btn: 'bg-bx-blue hover:brightness-110 text-white text-sm font-bold rounded-btn px-3 py-1.5',
    dotFilled: 'bg-bx-blue',
    frozenBg: 'bg-blue-50/50',
  },
}

// Funnel widths per phase (visually narrowing)
const FUNNEL_WIDTHS = ['w-full', 'max-w-[90%]', 'max-w-[78%]']

// ---------------------------------------------------------------------------
// Phase status badge (inline)
// ---------------------------------------------------------------------------

function PhaseBadge({ status, notCreated }) {
  if (notCreated) {
    return <span className="badge badge-empty">Not started</span>
  }
  if (status === 'frozen') {
    return (
      <span className="badge badge-phase-c inline-flex items-center gap-1">
        <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
            d="M12 15v2m-6 4h12a2 2 0 002-2v-6a2 2 0 00-2-2H6a2 2 0 00-2 2v6a2 2 0 002 2zm10-10V7a4 4 0 00-8 0v4h8z" />
        </svg>
        Frozen
      </span>
    )
  }
  return <span className="badge badge-active">Active</span>
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
            border-2 border-dashed border-bx-light-border rounded-card p-5
            flex items-center justify-between gap-4
            bg-bx-light-warm
          `}
        >
          <div className="flex items-center gap-3">
            <div className="w-8 h-8 rounded-lg bg-gray-100 flex items-center justify-center">
              <span className="text-sm font-bold text-gray-400">{phaseType.short}</span>
            </div>
            <div>
              <p className="text-sm font-semibold text-gray-400">{phase.label}</p>
              <p className="text-sm text-gray-400">{phaseType.label}</p>
            </div>
          </div>
          <button
            className={`flex items-center gap-1.5 px-3 py-1.5 text-sm font-semibold rounded-lg transition-colors ${colors.btn}`}
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
          rounded-card border shadow-card overflow-hidden
          transition-all duration-200 hover:shadow-card-hover
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
                  <p className="text-sm text-gray-500">{phaseType.label}</p>
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
                  className={`flex items-center gap-1 px-2.5 py-1 text-sm font-medium rounded-lg transition-colors ${colors.btn}`}
                  onClick={() => navigate(`/project/${projectId}/phase/${phase.id}`)}
                >
                  Open
                  <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
                  </svg>
                </button>
              </div>
            </div>

            {/* Card body — stats in dark inset */}
            <div className="px-4 py-3">
              <div className="dark-inset flex items-center gap-4 flex-wrap">
                <div className="text-center px-2">
                  <p className="stat-value text-bx-mint">{stats.total_molecules ?? 0}</p>
                  <p className="stat-label">molecules</p>
                </div>
                <div className="text-white/[.07] text-lg font-light">|</div>
                <div className="text-center px-2">
                  <p className="stat-value text-bx-cyan">{stats.bookmarked ?? 0}</p>
                  <p className="stat-label">bookmarked</p>
                </div>
                <div className="text-white/[.07] text-lg font-light">|</div>
                <div className="text-center px-2">
                  <p className="stat-value text-bx-sub">{stats.runs_completed ?? 0}</p>
                  <p className="stat-label">runs</p>
                </div>

                {/* Last run type */}
                {stats.last_run_type && (
                  <>
                    <div className="text-white/[.07] text-lg font-light">|</div>
                    <div className="text-center px-2">
                      <p className="text-sm font-semibold text-bx-sub capitalize">{stats.last_run_type}</p>
                      <p className="stat-label">last run</p>
                    </div>
                  </>
                )}

                {/* Running indicator */}
                {(stats.runs_running ?? 0) > 0 && (
                  <div className="flex items-center gap-1.5 badge badge-running ml-auto">
                    <BindXLogo variant="loading" size={20} />
                    {stats.runs_running} run{stats.runs_running !== 1 ? 's' : ''} in progress
                  </div>
                )}
              </div>
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
          <span className="badge badge-empty">
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
// Campaign Strategy Brief — context for AI agents
// ---------------------------------------------------------------------------

function CampaignStrategyBrief({ campaign, onSave }) {
  const { addToast } = useToast()
  const [editing, setEditing] = useState(false)
  const [notes, setNotes] = useState(campaign?.strategy_notes || '')
  const [saving, setSaving] = useState(false)

  const hasContent = !!campaign?.strategy_notes?.trim()

  const handleSave = async () => {
    if (!campaign?.id) return
    setSaving(true)
    try {
      await onSave(campaign.id, { strategy_notes: notes })
      addToast('Strategy brief saved', 'success')
      setEditing(false)
    } catch (err) {
      addToast(err.userMessage || 'Failed to save', 'error')
    } finally {
      setSaving(false)
    }
  }

  if (!editing && !hasContent) {
    // Empty state — invitation to write
    return (
      <div className="bg-gradient-to-br from-violet-50 to-indigo-50 rounded-xl border border-violet-200 px-5 py-5">
        <div className="flex items-start gap-3">
          <div className="w-9 h-9 rounded-lg bg-violet-100 flex items-center justify-center shrink-0 mt-0.5">
            <svg className="w-5 h-5 text-violet-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
                d="M9.813 15.904L9 18.75l-.813-2.846a4.5 4.5 0 00-3.09-3.09L2.25 12l2.846-.813a4.5 4.5 0 003.09-3.09L9 5.25l.813 2.846a4.5 4.5 0 003.09 3.09L15.75 12l-2.846.813a4.5 4.5 0 00-3.09 3.09zM18.259 8.715L18 9.75l-.259-1.035a3.375 3.375 0 00-2.455-2.456L14.25 6l1.036-.259a3.375 3.375 0 002.455-2.456L18 2.25l.259 1.035a3.375 3.375 0 002.455 2.456L21.75 6l-1.036.259a3.375 3.375 0 00-2.455 2.456zM16.894 20.567L16.5 21.75l-.394-1.183a2.25 2.25 0 00-1.423-1.423L13.5 18.75l1.183-.394a2.25 2.25 0 001.423-1.423l.394-1.183.394 1.183a2.25 2.25 0 001.423 1.423l1.183.394-1.183.394a2.25 2.25 0 00-1.423 1.423z" />
            </svg>
          </div>
          <div className="flex-1">
            <h3 className="text-sm font-semibold text-violet-900 mb-1">AI Strategy Brief</h3>
            <p className="text-sm text-violet-700 leading-relaxed mb-3">
              Guide the AI agents by sharing your expertise on this target. What do you know about the binding site?
              What properties matter most? Any known SAR, selectivity concerns, or chemical series to explore?
            </p>
            <button
              onClick={() => setEditing(true)}
              className="inline-flex items-center gap-1.5 px-3 py-1.5 text-sm font-semibold bg-violet-600 text-white rounded-lg hover:bg-violet-700 transition-colors"
            >
              <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 4v16m8-8H4" />
              </svg>
              Write Strategy Brief
            </button>
          </div>
        </div>
      </div>
    )
  }

  if (editing) {
    return (
      <div className="bg-white rounded-xl border border-violet-200 shadow-sm overflow-hidden">
        <div className="px-5 py-3 bg-violet-50 border-b border-violet-100">
          <div className="flex items-center gap-2">
            <svg className="w-4 h-4 text-violet-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
                d="M9.813 15.904L9 18.75l-.813-2.846a4.5 4.5 0 00-3.09-3.09L2.25 12l2.846-.813a4.5 4.5 0 003.09-3.09L9 5.25l.813 2.846a4.5 4.5 0 003.09 3.09L15.75 12l-2.846.813a4.5 4.5 0 00-3.09 3.09z" />
            </svg>
            <h3 className="text-sm font-semibold text-violet-900">AI Strategy Brief</h3>
            <span className="text-sm text-violet-500 ml-1">This context will be used by AI agents during runs</span>
          </div>
        </div>
        <div className="px-5 py-4">
          <textarea
            value={notes}
            onChange={e => setNotes(e.target.value)}
            rows={6}
            placeholder={"Share your knowledge and priorities for this campaign. For example:\n\n\u2022 Known SAR: Hinge-binding motif is critical, prefer aminopyrimidine scaffolds\n\u2022 Selectivity: Must avoid CDK2 off-target (structurally similar pocket)\n\u2022 ADMET priorities: Oral bioavailability important, avoid CYP3A4 inhibition\n\u2022 Chemical space: Explore fragments from FBDD screen (MW < 300)\n\u2022 Constraints: LogP < 5, no reactive warheads unless covalent strategy"}
            className="w-full text-sm px-4 py-3 rounded-lg border border-gray-200 focus:outline-none focus:ring-2 focus:ring-violet-200 focus:border-violet-400 resize-none leading-relaxed"
          />
          <div className="flex items-center justify-between mt-3">
            <p className="text-sm text-gray-400">
              The more context you provide, the better the AI agents can optimize compound selection and scoring.
            </p>
            <div className="flex gap-2">
              <button
                onClick={() => { setNotes(campaign?.strategy_notes || ''); setEditing(false) }}
                className="px-3 py-1.5 text-sm text-gray-500 border border-gray-200 rounded-lg hover:bg-gray-50 transition-colors"
              >
                Cancel
              </button>
              <button
                onClick={handleSave}
                disabled={saving}
                className={`px-4 py-1.5 text-sm font-semibold rounded-lg transition-colors ${
                  saving ? 'bg-gray-200 text-gray-400' : 'bg-violet-600 text-white hover:bg-violet-700'
                }`}
              >
                {saving ? 'Saving...' : 'Save Brief'}
              </button>
            </div>
          </div>
        </div>
      </div>
    )
  }

  // Display mode — has content
  return (
    <div className="card overflow-hidden">
      <div className="px-5 py-3 bg-violet-50/50 border-b border-violet-50 flex items-center justify-between">
        <div className="flex items-center gap-2">
          <svg className="w-4 h-4 text-violet-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
              d="M9.813 15.904L9 18.75l-.813-2.846a4.5 4.5 0 00-3.09-3.09L2.25 12l2.846-.813a4.5 4.5 0 003.09-3.09L9 5.25l.813 2.846a4.5 4.5 0 003.09 3.09L15.75 12l-2.846.813a4.5 4.5 0 00-3.09 3.09z" />
          </svg>
          <h3 className="text-sm font-semibold text-gray-700">AI Strategy Brief</h3>
          <span className="text-sm text-violet-400 bg-violet-50 px-2 py-0.5 rounded">Used by AI agents</span>
        </div>
        <button
          onClick={() => { setNotes(campaign.strategy_notes || ''); setEditing(true) }}
          className="text-sm text-violet-600 hover:text-violet-800 transition-colors flex items-center gap-1"
        >
          <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
              d="M11 5H6a2 2 0 00-2 2v11a2 2 0 002 2h11a2 2 0 002-2v-5m-1.414-9.414a2 2 0 112.828 2.828L11.828 15H9v-2.828l8.586-8.586z" />
          </svg>
          Edit
        </button>
      </div>
      <div className="px-5 py-4">
        <p className="text-sm text-gray-700 leading-relaxed whitespace-pre-line">{campaign.strategy_notes}</p>
      </div>
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
          <label className="block text-sm font-semibold text-gray-500 mb-1">Project Name</label>
          <input
            value={name}
            onChange={e => setName(e.target.value)}
            className="w-full text-sm px-3 py-2 rounded-lg border border-gray-200 focus:outline-none focus:ring-2 focus:ring-bx-mint/20 focus:border-bx-mint"
          />
        </div>
        <div>
          <label className="block text-sm font-semibold text-gray-500 mb-1">Description</label>
          <textarea
            value={description}
            onChange={e => setDescription(e.target.value)}
            rows={2}
            className="w-full text-sm px-3 py-2 rounded-lg border border-gray-200 resize-none focus:outline-none focus:ring-2 focus:ring-bx-mint/20 focus:border-bx-mint"
          />
        </div>
      </div>
      <div className="flex gap-2">
        <button
          onClick={() => onSave({ name, description })}
          className="btn-secondary-dark"
        >
          Save
        </button>
        <button
          onClick={onCancel}
          className="px-4 py-1.5 bg-gray-100 text-gray-600 text-sm font-medium rounded-lg hover:bg-gray-200 transition-colors"
        >
          Cancel
        </button>
      </div>
    </div>
  )
}

// ---------------------------------------------------------------------------
// Target Summary — rich display of all target data
// ---------------------------------------------------------------------------

function TargetSummary({ project, onEdit }) {
  const tp = project.target_preview || {}
  const uniprot = tp.uniprot || {}
  const chembl = tp.chembl || {}
  const pubchem = tp.pubchem || {}
  const structure = tp.structure || {}
  const pockets = project.pockets_detected || tp.pockets || []
  const selectedPocketIdx = tp.selected_pocket_index ?? 0
  const selectedPocket = pockets[selectedPocketIdx] || null
  const pdbUrl = structure.download_url || null

  const [funcExpanded, setFuncExpanded] = useState(false)

  return (
    <div className="space-y-4">
      {/* Header bar */}
      <div className="card overflow-hidden">
        <div className="bg-gradient-to-r from-bx-surface to-bx-elevated px-5 py-4">
          <div className="flex items-start justify-between">
            <div>
              <h3 className="text-white font-bold text-base">{project.target_name}</h3>
              <div className="flex items-center gap-2 mt-0.5 flex-wrap">
                {uniprot.gene && (
                  <span className="text-blue-200 text-sm font-mono font-semibold">{uniprot.gene}</span>
                )}
                {uniprot.organism && (
                  <span className="text-blue-300 text-sm italic">{uniprot.organism}</span>
                )}
                {project.target_input_value && (
                  <span className="text-sm font-mono bg-white/10 text-blue-100 px-2 py-0.5 rounded">
                    {project.target_input_value}
                  </span>
                )}
                {uniprot.seqLen > 0 && (
                  <span className="text-sm text-blue-300">{uniprot.seqLen} aa</span>
                )}
              </div>
            </div>
            <button
              onClick={onEdit}
              className="flex items-center gap-1.5 text-sm text-blue-200 hover:text-white transition-colors bg-white/10 px-3 py-1.5 rounded-lg"
            >
              <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                  d="M11 5H6a2 2 0 00-2 2v11a2 2 0 002 2h11a2 2 0 002-2v-5m-1.414-9.414a2 2 0 112.828 2.828L11.828 15H9v-2.828l8.586-8.586z" />
              </svg>
              Edit target
            </button>
          </div>
        </div>

        {/* Quick stats row */}
        <div className="px-5 py-3 flex items-center gap-6 flex-wrap text-sm">
          {project.target_pdb_id && (
            <span className="flex items-center gap-1.5">
              <span className="text-sm text-gray-400">PDB</span>
              <span className="font-mono font-semibold text-gray-700">{project.target_pdb_id}</span>
            </span>
          )}
          {project.structure_method && (
            <span className="text-gray-500">
              <span className="font-medium text-gray-700">{project.structure_method}</span>
              {project.structure_resolution && ` at ${project.structure_resolution.toFixed(2)} \u00C5`}
            </span>
          )}
          {pockets.length > 0 && (
            <span className="text-gray-500">
              <span className="font-medium text-gray-700">{pockets.length}</span> pocket{pockets.length > 1 ? 's' : ''} detected
            </span>
          )}
          <span className="text-gray-500">
            <span className="font-medium text-gray-700">{(chembl.n_actives || 0).toLocaleString()}</span> ChEMBL actives
          </span>
          <span className="text-gray-500">
            <span className="font-medium text-gray-700">{(project.pubchem_compounds_count || pubchem.n_compounds || 0).toLocaleString()}</span> PubChem compounds
          </span>
        </div>
      </div>

      {/* Two-column: 3D viewer + info panels */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
        {/* Left: 3D Viewer + Domains/Keywords underneath */}
        <div className="space-y-4">
          {pdbUrl ? (
            <ProteinViewer
              pdbUrl={pdbUrl}
              selectedPocket={selectedPocket}
              uniprotFeatures={uniprot.activeSites || uniprot.bindingSites || uniprot.domains ? {
                activeSites: uniprot.activeSites,
                bindingSites: uniprot.bindingSites,
                domains: uniprot.domains,
              } : null}
              height={400}
            />
          ) : (
            <div className="h-[400px] bg-[#0f1923] rounded-xl flex items-center justify-center border border-gray-200">
              <p className="text-gray-500 text-sm">No 3D structure available</p>
            </div>
          )}

        </div>

        {/* Right: Function + Pockets + ChEMBL + Diseases */}
        <div className="space-y-4">
          {/* Function — truncated with expand */}
          {uniprot.function && (
            <div className="card px-5 py-4">
              <p className="text-sm font-semibold text-gray-400 uppercase tracking-wider mb-1.5">Function</p>
              <div className="relative">
                <p className={`text-sm text-gray-700 leading-relaxed ${!funcExpanded ? 'line-clamp-4' : ''}`}>
                  {uniprot.function}
                </p>
                {uniprot.function.length > 300 && (
                  <button
                    onClick={() => setFuncExpanded(v => !v)}
                    className="mt-1 text-sm text-bx-light-text font-medium hover:underline"
                  >
                    {funcExpanded ? 'Show less' : 'Read more...'}
                  </button>
                )}
              </div>
            </div>
          )}

          {/* Pockets */}
          {pockets.length > 0 && (
            <div className="card px-5 py-4">
              <p className="text-sm font-semibold text-gray-400 uppercase tracking-wider mb-2">Binding Pockets</p>
              <div className="space-y-2">
                {pockets.map((p, i) => {
                  const prob = Math.round((p.probability || 0) * 100)
                  const residues = Array.isArray(p.residues)
                    ? (typeof p.residues[0] === 'string' && p.residues[0].includes(' ') ? p.residues[0].split(' ') : p.residues)
                    : []
                  const isSelected = i === selectedPocketIdx
                  const center = Array.isArray(p.center) ? p.center : null
                  return (
                    <div key={i} className={`rounded-lg p-2.5 ${isSelected ? 'bg-amber-50 border border-amber-200' : 'bg-gray-50'}`}>
                      <div className="flex items-center justify-between">
                        <div className="flex items-center gap-2">
                          <span className={`w-5 h-5 rounded-full text-[10px] font-bold flex items-center justify-center flex-shrink-0 ${isSelected ? 'bg-amber-400 text-white' : 'bg-gray-200 text-gray-500'}`}>
                            {i + 1}
                          </span>
                          <div>
                            <span className="text-sm font-medium text-gray-700">Pocket #{i + 1}</span>
                            {isSelected && <span className="text-[10px] text-amber-500 font-semibold ml-1.5">selected</span>}
                          </div>
                        </div>
                        <span className={`text-sm font-bold tabular-nums ${prob >= 80 ? 'text-green-600' : prob >= 50 ? 'text-amber-600' : 'text-red-500'}`}>
                          {prob}%
                        </span>
                      </div>
                      <div className="mt-1.5 ml-7 flex flex-wrap gap-x-3 gap-y-0.5 text-[10px] text-gray-400">
                        {p.method && <span>Method: <span className="text-gray-600 font-medium">{p.method}</span></span>}
                        <span>Residues: <span className="text-gray-600 font-medium">{residues.length}</span></span>
                        {p.p2rank_score != null && <span>P2Rank: <span className="text-gray-600 font-medium tabular-nums">{p.p2rank_score.toFixed(2)}</span></span>}
                        {p.volume > 0 && <span>Volume: <span className="text-gray-600 font-medium tabular-nums">{Math.round(p.volume)} A³</span></span>}
                        {center && (
                          <span>Center: <span className="text-gray-600 font-medium tabular-nums">
                            ({center[0].toFixed(1)}, {center[1].toFixed(1)}, {center[2].toFixed(1)})
                          </span></span>
                        )}
                      </div>
                      {residues.length > 0 && (
                        <p className="mt-1 ml-7 text-[9px] text-gray-400 leading-relaxed break-all">
                          {residues.slice(0, 15).join(', ')}{residues.length > 15 ? `, +${residues.length - 15} more` : ''}
                        </p>
                      )}
                    </div>
                  )
                })}
              </div>
            </div>
          )}

          {/* ChEMBL info — compact */}
          <div className="bg-gray-50 border border-gray-100 rounded-lg px-4 py-2.5 flex items-center justify-between">
            <p className="text-xs font-semibold text-gray-500 uppercase tracking-wider">ChEMBL</p>
            <div className="flex items-center gap-3 text-xs">
              <span className="text-gray-600"><strong className="text-gray-700">{(chembl.n_actives || 0).toLocaleString()}</strong> actives</span>
              {chembl.n_with_ic50 > 0 && (
                <span className="text-gray-500"><strong className="text-gray-600">{chembl.n_with_ic50?.toLocaleString()}</strong> IC50</span>
              )}
              {chembl.target_chembl_id && (
                <span className="text-gray-400 font-mono text-[10px]">{chembl.target_chembl_id}</span>
              )}
            </div>
          </div>
          {/* PubChem info — compact */}
          <div className="bg-blue-50 border border-blue-100 rounded-lg px-4 py-2.5 flex items-center justify-between">
            <p className="text-xs font-semibold text-blue-500 uppercase tracking-wider">PubChem</p>
            <div className="flex items-center gap-3 text-xs">
              <span className="text-blue-600"><strong className="text-blue-700">{(project.pubchem_compounds_count || pubchem.n_compounds || 0).toLocaleString()}</strong> compounds</span>
              {pubchem.gene_name && (
                <span className="text-blue-400 font-mono text-[10px]">{pubchem.gene_name}</span>
              )}
            </div>
          </div>

          {/* Details — collapsible: Domains + Keywords + Diseases */}
          {(uniprot.domains?.length > 0 || uniprot.keywords?.length > 0 || uniprot.diseases?.length > 0) && (
            <details className="card overflow-hidden group">
              <summary className="px-5 py-3 flex items-center gap-2 cursor-pointer select-none hover:bg-gray-50 transition-colors">
                <svg className="w-3.5 h-3.5 text-gray-400 transition-transform group-open:rotate-90" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
                </svg>
                <span className="text-sm font-semibold text-gray-400 uppercase tracking-wider">Details</span>
              </summary>
              <div className="px-5 pb-4 space-y-3">
                {uniprot.domains?.length > 0 && (
                  <div>
                    <p className="text-xs font-semibold text-gray-400 uppercase tracking-wider mb-1.5">Protein Domains</p>
                    <div className="flex flex-wrap gap-1.5">
                      {uniprot.domains.map((d, i) => (
                        <span key={i} className="text-sm bg-indigo-50 text-indigo-700 px-2 py-0.5 rounded-md">
                          {d.name} ({d.start}-{d.end})
                        </span>
                      ))}
                    </div>
                  </div>
                )}
                {uniprot.keywords?.length > 0 && (
                  <div>
                    <p className="text-xs font-semibold text-gray-400 uppercase tracking-wider mb-1.5">Keywords</p>
                    <div className="flex flex-wrap gap-1.5">
                      {uniprot.keywords.slice(0, 15).map((k, i) => (
                        <span key={i} className="text-sm bg-blue-50 text-blue-700 px-2 py-0.5 rounded-md">{k}</span>
                      ))}
                    </div>
                  </div>
                )}
                {uniprot.diseases?.length > 0 && (
                  <div>
                    <p className="text-xs font-semibold text-gray-400 uppercase tracking-wider mb-1.5">Associated Diseases</p>
                    <div className="flex flex-wrap gap-1.5">
                      {uniprot.diseases.slice(0, 8).map((d, i) => (
                        <span key={i} className="text-sm bg-red-50 text-red-700 px-2 py-0.5 rounded-md">{d}</span>
                      ))}
                    </div>
                  </div>
                )}
              </div>
            </details>
          )}

          {/* Receptor Preparation Report */}
          <ReceptorPrepReport report={project.receptor_prep_report} config={project.receptor_prep_config} projectId={project.id} hasTarget={!!project.target_preview} />
        </div>
      </div>
    </div>
  )
}

// ---------------------------------------------------------------------------
// ReceptorPrepReport — detailed preparation pipeline status
// ---------------------------------------------------------------------------

const PREP_STEPS = [
  { key: 'P2_cleaning',       label: 'Smart Cleaning',         desc: 'Remove waters, buffers; keep cofactors & metals' },
  { key: 'P3_nonstandard',    label: 'Non-standard Residues',  desc: 'Convert modified residues (MSE→MET, SEP→SER…)' },
  { key: 'P4_altloc',         label: 'Alternate Conformations', desc: 'Keep highest-occupancy conformation only' },
  { key: 'P5_missing_atoms',  label: 'Missing Atoms',          desc: 'Rebuild incomplete side chains' },
  { key: 'P6_pka_analysis',   label: 'pKa Analysis (PropKa)',  desc: 'Calculate pKa of ionizable residues' },
  { key: 'P7_hydrogens',      label: 'Add Hydrogens',          desc: 'Protonate at pH 7.4 with PDBFixer' },
  { key: 'P8_minimization',   label: 'Energy Minimization',    desc: 'AMBER14 + Cα restraints (AlphaFold structures)' },
  { key: 'P9_conversion',     label: 'PDBQT Conversion',       desc: 'Convert to docking format via Meeko/obabel' },
]

function ReceptorPrepReport({ report, config, projectId, hasTarget }) {
  const [expanded, setExpanded] = useState(false)
  const [preparing, setPreparing] = useState(false)
  const [error, setError] = useState(null)
  const [showAdvanced, setShowAdvanced] = useState(false)
  const [ph, setPh] = useState(7.4)
  const [keepMetals, setKeepMetals] = useState(true)
  const [keepCofactors, setKeepCofactors] = useState(true)
  const [minimize, setMinimize] = useState('auto') // 'auto' | 'yes' | 'no'

  // No target configured — don't show anything
  if (!hasTarget) return null

  // Not yet prepared and not in progress — show "Prepare" button
  if (!report && !config) {
    const handlePrepare = async () => {
      setPreparing(true)
      setError(null)
      try {
        const { v9PrepareReceptor } = await import('../api')
        await v9PrepareReceptor(projectId, {
          ph,
          keep_metals: keepMetals,
          keep_cofactors: keepCofactors,
          minimize: minimize === 'auto' ? null : minimize === 'yes',
        })
        window.location.reload()
      } catch (err) {
        setError(err.userMessage || err.message || 'Preparation failed')
        setPreparing(false)
      }
    }
    return (
      <div className="card px-5 py-4">
        <p className="text-sm font-semibold text-gray-400 uppercase tracking-wider mb-2">Receptor Preparation</p>
        <p className="text-xs text-gray-500 mb-3">
          7-step pipeline: cleaning, non-standard residues, protonation, energy minimization, PDBQT conversion.
        </p>

        {/* Advanced settings — collapsible */}
        <details open={showAdvanced} onToggle={e => setShowAdvanced(e.target.open)} className="mb-3">
          <summary className="text-[11px] font-medium text-gray-400 cursor-pointer hover:text-gray-600 select-none">
            Advanced settings
          </summary>
          <div className="mt-2 space-y-2.5 pl-1">
            {/* pH */}
            <div className="flex items-center justify-between gap-3">
              <div>
                <span className="text-xs font-medium text-gray-600">Protonation pH</span>
                <p className="text-[10px] text-gray-400">Standard physiological = 7.4</p>
              </div>
              <input
                type="number"
                step="0.1"
                min="0"
                max="14"
                value={ph}
                onChange={e => setPh(parseFloat(e.target.value) || 7.4)}
                className="w-16 text-sm text-center px-2 py-1 rounded border border-gray-200 focus:outline-none focus:ring-1 focus:ring-bx-mint/40 focus:border-bx-mint"
              />
            </div>
            {/* Keep metals */}
            <label className="flex items-center justify-between gap-3 cursor-pointer">
              <div>
                <span className="text-xs font-medium text-gray-600">Keep catalytic metals</span>
                <p className="text-[10px] text-gray-400">ZN, MG, FE, MN, CA, CO…</p>
              </div>
              <input type="checkbox" checked={keepMetals} onChange={e => setKeepMetals(e.target.checked)}
                className="w-4 h-4 rounded border-gray-300 text-bx-mint focus:ring-bx-mint/30" />
            </label>
            {/* Keep cofactors */}
            <label className="flex items-center justify-between gap-3 cursor-pointer">
              <div>
                <span className="text-xs font-medium text-gray-600">Keep cofactors</span>
                <p className="text-[10px] text-gray-400">NAD, FAD, HEM, ATP, PLP…</p>
              </div>
              <input type="checkbox" checked={keepCofactors} onChange={e => setKeepCofactors(e.target.checked)}
                className="w-4 h-4 rounded border-gray-300 text-bx-mint focus:ring-bx-mint/30" />
            </label>
            {/* Minimization */}
            <div className="flex items-center justify-between gap-3">
              <div>
                <span className="text-xs font-medium text-gray-600">Energy minimization</span>
                <p className="text-[10px] text-gray-400">Auto = on for AlphaFold, off for PDB</p>
              </div>
              <select
                value={minimize}
                onChange={e => setMinimize(e.target.value)}
                className="text-xs px-2 py-1 rounded border border-gray-200 focus:outline-none focus:ring-1 focus:ring-bx-mint/40 focus:border-bx-mint bg-white"
              >
                <option value="auto">Auto</option>
                <option value="yes">Force on</option>
                <option value="no">Force off</option>
              </select>
            </div>
          </div>
        </details>

        {error && <p className="text-xs text-red-500 mb-2">{error}</p>}
        <button
          onClick={handlePrepare}
          disabled={preparing}
          className="btn-primary text-sm"
        >
          {preparing ? (
            <span className="flex items-center gap-2">
              <svg className="w-3.5 h-3.5 animate-spin" fill="none" viewBox="0 0 24 24">
                <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
                <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
              </svg>
              Preparing…
            </span>
          ) : 'Prepare Receptor'}
        </button>
      </div>
    )
  }

  // Config set but no report yet — preparation in progress
  if (!report) {
    return (
      <div className="card px-5 py-4">
        <div className="flex items-center gap-2">
          <svg className="w-4 h-4 text-amber-500 animate-pulse" fill="currentColor" viewBox="0 0 20 20">
            <path fillRule="evenodd" d="M10 18a8 8 0 100-16 8 8 0 000 16zm1-12a1 1 0 10-2 0v4a1 1 0 00.293.707l2.828 2.829a1 1 0 101.415-1.415L11 9.586V6z" clipRule="evenodd" />
          </svg>
          <p className="text-sm font-semibold text-amber-600">Receptor preparation in progress…</p>
        </div>
        <p className="text-xs text-gray-400 mt-1">This may take 1-2 minutes. Refresh the page to check.</p>
      </div>
    )
  }

  const stepsCompleted = report.steps_completed || []
  const stepsSkipped = report.steps_skipped || []
  const conversionMethod = report.conversion_method || 'unknown'
  const shiftedPka = report.shifted_pka_residues || []
  const nonstdResidues = report.nonstandard_residues || []
  const missingAtoms = report.missing_atoms_added ?? 0
  const minimizationConverged = report.minimization_converged
  const cleaning = report.cleaning || {}
  const removedHetCount = cleaning.hetatm_removed ?? 0
  const keptHet = cleaning.het_kept || []
  // Separate metals vs cofactors from het_kept list
  const metalNames = new Set(['ZN','MG','FE','FE2','MN','CA','CO','NI','CU','CU1','MO','SE','W','CD'])
  const keptMetals = keptHet.filter(h => metalNames.has(h))
  const keptCofactors = keptHet.filter(h => !metalNames.has(h))

  const totalSteps = PREP_STEPS.length
  const completedCount = stepsCompleted.length
  const allGood = completedCount >= 5 // at least 5/7 steps = good

  // Determine step status — backend sends tuples as arrays ["P4_altloc", "reason"]
  // and combined keys like "P6+P7_protonation" for grouped skips
  function stepStatus(key) {
    if (stepsCompleted.includes(key)) return 'done'
    const skip = stepsSkipped.find(s => {
      const sKey = Array.isArray(s) ? s[0] : (typeof s === 'string' ? s : s?.step)
      return sKey === key || (sKey && sKey.includes(key.replace(/^P\d+_/, '')))
    })
    if (skip) return 'skipped'
    return 'pending'
  }
  function skipReason(key) {
    const skip = stepsSkipped.find(s => {
      const sKey = Array.isArray(s) ? s[0] : (typeof s === 'string' ? s : s?.step)
      return sKey === key || (sKey && sKey.includes(key.replace(/^P\d+_/, '')))
    })
    if (!skip) return null
    if (Array.isArray(skip)) return skip[1] || 'Skipped'
    return typeof skip === 'string' ? 'Dependency not available' : skip.reason
  }

  return (
    <details className="card overflow-hidden group" open={expanded} onToggle={e => setExpanded(e.target.open)}>
      <summary className="px-5 py-3 flex items-center gap-2 cursor-pointer select-none hover:bg-gray-50 transition-colors">
        <svg className="w-3.5 h-3.5 text-gray-400 transition-transform group-open:rotate-90" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
        </svg>
        <span className="text-sm font-semibold text-gray-400 uppercase tracking-wider">Receptor Preparation</span>
        <span className={`ml-auto text-xs font-bold px-2 py-0.5 rounded-full ${allGood ? 'bg-green-100 text-green-700' : 'bg-amber-100 text-amber-700'}`}>
          {completedCount}/{totalSteps} steps
        </span>
      </summary>

      <div className="px-5 pb-4 space-y-3">
        {/* Pipeline steps */}
        <div className="space-y-1.5">
          {PREP_STEPS.map(step => {
            const status = stepStatus(step.key)
            const reason = skipReason(step.key)
            return (
              <div key={step.key} className="flex items-start gap-2">
                {status === 'done' ? (
                  <svg className="w-4 h-4 text-green-500 mt-0.5 shrink-0" fill="currentColor" viewBox="0 0 20 20">
                    <path fillRule="evenodd" d="M10 18a8 8 0 100-16 8 8 0 000 16zm3.707-9.293a1 1 0 00-1.414-1.414L9 10.586 7.707 9.293a1 1 0 00-1.414 1.414l2 2a1 1 0 001.414 0l4-4z" clipRule="evenodd" />
                  </svg>
                ) : status === 'skipped' ? (
                  <svg className="w-4 h-4 text-amber-400 mt-0.5 shrink-0" fill="currentColor" viewBox="0 0 20 20">
                    <path fillRule="evenodd" d="M10 18a8 8 0 100-16 8 8 0 000 16zm1-11a1 1 0 10-2 0v3.586L7.707 9.293a1 1 0 00-1.414 1.414l3 3a1 1 0 001.414 0l3-3a1 1 0 00-1.414-1.414L11 10.586V7z" clipRule="evenodd" />
                  </svg>
                ) : (
                  <div className="w-4 h-4 rounded-full border-2 border-gray-200 mt-0.5 shrink-0" />
                )}
                <div className="flex-1 min-w-0">
                  <div className="flex items-center gap-2">
                    <span className={`text-sm font-medium ${status === 'done' ? 'text-gray-700' : status === 'skipped' ? 'text-amber-600' : 'text-gray-400'}`}>
                      {step.label}
                    </span>
                    {status === 'skipped' && reason && (
                      <span className="text-[10px] text-amber-500 bg-amber-50 px-1.5 py-0.5 rounded">{reason}</span>
                    )}
                  </div>
                  <p className="text-[10px] text-gray-400">{step.desc}</p>
                </div>
              </div>
            )
          })}
        </div>

        {/* Details section — only when we have interesting data */}
        {(removedHetCount > 0 || keptMetals.length > 0 || nonstdResidues.length > 0 || missingAtoms > 0 || shiftedPka.length > 0 || conversionMethod !== 'unknown') && (
          <div className="border-t border-gray-100 pt-3 space-y-2">
            {/* Cleaning stats */}
            {(removedHetCount > 0 || keptMetals.length > 0 || keptCofactors.length > 0) && (
              <div className="flex flex-wrap gap-x-4 gap-y-1 text-[11px]">
                {removedHetCount > 0 && (
                  <span className="text-gray-500">
                    <strong className="text-gray-700">{removedHetCount}</strong> HETATM removed
                  </span>
                )}
                {keptMetals.length > 0 && (
                  <span className="text-gray-500">
                    Metals kept: <strong className="text-gray-700">{keptMetals.join(', ')}</strong>
                  </span>
                )}
                {keptCofactors.length > 0 && (
                  <span className="text-gray-500">
                    Cofactors: <strong className="text-gray-700">{keptCofactors.join(', ')}</strong>
                  </span>
                )}
              </div>
            )}

            {/* Non-standard residues */}
            {nonstdResidues.length > 0 && (
              <div className="text-[11px] text-gray-500">
                <span className="font-medium text-gray-600">Non-standard residues replaced:</span>{' '}
                {nonstdResidues.map((r, i) => (
                  <span key={i}>
                    {i > 0 && ', '}
                    <span className="font-mono text-gray-700">{typeof r === 'string' ? r : `${r.from}→${r.to}`}</span>
                  </span>
                ))}
              </div>
            )}

            {/* Missing atoms */}
            {missingAtoms > 0 && (
              <div className="text-[11px] text-gray-500">
                <strong className="text-gray-700">{missingAtoms}</strong> missing atoms reconstructed
              </div>
            )}

            {/* Shifted pKa residues */}
            {shiftedPka.length > 0 && (
              <div>
                <p className="text-[11px] font-medium text-amber-600 mb-1">Shifted pKa residues (>1 unit from model):</p>
                <div className="flex flex-wrap gap-1">
                  {shiftedPka.map((r, i) => (
                    <span key={i} className="text-[10px] font-mono bg-amber-50 text-amber-700 px-1.5 py-0.5 rounded">
                      {typeof r === 'string' ? r : `${r.residue} ${r.chain || ''}${r.resnum || ''} pKa=${r.pka?.toFixed(1) ?? '?'}`}
                    </span>
                  ))}
                </div>
              </div>
            )}

            {/* Minimization — only show if it was actually attempted */}
            {stepsCompleted.includes('P8_minimization') && (
              <div className="text-[11px] text-gray-500">
                Minimization: <strong className={minimizationConverged ? 'text-green-600' : 'text-amber-600'}>
                  {minimizationConverged ? 'converged' : 'did not converge (500 iterations)'}
                </strong>
              </div>
            )}

            {/* Conversion method */}
            <div className="text-[11px] text-gray-500">
              PDBQT conversion: <strong className="text-gray-700 font-mono">{conversionMethod}</strong>
            </div>

            {/* Config used */}
            {config && (
              <div className="text-[10px] text-gray-400 mt-1">
                pH {config.ph ?? 7.4} · metals {config.keep_metals !== false ? 'on' : 'off'} · cofactors {config.keep_cofactors !== false ? 'on' : 'off'} · minimization {config.minimize == null ? 'auto' : config.minimize ? 'on' : 'off'}
              </div>
            )}
          </div>
        )}

        {/* Re-prepare with different settings */}
        <div className="border-t border-gray-100 pt-2 mt-1">
          {!showAdvanced ? (
            <button
              onClick={() => {
                // Pre-fill with current config values
                setPh(config?.ph ?? 7.4)
                setKeepMetals(config?.keep_metals !== false)
                setKeepCofactors(config?.keep_cofactors !== false)
                setMinimize(config?.minimize == null ? 'auto' : config?.minimize ? 'yes' : 'no')
                setShowAdvanced(true)
              }}
              className="text-[11px] text-gray-400 hover:text-bx-light-text transition-colors"
            >
              Re-prepare with different settings →
            </button>
          ) : (
            <div className="mt-2 space-y-2.5">
              <p className="text-xs font-medium text-gray-600 mb-2">Preparation settings</p>
              {/* pH */}
              <div className="flex items-center justify-between gap-3">
                <div>
                  <span className="text-xs font-medium text-gray-600">Protonation pH</span>
                  <p className="text-[10px] text-gray-400">Standard physiological = 7.4</p>
                </div>
                <input
                  type="number"
                  step="0.1"
                  min="0"
                  max="14"
                  value={ph}
                  onChange={e => setPh(parseFloat(e.target.value) || 7.4)}
                  className="w-16 text-sm text-center px-2 py-1 rounded border border-gray-200 focus:outline-none focus:ring-1 focus:ring-bx-mint/40 focus:border-bx-mint"
                />
              </div>
              {/* Keep metals */}
              <label className="flex items-center justify-between gap-3 cursor-pointer">
                <div>
                  <span className="text-xs font-medium text-gray-600">Keep catalytic metals</span>
                  <p className="text-[10px] text-gray-400">ZN, MG, FE, MN, CA, CO…</p>
                </div>
                <input type="checkbox" checked={keepMetals} onChange={e => setKeepMetals(e.target.checked)}
                  className="w-4 h-4 rounded border-gray-300 text-bx-mint focus:ring-bx-mint/30" />
              </label>
              {/* Keep cofactors */}
              <label className="flex items-center justify-between gap-3 cursor-pointer">
                <div>
                  <span className="text-xs font-medium text-gray-600">Keep cofactors</span>
                  <p className="text-[10px] text-gray-400">NAD, FAD, HEM, ATP, PLP…</p>
                </div>
                <input type="checkbox" checked={keepCofactors} onChange={e => setKeepCofactors(e.target.checked)}
                  className="w-4 h-4 rounded border-gray-300 text-bx-mint focus:ring-bx-mint/30" />
              </label>
              {/* Minimization */}
              <div className="flex items-center justify-between gap-3">
                <div>
                  <span className="text-xs font-medium text-gray-600">Energy minimization</span>
                  <p className="text-[10px] text-gray-400">Auto = on for AlphaFold, off for PDB</p>
                </div>
                <select
                  value={minimize}
                  onChange={e => setMinimize(e.target.value)}
                  className="text-xs px-2 py-1 rounded border border-gray-200 focus:outline-none focus:ring-1 focus:ring-bx-mint/40 focus:border-bx-mint bg-white"
                >
                  <option value="auto">Auto</option>
                  <option value="yes">Force on</option>
                  <option value="no">Force off</option>
                </select>
              </div>
              {/* Action buttons */}
              <div className="flex items-center gap-2 pt-1">
                <button
                  onClick={async () => {
                    setPreparing(true)
                    setError(null)
                    try {
                      const { v9PrepareReceptor } = await import('../api')
                      await v9PrepareReceptor(projectId, {
                        force: true,
                        ph,
                        keep_metals: keepMetals,
                        keep_cofactors: keepCofactors,
                        minimize: minimize === 'auto' ? null : minimize === 'yes',
                      })
                      window.location.reload()
                    } catch (err) {
                      setError(err.userMessage || err.message || 'Re-preparation failed')
                      setPreparing(false)
                    }
                  }}
                  disabled={preparing}
                  className="btn-primary text-xs py-1.5 px-3"
                >
                  {preparing ? 'Re-preparing…' : 'Re-prepare'}
                </button>
                <button
                  onClick={() => setShowAdvanced(false)}
                  className="text-xs text-gray-400 hover:text-gray-600 transition-colors"
                >
                  Cancel
                </button>
              </div>
            </div>
          )}
          {error && <p className="text-[10px] text-red-500 mt-1">{error}</p>}
        </div>
      </div>
    </details>
  )
}

// ---------------------------------------------------------------------------
// ProjectHome
// ---------------------------------------------------------------------------

export default function ProjectHome() {
  const { projectId } = useParams()
  const navigate = useNavigate()
  const { projects, selectProject, updateProject, deleteProject, createCampaign, updateCampaign, createPhase, createRun } = useWorkspace()
  const { addToast } = useToast()

  const [editing, setEditing] = useState(false)
  const [showPhaseCreator, setShowPhaseCreator] = useState(false)
  const [creatingForPhase, setCreatingForPhase] = useState(null)
  const [campaignCreating, setCampaignCreating] = useState(false)
  const [showAgentPanel, setShowAgentPanel] = useState(false)

  useEffect(() => {
    if (projectId) selectProject(projectId)
  }, [projectId, selectProject])

  const project = projects.find(p => p.id === projectId)

  const campaign = project?.campaigns?.[0] || null
  const dbPhases = campaign?.phases || []
  const hasTarget = !!project?.target_preview

  // Campaign KPIs aggregated from all phases (must be before early return — Rules of Hooks)
  const campaignKpi = useMemo(() => {
    if (!dbPhases.length) return null
    return dbPhases.reduce((acc, p) => {
      const s = p.stats || {}
      acc.totalMolecules += s.total_molecules ?? 0
      acc.totalBookmarked += s.bookmarked ?? 0
      acc.totalRuns += (s.runs_completed ?? 0) + (s.runs_running ?? 0)
      return acc
    }, { totalMolecules: 0, totalBookmarked: 0, totalRuns: 0 })
  }, [dbPhases])

  if (!project) {
    return (
      <div className="text-center py-20">
        <p className="text-gray-400 text-sm">Project not found.</p>
        <Link to="/" className="text-bx-light-text text-sm hover:underline mt-2 inline-block">
          Back to projects
        </Link>
      </div>
    )
  }

  // Build funnel: always show 3 slots (A/B/C), merge with existing DB phases
  const FUNNEL_TEMPLATE = ['hit_discovery', 'hit_to_lead', 'lead_optimization']
  const phases = FUNNEL_TEMPLATE.map(type => {
    const existing = dbPhases.find(p => p.type === type)
    if (existing) return existing
    return { id: `placeholder-${type}`, type, label: PHASE_TYPES[type]?.label || type, status: 'not_created' }
  })

  return (
    <div className="space-y-5">
      {/* Back navigation */}
      <div>
        <Link
          to="/"
          className="inline-flex items-center gap-1.5 text-sm text-gray-400 hover:text-bx-light-text transition-colors"
        >
          <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 19l-7-7 7-7" />
          </svg>
          All Projects
        </Link>
      </div>

      {/* ---- Project header ---- */}
      <div className="card px-6 py-5">
        <div className="flex items-start justify-between gap-4">
          <div className="flex-1 min-w-0">
            <div className="flex items-center gap-3 flex-wrap mb-1.5">
              <h1 className="text-[1.6rem] font-bold text-bx-light-text leading-tight font-serif">{project.name}</h1>
              <span className={`badge ${project.status === 'active' ? 'badge-active' : 'badge-empty'}`}>
                {project.status === 'active' ? 'Active' : project.status || 'Active'}
              </span>
            </div>

            {project.description && (
              <p className="text-bx-light-text2 text-[.78rem] leading-relaxed mb-3 max-w-2xl">
                {project.description}
              </p>
            )}

            {/* Target metadata row — design system tags */}
            <div className="flex items-center flex-wrap gap-2 text-sm">
              {project.target_name && (
                <span className="tag tag-gene-l">{project.target_name}</span>
              )}
              {project.target_input_value && (
                <span className="tag tag-uni-l">{project.target_input_value}</span>
              )}
              {project.target_pdb_id && (
                <span className="tag tag-pdb-l">PDB: {project.target_pdb_id}</span>
              )}
              {project.updated_at && (
                <span className="text-bx-light-muted text-sm font-mono ml-auto">
                  Updated {fmtDate(project.updated_at)}
                </span>
              )}
            </div>
          </div>

          {/* Action buttons */}
          <div className="flex items-center gap-2 shrink-0">
            <button
              className="btn-secondary-light"
              onClick={() => setEditing(e => !e)}
            >
              <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                  d="M11 5H6a2 2 0 00-2 2v11a2 2 0 002 2h11a2 2 0 002-2v-5m-1.414-9.414a2 2 0 112.828 2.828L11.828 15H9v-2.828l8.586-8.586z" />
              </svg>
              Edit
            </button>
            <button
              className="btn-secondary-light"
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
                  addToast(err.userMessage || 'Failed to update project', 'error')
                }
                setEditing(false)
              }}
              onCancel={() => setEditing(false)}
            />
          </div>
        )}
      </div>

      {/* ---- Target Setup Gate ---- */}
      {!hasTarget && (
        <div className="card px-6 py-10 text-center border-dashed">
          <div className="mb-3 mx-auto w-fit">
            <BindXLogo variant="idle" size={48} />
          </div>
          <h3 className="text-lg font-bold text-bx-light-text mb-1 font-serif">Set Up Your Target</h3>
          <p className="text-[.78rem] text-bx-light-text2 mb-5 max-w-md mx-auto">
            Before creating campaigns, you need to configure your protein target, select a 3D structure, and choose a binding pocket.
          </p>
          <button
            onClick={() => navigate(`/project/${projectId}/target-setup`)}
            className="btn-primary text-sm py-2.5 px-5"
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                d="M9.75 3.104v5.714a2.25 2.25 0 01-.659 1.591L5 14.5M9.75 3.104c-.251.023-.501.05-.75.082m.75-.082a24.301 24.301 0 014.5 0m0 0v5.714c0 .597.237 1.17.659 1.591L19 14.5M14.25 3.104c.251.023.501.05.75.082M19 14.5l-1.47 4.9a2.25 2.25 0 01-2.156 1.6H8.626a2.25 2.25 0 01-2.156-1.6L5 14.5m14 0H5" />
            </svg>
            Configure Target
          </button>
        </div>
      )}

      {/* ---- Target Summary (when configured) ---- */}
      {hasTarget && <TargetSummary project={project} onEdit={() => navigate(`/project/${projectId}/target-setup`)} />}

      {/* ---- Campaign section ---- */}
      {campaign ? (
        <>
          {/* Campaign header */}
          <div className="flex items-center justify-between">
            <div>
              <p className="text-sm text-gray-400 uppercase tracking-wider font-semibold mb-0.5">Campaign</p>
              <h2 className="text-lg font-bold text-gray-800">{campaign.name}</h2>
            </div>
            <div className="flex items-center gap-4">
              {campaignKpi && (
                <div className="flex items-center gap-3 text-sm">
                  <span className="text-gray-500">
                    <strong className="text-gray-700">{campaignKpi.totalMolecules}</strong> molecules
                  </span>
                  <span className="text-gray-300">|</span>
                  <span className="text-gray-500">
                    <strong className="text-yellow-600">{campaignKpi.totalBookmarked}</strong> bookmarked
                  </span>
                  <span className="text-gray-300">|</span>
                  <span className="text-gray-500">
                    <strong className="text-gray-700">{campaignKpi.totalRuns}</strong> runs
                  </span>
                </div>
              )}
              <span className="text-sm text-gray-400">
                {dbPhases.length}/{FUNNEL_TEMPLATE.length} phases
              </span>
              <button
                onClick={() => setShowAgentPanel(true)}
                className="flex items-center gap-1.5 px-3 py-1.5 rounded-lg text-sm font-medium bg-gradient-to-r from-bx-mint/10 to-emerald-50 text-emerald-700 border border-emerald-200 hover:border-emerald-300 hover:bg-emerald-50 transition-colors"
              >
                <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                    d="M9.75 17L9 20l-1 1h8l-1-1-.75-3M3 13h18M5 17H3a2 2 0 01-2-2V5a2 2 0 012-2h14a2 2 0 012 2v10a2 2 0 01-2 2h-2" />
                </svg>
                AI Agent
              </button>
            </div>
          </div>

          {/* AI Strategy Brief */}
          <CampaignStrategyBrief campaign={campaign} onSave={updateCampaign} />

          {/* Drug Discovery Funnel */}
          <div className="card px-6 py-6">
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
                : 'bg-bx-surface text-white hover:bg-bx-elevated'
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
                <BindXLogo variant="loading" size={16} />
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

      {/* AI Agent panel */}
      <CampaignAgentPanel
        isOpen={showAgentPanel}
        onClose={() => setShowAgentPanel(false)}
        project={project}
        campaign={campaign}
      />

      {/* Phase creator modal */}
      {showPhaseCreator && (
        <PhaseCreator
          campaignId={campaign?.id}
          existingPhases={dbPhases}
          onCreatePhase={async (newPhaseConfig) => {
            try {
              const created = await createPhase(newPhaseConfig.campaign_id, { type: newPhaseConfig.type })
              addToast('Phase created', 'success')

              // Auto-import bookmarks from previous phase if available
              if (created?.id && dbPhases.length > 0) {
                const prevPhase = dbPhases[dbPhases.length - 1]
                const prevBookmarked = prevPhase.stats?.bookmarked ?? 0
                if (prevBookmarked > 0) {
                  try {
                    await createRun(created.id, {
                      type: 'import',
                      config: {
                        source: 'phase_selection',
                        selection: 'bookmarked',
                        source_phase_id: prevPhase.id,
                      },
                    })
                    addToast(`Importing ${prevBookmarked} bookmarked molecules from ${prevPhase.label || 'previous phase'}...`, 'success')
                  } catch (_) {
                    // Non-blocking — phase is created even if auto-import fails
                  }
                }
              }

              // Auto-navigate to the new phase
              if (created?.id) {
                navigate(`/project/${projectId}/phase/${created.id}`)
              }
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
