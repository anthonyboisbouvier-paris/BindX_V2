import React, { useMemo, useEffect } from 'react'
import { useParams, useSearchParams, useNavigate, Link } from 'react-router-dom'
import { useWorkspace } from '../contexts/WorkspaceContext.jsx'
import { flattenMoleculeProperties, PHASE_TYPES } from '../lib/columns.js'
import ComparisonView from '../components/ComparisonView.jsx'

export default function ComparisonPage() {
  const { projectId, phaseId } = useParams()
  const [searchParams] = useSearchParams()
  const navigate = useNavigate()
  const {
    currentProject,
    currentPhase,
    phaseMolecules,
    selectProject,
    selectPhase,
    toggleBookmark,
  } = useWorkspace()

  // Sync URL params to workspace state
  useEffect(() => { if (projectId) selectProject(projectId) }, [projectId]) // eslint-disable-line react-hooks/exhaustive-deps
  useEffect(() => { if (phaseId) selectPhase(phaseId) }, [phaseId]) // eslint-disable-line react-hooks/exhaustive-deps

  // Parse molecule IDs from query params
  const moleculeIds = useMemo(() => {
    const ids = searchParams.get('ids')
    return ids ? ids.split(',').filter(Boolean) : []
  }, [searchParams])

  // Get flat molecules for selected IDs
  const selectedMolecules = useMemo(() => {
    if (!phaseMolecules.length || !moleculeIds.length) return []
    const idSet = new Set(moleculeIds)
    return phaseMolecules
      .filter(m => idSet.has(String(m.id)))
      .map(flattenMoleculeProperties)
  }, [phaseMolecules, moleculeIds])

  const phaseTypeMeta = currentPhase ? (PHASE_TYPES[currentPhase.type] || {}) : {}

  const handleBack = () => {
    navigate(`/project/${projectId}/phase/${phaseId}`)
  }

  if (!currentPhase) {
    return (
      <div className="flex items-center justify-center h-64">
        <div className="w-8 h-8 border-2 border-bx-accent border-t-transparent rounded-full animate-spin" />
      </div>
    )
  }

  return (
    <div className="space-y-4 pb-8">
      {/* Breadcrumb */}
      <nav className="flex items-center gap-1.5 text-sm text-gray-400">
        <Link to="/" className="hover:text-bx-mint transition-colors">Projects</Link>
        <svg className="w-3 h-3 text-gray-300" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
        </svg>
        <Link to={`/project/${projectId}`} className="hover:text-bx-light-text transition-colors truncate max-w-[120px]">
          {currentProject?.name || 'Project'}
        </Link>
        <svg className="w-3 h-3 text-gray-300" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
        </svg>
        <Link to={`/project/${projectId}/phase/${phaseId}`} className="hover:text-bx-light-text transition-colors">
          Phase {phaseTypeMeta.short || '?'}
        </Link>
        <svg className="w-3 h-3 text-gray-300" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
        </svg>
        <span className="text-gray-600 font-medium">Comparison ({selectedMolecules.length})</span>
      </nav>

      {/* Title */}
      <div className="card px-5 py-4">
        <div className="flex items-center justify-between">
          <div>
            <h1 className="text-xl font-bold text-bx-light-text">Molecule Comparison</h1>
            <p className="text-sm text-gray-400 mt-0.5">{selectedMolecules.length} molecules side by side</p>
          </div>
          <button
            onClick={handleBack}
            className="flex items-center gap-1.5 px-3 py-2 rounded-lg text-sm font-medium border border-gray-200 text-gray-600 hover:bg-gray-50 transition-colors"
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10 19l-7-7m0 0l7-7m-7 7h18" />
            </svg>
            Dashboard
          </button>
        </div>
      </div>

      {selectedMolecules.length < 2 ? (
        <div className="card p-12 text-center">
          <p className="text-gray-400">No molecules selected for comparison. Go back to the dashboard and select 2-6 molecules.</p>
        </div>
      ) : (
        <ComparisonView
          molecules={selectedMolecules}
          onBack={handleBack}
          onBookmark={(molId) => toggleBookmark(molId)}
        />
      )}
    </div>
  )
}
