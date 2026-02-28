import React, { useState, useCallback, useEffect } from 'react'
import { useParams, useNavigate } from 'react-router-dom'
import { useWorkspace } from '../contexts/WorkspaceContext.jsx'
import { useToast } from '../contexts/ToastContext.jsx'
import { previewTarget, detectPockets } from '../api.js'
import ProteinViewer from '../components/ProteinViewer.jsx'

// ---------------------------------------------------------------------------
// UniProt public API helper
// ---------------------------------------------------------------------------

async function fetchUniProtInfo(uniprotId) {
  const res = await fetch(`https://rest.uniprot.org/uniprotkb/${uniprotId}.json`)
  if (!res.ok) {
    if (res.status === 400 || res.status === 404) {
      throw new Error(`UniProt ID "${uniprotId}" not found. Please check the ID and try again (e.g. P00533 for EGFR).`)
    }
    throw new Error(`UniProt service error (${res.status}). Please try again later.`)
  }
  const data = await res.json()

  const name = data.proteinDescription?.recommendedName?.fullName?.value
    || data.proteinDescription?.submittedName?.[0]?.fullName?.value
    || uniprotId
  const gene = data.genes?.[0]?.geneName?.value || ''
  const organism = data.organism?.scientificName || ''
  const seqLen = data.sequence?.length || 0

  const funcComment = data.comments?.find(c => c.commentType === 'FUNCTION')
  const funcText = funcComment?.texts?.[0]?.value || ''

  const diseases = (data.comments || [])
    .filter(c => c.commentType === 'DISEASE')
    .map(c => c.disease?.diseaseId || c.texts?.[0]?.value)
    .filter(Boolean)

  const keywords = (data.keywords || []).map(k => k.name)

  const features = data.features || []
  const activeSites = features
    .filter(f => f.type === 'Active site')
    .map(f => ({ position: f.location?.start?.value }))
  const bindingSites = features
    .filter(f => f.type === 'Binding site')
    .map(f => ({ start: f.location?.start?.value, end: f.location?.end?.value }))
  const domains = features
    .filter(f => f.type === 'Domain')
    .map(f => ({ name: f.description, start: f.location?.start?.value, end: f.location?.end?.value }))

  return {
    name, gene, organism, seqLen,
    function: funcText,
    diseases, keywords,
    activeSites, bindingSites, domains,
  }
}

// ---------------------------------------------------------------------------
// Step indicator
// ---------------------------------------------------------------------------

function StepIndicator({ current }) {
  const steps = [
    { key: 'input', label: 'Input Target' },
    { key: 'preview', label: 'Review & Pockets' },
    { key: 'validate', label: 'Validate' },
  ]
  return (
    <div className="flex items-center gap-2 mb-6">
      {steps.map((s, i) => {
        const idx = steps.findIndex(x => x.key === current)
        const done = i < idx
        const active = i === idx
        return (
          <div key={s.key} className="flex items-center gap-2">
            {i > 0 && <div className={`w-8 h-px ${done || active ? 'bg-[#0f131d]' : 'bg-gray-200'}`} />}
            <div className="flex items-center gap-1.5">
              <span className={`w-6 h-6 rounded-full text-xs font-bold flex items-center justify-center ${
                done ? 'bg-[#00e6a0] text-white'
                  : active ? 'bg-[#0f131d] text-white'
                  : 'bg-gray-200 text-gray-400'
              }`}>
                {done ? '\u2713' : i + 1}
              </span>
              <span className={`text-xs font-medium ${active ? 'text-[#0f131d]' : 'text-gray-400'}`}>
                {s.label}
              </span>
            </div>
          </div>
        )
      })}
    </div>
  )
}

// ---------------------------------------------------------------------------
// Pocket card
// ---------------------------------------------------------------------------

function PocketCard({ pocket, index, selected, onSelect }) {
  const prob = Math.round((pocket.probability || 0) * 100)
  const residues = Array.isArray(pocket.residues)
    ? (typeof pocket.residues[0] === 'string' && pocket.residues[0].includes(' ')
      ? pocket.residues[0].split(' ')
      : pocket.residues)
    : []

  return (
    <button
      onClick={() => onSelect(index)}
      className={`w-full text-left p-3 rounded-lg border-2 transition-all ${
        selected
          ? 'border-amber-400 bg-amber-50 shadow-sm'
          : 'border-gray-200 bg-white hover:border-gray-300'
      }`}
    >
      <div className="flex items-center justify-between mb-2">
        <div className="flex items-center gap-2">
          <span className={`w-5 h-5 rounded-full text-xs font-bold flex items-center justify-center ${
            selected ? 'bg-amber-400 text-white' : 'bg-gray-200 text-gray-500'
          }`}>
            {index + 1}
          </span>
          <span className="text-sm font-semibold text-gray-700">
            Pocket #{index + 1}
          </span>
          {pocket.method && (
            <span className="text-xs text-gray-400 bg-gray-100 px-1.5 py-0.5 rounded">
              {pocket.method}
            </span>
          )}
        </div>
        <span className={`text-sm font-bold ${prob >= 80 ? 'text-green-600' : prob >= 50 ? 'text-amber-600' : 'text-red-500'}`}>
          {prob}%
        </span>
      </div>
      <p className="text-xs text-gray-500 mb-1">{pocket.explanation}</p>
      <p className="text-xs text-gray-400">{residues.length} residues</p>
    </button>
  )
}

// ---------------------------------------------------------------------------
// TargetSetup — main component
// ---------------------------------------------------------------------------

export default function TargetSetup() {
  const { projectId } = useParams()
  const navigate = useNavigate()
  const { projects, updateProject, refreshProjects } = useWorkspace()
  const { addToast } = useToast()

  const project = projects.find(p => p.id === projectId)

  // Input state
  const [inputMode, setInputMode] = useState('uniprot')
  const [uniprotId, setUniprotId] = useState('')
  const [fastaInput, setFastaInput] = useState('')

  // Loading + error
  const [loading, setLoading] = useState(false)
  const [loadingStep, setLoadingStep] = useState('')
  const [error, setError] = useState(null)

  // Results
  const [uniprotInfo, setUniprotInfo] = useState(null)
  const [previewData, setPreviewData] = useState(null)
  const [selectedPocketIdx, setSelectedPocketIdx] = useState(0)
  const [selectedStructureIdx, setSelectedStructureIdx] = useState(0)

  // Pockets cache per structure source (keyed by download_url)
  const [pocketsCache, setPocketsCache] = useState({})
  const [pocketsLoading, setPocketsLoading] = useState(false)

  // Saving
  const [saving, setSaving] = useState(false)
  const [initialized, setInitialized] = useState(false)

  // ---------------------------------------------------------------------------
  // Pre-fill from saved project data
  // ---------------------------------------------------------------------------

  useEffect(() => {
    if (initialized || !project) return
    const tp = project.target_preview
    if (!tp) { setInitialized(true); return }

    // Restore state from saved target_preview
    if (project.target_input_value) {
      setUniprotId(project.target_input_value)
    }
    if (tp.uniprot) {
      setUniprotInfo(tp.uniprot)
    }

    // Reconstruct previewData from saved target_preview
    // Prefer saved structures list; fallback to single structure + alphafold_structure (old format)
    let structures = tp.structures || []
    if (structures.length === 0) {
      if (tp.structure) structures.push(tp.structure)
      if (tp.alphafold_structure) structures.push(tp.alphafold_structure)
    }

    setPreviewData({
      protein_name: project.target_name,
      structures: structures.length > 0 ? structures : undefined,
      structure: tp.structure || structures[0] || null,
      pockets: tp.pockets || project.pockets_detected || [],
      chembl_info: tp.chembl || null,
    })
    setSelectedPocketIdx(tp.selected_pocket_index ?? 0)
    setSelectedStructureIdx(0)
    setInitialized(true)
  }, [project, initialized])

  // Step tracking
  const step = previewData ? 'preview' : 'input'

  // ---------------------------------------------------------------------------
  // Search handler
  // ---------------------------------------------------------------------------

  const handleSearch = useCallback(async () => {
    const id = uniprotId.trim().toUpperCase()
    if (!id) return

    setLoading(true)
    setError(null)
    setPreviewData(null)
    setUniprotInfo(null)

    try {
      // Step 1: Fetch UniProt info
      setLoadingStep('Fetching protein data from UniProt...')
      const info = await fetchUniProtInfo(id)
      setUniprotInfo(info)

      // Step 2: Backend — structure + pocket detection + ChEMBL
      setLoadingStep('Resolving structure & detecting pockets (P2Rank)...')
      const preview = await previewTarget(id)

      // Backend already returns all structure options (PDB + AlphaFold + ESMFold)
      // No need to fetch AlphaFold separately — it's included in preview.structures
      setPreviewData(preview)
      setSelectedPocketIdx(0)
      setSelectedStructureIdx(0)

      // Cache pockets for the initial (PDB) structure
      const firstUrl = preview.structures?.[0]?.download_url
      if (firstUrl && preview.pockets?.length > 0) {
        setPocketsCache({ [firstUrl]: preview.pockets })
      } else {
        setPocketsCache({})
      }

    } catch (err) {
      setError(err.userMessage || err.message || 'Failed to fetch target data')
    } finally {
      setLoading(false)
      setLoadingStep('')
    }
  }, [uniprotId])

  // ---------------------------------------------------------------------------
  // Validate & save target
  // ---------------------------------------------------------------------------

  const handleValidate = useCallback(async () => {
    if (!previewData || !projectId) return

    // Use pockets for the currently selected structure
    const currentPdbUrl = previewData.structures?.[selectedStructureIdx]?.download_url
    const currentPockets = (currentPdbUrl && pocketsCache[currentPdbUrl])
      ? pocketsCache[currentPdbUrl]
      : (selectedStructureIdx === 0 ? (previewData.pockets || []) : [])

    setSaving(true)
    try {
      const selectedStructure = previewData.structures?.[selectedStructureIdx] || previewData.structure

      await updateProject(projectId, {
        target_input_type: 'uniprot',
        target_input_value: uniprotId.trim().toUpperCase(),
        target_name: uniprotInfo?.name || previewData.protein_name,
        target_pdb_id: selectedStructure?.pdb_id || null,
        structure_source: selectedStructure?.source || null,
        structure_resolution: selectedStructure?.resolution || null,
        structure_method: selectedStructure?.method || null,
        cocrystal_ligand: selectedStructure?.ligand_name || null,
        target_preview: {
          uniprot: uniprotInfo,
          chembl: previewData.chembl_info,
          structure: selectedStructure,
          structures: previewData.structures || [],
          pockets: currentPockets,
          selected_pocket_index: selectedPocketIdx,
        },
        pockets_detected: currentPockets,
        chembl_actives_count: previewData.chembl_info?.n_actives || null,
        chembl_median_ic50: null,
      })

      await refreshProjects()
      addToast('Target validated and saved', 'success')
      navigate(`/project/${projectId}`)
    } catch (err) {
      addToast(err.userMessage || err.message || 'Failed to save target', 'error')
    } finally {
      setSaving(false)
    }
  }, [previewData, projectId, selectedStructureIdx, selectedPocketIdx, uniprotId, uniprotInfo, updateProject, refreshProjects, addToast, navigate, pocketsCache])

  // ---------------------------------------------------------------------------
  // Selected structure for 3D viewer
  // ---------------------------------------------------------------------------

  const selectedStructure = previewData?.structures?.[selectedStructureIdx] || previewData?.structure
  const pdbUrl = selectedStructure?.download_url || null

  // Pockets for the currently selected structure (from cache or initial data)
  const pocketsForStructure = pdbUrl && pocketsCache[pdbUrl]
    ? pocketsCache[pdbUrl]
    : (selectedStructureIdx === 0 ? (previewData?.pockets || []) : [])
  const selectedPocket = pocketsForStructure[selectedPocketIdx] || null

  // ---------------------------------------------------------------------------
  // Auto-detect pockets when switching to a structure not yet in cache
  // ---------------------------------------------------------------------------

  const handleStructureChange = useCallback(async (url, structure, structIdx) => {
    if (!url || !previewData) return
    // Index 0 = initial structure, pockets already loaded from preview-target
    if (structIdx === 0 && previewData.pockets?.length > 0) return
    // Already cached
    if (pocketsCache[url]) return

    setPocketsLoading(true)
    setSelectedPocketIdx(0)
    try {
      const result = await detectPockets(url, structure?.label || 'structure', structure?.ligand_id)
      setPocketsCache(prev => ({ ...prev, [url]: result.pockets || [] }))
    } catch (err) {
      console.warn('[TargetSetup] Pocket detection failed:', err.message)
      setPocketsCache(prev => ({ ...prev, [url]: [] }))
    } finally {
      setPocketsLoading(false)
    }
  }, [previewData, pocketsCache])

  // ---------------------------------------------------------------------------
  // Render
  // ---------------------------------------------------------------------------

  return (
    <div className="space-y-5">
      {/* Back link */}
      <button
        onClick={() => navigate(`/project/${projectId}`)}
        className="inline-flex items-center gap-1.5 text-xs text-gray-400 hover:text-[#0f131d] transition-colors"
      >
        <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 19l-7-7 7-7" />
        </svg>
        Back to project
      </button>

      <StepIndicator current={step} />

      {/* ========== STEP 1: Input ========== */}
      <div className="bg-white rounded-xl border border-gray-100 shadow-sm px-6 py-5">
        <h2 className="text-lg font-bold text-[#0f131d] mb-1">Target Setup</h2>
        <p className="text-sm text-gray-500 mb-5">
          Enter a UniProt ID to fetch protein data, resolve 3D structure, and detect binding pockets.
        </p>

        {/* Mode tabs */}
        <div className="flex gap-2 mb-4">
          <button
            onClick={() => setInputMode('uniprot')}
            className={`px-3 py-1.5 text-xs font-semibold rounded-lg transition-colors ${
              inputMode === 'uniprot'
                ? 'bg-[#0f131d] text-white'
                : 'bg-gray-100 text-gray-500 hover:bg-gray-200'
            }`}
          >
            UniProt ID
          </button>
          <button
            onClick={() => setInputMode('fasta')}
            className={`px-3 py-1.5 text-xs font-semibold rounded-lg transition-colors ${
              inputMode === 'fasta'
                ? 'bg-[#0f131d] text-white'
                : 'bg-gray-100 text-gray-500 hover:bg-gray-200'
            }`}
          >
            FASTA Sequence
          </button>
        </div>

        {inputMode === 'uniprot' ? (
          <div className="flex gap-3">
            <input
              type="text"
              value={uniprotId}
              onChange={e => setUniprotId(e.target.value)}
              onKeyDown={e => e.key === 'Enter' && handleSearch()}
              placeholder="e.g. P00533 (EGFR)"
              className="flex-1 px-4 py-2.5 text-sm rounded-lg border border-gray-200 focus:outline-none focus:ring-2 focus:ring-[#0f131d]/20 focus:border-[#0f131d] font-mono"
              disabled={loading}
            />
            <button
              onClick={handleSearch}
              disabled={loading || !uniprotId.trim()}
              className={`px-5 py-2.5 text-sm font-semibold rounded-lg transition-colors flex items-center gap-2 ${
                loading || !uniprotId.trim()
                  ? 'bg-gray-200 text-gray-400 cursor-not-allowed'
                  : 'bg-[#0f131d] text-white hover:bg-[#141925]'
              }`}
            >
              {loading ? (
                <>
                  <svg className="w-4 h-4 animate-spin" fill="none" viewBox="0 0 24 24">
                    <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
                    <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
                  </svg>
                  Searching...
                </>
              ) : (previewData ? 'Re-search' : 'Search')}
            </button>
          </div>
        ) : (
          <div className="space-y-3">
            <textarea
              value={fastaInput}
              onChange={e => setFastaInput(e.target.value)}
              placeholder=">protein_name&#10;MTEYKLVVVGAGGVGKSALTIQLIQNHFVDE..."
              rows={5}
              className="w-full px-4 py-3 text-sm rounded-lg border border-gray-200 focus:outline-none focus:ring-2 focus:ring-[#0f131d]/20 focus:border-[#0f131d] font-mono resize-none"
              disabled={true}
            />
            <div className="bg-amber-50 border border-amber-100 rounded-lg px-4 py-3 space-y-2">
              <p className="text-xs text-amber-700 font-medium">Structure prediction from FASTA sequence</p>
              <div className="grid grid-cols-1 sm:grid-cols-2 gap-2 text-xs text-amber-600">
                <div className="flex items-start gap-2">
                  <span className="text-amber-400 mt-0.5">&#9679;</span>
                  <span><strong>ESMFold</strong> — Fast prediction (~2 min). Requires GPU backend. <em>Coming soon.</em></span>
                </div>
                <div className="flex items-start gap-2">
                  <span className="text-amber-400 mt-0.5">&#9679;</span>
                  <span><strong>AlphaFold</strong> — Pre-computed structures are auto-fetched when using UniProt ID input.</span>
                </div>
              </div>
            </div>
          </div>
        )}

        {/* Loading step indicator */}
        {loading && loadingStep && (
          <div className="mt-4 flex items-center gap-2 text-sm text-blue-600">
            <svg className="w-4 h-4 animate-spin" fill="none" viewBox="0 0 24 24">
              <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
              <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
            </svg>
            {loadingStep}
          </div>
        )}

        {/* Error */}
        {error && (
          <div className="mt-4 px-4 py-3 bg-red-50 border border-red-100 rounded-lg text-sm text-red-600">
            {error}
          </div>
        )}
      </div>

      {/* ========== STEP 2: Preview ========== */}
      {previewData && (
        <>
          {/* UniProt protein info */}
          {uniprotInfo && (
            <div className="bg-white rounded-xl border border-gray-100 shadow-sm overflow-hidden">
              <div className="bg-gradient-to-r from-[#0f131d] to-[#141925] px-5 py-4">
                <div className="flex items-center gap-3">
                  <div>
                    <h3 className="text-white font-bold text-base">{uniprotInfo.name}</h3>
                    <div className="flex items-center gap-2 mt-0.5">
                      {uniprotInfo.gene && (
                        <span className="text-blue-200 text-sm font-mono">{uniprotInfo.gene}</span>
                      )}
                      <span className="text-blue-300 text-xs italic">{uniprotInfo.organism}</span>
                      <span className="text-xs font-mono bg-white/10 text-blue-100 px-2 py-0.5 rounded">
                        {uniprotId.toUpperCase()}
                      </span>
                      <span className="text-xs text-blue-300">{uniprotInfo.seqLen} aa</span>
                    </div>
                  </div>
                </div>
              </div>

              <div className="px-5 py-4 space-y-3">
                {uniprotInfo.function && (
                  <div>
                    <p className="text-xs font-semibold text-gray-500 uppercase tracking-wider mb-1">Function</p>
                    <p className="text-sm text-gray-700 leading-relaxed line-clamp-3">{uniprotInfo.function}</p>
                  </div>
                )}
                {uniprotInfo.diseases?.length > 0 && (
                  <div>
                    <p className="text-xs font-semibold text-gray-500 uppercase tracking-wider mb-1">Associated Diseases</p>
                    <div className="flex flex-wrap gap-1">
                      {uniprotInfo.diseases.slice(0, 5).map((d, i) => (
                        <span key={i} className="text-xs bg-red-50 text-red-700 px-2 py-0.5 rounded-md">{d}</span>
                      ))}
                    </div>
                  </div>
                )}
                {uniprotInfo.keywords?.length > 0 && (
                  <div>
                    <p className="text-xs font-semibold text-gray-500 uppercase tracking-wider mb-1">Keywords</p>
                    <div className="flex flex-wrap gap-1">
                      {uniprotInfo.keywords.slice(0, 10).map((k, i) => (
                        <span key={i} className="text-xs bg-blue-50 text-blue-700 px-2 py-0.5 rounded-md">{k}</span>
                      ))}
                    </div>
                  </div>
                )}
                {uniprotInfo.domains?.length > 0 && (
                  <div>
                    <p className="text-xs font-semibold text-gray-500 uppercase tracking-wider mb-1">Domains</p>
                    <div className="flex flex-wrap gap-1">
                      {uniprotInfo.domains.map((d, i) => (
                        <span key={i} className="text-xs bg-indigo-50 text-indigo-700 px-2 py-0.5 rounded-md">
                          {d.name} ({d.start}-{d.end})
                        </span>
                      ))}
                    </div>
                  </div>
                )}
              </div>
            </div>
          )}

          {/* Structure + 3D viewer + Pockets — two-column layout */}
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-5">
            {/* Left: 3D viewer */}
            <div>
              <h3 className="text-sm font-semibold text-gray-700 mb-3 flex items-center gap-2">
                <svg className="w-4 h-4 text-gray-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8} d="M14 10l-2 1m0 0l-2-1m2 1v2.5M20 7l-2 1m2-1l-2-1m2 1v2.5M14 4l-2-1-2 1M4 7l2-1M4 7l2 1M4 7v2.5M12 21l-2-1m2 1l2-1m-2 1v-2.5M6 18l-2-1v-2.5M18 18l2-1v-2.5" />
                </svg>
                3D Structure
                {selectedStructure?.pdb_id && (
                  <span className="text-xs font-mono bg-blue-50 text-blue-700 px-2 py-0.5 rounded">
                    {selectedStructure.pdb_id}
                  </span>
                )}
                {selectedStructure?.source === 'alphafold' && (
                  <span className="text-xs bg-purple-50 text-purple-700 px-2 py-0.5 rounded">AlphaFold</span>
                )}
              </h3>
              {pdbUrl ? (
                <ProteinViewer
                  pdbUrl={pdbUrl}
                  selectedPocket={selectedPocket}
                  uniprotFeatures={uniprotInfo ? {
                    activeSites: uniprotInfo.activeSites,
                    bindingSites: uniprotInfo.bindingSites,
                    domains: uniprotInfo.domains,
                  } : null}
                  height={420}
                />
              ) : (
                <div className="h-[420px] bg-[#0f1923] rounded-xl flex items-center justify-center">
                  <p className="text-gray-500 text-sm">No structure available for 3D view</p>
                </div>
              )}
            </div>

            {/* Right: Structure selection + Pockets */}
            <div className="space-y-5">
              {/* Structure selection — show available structures (hide ESMFold until supported) */}
              {previewData.structures?.filter(s => s.source !== 'esmfold').length > 0 && (
                <div>
                  <h3 className="text-sm font-semibold text-gray-700 mb-2">Available Structures</h3>
                  <div className="space-y-2">
                    {previewData.structures.filter(s => s.source !== 'esmfold').map((s, i) => {
                      // Find the real index in the full structures array for selection
                      const realIdx = previewData.structures.indexOf(s)
                      return (
                      <button
                        key={i}
                        onClick={() => {
                          setSelectedStructureIdx(realIdx)
                          const struct = previewData.structures[realIdx]
                          handleStructureChange(struct?.download_url, struct, realIdx)
                        }}
                        className={`w-full text-left p-3 rounded-lg border-2 transition-all ${
                          selectedStructureIdx === realIdx
                            ? 'border-blue-400 bg-blue-50'
                            : 'border-gray-200 bg-white hover:border-gray-300'
                        }`}
                      >
                        <div className="flex items-center justify-between">
                          <div className="flex items-center gap-2">
                            <span className="text-sm font-semibold text-gray-700">
                              {s.pdb_id || s.label}
                            </span>
                            <span className={`text-xs px-1.5 py-0.5 rounded ${
                              s.source === 'alphafold'
                                ? 'bg-purple-100 text-purple-700'
                                : 'bg-gray-100 text-gray-500'
                            }`}>
                              {s.source === 'alphafold' ? 'AlphaFold' : (s.method || 'Experimental')}
                            </span>
                          </div>
                          {s.resolution && (
                            <span className="text-xs text-gray-500">
                              {s.resolution.toFixed(2)} \u00C5
                            </span>
                          )}
                          {s.source === 'alphafold' && s.confidence && (
                            <span className="text-xs text-purple-500">
                              pLDDT: {s.confidence.toFixed(1)}
                            </span>
                          )}
                        </div>
                        {s.ligand_name && (
                          <p className="text-xs text-amber-600 mt-1">Co-crystal: {s.ligand_name}</p>
                        )}
                        {s.recommended && (
                          <span className="text-xs text-green-600 font-medium">Recommended</span>
                        )}
                        {s.source === 'alphafold' && (
                          <p className="text-xs text-purple-500 mt-1">Full-length predicted structure (may look different from PDB which only covers crystallized domains)</p>
                        )}
                      </button>
                    )})}
                  </div>
                </div>
              )}

              {/* Pockets */}
              <div>
                <h3 className="text-sm font-semibold text-gray-700 mb-2 flex items-center gap-2">
                  Detected Pockets
                  {pocketsLoading ? (
                    <span className="text-xs text-gray-400 font-normal flex items-center gap-1">
                      <span className="w-3 h-3 border border-blue-400 border-t-transparent rounded-full animate-spin" />
                      Detecting pockets...
                    </span>
                  ) : (
                    <span className="text-xs text-gray-400 font-normal">
                      ({pocketsForStructure.length} found via P2Rank)
                    </span>
                  )}
                </h3>
                {pocketsLoading ? (
                  <div className="text-center py-6 bg-gray-50 rounded-lg">
                    <div className="w-8 h-8 border-2 border-blue-400 border-t-transparent rounded-full animate-spin mx-auto mb-2" />
                    <p className="text-gray-400 text-sm">Running P2Rank on this structure...</p>
                  </div>
                ) : pocketsForStructure.length > 0 ? (
                  <div className="space-y-2">
                    {pocketsForStructure.map((p, i) => (
                      <PocketCard
                        key={i}
                        pocket={p}
                        index={i}
                        selected={selectedPocketIdx === i}
                        onSelect={setSelectedPocketIdx}
                      />
                    ))}
                  </div>
                ) : (
                  <div className="text-center py-6 bg-gray-50 rounded-lg">
                    <p className="text-gray-400 text-sm">No pockets detected on this structure</p>
                  </div>
                )}
              </div>

              {/* ChEMBL info */}
              {previewData.chembl_info?.has_data && (
                <div className="bg-green-50 border border-green-100 rounded-lg p-4">
                  <h3 className="text-sm font-semibold text-green-800 mb-2">ChEMBL Data Available</h3>
                  <div className="grid grid-cols-2 gap-3">
                    <div className="text-center">
                      <p className="text-xl font-bold text-green-700">
                        {previewData.chembl_info.n_actives?.toLocaleString()}
                      </p>
                      <p className="text-xs text-green-600">Active compounds</p>
                    </div>
                    <div className="text-center">
                      <p className="text-xl font-bold text-green-700">
                        {previewData.chembl_info.n_with_ic50?.toLocaleString()}
                      </p>
                      <p className="text-xs text-green-600">With IC50 data</p>
                    </div>
                  </div>
                  <p className="text-xs text-green-600 mt-2">
                    Target: {previewData.chembl_info.target_chembl_id}
                  </p>
                </div>
              )}
            </div>
          </div>

          {/* Validate button */}
          <div className="flex items-center justify-end gap-3 pt-2">
            <button
              onClick={() => {
                setPreviewData(null)
                setUniprotInfo(null)
              }}
              className="px-4 py-2.5 text-sm text-gray-500 border border-gray-200 rounded-lg hover:bg-gray-50 transition-colors"
            >
              Reset
            </button>
            <button
              onClick={handleValidate}
              disabled={saving || pocketsForStructure.length === 0}
              className={`px-6 py-2.5 text-sm font-semibold rounded-lg transition-colors flex items-center gap-2 ${
                saving || pocketsForStructure.length === 0
                  ? 'bg-gray-200 text-gray-400 cursor-not-allowed'
                  : 'bg-[#00e6a0] text-white hover:bg-[#00c98b]'
              }`}
            >
              {saving ? (
                <>
                  <svg className="w-4 h-4 animate-spin" fill="none" viewBox="0 0 24 24">
                    <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
                    <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
                  </svg>
                  Saving...
                </>
              ) : (
                <>
                  <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 13l4 4L19 7" />
                  </svg>
                  Validate Target
                </>
              )}
            </button>
          </div>
        </>
      )}
    </div>
  )
}
