import React, { useState, useCallback, useEffect, useMemo } from 'react'
import { useParams, useNavigate } from 'react-router-dom'
import { useWorkspace } from '../contexts/WorkspaceContext.jsx'
import { useToast } from '../contexts/ToastContext.jsx'
import { previewTarget, detectPockets } from '../api.js'
import ProteinViewer from '../components/ProteinViewer.jsx'
import ViewerDiagnostic from '../components/ViewerDiagnostic.jsx'
import BindXLogo from '../components/BindXLogo.jsx'

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
// RCSB PDB API helper — fetch structure info from a PDB code
// ---------------------------------------------------------------------------

async function fetchPDBInfo(pdbId) {
  const id = pdbId.trim().toUpperCase()
  const res = await fetch(`https://data.rcsb.org/rest/v1/core/entry/${id}`)
  if (!res.ok) {
    if (res.status === 404) throw new Error(`PDB ID "${id}" not found. Please check the code (e.g. 1M17, 4HJO).`)
    throw new Error(`RCSB PDB service error (${res.status}).`)
  }
  const data = await res.json()

  const title = data.struct?.title || id
  const method = data.exptl?.[0]?.method || null
  const resolution = data.rcsb_entry_info?.resolution_combined?.[0] || null

  // Get ligand IDs (non-polymer entities that are not water)
  const ligands = (data.rcsb_entry_info?.nonpolymer_bound_components || [])
    .filter(l => l !== 'HOH')

  // Try to get UniProt mapping
  let uniprotId = null
  try {
    const mapRes = await fetch(`https://data.rcsb.org/rest/v1/core/polymer_entity/${id}/1`)
    if (mapRes.ok) {
      const mapData = await mapRes.json()
      const uniprotRef = mapData.rcsb_polymer_entity_container_identifiers?.reference_sequence_identifiers
        ?.find(r => r.database_name === 'UniProt')
      uniprotId = uniprotRef?.database_accession || null
    }
  } catch (_) { /* ignore */ }

  return {
    pdb_id: id,
    title,
    method,
    resolution,
    ligands,
    uniprotId,
    download_url: `https://files.rcsb.org/download/${id}.pdb`,
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
            {i > 0 && <div className={`w-8 h-px ${done || active ? 'bg-bx-surface' : 'bg-gray-200'}`} />}
            <div className="flex items-center gap-1.5">
              <span className={`w-6 h-6 rounded-full text-sm font-bold flex items-center justify-center ${
                done ? 'bg-bx-mint text-white'
                  : active ? 'bg-bx-surface text-white'
                  : 'bg-gray-200 text-gray-400'
              }`}>
                {done ? '\u2713' : i + 1}
              </span>
              <span className={`text-sm font-medium ${active ? 'text-bx-light-text' : 'text-gray-400'}`}>
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
          <span className={`w-5 h-5 rounded-full text-sm font-bold flex items-center justify-center ${
            selected ? 'bg-amber-400 text-white' : 'bg-gray-200 text-gray-500'
          }`}>
            {index + 1}
          </span>
          <span className="text-sm font-semibold text-gray-700">
            Pocket #{index + 1}
          </span>
          {pocket.method && (
            <span className="text-sm text-gray-400 bg-gray-100 px-1.5 py-0.5 rounded">
              {pocket.method}
            </span>
          )}
        </div>
        <span className={`text-sm font-bold ${prob >= 80 ? 'text-green-600' : prob >= 50 ? 'text-amber-600' : 'text-red-500'}`}>
          {prob}%
        </span>
      </div>
      <p className="text-sm text-gray-500 mb-1">{pocket.explanation}</p>
      <p className="text-sm text-gray-400">{residues.length} residues</p>
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
  const [pdbCode, setPdbCode] = useState('')
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

  // Live diagnostic from ProteinViewer
  const [viewerDiag, setViewerDiag] = useState(null)

  // Memoize uniprotFeatures to avoid infinite re-render loop
  // (new object reference on every render triggers ProteinViewer useEffect → removeAllSurfaces → addSurface async killed)
  const uniprotFeaturesMemo = useMemo(() => {
    if (!uniprotInfo) return null
    return {
      activeSites: uniprotInfo.activeSites,
      bindingSites: uniprotInfo.bindingSites,
      domains: uniprotInfo.domains,
    }
  }, [uniprotInfo])

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
    if (project.target_input_type === 'pdb') {
      setInputMode('pdb')
      setPdbCode(project.target_input_value || '')
    } else if (project.target_input_value) {
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
  // Search by PDB code
  // ---------------------------------------------------------------------------

  const handleSearchPDB = useCallback(async () => {
    const code = pdbCode.trim().toUpperCase()
    if (!code) return

    setLoading(true)
    setError(null)
    setPreviewData(null)
    setUniprotInfo(null)

    try {
      // Step 1: Fetch PDB info from RCSB
      setLoadingStep('Fetching structure info from RCSB PDB...')
      const pdbInfo = await fetchPDBInfo(code)

      // Step 2: If we found a UniProt mapping, fetch UniProt info
      if (pdbInfo.uniprotId) {
        setLoadingStep('Fetching protein data from UniProt...')
        try {
          const info = await fetchUniProtInfo(pdbInfo.uniprotId)
          setUniprotInfo(info)
          setUniprotId(pdbInfo.uniprotId)
        } catch (_) { /* non-critical */ }
      }

      // Step 3: Detect pockets via backend
      setLoadingStep('Detecting binding pockets (P2Rank)...')
      const ligandId = pdbInfo.ligands?.[0] || null
      const pocketsResult = await detectPockets(pdbInfo.download_url, pdbInfo.pdb_id, ligandId)
      const pockets = pocketsResult.pockets || []

      // Build structure entry
      const structure = {
        pdb_id: pdbInfo.pdb_id,
        source: 'pdb',
        method: pdbInfo.method,
        resolution: pdbInfo.resolution,
        download_url: pdbInfo.download_url,
        label: `${pdbInfo.pdb_id} — ${pdbInfo.title}`,
        ligand_name: ligandId,
        ligand_id: ligandId,
        recommended: true,
      }

      setPreviewData({
        protein_name: pdbInfo.title,
        structures: [structure],
        structure,
        pockets,
        chembl_info: null,
      })
      setSelectedPocketIdx(0)
      setSelectedStructureIdx(0)
      setPocketsCache({ [pdbInfo.download_url]: pockets })

    } catch (err) {
      setError(err.message || 'Failed to fetch PDB data')
    } finally {
      setLoading(false)
      setLoadingStep('')
    }
  }, [pdbCode])

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

      const inputType = inputMode === 'pdb' ? 'pdb' : 'uniprot'
      const inputValue = inputMode === 'pdb' ? pdbCode.trim().toUpperCase() : uniprotId.trim().toUpperCase()

      await updateProject(projectId, {
        target_input_type: inputType,
        target_input_value: inputValue,
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
  }, [previewData, projectId, selectedStructureIdx, selectedPocketIdx, uniprotId, pdbCode, inputMode, uniprotInfo, updateProject, refreshProjects, addToast, navigate, pocketsCache])

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
        className="inline-flex items-center gap-1.5 text-sm text-gray-400 hover:text-bx-light-text transition-colors"
      >
        <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 19l-7-7 7-7" />
        </svg>
        Back to project
      </button>

      <StepIndicator current={step} />

      {/* ========== STEP 1: Input ========== */}
      <div className="card px-6 py-5">
        <h2 className="text-lg font-bold text-bx-light-text mb-1">Target Setup</h2>
        <p className="text-sm text-gray-500 mb-5">
          Enter a UniProt ID, PDB code, or FASTA sequence to set up your target.
        </p>

        {/* Mode tabs */}
        <div className="flex gap-2 mb-4">
          <button
            onClick={() => setInputMode('uniprot')}
            className={`px-3 py-1.5 text-sm font-semibold rounded-lg transition-colors ${
              inputMode === 'uniprot'
                ? 'bg-bx-surface text-white'
                : 'bg-gray-100 text-gray-500 hover:bg-gray-200'
            }`}
          >
            UniProt ID
          </button>
          <button
            onClick={() => setInputMode('pdb')}
            className={`px-3 py-1.5 text-sm font-semibold rounded-lg transition-colors ${
              inputMode === 'pdb'
                ? 'bg-bx-surface text-white'
                : 'bg-gray-100 text-gray-500 hover:bg-gray-200'
            }`}
          >
            PDB Code
          </button>
          <button
            onClick={() => setInputMode('fasta')}
            className={`px-3 py-1.5 text-sm font-semibold rounded-lg transition-colors ${
              inputMode === 'fasta'
                ? 'bg-bx-surface text-white'
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
              className="flex-1 px-4 py-2.5 text-sm rounded-lg border border-gray-200 focus:outline-none focus:ring-2 focus:ring-bx-mint/20 focus:border-bx-surface font-mono"
              disabled={loading}
            />
            <button
              onClick={handleSearch}
              disabled={loading || !uniprotId.trim()}
              className={`px-5 py-2.5 text-sm font-semibold rounded-lg transition-colors flex items-center gap-2 ${
                loading || !uniprotId.trim()
                  ? 'bg-gray-200 text-gray-400 cursor-not-allowed'
                  : 'bg-bx-surface text-white hover:bg-bx-elevated'
              }`}
            >
              {loading ? (
                <>
                  <BindXLogo variant="loading" size={16} />
                  Searching...
                </>
              ) : (previewData ? 'Re-search' : 'Search')}
            </button>
          </div>
        ) : inputMode === 'pdb' ? (
          <div className="flex gap-3">
            <input
              type="text"
              value={pdbCode}
              onChange={e => setPdbCode(e.target.value)}
              onKeyDown={e => e.key === 'Enter' && handleSearchPDB()}
              placeholder="e.g. 1M17, 4HJO, 6LU7"
              className="flex-1 px-4 py-2.5 text-sm rounded-lg border border-gray-200 focus:outline-none focus:ring-2 focus:ring-bx-mint/20 focus:border-bx-surface font-mono"
              disabled={loading}
            />
            <button
              onClick={handleSearchPDB}
              disabled={loading || !pdbCode.trim()}
              className={`px-5 py-2.5 text-sm font-semibold rounded-lg transition-colors flex items-center gap-2 ${
                loading || !pdbCode.trim()
                  ? 'bg-gray-200 text-gray-400 cursor-not-allowed'
                  : 'bg-bx-surface text-white hover:bg-bx-elevated'
              }`}
            >
              {loading ? (
                <>
                  <BindXLogo variant="loading" size={16} />
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
              className="w-full px-4 py-3 text-sm rounded-lg border border-gray-200 focus:outline-none focus:ring-2 focus:ring-bx-mint/20 focus:border-bx-surface font-mono resize-none"
              disabled={true}
            />
            <div className="bg-amber-50 border border-amber-100 rounded-lg px-4 py-3 space-y-2">
              <p className="text-sm text-amber-700 font-medium">Structure prediction from FASTA sequence</p>
              <div className="grid grid-cols-1 sm:grid-cols-2 gap-2 text-sm text-amber-600">
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
          <div className="mt-4 flex items-center gap-3 text-sm text-blue-600">
            <BindXLogo variant="loading" size={24} />
            <span>{loadingStep}</span>
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
            <div className="card overflow-hidden">
              <div className="bg-gradient-to-r from-bx-surface to-bx-elevated px-5 py-4">
                <div className="flex items-center gap-3">
                  <div>
                    <h3 className="text-white font-bold text-base">{uniprotInfo.name}</h3>
                    <div className="flex items-center gap-2 mt-0.5">
                      {uniprotInfo.gene && (
                        <span className="text-blue-200 text-sm font-mono">{uniprotInfo.gene}</span>
                      )}
                      <span className="text-blue-300 text-sm italic">{uniprotInfo.organism}</span>
                      <span className="text-sm font-mono bg-white/10 text-blue-100 px-2 py-0.5 rounded">
                        {uniprotId.toUpperCase()}
                      </span>
                      <span className="text-sm text-blue-300">{uniprotInfo.seqLen} aa</span>
                    </div>
                  </div>
                </div>
              </div>

              <div className="px-5 py-4 space-y-3">
                {uniprotInfo.function && (
                  <div>
                    <p className="text-sm font-semibold text-gray-500 uppercase tracking-wider mb-1">Function</p>
                    <p className="text-sm text-gray-700 leading-relaxed line-clamp-3">{uniprotInfo.function}</p>
                  </div>
                )}
                {uniprotInfo.diseases?.length > 0 && (
                  <div>
                    <p className="text-sm font-semibold text-gray-500 uppercase tracking-wider mb-1">Associated Diseases</p>
                    <div className="flex flex-wrap gap-1">
                      {uniprotInfo.diseases.slice(0, 5).map((d, i) => (
                        <span key={i} className="text-sm bg-red-50 text-red-700 px-2 py-0.5 rounded-md">{d}</span>
                      ))}
                    </div>
                  </div>
                )}
                {uniprotInfo.keywords?.length > 0 && (
                  <div>
                    <p className="text-sm font-semibold text-gray-500 uppercase tracking-wider mb-1">Keywords</p>
                    <div className="flex flex-wrap gap-1">
                      {uniprotInfo.keywords.slice(0, 10).map((k, i) => (
                        <span key={i} className="text-sm bg-blue-50 text-blue-700 px-2 py-0.5 rounded-md">{k}</span>
                      ))}
                    </div>
                  </div>
                )}
                {uniprotInfo.domains?.length > 0 && (
                  <div>
                    <p className="text-sm font-semibold text-gray-500 uppercase tracking-wider mb-1">Domains</p>
                    <div className="flex flex-wrap gap-1">
                      {uniprotInfo.domains.map((d, i) => (
                        <span key={i} className="text-sm bg-indigo-50 text-indigo-700 px-2 py-0.5 rounded-md">
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
                  <span className="text-sm font-mono bg-blue-50 text-blue-700 px-2 py-0.5 rounded">
                    {selectedStructure.pdb_id}
                  </span>
                )}
                {selectedStructure?.source === 'alphafold' && (
                  <span className="text-sm bg-purple-50 text-purple-700 px-2 py-0.5 rounded">AlphaFold</span>
                )}
              </h3>
              {pdbUrl ? (
                <ProteinViewer
                  pdbUrl={pdbUrl}
                  selectedPocket={selectedPocket}
                  uniprotFeatures={uniprotFeaturesMemo}
                  height={420}
                  onDiagnostic={setViewerDiag}
                />
              ) : (
                <div className="h-[420px] bg-[#0f1923] rounded-xl flex items-center justify-center">
                  <p className="text-gray-500 text-sm">No structure available for 3D view</p>
                </div>
              )}

              {/* Viewer Diagnostic — below 3D viewer */}
              <ViewerDiagnostic
                pdbUrl={pdbUrl}
                selectedPocket={selectedPocket}
                uniprotFeatures={uniprotFeaturesMemo}
                selectedStructure={selectedStructure}
                live={viewerDiag}
              />
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
                            <span className={`text-sm px-1.5 py-0.5 rounded ${
                              s.source === 'alphafold'
                                ? 'bg-purple-100 text-purple-700'
                                : 'bg-gray-100 text-gray-500'
                            }`}>
                              {s.source === 'alphafold' ? 'AlphaFold' : (s.method || 'Experimental')}
                            </span>
                          </div>
                          {s.resolution && (
                            <span className="text-sm text-gray-500">
                              {s.resolution.toFixed(2)} \u00C5
                            </span>
                          )}
                          {s.source === 'alphafold' && s.confidence && (
                            <span className="text-sm text-purple-500">
                              pLDDT: {s.confidence.toFixed(1)}
                            </span>
                          )}
                        </div>
                        {s.ligand_name && (
                          <p className="text-sm text-amber-600 mt-1">Co-crystal: {s.ligand_name}</p>
                        )}
                        {s.recommended && (
                          <span className="text-sm text-green-600 font-medium">Recommended</span>
                        )}
                        {s.source === 'alphafold' && (
                          <p className="text-sm text-purple-500 mt-1">Full-length predicted structure (may look different from PDB which only covers crystallized domains)</p>
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
                    <span className="text-sm text-gray-400 font-normal flex items-center gap-1">
                      <BindXLogo variant="loading" size={12} />
                      Detecting pockets...
                    </span>
                  ) : (
                    <span className="text-sm text-gray-400 font-normal">
                      ({pocketsForStructure.length} found via P2Rank)
                    </span>
                  )}
                </h3>
                {pocketsLoading ? (
                  <div className="flex flex-col items-center py-6 bg-gray-50 rounded-lg">
                    <BindXLogo variant="loading" size={48} label="Running P2Rank on this structure..." />
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
                      <p className="text-sm text-green-600">Active compounds</p>
                    </div>
                    <div className="text-center">
                      <p className="text-xl font-bold text-green-700">
                        {previewData.chembl_info.n_with_ic50?.toLocaleString()}
                      </p>
                      <p className="text-sm text-green-600">With IC50 data</p>
                    </div>
                  </div>
                  <p className="text-sm text-green-600 mt-2">
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
                  : 'bg-bx-mint text-white hover:bg-bx-mint-dim'
              }`}
            >
              {saving ? (
                <>
                  <BindXLogo variant="loading" size={16} />
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
