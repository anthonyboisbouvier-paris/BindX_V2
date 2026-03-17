import React, { useState, useCallback, useEffect, useMemo } from 'react'
import { useParams, useNavigate } from 'react-router-dom'
import { useWorkspace } from '../contexts/WorkspaceContext.jsx'
import { useToast } from '../contexts/ToastContext.jsx'
import { previewTarget, detectPockets, v9PredictStructure } from '../api.js'
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
  const seqValue = data.sequence?.value || ''

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

  // 1-letter → 3-letter amino acid code mapping
  const AA_MAP = {
    A:'ALA', R:'ARG', N:'ASN', D:'ASP', C:'CYS', E:'GLU', Q:'GLN', G:'GLY',
    H:'HIS', I:'ILE', L:'LEU', K:'LYS', M:'MET', F:'PHE', P:'PRO', S:'SER',
    T:'THR', W:'TRP', Y:'TYR', V:'VAL',
  }
  const getAA = (pos) => {
    if (!seqValue || pos < 1 || pos > seqValue.length) return null
    return AA_MAP[seqValue[pos - 1]] || null
  }

  // Extract individual functional residues for interaction analysis
  const functionalResidues = []
  const seen = new Set()

  // Active sites → importance: "key"
  activeSites.forEach(as => {
    if (as.position && !seen.has(as.position)) {
      seen.add(as.position)
      const aa = getAA(as.position)
      functionalResidues.push({
        number: as.position, aa, role: 'active_site',
        importance: 'key', description: 'Active site',
        source: 'uniprot', enabled: true,
      })
    }
  })

  // Binding sites (ranges) → each position, importance: "key"
  // Skip ranges > 20 residues (likely entire domains, not specific binding sites)
  bindingSites.forEach(bs => {
    if (!bs.start || !bs.end) return
    if (bs.end - bs.start > 20) return
    for (let pos = bs.start; pos <= bs.end; pos++) {
      if (!seen.has(pos)) {
        seen.add(pos)
        const aa = getAA(pos)
        functionalResidues.push({
          number: pos, aa, role: 'binding',
          importance: 'key',
          description: `Binding site (${bs.start}-${bs.end})`,
          source: 'uniprot', enabled: true,
        })
      }
    }
  })

  // Site features (metal binding, nucleotide binding) → importance: "normal"
  // Expand ranges like binding sites (with same >20 safety)
  const siteFeatures = features.filter(f => f.type === 'Site' || f.type === 'Metal binding')
  siteFeatures.forEach(f => {
    const start = f.location?.start?.value
    const end = f.location?.end?.value || start
    if (!start) return
    if (end - start > 20) return
    const role = f.type === 'Metal binding' ? 'metal_binding' : 'site'
    for (let pos = start; pos <= end; pos++) {
      if (!seen.has(pos)) {
        seen.add(pos)
        const aa = getAA(pos)
        functionalResidues.push({
          number: pos, aa, role,
          importance: 'normal',
          description: f.description || f.type,
          source: 'uniprot', enabled: true,
        })
      }
    }
  })

  return {
    name, gene, organism, seqLen, seqValue,
    function: funcText,
    diseases, keywords,
    activeSites, bindingSites, domains,
    functionalResidues,
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
// Functional Residues Editor
// ---------------------------------------------------------------------------

const ROLE_BADGES = {
  active_site:    { bg: 'bg-red-100',    text: 'text-red-700',    label: 'Active Site' },
  binding:        { bg: 'bg-blue-100',   text: 'text-blue-700',   label: 'Binding' },
  metal_binding:  { bg: 'bg-amber-100',  text: 'text-amber-700',  label: 'Metal' },
  site:           { bg: 'bg-gray-100',   text: 'text-gray-600',   label: 'Site' },
  pocket:         { bg: 'bg-purple-100', text: 'text-purple-700', label: 'Pocket' },
  custom:         { bg: 'bg-teal-100',   text: 'text-teal-700',   label: 'Custom' },
}

// 1-letter → 3-letter amino acid code mapping (duplicated for standalone use)
const AA_MAP_EDITOR = {
  A:'ALA', R:'ARG', N:'ASN', D:'ASP', C:'CYS', E:'GLU', Q:'GLN', G:'GLY',
  H:'HIS', I:'ILE', L:'LEU', K:'LYS', M:'MET', F:'PHE', P:'PRO', S:'SER',
  T:'THR', W:'TRP', Y:'TYR', V:'VAL',
}

function FunctionalResiduesEditor({ residues, onChange, pockets, selectedPocketIdx, onSelectResidue, selectedResidueNum, seqValue }) {
  const [showAdd, setShowAdd] = useState(false)
  const [newNumber, setNewNumber] = useState('')
  const [newDescription, setNewDescription] = useState('')
  const [newRole, setNewRole] = useState('custom')

  // Resolve AA code from sequence for custom residues
  const resolveAA = (pos) => {
    if (!seqValue || pos < 1 || pos > seqValue.length) return null
    return AA_MAP_EDITOR[seqValue[pos - 1]] || null
  }

  const enabledCount = residues.filter(r => r.enabled).length

  const handleToggle = (idx) => {
    const updated = [...residues]
    updated[idx] = { ...updated[idx], enabled: !updated[idx].enabled }
    onChange(updated)
  }

  const handleRemove = (idx) => {
    onChange(residues.filter((_, i) => i !== idx))
  }

  const handleAdd = () => {
    const num = parseInt(newNumber, 10)
    if (isNaN(num) || num <= 0) return
    if (residues.some(r => r.number === num)) return
    const aa = resolveAA(num)
    onChange([...residues, {
      number: num, aa, role: newRole,
      importance: newRole === 'active_site' || newRole === 'binding' ? 'key' : 'normal',
      description: newDescription || 'User-defined',
      source: 'user', enabled: true,
    }])
    setNewNumber('')
    setNewDescription('')
    setShowAdd(false)
  }

  const handleImportFromPocket = () => {
    const pocket = pockets?.[selectedPocketIdx]
    if (!pocket) return
    const pocketResidues = Array.isArray(pocket.residues)
      ? (typeof pocket.residues[0] === 'string' && pocket.residues[0].includes(' ')
          ? pocket.residues[0].split(' ')
          : pocket.residues)
      : []
    const existingNums = new Set(residues.map(r => r.number))
    const newResidues = []
    pocketResidues.forEach(r => {
      const str = String(r).trim()
      let aa = null, num = null

      // Format 1: "ALA_A_123" (co-crystallized: resName_chain_resNum)
      const m1 = str.match(/^([A-Z]{3})_[A-Z]_(\d+)$/i)
      if (m1) { aa = m1[1].toUpperCase(); num = parseInt(m1[2], 10) }

      // Format 2: "ALA_123_A" (P2Rank: resName_resNum_chain)
      if (!num) {
        const m2 = str.match(/^([A-Z]{3})_(\d+)_[A-Z]$/i)
        if (m2) { aa = m2[1].toUpperCase(); num = parseInt(m2[2], 10) }
      }

      // Format 3: "ALA123" (compact)
      if (!num) {
        const m3 = str.match(/^([A-Z]{3})(\d+)$/i)
        if (m3) { aa = m3[1].toUpperCase(); num = parseInt(m3[2], 10) }
      }

      // Format 4: just a number "123"
      if (!num) {
        const m4 = str.match(/^(\d+)$/)
        if (m4) { num = parseInt(m4[1], 10) }
      }

      // Fallback: extract any number from the string
      if (!num) {
        const mFall = str.match(/(\d+)/)
        if (mFall) { num = parseInt(mFall[1], 10) }
        const mAA = str.match(/([A-Z]{3})/i)
        if (mAA) { aa = mAA[1].toUpperCase() }
      }

      if (num && !existingNums.has(num)) {
        // If no AA from parsing, try resolving from sequence
        if (!aa) aa = resolveAA(num)
        existingNums.add(num)
        newResidues.push({
          number: num, aa, role: 'pocket',
          importance: 'normal',
          description: `From pocket #${selectedPocketIdx + 1}`,
          source: 'pocket', enabled: true,
        })
      }
    })
    if (newResidues.length > 0) onChange([...residues, ...newResidues])
  }

  return (
    <div className="card px-5 py-4">
      <div className="flex items-center justify-between mb-3">
        <div className="flex items-center gap-2">
          <svg className="w-4 h-4 text-lime-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12l2 2 4-4m5.618-4.016A11.955 11.955 0 0112 2.944a11.955 11.955 0 01-8.618 3.04A12.02 12.02 0 003 9c0 5.591 3.824 10.29 9 11.622 5.176-1.332 9-6.03 9-11.622 0-1.042-.133-2.052-.382-3.016z" />
          </svg>
          <h3 className="text-sm font-semibold text-gray-700">Key Binding Residues</h3>
          <span className="text-xs bg-lime-100 text-lime-700 px-2 py-0.5 rounded-full font-medium">
            {enabledCount} active
          </span>
        </div>
        <div className="flex items-center gap-2">
          {pockets?.length > 0 && (
            <button
              onClick={handleImportFromPocket}
              className="text-xs text-purple-600 hover:text-purple-700 font-medium px-2 py-1 rounded hover:bg-purple-50 transition-colors"
            >
              Import from pocket #{selectedPocketIdx + 1}
            </button>
          )}
          <button
            onClick={() => setShowAdd(!showAdd)}
            className="text-xs text-bx-light-text hover:text-bx-surface font-medium px-2 py-1 rounded hover:bg-blue-50 transition-colors"
          >
            + Add custom
          </button>
        </div>
      </div>

      {/* Add custom residue form */}
      {showAdd && (
        <div className="flex items-end gap-2 mb-3 p-3 bg-gray-50 rounded-lg border border-gray-100">
          <div className="flex-shrink-0">
            <label className="text-[10px] text-gray-400 uppercase tracking-wider block mb-1">Residue #</label>
            <input
              type="number" min={1} value={newNumber}
              onChange={e => setNewNumber(e.target.value)}
              className="w-20 px-2 py-1.5 text-sm rounded border border-gray-200 focus:outline-none focus:ring-1 focus:ring-bx-mint font-mono"
              placeholder="793"
            />
          </div>
          <div className="flex-1">
            <label className="text-[10px] text-gray-400 uppercase tracking-wider block mb-1">Description</label>
            <input
              type="text" value={newDescription}
              onChange={e => setNewDescription(e.target.value)}
              className="w-full px-2 py-1.5 text-sm rounded border border-gray-200 focus:outline-none focus:ring-1 focus:ring-bx-mint"
              placeholder="e.g. Hinge region"
            />
          </div>
          <div className="flex-shrink-0">
            <label className="text-[10px] text-gray-400 uppercase tracking-wider block mb-1">Role</label>
            <select
              value={newRole} onChange={e => setNewRole(e.target.value)}
              className="px-2 py-1.5 text-sm rounded border border-gray-200 focus:outline-none focus:ring-1 focus:ring-bx-mint"
            >
              <option value="active_site">Active Site</option>
              <option value="binding">Binding</option>
              <option value="metal_binding">Metal</option>
              <option value="custom">Custom</option>
            </select>
          </div>
          <button
            onClick={handleAdd} disabled={!newNumber}
            className="px-3 py-1.5 text-sm font-semibold rounded bg-bx-surface text-white hover:bg-bx-elevated disabled:bg-gray-200 disabled:text-gray-400 transition-colors"
          >
            Add
          </button>
        </div>
      )}

      {/* Residues table */}
      {residues.length > 0 ? (
        <div className="max-h-[200px] overflow-y-auto scrollbar-thin">
          <table className="w-full text-sm">
            <thead>
              <tr className="text-[10px] text-gray-400 uppercase tracking-wider border-b border-gray-100">
                <th className="text-left py-1.5 w-8"></th>
                <th className="text-left py-1.5 w-24">Residue</th>
                <th className="text-left py-1.5">Role</th>
                <th className="text-center py-1.5 w-8">Key</th>
                <th className="text-left py-1.5">Description</th>
                <th className="w-8"></th>
              </tr>
            </thead>
            <tbody className="divide-y divide-gray-50">
              {residues.map((r, i) => {
                const badge = ROLE_BADGES[r.role] || ROLE_BADGES.custom
                const isSelected = selectedResidueNum === r.number
                return (
                  <tr
                    key={`${r.number}-${i}`}
                    className={`${r.enabled ? '' : 'opacity-40'} ${isSelected ? 'bg-cyan-50 ring-1 ring-cyan-300' : 'hover:bg-gray-50'} transition-colors cursor-pointer`}
                    onClick={() => onSelectResidue?.(isSelected ? null : r)}
                  >
                    <td className="py-1.5">
                      <input type="checkbox" checked={r.enabled} onChange={(e) => { e.stopPropagation(); handleToggle(i) }}
                        className="w-3.5 h-3.5 rounded border-gray-300 text-bx-mint focus:ring-bx-mint cursor-pointer" />
                    </td>
                    <td className="py-1.5 font-mono font-semibold text-gray-700">
                      {r.aa && <span className="text-gray-400 font-normal mr-1">{r.aa}</span>}
                      {r.number}
                    </td>
                    <td className="py-1.5">
                      <span className={`inline-flex items-center px-1.5 py-0.5 rounded-full text-[10px] font-semibold ${badge.bg} ${badge.text}`}>
                        {badge.label}
                      </span>
                    </td>
                    <td className="py-1.5 text-center">
                      {r.importance === 'key' ? (
                        <svg className="w-3.5 h-3.5 text-amber-400 mx-auto" fill="currentColor" viewBox="0 0 20 20">
                          <path d="M9.049 2.927c.3-.921 1.603-.921 1.902 0l1.07 3.292a1 1 0 00.95.69h3.462c.969 0 1.371 1.24.588 1.81l-2.8 2.034a1 1 0 00-.364 1.118l1.07 3.292c.3.921-.755 1.688-1.54 1.118l-2.8-2.034a1 1 0 00-1.175 0l-2.8 2.034c-.784.57-1.838-.197-1.539-1.118l1.07-3.292a1 1 0 00-.364-1.118L2.98 8.72c-.783-.57-.38-1.81.588-1.81h3.461a1 1 0 00.951-.69l1.07-3.292z" />
                        </svg>
                      ) : (
                        <span className="text-gray-300">-</span>
                      )}
                    </td>
                    <td className="py-1.5 text-xs text-gray-500 truncate max-w-[140px]">{r.description}</td>
                    <td className="py-1.5">
                      {r.source === 'user' || r.source === 'pocket' ? (
                        <button onClick={() => handleRemove(i)} className="text-gray-300 hover:text-red-400 transition-colors">
                          <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
                          </svg>
                        </button>
                      ) : null}
                    </td>
                  </tr>
                )
              })}
            </tbody>
          </table>
        </div>
      ) : (
        <p className="text-xs text-gray-400 text-center py-4">
          No functional residues detected. Add custom residues or import from a pocket.
        </p>
      )}
    </div>
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

  // Functional residues (auto-detected + user-editable)
  const [functionalResidues, setFunctionalResidues] = useState([])
  // Highlighted residue for 3D viewer (click from editor)
  const [highlightedResidue, setHighlightedResidue] = useState(null)

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
      // Restore functional residues from saved data or re-extract from UniProt info
      if (tp.functional_residues?.residues?.length > 0) {
        setFunctionalResidues(tp.functional_residues.residues)
      } else if (tp.uniprot?.functionalResidues?.length > 0) {
        setFunctionalResidues(tp.uniprot.functionalResidues)
      }
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
      setFunctionalResidues(info.functionalResidues || [])

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
          setFunctionalResidues(info.functionalResidues || [])
        } catch (_) { /* non-critical */ }
      }

      // Step 3: Detect pockets via backend
      setLoadingStep('Detecting binding pockets (P2Rank)...')
      const ligandId = pdbInfo.ligands?.[0] || null
      const pocketsResult = await detectPockets(pdbInfo.download_url, pdbInfo.pdb_id, ligandId, pdbInfo.source || '')
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
          pubchem: previewData.pubchem_info,
          structure: selectedStructure,
          structures: previewData.structures || [],
          pockets: currentPockets,
          selected_pocket_index: selectedPocketIdx,
          functional_residues: {
            source: inputMode,
            residues: functionalResidues.filter(r => r.enabled),
            key_hbond_residues: functionalResidues
              .filter(r => r.enabled && r.importance === 'key')
              .map(r => r.number),
            all_residue_numbers: functionalResidues
              .filter(r => r.enabled).map(r => r.number),
            auto_detected: true,
          },
        },
        pockets_detected: currentPockets,
        chembl_actives_count: previewData.chembl_info?.n_actives || null,
        chembl_median_ic50: null,
        pubchem_compounds_count: previewData.pubchem_info?.n_compounds || null,
      })

      await refreshProjects()

      // Auto-trigger receptor preparation (7-step pipeline, project-level)
      try {
        const { v9PrepareReceptor } = await import('../api')
        await v9PrepareReceptor(projectId)
      } catch (prepErr) {
        console.warn('Receptor preparation trigger failed (will retry at docking):', prepErr)
      }

      addToast('Target validated and saved', 'success')
      navigate(`/project/${projectId}`)
    } catch (err) {
      addToast(err.userMessage || err.message || 'Failed to save target', 'error')
    } finally {
      setSaving(false)
    }
  }, [previewData, projectId, selectedStructureIdx, selectedPocketIdx, uniprotId, pdbCode, inputMode, uniprotInfo, functionalResidues, updateProject, refreshProjects, addToast, navigate, pocketsCache])

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
      const result = await detectPockets(url, structure?.label || 'structure', structure?.ligand_id, structure?.source || '')
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
              onClick={async () => {
                try {
                  const res = await fetch('https://rest.uniprot.org/uniprotkb/search?query=(reviewed:true)+AND+(organism_id:9606)+AND+(length:[200+TO+1000])&format=json&size=500&fields=accession')
                  if (!res.ok) throw new Error('UniProt search failed')
                  const data = await res.json()
                  const entries = data.results || []
                  if (entries.length > 0) {
                    const random = entries[Math.floor(Math.random() * entries.length)]
                    setUniprotId(random.primaryAccession)
                  }
                } catch (err) {
                  console.warn('Random UniProt failed:', err)
                  // Fallback: curated list of well-known drug targets
                  const targets = ['P00533','P04637','P10275','P23458','P00519','P15056','P07900','P35354','P08684','P11511','P11309','P38398','P29597','P29274','P14416','Q16539','P35968','P00918','P21802','P09874']
                  setUniprotId(targets[Math.floor(Math.random() * targets.length)])
                }
              }}
              disabled={loading}
              className="px-3 py-2.5 text-sm font-medium rounded-lg transition-colors border border-gray-200 text-gray-500 hover:text-bx-surface hover:border-bx-surface hover:bg-blue-50"
              title="Load a random human protein target"
            >
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
              </svg>
            </button>
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
              onClick={() => {
                const pdbTargets = ['1M17','4HJO','6LU7','3ERT','1OHR','2HYY','3NYA','4MQT','5HT1','6GCZ','1XKK','3PXF','4DKL','5C1M','2RGP','3QKK','1ERE','4ASD','5UIG','3EQM','1YCR','4JXS','2VT4','3K5V','5T35','6HD4','1GOS','4WKQ','2W96','3G0E']
                setPdbCode(pdbTargets[Math.floor(Math.random() * pdbTargets.length)])
              }}
              disabled={loading}
              className="px-3 py-2.5 text-sm font-medium rounded-lg transition-colors border border-gray-200 text-gray-500 hover:text-bx-surface hover:border-bx-surface hover:bg-blue-50"
              title="Load a random PDB structure"
            >
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
              </svg>
            </button>
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
              disabled={loading}
            />
            <div className="flex items-center gap-3">
              <button
                onClick={async () => {
                  try {
                    const targets = ['P00533','P04637','P10275','P23458','P00519','P15056','P07900','P35354','P08684','P11511','P11309','P38398','P29597','P29274','P14416','Q16539','P35968','P00918','P21802','P09874']
                    const acc = targets[Math.floor(Math.random() * targets.length)]
                    const res = await fetch(`https://rest.uniprot.org/uniprotkb/${acc}.fasta`)
                    if (res.ok) {
                      const fasta = await res.text()
                      setFastaInput(fasta.trim())
                    }
                  } catch (err) {
                    console.warn('Random FASTA failed:', err)
                  }
                }}
                disabled={loading}
                className="px-3 py-2.5 text-sm font-medium rounded-lg transition-colors border border-gray-200 text-gray-500 hover:text-bx-surface hover:border-bx-surface hover:bg-blue-50"
                title="Load a random FASTA sequence"
              >
                <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
                </svg>
              </button>
              <button
                onClick={async () => {
                  if (!fastaInput.trim()) return
                  setLoading(true)
                  setError(null)
                  setLoadingStep('Predicting structure with ESMFold...')
                  try {
                    const result = await v9PredictStructure(fastaInput)
                    // Create a blob URL from the PDB text
                    const blob = new Blob([result.pdb_text], { type: 'chemical/x-pdb' })
                    const pdbBlobUrl = URL.createObjectURL(blob)
                    setPreviewData({
                      structures: [{
                        pdb_id: 'esmfold_prediction',
                        source: 'esmfold',
                        method: result.method,
                        resolution: null,
                        download_url: pdbBlobUrl,
                      }],
                      uniprot_id: null,
                    })
                    setSelectedStructureIdx(0)
                    addToast(`Structure predicted: ${result.sequence_length} residues`, 'success')
                  } catch (e) {
                    setError(e.response?.data?.detail || e.message || 'Prediction failed')
                  } finally {
                    setLoading(false)
                    setLoadingStep(null)
                  }
                }}
                disabled={loading || !fastaInput.trim()}
                className={`px-5 py-2.5 text-sm font-semibold rounded-lg transition-colors flex items-center gap-2 ${
                  loading || !fastaInput.trim()
                    ? 'bg-gray-200 text-gray-400 cursor-not-allowed'
                    : 'bg-purple-600 text-white hover:bg-purple-700'
                }`}
              >
                {loading ? (
                  <>
                    <BindXLogo variant="loading" size={16} />
                    Predicting...
                  </>
                ) : 'Predict with ESMFold'}
              </button>
              <span className="text-xs text-gray-400">~1-3 min for typical proteins</span>
            </div>
            <div className="bg-blue-50 border border-blue-100 rounded-lg px-4 py-3 space-y-2">
              <p className="text-sm text-blue-700 font-medium">Structure prediction from FASTA sequence</p>
              <div className="grid grid-cols-1 sm:grid-cols-2 gap-2 text-sm text-blue-600">
                <div className="flex items-start gap-2">
                  <span className="text-blue-400 mt-0.5">&#9679;</span>
                  <span><strong>ESMFold</strong> — Fast prediction (~2 min). Uses Meta's ESMFold API. Max 2000 residues.</span>
                </div>
                <div className="flex items-start gap-2">
                  <span className="text-blue-400 mt-0.5">&#9679;</span>
                  <span><strong>AlphaFold</strong> — Pre-computed structures are auto-fetched when using UniProt ID input.</span>
                </div>
              </div>
            </div>
          </div>
        )}

        {/* Loading step indicator — prominent centered display */}
        {loading && loadingStep && (
          <div className="mt-8 flex flex-col items-center justify-center gap-4 py-12">
            <BindXLogo variant="loading" size={80} />
            <p className="text-sm font-medium text-gray-600 animate-pulse">{loadingStep}</p>
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
                {/* ChEMBL + PubChem — inline in protein header */}
                {(previewData?.chembl_info?.has_data || previewData?.pubchem_info?.has_data) && (
                  <div>
                    <p className="text-sm font-semibold text-gray-500 uppercase tracking-wider mb-1">Known Activity Data</p>
                  <div className="flex flex-wrap gap-3">
                    {previewData.chembl_info?.has_data && (
                      <div className="flex items-center gap-2 bg-green-50 border border-green-100 rounded-lg px-3 py-2">
                        <div className="w-2 h-2 rounded-full bg-green-500 flex-shrink-0" />
                        <div className="text-sm">
                          <span className="font-bold text-green-700">{previewData.chembl_info.n_actives?.toLocaleString()}</span>
                          <span className="text-green-600 ml-1">ChEMBL actives</span>
                          {previewData.chembl_info.n_with_ic50 > 0 && (
                            <span className="text-green-500 ml-1">({previewData.chembl_info.n_with_ic50?.toLocaleString()} IC50)</span>
                          )}
                        </div>
                      </div>
                    )}
                    {previewData.pubchem_info?.has_data && (
                      <div className="flex items-center gap-2 bg-blue-50 border border-blue-100 rounded-lg px-3 py-2">
                        <div className="w-2 h-2 rounded-full bg-blue-500 flex-shrink-0" />
                        <div className="text-sm">
                          <span className="font-bold text-blue-700">{previewData.pubchem_info.n_compounds?.toLocaleString()}</span>
                          <span className="text-blue-600 ml-1">PubChem compounds</span>
                        </div>
                      </div>
                    )}
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
                  highlightedResidue={highlightedResidue}
                  functionalResidues={functionalResidues}
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
              {/* Structure selection — show available structures */}
              {previewData.structures?.length > 0 && (
                <div>
                  <h3 className="text-sm font-semibold text-gray-700 mb-2">
                    Available Structures
                    <span className="ml-1.5 text-xs font-normal text-gray-400">
                      ({previewData.structures.length})
                    </span>
                  </h3>
                  <div className="space-y-2 max-h-[220px] overflow-y-auto pr-1 scrollbar-thin">
                    {previewData.structures.map((s, i) => {
                      const realIdx = i
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
                                : s.source === 'esmfold'
                                ? 'bg-green-100 text-green-700'
                                : 'bg-gray-100 text-gray-500'
                            }`}>
                              {s.source === 'alphafold' ? 'AlphaFold' : s.source === 'esmfold' ? 'ESMFold' : (s.method || 'Experimental')}
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
                  <div className="space-y-2 max-h-[220px] overflow-y-auto pr-1 scrollbar-thin">
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

              {/* Functional Residues Editor — below pockets */}
              <FunctionalResiduesEditor
                residues={functionalResidues}
                onChange={setFunctionalResidues}
                pockets={pocketsForStructure}
                selectedPocketIdx={selectedPocketIdx}
                onSelectResidue={(r) => setHighlightedResidue(r)}
                selectedResidueNum={highlightedResidue?.number || null}
                seqValue={uniprotInfo?.seqValue || ''}
              />

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
