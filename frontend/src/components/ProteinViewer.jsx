import React, { useEffect, useRef, useState, useCallback } from 'react'
import BindXLogo from './BindXLogo.jsx'

// 3Dmol.js CDN loader
function load3Dmol() {
  return new Promise((resolve, reject) => {
    if (window.$3Dmol) { resolve(window.$3Dmol); return }
    const existing = document.getElementById('3dmol-script')
    if (existing) {
      existing.addEventListener('load', () => resolve(window.$3Dmol))
      existing.addEventListener('error', reject)
      return
    }
    const script = document.createElement('script')
    script.id = '3dmol-script'
    script.src = 'https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js'
    script.onload = () => window.$3Dmol ? resolve(window.$3Dmol) : reject(new Error('3Dmol unavailable'))
    script.onerror = () => reject(new Error('Failed to load 3Dmol.js'))
    document.head.appendChild(script)
  })
}

// ---------------------------------------------------------------------------
// Rendering styles
// ---------------------------------------------------------------------------
const STYLES = ['cartoon', 'surface', 'tube', 'ball+stick', 'stick', 'line', 'sphere']
const STYLE_LABELS = {
  cartoon: 'Cartoon', surface: 'Surface', tube: 'Tube',
  'ball+stick': 'Ball+Stick', stick: 'Stick', line: 'Line', sphere: 'Sphere',
}

// Pocket color presets
const POCKET_COLORS = [
  { value: '#f59e0b', label: 'Amber' },
  { value: '#ef4444', label: 'Red' },
  { value: '#00e6a0', label: 'Green' },
  { value: '#3b82f6', label: 'Blue' },
  { value: '#a855f7', label: 'Purple' },
  { value: '#ec4899', label: 'Pink' },
  { value: '#06b6d4', label: 'Cyan' },
  { value: '#ffffff', label: 'White' },
]

// ---------------------------------------------------------------------------
// Color schemes — organized by category
// ---------------------------------------------------------------------------
const COLOR_GROUPS = [
  {
    label: 'Protein',
    options: [
      { value: 'default',    label: 'Blue (Default)' },
      { value: 'ss',         label: 'Secondary Structure' },
      { value: 'chain',      label: 'Chain' },
      { value: 'spectrum',   label: 'Rainbow (N→C)' },
      { value: 'bfactor',    label: 'B-factor / pLDDT' },
      { value: 'hydro',      label: 'Hydrophobicity' },
    ],
  },
  {
    label: 'Residue',
    options: [
      { value: 'amino',      label: 'Amino Acid Type' },
      { value: 'shapely',    label: 'Shapely' },
      { value: 'charge',     label: 'Charge (+/-/0)' },
    ],
  },
  {
    label: 'Atom',
    options: [
      { value: 'Jmol',       label: 'Jmol (CPK)' },
      { value: 'rasmol',     label: 'RasMol' },
    ],
  },
  {
    label: 'Custom',
    options: [
      { value: 'custom',     label: 'Pick Color...' },
    ],
  },
]

// Background presets
const BG_PRESETS = [
  { value: '#0f1923', label: 'Dark Navy' },
  { value: '#000000', label: 'Black' },
  { value: '#1a1a2e', label: 'Deep Blue' },
  { value: '#ffffff', label: 'White' },
  { value: '#f0f0f0', label: 'Light Gray' },
]

// Hydrophobicity scale (Kyte-Doolittle) mapped to color
const HYDRO_SCALE = {
  ILE: 4.5, VAL: 4.2, LEU: 3.8, PHE: 2.8, CYS: 2.5,
  MET: 1.9, ALA: 1.8, GLY: -0.4, THR: -0.7, SER: -0.8,
  TRP: -0.9, TYR: -1.3, PRO: -1.6, HIS: -3.2, GLU: -3.5,
  GLN: -3.5, ASP: -3.5, ASN: -3.5, LYS: -3.9, ARG: -4.5,
}

function getColorSpec(mode, customColor) {
  switch (mode) {
    case 'ss':       return { colorscheme: 'ssJmol' }
    case 'chain':    return { colorscheme: 'chain' }
    case 'spectrum': return { color: 'spectrum' }
    case 'amino':    return { colorscheme: 'amino' }
    case 'shapely':  return { colorscheme: 'shapely' }
    case 'Jmol':     return { colorscheme: 'Jmol' }
    case 'rasmol':   return { colorscheme: 'rasmol' }
    case 'bfactor':  return { colorscheme: { prop: 'b', gradient: 'rwb', min: 0, max: 100 } }
    case 'custom':   return { color: customColor || '#4a9eff' }
    case 'hydro':    return null
    case 'charge':   return null
    default:         return { color: '#4a9eff' }
  }
}

function hydroColor(resName) {
  const val = HYDRO_SCALE[resName] ?? 0
  const norm = (val + 4.5) / 9.0
  const r = Math.round(norm > 0.5 ? 255 : 255 * (norm * 2))
  const b = Math.round(norm < 0.5 ? 255 : 255 * ((1 - norm) * 2))
  const g = Math.round(norm > 0.5 ? 255 * ((1 - norm) * 2) : 255 * (norm * 2))
  return `rgb(${r},${g},${b})`
}

const CHARGE_MAP = {
  ARG: '#3b82f6', LYS: '#3b82f6', HIS: '#93c5fd',
  ASP: '#ef4444', GLU: '#ef4444',
}

// Build a style spec for a given rendering style + colorSpec
function buildStyleSpec(renderStyle, colorSpec) {
  switch (renderStyle) {
    case 'surface':
      return { cartoon: { ...(colorSpec || {}), opacity: 0.3 } }
    case 'ball+stick':
      return {
        stick: { ...(colorSpec || {}), radius: 0.12 },
        sphere: { ...(colorSpec || {}), radius: 0.3 },
      }
    case 'stick':
      return { stick: { ...(colorSpec || {}), radius: 0.15 } }
    case 'line':
      return { line: { ...(colorSpec || {}) } }
    case 'sphere':
      return { sphere: { ...(colorSpec || {}) } }
    case 'tube':
      return { cartoon: { ...(colorSpec || {}), tubes: true, thickness: 0.4 } }
    default: // cartoon
      return { cartoon: { ...(colorSpec || {}) } }
  }
}

// SVG icon components
function CameraIcon() {
  return (
    <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M3 9a2 2 0 012-2h.93a2 2 0 001.664-.89l.812-1.22A2 2 0 0110.07 4h3.86a2 2 0 011.664.89l.812 1.22A2 2 0 0018.07 7H19a2 2 0 012 2v9a2 2 0 01-2 2H5a2 2 0 01-2-2V9z" />
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 13a3 3 0 11-6 0 3 3 0 016 0z" />
    </svg>
  )
}

function CrosshairIcon() {
  return (
    <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <circle cx="12" cy="12" r="10" strokeWidth={2} />
      <path strokeLinecap="round" strokeWidth={2} d="M12 2v4m0 12v4M2 12h4m12 0h4" />
    </svg>
  )
}

function TagIcon() {
  return (
    <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M7 7h.01M7 3h5c.512 0 1.024.195 1.414.586l7 7a2 2 0 010 2.828l-7 7a2 2 0 01-2.828 0l-7-7A2 2 0 013 12V7a4 4 0 014-4z" />
    </svg>
  )
}

function RulerIcon() {
  return (
    <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M3 21l9-9m0 0l9-9m-9 9l-3-3m3 3l3-3m-3 3l-3 3m3-3l3 3" />
    </svg>
  )
}

function ScissorsIcon() {
  return (
    <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 6h16M4 12h16M4 18h16" />
    </svg>
  )
}

// Zone color presets
const ZONE_COLORS = ['#10b981', '#f43f5e', '#a855f7', '#06b6d4', '#f59e0b', '#6366f1', '#ec4899', '#14b8a6']

// Parse residue input: "745, 790, 855-870" → sorted unique array of ints
function parseResidueInput(input) {
  const residues = []
  const parts = input.split(/[,;\s]+/).filter(Boolean)
  parts.forEach(part => {
    const range = part.match(/^(\d+)\s*-\s*(\d+)$/)
    if (range) {
      const start = parseInt(range[1]), end = parseInt(range[2])
      for (let i = Math.min(start, end); i <= Math.max(start, end); i++) residues.push(i)
    } else {
      const n = parseInt(part)
      if (!isNaN(n)) residues.push(n)
    }
  })
  return [...new Set(residues)].sort((a, b) => a - b)
}

// Compact residue array to string: [1,2,3,5,7,8,9] → "1-3, 5, 7-9"
function residuesToString(res) {
  if (!res || res.length === 0) return ''
  const parts = []
  let i = 0
  while (i < res.length) {
    let j = i
    while (j + 1 < res.length && res[j + 1] === res[j] + 1) j++
    parts.push(j > i ? `${res[i]}-${res[j]}` : `${res[i]}`)
    i = j + 1
  }
  return parts.join(', ')
}

// Common solvent/ion/sugar residue names to exclude from ligand display
const SOLVENT_RESNAMES = new Set([
  // Water
  'HOH','WAT','H2O','DOD','SOL',
  // Common crystallization additives & buffers
  'SO4','PO4','GOL','EDO','ACE','PEG','DMS','BME','FMT','ACT',
  'TRS','MPD','PG4','EPE','MES','CIT','TAR','SUC','MLI','IMD',
  // Ions
  'CL','NA','MG','ZN','CA','K','MN','FE','CU','CO','NI',
  'IOD','BR','NH4','CD','HG','PT','AU','SR','BA','CS','RB','LI',
  // Glycosylation sugars (very common in crystal structures)
  'NAG','MAN','GAL','BMA','FUC','BGC','NDG','SIA','GLC','XYP',
  'A2G','NGA','RAM','RIB','BDP','GLA','ARA',
  // Modified residues / amino acid variants
  'SEP','TPO','PTR','MSE','CSO','OCS','KCX','LLP','CME','CSD',
  // Nucleotides & cofactors (not drug-like)
  'UNL',
])

// ---------------------------------------------------------------------------
// ZoneCard — isolated sub-component to avoid re-rendering the entire viewer on keystroke
// ---------------------------------------------------------------------------
function ZoneCard({ zone, isDarkBg, isEditing, onUpdate, onToggleVisibility, onToggleEdit, onDelete }) {
  const [localInput, setLocalInput] = useState(residuesToString(zone.residues))
  const prevResiduesRef = useRef(zone.residues)

  // Sync local input when residues change externally (click-to-add)
  useEffect(() => {
    if (prevResiduesRef.current !== zone.residues) {
      setLocalInput(residuesToString(zone.residues))
      prevResiduesRef.current = zone.residues
    }
  }, [zone.residues])

  const commitResidues = () => {
    const parsed = parseResidueInput(localInput)
    onUpdate(zone.id, { residues: parsed })
  }

  return (
    <div className={`rounded-lg border p-2 space-y-1.5 ${
      isEditing
        ? isDarkBg ? 'border-emerald-500/50 bg-emerald-900/20' : 'border-emerald-400 bg-emerald-50'
        : isDarkBg ? 'border-gray-600' : 'border-gray-200'
    }`}>
      {/* Zone header */}
      <div className="flex items-center gap-1.5">
        <input
          type="color"
          value={zone.color}
          onChange={(e) => onUpdate(zone.id, { color: e.target.value })}
          className="w-5 h-5 rounded border-0 cursor-pointer flex-shrink-0"
          title="Zone color"
        />
        <input
          type="text"
          value={zone.name}
          onChange={(e) => onUpdate(zone.id, { name: e.target.value })}
          className={`flex-1 text-xs px-1.5 py-0.5 rounded border bg-transparent outline-none ${
            isDarkBg ? 'border-gray-600 text-gray-200 focus:border-gray-400' : 'border-gray-300 text-gray-700 focus:border-gray-500'
          }`}
          placeholder="Zone name"
        />
        <button
          onClick={() => onToggleVisibility(zone.id)}
          title={zone.visible ? 'Hide zone' : 'Show zone'}
          className={`p-0.5 rounded text-xs ${
            zone.visible
              ? isDarkBg ? 'text-emerald-400 hover:bg-white/10' : 'text-emerald-600 hover:bg-gray-100'
              : isDarkBg ? 'text-gray-600 hover:bg-white/10' : 'text-gray-400 hover:bg-gray-100'
          }`}
        >
          <svg className="w-3.5 h-3.5" fill={zone.visible ? 'currentColor' : 'none'} stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 12a3 3 0 11-6 0 3 3 0 016 0z" />
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M2.458 12C3.732 7.943 7.523 5 12 5c4.478 0 8.268 2.943 9.542 7-1.274 4.057-5.064 7-9.542 7-4.477 0-8.268-2.943-9.542-7z" />
          </svg>
        </button>
        <button
          onClick={() => onToggleEdit(zone.id)}
          title={isEditing ? 'Stop click-to-add' : 'Click residues on structure to add/remove'}
          className={`p-0.5 rounded text-xs ${
            isEditing
              ? 'bg-emerald-600 text-white'
              : isDarkBg ? 'text-gray-400 hover:text-white hover:bg-white/10' : 'text-gray-500 hover:bg-gray-100'
          }`}
        >
          <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15.232 5.232l3.536 3.536m-2.036-5.036a2.5 2.5 0 113.536 3.536L6.5 21.036H3v-3.572L16.732 3.732z" />
          </svg>
        </button>
        <button
          onClick={() => onDelete(zone.id)}
          title="Delete zone"
          className={`p-0.5 rounded text-xs ${isDarkBg ? 'text-red-400 hover:bg-red-900/30' : 'text-red-500 hover:bg-red-50'}`}
        >
          <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 7l-.867 12.142A2 2 0 0116.138 21H7.862a2 2 0 01-1.995-1.858L5 7m5 4v6m4-6v6m1-10V4a1 1 0 00-1-1h-4a1 1 0 00-1 1v3M4 7h16" />
          </svg>
        </button>
      </div>

      {/* Residue text input — local state, commits on blur/Enter */}
      <input
        type="text"
        value={localInput}
        onChange={(e) => setLocalInput(e.target.value)}
        onBlur={commitResidues}
        onKeyDown={(e) => { if (e.key === 'Enter') commitResidues() }}
        className={`w-full text-xs px-1.5 py-1 rounded border bg-transparent outline-none font-mono ${
          isDarkBg ? 'border-gray-600 text-gray-300 focus:border-emerald-500 placeholder-gray-600' : 'border-gray-300 text-gray-700 focus:border-emerald-400 placeholder-gray-400'
        }`}
        placeholder="745, 790, 855-870"
      />

      {/* Residue count */}
      <div className={`text-xs ${isDarkBg ? 'text-gray-500' : 'text-gray-400'}`}>
        {zone.residues.length} residue{zone.residues.length !== 1 ? 's' : ''}
        {isEditing && (
          <span className="ml-2 text-emerald-400 animate-pulse">Click on structure to add/remove</span>
        )}
      </div>
    </div>
  )
}

// Props: pdbUrl, pdbData, selectedPocket, uniprotFeatures, height, onDiagnostic
export default function ProteinViewer({ pdbUrl, pdbData, selectedPocket, uniprotFeatures, height = 480, onDiagnostic }) {
  const wrapperRef = useRef(null)
  const containerRef = useRef(null)
  const viewerRef = useRef(null)
  const settingsRef = useRef(null)
  const settingsBtnRef = useRef(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState(null)
  const [ready, setReady] = useState(false)
  const [style, setStyle] = useState('cartoon')
  const [colorMode, setColorMode] = useState('ss')
  const [customColor, setCustomColor] = useState('#4a9eff')
  const [bgColor, setBgColor] = useState('#0f1923')
  const [pocketColor, setPocketColor] = useState('#f59e0b')
  const [ligandStyle, setLigandStyle] = useState('ball+stick') // 'ball+stick' | 'stick' | 'sphere' | 'line' | 'off'
  const [pocketRadius, setPocketRadius] = useState(5)
  const [pocketOpacity, setPocketOpacity] = useState(0.75)
  const [pocketStyle, setPocketStyle] = useState('surface')
  const [spinning, setSpinning] = useState(false)
  const [showDomains, setShowDomains] = useState(false)
  const [showActiveSites, setShowActiveSites] = useState(true)
  const [showBindingSites, setShowBindingSites] = useState(true)
  const [isFullscreen, setIsFullscreen] = useState(false)
  const [showSettings, setShowSettings] = useState(false)

  // Expert tools state
  const [labelMode, setLabelMode] = useState(false) // click → show residue label
  const [distanceMode, setDistanceMode] = useState(false) // click 2 atoms → distance
  const [distanceAtom1, setDistanceAtom1] = useState(null) // first atom for distance measurement
  const [selectedResidues, setSelectedResidues] = useState([]) // multi-select
  const [showSlab, setShowSlab] = useState(false)
  const [slabNear, setSlabNear] = useState(-50)
  const [slabFar, setSlabFar] = useState(50)

  // Custom zones state
  const [customZones, setCustomZones] = useState([]) // [{ id, name, color, residues: number[], visible: true }]
  const [showZonesPanel, setShowZonesPanel] = useState(false)
  const [zoneEditMode, setZoneEditMode] = useState(null) // zone id being click-edited, or null
  const zoneEditModeRef = useRef(null)
  const customZonesRef = useRef([])
  const zonesPanelRef = useRef(null)
  const zonesBtnRef = useRef(null)
  const nextZoneId = useRef(1)

  // 3Dmol library reference (needed for SurfaceType constants in addSurface)
  const mol3dRef = useRef(null)

  // Stable ref for onDiagnostic callback — avoids useEffect re-trigger on parent re-render
  const onDiagnosticRef = useRef(onDiagnostic)
  useEffect(() => { onDiagnosticRef.current = onDiagnostic }, [onDiagnostic])

  // Track labels and distance shapes for cleanup
  const labelsRef = useRef([])
  const distanceLinesRef = useRef([])

  // Refs for click handler — avoids stale closure in 3Dmol.js setClickable
  const labelModeRef = useRef(false)
  const distanceModeRef = useRef(false)
  const distanceAtom1Ref = useRef(null)

  // ---------------------------------------------------------------------------
  // Click-outside handler for settings panel
  // ---------------------------------------------------------------------------
  useEffect(() => {
    if (!showSettings) return
    const handler = (e) => {
      if (
        settingsRef.current && !settingsRef.current.contains(e.target) &&
        settingsBtnRef.current && !settingsBtnRef.current.contains(e.target)
      ) {
        setShowSettings(false)
      }
    }
    document.addEventListener('mousedown', handler)
    return () => document.removeEventListener('mousedown', handler)
  }, [showSettings])

  // Click-outside handler for zones panel
  useEffect(() => {
    if (!showZonesPanel) return
    const handler = (e) => {
      if (
        zonesPanelRef.current && !zonesPanelRef.current.contains(e.target) &&
        zonesBtnRef.current && !zonesBtnRef.current.contains(e.target)
      ) {
        setShowZonesPanel(false)
      }
    }
    document.addEventListener('mousedown', handler)
    return () => document.removeEventListener('mousedown', handler)
  }, [showZonesPanel])

  // Keep refs in sync for click handler
  useEffect(() => { zoneEditModeRef.current = zoneEditMode }, [zoneEditMode])
  useEffect(() => { customZonesRef.current = customZones }, [customZones])

  // ---------------------------------------------------------------------------
  // Load PDB data
  // ---------------------------------------------------------------------------
  const initViewer = useCallback(async () => {
    if (!containerRef.current) return
    if (!pdbUrl && !pdbData) return

    setLoading(true)
    setError(null)
    setReady(false)

    try {
      const $3Dmol = await load3Dmol()
      mol3dRef.current = $3Dmol

      if (viewerRef.current) {
        try { viewerRef.current.clear() } catch (_) {}
      }
      containerRef.current.innerHTML = ''

      const viewer = $3Dmol.createViewer(containerRef.current, {
        backgroundColor: bgColor,
      })
      viewerRef.current = viewer

      let data = pdbData
      if (!data && pdbUrl) {
        const res = await fetch(pdbUrl)
        if (!res.ok) throw new Error(`Failed to fetch PDB: ${res.status}`)
        data = await res.text()
      }

      viewer.addModel(data, 'pdb')
      viewer.setStyle({}, { cartoon: { colorscheme: 'ssJmol' } })
      viewer.zoomTo()
      viewer.zoom(0.85)
      viewer.render()

      // Reset expert tool state
      labelsRef.current = []
      distanceLinesRef.current = []
      setDistanceAtom1(null)
      setSelectedResidues([])

      setReady(true)
      setLoading(false)
    } catch (err) {
      setError(err.message)
      setLoading(false)
    }
  }, [pdbUrl, pdbData])

  useEffect(() => { initViewer() }, [initViewer])

  // ---------------------------------------------------------------------------
  // Background color change
  // ---------------------------------------------------------------------------
  useEffect(() => {
    const viewer = viewerRef.current
    if (!viewer || !ready) return
    viewer.setBackgroundColor(bgColor)
    viewer.render()
  }, [bgColor, ready])

  // ---------------------------------------------------------------------------
  // Spin toggle
  // ---------------------------------------------------------------------------
  useEffect(() => {
    const viewer = viewerRef.current
    if (!viewer || !ready) return
    if (spinning) {
      viewer.spin('y', 1)
    } else {
      viewer.spin(false)
    }
  }, [spinning, ready])

  // ---------------------------------------------------------------------------
  // Slab / Clip
  // ---------------------------------------------------------------------------
  useEffect(() => {
    const viewer = viewerRef.current
    if (!viewer || !ready) return
    if (showSlab) {
      viewer.setSlab(slabNear, slabFar)
    } else {
      viewer.setSlab(-999, 999)
    }
    viewer.render()
  }, [showSlab, slabNear, slabFar, ready])

  // ---------------------------------------------------------------------------
  // Keep refs in sync with state (for click handler closure)
  // ---------------------------------------------------------------------------
  useEffect(() => { labelModeRef.current = labelMode }, [labelMode])
  useEffect(() => { distanceModeRef.current = distanceMode }, [distanceMode])
  useEffect(() => { distanceAtom1Ref.current = distanceAtom1 }, [distanceAtom1])

  // ---------------------------------------------------------------------------
  // Click handler for labels, distance, multi-select
  // Uses refs to avoid stale closure — setClickable is called ONCE
  // ---------------------------------------------------------------------------
  useEffect(() => {
    const viewer = viewerRef.current
    if (!viewer || !ready) return

    const handleClick = (atom) => {
      if (!atom) return

      const resLabel = `${atom.resn} ${atom.resi}${atom.chain ? ':' + atom.chain : ''}`

      // Label mode (read from ref, not state)
      if (labelModeRef.current) {
        viewer.addLabel(resLabel, {
          position: { x: atom.x, y: atom.y, z: atom.z },
          fontSize: 12,
          fontColor: 'white',
          backgroundColor: 'rgba(0,0,0,0.7)',
          backgroundOpacity: 0.7,
          borderRadius: 4,
          padding: 4,
        })
        viewer.render()
        return
      }

      // Distance measurement mode (read from refs)
      if (distanceModeRef.current) {
        const atom1 = distanceAtom1Ref.current
        if (!atom1) {
          // First atom — store in ref AND state
          const atomData = { x: atom.x, y: atom.y, z: atom.z, label: resLabel, atom: atom.atom }
          distanceAtom1Ref.current = atomData
          setDistanceAtom1(atomData)
          viewer.addSphere({
            center: { x: atom.x, y: atom.y, z: atom.z },
            radius: 0.5,
            color: '#00e6a0',
            opacity: 0.8,
          })
          viewer.addLabel(`${resLabel} (${atom.atom})`, {
            position: { x: atom.x, y: atom.y, z: atom.z },
            fontSize: 11,
            fontColor: '#00e6a0',
            backgroundColor: 'rgba(0,0,0,0.7)',
            backgroundOpacity: 0.7,
          })
          viewer.render()
        } else {
          // Second atom — compute distance and draw
          const dx = atom.x - atom1.x
          const dy = atom.y - atom1.y
          const dz = atom.z - atom1.z
          const dist = Math.sqrt(dx * dx + dy * dy + dz * dz).toFixed(1)

          viewer.addCylinder({
            start: { x: atom1.x, y: atom1.y, z: atom1.z },
            end: { x: atom.x, y: atom.y, z: atom.z },
            radius: 0.07,
            color: '#00e6a0',
            fromCap: true,
            toCap: true,
          })

          const mid = {
            x: (atom.x + atom1.x) / 2,
            y: (atom.y + atom1.y) / 2,
            z: (atom.z + atom1.z) / 2,
          }
          viewer.addLabel(`${dist} \u00C5`, {
            position: mid,
            fontSize: 13,
            fontColor: '#00e6a0',
            backgroundColor: 'rgba(0,0,0,0.8)',
            backgroundOpacity: 0.8,
            borderRadius: 4,
            padding: 5,
          })

          viewer.addSphere({
            center: { x: atom.x, y: atom.y, z: atom.z },
            radius: 0.5,
            color: '#00e6a0',
            opacity: 0.8,
          })
          viewer.addLabel(`${resLabel} (${atom.atom})`, {
            position: { x: atom.x, y: atom.y, z: atom.z },
            fontSize: 11,
            fontColor: '#00e6a0',
            backgroundColor: 'rgba(0,0,0,0.7)',
            backgroundOpacity: 0.7,
          })

          distanceLinesRef.current.push({ start: atom1, end: atom, dist })
          distanceAtom1Ref.current = null
          setDistanceAtom1(null)
          viewer.render()
        }
        return
      }

      // Zone edit mode — click to add/remove residues from the active zone
      if (zoneEditModeRef.current) {
        const zoneId = zoneEditModeRef.current
        const resi = atom.resi
        setCustomZones(prev => prev.map(z => {
          if (z.id !== zoneId) return z
          const has = z.residues.includes(resi)
          return { ...z, residues: has ? z.residues.filter(r => r !== resi) : [...z.residues, resi].sort((a, b) => a - b) }
        }))
        return
      }

      // Multi-select residue mode
      const resi = atom.resi
      const chain = atom.chain || ''
      const key = `${chain}_${resi}`
      setSelectedResidues(prev => {
        const exists = prev.find(r => r.key === key)
        if (exists) return prev.filter(r => r.key !== key)
        return [...prev, { key, resi, chain, resn: atom.resn }]
      })
    }

    viewer.setClickable({}, true, handleClick)
    return () => {
      try { viewer.setClickable({}, false) } catch (_) {}
    }
  }, [ready]) // Only depends on ready — refs handle the rest

  // ---------------------------------------------------------------------------
  // Re-apply styles when anything changes
  // ---------------------------------------------------------------------------
  useEffect(() => {
    const viewer = viewerRef.current
    if (!viewer || !ready) return

    const models = viewer.getModelList()
    const protein = models?.[0]
    if (!protein) return

    // Diagnostic data — collected during rendering
    const diag = {
      structureAtoms: 0,
      ligandAtomsFound: 0,
      pocketResiduesRaw: 0,
      pocketResiduesParsed: 0,
      pocketAtomsMatched: 0,
      pocketSurfaceAttempted: false,
      pocketRenderMode: pocketStyle,
      activeSitesResolved: 0,
      bindingSitesResolved: 0,
      domainsResolved: 0,
      warnings: [],
    }

    // Count total structure atoms
    try { diag.structureAtoms = viewer.selectedAtoms({ model: protein }).length } catch (_) {}

    try { viewer.removeAllSurfaces() } catch (_) {}
    try { viewer.removeAllShapes() } catch (_) {}
    try { viewer.removeAllLabels() } catch (_) {}

    // Helper: addSurface is async in 3Dmol.js — pass callback to re-render when ready
    // Signature: addSurface(type, data, atomsel, allsel, focus, surfacecallback)
    const doAddSurface = (opts, atomsel) => {
      const ST = mol3dRef.current?.SurfaceType?.VDW ?? 1
      viewer.addSurface(ST, opts, atomsel || { model: protein }, undefined, undefined, () => {
        viewer.render()
      })
    }
    labelsRef.current = []

    const colorSpec = getColorSpec(colorMode, customColor)

    // Per-residue coloring modes (hydrophobicity, charge)
    if (colorMode === 'hydro' || colorMode === 'charge') {
      const baseColor = '#e2e8f0'
      const baseSpec = { color: baseColor }
      viewer.setStyle({ model: protein }, buildStyleSpec(style, baseSpec))

      if (colorMode === 'hydro') {
        Object.entries(HYDRO_SCALE).forEach(([resName]) => {
          const c = hydroColor(resName)
          viewer.addStyle({ resn: resName, model: protein }, buildStyleSpec(style, { color: c }))
        })
      } else {
        Object.entries(CHARGE_MAP).forEach(([resName, c]) => {
          viewer.addStyle({ resn: resName, model: protein }, buildStyleSpec(style, { color: c }))
        })
      }

      if (style === 'surface') {
        viewer.setStyle({ model: protein }, { cartoon: { color: '#888', opacity: 0.15 } })
        doAddSurface({
          opacity: 0.9,
          colorscheme: colorMode === 'hydro' ? 'amino' : 'charge',
        })
      }
    } else {
      // Standard coloring modes
      const cSpec = colorMode === 'default' ? { color: '#4a9eff' } : colorSpec

      if (style === 'surface') {
        viewer.setStyle({ model: protein }, { cartoon: { ...(cSpec || {}), opacity: 0.15 } })
        let surfOpts = { opacity: 0.9 }
        if (colorMode === 'bfactor') {
          surfOpts.colorscheme = { prop: 'b', gradient: 'rwb', min: 0, max: 100 }
        } else if (colorMode === 'default') {
          surfOpts.color = '#4a9eff'
        } else if (colorMode === 'custom') {
          surfOpts.color = customColor || '#4a9eff'
        } else if (colorMode === 'spectrum') {
          surfOpts.color = 'spectrum'
        } else if (cSpec?.colorscheme) {
          surfOpts.colorscheme = cSpec.colorscheme
        } else if (cSpec?.color) {
          surfOpts.color = cSpec.color
        }
        doAddSurface(surfOpts)
      } else {
        viewer.setStyle({ model: protein }, buildStyleSpec(style, cSpec))
      }
    }

    // Ligand display — show co-crystallized small molecules
    // Step 1: Always clear ALL HETATM atoms first
    viewer.setStyle({ hetflag: true, model: protein }, {})

    // Count ligand atoms for diagnostic
    try {
      const ligAtoms = viewer.selectedAtoms({ hetflag: true, not: { resn: [...SOLVENT_RESNAMES] }, model: protein })
      diag.ligandAtomsFound = ligAtoms.length
      if (ligAtoms.length === 0) diag.warnings.push('No ligand (HETATM) atoms found in this structure.')
    } catch (_) {}

    // Step 2: Re-apply ligand style only for non-solvent HETATM (using setStyle, not addStyle)
    if (ligandStyle !== 'off') {
      const ligSel = { hetflag: true, not: { resn: [...SOLVENT_RESNAMES] }, model: protein }
      const ligSpec = (() => {
        switch (ligandStyle) {
          case 'stick':   return { stick: { colorscheme: 'default', radius: 0.15, opacity: 1.0 } }
          case 'sphere':  return { sphere: { colorscheme: 'default', scale: 1.0, opacity: 0.9 } }
          case 'line':    return { line: { colorscheme: 'default', linewidth: 2 } }
          default:        return { stick: { colorscheme: 'default', radius: 0.2, opacity: 1.0 }, sphere: { colorscheme: 'default', radius: 0.4, opacity: 0.8 } }
        }
      })()
      viewer.setStyle(ligSel, ligSpec)
    }

    // Pocket overlay
    if (selectedPocket) {
      const residues = selectedPocket.residues || []
      diag.pocketResiduesRaw = residues.length
      diag.pocketRenderMode = pocketStyle
      diag.pocketMethod = selectedPocket.method || 'unknown'
      diag.pocketRawSamples = residues.slice(0, 5)
      diag.pocketHasCenter = !!(selectedPocket.center && selectedPocket.center.length === 3)
      diag.pocketCenter = selectedPocket.center || null

      // Parse residue labels — supports multiple formats
      const parsedResidues = residues.map(r => {
        const parts = r.split('_')
        if (parts.length >= 3) {
          return { resi: parseInt(parts[parts.length - 1]), chain: parts[parts.length - 2], raw: r, format: '3-part' }
        } else if (parts.length === 2) {
          const num = parseInt(parts[1])
          if (!isNaN(num)) return { resi: num, chain: null, raw: r, format: '2-part-num' }
          const match = parts[1].match(/(\d+)$/)
          if (match) return { resi: parseInt(match[1]), chain: parts[0].length === 1 ? parts[0] : null, raw: r, format: '2-part-legacy' }
        }
        const fallback = r.match(/(\d+)$/)
        if (fallback) return { resi: parseInt(fallback[1]), chain: null, raw: r, format: 'fallback' }
        return null
      }).filter(p => p != null && !isNaN(p.resi))

      diag.pocketResiduesParsed = parsedResidues.length
      diag.pocketParseFormat = parsedResidues[0]?.format || 'none'
      diag.pocketParsedSamples = parsedResidues.slice(0, 5).map(p => `resi=${p.resi} chain=${p.chain || '?'} (${p.format})`)
      if (parsedResidues.length === 0 && residues.length > 0) {
        diag.warnings.push(`All ${residues.length} pocket residues failed to parse. First raw value: "${residues[0]}". Pocket cannot render in Surface or Stick mode.`)
      }

      // Build per-chain selectors
      const chainGroups = {}
      parsedResidues.forEach(({ resi, chain }) => {
        const key = chain || '_all'
        if (!chainGroups[key]) chainGroups[key] = []
        chainGroups[key].push(resi)
      })
      diag.pocketChains = Object.keys(chainGroups)

      // Inspect the model: what chains and residue range does it have?
      try {
        const allAtoms = viewer.selectedAtoms({ model: protein })
        const modelChains = new Set()
        let minResi = Infinity, maxResi = -Infinity
        allAtoms.forEach(a => {
          if (a.chain) modelChains.add(a.chain)
          if (a.resi != null) {
            if (a.resi < minResi) minResi = a.resi
            if (a.resi > maxResi) maxResi = a.resi
          }
        })
        diag.modelChains = [...modelChains].sort()
        diag.modelResidueRange = minResi <= maxResi ? [minResi, maxResi] : null
      } catch (_) {
        diag.modelChains = []
        diag.modelResidueRange = null
      }

      // Check how many atoms actually match in the model — per chain group
      let totalMatched = 0
      diag.pocketChainMatches = {}
      Object.entries(chainGroups).forEach(([chain, resis]) => {
        const sel = chain === '_all'
          ? { resi: resis, model: protein }
          : { resi: resis, chain, model: protein }
        try {
          const matched = viewer.selectedAtoms(sel).length
          diag.pocketChainMatches[chain] = { residues: resis.length, atoms: matched }
          totalMatched += matched
        } catch (err) {
          diag.pocketChainMatches[chain] = { residues: resis.length, atoms: 0, error: err.message }
        }
      })
      diag.pocketAtomsMatched = totalMatched

      // Detailed warnings based on real data
      if (totalMatched === 0 && parsedResidues.length > 0) {
        const pocketChainKeys = Object.keys(chainGroups).filter(c => c !== '_all')
        const chainMismatch = pocketChainKeys.length > 0 && diag.modelChains.length > 0
          && !pocketChainKeys.some(c => diag.modelChains.includes(c))
        if (chainMismatch) {
          diag.warnings.push(
            `Chain mismatch: pocket references chain(s) [${pocketChainKeys.join(', ')}] but model only has chain(s) [${diag.modelChains.join(', ')}].`
          )
        }
        const pocketResis = parsedResidues.map(p => p.resi)
        const pocketMin = Math.min(...pocketResis), pocketMax = Math.max(...pocketResis)
        if (diag.modelResidueRange) {
          const [mMin, mMax] = diag.modelResidueRange
          const overlap = pocketMin <= mMax && pocketMax >= mMin
          if (!overlap) {
            diag.warnings.push(
              `Residue range mismatch: pocket residues span ${pocketMin}–${pocketMax} but model residues span ${mMin}–${mMax}. No overlap.`
            )
          }
        }
        if (diag.warnings.length === 0 || (diag.warnings.length > 0 && !chainMismatch)) {
          diag.warnings.push(
            `0 atoms matched for ${parsedResidues.length} parsed residues. Surface and Stick modes cannot render. Sphere mode works if center coordinates are available.`
          )
        }
      }

      // Surface-specific warnings
      if (pocketStyle === 'surface' && totalMatched > 0 && totalMatched < 20) {
        diag.warnings.push(`SES surface: only ${totalMatched} atoms selected. The surface may appear very small or incomplete.`)
      }

      const applyPocketStyle = (sel) => {
        if (pocketStyle === 'sphere+stick') {
          const center = selectedPocket.center
          if (center && center.length === 3) {
            viewer.addSphere({
              center: { x: center[0], y: center[1], z: center[2] },
              radius: pocketRadius,
              color: pocketColor,
              opacity: pocketOpacity,
              wireframe: false,
            })
            diag.pocketSphereRendered = true
          } else {
            diag.warnings.push('Sphere mode: no center coordinates available. Sphere cannot render.')
            diag.pocketSphereRendered = false
          }
        } else if (pocketStyle === 'surface') {
          diag.pocketSurfaceAttempted = true
          diag.pocketSurfaceType = 'VDW'
          // Sticks as anchor + VDW surface on pocket residues
          viewer.addStyle(sel, { stick: { color: pocketColor, radius: 0.1, opacity: 0.4 } })
          doAddSurface({ opacity: pocketOpacity, color: pocketColor }, sel)
        } else {
          viewer.addStyle(sel, {
            stick: { color: pocketColor, radius: 0.18 },
            sphere: { color: pocketColor, radius: 0.35 },
          })
          diag.pocketStickRendered = true
        }
      }

      // Apply per chain — always attempt rendering
      Object.entries(chainGroups).forEach(([chain, resis]) => {
        const sel = chain === '_all'
          ? { resi: resis, model: protein }
          : { resi: resis, chain, model: protein }
        applyPocketStyle(sel)
      })
    }

    // Selected residues highlight (multi-select)
    if (selectedResidues.length > 0) {
      const resis = selectedResidues.map(r => r.resi)
      viewer.addStyle(
        { resi: resis, model: protein },
        { stick: { color: '#06b6d4', radius: 0.16, opacity: 0.9 } }
      )
      // Add labels for each selected residue
      selectedResidues.forEach(r => {
        viewer.addLabel(`${r.resn} ${r.resi}`, {
          position: undefined,
          sel: { resi: r.resi, atom: 'CA', model: protein },
          fontSize: 11,
          fontColor: '#06b6d4',
          backgroundColor: 'rgba(0,0,0,0.7)',
          backgroundOpacity: 0.7,
          borderRadius: 3,
        })
      })
    }

    // Helper: build annotation style matching the current protein style with a custom color
    const buildAnnotationStyle = (color) => {
      switch (style) {
        case 'cartoon':    return { cartoon: { color } }
        case 'tube':       return { cartoon: { color, tubes: true, thickness: 0.4 } }
        case 'ball+stick': return { stick: { color, radius: 0.18 }, sphere: { color, radius: 0.4 } }
        case 'stick':      return { stick: { color, radius: 0.18 } }
        case 'line':       return { line: { color, linewidth: 3 } }
        case 'sphere':     return { sphere: { color } }
        case 'surface':    return { stick: { color, radius: 0.18 } }
        default:           return { cartoon: { color } }
      }
    }

    // UniProt annotations — follow the protein style with annotation color
    if (uniprotFeatures) {
      if (showActiveSites && uniprotFeatures.activeSites) {
        const resis = uniprotFeatures.activeSites.map(f => parseInt(f.position)).filter(n => !isNaN(n))
        if (resis.length > 0) {
          try {
            const matched = viewer.selectedAtoms({ resi: resis, model: protein }).length
            diag.activeSitesResolved = matched
            if (matched === 0) diag.warnings.push(`Active sites: ${resis.length} residues annotated but 0 atoms found in the model. Residue numbering may differ between UniProt and this structure.`)
          } catch (_) {}
          viewer.addStyle({ resi: resis, model: protein }, buildAnnotationStyle('#ef4444'))
          if (style === 'surface') {
            doAddSurface({ color: '#ef4444', opacity: 0.7 }, { resi: resis, model: protein })
          }
        }
      }
      if (showBindingSites && uniprotFeatures.bindingSites) {
        const resis = uniprotFeatures.bindingSites.flatMap(f => {
          if (f.start && f.end) {
            const arr = []
            for (let i = parseInt(f.start); i <= parseInt(f.end); i++) arr.push(i)
            return arr
          }
          return [parseInt(f.position)].filter(n => !isNaN(n))
        })
        if (resis.length > 0) {
          try {
            const matched = viewer.selectedAtoms({ resi: resis, model: protein }).length
            diag.bindingSitesResolved = matched
            if (matched === 0) diag.warnings.push(`Binding sites: ${resis.length} residues annotated but 0 atoms found in the model.`)
          } catch (_) {}
          viewer.addStyle({ resi: resis, model: protein }, buildAnnotationStyle('#f97316'))
          if (style === 'surface') {
            doAddSurface({ color: '#f97316', opacity: 0.7 }, { resi: resis, model: protein })
          }
        }
      }
      if (showDomains && uniprotFeatures.domains) {
        const domainColors = ['#6366f1', '#8b5cf6', '#ec4899', '#14b8a6', '#f59e0b', '#06b6d4']
        let domainsMatched = 0
        uniprotFeatures.domains.forEach((d, i) => {
          const start = parseInt(d.start), end = parseInt(d.end)
          if (!isNaN(start) && !isNaN(end)) {
            const resis = []
            for (let j = start; j <= end; j++) resis.push(j)
            try { domainsMatched += viewer.selectedAtoms({ resi: resis, model: protein }).length } catch (_) {}
            const dColor = domainColors[i % domainColors.length]
            viewer.addStyle({ resi: resis, model: protein }, buildAnnotationStyle(dColor))
            if (style === 'surface') {
              doAddSurface({ color: dColor, opacity: 0.7 }, { resi: resis })
            }
          }
        })
        diag.domainsResolved = domainsMatched
      }
    }

    // Custom zones — addStyle to overlay on top of base protein style
    customZones.forEach(zone => {
      if (!zone.visible || zone.residues.length === 0) return
      const zoneSel = { resi: zone.residues, model: protein }
      viewer.addStyle(zoneSel, buildAnnotationStyle(zone.color))
      if (style === 'surface') {
        doAddSurface({ color: zone.color, opacity: 0.7 }, { resi: zone.residues, model: protein })
      }
    })

    viewer.render()

    // Report diagnostic to parent
    if (onDiagnosticRef.current) onDiagnosticRef.current(diag)

  }, [ready, style, colorMode, customColor, pocketColor, pocketRadius, pocketOpacity, pocketStyle, selectedPocket, selectedResidues, uniprotFeatures, showDomains, showActiveSites, showBindingSites, ligandStyle, customZones])

  // ---------------------------------------------------------------------------
  // Actions
  // ---------------------------------------------------------------------------
  const handleReset = () => {
    const viewer = viewerRef.current
    if (!viewer) return
    viewer.zoomTo()
    viewer.zoom(0.85)
    viewer.render()
  }

  // Screenshot / Export PNG
  const handleScreenshot = useCallback(() => {
    const viewer = viewerRef.current
    if (!viewer) return
    try {
      const canvas = containerRef.current?.querySelector('canvas')
      if (!canvas) return
      // Force a render to ensure the canvas is up to date
      viewer.render()
      const dataUrl = canvas.toDataURL('image/png')
      const link = document.createElement('a')
      link.download = 'protein-structure.png'
      link.href = dataUrl
      link.click()
    } catch (err) {
      console.warn('Screenshot failed:', err)
    }
  }, [])

  // Zoom to pocket
  const handleZoomToPocket = useCallback(() => {
    const viewer = viewerRef.current
    if (!viewer || !selectedPocket) return
    const residues = selectedPocket.residues || []
    const resSels = residues.map(r => {
      const parts = r.split('_')
      if (parts.length >= 3) return parseInt(parts[parts.length - 1])
      if (parts.length === 2) {
        const n = parseInt(parts[1])
        if (!isNaN(n)) return n
        const m = parts[1].match(/(\d+)$/)
        return m ? parseInt(m[1]) : null
      }
      const fb = r.match(/(\d+)$/)
      return fb ? parseInt(fb[1]) : null
    }).filter(n => n != null && !isNaN(n))
    if (resSels.length > 0) {
      viewer.zoomTo({ resi: resSels })
      viewer.zoom(0.9)
    } else if (selectedPocket.center) {
      const c = selectedPocket.center
      viewer.zoomTo({ x: c[0], y: c[1], z: c[2] }, 12)
    }
    viewer.render()
  }, [selectedPocket])

  // Clear all labels
  const clearLabels = useCallback(() => {
    const viewer = viewerRef.current
    if (!viewer) return
    try { viewer.removeAllLabels() } catch (_) {}
    labelsRef.current = []
    viewer.render()
  }, [])

  // Clear distance measurements
  const clearDistances = useCallback(() => {
    const viewer = viewerRef.current
    if (!viewer) return
    distanceLinesRef.current = []
    setDistanceAtom1(null)
    // Re-render will clear shapes
    try { viewer.removeAllShapes() } catch (_) {}
    viewer.render()
  }, [])

  // Clear selected residues
  const clearSelectedResidues = useCallback(() => {
    setSelectedResidues([])
  }, [])

  // Fullscreen toggle
  const toggleFullscreen = useCallback(() => {
    const el = wrapperRef.current
    if (!el) return
    if (!document.fullscreenElement) {
      el.requestFullscreen().catch(() => {})
    } else {
      document.exitFullscreen().catch(() => {})
    }
  }, [])

  useEffect(() => {
    const handler = () => {
      const isFull = !!document.fullscreenElement
      setIsFullscreen(isFull)
      setTimeout(() => {
        const viewer = viewerRef.current
        if (viewer) {
          viewer.resize()
          viewer.render()
        }
      }, 100)
    }
    document.addEventListener('fullscreenchange', handler)
    return () => document.removeEventListener('fullscreenchange', handler)
  }, [])

  // Disable other modes when one is activated
  const toggleLabelMode = useCallback(() => {
    setLabelMode(v => !v)
    setDistanceMode(false)
    setDistanceAtom1(null)
  }, [])

  const toggleDistanceMode = useCallback(() => {
    setDistanceMode(v => !v)
    setLabelMode(false)
    setDistanceAtom1(null)
  }, [])

  // --- Custom zones helpers ---
  const addZone = useCallback(() => {
    const id = nextZoneId.current++
    const color = ZONE_COLORS[(id - 1) % ZONE_COLORS.length]
    setCustomZones(prev => [...prev, { id, name: `Zone ${id}`, color, residues: [], visible: true }])
    return id
  }, [])

  const updateZone = useCallback((zoneId, updates) => {
    setCustomZones(prev => prev.map(z => z.id === zoneId ? { ...z, ...updates } : z))
  }, [])

  const deleteZone = useCallback((zoneId) => {
    setCustomZones(prev => prev.filter(z => z.id !== zoneId))
    if (zoneEditMode === zoneId) setZoneEditMode(null)
  }, [zoneEditMode])

  const toggleZoneVisibility = useCallback((zoneId) => {
    setCustomZones(prev => prev.map(z => z.id === zoneId ? { ...z, visible: !z.visible } : z))
  }, [])

  const toggleZoneEditMode = useCallback((zoneId) => {
    setZoneEditMode(prev => prev === zoneId ? null : zoneId)
    // Disable other tools when entering zone edit mode
    if (zoneEditMode !== zoneId) {
      setLabelMode(false)
      setDistanceMode(false)
      setDistanceAtom1(null)
    }
  }, [zoneEditMode])

  const isDarkBg = bgColor === '#0f1923' || bgColor === '#000000' || bgColor === '#1a1a2e'

  // Active tool indicator
  const activeToolName = zoneEditMode
    ? `Zone edit: ${customZones.find(z => z.id === zoneEditMode)?.name || 'Zone'} (click to add/remove residues)`
    : labelMode ? 'Label' : distanceMode ? (distanceAtom1 ? 'Distance (click 2nd atom)' : 'Distance (click 1st atom)') : null

  return (
    <div ref={wrapperRef} className={`rounded-xl border border-gray-200 overflow-hidden ${isFullscreen ? 'flex flex-col' : ''}`} style={{ backgroundColor: bgColor }}>
      {/* Toolbar — Row 1: Style pills */}
      <div className={`flex items-center gap-2 px-3 py-1.5 border-b ${isDarkBg ? 'bg-[#1a2d42] border-gray-700/50' : 'bg-gray-100 border-gray-200'}`}>
        <div className="flex items-center gap-0.5 bg-black/20 rounded-lg p-0.5">
          {STYLES.map(s => (
            <button
              key={s}
              onClick={() => setStyle(s)}
              className={`px-2 py-1 text-[11px] rounded-md transition-all font-medium whitespace-nowrap ${
                style === s
                  ? 'bg-blue-600 text-white shadow-sm'
                  : isDarkBg
                    ? 'text-gray-400 hover:text-white hover:bg-white/10'
                    : 'text-gray-500 hover:text-gray-800 hover:bg-white/60'
              }`}
            >
              {STYLE_LABELS[s]}
            </button>
          ))}
        </div>

        <div className="w-px h-5 bg-gray-600/40" />

        <select
          value={colorMode}
          onChange={e => setColorMode(e.target.value)}
          className={`text-[11px] border rounded-md px-1.5 py-1 max-w-[130px] ${
            isDarkBg
              ? 'bg-[#0f1923] text-gray-300 border-gray-600'
              : 'bg-white text-gray-700 border-gray-300'
          }`}
        >
          {COLOR_GROUPS.map(group => (
            <optgroup key={group.label} label={group.label}>
              {group.options.map(opt => (
                <option key={opt.value} value={opt.value}>{opt.label}</option>
              ))}
            </optgroup>
          ))}
        </select>

        {colorMode === 'custom' && (
          <input
            type="color"
            value={customColor}
            onChange={e => setCustomColor(e.target.value)}
            className="w-5 h-5 rounded border-0 cursor-pointer"
            title="Pick protein color"
          />
        )}
      </div>

      {/* Toolbar — Row 2: Tool icons */}
      <div className={`flex items-center gap-0.5 px-3 py-1 border-b ${isDarkBg ? 'bg-[#162233] border-gray-700/50' : 'bg-gray-50 border-gray-200'}`}>
          <button onClick={handleScreenshot} title="Screenshot"
            className={`p-1.5 rounded-md transition-colors ${isDarkBg ? 'text-gray-400 hover:text-white hover:bg-white/10' : 'text-gray-500 hover:text-gray-800 hover:bg-gray-200'}`}>
            <CameraIcon />
          </button>

          {selectedPocket && (
            <button onClick={handleZoomToPocket} title="Zoom to pocket"
              className={`p-1.5 rounded-md transition-colors ${isDarkBg ? 'text-gray-400 hover:text-white hover:bg-white/10' : 'text-gray-500 hover:text-gray-800 hover:bg-gray-200'}`}>
              <CrosshairIcon />
            </button>
          )}

          <div className="w-px h-4 bg-gray-600/40" />

          <button onClick={toggleLabelMode} title={labelMode ? 'Disable label mode' : 'Label mode'}
            className={`p-1.5 rounded-md transition-colors ${labelMode ? 'bg-cyan-600 text-white' : isDarkBg ? 'text-gray-400 hover:text-white hover:bg-white/10' : 'text-gray-500 hover:text-gray-800 hover:bg-gray-200'}`}>
            <TagIcon />
          </button>
          <button onClick={toggleDistanceMode} title={distanceMode ? 'Disable distance mode' : 'Distance measurement'}
            className={`p-1.5 rounded-md transition-colors ${distanceMode ? 'bg-green-600 text-white' : isDarkBg ? 'text-gray-400 hover:text-white hover:bg-white/10' : 'text-gray-500 hover:text-gray-800 hover:bg-gray-200'}`}>
            <RulerIcon />
          </button>
          <button onClick={() => setShowSlab(v => !v)} title={showSlab ? 'Disable clipping' : 'Slab / Clip'}
            className={`p-1.5 rounded-md transition-colors ${showSlab ? 'bg-purple-600 text-white' : isDarkBg ? 'text-gray-400 hover:text-white hover:bg-white/10' : 'text-gray-500 hover:text-gray-800 hover:bg-gray-200'}`}>
            <ScissorsIcon />
          </button>

          {/* Custom Zones */}
          <div className="relative">
            <button ref={zonesBtnRef} onClick={() => setShowZonesPanel(v => !v)} title="Custom Zones"
              className={`p-1.5 rounded-md transition-colors ${showZonesPanel || customZones.length > 0 ? 'bg-emerald-600 text-white' : isDarkBg ? 'text-gray-400 hover:text-white hover:bg-white/10' : 'text-gray-500 hover:text-gray-800 hover:bg-gray-200'}`}>
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M17.657 18.657A8 8 0 016.343 7.343S7 9 9 10c0-2 .5-5 2.986-7C14 5 16.09 5.777 17.656 7.343A7.975 7.975 0 0120 13a7.975 7.975 0 01-2.343 5.657z" />
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9.879 16.121A3 3 0 1012.015 11L11 14H9c0 .768.293 1.536.879 2.121z" />
              </svg>
            </button>
            {showZonesPanel && (
              <div ref={zonesPanelRef} className={`absolute left-0 top-full mt-1 z-30 rounded-lg shadow-lg border p-3 min-w-[300px] max-w-[380px] space-y-2 ${
                isDarkBg ? 'bg-[#1a2d42] border-gray-600' : 'bg-white border-gray-200'
              }`}>
                <div className="flex items-center justify-between mb-1">
                  <p className={`text-xs font-semibold ${isDarkBg ? 'text-gray-300' : 'text-gray-600'}`}>Custom Zones</p>
                  <button onClick={addZone} className="text-xs px-2 py-0.5 rounded bg-emerald-600 text-white hover:bg-emerald-700 transition-colors">+ New Zone</button>
                </div>
                {customZones.length === 0 && (
                  <p className={`text-xs py-2 ${isDarkBg ? 'text-gray-500' : 'text-gray-400'}`}>
                    No zones defined. Create one to tag and color residue groups.
                  </p>
                )}
                {customZones.map(zone => (
                  <ZoneCard key={zone.id} zone={zone} isDarkBg={isDarkBg} isEditing={zoneEditMode === zone.id}
                    onUpdate={updateZone} onToggleVisibility={toggleZoneVisibility} onToggleEdit={toggleZoneEditMode} onDelete={deleteZone} />
                ))}
              </div>
            )}
          </div>

          <div className="w-px h-4 bg-gray-600/40" />

          <button onClick={() => setSpinning(v => !v)} title={spinning ? 'Stop spin' : 'Auto-rotate'}
            className={`p-1.5 rounded-md transition-colors ${spinning ? 'bg-blue-600 text-white' : isDarkBg ? 'text-gray-400 hover:text-white hover:bg-white/10' : 'text-gray-500 hover:text-gray-800 hover:bg-gray-200'}`}>
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
            </svg>
          </button>

          {/* Settings gear */}
          <div className="relative">
            <button ref={settingsBtnRef} onClick={() => setShowSettings(v => !v)} title="Viewer settings"
              className={`p-1.5 rounded-md transition-colors ${showSettings ? 'bg-blue-600 text-white' : isDarkBg ? 'text-gray-400 hover:text-white hover:bg-white/10' : 'text-gray-500 hover:text-gray-800 hover:bg-gray-200'}`}>
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10.325 4.317c.426-1.756 2.924-1.756 3.35 0a1.724 1.724 0 002.573 1.066c1.543-.94 3.31.826 2.37 2.37a1.724 1.724 0 001.066 2.573c1.756.426 1.756 2.924 0 3.35a1.724 1.724 0 00-1.066 2.573c.94 1.543-.826 3.31-2.37 2.37a1.724 1.724 0 00-2.573 1.066c-.426 1.756-2.924 1.756-3.35 0a1.724 1.724 0 00-2.573-1.066c-1.543.94-3.31-.826-2.37-2.37a1.724 1.724 0 00-1.066-2.573c-1.756-.426-1.756-2.924 0-3.35a1.724 1.724 0 001.066-2.573c-.94-1.543.826-3.31 2.37-2.37.996.608 2.296.07 2.572-1.065z" />
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 12a3 3 0 11-6 0 3 3 0 016 0z" />
              </svg>
            </button>
            {showSettings && (
              <div ref={settingsRef} className={`absolute right-0 top-full mt-1 z-30 rounded-lg shadow-lg border p-3 min-w-[240px] space-y-3 ${
                isDarkBg ? 'bg-[#1a2d42] border-gray-600' : 'bg-white border-gray-200'
              }`}>
                {/* Background */}
                <div>
                  <p className={`text-xs font-semibold mb-1.5 ${isDarkBg ? 'text-gray-300' : 'text-gray-600'}`}>Background</p>
                  <div className="flex items-center gap-1.5 flex-wrap">
                    {BG_PRESETS.map(bg => (
                      <button key={bg.value} onClick={() => setBgColor(bg.value)} title={bg.label}
                        className={`w-6 h-6 rounded-full border-2 transition-all ${bgColor === bg.value ? 'border-blue-500 scale-110' : 'border-gray-400 hover:border-gray-300'}`}
                        style={{ backgroundColor: bg.value }} />
                    ))}
                  </div>
                </div>
                {/* Pocket options */}
                {selectedPocket && (
                  <div className={`pt-2 border-t ${isDarkBg ? 'border-gray-600' : 'border-gray-200'}`}>
                    <p className={`text-xs font-semibold mb-1.5 ${isDarkBg ? 'text-gray-300' : 'text-gray-600'}`}>Pocket Display</p>
                    <div className="flex items-center gap-1.5 mb-2">
                      {POCKET_COLORS.map(pc => (
                        <button key={pc.value} onClick={() => setPocketColor(pc.value)} title={pc.label}
                          className={`w-5 h-5 rounded-full border-2 transition-all ${pocketColor === pc.value ? 'border-blue-500 scale-110' : 'border-gray-500/50 hover:border-gray-400'}`}
                          style={{ backgroundColor: pc.value }} />
                      ))}
                    </div>
                    <div className="flex items-center gap-1 mb-2">
                      {[
                        { v: 'sphere+stick', l: 'Sphere' },
                        { v: 'surface', l: 'Surface' },
                        { v: 'stick', l: 'Ball+Stick' },
                      ].map(ps => (
                        <button key={ps.v} onClick={() => setPocketStyle(ps.v)}
                          className={`px-2 py-0.5 text-xs rounded transition-colors ${
                            pocketStyle === ps.v ? 'bg-amber-500 text-white'
                              : isDarkBg ? 'text-gray-400 hover:text-white hover:bg-white/10' : 'text-gray-500 hover:text-gray-700 hover:bg-gray-100'
                          }`}>
                          {ps.l}
                        </button>
                      ))}
                    </div>
                    {pocketStyle === 'sphere+stick' && (
                      <div className="flex items-center gap-2 mb-1.5">
                        <span className={`text-xs w-12 ${isDarkBg ? 'text-gray-400' : 'text-gray-500'}`}>Radius</span>
                        <input type="range" min="2" max="12" step="0.5" value={pocketRadius}
                          onChange={e => setPocketRadius(parseFloat(e.target.value))} className="flex-1 h-1 accent-amber-500" />
                        <span className={`text-xs w-6 text-right ${isDarkBg ? 'text-gray-400' : 'text-gray-500'}`}>{pocketRadius}</span>
                      </div>
                    )}
                    <div className="flex items-center gap-2">
                      <span className={`text-xs w-12 ${isDarkBg ? 'text-gray-400' : 'text-gray-500'}`}>Opacity</span>
                      <input type="range" min="0.05" max="0.8" step="0.05" value={pocketOpacity}
                        onChange={e => setPocketOpacity(parseFloat(e.target.value))} className="flex-1 h-1 accent-amber-500" />
                      <span className={`text-xs w-6 text-right ${isDarkBg ? 'text-gray-400' : 'text-gray-500'}`}>{Math.round(pocketOpacity * 100)}%</span>
                    </div>
                  </div>
                )}
                {/* Clear buttons */}
                <div className={`pt-2 border-t flex flex-wrap gap-1.5 ${isDarkBg ? 'border-gray-600' : 'border-gray-200'}`}>
                  <button onClick={clearLabels}
                    className={`px-2 py-0.5 text-xs rounded border transition-colors ${isDarkBg ? 'border-gray-600 text-gray-400 hover:text-white hover:bg-white/10' : 'border-gray-300 text-gray-500 hover:bg-gray-100'}`}>
                    Clear labels
                  </button>
                  <button onClick={clearDistances}
                    className={`px-2 py-0.5 text-xs rounded border transition-colors ${isDarkBg ? 'border-gray-600 text-gray-400 hover:text-white hover:bg-white/10' : 'border-gray-300 text-gray-500 hover:bg-gray-100'}`}>
                    Clear distances
                  </button>
                  {selectedResidues.length > 0 && (
                    <button onClick={clearSelectedResidues}
                      className={`px-2 py-0.5 text-xs rounded border transition-colors ${isDarkBg ? 'border-gray-600 text-gray-400 hover:text-white hover:bg-white/10' : 'border-gray-300 text-gray-500 hover:bg-gray-100'}`}>
                      Clear selection ({selectedResidues.length})
                    </button>
                  )}
                </div>
              </div>
            )}
          </div>

          <button onClick={handleReset} title="Reset view"
            className={`p-1.5 rounded-md transition-colors ${isDarkBg ? 'text-gray-400 hover:text-white hover:bg-white/10' : 'text-gray-500 hover:text-gray-800 hover:bg-gray-200'}`}>
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0zM10 7v3m0 0v3m0-3h3m-3 0H7" />
            </svg>
          </button>

          <button onClick={toggleFullscreen} title={isFullscreen ? 'Exit fullscreen' : 'Fullscreen'}
            className={`p-1.5 rounded-md transition-colors ${isDarkBg ? 'text-gray-400 hover:text-white hover:bg-white/10' : 'text-gray-500 hover:text-gray-800 hover:bg-gray-200'}`}>
            {isFullscreen ? (
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
              </svg>
            ) : (
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 8V4m0 0h4M4 4l5 5m11-1V4m0 0h-4m4 0l-5 5M4 16v4m0 0h4m-4 0l5-5m11 5l-5-5m5 5v-4m0 4h-4" />
              </svg>
            )}
          </button>
      </div>

      {/* Active tool indicator */}
      {activeToolName && (
        <div className={`flex items-center gap-2 px-3 py-1 text-xs border-b ${
          isDarkBg ? 'bg-[#1a2d42]/80 border-gray-700' : 'bg-blue-50 border-gray-200'
        }`}>
          <span className="w-2 h-2 rounded-full bg-green-400 animate-pulse" />
          <span className={isDarkBg ? 'text-gray-300' : 'text-gray-600'}>
            Tool: <strong>{activeToolName}</strong>
          </span>
          <button
            onClick={() => { setLabelMode(false); setDistanceMode(false); setDistanceAtom1(null); setZoneEditMode(null) }}
            className={`ml-auto text-xs px-1.5 py-0.5 rounded ${isDarkBg ? 'text-gray-400 hover:text-white hover:bg-white/10' : 'text-gray-500 hover:bg-gray-200'}`}
          >
            Cancel
          </button>
        </div>
      )}

      {/* Slab / Clip inline controls */}
      {showSlab && (
        <div className={`flex items-center gap-4 px-3 py-1.5 border-b text-xs ${
          isDarkBg ? 'bg-purple-900/20 border-gray-700' : 'bg-purple-50 border-gray-200'
        }`}>
          <span className={`font-semibold ${isDarkBg ? 'text-purple-300' : 'text-purple-600'}`}>Clip</span>
          <div className="flex items-center gap-2 flex-1 max-w-[200px]">
            <span className={isDarkBg ? 'text-gray-400' : 'text-gray-500'}>Near</span>
            <input
              type="range" min="-100" max="0" step="1"
              value={slabNear}
              onChange={e => setSlabNear(parseInt(e.target.value))}
              className="flex-1 h-1 accent-purple-500"
            />
          </div>
          <div className="flex items-center gap-2 flex-1 max-w-[200px]">
            <span className={isDarkBg ? 'text-gray-400' : 'text-gray-500'}>Far</span>
            <input
              type="range" min="0" max="100" step="1"
              value={slabFar}
              onChange={e => setSlabFar(parseInt(e.target.value))}
              className="flex-1 h-1 accent-purple-500"
            />
          </div>
          <button
            onClick={() => { setSlabNear(-50); setSlabFar(50) }}
            className={`px-1.5 py-0.5 rounded ${isDarkBg ? 'text-gray-400 hover:text-white hover:bg-white/10' : 'text-gray-500 hover:bg-gray-200'}`}
          >
            Reset
          </button>
          <button
            onClick={() => setShowSlab(false)}
            className={`px-1.5 py-0.5 rounded ${isDarkBg ? 'text-gray-400 hover:text-white hover:bg-white/10' : 'text-gray-500 hover:bg-gray-200'}`}
          >
            Close
          </button>
        </div>
      )}

      {/* Viewer */}
      <div className={`relative ${isFullscreen ? 'flex-1' : ''}`} style={isFullscreen ? {} : { height }}>
        {loading && (
          <div className="absolute inset-0 flex flex-col items-center justify-center z-10">
            <BindXLogo variant="loading" size={56} label="Loading 3D structure..." />
          </div>
        )}
        {error && (
          <div className="absolute inset-0 flex flex-col items-center justify-center z-10">
            <BindXLogo variant="error" size={48} />
            <p className="text-red-400 text-sm mt-2 mb-1">Failed to load structure</p>
            <p className="text-gray-500 text-xs">{error}</p>
            <button onClick={initViewer} className="mt-3 px-3 py-1 text-xs bg-blue-600 text-white rounded hover:bg-blue-700">Retry</button>
          </div>
        )}
        <div
          ref={containerRef}
          className="w-full h-full"
          style={{
            visibility: loading || error ? 'hidden' : 'visible',
            cursor: labelMode ? 'crosshair' : distanceMode ? 'crosshair' : 'default',
          }}
        />

        {/* Color legend overlay */}
        {ready && colorMode === 'bfactor' && (
          <div className={`absolute bottom-2 left-2 rounded px-2 py-1 text-xs flex items-center gap-2 pointer-events-none ${isDarkBg ? 'bg-black/60 text-gray-300' : 'bg-white/80 text-gray-600'}`}>
            <span>pLDDT / B-factor:</span>
            <span className="flex items-center gap-0.5">
              <span className="w-3 h-3 rounded" style={{ background: '#2563eb' }} /> Low
              <span className="w-3 h-3 rounded mx-0.5" style={{ background: '#ffffff' }} />
              <span className="w-3 h-3 rounded" style={{ background: '#dc2626' }} /> High
            </span>
          </div>
        )}
        {ready && colorMode === 'hydro' && (
          <div className={`absolute bottom-2 left-2 rounded px-2 py-1 text-xs flex items-center gap-2 pointer-events-none ${isDarkBg ? 'bg-black/60 text-gray-300' : 'bg-white/80 text-gray-600'}`}>
            <span>Hydrophobicity:</span>
            <span className="flex items-center gap-0.5">
              <span className="w-3 h-3 rounded" style={{ background: '#3b82f6' }} /> Hydrophilic
              <span className="w-3 h-3 rounded mx-0.5" style={{ background: '#ffffff' }} />
              <span className="w-3 h-3 rounded" style={{ background: '#ef4444' }} /> Hydrophobic
            </span>
          </div>
        )}
        {ready && colorMode === 'charge' && (
          <div className={`absolute bottom-2 left-2 rounded px-2 py-1 text-xs flex items-center gap-2 pointer-events-none ${isDarkBg ? 'bg-black/60 text-gray-300' : 'bg-white/80 text-gray-600'}`}>
            <span className="flex items-center gap-1">
              <span className="w-3 h-3 rounded" style={{ background: '#3b82f6' }} /> + (Arg, Lys)
              <span className="w-3 h-3 rounded" style={{ background: '#ef4444' }} /> - (Asp, Glu)
              <span className="w-3 h-3 rounded" style={{ background: '#e2e8f0' }} /> Neutral
            </span>
          </div>
        )}

        {/* Selected residues chip (bottom-right) */}
        {ready && selectedResidues.length > 0 && (
          <div className={`absolute bottom-2 right-2 rounded-lg px-2.5 py-1.5 text-xs max-w-[200px] ${isDarkBg ? 'bg-black/70 text-gray-300' : 'bg-white/90 text-gray-600 border border-gray-200'}`}>
            <div className="flex items-center justify-between mb-1">
              <span className="font-semibold text-cyan-400">Selected ({selectedResidues.length})</span>
              <button onClick={clearSelectedResidues} className="text-gray-500 hover:text-white ml-2">x</button>
            </div>
            <div className="flex flex-wrap gap-1">
              {selectedResidues.slice(0, 8).map(r => (
                <span key={r.key} className="bg-cyan-900/40 px-1.5 py-0.5 rounded text-cyan-300">{r.resn}{r.resi}</span>
              ))}
              {selectedResidues.length > 8 && <span className="text-gray-500">+{selectedResidues.length - 8}</span>}
            </div>
          </div>
        )}
      </div>

      {/* Legend bar */}
      <div className={`flex items-center gap-3 px-3 py-2 border-t text-xs flex-wrap ${
        isDarkBg ? 'bg-[#1a2d42] border-gray-700 text-gray-400' : 'bg-gray-100 border-gray-200 text-gray-500'
      }`}>
        <span className="flex items-center gap-1.5">
          <span className="w-2.5 h-2.5 rounded-full bg-blue-400" /> Protein
        </span>
        {selectedPocket && (
          <span className="flex items-center gap-1.5">
            <span className="w-2.5 h-2.5 rounded-full" style={{ backgroundColor: pocketColor }} /> Pocket
          </span>
        )}
        <span className="flex items-center gap-1.5">
          <span className="w-2.5 h-2.5 rounded-full bg-emerald-400" style={{ opacity: ligandStyle !== 'off' ? 1 : 0.3 }} />
          <select
            value={ligandStyle}
            onChange={(e) => setLigandStyle(e.target.value)}
            className={`bg-transparent border rounded px-1 py-0 text-xs cursor-pointer outline-none ${
              ligandStyle !== 'off'
                ? 'border-emerald-500/50 text-emerald-300'
                : `border-gray-600 ${isDarkBg ? 'text-gray-500' : 'text-gray-400'}`
            }`}
          >
            <option value="off" className="bg-gray-800 text-gray-300">Ligand: Off</option>
            <option value="ball+stick" className="bg-gray-800 text-gray-300">Ligand: Ball+Stick</option>
            <option value="stick" className="bg-gray-800 text-gray-300">Ligand: Stick</option>
            <option value="sphere" className="bg-gray-800 text-gray-300">Ligand: Sphere</option>
            <option value="line" className="bg-gray-800 text-gray-300">Ligand: Line</option>
          </select>
        </span>
        {selectedResidues.length > 0 && (
          <span className="flex items-center gap-1.5">
            <span className="w-2.5 h-2.5 rounded-full bg-cyan-400" /> Selected
          </span>
        )}
        {uniprotFeatures?.activeSites?.length > 0 && (
          <button
            onClick={() => setShowActiveSites(v => !v)}
            className={`flex items-center gap-1.5 px-2 py-0.5 rounded border transition-colors ${
              showActiveSites ? 'bg-red-900/30 border-red-500/50 text-red-300' : `border-gray-600 ${isDarkBg ? 'text-gray-500 hover:text-gray-300' : 'text-gray-400 hover:text-gray-600'}`
            }`}
          >
            <span className="w-2.5 h-2.5 rounded-full bg-red-500" style={{ opacity: showActiveSites ? 1 : 0.3 }} />
            Active Sites
          </button>
        )}
        {uniprotFeatures?.bindingSites?.length > 0 && (
          <button
            onClick={() => setShowBindingSites(v => !v)}
            className={`flex items-center gap-1.5 px-2 py-0.5 rounded border transition-colors ${
              showBindingSites ? 'bg-orange-900/30 border-orange-500/50 text-orange-300' : `border-gray-600 ${isDarkBg ? 'text-gray-500 hover:text-gray-300' : 'text-gray-400 hover:text-gray-600'}`
            }`}
          >
            <span className="w-2.5 h-2.5 rounded-full bg-orange-500" style={{ opacity: showBindingSites ? 1 : 0.3 }} />
            Binding Sites
          </button>
        )}
        {uniprotFeatures?.domains?.length > 0 && (
          <button
            onClick={() => setShowDomains(v => !v)}
            className={`flex items-center gap-1.5 px-2 py-0.5 rounded border transition-colors ${
              showDomains ? 'bg-indigo-900/30 border-indigo-500/50 text-indigo-300' : `border-gray-600 ${isDarkBg ? 'text-gray-500 hover:text-gray-300' : 'text-gray-400 hover:text-gray-600'}`
            }`}
          >
            <span className="w-2.5 h-2.5 rounded-full bg-indigo-500" style={{ opacity: showDomains ? 1 : 0.3 }} />
            Domains
          </button>
        )}
        {customZones.filter(z => z.visible && z.residues.length > 0).map(zone => (
          <span key={zone.id} className="flex items-center gap-1.5">
            <span className="w-2.5 h-2.5 rounded-full" style={{ backgroundColor: zone.color }} />
            <span style={{ color: zone.color }}>{zone.name}</span>
          </span>
        ))}
        <span className={`ml-auto hidden sm:inline ${isDarkBg ? 'text-gray-500' : 'text-gray-400'}`}>
          {zoneEditMode
            ? 'Click residues to add/remove from zone'
            : labelMode || distanceMode
              ? 'Click atoms to interact'
              : 'Drag: rotate | Scroll: zoom | Right-click: pan | Click: select residue'}
        </span>
      </div>

    </div>
  )
}

