import React, { useEffect, useRef, useState, useCallback } from 'react'

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
        sphere: { ...(colorSpec || {}), radius: 0.3, opacity: 0.8 },
      }
    case 'stick':
      return { stick: { ...(colorSpec || {}), radius: 0.15 } }
    case 'line':
      return { line: { ...(colorSpec || {}) } }
    case 'sphere':
      return { sphere: { ...(colorSpec || {}), radius: 1.5 } }
    case 'tube':
      return { cartoon: { ...(colorSpec || {}), opacity: 0.85, tubes: true, thickness: 0.4 } }
    default: // cartoon
      return { cartoon: { ...(colorSpec || {}), opacity: 0.85 } }
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

// Props: pdbUrl, pdbData, selectedPocket, uniprotFeatures, height
export default function ProteinViewer({ pdbUrl, pdbData, selectedPocket, uniprotFeatures, height = 480 }) {
  const wrapperRef = useRef(null)
  const containerRef = useRef(null)
  const viewerRef = useRef(null)
  const settingsRef = useRef(null)
  const settingsBtnRef = useRef(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState(null)
  const [ready, setReady] = useState(false)
  const [style, setStyle] = useState('cartoon')
  const [colorMode, setColorMode] = useState('default')
  const [customColor, setCustomColor] = useState('#4a9eff')
  const [bgColor, setBgColor] = useState('#0f1923')
  const [pocketColor, setPocketColor] = useState('#f59e0b')
  const [pocketRadius, setPocketRadius] = useState(5)
  const [pocketOpacity, setPocketOpacity] = useState(0.25)
  const [pocketStyle, setPocketStyle] = useState('sphere+stick')
  const [spinning, setSpinning] = useState(false)
  const [showDomains, setShowDomains] = useState(false)
  const [showActiveSites, setShowActiveSites] = useState(false)
  const [showBindingSites, setShowBindingSites] = useState(false)
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
      viewer.setStyle({}, { cartoon: { color: '#4a9eff', opacity: 0.85 } })
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

    try { viewer.removeAllSurfaces() } catch (_) {}
    try { viewer.removeAllShapes() } catch (_) {}
    try { viewer.removeAllLabels() } catch (_) {}
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
        // Surface mode: thin cartoon inside + opaque surface on top
        viewer.setStyle({ model: protein }, { cartoon: { color: '#888', opacity: 0.08 } })
        const surfaceType = window.$3Dmol?.SurfaceType?.VDW ?? 1
        viewer.addSurface(surfaceType, {
          opacity: 0.82,
          colorscheme: colorMode === 'hydro' ? 'amino' : 'charge',
        }, { model: protein })
      }
    } else {
      // Standard coloring modes
      const cSpec = colorMode === 'default' ? { color: '#4a9eff' } : colorSpec

      if (style === 'surface') {
        // Surface mode: thin cartoon inside + opaque colored surface
        viewer.setStyle({ model: protein }, { cartoon: { ...(cSpec || {}), opacity: 0.08 } })
        const surfaceType = window.$3Dmol?.SurfaceType?.VDW ?? 1
        let surfOpts = { opacity: 0.82 }
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
        viewer.addSurface(surfaceType, surfOpts, { model: protein })
      } else {
        viewer.setStyle({ model: protein }, buildStyleSpec(style, cSpec))
      }
    }

    // Pocket overlay — harmonized sphere+stick vs surface
    if (selectedPocket) {
      const residues = selectedPocket.residues || []
      const resSels = residues.map(r => {
        const parts = r.split('_')
        return parts.length >= 2 ? parseInt(parts[parts.length - 1]) : null
      }).filter(n => n != null)

      if (resSels.length > 0) {
        if (pocketStyle === 'sphere+stick') {
          // Show residue sticks AND a translucent sphere at center
          viewer.addStyle(
            { resi: resSels, model: protein },
            { stick: { color: pocketColor, radius: 0.14, opacity: 0.85 } }
          )
          const center = selectedPocket.center
          if (center && center.length === 3) {
            viewer.addSphere({
              center: { x: center[0], y: center[1], z: center[2] },
              radius: pocketRadius,
              color: pocketColor,
              opacity: pocketOpacity,
              wireframe: false,
            })
          }
        } else if (pocketStyle === 'surface') {
          // Surface covering the pocket residues
          viewer.addStyle(
            { resi: resSels, model: protein },
            { stick: { color: pocketColor, radius: 0.1, opacity: 0.5 } }
          )
          const surfType = window.$3Dmol?.SurfaceType?.VDW ?? 1
          viewer.addSurface(surfType, {
            opacity: pocketOpacity + 0.15,
            color: pocketColor,
          }, { resi: resSels, model: protein })
        } else {
          // stick (ball+stick)
          viewer.addStyle(
            { resi: resSels, model: protein },
            {
              stick: { color: pocketColor, radius: 0.18, opacity: 0.9 },
              sphere: { color: pocketColor, radius: 0.35, opacity: 0.6 },
            }
          )
        }
      }
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

    // UniProt annotations
    if (uniprotFeatures) {
      if (showActiveSites && uniprotFeatures.activeSites) {
        const resis = uniprotFeatures.activeSites.map(f => parseInt(f.position)).filter(n => !isNaN(n))
        if (resis.length > 0) {
          viewer.addStyle({ resi: resis, model: protein }, {
            stick: { color: '#ef4444', radius: 0.18, opacity: 0.95 },
            sphere: { color: '#ef4444', radius: 0.4, opacity: 0.7 },
          })
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
          viewer.addStyle({ resi: resis, model: protein }, {
            stick: { color: '#f97316', radius: 0.18, opacity: 0.95 },
          })
        }
      }
      if (showDomains && uniprotFeatures.domains) {
        const domainColors = ['#6366f1', '#8b5cf6', '#ec4899', '#14b8a6', '#f59e0b', '#06b6d4']
        uniprotFeatures.domains.forEach((d, i) => {
          const start = parseInt(d.start), end = parseInt(d.end)
          if (!isNaN(start) && !isNaN(end)) {
            const resis = []
            for (let j = start; j <= end; j++) resis.push(j)
            viewer.addStyle({ resi: resis, model: protein }, {
              cartoon: { color: domainColors[i % domainColors.length], opacity: 0.95 },
            })
          }
        })
      }
    }

    viewer.render()
  }, [ready, style, colorMode, customColor, pocketColor, pocketRadius, pocketOpacity, pocketStyle, selectedPocket, selectedResidues, uniprotFeatures, showDomains, showActiveSites, showBindingSites])

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
      return parts.length >= 2 ? parseInt(parts[parts.length - 1]) : null
    }).filter(n => n != null)
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

  const isDarkBg = bgColor === '#0f1923' || bgColor === '#000000' || bgColor === '#1a1a2e'

  // Active tool indicator
  const activeToolName = labelMode ? 'Label' : distanceMode ? (distanceAtom1 ? 'Distance (click 2nd atom)' : 'Distance (click 1st atom)') : null

  return (
    <div ref={wrapperRef} className={`rounded-xl border border-gray-200 overflow-hidden ${isFullscreen ? 'flex flex-col' : ''}`} style={{ backgroundColor: bgColor }}>
      {/* Toolbar */}
      <div className={`flex items-center justify-between px-3 py-2 border-b ${isDarkBg ? 'bg-[#1a2d42] border-gray-700' : 'bg-gray-100 border-gray-200'}`}>
        {/* Style buttons */}
        <div className="flex items-center gap-0.5 flex-wrap">
          {STYLES.map(s => (
            <button
              key={s}
              onClick={() => setStyle(s)}
              className={`px-2 py-1 text-xs rounded transition-colors ${
                style === s
                  ? 'bg-blue-600 text-white'
                  : isDarkBg
                    ? 'text-gray-400 hover:text-white hover:bg-white/10'
                    : 'text-gray-500 hover:text-gray-800 hover:bg-gray-200'
              }`}
            >
              {STYLE_LABELS[s]}
            </button>
          ))}
        </div>

        {/* Right controls */}
        <div className="flex items-center gap-1">
          {/* Color mode dropdown */}
          <select
            value={colorMode}
            onChange={e => setColorMode(e.target.value)}
            className={`text-xs border rounded px-2 py-1 max-w-[140px] ${
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
              className="w-6 h-6 rounded border-0 cursor-pointer"
              title="Pick protein color"
            />
          )}

          <div className="w-px h-5 bg-gray-600 mx-0.5" />

          {/* Screenshot */}
          <button
            onClick={handleScreenshot}
            title="Save current view as PNG image"
            className={`p-1 rounded transition-colors ${isDarkBg ? 'text-gray-400 hover:text-white hover:bg-white/10' : 'text-gray-500 hover:text-gray-800 hover:bg-gray-200'}`}
          >
            <CameraIcon />
          </button>

          {/* Zoom to pocket */}
          {selectedPocket && (
            <button
              onClick={handleZoomToPocket}
              title="Zoom to selected pocket — centers the view on the binding site"
              className={`p-1 rounded transition-colors ${isDarkBg ? 'text-gray-400 hover:text-white hover:bg-white/10' : 'text-gray-500 hover:text-gray-800 hover:bg-gray-200'}`}
            >
              <CrosshairIcon />
            </button>
          )}

          {/* Label toggle */}
          <button
            onClick={toggleLabelMode}
            title={labelMode ? 'Disable label mode' : 'Label mode — click any atom to display its residue name (e.g. ALA 42:A)'}
            className={`p-1 rounded transition-colors ${
              labelMode
                ? 'bg-cyan-600 text-white'
                : isDarkBg
                  ? 'text-gray-400 hover:text-white hover:bg-white/10'
                  : 'text-gray-500 hover:text-gray-800 hover:bg-gray-200'
            }`}
          >
            <TagIcon />
          </button>

          {/* Distance measurement toggle */}
          <button
            onClick={toggleDistanceMode}
            title={distanceMode ? 'Disable distance mode' : 'Distance measurement — click two atoms to measure the distance in Angstroms'}
            className={`p-1 rounded transition-colors ${
              distanceMode
                ? 'bg-green-600 text-white'
                : isDarkBg
                  ? 'text-gray-400 hover:text-white hover:bg-white/10'
                  : 'text-gray-500 hover:text-gray-800 hover:bg-gray-200'
            }`}
          >
            <RulerIcon />
          </button>

          {/* Slab / Clip toggle */}
          <button
            onClick={() => setShowSlab(v => !v)}
            title={showSlab ? 'Disable clipping' : 'Slab / Clip — slice the structure with near/far planes to see inside'}
            className={`p-1 rounded transition-colors ${
              showSlab
                ? 'bg-purple-600 text-white'
                : isDarkBg
                  ? 'text-gray-400 hover:text-white hover:bg-white/10'
                  : 'text-gray-500 hover:text-gray-800 hover:bg-gray-200'
            }`}
          >
            <ScissorsIcon />
          </button>

          <div className="w-px h-5 bg-gray-600 mx-0.5" />

          {/* Spin */}
          <button
            onClick={() => setSpinning(v => !v)}
            title={spinning ? 'Stop spin' : 'Auto-rotate'}
            className={`p-1 rounded transition-colors ${
              spinning
                ? 'bg-blue-600 text-white'
                : isDarkBg
                  ? 'text-gray-400 hover:text-white hover:bg-white/10'
                  : 'text-gray-500 hover:text-gray-800 hover:bg-gray-200'
            }`}
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
            </svg>
          </button>

          {/* Settings gear */}
          <div className="relative">
            <button
              ref={settingsBtnRef}
              onClick={() => setShowSettings(v => !v)}
              title="Viewer settings"
              className={`p-1 rounded transition-colors ${
                showSettings
                  ? 'bg-blue-600 text-white'
                  : isDarkBg
                    ? 'text-gray-400 hover:text-white hover:bg-white/10'
                    : 'text-gray-500 hover:text-gray-800 hover:bg-gray-200'
              }`}
            >
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
                      <button
                        key={bg.value}
                        onClick={() => setBgColor(bg.value)}
                        title={bg.label}
                        className={`w-6 h-6 rounded-full border-2 transition-all ${
                          bgColor === bg.value ? 'border-blue-500 scale-110' : 'border-gray-400 hover:border-gray-300'
                        }`}
                        style={{ backgroundColor: bg.value }}
                      />
                    ))}
                  </div>
                </div>

                {/* Pocket options */}
                {selectedPocket && (
                  <div className={`pt-2 border-t ${isDarkBg ? 'border-gray-600' : 'border-gray-200'}`}>
                    <p className={`text-xs font-semibold mb-1.5 ${isDarkBg ? 'text-gray-300' : 'text-gray-600'}`}>Pocket Display</p>

                    {/* Color */}
                    <div className="flex items-center gap-1.5 mb-2">
                      {POCKET_COLORS.map(pc => (
                        <button
                          key={pc.value}
                          onClick={() => setPocketColor(pc.value)}
                          title={pc.label}
                          className={`w-5 h-5 rounded-full border-2 transition-all ${
                            pocketColor === pc.value ? 'border-blue-500 scale-110' : 'border-gray-500/50 hover:border-gray-400'
                          }`}
                          style={{ backgroundColor: pc.value }}
                        />
                      ))}
                    </div>

                    {/* Style */}
                    <div className="flex items-center gap-1 mb-2">
                      {[
                        { v: 'sphere+stick', l: 'Sphere' },
                        { v: 'surface', l: 'Surface' },
                        { v: 'stick', l: 'Ball+Stick' },
                      ].map(ps => (
                        <button
                          key={ps.v}
                          onClick={() => setPocketStyle(ps.v)}
                          className={`px-2 py-0.5 text-xs rounded transition-colors ${
                            pocketStyle === ps.v
                              ? 'bg-amber-500 text-white'
                              : isDarkBg
                                ? 'text-gray-400 hover:text-white hover:bg-white/10'
                                : 'text-gray-500 hover:text-gray-700 hover:bg-gray-100'
                          }`}
                        >
                          {ps.l}
                        </button>
                      ))}
                    </div>

                    {/* Radius (only for sphere mode) */}
                    {pocketStyle === 'sphere+stick' && (
                      <div className="flex items-center gap-2 mb-1.5">
                        <span className={`text-xs w-12 ${isDarkBg ? 'text-gray-400' : 'text-gray-500'}`}>Radius</span>
                        <input
                          type="range" min="2" max="12" step="0.5"
                          value={pocketRadius}
                          onChange={e => setPocketRadius(parseFloat(e.target.value))}
                          className="flex-1 h-1 accent-amber-500"
                        />
                        <span className={`text-xs w-6 text-right ${isDarkBg ? 'text-gray-400' : 'text-gray-500'}`}>{pocketRadius}</span>
                      </div>
                    )}

                    {/* Opacity */}
                    <div className="flex items-center gap-2">
                      <span className={`text-xs w-12 ${isDarkBg ? 'text-gray-400' : 'text-gray-500'}`}>Opacity</span>
                      <input
                        type="range" min="0.05" max="0.8" step="0.05"
                        value={pocketOpacity}
                        onChange={e => setPocketOpacity(parseFloat(e.target.value))}
                        className="flex-1 h-1 accent-amber-500"
                      />
                      <span className={`text-xs w-6 text-right ${isDarkBg ? 'text-gray-400' : 'text-gray-500'}`}>{Math.round(pocketOpacity * 100)}%</span>
                    </div>
                  </div>
                )}

                {/* Clear buttons */}
                <div className={`pt-2 border-t flex flex-wrap gap-1.5 ${isDarkBg ? 'border-gray-600' : 'border-gray-200'}`}>
                  <button
                    onClick={clearLabels}
                    className={`px-2 py-0.5 text-xs rounded border transition-colors ${isDarkBg ? 'border-gray-600 text-gray-400 hover:text-white hover:bg-white/10' : 'border-gray-300 text-gray-500 hover:bg-gray-100'}`}
                  >
                    Clear labels
                  </button>
                  <button
                    onClick={clearDistances}
                    className={`px-2 py-0.5 text-xs rounded border transition-colors ${isDarkBg ? 'border-gray-600 text-gray-400 hover:text-white hover:bg-white/10' : 'border-gray-300 text-gray-500 hover:bg-gray-100'}`}
                  >
                    Clear distances
                  </button>
                  {selectedResidues.length > 0 && (
                    <button
                      onClick={clearSelectedResidues}
                      className={`px-2 py-0.5 text-xs rounded border transition-colors ${isDarkBg ? 'border-gray-600 text-gray-400 hover:text-white hover:bg-white/10' : 'border-gray-300 text-gray-500 hover:bg-gray-100'}`}
                    >
                      Clear selection ({selectedResidues.length})
                    </button>
                  )}
                </div>
              </div>
            )}
          </div>

          {/* Reset view */}
          <button
            onClick={handleReset}
            title="Reset view"
            className={`p-1 rounded transition-colors ${isDarkBg ? 'text-gray-400 hover:text-white hover:bg-white/10' : 'text-gray-500 hover:text-gray-800 hover:bg-gray-200'}`}
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0zM10 7v3m0 0v3m0-3h3m-3 0H7" />
            </svg>
          </button>

          {/* Fullscreen */}
          <button
            onClick={toggleFullscreen}
            title={isFullscreen ? 'Exit fullscreen' : 'Fullscreen'}
            className={`p-1 rounded transition-colors ${isDarkBg ? 'text-gray-400 hover:text-white hover:bg-white/10' : 'text-gray-500 hover:text-gray-800 hover:bg-gray-200'}`}
          >
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
            onClick={() => { setLabelMode(false); setDistanceMode(false); setDistanceAtom1(null) }}
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
            <div className="w-10 h-10 border-2 border-blue-500 border-t-transparent rounded-full animate-spin mb-3" />
            <p className={`text-sm ${isDarkBg ? 'text-gray-400' : 'text-gray-500'}`}>Loading 3D structure...</p>
          </div>
        )}
        {error && (
          <div className="absolute inset-0 flex flex-col items-center justify-center z-10">
            <p className="text-red-400 text-sm mb-2">Failed to load structure</p>
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
        <span className={`ml-auto hidden sm:inline ${isDarkBg ? 'text-gray-500' : 'text-gray-400'}`}>
          {labelMode || distanceMode
            ? 'Click atoms to interact'
            : 'Drag: rotate | Scroll: zoom | Right-click: pan | Click: select residue'}
        </span>
      </div>
    </div>
  )
}
