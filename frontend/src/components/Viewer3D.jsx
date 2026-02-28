import React, { useEffect, useRef, useState, useCallback } from 'react'
import { getProteinUrl, getPoseUrl, getReportUrl } from '../api.js'
import BindXLogo from './BindXLogo.jsx'

// --------------------------------------------------
// 3Dmol loader — CDN script tag approach
// --------------------------------------------------
function load3Dmol() {
  return new Promise((resolve, reject) => {
    if (window.$3Dmol) {
      resolve(window.$3Dmol)
      return
    }
    const existing = document.getElementById('3dmol-script')
    if (existing) {
      existing.addEventListener('load', () => resolve(window.$3Dmol))
      existing.addEventListener('error', reject)
      return
    }
    const script = document.createElement('script')
    script.id = '3dmol-script'
    script.src = 'https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js'
    script.onload = () => {
      if (window.$3Dmol) resolve(window.$3Dmol)
      else reject(new Error('3Dmol unavailable after loading'))
    }
    script.onerror = () => reject(new Error('Failed to load 3Dmol.js'))
    document.head.appendChild(script)
  })
}

// --------------------------------------------------
// Color mode helpers
// --------------------------------------------------
const HYDRO_SCALE = {
  ILE: 4.5, VAL: 4.2, LEU: 3.8, PHE: 2.8, CYS: 2.5, MET: 1.9, ALA: 1.8,
  GLY: -0.4, THR: -0.7, SER: -0.8, TRP: -0.9, TYR: -1.3, PRO: -1.6,
  HIS: -3.2, GLU: -3.5, GLN: -3.5, ASP: -3.5, ASN: -3.5, LYS: -3.9, ARG: -4.5,
}

function hydroColor(atom) {
  const val = HYDRO_SCALE[atom.resn] ?? 0
  const t = (val + 4.5) / 9.0
  if (t < 0.5) {
    const s = t * 2
    return `rgb(${Math.round(s * 255)},${Math.round(s * 255)},255)`
  }
  const s = (t - 0.5) * 2
  return `rgb(255,${Math.round((1 - s) * 255)},${Math.round((1 - s) * 255)})`
}

function getColorSpec(mode) {
  switch (mode) {
    case 'ss':      return { colorscheme: 'ssJmol' }
    case 'hydro':   return { colorfunc: hydroColor }
    case 'bfactor': return { colorscheme: { prop: 'b', gradient: 'rwb', min: 0, max: 100 } }
    case 'chain':   return { colorscheme: 'chain' }
    default:        return { color: '#4a9eff' }
  }
}

// Domain palette — cycles through 6 distinct colours
const DOMAIN_COLORS = ['#6366f1', '#8b5cf6', '#ec4899', '#14b8a6', '#f59e0b', '#06b6d4']

// Human-readable label for each color mode
const COLOR_MODE_LABELS = {
  default:  'Default (blue)',
  ss:       'Secondary Structure',
  hydro:    'Hydrophobicity',
  bfactor:  'B-factor / pLDDT',
  chain:    'Chain',
}

// --------------------------------------------------
// Simplified toolbar — used in V3 dashboard flow
// Props: onBack, onPrev, onNext, jobId, onReset
// --------------------------------------------------
function SimplifiedToolbar({ onBack, onPrev, onNext, jobId, onReset }) {
  return (
    <div className="flex items-center justify-between px-4 py-3 bg-bx-surface border-b border-bx-bg">
      {/* Left: back */}
      <button
        onClick={onBack}
        className="flex items-center gap-1.5 px-3 py-1.5 text-white/80 hover:text-white hover:bg-white/10 rounded-lg text-sm font-medium transition-colors"
      >
        <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10 19l-7-7m0 0l7-7m-7 7h18" />
        </svg>
        Back
      </button>

      {/* Center: title */}
      <h3 className="text-white font-semibold text-sm flex items-center gap-2">
        <svg className="w-4 h-4 text-bx-mint" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 12a3 3 0 11-6 0 3 3 0 016 0z" />
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M2.458 12C3.732 7.943 7.523 5 12 5c4.478 0 8.268 2.943 9.542 7-1.274 4.057-5.064 7-9.542 7-4.477 0-8.268-2.943-9.542-7z" />
        </svg>
        3D Viewer
      </h3>

      {/* Right: navigation + PDF */}
      <div className="flex items-center gap-1.5">
        <button
          onClick={onPrev}
          disabled={!onPrev}
          className="flex items-center gap-1 px-3 py-1.5 text-white/80 hover:text-white hover:bg-white/10 rounded-lg text-sm font-medium transition-colors disabled:opacity-30 disabled:cursor-not-allowed"
        >
          <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 19l-7-7 7-7" />
          </svg>
          Previous
        </button>
        <button
          onClick={onNext}
          disabled={!onNext}
          className="flex items-center gap-1 px-3 py-1.5 text-white/80 hover:text-white hover:bg-white/10 rounded-lg text-sm font-medium transition-colors disabled:opacity-30 disabled:cursor-not-allowed"
        >
          Next
          <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
          </svg>
        </button>
        {jobId && (
          <a
            href={getReportUrl(jobId)}
            target="_blank"
            rel="noopener noreferrer"
            download
            title="Download PDF report"
            className="flex items-center gap-1 px-3 py-1.5 bg-bx-mint hover:bg-bx-mint-dim text-white rounded-lg text-sm font-medium transition-colors"
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                d="M12 10v6m0 0l-3-3m3 3l3-3m2 8H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
            </svg>
            PDF
          </a>
        )}
        {/* Reset view */}
        <button
          onClick={onReset}
          title="Reset view"
          className="p-1.5 text-white/70 hover:text-white hover:bg-white/10 rounded transition-colors"
        >
          <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 8V4m0 0h4M4 4l5 5m11-1V4m0 0h-4m4 0l-5 5M4 16v4m0 0h4m-4 0l5-5m11 5l-5-5m5 5v-4m0 4h-4" />
          </svg>
        </button>
      </div>
    </div>
  )
}

// --------------------------------------------------
// Full toolbar (V2 style) — with color mode dropdown
// Props: topResults, selectedPoseIndex, onPoseSelect,
//        currentStyle, onStyleChange,
//        colorMode, onColorChange,
//        onReset
// --------------------------------------------------
function FullToolbar({
  topResults,
  selectedPoseIndex,
  onPoseSelect,
  currentStyle,
  onStyleChange,
  colorMode,
  onColorChange,
  onReset,
}) {
  return (
    <div className="flex items-center justify-between px-4 py-3 bg-bx-surface border-b border-bx-bg">
      <h3 className="text-white font-semibold text-sm flex items-center gap-2">
        <svg className="w-4 h-4 text-bx-mint" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 12a3 3 0 11-6 0 3 3 0 016 0z" />
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M2.458 12C3.732 7.943 7.523 5 12 5c4.478 0 8.268 2.943 9.542 7-1.274 4.057-5.064 7-9.542 7-4.477 0-8.268-2.943-9.542-7z" />
        </svg>
        3D Viewer
      </h3>

      <div className="flex items-center gap-2 flex-wrap justify-end">
        {/* Pose selector */}
        {topResults.length > 0 && (
          <select
            value={selectedPoseIndex ?? ''}
            onChange={(e) => onPoseSelect(e.target.value !== '' ? parseInt(e.target.value) : null)}
            className="text-sm bg-bx-elevated text-white border border-white/20 rounded px-2 py-1 focus:outline-none focus:ring-1 focus:ring-bx-mint"
          >
            <option value="">-- Select a ligand --</option>
            {topResults.map((r, i) => (
              <option key={i} value={i}>
                #{i + 1} {r.name || r.ligand_name || `Ligand ${i + 1}`}
              </option>
            ))}
          </select>
        )}

        {/* Style buttons */}
        <div className="flex rounded overflow-hidden border border-white/20">
          {['cartoon', 'surface', 'stick'].map((style) => (
            <button
              key={style}
              onClick={() => onStyleChange(style)}
              className={`px-2 py-1 text-sm font-medium transition-colors ${
                currentStyle === style
                  ? 'bg-bx-mint text-white'
                  : 'text-white/70 hover:text-white hover:bg-white/10'
              }`}
            >
              {style === 'cartoon' ? 'Ribbon' : style === 'surface' ? 'Surface' : 'Stick'}
            </button>
          ))}
        </div>

        {/* Color mode dropdown */}
        <select
          value={colorMode}
          onChange={(e) => onColorChange(e.target.value)}
          title="Color mode"
          className="text-sm bg-bx-elevated text-white border border-white/20 rounded px-2 py-1 focus:outline-none focus:ring-1 focus:ring-bx-mint"
        >
          {Object.entries(COLOR_MODE_LABELS).map(([value, label]) => (
            <option key={value} value={value}>{label}</option>
          ))}
        </select>

        {/* Reset view */}
        <button
          onClick={onReset}
          title="Reset view"
          className="p-1.5 text-white/70 hover:text-white hover:bg-white/10 rounded transition-colors"
        >
          <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 8V4m0 0h4M4 4l5 5m11-1V4m0 0h-4m4 0l-5 5M4 16v4m0 0h4m-4 0l5-5m11 5l-5-5m5 5v-4m0 4h-4" />
          </svg>
        </button>
      </div>
    </div>
  )
}

// --------------------------------------------------
// Viewer3D
// Props:
//   jobId, results, selectedPoseIndex, onPoseSelect — V2 full mode
//   simplified — boolean, enables simplified toolbar
//   onBack, onPrev, onNext — simplified mode callbacks
//   pocketCenter, pocketResidues — pocket overlay
//   uniprotFeatures — { activeSites, bindingSites, domains }
//     activeSites:  [{ position }]  or [{ residues:[...] }]
//     bindingSites: [{ position }]  or [{ residues:[...] }]
//     domains:      [{ start, end, name }]
// --------------------------------------------------
export default function Viewer3D({
  jobId,
  results,
  selectedPoseIndex,
  onPoseSelect,
  simplified = false,
  onBack,
  onPrev,
  onNext,
  pocketCenter = null,
  pocketResidues = null,
  uniprotFeatures = null,
}) {
  const containerRef = useRef(null)
  const viewerRef = useRef(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState(null)
  const [viewerReady, setViewerReady] = useState(false)
  const [currentStyle, setCurrentStyle] = useState('cartoon')
  const [colorMode, setColorMode] = useState('default')
  const [ligandLoaded, setLigandLoaded] = useState(false)
  const [noPose, setNoPose] = useState(false)
  const [showPocket, setShowPocket] = useState(true)
  const [annotations, setAnnotations] = useState({
    activeSites: false,
    bindingSites: false,
    domains: false,
  })

  const surfaceActiveRef = useRef(false)

  const topResults = results?.results?.slice(0, 10) || []

  // --------------------------------------------------
  // Helpers: extract residue numbers from UniProt feature arrays
  // --------------------------------------------------
  function extractResidues(featureArr) {
    if (!featureArr || !Array.isArray(featureArr)) return []
    const nums = []
    featureArr.forEach((f) => {
      // Flat position
      if (f.position != null) {
        const n = parseInt(f.position)
        if (!isNaN(n)) nums.push(n)
      }
      // Nested residues array
      if (Array.isArray(f.residues)) {
        f.residues.forEach((r) => {
          const n = parseInt(r?.position ?? r)
          if (!isNaN(n)) nums.push(n)
        })
      }
      // start/end range (domains)
      if (f.start != null && f.end != null) {
        const start = parseInt(f.start)
        const end = parseInt(f.end)
        if (!isNaN(start) && !isNaN(end)) {
          for (let i = start; i <= end; i++) nums.push(i)
        }
      }
    })
    return [...new Set(nums)]
  }

  // --------------------------------------------------
  // initViewer — load protein, set up viewer
  // --------------------------------------------------
  const initViewer = useCallback(async () => {
    if (!containerRef.current || !jobId) return

    setLoading(true)
    setError(null)
    setLigandLoaded(false)
    setNoPose(false)

    try {
      const $3Dmol = await load3Dmol()

      // Clear any existing viewer
      if (viewerRef.current) {
        try { viewerRef.current.clear() } catch (_) {}
      }
      containerRef.current.innerHTML = ''

      // Create viewer
      const viewer = $3Dmol.createViewer(containerRef.current, {
        backgroundColor: '#0f1923',
        id: 'dockit-viewer',
      })
      viewerRef.current = viewer

      // Load protein structure
      const proteinUrl = getProteinUrl(jobId)
      const response = await fetch(proteinUrl)
      if (!response.ok) throw new Error(`Error loading protein: ${response.status}`)
      const pdbData = await response.text()

      viewer.addModel(pdbData, 'pdb')
      viewer.setStyle({}, { cartoon: { color: '#4a9eff', opacity: 0.85 } })
      viewer.zoomTo()
      viewer.render()

      setViewerReady(true)
      setLoading(false)
    } catch (err) {
      console.error('[Viewer3D] Init error:', err)
      setError(err.message || 'Unable to load the 3D viewer.')
      setLoading(false)
    }
  }, [jobId])

  // --------------------------------------------------
  // loadPose — fetch and display a docking pose
  // --------------------------------------------------
  const loadPose = useCallback(async (poseIndex) => {
    const viewer = viewerRef.current
    if (!viewer || poseIndex === null || poseIndex === undefined) return

    try {
      // Remove previous ligand models (all non-protein models)
      const models = viewer.getModelList()
      if (models && models.length > 1) {
        for (let i = 1; i < models.length; i++) {
          viewer.removeModel(models[i])
        }
      }

      const poseUrl = getPoseUrl(jobId, poseIndex)
      let poseLoaded = false
      try {
        const response = await fetch(poseUrl)
        if (response.ok) {
          const poseData = await response.text()
          const isSdf = /V2000|V3000|\$\$\$\$/.test(poseData)
          const isPdbqt = /BRANCH|ENDROOT|TORSDOF/.test(poseData)
          const format = isSdf ? 'sdf' : isPdbqt ? 'pdbqt' : 'pdb'
          viewer.addModel(poseData, format)
          poseLoaded = true
        }
      } catch (fetchErr) {
        console.warn(`[Viewer3D] Pose fetch failed for index ${poseIndex}:`, fetchErr)
      }

      // No pose file available — show "not docked" banner, do NOT render fake ligand
      if (!poseLoaded) {
        console.info(`[Viewer3D] No docking pose for index ${poseIndex} — showing banner`)
        setNoPose(true)
        setLigandLoaded(false)
        viewer.render()
        return
      }

      setNoPose(false)

      // Last model is the ligand
      const models2 = viewer.getModelList()
      const ligandModel = models2[models2.length - 1]
      viewer.setStyle({ model: ligandModel }, {
        stick: { color: '#00e6a0', radius: 0.2 },
        sphere: { color: '#00e6a0', radius: 0.35, opacity: 0.8 },
      })

      // Semi-transparent surface on protein when ligand is loaded (simplified mode default)
      if (simplified) {
        const proteinModel = models2[0]
        if (proteinModel) {
          viewer.setStyle({ model: proteinModel }, {
            cartoon: { color: '#4a9eff', opacity: 0.6 },
          })
          viewer.addSurface(
            window.$3Dmol?.SurfaceType?.VDW || 1,
            { opacity: 0.25, color: '#4a9eff' },
            { model: proteinModel }
          )
        }
      }

      // Zoom to ligand region (not the whole protein) so the pose is visible
      try {
        viewer.zoomTo({ model: ligandModel })
        viewer.zoom(0.7) // zoom out slightly to show surrounding pocket
      } catch (_) {
        viewer.zoomTo()
      }
      viewer.render()
      setLigandLoaded(true)
    } catch (err) {
      console.error('[Viewer3D] Pose load error:', err)
    }
  }, [jobId, simplified])

  // --------------------------------------------------
  // applyAllStyles — unified callback
  // Sets base protein style (cartoon/surface/stick x colorMode),
  // then overlays pocket highlight, then UniProt annotations.
  // This replaces the old applyStyle + pocket useEffect pair so
  // changing protein style no longer resets pocket residue sticks.
  // --------------------------------------------------
  const applyAllStyles = useCallback(() => {
    const viewer = viewerRef.current
    if (!viewer || !viewerReady) return

    const models = viewer.getModelList()
    const proteinModel = models?.[0]
    if (!proteinModel) return

    // 1. Remove any existing surfaces and shapes before rebuilding
    if (surfaceActiveRef.current) {
      try { viewer.removeAllSurfaces() } catch (_) {}
      surfaceActiveRef.current = false
    }
    try { viewer.removeAllShapes() } catch (_) {}

    const colorSpec = getColorSpec(colorMode)

    // 2. Base protein style
    if (currentStyle === 'surface') {
      // Coloured cartoon underlay (opacity 0.3), then translucent white surface
      viewer.setStyle({ model: proteinModel }, {
        cartoon: { ...colorSpec, opacity: 0.3 },
      })
      const $3Dmol = window.$3Dmol
      const surfaceType = $3Dmol?.SurfaceType?.VDW ?? 1
      viewer.addSurface(
        surfaceType,
        { opacity: 0.35, color: '#e2e8f0' },
        { model: proteinModel }
      )
      surfaceActiveRef.current = true
    } else if (currentStyle === 'cartoon') {
      viewer.setStyle({ model: proteinModel }, {
        cartoon: { ...colorSpec, opacity: 0.85 },
      })
    } else {
      // stick
      const stickSpec = colorMode === 'default'
        ? { colorscheme: 'default', radius: 0.15 }
        : { ...colorSpec, radius: 0.15 }
      viewer.setStyle({ model: proteinModel }, { stick: stickSpec })
    }

    // 3. Pocket overlay
    if (showPocket && pocketCenter && Array.isArray(pocketCenter) && pocketCenter.length === 3) {
      const [x, y, z] = pocketCenter
      viewer.addSphere({
        center: { x, y, z },
        radius: 6,
        color: '#f59e0b',
        opacity: 0.25,
        wireframe: false,
      })

      if (pocketResidues && pocketResidues.length > 0) {
        const resSels = pocketResidues.map((r) => {
          const parts = r.split('_')
          return parts.length >= 2 ? parseInt(parts[parts.length - 1]) : null
        }).filter((n) => n != null)

        if (resSels.length > 0) {
          viewer.addStyle(
            { resi: resSels, model: proteinModel },
            { stick: { color: '#f59e0b', radius: 0.12, opacity: 0.7 } }
          )
        }
      }
    }

    // 4. UniProt annotations
    if (uniprotFeatures) {
      // Active Sites — red spheres + sticks
      if (annotations.activeSites) {
        const residues = extractResidues(uniprotFeatures.activeSites)
        if (residues.length > 0) {
          viewer.addStyle(
            { resi: residues, model: proteinModel },
            {
              stick: { color: '#ef4444', radius: 0.18, opacity: 0.95 },
              sphere: { color: '#ef4444', radius: 0.4, opacity: 0.7 },
            }
          )
        }
      }

      // Binding Sites — orange sticks
      if (annotations.bindingSites) {
        const residues = extractResidues(uniprotFeatures.bindingSites)
        if (residues.length > 0) {
          viewer.addStyle(
            { resi: residues, model: proteinModel },
            { stick: { color: '#f97316', radius: 0.18, opacity: 0.95 } }
          )
        }
      }

      // Domains — per-domain coloured cartoon regions
      if (annotations.domains && Array.isArray(uniprotFeatures.domains)) {
        uniprotFeatures.domains.forEach((domain, idx) => {
          const color = DOMAIN_COLORS[idx % DOMAIN_COLORS.length]
          const start = parseInt(domain.start)
          const end = parseInt(domain.end)
          if (!isNaN(start) && !isNaN(end)) {
            const resiRange = []
            for (let i = start; i <= end; i++) resiRange.push(i)
            viewer.addStyle(
              { resi: resiRange, model: proteinModel },
              { cartoon: { color, opacity: 0.95 } }
            )
          }
        })
      }
    }

    viewer.render()
  }, [
    viewerReady,
    currentStyle,
    colorMode,
    showPocket,
    pocketCenter,
    pocketResidues,
    uniprotFeatures,
    annotations,
  ])

  // --------------------------------------------------
  // Effects
  // --------------------------------------------------

  // Initialize viewer when jobId changes
  useEffect(() => {
    initViewer()
  }, [initViewer])

  // Re-apply all styles whenever any relevant state changes
  useEffect(() => {
    applyAllStyles()
  }, [applyAllStyles])

  // Load pose when selection changes
  useEffect(() => {
    if (viewerReady && selectedPoseIndex !== null && selectedPoseIndex !== undefined) {
      loadPose(selectedPoseIndex)
    }
  }, [selectedPoseIndex, viewerReady, loadPose])

  // --------------------------------------------------
  // Handlers
  // --------------------------------------------------
  const handleStyleChange = (style) => {
    setCurrentStyle(style)
    // applyAllStyles triggers automatically via the effect above
  }

  const handleColorChange = (mode) => {
    setColorMode(mode)
    // applyAllStyles triggers automatically via the effect above
  }

  const toggleAnnotation = (key) => {
    setAnnotations((prev) => ({ ...prev, [key]: !prev[key] }))
    // applyAllStyles triggers automatically via the effect above
  }

  const handleReset = () => {
    const viewer = viewerRef.current
    if (!viewer) return
    viewer.zoomTo()
    viewer.render()
  }

  // --------------------------------------------------
  // Annotation toggle button — reusable sub-component
  // --------------------------------------------------
  function AnnotationToggle({ label, colorDot, active, onToggle }) {
    return (
      <button
        onClick={onToggle}
        className={`flex items-center gap-1.5 px-2 py-0.5 rounded border transition-colors ${
          active
            ? 'border-opacity-60 text-gray-700'
            : 'bg-white border-gray-200 text-gray-400 hover:border-gray-300'
        }`}
        style={active ? { backgroundColor: colorDot + '18', borderColor: colorDot } : {}}
      >
        <span
          className="w-3 h-3 rounded-full inline-block flex-shrink-0"
          style={{ backgroundColor: colorDot, opacity: active ? 1 : 0.3 }}
        />
        {label} {active ? 'ON' : 'OFF'}
      </button>
    )
  }

  // Determine whether to show annotation toggles in the legend
  const hasActiveSites = uniprotFeatures?.activeSites?.length > 0
  const hasBindingSites = uniprotFeatures?.bindingSites?.length > 0
  const hasDomains = uniprotFeatures?.domains?.length > 0

  // --------------------------------------------------
  // Render
  // --------------------------------------------------
  return (
    <div className="card overflow-hidden">
      {/* Toolbar */}
      {simplified ? (
        <SimplifiedToolbar
          onBack={onBack}
          onPrev={onPrev}
          onNext={onNext}
          jobId={jobId}
          onReset={handleReset}
        />
      ) : (
        <FullToolbar
          topResults={topResults}
          selectedPoseIndex={selectedPoseIndex}
          onPoseSelect={onPoseSelect}
          currentStyle={currentStyle}
          onStyleChange={handleStyleChange}
          colorMode={colorMode}
          onColorChange={handleColorChange}
          onReset={handleReset}
        />
      )}

      {/* Viewer container */}
      <div className="relative bg-gray-900" style={{ height: '420px' }}>
        {loading && (
          <div className="absolute inset-0 flex flex-col items-center justify-center z-10 bg-gray-900">
            <BindXLogo variant="loading" size={48} label="Loading 3D structure..." />
          </div>
        )}

        {error && (
          <div className="absolute inset-0 flex flex-col items-center justify-center z-10 bg-gray-900">
            <BindXLogo variant="error" size={48} />
            <p className="text-red-400 text-sm font-medium mb-1 mt-2">Visualization unavailable</p>
            <p className="text-white/40 text-sm text-center px-8">{error}</p>
            <button
              onClick={initViewer}
              className="mt-4 px-4 py-2 text-sm bg-bx-surface text-white rounded-lg hover:bg-bx-elevated transition-colors"
            >
              Retry
            </button>
          </div>
        )}

        <div
          ref={containerRef}
          className="w-full h-full viewer-container"
          style={{ visibility: loading || error ? 'hidden' : 'visible' }}
        />

        {/* "Not docked" banner — shown when selected molecule has no pose */}
        {noPose && !loading && !error && (
          <div className="absolute inset-x-0 bottom-0 flex items-center justify-center z-10 pointer-events-none pb-6">
            <div className="bg-gray-800/90 backdrop-blur-sm border border-gray-600 rounded-xl px-6 py-3 flex items-center gap-3 shadow-lg pointer-events-auto">
              <svg className="w-5 h-5 text-amber-400 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                  d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-2.5L13.732 4.5c-.77-.833-2.694-.833-3.464 0L3.34 16.5c-.77.833.192 2.5 1.732 2.5z" />
              </svg>
              <div>
                <p className="text-white text-sm font-semibold">No docking pose available</p>
                <p className="text-gray-400 text-sm">This molecule was not docked — only real docked poses are shown</p>
              </div>
            </div>
          </div>
        )}
      </div>

      {/* Legend + annotation toggles */}
      <div className="flex items-center gap-4 px-4 py-2 bg-gray-50 border-t border-gray-100 text-sm text-gray-500 flex-wrap">
        {/* Static legend items */}
        <span className="flex items-center gap-1.5">
          <span className="w-3 h-3 rounded-full bg-blue-400 inline-block" />
          Protein
        </span>
        <span className="flex items-center gap-1.5">
          <span className="w-3 h-3 rounded-full bg-bx-mint inline-block" />
          Ligand
        </span>

        {/* Pocket toggle */}
        {pocketCenter && (
          <button
            onClick={() => setShowPocket((v) => !v)}
            className={`flex items-center gap-1.5 px-2 py-0.5 rounded border transition-colors ${
              showPocket
                ? 'bg-amber-50 border-amber-300 text-amber-700'
                : 'bg-white border-gray-200 text-gray-400 hover:border-gray-300'
            }`}
          >
            <span className="w-3 h-3 rounded-full bg-amber-400 inline-block" style={{ opacity: showPocket ? 1 : 0.3 }} />
            Pocket {showPocket ? 'ON' : 'OFF'}
          </button>
        )}

        {/* UniProt annotation toggles */}
        {hasActiveSites && (
          <AnnotationToggle
            label="Active Sites"
            colorDot="#ef4444"
            active={annotations.activeSites}
            onToggle={() => toggleAnnotation('activeSites')}
          />
        )}
        {hasBindingSites && (
          <AnnotationToggle
            label="Binding Sites"
            colorDot="#f97316"
            active={annotations.bindingSites}
            onToggle={() => toggleAnnotation('bindingSites')}
          />
        )}
        {hasDomains && (
          <AnnotationToggle
            label="Domains"
            colorDot="#6366f1"
            active={annotations.domains}
            onToggle={() => toggleAnnotation('domains')}
          />
        )}

        <span className="ml-auto text-gray-400 hidden sm:inline">
          Rotate | Scroll: zoom | Right: pan
        </span>
      </div>
    </div>
  )
}
