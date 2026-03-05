import React, { useEffect, useRef, useState, useCallback } from 'react'
import { v9SmilesToMolblock } from '../api'

// Module-level cache: SMILES -> molblock
const molblockCache = new Map()

// Shared 3Dmol CDN loader (same pattern as MiniProteinViewer)
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

const LIGAND_STYLES = [
  { value: 'ball+stick', label: 'Ball+Stick' },
  { value: 'stick', label: 'Stick' },
  { value: 'sphere', label: 'Sphere' },
  { value: 'surface', label: 'Surface' },
  { value: 'line', label: 'Line' },
]

function getStyleSpec(mode) {
  switch (mode) {
    case 'stick':      return { stick: { radius: 0.15, colorscheme: 'default' } }
    case 'sphere':     return { sphere: { scale: 0.6, colorscheme: 'default' } }
    case 'line':       return { line: { linewidth: 2, colorscheme: 'default' } }
    case 'surface':    return { stick: { radius: 0.1, colorscheme: 'default', opacity: 0.4 } }
    default:           return { stick: { radius: 0.15, colorscheme: 'default' }, sphere: { scale: 0.25, colorscheme: 'default' } }
  }
}

// ---------------------------------------------------------------------------
// LigandViewer3D — 3D molecule viewer from SMILES with display controls
// ---------------------------------------------------------------------------
export default function LigandViewer3D({ smiles, height = 280 }) {
  const containerRef = useRef(null)
  const viewerRef = useRef(null)
  const animRef = useRef(null)
  const spinRef = useRef(true)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState(null)
  const [displayStyle, setDisplayStyle] = useState('ball+stick')
  const [spinning, setSpinning] = useState(true)
  const [ready, setReady] = useState(false)

  useEffect(() => {
    if (!smiles || !containerRef.current) return
    let cancelled = false

    async function init() {
      setLoading(true)
      setError(null)
      setReady(false)

      try {
        // Fetch molblock (with cache)
        let molblock = molblockCache.get(smiles)
        if (!molblock) {
          const data = await v9SmilesToMolblock(smiles)
          molblock = data.molblock
          molblockCache.set(smiles, molblock)
        }
        if (cancelled) return

        // Load 3Dmol
        const $3Dmol = await load3Dmol()
        if (cancelled) return

        // Clear container
        containerRef.current.innerHTML = ''

        const viewer = $3Dmol.createViewer(containerRef.current, {
          backgroundColor: '#0e1628',
          antialias: true,
        })
        viewerRef.current = viewer

        viewer.addModel(molblock, 'sdf')
        viewer.setStyle({}, getStyleSpec('ball+stick'))
        viewer.zoomTo()
        viewer.render()
        setLoading(false)
        setReady(true)

        // Auto-rotation loop
        function spin() {
          if (cancelled) return
          if (spinRef.current) {
            viewer.rotate(0.3, 'y')
            viewer.render()
          }
          animRef.current = requestAnimationFrame(spin)
        }
        animRef.current = requestAnimationFrame(spin)
      } catch (err) {
        if (!cancelled) {
          console.warn('[LigandViewer3D] Load failed:', err.message)
          setError(err.userMessage || err.message || '3D generation failed')
          setLoading(false)
        }
      }
    }

    init()

    return () => {
      cancelled = true
      if (animRef.current) cancelAnimationFrame(animRef.current)
      if (viewerRef.current) {
        try { viewerRef.current.clear() } catch {}
        viewerRef.current = null
      }
    }
  }, [smiles])

  // Update style when changed
  useEffect(() => {
    if (!ready || !viewerRef.current) return
    const viewer = viewerRef.current
    viewer.setStyle({}, getStyleSpec(displayStyle))
    // Add surface if needed
    viewer.removeAllSurfaces()
    if (displayStyle === 'surface') {
      try { viewer.addSurface(window.$3Dmol.SurfaceType.VDW, { opacity: 0.65, colorscheme: 'default' }) } catch {}
    }
    viewer.render()
  }, [displayStyle, ready])

  // Sync spin ref
  useEffect(() => {
    spinRef.current = spinning
  }, [spinning])

  const handleScreenshot = useCallback(() => {
    if (!viewerRef.current) return
    const dataUrl = viewerRef.current.pngURI()
    const link = document.createElement('a')
    link.download = 'ligand_3d.png'
    link.href = dataUrl
    link.click()
  }, [])

  if (!smiles) return null

  return (
    <div
      className="relative overflow-hidden rounded-xl bg-[#0e1628] border border-bx-surface/30"
      style={{ height, minHeight: height }}
    >
      {/* Header bar with controls */}
      <div className="absolute top-0 left-0 right-0 z-10 px-3 py-1.5 border-b border-white/5 flex items-center gap-2 bg-[#0e1628]/80 backdrop-blur-sm">
        <div className="w-2 h-2 rounded-full bg-red-500" />
        <div className="w-2 h-2 rounded-full bg-amber-400" />
        <div className="w-2 h-2 rounded-full bg-green-400" />
        <span className="ml-1 text-[10px] font-semibold text-white/50 uppercase tracking-wide">
          3D Ligand
        </span>

        <div className="ml-auto flex items-center gap-1">
          {/* Style pills */}
          <div className="flex items-center gap-0.5 bg-black/30 rounded-md p-0.5">
            {LIGAND_STYLES.map(s => (
              <button
                key={s.value}
                onClick={() => setDisplayStyle(s.value)}
                className={`px-1.5 py-0.5 text-[9px] rounded transition-all font-medium whitespace-nowrap ${
                  displayStyle === s.value
                    ? 'bg-bx-mint text-white shadow-sm'
                    : 'text-gray-400 hover:text-white hover:bg-white/10'
                }`}
              >
                {s.label}
              </button>
            ))}
          </div>

          {/* Spin toggle */}
          <button
            onClick={() => setSpinning(v => !v)}
            title={spinning ? 'Stop rotation' : 'Auto-rotate'}
            className={`p-1 rounded-md transition-colors ${spinning ? 'bg-blue-600 text-white' : 'text-gray-400 hover:text-white hover:bg-white/10'}`}
          >
            <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
            </svg>
          </button>

          {/* Screenshot */}
          <button
            onClick={handleScreenshot}
            title="Screenshot"
            className="p-1 rounded-md text-gray-400 hover:text-white hover:bg-white/10 transition-colors"
          >
            <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M3 9a2 2 0 012-2h.93a2 2 0 001.664-.89l.812-1.22A2 2 0 0110.07 4h3.86a2 2 0 011.664.89l.812 1.22A2 2 0 0018.07 7H19a2 2 0 012 2v9a2 2 0 01-2 2H5a2 2 0 01-2-2V9z" />
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 13a3 3 0 11-6 0 3 3 0 016 0z" />
            </svg>
          </button>
        </div>
      </div>

      {/* Loading spinner */}
      {loading && (
        <div className="absolute inset-0 flex items-center justify-center z-20">
          <div className="flex flex-col items-center gap-2">
            <div className="w-6 h-6 border-2 border-bx-mint/30 border-t-bx-mint rounded-full animate-spin" />
            <span className="text-[10px] text-white/40">Generating 3D...</span>
          </div>
        </div>
      )}

      {/* Error state */}
      {error && !loading && (
        <div className="absolute inset-0 flex flex-col items-center justify-center z-20 px-4">
          <svg className="w-8 h-8 text-red-400/60 mb-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
              d="M12 9v3.75m9-.75a9 9 0 11-18 0 9 9 0 0118 0zm-9 3.75h.008v.008H12v-.008z" />
          </svg>
          <p className="text-[10px] text-white/40 text-center">{error}</p>
        </div>
      )}

      {/* 3Dmol container */}
      <div
        ref={containerRef}
        style={{ width: '100%', height: '100%', position: 'relative' }}
      />
    </div>
  )
}
