import React, { useEffect, useRef, useState } from 'react'
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

// ---------------------------------------------------------------------------
// LigandViewer3D — 3D stick+sphere molecule viewer from SMILES
// Props:
//   smiles  — SMILES string to render
//   height  — viewer height (default 280)
// ---------------------------------------------------------------------------
export default function LigandViewer3D({ smiles, height = 280 }) {
  const containerRef = useRef(null)
  const viewerRef = useRef(null)
  const animRef = useRef(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState(null)

  useEffect(() => {
    if (!smiles || !containerRef.current) return
    let cancelled = false

    async function init() {
      setLoading(true)
      setError(null)

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
        viewer.setStyle({}, {
          stick: { radius: 0.15, colorscheme: 'default' },
          sphere: { scale: 0.25, colorscheme: 'default' },
        })
        viewer.zoomTo()
        viewer.render()
        setLoading(false)

        // Slow auto-rotation
        function spin() {
          if (cancelled) return
          viewer.rotate(0.3, 'y')
          viewer.render()
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

  if (!smiles) return null

  return (
    <div
      className="relative overflow-hidden rounded-xl bg-[#0e1628] border border-bx-surface/30"
      style={{ height, minHeight: height }}
    >
      {/* Header bar */}
      <div className="absolute top-0 left-0 right-0 z-10 px-3 py-2 border-b border-white/5 flex items-center gap-2 bg-[#0e1628]/80 backdrop-blur-sm">
        <div className="w-2 h-2 rounded-full bg-red-500" />
        <div className="w-2 h-2 rounded-full bg-amber-400" />
        <div className="w-2 h-2 rounded-full bg-green-400" />
        <span className="ml-1 text-[10px] font-semibold text-white/50 uppercase tracking-wide">
          3D Ligand
        </span>
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
