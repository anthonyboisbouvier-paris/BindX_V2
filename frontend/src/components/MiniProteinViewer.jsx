import React, { useEffect, useRef, useState } from 'react'

// Shared 3Dmol CDN loader (same as ProteinViewer)
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
// MiniProteinViewer — tiny 3D spinning protein for project cards
// Props:
//   pdbUrl    — URL to fetch PDB file
//   width     — viewer width (default 100%)
//   height    — viewer height (default 120)
// ---------------------------------------------------------------------------
export default function MiniProteinViewer({ pdbUrl, width = '100%', height = 120 }) {
  const containerRef = useRef(null)
  const viewerRef = useRef(null)
  const animRef = useRef(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState(false)

  useEffect(() => {
    if (!pdbUrl || !containerRef.current) return
    let cancelled = false

    async function init() {
      try {
        const $3Dmol = await load3Dmol()
        if (cancelled) return

        const res = await fetch(pdbUrl)
        if (!res.ok) throw new Error(`HTTP ${res.status}`)
        const pdbData = await res.text()
        if (cancelled) return

        // Clear container
        containerRef.current.innerHTML = ''

        const viewer = $3Dmol.createViewer(containerRef.current, {
          backgroundColor: 'white',
          antialias: true,
          cartoonQuality: 4,
        })
        viewerRef.current = viewer

        viewer.addModel(pdbData, 'pdb')
        viewer.setStyle({}, {
          cartoon: {
            color: 'spectrum',
            opacity: 0.95,
          },
        })
        viewer.zoomTo()
        viewer.render()
        setLoading(false)

        // Slow auto-rotation
        let angle = 0
        function spin() {
          if (cancelled) return
          angle += 0.3
          viewer.rotate(0.3, 'y')
          viewer.render()
          animRef.current = requestAnimationFrame(spin)
        }
        animRef.current = requestAnimationFrame(spin)
      } catch (err) {
        if (!cancelled) {
          console.warn('[MiniProteinViewer] Load failed:', err.message)
          setError(true)
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
  }, [pdbUrl])

  if (!pdbUrl) return null
  if (error) return null

  return (
    <div
      className="relative overflow-hidden rounded-lg"
      style={{ width, height, minHeight: height }}
    >
      {loading && (
        <div className="absolute inset-0 flex items-center justify-center bg-gray-50">
          <div className="w-4 h-4 border-2 border-bx-mint/30 border-t-bx-mint rounded-full animate-spin" />
        </div>
      )}
      <div
        ref={containerRef}
        style={{ width: '100%', height: '100%', position: 'relative' }}
      />
    </div>
  )
}
