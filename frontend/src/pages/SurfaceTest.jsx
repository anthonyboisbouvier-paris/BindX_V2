import React, { useEffect, useRef, useState } from 'react'

const PDB_ID = '1ycr' // small protein, MDM2 + p53 peptide

export default function SurfaceTest() {
  const containerRef = useRef(null)
  const [log, setLog] = useState([])
  const [viewer, setViewer] = useState(null)

  const addLog = (msg) => setLog((prev) => [...prev, `${new Date().toLocaleTimeString()} — ${msg}`])

  // Load 3Dmol.js
  useEffect(() => {
    const script = document.createElement('script')
    script.src = 'https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js'
    script.onload = () => {
      addLog('3Dmol.js loaded OK')
      initViewer()
    }
    script.onerror = () => addLog('ERROR: Failed to load 3Dmol.js')
    document.head.appendChild(script)
  }, [])

  const initViewer = async () => {
    const $3Dmol = window.$3Dmol
    if (!$3Dmol) { addLog('ERROR: $3Dmol not available'); return }

    addLog(`3Dmol version: ${$3Dmol.version || 'unknown'}`)
    addLog(`SurfaceType enum: ${JSON.stringify($3Dmol.SurfaceType)}`)

    const v = $3Dmol.createViewer(containerRef.current, {
      backgroundColor: '#1a1a2e',
      antialias: true,
    })

    addLog(`Fetching PDB ${PDB_ID}...`)
    try {
      const resp = await fetch(`https://files.rcsb.org/view/${PDB_ID}.pdb`)
      const pdbData = await resp.text()
      addLog(`PDB fetched: ${pdbData.length} chars`)

      const model = v.addModel(pdbData, 'pdb')
      addLog(`Model added — atoms: ${model.selectedAtoms({}).length}`)

      // Cartoon base style
      v.setStyle({ model }, { cartoon: { color: '#4a9eff', opacity: 0.3 } })
      v.zoomTo()
      v.render()
      addLog('Cartoon rendered OK')

      setViewer({ v, model, $3Dmol })
    } catch (err) {
      addLog(`ERROR fetching PDB: ${err.message}`)
    }
  }

  // === SECTION 1: Surface type tests ===

  const testProteinSurface = () => {
    if (!viewer) return
    const { v, model, $3Dmol } = viewer
    v.removeAllSurfaces()
    addLog('--- PROTEIN surface (VDW + model) ---')
    v.addSurface($3Dmol.SurfaceType.VDW,
      { opacity: 0.85, color: '#4a9eff' },
      { model },
      undefined, undefined,
      () => { addLog('OK — callback fired'); v.render() }
    )
  }

  const testSES = () => {
    if (!viewer) return
    const { v, model, $3Dmol } = viewer
    v.removeAllSurfaces()
    addLog('--- SES surface ---')
    const surfType = $3Dmol.SurfaceType.SES
    addLog(`SurfaceType.SES = ${surfType}`)
    if (surfType === undefined) {
      addLog('SES UNDEFINED — falling back to MS')
      v.addSurface($3Dmol.SurfaceType.MS,
        { opacity: 0.85, color: '#00e6a0' },
        { model }, undefined, undefined,
        () => { addLog('MS callback fired'); v.render() }
      )
    } else {
      v.addSurface(surfType,
        { opacity: 0.85, color: '#00e6a0' },
        { model }, undefined, undefined,
        () => { addLog('SES callback fired'); v.render() }
      )
    }
  }

  const testNoModel = () => {
    if (!viewer) return
    const { v, $3Dmol } = viewer
    v.removeAllSurfaces()
    addLog('--- Surface WITHOUT model selector (atomsel={}) ---')
    v.addSurface($3Dmol.SurfaceType.VDW,
      { opacity: 0.85, color: '#f59e0b' },
      {},
      undefined, undefined,
      () => { addLog('No-model callback fired'); v.render() }
    )
  }

  const testNumericType = () => {
    if (!viewer) return
    const { v, model } = viewer
    v.removeAllSurfaces()
    addLog('--- Numeric type=1 (VDW) ---')
    v.addSurface(1,
      { opacity: 0.85, color: '#8b5cf6' },
      { model }, undefined, undefined,
      () => { addLog('Numeric VDW callback fired'); v.render() }
    )
  }

  const testPocketSurface = () => {
    if (!viewer) return
    const { v, model, $3Dmol } = viewer
    v.removeAllSurfaces()
    addLog('--- POCKET surface (resi subset) ---')
    const pocketResis = [25, 26, 27, 28, 29, 50, 51, 52, 53, 54]
    addLog(`Pocket residues: ${pocketResis.join(',')}`)
    const sel = { resi: pocketResis, model }
    v.addStyle(sel, { stick: { color: '#f43f5e', radius: 0.12 } })
    v.addSurface($3Dmol.SurfaceType.VDW,
      { opacity: 0.7, color: '#f43f5e' },
      sel, undefined, undefined,
      () => { addLog('Pocket surface callback fired'); v.render() }
    )
    v.render()
  }

  // === SECTION 2: Simulating ProteinViewer conditions ===

  const testCartoonUnderSurface = () => {
    if (!viewer) return
    const { v, model, $3Dmol } = viewer
    v.removeAllSurfaces()
    addLog('--- Cartoon (0.15 opacity) + Surface overlay ---')
    // This is exactly what ProteinViewer does
    v.setStyle({ model }, { cartoon: { color: '#888', opacity: 0.15 } })
    v.addSurface($3Dmol.SurfaceType.VDW,
      { opacity: 0.9, color: '#4a9eff' },
      { model }, undefined, undefined,
      () => { addLog('Surface over cartoon callback fired'); v.render() }
    )
    v.render()
  }

  const testCartoonThenPocketSurface = () => {
    if (!viewer) return
    const { v, model, $3Dmol } = viewer
    v.removeAllSurfaces()
    addLog('--- Cartoon + Pocket sticks + Pocket surface ---')
    // Protein as cartoon
    v.setStyle({ model }, { cartoon: { color: '#4a9eff' } })
    // Pocket sticks + surface
    const pocketResis = [25, 26, 27, 28, 29, 50, 51, 52, 53, 54]
    const sel = { resi: pocketResis, model }
    v.addStyle(sel, { stick: { color: '#f43f5e', radius: 0.1, opacity: 0.4 } })
    v.addSurface($3Dmol.SurfaceType.VDW,
      { opacity: 0.7, color: '#f43f5e' },
      sel, undefined, undefined,
      () => { addLog('Pocket surface callback fired'); v.render() }
    )
    v.render()
  }

  const testMultipleSurfaces = () => {
    if (!viewer) return
    const { v, model, $3Dmol } = viewer
    v.removeAllSurfaces()
    addLog('--- Multiple surfaces (protein + pocket) ---')
    v.setStyle({ model }, { cartoon: { color: '#888', opacity: 0.15 } })
    // Full protein surface (transparent)
    v.addSurface($3Dmol.SurfaceType.VDW,
      { opacity: 0.3, color: '#4a9eff' },
      { model }, undefined, undefined,
      () => { addLog('Protein surface callback'); v.render() }
    )
    // Pocket surface (opaque)
    const pocketResis = [25, 26, 27, 28, 29, 50, 51, 52, 53, 54]
    v.addSurface($3Dmol.SurfaceType.VDW,
      { opacity: 0.85, color: '#f43f5e' },
      { resi: pocketResis, model }, undefined, undefined,
      () => { addLog('Pocket surface callback'); v.render() }
    )
    v.render()
  }

  const testColorschemeAminoSurface = () => {
    if (!viewer) return
    const { v, model, $3Dmol } = viewer
    v.removeAllSurfaces()
    addLog('--- Surface with colorscheme: "amino" (hydrophobicity) ---')
    v.setStyle({ model }, { cartoon: { color: '#888', opacity: 0.15 } })
    v.addSurface($3Dmol.SurfaceType.VDW,
      { opacity: 0.9, colorscheme: 'amino' },
      { model }, undefined, undefined,
      () => { addLog('Amino colorscheme surface callback'); v.render() }
    )
    v.render()
  }

  const testSpectrumSurface = () => {
    if (!viewer) return
    const { v, model, $3Dmol } = viewer
    v.removeAllSurfaces()
    addLog('--- Surface with color: "spectrum" ---')
    v.setStyle({ model }, { cartoon: { color: '#888', opacity: 0.15 } })
    v.addSurface($3Dmol.SurfaceType.VDW,
      { opacity: 0.9, color: 'spectrum' },
      { model }, undefined, undefined,
      () => { addLog('Spectrum surface callback'); v.render() }
    )
    v.render()
  }

  const testBfactorSurface = () => {
    if (!viewer) return
    const { v, model, $3Dmol } = viewer
    v.removeAllSurfaces()
    addLog('--- Surface with bfactor colorscheme ---')
    v.setStyle({ model }, { cartoon: { color: '#888', opacity: 0.15 } })
    v.addSurface($3Dmol.SurfaceType.VDW,
      { opacity: 0.9, colorscheme: { prop: 'b', gradient: 'rwb', min: 0, max: 100 } },
      { model }, undefined, undefined,
      () => { addLog('Bfactor surface callback'); v.render() }
    )
    v.render()
  }

  // === SECTION 3: Edge cases / debugging ===

  const testPocketNoModelSelector = () => {
    if (!viewer) return
    const { v, model, $3Dmol } = viewer
    v.removeAllSurfaces()
    addLog('--- Pocket surface WITHOUT model in selector ---')
    v.setStyle({ model }, { cartoon: { color: '#4a9eff' } })
    const pocketResis = [25, 26, 27, 28, 29, 50, 51, 52, 53, 54]
    // Intentionally omit model from surface selector
    v.addStyle({ resi: pocketResis, model }, { stick: { color: '#f43f5e', radius: 0.1 } })
    v.addSurface($3Dmol.SurfaceType.VDW,
      { opacity: 0.7, color: '#f43f5e' },
      { resi: pocketResis },  // NO model
      undefined, undefined,
      () => { addLog('Pocket (no model) callback fired'); v.render() }
    )
    v.render()
    addLog('Compare this with the "Pocket Surface" button above')
  }

  const testLowOpacity = () => {
    if (!viewer) return
    const { v, model, $3Dmol } = viewer
    v.removeAllSurfaces()
    addLog('--- Surface opacity=0.08 (very transparent, like old cartoon) ---')
    v.setStyle({ model }, { cartoon: { color: '#4a9eff', opacity: 0.08 } })
    v.addSurface($3Dmol.SurfaceType.VDW,
      { opacity: 0.08, color: '#4a9eff' },
      { model }, undefined, undefined,
      () => { addLog('Low opacity surface callback'); v.render() }
    )
    v.render()
    addLog('Should be barely visible — this shows what happens with wrong opacity')
  }

  const resetCartoon = () => {
    if (!viewer) return
    const { v, model } = viewer
    v.removeAllSurfaces()
    v.removeAllShapes()
    v.removeAllLabels()
    v.setStyle({ model }, { cartoon: { color: '#4a9eff' } })
    v.render()
    addLog('Reset to cartoon')
  }

  const clearSurfaces = () => {
    if (!viewer) return
    viewer.v.removeAllSurfaces()
    viewer.v.render()
    addLog('All surfaces cleared')
  }

  const btnClass = (color) =>
    `px-3 py-1.5 rounded text-xs font-medium transition-colors ${color}`

  return (
    <div className="min-h-screen bg-[#0f131d] text-white p-6">
      <h1 className="text-2xl font-bold mb-2">3Dmol Surface Test</h1>
      <p className="text-gray-400 mb-4">Page de test isolée — PDB: {PDB_ID} (MDM2)</p>

      <div className="flex gap-4">
        {/* 3D Viewer */}
        <div
          ref={containerRef}
          style={{ width: 600, height: 550, position: 'relative' }}
          className="border border-gray-600 rounded-lg flex-shrink-0"
        />

        {/* Controls + Log */}
        <div className="flex-1 flex flex-col gap-3 min-w-0">
          {/* Section 1: Basic surface types */}
          <div>
            <div className="text-xs text-gray-500 uppercase tracking-wider mb-1">Surface Types</div>
            <div className="flex flex-wrap gap-1.5">
              <button onClick={testProteinSurface} className={btnClass('bg-blue-600 hover:bg-blue-500')}>
                VDW + model
              </button>
              <button onClick={testSES} className={btnClass('bg-emerald-600 hover:bg-emerald-500')}>
                SES / MS
              </button>
              <button onClick={testNoModel} className={btnClass('bg-amber-600 hover:bg-amber-500')}>
                Sans model
              </button>
              <button onClick={testNumericType} className={btnClass('bg-violet-600 hover:bg-violet-500')}>
                Numeric type=1
              </button>
              <button onClick={testPocketSurface} className={btnClass('bg-rose-600 hover:bg-rose-500')}>
                Pocket (resi subset)
              </button>
            </div>
          </div>

          {/* Section 2: ProteinViewer simulation */}
          <div>
            <div className="text-xs text-gray-500 uppercase tracking-wider mb-1">Simulation ProteinViewer</div>
            <div className="flex flex-wrap gap-1.5">
              <button onClick={testCartoonUnderSurface} className={btnClass('bg-sky-700 hover:bg-sky-600')}>
                Cartoon + Surface
              </button>
              <button onClick={testCartoonThenPocketSurface} className={btnClass('bg-pink-700 hover:bg-pink-600')}>
                Cartoon + Pocket Surf
              </button>
              <button onClick={testMultipleSurfaces} className={btnClass('bg-indigo-700 hover:bg-indigo-600')}>
                Protein + Pocket Surf
              </button>
              <button onClick={testColorschemeAminoSurface} className={btnClass('bg-teal-700 hover:bg-teal-600')}>
                Amino colorscheme
              </button>
              <button onClick={testSpectrumSurface} className={btnClass('bg-cyan-700 hover:bg-cyan-600')}>
                Spectrum color
              </button>
              <button onClick={testBfactorSurface} className={btnClass('bg-orange-700 hover:bg-orange-600')}>
                Bfactor gradient
              </button>
            </div>
          </div>

          {/* Section 3: Edge cases */}
          <div>
            <div className="text-xs text-gray-500 uppercase tracking-wider mb-1">Edge Cases / Debug</div>
            <div className="flex flex-wrap gap-1.5">
              <button onClick={testPocketNoModelSelector} className={btnClass('bg-red-800 hover:bg-red-700')}>
                Pocket SANS model
              </button>
              <button onClick={testLowOpacity} className={btnClass('bg-gray-700 hover:bg-gray-600')}>
                Opacity 0.08
              </button>
              <button onClick={resetCartoon} className={btnClass('bg-gray-600 hover:bg-gray-500')}>
                Reset Cartoon
              </button>
              <button onClick={clearSurfaces} className={btnClass('bg-gray-600 hover:bg-gray-500')}>
                Clear Surfaces
              </button>
            </div>
          </div>

          {/* Log */}
          <div className="flex-1 bg-[#141925] rounded-lg p-3 overflow-y-auto max-h-[380px] font-mono text-xs">
            {log.map((l, i) => (
              <div key={i} className={l.includes('ERROR') ? 'text-red-400' : l.includes('OK') ? 'text-green-400' : 'text-gray-300'}>
                {l}
              </div>
            ))}
            {log.length === 0 && <div className="text-gray-500">Loading...</div>}
          </div>
        </div>
      </div>
    </div>
  )
}
