import React from 'react'

// ---------------------------------------------------------------------------
// Status icons
// ---------------------------------------------------------------------------
function Ok() { return <span className="text-emerald-500 font-bold shrink-0 text-sm leading-none">&#10003;</span> }
function Fail() { return <span className="text-red-400 shrink-0 text-sm leading-none">&#10007;</span> }
function Warning() { return <span className="text-amber-500 shrink-0 text-sm leading-none">&#9888;</span> }

// ---------------------------------------------------------------------------
// Row — one diagnostic line
// ---------------------------------------------------------------------------
function Row({ icon, children }) {
  return <div className="flex items-start gap-1.5 py-0.5">{icon}{children}</div>
}

// ---------------------------------------------------------------------------
// ViewerDiagnostic — live data from renderer + static metadata
// ---------------------------------------------------------------------------
export default function ViewerDiagnostic({ pdbUrl, selectedPocket, uniprotFeatures, selectedStructure, live }) {
  const activeSites  = uniprotFeatures?.activeSites  || []
  const bindingSites = uniprotFeatures?.bindingSites || []
  const domains      = uniprotFeatures?.domains      || []

  const isAlphaFold = pdbUrl?.includes('alphafold') || selectedStructure?.source === 'alphafold'
  const structureLabel = isAlphaFold ? 'AlphaFold' : 'RCSB PDB'
  const pdbId = selectedStructure?.pdb_id
    || (pdbUrl ? (pdbUrl.match(/([A-Z0-9]{4})\.pdb/i)?.[1] || pdbUrl.match(/\/([A-Z0-9]{4})\b/i)?.[1]) : null)

  const hasWarnings = live?.warnings?.length > 0

  return (
    <div className="card overflow-hidden">
      {/* Header */}
      <div className="px-4 py-2.5 border-b border-bx-light-border bg-bx-light-warm flex items-center justify-between">
        <h3 className="text-sm font-semibold text-bx-light-text flex items-center gap-2">
          <svg className="w-4 h-4 text-gray-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8} d="M9.75 3.104v5.714a2.25 2.25 0 01-.659 1.591L5 14.5M9.75 3.104c-.251.023-.501.05-.75.082m.75-.082a24.301 24.301 0 014.5 0m0 0v5.714c0 .597.237 1.17.659 1.591L19 14.5M14.25 3.104c.251.023.501.05.75.082M19 14.5l-1.47 4.9a2.25 2.25 0 01-2.156 1.6H8.626a2.25 2.25 0 01-2.156-1.6L5 14.5m14 0H5" />
          </svg>
          Viewer Diagnostic
        </h3>
        {hasWarnings && (
          <span className="text-sm bg-amber-50 text-amber-700 px-2 py-0.5 rounded-full font-medium">
            {live.warnings.length} warning{live.warnings.length > 1 ? 's' : ''}
          </span>
        )}
      </div>

      {/* Scrollable content */}
      <div className="px-4 py-3 space-y-3 text-sm text-gray-500 leading-relaxed overflow-y-auto" style={{ maxHeight: 360 }}>

        {/* ── Warnings (top, prominent) ── */}
        {hasWarnings && (
          <div className="space-y-1.5">
            {live.warnings.map((w, i) => (
              <div key={i} className="flex items-start gap-2 bg-amber-50 border border-amber-200 rounded-md px-3 py-2 text-amber-800">
                <Warning />
                <span>{w}</span>
              </div>
            ))}
          </div>
        )}

        {/* ── Structure ── */}
        <section>
          <h4 className="text-gray-700 font-semibold mb-1 text-sm uppercase tracking-wide">Structure</h4>
          <div className="space-y-0.5 ml-0.5">
            <Row icon={<Ok />}>
              <span>
                <span className="text-gray-700 font-medium">{structureLabel}</span>
                {pdbId && <> &mdash; <span className="font-mono text-gray-600">{pdbId}</span></>}
                {selectedStructure?.resolution && <> &mdash; {selectedStructure.resolution.toFixed(2)} &#8491;</>}
                {selectedStructure?.method && ` (${selectedStructure.method})`}
              </span>
            </Row>
            {live?.structureAtoms > 0 && (
              <Row icon={<Ok />}>
                <span><span className="font-mono text-gray-700">{live.structureAtoms.toLocaleString()}</span> atoms loaded in model</span>
              </Row>
            )}
            <p className="text-gray-400 mt-0.5">
              {isAlphaFold
                ? 'Predicted structure (DeepMind). Confidence = pLDDT in B-factor. No co-crystallized ligand. Disordered regions (pLDDT < 50) should be interpreted with caution.'
                : 'Experimental structure. Co-crystallized ligands may be present as HETATM records.'}
            </p>
          </div>
        </section>

        {/* ── Pocket ── */}
        <section>
          <h4 className="text-gray-700 font-semibold mb-1 text-sm uppercase tracking-wide">Pocket</h4>
          <div className="space-y-0.5 ml-0.5">
            {selectedPocket ? (<>
              {/* Detection method */}
              <Row icon={live?.pocketAtomsMatched > 0 ? <Ok /> : <Fail />}>
                <span>
                  <span className="text-gray-700 font-medium">
                    {selectedPocket.method === 'co-crystallized_ligand' ? 'Co-crystallized ligand' : selectedPocket.method === 'p2rank' ? 'P2Rank ML' : selectedPocket.method || '?'}
                  </span>
                  {selectedPocket.probability != null && <> &mdash; <span className="font-mono text-gray-700">{Math.round(selectedPocket.probability * 100)}%</span> confidence</>}
                  {selectedPocket.ligand_id && <> &mdash; ligand <span className="font-mono text-amber-700">{selectedPocket.ligand_id}</span></>}
                </span>
              </Row>

              {/* Live atom matching */}
              {live && (
                <Row icon={live.pocketAtomsMatched > 0 ? <Ok /> : <Fail />}>
                  <span>
                    <span className="font-mono text-gray-700">{live.pocketResiduesParsed}</span>/{live.pocketResiduesRaw} residues parsed
                    {' '}&rarr; <span className={`font-mono ${live.pocketAtomsMatched > 0 ? 'text-gray-700' : 'text-red-600 font-medium'}`}>{live.pocketAtomsMatched}</span> atoms matched in model
                    {live.pocketAtomsMatched === 0 && live.pocketResiduesParsed > 0 && (
                      <span className="text-red-500 font-medium"> &mdash; Surface and Stick modes will not render (no matching atoms). Sphere mode still works if center is available.</span>
                    )}
                  </span>
                </Row>
              )}

              {/* Surface rendering status */}
              {live?.pocketSurfaceAttempted && (
                <Row icon={live.pocketAtomsMatched > 20 ? <Ok /> : live.pocketAtomsMatched > 0 ? <Warning /> : <Fail />}>
                  <span>
                    SES surface: <span className="font-mono text-gray-700">{live.pocketAtomsMatched}</span> atoms in selection
                    {live.pocketAtomsMatched === 0 && ' — surface cannot render (0 atoms)'}
                    {live.pocketAtomsMatched > 0 && live.pocketAtomsMatched < 20 && ' — very few atoms, surface may appear small or incomplete'}
                    {live.pocketAtomsMatched >= 20 && ' — OK'}
                  </span>
                </Row>
              )}

              {/* Sphere mode status */}
              {live?.pocketRenderMode === 'sphere+stick' && (
                <Row icon={selectedPocket.center ? <Ok /> : <Fail />}>
                  <span>
                    Sphere mode: {selectedPocket.center
                      ? 'uses center coordinates (independent of residue matching)'
                      : 'no center coordinates available — sphere cannot render'}
                  </span>
                </Row>
              )}

              {selectedPocket.center && (
                <Row icon={<Ok />}>
                  <span>Center: <span className="font-mono text-gray-600">({selectedPocket.center.map(c => c.toFixed(1)).join(', ')})</span> &#8491;</span>
                </Row>
              )}

              {/* Method explanation */}
              <p className="text-gray-400 mt-1">
                {selectedPocket.method === 'co-crystallized_ligand'
                  ? `Ligand ${selectedPocket.ligand_id || '?'} is co-crystallized in the PDB. All protein residues within 6 \u212B of any ligand atom form the pocket. Experimentally validated binding site.`
                  : selectedPocket.method === 'p2rank'
                    ? 'P2Rank: ML model (random forest) trained on known holo structures. Predicts pockets by analyzing surface geometry, curvature, hydrophobicity, and solvent accessibility.'
                    : `Detection method: ${selectedPocket.method || 'unknown'}.`
                }
              </p>

              {/* Rendering mode explanation */}
              <p className="text-gray-400 border-l-2 border-gray-200 pl-2 mt-1">
                {live?.pocketRenderMode === 'surface'
                  ? 'Surface mode: SES (Solvent Excluded Surface, 1.4 \u212B probe) computed around pocket residues. Standard in drug discovery. Requires enough atoms in the selection to form a coherent surface.'
                  : live?.pocketRenderMode === 'sphere+stick'
                    ? 'Sphere mode: geometric sphere at the pocket barycenter. Does not depend on individual residue matching — only needs center coordinates.'
                    : 'Stick mode: ball+stick on pocket-lining residues. Each residue is rendered individually.'}
              </p>
            </>) : (
              <Row icon={<Fail />}>
                <span>No pocket selected. Choose one from the list above.</span>
              </Row>
            )}
          </div>
        </section>

        {/* ── Ligand ── */}
        <section>
          <h4 className="text-gray-700 font-semibold mb-1 text-sm uppercase tracking-wide">Ligand</h4>
          <div className="space-y-0.5 ml-0.5">
            {live ? (
              <Row icon={live.ligandAtomsFound > 0 ? <Ok /> : (isAlphaFold ? <Fail /> : <Warning />)}>
                <span>
                  {live.ligandAtomsFound > 0
                    ? <><span className="font-mono text-gray-700">{live.ligandAtomsFound}</span> ligand atoms detected (HETATM, excluding solvents)</>
                    : isAlphaFold
                      ? 'No ligand — AlphaFold structures are apo (ligand-free predictions)'
                      : 'No ligand atoms found. Structure may be apo or ligand was stripped during deposition.'
                  }
                </span>
              </Row>
            ) : (
              <p className="text-gray-400">Ligands auto-detected from HETATM records. Solvents (HOH, SO4, GOL, etc.) filtered out.</p>
            )}
          </div>
        </section>

        {/* ── UniProt Annotations ── */}
        <section>
          <h4 className="text-gray-700 font-semibold mb-1 text-sm uppercase tracking-wide">UniProt Annotations</h4>
          <div className="space-y-0.5 ml-0.5">
            {uniprotFeatures ? (<>
              <Row icon={activeSites.length > 0 ? (live?.activeSitesResolved > 0 ? <Ok /> : <Warning />) : <Fail />}>
                <span>
                  Active Sites: <span className="text-gray-700 font-medium">{activeSites.length}</span> annotated
                  {live?.activeSitesResolved > 0 && <> &mdash; <span className="font-mono text-gray-700">{live.activeSitesResolved}</span> atoms resolved in model</>}
                  {activeSites.length > 0 && live?.activeSitesResolved === 0 && <span className="text-amber-600"> &mdash; 0 atoms matched (numbering mismatch?)</span>}
                  {activeSites.length === 0 && <span className="text-gray-400"> (non-enzymatic or not annotated)</span>}
                </span>
              </Row>
              <Row icon={bindingSites.length > 0 ? (live?.bindingSitesResolved > 0 ? <Ok /> : <Warning />) : <Fail />}>
                <span>
                  Binding Sites: <span className="text-gray-700 font-medium">{bindingSites.length}</span> annotated
                  {live?.bindingSitesResolved > 0 && <> &mdash; <span className="font-mono text-gray-700">{live.bindingSitesResolved}</span> atoms resolved</>}
                  {bindingSites.length > 0 && live?.bindingSitesResolved === 0 && <span className="text-amber-600"> &mdash; 0 atoms matched</span>}
                  {bindingSites.length === 0 && <span className="text-gray-400"> (not annotated)</span>}
                </span>
              </Row>
              <Row icon={domains.length > 0 ? <Ok /> : <Fail />}>
                <span>
                  Domains: <span className="text-gray-700 font-medium">{domains.length}</span>
                  {domains.length > 0 && (<>
                    {' '}&mdash;{' '}
                    {domains.map((d, i) => (
                      <span key={i}>
                        {i > 0 && ', '}
                        <span className="text-gray-600">{d.description || d.type || '?'}</span>
                        <span className="text-gray-400"> ({d.start}&ndash;{d.end})</span>
                      </span>
                    ))}
                  </>)}
                  {domains.length === 0 && <span className="text-gray-400"> (not annotated)</span>}
                </span>
              </Row>
              <p className="text-gray-400 mt-0.5">Source: UniProt experimental annotations. Domains via InterPro / Pfam.</p>
            </>) : (
              <Row icon={<Fail />}>
                <span>No UniProt data. Annotations require target configuration via UniProt ID.</span>
              </Row>
            )}
          </div>
        </section>
      </div>
    </div>
  )
}
