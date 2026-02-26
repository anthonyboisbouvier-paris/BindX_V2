import React, { useState } from 'react'
import ADMETRadar from './ADMETRadar.jsx'
import SynthesisTree from './SynthesisTree.jsx'
import InfoTip from './InfoTip.jsx'
import InteractionView from './InteractionView.jsx'
import InteractionDiagram from './InteractionDiagram.jsx'

function PropRow({ label, value, unit = '', highlight = false, tooltip }) {
  return (
    <div className={`flex items-center justify-between py-1.5 px-3 rounded-lg ${highlight ? 'bg-dockit-green/10' : 'hover:bg-gray-50'}`}>
      <span className="text-xs font-medium text-gray-500 flex items-center">
        {label}
        {tooltip && <InfoTip text={tooltip} />}
      </span>
      <span className={`text-sm font-semibold ${highlight ? 'text-dockit-green' : 'text-gray-800'}`}>
        {value !== null && value !== undefined ? `${typeof value === 'number' ? Number(value).toFixed(2) : value}${unit}` : 'N/A'}
      </span>
    </div>
  )
}

function LipinskiIndicator({ value, optimal }) {
  // Check if value respects Lipinski rule
  const ok = optimal(value)
  return (
    <span className={`inline-flex items-center gap-1 text-xs ${ok ? 'text-green-600' : 'text-red-500'}`}>
      {ok ? (
        <svg className="w-3 h-3" fill="currentColor" viewBox="0 0 20 20">
          <path fillRule="evenodd" d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z" clipRule="evenodd" />
        </svg>
      ) : (
        <svg className="w-3 h-3" fill="currentColor" viewBox="0 0 20 20">
          <path fillRule="evenodd" d="M4.293 4.293a1 1 0 011.414 0L10 8.586l4.293-4.293a1 1 0 111.414 1.414L11.414 10l4.293 4.293a1 1 0 01-1.414 1.414L10 11.414l-4.293 4.293a1 1 0 01-1.414-1.414L8.586 10 4.293 5.707a1 1 0 010-1.414z" clipRule="evenodd" />
        </svg>
      )}
    </span>
  )
}

function TabButton({ label, active, onClick }) {
  return (
    <button
      type="button"
      onClick={onClick}
      className={`px-3 py-1.5 text-xs font-semibold transition-all duration-150 focus:outline-none ${
        active
          ? 'text-dockit-blue border-b-2 border-dockit-blue'
          : 'text-gray-400 hover:text-gray-600 border-b-2 border-transparent'
      }`}
    >
      {label}
    </button>
  )
}

export default function MoleculeCard({ molecule, rank }) {
  const [activeTab, setActiveTab] = useState('docking')
  const [interactionView, setInteractionView] = useState('table') // 'table' | 'diagram'

  if (!molecule) {
    return (
      <div className="card p-6">
        <div className="flex flex-col items-center justify-center h-40 text-gray-300">
          <svg className="w-12 h-12 mb-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
          </svg>
          <p className="text-sm">Select a ligand from the table</p>
        </div>
      </div>
    )
  }

  const name = molecule.name || molecule.molecule_name || `Ligand ${rank + 1}`
  const hasSvg = Boolean(molecule.svg || molecule.svg_2d || molecule.structure_svg)
  const svgContent = molecule.svg || molecule.svg_2d || molecule.structure_svg

  // V2 data availability
  const hasAdmet = Boolean(molecule.admet)
  const hasSynthesis = Boolean(molecule.synthesis_route)
  // V5bis: interaction data
  const hasInteractions = Boolean(molecule.interactions)
  const hasV2 = hasAdmet || hasSynthesis || hasInteractions

  // Lipinski rules
  const lipinski = [
    { label: 'MW <= 500 Da',    check: (m) => (m.mw || 0) <= 500 },
    { label: 'LogP <= 5',       check: (m) => (m.logp || 0) <= 5 },
    { label: 'HBD <= 5',        check: (m) => (m.hbd || 0) <= 5 },
    { label: 'HBA <= 10',       check: (m) => (m.hba || 0) <= 10 },
  ]
  const lipinskiPass = lipinski.filter(r => r.check(molecule)).length

  return (
    <div className="card overflow-hidden">
      {/* Card header */}
      <div className="px-4 py-3 bg-dockit-blue flex items-center justify-between">
        <div className="flex items-center gap-2">
          <span className="w-7 h-7 bg-dockit-green text-white text-xs font-bold rounded-full flex items-center justify-center flex-shrink-0">
            {rank + 1}
          </span>
          <div>
            <div className="flex items-center gap-1.5">
              <h3 className="text-white font-semibold text-sm truncate max-w-[160px]" title={name}>
                {name}
              </h3>
              {/* Docking status badge */}
              {molecule.docking_status === 'docked' ? (
                <span className="flex-shrink-0 px-1.5 py-0.5 text-[10px] font-semibold rounded bg-green-100 text-green-700 uppercase">
                  {molecule.docking_engine || 'Docked'}
                </span>
              ) : molecule.docking_status === 'failed' ? (
                <span className="flex-shrink-0 px-1.5 py-0.5 text-[10px] font-semibold rounded bg-red-100 text-red-700">
                  Failed
                </span>
              ) : (
                <span className="flex-shrink-0 px-1.5 py-0.5 text-[10px] font-semibold rounded bg-gray-100 text-gray-500">
                  Not docked
                </span>
              )}
            </div>
            {molecule.chembl_id && (
              <p className="text-white/50 text-xs font-mono">{molecule.chembl_id}</p>
            )}
          </div>
        </div>
        {molecule.composite_score !== undefined && (
          <div className="text-right">
            <div className="text-dockit-green font-bold text-lg">
              {Number(molecule.composite_score).toFixed(3)}
            </div>
            <div className="text-white/50 text-xs flex items-center justify-end">
              composite score
              <InfoTip text="Combined score (0-1) from binding affinity, drug-likeness (QED), and ADMET safety." />
            </div>
          </div>
        )}
      </div>

      {/* Tab navigation — only when V2 data is available */}
      {hasV2 && (
        <div className="flex items-center gap-0 px-4 pt-2 border-b border-gray-100 bg-white">
          <TabButton label="Docking" active={activeTab === 'docking'} onClick={() => setActiveTab('docking')} />
          {hasAdmet && (
            <TabButton label="ADMET" active={activeTab === 'admet'} onClick={() => setActiveTab('admet')} />
          )}
          {hasSynthesis && (
            <TabButton label="Synthesis" active={activeTab === 'synthesis'} onClick={() => setActiveTab('synthesis')} />
          )}
          {hasInteractions && (
            <TabButton label="Interactions" active={activeTab === 'interactions'} onClick={() => setActiveTab('interactions')} />
          )}
        </div>
      )}

      {/* Tab content */}

      {/* Tab: Docking (existing content — always visible when activeTab=docking OR no V2) */}
      {(!hasV2 || activeTab === 'docking') && (
        <>
          {/* 2D Structure */}
          <div className="p-4 border-b border-gray-100">
            <p className="text-xs font-medium text-gray-400 uppercase tracking-wide mb-2">2D Structure</p>
            <div className="bg-white rounded-lg border border-gray-100 flex items-center justify-center" style={{ minHeight: '160px' }}>
              {hasSvg ? (
                <div
                  className="w-full"
                  dangerouslySetInnerHTML={{ __html: svgContent }}
                  style={{ maxHeight: '200px', overflow: 'hidden' }}
                />
              ) : (
                <div className="flex flex-col items-center text-gray-200 py-8">
                  <svg className="w-10 h-10 mb-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M4 16l4.586-4.586a2 2 0 012.828 0L16 16m-2-2l1.586-1.586a2 2 0 012.828 0L20 14m-6-6h.01M6 20h12a2 2 0 002-2V6a2 2 0 00-2-2H6a2 2 0 00-2 2v12a2 2 0 002 2z" />
                  </svg>
                  <p className="text-xs">2D structure not available</p>
                </div>
              )}
            </div>

            {/* SMILES */}
            {molecule.smiles && (
              <div className="mt-2 p-2 bg-gray-50 rounded text-xs font-mono text-gray-400 truncate" title={molecule.smiles}>
                {molecule.smiles}
              </div>
            )}
          </div>

          {/* Properties */}
          <div className="p-4">
            <p className="text-xs font-medium text-gray-400 uppercase tracking-wide mb-3">Physicochemical Properties</p>
            <div className="space-y-0.5">
              <PropRow
                label="Binding Affinity"
                value={molecule.affinity}
                unit=" kcal/mol"
                highlight
                tooltip="Predicted binding energy (kcal/mol). More negative = stronger binding to the protein."
              />
              <PropRow
                label="Molecular Weight"
                value={molecule.mw ? Math.round(molecule.mw) : null}
                unit=" Da"
                tooltip="Molecular weight in Daltons. Lipinski rule: should be below 500 Da."
              />
              <PropRow
                label="LogP"
                value={molecule.logp}
                tooltip="Lipophilicity coefficient. Optimal range 0-5 for oral bioavailability."
              />
              <PropRow
                label="TPSA"
                value={molecule.tpsa}
                unit=" A²"
                tooltip="Topological Polar Surface Area. Values below 140 A² indicate good oral absorption."
              />
              <PropRow
                label="QED"
                value={molecule.qed}
                tooltip="Quantitative Estimate of Drug-likeness (0-1). Higher values indicate more drug-like properties."
              />
              {molecule.hbd !== undefined && (
                <PropRow
                  label="HBD (H-bond donors)"
                  value={molecule.hbd}
                  tooltip="Number of hydrogen bond donors. Lipinski rule: should be 5 or fewer."
                />
              )}
              {molecule.hba !== undefined && (
                <PropRow
                  label="HBA (H-bond acceptors)"
                  value={molecule.hba}
                  tooltip="Number of hydrogen bond acceptors. Lipinski rule: should be 10 or fewer."
                />
              )}
              {molecule.rotatable_bonds !== undefined && (
                <PropRow
                  label="Rotatable Bonds"
                  value={molecule.rotatable_bonds}
                  tooltip="Number of rotatable bonds. Fewer bonds generally correlate with better oral absorption."
                />
              )}
            </div>

            {/* Lipinski rule of 5 */}
            <div className="mt-4 p-3 bg-gray-50 rounded-lg">
              <div className="flex items-center justify-between mb-2">
                <p className="text-xs font-semibold text-gray-600 flex items-center">
                  Lipinski Rule of 5
                  <InfoTip text="The Rule of 5 predicts oral bioavailability. Molecules passing all 4 criteria are more likely to be orally active drugs." />
                </p>
                <span className={`text-xs font-bold px-2 py-0.5 rounded-full ${
                  lipinskiPass === 4 ? 'bg-green-100 text-green-700' :
                  lipinskiPass >= 3 ? 'bg-yellow-100 text-yellow-700' :
                  'bg-red-100 text-red-700'
                }`}>
                  {lipinskiPass}/4 criteria
                </span>
              </div>
              <div className="grid grid-cols-2 gap-x-4 gap-y-1">
                {lipinski.map((rule) => (
                  <div key={rule.label} className="flex items-center gap-1.5">
                    <LipinskiIndicator value={molecule} optimal={rule.check} />
                    <span className="text-xs text-gray-500">{rule.label}</span>
                  </div>
                ))}
              </div>
            </div>

            {/* Scoring Details */}
            {(molecule.vina_score || molecule.cnn_score > 0 || molecule.consensus_rank) && (
              <div className="mt-4 border-t border-gray-100 pt-4">
                <p className="text-xs font-medium text-gray-400 uppercase tracking-wide mb-3 flex items-center gap-1">
                  Scoring Details
                  <InfoTip text="Detailed scoring breakdown from the docking engine including Vina binding energy, CNN confidence scores, and consensus ranking." />
                </p>
                <div className="space-y-0.5">
                  {molecule.docking_engine && (
                    <div className="flex items-center justify-between py-1.5 px-3 rounded-lg hover:bg-gray-50">
                      <span className="text-xs font-medium text-gray-500">Engine</span>
                      <span className="text-xs font-semibold px-2 py-0.5 rounded-full bg-dockit-blue/10 text-dockit-blue">
                        {molecule.docking_engine}
                      </span>
                    </div>
                  )}
                  {molecule.vina_score != null && (
                    <PropRow label="Vina Score" value={molecule.vina_score} unit=" kcal/mol"
                      tooltip="AutoDock Vina binding energy. More negative = stronger predicted binding." />
                  )}
                  {molecule.cnn_score > 0 && (
                    <PropRow label="CNN Score" value={molecule.cnn_score}
                      tooltip="GNINA CNN confidence score (0-1). Higher = more confident the pose is correct." />
                  )}
                  {molecule.cnn_affinity > 0 && (
                    <PropRow label="CNN Affinity" value={molecule.cnn_affinity} unit=" pK"
                      tooltip="GNINA CNN predicted binding affinity in pK units. Higher = stronger predicted binding." />
                  )}
                  {molecule.consensus_rank != null && (
                    <PropRow label="Consensus Rank" value={molecule.consensus_rank}
                      tooltip="Overall ranking across multiple scoring functions. Lower = better." />
                  )}
                </div>
              </div>
            )}

            {/* V12: Pose Quality Panel */}
            {molecule.pose_quality && (
              <div className="mt-4 p-3 bg-blue-50/50 rounded-lg">
                <h4 className="text-xs font-semibold text-[#1e3a5f] mb-2">Pose Quality</h4>
                <div className="grid grid-cols-2 gap-2 text-xs">
                  <div className="flex justify-between">
                    <span className="text-gray-500">Contacts &lt;4A</span>
                    <span className={`font-medium ${
                      molecule.pose_quality.n_contacts_4A >= 5 ? 'text-green-600' :
                      molecule.pose_quality.n_contacts_4A >= 3 ? 'text-yellow-600' : 'text-red-500'
                    }`}>{molecule.pose_quality.n_contacts_4A}</span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-gray-500">H-bonds</span>
                    <span className="font-medium text-[#1e3a5f]">{molecule.pose_quality.n_hbonds}</span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-gray-500">Clashes</span>
                    <span className={`font-medium ${molecule.pose_quality.has_clashes ? 'text-red-500' : 'text-green-600'}`}>
                      {molecule.pose_quality.n_clashes}{!molecule.pose_quality.has_clashes && ' \u2713'}
                    </span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-gray-500">Quality</span>
                    <span className={`font-medium ${
                      molecule.pose_quality.interaction_quality >= 0.6 ? 'text-green-600' :
                      molecule.pose_quality.interaction_quality >= 0.3 ? 'text-yellow-600' : 'text-red-500'
                    }`}>{(molecule.pose_quality.interaction_quality * 100).toFixed(0)}%</span>
                  </div>
                </div>
                {molecule.pose_quality.key_residue_distances && Object.keys(molecule.pose_quality.key_residue_distances).length > 0 && (
                  <div className="mt-2 pt-2 border-t border-blue-100">
                    <span className="text-[10px] font-medium text-gray-500 uppercase tracking-wide">Key Residues</span>
                    <div className="mt-1 space-y-0.5">
                      {Object.entries(molecule.pose_quality.key_residue_distances).map(([res, dist]) => (
                        <div key={res} className="flex items-center justify-between text-xs">
                          <span className="font-mono text-gray-600">{res}</span>
                          <span className={`font-medium ${
                            dist < 3.5 ? 'text-green-600' : dist < 5.0 ? 'text-yellow-600' : 'text-red-500'
                          }`}>{dist}A {dist < 3.5 ? '\u25cf' : dist < 5.0 ? '\u25cf' : '\u25cf'}</span>
                        </div>
                      ))}
                    </div>
                  </div>
                )}
              </div>
            )}

            {/* Binding Interactions Summary */}
            {hasInteractions && molecule.interactions && (
              <div className="mt-4 border-t border-gray-100 pt-4">
                <p className="text-xs font-medium text-gray-400 uppercase tracking-wide mb-2 flex items-center gap-1">
                  Binding Interactions
                  <InfoTip text="Summary of protein-ligand interactions detected in the docking pose." />
                </p>
                <div className="flex items-center gap-3 flex-wrap">
                  {molecule.interactions.key_hbonds != null && (
                    <div className="flex items-center gap-1.5 px-2.5 py-1.5 bg-blue-50 rounded-lg">
                      <svg className="w-3.5 h-3.5 text-blue-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13.828 10.172a4 4 0 00-5.656 0l-4 4a4 4 0 105.656 5.656l1.102-1.101" />
                      </svg>
                      <span className="text-xs font-semibold text-blue-700">{molecule.interactions.key_hbonds}</span>
                      <span className="text-xs text-blue-500">H-bonds</span>
                    </div>
                  )}
                  {molecule.interactions.functional_contacts != null && (
                    <div className="flex items-center gap-1.5 px-2.5 py-1.5 bg-purple-50 rounded-lg">
                      <span className="text-xs font-semibold text-purple-700">{molecule.interactions.functional_contacts}</span>
                      <span className="text-xs text-purple-500">contacts</span>
                    </div>
                  )}
                  {molecule.interactions.interaction_quality != null && (
                    <div className="flex items-center gap-1.5 px-2.5 py-1.5 bg-green-50 rounded-lg">
                      <span className="text-xs font-semibold text-green-700">{Math.round(molecule.interactions.interaction_quality * 100)}%</span>
                      <span className="text-xs text-green-500">quality</span>
                    </div>
                  )}
                </div>
                {molecule.interactions.functional_residues && molecule.interactions.functional_residues.length > 0 && (
                  <p className="mt-2 text-xs text-gray-500">
                    <span className="font-medium text-gray-600">Key residues: </span>
                    {molecule.interactions.functional_residues.slice(0, 5).join(', ')}
                    {molecule.interactions.functional_residues.length > 5 && (
                      <span className="text-gray-400"> +{molecule.interactions.functional_residues.length - 5} more</span>
                    )}
                  </p>
                )}
                <button
                  onClick={() => setActiveTab('interactions')}
                  className="mt-2 text-xs font-medium text-dockit-blue hover:underline"
                >
                  View all interactions →
                </button>
              </div>
            )}
          </div>
        </>
      )}

      {/* Tab: ADMET */}
      {hasV2 && activeTab === 'admet' && (
        <div className="p-4">
          <ADMETRadar admet={molecule.admet} />
        </div>
      )}

      {/* Tab: Synthesis */}
      {hasV2 && activeTab === 'synthesis' && (
        <div className="p-4">
          <SynthesisTree synthesisRoute={molecule.synthesis_route} />
        </div>
      )}

      {/* Tab: Interactions (V5bis + V6.2 diagram toggle) */}
      {hasV2 && activeTab === 'interactions' && (
        <div className="p-4 space-y-3">
          {/* View toggle */}
          <div className="flex items-center gap-1 p-0.5 bg-gray-100 rounded-lg w-fit">
            <button
              onClick={() => setInteractionView('table')}
              className={`px-3 py-1 rounded-md text-xs font-semibold transition-colors ${
                interactionView === 'table'
                  ? 'bg-white text-dockit-blue shadow-sm'
                  : 'text-gray-500 hover:text-gray-700'
              }`}
            >
              Table
            </button>
            <button
              onClick={() => setInteractionView('diagram')}
              className={`flex items-center gap-1 px-3 py-1 rounded-md text-xs font-semibold transition-colors ${
                interactionView === 'diagram'
                  ? 'bg-white text-dockit-blue shadow-sm'
                  : 'text-gray-500 hover:text-gray-700'
              }`}
            >
              <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                  d="M13.828 10.172a4 4 0 00-5.656 0l-4 4a4 4 0 105.656 5.656l1.102-1.101m-.758-4.899a4 4 0 005.656 0l4-4a4 4 0 00-5.656-5.656l-1.1 1.1" />
              </svg>
              Diagram
            </button>
          </div>
          {/* Content */}
          {interactionView === 'table' ? (
            <InteractionView interactions={molecule.interactions} />
          ) : (
            <InteractionDiagram interactions={molecule.interactions} />
          )}
        </div>
      )}

      {/* V6.1: Quick action CTAs */}
      {hasV2 && (
        <div className="flex items-center gap-2 px-4 py-2.5 border-t border-gray-100 bg-gray-50/50">
          {hasInteractions && activeTab !== 'interactions' && (
            <button
              onClick={() => setActiveTab('interactions')}
              className="flex items-center gap-1 px-2.5 py-1.5 text-xs font-medium text-dockit-blue hover:bg-blue-50 rounded-md transition-colors"
            >
              <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13.828 10.172a4 4 0 00-5.656 0l-4 4a4 4 0 105.656 5.656l1.102-1.101" />
              </svg>
              View Interactions
            </button>
          )}
          {hasSynthesis && activeTab !== 'synthesis' && (
            <button
              onClick={() => setActiveTab('synthesis')}
              className="flex items-center gap-1 px-2.5 py-1.5 text-xs font-medium text-dockit-blue hover:bg-blue-50 rounded-md transition-colors"
            >
              <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
              </svg>
              View Synthesis Route
            </button>
          )}
        </div>
      )}
    </div>
  )
}
