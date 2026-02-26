import React from 'react'

// ---------------------------------------------------------------------------
// Gradient progress bar (red -> yellow -> green)
// ---------------------------------------------------------------------------

function GradientBar({ value, label }) {
  const pct = Math.round((value || 0) * 100)
  // Color stops: 0% red, 50% yellow, 100% green
  const barColor = pct >= 80 ? '#22c55e' : pct >= 55 ? '#f59e0b' : '#ef4444'

  return (
    <div className="flex items-center gap-3">
      <span className="text-xs text-gray-500 w-24 shrink-0">{label}</span>
      <div className="flex-1 h-2 bg-gray-100 rounded-full overflow-hidden">
        <div
          className="h-full rounded-full transition-all duration-500"
          style={{
            width: `${pct}%`,
            background: `linear-gradient(90deg, #ef4444 0%, #f59e0b 50%, #22c55e 100%)`,
            clipPath: `inset(0 ${100 - pct}% 0 0)`,
          }}
        />
      </div>
      <span
        className="text-xs font-semibold w-9 text-right shrink-0"
        style={{ color: barColor }}
      >
        {pct}%
      </span>
    </div>
  )
}

// ---------------------------------------------------------------------------
// TargetInfoCard
// ---------------------------------------------------------------------------
// Props: { assessment } — from MOCK_TARGET_ASSESSMENT[projectId]

export default function TargetInfoCard({ assessment, project }) {
  if (!assessment) return null

  const {
    druggability_score,
    evidence_score,
    confidence,
    disease_context,
    known_drugs,
    protein_family,
    gene_name,
    organism,
    sequence_length,
    pdb_structures,
    chembl_compounds,
    binding_site,
  } = assessment

  return (
    <div className="bg-white rounded-xl border border-gray-100 shadow-sm overflow-hidden">
      {/* Header */}
      <div className="bg-gradient-to-r from-[#1e3a5f] to-[#1e4a7f] px-5 py-4">
        <div className="flex items-start justify-between gap-4">
          <div>
            <div className="flex items-center gap-2 mb-1">
              <span className="text-white font-bold text-base">{gene_name || project?.target_name}</span>
              <span className="text-blue-200 text-sm">—</span>
              <span className="text-blue-200 text-sm">{protein_family}</span>
            </div>
            <div className="flex items-center gap-2 flex-wrap">
              {project?.uniprot_id && (
                <span className="text-xs font-mono bg-white/10 text-blue-100 px-2 py-0.5 rounded">
                  {project.uniprot_id}
                </span>
              )}
              {project?.target_pdb_id && (
                <span className="text-xs font-mono bg-white/10 text-blue-100 px-2 py-0.5 rounded">
                  PDB: {project.target_pdb_id}
                </span>
              )}
              {organism && (
                <span className="text-xs text-blue-300 italic">{organism}</span>
              )}
            </div>
          </div>

          {/* Druggability score circle */}
          <div className="shrink-0 text-center">
            <div className="w-14 h-14 rounded-full bg-white/10 border-2 border-[#22c55e] flex items-center justify-center">
              <span className="text-[#22c55e] font-bold text-base">
                {Math.round((druggability_score || 0) * 100)}
              </span>
            </div>
            <p className="text-blue-300 text-xs mt-1">Drug.</p>
          </div>
        </div>
      </div>

      {/* Body */}
      <div className="px-5 py-4 space-y-4">
        {/* Score bars */}
        <div className="space-y-2">
          <GradientBar value={druggability_score} label="Druggability" />
          <GradientBar value={evidence_score} label="Evidence" />
          <GradientBar value={confidence} label="Confidence" />
        </div>

        <div className="border-t border-gray-50" />

        {/* Quick stats */}
        <div className="grid grid-cols-3 gap-3">
          <div className="bg-gray-50 rounded-lg px-3 py-2 text-center">
            <p className="text-base font-bold text-[#1e3a5f]">{sequence_length?.toLocaleString()}</p>
            <p className="text-xs text-gray-400 mt-0.5">Amino acids</p>
          </div>
          <div className="bg-gray-50 rounded-lg px-3 py-2 text-center">
            <p className="text-base font-bold text-[#1e3a5f]">{pdb_structures}</p>
            <p className="text-xs text-gray-400 mt-0.5">PDB structures</p>
          </div>
          <div className="bg-gray-50 rounded-lg px-3 py-2 text-center">
            <p className="text-base font-bold text-[#1e3a5f]">{chembl_compounds?.toLocaleString()}</p>
            <p className="text-xs text-gray-400 mt-0.5">ChEMBL cpds</p>
          </div>
        </div>

        <div className="border-t border-gray-50" />

        {/* Disease + known drugs */}
        <div className="space-y-1.5">
          {disease_context && (
            <div className="flex items-start gap-2">
              <span className="text-xs font-medium text-gray-400 w-20 shrink-0 pt-0.5">Disease</span>
              <span className="text-sm text-gray-700">{disease_context}</span>
            </div>
          )}
          {known_drugs && known_drugs.length > 0 && (
            <div className="flex items-start gap-2">
              <span className="text-xs font-medium text-gray-400 w-20 shrink-0 pt-0.5">Known drugs</span>
              <div className="flex flex-wrap gap-1">
                {known_drugs.map(d => (
                  <span key={d} className="text-xs bg-blue-50 text-blue-700 px-2 py-0.5 rounded-md font-medium">
                    {d}
                  </span>
                ))}
              </div>
            </div>
          )}
        </div>

        {/* Binding site */}
        {binding_site && (
          <>
            <div className="border-t border-gray-50" />
            <div>
              <p className="text-xs font-semibold text-gray-500 uppercase tracking-wider mb-2">Binding Site</p>
              <div className="flex items-start gap-2 mb-2">
                <span className="text-xs font-medium text-gray-400 w-20 shrink-0 pt-0.5">Residues</span>
                <div className="flex flex-wrap gap-1">
                  {(binding_site.residues || []).map(r => (
                    <span key={r} className="text-xs font-mono bg-amber-50 text-amber-700 px-1.5 py-0.5 rounded">
                      {r}
                    </span>
                  ))}
                </div>
              </div>
              <div className="flex items-center gap-4">
                {binding_site.volume && (
                  <span className="text-xs text-gray-500">
                    Volume: <span className="font-semibold text-gray-700">{binding_site.volume} A3</span>
                  </span>
                )}
                {binding_site.druggability_probability && (
                  <span className="text-xs text-gray-500">
                    Pocket drug.: <span className="font-semibold text-[#22c55e]">
                      {Math.round(binding_site.druggability_probability * 100)}%
                    </span>
                  </span>
                )}
              </div>
            </div>
          </>
        )}
      </div>
    </div>
  )
}
