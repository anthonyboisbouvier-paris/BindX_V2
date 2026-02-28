import React from 'react'
import InfoTip from './InfoTip.jsx'

// --------------------------------------------------
// Traffic-light status dot
// --------------------------------------------------
function StatusDot({ status }) {
  if (status === 'safe') {
    return (
      <span className="inline-flex items-center gap-1.5 text-green-600 font-semibold text-sm">
        <span className="w-2.5 h-2.5 rounded-full bg-green-500 flex-shrink-0" />
        Safe
      </span>
    )
  }
  if (status === 'risk') {
    return (
      <span className="inline-flex items-center gap-1.5 text-red-600 font-semibold text-sm">
        <span className="w-2.5 h-2.5 rounded-full bg-red-500 flex-shrink-0" />
        Risk
      </span>
    )
  }
  return (
    <span className="inline-flex items-center gap-1.5 text-yellow-600 font-semibold text-sm">
      <span className="w-2.5 h-2.5 rounded-full bg-yellow-400 flex-shrink-0" />
      Unknown
    </span>
  )
}

// --------------------------------------------------
// Selectivity badge
// --------------------------------------------------
function SelectivityBadge({ nSafe, nTotal }) {
  const ratio = nTotal > 0 ? nSafe / nTotal : 0
  let colorClass = 'bg-red-100 text-red-700 border-red-200'
  if (ratio >= 0.9) colorClass = 'bg-green-100 text-green-700 border-green-200'
  else if (ratio >= 0.7) colorClass = 'bg-yellow-100 text-yellow-700 border-yellow-200'

  return (
    <span className={`inline-flex items-center gap-1.5 px-3 py-1 rounded-full text-sm font-semibold border ${colorClass}`}>
      <span className={`w-2 h-2 rounded-full ${ratio >= 0.9 ? 'bg-green-500' : ratio >= 0.7 ? 'bg-yellow-400' : 'bg-red-500'}`} />
      {nSafe}/{nTotal} anti-targets clear
    </span>
  )
}

// --------------------------------------------------
// SafetyReport component
// --------------------------------------------------
export default function SafetyReport({ offTargetResults, moleculeName }) {
  if (!offTargetResults) {
    return (
      <div className="card p-6 text-center text-gray-400 text-sm">
        No off-target safety data available.
      </div>
    )
  }

  const { results = {}, selectivity_score, n_safe = 0, n_total = 0, warnings = [] } = offTargetResults

  const entries = Object.entries(results)
  const riskEntries = entries.filter(([, v]) => v.status === 'risk')

  return (
    <div className="card overflow-hidden">
      {/* Header */}
      <div className="bg-bx-surface px-5 py-4 text-white">
        <div className="flex flex-col sm:flex-row sm:items-center sm:justify-between gap-2">
          <div>
            <h3 className="font-bold text-base">
              Off-Target Safety Report
              {moleculeName && (
                <span className="text-bx-mint ml-2 font-mono">{moleculeName}</span>
              )}
            </h3>
            <p className="text-white/60 text-sm mt-0.5">
              Cross-reactivity analysis against known anti-targets
              <InfoTip text="Anti-targets are proteins that, if inadvertently bound by the drug candidate, could cause adverse effects such as cardiac toxicity (hERG) or drug-drug interactions (CYPs)." />
            </p>
          </div>
          <SelectivityBadge nSafe={n_safe} nTotal={n_total} />
        </div>
        {selectivity_score !== undefined && (
          <div className="mt-2 text-white/70 text-sm font-mono">
            Selectivity score: {Number(selectivity_score).toFixed(2)}
            <InfoTip text="Composite selectivity score (0-1). A score close to 1 means the molecule is highly selective for the intended target and shows minimal binding to anti-targets." />
          </div>
        )}
      </div>

      {/* Table */}
      <div className="overflow-x-auto">
        <table className="w-full text-sm">
          <thead>
            <tr className="bg-gray-50 border-b border-gray-100">
              <th className="text-left px-4 py-3 text-sm font-semibold text-gray-500 uppercase tracking-wide">
                Anti-target
                <InfoTip text="Protein that may be unintentionally bound, causing off-target side effects." />
              </th>
              <th className="text-right px-4 py-3 text-sm font-semibold text-gray-500 uppercase tracking-wide">
                Score (kcal/mol)
                <InfoTip text="Predicted binding affinity to this anti-target. More negative = stronger binding (higher risk)." />
              </th>
              <th className="text-right px-4 py-3 text-sm font-semibold text-gray-500 uppercase tracking-wide">
                Threshold
                <InfoTip text="Risk threshold for this anti-target. Scores below this threshold indicate potential off-target activity." />
              </th>
              <th className="text-center px-4 py-3 text-sm font-semibold text-gray-500 uppercase tracking-wide">
                Status
              </th>
              <th className="text-left px-4 py-3 text-sm font-semibold text-gray-500 uppercase tracking-wide hidden sm:table-cell">
                Risk
              </th>
            </tr>
          </thead>
          <tbody className="divide-y divide-gray-50">
            {entries.length === 0 && (
              <tr>
                <td colSpan={5} className="px-4 py-8 text-center text-gray-400 text-sm italic">
                  No anti-target data available.
                </td>
              </tr>
            )}
            {entries.map(([target, data]) => (
              <tr
                key={target}
                className={`hover:bg-gray-50 transition-colors ${data.status === 'risk' ? 'bg-red-50/40' : ''}`}
              >
                <td className="px-4 py-3 font-medium text-gray-800 text-sm">{target}</td>
                <td className="px-4 py-3 text-right font-mono text-sm text-gray-600">
                  {data.score !== undefined ? Number(data.score).toFixed(1) : 'N/A'}
                </td>
                <td className="px-4 py-3 text-right font-mono text-sm text-gray-400">
                  {data.threshold !== undefined ? Number(data.threshold).toFixed(1) : 'N/A'}
                </td>
                <td className="px-4 py-3 text-center">
                  <StatusDot status={data.status} />
                </td>
                <td className="px-4 py-3 text-sm text-gray-500 hidden sm:table-cell">
                  {data.risk_description || '—'}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      {/* Warning box */}
      {riskEntries.length > 0 && (
        <div className="mx-5 my-4 p-4 bg-red-50 border border-red-200 rounded-lg">
          <div className="flex items-start gap-2">
            <svg className="w-4 h-4 text-red-500 flex-shrink-0 mt-0.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
            </svg>
            <div>
              <p className="text-sm font-semibold text-red-700 mb-1">
                Off-target activity detected ({riskEntries.length} target{riskEntries.length > 1 ? 's' : ''})
              </p>
              {warnings.map((w, i) => (
                <p key={i} className="text-sm text-red-600 leading-relaxed mb-1">{w}</p>
              ))}
              {warnings.length === 0 && riskEntries.map(([target, data]) => (
                <p key={target} className="text-sm text-red-600 leading-relaxed">
                  {target}: predicted binding ({Number(data.score).toFixed(1)} kcal/mol) exceeds risk threshold ({Number(data.threshold).toFixed(1)} kcal/mol). Risk: {data.risk_description}.
                </p>
              ))}
              <p className="text-sm text-red-500 mt-2 font-medium">
                Recommendation: Consider structural modifications to reduce binding to flagged anti-targets before further development.
              </p>
            </div>
          </div>
        </div>
      )}

      {/* SEA Broad Screening section */}
      {(offTargetResults.sea_results || offTargetResults.combined_selectivity !== undefined || offTargetResults.tier1_hits) && (
        <div className="mx-5 my-4 border border-blue-200 rounded-lg overflow-hidden">
          {/* SEA section header */}
          <div className="bg-blue-50 px-4 py-2.5 flex items-center justify-between gap-3 border-b border-blue-200">
            <div className="flex items-center gap-2">
              <svg className="w-4 h-4 text-blue-600 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
              </svg>
              <span className="text-sm font-semibold text-blue-800">
                Broad Screening (SEA): 3000+ targets
              </span>
              <InfoTip text="Similarity Ensemble Approach (SEA) screening against 3000+ protein targets. Identifies potential off-target activity based on chemical similarity to known ligands." />
            </div>
            {offTargetResults.combined_selectivity !== undefined && (() => {
              const cs = Number(offTargetResults.combined_selectivity)
              let badgeClass = 'bg-red-100 text-red-700 border-red-300'
              if (cs >= 0.8) badgeClass = 'bg-green-100 text-green-700 border-green-300'
              else if (cs >= 0.6) badgeClass = 'bg-yellow-100 text-yellow-700 border-yellow-300'
              return (
                <span className={`inline-flex items-center gap-1.5 px-2.5 py-0.5 rounded-full text-sm font-bold border ${badgeClass} flex-shrink-0`}>
                  <span className={`w-1.5 h-1.5 rounded-full ${cs >= 0.8 ? 'bg-green-500' : cs >= 0.6 ? 'bg-yellow-400' : 'bg-red-500'}`} />
                  Combined selectivity: {cs.toFixed(2)}
                </span>
              )
            })()}
          </div>

          <div className="px-4 py-3 bg-white space-y-3">
            {/* Screened count */}
            {offTargetResults.sea_results?.n_targets_screened !== undefined && (
              <p className="text-sm text-gray-500">
                <span className="font-semibold text-gray-700">{offTargetResults.sea_results.n_targets_screened.toLocaleString()}</span> targets screened
                {offTargetResults.sea_results.targets_hit?.length > 0 && (
                  <span> &mdash; <span className="font-semibold text-red-600">{offTargetResults.sea_results.targets_hit.length}</span> hits above threshold</span>
                )}
              </p>
            )}

            {/* Dangerous targets hit with probability bars */}
            {offTargetResults.sea_results?.targets_hit?.length > 0 && (
              <div className="space-y-1.5">
                <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wide">Predicted off-target hits</p>
                {offTargetResults.sea_results.targets_hit.map((hit, idx) => {
                  const prob = hit.probability ?? hit.p_value ?? hit.score ?? 0
                  const pct = Math.round(Math.min(1, Math.max(0, prob)) * 100)
                  let barColor = 'bg-red-400'
                  let textColor = 'text-red-600'
                  if (prob < 0.5) { barColor = 'bg-yellow-400'; textColor = 'text-yellow-700' }
                  if (prob < 0.3) { barColor = 'bg-orange-400'; textColor = 'text-orange-600' }
                  const targetName = hit.target_name || hit.target || hit.chembl_id || `Target ${idx + 1}`
                  return (
                    <div key={idx} className="flex items-center gap-2">
                      <span className="text-sm text-gray-700 w-36 truncate flex-shrink-0 font-medium" title={targetName}>
                        {targetName}
                      </span>
                      <div className="flex-1 h-2 bg-gray-100 rounded-full overflow-hidden">
                        <div
                          className={`h-full rounded-full ${barColor} transition-all`}
                          style={{ width: `${pct}%` }}
                        />
                      </div>
                      <span className={`text-sm font-semibold ${textColor} w-9 text-right flex-shrink-0`}>
                        {pct}%
                      </span>
                    </div>
                  )
                })}
              </div>
            )}

            {/* Tier 1 hits list */}
            {offTargetResults.tier1_hits?.length > 0 && (
              <div>
                <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wide mb-1.5">Tier 1 priority hits</p>
                <div className="flex flex-wrap gap-1.5">
                  {offTargetResults.tier1_hits.map((hit, idx) => {
                    const name = typeof hit === 'string' ? hit : (hit.target_name || hit.target || String(hit))
                    return (
                      <span
                        key={idx}
                        className="inline-flex items-center gap-1 px-2 py-0.5 rounded text-[10px] font-semibold bg-red-50 border border-red-200 text-red-700"
                      >
                        <span className="w-1.5 h-1.5 rounded-full bg-red-500 flex-shrink-0" />
                        {name}
                      </span>
                    )
                  })}
                </div>
              </div>
            )}

            {/* No hits message */}
            {(!offTargetResults.sea_results?.targets_hit || offTargetResults.sea_results.targets_hit.length === 0) &&
              !offTargetResults.tier1_hits?.length && (
              <p className="text-sm text-green-600 font-medium flex items-center gap-1.5">
                <svg className="w-3.5 h-3.5" fill="currentColor" viewBox="0 0 20 20">
                  <path fillRule="evenodd" d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z" clipRule="evenodd" />
                </svg>
                No significant SEA hits detected above threshold.
              </p>
            )}
          </div>
        </div>
      )}

      {/* Footer disclaimer */}
      <div className="px-5 py-3 bg-gray-50 border-t border-gray-100">
        <p className="text-sm text-gray-400 italic leading-relaxed">
          Computational predictions. Off-target binding should be confirmed experimentally (e.g., selectivity panels, hERG patch-clamp assays).
        </p>
      </div>
    </div>
  )
}
