import React, { useState, useMemo } from 'react'
import InfoTip, { TIPS } from './InfoTip.jsx'

// ---------------------------------------------------------------------------
// Enamine REAL Space — approximate property distributions (UI estimation only)
// Source: Gorgulla et al., VirtualFlow 2.0, bioRxiv 2023
// ---------------------------------------------------------------------------
const TOTAL_MOLECULES = 69e9
const AVG_TRANCHE_SIZE = 5600 // paper: "average of 5,600 molecules per occupied tranche"

const MW_CDF  = [[100, 0.01], [200, 0.04], [300, 0.14], [350, 0.30], [400, 0.50], [450, 0.70], [500, 0.82], [600, 0.93], [800, 0.99], [1000, 1.0]]
const LOGP_CDF = [[0, 0.04], [1, 0.10], [2, 0.25], [3, 0.48], [4, 0.68], [5, 0.83], [6, 0.92], [8, 0.98], [10, 1.0]]
const HBD_CDF = [[0, 0.08], [1, 0.25], [2, 0.48], [3, 0.70], [4, 0.85], [5, 0.93], [8, 0.99], [15, 1.0]]
const HBA_CDF = [[0, 0.01], [2, 0.06], [4, 0.22], [6, 0.50], [8, 0.78], [10, 0.92], [12, 0.97], [15, 0.99], [20, 1.0]]

// Lipinski baseline for cost scaling (MW≤500, logP≤5, HBD≤5, HBA≤10)
const LIPINSKI_FRAC = 0.58

function interpolateCDF(value, cdf) {
  if (value <= cdf[0][0]) return cdf[0][1]
  if (value >= cdf[cdf.length - 1][0]) return cdf[cdf.length - 1][1]
  for (let i = 0; i < cdf.length - 1; i++) {
    const [t1, f1] = cdf[i]
    const [t2, f2] = cdf[i + 1]
    if (value >= t1 && value <= t2) {
      return f1 + (f2 - f1) * (value - t1) / (t2 - t1)
    }
  }
  return 1.0
}

function estimateFilteredFraction(mw_max, logp_max, hbd_max, hba_max) {
  return interpolateCDF(mw_max, MW_CDF)
    * interpolateCDF(logp_max, LOGP_CDF)
    * interpolateCDF(hbd_max, HBD_CDF)
    * interpolateCDF(hba_max, HBA_CDF)
}

function formatBillions(n) {
  if (n >= 1e9) return `${(n / 1e9).toFixed(1)}B`
  if (n >= 1e6) return `${(n / 1e6).toFixed(0)}M`
  if (n >= 1e3) return `${(n / 1e3).toFixed(0)}K`
  return String(Math.round(n))
}

function formatCount(n) {
  if (n >= 1e9) return `${(n / 1e9).toFixed(1)}B`
  if (n >= 1e6) return `${(n / 1e6).toFixed(1)}M`
  if (n >= 1e3) return `${(n / 1e3).toFixed(0)}K`
  return String(Math.round(n))
}

// ---------------------------------------------------------------------------
// Strategy presets
// Paper-validated: reps=1, exhaust_pre=1, exhaust_pri=8 for all.
// Differentiation = coverage only. Costs calibrated for Lipinski filters.
// ---------------------------------------------------------------------------
const STRATEGY_PRESETS = {
  lean: {
    label: 'Lean',
    subtitle: 'Quick pocket validation',
    duration: '~3–6 h',
    top_n_results: 5000,
    cost_low: 100, cost_high: 200,
    budget_cap_usd: 250,
    exhaustiveness_prescreen: 1, exhaustiveness_primary: 8,
    reps_per_tranche: 1, tranche_pct: 0.05,
  },
  balanced: {
    label: 'Balanced',
    subtitle: 'Best accuracy/cost ratio',
    duration: '~12–24 h',
    top_n_results: 10000,
    cost_low: 300, cost_high: 600,
    budget_cap_usd: 700,
    exhaustiveness_prescreen: 1, exhaustiveness_primary: 8,
    reps_per_tranche: 1, tranche_pct: 0.2,
  },
  aggressive: {
    label: 'Aggressive',
    subtitle: 'Maximum chemical diversity',
    duration: '~24–48 h',
    top_n_results: 50000,
    cost_low: 1500, cost_high: 3000,
    budget_cap_usd: 3500,
    exhaustiveness_prescreen: 1, exhaustiveness_primary: 8,
    reps_per_tranche: 1, tranche_pct: 1.0,
  },
}

// ---------------------------------------------------------------------------
// Molecular filter presets
// ---------------------------------------------------------------------------
const FILTER_PRESETS = {
  no_filter: { label: 'No filter', mw_max: 1000, logp_max: 10, hbd_max: 15, hba_max: 20 },
  lipinski: { label: 'Lipinski (Ro5)', mw_max: 500, logp_max: 5, hbd_max: 5, hba_max: 10 },
  beyond_ro5: { label: 'Beyond Ro5', mw_max: 800, logp_max: 6, hbd_max: 6, hba_max: 15 },
  fragment: { label: 'Fragment (Ro3)', mw_max: 300, logp_max: 3, hbd_max: 3, hba_max: 3 },
}

// ---------------------------------------------------------------------------
// AFVSConfigForm
// ---------------------------------------------------------------------------
export default function AFVSConfigForm({ config, onChange, project, campaign }) {
  const [showAdvanced, setShowAdvanced] = useState(false)
  const set = (key, val) => onChange({ ...config, [key]: val })

  const applyPreset = (presetKey) => {
    const preset = STRATEGY_PRESETS[presetKey]
    onChange({
      ...config,
      strategy: presetKey,
      exhaustiveness_prescreen: preset.exhaustiveness_prescreen,
      exhaustiveness_primary: preset.exhaustiveness_primary,
      reps_per_tranche: preset.reps_per_tranche,
      tranche_pct: preset.tranche_pct,
      top_n_results: preset.top_n_results,
      budget_cap_usd: preset.budget_cap_usd,
    })
  }

  const applyFilterPreset = (key) => {
    const fp = FILTER_PRESETS[key]
    onChange({ ...config, mw_max: fp.mw_max, logp_max: fp.logp_max, hbd_max: fp.hbd_max, hba_max: fp.hba_max })
  }

  const pocketName = campaign?.pocket_config?.name || 'Default pocket'
  const targetName = project?.target_name || project?.target_pdb_id || 'Target'
  const activePreset = STRATEGY_PRESETS[config.strategy]

  // --- Filter impact ---
  const filterEstimate = useMemo(() => {
    const frac = estimateFilteredFraction(
      config.mw_max ?? 500, config.logp_max ?? 5,
      config.hbd_max ?? 5, config.hba_max ?? 10,
    )
    const remaining = TOTAL_MOLECULES * frac
    const tranches = Math.round(remaining / AVG_TRANCHE_SIZE)
    return { frac, remaining, pct: frac * 100, tranches }
  }, [config.mw_max, config.logp_max, config.hbd_max, config.hba_max])

  // --- Screening estimates ---
  const screeningEstimate = useMemo(() => {
    const reps = config.reps_per_tranche ?? 1
    const tranchePct = config.tranche_pct ?? 0.5

    const prescreenDockings = filterEstimate.tranches * reps
    const primaryTranches = Math.round(filterEstimate.tranches * tranchePct / 100)
    const primaryMolecules = primaryTranches * AVG_TRANCHE_SIZE

    const costMultiplier = LIPINSKI_FRAC > 0 ? filterEstimate.frac / LIPINSKI_FRAC : 1
    const costLow = activePreset ? Math.round(activePreset.cost_low * costMultiplier) : 0
    const costHigh = activePreset ? Math.round(activePreset.cost_high * costMultiplier) : 0
    const costLabel = `$${costLow.toLocaleString()}–${costHigh.toLocaleString()}`

    return { prescreenDockings, primaryTranches, primaryMolecules, costMultiplier, costLow, costHigh, costLabel }
  }, [filterEstimate, config.reps_per_tranche, config.tranche_pct, activePreset])

  // Budget warning
  const budgetBelowPreset = activePreset && config.budget_cap_usd < screeningEstimate.costLow * 0.5

  // Helper: compute primary molecules for a preset with current filters
  const presetPrimaryMol = (preset) => {
    const t = Math.round(filterEstimate.tranches * preset.tranche_pct / 100)
    return t * AVG_TRANCHE_SIZE
  }

  // Helper: compute cost label for a preset with current filters
  const presetCostLabel = (preset) => {
    const m = LIPINSKI_FRAC > 0 ? filterEstimate.frac / LIPINSKI_FRAC : 1
    return `$${Math.round(preset.cost_low * m).toLocaleString()}–${Math.round(preset.cost_high * m).toLocaleString()}`
  }

  return (
    <div className="space-y-5">
      {/* --- Header + pipeline visual --- */}
      <div className="p-4 rounded-xl border border-cyan-200 bg-cyan-50 space-y-3">
        <div>
          <p className="text-sm font-semibold text-cyan-700">Ultra-Large Virtual Screening</p>
          <p className="text-xs text-cyan-600 mt-1">
            Screen billions of commercially available molecules against your target.
            The best-scoring hits are automatically imported into your Phase A dashboard.
          </p>
        </div>

        {/* Pipeline visual */}
        <div className="flex items-center gap-0 text-[10px]">
          {[
            { label: 'Filter', desc: 'Drug-like rules', sub: `69B → ${formatBillions(filterEstimate.remaining)}`, icon: filterIcon },
            { label: 'Prescreen', desc: '1 rep/tranche', sub: `${formatCount(filterEstimate.tranches)} tranches`, icon: prescreenIcon },
            { label: 'Primary', desc: `Top ${config.tranche_pct ?? 0.5}% re-docked`, sub: `${formatBillions(screeningEstimate.primaryMolecules)} molecules`, icon: dockingIcon },
          ].map((s, i) => (
            <React.Fragment key={s.label}>
              {i > 0 && (
                <div className="flex-shrink-0 px-0.5">
                  <svg className="w-4 h-4 text-cyan-300" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
                  </svg>
                </div>
              )}
              <div className="flex items-center gap-1.5 px-2 py-1.5 rounded-lg bg-cyan-100/60 flex-shrink-0">
                <span className="text-cyan-500">{s.icon}</span>
                <div>
                  <p className="font-semibold text-cyan-700 leading-none">{s.label}</p>
                  <p className="text-cyan-500 leading-none mt-0.5">{s.desc}</p>
                  <p className="text-cyan-400 leading-none mt-0.5 font-medium">{s.sub}</p>
                </div>
              </div>
            </React.Fragment>
          ))}
        </div>

        {/* Target / Pocket */}
        <div className="flex gap-4 text-xs text-cyan-700">
          <span><strong>Target:</strong> {targetName}</span>
          <span><strong>Pocket:</strong> {pocketName}</span>
        </div>
      </div>

      {/* --- Strategy presets --- */}
      <div>
        <p className="text-sm font-semibold text-gray-400 uppercase tracking-wide mb-2">Strategy</p>
        <div className="grid grid-cols-3 gap-2">
          {Object.entries(STRATEGY_PRESETS).map(([key, preset]) => {
            const active = config.strategy === key
            const mol = presetPrimaryMol(preset)
            const cost = presetCostLabel(preset)

            return (
              <button
                key={key}
                onClick={() => applyPreset(key)}
                className={`p-3 rounded-xl border-2 text-left transition-all ${
                  active ? 'border-cyan-400 bg-cyan-50 shadow-sm' : 'border-gray-100 hover:border-cyan-200'
                }`}
              >
                {/* Title + badge */}
                <div className="flex items-center gap-1.5 mb-0.5">
                  <p className={`text-sm font-bold ${active ? 'text-cyan-700' : 'text-gray-700'}`}>
                    {preset.label}
                  </p>
                  {key === 'balanced' && (
                    <span className="text-[7px] font-bold uppercase tracking-wider text-white bg-cyan-500 px-1.5 py-0.5 rounded-full leading-none">
                      Recommended
                    </span>
                  )}
                </div>
                <p className="text-[10px] text-gray-400 mb-2">{preset.subtitle}</p>

                {/* 3-metric row — centered, aligned */}
                <div className="grid grid-cols-3 divide-x divide-gray-100 bg-white rounded-lg py-1.5">
                  <MetricBlock label="Coverage" value={`${preset.tranche_pct}% (${formatBillions(mol)})`} active={active} />
                  <MetricBlock label="Duration" value={preset.duration} active={active} />
                  <MetricBlock label="Est. cost" value={cost} active={active} highlight />
                </div>

                {/* Sub info */}
                <p className="text-[9px] text-gray-400 mt-1.5 leading-tight">
                  {formatBillions(mol)} molecules re-docked · Top {preset.top_n_results.toLocaleString()} imported
                </p>
              </button>
            )
          })}
        </div>
      </div>

      {/* --- Molecular Filters --- */}
      <div className="space-y-3">
        <div className="flex items-center justify-between">
          <p className="text-sm font-semibold text-gray-400 uppercase tracking-wide">
            Molecular Filters
          </p>
          <InfoTip text={TIPS.afvs_filters} size="xs" />
        </div>

        {/* Filter presets */}
        <div className="flex gap-1.5">
          {Object.entries(FILTER_PRESETS).map(([key, fp]) => {
            const isActive = config.mw_max === fp.mw_max && config.logp_max === fp.logp_max
              && config.hbd_max === fp.hbd_max && config.hba_max === fp.hba_max
            return (
              <button key={key} onClick={() => applyFilterPreset(key)}
                className={`text-[10px] font-medium px-2.5 py-1 rounded-full border transition-all ${
                  isActive
                    ? 'border-cyan-400 bg-cyan-50 text-cyan-700'
                    : 'border-gray-200 text-gray-500 hover:border-cyan-200 hover:text-cyan-600'
                }`}
              >
                {fp.label}
              </button>
            )
          })}
        </div>

        <div className="grid grid-cols-2 gap-3">
          <NumberInput label="MW max (Da)" value={config.mw_max ?? 500} onChange={v => set('mw_max', v)} min={100} max={1000} tip={TIPS.MW} />
          <NumberInput label="logP max" value={config.logp_max ?? 5} onChange={v => set('logp_max', v)} min={0} max={10} step={0.5} tip={TIPS.logP} />
          <NumberInput label="HBD max" value={config.hbd_max ?? 5} onChange={v => set('hbd_max', v)} min={0} max={15} tip={TIPS.HBD} />
          <NumberInput label="HBA max" value={config.hba_max ?? 10} onChange={v => set('hba_max', v)} min={0} max={20} tip={TIPS.HBA} />
        </div>

        {/* Filter impact bar */}
        <div className="bg-white border border-gray-200 rounded-xl p-3 space-y-1.5">
          <div className="flex items-center justify-between text-[10px]">
            <span className="text-gray-500 font-medium">Estimated space after filters</span>
            <span className="font-semibold text-gray-700 tabular-nums">
              {formatBillions(filterEstimate.remaining)} / 69B ({filterEstimate.pct.toFixed(0)}%)
            </span>
          </div>
          <div className="w-full bg-gray-100 rounded-full h-2 overflow-hidden">
            <div
              className="h-full bg-cyan-400 rounded-full transition-all duration-300"
              style={{ width: `${Math.max(2, filterEstimate.pct)}%` }}
            />
          </div>
          <div className="flex items-center justify-between text-[10px] text-gray-400">
            <span>~{formatCount(filterEstimate.tranches)} tranches · {formatCount(screeningEstimate.prescreenDockings)} prescreen dockings</span>
            {filterEstimate.frac < 0.95 && (
              <span className="text-green-600 font-medium">
                ~{Math.round((1 - filterEstimate.frac) * 100)}% cost reduction vs. no filter
              </span>
            )}
          </div>
          <p className="text-[9px] text-gray-400">
            Approximate — based on published Enamine REAL Space distributions (~5.6K mol/tranche).
          </p>
        </div>
      </div>

      {/* --- Results & Budget --- */}
      <div className="space-y-3">
        <p className="text-sm font-semibold text-gray-400 uppercase tracking-wide">Results & Budget</p>
        <div className="grid grid-cols-2 gap-3">
          <div>
            <label className="text-[10px] text-gray-500 font-medium mb-1 flex items-center gap-1">
              Top molecules to import
              <InfoTip text={TIPS.afvs_top_n} size="xs" />
            </label>
            <input type="number" min={100} max={100000} step={1000}
              value={config.top_n_results ?? 10000}
              onChange={e => set('top_n_results', Math.min(100000, parseInt(e.target.value) || 10000))}
              className="w-full text-sm border border-gray-200 rounded-lg px-3 py-2 bg-white tabular-nums
                         focus:outline-none focus:ring-2 focus:ring-cyan-300/30 focus:border-cyan-400" />
            <p className="text-[9px] text-gray-400 mt-0.5">
              Best {(config.top_n_results ?? 10000).toLocaleString()} from {formatBillions(screeningEstimate.primaryMolecules)} re-docked
            </p>
          </div>
          <div>
            <label className="text-[10px] text-gray-500 font-medium mb-1 flex items-center gap-1">
              Budget cap (USD)
              <InfoTip text={TIPS.afvs_budget} size="xs" />
            </label>
            <div className="relative">
              <span className="absolute left-3 top-1/2 -translate-y-1/2 text-sm text-gray-400 font-medium">$</span>
              <input type="number" min={10} max={50000} step={100}
                value={config.budget_cap_usd ?? 1500}
                onChange={e => set('budget_cap_usd', parseFloat(e.target.value) || 1500)}
                className="w-full text-sm border border-gray-200 rounded-lg pl-7 pr-3 py-2 bg-white tabular-nums
                           focus:outline-none focus:ring-2 focus:ring-cyan-300/30 focus:border-cyan-400" />
            </div>
          </div>
        </div>

        {/* Budget warning */}
        {budgetBelowPreset && (
          <div className="flex items-start gap-1.5 text-[10px] text-amber-600 bg-amber-50 border border-amber-200 rounded-lg px-3 py-2">
            <svg className="w-3.5 h-3.5 flex-shrink-0 mt-0.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
            </svg>
            <span>
              Budget is below the estimated cost for <strong>{activePreset.label}</strong> ({screeningEstimate.costLabel}).
              The run may stop early with partial results.
            </span>
          </div>
        )}
      </div>

      {/* --- Cost estimation --- */}
      {activePreset && (
        <div className="bg-cyan-50/50 border border-cyan-200 rounded-xl px-4 py-3 space-y-1">
          <div className="flex items-center justify-between">
            <div className="text-sm text-cyan-700">
              <span className="font-medium">Estimated cost:</span>{' '}
              <span className="font-semibold">{screeningEstimate.costLabel}</span>
            </div>
            <span className="text-[8px] font-medium text-cyan-400 bg-cyan-100 px-2 py-0.5 rounded-full uppercase tracking-wider">
              Uncalibrated
            </span>
          </div>
          <p className="text-[10px] text-cyan-500">
            Based on ~$20–30K per 1B dockings (cloud CPU, Spot pricing).
            Will be refined after your first AFVS run.
            {screeningEstimate.costMultiplier > 1.1 && (
              <> Wider filters increase cost (~{Math.round(screeningEstimate.costMultiplier * 100 - 100)}% vs. Lipinski).</>
            )}
            {screeningEstimate.costMultiplier < 0.9 && (
              <> Stricter filters reduce cost (~{Math.round((1 - screeningEstimate.costMultiplier) * 100)}% cheaper than Lipinski).</>
            )}
          </p>
        </div>
      )}

      {/* --- Billing info --- */}
      <div className="bg-amber-50 border border-amber-200 rounded-xl px-4 py-3 text-sm text-amber-700 flex items-start gap-2">
        <svg className="w-5 h-5 flex-shrink-0 mt-0.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
        </svg>
        <span>
          You will be charged based on actual compute time.
          The run <strong>stops automatically</strong> when your budget cap is reached — any results found so far are still imported.
        </span>
      </div>

      {/* --- Advanced Settings (collapsed) --- */}
      <div className="border border-gray-200 rounded-xl overflow-hidden">
        <button
          onClick={() => setShowAdvanced(v => !v)}
          className="w-full flex items-center justify-between px-4 py-2.5 text-sm font-medium text-gray-500 hover:bg-gray-50 transition-colors"
        >
          <span>Advanced Settings</span>
          <svg className={`w-4 h-4 transition-transform ${showAdvanced ? 'rotate-180' : ''}`} fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
          </svg>
        </button>

        {showAdvanced && (
          <div className="px-4 pb-4 space-y-5 border-t border-gray-100">
            {/* Docking scenario */}
            <div className="pt-3">
              <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wide mb-2">Docking Engine</p>
              <div className="flex gap-2">
                {[
                  { id: 'smina_vinardo', label: 'Smina / Vinardo', desc: 'Standard — fast and reliable', recommended: true },
                  { id: 'qvina2', label: 'QVina2', desc: 'Thorough — better for flexible pockets' },
                ].map(sc => (
                  <button key={sc.id}
                    onClick={() => set('docking_scenario', sc.id)}
                    className={`flex-1 p-2.5 rounded-xl border-2 text-left transition-all ${
                      config.docking_scenario === sc.id
                        ? 'border-cyan-400 bg-cyan-50'
                        : 'border-gray-100 hover:border-gray-200'
                    }`}
                  >
                    <div className="flex items-center gap-1.5">
                      <p className={`text-sm font-semibold ${config.docking_scenario === sc.id ? 'text-cyan-700' : 'text-gray-700'}`}>{sc.label}</p>
                      {sc.recommended && (
                        <span className="text-[8px] font-bold text-gray-400 uppercase">default</span>
                      )}
                    </div>
                    <p className="text-[10px] text-gray-400 mt-0.5">{sc.desc}</p>
                  </button>
                ))}
              </div>
            </div>

            {/* Docking parameters */}
            <div className="space-y-3">
              <p className="text-[10px] font-semibold text-gray-400 uppercase tracking-wide">Docking Parameters</p>
              <p className="text-[9px] text-gray-400 -mt-1">
                Set automatically by the strategy. Only adjust if you know what you're doing.
              </p>

              <SliderRow label="Prescreen thoroughness" min={1} max={8}
                value={config.exhaustiveness_prescreen ?? 1} onChange={v => set('exhaustiveness_prescreen', v)}
                leftLabel="Fast" rightLabel="Thorough" tip={TIPS.afvs_exhaustiveness_prescreen} />

              <SliderRow label="Primary thoroughness" min={1} max={16}
                value={config.exhaustiveness_primary ?? 8} onChange={v => set('exhaustiveness_primary', v)}
                leftLabel="Fast" rightLabel="Thorough" tip={TIPS.afvs_exhaustiveness_primary} />

              <SliderRow label="Representatives / tranche" min={1} max={3}
                value={config.reps_per_tranche ?? 1} onChange={v => set('reps_per_tranche', v)}
                leftLabel="1" rightLabel="3" tip={TIPS.afvs_reps} />

              <SliderRow label={`Coverage — ${formatBillions(screeningEstimate.primaryMolecules)} molecules`} min={0.01} max={2} step={0.01} unit="%"
                value={config.tranche_pct ?? 0.5} onChange={v => set('tranche_pct', v)}
                leftLabel="Narrow" rightLabel="Wide" tip={TIPS.afvs_coverage} />
            </div>
          </div>
        )}
      </div>
    </div>
  )
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function MetricBlock({ label, value, active, highlight }) {
  return (
    <div className="text-center px-1">
      <p className={`text-[11px] font-bold tabular-nums leading-tight ${
        active ? (highlight ? 'text-cyan-600' : 'text-cyan-700') : 'text-gray-700'
      }`}>
        {value}
      </p>
      <p className="text-[8px] text-gray-400 uppercase tracking-wide mt-0.5">{label}</p>
    </div>
  )
}

function SliderRow({ label, min, max, step = 1, value, onChange, unit = '', leftLabel, rightLabel, tip }) {
  return (
    <div>
      <div className="flex justify-between text-[10px] mb-1">
        <span className="text-gray-500 font-medium flex items-center gap-1">
          {label}
          {tip && <InfoTip text={tip} size="xs" />}
        </span>
        <span className="font-semibold text-gray-700 tabular-nums">{value}{unit}</span>
      </div>
      <input
        type="range" min={min} max={max} step={step} value={value}
        onChange={e => onChange(parseFloat(e.target.value))}
        className="w-full h-1.5 bg-gray-200 rounded-full appearance-none cursor-pointer accent-cyan-500"
      />
      {(leftLabel || rightLabel) && (
        <div className="flex justify-between text-[9px] text-gray-300 mt-0.5">
          <span>{leftLabel}</span>
          <span>{rightLabel}</span>
        </div>
      )}
    </div>
  )
}

function NumberInput({ label, value, onChange, min, max, step = 1, tip }) {
  return (
    <div>
      <label className="text-[10px] text-gray-500 font-medium mb-1 flex items-center gap-1">
        {label}
        {tip && <InfoTip text={tip} size="xs" />}
      </label>
      <input type="number" min={min} max={max} step={step}
        value={value}
        onChange={e => onChange(parseFloat(e.target.value) || min)}
        className="w-full text-sm border border-gray-200 rounded-lg px-3 py-2 bg-white tabular-nums
                   focus:outline-none focus:ring-2 focus:ring-cyan-300/30 focus:border-cyan-400" />
    </div>
  )
}

// ---------------------------------------------------------------------------
// Inline SVG icons
// ---------------------------------------------------------------------------
const filterIcon = <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M3 4a1 1 0 011-1h16a1 1 0 011 1v2.586a1 1 0 01-.293.707l-6.414 6.414a1 1 0 00-.293.707V17l-4 4v-6.586a1 1 0 00-.293-.707L3.293 7.293A1 1 0 013 6.586V4z" /></svg>
const prescreenIcon = <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 10V3L4 14h7v7l9-11h-7z" /></svg>
const dockingIcon = <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12l2 2 4-4m5.618-4.016A11.955 11.955 0 0112 2.944a11.955 11.955 0 01-8.618 3.04A12.02 12.02 0 003 9c0 5.591 3.824 10.29 9 11.622 5.176-1.332 9-6.03 9-11.622 0-1.042-.133-2.052-.382-3.016z" /></svg>
