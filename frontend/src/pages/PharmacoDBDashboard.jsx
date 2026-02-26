import React, { useState, useEffect, useCallback, useRef } from 'react'
import { supabase } from '../lib/supabase.js'
import {
  AreaChart, Area, BarChart, Bar, PieChart, Pie, Cell,
  ScatterChart, Scatter, XAxis, YAxis, ZAxis, CartesianGrid, Tooltip,
  Legend, ResponsiveContainer, ReferenceLine, LabelList,
} from 'recharts'

// ─── Color palette ───────────────────────────────────────────────────────────
const C = {
  blue:    '#1e3a5f',
  blueMid: '#2d5a8e',
  blueLight:'#3b82f6',
  green:   '#22c55e',
  greenDark:'#16a34a',
  teal:    '#14b8a6',
  emerald: '#10b981',
  amber:   '#f59e0b',
  violet:  '#8b5cf6',
  rose:    '#f43f5e',
  sky:     '#0ea5e9',
  indigo:  '#6366f1',
  orange:  '#f97316',
  cyan:    '#06b6d4',
}

const FAMILY_COLORS = [
  C.blueLight, C.teal, C.emerald, C.amber, C.violet,
  C.rose, C.sky, C.indigo, C.orange, C.cyan,
  '#a78bfa', '#34d399', '#fbbf24', '#60a5fa', '#f472b6',
]

const PHASE_COLORS = {
  '-1': '#ef4444', '0': '#94a3b8', '1': '#f59e0b',
  '2': '#3b82f6', '3': '#8b5cf6', '4': '#22c55e',
}

// ─── Custom hook: useRPC ──────────────────────────────────────────────────────
function useRPC(fnName, options = {}) {
  const { params = {}, enabled = true } = options
  const [data, setData]     = useState(null)
  const [loading, setLoading] = useState(true)
  const [error, setError]   = useState(null)
  const [lastFetch, setLastFetch] = useState(null)

  const fetch = useCallback(async () => {
    if (!enabled) return
    setLoading(true)
    try {
      const { data: result, error: err } = await supabase.rpc(fnName, params)
      if (err) throw err
      setData(result)
      setLastFetch(new Date())
      setError(null)
    } catch (e) {
      setError(e.message || 'RPC error')
    } finally {
      setLoading(false)
    }
  }, [fnName, JSON.stringify(params), enabled])

  useEffect(() => { fetch() }, [fetch])
  return { data, loading, error, refetch: fetch, lastFetch }
}

// ─── Animated counter ─────────────────────────────────────────────────────────
function AnimatedNumber({ value, duration = 1200 }) {
  const [displayed, setDisplayed] = useState(0)
  const rafRef = useRef(null)

  useEffect(() => {
    if (value == null || isNaN(value)) return
    const start = performance.now()
    const from  = displayed
    const to    = Number(value)

    const step = (now) => {
      const elapsed = now - start
      const progress = Math.min(elapsed / duration, 1)
      const eased = 1 - Math.pow(1 - progress, 3) // ease-out cubic
      setDisplayed(Math.round(from + (to - from) * eased))
      if (progress < 1) rafRef.current = requestAnimationFrame(step)
    }
    rafRef.current = requestAnimationFrame(step)
    return () => cancelAnimationFrame(rafRef.current)
  }, [value])

  return <>{displayed.toLocaleString()}</>
}

// ─── Skeleton shimmer ─────────────────────────────────────────────────────────
function Skeleton({ className = '' }) {
  return (
    <div
      className={`rounded-lg animate-pulse bg-gradient-to-r from-slate-200 via-slate-100 to-slate-200 ${className}`}
      style={{ backgroundSize: '200% 100%', animation: 'shimmer 1.5s infinite' }}
    />
  )
}

function ChartSkeleton({ h = 200 }) {
  return (
    <div className="space-y-3">
      <Skeleton className="h-4 w-1/3" />
      <Skeleton className={`w-full`} style={{ height: h }} />
    </div>
  )
}

// ─── Section wrapper ──────────────────────────────────────────────────────────
function Section({ title, subtitle, children, id }) {
  return (
    <section id={id} className="mb-12">
      <div className="mb-6">
        <h2 className="text-2xl font-bold text-slate-800 flex items-center gap-3">
          {title}
        </h2>
        {subtitle && <p className="text-slate-500 text-sm mt-1">{subtitle}</p>}
      </div>
      {children}
    </section>
  )
}

// ─── Card ─────────────────────────────────────────────────────────────────────
function Card({ children, className = '', style }) {
  return (
    <div
      className={`bg-white rounded-2xl shadow-card border border-slate-100 ${className}`}
      style={style}
    >
      {children}
    </div>
  )
}

// ─── Custom tooltip ───────────────────────────────────────────────────────────
function CustomTooltip({ active, payload, label, formatter }) {
  if (!active || !payload?.length) return null
  return (
    <div className="bg-slate-900 text-white px-3 py-2 rounded-lg shadow-xl text-xs border border-slate-700">
      {label != null && <div className="font-semibold text-slate-300 mb-1">{label}</div>}
      {payload.map((p, i) => (
        <div key={i} className="flex items-center gap-2">
          <span className="w-2 h-2 rounded-full inline-block" style={{ background: p.color || p.fill }} />
          <span className="text-slate-300">{p.name}:</span>
          <span className="font-bold">
            {formatter ? formatter(p.value) : p.value?.toLocaleString()}
          </span>
        </div>
      ))}
    </div>
  )
}

// ─── KPI Card ─────────────────────────────────────────────────────────────────
function KpiCard({ label, value, icon, gradient, sub }) {
  return (
    <Card className="overflow-hidden">
      <div
        className="h-1.5 w-full"
        style={{ background: gradient }}
      />
      <div className="p-5">
        <div className="flex items-start justify-between">
          <div>
            <p className="text-xs font-medium text-slate-500 uppercase tracking-wider mb-1">{label}</p>
            <p className="text-3xl font-extrabold text-slate-800 leading-none">
              {value == null ? <Skeleton className="h-9 w-28" /> : <AnimatedNumber value={value} />}
            </p>
            {sub && <p className="text-xs text-slate-400 mt-2">{sub}</p>}
          </div>
          <div
            className="w-12 h-12 rounded-xl flex items-center justify-center text-2xl flex-shrink-0"
            style={{ background: gradient, opacity: 0.15 + 0.85 }}
          >
            <span style={{ filter: 'drop-shadow(0 2px 4px rgba(0,0,0,0.2))' }}>{icon}</span>
          </div>
        </div>
      </div>
    </Card>
  )
}

// ─── "No data" state ─────────────────────────────────────────────────────────
function NoData({ message = 'No data available yet' }) {
  return (
    <div className="flex flex-col items-center justify-center h-40 text-slate-400">
      <svg className="w-10 h-10 mb-3 opacity-30" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M20 13V6a2 2 0 00-2-2H6a2 2 0 00-2 2v7m16 0v5a2 2 0 01-2 2H6a2 2 0 01-2-2v-5m16 0h-2.586a1 1 0 00-.707.293l-2.414 2.414a1 1 0 01-.707.293h-3.172a1 1 0 01-.707-.293l-2.414-2.414A1 1 0 006.586 13H4" />
      </svg>
      <p className="text-sm">{message}</p>
    </div>
  )
}

// ─── Section 1: KPI Cards ─────────────────────────────────────────────────────
function KpiSection({ data, loading }) {
  const kpis = [
    {
      label: 'Compounds',
      value: data?.compounds,
      icon: '⬡',
      gradient: 'linear-gradient(135deg, #1e3a5f, #3b82f6)',
      sub: `${(data?.compounds_with_activity || 0).toLocaleString()} with activity data`,
    },
    {
      label: 'Targets',
      value: data?.targets,
      icon: '🎯',
      gradient: 'linear-gradient(135deg, #14b8a6, #10b981)',
      sub: `${(data?.druggable_targets || 0).toLocaleString()} druggable`,
    },
    {
      label: 'Bioactivities',
      value: data?.bioactivities,
      icon: '📊',
      gradient: 'linear-gradient(135deg, #8b5cf6, #6366f1)',
      sub: `Avg pChEMBL: ${data?.avg_pchembl?.toFixed(2) || '—'}`,
    },
    {
      label: 'Approved Drugs',
      value: data?.approved_drugs,
      icon: '💊',
      gradient: 'linear-gradient(135deg, #22c55e, #16a34a)',
      sub: `${(data?.drug_mechanisms || 0).toLocaleString()} known mechanisms`,
    },
    {
      label: 'Assays',
      value: data?.assays,
      icon: '🔬',
      gradient: 'linear-gradient(135deg, #f59e0b, #f97316)',
      sub: `${(data?.drug_indications || 0).toLocaleString()} drug indications`,
    },
    {
      label: 'Database Size',
      value: data?.db_size ? parseInt(data.db_size) : null,
      icon: '🗄️',
      gradient: 'linear-gradient(135deg, #f43f5e, #e11d48)',
      sub: data?.db_size || 'Calculating...',
    },
  ]

  return (
    <div className="grid grid-cols-2 md:grid-cols-3 xl:grid-cols-6 gap-4">
      {loading
        ? Array.from({ length: 6 }).map((_, i) => (
            <Card key={i} className="p-5">
              <Skeleton className="h-4 w-20 mb-3" />
              <Skeleton className="h-9 w-28 mb-2" />
              <Skeleton className="h-3 w-24" />
            </Card>
          ))
        : kpis.map((kpi) => <KpiCard key={kpi.label} {...kpi} />)
      }
    </div>
  )
}

// ─── Section 2: Chemical Properties ──────────────────────────────────────────
function ChemicalPropertiesSection() {
  const mw   = useRPC('pharmaco_mw_distribution')
  const logp = useRPC('pharmaco_logp_distribution')

  const mwData   = (mw.data   || []).map(r => ({ name: `${r.mw_low}–${r.mw_high}`, count: r.count, bin: r.bin }))
  const logpData = (logp.data || []).map(r => ({ name: `${Number(r.logp_low).toFixed(1)}–${Number(r.logp_high).toFixed(1)}`, count: r.count }))

  return (
    <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
      {/* MW */}
      <Card className="p-6">
        <h3 className="font-semibold text-slate-700 mb-1">Molecular Weight Distribution</h3>
        <p className="text-xs text-slate-400 mb-5">Number of compounds per MW range (Da)</p>
        {mw.loading ? <ChartSkeleton h={220} /> : mwData.length === 0 ? <NoData /> : (
          <ResponsiveContainer width="100%" height={220}>
            <AreaChart data={mwData} margin={{ top: 5, right: 10, left: 0, bottom: 5 }}>
              <defs>
                <linearGradient id="mwGrad" x1="0" y1="0" x2="0" y2="1">
                  <stop offset="5%"  stopColor={C.blueLight} stopOpacity={0.4} />
                  <stop offset="95%" stopColor={C.blueLight} stopOpacity={0.02} />
                </linearGradient>
              </defs>
              <CartesianGrid strokeDasharray="3 3" stroke="#f1f5f9" />
              <XAxis dataKey="name" tick={{ fontSize: 10, fill: '#94a3b8' }} interval={2} />
              <YAxis tick={{ fontSize: 10, fill: '#94a3b8' }} tickFormatter={v => v >= 1000 ? `${(v/1000).toFixed(0)}k` : v} />
              <Tooltip content={<CustomTooltip formatter={v => v.toLocaleString() + ' compounds'} />} />
              <Area
                type="monotone" dataKey="count" name="Compounds"
                stroke={C.blueLight} strokeWidth={2}
                fill="url(#mwGrad)"
                dot={false} activeDot={{ r: 4, fill: C.blueLight }}
              />
              <ReferenceLine x="450–500" stroke={C.amber} strokeDasharray="4 2" label={{ value: 'Ro5 limit', fontSize: 9, fill: C.amber }} />
            </AreaChart>
          </ResponsiveContainer>
        )}
      </Card>

      {/* LogP */}
      <Card className="p-6">
        <h3 className="font-semibold text-slate-700 mb-1">LogP Distribution</h3>
        <p className="text-xs text-slate-400 mb-5">Lipophilicity profile across compounds</p>
        {logp.loading ? <ChartSkeleton h={220} /> : logpData.length === 0 ? <NoData /> : (
          <ResponsiveContainer width="100%" height={220}>
            <AreaChart data={logpData} margin={{ top: 5, right: 10, left: 0, bottom: 5 }}>
              <defs>
                <linearGradient id="logpGrad" x1="0" y1="0" x2="0" y2="1">
                  <stop offset="5%"  stopColor={C.teal} stopOpacity={0.4} />
                  <stop offset="95%" stopColor={C.teal} stopOpacity={0.02} />
                </linearGradient>
              </defs>
              <CartesianGrid strokeDasharray="3 3" stroke="#f1f5f9" />
              <XAxis dataKey="name" tick={{ fontSize: 10, fill: '#94a3b8' }} interval={2} />
              <YAxis tick={{ fontSize: 10, fill: '#94a3b8' }} tickFormatter={v => v >= 1000 ? `${(v/1000).toFixed(0)}k` : v} />
              <Tooltip content={<CustomTooltip formatter={v => v.toLocaleString() + ' compounds'} />} />
              <Area
                type="monotone" dataKey="count" name="Compounds"
                stroke={C.teal} strokeWidth={2}
                fill="url(#logpGrad)"
                dot={false} activeDot={{ r: 4, fill: C.teal }}
              />
              <ReferenceLine x="4.5–5.0" stroke={C.amber} strokeDasharray="4 2" label={{ value: 'Ro5 limit', fontSize: 9, fill: C.amber }} />
            </AreaChart>
          </ResponsiveContainer>
        )}
      </Card>
    </div>
  )
}

// ─── Section 3: Drug-Likeness ─────────────────────────────────────────────────
const RO5_COLORS = [C.green, C.emerald, C.amber, C.orange, C.rose, '#7f1d1d']

function DrugLikenessSection() {
  const { data, loading } = useRPC('pharmaco_druglikeness')

  const ro5Data = (data?.ro5_distribution || []).map((r, i) => ({
    name: r.violations === 0 ? 'Fully compliant' : `${r.violations} violation${r.violations > 1 ? 's' : ''}`,
    value: r.count,
    violations: r.violations,
  }))

  const qedData = (data?.qed_distribution || []).map(r => ({
    name: `${Number(r.qed_low).toFixed(1)}–${Number(r.qed_high).toFixed(1)}`,
    count: r.count,
  }))

  const avgQED     = data?.avg_qed ? Number(data.avg_qed).toFixed(3) : null
  const pctDrugLike = data?.pct_drug_like ? Number(data.pct_drug_like).toFixed(1) : null

  const RADIAN = Math.PI / 180
  const renderCustomPieLabel = ({ cx, cy, midAngle, innerRadius, outerRadius, percent, name }) => {
    if (percent < 0.04) return null
    const radius = innerRadius + (outerRadius - innerRadius) * 0.55
    const x = cx + radius * Math.cos(-midAngle * RADIAN)
    const y = cy + radius * Math.sin(-midAngle * RADIAN)
    return (
      <text x={x} y={y} fill="white" textAnchor="middle" dominantBaseline="central" fontSize={10} fontWeight={600}>
        {`${(percent * 100).toFixed(0)}%`}
      </text>
    )
  }

  return (
    <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
      {/* Ro5 Pie */}
      <Card className="p-6">
        <h3 className="font-semibold text-slate-700 mb-1">Lipinski Ro5 Violations</h3>
        <p className="text-xs text-slate-400 mb-3">Drug-likeness compliance</p>
        {loading ? <ChartSkeleton h={220} /> : ro5Data.length === 0 ? <NoData /> : (
          <>
            <ResponsiveContainer width="100%" height={180}>
              <PieChart>
                <Pie
                  data={ro5Data} cx="50%" cy="50%"
                  innerRadius={45} outerRadius={80}
                  paddingAngle={2}
                  dataKey="value"
                  labelLine={false}
                  label={renderCustomPieLabel}
                >
                  {ro5Data.map((_, i) => (
                    <Cell key={i} fill={RO5_COLORS[i % RO5_COLORS.length]} />
                  ))}
                </Pie>
                <Tooltip content={<CustomTooltip formatter={v => v.toLocaleString() + ' compounds'} />} />
              </PieChart>
            </ResponsiveContainer>
            <div className="mt-2 space-y-1">
              {ro5Data.map((r, i) => (
                <div key={i} className="flex items-center justify-between text-xs">
                  <div className="flex items-center gap-2">
                    <span className="w-2.5 h-2.5 rounded-sm" style={{ background: RO5_COLORS[i % RO5_COLORS.length] }} />
                    <span className="text-slate-600">{r.name}</span>
                  </div>
                  <span className="font-semibold text-slate-700">{r.value.toLocaleString()}</span>
                </div>
              ))}
            </div>
          </>
        )}
      </Card>

      {/* QED Distribution */}
      <Card className="p-6">
        <h3 className="font-semibold text-slate-700 mb-1">QED Score Distribution</h3>
        <p className="text-xs text-slate-400 mb-5">Quantitative Estimate of Drug-likeness</p>
        {loading ? <ChartSkeleton h={220} /> : qedData.length === 0 ? <NoData /> : (
          <ResponsiveContainer width="100%" height={220}>
            <AreaChart data={qedData} margin={{ top: 5, right: 10, left: 0, bottom: 5 }}>
              <defs>
                <linearGradient id="qedGrad" x1="0" y1="0" x2="0" y2="1">
                  <stop offset="5%"  stopColor={C.violet} stopOpacity={0.45} />
                  <stop offset="95%" stopColor={C.violet} stopOpacity={0.02} />
                </linearGradient>
              </defs>
              <CartesianGrid strokeDasharray="3 3" stroke="#f1f5f9" />
              <XAxis dataKey="name" tick={{ fontSize: 10, fill: '#94a3b8' }} />
              <YAxis tick={{ fontSize: 10, fill: '#94a3b8' }} tickFormatter={v => v >= 1000 ? `${(v/1000).toFixed(0)}k` : v} />
              <Tooltip content={<CustomTooltip formatter={v => v.toLocaleString() + ' compounds'} />} />
              <Area
                type="monotone" dataKey="count" name="Compounds"
                stroke={C.violet} strokeWidth={2}
                fill="url(#qedGrad)"
                dot={false} activeDot={{ r: 4, fill: C.violet }}
              />
            </AreaChart>
          </ResponsiveContainer>
        )}
      </Card>

      {/* Summary stats */}
      <Card className="p-6 flex flex-col justify-between" style={{ background: 'linear-gradient(135deg, #1e3a5f 0%, #2d5a8e 100%)' }}>
        <div>
          <h3 className="font-semibold text-white mb-1">Drug-Likeness Summary</h3>
          <p className="text-xs text-blue-200 mb-6">Key drug-like property metrics</p>
        </div>
        {loading ? (
          <div className="space-y-4">
            <Skeleton className="h-16 w-full opacity-20" />
            <Skeleton className="h-16 w-full opacity-20" />
          </div>
        ) : (
          <div className="space-y-5">
            <div className="bg-white/10 rounded-xl p-4 border border-white/20">
              <p className="text-blue-200 text-xs uppercase tracking-wider mb-1">Average QED Score</p>
              <p className="text-4xl font-black text-white">{avgQED ?? '—'}</p>
              <p className="text-blue-300 text-xs mt-1">0 = least drug-like, 1 = most drug-like</p>
            </div>
            <div className="bg-white/10 rounded-xl p-4 border border-white/20">
              <p className="text-blue-200 text-xs uppercase tracking-wider mb-1">Drug-Like Compounds</p>
              <p className="text-4xl font-black text-green-400">{pctDrugLike != null ? `${pctDrugLike}%` : '—'}</p>
              <p className="text-blue-300 text-xs mt-1">Pass all Lipinski Ro5 criteria</p>
              {pctDrugLike != null && (
                <div className="mt-3 h-2 bg-white/10 rounded-full overflow-hidden">
                  <div
                    className="h-full rounded-full transition-all duration-1000"
                    style={{ width: `${pctDrugLike}%`, background: 'linear-gradient(90deg, #22c55e, #10b981)' }}
                  />
                </div>
              )}
            </div>
          </div>
        )}
      </Card>
    </div>
  )
}

// ─── Section 4: Target Landscape ──────────────────────────────────────────────
function TargetLandscapeSection() {
  const families = useRPC('pharmaco_target_families')
  const types    = useRPC('pharmaco_target_types')

  const familyData = (families.data || [])
    .slice(0, 15)
    .map(r => ({ name: r.family || 'Unknown', count: r.count }))
    .reverse()

  const typeData = (types.data || [])
    .map((r, i) => ({ name: r.type || 'Unknown', value: r.count, fill: FAMILY_COLORS[i % FAMILY_COLORS.length] }))

  return (
    <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
      {/* Target Families */}
      <Card className="p-6">
        <h3 className="font-semibold text-slate-700 mb-1">Target Protein Families</h3>
        <p className="text-xs text-slate-400 mb-4">Top 15 protein families by target count</p>
        {families.loading ? <ChartSkeleton h={380} /> : familyData.length === 0 ? <NoData /> : (
          <ResponsiveContainer width="100%" height={380}>
            <BarChart data={familyData} layout="vertical" margin={{ top: 0, right: 40, left: 0, bottom: 0 }}>
              <CartesianGrid strokeDasharray="3 3" stroke="#f1f5f9" horizontal={false} />
              <XAxis type="number" tick={{ fontSize: 10, fill: '#94a3b8' }} tickFormatter={v => v.toLocaleString()} />
              <YAxis
                type="category" dataKey="name" width={140}
                tick={{ fontSize: 10, fill: '#64748b' }}
                tickFormatter={v => v.length > 20 ? v.slice(0, 20) + '…' : v}
              />
              <Tooltip content={<CustomTooltip formatter={v => v.toLocaleString() + ' targets'} />} />
              <Bar dataKey="count" name="Targets" radius={[0, 4, 4, 0]}>
                {familyData.map((_, i) => (
                  <Cell key={i} fill={FAMILY_COLORS[i % FAMILY_COLORS.length]} fillOpacity={0.85} />
                ))}
                <LabelList dataKey="count" position="right" style={{ fontSize: 10, fill: '#94a3b8' }} formatter={v => v.toLocaleString()} />
              </Bar>
            </BarChart>
          </ResponsiveContainer>
        )}
      </Card>

      {/* Target Types */}
      <Card className="p-6">
        <h3 className="font-semibold text-slate-700 mb-1">Target Types</h3>
        <p className="text-xs text-slate-400 mb-4">Distribution of target categories</p>
        {types.loading ? <ChartSkeleton h={380} /> : typeData.length === 0 ? <NoData /> : (
          <>
            <ResponsiveContainer width="100%" height={260}>
              <PieChart>
                <Pie
                  data={typeData} cx="50%" cy="50%"
                  outerRadius={110} innerRadius={55}
                  paddingAngle={3}
                  dataKey="value"
                  label={({ name, percent }) => percent > 0.04 ? `${(percent * 100).toFixed(0)}%` : ''}
                  labelLine={{ stroke: '#cbd5e1', strokeWidth: 1 }}
                >
                  {typeData.map((entry, i) => (
                    <Cell key={i} fill={entry.fill} />
                  ))}
                </Pie>
                <Tooltip content={<CustomTooltip formatter={v => v.toLocaleString() + ' targets'} />} />
              </PieChart>
            </ResponsiveContainer>
            <div className="mt-4 grid grid-cols-2 gap-2">
              {typeData.map((t, i) => (
                <div key={i} className="flex items-center gap-2 text-xs">
                  <span className="w-2.5 h-2.5 rounded-full flex-shrink-0" style={{ background: t.fill }} />
                  <span className="text-slate-600 truncate">{t.name}</span>
                  <span className="ml-auto font-semibold text-slate-700">{t.value.toLocaleString()}</span>
                </div>
              ))}
            </div>
          </>
        )}
      </Card>
    </div>
  )
}

// ─── Section 5: Bioactivity Analysis ─────────────────────────────────────────
function BioactivitySection() {
  const pchembl  = useRPC('pharmaco_pchembl_distribution')
  const actTypes = useRPC('pharmaco_activity_types')

  // Color each pChEMBL bin from gray (inactive) to vivid green (potent)
  const pchemblData = (pchembl.data || []).map((r, i, arr) => {
    const t = i / Math.max(arr.length - 1, 1)
    const r_ = Math.round(148 - t * 148)
    const g_ = Math.round(163 + t * (197 - 163))
    const b_ = Math.round(184 - t * 90)
    return {
      name: `${Number(r.pchembl_low).toFixed(1)}–${Number(r.pchembl_high).toFixed(1)}`,
      count: r.count,
      fill: `rgb(${r_},${g_},${b_})`,
    }
  })

  const actData = (actTypes.data || [])
    .slice(0, 15)
    .map(r => ({ name: r.type || 'Unknown', count: r.count, avg: Number(r.avg_pchembl || 0).toFixed(2) }))
    .reverse()

  return (
    <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
      {/* pChEMBL */}
      <Card className="p-6">
        <h3 className="font-semibold text-slate-700 mb-1">Potency Distribution (pChEMBL)</h3>
        <p className="text-xs text-slate-400 mb-5">Higher = more potent. &gt;6 considered active.</p>
        {pchembl.loading ? <ChartSkeleton h={260} /> : pchemblData.length === 0 ? <NoData /> : (
          <ResponsiveContainer width="100%" height={260}>
            <BarChart data={pchemblData} margin={{ top: 5, right: 10, left: 0, bottom: 5 }}>
              <CartesianGrid strokeDasharray="3 3" stroke="#f1f5f9" />
              <XAxis dataKey="name" tick={{ fontSize: 9, fill: '#94a3b8' }} interval={1} />
              <YAxis tick={{ fontSize: 10, fill: '#94a3b8' }} tickFormatter={v => v >= 1000 ? `${(v/1000).toFixed(0)}k` : v} />
              <Tooltip content={<CustomTooltip formatter={v => v.toLocaleString() + ' assays'} />} />
              <ReferenceLine x="5.5–6.0" stroke={C.amber} strokeDasharray="4 2" label={{ value: 'Active threshold', fontSize: 9, fill: C.amber }} />
              <Bar dataKey="count" name="Bioactivities" radius={[3, 3, 0, 0]}>
                {pchemblData.map((d, i) => (
                  <Cell key={i} fill={d.fill} />
                ))}
              </Bar>
            </BarChart>
          </ResponsiveContainer>
        )}
      </Card>

      {/* Activity types */}
      <Card className="p-6">
        <h3 className="font-semibold text-slate-700 mb-1">Activity Measurement Types</h3>
        <p className="text-xs text-slate-400 mb-4">Top 15 assay measurement types</p>
        {actTypes.loading ? <ChartSkeleton h={260} /> : actData.length === 0 ? <NoData /> : (
          <ResponsiveContainer width="100%" height={340}>
            <BarChart data={actData} layout="vertical" margin={{ top: 0, right: 60, left: 0, bottom: 0 }}>
              <CartesianGrid strokeDasharray="3 3" stroke="#f1f5f9" horizontal={false} />
              <XAxis type="number" tick={{ fontSize: 10, fill: '#94a3b8' }} tickFormatter={v => v >= 1000 ? `${(v/1000).toFixed(0)}k` : v} />
              <YAxis type="category" dataKey="name" width={60} tick={{ fontSize: 10, fill: '#64748b' }} />
              <Tooltip
                content={({ active, payload }) =>
                  active && payload?.length ? (
                    <div className="bg-slate-900 text-white px-3 py-2 rounded-lg shadow-xl text-xs">
                      <div className="font-semibold">{payload[0]?.payload?.name}</div>
                      <div>Count: <b>{payload[0]?.value?.toLocaleString()}</b></div>
                      <div>Avg pChEMBL: <b>{payload[0]?.payload?.avg}</b></div>
                    </div>
                  ) : null
                }
              />
              <Bar dataKey="count" name="Bioactivities" radius={[0, 4, 4, 0]} fill={C.emerald} fillOpacity={0.85}>
                <LabelList dataKey="count" position="right" style={{ fontSize: 9, fill: '#94a3b8' }} formatter={v => v >= 1000 ? `${(v/1000).toFixed(1)}k` : v} />
              </Bar>
            </BarChart>
          </ResponsiveContainer>
        )}
      </Card>
    </div>
  )
}

// ─── Section 6: Drug Pipeline ─────────────────────────────────────────────────
function DrugPipelineSection() {
  const phases = useRPC('pharmaco_clinical_phases')
  const mechs  = useRPC('pharmaco_mechanisms')

  const phaseData = (phases.data || []).map(r => ({
    name: r.label || `Phase ${r.phase}`,
    count: r.count,
    phase: String(r.phase),
    fill: PHASE_COLORS[String(r.phase)] || '#94a3b8',
  }))

  const mechData = (mechs.data || [])
    .slice(0, 15)
    .map(r => ({ name: r.action || 'Unknown', count: r.count }))
    .reverse()

  const maxPhase = Math.max(...phaseData.map(p => p.count), 1)

  return (
    <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
      {/* Funnel (custom) */}
      <Card className="p-6">
        <h3 className="font-semibold text-slate-700 mb-1">Clinical Development Pipeline</h3>
        <p className="text-xs text-slate-400 mb-5">Drugs by clinical phase</p>
        {phases.loading ? <ChartSkeleton h={260} /> : phaseData.length === 0 ? <NoData /> : (
          <div className="space-y-3">
            {phaseData.map((p) => {
              const pct = (p.count / maxPhase) * 100
              return (
                <div key={p.phase} className="group">
                  <div className="flex items-center justify-between mb-1.5">
                    <span className="text-sm font-medium text-slate-700">{p.name}</span>
                    <span className="text-sm font-bold" style={{ color: p.fill }}>
                      {p.count.toLocaleString()}
                    </span>
                  </div>
                  <div className="h-8 bg-slate-100 rounded-lg overflow-hidden">
                    <div
                      className="h-full rounded-lg flex items-center px-3 transition-all duration-700"
                      style={{ width: `${Math.max(pct, 4)}%`, background: p.fill, opacity: 0.9 }}
                    >
                      {pct > 20 && (
                        <span className="text-white text-xs font-semibold">
                          {pct.toFixed(0)}%
                        </span>
                      )}
                    </div>
                  </div>
                </div>
              )
            })}
          </div>
        )}
      </Card>

      {/* Mechanisms */}
      <Card className="p-6">
        <h3 className="font-semibold text-slate-700 mb-1">Mechanisms of Action</h3>
        <p className="text-xs text-slate-400 mb-4">Top 15 drug action types</p>
        {mechs.loading ? <ChartSkeleton h={340} /> : mechData.length === 0 ? <NoData /> : (
          <ResponsiveContainer width="100%" height={340}>
            <BarChart data={mechData} layout="vertical" margin={{ top: 0, right: 50, left: 0, bottom: 0 }}>
              <CartesianGrid strokeDasharray="3 3" stroke="#f1f5f9" horizontal={false} />
              <XAxis type="number" tick={{ fontSize: 10, fill: '#94a3b8' }} tickFormatter={v => v.toLocaleString()} />
              <YAxis
                type="category" dataKey="name" width={100}
                tick={{ fontSize: 10, fill: '#64748b' }}
                tickFormatter={v => v.length > 16 ? v.slice(0, 16) + '…' : v}
              />
              <Tooltip content={<CustomTooltip formatter={v => v.toLocaleString() + ' drugs'} />} />
              <Bar dataKey="count" name="Drugs" radius={[0, 4, 4, 0]}>
                {mechData.map((_, i) => (
                  <Cell key={i} fill={FAMILY_COLORS[i % FAMILY_COLORS.length]} fillOpacity={0.85} />
                ))}
                <LabelList dataKey="count" position="right" style={{ fontSize: 9, fill: '#94a3b8' }} formatter={v => v.toLocaleString()} />
              </Bar>
            </BarChart>
          </ResponsiveContainer>
        )}
      </Card>
    </div>
  )
}

// ─── Section 7: Chemical Space ────────────────────────────────────────────────
function ChemicalSpaceSection() {
  const { data, loading } = useRPC('pharmaco_chemical_space')

  // Color QED: 0=red, 0.5=amber, 1=green
  function qedColor(qed) {
    if (qed == null) return '#94a3b8'
    const t = Math.min(Math.max(Number(qed), 0), 1)
    if (t < 0.5) {
      // red → amber
      const f = t * 2
      const r = 244, g = Math.round(f * 159), b = Math.round(f * 30)
      return `rgb(${r},${g},${b})`
    } else {
      // amber → green
      const f = (t - 0.5) * 2
      const r = Math.round(244 - f * 210), g = Math.round(159 + f * 36), b = Math.round(30 - f * 30 + f * 94)
      return `rgb(${r},${g},${b})`
    }
  }

  const points = (data || []).map(r => ({
    mw:   Number(r.mw   || 0),
    logp: Number(r.logp || 0),
    qed:  Number(r.qed  || 0),
    phase: String(r.phase ?? ''),
    color: qedColor(r.qed),
  }))

  const CustomDot = (props) => {
    const { cx, cy, payload } = props
    if (cx == null || cy == null) return null
    return (
      <circle
        cx={cx} cy={cy} r={3.5}
        fill={payload.color}
        fillOpacity={0.7}
        stroke="white"
        strokeWidth={0.5}
      />
    )
  }

  const CustomScatterTooltip = ({ active, payload }) => {
    if (!active || !payload?.length) return null
    const d = payload[0]?.payload
    return (
      <div className="bg-slate-900 text-white px-3 py-2 rounded-lg shadow-xl text-xs border border-slate-700">
        <div className="grid grid-cols-2 gap-x-3 gap-y-1">
          <span className="text-slate-400">MW</span><span className="font-bold">{d.mw.toFixed(1)} Da</span>
          <span className="text-slate-400">LogP</span><span className="font-bold">{d.logp.toFixed(2)}</span>
          <span className="text-slate-400">QED</span><span className="font-bold">{d.qed.toFixed(3)}</span>
          {d.phase && <><span className="text-slate-400">Phase</span><span className="font-bold">{d.phase}</span></>}
        </div>
      </div>
    )
  }

  // QED gradient legend
  const legendStops = Array.from({ length: 10 }, (_, i) => ({
    qed: i / 9,
    color: qedColor(i / 9),
    label: (i / 9).toFixed(1),
  }))

  return (
    <Card className="p-6">
      <div className="flex items-start justify-between mb-2">
        <div>
          <h3 className="font-semibold text-slate-700">Chemical Space Explorer</h3>
          <p className="text-xs text-slate-400 mt-1">
            MW vs LogP for 2,000 sampled compounds. Color encodes QED score (green = drug-like, red = non-drug-like).
            Dashed lines mark Lipinski Ro5 limits (MW=500, LogP=5).
          </p>
        </div>
        <div className="flex-shrink-0 ml-6">
          <p className="text-xs text-slate-400 mb-1 text-center">QED Score</p>
          <div className="flex items-center gap-1">
            <span className="text-xs text-red-500 font-medium">0</span>
            <div
              className="h-3 w-32 rounded-full"
              style={{
                background: 'linear-gradient(to right, rgb(244,0,0), rgb(244,159,30), rgb(34,197,94))',
              }}
            />
            <span className="text-xs text-green-600 font-medium">1</span>
          </div>
        </div>
      </div>
      {loading ? (
        <Skeleton className="w-full h-96 mt-4" />
      ) : points.length === 0 ? (
        <NoData message="No chemical space data yet" />
      ) : (
        <ResponsiveContainer width="100%" height={420}>
          <ScatterChart margin={{ top: 10, right: 30, left: 10, bottom: 30 }}>
            <CartesianGrid strokeDasharray="3 3" stroke="#f1f5f9" />
            <XAxis
              type="number" dataKey="mw" name="MW" domain={[0, 900]}
              tick={{ fontSize: 11, fill: '#94a3b8' }}
              label={{ value: 'Molecular Weight (Da)', position: 'insideBottom', offset: -15, fontSize: 12, fill: '#64748b' }}
            />
            <YAxis
              type="number" dataKey="logp" name="LogP" domain={[-5, 12]}
              tick={{ fontSize: 11, fill: '#94a3b8' }}
              label={{ value: 'LogP', angle: -90, position: 'insideLeft', offset: 10, fontSize: 12, fill: '#64748b' }}
            />
            <ZAxis range={[25, 25]} />
            <Tooltip content={<CustomScatterTooltip />} />
            {/* Lipinski reference lines */}
            <ReferenceLine
              x={500} stroke={C.amber} strokeDasharray="6 3" strokeWidth={1.5}
              label={{ value: 'MW 500', position: 'top', fontSize: 10, fill: C.amber }}
            />
            <ReferenceLine
              y={5} stroke={C.amber} strokeDasharray="6 3" strokeWidth={1.5}
              label={{ value: 'LogP 5', position: 'right', fontSize: 10, fill: C.amber }}
            />
            {/* Render as individual scatter points with custom dots */}
            <Scatter
              data={points}
              shape={<CustomDot />}
              isAnimationActive={false}
            />
          </ScatterChart>
        </ResponsiveContainer>
      )}
      {!loading && points.length > 0 && (
        <p className="text-xs text-slate-400 text-right mt-1">{points.length.toLocaleString()} compounds shown</p>
      )}
    </Card>
  )
}

// ─── Section 8: Ingestion Status ──────────────────────────────────────────────
const STATUS_STYLES = {
  completed: { bg: 'bg-green-100', text: 'text-green-700', dot: 'bg-green-500', label: 'Completed' },
  running:   { bg: 'bg-blue-100',  text: 'text-blue-700',  dot: 'bg-blue-500 animate-pulse', label: 'Running' },
  pending:   { bg: 'bg-amber-100', text: 'text-amber-700', dot: 'bg-amber-400', label: 'Pending' },
  failed:    { bg: 'bg-red-100',   text: 'text-red-700',   dot: 'bg-red-500', label: 'Failed' },
}

function IngestionStatusSection() {
  const { data, loading } = useRPC('pharmaco_ingestion_status')

  const rows = data || []

  function fmtDate(d) {
    if (!d) return '—'
    return new Date(d).toLocaleString()
  }

  function statusBadge(status) {
    const s = STATUS_STYLES[status?.toLowerCase()] || { bg: 'bg-slate-100', text: 'text-slate-600', dot: 'bg-slate-400', label: status }
    return (
      <span className={`inline-flex items-center gap-1.5 px-2.5 py-1 rounded-full text-xs font-medium ${s.bg} ${s.text}`}>
        <span className={`w-1.5 h-1.5 rounded-full ${s.dot}`} />
        {s.label}
      </span>
    )
  }

  return (
    <Card className="overflow-hidden">
      <div className="px-6 py-4 border-b border-slate-100 flex items-center justify-between">
        <div>
          <h3 className="font-semibold text-slate-700">Data Ingestion Status</h3>
          <p className="text-xs text-slate-400 mt-0.5">Pipeline steps and import progress</p>
        </div>
        {rows.some(r => r.status?.toLowerCase() === 'running') && (
          <span className="inline-flex items-center gap-2 text-xs text-blue-600 bg-blue-50 px-3 py-1.5 rounded-full border border-blue-200">
            <span className="w-2 h-2 rounded-full bg-blue-500 animate-pulse" />
            Ingestion in progress
          </span>
        )}
      </div>
      {loading ? (
        <div className="p-6 space-y-3">
          {Array.from({ length: 5 }).map((_, i) => (
            <Skeleton key={i} className="h-10 w-full" />
          ))}
        </div>
      ) : rows.length === 0 ? (
        <div className="p-6">
          <NoData message="No ingestion records found" />
        </div>
      ) : (
        <div className="overflow-x-auto">
          <table className="w-full text-sm">
            <thead>
              <tr className="bg-slate-50 border-b border-slate-100">
                <th className="text-left px-6 py-3 text-xs font-semibold text-slate-500 uppercase tracking-wider">Source</th>
                <th className="text-left px-6 py-3 text-xs font-semibold text-slate-500 uppercase tracking-wider">Step</th>
                <th className="text-left px-6 py-3 text-xs font-semibold text-slate-500 uppercase tracking-wider">Status</th>
                <th className="text-right px-6 py-3 text-xs font-semibold text-slate-500 uppercase tracking-wider">Rows Inserted</th>
                <th className="text-left px-6 py-3 text-xs font-semibold text-slate-500 uppercase tracking-wider">Started</th>
                <th className="text-left px-6 py-3 text-xs font-semibold text-slate-500 uppercase tracking-wider">Completed</th>
              </tr>
            </thead>
            <tbody>
              {rows.map((r, i) => (
                <tr key={i} className={`border-b border-slate-50 hover:bg-slate-50 transition-colors ${i % 2 === 0 ? '' : 'bg-slate-50/40'}`}>
                  <td className="px-6 py-3 font-medium text-slate-700">{r.source}</td>
                  <td className="px-6 py-3 text-slate-600">{r.step}</td>
                  <td className="px-6 py-3">{statusBadge(r.status)}</td>
                  <td className="px-6 py-3 text-right font-mono text-slate-700">
                    {r.rows_inserted != null ? r.rows_inserted.toLocaleString() : '—'}
                  </td>
                  <td className="px-6 py-3 text-xs text-slate-500">{fmtDate(r.started_at)}</td>
                  <td className="px-6 py-3 text-xs text-slate-500">{fmtDate(r.completed_at)}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
    </Card>
  )
}

// ─── Scroll spy nav ───────────────────────────────────────────────────────────
const NAV_SECTIONS = [
  { id: 'kpi',        label: 'Overview' },
  { id: 'chemical',   label: 'Chemical Props' },
  { id: 'druglike',   label: 'Drug-Likeness' },
  { id: 'targets',    label: 'Targets' },
  { id: 'bioactivity',label: 'Bioactivity' },
  { id: 'pipeline',   label: 'Pipeline' },
  { id: 'space',      label: 'Chem. Space' },
  { id: 'ingestion',  label: 'Ingestion' },
]

function FloatingNav() {
  const [active, setActive] = useState('kpi')

  useEffect(() => {
    const observer = new IntersectionObserver(
      entries => {
        entries.forEach(e => { if (e.isIntersecting) setActive(e.target.id) })
      },
      { threshold: 0.3, rootMargin: '-20% 0px -60% 0px' }
    )
    NAV_SECTIONS.forEach(s => {
      const el = document.getElementById(s.id)
      if (el) observer.observe(el)
    })
    return () => observer.disconnect()
  }, [])

  const scrollTo = (id) => {
    document.getElementById(id)?.scrollIntoView({ behavior: 'smooth', block: 'start' })
  }

  return (
    <nav className="fixed left-4 top-1/2 -translate-y-1/2 z-50 hidden xl:flex flex-col gap-1">
      {NAV_SECTIONS.map(s => (
        <button
          key={s.id}
          onClick={() => scrollTo(s.id)}
          className={`flex items-center gap-2 px-3 py-1.5 rounded-lg text-xs font-medium transition-all ${
            active === s.id
              ? 'bg-[#1e3a5f] text-white shadow-md'
              : 'bg-white/80 text-slate-500 hover:bg-white hover:text-slate-700 shadow-sm'
          }`}
        >
          <span className={`w-1.5 h-1.5 rounded-full ${active === s.id ? 'bg-green-400' : 'bg-slate-300'}`} />
          {s.label}
        </button>
      ))}
    </nav>
  )
}

// ─── Main Dashboard ───────────────────────────────────────────────────────────
export default function PharmacoDBDashboard() {
  const overview = useRPC('pharmaco_overview')
  const [lastRefresh, setLastRefresh] = useState(new Date())
  const [refreshKey, setRefreshKey] = useState(0)

  // Auto-refresh every 60s
  useEffect(() => {
    const id = setInterval(() => {
      setRefreshKey(k => k + 1)
      setLastRefresh(new Date())
    }, 60_000)
    return () => clearInterval(id)
  }, [])

  const handleRefresh = () => {
    setRefreshKey(k => k + 1)
    setLastRefresh(new Date())
    overview.refetch()
  }

  const isIngesting = overview.data &&
    (overview.data.bioactivities === 0 || overview.data.compounds === 0)

  return (
    <div className="min-h-screen bg-slate-50">
      <FloatingNav />

      {/* ── Header ── */}
      <header
        style={{ background: 'linear-gradient(135deg, #0f2540 0%, #1e3a5f 40%, #1a4a7a 70%, #0f3d2e 100%)' }}
        className="relative overflow-hidden"
      >
        {/* Background decoration */}
        <div className="absolute inset-0 pointer-events-none">
          <div className="absolute top-0 right-0 w-96 h-96 rounded-full opacity-5" style={{ background: C.green, transform: 'translate(30%, -40%)' }} />
          <div className="absolute bottom-0 left-0 w-64 h-64 rounded-full opacity-5" style={{ background: C.blueLight, transform: 'translate(-30%, 40%)' }} />
          {/* Grid dots */}
          <svg className="absolute inset-0 w-full h-full opacity-5" xmlns="http://www.w3.org/2000/svg">
            <defs>
              <pattern id="grid" width="40" height="40" patternUnits="userSpaceOnUse">
                <circle cx="1" cy="1" r="1" fill="white" />
              </pattern>
            </defs>
            <rect width="100%" height="100%" fill="url(#grid)" />
          </svg>
        </div>

        <div className="relative max-w-screen-2xl mx-auto px-8 xl:px-24 py-10">
          <div className="flex items-start justify-between flex-wrap gap-6">
            <div>
              <div className="flex items-center gap-3 mb-2">
                <div className="w-10 h-10 rounded-xl bg-green-500/20 border border-green-500/30 flex items-center justify-center">
                  <svg className="w-5 h-5 text-green-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 3H5a2 2 0 00-2 2v4m6-6h10a2 2 0 012 2v4M9 3v18m0 0h10a2 2 0 002-2V9M9 21H5a2 2 0 01-2-2V9m0 0h18" />
                  </svg>
                </div>
                <span className="text-green-400 text-sm font-semibold tracking-wider uppercase">PharmacoDB</span>
              </div>
              <h1 className="text-4xl xl:text-5xl font-black text-white leading-tight">
                Pharmacological
                <br />
                <span style={{ WebkitTextFillColor: 'transparent', WebkitBackgroundClip: 'text', backgroundClip: 'text', background: 'linear-gradient(90deg, #22c55e, #14b8a6)' }}>
                  Database Dashboard
                </span>
              </h1>
              <p className="text-blue-200 mt-3 text-sm max-w-lg">
                Comprehensive pharmacological data covering compounds, targets, bioactivities,
                and drug mechanisms sourced from ChEMBL, UniProt, and major pharmacological databases.
              </p>
            </div>

            {/* Controls */}
            <div className="flex flex-col items-end gap-3">
              {isIngesting && (
                <div className="flex items-center gap-2 bg-amber-400/10 border border-amber-400/30 text-amber-300 px-4 py-2 rounded-xl text-sm">
                  <span className="w-2 h-2 rounded-full bg-amber-400 animate-pulse" />
                  Data ingestion in progress
                </div>
              )}
              <button
                onClick={handleRefresh}
                className="flex items-center gap-2 bg-white/10 hover:bg-white/20 border border-white/20 text-white px-4 py-2.5 rounded-xl text-sm font-medium transition-all"
              >
                <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
                </svg>
                Refresh
              </button>
              <p className="text-blue-300 text-xs">
                Last updated: {lastRefresh.toLocaleTimeString()} · Auto-refresh: 60s
              </p>
            </div>
          </div>

          {/* Quick stats bar */}
          {!overview.loading && overview.data && (
            <div className="mt-8 flex flex-wrap gap-6">
              {[
                { label: 'Compounds', value: overview.data.compounds, color: '#60a5fa' },
                { label: 'Targets',   value: overview.data.targets,   color: '#34d399' },
                { label: 'Bioactivities', value: overview.data.bioactivities, color: '#a78bfa' },
                { label: 'Approved Drugs', value: overview.data.approved_drugs, color: '#fbbf24' },
              ].map(s => (
                <div key={s.label} className="flex items-center gap-2">
                  <span className="w-2 h-2 rounded-full" style={{ background: s.color }} />
                  <span className="text-white font-bold text-sm">{s.value?.toLocaleString() ?? '—'}</span>
                  <span className="text-blue-300 text-sm">{s.label}</span>
                </div>
              ))}
            </div>
          )}
        </div>
      </header>

      {/* ── Content ── */}
      <main className="max-w-screen-2xl mx-auto px-4 xl:px-24 py-10 space-y-12">

        {/* Section 1: KPIs */}
        <div id="kpi" className="scroll-mt-6">
          <Section title="Overview" subtitle="Key database statistics at a glance">
            <KpiSection data={overview.data} loading={overview.loading} />
          </Section>
        </div>

        {/* Section 2: Chemical Properties */}
        <div id="chemical" className="scroll-mt-6">
          <Section
            title="Chemical Property Distributions"
            subtitle="Molecular weight and lipophilicity profiles across all compounds"
          >
            <ChemicalPropertiesSection key={`chem-${refreshKey}`} />
          </Section>
        </div>

        {/* Section 3: Drug-Likeness */}
        <div id="druglike" className="scroll-mt-6">
          <Section
            title="Drug-Likeness Analysis"
            subtitle="Lipinski Ro5 compliance and Quantitative Estimate of Drug-likeness (QED)"
          >
            <DrugLikenessSection key={`dl-${refreshKey}`} />
          </Section>
        </div>

        {/* Section 4: Target Landscape */}
        <div id="targets" className="scroll-mt-6">
          <Section
            title="Target Landscape"
            subtitle="Protein families and target classifications across the druggable genome"
          >
            <TargetLandscapeSection key={`tl-${refreshKey}`} />
          </Section>
        </div>

        {/* Section 5: Bioactivity */}
        <div id="bioactivity" className="scroll-mt-6">
          <Section
            title="Bioactivity Analysis"
            subtitle="Potency distribution and assay type breakdown for all bioactivity measurements"
          >
            <BioactivitySection key={`bio-${refreshKey}`} />
          </Section>
        </div>

        {/* Section 6: Drug Pipeline */}
        <div id="pipeline" className="scroll-mt-6">
          <Section
            title="Drug Pipeline and Mechanisms"
            subtitle="Clinical development phases and mechanisms of action for approved and investigational drugs"
          >
            <DrugPipelineSection key={`pipe-${refreshKey}`} />
          </Section>
        </div>

        {/* Section 7: Chemical Space */}
        <div id="space" className="scroll-mt-6">
          <Section
            title="Chemical Space Explorer"
            subtitle="Interactive 2D projection of the chemical space — MW vs LogP colored by QED drug-likeness score"
          >
            <ChemicalSpaceSection key={`space-${refreshKey}`} />
          </Section>
        </div>

        {/* Section 8: Ingestion */}
        <div id="ingestion" className="scroll-mt-6">
          <Section
            title="Data Ingestion Status"
            subtitle="Real-time monitoring of the database population pipeline"
          >
            <IngestionStatusSection key={`ing-${refreshKey}`} />
          </Section>
        </div>

        {/* Footer */}
        <footer className="border-t border-slate-200 pt-8 pb-4">
          <div className="flex items-center justify-between text-xs text-slate-400">
            <span>PharmacoDB — BindX V9 · Data sourced from ChEMBL, UniProt, PubChem, BindingDB</span>
            <span>Auto-refresh every 60s</span>
          </div>
        </footer>

      </main>
    </div>
  )
}
