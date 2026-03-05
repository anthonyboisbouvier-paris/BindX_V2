import React, { useState, useRef, useEffect, useCallback } from 'react'
import BindXLogo from './BindXLogo.jsx'
import { queryAgent } from '../api.js'

// ---------------------------------------------------------------------------
// Chat message components
// ---------------------------------------------------------------------------

function UserMessage({ text }) {
  return (
    <div className="flex justify-end">
      <div className="max-w-[85%] bg-bx-surface text-white px-4 py-2.5 rounded-2xl rounded-br-sm text-sm leading-relaxed">
        {text}
      </div>
    </div>
  )
}

// ---------------------------------------------------------------------------
// Structured agent response renderer
// ---------------------------------------------------------------------------

const PRIORITY_COLORS = {
  high: 'bg-red-100 text-red-700',
  medium: 'bg-amber-100 text-amber-700',
  low: 'bg-blue-100 text-blue-700',
}

const QUALITY_COLORS = {
  good: 'bg-emerald-100 text-emerald-700',
  moderate: 'bg-amber-100 text-amber-700',
  poor: 'bg-red-100 text-red-700',
}

function ConfidenceBadge({ value }) {
  if (value == null) return null
  const pct = Math.round(value * 100)
  const color = pct >= 80 ? 'text-emerald-600' : pct >= 60 ? 'text-amber-600' : 'text-red-500'
  return (
    <div className="flex items-center gap-2">
      <div className="flex-1 h-1.5 bg-gray-200 rounded-full overflow-hidden">
        <div className={`h-full rounded-full ${pct >= 80 ? 'bg-emerald-400' : pct >= 60 ? 'bg-amber-400' : 'bg-red-400'}`}
          style={{ width: `${pct}%` }} />
      </div>
      <span className={`text-xs font-bold tabular-nums ${color}`}>{pct}%</span>
    </div>
  )
}

function SectionLabel({ children }) {
  return <h4 className="text-[10px] font-bold uppercase tracking-wider text-gray-400 mt-3 mb-1.5 first:mt-0">{children}</h4>
}

function StructuredAgentResponse({ data }) {
  if (!data || typeof data !== 'object') return <span>{String(data)}</span>

  // Metadata fields to exclude from display
  const META_KEYS = new Set(['agent_name', 'model', 'version', 'timestamp', 'input_hash', 'available'])

  // Normalize field names — handle both API naming conventions
  const summary = data.summary || (typeof data.analysis === 'string' ? data.analysis : null)
    || (typeof data.recommendation === 'string' ? data.recommendation : null)
  const runQuality = data.run_quality
  const hitRate = data.hit_rate_assessment
  const topCandidates = data.top_candidates
  const chemicalSeries = data.chemical_series
  const propertyAlerts = data.property_alerts
  const nextSteps = data.recommended_next_steps || data.next_steps
  const comparableCompounds = data.comparable_clinical_compounds || data.comparable_compounds
  const keyPapers = data.key_papers || data.literature_references
  const strengths = data.strengths
  const weaknesses = data.weaknesses
  const screeningStrategy = data.screening_strategy
  const riskMitigation = data.risk_mitigation
  const confidenceRationale = data.confidence_rationale

  // Parse confidence — handle "85%" string, 0.85 number, or 85 number
  let confidencePct = null
  if (data.confidence != null) {
    if (typeof data.confidence === 'number') {
      confidencePct = data.confidence <= 1 ? Math.round(data.confidence * 100) : Math.round(data.confidence)
    } else if (typeof data.confidence === 'string') {
      const match = data.confidence.match(/(\d+)/)
      if (match) confidencePct = parseInt(match[1])
    }
  }

  // Known keys to exclude from "rest" catch-all
  const KNOWN_KEYS = new Set([
    'summary', 'analysis', 'recommendation', 'confidence', 'confidence_rationale',
    'run_quality', 'hit_rate_assessment', 'top_candidates', 'chemical_series',
    'property_alerts', 'recommended_next_steps', 'next_steps',
    'comparable_clinical_compounds', 'comparable_compounds',
    'key_papers', 'literature_references', 'strengths', 'weaknesses',
    'screening_strategy', 'risk_mitigation',
    ...META_KEYS
  ])

  // Remaining unknown fields
  const restEntries = Object.entries(data)
    .filter(([k, v]) => !KNOWN_KEYS.has(k) && v != null && v !== '' && v !== true && v !== false
      && !(Array.isArray(v) && v.length === 0) && !(typeof v === 'object' && !Array.isArray(v) && Object.keys(v).length === 0))

  // Helper: render a string array as bullet list
  const BulletList = ({ items, colorClass = 'text-gray-600' }) => (
    <div className="space-y-1">
      {items.map((item, i) => (
        <div key={i} className="flex gap-1.5 items-start">
          <span className={`mt-1 w-1.5 h-1.5 rounded-full flex-shrink-0 ${colorClass === 'text-emerald-600' ? 'bg-emerald-400' : colorClass === 'text-amber-600' ? 'bg-amber-400' : 'bg-gray-400'}`} />
          <span className={`text-[11px] leading-snug ${colorClass}`}>{typeof item === 'string' ? item : JSON.stringify(item)}</span>
        </div>
      ))}
    </div>
  )

  return (
    <div className="space-y-0.5">
      {/* Summary */}
      {summary && <p className="text-[13px] text-gray-700 leading-relaxed">{summary}</p>}

      {/* Quality + Confidence row */}
      {(runQuality || confidencePct != null) && (
        <div className="flex items-center gap-2 flex-wrap mt-2">
          {runQuality && (
            <span className={`text-[10px] font-bold px-2 py-0.5 rounded-full ${QUALITY_COLORS[runQuality] || 'bg-gray-100 text-gray-600'}`}>
              Quality: {runQuality}
            </span>
          )}
          {confidencePct != null && (
            <div className="flex items-center gap-1.5 flex-1 min-w-[120px]">
              <span className="text-[10px] text-gray-400 font-medium">Confidence</span>
              <div className="flex-1"><ConfidenceBadge value={confidencePct / 100} /></div>
            </div>
          )}
        </div>
      )}
      {(confidenceRationale || (typeof data.confidence === 'string' && data.confidence.length > 5)) && (
        <p className="text-[11px] text-gray-400 italic mt-0.5">
          {confidenceRationale || (typeof data.confidence === 'string' && !data.confidence.match(/^\d+%?$/) ? data.confidence : null)}
        </p>
      )}

      {/* Hit rate */}
      {hitRate && (
        <>
          <SectionLabel>Hit Rate</SectionLabel>
          <p className="text-xs text-gray-600">{hitRate}</p>
        </>
      )}

      {/* Strengths */}
      {Array.isArray(strengths) && strengths.length > 0 && (
        <>
          <SectionLabel>Strengths</SectionLabel>
          <BulletList items={strengths} colorClass="text-emerald-600" />
        </>
      )}

      {/* Weaknesses */}
      {Array.isArray(weaknesses) && weaknesses.length > 0 && (
        <>
          <SectionLabel>Weaknesses</SectionLabel>
          <BulletList items={weaknesses} colorClass="text-amber-600" />
        </>
      )}

      {/* Top candidates */}
      {topCandidates?.length > 0 && (
        <>
          <SectionLabel>Top Candidates</SectionLabel>
          <div className="space-y-1">
            {topCandidates.map((c, i) => (
              <div key={i} className="flex gap-2 items-start">
                <span className="w-5 h-5 rounded-full bg-bx-mint/10 text-bx-mint text-[10px] font-bold flex items-center justify-center flex-shrink-0 mt-0.5">
                  {i + 1}
                </span>
                <div className="min-w-0">
                  <span className="text-xs font-semibold text-gray-700">{typeof c === 'string' ? c : c.name}</span>
                  {c.rationale && <p className="text-[11px] text-gray-500 leading-snug">{c.rationale}</p>}
                </div>
              </div>
            ))}
          </div>
        </>
      )}

      {/* Chemical series */}
      {chemicalSeries?.length > 0 && (
        <>
          <SectionLabel>Chemical Series</SectionLabel>
          <div className="space-y-1.5">
            {chemicalSeries.map((s, i) => (
              <div key={i} className="bg-white border border-gray-100 rounded-lg px-2.5 py-2">
                <div className="flex items-center justify-between">
                  <span className="text-xs font-semibold text-gray-700">{typeof s === 'string' ? s : s.name}</span>
                  {s.n_members && <span className="text-[10px] text-gray-400">{s.n_members} mol</span>}
                </div>
                {s.avg_affinity && <span className="text-[10px] text-gray-500">Avg: {s.avg_affinity} kcal/mol</span>}
                {s.strengths?.length > 0 && (
                  <div className="flex flex-wrap gap-1 mt-1">
                    {s.strengths.map((st, j) => (
                      <span key={j} className="text-[9px] px-1.5 py-0.5 rounded bg-emerald-50 text-emerald-600">{st}</span>
                    ))}
                  </div>
                )}
                {s.concerns?.length > 0 && (
                  <div className="flex flex-wrap gap-1 mt-1">
                    {s.concerns.map((c, j) => (
                      <span key={j} className="text-[9px] px-1.5 py-0.5 rounded bg-amber-50 text-amber-600">{c}</span>
                    ))}
                  </div>
                )}
              </div>
            ))}
          </div>
        </>
      )}

      {/* Property alerts */}
      {propertyAlerts?.length > 0 && (
        <>
          <SectionLabel>Alerts</SectionLabel>
          <BulletList items={propertyAlerts.map(a => typeof a === 'string' ? a : a.message || JSON.stringify(a))} colorClass="text-amber-600" />
        </>
      )}

      {/* Screening strategy */}
      {screeningStrategy && typeof screeningStrategy === 'object' && (
        <>
          <SectionLabel>Screening Strategy</SectionLabel>
          <div className="bg-white border border-gray-100 rounded-lg px-2.5 py-2 space-y-1">
            {screeningStrategy.recommended_mode && (
              <div className="flex items-center gap-2">
                <span className="text-[10px] text-gray-400">Mode:</span>
                <span className="text-[10px] font-semibold text-gray-700">{screeningStrategy.recommended_mode}</span>
              </div>
            )}
            {screeningStrategy.expected_hit_rate && (
              <div className="flex items-center gap-2">
                <span className="text-[10px] text-gray-400">Expected hit rate:</span>
                <span className="text-[10px] font-semibold text-gray-700">{screeningStrategy.expected_hit_rate}</span>
              </div>
            )}
            {screeningStrategy.library_focus && (
              <p className="text-[11px] text-gray-600 leading-snug">{screeningStrategy.library_focus}</p>
            )}
            {screeningStrategy.key_considerations?.length > 0 && (
              <BulletList items={screeningStrategy.key_considerations} />
            )}
          </div>
        </>
      )}

      {/* Risk mitigation */}
      {Array.isArray(riskMitigation) && riskMitigation.length > 0 && (
        <>
          <SectionLabel>Risk Mitigation</SectionLabel>
          <BulletList items={riskMitigation} />
        </>
      )}

      {/* Next steps — handle both object array [{action, priority}] and string array */}
      {Array.isArray(nextSteps) && nextSteps.length > 0 && (
        <>
          <SectionLabel>Next Steps</SectionLabel>
          <div className="space-y-1.5">
            {nextSteps.map((s, i) => {
              if (typeof s === 'string') {
                return (
                  <div key={i} className="flex gap-2 items-start">
                    <span className="w-5 h-5 rounded-full bg-blue-50 text-blue-600 text-[10px] font-bold flex items-center justify-center flex-shrink-0 mt-0.5">
                      {i + 1}
                    </span>
                    <span className="text-[11px] text-gray-600 leading-snug">{s}</span>
                  </div>
                )
              }
              return (
                <div key={i} className="bg-white border border-gray-100 rounded-lg px-2.5 py-2">
                  <div className="flex items-center gap-1.5 mb-0.5">
                    {s.priority && (
                      <span className={`text-[9px] font-bold px-1.5 py-0.5 rounded-full ${PRIORITY_COLORS[s.priority] || 'bg-gray-100 text-gray-500'}`}>
                        {s.priority}
                      </span>
                    )}
                    <span className="text-xs font-medium text-gray-700">{s.action || s.step || JSON.stringify(s)}</span>
                  </div>
                  {s.rationale && <p className="text-[11px] text-gray-400 leading-snug">{s.rationale}</p>}
                </div>
              )
            })}
          </div>
        </>
      )}

      {/* Comparable compounds */}
      {Array.isArray(comparableCompounds) && comparableCompounds.length > 0 && (
        <>
          <SectionLabel>Comparable Compounds</SectionLabel>
          <div className="flex flex-wrap gap-1">
            {comparableCompounds.map((c, i) => (
              <span key={i} className="text-[10px] px-2 py-0.5 rounded-full bg-blue-50 text-blue-600 font-medium">
                {typeof c === 'string' ? c : c.name || JSON.stringify(c)}
              </span>
            ))}
          </div>
        </>
      )}

      {/* Papers / references */}
      {Array.isArray(keyPapers) && keyPapers.length > 0 && (
        <>
          <SectionLabel>Key Papers</SectionLabel>
          {keyPapers.map((p, i) => (
            <div key={i} className="text-[11px] text-gray-500">
              {typeof p === 'string' ? p : (
                <>
                  {p.title}{p.pmid && <span className="text-gray-400 ml-1">(PMID: {p.pmid})</span>}
                  {p.finding && <p className="text-[10px] text-gray-400 italic">{p.finding}</p>}
                </>
              )}
            </div>
          ))}
        </>
      )}

      {/* Nested analysis object */}
      {typeof data.analysis === 'object' && data.analysis !== null && (
        <StructuredAgentResponse data={data.analysis} />
      )}

      {/* Unknown fields (excluding metadata) */}
      {restEntries.length > 0 && restEntries.map(([key, value]) => (
        <div key={key}>
          <SectionLabel>{key.replace(/_/g, ' ').replace(/\b\w/g, c => c.toUpperCase())}</SectionLabel>
          {Array.isArray(value) ? (
            <BulletList items={value.map(v => typeof v === 'string' ? v : JSON.stringify(v))} />
          ) : typeof value === 'object' ? (
            <div className="bg-white border border-gray-100 rounded-lg px-2.5 py-2 text-[11px] text-gray-600">
              {Object.entries(value).map(([k, v]) => (
                <div key={k} className="flex gap-2">
                  <span className="text-gray-400">{k.replace(/_/g, ' ')}:</span>
                  <span className="font-medium">{typeof v === 'string' ? v : JSON.stringify(v)}</span>
                </div>
              ))}
            </div>
          ) : (
            <p className="text-xs text-gray-600">{String(value)}</p>
          )}
        </div>
      ))}
    </div>
  )
}

function AgentMessage({ text, agent, timestamp }) {
  const isStructured = text && typeof text === 'object'

  return (
    <div className="flex gap-2.5">
      <div className="w-7 h-7 rounded-full bg-gradient-to-br from-bx-mint to-emerald-400 flex items-center justify-center flex-shrink-0 mt-0.5">
        <svg className="w-3.5 h-3.5 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
            d="M9.75 17L9 20l-1 1h8l-1-1-.75-3M3 13h18M5 17H3a2 2 0 01-2-2V5a2 2 0 012-2h14a2 2 0 012 2v10a2 2 0 01-2 2h-2" />
        </svg>
      </div>
      <div className="max-w-[85%] flex-1 min-w-0">
        <div className="flex items-center gap-2 mb-1">
          <span className="text-[10px] font-semibold text-bx-mint uppercase tracking-wider">{agent || 'BindX Agent'}</span>
          {timestamp && <span className="text-[10px] text-gray-400">{new Date(timestamp).toLocaleTimeString()}</span>}
        </div>
        <div className="bg-gray-50 border border-gray-100 px-4 py-3 rounded-2xl rounded-tl-sm">
          {isStructured ? <StructuredAgentResponse data={text} /> : (
            <p className="text-sm text-gray-700 leading-relaxed whitespace-pre-wrap">{text}</p>
          )}
        </div>
      </div>
    </div>
  )
}

function TypingIndicator() {
  return (
    <div className="flex gap-2.5">
      <div className="w-7 h-7 rounded-full bg-gradient-to-br from-bx-mint to-emerald-400 flex items-center justify-center flex-shrink-0">
        <BindXLogo variant="loading" size={20} />
      </div>
      <div className="bg-gray-50 border border-gray-100 px-4 py-3 rounded-2xl rounded-tl-sm flex items-center gap-1.5">
        <span className="w-1.5 h-1.5 bg-gray-400 rounded-full animate-bounce" style={{ animationDelay: '0ms' }} />
        <span className="w-1.5 h-1.5 bg-gray-400 rounded-full animate-bounce" style={{ animationDelay: '150ms' }} />
        <span className="w-1.5 h-1.5 bg-gray-400 rounded-full animate-bounce" style={{ animationDelay: '300ms' }} />
      </div>
    </div>
  )
}

// ---------------------------------------------------------------------------
// Agent type picker
// ---------------------------------------------------------------------------
const AGENTS = [
  { key: 'run_analysis', label: 'Run Analysis', desc: 'Hit selection, clustering, scoring trends' },
  { key: 'candidate', label: 'Molecule Review', desc: 'ADMET, safety, synthesis assessment' },
  { key: 'target', label: 'Target Strategy', desc: 'Druggability, screening recommendations' },
  { key: 'optimization', label: 'Lead Optimization', desc: 'SAR, modifications, multi-parameter' },
]

// ---------------------------------------------------------------------------
// Build campaign context from workspace data
// ---------------------------------------------------------------------------
function quartiles(arr) {
  if (!arr.length) return null
  const s = [...arr].sort((a, b) => a - b)
  const q = (p) => { const i = p * (s.length - 1); const lo = Math.floor(i); return lo === i ? s[lo] : s[lo] + (s[lo + 1] - s[lo]) * (i - lo) }
  return { min: s[0], q25: +q(0.25).toFixed(3), median: +q(0.5).toFixed(3), q75: +q(0.75).toFixed(3), max: s[s.length - 1] }
}

function buildContext(project, campaign, phase, molecules) {
  const ctx = {
    project_name: project?.name || 'Unknown',
    target: project?.target_input_value || '',
  }

  if (campaign) {
    ctx.campaign_name = campaign.name
    ctx.pocket = campaign.pocket_label || ''
  }

  if (phase) {
    ctx.phase = phase.label
    ctx.phase_type = phase.type
    ctx.molecule_count = molecules?.length || 0
  }

  if (molecules?.length > 0) {
    const scores = molecules.filter(m => m.docking_score != null).map(m => m.docking_score)
    const composites = molecules.filter(m => m.composite_score != null).map(m => m.composite_score)
    const qeds = molecules.filter(m => m.QED != null).map(m => m.QED)
    const logPs = molecules.filter(m => m.logP != null).map(m => m.logP)
    const mws = molecules.filter(m => m.MW != null).map(m => m.MW)
    const bookmarked = molecules.filter(m => m.bookmarked)

    ctx.molecules_summary = {
      total: molecules.length,
      bookmarked: bookmarked.length,
      with_docking: scores.length,
      with_composite: composites.length,
      avg_docking: scores.length ? +(scores.reduce((a, b) => a + b, 0) / scores.length).toFixed(2) : null,
      avg_composite: composites.length ? +(composites.reduce((a, b) => a + b, 0) / composites.length).toFixed(2) : null,
      best_docking: scores.length ? +Math.min(...scores).toFixed(2) : null,
      best_composite: composites.length ? +Math.max(...composites).toFixed(2) : null,
    }

    // Score distributions (quartiles)
    ctx.distributions = {}
    if (scores.length) ctx.distributions.docking_score = quartiles(scores)
    if (composites.length) ctx.distributions.composite_score = quartiles(composites)
    if (qeds.length) ctx.distributions.QED = quartiles(qeds)
    if (logPs.length) ctx.distributions.logP = quartiles(logPs)
    if (mws.length) ctx.distributions.MW = quartiles(mws)

    // Scaffold diversity
    const scaffolds = new Set(molecules.map(m => m.scaffold).filter(Boolean))
    const clusters = new Set(molecules.map(m => m.cluster_id).filter(v => v != null))
    ctx.chemical_diversity = {
      unique_scaffolds: scaffolds.size,
      n_clusters: clusters.size,
      ai_generated: molecules.filter(m => m.ai_generated).length,
      max_generation_level: Math.max(0, ...molecules.map(m => m.generation_level || 0)),
    }

    // Safety flags summary
    const safetyRed = molecules.filter(m => m.safety_color_code === 'red').length
    const safetyYellow = molecules.filter(m => m.safety_color_code === 'yellow').length
    const safetyGreen = molecules.filter(m => m.safety_color_code === 'green').length
    const painsAlerts = molecules.filter(m => m.pains_alert).length
    const lipinski = molecules.filter(m => m.lipinski_pass).length
    ctx.safety_overview = {
      red: safetyRed,
      yellow: safetyYellow,
      green: safetyGreen,
      pains_alerts: painsAlerts,
      lipinski_pass: lipinski,
      lipinski_total: molecules.filter(m => m.lipinski_pass != null).length,
    }

    // Off-target summary
    const withSelectivity = molecules.filter(m => m.selectivity_score != null)
    if (withSelectivity.length) {
      ctx.selectivity_overview = {
        count: withSelectivity.length,
        avg_selectivity: +(withSelectivity.reduce((s, m) => s + m.selectivity_score, 0) / withSelectivity.length).toFixed(3),
        avg_off_target_hits: +(withSelectivity.reduce((s, m) => s + (m.off_target_hits || 0), 0) / withSelectivity.length).toFixed(1),
      }
    }

    // Top 10 molecules by composite score (more context)
    ctx.top_molecules = [...molecules]
      .filter(m => m.composite_score != null)
      .sort((a, b) => (b.composite_score || 0) - (a.composite_score || 0))
      .slice(0, 10)
      .map(m => ({
        name: m.name || m.id,
        smiles: m.smiles,
        docking_score: m.docking_score,
        composite_score: m.composite_score,
        QED: m.QED,
        logP: m.logP,
        MW: m.MW,
        selectivity_score: m.selectivity_score,
        safety_color_code: m.safety_color_code,
        bookmarked: m.bookmarked,
        generation_level: m.generation_level || 0,
      }))
  }

  return ctx
}

// ---------------------------------------------------------------------------
// CampaignAgentPanel — Chat interface for AI agent
// ---------------------------------------------------------------------------
// Props:
//   isOpen     — boolean
//   onClose    — () => void
//   project    — current project object
//   campaign   — current campaign object
//   phase      — current phase object
//   molecules  — flat molecule array

export default function CampaignAgentPanel({ isOpen, onClose, project, campaign, phase, molecules }) {
  const [messages, setMessages] = useState([])
  const [input, setInput] = useState('')
  const [loading, setLoading] = useState(false)
  const [selectedAgent, setSelectedAgent] = useState('run_analysis')
  const scrollRef = useRef(null)
  const inputRef = useRef(null)

  // Auto-scroll on new messages
  useEffect(() => {
    if (scrollRef.current) {
      scrollRef.current.scrollTop = scrollRef.current.scrollHeight
    }
  }, [messages, loading])

  // Focus input when panel opens
  useEffect(() => {
    if (isOpen && inputRef.current) {
      setTimeout(() => inputRef.current?.focus(), 100)
    }
  }, [isOpen])

  const handleSend = useCallback(async () => {
    const q = input.trim()
    if (!q || loading) return

    const userMsg = { role: 'user', text: q, timestamp: new Date().toISOString() }
    setMessages(prev => [...prev, userMsg])
    setInput('')
    setLoading(true)

    try {
      const context = buildContext(project, campaign, phase, molecules)
      context.user_question = q

      const result = await queryAgent(selectedAgent, context, project?.id)

      // Handle unavailable agent (missing API key, etc.)
      if (result.available === false) {
        const agentMsg = {
          role: 'agent',
          agent: 'System',
          text: `The AI Agent is not available yet. ${result.fallback || 'The backend API key may not be configured.'}. This feature requires an OpenAI API key to be set on the server.`,
          timestamp: new Date().toISOString(),
          isError: true,
        }
        setMessages(prev => [...prev, agentMsg])
      } else {
        // Extract text — prefer string fields, fall back to formatting the whole object
        const rawText = (typeof result.analysis === 'string' && result.analysis)
          || (typeof result.recommendation === 'string' && result.recommendation)
          || (typeof result.summary === 'string' && result.summary)
          || null
        const agentMsg = {
          role: 'agent',
          agent: AGENTS.find(a => a.key === selectedAgent)?.label || selectedAgent,
          text: rawText || result,
          timestamp: new Date().toISOString(),
        }
        setMessages(prev => [...prev, agentMsg])
      }
    } catch (err) {
      const errorMsg = {
        role: 'agent',
        agent: 'System',
        text: `Agent unavailable: ${err.userMessage || err.message || 'Connection error'}. The agent backend may not be running or the API key may not be configured.`,
        timestamp: new Date().toISOString(),
        isError: true,
      }
      setMessages(prev => [...prev, errorMsg])
    } finally {
      setLoading(false)
    }
  }, [input, loading, selectedAgent, project, campaign, phase, molecules])

  if (!isOpen) return null

  return (
    <div className="fixed inset-y-0 right-0 z-40 flex" onClick={onClose}>
      {/* Backdrop */}
      <div className="absolute inset-0 bg-black/20" />

      {/* Panel */}
      <div
        className="relative ml-auto w-[420px] max-w-[90vw] bg-white shadow-2xl flex flex-col border-l border-gray-200"
        onClick={e => e.stopPropagation()}
      >
        {/* Header */}
        <div className="flex items-center justify-between px-5 py-3.5 border-b border-gray-100 bg-gradient-to-r from-gray-50 to-white flex-shrink-0">
          <div className="flex items-center gap-2.5">
            <div className="w-8 h-8 rounded-xl bg-gradient-to-br from-bx-mint to-emerald-400 flex items-center justify-center">
              <svg className="w-4 h-4 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                  d="M9.75 17L9 20l-1 1h8l-1-1-.75-3M3 13h18M5 17H3a2 2 0 01-2-2V5a2 2 0 012-2h14a2 2 0 012 2v10a2 2 0 01-2 2h-2" />
              </svg>
            </div>
            <div>
              <h3 className="text-sm font-bold text-gray-800">BindX AI Agent</h3>
              <p className="text-[10px] text-gray-400">{project?.name || 'Campaign'} analysis</p>
            </div>
          </div>
          <button
            onClick={onClose}
            className="p-1.5 rounded-lg hover:bg-gray-100 text-gray-400 hover:text-gray-700 transition-colors"
          >
            <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
            </svg>
          </button>
        </div>

        {/* Agent selector */}
        <div className="px-5 py-2.5 border-b border-gray-100 bg-gray-50/50">
          <div className="flex gap-1.5 overflow-x-auto" style={{ scrollbarWidth: 'none' }}>
            {AGENTS.map(agent => (
              <button
                key={agent.key}
                onClick={() => setSelectedAgent(agent.key)}
                className={`flex-shrink-0 px-2.5 py-1.5 rounded-lg text-[11px] font-semibold transition-all ${
                  selectedAgent === agent.key
                    ? 'bg-bx-surface text-white shadow-sm'
                    : 'bg-white border border-gray-200 text-gray-600 hover:border-gray-300 hover:bg-gray-50'
                }`}
                title={agent.desc}
              >
                {agent.label}
              </button>
            ))}
          </div>
        </div>

        {/* Context summary */}
        {molecules?.length > 0 && (
          <div className="px-5 py-2 border-b border-gray-100 bg-blue-50/30">
            <p className="text-[10px] text-gray-500">
              Context: <span className="font-semibold text-gray-700">{molecules.length}</span> molecules
              {phase && <> in <span className="font-semibold text-gray-700">{phase.label}</span></>}
              {molecules.filter(m => m.bookmarked).length > 0 && (
                <> ({molecules.filter(m => m.bookmarked).length} bookmarked)</>
              )}
            </p>
          </div>
        )}

        {/* Messages area */}
        <div ref={scrollRef} className="flex-1 overflow-y-auto px-5 py-4 space-y-4"
          style={{ scrollbarWidth: 'thin', scrollbarColor: '#e5e7eb transparent' }}>

          {messages.length === 0 && (
            <div className="flex flex-col items-center justify-center py-12 text-center">
              <div className="w-14 h-14 rounded-2xl bg-gradient-to-br from-bx-mint/10 to-emerald-50 flex items-center justify-center mb-4">
                <svg className="w-7 h-7 text-bx-mint" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                    d="M8 10h.01M12 10h.01M16 10h.01M9 16H5a2 2 0 01-2-2V6a2 2 0 012-2h14a2 2 0 012 2v8a2 2 0 01-2 2h-5l-5 5v-5z" />
                </svg>
              </div>
              <p className="text-sm font-semibold text-gray-700 mb-1">Ask the AI Agent</p>
              <p className="text-xs text-gray-400 max-w-[240px] leading-relaxed">
                Analyze your campaign data, get recommendations for next runs, identify scaffolds, or assess molecule safety.
              </p>
              <div className="mt-4 space-y-1.5 w-full max-w-[280px]">
                {[
                  'What are the best hits so far?',
                  'Which molecules should I prioritize?',
                  'Recommend the next analysis to run',
                  'Are there any safety red flags?',
                ].map(suggestion => (
                  <button
                    key={suggestion}
                    onClick={() => { setInput(suggestion); inputRef.current?.focus() }}
                    className="w-full text-left px-3 py-2 rounded-lg border border-gray-200 text-xs text-gray-600 hover:border-bx-mint/40 hover:bg-emerald-50/30 hover:text-gray-800 transition-all"
                  >
                    {suggestion}
                  </button>
                ))}
              </div>
            </div>
          )}

          {messages.map((msg, i) => (
            msg.role === 'user'
              ? <UserMessage key={i} text={msg.text} />
              : <AgentMessage key={i} text={msg.text} agent={msg.agent} timestamp={msg.timestamp} />
          ))}

          {loading && <TypingIndicator />}
        </div>

        {/* Input area */}
        <div className="px-5 py-3.5 border-t border-gray-100 bg-white flex-shrink-0">
          <div className="flex items-center gap-2">
            <input
              ref={inputRef}
              type="text"
              value={input}
              onChange={e => setInput(e.target.value)}
              onKeyDown={e => e.key === 'Enter' && !e.shiftKey && handleSend()}
              placeholder={`Ask ${AGENTS.find(a => a.key === selectedAgent)?.label || 'agent'}...`}
              className="flex-1 text-sm px-3.5 py-2.5 rounded-xl border border-gray-200 bg-gray-50
                         focus:outline-none focus:ring-2 focus:ring-bx-mint/20 focus:border-bx-mint
                         transition-colors placeholder-gray-400"
              disabled={loading}
            />
            <button
              onClick={handleSend}
              disabled={loading || !input.trim()}
              className="p-2.5 bg-bx-surface text-white rounded-xl hover:bg-bx-elevated
                         disabled:opacity-40 disabled:cursor-not-allowed transition-colors flex-shrink-0"
            >
              {loading ? (
                <BindXLogo variant="loading" size={20} />
              ) : (
                <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
                    d="M12 19l9 2-9-18-9 18 9-2zm0 0v-8" />
                </svg>
              )}
            </button>
          </div>
        </div>
      </div>
    </div>
  )
}
