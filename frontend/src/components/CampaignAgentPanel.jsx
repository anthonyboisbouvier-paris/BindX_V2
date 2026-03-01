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

function AgentMessage({ text, agent, timestamp }) {
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
        <div className="bg-gray-50 border border-gray-100 px-4 py-2.5 rounded-2xl rounded-tl-sm text-sm text-gray-700 leading-relaxed whitespace-pre-wrap">
          {text}
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
    // Summary stats
    const scores = molecules.filter(m => m.docking_score != null).map(m => m.docking_score)
    const composites = molecules.filter(m => m.composite_score != null).map(m => m.composite_score)
    const bookmarked = molecules.filter(m => m.bookmarked)

    ctx.molecules_summary = {
      total: molecules.length,
      bookmarked: bookmarked.length,
      with_docking: scores.length,
      with_composite: composites.length,
      avg_docking: scores.length ? (scores.reduce((a, b) => a + b, 0) / scores.length).toFixed(2) : null,
      avg_composite: composites.length ? (composites.reduce((a, b) => a + b, 0) / composites.length).toFixed(2) : null,
      best_docking: scores.length ? Math.min(...scores).toFixed(2) : null,
      best_composite: composites.length ? Math.max(...composites).toFixed(2) : null,
    }

    // Top 5 molecules by composite score
    ctx.top_molecules = [...molecules]
      .filter(m => m.composite_score != null)
      .sort((a, b) => (b.composite_score || 0) - (a.composite_score || 0))
      .slice(0, 5)
      .map(m => ({
        name: m.name || m.id,
        smiles: m.smiles,
        docking_score: m.docking_score,
        composite_score: m.composite_score,
        QED: m.QED,
        logP: m.logP,
        bookmarked: m.bookmarked,
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
        const agentMsg = {
          role: 'agent',
          agent: AGENTS.find(a => a.key === selectedAgent)?.label || selectedAgent,
          text: result.analysis || result.recommendation || result.summary || JSON.stringify(result, null, 2),
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
