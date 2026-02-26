import React, { useState, useRef, useEffect } from 'react'
import Badge from './Badge.jsx'

// ---------------------------------------------------------------------------
// Relative time
// ---------------------------------------------------------------------------

function relativeTime(iso) {
  if (!iso) return ''
  const diff = Date.now() - new Date(iso).getTime()
  const mins = Math.floor(diff / 60000)
  if (mins < 1) return 'just now'
  if (mins < 60) return `${mins}m ago`
  const hrs = Math.floor(mins / 60)
  if (hrs < 24) return `${hrs}h ago`
  return `${Math.floor(hrs / 24)}d ago`
}

// ---------------------------------------------------------------------------
// Priority badge
// ---------------------------------------------------------------------------

function PriorityBadge({ priority }) {
  const map = { high: 'red', medium: 'yellow', low: 'green' }
  return (
    <Badge variant={map[priority] || 'gray'} size="sm">
      {priority?.toUpperCase()}
    </Badge>
  )
}

// ---------------------------------------------------------------------------
// Type icon
// ---------------------------------------------------------------------------

function TypeIcon({ type }) {
  if (type === 'recommendation') {
    return (
      <svg className="w-4 h-4 text-amber-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
          d="M9.663 17h4.673M12 3v1m6.364 1.636l-.707.707M21 12h-1M4 12H3m3.343-5.657l-.707-.707m2.828 9.9a5 5 0 117.072 0l-.548.547A3.374 3.374 0 0014 18.469V19a2 2 0 11-4 0v-.531c0-.895-.356-1.754-.988-2.386l-.548-.547z" />
      </svg>
    )
  }
  if (type === 'alert') {
    return (
      <svg className="w-4 h-4 text-red-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
          d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
      </svg>
    )
  }
  // observation
  return (
    <svg className="w-4 h-4 text-blue-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
        d="M15 12a3 3 0 11-6 0 3 3 0 016 0z" />
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
        d="M2.458 12C3.732 7.943 7.523 5 12 5c4.478 0 8.268 2.943 9.542 7-1.274 4.057-5.064 7-9.542 7-4.477 0-8.268-2.943-9.542-7z" />
    </svg>
  )
}

// ---------------------------------------------------------------------------
// Single insight card
// ---------------------------------------------------------------------------

function InsightCard({ insight, onAction }) {
  const [expanded, setExpanded] = useState(insight.priority === 'high')

  const borderColor = {
    high: 'border-l-red-400',
    medium: 'border-l-amber-400',
    low: 'border-l-green-400',
  }[insight.priority] || 'border-l-gray-300'

  return (
    <div className={`border-l-2 ${borderColor} pl-3 py-2`}>
      <div
        className="flex items-start gap-2 cursor-pointer select-none"
        onClick={() => setExpanded(e => !e)}
      >
        <TypeIcon type={insight.type} />
        <div className="flex-1 min-w-0">
          <div className="flex items-center gap-2 flex-wrap mb-0.5">
            <span className="text-sm font-semibold text-gray-800 leading-tight">{insight.title}</span>
          </div>
          <div className="flex items-center gap-2">
            <PriorityBadge priority={insight.priority} />
            <span className="text-xs text-gray-400">{relativeTime(insight.timestamp)}</span>
          </div>
        </div>
        <svg
          className={`w-4 h-4 text-gray-300 shrink-0 transition-transform mt-0.5 ${expanded ? 'rotate-180' : ''}`}
          fill="none" stroke="currentColor" viewBox="0 0 24 24"
        >
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
        </svg>
      </div>

      {expanded && insight.body && (
        <div className="mt-2 space-y-2">
          <p className="text-xs text-gray-500 leading-relaxed">{insight.body}</p>
          {insight.action && (
            <button
              className="inline-flex items-center gap-1.5 px-3 py-1.5 bg-[#1e3a5f] text-white text-xs font-medium rounded-lg hover:bg-[#1e4a7f] transition-colors"
              onClick={(e) => {
                e.stopPropagation()
                if (onAction) onAction(insight.action)
              }}
            >
              <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M14.752 11.168l-3.197-2.132A1 1 0 0010 9.87v4.263a1 1 0 001.555.832l3.197-2.132a1 1 0 000-1.664z" />
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
              </svg>
              {insight.action.label}
            </button>
          )}
        </div>
      )}
    </div>
  )
}

// ---------------------------------------------------------------------------
// CampaignAgentPanel
// ---------------------------------------------------------------------------
// Props:
//   insights   — array of insight objects (from MOCK_AGENT_INSIGHTS)
//   onAskAgent — (question: string) => void

export default function CampaignAgentPanel({ insights: initialInsights = [], onAskAgent }) {
  const [insights, setInsights] = useState(initialInsights)
  const [question, setQuestion] = useState('')
  const [loading, setLoading] = useState(false)
  const inputRef = useRef(null)

  // Sync if prop changes
  useEffect(() => {
    setInsights(initialInsights)
  }, [initialInsights])

  const handleAsk = async () => {
    const q = question.trim()
    if (!q) return

    setLoading(true)
    const questionInsight = {
      id: `user-q-${Date.now()}`,
      type: 'observation',
      priority: 'medium',
      timestamp: new Date().toISOString(),
      title: q,
      body: 'Analyzing your question against campaign data... Based on current Phase A and B data, the key observation is that Osimertinib-v1 stands out with the highest composite score (91.3) and CNN affinity (8.0). Recommend prioritizing for full ADMET run.',
    }

    // Simulate async response
    setTimeout(() => {
      setInsights(prev => [questionInsight, ...prev])
      setLoading(false)
      setQuestion('')
    }, 800)

    if (onAskAgent) onAskAgent(q)
  }

  const handleAction = (action) => {
    console.log('[CampaignAgentPanel] Action triggered:', action)
  }

  return (
    <div className="bg-white rounded-xl border border-gray-100 shadow-sm overflow-hidden">
      {/* Header */}
      <div className="px-5 py-3.5 border-b border-gray-50 flex items-center gap-2">
        <div className="w-6 h-6 rounded-full bg-[#1e3a5f] flex items-center justify-center shrink-0">
          <svg className="w-3.5 h-3.5 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
              d="M9.75 17L9 20l-1 1h8l-1-1-.75-3M3 13h18M5 17H3a2 2 0 01-2-2V5a2 2 0 012-2h14a2 2 0 012 2v10a2 2 0 01-2 2h-2" />
          </svg>
        </div>
        <span className="text-sm font-semibold text-gray-700">AI Agent Insights</span>
        <span className="ml-auto text-xs text-gray-400 bg-gray-50 px-2 py-0.5 rounded-full">
          {insights.length} insight{insights.length !== 1 ? 's' : ''}
        </span>
      </div>

      {/* Insight list */}
      <div className="px-5 py-3 space-y-3 max-h-80 overflow-y-auto">
        {insights.length === 0 && (
          <p className="text-sm text-gray-400 text-center py-4">No insights yet.</p>
        )}
        {insights.map(ins => (
          <InsightCard key={ins.id} insight={ins} onAction={handleAction} />
        ))}
      </div>

      {/* Agent chat input */}
      <div className="px-5 py-3.5 border-t border-gray-50">
        <div className="flex items-center gap-2">
          <input
            ref={inputRef}
            type="text"
            value={question}
            onChange={e => setQuestion(e.target.value)}
            onKeyDown={e => e.key === 'Enter' && handleAsk()}
            placeholder="Ask the AI agent about this campaign..."
            className="flex-1 text-sm px-3 py-2 rounded-lg border border-gray-200 bg-gray-50 focus:outline-none focus:ring-2 focus:ring-[#1e3a5f]/20 focus:border-[#1e3a5f] transition-colors placeholder-gray-400"
            disabled={loading}
          />
          <button
            onClick={handleAsk}
            disabled={loading || !question.trim()}
            className="px-3 py-2 bg-[#1e3a5f] text-white text-xs font-semibold rounded-lg hover:bg-[#1e4a7f] disabled:opacity-40 transition-colors flex items-center gap-1.5 shrink-0"
          >
            {loading ? (
              <svg className="w-3.5 h-3.5 animate-spin" fill="none" viewBox="0 0 24 24">
                <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
                <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
              </svg>
            ) : (
              <>
                Ask
                <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M14 5l7 7m0 0l-7 7m7-7H3" />
                </svg>
              </>
            )}
          </button>
        </div>
      </div>
    </div>
  )
}
