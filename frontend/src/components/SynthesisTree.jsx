import React, { useState, useRef, useLayoutEffect, useEffect, useCallback } from 'react'
import { createPortal } from 'react-dom'
import InfoTip from './InfoTip.jsx'

// --------------------------------------------------
// smiles-drawer lazy-load (same pattern as MoleculeTable)
// --------------------------------------------------
let _sdPromise = null
function getSmilesDrawer() {
  if (!_sdPromise) {
    _sdPromise = import('smiles-drawer').then(mod => mod.default || mod)
  }
  return _sdPromise
}
const _treeSvgCache = new Map()

// --------------------------------------------------
// Color helpers
// --------------------------------------------------
function confidenceColor(conf) {
  if (conf === null || conf === undefined) return { bg: 'bg-gray-100', text: 'text-gray-500', border: 'border-gray-200' }
  if (conf >= 0.75) return { bg: 'bg-green-100', text: 'text-green-700', border: 'border-green-200' }
  if (conf >= 0.5)  return { bg: 'bg-yellow-100', text: 'text-yellow-700', border: 'border-yellow-200' }
  return { bg: 'bg-red-100', text: 'text-red-600', border: 'border-red-200' }
}

function costColor(cost) {
  if (!cost) return 'text-gray-400'
  const c = cost.toLowerCase()
  if (c === 'low')    return 'text-green-600'
  if (c === 'medium') return 'text-yellow-600'
  return 'text-red-600'
}

// --------------------------------------------------
// Summary bar at top of the card
// --------------------------------------------------
function SummaryBar({ nSteps, confidence, estimatedCost, allReagentsAvailable }) {
  const confColor = confidenceColor(confidence)
  return (
    <div className="flex flex-wrap items-center gap-3 px-4 py-3 bg-gray-50 border-b border-gray-100">
      {nSteps !== undefined && nSteps !== null && (
        <span className="flex items-center gap-1.5 text-sm font-medium text-gray-700">
          <svg className="w-4 h-4 text-purple-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 6h16M4 12h16M4 18h7" />
          </svg>
          <span>{nSteps} step{nSteps !== 1 ? 's' : ''}</span>
          <InfoTip text="Number of chemical reactions required to synthesize this molecule from starting materials." />
        </span>
      )}
      {confidence !== undefined && confidence !== null && (
        <span className={`flex items-center gap-1 px-2 py-0.5 rounded-full text-sm font-semibold ${confColor.bg} ${confColor.text} border ${confColor.border}`}>
          {Math.round(confidence * 100)}% confidence
          <InfoTip text="Confidence that this synthesis route is feasible based on known chemical reactions." />
        </span>
      )}
      {estimatedCost && (
        <span className={`flex items-center gap-1 text-sm font-medium ${costColor(estimatedCost)}`}>
          <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
              d="M12 8c-1.657 0-3 .895-3 2s1.343 2 3 2 3 .895 3 2-1.343 2-3 2m0-8c1.11 0 2.08.402 2.599 1M12 8V7m0 1v8m0 0v1m0-1c-1.11 0-2.08-.402-2.599-1M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
          </svg>
          Cost: {estimatedCost}
        </span>
      )}
      {allReagentsAvailable != null && (
        <span className={`flex items-center gap-1 text-sm font-medium ${allReagentsAvailable ? 'text-green-600' : 'text-orange-500'}`}>
          {allReagentsAvailable ? (
            <svg className="w-3.5 h-3.5" fill="currentColor" viewBox="0 0 20 20">
              <path fillRule="evenodd" d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z" clipRule="evenodd" />
            </svg>
          ) : (
            <svg className="w-3.5 h-3.5" fill="currentColor" viewBox="0 0 20 20">
              <path fillRule="evenodd" d="M8.257 3.099c.765-1.36 2.722-1.36 3.486 0l5.58 9.92c.75 1.334-.213 2.98-1.742 2.98H4.42c-1.53 0-2.493-1.646-1.743-2.98l5.58-9.92zM11 13a1 1 0 11-2 0 1 1 0 012 0zm-1-8a1 1 0 00-1 1v3a1 1 0 002 0V6a1 1 0 00-1-1z" clipRule="evenodd" />
            </svg>
          )}
          {allReagentsAvailable ? 'All reagents available' : 'Some reagents missing'}
        </span>
      )}
    </div>
  )
}

// --------------------------------------------------
// NodeSmilesImage — compact 2D structure for tree nodes
// --------------------------------------------------
const IMG_W = 100, IMG_H = 65

function NodeSmilesImage({ smiles }) {
  const [svgHtml, setSvgHtml] = useState(() => _treeSvgCache.get(smiles) || null)
  const [failed, setFailed] = useState(false)

  useEffect(() => {
    if (!smiles) return
    if (_treeSvgCache.has(smiles)) { setSvgHtml(_treeSvgCache.get(smiles)); return }
    let cancelled = false
    getSmilesDrawer().then(SD => {
      if (cancelled) return
      try {
        SD.parse(smiles, (tree) => {
          if (cancelled) return
          const svgEl = document.createElementNS('http://www.w3.org/2000/svg', 'svg')
          svgEl.setAttribute('width', String(IMG_W))
          svgEl.setAttribute('height', String(IMG_H))
          document.body.appendChild(svgEl)
          const drawer = new SD.SvgDrawer({ width: IMG_W, height: IMG_H })
          drawer.draw(tree, svgEl, 'light')
          const html = svgEl.outerHTML
          document.body.removeChild(svgEl)
          _treeSvgCache.set(smiles, html)
          if (!cancelled) setSvgHtml(html)
        }, () => { if (!cancelled) setFailed(true) })
      } catch { if (!cancelled) setFailed(true) }
    }).catch(() => { if (!cancelled) setFailed(true) })
    return () => { cancelled = true }
  }, [smiles])

  if (!smiles) return null

  if (failed) {
    const t = smiles.length > 20 ? smiles.slice(0, 20) + '…' : smiles
    return (
      <div className="flex items-center justify-center bg-gray-50 rounded" style={{ width: IMG_W, height: IMG_H }}>
        <span className="font-mono text-[8px] text-gray-400 px-1 text-center break-all leading-tight" title={smiles}>{t}</span>
      </div>
    )
  }

  if (!svgHtml) {
    return (
      <div className="flex items-center justify-center" style={{ width: IMG_W, height: IMG_H }}>
        <div className="w-4 h-4 border-2 border-gray-200 border-t-purple-400 rounded-full animate-spin" />
      </div>
    )
  }

  return (
    <div
      className="flex items-center justify-center"
      style={{ width: IMG_W, height: IMG_H }}
      dangerouslySetInnerHTML={{ __html: svgHtml }}
    />
  )
}

// --------------------------------------------------
// Pure-CSS connector tree (no JS measurement needed)
// Uses CSS grid for perfect horizontal alignment
// --------------------------------------------------

// Card width constant
const CARD_W = 150

function TreeNode({ node, depth }) {
  const [expanded, setExpanded] = useState(true)
  if (!node) return null

  const isTarget  = node.type === 'target'   || depth === 0
  const isReagent = node.type === 'reagent'  || node.type === 'commercial'
  const hasAvailability = node.available != null  // null/undefined = no data (RDKit)
  const available = node.available === true
  const hasChildren = Array.isArray(node.children) && node.children.length > 0
  const canToggle = hasChildren && !isTarget

  const typeLabel = isTarget ? 'TARGET' : isReagent ? 'REAGENT' : 'INTERMEDIATE'
  const smilesText = node.smiles || ''
  const truncSmiles = smilesText.length > 18 ? smilesText.slice(0, 18) + '…' : smilesText

  // Card styling — neutral gray for reagents without availability data
  let borderColor, labelColor, accentEl
  if (isTarget) {
    borderColor = 'border-purple-400'
    labelColor = 'text-purple-600'
    accentEl = <div className="h-0.5 bg-gradient-to-r from-purple-500 to-purple-700" />
  } else if (isReagent && hasAvailability && available) {
    borderColor = 'border-green-300'
    labelColor = 'text-green-600'
    accentEl = <div className="absolute left-0 top-0 bottom-0 w-0.5 bg-green-400" />
  } else if (isReagent && hasAvailability) {
    borderColor = 'border-orange-300'
    labelColor = 'text-orange-500'
    accentEl = <div className="absolute left-0 top-0 bottom-0 w-0.5 bg-orange-400" />
  } else if (isReagent) {
    // No availability data (RDKit heuristic) — neutral styling
    borderColor = 'border-gray-300'
    labelColor = 'text-gray-500'
    accentEl = <div className="absolute left-0 top-0 bottom-0 w-0.5 bg-gray-300" />
  } else {
    borderColor = 'border-gray-200'
    labelColor = 'text-gray-400'
    accentEl = null
  }

  return (
    <div className="flex flex-col items-center" style={{ minWidth: CARD_W }}>
      {/* ---- Card ---- */}
      <div
        className={`relative bg-white rounded-md shadow-sm border-2 overflow-hidden transition-shadow ${borderColor} ${canToggle ? 'cursor-pointer hover:shadow-md' : ''}`}
        style={{ width: CARD_W }}
        onClick={() => canToggle && setExpanded(e => !e)}
      >
        {accentEl}
        <div className={`px-2 py-1.5 ${isReagent ? 'pl-3' : ''}`}>
          {/* Type label + chevron */}
          <div className="flex items-center justify-between mb-1">
            <span className={`text-[9px] font-bold uppercase tracking-wider ${labelColor}`}>{typeLabel}</span>
            {canToggle && (
              <svg className={`w-3 h-3 text-gray-400 transition-transform ${expanded ? '' : '-rotate-90'}`} fill="currentColor" viewBox="0 0 20 20">
                <path fillRule="evenodd" d="M5.293 7.293a1 1 0 011.414 0L10 10.586l3.293-3.293a1 1 0 111.414 1.414l-4 4a1 1 0 01-1.414 0l-4-4a1 1 0 010-1.414z" clipRule="evenodd" />
              </svg>
            )}
          </div>

          {/* 2D structure */}
          <div className="flex justify-center mb-1">
            <NodeSmilesImage smiles={node.smiles} />
          </div>

          {/* SMILES mono */}
          {smilesText && (
            <p className="font-mono text-[8px] text-gray-400 text-center truncate" title={smilesText}>{truncSmiles}</p>
          )}

          {/* Supplier pill — only when availability data exists */}
          {node.supplier && node.supplier !== 'N/A' && (
            <div className="mt-1 flex justify-center">
              <span className="inline-flex items-center gap-0.5 px-1.5 py-0.5 bg-blue-50 text-blue-600 border border-blue-200 rounded-full text-[9px] font-medium truncate max-w-[130px]">
                <svg className="w-2 h-2 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 21V5a2 2 0 00-2-2H7a2 2 0 00-2 2v16m14 0h2m-2 0h-5m-9 0H3m2 0h5M9 7h1m-1 4h1m4-4h1m-1 4h1m-5 10v-5a1 1 0 011-1h2a1 1 0 011 1v5m-4 0h4" />
                </svg>
                {node.supplier}
              </span>
            </div>
          )}

          {/* Availability — only show when data exists (not for RDKit heuristic) */}
          {isReagent && hasAvailability && (
            <div className="mt-1 flex justify-center">
              {available ? (
                <span className="inline-flex items-center gap-0.5 text-[9px] font-semibold text-green-600">
                  <svg className="w-2.5 h-2.5" fill="currentColor" viewBox="0 0 20 20">
                    <path fillRule="evenodd" d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z" clipRule="evenodd" />
                  </svg>
                  Available
                </span>
              ) : (
                <span className="inline-flex items-center gap-0.5 text-[9px] font-semibold text-orange-500">
                  <svg className="w-2.5 h-2.5" fill="currentColor" viewBox="0 0 20 20">
                    <path fillRule="evenodd" d="M8.257 3.099c.765-1.36 2.722-1.36 3.486 0l5.58 9.92c.75 1.334-.213 2.98-1.742 2.98H4.42c-1.53 0-2.493-1.646-1.743-2.98l5.58-9.92zM11 13a1 1 0 11-2 0 1 1 0 012 0zm-1-8a1 1 0 00-1 1v3a1 1 0 002 0V6a1 1 0 00-1-1z" clipRule="evenodd" />
                  </svg>
                  Not available
                </span>
              )}
            </div>
          )}
        </div>
      </div>

      {/* ---- Connector + children ---- */}
      {hasChildren && expanded && (
        <div className="flex flex-col items-center w-full">
          {/* Vertical line from card */}
          <div className="w-px h-3 bg-gray-300" />

          {/* Reaction pill */}
          {node.reaction && (
            <>
              <div className="px-2 py-0.5 rounded-full text-[9px] text-purple-600 bg-purple-50 border border-purple-200 font-medium whitespace-nowrap max-w-[160px] truncate" title={node.reaction}>
                {node.reaction}
              </div>
              <div className="w-px h-2 bg-gray-300" />
            </>
          )}

          {/* Arrow */}
          <svg className="w-2.5 h-2 text-gray-400" viewBox="0 0 10 8" fill="currentColor">
            <path d="M5 8L0 0h10z" />
          </svg>

          {/* Children layout */}
          {node.children.length === 1 ? (
            <>
              <div className="w-px h-2 bg-gray-300" />
              <TreeNode node={node.children[0]} depth={depth + 1} />
            </>
          ) : (
            <ForkLayout>
              {node.children.map((child, i) => (
                <TreeNode key={i} node={child} depth={depth + 1} />
              ))}
            </ForkLayout>
          )}
        </div>
      )}
    </div>
  )
}

// --------------------------------------------------
// ForkLayout — uses ref measurement for horizontal bar
// --------------------------------------------------
function ForkLayout({ children }) {
  const containerRef = useRef(null)
  const [bar, setBar] = useState(null)
  const childCount = React.Children.count(children)

  const measure = useCallback(() => {
    if (!containerRef.current || childCount < 2) return
    const cols = containerRef.current.querySelectorAll('[data-fork-col]')
    if (cols.length < 2) return
    const parentRect = containerRef.current.getBoundingClientRect()
    const first = cols[0].getBoundingClientRect()
    const last = cols[cols.length - 1].getBoundingClientRect()
    const l = first.left + first.width / 2 - parentRect.left
    const r = last.left + last.width / 2 - parentRect.left
    setBar({ left: l, width: r - l })
  }, [childCount])

  // Measure after layout + on resize
  useLayoutEffect(() => {
    measure()
  }, [measure, children])

  useEffect(() => {
    if (!containerRef.current) return
    const ro = new ResizeObserver(measure)
    ro.observe(containerRef.current)
    return () => ro.disconnect()
  }, [measure])

  return (
    <div className="relative pt-px" ref={containerRef}>
      {/* Horizontal bar */}
      {bar && bar.width > 0 && (
        <div
          className="absolute top-0 h-px bg-gray-300"
          style={{ left: bar.left, width: bar.width }}
        />
      )}
      {/* Children columns */}
      <div className="flex items-start justify-center gap-3">
        {React.Children.map(children, (child, i) => (
          <div key={i} data-fork-col className="flex flex-col items-center">
            <div className="w-px h-3 bg-gray-300" />
            {child}
          </div>
        ))}
      </div>
    </div>
  )
}

// --------------------------------------------------
// Fullscreen modal for tree view
// --------------------------------------------------
function FullscreenTreeModal({ open, onClose, children }) {
  useEffect(() => {
    if (!open) return
    const handler = (e) => e.key === 'Escape' && onClose()
    window.addEventListener('keydown', handler)
    return () => window.removeEventListener('keydown', handler)
  }, [open, onClose])

  if (!open) return null

  return createPortal(
    <div className="fixed inset-0 z-[9999] flex flex-col bg-white">
      {/* Top bar */}
      <div className="flex items-center justify-between px-4 py-2 border-b border-gray-200 bg-gray-50 flex-shrink-0">
        <span className="text-sm font-semibold text-gray-700 flex items-center gap-2">
          <svg className="w-4 h-4 text-purple-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 11H5m14 0a2 2 0 012 2v6a2 2 0 01-2 2H5a2 2 0 01-2-2v-6a2 2 0 012-2m14 0V9a2 2 0 00-2-2M5 11V9a2 2 0 012-2m0 0V5a2 2 0 012-2h6a2 2 0 012 2v2M7 7h10" />
          </svg>
          Retrosynthesis Route
        </span>
        <button
          onClick={onClose}
          className="flex items-center gap-1 px-3 py-1.5 text-sm text-gray-600 hover:text-gray-900 hover:bg-gray-100 rounded-md transition-colors"
        >
          <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
          </svg>
          Close
          <kbd className="ml-1 text-[10px] px-1 py-0.5 bg-gray-200 rounded text-gray-500">Esc</kbd>
        </button>
      </div>
      {/* Scrollable tree area */}
      <div className="flex-1 overflow-auto p-6">
        {children}
      </div>
    </div>,
    document.body
  )
}

// --------------------------------------------------
// Linear step timeline (used when no tree data)
// --------------------------------------------------
function StepTimeline({ steps }) {
  if (!steps || !steps.length) return null
  return (
    <div className="space-y-3 px-4 py-3">
      {steps.map((step, i) => {
        const confColor = confidenceColor(step.confidence)
        const reactants = step.reactant_names || step.reactants || []
        return (
          <div key={i} className="flex gap-3">
            {/* Step number + connector */}
            <div className="flex flex-col items-center flex-shrink-0">
              <div className="w-7 h-7 rounded-full bg-purple-600 text-white text-sm font-bold flex items-center justify-center shadow-sm">
                {i + 1}
              </div>
              {i < steps.length - 1 && (
                <div className="w-px flex-1 bg-gray-200 my-1 min-h-[1.5rem]" />
              )}
            </div>

            {/* Step card */}
            <div className="flex-1 pb-3">
              <div className="bg-white border border-gray-200 rounded-lg p-3 shadow-sm">
                <div className="flex items-start justify-between gap-2 mb-2">
                  <span className="text-sm font-semibold text-gray-800 truncate">
                    {step.reaction || `Step ${i + 1}`}
                  </span>
                  {step.confidence !== undefined && step.confidence !== null && (
                    <span className={`flex-shrink-0 text-sm px-1.5 py-0.5 rounded-full font-medium ${confColor.bg} ${confColor.text}`}>
                      {Math.round(step.confidence * 100)}%
                    </span>
                  )}
                </div>

                {reactants.length > 0 && (
                  <div className="mb-1.5">
                    <p className="text-[10px] text-gray-400 uppercase tracking-wide mb-1">Reagents</p>
                    <div className="flex flex-wrap gap-1">
                      {reactants.map((r, ri) => (
                        <span
                          key={ri}
                          className="px-1.5 py-0.5 bg-green-50 text-green-700 border border-green-200 rounded text-[10px] font-mono truncate max-w-[140px]"
                          title={r}
                        >
                          {r}
                        </span>
                      ))}
                    </div>
                  </div>
                )}

                {step.conditions && (
                  <p className="text-[10px] text-gray-400 italic truncate" title={step.conditions}>
                    {step.conditions}
                  </p>
                )}
              </div>
            </div>
          </div>
        )
      })}
    </div>
  )
}

// --------------------------------------------------
// SynthesisTree — visual retrosynthesis route component
// Props: { synthesisRoute } — object with n_steps, confidence, estimated_cost,
//   all_reagents_available, steps[], tree{}, all_reagents[]
// --------------------------------------------------
export default function SynthesisTree({ synthesisRoute }) {
  const [viewMode, setViewMode] = useState('tree') // 'tree' | 'steps'
  const [fullscreen, setFullscreen] = useState(false)

  if (!synthesisRoute) {
    return (
      <div className="card p-6 flex flex-col items-center justify-center h-48 text-gray-300">
        <svg className="w-10 h-10 mb-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
            d="M19 11H5m14 0a2 2 0 012 2v6a2 2 0 01-2 2H5a2 2 0 01-2-2v-6a2 2 0 012-2m14 0V9a2 2 0 00-2-2M5 11V9a2 2 0 012-2m0 0V5a2 2 0 012-2h6a2 2 0 012 2v2M7 7h10" />
        </svg>
        <p className="text-sm">Synthesis route not available</p>
      </div>
    )
  }

  const hasTree  = Boolean(synthesisRoute.tree)
  const hasSteps = Array.isArray(synthesisRoute.steps) && synthesisRoute.steps.length > 0

  const treeContent = hasTree ? (
    <div className="flex justify-center py-2 min-w-max">
      <TreeNode node={synthesisRoute.tree} depth={0} />
    </div>
  ) : null

  return (
    <div className="card overflow-hidden">
      {/* Header */}
      <div className="px-4 py-3 bg-bx-surface flex items-center justify-between">
        <h3 className="text-white font-semibold text-sm flex items-center gap-2">
          <svg className="w-4 h-4 text-bx-mint" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
              d="M19 11H5m14 0a2 2 0 012 2v6a2 2 0 01-2 2H5a2 2 0 01-2-2v-6a2 2 0 012-2m14 0V9a2 2 0 00-2-2M5 11V9a2 2 0 012-2m0 0V5a2 2 0 012-2h6a2 2 0 012 2v2M7 7h10" />
          </svg>
          Retrosynthesis Route
          <InfoTip text="AI-planned synthesis pathway showing how to build this molecule step-by-step from commercially available reagents." />
        </h3>

        <div className="flex items-center gap-2">
          {/* Toggle tree / steps only when both are available */}
          {hasTree && hasSteps && (
            <div className="flex rounded-md overflow-hidden border border-white/20">
              {[{ key: 'tree', label: 'Tree' }, { key: 'steps', label: 'Steps' }].map(({ key, label }) => (
                <button
                  key={key}
                  onClick={() => setViewMode(key)}
                  className={`px-3 py-1 text-sm font-medium transition-colors ${
                    viewMode === key
                      ? 'bg-white text-bx-light-text'
                      : 'text-white/70 hover:text-white hover:bg-white/10'
                  }`}
                >
                  {label}
                </button>
              ))}
            </div>
          )}

          {/* Fullscreen button */}
          {hasTree && viewMode === 'tree' && (
            <button
              onClick={() => setFullscreen(true)}
              className="p-1.5 text-white/60 hover:text-white hover:bg-white/10 rounded transition-colors"
              title="Fullscreen"
            >
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 8V4m0 0h4M4 4l5 5m11-1V4m0 0h-4m4 0l-5 5M4 16v4m0 0h4m-4 0l5-5m11 5l-5-5m5 5v-4m0 4h-4" />
              </svg>
            </button>
          )}
        </div>
      </div>

      {/* Summary bar */}
      <SummaryBar
        nSteps={synthesisRoute.n_steps}
        confidence={synthesisRoute.confidence}
        estimatedCost={synthesisRoute.estimated_cost}
        allReagentsAvailable={synthesisRoute.all_reagents_available}
      />

      {/* Main content */}
      <div className="p-4 overflow-x-auto">
        {viewMode === 'tree' && treeContent ? treeContent : hasSteps ? (
          <StepTimeline steps={synthesisRoute.steps} />
        ) : (
          <div className="py-6 text-center text-gray-400 text-sm">
            No synthesis data available.
          </div>
        )}
      </div>

      {/* All reagents list */}
      {synthesisRoute.all_reagents && synthesisRoute.all_reagents.length > 0 && (
        <div className="px-4 pb-4">
          <p className="text-sm font-medium text-gray-400 uppercase tracking-wide mb-2 flex items-center">
            All Reagents
            <InfoTip text="Complete list of starting materials and reagents needed for the synthesis." />
          </p>
          <div className="flex flex-wrap gap-1.5">
            {synthesisRoute.all_reagents.map((r, i) => (
              <span
                key={i}
                className="px-2 py-0.5 bg-gray-50 border border-gray-200 rounded text-sm text-gray-600 font-mono truncate max-w-[160px]"
                title={r}
              >
                {r}
              </span>
            ))}
          </div>
        </div>
      )}

      {/* Cost breakdown section */}
      {synthesisRoute.cost_estimate && (
        <div className="px-4 pb-4">
          <p className="text-sm font-medium text-gray-400 uppercase tracking-wide mb-2 flex items-center">
            Cost Estimate
            <InfoTip text="Estimated synthesis cost in USD, broken down by reagent and labor costs." />
          </p>
          <div className="bg-gray-50 border border-gray-200 rounded-lg p-3 space-y-1.5">
            {synthesisRoute.cost_estimate.total_cost_usd !== undefined && (
              <div className="flex items-center justify-between">
                <span className="text-sm font-semibold text-gray-700">Total Cost</span>
                <span className="text-sm font-bold text-bx-light-text font-mono">
                  ${Number(synthesisRoute.cost_estimate.total_cost_usd).toLocaleString('en-US', { minimumFractionDigits: 2, maximumFractionDigits: 2 })}
                </span>
              </div>
            )}
            {synthesisRoute.cost_estimate.reagent_cost !== undefined && (
              <div className="flex items-center justify-between">
                <span className="text-sm text-gray-500">Reagents</span>
                <span className="text-sm font-medium text-gray-700 font-mono">
                  ${Number(synthesisRoute.cost_estimate.reagent_cost).toLocaleString('en-US', { minimumFractionDigits: 2, maximumFractionDigits: 2 })}
                </span>
              </div>
            )}
            {synthesisRoute.cost_estimate.labor_cost !== undefined && (
              <div className="flex items-center justify-between">
                <span className="text-sm text-gray-500">Labor</span>
                <span className="text-sm font-medium text-gray-700 font-mono">
                  ${Number(synthesisRoute.cost_estimate.labor_cost).toLocaleString('en-US', { minimumFractionDigits: 2, maximumFractionDigits: 2 })}
                </span>
              </div>
            )}
          </div>
        </div>
      )}

      {/* Reagent availability table */}
      {Array.isArray(synthesisRoute.reagent_availability) && synthesisRoute.reagent_availability.length > 0 && (
        <div className="px-4 pb-4">
          <p className="text-sm font-medium text-gray-400 uppercase tracking-wide mb-2 flex items-center">
            Reagent Availability
            <InfoTip text="Commercial availability and supplier information for each required reagent." />
          </p>
          <div className="border border-gray-200 rounded-lg overflow-hidden">
            <table className="w-full text-sm">
              <thead>
                <tr className="bg-gray-50 border-b border-gray-200">
                  <th className="text-left px-3 py-2 font-semibold text-gray-500 uppercase tracking-wide text-[10px]">Reagent</th>
                  <th className="text-left px-3 py-2 font-semibold text-gray-500 uppercase tracking-wide text-[10px] hidden sm:table-cell">Supplier</th>
                  <th className="text-left px-3 py-2 font-semibold text-gray-500 uppercase tracking-wide text-[10px] hidden sm:table-cell">Catalog ID</th>
                  <th className="text-center px-3 py-2 font-semibold text-gray-500 uppercase tracking-wide text-[10px]">Status</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-gray-100">
                {synthesisRoute.reagent_availability.map((reagent, idx) => {
                  const available = reagent.available !== false
                  const name = reagent.name || reagent.smiles || `Reagent ${idx + 1}`
                  return (
                    <tr key={idx} className={available ? '' : 'bg-orange-50/40'}>
                      <td className="px-3 py-2 font-medium text-gray-700 truncate max-w-[120px]" title={name}>
                        {name}
                      </td>
                      <td className="px-3 py-2 text-gray-500 hidden sm:table-cell truncate max-w-[100px]">
                        {reagent.supplier || '—'}
                      </td>
                      <td className="px-3 py-2 font-mono text-gray-400 hidden sm:table-cell">
                        {reagent.catalog_id || '—'}
                      </td>
                      <td className="px-3 py-2 text-center">
                        {available ? (
                          <span className="inline-flex items-center gap-1 px-1.5 py-0.5 rounded text-[10px] font-semibold bg-green-100 text-green-700 border border-green-200">
                            <svg className="w-2.5 h-2.5" fill="currentColor" viewBox="0 0 20 20">
                              <path fillRule="evenodd" d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z" clipRule="evenodd" />
                            </svg>
                            Available
                          </span>
                        ) : (
                          <span className="inline-flex items-center gap-1 px-1.5 py-0.5 rounded text-[10px] font-semibold bg-red-100 text-red-700 border border-red-200">
                            <svg className="w-2.5 h-2.5" fill="currentColor" viewBox="0 0 20 20">
                              <path fillRule="evenodd" d="M4.293 4.293a1 1 0 011.414 0L10 8.586l4.293-4.293a1 1 0 111.414 1.414L11.414 10l4.293 4.293a1 1 0 01-1.414 1.414L10 11.414l-4.293 4.293a1 1 0 01-1.414-1.414L8.586 10 4.293 5.707a1 1 0 010-1.414z" clipRule="evenodd" />
                            </svg>
                            Not Available
                          </span>
                        )}
                      </td>
                    </tr>
                  )
                })}
              </tbody>
            </table>
          </div>
        </div>
      )}

      {/* Legend */}
      <div className="px-4 py-2 bg-gray-50 border-t border-gray-100 flex flex-wrap gap-4 text-sm text-gray-400">
        <span className="flex items-center gap-1.5">
          <span className="w-3 h-3 rounded border-2 border-purple-400 bg-white inline-block relative overflow-hidden">
            <span className="absolute top-0 left-0 right-0 h-0.5 bg-gradient-to-r from-purple-500 to-purple-700" />
          </span>
          Target
        </span>
        <span className="flex items-center gap-1.5">
          <span className="w-3 h-3 rounded border border-gray-200 bg-white inline-block" />
          Intermediate
        </span>
        {synthesisRoute.all_reagents_available != null ? (
          <>
            <span className="flex items-center gap-1.5">
              <span className="w-3 h-3 rounded border-2 border-green-300 bg-white inline-block relative overflow-hidden">
                <span className="absolute left-0 top-0 bottom-0 w-0.5 bg-green-400" />
              </span>
              Available
            </span>
            <span className="flex items-center gap-1.5">
              <span className="w-3 h-3 rounded border-2 border-orange-300 bg-white inline-block relative overflow-hidden">
                <span className="absolute left-0 top-0 bottom-0 w-0.5 bg-orange-400" />
              </span>
              Missing
            </span>
          </>
        ) : (
          <span className="flex items-center gap-1.5">
            <span className="w-3 h-3 rounded border-2 border-gray-300 bg-white inline-block relative overflow-hidden">
              <span className="absolute left-0 top-0 bottom-0 w-0.5 bg-gray-300" />
            </span>
            Reagent
          </span>
        )}
      </div>

      {/* Fullscreen modal */}
      <FullscreenTreeModal open={fullscreen} onClose={() => setFullscreen(false)}>
        <div className="flex justify-center min-w-max">
          {hasTree && <TreeNode node={synthesisRoute.tree} depth={0} />}
        </div>
      </FullscreenTreeModal>
    </div>
  )
}
