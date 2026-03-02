import React, { useState, useCallback, useEffect, useRef } from 'react'

// ---------------------------------------------------------------------------
// Operator definitions per column type
// ---------------------------------------------------------------------------
const OPERATORS_BY_TYPE = {
  number: [
    { value: '>', label: '>' },
    { value: '>=', label: '>=' },
    { value: '<', label: '<' },
    { value: '<=', label: '<=' },
    { value: '=', label: '=' },
    { value: 'between', label: 'between' },
  ],
  text: [
    { value: 'contains', label: 'contains' },
    { value: 'not_contains', label: 'not contains' },
    { value: 'equals', label: 'equals' },
    { value: 'starts_with', label: 'starts with' },
  ],
  smiles: [
    { value: 'contains', label: 'contains' },
    { value: 'not_contains', label: 'not contains' },
    { value: 'equals', label: 'equals' },
    { value: 'starts_with', label: 'starts with' },
  ],
  boolean: [
    { value: 'is_yes', label: 'is Yes' },
    { value: 'is_no', label: 'is No' },
  ],
}

// Fallback for unrecognized types
const DEFAULT_OPERATORS = OPERATORS_BY_TYPE.text

function getOperators(type) {
  return OPERATORS_BY_TYPE[type] || DEFAULT_OPERATORS
}

// ---------------------------------------------------------------------------
// Build a filter function from a single rule
// ---------------------------------------------------------------------------
function buildRuleFn(rule) {
  const { column, operator, value, valueTo, colType } = rule
  if (!column || !operator) return null

  // Boolean rules don't need a value input
  if (colType === 'boolean') {
    return (mol) => {
      if (!mol) return false
      const v = mol[column]
      if (operator === 'is_yes') return !!v
      if (operator === 'is_no') return !v
      return true
    }
  }

  // Number rules
  if (colType === 'number') {
    const num = parseFloat(value)
    if (operator === 'between') {
      const numTo = parseFloat(valueTo)
      if (isNaN(num) || isNaN(numTo)) return null
      return (mol) => {
        if (!mol) return false
        const v = parseFloat(mol[column])
        if (isNaN(v)) return false
        return v >= Math.min(num, numTo) && v <= Math.max(num, numTo)
      }
    }
    if (isNaN(num)) return null
    return (mol) => {
      if (!mol) return false
      const v = parseFloat(mol[column])
      if (isNaN(v)) return false
      switch (operator) {
        case '>':  return v > num
        case '>=': return v >= num
        case '<':  return v < num
        case '<=': return v <= num
        case '=':  return v === num
        default:   return true
      }
    }
  }

  // Text / SMILES rules
  const str = (value || '').toLowerCase()
  if (!str && operator !== 'equals') return null

  return (mol) => {
    const v = (mol[column] ?? '').toString().toLowerCase()
    switch (operator) {
      case 'contains':     return v.includes(str)
      case 'not_contains': return !v.includes(str)
      case 'equals':       return v === str
      case 'starts_with':  return v.startsWith(str)
      default:             return true
    }
  }
}

// ---------------------------------------------------------------------------
// Build composed filter from all rules + mode
// ---------------------------------------------------------------------------
function buildFilterFn(rules, mode) {
  const fns = rules.map(buildRuleFn).filter(Boolean)
  if (fns.length === 0) return null

  if (mode === 'OR') {
    return (mol) => fns.some(fn => fn(mol))
  }
  return (mol) => fns.every(fn => fn(mol))
}

// ---------------------------------------------------------------------------
// Create a fresh empty rule
// ---------------------------------------------------------------------------
let ruleIdCounter = 0
function createEmptyRule() {
  return { id: ++ruleIdCounter, column: '', operator: '', value: '', valueTo: '', colType: 'text' }
}

// ---------------------------------------------------------------------------
// AdvancedFilterBuilder
// ---------------------------------------------------------------------------
export default function AdvancedFilterBuilder({ columns, onFilterChange, isOpen, onToggle }) {
  const [mode, setMode] = useState('AND') // 'AND' | 'OR'
  const [rules, setRules] = useState(() => [createEmptyRule()])

  // Filterable columns — only number, text, boolean, smiles
  const filterableColumns = (columns || []).filter(c =>
    ['number', 'text', 'boolean', 'smiles'].includes(c.type)
  )

  // Emit filter whenever rules or mode change
  const prevFilterRef = useRef(null)
  useEffect(() => {
    const fn = buildFilterFn(rules, mode)
    // Avoid calling if both null
    if (fn === null && prevFilterRef.current === null) return
    prevFilterRef.current = fn
    onFilterChange(fn)
  }, [rules, mode, onFilterChange])

  const handleColumnChange = useCallback((ruleId, colKey) => {
    setRules(prev => prev.map(r => {
      if (r.id !== ruleId) return r
      const colDef = filterableColumns.find(c => c.key === colKey)
      const colType = colDef?.type || 'text'
      const operators = getOperators(colType)
      return {
        ...r,
        column: colKey,
        colType,
        operator: operators[0]?.value || '',
        value: '',
        valueTo: '',
      }
    }))
  }, [filterableColumns])

  const handleOperatorChange = useCallback((ruleId, op) => {
    setRules(prev => prev.map(r => {
      if (r.id !== ruleId) return r
      return { ...r, operator: op, value: r.value, valueTo: '' }
    }))
  }, [])

  const handleValueChange = useCallback((ruleId, val) => {
    setRules(prev => prev.map(r => r.id === ruleId ? { ...r, value: val } : r))
  }, [])

  const handleValueToChange = useCallback((ruleId, val) => {
    setRules(prev => prev.map(r => r.id === ruleId ? { ...r, valueTo: val } : r))
  }, [])

  const handleDeleteRule = useCallback((ruleId) => {
    setRules(prev => {
      const next = prev.filter(r => r.id !== ruleId)
      return next.length === 0 ? [createEmptyRule()] : next
    })
  }, [])

  const handleAddRule = useCallback(() => {
    setRules(prev => [...prev, createEmptyRule()])
  }, [])

  const handleClearAll = useCallback(() => {
    setRules([createEmptyRule()])
  }, [])

  if (!isOpen) return null

  return (
    <div className="bg-white border border-gray-200 rounded-xl p-4 space-y-3 shadow-xs">
      {/* Header row: title + AND/OR toggle + clear */}
      <div className="flex items-center justify-between">
        <div className="flex items-center gap-3">
          <span className="text-xs font-semibold text-gray-700">Advanced Filter</span>

          {/* AND / OR toggle */}
          <div className="flex rounded-lg overflow-hidden border border-gray-200">
            <button
              onClick={() => setMode('AND')}
              className={`px-2.5 py-1 text-xs font-medium transition-colors ${
                mode === 'AND'
                  ? 'bg-bx-surface text-white'
                  : 'bg-white text-gray-500 hover:bg-gray-50'
              }`}
            >
              AND
            </button>
            <button
              onClick={() => setMode('OR')}
              className={`px-2.5 py-1 text-xs font-medium transition-colors ${
                mode === 'OR'
                  ? 'bg-bx-surface text-white'
                  : 'bg-white text-gray-500 hover:bg-gray-50'
              }`}
            >
              OR
            </button>
          </div>
        </div>

        <button
          onClick={handleClearAll}
          className="text-xs text-gray-400 hover:text-bx-red transition-colors"
        >
          Clear all
        </button>
      </div>

      {/* Rules */}
      <div className="space-y-2">
        {rules.map((rule) => {
          const operators = getOperators(rule.colType)
          const needsValue = rule.colType !== 'boolean'
          const isBetween = rule.operator === 'between'

          return (
            <div key={rule.id} className="flex items-center gap-2">
              {/* Column select */}
              <select
                value={rule.column}
                onChange={(e) => handleColumnChange(rule.id, e.target.value)}
                className="text-xs border border-gray-200 rounded-lg px-2 py-1.5 bg-white text-gray-700 min-w-[120px] focus:outline-none focus:ring-1 focus:ring-bx-cyan/40"
              >
                <option value="">Column...</option>
                {filterableColumns.map(col => (
                  <option key={col.key} value={col.key}>{col.label}</option>
                ))}
              </select>

              {/* Operator select */}
              <select
                value={rule.operator}
                onChange={(e) => handleOperatorChange(rule.id, e.target.value)}
                className="text-xs border border-gray-200 rounded-lg px-2 py-1.5 bg-white text-gray-700 min-w-[100px] focus:outline-none focus:ring-1 focus:ring-bx-cyan/40"
                disabled={!rule.column}
              >
                {operators.map(op => (
                  <option key={op.value} value={op.value}>{op.label}</option>
                ))}
              </select>

              {/* Value input */}
              {needsValue && (
                <>
                  <input
                    type={rule.colType === 'number' ? 'number' : 'text'}
                    value={rule.value}
                    onChange={(e) => handleValueChange(rule.id, e.target.value)}
                    placeholder={isBetween ? 'min' : 'value'}
                    className="text-xs border border-gray-200 rounded-lg px-2 py-1.5 bg-white text-gray-700 w-[80px] focus:outline-none focus:ring-1 focus:ring-bx-cyan/40"
                    disabled={!rule.column}
                  />
                  {isBetween && (
                    <>
                      <span className="text-xs text-gray-400">to</span>
                      <input
                        type="number"
                        value={rule.valueTo}
                        onChange={(e) => handleValueToChange(rule.id, e.target.value)}
                        placeholder="max"
                        className="text-xs border border-gray-200 rounded-lg px-2 py-1.5 bg-white text-gray-700 w-[80px] focus:outline-none focus:ring-1 focus:ring-bx-cyan/40"
                      />
                    </>
                  )}
                </>
              )}

              {/* Delete button */}
              <button
                onClick={() => handleDeleteRule(rule.id)}
                className="flex-shrink-0 p-1 rounded-md text-gray-300 hover:text-bx-red hover:bg-red-50 transition-colors"
                title="Remove rule"
              >
                <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
                </svg>
              </button>
            </div>
          )
        })}
      </div>

      {/* Add rule button */}
      <button
        onClick={handleAddRule}
        className="flex items-center gap-1.5 text-xs text-bx-cyan hover:text-bx-blue transition-colors font-medium"
      >
        <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 4v16m8-8H4" />
        </svg>
        Add rule
      </button>
    </div>
  )
}
