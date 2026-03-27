import React from 'react'

const MAX_PROPERTIES = 3

/**
 * Chips/pills selector for multi-property SAR comparison.
 * Max 3 selected simultaneously. First = primary (full analysis), others = comparison.
 */
export default function MultiPropertySelector({ availableProperties, selectedKeys, onSelectionChange }) {
  if (!availableProperties?.length) return null

  const toggle = (key) => {
    if (selectedKeys.includes(key)) {
      // Deselect — but always keep at least 1
      if (selectedKeys.length <= 1) return
      onSelectionChange(selectedKeys.filter(k => k !== key))
    } else {
      // Select — respect max
      if (selectedKeys.length >= MAX_PROPERTIES) {
        // Replace last selected
        onSelectionChange([...selectedKeys.slice(0, MAX_PROPERTIES - 1), key])
      } else {
        onSelectionChange([...selectedKeys, key])
      }
    }
  }

  return (
    <div className="card px-5 py-3">
      <div className="flex items-center gap-3 mb-2">
        <h2 className="text-xs font-semibold text-gray-500 uppercase tracking-wide">Compare Properties</h2>
        <span className="text-[10px] text-gray-400">
          Select up to {MAX_PROPERTIES} properties to compare (first = primary)
        </span>
      </div>
      <div className="flex flex-wrap gap-2">
        {availableProperties.map((p) => {
          const isSelected = selectedKeys.includes(p.key)
          const isPrimary = selectedKeys[0] === p.key
          const idx = selectedKeys.indexOf(p.key)
          return (
            <button
              key={p.key}
              onClick={() => toggle(p.key)}
              className={`
                relative flex items-center gap-1.5 px-3 py-1.5 rounded-full text-xs font-medium
                border transition-all cursor-pointer
                ${isPrimary
                  ? 'border-amber-400 bg-amber-50 text-amber-700 ring-2 ring-amber-200'
                  : isSelected
                    ? 'border-blue-300 bg-blue-50 text-blue-700 ring-1 ring-blue-200'
                    : p.auto_eligible
                      ? 'border-gray-200 bg-white text-gray-600 hover:border-gray-300 hover:bg-gray-50'
                      : 'border-gray-100 bg-gray-50 text-gray-400 hover:border-gray-200'
                }
              `}
              title={`${p.label}: ${p.n_values} values, ${Math.round(p.coverage * 100)}% coverage${p.auto_eligible ? '' : ' (low coverage)'}`}
            >
              {isSelected && (
                <span className={`inline-flex items-center justify-center w-4 h-4 rounded-full text-[9px] font-bold text-white ${isPrimary ? 'bg-amber-500' : 'bg-blue-500'}`}>
                  {idx + 1}
                </span>
              )}
              {p.label}
              <span className={`text-[10px] ${isSelected ? 'opacity-75' : 'text-gray-400'}`}>
                {Math.round(p.coverage * 100)}%
              </span>
            </button>
          )
        })}
      </div>
      {selectedKeys.length > 1 && (
        <div className="flex items-center gap-2 mt-2 text-[10px] text-gray-400">
          <span className="inline-block w-3 h-3 rounded-full bg-amber-400" /> Primary
          <span className="inline-block w-3 h-3 rounded-full bg-blue-400 ml-2" /> Comparison
        </div>
      )}
    </div>
  )
}
