import React from 'react'
import InfoTip from './InfoTip.jsx'

// --------------------------------------------------
// Interaction type config: color + label
// --------------------------------------------------
const INTERACTION_TYPES = {
  HBDonor:     { label: 'H-Bond Donor',    bg: 'bg-blue-100',   text: 'text-blue-700',   dot: '#3b82f6' },
  HBAcceptor:  { label: 'H-Bond Acceptor', bg: 'bg-red-100',    text: 'text-red-700',    dot: '#ef4444' },
  Hydrophobic: { label: 'Hydrophobic',     bg: 'bg-gray-100',   text: 'text-gray-600',   dot: '#6b7280' },
  PiStacking:  { label: 'Pi Stacking',     bg: 'bg-purple-100', text: 'text-purple-700', dot: '#8b5cf6' },
  CationPi:    { label: 'Cation-Pi',       bg: 'bg-orange-100', text: 'text-orange-700', dot: '#f97316' },
  VDW:         { label: 'Van der Waals',   bg: 'bg-gray-50',    text: 'text-gray-400',   dot: '#d1d5db' },
}

function getTypeConfig(type) {
  return INTERACTION_TYPES[type] || {
    label: type,
    bg: 'bg-gray-100',
    text: 'text-gray-600',
    dot: '#9ca3af',
  }
}

// --------------------------------------------------
// Interaction quality bar
// --------------------------------------------------
function QualityBar({ value }) {
  const pct = Math.round((value || 0) * 100)
  let colorClass = 'bg-red-500'
  let labelColor = 'text-red-600'
  let qualityLabel = 'Low'
  if (value >= 0.6) {
    colorClass = 'bg-bx-mint'
    labelColor = 'text-green-700'
    qualityLabel = 'Good'
  } else if (value >= 0.3) {
    colorClass = 'bg-yellow-400'
    labelColor = 'text-yellow-700'
    qualityLabel = 'Moderate'
  }

  return (
    <div className="flex items-center gap-3">
      <div className="flex-1 h-2.5 bg-gray-100 rounded-full overflow-hidden">
        <div
          className={`h-full rounded-full transition-all duration-500 ${colorClass}`}
          style={{ width: `${pct}%` }}
        />
      </div>
      <span className={`text-sm font-bold w-16 text-right ${labelColor}`}>
        {qualityLabel} ({pct}%)
      </span>
    </div>
  )
}

// --------------------------------------------------
// Type badge
// --------------------------------------------------
function TypeBadge({ type }) {
  const cfg = getTypeConfig(type)
  return (
    <span className={`inline-flex items-center gap-1 px-2 py-0.5 rounded-full text-sm font-semibold ${cfg.bg} ${cfg.text}`}>
      <span
        className="w-1.5 h-1.5 rounded-full flex-shrink-0"
        style={{ backgroundColor: cfg.dot }}
      />
      {cfg.label}
    </span>
  )
}

// --------------------------------------------------
// Functional checkmark
// --------------------------------------------------
function FunctionalIcon({ is_functional }) {
  if (is_functional) {
    return (
      <svg className="w-4 h-4 text-bx-mint" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M5 13l4 4L19 7" />
      </svg>
    )
  }
  return <span className="text-gray-300 text-sm leading-none">&#8212;</span>
}

// --------------------------------------------------
// Main InteractionView
// --------------------------------------------------
export default function InteractionView({ interactions }) {
  if (!interactions) {
    return (
      <div className="flex flex-col items-center justify-center py-10 text-gray-300">
        <svg className="w-10 h-10 mb-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
            d="M9.172 16.172a4 4 0 015.656 0M9 10h.01M15 10h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
        </svg>
        <p className="text-sm">Interaction data not available</p>
      </div>
    )
  }

  const {
    interactions: rows = [],
    functional_contacts = 0,
    total_functional = 0,
    interaction_quality = 0,
    key_hbonds = 0,
    method = 'mock',
    summary = '',
  } = interactions

  const isMock = method === 'mock'

  return (
    <div className="space-y-4">
      {/* Header */}
      <div className="flex items-center gap-2">
        <h4 className="text-sm font-semibold text-gray-700 uppercase tracking-wide">
          Protein-Ligand Interactions
        </h4>
        <InfoTip text="Interaction analysis reveals which protein residues make direct contact with the ligand. Functional contacts with key residues (catalytic site, binding motifs) are especially important for predicting biological activity." />
      </div>

      {/* Quality bar */}
      <div className="p-3 bg-gray-50 rounded-lg space-y-2">
        <div className="flex items-center justify-between text-sm text-gray-500 mb-1">
          <span className="font-medium">Interaction Quality</span>
          <InfoTip text="Overall quality score based on the fraction of key functional residues contacted and the number of hydrogen bonds formed." />
        </div>
        <QualityBar value={interaction_quality} />
        <p className="text-sm text-gray-500 mt-1">
          {summary || `${functional_contacts}/${total_functional} functional residues contacted — ${key_hbonds} key H-bond${key_hbonds !== 1 ? 's' : ''}`}
        </p>
      </div>

      {/* Interaction table */}
      {rows.length > 0 ? (
        <div className="overflow-hidden rounded-lg border border-gray-100">
          <table className="w-full text-sm">
            <thead>
              <tr className="bg-gray-50 border-b border-gray-100">
                <th className="text-left px-3 py-2 font-semibold text-gray-500 uppercase tracking-wide">
                  Residue
                </th>
                <th className="text-left px-3 py-2 font-semibold text-gray-500 uppercase tracking-wide">
                  Type
                </th>
                <th className="text-center px-3 py-2 font-semibold text-gray-500 uppercase tracking-wide">
                  Functional
                  <InfoTip text="Residue is part of the known functional/catalytic site of the protein." />
                </th>
              </tr>
            </thead>
            <tbody className="divide-y divide-gray-50">
              {rows.map((row, idx) => (
                <tr
                  key={`${row.residue}-${idx}`}
                  className={`hover:bg-gray-50 transition-colors ${row.is_functional ? 'bg-green-50/30' : ''}`}
                >
                  <td className="px-3 py-2 font-mono font-semibold text-bx-light-text">
                    {row.residue}
                    {row.residue_number !== undefined && row.residue_number !== null && (
                      <span className="text-gray-400 font-normal ml-1 text-[10px]">
                        #{row.residue_number}
                      </span>
                    )}
                  </td>
                  <td className="px-3 py-2">
                    <TypeBadge type={row.type} />
                  </td>
                  <td className="px-3 py-2 text-center">
                    <FunctionalIcon is_functional={row.is_functional} />
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      ) : (
        <p className="text-sm text-gray-400 text-center py-4">No interaction contacts recorded</p>
      )}

      {/* Method badge */}
      <div className="flex justify-end">
        <span className={`inline-flex items-center gap-1.5 px-2.5 py-1 rounded-full text-sm font-medium border ${
          isMock
            ? 'bg-yellow-50 text-yellow-700 border-yellow-200'
            : 'bg-green-50 text-green-700 border-green-200'
        }`}>
          <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
              d="M9 3H5a2 2 0 00-2 2v4m6-6h10a2 2 0 012 2v4M9 3v18m0 0h10a2 2 0 002-2V9M9 21H5a2 2 0 01-2-2V9m0 0h18" />
          </svg>
          Analysis: {isMock ? 'mock (estimated)' : 'ProLIF'}
        </span>
      </div>
    </div>
  )
}
