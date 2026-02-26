import React, { useMemo } from 'react'

// --------------------------------------------------
// Interaction type visual config
// --------------------------------------------------
const TYPE_CONFIG = {
  Hbond:       { label: 'H-Bond',       stroke: '#3b82f6', dash: '',      width: 2.5 },
  HBDonor:     { label: 'H-Bond',       stroke: '#3b82f6', dash: '',      width: 2.5 },
  HBAcceptor:  { label: 'H-Bond',       stroke: '#3b82f6', dash: '',      width: 2.5 },
  Hydrophobic: { label: 'Hydrophobic',  stroke: '#9ca3af', dash: '5,4',   width: 1.5 },
  PiStacking:  { label: 'Pi-Stacking',  stroke: '#8b5cf6', dash: '',      width: 2 },
  CationPi:    { label: 'Cation-Pi',    stroke: '#f97316', dash: '',      width: 2 },
  SaltBridge:  { label: 'Salt Bridge',  stroke: '#0d9488', dash: '3,2',   width: 2 },
  VDW:         { label: 'Van der Waals',stroke: '#d1d5db', dash: '2,3',   width: 1 },
}

function getTypeConfig(type) {
  return TYPE_CONFIG[type] || { label: type || 'Unknown', stroke: '#9ca3af', dash: '', width: 1.5 }
}

// --------------------------------------------------
// Build unique type list for legend (in order of priority)
// --------------------------------------------------
const TYPE_ORDER = ['Hbond','HBDonor','HBAcceptor','Hydrophobic','PiStacking','CationPi','SaltBridge','VDW']

// --------------------------------------------------
// SVG dimensions
// --------------------------------------------------
const CX     = 200
const CY     = 200
const RADIUS = 130   // residue positions
const VIEWBOX_SIZE = 400
const LIGAND_R = 28
const RES_R    = 20

// --------------------------------------------------
// Evenly distribute residues on a circle
// --------------------------------------------------
function residuePositions(count) {
  return Array.from({ length: count }, (_, i) => {
    const angle = (2 * Math.PI * i) / count - Math.PI / 2
    return {
      x: CX + RADIUS * Math.cos(angle),
      y: CY + RADIUS * Math.sin(angle),
      angle,
    }
  })
}

// --------------------------------------------------
// Shorten residue label if too long
// --------------------------------------------------
function shortLabel(residue) {
  if (!residue) return '?'
  // Remove chain prefix like "A:" if present
  const clean = residue.replace(/^[A-Z]:/, '').trim()
  return clean.length > 7 ? clean.slice(0, 7) : clean
}

// --------------------------------------------------
// Curved path from ligand center to residue
// --------------------------------------------------
function pathTo(rx, ry) {
  // Straight line slightly offset from center dot
  const dx = rx - CX
  const dy = ry - CY
  const len = Math.sqrt(dx * dx + dy * dy) || 1
  // Start on ligand circle edge
  const sx = CX + (LIGAND_R + 2) * (dx / len)
  const sy = CY + (LIGAND_R + 2) * (dy / len)
  // End on residue circle edge
  const ex = rx - (RES_R + 2) * (dx / len)
  const ey = ry - (RES_R + 2) * (dy / len)
  return `M ${sx} ${sy} L ${ex} ${ey}`
}

// --------------------------------------------------
// Legend item
// --------------------------------------------------
function LegendItem({ type, config }) {
  return (
    <div className="flex items-center gap-1.5">
      <svg width={22} height={10}>
        <line
          x1={0} y1={5} x2={22} y2={5}
          stroke={config.stroke}
          strokeWidth={config.width}
          strokeDasharray={config.dash}
          strokeLinecap="round"
        />
      </svg>
      <span className="text-[10px] text-gray-600">{config.label}</span>
    </div>
  )
}

// --------------------------------------------------
// Main InteractionDiagram component
// --------------------------------------------------
export default function InteractionDiagram({ interactions }) {
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

  // Normalise contacts: support both `interactions.interactions` (array) and
  // `interactions.residue_contacts` (dict) and `interactions.contacts` (array)
  const rawContacts = useMemo(() => {
    if (Array.isArray(interactions.contacts)) {
      return interactions.contacts
    }
    if (Array.isArray(interactions.interactions)) {
      return interactions.interactions.map((row) => ({
        residue: row.residue,
        type: row.type,
        functional: row.is_functional ?? row.functional ?? false,
      }))
    }
    if (interactions.residue_contacts && typeof interactions.residue_contacts === 'object') {
      return Object.entries(interactions.residue_contacts).map(([res, info]) => ({
        residue: res,
        type: typeof info === 'string' ? info : (info.type || 'Unknown'),
        functional: info.functional ?? info.is_functional ?? false,
      }))
    }
    return []
  }, [interactions])

  // De-duplicate residues, keeping highest-priority type per residue
  const contacts = useMemo(() => {
    const map = new Map()
    rawContacts.forEach((c) => {
      const existing = map.get(c.residue)
      if (!existing) {
        map.set(c.residue, c)
      } else {
        // Keep the type with higher priority in TYPE_ORDER
        const currIdx = TYPE_ORDER.indexOf(existing.type)
        const newIdx  = TYPE_ORDER.indexOf(c.type)
        if (newIdx < currIdx) map.set(c.residue, c)
        // Preserve functional flag
        if (c.functional) map.set(c.residue, { ...map.get(c.residue), functional: true })
      }
    })
    return Array.from(map.values())
  }, [rawContacts])

  const positions = useMemo(() => residuePositions(contacts.length), [contacts.length])

  // Unique interaction types present (for legend)
  const presentTypes = useMemo(() => {
    const seen = new Set()
    contacts.forEach((c) => seen.add(c.type))
    return TYPE_ORDER.filter((t) => seen.has(t)).concat(
      [...seen].filter((t) => !TYPE_ORDER.includes(t))
    )
  }, [contacts])

  const ligandName = interactions.ligand_name || 'Ligand'

  if (contacts.length === 0) {
    return (
      <p className="text-xs text-gray-400 text-center py-6">No interaction contacts to display</p>
    )
  }

  return (
    <div className="space-y-3">
      {/* SVG diagram */}
      <div className="bg-white rounded-xl border border-gray-100 overflow-hidden">
        <svg
          viewBox={`0 0 ${VIEWBOX_SIZE} ${VIEWBOX_SIZE}`}
          className="w-full max-w-xs mx-auto"
          style={{ display: 'block' }}
        >
          {/* Background */}
          <rect width={VIEWBOX_SIZE} height={VIEWBOX_SIZE} fill="#f9fafb" />

          {/* Title */}
          <text
            x={CX} y={18}
            textAnchor="middle"
            fontSize={10}
            fill="#6b7280"
            fontWeight={600}
          >
            Protein-Ligand Interaction Map
          </text>

          {/* Lines from ligand to residues */}
          {contacts.map((contact, idx) => {
            const pos = positions[idx]
            const cfg = getTypeConfig(contact.type)
            return (
              <path
                key={`line-${idx}`}
                d={pathTo(pos.x, pos.y)}
                stroke={cfg.stroke}
                strokeWidth={cfg.width}
                strokeDasharray={cfg.dash}
                strokeLinecap="round"
                fill="none"
                opacity={0.8}
              />
            )
          })}

          {/* Residue circles */}
          {contacts.map((contact, idx) => {
            const pos = positions[idx]
            const cfg = getTypeConfig(contact.type)
            const isFunctional = contact.functional === true
            const label = shortLabel(contact.residue)

            // Label position: push outward from center
            const dx = pos.x - CX
            const dy = pos.y - CY
            const len = Math.sqrt(dx * dx + dy * dy) || 1
            const lx = pos.x + (RES_R + 14) * (dx / len)
            const ly = pos.y + (RES_R + 14) * (dy / len)

            return (
              <g key={`res-${idx}`}>
                {/* Functional highlight ring */}
                {isFunctional && (
                  <circle
                    cx={pos.x} cy={pos.y}
                    r={RES_R + 5}
                    fill="none"
                    stroke="#22c55e"
                    strokeWidth={2}
                    strokeDasharray="3,2"
                    opacity={0.8}
                  />
                )}
                {/* Residue circle */}
                <circle
                  cx={pos.x} cy={pos.y}
                  r={RES_R}
                  fill="#fff"
                  stroke={cfg.stroke}
                  strokeWidth={2}
                />
                {/* Residue name inside circle */}
                <text
                  x={pos.x} y={pos.y + 3.5}
                  textAnchor="middle"
                  fontSize={label.length > 5 ? 7 : 8}
                  fontWeight={600}
                  fill={cfg.stroke}
                >
                  {label}
                </text>
                {/* Label outside the circle */}
                <text
                  x={lx} y={ly + 3}
                  textAnchor="middle"
                  fontSize={8}
                  fill="#374151"
                >
                  {contact.residue || ''}
                </text>
              </g>
            )
          })}

          {/* Central ligand node */}
          <circle
            cx={CX} cy={CY}
            r={LIGAND_R}
            fill="#1e3a5f"
            stroke="#fff"
            strokeWidth={2}
          />
          <text
            x={CX} y={CY - 4}
            textAnchor="middle"
            fontSize={8}
            fontWeight={700}
            fill="#fff"
          >
            {ligandName.length > 8 ? ligandName.slice(0, 8) : ligandName}
          </text>
          <text
            x={CX} y={CY + 8}
            textAnchor="middle"
            fontSize={7}
            fill="#93c5fd"
          >
            (ligand)
          </text>

          {/* Functional legend note at bottom */}
          {contacts.some((c) => c.functional) && (
            <g>
              <circle cx={12} cy={VIEWBOX_SIZE - 12} r={6} fill="none" stroke="#22c55e" strokeWidth={1.5} strokeDasharray="2,1.5" />
              <text x={22} y={VIEWBOX_SIZE - 8} fontSize={8} fill="#16a34a">
                = functional residue
              </text>
            </g>
          )}
        </svg>
      </div>

      {/* Legend */}
      {presentTypes.length > 0 && (
        <div className="flex flex-wrap gap-x-4 gap-y-1 px-1">
          {presentTypes.map((type) => (
            <LegendItem key={type} type={type} config={getTypeConfig(type)} />
          ))}
          {contacts.some((c) => c.functional) && (
            <div className="flex items-center gap-1.5">
              <svg width={14} height={14}>
                <circle cx={7} cy={7} r={6} fill="none" stroke="#22c55e" strokeWidth={1.5} strokeDasharray="2,1.5" />
              </svg>
              <span className="text-[10px] text-gray-600">Functional residue</span>
            </div>
          )}
        </div>
      )}
    </div>
  )
}
