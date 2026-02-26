import React from 'react'

const VARIANT_CLASSES = {
  green:  'bg-green-100 text-green-700',
  yellow: 'bg-yellow-100 text-yellow-700',
  red:    'bg-red-100 text-red-700',
  blue:   'bg-blue-100 text-blue-700',
  purple: 'bg-purple-100 text-purple-700',
  gray:   'bg-gray-100 text-gray-600',
}

const SIZE_CLASSES = {
  sm: 'px-2 py-0.5 text-xs',
  md: 'px-3 py-1 text-sm',
}

/**
 * Reusable badge component.
 * Props:
 *   variant — green | yellow | red | blue | purple | gray (default: gray)
 *   size    — sm | md (default: sm)
 *   children — badge label content
 */
export default function Badge({ variant = 'gray', size = 'sm', children }) {
  const variantClass = VARIANT_CLASSES[variant] ?? VARIANT_CLASSES.gray
  const sizeClass = SIZE_CLASSES[size] ?? SIZE_CLASSES.sm

  return (
    <span className={`inline-block rounded-full font-semibold leading-none ${sizeClass} ${variantClass}`}>
      {children}
    </span>
  )
}
