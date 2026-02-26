import React, { createContext, useContext, useState, useCallback, useRef } from 'react'

const ToastContext = createContext(null)

let _id = 0
function nextId() { return ++_id }

/**
 * ToastProvider — wrap the app to enable useToast() anywhere.
 * Renders the Toast stack at bottom-right of the viewport.
 */
export function ToastProvider({ children }) {
  const [toasts, setToasts] = useState([])
  const timerRefs = useRef({})

  const dismiss = useCallback((id) => {
    setToasts(prev => prev.map(t => t.id === id ? { ...t, leaving: true } : t))
    // Remove after animation
    setTimeout(() => {
      setToasts(prev => prev.filter(t => t.id !== id))
      if (timerRefs.current[id]) {
        clearTimeout(timerRefs.current[id])
        delete timerRefs.current[id]
      }
    }, 300)
  }, [])

  const addToast = useCallback((message, type = 'info', duration = 4000) => {
    const id = nextId()
    setToasts(prev => {
      const next = [...prev, { id, message, type, leaving: false }]
      // Cap at 3 visible
      if (next.length > 3) {
        const overflow = next[0]
        setTimeout(() => dismiss(overflow.id), 0)
      }
      return next
    })
    if (duration > 0) {
      timerRefs.current[id] = setTimeout(() => dismiss(id), duration)
    }
    return id
  }, [dismiss])

  return (
    <ToastContext.Provider value={{ addToast, dismiss }}>
      {children}
      <ToastStack toasts={toasts} onDismiss={dismiss} />
    </ToastContext.Provider>
  )
}

export function useToast() {
  const ctx = useContext(ToastContext)
  if (!ctx) throw new Error('useToast must be used within ToastProvider')
  return ctx
}

// ---------------------------------------------------------------------------
// Internal stack renderer
// ---------------------------------------------------------------------------
const TYPE_CONFIG = {
  success: {
    border: 'border-l-green-500',
    icon: (
      <svg className="w-4 h-4 text-green-500 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M5 13l4 4L19 7" />
      </svg>
    ),
    label: 'text-green-700',
  },
  error: {
    border: 'border-l-red-500',
    icon: (
      <svg className="w-4 h-4 text-red-500 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2.5} d="M6 18L18 6M6 6l12 12" />
      </svg>
    ),
    label: 'text-red-700',
  },
  warning: {
    border: 'border-l-amber-500',
    icon: (
      <svg className="w-4 h-4 text-amber-500 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
          d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z" />
      </svg>
    ),
    label: 'text-amber-700',
  },
  info: {
    border: 'border-l-blue-500',
    icon: (
      <svg className="w-4 h-4 text-blue-500 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
          d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
      </svg>
    ),
    label: 'text-blue-700',
  },
}

function ToastStack({ toasts, onDismiss }) {
  if (toasts.length === 0) return null
  return (
    <div className="fixed bottom-5 right-5 z-[9999] flex flex-col gap-2 items-end pointer-events-none">
      {toasts.map(toast => (
        <ToastItem key={toast.id} toast={toast} onDismiss={onDismiss} />
      ))}
    </div>
  )
}

function ToastItem({ toast, onDismiss }) {
  const cfg = TYPE_CONFIG[toast.type] || TYPE_CONFIG.info

  return (
    <div
      className={`
        pointer-events-auto
        flex items-start gap-3 min-w-[280px] max-w-sm
        bg-white rounded-xl shadow-lg border border-gray-100 border-l-4
        ${cfg.border}
        px-4 py-3
        transition-all duration-300
        ${toast.leaving
          ? 'opacity-0 translate-x-4'
          : 'opacity-100 translate-x-0'
        }
      `}
      style={{
        animation: toast.leaving ? undefined : 'toast-slide-in 0.25s ease-out',
      }}
    >
      <div className="mt-0.5">{cfg.icon}</div>
      <p className={`flex-1 text-sm font-medium ${cfg.label}`}>{toast.message}</p>
      <button
        onClick={() => onDismiss(toast.id)}
        className="flex-shrink-0 p-0.5 rounded hover:bg-gray-100 text-gray-400 hover:text-gray-600 transition-colors"
        aria-label="Dismiss"
      >
        <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
        </svg>
      </button>
    </div>
  )
}
