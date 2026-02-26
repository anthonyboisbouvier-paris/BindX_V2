/**
 * Toast.jsx — Re-exports the toast primitives from ToastContext for convenience.
 * Usage:
 *   import { useToast } from '../contexts/ToastContext'
 *   const { addToast } = useToast()
 *   addToast('Run queued successfully', 'success')
 *   addToast('An error occurred', 'error')
 *   addToast('Check your inputs', 'warning')
 *   addToast('Analysis started', 'info')
 *
 * The ToastProvider must wrap the component tree (e.g. in main.jsx or App.jsx).
 */
export { ToastProvider, useToast } from '../contexts/ToastContext.jsx'
