import { useEffect, useState, useRef, useCallback } from 'react'
import { useLocation } from 'react-router-dom'
import useWorkspaceStore from '../stores/workspaceStore.js'

/**
 * NavigationProgress — top loading bar visible during route changes and data loading.
 * Stays visible while Zustand store reports loading state (projects or molecules).
 * 3px gradient bar with glow effect for clear visual feedback.
 */
export default function NavigationProgress() {
  const location = useLocation()
  const [progress, setProgress] = useState(0)
  const [visible, setVisible] = useState(false)
  const timerRef = useRef(null)
  const trickleRef = useRef(null)
  const prevPath = useRef(location.pathname)

  // Subscribe to store loading states
  const storeLoading = useWorkspaceStore(s => s.loading)
  const molLoading = useWorkspaceStore(s => s.moleculesLoading)
  const isLoading = storeLoading || molLoading

  // Trickle: slowly increment progress while waiting
  const startTrickle = useCallback(() => {
    if (trickleRef.current) return
    trickleRef.current = setInterval(() => {
      setProgress(p => {
        if (p >= 90) return p
        // Slow down as we approach 90%
        const inc = p < 50 ? 3 : p < 70 ? 2 : 0.5
        return Math.min(p + inc, 90)
      })
    }, 300)
  }, [])

  const stopTrickle = useCallback(() => {
    if (trickleRef.current) {
      clearInterval(trickleRef.current)
      trickleRef.current = null
    }
  }, [])

  // Route change → start the bar
  useEffect(() => {
    if (location.pathname === prevPath.current) return
    prevPath.current = location.pathname

    clearTimeout(timerRef.current)
    stopTrickle()
    setVisible(true)
    setProgress(20)
    startTrickle()
  }, [location.pathname, startTrickle, stopTrickle])

  // Store loading state → keep bar alive or start it
  useEffect(() => {
    if (isLoading) {
      if (!visible) {
        setVisible(true)
        setProgress(15)
      }
      startTrickle()
    } else if (visible) {
      // Loading finished → snap to 100% and fade out
      stopTrickle()
      setProgress(100)
      timerRef.current = setTimeout(() => {
        setVisible(false)
        setProgress(0)
      }, 400)
    }
  }, [isLoading, visible, startTrickle, stopTrickle])

  // Cleanup
  useEffect(() => {
    return () => {
      clearTimeout(timerRef.current)
      stopTrickle()
    }
  }, [stopTrickle])

  if (!visible && progress === 0) return null

  return (
    <div
      className="fixed top-0 left-0 right-0 z-[9999] pointer-events-none"
      style={{ opacity: visible ? 1 : 0, transition: 'opacity 0.4s ease' }}
    >
      {/* Main bar */}
      <div
        className="h-[4px] bg-gradient-to-r from-emerald-400 via-cyan-400 to-emerald-400"
        style={{
          width: `${progress}%`,
          transition: progress === 0 ? 'none' : progress === 100 ? 'width 0.2s ease' : 'width 0.4s ease',
          boxShadow: '0 0 12px rgba(16,185,129,0.7), 0 0 4px rgba(16,185,129,0.4)',
        }}
      />
      {/* Pulsing glow underneath */}
      <div
        className="h-[3px] -mt-[2px] bg-gradient-to-r from-emerald-400/0 via-cyan-300/50 to-emerald-400/0 animate-pulse"
        style={{
          width: `${progress}%`,
          transition: progress === 0 ? 'none' : 'width 0.4s ease',
        }}
      />
    </div>
  )
}
