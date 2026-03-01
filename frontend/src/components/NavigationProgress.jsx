import { useEffect, useState, useRef } from 'react'
import { useLocation } from 'react-router-dom'

/**
 * NavigationProgress — thin gradient bar at the top of the viewport.
 * Animates 0 → 30 → 60 on location change, then snaps to 100% and fades out.
 */
export default function NavigationProgress() {
  const location = useLocation()
  const [progress, setProgress] = useState(0)
  const [visible, setVisible] = useState(false)
  const timerRef = useRef(null)

  useEffect(() => {
    // Start animation
    setVisible(true)
    setProgress(30)

    timerRef.current = setTimeout(() => setProgress(60), 150)

    return () => clearTimeout(timerRef.current)
  }, [location.pathname])

  // When the page finishes loading (Suspense resolved), snap to 100
  useEffect(() => {
    if (progress >= 60) {
      const t = setTimeout(() => {
        setProgress(100)
        setTimeout(() => {
          setVisible(false)
          setProgress(0)
        }, 300)
      }, 200)
      return () => clearTimeout(t)
    }
  }, [progress])

  if (!visible && progress === 0) return null

  return (
    <div
      className="fixed top-0 left-0 right-0 z-[9999] h-[2px] pointer-events-none"
      style={{ opacity: visible ? 1 : 0, transition: 'opacity 0.3s ease' }}
    >
      <div
        className="h-full bg-gradient-to-r from-bx-mint via-bx-cyan to-bx-mint"
        style={{
          width: `${progress}%`,
          transition: progress === 0 ? 'none' : 'width 0.4s ease',
        }}
      />
    </div>
  )
}
