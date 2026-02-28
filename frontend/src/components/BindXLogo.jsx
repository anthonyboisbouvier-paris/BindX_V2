import React, { useEffect, useRef, useCallback } from 'react'

// ---------------------------------------------------------------------------
// BindX Logo — Animated component with multiple variants
// Based on: frontend/src/assets/bindx-anims-final.html
//
// Variants:
//   splash    — Convergence orbitale (app launch, page load)
//   loading   — Canvas orbits with physics (computation, import)
//   screening — Particles flowing through logo (pipeline, screening)
//   success   — Merge + sparkles + checkmark (run completed)
//   error     — Red version + cross (failure, timeout)
//   phase     — Chromatic shift (phase transition A→B→C)
//   idle      — Static breathing logo (default)
//
// Props:
//   variant: string (default: 'idle')
//   size: number (default: 80) — width/height in px
//   className: string
//   label: string — optional text below the animation
// ---------------------------------------------------------------------------

const V = {
  bg: '#0f131d', s2: '#141925',
  mint: '#00e6a0', cyan: '#06b6d4', blue: '#3b82f6', violet: '#8b5cf6', red: '#f43f5e',
}

// ---- SVG base logo (3 circles + dot) ----
function LogoSVG({ size, r1Props, r2Props, r3Props, dotProps, children, className = '' }) {
  return (
    <svg
      className={className}
      width={size}
      height={size}
      viewBox="0 0 56 56"
      fill="none"
      style={{ overflow: 'visible' }}
    >
      {children}
      <g {...(r1Props || {})}><circle cx="18" cy="28" r="15" stroke={V.mint} strokeWidth="1.5" fill="none" /></g>
      <g {...(r2Props || {})}><circle cx="28" cy="28" r="12" stroke={V.cyan} strokeWidth="2" fill="none" /></g>
      <g {...(r3Props || {})}><circle cx="38" cy="28" r="9" stroke={V.blue} strokeWidth="2.5" fill="none" /></g>
      <circle cx="38" cy="28" r="4" fill={V.blue} {...(dotProps || {})} />
    </svg>
  )
}

// ============================================================================
// IDLE — Static with subtle breathing
// ============================================================================
function IdleLogo({ size }) {
  return (
    <div className="bindx-idle" style={{ width: size, height: size }}>
      <LogoSVG
        size={size}
        r1Props={{ style: { opacity: 0.45 } }}
        r2Props={{ style: { opacity: 0.7 } }}
        r3Props={{ style: { opacity: 1 } }}
        dotProps={{ opacity: 0.85 }}
      />
      <style>{`
        .bindx-idle svg g:first-of-type { animation: bx-breathe 4s ease-in-out infinite; }
        .bindx-idle svg g:nth-of-type(2) { animation: bx-breathe 4s ease-in-out infinite 0.3s; }
        .bindx-idle svg g:nth-of-type(3) { animation: bx-breathe 4s ease-in-out infinite 0.6s; }
        @keyframes bx-breathe { 0%,100% { transform: scale(1); } 50% { transform: scale(1.03); } }
      `}</style>
    </div>
  )
}

// ============================================================================
// SPLASH — Convergence orbitale
// ============================================================================
function SplashLogo({ size }) {
  return (
    <div className="bindx-splash" style={{ width: size, height: size, position: 'relative' }}>
      <svg width={size} height={size} viewBox="0 0 56 56" fill="none" style={{ overflow: 'visible' }}>
        <circle className="bx-sp-glow" cx="38" cy="28" fill={V.blue} r="8" opacity="0" />
        <g className="bx-sp-r1" style={{ transformOrigin: '18px 28px' }}>
          <circle cx="18" cy="28" r="15" stroke={V.mint} strokeWidth="1.5" fill="none" />
        </g>
        <g className="bx-sp-r2" style={{ transformOrigin: '28px 28px' }}>
          <circle cx="28" cy="28" r="12" stroke={V.cyan} strokeWidth="2" fill="none" />
        </g>
        <g className="bx-sp-r3" style={{ transformOrigin: '38px 28px' }}>
          <circle cx="38" cy="28" r="9" stroke={V.blue} strokeWidth="2.5" fill="none" />
        </g>
        <circle className="bx-sp-dot" cx="38" cy="28" fill={V.blue} r="0" />
      </svg>
      <style>{`
        .bx-sp-r1{animation:bx-sp1 1.6s cubic-bezier(.22,1,.36,1) forwards,bx-idle1 6s ease-in-out 2.2s infinite}
        .bx-sp-r2{animation:bx-sp2 1.4s cubic-bezier(.22,1,.36,1) .12s forwards,bx-idle2 6s ease-in-out 2.4s infinite}
        .bx-sp-r3{animation:bx-sp3 1.5s cubic-bezier(.22,1,.36,1) .22s forwards,bx-idle3 6s ease-in-out 2.6s infinite}
        .bx-sp-dot{animation:bx-spdot .45s cubic-bezier(.34,1.56,.64,1) 1.1s forwards}
        .bx-sp-glow{animation:bx-spglow .8s ease 1.1s forwards}
        @keyframes bx-sp1{0%{transform:translate(-70px,-25px) scale(1.25) rotate(-8deg);opacity:0}35%{opacity:.3}100%{transform:translate(0,0) scale(1) rotate(0);opacity:.45}}
        @keyframes bx-sp2{0%{transform:translate(5px,-55px) scale(.65) rotate(5deg);opacity:0}35%{opacity:.5}100%{transform:translate(0,0) scale(1) rotate(0);opacity:.7}}
        @keyframes bx-sp3{0%{transform:translate(55px,35px) scale(1.3) rotate(10deg);opacity:0}35%{opacity:.6}100%{transform:translate(0,0) scale(1) rotate(0);opacity:1}}
        @keyframes bx-spdot{0%{r:0;opacity:0}65%{r:5.5;opacity:1}100%{r:4;opacity:.85}}
        @keyframes bx-spglow{0%{r:4;opacity:0}40%{opacity:.2}100%{r:28;opacity:0}}
        @keyframes bx-idle1{0%,100%{transform:translate(0,0) rotate(0)}50%{transform:translate(-2.5px,-1.5px) rotate(-.5deg)}}
        @keyframes bx-idle2{0%,100%{transform:translate(0,0)}50%{transform:translate(.5px,1.5px)}}
        @keyframes bx-idle3{0%,100%{transform:translate(0,0) rotate(0)}50%{transform:translate(2px,-1px) rotate(.3deg)}}
      `}</style>
    </div>
  )
}

// ============================================================================
// LOADING — Canvas physics orbits
// ============================================================================
function LoadingLogo({ size }) {
  const canvasRef = useRef(null)
  const rafRef = useRef(null)

  useEffect(() => {
    const canvas = canvasRef.current
    if (!canvas) return
    const ctx = canvas.getContext('2d')
    const W = size, H = size
    const dpr = window.devicePixelRatio || 1
    canvas.width = W * dpr
    canvas.height = H * dpr
    canvas.style.width = W + 'px'
    canvas.style.height = H + 'px'
    ctx.scale(dpr, dpr)

    const anchor = { x: W / 2, y: H / 2 }
    const DOT_R = W * 0.025
    const dot = { x: anchor.x, y: anchor.y, vx: 0, vy: 0 }
    const SPRING = 0.03, DAMPING = 0.93, CENTRIFUGAL = 0.06

    const rings = [
      { radius: W * 0.225, strokeW: W * 0.008, color: [0, 230, 160], opacity: 0.45, angle: 0, mass: 3.2,
        speedFn: t => 0.4 + 0.8 * (0.5 + 0.5 * Math.sin(t * 0.8)) },
      { radius: W * 0.175, strokeW: W * 0.011, color: [6, 182, 212], opacity: 0.7, angle: Math.PI * 0.7, mass: 2.4,
        speedFn: t => -(0.5 + 1.0 * (0.5 + 0.5 * Math.sin(t * 0.6 + 1))) },
      { radius: W * 0.131, strokeW: W * 0.014, color: [59, 130, 246], opacity: 0.9, angle: Math.PI * 1.4, mass: 1.7,
        speedFn: t => 0.7 + 1.4 * (0.5 + 0.5 * Math.sin(t * 1.0 + 2)) },
    ]

    let lastT = performance.now(), elapsed = 0

    function tick(now) {
      const dt = Math.min((now - lastT) / 1000, 0.05)
      lastT = now
      elapsed += dt
      ctx.clearRect(0, 0, W, H)

      let fx = 0, fy = 0
      rings.forEach(r => {
        const spd = r.speedFn(elapsed)
        r.angle += spd * dt
        const orbitDist = r.radius - DOT_R
        const cx = dot.x + orbitDist * Math.cos(r.angle)
        const cy = dot.y + orbitDist * Math.sin(r.angle)
        const dirX = cx - dot.x, dirY = cy - dot.y
        const dist = Math.sqrt(dirX * dirX + dirY * dirY) || 1
        const force = r.mass * spd * spd * CENTRIFUGAL
        fx += (dirX / dist) * force
        fy += (dirY / dist) * force
        r._cx = cx
        r._cy = cy
      })

      fx += (anchor.x - dot.x) * SPRING
      fy += (anchor.y - dot.y) * SPRING
      dot.vx = (dot.vx + fx * dt) * DAMPING
      dot.vy = (dot.vy + fy * dt) * DAMPING
      dot.x += dot.vx
      dot.y += dot.vy

      rings.forEach(r => {
        const orbitDist = r.radius - DOT_R
        r._cx = dot.x + orbitDist * Math.cos(r.angle)
        r._cy = dot.y + orbitDist * Math.sin(r.angle)
      })

      // Glow
      const glR = W * 0.1 + Math.sin(now / 500) * W * 0.019
      const grd = ctx.createRadialGradient(dot.x, dot.y, 0, dot.x, dot.y, glR)
      grd.addColorStop(0, 'rgba(59,130,246,0.08)')
      grd.addColorStop(1, 'rgba(59,130,246,0)')
      ctx.fillStyle = grd
      ctx.beginPath()
      ctx.arc(dot.x, dot.y, glR, 0, Math.PI * 2)
      ctx.fill()

      // Rings
      rings.forEach(r => {
        ctx.beginPath()
        ctx.arc(r._cx, r._cy, r.radius, 0, Math.PI * 2)
        ctx.strokeStyle = `rgba(${r.color[0]},${r.color[1]},${r.color[2]},${r.opacity})`
        ctx.lineWidth = r.strokeW
        ctx.stroke()
      })

      // Dot
      ctx.beginPath()
      ctx.arc(dot.x, dot.y, DOT_R, 0, Math.PI * 2)
      ctx.fillStyle = V.blue
      ctx.fill()

      rafRef.current = requestAnimationFrame(tick)
    }
    rafRef.current = requestAnimationFrame(tick)

    return () => { if (rafRef.current) cancelAnimationFrame(rafRef.current) }
  }, [size])

  return (
    <div style={{ display: 'flex', flexDirection: 'column', alignItems: 'center', gap: size > 40 ? 8 : 0 }}>
      <canvas ref={canvasRef} style={{ display: 'block' }} />
      {size > 40 && (
        <div style={{ display: 'flex', gap: 6, alignItems: 'center' }}>
          <span className="bx-ldot" style={{ background: V.mint, animationDelay: '0s' }} />
          <span className="bx-ldot" style={{ background: V.cyan, animationDelay: '0.25s' }} />
          <span className="bx-ldot" style={{ background: V.blue, animationDelay: '0.5s' }} />
        </div>
      )}
      <style>{`
        .bx-ldot{width:5px;height:5px;border-radius:50%;animation:bx-ldpulse 1.8s ease-in-out infinite}
        @keyframes bx-ldpulse{0%,100%{opacity:.2;transform:scale(.7)}50%{opacity:1;transform:scale(1)}}
      `}</style>
    </div>
  )
}

// ============================================================================
// SUCCESS — Merge + sparkles + checkmark
// ============================================================================
function SuccessLogo({ size }) {
  const sparkles = [
    { sx: '-55px', sy: '-38px', cls: 'mint', d: 0 },
    { sx: '50px', sy: '-32px', cls: 'cyan', d: 0.04 },
    { sx: '58px', sy: '18px', cls: 'blue', d: 0.08 },
    { sx: '-45px', sy: '42px', cls: 'mint', d: 0.12 },
    { sx: '22px', sy: '-52px', cls: 'violet', d: 0.06 },
    { sx: '-58px', sy: '12px', cls: 'cyan', d: 0.10 },
    { sx: '38px', sy: '44px', cls: 'blue', d: 0.14 },
    { sx: '-18px', sy: '-50px', cls: 'mint', d: 0.03 },
  ]
  const colors = { mint: V.mint, cyan: V.cyan, blue: V.blue, violet: V.violet }
  const scale = size / 130

  return (
    <div className="bindx-ok" style={{ position: 'relative', display: 'flex', alignItems: 'center', justifyContent: 'center', width: size, height: size }}>
      {sparkles.map((s, i) => (
        <span key={i} className="bx-spark" style={{
          '--sx': `calc(${s.sx} * ${scale})`, '--sy': `calc(${s.sy} * ${scale})`,
          left: '50%', top: '50%', background: colors[s.cls],
          animation: `bx-spark 4s ease infinite ${s.d}s`,
          width: 3 * scale, height: 3 * scale,
        }} />
      ))}
      <svg width={size} height={size} viewBox="0 0 56 56" fill="none" style={{ overflow: 'visible' }}>
        <circle className="bx-ok-dbg" cx="38" cy="28" fill={V.mint} r="0" opacity="0" />
        <g className="bx-ok-r1" style={{ transformOrigin: '18px 28px' }}>
          <circle cx="18" cy="28" r="15" stroke={V.mint} strokeWidth="1.5" fill="none" opacity=".45" />
        </g>
        <g className="bx-ok-r2" style={{ transformOrigin: '28px 28px' }}>
          <circle cx="28" cy="28" r="12" stroke={V.cyan} strokeWidth="2" fill="none" opacity=".7" />
        </g>
        <g className="bx-ok-r3" style={{ transformOrigin: '38px 28px' }}>
          <circle cx="38" cy="28" r="9" stroke={V.blue} strokeWidth="2.5" fill="none" />
        </g>
        <circle className="bx-ok-dot" cx="38" cy="28" fill={V.mint} r="0" />
        <path className="bx-ok-chk" d="M34.5 28.5L37 31L42 25.5" stroke="#0f131d" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" fill="none" />
      </svg>
      <style>{`
        .bx-spark{position:absolute;border-radius:50%;pointer-events:none;opacity:0}
        @keyframes bx-spark{0%,52%{opacity:0;transform:translate(0,0) scale(0)}58%{opacity:.9;transform:translate(var(--sx),var(--sy)) scale(1.2)}76%{opacity:.2;transform:translate(calc(var(--sx)*1.6),calc(var(--sy)*1.6)) scale(.4)}90%,100%{opacity:0;transform:translate(calc(var(--sx)*2),calc(var(--sy)*2)) scale(0)}}
        .bx-ok-r1{animation:bx-okr1 4s ease infinite}.bx-ok-r2{animation:bx-okr2 4s ease infinite}.bx-ok-r3{animation:bx-okr3 4s ease infinite}
        .bx-ok-dbg{animation:bx-okdbg 4s ease infinite}.bx-ok-dot{animation:bx-okdot 4s ease infinite}
        .bx-ok-chk{stroke-dasharray:22;stroke-dashoffset:22;animation:bx-okchk 4s ease infinite}
        @keyframes bx-okr1{0%,28%{transform:translate(0,0);opacity:.45}48%{transform:translate(10px,0) scale(.82);opacity:.8}55%{transform:translate(12px,0) scale(.78);opacity:.9}68%{transform:translate(2px,0);opacity:.55}85%,100%{transform:translate(0,0);opacity:.45}}
        @keyframes bx-okr2{0%,28%{transform:translate(0,0);opacity:.7}48%{transform:translate(5px,0) scale(.88);opacity:.85}55%{transform:translate(6px,0) scale(.85);opacity:1}68%{transform:translate(1px,0);opacity:.75}85%,100%{transform:translate(0,0);opacity:.7}}
        @keyframes bx-okr3{0%,28%{transform:translate(0,0);opacity:1}48%{transform:translate(-1px,0) scale(.94)}55%{transform:translate(-1px,0) scale(.92)}68%{transform:translate(0,0)}85%,100%{transform:translate(0,0);opacity:1}}
        @keyframes bx-okdbg{0%,50%{r:0;opacity:0}57%{r:8;opacity:.2}66%{r:22;opacity:.04}78%,100%{r:28;opacity:0}}
        @keyframes bx-okdot{0%,50%{r:0;opacity:0}56%{r:6.5;opacity:1}64%{r:5.5;opacity:.92}100%{r:5.5;opacity:.92}}
        @keyframes bx-okchk{0%,58%{stroke-dashoffset:22;opacity:0}70%{stroke-dashoffset:0;opacity:1}100%{stroke-dashoffset:0;opacity:1}}
      `}</style>
    </div>
  )
}

// ============================================================================
// ERROR — Red version + cross
// ============================================================================
function ErrorLogo({ size }) {
  const sparkles = [
    { sx: '-50px', sy: '-35px', d: 0 },
    { sx: '48px', sy: '-28px', d: 0.05 },
    { sx: '52px', sy: '20px', d: 0.10 },
    { sx: '-42px', sy: '38px', d: 0.13 },
    { sx: '18px', sy: '-48px', d: 0.07 },
    { sx: '-52px', sy: '10px', d: 0.11 },
  ]
  const scale = size / 130

  return (
    <div className="bindx-err" style={{ position: 'relative', display: 'flex', alignItems: 'center', justifyContent: 'center', width: size, height: size }}>
      {sparkles.map((s, i) => (
        <span key={i} className="bx-spark-err" style={{
          '--sx': `calc(${s.sx} * ${scale})`, '--sy': `calc(${s.sy} * ${scale})`,
          left: '50%', top: '50%', background: i % 2 === 0 ? V.red : '#fb7185',
          animation: `bx-spark 4s ease infinite ${s.d}s`,
          width: 3 * scale, height: 3 * scale,
        }} />
      ))}
      <svg width={size} height={size} viewBox="0 0 56 56" fill="none" style={{ overflow: 'visible' }}>
        <circle className="bx-er-dbg" cx="38" cy="28" fill={V.red} r="0" opacity="0" />
        <g className="bx-er-r1" style={{ transformOrigin: '18px 28px' }}>
          <circle cx="18" cy="28" r="15" stroke={V.red} strokeWidth="1.5" fill="none" opacity=".45" />
        </g>
        <g className="bx-er-r2" style={{ transformOrigin: '28px 28px' }}>
          <circle cx="28" cy="28" r="12" stroke={V.red} strokeWidth="2" fill="none" opacity=".7" />
        </g>
        <g className="bx-er-r3" style={{ transformOrigin: '38px 28px' }}>
          <circle cx="38" cy="28" r="9" stroke={V.red} strokeWidth="2.5" fill="none" />
        </g>
        <circle className="bx-er-dot" cx="38" cy="28" fill={V.red} r="0" />
        <line className="bx-er-x1" x1="35.5" y1="25.5" x2="40.5" y2="30.5" stroke="#0f131d" strokeWidth="2" strokeLinecap="round" />
        <line className="bx-er-x2" x1="40.5" y1="25.5" x2="35.5" y2="30.5" stroke="#0f131d" strokeWidth="2" strokeLinecap="round" />
      </svg>
      <style>{`
        .bx-spark-err{position:absolute;border-radius:50%;pointer-events:none;opacity:0}
        @keyframes bx-spark{0%,52%{opacity:0;transform:translate(0,0) scale(0)}58%{opacity:.9;transform:translate(var(--sx),var(--sy)) scale(1.2)}76%{opacity:.2;transform:translate(calc(var(--sx)*1.6),calc(var(--sy)*1.6)) scale(.4)}90%,100%{opacity:0;transform:translate(calc(var(--sx)*2),calc(var(--sy)*2)) scale(0)}}
        .bx-er-r1{animation:bx-err1 4s ease infinite}.bx-er-r2{animation:bx-err2 4s ease infinite}.bx-er-r3{animation:bx-err3 4s ease infinite}
        .bx-er-dbg{animation:bx-erdbg 4s ease infinite}.bx-er-dot{animation:bx-erdot 4s ease infinite}
        .bx-er-x1,.bx-er-x2{stroke-dasharray:14;stroke-dashoffset:14;animation:bx-erx 4s ease infinite}
        .bx-er-x2{animation-delay:.06s}
        @keyframes bx-err1{0%,28%{transform:translate(0,0);opacity:.45}48%{transform:translate(10px,0) scale(.82);opacity:.8}55%{transform:translate(12px,0) scale(.78);opacity:.9}68%{transform:translate(2px,0);opacity:.55}85%,100%{transform:translate(0,0);opacity:.45}}
        @keyframes bx-err2{0%,28%{transform:translate(0,0);opacity:.7}48%{transform:translate(5px,0) scale(.88);opacity:.85}55%{transform:translate(6px,0) scale(.85);opacity:1}68%{transform:translate(1px,0);opacity:.75}85%,100%{transform:translate(0,0);opacity:.7}}
        @keyframes bx-err3{0%,28%{transform:translate(0,0);opacity:1}48%{transform:translate(-1px,0) scale(.94)}55%{transform:translate(-1px,0) scale(.92)}68%{transform:translate(0,0)}85%,100%{transform:translate(0,0);opacity:1}}
        @keyframes bx-erdbg{0%,50%{r:0;opacity:0}57%{r:8;opacity:.2}66%{r:22;opacity:.04}78%,100%{r:28;opacity:0}}
        @keyframes bx-erdot{0%,50%{r:0;opacity:0}56%{r:6.5;opacity:1}64%{r:5.5;opacity:.92}100%{r:5.5;opacity:.92}}
        @keyframes bx-erx{0%,58%{stroke-dashoffset:14;opacity:0}70%{stroke-dashoffset:0;opacity:1}100%{stroke-dashoffset:0;opacity:1}}
      `}</style>
    </div>
  )
}

// ============================================================================
// PHASE — Chromatic shift transition
// ============================================================================
function PhaseLogo({ size }) {
  return (
    <div className="bindx-phase" style={{ width: size, height: size }}>
      <svg width={size} height={size} viewBox="0 0 56 56" fill="none" style={{ overflow: 'visible' }}>
        <g className="bx-ph-r1" style={{ transformOrigin: '18px 28px' }}>
          <circle cx="18" cy="28" r="15" stroke={V.mint} strokeWidth="1.5" fill="none" opacity=".45" />
        </g>
        <g className="bx-ph-r2" style={{ transformOrigin: '28px 28px' }}>
          <circle cx="28" cy="28" r="12" stroke={V.cyan} strokeWidth="2" fill="none" opacity=".7" />
        </g>
        <g className="bx-ph-r3" style={{ transformOrigin: '38px 28px' }}>
          <circle cx="38" cy="28" r="9" stroke={V.blue} strokeWidth="2.5" fill="none" />
        </g>
        <circle className="bx-ph-dot" cx="38" cy="28" fill={V.blue} r="4" opacity=".85" />
      </svg>
      <style>{`
        .bx-ph-r1{animation:bx-phr1 6s ease-in-out infinite}.bx-ph-r2{animation:bx-phr2 6s ease-in-out infinite}
        .bx-ph-r3{animation:bx-phr3 6s ease-in-out infinite}.bx-ph-dot{animation:bx-phdot 6s ease-in-out infinite}
        @keyframes bx-phr1{0%,18%{stroke:${V.mint};opacity:.45;transform:translate(0,0)}28%{opacity:.2;transform:translate(-10px,0) scale(.85)}40%{opacity:0;transform:translate(-18px,0) scale(.65)}52%{opacity:0;transform:translate(18px,0) scale(.65);stroke:${V.cyan}}62%{opacity:.2;transform:translate(10px,0) scale(.85)}78%,100%{stroke:${V.cyan};opacity:.45;transform:translate(0,0)}}
        @keyframes bx-phr2{0%,18%{stroke:${V.cyan};opacity:.7}38%,62%{opacity:.45}78%,100%{stroke:${V.blue};opacity:.7}}
        @keyframes bx-phr3{0%,18%{stroke:${V.blue};opacity:1}38%,62%{opacity:.65}78%,100%{stroke:${V.violet};opacity:1}}
        @keyframes bx-phdot{0%,18%{fill:${V.blue};opacity:.85}38%{opacity:.25}62%{opacity:.25}78%,100%{fill:${V.violet};opacity:.85}}
      `}</style>
    </div>
  )
}

// ============================================================================
// STATIC — No animation, size-adaptive strokes
// ============================================================================
function StaticLogo({ size }) {
  // Adaptive strokes for legibility at small sizes (per design system reference)
  const sw1 = size <= 18 ? 3 : size <= 26 ? 2 : size <= 32 ? 1.8 : 1.5
  const sw2 = size <= 18 ? 3.5 : size <= 26 ? 2.5 : size <= 32 ? 2.2 : 2
  const sw3 = size <= 18 ? 4 : size <= 26 ? 3 : size <= 32 ? 2.8 : 2.5
  const dotR = size <= 18 ? 5.5 : size <= 26 ? 4.5 : 4
  const dotOp = size <= 26 ? 0.5 : 0.45
  return (
    <svg width={size} height={size} viewBox="0 0 56 56" fill="none" style={{ overflow: 'visible' }}>
      <circle cx="18" cy="28" r="15" stroke={V.mint} strokeWidth={sw1} fill="none" opacity={dotOp} />
      <circle cx="28" cy="28" r="12" stroke={V.cyan} strokeWidth={sw2} fill="none" opacity=".7" />
      <circle cx="38" cy="28" r="9" stroke={V.blue} strokeWidth={sw3} fill="none" />
      <circle cx="38" cy="28" r={dotR} fill={V.blue} opacity=".85" />
    </svg>
  )
}

// ============================================================================
// Main export
// ============================================================================
export default function BindXLogo({ variant = 'idle', size = 80, className = '', label }) {
  const Component = {
    static: StaticLogo,
    idle: IdleLogo,
    splash: SplashLogo,
    loading: LoadingLogo,
    success: SuccessLogo,
    error: ErrorLogo,
    phase: PhaseLogo,
  }[variant] || IdleLogo

  return (
    <div className={`flex flex-col items-center gap-2 ${className}`}>
      <Component size={size} />
      {label && (
        <p className="text-sm text-gray-400 animate-pulse">{label}</p>
      )}
    </div>
  )
}
