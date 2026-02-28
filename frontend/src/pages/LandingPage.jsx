import React, { useEffect, useRef } from 'react'
import { Link } from 'react-router-dom'
import './landing.css'

/* ═══════════════════════════════════════════════
   SVG ICONS (reusable)
   ═══════════════════════════════════════════════ */
const CheckIcon = () => (
  <svg viewBox="0 0 14 14" fill="none"><path d="M3 7l3 3 5-5" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round" /></svg>
)

const LogoSvg = () => (
  <svg viewBox="0 0 26 26" fill="none">
    <defs><linearGradient id="nlg" x1="0" y1="0" x2="26" y2="26"><stop stopColor="#00e6a0" /><stop offset="1" stopColor="#06b6d4" /></linearGradient></defs>
    <path d="M8.5 5.5L13 3L17.5 5.5V10.5L13 13L8.5 10.5Z" stroke="url(#nlg)" strokeWidth="1.4" fill="none" />
    <path d="M8.5 15.5L13 13L17.5 15.5V20.5L13 23L8.5 20.5Z" stroke="url(#nlg)" strokeWidth="1.4" fill="none" />
    <circle cx="13" cy="13" r="2.2" fill="url(#nlg)" />
    <circle cx="8.5" cy="5.5" r="1.1" fill="url(#nlg)" opacity=".5" />
    <circle cx="17.5" cy="5.5" r="1.1" fill="url(#nlg)" opacity=".5" />
    <circle cx="8.5" cy="20.5" r="1.1" fill="url(#nlg)" opacity=".5" />
    <circle cx="17.5" cy="20.5" r="1.1" fill="url(#nlg)" opacity=".5" />
  </svg>
)

const FooterLogoSvg = () => (
  <svg viewBox="0 0 18 18" fill="none">
    <defs><linearGradient id="flg2" x1="0" y1="0" x2="18" y2="18"><stop stopColor="#00e6a0" /><stop offset="1" stopColor="#06b6d4" /></linearGradient></defs>
    <path d="M6 4L9 2.5L12 4V7.5L9 9L6 7.5Z" stroke="url(#flg2)" strokeWidth="1" fill="none" />
    <path d="M6 10.5L9 9L12 10.5V14L9 15.5L6 14Z" stroke="url(#flg2)" strokeWidth="1" fill="none" />
    <circle cx="9" cy="9" r="1.5" fill="url(#flg2)" />
  </svg>
)

/* ═══════════════════════════════════════════════
   MOLECULAR NETWORK CANVAS (background)
   ═══════════════════════════════════════════════ */
function useMolecularNetwork(canvasRef) {
  useEffect(() => {
    const c = canvasRef.current
    if (!c) return
    const ctx = c.getContext('2d')
    let W, H, animId
    const pts = []
    const mouse = { x: null, y: null }

    class P {
      constructor() {
        this.x = Math.random() * W
        this.y = Math.random() * H
        this.vx = (Math.random() - 0.5) * 0.35
        this.vy = (Math.random() - 0.5) * 0.35
        this.r = Math.random() * 2 + 0.6
        this.a = Math.random() * 0.35 + 0.1
        this.c = Math.random() > 0.6 ? '0,230,160' : '59,130,246'
      }
      update() {
        if (mouse.x !== null) {
          const dx = mouse.x - this.x, dy = mouse.y - this.y
          const d = Math.sqrt(dx * dx + dy * dy)
          if (d < 200) { const f = 0.0004 * (200 - d); this.vx += dx * f; this.vy += dy * f }
        }
        this.vx *= 0.995; this.vy *= 0.995
        this.x += this.vx; this.y += this.vy
        if (this.x < 0 || this.x > W) this.vx *= -1
        if (this.y < 0 || this.y > H) this.vy *= -1
      }
      draw() {
        ctx.beginPath()
        ctx.arc(this.x, this.y, this.r, 0, Math.PI * 2)
        ctx.fillStyle = `rgba(${this.c},${this.a})`
        ctx.fill()
      }
    }

    function resize() {
      W = c.width = window.innerWidth
      H = c.height = window.innerHeight
    }

    function loop() {
      ctx.clearRect(0, 0, W, H)
      pts.forEach(p => { p.update(); p.draw() })
      for (let i = 0; i < pts.length; i++) {
        for (let j = i + 1; j < pts.length; j++) {
          const dx = pts[i].x - pts[j].x, dy = pts[i].y - pts[j].y
          const d = Math.sqrt(dx * dx + dy * dy)
          if (d < 200) {
            ctx.beginPath()
            ctx.moveTo(pts[i].x, pts[i].y)
            ctx.lineTo(pts[j].x, pts[j].y)
            ctx.strokeStyle = `rgba(0,230,160,${0.07 * (1 - d / 200)})`
            ctx.lineWidth = 0.6
            ctx.stroke()
          }
        }
      }
      animId = requestAnimationFrame(loop)
    }

    function onMouseMove(e) { mouse.x = e.clientX; mouse.y = e.clientY }

    resize()
    window.addEventListener('resize', resize)
    window.addEventListener('mousemove', onMouseMove)
    for (let i = 0; i < 70; i++) pts.push(new P())
    loop()

    return () => {
      cancelAnimationFrame(animId)
      window.removeEventListener('resize', resize)
      window.removeEventListener('mousemove', onMouseMove)
    }
  }, [canvasRef])
}

/* ═══════════════════════════════════════════════
   FOOTER PARTICLE CANVAS
   ═══════════════════════════════════════════════ */
function useFooterParticles(canvasRef, containerRef) {
  useEffect(() => {
    const fc = canvasRef.current
    const ft = containerRef.current
    if (!fc || !ft) return
    const ctx = fc.getContext('2d')
    let fw, fh, animId
    const fp = []

    function resize() {
      fw = fc.width = ft.offsetWidth
      fh = fc.height = ft.offsetHeight
    }

    resize()
    window.addEventListener('resize', resize)
    for (let i = 0; i < 45; i++) {
      fp.push({ x: Math.random() * fw, y: Math.random() * fh, vx: (Math.random() - 0.5) * 0.2, vy: (Math.random() - 0.5) * 0.2, r: Math.random() * 1.5 + 0.4 })
    }

    function loop() {
      ctx.clearRect(0, 0, fw, fh)
      fp.forEach(p => {
        p.x += p.vx; p.y += p.vy
        if (p.x < 0 || p.x > fw) p.vx *= -1
        if (p.y < 0 || p.y > fh) p.vy *= -1
        ctx.beginPath()
        ctx.arc(p.x, p.y, p.r, 0, Math.PI * 2)
        ctx.fillStyle = 'rgba(0,230,160,.12)'
        ctx.fill()
      })
      for (let i = 0; i < fp.length; i++) {
        for (let j = i + 1; j < fp.length; j++) {
          const dx = fp[i].x - fp[j].x, dy = fp[i].y - fp[j].y
          const d = Math.sqrt(dx * dx + dy * dy)
          if (d < 140) {
            ctx.beginPath()
            ctx.moveTo(fp[i].x, fp[i].y)
            ctx.lineTo(fp[j].x, fp[j].y)
            ctx.strokeStyle = `rgba(0,230,160,${0.1 * (1 - d / 140)})`
            ctx.lineWidth = 0.6
            ctx.stroke()
          }
        }
      }
      animId = requestAnimationFrame(loop)
    }

    loop()
    return () => {
      cancelAnimationFrame(animId)
      window.removeEventListener('resize', resize)
    }
  }, [canvasRef, containerRef])
}

/* ═══════════════════════════════════════════════
   SCROLL REVEAL HOOK
   ═══════════════════════════════════════════════ */
function useScrollReveal(rootRef) {
  useEffect(() => {
    const root = rootRef.current
    if (!root) return
    const obs = new IntersectionObserver(entries => {
      entries.forEach(e => {
        if (e.isIntersecting) { e.target.classList.add('vis'); obs.unobserve(e.target) }
      })
    }, { threshold: 0.1 })

    root.querySelectorAll('.phase, .fh-row, .ba-wrap, .arch-tree').forEach(el => obs.observe(el))
    return () => obs.disconnect()
  }, [rootRef])
}

/* ═══════════════════════════════════════════════
   LANDING PAGE COMPONENT
   ═══════════════════════════════════════════════ */
export default function LandingPage() {
  const molCanvasRef = useRef(null)
  const footerCanvasRef = useRef(null)
  const footerRef = useRef(null)
  const rootRef = useRef(null)

  useMolecularNetwork(molCanvasRef)
  useFooterParticles(footerCanvasRef, footerRef)
  useScrollReveal(rootRef)

  return (
    <div className="landing-root" ref={rootRef}>
      {/* Background molecular network */}
      <canvas ref={molCanvasRef} className="mol-bg" />

      {/* ═══ NAV ═══ */}
      <nav className="lnav">
        <a className="logo" href="#top">
          <LogoSvg />
          <span className="logo-t">Bind<span style={{ color: 'var(--mint)' }}>X</span></span>
        </a>
        <div className="nr">
          <a href="#how">How it works</a>
          <a href="#compute">Computations</a>
          <a href="#libraries">Libraries</a>
          <a href="#pricing">Pricing</a>
          <a href="#early-access">Early Access</a>
          <Link to="/login" className="btn btn-o">Sign in</Link>
          <Link to="/register" className="btn btn-m">Start free →</Link>
        </div>
      </nav>

      {/* ═══ HERO ═══ */}
      <section className="hero">
        <div className="hero-tag"><span>Free at launch · In Silico Drug Discovery</span></div>
        <h1>The complete in silico pipeline.<br /><span className="a">In one place.</span></h1>
        <p className="hero-p">Docking, ADMET, safety, retrosynthesis, AI generation, and multi-objective scoring — unified in a <strong>single dashboard</strong>. Paste a target. Run your computations. Get candidates. No scripting. No file juggling.</p>
        <p className="hero-p2">From target identification to synthesis-ready candidates — one workflow, one dashboard, zero scripting. And it's free.</p>
        <div className="hero-btns">
          <Link to="/register" className="btn btn-m" style={{ padding: '.6rem 1.8rem', fontSize: '.85rem' }}>Start a project — it's free</Link>
          <a href="#how" className="btn btn-o" style={{ padding: '.6rem 1.8rem', fontSize: '.85rem' }}>See how it works ↓</a>
        </div>
        <div className="hstats">
          <div className="hs"><div className="n">1M+ <span className="u">molecules</span></div><div className="l">screening capacity</div></div>
          <div className="hs"><div className="n">9 <span className="u">types</span></div><div className="l">of computation</div></div>
          <div className="hs"><div className="n">4 <span className="u">engines</span></div><div className="l">for docking</div></div>
          <div className="hs"><div className="n">37 <span className="u">properties</span></div><div className="l">per molecule</div></div>
        </div>
      </section>

      {/* ═══ BEFORE / AFTER ═══ */}
      <section className="sec sec-alt" id="why">
        <div className="inner">
          <div className="center" style={{ marginBottom: '3.5rem' }}>
            <div className="ey">Why BindX</div>
            <h2 className="sh">The in silico workflow, <span className="a">simplified.</span></h2>
            <p className="ss">What used to require 6 tools, 4 file formats, and days of scripting — now runs in one place.</p>
          </div>
          <div className="ba-wrap">
            <div className="ba-col ba-before">
              <div className="ba-badge" style={{ background: 'rgba(244,63,94,.08)', color: 'var(--rose)', borderColor: 'rgba(244,63,94,.15)' }}>Without BindX</div>
              <div className="ba-items">
                <div className="ba-item"><div className="ba-dot" style={{ background: 'var(--rose)' }} /><div className="ba-txt"><strong>Fragmented tools</strong><span>AutoDock for docking · ADMETlab for ADMET · RDKit for descriptors · custom scripts for scoring</span></div></div>
                <div className="ba-item"><div className="ba-dot" style={{ background: 'var(--rose)' }} /><div className="ba-txt"><strong>Constant file conversion</strong><span>SDF → SMILES → CSV → SDF at every step of the pipeline</span></div></div>
                <div className="ba-item"><div className="ba-dot" style={{ background: 'var(--rose)' }} /><div className="ba-txt"><strong>Missing capabilities</strong><span>Safety profiling, retrosynthesis, AI generation? Separate expensive products — or not available at all</span></div></div>
                <div className="ba-item"><div className="ba-dot" style={{ background: 'var(--rose)' }} /><div className="ba-txt"><strong>Enterprise pricing</strong><span>$7,500/yr academic · $50,000+ commercial license</span></div></div>
              </div>
            </div>
            <div className="ba-divider">
              <svg viewBox="0 0 24 60" fill="none" width="24"><line x1="12" y1="0" x2="12" y2="60" stroke="var(--border)" strokeWidth="1" /><polygon points="6,28 18,28 12,36" fill="var(--mint)" opacity=".6" /></svg>
            </div>
            <div className="ba-col ba-after">
              <div className="ba-badge" style={{ background: 'var(--mint-d)', color: 'var(--mint)', borderColor: 'rgba(0,230,160,.15)' }}>With BindX</div>
              <div className="ba-items">
                <div className="ba-item"><div className="ba-dot" style={{ background: 'var(--mint)' }} /><div className="ba-txt"><strong>One unified platform</strong><span>Docking, ADMET, safety, retrosynthesis, AI generation, scoring — all in one workflow</span></div></div>
                <div className="ba-item"><div className="ba-dot" style={{ background: 'var(--mint)' }} /><div className="ba-txt"><strong>One cumulative dashboard</strong><span>Every run enriches the same table. No files. No imports. Nothing gets lost</span></div></div>
                <div className="ba-item"><div className="ba-dot" style={{ background: 'var(--mint)' }} /><div className="ba-txt"><strong>No scripting required</strong><span>Select computations with checkboxes. Click Run. AI agents handle the configuration</span></div></div>
                <div className="ba-item"><div className="ba-dot" style={{ background: 'var(--mint)' }} /><div className="ba-txt"><strong>Free at launch</strong><span>Full platform access · GPU docking: €1 / 100 molecules</span></div></div>
              </div>
            </div>
          </div>
        </div>
      </section>

      {/* ═══ HOW IT WORKS ═══ */}
      <section className="sec" id="how">
        <div className="inner">
          <div className="center" style={{ marginBottom: '1rem' }}>
            <div className="ey">How It Works</div>
            <h2 className="sh">Five clicks from target to <span className="a">candidates.</span></h2>
            <p className="ss">Every project follows the same logical flow. No manual configuration needed.</p>
          </div>
          <div className="wf-steps">
            {[
              { color: 'var(--mint)',   n: 'Step 1', t: 'Define target',     d: 'UniProt ID, PDB ID, or FASTA. Auto-resolve structure & pockets.' },
              { color: 'var(--cyan)',   n: 'Step 2', t: 'Import compounds',  d: 'Public databases, SDF, or SMILES. Up to 1M+ molecules per library.' },
              { color: 'var(--blue)',   n: 'Step 3', t: 'Run calculations',  d: 'Docking, ADMET, scoring, confidence. Select with checkboxes.' },
              { color: 'var(--violet)', n: 'Step 4', t: 'Optimize leads',    d: 'AI generation, safety profiling, retrosynthesis, Pareto optimization.' },
              { color: 'var(--amber)',  n: 'Step 5', t: 'Export report',     d: 'PDF, CSV, or SDF. Full campaign data in one click.' },
            ].map(s => (
              <div className="wf-step" key={s.n}>
                <div className="wf-dot" style={{ background: s.color }} />
                <div className="wf-num">{s.n}</div>
                <div className="wf-t">{s.t}</div>
                <div className="wf-d">{s.d}</div>
              </div>
            ))}
          </div>
        </div>
      </section>

      {/* ═══ ARCHITECTURE ═══ */}
      <section className="sec sec-alt">
        <div className="inner">
          <div className="arch-row">
            <div className="arch-text">
              <div className="ey">Architecture V9</div>
              <h3>Modular by <span className="a">design.</span></h3>
              <p>BindX replaces the monolithic pipeline with a <strong>4-level hierarchy</strong> that maps to how computational chemists actually work: define a target, select a pocket, run phases, iterate.</p>
              <p>The <strong>Phase Dashboard</strong> is the central workspace. Every run — import, calculation, or generation — adds columns (properties) or rows (compounds) to a single, cumulative table. Nothing gets overwritten. Everything is traceable.</p>
              <p>Freeze a phase when you're satisfied. Advance to the next. Go back if you need to. The whole campaign history is preserved.</p>
            </div>
            <div className="arch-tree">
              <div className="at-line"><span className="at-icon">📁</span><span className="at-key">Project</span> <span className="at-dim">—</span> <span className="at-val">EGFR Inhibitors</span></div>
              <div className="at-line at-indent"><span className="at-icon">🎯</span><span className="at-val at-dim">Target: P00533 — PDB 1M17 — 5 pockets — 3400 ChEMBL actives</span></div>
              <div className="at-line at-indent"><span className="at-icon">🧫</span><span className="at-key" style={{ color: 'var(--cyan)' }}>Campaign</span> <span className="at-dim">—</span> <span className="at-val">ATP Pocket (druggability 0.92)</span></div>
              <div className="at-line at-indent2"><span className="at-icon">🔬</span><span className="at-key" style={{ color: 'var(--mint)' }}>Phase A</span> <span className="at-dim">—</span> <span className="at-val">Hit Discovery</span> <span className="at-tag at-tag-active">active</span></div>
              <div className="at-line at-indent2" style={{ paddingLeft: '4.8rem' }}><span className="at-dim">▸</span> <span className="at-key" style={{ color: 'var(--amber)' }}>Run Import</span> <span className="at-dim">— 20,000 mols from SDF</span></div>
              <div className="at-line at-indent2" style={{ paddingLeft: '4.8rem' }}><span className="at-dim">▸</span> <span className="at-key" style={{ color: 'var(--amber)' }}>Run Calcul</span> <span className="at-dim">— Docking GNINA GPU + ADMET + Confidence</span></div>
              <div className="at-line at-indent2" style={{ paddingLeft: '4.8rem' }}><span className="at-dim">▸</span> <span className="at-key" style={{ color: 'var(--amber)' }}>Run Calcul</span> <span className="at-dim">— Scoring (custom weights)</span></div>
              <div className="at-sep" />
              <div className="at-line at-indent2"><span className="at-icon">🧪</span><span className="at-key" style={{ color: 'var(--blue)' }}>Phase B</span> <span className="at-dim">—</span> <span className="at-val">Hit-to-Lead</span> <span className="at-tag at-tag-frozen">frozen</span></div>
              <div className="at-line at-indent2" style={{ paddingLeft: '4.8rem' }}><span className="at-dim">▸</span> <span className="at-key" style={{ color: 'var(--amber)' }}>Run Generation</span> <span className="at-dim">— 3 iter × 50 mols (batch)</span></div>
              <div className="at-sep" />
              <div className="at-line at-indent2"><span className="at-icon">💊</span><span className="at-key" style={{ color: 'var(--violet)' }}>Phase C</span> <span className="at-dim">—</span> <span className="at-val">Lead Optimization</span> <span className="at-tag at-tag-done">completed</span></div>
              <div className="at-line at-indent2" style={{ paddingLeft: '4.8rem' }}><span className="at-dim">▸</span> <span className="at-key" style={{ color: 'var(--amber)' }}>Run Calcul</span> <span className="at-dim">— Off-target + Retrosynthesis + Safety</span></div>
            </div>
          </div>
        </div>
      </section>

      {/* ═══ PHASES ═══ */}
      <section className="sec">
        <div className="inner">
          <div className="center" style={{ marginBottom: '3rem' }}>
            <div className="ey">Discovery Funnel</div>
            <h2 className="sh">Screen → Generate → <span className="a">Optimize.</span></h2>
            <p className="ss">Three phases that mirror how drug discovery actually works — from broad screening to synthesis-ready candidates.</p>
          </div>
          <div className="phases">

            {/* Phase A */}
            <div className="phase">
              <div className="tline tline-m" />
              <div className="ph-label"><span className="d" style={{ background: 'var(--mint)' }} /> Phase A</div>
              <div className="ph-icon"><svg viewBox="0 0 34 34" fill="none"><circle cx="17" cy="17" r="14" stroke="var(--mint)" strokeWidth="1.2" opacity=".25" /><circle cx="17" cy="17" r="8" stroke="var(--mint)" strokeWidth="1.2" opacity=".45" /><circle cx="17" cy="17" r="2.5" fill="var(--mint)" /></svg></div>
              <div className="ph-title">Hit Discovery</div>
              <div className="ph-desc">Screen a compound library — from 10,000 to <strong style={{ color: 'var(--text)' }}>1,000,000+ molecules</strong> — against your target pocket. Dock, compute ADMET, score, and filter to your top hits.</div>
              <div className="ph-detail">
                <strong>Input:</strong> Public libraries (ZINC, Enamine REAL, ChEMBL), custom SDF/SMILES, or curated subsets<br />
                <strong>Output:</strong> Ranked dashboard with 15+ properties per molecule<br />
                <strong>Then:</strong> Bookmark top 50–200 hits → freeze → advance
              </div>
              <div className="ph-metrics">
                <div className="pm"><div className="v">1M+</div><div className="l">screened</div></div>
                <div className="pm"><div className="v">~200</div><div className="l">hits</div></div>
                <div className="pm"><div className="v">15+</div><div className="l">props/mol</div></div>
              </div>
            </div>

            {/* Phase B */}
            <div className="phase">
              <div className="tline tline-b" />
              <div className="ph-label"><span className="d" style={{ background: 'var(--blue)' }} /> Phase B</div>
              <div className="ph-icon"><svg viewBox="0 0 34 34" fill="none"><path d="M9 25L17 7L25 25" stroke="var(--blue)" strokeWidth="1.2" opacity=".4" /><circle cx="17" cy="7" r="2.5" fill="var(--blue)" /><circle cx="9" cy="25" r="2" fill="var(--blue)" opacity=".5" /><circle cx="25" cy="25" r="2" fill="var(--blue)" opacity=".5" /></svg></div>
              <div className="ph-title">Hit-to-Lead</div>
              <div className="ph-desc">Generate novel analogs with AI — batch scaffold hopping or per-molecule R-group editing. Every generated compound is auto-scored and parent-linked.</div>
              <div className="ph-detail">
                <strong>Input:</strong> 50–200 frozen hits from Phase A<br />
                <strong>Output:</strong> 500–2,000 novel analogs, all pre-scored<br />
                <strong>Then:</strong> Pareto multi-objective → select ~30 leads
              </div>
              <div className="ph-metrics">
                <div className="pm"><div className="v">2k</div><div className="l">generated</div></div>
                <div className="pm"><div className="v">~30</div><div className="l">leads</div></div>
                <div className="pm"><div className="v">25+</div><div className="l">props/mol</div></div>
              </div>
            </div>

            {/* Phase C */}
            <div className="phase">
              <div className="tline tline-v" />
              <div className="ph-label"><span className="d" style={{ background: 'var(--violet)' }} /> Phase C</div>
              <div className="ph-icon"><svg viewBox="0 0 34 34" fill="none"><rect x="8" y="10" width="18" height="16" rx="3" stroke="var(--violet)" strokeWidth="1.2" opacity=".4" /><circle cx="14" cy="19" r="2" fill="var(--violet)" opacity=".5" /><circle cx="20" cy="19" r="2" fill="var(--violet)" opacity=".5" /><circle cx="17" cy="24" r="2" fill="var(--violet)" opacity=".5" /><path d="M13 10V7h8v3" stroke="var(--violet)" strokeWidth="1.2" strokeLinecap="round" opacity=".4" /></svg></div>
              <div className="ph-title">Lead Optimization</div>
              <div className="ph-desc">Full safety profiling, synthesis feasibility, off-target selectivity. This is where you decide what's worth synthesizing.</div>
              <div className="ph-detail">
                <strong>Input:</strong> 15–50 leads from Phase B<br />
                <strong>Output:</strong> Complete safety + retrosynthesis + selectivity profiles<br />
                <strong>Then:</strong> 3–10 candidates ready for synthesis → PDF report
              </div>
              <div className="ph-metrics">
                <div className="pm"><div className="v">37</div><div className="l">total props</div></div>
                <div className="pm"><div className="v">3–10</div><div className="l">candidates</div></div>
                <div className="pm"><div className="v">PDF</div><div className="l">report</div></div>
              </div>
            </div>

          </div>
        </div>
      </section>

      {/* ═══ COMPUTATIONS ═══ */}
      <section className="sec" id="compute">
        <div className="inner">
          <div className="center" style={{ marginBottom: '3rem' }}>
            <div className="ey">Computations</div>
            <h2 className="sh">Every analysis you need. <span className="a">One checkbox each.</span></h2>
            <p className="ss">Select computations, click Run. BindX handles execution order, resource allocation, and data merging automatically.</p>
          </div>
          <div className="comp-grid">
            {/* Molecular Docking */}
            <div className="cc cm">
              <div className="cc-top">
                <span className="cc-name"><svg viewBox="0 0 15 15" fill="none"><circle cx="7.5" cy="7.5" r="5.5" stroke="var(--mint)" strokeWidth="1" /><circle cx="7.5" cy="7.5" r="2" fill="var(--mint)" /></svg> Molecular Docking</span>
                <span className="cc-count cnt-m">4 outputs</span>
              </div>
              <div className="cc-desc">4 engine choices. Dock compounds into target binding pocket.</div>
              <div className="cc-outs"><span className="cc-out">docking_score</span><span className="cc-out">CNNscore</span><span className="cc-out">CNNaffinity</span><span className="cc-out">poses</span></div>
            </div>
            {/* ADMET */}
            <div className="cc cm">
              <div className="cc-top">
                <span className="cc-name"><svg viewBox="0 0 15 15" fill="none"><rect x="2" y="2" width="11" height="11" rx="2" stroke="var(--mint)" strokeWidth="1" /></svg> ADMET</span>
                <span className="cc-count cnt-m">8 outputs</span>
              </div>
              <div className="cc-desc">Full pharmacokinetic profile. Absorption to toxicity.</div>
              <div className="cc-outs"><span className="cc-out">logP</span><span className="cc-out">solubility</span><span className="cc-out">BBB</span><span className="cc-out">hERG</span><span className="cc-out">metab_stab</span><span className="cc-out">oral_bioavail</span><span className="cc-out">PPB</span><span className="cc-out">CYP</span></div>
            </div>
            {/* Scoring */}
            <div className="cc cb">
              <div className="cc-top">
                <span className="cc-name"><svg viewBox="0 0 15 15" fill="none"><polyline points="2,12 5,6 8,9 11,3 14,7" stroke="var(--blue)" strokeWidth="1" fill="none" strokeLinecap="round" /></svg> Scoring</span>
                <span className="cc-count cnt-b">1 output</span>
              </div>
              <div className="cc-desc">Weighted composite with visual editor.</div>
              <div className="cc-outs"><span className="cc-out">composite_score</span></div>
            </div>
            {/* Enrichment */}
            <div className="cc cb">
              <div className="cc-top">
                <span className="cc-name"><svg viewBox="0 0 15 15" fill="none"><circle cx="5" cy="10" r="2.5" stroke="var(--blue)" strokeWidth="1" /><circle cx="10" cy="5" r="2.5" stroke="var(--blue)" strokeWidth="1" /></svg> Enrichment</span>
                <span className="cc-count cnt-b">4 outputs</span>
              </div>
              <div className="cc-desc">ProLIF interaction fingerprints + clustering.</div>
              <div className="cc-outs"><span className="cc-out">interactions</span><span className="cc-out">inter_count</span><span className="cc-out">cluster_id</span><span className="cc-out">scaffold</span></div>
            </div>
            {/* Clustering */}
            <div className="cc cb">
              <div className="cc-top">
                <span className="cc-name"><svg viewBox="0 0 15 15" fill="none"><circle cx="4" cy="11" r="1.5" stroke="var(--blue)" strokeWidth="1" /><circle cx="11" cy="4" r="1.5" stroke="var(--blue)" strokeWidth="1" /><circle cx="11" cy="11" r="1.5" stroke="var(--blue)" strokeWidth="1" /></svg> Clustering</span>
                <span className="cc-count cnt-b">3 outputs</span>
              </div>
              <div className="cc-desc">Diversity, scaffolds, Tanimoto similarity.</div>
              <div className="cc-outs"><span className="cc-out">cluster_id</span><span className="cc-out">scaffold</span><span className="cc-out">tanimoto</span></div>
            </div>
            {/* Off-target */}
            <div className="cc ca">
              <div className="cc-top">
                <span className="cc-name"><svg viewBox="0 0 15 15" fill="none"><circle cx="7.5" cy="7.5" r="5.5" stroke="var(--amber)" strokeWidth="1" /><circle cx="4.5" cy="7.5" r="1.2" fill="var(--amber)" opacity=".4" /><circle cx="10.5" cy="7.5" r="1.2" fill="var(--amber)" opacity=".4" /></svg> Off-target</span>
                <span className="cc-count cnt-a">3 outputs</span>
              </div>
              <div className="cc-desc">Selectivity profiling against related targets.</div>
              <div className="cc-outs"><span className="cc-out">selectivity</span><span className="cc-out">off_target_hits</span><span className="cc-out">ratio</span></div>
            </div>
            {/* Safety */}
            <div className="cc cr">
              <div className="cc-top">
                <span className="cc-name"><svg viewBox="0 0 15 15" fill="none"><path d="M7.5 2L9.5 6H5.5L7.5 2Z" stroke="var(--rose)" strokeWidth="1" fill="none" /><circle cx="7.5" cy="11" r="1" fill="var(--rose)" /><line x1="7.5" y1="8" x2="7.5" y2="10" stroke="var(--rose)" strokeWidth="1" strokeLinecap="round" /></svg> Safety</span>
                <span className="cc-count cnt-r">6 outputs</span>
              </div>
              <div className="cc-desc">Full toxicity: hERG, AMES, hepatotox, carcinogenicity.</div>
              <div className="cc-outs"><span className="cc-out">hERG</span><span className="cc-out">ames</span><span className="cc-out">hepatotox</span><span className="cc-out">skin_sens</span><span className="cc-out">carcino</span><span className="cc-out">color_code</span></div>
            </div>
            {/* Retrosynthesis */}
            <div className="cc cv">
              <div className="cc-top">
                <span className="cc-name"><svg viewBox="0 0 15 15" fill="none"><path d="M7.5 2v3M7.5 5L4 9M7.5 5L11 9M4 9v3M11 9v3" stroke="var(--violet)" strokeWidth="1" strokeLinecap="round" /></svg> Retrosynthesis</span>
                <span className="cc-count cnt-v">4 outputs</span>
              </div>
              <div className="cc-desc">Route, cost, reagent availability.</div>
              <div className="cc-outs"><span className="cc-out">n_steps</span><span className="cc-out">confidence</span><span className="cc-out">cost</span><span className="cc-out">reagents</span></div>
            </div>
            {/* Confidence */}
            <div className="cc cb">
              <div className="cc-top">
                <span className="cc-name"><svg viewBox="0 0 15 15" fill="none"><path d="M7.5 2v9M4.5 5l3-3 3 3" stroke="var(--blue)" strokeWidth="1" strokeLinecap="round" strokeLinejoin="round" /></svg> Confidence</span>
                <span className="cc-count cnt-b">4 outputs</span>
              </div>
              <div className="cc-desc">PAINS, applicability domain, convergence.</div>
              <div className="cc-outs"><span className="cc-out">confidence</span><span className="cc-out">pains</span><span className="cc-out">AD</span><span className="cc-out">flags</span></div>
            </div>
          </div>
        </div>
      </section>

      {/* ═══ ENGINES ═══ */}
      <section className="sec sec-alt" id="engines">
        <div className="inner">
          <div className="center" style={{ marginBottom: '3rem' }}>
            <div className="ey">Docking Engines</div>
            <h2 className="sh">4 docking engines. <span className="a">From classical to deep learning.</span></h2>
            <p className="ss">Fast CPU screening or cutting-edge diffusion-based docking. Two engines are free. Two use on-demand GPU.</p>
          </div>
          <div className="eng-grid">
            {/* Vina */}
            <div className="eng">
              <div className="badge" style={{ background: 'var(--mint-d)', color: 'var(--mint)' }}>Fast</div>
              <div className="eng-icon"><svg viewBox="0 0 30 30" fill="none"><rect x="5" y="8" width="20" height="14" rx="2.5" stroke="var(--mint)" strokeWidth="1.2" /><line x1="10" y1="12" x2="20" y2="12" stroke="var(--mint)" strokeWidth=".8" opacity=".3" /><line x1="10" y1="15" x2="18" y2="15" stroke="var(--mint)" strokeWidth=".8" opacity=".3" /><line x1="10" y1="18" x2="16" y2="18" stroke="var(--mint)" strokeWidth=".8" opacity=".3" /></svg></div>
              <div className="eng-n">AutoDock Vina</div>
              <div className="eng-sub">CPU · Empirical</div>
              <div className="eng-d">Gold standard for HTVS. Reliable scoring, high throughput.</div>
              <div className="eng-spec">Best for: initial screening<br />Free — CPU only</div>
            </div>
            {/* GNINA CPU */}
            <div className="eng">
              <div className="badge" style={{ background: 'var(--blue-d)', color: 'var(--blue)' }}>CNN</div>
              <div className="eng-icon"><svg viewBox="0 0 30 30" fill="none"><circle cx="10" cy="10" r="2.5" stroke="var(--blue)" strokeWidth="1" /><circle cx="20" cy="10" r="2.5" stroke="var(--blue)" strokeWidth="1" /><circle cx="10" cy="20" r="2.5" stroke="var(--blue)" strokeWidth="1" /><circle cx="20" cy="20" r="2.5" stroke="var(--blue)" strokeWidth="1" /><line x1="12.5" y1="10" x2="17.5" y2="10" stroke="var(--blue)" strokeWidth=".8" /><line x1="10" y1="12.5" x2="10" y2="17.5" stroke="var(--blue)" strokeWidth=".8" /></svg></div>
              <div className="eng-n">GNINA CPU</div>
              <div className="eng-sub">CPU · CNN-Scored</div>
              <div className="eng-d">Vina sampling + CNN scoring. Better pose ranking.</div>
              <div className="eng-spec">Best for: accuracy on CPU<br />Free — CPU only</div>
            </div>
            {/* GNINA GPU */}
            <div className="eng">
              <div className="badge" style={{ background: 'var(--amber-d)', color: 'var(--amber)' }}>GPU</div>
              <div className="eng-icon"><svg viewBox="0 0 30 30" fill="none"><rect x="7" y="7" width="16" height="16" rx="2.5" stroke="var(--amber)" strokeWidth="1.2" /><rect x="10" y="10" width="4" height="4" rx=".8" fill="var(--amber)" opacity=".3" /><rect x="16" y="10" width="4" height="4" rx=".8" fill="var(--amber)" opacity=".3" /><rect x="10" y="16" width="4" height="4" rx=".8" fill="var(--amber)" opacity=".3" /><rect x="16" y="16" width="4" height="4" rx=".8" fill="var(--amber)" opacity=".3" /></svg></div>
              <div className="eng-n">GNINA GPU</div>
              <div className="eng-sub">RunPod · Serverless</div>
              <div className="eng-d">Full GPU acceleration. Maximum throughput for large libraries.</div>
              <div className="eng-spec">Best for: 10k+ compounds<br />€1 per 100 molecules</div>
            </div>
            {/* DiffDock */}
            <div className="eng">
              <div className="badge" style={{ background: 'var(--violet-d)', color: 'var(--violet)' }}>AI</div>
              <div className="eng-icon"><svg viewBox="0 0 30 30" fill="none"><circle cx="15" cy="15" r="8" stroke="var(--violet)" strokeWidth="1.2" strokeDasharray="2.5 1.5" /><circle cx="15" cy="15" r="3" fill="var(--violet)" opacity=".25" /><circle cx="15" cy="15" r="1.2" fill="var(--violet)" /></svg></div>
              <div className="eng-n">DiffDock</div>
              <div className="eng-sub">Diffusion · Deep Learning</div>
              <div className="eng-d">Generative diffusion model. State-of-the-art for flexible targets.</div>
              <div className="eng-spec">Best for: challenging sites<br />€1 per 100 molecules</div>
            </div>
          </div>
        </div>
      </section>

      {/* ═══ COMPOUND LIBRARIES ═══ */}
      <section className="sec" id="libraries">
        <div className="inner">
          <div className="fh-row">
            <div className="fh-text">
              <div className="ey">Compound Libraries</div>
              <h3>Connected to every <span className="a">public database.</span></h3>
              <p>Don't waste time downloading, converting, and filtering compound libraries manually. BindX is <strong>pre-connected to all major public compound databases</strong> — and gives you a smart interface to build the perfect screening set in minutes.</p>
              <p>Filter by molecular weight, logP, Lipinski compliance, scaffold diversity, or target class. Preview your library statistics before committing to a screen. Or upload your own SDF/SMILES files.</p>
              <div className="fh-tags">
                {['ZINC 22', 'Enamine REAL', 'ChEMBL 34', 'MolPort', 'PubChem', 'DrugBank', 'Custom SDF', 'SMILES upload'].map(t => (
                  <span className="fh-tag" key={t}>{t}</span>
                ))}
              </div>
            </div>
            <div className="fh-visual">
              <div className="fh-terminal">
                <div className="ft-line"><span className="ft-p">→</span><span className="ft-s">Source:</span> <span className="ft-v">Enamine REAL Diversity</span></div>
                <div className="ft-line"><span className="ft-p">→</span><span className="ft-s">Filters:</span> <span className="ft-d">MW 250–500 · logP ≤ 4 · Lipinski ✓</span></div>
                <div className="ft-line"><span className="ft-p">→</span><span className="ft-s">Diversity:</span> <span className="ft-d">Tanimoto cutoff 0.6 · max 50k</span></div>
                <div className="ft-line"><span className="ft-p">→</span><span className="ft-s">Preview:</span> <span className="ft-h">48,721 compounds</span> <span className="ft-d">· 312 scaffolds</span></div>
                <div className="ft-line"><span className="ft-p">→</span><span className="ft-s">Stats:</span> <span className="ft-d">MW avg 387 · logP avg 2.1 · QED avg 0.72</span></div>
                <div className="ft-line"><span className="ft-p" style={{ color: 'var(--mint)' }}>✓</span><span style={{ color: 'var(--mint)' }}>Library ready — import to Phase A</span></div>
              </div>
            </div>
          </div>
        </div>
      </section>

      {/* ═══ AI AGENTS ═══ */}
      <section className="sec sec-alt">
        <div className="inner">
          <div className="fh-row rev">
            <div className="fh-text">
              <div className="ey">AI-Assisted Pipeline</div>
              <h3>Guided by AI. <span className="a">Controlled by you.</span></h3>
              <p>Every step of the pipeline is assisted by <strong>specialized AI agents</strong> — target analysis, pocket selection, library curation, docking configuration, lead selection, generation strategy. They analyze your data in real-time and suggest optimal next steps.</p>
              <p>But you're always in control. Follow the AI recommendations end-to-end for a fast, automated campaign. Or override any decision — every parameter, every threshold, every selection is yours to change. <strong>100% autopilot or 100% manual. Or anything in between.</strong></p>
              <div className="fh-tags">
                {['Target agent', 'Library agent', 'Docking agent', 'Selection agent', 'Generation agent', 'Safety agent'].map(t => (
                  <span className="fh-tag" key={t}>{t}</span>
                ))}
              </div>
            </div>
            <div className="fh-visual">
              <div className="slider-visual">
                <div style={{ fontFamily: 'var(--mono)', fontSize: '.6rem', color: 'var(--dim)', textAlign: 'center', marginBottom: '1.2rem', letterSpacing: '.1em', textTransform: 'uppercase' }}>Autonomy Level</div>
                <div className="slider-bar">
                  <div className="slider-label" style={{ color: 'var(--violet)' }}>Full AI</div>
                  <div className="slider-track"><div className="slider-thumb" /></div>
                  <div className="slider-label" style={{ color: 'var(--mint)' }}>Full manual</div>
                </div>
                <div className="slider-row" style={{ marginTop: '1.5rem' }}>
                  <div className="slider-opt"><div className="v" style={{ color: 'var(--violet)' }}>Autopilot</div><div className="l">AI picks pocket, library,<br />engines, thresholds</div></div>
                  <div className="slider-opt"><div className="v" style={{ color: 'var(--blue)' }}>Guided</div><div className="l">AI recommends,<br />you approve each step</div></div>
                  <div className="slider-opt"><div className="v" style={{ color: 'var(--mint)' }}>Expert</div><div className="l">Full control on every<br />parameter and selection</div></div>
                </div>
                <div style={{ marginTop: '1.5rem', paddingTop: '1rem', borderTop: '1px solid var(--border)', fontSize: '.72rem', color: 'var(--dim)', lineHeight: 1.5, textAlign: 'center' }}>
                  "Dock the top pocket with GNINA GPU, exhaustiveness 12, then run ADMET + scoring on the top 500 by CNNaffinity."<br />
                  <span style={{ color: 'var(--mint)', fontWeight: 600 }}>— AI suggestion you can accept, edit, or ignore</span>
                </div>
              </div>
            </div>
          </div>
        </div>
      </section>

      {/* ═══ EXPERIMENTAL RESULTS ═══ */}
      <section className="sec">
        <div className="inner">
          <div className="fh-row">
            <div className="fh-text">
              <div className="ey">Experimental Results</div>
              <h3>Close the loop between <span className="a">prediction and reality.</span></h3>
              <p>When your candidates reach the lab, the story doesn't end. BindX includes a dedicated space to <strong>upload your experimental data</strong> — IC₅₀, Ki, selectivity ratios, ADMET assay results — and link them back to the molecules in your dashboard.</p>
              <p>This does two things: it lets you <strong>compare in silico predictions vs. real-world outcomes</strong> directly in the interface. And over time, this data is used to <strong>recalibrate scoring models</strong> for your target — making each new round of optimization more accurate than the last.</p>
              <div className="fh-tags">
                {['IC₅₀ / Ki upload', 'Selectivity data', 'ADMET assays', 'Prediction vs. reality', 'Model recalibration'].map(t => (
                  <span className="fh-tag" key={t}>{t}</span>
                ))}
              </div>
            </div>
            <div className="fh-visual">
              <div className="lab-flow">
                <div className="lab-step"><span className="num">01</span><span className="txt"><strong>Discover & optimize</strong> — run the full pipeline in BindX</span></div>
                <div className="lab-arrow">↓</div>
                <div className="lab-step"><span className="num">02</span><span className="txt"><strong>Synthesize & assay</strong> — test your top candidates at the bench</span></div>
                <div className="lab-arrow">↓</div>
                <div className="lab-step"><span className="num">03</span><span className="txt"><strong>Upload results</strong> — paste IC₅₀, Ki, selectivity into the Experimental Results panel</span></div>
                <div className="lab-arrow">↓</div>
                <div className="lab-step"><span className="num">04</span><span className="txt"><strong>Compare</strong> — see where predictions matched or diverged, per molecule</span></div>
                <div className="lab-arrow">↓</div>
                <div className="lab-step" style={{ borderColor: 'rgba(0,230,160,.15)', background: 'var(--mint-dd)' }}><span className="num" style={{ color: 'var(--mint)' }}>05</span><span className="txt"><strong style={{ color: 'var(--mint)' }}>Recalibrate</strong> — BindX adjusts its scoring models using your real data for the next round</span></div>
              </div>
            </div>
          </div>
        </div>
      </section>

      {/* ═══ PRICING ═══ */}
      <section className="sec" id="pricing">
        <div className="inner">
          <div className="center" style={{ marginBottom: '3rem' }}>
            <div className="ey">Pricing</div>
            <h2 className="sh">Free to use. <span className="a">Seriously.</span></h2>
            <p className="ss">Other platforms charge $7,500–200,000/year. BindX is free — you only pay for GPU when you need it.</p>
          </div>
          <div className="pricing-row">
            {/* Free tier */}
            <div className="price-card featured">
              <div className="price-inner">
                <div className="price-label">Platform Access</div>
                <div className="price-amount"><span className="a">Free</span></div>
                <div className="price-note">Not a trial. Not a demo. The full platform.</div>
                <ul className="price-list">
                  <li><CheckIcon /> Unlimited projects & campaigns</li>
                  <li><CheckIcon /> All 9 computation types — 37 properties</li>
                  <li><CheckIcon /> Vina + GNINA CPU docking — free</li>
                  <li><CheckIcon /> ADMET, safety, retrosynthesis, confidence</li>
                  <li><CheckIcon /> AI-powered molecular generation</li>
                  <li><CheckIcon /> Pareto optimization & 3D viewer</li>
                  <li><CheckIcon /> CSV / SDF / PDF export</li>
                </ul>
                <Link to="/register" className="btn btn-m" style={{ width: '100%', justifyContent: 'center', padding: '.6rem' }}>Start free →</Link>
              </div>
            </div>
            {/* GPU tier */}
            <div className="price-card">
              <div className="price-inner">
                <div className="price-label">GPU Compute</div>
                <div className="price-amount">€1 <span style={{ fontSize: '1rem', color: 'var(--sub)', fontWeight: 400, fontFamily: 'var(--body)' }}>/ 100 molecules</span></div>
                <div className="price-note">Pay only when you use GPU-accelerated docking.</div>
                <div className="gpu-note">
                  <p><strong>When do you need GPU?</strong> Only for GNINA GPU and DiffDock engines. CPU-based docking (Vina, GNINA CPU) is free and works for most screening campaigns. GPU is recommended for libraries above 10,000 compounds or when using deep learning docking.</p>
                </div>
                <ul className="price-list">
                  <li><CheckIcon /> GNINA GPU — serverless via RunPod</li>
                  <li><CheckIcon /> DiffDock — diffusion-based docking</li>
                  <li><CheckIcon /> No subscription — pay per use</li>
                  <li><CheckIcon /> Example: 20k compounds = €200</li>
                </ul>
                <a href="#compute" className="btn btn-o" style={{ width: '100%', justifyContent: 'center', padding: '.6rem' }}>Learn more</a>
              </div>
            </div>
          </div>
        </div>
      </section>

      {/* ═══ EARLY ACCESS ═══ */}
      <section className="sec sec-alt" id="early-access">
        <div className="inner">
          <div className="ea-card">
            <div className="ea-inner">
              <div className="ey">Scientific Early Access Program</div>
              <h2 className="sh">Validate with us.<br />Get <span className="a">full access — free GPU included.</span></h2>
              <p className="ss" style={{ maxWidth: 600 }}>BindX predictions need real-world validation to improve. We're inviting computational chemists with active projects to use the full platform — including free GPU compute — in exchange for sharing anonymized experimental outcomes.</p>
              <div className="ea-grid">
                <div className="ea-col">
                  <h4><svg viewBox="0 0 16 16" fill="none" width="16"><path d="M3 8l4 4 6-6" stroke="var(--mint)" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round" /></svg> What you get</h4>
                  <ul>
                    <li><svg viewBox="0 0 14 14" fill="none"><path d="M3 7l3 3 5-5" stroke="var(--mint)" strokeWidth="1.3" strokeLinecap="round" strokeLinejoin="round" /></svg> Full platform access — no restrictions</li>
                    <li><svg viewBox="0 0 14 14" fill="none"><path d="M3 7l3 3 5-5" stroke="var(--mint)" strokeWidth="1.3" strokeLinecap="round" strokeLinejoin="round" /></svg> Free GPU compute (GNINA GPU + DiffDock)</li>
                    <li><svg viewBox="0 0 14 14" fill="none"><path d="M3 7l3 3 5-5" stroke="var(--mint)" strokeWidth="1.3" strokeLinecap="round" strokeLinejoin="round" /></svg> Priority support & feature requests</li>
                    <li><svg viewBox="0 0 14 14" fill="none"><path d="M3 7l3 3 5-5" stroke="var(--mint)" strokeWidth="1.3" strokeLinecap="round" strokeLinejoin="round" /></svg> Early access to future modules</li>
                  </ul>
                </div>
                <div className="ea-col">
                  <h4><svg viewBox="0 0 16 16" fill="none" width="16"><circle cx="8" cy="8" r="6" stroke="var(--blue)" strokeWidth="1.3" /><path d="M8 5v3l2 2" stroke="var(--blue)" strokeWidth="1.3" strokeLinecap="round" /></svg> What we ask</h4>
                  <ul>
                    <li><svg viewBox="0 0 14 14" fill="none"><circle cx="7" cy="7" r="2" fill="var(--blue)" /></svg> Share anonymized experimental results (IC₅₀, Ki, selectivity) when compounds are tested</li>
                    <li><svg viewBox="0 0 14 14" fill="none"><circle cx="7" cy="7" r="2" fill="var(--blue)" /></svg> Flag where BindX predictions matched or diverged from reality</li>
                    <li><svg viewBox="0 0 14 14" fill="none"><circle cx="7" cy="7" r="2" fill="var(--blue)" /></svg> Your data stays yours — we only publish aggregate prediction benchmarks</li>
                    <li><svg viewBox="0 0 14 14" fill="none"><circle cx="7" cy="7" r="2" fill="var(--blue)" /></svg> Optional: workflow feedback and feature requests</li>
                  </ul>
                </div>
              </div>
              <div style={{ marginTop: '2rem', textAlign: 'center' }}>
                <Link to="/register" className="btn btn-m" style={{ padding: '.65rem 2rem', fontSize: '.85rem' }}>Apply for Early Access</Link>
              </div>
            </div>
          </div>
        </div>
      </section>

      {/* ═══ DIFFERENTIATION TABLE ═══ */}
      <section className="sec">
        <div className="inner">
          <div className="center" style={{ marginBottom: '2rem' }}>
            <div className="ey">Honest Comparison</div>
            <h2 className="sh">Where BindX <span className="a">fits</span></h2>
          </div>
          <table className="diff-table">
            <thead>
              <tr>
                <th>Capability</th>
                <th className="bindx">BindX</th>
                <th>Enterprise platforms</th>
                <th>Open-source tools</th>
              </tr>
            </thead>
            <tbody>
              <tr><td className="feat">Molecular docking</td><td className="yes">✓ 4 engines incl. DiffDock</td><td className="yes">✓ Glide, GOLD, etc.</td><td className="mid">∼ Vina only, manual setup</td></tr>
              <tr><td className="feat">ADMET profiling (8 props)</td><td className="yes">✓ Built-in</td><td className="yes">✓ Separate module</td><td className="mid">∼ ADMETlab, separate tool</td></tr>
              <tr><td className="feat">Safety profiling (6 props)</td><td className="yes">✓ Built-in</td><td className="mid">∼ Add-on</td><td className="no">✗ Not available</td></tr>
              <tr><td className="feat">Retrosynthesis</td><td className="yes">✓ Built-in</td><td className="mid">∼ Separate product</td><td className="no">✗ Not available</td></tr>
              <tr><td className="feat">AI-powered generation</td><td className="yes">✓ Batch + molecule mode</td><td className="mid">∼ Separate product</td><td className="mid">∼ REINVENT, manual</td></tr>
              <tr><td className="feat">Pareto optimization</td><td className="yes">✓ Built-in</td><td className="mid">∼ Requires scripting</td><td className="no">✗ Not available</td></tr>
              <tr><td className="feat">Unified dashboard</td><td className="yes">✓ Cumulative, progressive</td><td className="mid">∼ LiveDesign-like</td><td className="no">✗ Separate files</td></tr>
              <tr><td className="feat">No scripting required</td><td className="yes">✓ Click-based workflow</td><td className="yes">✓ GUI-based</td><td className="no">✗ Command-line</td></tr>
              <tr><td className="feat">AI agents assist every step</td><td className="yes">✓ 6 specialized agents</td><td className="mid">∼ Basic wizards</td><td className="no">✗ Manual only</td></tr>
              <tr><td className="feat">Public library connectors</td><td className="yes">✓ ZINC, Enamine, ChEMBL+</td><td className="yes">✓ Vendor catalogs</td><td className="mid">∼ Manual download</td></tr>
              <tr><td className="feat">Price</td><td className="yes">Free (GPU: €1/100 mols)</td><td className="price">$7,500–200,000+/yr</td><td className="yes">Free (but fragmented)</td></tr>
            </tbody>
          </table>
        </div>
      </section>

      {/* ═══ CTA ═══ */}
      <section className="cta-sec">
        <div className="cta-i cta">
          <h2>Target → Hits → Leads → <span className="a">Candidates.</span></h2>
          <p>The full in silico pipeline in one place. Free to start. Built for real science.</p>
          <div className="cta-btns">
            <Link to="/register" className="btn btn-m" style={{ padding: '.7rem 2.2rem', fontSize: '.9rem' }}>Create your first project</Link>
            <a href="#early-access" className="btn btn-o" style={{ padding: '.7rem 2.2rem', fontSize: '.9rem' }}>Apply for Early Access</a>
          </div>
        </div>
      </section>

      {/* ═══ FOOTER ═══ */}
      <footer className="lfoot" ref={footerRef}>
        <canvas ref={footerCanvasRef} style={{ position: 'absolute', inset: 0, zIndex: 0, pointerEvents: 'none' }} />
        <div className="fl">
          <FooterLogoSvg />
          <span>BindX — In Silico Drug Discovery Platform</span>
        </div>
        <div className="fr">
          <a href="#">Documentation</a>
          <a href="#">API</a>
          <a href="#">Contact</a>
          <a href="#">Privacy</a>
        </div>
      </footer>
    </div>
  )
}
