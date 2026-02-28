import React, { useEffect } from 'react'
import { Outlet, NavLink, Link, useParams, useNavigate } from 'react-router-dom'
import { useWorkspace } from '../contexts/WorkspaceContext.jsx'
import { useAuth } from '../contexts/AuthContext'

// ---------------------------------------------------------------------------
// Brand logo SVG
// ---------------------------------------------------------------------------

function BindXNavLogo({ size = 32 }) {
  // Static 3-circles logo — stroke/dot sizes adapt for small sizes (per design system)
  const small = size <= 26
  const sw1 = small ? 2 : 1.5
  const sw2 = small ? 2.5 : 2
  const sw3 = small ? 3 : 2.5
  const dotR = small ? 4.5 : 4
  return (
    <svg width={size} height={size} viewBox="0 0 56 56" fill="none" style={{ overflow: 'visible' }}>
      <circle cx="18" cy="28" r="15" stroke="#00e6a0" strokeWidth={sw1} fill="none" opacity=".5" />
      <circle cx="28" cy="28" r="12" stroke="#06b6d4" strokeWidth={sw2} fill="none" opacity=".7" />
      <circle cx="38" cy="28" r="9" stroke="#3b82f6" strokeWidth={sw3} fill="none" />
      <circle cx="38" cy="28" r={dotR} fill="#3b82f6" opacity=".85" />
    </svg>
  )
}

// ---------------------------------------------------------------------------
// Icons
// ---------------------------------------------------------------------------

function IconFolder() {
  return (
    <svg className="w-4 h-4 shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
        d="M3 7v10a2 2 0 002 2h14a2 2 0 002-2V9a2 2 0 00-2-2h-6l-2-2H5a2 2 0 00-2 2z" />
    </svg>
  )
}

function IconArrowLeft() {
  return (
    <svg className="w-4 h-4 shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10 19l-7-7m0 0l7-7m-7 7h18" />
    </svg>
  )
}

function IconLock() {
  return (
    <svg className="w-3.5 h-3.5 shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
        d="M12 15v2m-6 4h12a2 2 0 002-2v-6a2 2 0 00-2-2H6a2 2 0 00-2 2v6a2 2 0 002 2zm10-10V7a4 4 0 00-8 0v4h8z" />
    </svg>
  )
}

function IconBook() {
  return (
    <svg className="w-4 h-4 shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
        d="M12 6.253v13m0-13C10.832 5.477 9.246 5 7.5 5S4.168 5.477 3 6.253v13C4.168 18.477 5.754 18 7.5 18s3.332.477 4.5 1.253m0-13C13.168 5.477 14.754 5 16.5 5c1.746 0 3.332.477 4.5 1.253v13C19.832 18.477 18.246 18 16.5 18c-1.746 0-3.332.477-4.5 1.253" />
    </svg>
  )
}

function IconFlask() {
  return (
    <svg className="w-4 h-4 shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
        d="M9.75 3.104v5.714a2.25 2.25 0 01-.659 1.591L5 14.5M9.75 3.104c-.251.023-.501.05-.75.082m.75-.082a24.301 24.301 0 014.5 0m0 0v5.714c0 .597.237 1.17.659 1.591L19 14.5M14.25 3.104c.251.023.501.05.75.082M19 14.5l-1.47 4.9a2.25 2.25 0 01-2.156 1.6H8.626a2.25 2.25 0 01-2.156-1.6L5 14.5m14 0H5" />
    </svg>
  )
}


// ---------------------------------------------------------------------------
// Phase status indicator
// ---------------------------------------------------------------------------

function PhaseStatusDot({ status, notCreated }) {
  if (notCreated) {
    return <span className="w-2 h-2 rounded-full bg-white/20 shrink-0" />
  }
  if (status === 'frozen') {
    return <IconLock />
  }
  if (status === 'active') {
    return <span className="w-2 h-2 rounded-full bg-bx-mint shrink-0 shadow-sm" style={{ boxShadow: '0 0 0 2px rgba(0,230,160,0.3)' }} />
  }
  // fallback
  return <span className="w-2 h-2 rounded-full bg-white/30 shrink-0" />
}

// ---------------------------------------------------------------------------
// Phase nav item
// ---------------------------------------------------------------------------

function PhaseNavItem({ phase, projectId, currentPhaseId }) {
  const notCreated = !phase.created_at
  const isFrozen = phase.status === 'frozen'
  const isCurrentPhase = phase.id === currentPhaseId

  const baseClass = 'flex items-center gap-2 px-2.5 py-1.5 rounded-[7px] text-[.68rem] transition-colors duration-150 w-full text-left'

  if (notCreated) {
    return (
      <div className={`${baseClass} text-bx-dim cursor-default`}>
        <PhaseStatusDot status={phase.status} notCreated={true} />
        <span className="flex-1 truncate">{phase.label}</span>
        <span className="text-white/20 text-[10px]">Not started</span>
      </div>
    )
  }

  return (
    <NavLink
      to={`/project/${projectId}/phase/${phase.id}`}
      className={({ isActive }) =>
        [
          baseClass,
          isFrozen ? 'opacity-70' : '',
          isActive
            ? 'bg-bx-mint/[.08] text-bx-mint font-semibold border border-bx-mint/10'
            : 'text-bx-sub hover:text-white hover:bg-white/5',
        ].join(' ')
      }
    >
      <PhaseStatusDot status={phase.status} notCreated={false} />
      <span className="flex-1 truncate">{phase.label}</span>
      {isFrozen && (
        <span className="text-white/30 text-[10px] font-mono">Frozen</span>
      )}
    </NavLink>
  )
}

// ---------------------------------------------------------------------------
// Project tree — campaign + phases
// ---------------------------------------------------------------------------

function ProjectTree({ project, projectId, currentPhaseId }) {
  const campaign = project.campaigns?.[0] || null

  return (
    <div className="mt-2">
      {/* Project name as nav link */}
      <NavLink
        to={`/project/${projectId}`}
        end
        className={({ isActive }) =>
          [
            'flex items-center gap-2 px-2.5 py-1.5 rounded-[7px] text-[.68rem] font-medium transition-colors duration-150',
            isActive
              ? 'bg-bx-mint/[.08] text-bx-mint font-semibold border border-bx-mint/10'
              : 'text-white/80 hover:text-white hover:bg-white/5',
          ].join(' ')
        }
      >
        <IconFolder />
        <span className="flex-1 truncate font-medium">{project.name}</span>
      </NavLink>

      {campaign ? (
        <div className="ml-3 border-l border-white/[.07] pl-2 mt-1">
          {/* Campaign header — not a link, just a label */}
          <div className="px-2 py-1.5 mb-0.5">
            <p className="text-bx-dim text-[10px] uppercase tracking-[.1em] font-bold truncate">
              {campaign.name}
            </p>
          </div>

          {/* Phase items */}
          <div className="space-y-0.5">
            {campaign.phases.map((phase) => (
              <PhaseNavItem
                key={phase.id}
                phase={phase}
                projectId={projectId}
                currentPhaseId={currentPhaseId}
              />
            ))}
          </div>
        </div>
      ) : (
        <div className="ml-3 border-l border-white/[.07] pl-2 mt-1">
          <p className="px-2 py-1.5 text-bx-dim text-xs italic">No campaigns yet</p>
        </div>
      )}
    </div>
  )
}

// ---------------------------------------------------------------------------
// Bottom nav links
// ---------------------------------------------------------------------------

function BottomNavLink({ to, icon, label }) {
  return (
    <NavLink
      to={to}
      className={({ isActive }) =>
        [
          'flex items-center gap-2 px-2.5 py-1.5 text-[.68rem] rounded-[7px] font-medium transition-colors duration-150',
          isActive
            ? 'bg-bx-mint/[.08] text-bx-mint font-semibold border border-bx-mint/10'
            : 'text-bx-sub hover:text-white hover:bg-white/5',
        ].join(' ')
      }
    >
      {icon}
      <span>{label}</span>
    </NavLink>
  )
}

// ---------------------------------------------------------------------------
// SidebarLayout — main export
// ---------------------------------------------------------------------------

function IconLogout() {
  return (
    <svg className="w-4 h-4 shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
        d="M15.75 9V5.25A2.25 2.25 0 0013.5 3h-6a2.25 2.25 0 00-2.25 2.25v13.5A2.25 2.25 0 007.5 21h6a2.25 2.25 0 002.25-2.25V15m3-3l3-3m0 0l-3-3m3 3H9" />
    </svg>
  )
}

export default function SidebarLayout() {
  const { projectId, phaseId } = useParams()
  const { projects, currentProject, selectProject } = useWorkspace()
  const { user, logout } = useAuth()
  const navigate = useNavigate()

  const handleLogout = async () => {
    await logout()
    navigate('/welcome')
  }

  // Sync workspace state with URL param
  useEffect(() => {
    if (projectId) {
      selectProject(projectId)
    }
  }, [projectId, selectProject])

  // Mobile phase tabs
  const mobilePhaseTabs = currentProject?.campaigns?.[0]?.phases ?? []

  return (
    <div className="flex min-h-screen">
      {/* Desktop sidebar */}
      <aside className="hidden md:flex md:flex-col w-[230px] bg-bx-bg shrink-0 min-h-screen">

        {/* Brand header */}
        <div className="px-4 py-4 border-b border-white/[.07]">
          <Link to="/" className="flex items-center gap-2.5 group">
            <BindXNavLogo size={26} />
            <div>
              <p className="text-white font-extrabold text-[1.05rem] leading-none tracking-tight group-hover:text-bx-mint transition-colors" style={{ letterSpacing: '-.02em' }}>
                Bind<span className="text-bx-mint">X</span>
              </p>
              <p className="text-bx-dim text-[10px] leading-none mt-0.5 font-mono">v9.0</p>
            </div>
          </Link>
        </div>

        {/* Navigation link — context-dependent */}
        <div className="px-3 pt-3 pb-1">
          <p className="px-2 mb-1 text-[10px] font-bold uppercase tracking-[.1em] text-bx-dim">Main</p>
          {projectId ? (
            <NavLink
              to="/"
              className="flex items-center gap-2 px-2.5 py-1.5 rounded-[7px] text-[.68rem] font-medium transition-colors duration-150 text-bx-sub hover:text-white hover:bg-white/5"
            >
              <IconArrowLeft />
              <span>All Projects</span>
            </NavLink>
          ) : (
            <NavLink
              to="/"
              className={({ isActive }) =>
                [
                  'flex items-center gap-2 px-2.5 py-1.5 rounded-[7px] text-[.68rem] font-medium transition-colors duration-150',
                  isActive
                    ? 'bg-bx-mint/[.08] text-bx-mint font-semibold border border-bx-mint/10'
                    : 'text-bx-sub hover:text-white hover:bg-white/5',
                ].join(' ')
              }
            >
              <IconFolder />
              <span>Projects</span>
            </NavLink>
          )}
        </div>

        {/* Divider + project tree (when a project is selected) */}
        {projectId && currentProject && (
          <div className="px-3 pb-2">
            <div className="h-px bg-white/[.07] my-2" />
            <ProjectTree
              project={currentProject}
              projectId={projectId}
              currentPhaseId={phaseId}
            />
          </div>
        )}

        {/* Placeholder when no project is selected */}
        {!projectId && (
          <nav className="flex-1 px-3 py-3" aria-label="Sidebar placeholder">
            <p className="px-3 text-bx-dim text-xs italic">Select a project to begin</p>
          </nav>
        )}

        {/* Spacer */}
        <div className="flex-1" />

        {/* Bottom links */}
        <div className="px-3 py-3 border-t border-white/[.07] space-y-0.5">
          <p className="px-2 mb-1 text-[10px] font-bold uppercase tracking-[.1em] text-bx-dim">Resources</p>
          <BottomNavLink to="/references" icon={<IconBook />} label="References" />
          <BottomNavLink to="/methodology" icon={<IconFlask />} label="Methodology" />
          <button
            onClick={handleLogout}
            className="flex items-center gap-2 px-2.5 py-1.5 text-[.68rem] rounded-[7px] font-medium transition-colors duration-150 text-bx-sub hover:text-white hover:bg-white/5 w-full"
          >
            <IconLogout />
            <span>Sign out</span>
          </button>
          {user?.email && (
            <p className="text-bx-dim text-[10px] px-3 pt-0.5 truncate font-mono">{user.email}</p>
          )}
          <p className="text-bx-dim text-[10px] px-3 pt-0.5 font-mono">BindX v9.0</p>
        </div>
      </aside>

      {/* Mobile horizontal nav */}
      <div className="md:hidden fixed top-0 inset-x-0 z-50 bg-bx-bg border-b border-white/[.07]">
        <div className="flex items-center justify-between px-4 h-12">
          <Link to="/" className="flex items-center gap-2">
            <BindXNavLogo size={24} />
            <span className="text-white font-bold text-sm">BindX</span>
          </Link>

          {/* Phase tabs if inside a project */}
          {projectId && mobilePhaseTabs.length > 0 && (
            <nav className="flex items-center gap-1">
              {mobilePhaseTabs
                .filter(p => p.created_at)
                .map(phase => (
                  <NavLink
                    key={phase.id}
                    to={`/project/${projectId}/phase/${phase.id}`}
                    className={({ isActive }) =>
                      [
                        'px-2.5 py-1 rounded text-xs font-medium transition-colors',
                        isActive
                          ? 'bg-white/15 text-white'
                          : 'text-white/50 hover:text-white',
                      ].join(' ')
                    }
                  >
                    {phase.label}
                  </NavLink>
                ))}
            </nav>
          )}

          {/* Context-dependent back link on mobile */}
          {projectId ? (
            <NavLink to="/" className="text-white/60 hover:text-white text-xs">
              All Projects
            </NavLink>
          ) : (
            <NavLink to="/welcome" className="text-white/60 hover:text-white text-xs">
              Home
            </NavLink>
          )}

          <button
            onClick={handleLogout}
            className="text-white/40 hover:text-white text-xs ml-2"
          >
            Sign out
          </button>
        </div>
      </div>

      {/* Content area */}
      <main className="flex-1 min-w-0 bg-bx-light-bg">
        {/* Mobile spacer */}
        <div className="md:hidden h-12" aria-hidden="true" />
        <div className="max-w-7xl mx-auto" style={{ padding: '1.8rem 2.2rem' }}>
          <Outlet />
        </div>
      </main>
    </div>
  )
}
