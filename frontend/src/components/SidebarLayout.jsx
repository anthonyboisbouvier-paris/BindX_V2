import React, { useEffect } from 'react'
import { Outlet, NavLink, Link, useParams, useNavigate } from 'react-router-dom'
import { useWorkspace } from '../contexts/WorkspaceContext.jsx'

// ---------------------------------------------------------------------------
// Brand logo SVG
// ---------------------------------------------------------------------------

function BindXLogo({ size = 32 }) {
  return (
    <svg width={size} height={size} viewBox="0 0 40 40" fill="none" xmlns="http://www.w3.org/2000/svg">
      <rect width="40" height="40" rx="8" fill="#1e3a5f" />
      <circle cx="20" cy="15" r="6" fill="none" stroke="#22c55e" strokeWidth="2.5" />
      <circle cx="10" cy="29" r="3.5" fill="none" stroke="#22c55e" strokeWidth="2" />
      <circle cx="30" cy="29" r="3.5" fill="none" stroke="#22c55e" strokeWidth="2" />
      <line x1="20" y1="21" x2="10" y2="29" stroke="#22c55e" strokeWidth="1.8" />
      <line x1="20" y1="21" x2="30" y2="29" stroke="#22c55e" strokeWidth="1.8" />
      <circle cx="20" cy="15" r="2.5" fill="#22c55e" />
      <circle cx="10" cy="29" r="2" fill="#22c55e" />
      <circle cx="30" cy="29" r="2" fill="#22c55e" />
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

function IconDatabase() {
  return (
    <svg className="w-4 h-4 shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.8}
        d="M4 7c0-1.657 3.582-3 8-3s8 1.343 8 3M4 7v5c0 1.657 3.582 3 8 3s8-1.343 8-3V7M4 7v5m16-5v5M4 12v5c0 1.657 3.582 3 8 3s8-1.343 8-3v-5" />
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
    return <span className="w-2 h-2 rounded-full bg-dockit-green shrink-0 shadow-sm" style={{ boxShadow: '0 0 0 2px rgba(34,197,94,0.3)' }} />
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

  const baseClass = 'flex items-center gap-2 px-2.5 py-1.5 rounded-md text-xs transition-colors duration-150 w-full text-left'

  if (notCreated) {
    return (
      <div className={`${baseClass} text-white/25 cursor-default`}>
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
            ? 'bg-white/10 text-white font-medium'
            : 'text-white/60 hover:text-white hover:bg-white/5',
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
            'flex items-center gap-2 px-3 py-2 rounded-lg text-sm transition-colors duration-150',
            isActive
              ? 'bg-white/10 text-white font-semibold'
              : 'text-white/80 hover:text-white hover:bg-white/5',
          ].join(' ')
        }
      >
        <IconFolder />
        <span className="flex-1 truncate font-medium">{project.name}</span>
      </NavLink>

      {campaign ? (
        <div className="ml-3 border-l border-white/10 pl-2 mt-1">
          {/* Campaign header — not a link, just a label */}
          <div className="px-2 py-1.5 mb-0.5">
            <p className="text-white/40 text-[10px] uppercase tracking-wider font-semibold truncate">
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
        <div className="ml-3 border-l border-white/10 pl-2 mt-1">
          <p className="px-2 py-1.5 text-white/25 text-xs italic">No campaigns yet</p>
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
          'flex items-center gap-2 px-3 py-1.5 text-xs rounded-md transition-colors duration-150',
          isActive
            ? 'bg-white/10 text-white font-medium'
            : 'text-white/40 hover:text-white hover:bg-white/5',
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

export default function SidebarLayout() {
  const { projectId, phaseId } = useParams()
  const { projects, currentProject, selectProject } = useWorkspace()

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
      <aside className="hidden md:flex md:flex-col w-60 bg-[#1e3a5f] shrink-0 min-h-screen">

        {/* Brand header */}
        <div className="px-4 py-4 border-b border-white/10">
          <Link to="/" className="flex items-center gap-2.5 group">
            <BindXLogo size={32} />
            <div>
              <p className="text-white font-bold text-base leading-none tracking-tight group-hover:text-dockit-green transition-colors">
                BindX
              </p>
              <p className="text-white/30 text-[10px] leading-none mt-0.5 font-mono">v9.0</p>
            </div>
          </Link>
        </div>

        {/* All Projects link */}
        <div className="px-3 pt-3 pb-1">
          <NavLink
            to="/"
            end
            className={({ isActive }) =>
              [
                'flex items-center gap-2 px-3 py-1.5 rounded-md text-sm transition-colors duration-150',
                isActive && !projectId
                  ? 'bg-white/10 text-white font-medium'
                  : 'text-white/55 hover:text-white hover:bg-white/5',
              ].join(' ')
            }
          >
            <IconArrowLeft />
            <span>All Projects</span>
          </NavLink>
        </div>

        {/* Divider + project tree (when a project is selected) */}
        {projectId && currentProject && (
          <div className="px-3 pb-2">
            <div className="h-px bg-white/10 my-2" />
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
            <p className="px-3 text-white/25 text-xs italic">Select a project to begin</p>
          </nav>
        )}

        {/* Spacer */}
        <div className="flex-1" />

        {/* Bottom links */}
        <div className="px-3 py-3 border-t border-white/10 space-y-0.5">
          <BottomNavLink to="/pharmacodb" icon={<IconDatabase />} label="PharmacoDB" />
          <BottomNavLink to="/references" icon={<IconBook />} label="References" />
          <BottomNavLink to="/methodology" icon={<IconFlask />} label="Methodology" />
          <p className="text-white/20 text-[10px] px-3 pt-1.5 font-mono">BindX v9.0</p>
        </div>
      </aside>

      {/* Mobile horizontal nav */}
      <div className="md:hidden fixed top-0 inset-x-0 z-50 bg-[#1e3a5f] border-b border-white/10">
        <div className="flex items-center justify-between px-4 h-12">
          <Link to="/" className="flex items-center gap-2">
            <BindXLogo size={24} />
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

          {/* Back to all projects link on mobile */}
          {!projectId && (
            <NavLink to="/" className="text-white/60 hover:text-white text-xs">
              All Projects
            </NavLink>
          )}
        </div>
      </div>

      {/* Content area */}
      <main className="flex-1 min-w-0 bg-dockit-gray">
        {/* Mobile spacer */}
        <div className="md:hidden h-12" aria-hidden="true" />
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
          <Outlet />
        </div>
      </main>
    </div>
  )
}
