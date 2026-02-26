import React from 'react'
import { Routes, Route, Navigate } from 'react-router-dom'
import { WorkspaceProvider } from './contexts/WorkspaceContext.jsx'
import { ToastProvider } from './contexts/ToastContext.jsx'

// Layout
import SidebarLayout from './components/SidebarLayout.jsx'

// Main pages
import ProjectListPage from './pages/ProjectListPage.jsx'
import ProjectHome from './pages/ProjectHome.jsx'

// Phase dashboard
import PhaseDashboard from './pages/PhaseDashboard.jsx'

// Static pages
import MethodologyPage from './components/MethodologyPage.jsx'

// PharmacoDB
import PharmacoDBDashboard from './pages/PharmacoDBDashboard.jsx'

// --------------------------------------------------
// App — V9 routing
// --------------------------------------------------
export default function App() {
  return (
    <WorkspaceProvider>
    <ToastProvider>
      <Routes>
        {/* TODO: Auth pages (Supabase Auth V9) */}
        <Route path="/login" element={<Navigate to="/" replace />} />
        <Route path="/register" element={<Navigate to="/" replace />} />

        {/* Main app — with sidebar */}
        <Route element={<SidebarLayout />}>
          {/* Home — Project list */}
          <Route path="/" element={<ProjectListPage />} />

          {/* Project overview */}
          <Route path="/project/:projectId" element={<ProjectHome />} />

          {/* Phase dashboard */}
          <Route
            path="/project/:projectId/phase/:phaseId"
            element={<PhaseDashboard />}
          />

          {/* Static reference pages */}
          <Route path="/methodology" element={<MethodologyPage />} />

          {/* PharmacoDB */}
          <Route path="/pharmacodb" element={<PharmacoDBDashboard />} />

          {/* Catch-all → home */}
          <Route path="*" element={<Navigate to="/" replace />} />
        </Route>
      </Routes>
    </ToastProvider>
    </WorkspaceProvider>
  )
}
