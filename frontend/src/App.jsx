import React from 'react'
import { Routes, Route, Navigate } from 'react-router-dom'
import { WorkspaceProvider } from './contexts/WorkspaceContext.jsx'
import { ToastProvider } from './contexts/ToastContext.jsx'
import { useAuth } from './contexts/AuthContext'

// Layout
import SidebarLayout from './components/SidebarLayout.jsx'
import ErrorBoundary from './components/ErrorBoundary.jsx'

// Auth pages
import LoginPage from './pages/LoginPage.jsx'
import RegisterPage from './pages/RegisterPage.jsx'
import LandingPage from './pages/LandingPage.jsx'

// Main pages
import ProjectListPage from './pages/ProjectListPage.jsx'
import ProjectHome from './pages/ProjectHome.jsx'

// Phase dashboard
import PhaseDashboard from './pages/PhaseDashboard.jsx'

// Target setup
import TargetSetup from './pages/TargetSetup.jsx'

// Static pages
import MethodologyPage from './components/MethodologyPage.jsx'

// --------------------------------------------------
// ProtectedRoute — redirects to /login if not authenticated
// --------------------------------------------------
function ProtectedRoute({ children }) {
  const { isAuthenticated, loading } = useAuth()

  if (loading) {
    return (
      <div className="min-h-screen bg-dockit-gray flex items-center justify-center">
        <div className="w-8 h-8 border-2 border-bx-mint border-t-transparent rounded-full animate-spin" />
      </div>
    )
  }

  if (!isAuthenticated) {
    return <Navigate to="/welcome" replace />
  }

  return children
}

// --------------------------------------------------
// App — V9 routing
// --------------------------------------------------
export default function App() {
  return (
    <ErrorBoundary>
    <ToastProvider>
    <WorkspaceProvider>
      <Routes>
        {/* Public pages (no sidebar, no auth required) */}
        <Route path="/welcome" element={<LandingPage />} />
        <Route path="/login" element={<LoginPage />} />
        <Route path="/register" element={<RegisterPage />} />

        {/* Main app — protected, with sidebar */}
        <Route element={
          <ProtectedRoute>
            <SidebarLayout />
          </ProtectedRoute>
        }>
          {/* Home — Project list */}
          <Route path="/" element={<ProjectListPage />} />

          {/* Project overview */}
          <Route path="/project/:projectId" element={<ProjectHome />} />

          {/* Target setup */}
          <Route path="/project/:projectId/target-setup" element={<TargetSetup />} />

          {/* Phase dashboard */}
          <Route
            path="/project/:projectId/phase/:phaseId"
            element={<PhaseDashboard />}
          />

          {/* Static reference pages */}
          <Route path="/methodology" element={<MethodologyPage />} />

          {/* Catch-all → home */}
          <Route path="*" element={<Navigate to="/" replace />} />
        </Route>
      </Routes>
    </WorkspaceProvider>
    </ToastProvider>
    </ErrorBoundary>
  )
}
