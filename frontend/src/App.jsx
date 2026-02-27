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

// Main pages
import ProjectListPage from './pages/ProjectListPage.jsx'
import ProjectHome from './pages/ProjectHome.jsx'

// Phase dashboard
import PhaseDashboard from './pages/PhaseDashboard.jsx'

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
        <div className="w-8 h-8 border-2 border-[#1e3a5f] border-t-transparent rounded-full animate-spin" />
      </div>
    )
  }

  if (!isAuthenticated) {
    return <Navigate to="/login" replace />
  }

  return children
}

// --------------------------------------------------
// App — V9 routing
// --------------------------------------------------
export default function App() {
  return (
    <ErrorBoundary>
    <WorkspaceProvider>
    <ToastProvider>
      <Routes>
        {/* Auth pages (no sidebar, no auth required) */}
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
    </ToastProvider>
    </WorkspaceProvider>
    </ErrorBoundary>
  )
}
