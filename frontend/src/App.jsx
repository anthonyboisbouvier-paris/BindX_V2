import React, { lazy, Suspense } from 'react'
import { Routes, Route, Navigate } from 'react-router-dom'
import { WorkspaceProvider } from './contexts/WorkspaceContext.jsx'
import { ToastProvider } from './contexts/ToastContext.jsx'
import { useAuth } from './contexts/AuthContext'

// Layout — always loaded (shell)
import SidebarLayout from './components/SidebarLayout.jsx'
import ErrorBoundary from './components/ErrorBoundary.jsx'
import NavigationProgress from './components/NavigationProgress.jsx'
import BindXLogo from './components/BindXLogo.jsx'

// Lazy-loaded pages (code-split per route)
const LandingPage = lazy(() => import('./pages/LandingPage.jsx'))
const LoginPage = lazy(() => import('./pages/LoginPage.jsx'))
const RegisterPage = lazy(() => import('./pages/RegisterPage.jsx'))
const ProjectListPage = lazy(() => import('./pages/ProjectListPage.jsx'))
const ProjectHome = lazy(() => import('./pages/ProjectHome.jsx'))
const PhaseDashboard = lazy(() => import('./pages/PhaseDashboard/index.jsx'))
const TargetSetup = lazy(() => import('./pages/TargetSetup.jsx'))
const MethodologyPage = lazy(() => import('./components/MethodologyPage.jsx'))
const SurfaceTest = lazy(() => import('./pages/SurfaceTest.jsx'))

// --------------------------------------------------
// ProtectedRoute — redirects to /login if not authenticated
// --------------------------------------------------
function ProtectedRoute({ children }) {
  const { isAuthenticated, loading } = useAuth()

  if (loading) {
    return (
      <div className="min-h-screen bg-bx-surface flex items-center justify-center">
        <BindXLogo variant="splash" size={120} />
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
    <NavigationProgress />
    <Suspense fallback={
      <div className="min-h-screen bg-bx-bg flex items-center justify-center">
        <BindXLogo variant="loading" size={48} />
      </div>
    }>
      <Routes>
        {/* Debug test page (no auth) */}
        <Route path="/surface-test" element={<SurfaceTest />} />

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
    </Suspense>
    </WorkspaceProvider>
    </ToastProvider>
    </ErrorBoundary>
  )
}
