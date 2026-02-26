import React, { createContext, useContext, useState, useEffect, useCallback } from 'react'
import apiClient from '../api'

const AuthContext = createContext(null)

export function AuthProvider({ children }) {
  const [user, setUser] = useState(null)
  const [token, setToken] = useState(() => localStorage.getItem('dockit_token'))
  const [loading, setLoading] = useState(true)

  // Set up axios interceptor for auth header
  useEffect(() => {
    const interceptorId = apiClient.interceptors.request.use((config) => {
      const currentToken = localStorage.getItem('dockit_token')
      if (currentToken) {
        config.headers.Authorization = `Bearer ${currentToken}`
      }
      return config
    })
    return () => apiClient.interceptors.request.eject(interceptorId)
  }, [])

  // Load user on mount if token exists
  useEffect(() => {
    if (token) {
      apiClient.get('/auth/me')
        .then(res => {
          setUser(res.data)
          setLoading(false)
        })
        .catch(() => {
          // Token invalid/expired
          localStorage.removeItem('dockit_token')
          setToken(null)
          setUser(null)
          setLoading(false)
        })
    } else {
      setLoading(false)
    }
  }, [token])

  const login = useCallback(async (email, password) => {
    const res = await apiClient.post('/auth/login', { email, password })
    const { user: userData, token: newToken } = res.data
    localStorage.setItem('dockit_token', newToken)
    setToken(newToken)
    setUser(userData)
    return userData
  }, [])

  const register = useCallback(async (email, username, password) => {
    const res = await apiClient.post('/auth/register', { email, username, password })
    const { user: userData, token: newToken } = res.data
    localStorage.setItem('dockit_token', newToken)
    setToken(newToken)
    setUser(userData)
    return userData
  }, [])

  const logout = useCallback(() => {
    localStorage.removeItem('dockit_token')
    setToken(null)
    setUser(null)
  }, [])

  return (
    <AuthContext.Provider value={{ user, token, loading, login, register, logout, isAuthenticated: !!user }}>
      {children}
    </AuthContext.Provider>
  )
}

export function useAuth() {
  const ctx = useContext(AuthContext)
  if (!ctx) throw new Error('useAuth must be used within AuthProvider')
  return ctx
}

export default AuthContext
