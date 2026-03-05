/**
 * BindX V9 — User Settings Zustand store.
 *
 * Persisted to localStorage under "bindx-user-settings".
 * Manages column color overrides and Pareto objective overrides.
 *
 * Empty maps on first load → identical behaviour to hardcoded defaults.
 */

import { create } from 'zustand'
import { persist } from 'zustand/middleware'

const useSettingsStore = create(
  persist(
    (set, get) => ({
      // ===== Column color overrides =====
      // { [colKey]: { mode: 'none'|'higher-better'|'lower-better'|'custom', intervals?: [...] } }
      // intervals: [{ min: number|null, max: number|null, color: string }]
      columnColorOverrides: {},

      setColumnColor: (colKey, config) =>
        set(s => ({
          columnColorOverrides: { ...s.columnColorOverrides, [colKey]: config },
        })),

      resetColumnColor: (colKey) =>
        set(s => {
          const next = { ...s.columnColorOverrides }
          delete next[colKey]
          return { columnColorOverrides: next }
        }),

      // ===== Pareto objective overrides =====
      // { [objKey]: { higher_is_better?: boolean, enabled?: boolean } }
      paretoOverrides: {},

      setParetoObjective: (objKey, config) =>
        set(s => ({
          paretoOverrides: { ...s.paretoOverrides, [objKey]: { ...s.paretoOverrides[objKey], ...config } },
        })),

      resetParetoObjective: (objKey) =>
        set(s => {
          const next = { ...s.paretoOverrides }
          delete next[objKey]
          return { paretoOverrides: next }
        }),

      // ===== Analytics charts =====
      // Array of declarative chart configs (JSON, no functions)
      analyticsCharts: [],

      addAnalyticsChart: (config) =>
        set(s => ({
          analyticsCharts: [...s.analyticsCharts, { ...config, id: `chart-${Date.now()}-${Math.random().toString(36).slice(2, 7)}` }],
        })),

      updateAnalyticsChart: (id, updates) =>
        set(s => ({
          analyticsCharts: s.analyticsCharts.map(c => c.id === id ? { ...c, ...updates } : c),
        })),

      removeAnalyticsChart: (id) =>
        set(s => ({
          analyticsCharts: s.analyticsCharts.filter(c => c.id !== id),
        })),

      resetAnalyticsCharts: () =>
        set({ analyticsCharts: [] }),

      // ===== Bulk reset =====
      resetAllSettings: () =>
        set({ columnColorOverrides: {}, paretoOverrides: {}, analyticsCharts: [] }),
    }),
    {
      name: 'bindx-user-settings',
    }
  )
)

export default useSettingsStore
