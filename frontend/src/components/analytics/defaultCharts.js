/**
 * BindX V9 — Default Analytics Charts
 * Only the correlation matrix is shown by default.
 * Use the AI Chart Advisor to generate additional charts.
 */
export function generateDefaultCharts(availableColumns) {
  const numericCount = availableColumns.filter(c => c.type === 'number').length
  if (numericCount >= 3) return [{ type: 'correlation' }]
  return []
}
