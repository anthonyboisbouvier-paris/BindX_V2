/**
 * BindX V9 — Analytics Chart Registry
 * Maps chart type strings to React components.
 */
import HistogramChart from './HistogramChart.jsx'
import CategoryBarChart from './CategoryBarChart.jsx'
import PieBreakdown from './PieBreakdown.jsx'
import ScatterPlot from './ScatterPlot.jsx'
import BoxPlotChart from './BoxPlotChart.jsx'
import CorrelationMatrix from './CorrelationMatrix.jsx'
import BubbleChart from './BubbleChart.jsx'
import TopNRanked from './TopNRanked.jsx'

const CHART_REGISTRY = {
  histogram:   { component: HistogramChart,    label: 'Histogram' },
  bar:         { component: CategoryBarChart,  label: 'Bar Chart' },
  pie:         { component: PieBreakdown,      label: 'Pie Chart' },
  scatter:     { component: ScatterPlot,       label: 'Scatter' },
  box:         { component: BoxPlotChart,      label: 'Box Plot' },
  correlation: { component: CorrelationMatrix, label: 'Correlation' },
  bubble:      { component: BubbleChart,       label: 'Bubble' },
  topn:        { component: TopNRanked,        label: 'Top N' },
}

export default CHART_REGISTRY

export const CHART_TYPES = Object.entries(CHART_REGISTRY).map(([type, { label }]) => ({ type, label }))
