"""
BindX pipeline -- AI Chart Advisor Agent.

Lightweight single-stage agent that generates analytics chart configurations
from a user's natural-language prompt and the available dataset columns.
Uses gpt-4o-mini for fast, low-cost responses.
"""

from __future__ import annotations

import logging

from pipeline.agents.base_agent import BaseAgent

logger = logging.getLogger(__name__)

SYSTEM_PROMPT = """\
You are an expert data-visualization assistant for a drug-discovery platform called BindX.
The user describes what they want to see, and you return the best chart configurations.

## AVAILABLE CHART TYPES

Each chart is a JSON object with a `type` field and type-specific keys:

1. **scatter** — 2D scatter plot
   `{ "type": "scatter", "xKey": "<col>", "yKey": "<col>", "colorKey": "<col>" (optional) }`

2. **histogram** — Distribution of a numeric column
   `{ "type": "histogram", "key": "<col>", "bins": 20 }`

3. **bar** — Category counts (categorical column)
   `{ "type": "bar", "key": "<col>" }`

4. **pie** — Category proportions (categorical column)
   `{ "type": "pie", "key": "<col>" }`

5. **box** — Box plot of a metric grouped by a category
   `{ "type": "box", "metricKey": "<col>", "groupKey": "<col>" }`

6. **correlation** — Correlation matrix (all numeric columns or a subset)
   `{ "type": "correlation", "corrColumns": ["<col>", ...] }` (corrColumns optional)

7. **bubble** — Bubble chart (3 numeric axes + optional color)
   `{ "type": "bubble", "xKey": "<col>", "yKey": "<col>", "sizeKey": "<col>", "colorKey": "<col>" (optional) }`

8. **topn** — Top N ranked molecules by a metric
   `{ "type": "topn", "key": "<col>", "n": 15 }`

## RULES

- Only use column keys that exist in the provided `columns` list.
- Return between 1 and 6 charts — pick only the most relevant ones for the user's question.
- For docking scores: lower (more negative) is better.
- For QED, composite_score, confidence: higher is better.
- Prefer scatter plots for relationships, histograms for distributions, bar/pie for categories.
- If the user asks about correlations, use a correlation chart with the relevant subset of columns.
- If the prompt is vague, pick a reasonable mix of 2-3 charts.

## OUTPUT FORMAT

Return ONLY valid JSON:
```json
{
  "charts": [ ... ],
  "explanation": "Brief explanation of why these charts answer the user's question."
}
```
"""


class ChartAdvisorAgent(BaseAgent):
    """Generates analytics chart configs from a natural-language prompt."""

    def __init__(self):
        super().__init__(
            name="Chart Advisor",
            system_prompt=SYSTEM_PROMPT,
            agent_type="chart_advisor",
            model="gpt-4o-mini",
            max_retries=2,
            timeout=30.0,
            enable_research=False,
            two_stage=False,
        )

    async def query(self, context: dict) -> dict:
        """Override to return a simpler response shape for chart configs."""
        import time

        input_hash = self._compute_input_hash(context)
        timestamp = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

        client = self._get_client()
        if client is None:
            return {
                "available": False,
                "agent_name": self.name,
                "error": "OpenAI API key not configured",
                "charts": [],
                "explanation": "",
            }

        # Build a focused user message
        prompt = context.get("prompt", "")
        columns = context.get("columns", [])
        molecule_count = context.get("molecule_count", 0)

        user_message = (
            f"User prompt: {prompt}\n\n"
            f"Dataset: {molecule_count} molecules\n\n"
            f"Available columns:\n"
        )
        for col in columns:
            user_message += f"- {col['key']} ({col.get('type', 'string')}) — {col.get('label', col['key'])}"
            if col.get('group'):
                user_message += f" [group: {col['group']}]"
            user_message += "\n"

        try:
            messages = [
                {"role": "system", "content": self.system_prompt},
                {"role": "user", "content": user_message},
            ]
            result = self._call_llm(client, messages, temperature=0.3, max_tokens=1500)

            charts = result.get("charts", [])
            explanation = result.get("explanation", "")

            # Validate that all chart keys reference real columns
            col_keys = {c["key"] for c in columns}
            validated = []
            for chart in charts:
                valid = True
                for k in ("key", "xKey", "yKey", "sizeKey", "colorKey", "metricKey", "groupKey"):
                    val = chart.get(k)
                    if val and val not in col_keys:
                        valid = False
                        break
                corr_cols = chart.get("corrColumns")
                if corr_cols:
                    chart["corrColumns"] = [c for c in corr_cols if c in col_keys]
                if valid:
                    validated.append(chart)

            return {
                "available": True,
                "agent_name": self.name,
                "charts": validated[:6],
                "explanation": explanation,
            }

        except Exception as exc:
            logger.warning("ChartAdvisor query failed: %s", exc)
            return {
                "available": False,
                "agent_name": self.name,
                "error": str(exc),
                "charts": [],
                "explanation": "",
            }
