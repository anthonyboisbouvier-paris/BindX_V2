"""
DockIt pipeline -- Agent Computational Tools.

Provides a ToolRunner that executes pipeline functions in parallel,
injecting factual computed results into agent context before LLM calls.
This grounds LLM reasoning with verified numbers from RDKit/ADMET/scoring.
"""

from __future__ import annotations

import copy
import logging
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from typing import Any, Callable, Optional

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class AgentTool:
    """Descriptor for a single computational tool available to an agent."""

    name: str
    description: str
    fn: Callable[..., Any]
    extract_args: Callable[[dict], dict]
    timeout_sec: int = 30


class ToolRunner:
    """Execute a list of AgentTools in parallel and collect results.

    Each tool runs independently — failures in one tool do not block others.
    Returns a dict mapping tool names to their results, plus ``_meta`` stats.
    """

    def __init__(self, tools: list[AgentTool], max_workers: int = 4):
        self.tools = tools
        self.max_workers = max_workers

    def run(self, context: dict) -> dict:
        """Run all tools in parallel against the given context.

        Parameters
        ----------
        context : dict
            Agent context from which tool arguments are extracted.

        Returns
        -------
        dict
            ``{tool_name: result, ..., _meta: {tools_run, tools_succeeded, elapsed_sec}}``
        """
        results: dict[str, Any] = {}
        t0 = time.monotonic()
        tools_run = 0
        tools_succeeded = 0

        with ThreadPoolExecutor(max_workers=self.max_workers) as pool:
            futures = {}
            for tool in self.tools:
                try:
                    kwargs = tool.extract_args(context)
                except Exception as exc:
                    logger.debug("Tool '%s' arg extraction failed: %s", tool.name, exc)
                    continue
                if kwargs is None:
                    continue
                futures[pool.submit(self._run_single, tool, kwargs)] = tool.name
                tools_run += 1

            for future in as_completed(futures):
                name = futures[future]
                try:
                    result = future.result(timeout=60)
                    if result is not None:
                        results[name] = result
                        tools_succeeded += 1
                except Exception as exc:
                    logger.warning("Tool '%s' failed: %s", name, exc)

        elapsed = round(time.monotonic() - t0, 2)
        results["_meta"] = {
            "tools_run": tools_run,
            "tools_succeeded": tools_succeeded,
            "elapsed_sec": elapsed,
        }
        logger.info(
            "ToolRunner: %d/%d tools succeeded in %.1fs",
            tools_succeeded, tools_run, elapsed,
        )
        return results

    @staticmethod
    def _run_single(tool: AgentTool, kwargs: dict) -> Any:
        """Execute a single tool with timeout."""
        t0 = time.monotonic()
        try:
            result = tool.fn(**kwargs)
            elapsed = round(time.monotonic() - t0, 2)
            logger.debug("Tool '%s' completed in %.1fs", tool.name, elapsed)
            return result
        except Exception as exc:
            logger.warning("Tool '%s' raised: %s", tool.name, exc)
            return None


# ---------------------------------------------------------------------------
# Shared tool functions — reused by Candidate + Optimization agents
# ---------------------------------------------------------------------------

def tool_compute_properties(smiles: str) -> dict:
    """Compute full molecular property profile including drug-likeness rules.

    Wraps scoring.compute_properties + PAINS + SA score + Lipinski/Veber counts.
    """
    from pipeline.scoring import compute_properties, compute_pains_alert, compute_sa_score

    props = compute_properties(smiles)
    if not props or props.get("MW") is None:
        return {"error": "Could not compute properties", "smiles": smiles}

    pains = compute_pains_alert(smiles)
    sa_score = compute_sa_score(smiles)

    # Lipinski violations
    lipinski_violations = 0
    if (props.get("MW") or 0) > 500:
        lipinski_violations += 1
    if (props.get("logP") or 0) > 5:
        lipinski_violations += 1
    if (props.get("hbd") or 0) > 5:
        lipinski_violations += 1
    if (props.get("hba") or 0) > 10:
        lipinski_violations += 1

    # Veber violations
    veber_violations = 0
    if (props.get("tpsa") or 0) > 140:
        veber_violations += 1
    if (props.get("rotatable_bonds") or 0) > 10:
        veber_violations += 1

    return {
        **props,
        "smiles": smiles,
        "pains_alert": pains,
        "sa_score": sa_score,
        "lipinski_violations": lipinski_violations,
        "veber_violations": veber_violations,
        "lipinski_pass": lipinski_violations <= 1,
        "veber_pass": veber_violations == 0,
    }


def tool_admet_profile(smiles: str) -> dict:
    """Full ADMET profile with hERG specialized + applicability domain check.

    Wraps admet.predict_admet + compute_admet_composite + predict_herg_specialized
    + check_applicability_domain.
    """
    from pipeline.admet import (
        predict_admet,
        predict_herg_specialized,
        check_applicability_domain,
    )

    admet_list = predict_admet([smiles])
    admet = admet_list[0] if admet_list else {}

    herg = predict_herg_specialized(smiles)
    domain = check_applicability_domain(smiles)

    # Remove verbose sub-dicts for LLM context, keep summary
    return {
        "smiles": smiles,
        "composite_score": admet.get("composite_score"),
        "color_code": admet.get("color_code"),
        "flags": admet.get("flags", []),
        "absorption": admet.get("absorption"),
        "metabolism": admet.get("metabolism"),
        "toxicity": admet.get("toxicity"),
        "herg": herg,
        "applicability_domain": domain,
    }


def tool_scaffold_analysis(smiles: str) -> dict:
    """Scaffold decomposition without SVG (too large for LLM context).

    Wraps scaffold_analysis.analyze_scaffold, strips annotated_svg.
    """
    from pipeline.scaffold_analysis import analyze_scaffold

    result = analyze_scaffold(smiles)
    # Remove SVG — too large for LLM context
    result.pop("annotated_svg", None)
    return result


def tool_synthesis_assessment(smiles: str) -> dict:
    """Retrosynthetic route + cost estimation.

    Wraps retrosynthesis.plan_synthesis + estimate_synthesis_cost.
    """
    from pipeline.retrosynthesis import plan_synthesis, estimate_synthesis_cost

    route = plan_synthesis(smiles, max_depth=4, timeout_sec=30)
    cost = estimate_synthesis_cost(route)

    return {
        "smiles": smiles,
        "n_steps": route.get("n_steps"),
        "confidence": route.get("confidence"),
        "all_reagents_available": route.get("all_reagents_available"),
        "estimated_cost": route.get("estimated_cost"),
        "cost_detail": cost,
        "steps_summary": [
            {
                "reaction": s.get("reaction_type", s.get("reaction", "unknown")),
                "n_reactants": len(s.get("reactants", [])),
            }
            for s in route.get("steps", [])
        ],
    }


# ---------------------------------------------------------------------------
# Arg extractors — pull the right fields from context for each tool
# ---------------------------------------------------------------------------

def _extract_smiles(context: dict) -> Optional[dict]:
    """Extract SMILES from context."""
    smiles = context.get("smiles") or context.get("candidate_smiles")
    if not smiles:
        return None
    return {"smiles": smiles}


def _extract_uniprot(context: dict) -> Optional[dict]:
    """Extract UniProt ID from context."""
    uid = context.get("uniprot_id") or context.get("target_uniprot_id")
    if not uid:
        return None
    return {"uniprot_id": uid}


def _extract_smiles_and_uniprot(context: dict) -> Optional[dict]:
    """Extract both SMILES and UniProt ID."""
    smiles = context.get("smiles") or context.get("candidate_smiles")
    uid = context.get("uniprot_id") or context.get("target_uniprot_id")
    if not smiles or not uid:
        return None
    return {"smiles": smiles, "uniprot_id": uid}
