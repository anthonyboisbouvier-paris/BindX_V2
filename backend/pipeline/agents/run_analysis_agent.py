"""
DockIt pipeline -- Agent 2: Run Analysis & Lead Selection.

Mega-expert agent specializing in SAR analysis, chemical series identification,
hit rate benchmarking, and screening run quality assessment.
"""

from __future__ import annotations

import copy
import logging

import numpy as np

from pipeline.agents.base_agent import BaseAgent
from pipeline.agents.tools import AgentTool, ToolRunner, tool_compute_properties

logger = logging.getLogger(__name__)

SYSTEM_PROMPT = """You are a senior medicinal chemist with 20+ years of experience in hit-to-lead analysis, chemical series identification, and virtual screening campaign evaluation. You are an expert in SAR (structure-activity relationships), molecular fingerprints, and computational hit assessment.

## YOUR EXPERTISE
- Hit identification: distinguishing true hits from false positives
- Chemical series analysis: scaffold clustering, SAR trends, series prioritization
- Virtual screening benchmarking: expected hit rates by target class and method
- Property assessment: drug-like property distributions, PAINS/assay interference
- Campaign evaluation: enrichment factors, score distributions, method comparison

## INPUT DATA
You receive structured data from a completed virtual screening run including:
- Run parameters (mode, docking engine, ligand count, target)
- Top-scored molecules with properties (affinity, QED, LogP, MW, ADMET, CNN scores)
- Score distributions and summary statistics
- Hit classification counts (by score threshold)
- Research data (_research_data) if available: ChEMBL benchmarks, PubMed literature

## SCORING INTERPRETATION RUBRIC
Affinity scores (kcal/mol, more negative = better):
- **< -10.0**: Exceptional — sub-nanomolar predicted binding
- **-8.0 to -10.0**: Strong — low nanomolar range
- **-6.0 to -8.0**: Moderate — micromolar range, typical for initial hits
- **-4.0 to -6.0**: Weak — marginal binding, likely need optimization
- **> -4.0**: Very weak — likely non-specific

QED (Quantitative Estimate of Drug-likeness, 0-1):
- **> 0.7**: Excellent drug-likeness
- **0.5-0.7**: Good, acceptable for lead optimization
- **0.3-0.5**: Moderate, may need significant optimization
- **< 0.3**: Poor, likely challenging to develop

Hit rate benchmarks by target class:
- Kinases: 0.5-3% typical hit rate
- GPCRs: 0.1-1% typical hit rate
- Nuclear receptors: 1-5% typical hit rate
- Protein-protein interactions: <0.1% typical hit rate
- Novel targets (no known drugs): <0.5% expected

## RESEARCH INTEGRATION
When _research_data is provided:
- Compare hit rate against published benchmarks for similar targets
- Reference known active chemotypes from ChEMBL
- Cite relevant virtual screening methodology papers
- Identify if discovered scaffolds match known bioactive scaffolds

## QUALITY GUIDELINES
GOOD: "Run identified 8 hits (1.6% hit rate) from 500 screened molecules against EGFR kinase. This is within expected range for kinase targets (0.5-3%). The top hit (-9.2 kcal/mol, QED 0.68) shares a quinazoline scaffold with known EGFR inhibitors (erlotinib, gefitinib)."
BAD: "The run produced some good hits with reasonable scores."

## OUTPUT FORMAT
Return a JSON object with these exact keys:

{
  "summary": "3-5 sentence overview of run quality and key findings",
  "run_quality": "excellent" | "good" | "moderate" | "poor",
  "hit_rate_assessment": "detailed assessment of hit rate vs benchmarks for this target class",
  "chemical_series": [
    {
      "name": "descriptive series name (e.g., 'Quinazoline series')",
      "representative_smiles": "SMILES of best example",
      "n_members": 0,
      "avg_affinity": 0.0,
      "strengths": ["specific strengths of this series"],
      "concerns": ["specific concerns"]
    }
  ],
  "top_candidates": [
    {
      "name": "molecule name",
      "rationale": "specific reason this molecule stands out (cite scores)"
    }
  ],
  "property_alerts": ["list of concerning property trends with specifics"],
  "recommended_next_steps": [
    {
      "action": "specific description of recommended action",
      "priority": "high" | "medium" | "low",
      "rationale": "scientific justification with data references"
    }
  ],
  "confidence": 0.0-1.0,
  "confidence_rationale": "explanation of confidence level",
  "comparable_clinical_compounds": ["known drugs with similar scaffolds"],
  "key_papers": [{"title": "...", "pmid": "...", "finding": "..."}],
  "literature_references": [{"title": "...", "pmid": "...", "url": "...", "finding": "..."}]
}"""


_COMPUTED_DATA_PROMPT = """
## COMPUTED DATA USAGE (Run Analysis)
When _computed_data is present:
- **score_distribution**: Exact statistical summary (mean, median, std, percentiles) of affinity/QED/composite scores.
  Use these numbers directly — do NOT re-estimate distribution stats.
- **clustering**: Butina clustering results with cluster IDs and representatives.
  Map clusters to chemical_series in your output. Use cluster_count and sizes as FACT.
- **pareto_analysis**: Mathematically correct Pareto front — cite pareto_front_size and front members.
- **hard_cutoff_analysis**: Exact count of eliminated molecules and reasons. Report these verbatim.
- **known_actives_comparison**: Property comparison of hits vs known ChEMBL actives. Use deltas for assessment.
"""


# -- Tool functions specific to Run Analysis --

def _tool_score_distribution(molecules: list) -> dict:
    """Compute exact score statistics from the run results."""
    from pipeline.scoring import enrich_consensus_detail

    mols = copy.deepcopy(molecules)
    mols = enrich_consensus_detail(mols)

    def _stats(key):
        vals = [m.get(key) for m in mols if m.get(key) is not None]
        if not vals:
            return None
        arr = np.array(vals, dtype=float)
        return {
            "mean": round(float(np.mean(arr)), 3),
            "median": round(float(np.median(arr)), 3),
            "std": round(float(np.std(arr)), 3),
            "min": round(float(np.min(arr)), 3),
            "max": round(float(np.max(arr)), 3),
            "p25": round(float(np.percentile(arr, 25)), 3),
            "p75": round(float(np.percentile(arr, 75)), 3),
            "n": len(vals),
        }

    # Check consensus agreement
    n_consensus = sum(1 for m in mols if m.get("consensus_detail", {}).get("methods_agree", 0) >= 2)

    return {
        "affinity": _stats("affinity"),
        "composite_score": _stats("composite_score"),
        "qed": _stats("qed"),
        "logp": _stats("logP"),
        "mw": _stats("mw"),
        "n_molecules": len(mols),
        "n_consensus_agree": n_consensus,
    }


def _tool_clustering(molecules: list) -> dict:
    """Run Butina clustering on molecules."""
    from pipeline.scoring import cluster_results

    mols = copy.deepcopy(molecules)
    clustered = cluster_results(mols)

    # Summarize clusters
    clusters = {}
    for m in clustered:
        cid = m.get("cluster_id", -1)
        if cid not in clusters:
            clusters[cid] = {"members": [], "representative": None}
        clusters[cid]["members"].append(m.get("name", m.get("smiles", "?")))
        if m.get("is_representative"):
            clusters[cid]["representative"] = m.get("name", m.get("smiles"))

    summary = []
    for cid, data in sorted(clusters.items()):
        summary.append({
            "cluster_id": cid,
            "size": len(data["members"]),
            "representative": data["representative"],
            "members": data["members"][:5],  # Limit for context size
        })

    return {
        "cluster_count": len(clusters),
        "clusters": summary[:15],  # Limit top clusters
    }


def _tool_pareto_analysis(molecules: list) -> dict:
    """Run Pareto ranking on molecules."""
    from pipeline.scoring import pareto_ranking

    mols = copy.deepcopy(molecules)
    ranked = pareto_ranking(mols)

    front = [m for m in ranked if m.get("pareto_front")]
    return {
        "pareto_front_size": len(front),
        "total_molecules": len(ranked),
        "front_members": [
            {"name": m.get("name", m.get("smiles", "?")), "composite_score": m.get("composite_score"),
             "affinity": m.get("affinity")}
            for m in front[:10]
        ],
        "n_ranks": max((m.get("pareto_rank", 0) for m in ranked), default=0) + 1,
    }


def _tool_hard_cutoff_analysis(molecules: list) -> dict:
    """Apply hard cutoffs and report elimination stats."""
    from pipeline.scoring import apply_hard_cutoffs

    mols = copy.deepcopy(molecules)
    passed, eliminated = apply_hard_cutoffs(mols)

    # Count elimination reasons
    reasons = {}
    for m in eliminated:
        reason = m.get("elimination_reason", "unknown")
        reasons[reason] = reasons.get(reason, 0) + 1

    return {
        "total": len(mols),
        "passed": len(passed),
        "eliminated": len(eliminated),
        "elimination_reasons": reasons,
    }


def _tool_known_actives_comparison(uniprot_id: str, molecules: list) -> dict:
    """Compare hit properties against known active compounds from ChEMBL."""
    from pipeline.ligands import fetch_chembl_ligands

    known = fetch_chembl_ligands(uniprot_id, max_count=30)
    if not known:
        return {"message": "No known actives found in ChEMBL"}

    known_props = []
    for lig in known:
        smi = lig.get("smiles")
        if smi:
            props = tool_compute_properties(smi)
            if props and "error" not in props:
                known_props.append(props)

    if not known_props:
        return {"n_known": len(known), "profiles_computed": 0}

    def _mean(lst, key):
        vals = [x.get(key) for x in lst if x.get(key) is not None]
        return round(float(np.mean(vals)), 2) if vals else None

    hit_mw = _mean(molecules, "mw") or _mean(molecules, "MW")
    hit_logp = _mean(molecules, "logP")
    hit_qed = _mean(molecules, "qed")
    known_mw = _mean(known_props, "MW")
    known_logp = _mean(known_props, "logP")
    known_qed = _mean(known_props, "qed")

    return {
        "n_known_actives": len(known),
        "n_known_profiled": len(known_props),
        "known_mean_mw": known_mw,
        "known_mean_logp": known_logp,
        "known_mean_qed": known_qed,
        "hits_mean_mw": hit_mw,
        "hits_mean_logp": hit_logp,
        "hits_mean_qed": hit_qed,
        "delta_mw": round(hit_mw - known_mw, 2) if hit_mw and known_mw else None,
        "delta_logp": round(hit_logp - known_logp, 2) if hit_logp and known_logp else None,
        "delta_qed": round(hit_qed - known_qed, 2) if hit_qed and known_qed else None,
    }


class RunAnalysisAgent(BaseAgent):
    """Agent 2: Analyzes screening runs and recommends lead selection."""

    def __init__(self, model: str = "gpt-4o"):
        super().__init__(
            name="Run Analysis Agent",
            system_prompt=SYSTEM_PROMPT + _COMPUTED_DATA_PROMPT,
            agent_type="run_analysis",
            model=model,
        )

    def _compute_tools(self, context: dict) -> dict:
        """Run analysis-specific computational tools."""
        molecules = context.get("molecules")
        if not molecules or not isinstance(molecules, list):
            return {}

        uid = context.get("uniprot_id") or context.get("target_uniprot_id")

        tools = [
            AgentTool(
                name="score_distribution",
                description="Exact score statistics (mean/median/std/percentiles)",
                fn=_tool_score_distribution,
                extract_args=lambda ctx: {"molecules": ctx["molecules"]},
            ),
            AgentTool(
                name="clustering",
                description="Butina fingerprint clustering into chemical series",
                fn=_tool_clustering,
                extract_args=lambda ctx: {"molecules": ctx["molecules"]},
            ),
            AgentTool(
                name="pareto_analysis",
                description="Multi-objective Pareto front ranking",
                fn=_tool_pareto_analysis,
                extract_args=lambda ctx: {"molecules": ctx["molecules"]},
            ),
            AgentTool(
                name="hard_cutoff_analysis",
                description="Hard cutoff elimination counts and reasons",
                fn=_tool_hard_cutoff_analysis,
                extract_args=lambda ctx: {"molecules": ctx["molecules"]},
            ),
        ]

        if uid:
            tools.append(AgentTool(
                name="known_actives_comparison",
                description="Compare hit properties vs known ChEMBL actives",
                fn=_tool_known_actives_comparison,
                extract_args=lambda ctx: {
                    "uniprot_id": ctx.get("uniprot_id") or ctx.get("target_uniprot_id"),
                    "molecules": ctx["molecules"],
                },
                timeout_sec=15,
            ))

        return ToolRunner(tools).run(context)
