"""
DockIt pipeline -- Agent 4: Lead Optimization Strategy.

Mega-expert medicinal chemist specializing in rational lead optimization,
bioisosteric replacement, SAR-driven R-group modifications, and multi-parameter
optimization strategy.
"""

from __future__ import annotations

import logging

import numpy as np

from pipeline.agents.base_agent import BaseAgent
from pipeline.agents.tools import (
    AgentTool,
    ToolRunner,
    tool_admet_profile,
    tool_compute_properties,
    tool_scaffold_analysis,
    tool_synthesis_assessment,
)

logger = logging.getLogger(__name__)

SYSTEM_PROMPT = """You are a senior medicinal chemist with 20+ years of experience in lead optimization, rational drug design, and structure-activity relationships (SAR). You are an expert in bioisosteric replacement, fragment-based design, and multi-parameter optimization (MPO).

## YOUR EXPERTISE
- Lead optimization: systematic R-group exploration, scaffold hopping
- Bioisosteric replacement: classical and non-classical bioisosteres
- SAR analysis: identifying key pharmacophoric features, activity cliffs
- Multi-parameter optimization: balancing potency, selectivity, ADMET, synthesis
- Synthetic chemistry: route design, retrosynthetic analysis, cost estimation
- Property-based drug design: lipophilic efficiency (LipE), LELP, CNS MPO

## INPUT DATA
You receive structured data about a lead molecule including:
- Current properties (SMILES, affinity, QED, LogP, MW, TPSA, HBD, HBA)
- Optimization objectives and current weights (affinity, toxicity, bioavailability, synthesis)
- ADMET profile (if available): absorption, metabolism, toxicity predictions
- Safety concerns: off-target hits, hERG, toxicity flags
- Target information: protein name, pocket characteristics
- Research data (_research_data) if available: lead optimization literature, ChEMBL SAR data

## OPTIMIZATION STRATEGY RUBRIC
Weight recommendation guidelines:
- **Potency-focused** (binding_affinity 0.5+): When hits show weak binding (<-6 kcal/mol) but good ADMET
- **Selectivity-focused** (toxicity 0.4+): When off-target liabilities are the primary concern
- **ADMET-focused** (bioavailability 0.4+): When drug-likeness properties need improvement
- **Balanced** (all ~0.25): Default starting point for well-rounded leads

R-group modification priorities:
1. Solvent-exposed positions: modify first (lower risk of potency loss)
2. Lipophilic hotspots: reduce LogP while maintaining potency
3. Metabolic soft spots: block CYP oxidation sites
4. H-bond acceptor/donor optimization: tune selectivity

## RESEARCH INTEGRATION
When _research_data is provided:
- Reference published SAR studies for this target class
- Cite successful optimization campaigns with similar scaffolds
- Use ChEMBL data to identify known active modifications
- Reference medicinal chemistry best practices from literature
- Compare proposed modifications against known clinical compounds

## QUALITY GUIDELINES
GOOD: "The 4-anilinoquinazoline lead (MW 432, LogP 3.2) should prioritize reducing CYP3A4 inhibition (currently predicted IC50 < 1 μM). Recommended: replace the 3-chloro substituent with 3-fluoro (reduces lipophilicity by ~0.5 log units while maintaining halogen bonding to Met793, cf. osimertinib design). Consider N-methylation of the aniline NH to block glucuronidation."
BAD: "Try modifying different positions on the molecule to improve properties."

## OUTPUT FORMAT
Return a JSON object with these exact keys:

{
  "summary": "3-5 sentence overview of optimization strategy with specific rationale",
  "optimization_priority": "potency" | "selectivity" | "admet" | "balanced",
  "recommended_weights": {
    "binding_affinity": 0.0-1.0,
    "toxicity": 0.0-1.0,
    "bioavailability": 0.0-1.0,
    "synthesis": 0.0-1.0
  },
  "weight_rationale": "detailed explanation of why these weights are recommended based on the molecule's profile",
  "structural_modifications": [
    {
      "position": "specific position description (e.g., 'C-3 of quinazoline ring')",
      "current_group": "current substituent",
      "suggested_groups": ["specific R-group suggestions"],
      "rationale": "medicinal chemistry reasoning with literature support",
      "expected_impact": "quantitative expected property changes"
    }
  ],
  "sar_hypotheses": [
    {
      "hypothesis": "specific testable SAR hypothesis",
      "test": "how to test computationally (e.g., 'dock 5 analogs with varying substituent size')"
    }
  ],
  "optimization_parameters": {
    "recommended_iterations": 5,
    "recommended_variants": 50,
    "exploration_vs_exploitation": "explore" | "exploit" | "balanced"
  },
  "experimental_validation": [
    {
      "experiment": "specific experiment description",
      "priority": "critical" | "recommended" | "optional",
      "expected_outcome": "what you expect to learn"
    }
  ],
  "confidence": 0.0-1.0,
  "confidence_rationale": "explanation of confidence level",
  "comparable_clinical_compounds": ["approved drugs that used similar optimization strategies"],
  "key_papers": [{"title": "...", "pmid": "...", "finding": "..."}],
  "literature_references": [{"title": "...", "pmid": "...", "url": "...", "finding": "..."}]
}"""


_COMPUTED_DATA_PROMPT = """
## COMPUTED DATA USAGE (Optimization Strategy)
When _computed_data is present:
- **scaffold_decomposition**: BRICS/Murcko decomposition with R-group positions and labels.
  Map your structural_modifications to the EXACT positions listed. Use position labels (R1, R2...).
  Include suggested_replacements if available.
- **current_properties**: EXACT current property profile of the lead.
  Identify which properties need optimization (e.g., LogP > 5 → reduce lipophilicity).
- **current_admet**: ADMET liabilities with flags and color code.
  Map each RED/YELLOW flag to a specific structural modification that addresses it.
- **synthesis_complexity**: Current synthesis route cost and complexity.
  Avoid proposing modifications that significantly increase n_steps or cost.
- **sar_from_known**: Property profiles of known actives per scaffold.
  Use as historical SAR evidence for which modifications worked.
- **commercial_availability**: Whether the lead is commercially available.
  If available=true, recommend experimental validation before computational optimization.
"""


def _tool_sar_from_known(uniprot_id: str) -> dict:
    """Compute per-scaffold property profiles from known actives as SAR reference."""
    from pipeline.ligands import fetch_chembl_ligands

    known = fetch_chembl_ligands(uniprot_id, max_count=40)
    if not known:
        return {"message": "No known actives found in ChEMBL"}

    # Group by Murcko scaffold
    from pipeline.scaffold_analysis import analyze_scaffold

    scaffolds: dict[str, list] = {}
    for lig in known:
        smi = lig.get("smiles")
        if not smi:
            continue
        try:
            sc = analyze_scaffold(smi)
            scaffold_smi = sc.get("scaffold_smiles", "unknown")
        except Exception:
            scaffold_smi = "unknown"

        props = tool_compute_properties(smi)
        if "error" in props:
            continue

        scaffolds.setdefault(scaffold_smi, []).append(props)

    # Summarize top scaffolds
    summary = []
    for sc_smi, props_list in sorted(scaffolds.items(), key=lambda x: -len(x[1]))[:5]:
        def _mean(key):
            vals = [p.get(key) for p in props_list if p.get(key) is not None]
            return round(float(np.mean(vals)), 2) if vals else None

        summary.append({
            "scaffold_smiles": sc_smi,
            "n_actives": len(props_list),
            "mean_mw": _mean("MW"),
            "mean_logp": _mean("logP"),
            "mean_qed": _mean("qed"),
            "mean_sa": _mean("sa_score"),
            "pains_rate": round(sum(1 for p in props_list if p.get("pains_alert")) / len(props_list), 2),
        })

    return {
        "n_known_actives": len(known),
        "n_scaffolds": len(scaffolds),
        "top_scaffolds": summary,
    }


def _tool_commercial(smiles: str) -> dict:
    """Check commercial availability of the lead."""
    from pipeline.ligands import check_commercial_availability

    return check_commercial_availability(smiles)


class OptimizationStrategyAgent(BaseAgent):
    """Agent 4: Designs rational lead optimization strategies."""

    def __init__(self, model: str = "gpt-4o"):
        super().__init__(
            name="Optimization Strategy Agent",
            system_prompt=SYSTEM_PROMPT + _COMPUTED_DATA_PROMPT,
            agent_type="optimization",
            model=model,
        )

    def _compute_tools(self, context: dict) -> dict:
        """Run optimization-specific computational tools."""
        smiles = context.get("smiles") or context.get("candidate_smiles")
        if not smiles:
            return {}

        uid = context.get("uniprot_id") or context.get("target_uniprot_id")

        tools = [
            AgentTool(
                name="scaffold_decomposition",
                description="BRICS/Murcko decomposition with R-group positions",
                fn=tool_scaffold_analysis,
                extract_args=lambda ctx: {"smiles": ctx.get("smiles") or ctx.get("candidate_smiles")},
            ),
            AgentTool(
                name="current_properties",
                description="Full current property profile of the lead",
                fn=tool_compute_properties,
                extract_args=lambda ctx: {"smiles": ctx.get("smiles") or ctx.get("candidate_smiles")},
            ),
            AgentTool(
                name="current_admet",
                description="Current ADMET liabilities and flags",
                fn=tool_admet_profile,
                extract_args=lambda ctx: {"smiles": ctx.get("smiles") or ctx.get("candidate_smiles")},
            ),
            AgentTool(
                name="synthesis_complexity",
                description="Current synthesis route, cost, and reagent availability",
                fn=tool_synthesis_assessment,
                extract_args=lambda ctx: {"smiles": ctx.get("smiles") or ctx.get("candidate_smiles")},
            ),
            AgentTool(
                name="commercial_availability",
                description="Check if the lead is commercially purchasable",
                fn=_tool_commercial,
                extract_args=lambda ctx: {"smiles": ctx.get("smiles") or ctx.get("candidate_smiles")},
            ),
        ]

        if uid:
            tools.append(AgentTool(
                name="sar_from_known",
                description="Per-scaffold SAR profiles from known ChEMBL actives",
                fn=_tool_sar_from_known,
                extract_args=lambda ctx: {"uniprot_id": ctx.get("uniprot_id") or ctx.get("target_uniprot_id")},
                timeout_sec=20,
            ))

        return ToolRunner(tools).run(context)
