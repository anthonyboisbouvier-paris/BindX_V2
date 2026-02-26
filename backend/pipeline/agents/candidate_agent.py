"""
DockIt pipeline -- Agent 3: Candidate Scientific Evaluation.

Mega-expert pharmacologist with deep expertise in drug-likeness assessment,
ADMET profiling, safety pharmacology, and clinical compound comparison.
"""

from __future__ import annotations

import logging
from pathlib import Path
from tempfile import mkdtemp

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

SYSTEM_PROMPT = """You are a senior pharmacologist and drug discovery scientist with 20+ years of experience in candidate evaluation, ADMET profiling, and preclinical risk assessment. You hold expertise in medicinal chemistry, pharmacokinetics, and regulatory toxicology.

## YOUR EXPERTISE
- Drug-likeness: Lipinski, Veber, lead-likeness, Pfizer 3/75 rule, GSK 4/400
- ADMET profiling: absorption, distribution, metabolism, excretion, toxicity
- Safety pharmacology: hERG liability, CYP inhibition, reactive metabolites, mutagenicity
- Clinical benchmarking: comparison with approved drugs and clinical candidates
- Risk assessment: benefit-risk evaluation, therapeutic index estimation

## INPUT DATA
You receive structured data about a single molecule candidate including:
- Docking scores (affinity kcal/mol, composite score 0-1, CNN score)
- Physicochemical properties (MW, LogP, QED, TPSA, HBD, HBA, rotatable bonds)
- ADMET predictions (if available): bioavailability, permeability, CYP inhibition, clearance
- Safety flags: off-target hits, hERG risk, toxicity predictions
- Synthesis feasibility: route, steps, cost estimate, SA score
- Confidence metrics
- Research data (_research_data) if available: SAR literature, clinical trials, ChEMBL data

## SCORING INTERPRETATION RUBRIC
Composite Score (0-1):
- **0.8-1.0**: Outstanding candidate — top-tier across all metrics
- **0.6-0.8**: Strong candidate — advance with standard monitoring
- **0.4-0.6**: Moderate candidate — advance with caution, address liabilities
- **0.2-0.4**: Weak candidate — significant optimization needed
- **0.0-0.2**: Poor candidate — likely not viable without major changes

Physicochemical guidelines (Lipinski Ro5):
- MW < 500 Da (ideal: 300-450)
- LogP < 5 (ideal: 1-3)
- HBD <= 5, HBA <= 10
- TPSA 20-130 A² for oral bioavailability
- Rotatable bonds <= 10

ADMET color codes:
- Green: Low risk, favorable properties
- Yellow: Moderate risk, needs monitoring
- Red: High risk, likely liability

## RESEARCH INTEGRATION
When _research_data is provided:
- Compare candidate properties against known SAR for this target
- Reference clinical compounds with similar scaffolds
- Cite relevant ADMET/safety literature
- Identify if the compound falls within known pharmacophore models
- Flag if similar scaffolds have shown clinical failures

## QUALITY GUIDELINES
GOOD: "Candidate CHEMBL123 (MW 432, LogP 2.8, QED 0.72) shows strong binding (-8.9 kcal/mol) with favorable drug-likeness. However, TPSA of 142 A² exceeds the oral bioavailability threshold (130 A²), and the predicted hERG IC50 of 4.2 μM raises moderate cardiac liability concern. The quinazoline scaffold is well-validated in EGFR inhibitors (cf. erlotinib, MW 393)."
BAD: "This molecule has good binding and acceptable properties."

## OUTPUT FORMAT
Return a JSON object with these exact keys:

{
  "summary": "3-5 sentence executive assessment with specific data points",
  "assessment": "advance" | "caution" | "stop",
  "confidence": 0.0-1.0,
  "confidence_rationale": "explanation of confidence level",
  "strengths": [
    {"category": "binding" | "admet" | "selectivity" | "synthesis" | "novelty", "detail": "specific description with numbers"}
  ],
  "risks": [
    {"category": "binding" | "admet" | "selectivity" | "synthesis" | "toxicity", "severity": "high" | "medium" | "low", "detail": "specific risk description"}
  ],
  "drug_likeness": {
    "lipinski_violations": 0,
    "lead_likeness": true,
    "oral_bioavailability_prediction": "favorable" | "moderate" | "poor",
    "key_liabilities": ["specific property concerns"]
  },
  "validation_actions": [
    {
      "action": "specific experiment description",
      "priority": "critical" | "recommended" | "optional",
      "rationale": "scientific justification"
    }
  ],
  "comparable_drugs": ["list of known drugs with similar scaffolds or mechanism"],
  "comparable_clinical_compounds": ["approved/clinical drugs for comparison"],
  "key_papers": [{"title": "...", "pmid": "...", "finding": "..."}],
  "literature_references": [{"title": "...", "pmid": "...", "url": "...", "finding": "..."}]
}"""


_COMPUTED_DATA_PROMPT = """
## COMPUTED DATA USAGE (Candidate Evaluation)
When _computed_data is present:
- **molecular_properties**: EXACT MW, LogP, QED, TPSA, HBD, HBA, PAINS alert, SA score, Lipinski/Veber pass.
  Cite these numbers verbatim in your assessment. Do NOT estimate them.
- **admet_profile**: Full ADMET prediction with hERG IC50, applicability domain, flags, and color code.
  Report every RED flag. Use composite_score for overall ADMET quality.
- **off_target_screening**: SEA broad scan + 10-panel docking results with combined selectivity score.
  Flag any tier1_hits or low combined_selectivity (< 0.5).
- **scaffold_analysis**: Murcko/BRICS decomposition with R-group positions.
  Use scaffold_smiles for structural comparison.
- **synthesis_assessment**: Retrosynthetic route, step count, cost estimate, reagent availability.
  Report n_steps and cost_detail.total_cost_usd.
- **similar_approved_drugs**: Property comparison vs approved drugs targeting the same protein.
  Report deltas (candidate vs approved) for MW, LogP, QED.
"""


def _tool_off_target(smiles: str) -> dict:
    """Run combined off-target screening."""
    from pipeline.off_target import combined_off_target_screening

    work_dir = Path(mkdtemp(prefix="agent_offtarget_"))
    result = combined_off_target_screening(smiles, work_dir)
    return {
        "combined_selectivity": result.get("combined_selectivity"),
        "tier1_hits": result.get("tier1_hits", []),
        "tier2_safe_count": result.get("tier2_safe_count"),
        "sea_total_targets": result.get("sea_results", {}).get("total_targets_scanned", 0),
        "sea_hits": result.get("sea_results", {}).get("n_hits", 0),
    }


def _tool_similar_approved(smiles: str, uniprot_id: str) -> dict:
    """Compare candidate properties vs approved drugs for the same target."""
    from pipeline.ligands import fetch_chembl_ligands

    import numpy as np

    known = fetch_chembl_ligands(uniprot_id, max_count=20)
    if not known:
        return {"message": "No approved/active compounds found in ChEMBL"}

    candidate_props = tool_compute_properties(smiles)
    if "error" in candidate_props:
        return candidate_props

    known_props = []
    for lig in known:
        smi = lig.get("smiles")
        if smi:
            p = tool_compute_properties(smi)
            if "error" not in p:
                known_props.append(p)

    if not known_props:
        return {"n_approved": len(known), "profiles_computed": 0}

    def _mean(lst, key):
        vals = [x.get(key) for x in lst if x.get(key) is not None]
        return round(float(np.mean(vals)), 2) if vals else None

    approved_mw = _mean(known_props, "MW")
    approved_logp = _mean(known_props, "logP")
    approved_qed = _mean(known_props, "qed")

    return {
        "candidate": {
            "MW": candidate_props.get("MW"),
            "logP": candidate_props.get("logP"),
            "qed": candidate_props.get("qed"),
            "tpsa": candidate_props.get("tpsa"),
        },
        "approved_mean": {
            "MW": approved_mw,
            "logP": approved_logp,
            "qed": approved_qed,
        },
        "delta_mw": round(candidate_props.get("MW", 0) - (approved_mw or 0), 2) if approved_mw else None,
        "delta_logp": round(candidate_props.get("logP", 0) - (approved_logp or 0), 2) if approved_logp else None,
        "delta_qed": round(candidate_props.get("qed", 0) - (approved_qed or 0), 2) if approved_qed else None,
        "n_approved_profiled": len(known_props),
    }


class CandidateEvaluationAgent(BaseAgent):
    """Agent 3: In-depth scientific evaluation of individual candidates."""

    def __init__(self, model: str = "gpt-4o"):
        super().__init__(
            name="Candidate Evaluation Agent",
            system_prompt=SYSTEM_PROMPT + _COMPUTED_DATA_PROMPT,
            agent_type="candidate",
            model=model,
        )

    def _compute_tools(self, context: dict) -> dict:
        """Run candidate-specific computational tools."""
        smiles = context.get("smiles") or context.get("candidate_smiles")
        if not smiles:
            return {}

        uid = context.get("uniprot_id") or context.get("target_uniprot_id")

        tools = [
            AgentTool(
                name="molecular_properties",
                description="Full physicochemical profile (MW, LogP, QED, PAINS, SA, Lipinski, Veber)",
                fn=tool_compute_properties,
                extract_args=lambda ctx: {"smiles": ctx.get("smiles") or ctx.get("candidate_smiles")},
            ),
            AgentTool(
                name="admet_profile",
                description="ADMET prediction + hERG + applicability domain",
                fn=tool_admet_profile,
                extract_args=lambda ctx: {"smiles": ctx.get("smiles") or ctx.get("candidate_smiles")},
            ),
            AgentTool(
                name="off_target_screening",
                description="Combined SEA + 10-panel off-target screening",
                fn=_tool_off_target,
                extract_args=lambda ctx: {"smiles": ctx.get("smiles") or ctx.get("candidate_smiles")},
            ),
            AgentTool(
                name="scaffold_analysis",
                description="Murcko/BRICS scaffold decomposition with R-group positions",
                fn=tool_scaffold_analysis,
                extract_args=lambda ctx: {"smiles": ctx.get("smiles") or ctx.get("candidate_smiles")},
            ),
            AgentTool(
                name="synthesis_assessment",
                description="Retrosynthetic route + cost estimate",
                fn=tool_synthesis_assessment,
                extract_args=lambda ctx: {"smiles": ctx.get("smiles") or ctx.get("candidate_smiles")},
            ),
        ]

        if uid:
            tools.append(AgentTool(
                name="similar_approved_drugs",
                description="Property delta vs approved drugs for this target",
                fn=_tool_similar_approved,
                extract_args=lambda ctx: {
                    "smiles": ctx.get("smiles") or ctx.get("candidate_smiles"),
                    "uniprot_id": ctx.get("uniprot_id") or ctx.get("target_uniprot_id"),
                },
                timeout_sec=15,
            ))

        return ToolRunner(tools).run(context)
