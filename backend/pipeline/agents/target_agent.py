"""
DockIt pipeline -- Agent 1: Target Assessment & Campaign Strategy.

Mega-expert agent with 20+ years of computational drug discovery expertise.
Interprets assessment scores, recommends GO/CAUTION/NO-GO,
proposes screening strategy, and integrates literature findings.
"""

from __future__ import annotations

import logging

from pipeline.agents.base_agent import BaseAgent
from pipeline.agents.tools import AgentTool, ToolRunner, tool_compute_properties

logger = logging.getLogger(__name__)

SYSTEM_PROMPT = """You are a world-class computational drug discovery scientist with 20+ years of experience in target identification, validation, and campaign strategy. You hold deep expertise in structural biology, pharmacology, and translational medicine.

## YOUR EXPERTISE
- Target validation: genetic association, clinical evidence, mechanism of action
- Structural druggability: pocket analysis, allosteric sites, cryptic binding sites
- Competitive landscape: approved drugs, clinical pipeline, patent status
- Safety pharmacology: off-target liability, essential gene analysis, tissue expression
- Screening strategy: library design, hit rate prediction, cascade planning

## INPUT DATA
You receive structured data about a therapeutic target including:
- 5 assessment scores (evidence 0-1, druggability 0-1, novelty 0-1, safety 0-1, feasibility 0-1)
- Protein information (name, class, structure source, pocket data)
- Disease context and ChEMBL statistics
- Research data (_research_data) if available: PubMed papers, ChEMBL profile, UniProt info, clinical trials, drug-gene interactions

## SCORING INTERPRETATION RUBRIC
For each score dimension (0-1 scale):
- **0.8-1.0**: Excellent — strong evidence, high confidence
- **0.6-0.8**: Good — solid data, some gaps remain
- **0.4-0.6**: Moderate — mixed evidence, key uncertainties
- **0.2-0.4**: Weak — limited data, significant concerns
- **0.0-0.2**: Poor — very limited or contradictory evidence

## RESEARCH INTEGRATION
When _research_data is provided:
- Reference specific PubMed papers by PMID to support your claims
- Use ChEMBL bioactivity data to assess druggability and competitive landscape
- Cite clinical trials to evaluate translational potential
- Use UniProt functional data to inform mechanism-based reasoning
- Reference drug-gene interactions for safety context

## QUALITY GUIDELINES
GOOD example: "The EGFR kinase domain (PDB: 1M17, 1.8A resolution) shows a well-defined ATP-binding pocket with druggability score 0.82. ChEMBL reports 15,247 bioactivities with 4 approved drugs (erlotinib, gefitinib, osimertinib, lapatinib), confirming clinical tractability."
BAD example: "This target looks druggable and has some known drugs."

Be specific. Cite data. Quantify when possible. Avoid vague statements.

## OUTPUT FORMAT
Return a JSON object with these exact keys:

{
  "summary": "3-5 sentence executive summary with specific data points",
  "recommendation": "GO" | "CAUTION" | "NO-GO",
  "confidence": 0.0-1.0,
  "confidence_rationale": "Explain why confidence is at this level",
  "strengths": ["list of specific strengths with data references"],
  "weaknesses": ["list of specific weaknesses with data references"],
  "screening_strategy": {
    "recommended_mode": "rapid" | "standard" | "deep",
    "library_focus": "detailed description of recommended library composition",
    "expected_hit_rate": "low (<0.1%)" | "moderate (0.1-1%)" | "high (>1%)",
    "key_considerations": ["list of important strategic factors"]
  },
  "risk_mitigation": ["list of specific actions to mitigate identified risks"],
  "next_steps": ["ordered list of recommended next actions with rationale"],
  "comparable_clinical_compounds": ["list of approved/clinical drugs targeting this protein"],
  "key_papers": [{"title": "...", "pmid": "...", "finding": "key finding relevant to assessment"}],
  "literature_references": [{"title": "...", "pmid": "...", "url": "...", "finding": "..."}]
}"""


_COMPUTED_DATA_PROMPT = """
## COMPUTED DATA USAGE (Target Assessment)
When _computed_data is present:
- **target_assessment_scores**: Use these 5 scores as GROUND TRUTH. Do NOT re-estimate them.
  Report the exact composite_score and per-dimension scores. Base your GO/CAUTION/NO-GO on these.
- **known_ligands_profile**: This is the EXACT property distribution of known active ligands from ChEMBL.
  Use mean/median MW, LogP, QED as benchmarks for druggability assessment.
- **ligand_strategy**: This is the data-driven recommendation for library composition.
  Incorporate it into your screening_strategy.
"""


# -- Tool functions specific to Target Assessment --

def _tool_known_ligands_profile(uniprot_id: str) -> dict:
    """Fetch known ligands and compute their property distribution."""
    from pipeline.ligands import fetch_chembl_ligands

    ligands = fetch_chembl_ligands(uniprot_id, max_count=50)
    if not ligands:
        return {"n_ligands": 0, "message": "No known ligands found in ChEMBL"}

    profiles = []
    for lig in ligands:
        smi = lig.get("smiles")
        if not smi:
            continue
        props = tool_compute_properties(smi)
        if props and "error" not in props:
            profiles.append(props)

    if not profiles:
        return {"n_ligands": len(ligands), "profiles_computed": 0}

    import numpy as np

    def _stats(values):
        arr = [v for v in values if v is not None]
        if not arr:
            return None
        return {"mean": round(float(np.mean(arr)), 2), "median": round(float(np.median(arr)), 2),
                "std": round(float(np.std(arr)), 2), "min": round(float(np.min(arr)), 2),
                "max": round(float(np.max(arr)), 2)}

    return {
        "n_ligands": len(ligands),
        "profiles_computed": len(profiles),
        "mw": _stats([p.get("MW") for p in profiles]),
        "logp": _stats([p.get("logP") for p in profiles]),
        "qed": _stats([p.get("qed") for p in profiles]),
        "tpsa": _stats([p.get("tpsa") for p in profiles]),
        "hbd": _stats([p.get("hbd") for p in profiles]),
        "hba": _stats([p.get("hba") for p in profiles]),
        "pains_alert_count": sum(1 for p in profiles if p.get("pains_alert")),
        "lipinski_pass_rate": round(sum(1 for p in profiles if p.get("lipinski_pass")) / len(profiles), 2),
    }


def _tool_target_assessment_scores(uniprot_id: str, disease_context: str = None) -> dict:
    """Run the full target assessment pipeline and return scores."""
    from pipeline.target_assessment import assess_target

    result = assess_target(uniprot_id, disease_context=disease_context)
    # Return only the scores + recommendation, not the full verbose output
    return {
        "composite_score": result.get("composite_score"),
        "recommendation": result.get("recommendation"),
        "scores": result.get("scores"),
        "flags": result.get("flags", []),
        "rationale": result.get("rationale"),
    }


def _tool_ligand_strategy(uniprot_id: str) -> dict:
    """Data-driven ligand strategy recommendation."""
    from pipeline.ligands import auto_select_ligand_strategy

    return auto_select_ligand_strategy(uniprot_id)


class TargetAssessmentAgent(BaseAgent):
    """Agent 1: Interprets target assessment scores and recommends strategy."""

    def __init__(self, model: str = "gpt-4o"):
        super().__init__(
            name="Target Assessment Agent",
            system_prompt=SYSTEM_PROMPT + _COMPUTED_DATA_PROMPT,
            agent_type="target",
            model=model,
        )

    def _compute_tools(self, context: dict) -> dict:
        """Run target-specific computational tools."""
        uid = context.get("uniprot_id") or context.get("target_uniprot_id")
        if not uid:
            return {}

        disease = context.get("disease_context")

        tools = [
            AgentTool(
                name="known_ligands_profile",
                description="Property distribution of known ChEMBL active ligands",
                fn=_tool_known_ligands_profile,
                extract_args=lambda ctx: {"uniprot_id": uid},
                timeout_sec=15,
            ),
            AgentTool(
                name="target_assessment_scores",
                description="5-dimension target assessment scores with provenance",
                fn=_tool_target_assessment_scores,
                extract_args=lambda ctx: {"uniprot_id": uid, "disease_context": disease},
                timeout_sec=30,
            ),
            AgentTool(
                name="ligand_strategy",
                description="Data-driven ligand sourcing strategy recommendation",
                fn=_tool_ligand_strategy,
                extract_args=lambda ctx: {"uniprot_id": uid},
                timeout_sec=10,
            ),
        ]
        return ToolRunner(tools).run(context)
