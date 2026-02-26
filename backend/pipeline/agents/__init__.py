"""
DockIt pipeline -- AI Scientific Agents package.

4 specialized LLM agents (OpenAI GPT-4o) providing advisory analysis
at each stage of the drug discovery workflow.
Includes research tools for literature and database integration.
"""

from pipeline.agents.base_agent import BaseAgent
from pipeline.agents.target_agent import TargetAssessmentAgent
from pipeline.agents.run_analysis_agent import RunAnalysisAgent
from pipeline.agents.candidate_agent import CandidateEvaluationAgent
from pipeline.agents.optimization_agent import OptimizationStrategyAgent
from pipeline.agents.research import build_research_package

__all__ = [
    "BaseAgent",
    "TargetAssessmentAgent",
    "RunAnalysisAgent",
    "CandidateEvaluationAgent",
    "OptimizationStrategyAgent",
    "build_research_package",
]
