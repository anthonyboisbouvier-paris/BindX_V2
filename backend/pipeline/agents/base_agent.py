"""
DockIt pipeline -- Base LLM Agent v2.0.

Provides the foundation for all 4 AI scientific agents.
Uses OpenAI GPT-4o with JSON mode, two-stage reasoning,
research integration, retry logic, and graceful fallback.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import time
from typing import Optional

logger = logging.getLogger(__name__)

# Agent version for tracking
AGENT_VERSION = "2.1.0"


class BaseAgent:
    """Base class for all AI scientific agents.

    Parameters
    ----------
    name : str
        Human-readable agent name (e.g. "Target Assessment Agent").
    system_prompt : str
        System prompt defining the agent's expertise and output format.
    agent_type : str
        Agent type key for research package routing (target, run_analysis, candidate, optimization).
    model : str
        OpenAI model ID. Default: gpt-4o.
    max_retries : int
        Number of retry attempts on API failure.
    timeout : float
        Request timeout in seconds.
    enable_research : bool
        Whether to gather research data before LLM call.
    two_stage : bool
        Whether to use two-stage LLM reasoning.
    """

    def __init__(
        self,
        name: str,
        system_prompt: str,
        agent_type: str = "target",
        model: str = "gpt-4o",
        max_retries: int = 3,
        timeout: float = 90.0,
        enable_research: bool = True,
        two_stage: bool = True,
    ):
        self.name = name
        self.system_prompt = system_prompt
        self.agent_type = agent_type
        self.model = model
        self.max_retries = max_retries
        self.timeout = timeout
        self.enable_research = enable_research
        self.two_stage = two_stage

    def _get_client(self):
        """Lazy-initialize OpenAI client."""
        api_key = os.environ.get("OPENAI_API_KEY")
        if not api_key:
            return None
        try:
            from openai import OpenAI
            return OpenAI(
                api_key=api_key,
                timeout=self.timeout,
            )
        except ImportError:
            logger.warning("openai package not installed, agent %s unavailable", self.name)
            return None

    def _compute_input_hash(self, context: dict) -> str:
        """Compute deterministic hash of input context."""
        return hashlib.md5(
            json.dumps(context, sort_keys=True, default=str).encode()
        ).hexdigest()[:16]

    def _compute_tools(self, context: dict) -> dict:
        """Run computational tools to inject factual data before LLM call.

        Override in subclasses to register agent-specific tools.
        Returns {} by default (no tools).
        """
        return {}

    def _effective_system_prompt(self, enriched_context: dict) -> str:
        """Build system prompt, appending computed-data instructions if present."""
        prompt = self.system_prompt
        if enriched_context.get("_computed_data"):
            prompt += (
                "\n\n## COMPUTED DATA (GROUND TRUTH)\n"
                "The `_computed_data` field in the input contains FACTUAL COMPUTED RESULTS "
                "from validated RDKit/ADMET/scoring pipeline functions. These are NOT estimates — "
                "they are exact computed values.\n"
                "- NEVER contradict computed values. They override any prior assumptions.\n"
                "- ALWAYS cite specific numbers from _computed_data (e.g., 'MW = 432.5 Da').\n"
                "- Use computed flags and scores as ground truth for your assessment.\n"
                "- If _computed_data shows a PAINS alert or Lipinski violation, you MUST flag it.\n"
            )
        return prompt

    def _gather_research(self, context: dict) -> dict:
        """Gather research data from public APIs based on agent type."""
        if not self.enable_research:
            return {}

        try:
            from pipeline.agents.research import build_research_package
        except ImportError:
            logger.debug("Research module not available, skipping research gathering")
            return {}

        # Extract identifiers from context
        uniprot_id = context.get("uniprot_id") or context.get("target_uniprot_id")
        gene_name = context.get("gene_name") or context.get("target_gene")
        target_name = context.get("target_name") or context.get("protein_name") or context.get("target_protein")

        if not any([uniprot_id, gene_name, target_name]):
            logger.debug("No identifiers found in context for research gathering")
            return {}

        try:
            package = build_research_package(
                agent_type=self.agent_type,
                uniprot_id=uniprot_id,
                gene_name=gene_name,
                target_name=target_name,
            )
            logger.info(
                "Research package for agent '%s': %d papers, %d trials, chembl=%s",
                self.name,
                package.get("research_summary", {}).get("pubmed_papers_found", 0),
                package.get("research_summary", {}).get("clinical_trials_found", 0),
                package.get("research_summary", {}).get("has_chembl_profile", False),
            )
            return package
        except Exception as exc:
            logger.warning("Research gathering failed for agent '%s': %s", self.name, exc)
            return {}

    def _call_llm(self, client, messages: list, temperature: float = 0.2, max_tokens: int = 2000) -> dict:
        """Call the LLM with retry logic. Returns parsed JSON or raises."""
        for attempt in range(1, self.max_retries + 1):
            try:
                response = client.chat.completions.create(
                    model=self.model,
                    messages=messages,
                    response_format={"type": "json_object"},
                    temperature=temperature,
                    max_tokens=max_tokens,
                )
                raw = response.choices[0].message.content
                result = json.loads(raw)

                logger.info(
                    "Agent '%s' LLM call success (tokens: %d/%d, attempt %d)",
                    self.name,
                    response.usage.prompt_tokens if response.usage else 0,
                    response.usage.completion_tokens if response.usage else 0,
                    attempt,
                )
                return result

            except json.JSONDecodeError as jde:
                logger.warning("Agent '%s' invalid JSON (attempt %d): %s", self.name, attempt, jde)
                if attempt == self.max_retries:
                    raise
            except Exception as exc:
                logger.warning("Agent '%s' API error (attempt %d): %s", self.name, attempt, exc)
                if attempt == self.max_retries:
                    raise
                time.sleep(min(2 ** attempt, 10))

        raise RuntimeError("Max retries exceeded")

    async def query(self, context: dict) -> dict:
        """Query the LLM agent with structured context.

        Parameters
        ----------
        context : dict
            Structured input data (agent-specific schema).

        Returns
        -------
        dict
            Structured response with ``available``, ``agent_name``, ``analysis``,
            ``version``, ``timestamp``, ``input_hash``, ``research_summary``.
        """
        input_hash = self._compute_input_hash(context)
        timestamp = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

        client = self._get_client()
        if client is None:
            # Still run computational tools even without API key
            computed_data = self._compute_tools(context)
            resp = self._fallback_response(input_hash, timestamp, "OpenAI API key not configured or openai not installed")
            if computed_data.get("_meta", {}).get("tools_succeeded", 0) > 0:
                resp["_computed_data"] = computed_data
            return resp

        # Step 1: Gather research data
        research = self._gather_research(context)
        research_summary = research.get("research_summary", {})

        # Enrich context with research
        enriched_context = {**context}
        if research.get("research_available"):
            enriched_context["_research_data"] = {
                k: v for k, v in research.items()
                if k not in ("agent_type", "research_available", "research_summary")
            }

        # Step 2: Run computational tools (pipeline-based, no API key needed)
        computed_data = self._compute_tools(enriched_context)
        if computed_data.get("_meta", {}).get("tools_succeeded", 0) > 0:
            enriched_context["_computed_data"] = computed_data
            logger.info(
                "Agent '%s' computed_data: %d/%d tools succeeded",
                self.name,
                computed_data["_meta"]["tools_succeeded"],
                computed_data["_meta"]["tools_run"],
            )

        # Build effective system prompt (includes computed-data instructions if present)
        effective_prompt = self._effective_system_prompt(enriched_context)

        user_message = json.dumps(enriched_context, indent=2, default=str)

        try:
            if self.two_stage:
                # Stage 1: Analysis
                logger.info("Agent '%s' Stage 1: analyzing data (input_hash=%s)", self.name, input_hash)
                stage1_messages = [
                    {"role": "system", "content": effective_prompt},
                    {"role": "user", "content": user_message},
                    {"role": "user", "content": (
                        "STAGE 1: Analyze all the provided data including any research data (_research_data). "
                        "The _computed_data field (if present) contains FACTUAL COMPUTED RESULTS from the pipeline. "
                        "You MUST use these exact values as ground truth — do not estimate or contradict them. "
                        "Identify key insights, patterns, and critical findings. "
                        "Return your analysis as a JSON object with a key 'stage1_analysis' containing your detailed findings."
                    )},
                ]
                stage1_result = self._call_llm(client, stage1_messages, temperature=0.2, max_tokens=2000)

                # Stage 2: Expert recommendation
                logger.info("Agent '%s' Stage 2: generating recommendation", self.name)
                stage2_messages = [
                    {"role": "system", "content": effective_prompt},
                    {"role": "user", "content": user_message},
                    {"role": "assistant", "content": json.dumps(stage1_result, default=str)},
                    {"role": "user", "content": (
                        "STAGE 2: Based on your analysis above, provide your final expert recommendation. "
                        "Follow the exact JSON output schema specified in your system prompt. "
                        "Incorporate insights from the research literature data if available. "
                        "IMPORTANT: All numbers from _computed_data are VERIFIED — cite them exactly. "
                        "Include 'literature_references' (list of {title, pmid, url, finding}) "
                        "and 'key_papers' (list of {title, finding, pmid}) if research papers were provided."
                    )},
                ]
                analysis = self._call_llm(client, stage2_messages, temperature=0.2, max_tokens=4000)

            else:
                # Single-stage (fallback)
                logger.info("Agent '%s' single-stage query (input_hash=%s)", self.name, input_hash)
                messages = [
                    {"role": "system", "content": effective_prompt},
                    {"role": "user", "content": user_message},
                ]
                analysis = self._call_llm(client, messages, temperature=0.3, max_tokens=4000)

            # Inject literature references from research package if not already in analysis
            if research.get("pubmed_papers") and "literature_references" not in analysis:
                analysis["literature_references"] = [
                    {
                        "title": p.get("title", ""),
                        "pmid": p.get("pmid", ""),
                        "url": p.get("url", ""),
                        "finding": (p.get("abstract", "") or "")[:200],
                    }
                    for p in research["pubmed_papers"]
                ]

            return {
                "available": True,
                "agent_name": self.name,
                "model": self.model,
                "analysis": analysis,
                "version": AGENT_VERSION,
                "timestamp": timestamp,
                "input_hash": input_hash,
                "research_summary": research_summary,
            }

        except Exception as exc:
            logger.warning("Agent '%s' query failed: %s", self.name, exc)
            return self._fallback_response(input_hash, timestamp, f"Query failed: {exc}")

    def _fallback_response(self, input_hash: str, timestamp: str, reason: str) -> dict:
        """Generate graceful fallback when OpenAI is unavailable."""
        logger.info("Agent '%s' returning fallback: %s", self.name, reason)
        return {
            "available": False,
            "agent_name": self.name,
            "model": self.model,
            "analysis": None,
            "fallback": reason,
            "version": AGENT_VERSION,
            "timestamp": timestamp,
            "input_hash": input_hash,
            "research_summary": {},
        }

    def query_sync(self, context: dict) -> dict:
        """Synchronous wrapper for query(). Use in non-async contexts."""
        import asyncio
        return asyncio.run(self.query(context))
