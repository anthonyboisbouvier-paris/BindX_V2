"""
BindX — Legacy V8 endpoints still used by the V9 frontend.

4 endpoints:
- POST /api/preview-target   (TargetSetup page)
- POST /api/detect-pockets    (TargetSetup page)
- POST /api/molecule/analyze-scaffold (MoleculeTable)
- POST /api/agent/{agent_name}/query  (PhaseDashboard)

These will be migrated to proper V9 routers as the corresponding
V9 features are implemented.
"""

from __future__ import annotations

import asyncio
import importlib
import json
import logging
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Optional

import requests as http_requests
from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse

from models import AgentQuery, ScaffoldAnalysisResponse
from pipeline.structure import query_rcsb_pdb_multi

logger = logging.getLogger("dockit.api")

router = APIRouter()

# Thread pools for blocking I/O
_preview_pool = ThreadPoolExecutor(max_workers=4)
_assessment_pool = ThreadPoolExecutor(max_workers=2, thread_name_prefix="assessment")

AGENT_REGISTRY = {
    "target": "pipeline.agents.target_agent.TargetAssessmentAgent",
    "run_analysis": "pipeline.agents.run_analysis_agent.RunAnalysisAgent",
    "candidate": "pipeline.agents.candidate_agent.CandidateEvaluationAgent",
    "optimization": "pipeline.agents.optimization_agent.OptimizationStrategyAgent",
    "chart_advisor": "pipeline.agents.chart_advisor.ChartAdvisorAgent",
}


def shutdown_pools() -> None:
    """Gracefully shut down thread pools (called from app lifespan)."""
    _preview_pool.shutdown(wait=False)
    _assessment_pool.shutdown(wait=False)


# ---------------------------------------------------------------------------
# Helpers — protein / structure lookups
# ---------------------------------------------------------------------------

def _fetch_protein_name(uniprot_id: str) -> Optional[str]:
    """Fetch the recommended protein name from the UniProt REST API."""
    try:
        resp = http_requests.get(
            f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json",
            timeout=(5, 15),
        )
        resp.raise_for_status()
        data = resp.json()
        protein_desc = data.get("proteinDescription", {})
        rec_name = protein_desc.get("recommendedName", {})
        full_name = rec_name.get("fullName", {})
        name = full_name.get("value")
        if name:
            return name
        submitted = protein_desc.get("submittedName", [])
        if submitted and isinstance(submitted, list):
            return submitted[0].get("fullName", {}).get("value")
        return None
    except Exception as exc:
        logger.warning("Failed to fetch protein name for %s: %s", uniprot_id, exc)
        return None


def _fetch_alphafold_info(uniprot_id: str) -> Optional[dict]:
    """Fetch AlphaFold DB prediction info for a UniProt ID."""
    try:
        resp = http_requests.get(
            f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}",
            timeout=(5, 10),
        )
        if resp.status_code != 200:
            return None
        data = resp.json()
        entry = data[0] if isinstance(data, list) else data
        pdb_url = entry.get("pdbUrl")
        if not pdb_url:
            return None
        return {
            "pdb_url": pdb_url,
            "confidence": entry.get("globalMetricValue"),
        }
    except Exception:
        return None


def _fetch_chembl_info(uniprot_id: str) -> dict:
    """Fetch ChEMBL target info and activity counts for a UniProt accession."""
    result: dict = {
        "has_data": False,
        "n_actives": 0,
        "n_with_ic50": 0,
        "target_chembl_id": None,
        "description": "No ChEMBL data available for this target.",
    }
    try:
        target_resp = http_requests.get(
            "https://www.ebi.ac.uk/chembl/api/data/target.json",
            params={
                "target_components__accession": uniprot_id,
                "limit": 1,
                "format": "json",
            },
            timeout=(5, 15),
        )
        target_resp.raise_for_status()
        target_data = target_resp.json()
        targets = target_data.get("targets", [])
        if not targets:
            return result

        target_chembl_id = targets[0].get("target_chembl_id")
        if not target_chembl_id:
            return result
        result["target_chembl_id"] = target_chembl_id

        activity_resp = http_requests.get(
            "https://www.ebi.ac.uk/chembl/api/data/activity.json",
            params={
                "target_chembl_id": target_chembl_id,
                "limit": 0,
                "format": "json",
            },
            timeout=(5, 15),
        )
        activity_resp.raise_for_status()
        activity_data = activity_resp.json()
        total_activities = activity_data.get("page_meta", {}).get("total_count", 0)

        ic50_resp = http_requests.get(
            "https://www.ebi.ac.uk/chembl/api/data/activity.json",
            params={
                "target_chembl_id": target_chembl_id,
                "standard_type": "IC50",
                "limit": 0,
                "format": "json",
            },
            timeout=(5, 15),
        )
        ic50_resp.raise_for_status()
        ic50_data = ic50_resp.json()
        n_ic50 = ic50_data.get("page_meta", {}).get("total_count", 0)

        result["has_data"] = total_activities > 0
        result["n_actives"] = total_activities
        result["n_with_ic50"] = n_ic50

        if total_activities > 0:
            result["description"] = (
                f"{total_activities} known active compounds found in ChEMBL, "
                f"including {n_ic50} with IC50 data. "
                f"ChEMBL is recommended as the primary screening library for this target."
            )
        else:
            result["description"] = (
                "No known active compounds found in ChEMBL for this target. "
                "ZINC drug-like library or de novo generation is recommended."
            )
    except Exception as exc:
        logger.warning("ChEMBL lookup failed for %s: %s", uniprot_id, exc)
        result["description"] = f"ChEMBL lookup failed: {exc}"
    return result


def _fetch_pubchem_info(uniprot_id: str) -> dict:
    """Fetch PubChem bioassay compound count for a UniProt target.

    Resolves UniProt → gene name → PubChem assays → active CID count.
    """
    result: dict = {
        "has_data": False,
        "n_compounds": 0,
        "gene_name": None,
        "description": "No PubChem data available for this target.",
    }
    try:
        # Resolve gene name
        from pipeline.ligands import resolve_gene_name
        gene = resolve_gene_name(uniprot_id)
        if not gene:
            return result
        result["gene_name"] = gene

        # Query PubChem assays for this gene
        pubchem_base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        aids_resp = http_requests.get(
            f"{pubchem_base}/assay/target/genesymbol/{gene}/aids/JSON",
            timeout=(5, 15),
        )
        if aids_resp.status_code == 404:
            result["description"] = f"No PubChem assays found for gene {gene}."
            return result
        aids_resp.raise_for_status()
        aids_data = aids_resp.json()
        all_aids = aids_data.get("IdentifierList", {}).get("AID", [])
        if not all_aids:
            return result

        # Count unique active CIDs across first 10 assays
        seen_cids: set = set()
        for aid in all_aids[:10]:
            try:
                cid_resp = http_requests.get(
                    f"{pubchem_base}/assay/aid/{aid}/cids/JSON?cids_type=active",
                    timeout=(5, 10),
                )
                if cid_resp.status_code != 200:
                    continue
                data = cid_resp.json()
                # PubChem returns CIDs in two possible formats
                cids = data.get("IdentifierList", {}).get("CID", [])
                if not cids:
                    info_list = data.get("InformationList", {}).get("Information", [])
                    if info_list:
                        cids = info_list[0].get("CID", [])
                seen_cids.update(cids)
            except Exception:
                continue

        n = len(seen_cids)
        result["has_data"] = n > 0
        result["n_compounds"] = n
        if n > 0:
            result["description"] = (
                f"{n:,} active compounds found in PubChem bioassays for {gene} "
                f"(across {min(len(all_aids), 10)} assays)."
            )
    except Exception as exc:
        logger.warning("PubChem info lookup failed for %s: %s", uniprot_id, exc)
        result["description"] = f"PubChem lookup failed: {exc}"
    return result


# ---------------------------------------------------------------------------
# Helpers — structure building
# ---------------------------------------------------------------------------

def _build_structures_list(pdb_info, uniprot_id: str) -> list[dict]:
    """Return all available structure options for the preview, in priority order."""
    structures: list[dict] = []

    pdb_list: list[dict] = []
    if isinstance(pdb_info, list):
        pdb_list = pdb_info
    elif isinstance(pdb_info, dict):
        pdb_list = [pdb_info]

    for i, pdb in enumerate(pdb_list):
        pdb_id = pdb.get("pdb_id", "unknown")
        resolution = pdb.get("resolution")
        method = pdb.get("method", "UNKNOWN")
        ligand_id = pdb.get("ligand_id")
        ligand_name = pdb.get("ligand_name")
        source = "pdb_holo" if ligand_id else "pdb_experimental"
        label = "PDB Holo" if ligand_id else "PDB Experimental"
        confidence = 0.98 if ligand_id else 0.90
        desc = f"Experimental {method} structure ({pdb_id})"
        if resolution:
            desc += f" at {resolution} Å"
        if ligand_id:
            lig = f"{ligand_id} ({ligand_name})" if ligand_name else ligand_id
            desc += f" — co-crystallized ligand {lig}"
        structures.append({
            "source": source,
            "label": label,
            "pdb_id": pdb_id,
            "resolution": resolution,
            "method": method,
            "ligand_id": ligand_id,
            "ligand_name": ligand_name,
            "confidence": confidence,
            "recommended": i == 0,
            "explanation": desc,
            "download_url": pdb.get("download_url"),
        })

    af_info = _fetch_alphafold_info(uniprot_id)
    if af_info:
        structures.append({
            "source": "alphafold",
            "label": "AlphaFold",
            "pdb_id": None,
            "resolution": None,
            "method": "AI prediction",
            "ligand_id": None,
            "ligand_name": None,
            "confidence": af_info.get("confidence") or 0.85,
            "recommended": len(structures) == 0,
            "explanation": (
                f"AlphaFold DB predicted structure for {uniprot_id}. "
                "Full-length protein — may differ from PDB which only covers "
                "crystallized domains."
            ),
            "download_url": af_info["pdb_url"],
        })

    structures.append({
        "source": "esmfold",
        "label": "ESMFold",
        "pdb_id": None,
        "resolution": None,
        "method": "ESMFold (on-demand)",
        "ligand_id": None,
        "ligand_name": None,
        "confidence": 0.60,
        "recommended": len(structures) == 0,
        "explanation": (
            "On-demand ESMFold structure prediction. Lower quality than "
            "experimental structures, but always available."
        ),
        "download_url": None,
    })

    return structures


def _build_structure_info(pdb_info: Optional[dict], uniprot_id: str) -> dict:
    """Build the ``structure`` section of the preview response."""
    if pdb_info is not None:
        pdb_id = pdb_info.get("pdb_id", "unknown")
        resolution = pdb_info.get("resolution")
        method = pdb_info.get("method", "UNKNOWN")
        ligand_id = pdb_info.get("ligand_id")
        ligand_name = pdb_info.get("ligand_name")

        explanation_parts = [
            f"Experimental {method} structure from PDB ({pdb_id})"
        ]
        if resolution:
            explanation_parts[0] += f" at {resolution} A resolution"
        if ligand_id:
            lig_desc = ligand_id
            if ligand_name:
                lig_desc = f"{ligand_id} ({ligand_name})"
            explanation_parts.append(f"with co-crystallized ligand {lig_desc}")
        explanation_parts.append("This is the highest quality structure available.")

        return {
            "source": "pdb_experimental",
            "pdb_id": pdb_id,
            "resolution": resolution,
            "method": method,
            "ligand_id": ligand_id,
            "ligand_name": ligand_name,
            "explanation": ". ".join(explanation_parts) + ".",
        }

    af_info = _fetch_alphafold_info(uniprot_id)
    if af_info:
        return {
            "source": "alphafold",
            "pdb_id": None,
            "resolution": None,
            "method": "AI prediction",
            "ligand_id": None,
            "ligand_name": None,
            "download_url": af_info["pdb_url"],
            "explanation": (
                f"No experimental structure found in PDB for {uniprot_id}. "
                "AlphaFold DB has a predicted structure available. "
                "AI-predicted structures are suitable for pocket detection "
                "but may be less accurate for docking near flexible loops."
            ),
        }

    return {
        "source": "none",
        "pdb_id": None,
        "resolution": None,
        "method": None,
        "ligand_id": None,
        "ligand_name": None,
        "explanation": (
            f"No experimental or AlphaFold structure found for {uniprot_id}. "
            "ESMFold on-demand folding will be attempted when the pipeline runs. "
            "This may take longer and the predicted structure may have lower confidence."
        ),
    }


# ---------------------------------------------------------------------------
# Helpers — pocket detection
# ---------------------------------------------------------------------------

def _build_pockets_info(
    pdb_info: Optional[dict],
    uniprot_id: str,
) -> list[dict]:
    """Detect pockets if a PDB structure can be quickly obtained."""
    if pdb_info is None:
        return []
    download_url = pdb_info.get("download_url", "")
    ligand_id = pdb_info.get("ligand_id")
    return _detect_pockets_from_url(download_url, uniprot_id, ligand_id)


def _detect_pockets_from_url(
    download_url: str,
    label: str,
    ligand_id: Optional[str] = None,
    structure_source: str = "",
) -> list[dict]:
    """Download a PDB file from *download_url* and run pocket detection."""
    if not download_url:
        return []

    from pipeline.pockets import detect_pockets

    try:
        with tempfile.TemporaryDirectory(prefix="dockit_preview_") as tmp_dir:
            tmp_path = Path(tmp_dir)
            pdb_path = tmp_path / f"{label}.pdb"

            resp = http_requests.get(download_url, timeout=(5, 30))
            resp.raise_for_status()
            content = resp.text
            if "ATOM" not in content or len(content) < 100:
                logger.warning("Preview PDB download has no ATOM records")
                return []

            pdb_path.write_text(content)

            raw_pockets = detect_pockets(
                pdb_path=pdb_path,
                work_dir=tmp_path,
                ligand_id=ligand_id,
                structure_source=structure_source,
            )

            formatted: list[dict] = []
            for i, pocket in enumerate(raw_pockets[:5]):
                method = pocket.get("method", "unknown")
                probability = pocket.get("probability", 0.0)
                residues = pocket.get("residues", [])
                center = pocket.get("center", (0.0, 0.0, 0.0))

                if method == "co-crystallized_ligand":
                    lig = pocket.get("ligand_id", ligand_id or "unknown")
                    explanation = (
                        f"Pocket defined by co-crystallized ligand {lig}. "
                        "This is an experimentally validated binding site."
                    )
                elif method in ("p2rank", "p2rank_mock"):
                    explanation = (
                        f"Pocket detected by P2Rank with "
                        f"{probability:.0%} probability."
                    )
                elif method == "fpocket":
                    explanation = (
                        f"Pocket detected by fpocket geometric algorithm "
                        f"(druggability score: {probability:.2f})."
                    )
                else:
                    explanation = (
                        f"Pocket detected by {method} method "
                        f"(score: {probability:.2f})."
                    )

                formatted.append({
                    "rank": i + 1,
                    "method": method,
                    "probability": round(probability, 4),
                    "residues_count": len(residues) if isinstance(residues, list) else 0,
                    "residues": residues[:50] if isinstance(residues, list) else [],
                    "center": list(center) if isinstance(center, tuple) else center,
                    "selected": (i == 0),
                    "explanation": explanation,
                    "volume": round(pocket.get("volume", 0.0), 1),
                })

            return formatted

    except Exception as exc:
        logger.warning("Pocket detection failed during preview for %s: %s", label, exc)
        return []


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------

@router.post("/detect-pockets")
async def detect_pockets_endpoint(body: dict) -> dict:
    """Run P2Rank pocket detection on any PDB structure URL."""
    download_url = body.get("download_url", "").strip()
    if not download_url:
        raise HTTPException(400, "download_url is required")

    label = body.get("label", "structure")
    ligand_id = body.get("ligand_id")
    structure_source = body.get("structure_source", "")

    loop = asyncio.get_event_loop()
    try:
        pockets = await loop.run_in_executor(
            _preview_pool, _detect_pockets_from_url, download_url, label, ligand_id, structure_source,
        )
    except Exception as exc:
        logger.warning("Pocket detection endpoint failed: %s", exc)
        pockets = []

    return {"pockets": pockets}


@router.post("/preview-target")
async def preview_target(body: dict) -> dict:
    """Preview protein structure and target information before running the pipeline."""
    uniprot_id = body.get("uniprot_id", "").strip().upper()
    if not uniprot_id:
        raise HTTPException(status_code=400, detail="uniprot_id is required")

    loop = asyncio.get_event_loop()

    protein_name_future = loop.run_in_executor(
        _preview_pool, _fetch_protein_name, uniprot_id,
    )
    pdb_info_future = loop.run_in_executor(
        _preview_pool, query_rcsb_pdb_multi, uniprot_id,
    )
    chembl_future = loop.run_in_executor(
        _preview_pool, _fetch_chembl_info, uniprot_id,
    )
    pubchem_future = loop.run_in_executor(
        _preview_pool, _fetch_pubchem_info, uniprot_id,
    )

    protein_name, pdb_info, chembl_info, pubchem_info = await asyncio.gather(
        protein_name_future,
        pdb_info_future,
        chembl_future,
        pubchem_future,
        return_exceptions=True,
    )

    if isinstance(protein_name, BaseException):
        logger.warning("Protein name lookup raised: %s", protein_name)
        protein_name = None
    if isinstance(pdb_info, BaseException):
        logger.warning("RCSB PDB lookup raised: %s", pdb_info)
        pdb_info = []
    if isinstance(chembl_info, BaseException):
        logger.warning("ChEMBL lookup raised: %s", chembl_info)
        chembl_info = {
            "has_data": False,
            "n_actives": 0,
            "n_with_ic50": 0,
            "target_chembl_id": None,
            "description": f"ChEMBL lookup failed: {chembl_info}",
        }
    if isinstance(pubchem_info, BaseException):
        logger.warning("PubChem lookup raised: %s", pubchem_info)
        pubchem_info = {
            "has_data": False,
            "n_compounds": 0,
            "gene_name": None,
            "description": f"PubChem lookup failed: {pubchem_info}",
        }

    structures_list = _build_structures_list(pdb_info, uniprot_id)
    best_pdb = (
        pdb_info[0]
        if isinstance(pdb_info, list) and pdb_info
        else (pdb_info if isinstance(pdb_info, dict) else None)
    )
    structure_info = (
        structures_list[0]
        if structures_list
        else _build_structure_info(best_pdb, uniprot_id)
    )

    pockets_info: list[dict] = []
    try:
        pockets_info = await loop.run_in_executor(
            _preview_pool, _build_pockets_info, best_pdb, uniprot_id,
        )
    except Exception as exc:
        logger.warning("Pocket detection raised during preview: %s", exc)

    return {
        "uniprot_id": uniprot_id,
        "protein_name": protein_name or uniprot_id,
        "structure": structure_info,
        "structures": structures_list,
        "pockets": pockets_info,
        "chembl_info": chembl_info,
        "pubchem_info": pubchem_info,
    }


@router.post("/molecule/analyze-scaffold")
async def analyze_scaffold_endpoint(body: dict) -> ScaffoldAnalysisResponse:
    """Analyze a molecule's scaffold, R-groups, and return annotated SVG."""
    smiles = body.get("smiles", "")
    if not smiles:
        raise HTTPException(status_code=422, detail="smiles field is required")

    try:
        from pipeline.scaffold_analysis import analyze_scaffold
        result = analyze_scaffold(smiles)
        return ScaffoldAnalysisResponse(**result)
    except Exception as exc:
        logger.error("Scaffold analysis failed for %s: %s", smiles[:60], exc)
        raise HTTPException(status_code=500, detail=f"Scaffold analysis failed: {exc}")


@router.post("/agent/{agent_name}/query")
async def query_agent(agent_name: str, req: AgentQuery) -> JSONResponse:
    """Query a specific AI agent with arbitrary context."""
    if agent_name not in AGENT_REGISTRY:
        raise HTTPException(
            404,
            f"Unknown agent '{agent_name}'. Available: {list(AGENT_REGISTRY.keys())}",
        )

    enriched_context = dict(req.context)
    if req.project_id:
        try:
            from database import get_db
            from models import ProjectORM
            with get_db() as session:
                proj = session.query(ProjectORM).filter_by(id=req.project_id).first()
                if proj:
                    if proj.uniprot_id and "uniprot_id" not in enriched_context:
                        enriched_context["uniprot_id"] = proj.uniprot_id
                    if proj.target_preview_json:
                        try:
                            tp = (
                                json.loads(proj.target_preview_json)
                                if isinstance(proj.target_preview_json, str)
                                else proj.target_preview_json
                            )
                            if tp.get("protein_name") and "target_name" not in enriched_context:
                                enriched_context["target_name"] = tp["protein_name"]
                            assessment = tp.get("assessment", {})
                            if assessment.get("provenance", {}).get("gene_name"):
                                enriched_context["gene_name"] = assessment["provenance"]["gene_name"]
                        except (json.JSONDecodeError, TypeError):
                            pass
        except Exception as proj_exc:
            logger.debug("Could not enrich context from project %s: %s", req.project_id, proj_exc)

    module_path, class_name = AGENT_REGISTRY[agent_name].rsplit(".", 1)

    try:
        mod = importlib.import_module(module_path)
        agent_class = getattr(mod, class_name)
        agent = agent_class()

        loop = asyncio.get_running_loop()
        result = await loop.run_in_executor(
            _assessment_pool,
            lambda: agent.query_sync(enriched_context),
        )

        return JSONResponse(result)

    except Exception as exc:
        logger.error("Agent '%s' query failed: %s", agent_name, exc)
        return JSONResponse({
            "available": False,
            "agent_name": agent_name,
            "fallback": "Agent analysis failed. Please try again.",
        })
