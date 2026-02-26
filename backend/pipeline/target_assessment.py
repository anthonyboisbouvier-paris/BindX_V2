"""
DockIt pipeline -- Target Assessment Engine.

Evaluates the scientific quality of a therapeutic target BEFORE screening.
Produces 5 scores (evidence, druggability, novelty, safety, feasibility)
and an overall GO / CAUTION / NO-GO recommendation.

Reuses existing pipeline modules (ligands, pockets, structure, off_target)
plus new Open Targets API integration.
"""

from __future__ import annotations

import hashlib
import json
import logging
import time
from typing import Optional

import requests

logger = logging.getLogger(__name__)

HTTP_TIMEOUT = (10, 30)

# Default weights for score aggregation
DEFAULT_WEIGHTS = {
    "evidence": 0.25,
    "druggability": 0.30,
    "novelty": 0.15,
    "safety": 0.15,
    "feasibility": 0.15,
}


# ---------------------------------------------------------------------------
# Open Targets helpers
# ---------------------------------------------------------------------------

def _fetch_open_targets_association(uniprot_id: str, disease_context: Optional[str] = None) -> dict:
    """Fetch target-disease association data from Open Targets Platform.

    Parameters
    ----------
    uniprot_id : str
        UniProt accession (e.g. "P00533").
    disease_context : str, optional
        EFO disease ID or free-text disease name.

    Returns
    -------
    dict
        Keys: ``overall_score``, ``n_associations``, ``top_diseases``, ``has_data``.
    """
    try:
        # Open Targets GraphQL API
        query = """
        query TargetAssociations($ensemblId: String!) {
          target(ensemblId: $ensemblId) {
            id
            approvedSymbol
            associatedDiseases(page: {index: 0, size: 10}) {
              count
              rows {
                disease { id name }
                score
                datatypeScores { componentId score }
              }
            }
          }
        }
        """
        # First resolve UniProt to Ensembl ID via UniProt API
        ensembl_id = _resolve_ensembl_id(uniprot_id)
        if not ensembl_id:
            return {"overall_score": 0.0, "n_associations": 0, "top_diseases": [], "has_data": False}

        resp = requests.post(
            "https://api.platform.opentargets.org/api/v4/graphql",
            json={"query": query, "variables": {"ensemblId": ensembl_id}},
            timeout=HTTP_TIMEOUT,
        )
        if resp.status_code != 200:
            logger.warning("Open Targets API returned %d", resp.status_code)
            return {"overall_score": 0.0, "n_associations": 0, "top_diseases": [], "has_data": False}

        data = resp.json().get("data", {}).get("target", {})
        if not data:
            return {"overall_score": 0.0, "n_associations": 0, "top_diseases": [], "has_data": False}

        assoc = data.get("associatedDiseases", {})
        n_associations = assoc.get("count", 0)
        rows = assoc.get("rows", [])
        top_diseases = [
            {"name": r["disease"]["name"], "score": r["score"]}
            for r in rows[:5]
        ]

        # If disease_context provided, find matching association
        best_score = 0.0
        if disease_context and rows:
            disease_lower = disease_context.lower()
            for r in rows:
                if disease_lower in r["disease"]["name"].lower() or disease_lower in r["disease"].get("id", "").lower():
                    best_score = max(best_score, r["score"])
        elif rows:
            best_score = rows[0]["score"]

        return {
            "overall_score": best_score,
            "n_associations": n_associations,
            "top_diseases": top_diseases,
            "has_data": n_associations > 0,
        }

    except Exception as exc:
        logger.warning("Open Targets fetch failed for %s: %s", uniprot_id, exc)
        return {"overall_score": 0.0, "n_associations": 0, "top_diseases": [], "has_data": False}


def _resolve_ensembl_id(uniprot_id: str) -> Optional[str]:
    """Resolve UniProt accession to Ensembl gene ID via UniProt API."""
    try:
        resp = requests.get(
            f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json",
            timeout=(5, 15),
        )
        if resp.status_code != 200:
            return None
        data = resp.json()
        xrefs = data.get("uniProtKBCrossReferences", [])
        for xref in xrefs:
            if xref.get("database") == "Ensembl":
                props = xref.get("properties", [])
                for prop in props:
                    if prop.get("key") == "GeneId":
                        return prop.get("value")
        # Try gene-centric approach
        gene_name = data.get("genes", [{}])[0].get("geneName", {}).get("value")
        if gene_name:
            # Use Open Targets search to resolve
            search_resp = requests.get(
                f"https://api.platform.opentargets.org/api/v4/graphql",
                json={
                    "query": 'query Search($q: String!) { search(queryString: $q, entityNames: ["target"], page: {size: 1}) { hits { id } } }',
                    "variables": {"q": gene_name},
                },
                timeout=(5, 10),
            )
            if search_resp.status_code == 200:
                hits = search_resp.json().get("data", {}).get("search", {}).get("hits", [])
                if hits:
                    return hits[0]["id"]
        return None
    except Exception as exc:
        logger.debug("Ensembl ID resolution failed for %s: %s", uniprot_id, exc)
        return None


def _get_chembl_target_info(uniprot_id: str) -> dict:
    """Get ChEMBL target statistics using the existing ligands module pattern.

    Returns
    -------
    dict
        Keys: ``chembl_id``, ``n_actives``, ``n_assays``, ``target_type``,
        ``target_class``, ``has_data``.
    """
    from pipeline.ligands import CHEMBL_BASE, HTTP_TIMEOUT as CHEMBL_TIMEOUT

    result = {
        "chembl_id": None,
        "n_actives": 0,
        "n_assays": 0,
        "target_type": "unknown",
        "target_class": "unknown",
        "has_data": False,
    }

    try:
        # Resolve UniProt to ChEMBL target ID
        resp = requests.get(
            f"{CHEMBL_BASE}/target.json",
            params={"target_components__accession": uniprot_id, "limit": 1},
            headers={"Accept": "application/json"},
            timeout=CHEMBL_TIMEOUT,
        )
        if resp.status_code != 200:
            return result

        targets = resp.json().get("targets", [])
        if not targets:
            return result

        target = targets[0]
        chembl_id = target.get("target_chembl_id", "")
        result["chembl_id"] = chembl_id
        result["target_type"] = target.get("target_type", "unknown")

        # Get target class from target_components
        components = target.get("target_components", [])
        if components:
            descs = components[0].get("target_component_synonyms", [])
            for d in descs:
                if d.get("syn_type") == "PROTEIN_CLASS":
                    result["target_class"] = d.get("component_synonym", "unknown")
                    break

        # Count bioactive compounds
        activity_resp = requests.get(
            f"{CHEMBL_BASE}/activity.json",
            params={
                "target_chembl_id": chembl_id,
                "limit": 1,
                "pchembl_value__isnull": "false",
            },
            headers={"Accept": "application/json"},
            timeout=CHEMBL_TIMEOUT,
        )
        if activity_resp.status_code == 200:
            page_meta = activity_resp.json().get("page_meta", {})
            result["n_actives"] = page_meta.get("total_count", 0)
            result["has_data"] = result["n_actives"] > 0

        # Count assays
        assay_resp = requests.get(
            f"{CHEMBL_BASE}/assay.json",
            params={"target_chembl_id": chembl_id, "limit": 1},
            headers={"Accept": "application/json"},
            timeout=CHEMBL_TIMEOUT,
        )
        if assay_resp.status_code == 200:
            result["n_assays"] = assay_resp.json().get("page_meta", {}).get("total_count", 0)

    except Exception as exc:
        logger.warning("ChEMBL target info fetch failed for %s: %s", uniprot_id, exc)

    return result


# ---------------------------------------------------------------------------
# 5 Scoring Modules
# ---------------------------------------------------------------------------

def compute_evidence_score(uniprot_id: str, disease_context: Optional[str] = None) -> dict:
    """Score 1: Evidence — therapeutic validation of the target.

    Based on Open Targets association score + ChEMBL activity count.

    Returns
    -------
    dict
        ``score`` (0-1), ``details``, ``provenance``.
    """
    ot_data = _fetch_open_targets_association(uniprot_id, disease_context)
    chembl_info = _get_chembl_target_info(uniprot_id)

    # Open Targets component (0-0.6)
    ot_score = min(ot_data["overall_score"], 1.0) * 0.6

    # ChEMBL activity component (0-0.4)
    n_actives = chembl_info["n_actives"]
    if n_actives >= 500:
        chembl_component = 0.4
    elif n_actives >= 100:
        chembl_component = 0.3
    elif n_actives >= 10:
        chembl_component = 0.2
    elif n_actives > 0:
        chembl_component = 0.1
    else:
        chembl_component = 0.0

    score = round(min(ot_score + chembl_component, 1.0), 3)

    return {
        "score": score,
        "details": {
            "open_targets_score": ot_data["overall_score"],
            "n_disease_associations": ot_data["n_associations"],
            "top_diseases": ot_data["top_diseases"][:3],
            "chembl_actives": n_actives,
            "chembl_id": chembl_info["chembl_id"],
        },
        "provenance": {
            "open_targets": ot_data["has_data"],
            "chembl": chembl_info["has_data"],
        },
    }


def compute_druggability_score(preview_data: Optional[dict] = None) -> dict:
    """Score 2: Druggability — structural feasibility for small-molecule targeting.

    Based on pocket quality, structure source, and disorder fraction.

    Parameters
    ----------
    preview_data : dict, optional
        Existing target_preview_json data (structure, pockets, disorder).

    Returns
    -------
    dict
        ``score`` (0-1), ``details``, ``provenance``.
    """
    if not preview_data:
        return {"score": 0.3, "details": {"reason": "no_preview_data"}, "provenance": {}}

    score = 0.0
    details = {}

    # Structure quality component (0-0.35)
    structures = preview_data.get("structures", [])
    structure = preview_data.get("structure", {})
    if not structures and structure:
        structures = [structure]

    best_structure = structures[0] if structures else {}
    source = best_structure.get("source", "")
    resolution = best_structure.get("resolution")

    if source in ("pdb_holo", "pdb_experimental"):
        struct_score = 0.35
        if resolution and resolution <= 2.0:
            struct_score = 0.35
        elif resolution and resolution <= 3.0:
            struct_score = 0.30
    elif source == "alphafold":
        confidence = best_structure.get("confidence", 0.7)
        struct_score = 0.20 + (confidence * 0.10)
    elif source == "esmfold":
        struct_score = 0.15
    else:
        struct_score = 0.10

    score += struct_score
    details["structure_source"] = source
    details["structure_score"] = round(struct_score, 3)

    # Pocket quality component (0-0.40)
    pockets = preview_data.get("pockets", [])
    if pockets:
        top_pocket = pockets[0]
        prob = top_pocket.get("probability", 0) or top_pocket.get("druggability", 0)
        method = top_pocket.get("method", "")

        if method == "co-crystallized_ligand":
            pocket_score = 0.40
        elif prob >= 0.8:
            pocket_score = 0.35
        elif prob >= 0.6:
            pocket_score = 0.28
        elif prob >= 0.4:
            pocket_score = 0.18
        else:
            pocket_score = 0.10

        details["pocket_probability"] = prob
        details["pocket_method"] = method
    else:
        pocket_score = 0.05
        details["pocket_status"] = "none_detected"

    score += pocket_score
    details["pocket_score"] = round(pocket_score, 3)

    # Disorder penalty (0-0.25 base, reduced by disorder fraction)
    disorder = preview_data.get("disorder", {})
    fraction_disordered = disorder.get("fraction_disordered", 0.1)
    disorder_penalty = min(fraction_disordered * 0.5, 0.25)
    disorder_component = 0.25 - disorder_penalty
    score += disorder_component
    details["fraction_disordered"] = round(fraction_disordered, 3)
    details["disorder_component"] = round(disorder_component, 3)

    final_score = round(min(score, 1.0), 3)

    return {
        "score": final_score,
        "details": details,
        "provenance": {
            "structure_available": len(structures) > 0,
            "pockets_detected": len(pockets) > 0,
        },
    }


def compute_novelty_score(chembl_info: dict) -> dict:
    """Score 3: Novelty — how unexplored the target is (higher = more novel).

    Inverse of ChEMBL saturation.

    Returns
    -------
    dict
        ``score`` (0-1), ``details``.
    """
    n_actives = chembl_info.get("n_actives", 0)
    n_assays = chembl_info.get("n_assays", 0)

    # More actives = less novel (saturated target)
    if n_actives >= 1000:
        active_component = 0.1  # Heavily explored
    elif n_actives >= 500:
        active_component = 0.25
    elif n_actives >= 100:
        active_component = 0.45
    elif n_actives >= 10:
        active_component = 0.65
    elif n_actives > 0:
        active_component = 0.80
    else:
        active_component = 0.95  # Very novel / unexplored

    # Assay diversity penalty: many assays = well-studied
    if n_assays >= 500:
        assay_factor = 0.8
    elif n_assays >= 100:
        assay_factor = 0.85
    elif n_assays >= 10:
        assay_factor = 0.95
    else:
        assay_factor = 1.0

    score = round(min(active_component * assay_factor, 1.0), 3)

    return {
        "score": score,
        "details": {
            "n_actives": n_actives,
            "n_assays": n_assays,
            "saturation_level": (
                "high" if n_actives >= 500
                else "moderate" if n_actives >= 50
                else "low" if n_actives > 0
                else "unexplored"
            ),
        },
    }


def compute_safety_score(uniprot_id: str, protein_info: Optional[dict] = None) -> dict:
    """Score 4: Safety — risk of off-target liabilities (higher = MORE risk).

    Based on anti-target panel overlap + known dangerous target classification.

    Returns
    -------
    dict
        ``score`` (0-1 risk), ``details``, ``flags``.
    """
    from pipeline.off_target import OFF_TARGET_PANEL, DANGEROUS_TARGETS

    flags = []
    risk_score = 0.0

    # Check if target itself is in the dangerous targets panel
    target_is_dangerous = False
    for name, description in DANGEROUS_TARGETS.items():
        if name.upper() in (uniprot_id or "").upper():
            target_is_dangerous = True
            flags.append(f"Target is a known dangerous target: {name} ({description})")
            risk_score += 0.4
            break

    # Check if target is in the off-target panel (targeting an off-target)
    for panel_name, panel_info in OFF_TARGET_PANEL.items():
        if panel_info.get("uniprot", "").upper() == uniprot_id.upper():
            target_is_dangerous = True
            flags.append(f"Target is a common anti-target: {panel_name} ({panel_info['risk']})")
            risk_score += 0.5
            break

    # Protein class risk heuristics
    if protein_info:
        target_class = protein_info.get("target_class", "").lower()
        target_type = protein_info.get("target_type", "").lower()

        risky_classes = ["ion channel", "gpcr", "nuclear receptor", "transporter"]
        for rc in risky_classes:
            if rc in target_class or rc in target_type:
                risk_score += 0.15
                flags.append(f"Target class '{rc}' has higher off-target risk")
                break

    # Clamp
    risk_score = round(min(risk_score, 1.0), 3)

    return {
        "score": risk_score,
        "details": {
            "is_dangerous_target": target_is_dangerous,
            "n_flags": len(flags),
            "protein_class": protein_info.get("target_class", "unknown") if protein_info else "unknown",
        },
        "flags": flags,
    }


def compute_feasibility_score(chembl_info: dict, target_class: Optional[str] = None) -> dict:
    """Score 5: Feasibility — practical screening feasibility.

    Based on ChEMBL assay availability + target class heuristics.

    Returns
    -------
    dict
        ``score`` (0-1), ``details``.
    """
    n_assays = chembl_info.get("n_assays", 0)
    n_actives = chembl_info.get("n_actives", 0)
    target_type = chembl_info.get("target_type", "unknown")

    # Assay availability (0-0.4)
    if n_assays >= 100:
        assay_component = 0.4
    elif n_assays >= 20:
        assay_component = 0.3
    elif n_assays >= 5:
        assay_component = 0.2
    elif n_assays > 0:
        assay_component = 0.1
    else:
        assay_component = 0.05

    # Active compounds (screening material) (0-0.3)
    if n_actives >= 50:
        active_component = 0.3
    elif n_actives >= 10:
        active_component = 0.2
    elif n_actives > 0:
        active_component = 0.1
    else:
        active_component = 0.0

    # Target class heuristic (0-0.3)
    class_component = 0.15  # default neutral
    if target_class:
        tc = target_class.lower()
        # Well-established druggable classes
        if any(k in tc for k in ["kinase", "protease", "gpcr", "nuclear receptor"]):
            class_component = 0.3
        elif any(k in tc for k in ["enzyme", "transferase", "hydrolase"]):
            class_component = 0.25
        elif any(k in tc for k in ["channel", "transporter"]):
            class_component = 0.20
        elif "protein-protein" in tc:
            class_component = 0.05  # PPI targets are hard
    elif target_type:
        tt = target_type.lower()
        if "single protein" in tt:
            class_component = 0.20
        elif "protein complex" in tt:
            class_component = 0.10

    score = round(min(assay_component + active_component + class_component, 1.0), 3)

    return {
        "score": score,
        "details": {
            "n_assays": n_assays,
            "n_actives": n_actives,
            "target_type": target_type,
            "target_class": target_class or "unknown",
            "assay_component": round(assay_component, 3),
            "active_component": round(active_component, 3),
            "class_component": round(class_component, 3),
        },
    }


# ---------------------------------------------------------------------------
# Aggregation & Recommendation
# ---------------------------------------------------------------------------

def aggregate_scores(
    scores: dict[str, dict],
    weights: Optional[dict[str, float]] = None,
) -> dict:
    """Compute weighted composite score from 5 individual scores.

    Parameters
    ----------
    scores : dict
        Keys: evidence, druggability, novelty, safety, feasibility.
        Each value is a dict with at least ``score`` key (0-1).
    weights : dict, optional
        Custom weights. Defaults to DEFAULT_WEIGHTS.

    Returns
    -------
    dict
        ``composite_score``, ``weighted_scores``, ``weights_used``.
    """
    w = weights or DEFAULT_WEIGHTS.copy()

    # Normalise weights to sum to 1.0
    total_w = sum(w.values())
    if total_w > 0:
        w = {k: v / total_w for k, v in w.items()}

    # Safety is inverted: high safety_risk = bad, so we use (1 - safety_risk)
    safety_risk = scores.get("safety", {}).get("score", 0)
    adjusted_scores = {
        "evidence": scores.get("evidence", {}).get("score", 0),
        "druggability": scores.get("druggability", {}).get("score", 0),
        "novelty": scores.get("novelty", {}).get("score", 0),
        "safety": 1.0 - safety_risk,  # Invert: low risk = high safety score
        "feasibility": scores.get("feasibility", {}).get("score", 0),
    }

    composite = sum(adjusted_scores[k] * w.get(k, 0) for k in adjusted_scores)
    composite = round(min(max(composite, 0), 1.0), 3)

    weighted_scores = {k: round(adjusted_scores[k] * w.get(k, 0), 3) for k in adjusted_scores}

    return {
        "composite_score": composite,
        "adjusted_scores": adjusted_scores,
        "weighted_scores": weighted_scores,
        "weights_used": w,
    }


def generate_recommendation(aggregate: dict, scores: dict) -> dict:
    """Generate GO / CAUTION / NO-GO recommendation.

    Parameters
    ----------
    aggregate : dict
        Output of ``aggregate_scores()``.
    scores : dict
        Individual score dicts with ``flags`` lists.

    Returns
    -------
    dict
        ``recommendation``, ``rationale``, ``flags``, ``critical_flags``.
    """
    composite = aggregate["composite_score"]
    all_flags = []
    critical_flags = []

    # Collect flags from all scores
    for name, score_data in scores.items():
        flags = score_data.get("flags", [])
        all_flags.extend(flags)

    # Check critical conditions
    safety_risk = scores.get("safety", {}).get("score", 0)
    if safety_risk > 0.8:
        critical_flags.append(f"High safety risk ({safety_risk:.2f}): target has significant off-target concerns")

    druggability = scores.get("druggability", {}).get("score", 0)
    if druggability < 0.2:
        critical_flags.append(f"Very low druggability ({druggability:.2f}): poor structural basis for drug design")

    evidence = scores.get("evidence", {}).get("score", 0)
    if evidence < 0.1:
        critical_flags.append(f"Minimal evidence ({evidence:.2f}): target lacks therapeutic validation")

    # Determine recommendation
    has_critical = len(critical_flags) > 0

    if composite >= 0.65 and not has_critical:
        recommendation = "GO"
        rationale = (
            f"Strong composite score ({composite:.2f}/1.00). "
            f"Target shows good evidence, structural druggability, and acceptable safety profile. "
            f"Recommended to proceed with screening."
        )
    elif composite >= 0.40 or (composite >= 0.65 and has_critical):
        recommendation = "CAUTION"
        reasons = []
        if has_critical:
            reasons.append("critical flags detected")
        if composite < 0.65:
            reasons.append(f"moderate composite score ({composite:.2f})")
        if druggability < 0.4:
            reasons.append("limited structural druggability")
        rationale = (
            f"Proceed with caution ({composite:.2f}/1.00). "
            f"Issues: {'; '.join(reasons)}. "
            f"Consider additional validation before large-scale screening."
        )
    else:
        recommendation = "NO-GO"
        reasons = []
        if has_critical:
            reasons.extend(critical_flags)
        if composite < 0.40:
            reasons.append(f"low composite score ({composite:.2f})")
        rationale = (
            f"Not recommended for screening ({composite:.2f}/1.00). "
            f"Blocking issues: {'; '.join(reasons[:3])}. "
            f"Consider alternative targets or fundamental research first."
        )

    # Override: safety_risk > 0.8 → at least CAUTION
    if safety_risk > 0.8 and recommendation == "GO":
        recommendation = "CAUTION"
        rationale = f"Downgraded from GO due to high safety risk ({safety_risk:.2f}). " + rationale

    return {
        "recommendation": recommendation,
        "rationale": rationale,
        "flags": all_flags,
        "critical_flags": critical_flags,
    }


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def assess_target(
    uniprot_id: str,
    disease_context: Optional[str] = None,
    modality: str = "small_molecule",
    preview_data: Optional[dict] = None,
    weights: Optional[dict[str, float]] = None,
) -> dict:
    """Run full target assessment pipeline.

    Parameters
    ----------
    uniprot_id : str
        UniProt accession (e.g. "P00533").
    disease_context : str, optional
        Disease context (EFO ID or name) for evidence scoring.
    modality : str
        Drug modality: "small_molecule" (default), "biologic", "degrader".
    preview_data : dict, optional
        Existing target_preview_json for druggability scoring.
    weights : dict, optional
        Custom aggregation weights.

    Returns
    -------
    dict
        Complete assessment with scores, composite_score, recommendation,
        flags, rationale, provenance, and timing info.
    """
    start = time.time()
    logger.info("Starting target assessment for %s (disease=%s, modality=%s)",
                uniprot_id, disease_context, modality)

    # Compute input hash for idempotence
    input_hash = hashlib.md5(
        json.dumps({
            "uniprot_id": uniprot_id,
            "disease_context": disease_context,
            "modality": modality,
        }, sort_keys=True).encode()
    ).hexdigest()

    # Step 1: Gather ChEMBL target info (shared across multiple scores)
    chembl_info = _get_chembl_target_info(uniprot_id)

    # Step 2: Compute 5 scores
    scores = {}

    scores["evidence"] = compute_evidence_score(uniprot_id, disease_context)
    logger.info("Evidence score: %.3f", scores["evidence"]["score"])

    scores["druggability"] = compute_druggability_score(preview_data)
    logger.info("Druggability score: %.3f", scores["druggability"]["score"])

    scores["novelty"] = compute_novelty_score(chembl_info)
    logger.info("Novelty score: %.3f", scores["novelty"]["score"])

    scores["safety"] = compute_safety_score(uniprot_id, chembl_info)
    logger.info("Safety risk score: %.3f", scores["safety"]["score"])

    scores["feasibility"] = compute_feasibility_score(chembl_info, chembl_info.get("target_class"))
    logger.info("Feasibility score: %.3f", scores["feasibility"]["score"])

    # Step 3: Aggregate
    aggregate = aggregate_scores(scores, weights)
    logger.info("Composite score: %.3f", aggregate["composite_score"])

    # Step 4: Generate recommendation
    recommendation = generate_recommendation(aggregate, scores)
    logger.info("Recommendation: %s", recommendation["recommendation"])

    elapsed = round(time.time() - start, 2)

    return {
        "uniprot_id": uniprot_id,
        "disease_context": disease_context,
        "modality": modality,
        "scores": {
            name: {
                "score": data["score"],
                "details": data.get("details", {}),
            }
            for name, data in scores.items()
        },
        "composite_score": aggregate["composite_score"],
        "adjusted_scores": aggregate["adjusted_scores"],
        "weighted_scores": aggregate["weighted_scores"],
        "weights_used": aggregate["weights_used"],
        "recommendation": recommendation["recommendation"],
        "rationale": recommendation["rationale"],
        "flags": recommendation["flags"],
        "critical_flags": recommendation["critical_flags"],
        "provenance": {
            name: data.get("provenance", {})
            for name, data in scores.items()
        },
        "input_hash": input_hash,
        "elapsed_seconds": elapsed,
    }
