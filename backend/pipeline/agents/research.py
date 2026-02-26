"""
DockIt pipeline -- Research Tools for AI Agents.

Provides literature and database search functions using public APIs
(no API keys required) to enrich agent context with real scientific data.
"""

from __future__ import annotations

import logging
import xml.etree.ElementTree as ET
from typing import Optional

import requests

logger = logging.getLogger(__name__)

TIMEOUT = (5, 15)  # (connect, read) seconds


# ---------------------------------------------------------------------------
# PubMed (NCBI E-utilities)
# ---------------------------------------------------------------------------

def search_pubmed(query: str, max_results: int = 5) -> list[dict]:
    """Search PubMed and return paper metadata."""
    try:
        # Step 1: search for PMIDs
        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        search_resp = requests.get(search_url, params={
            "db": "pubmed",
            "term": query,
            "retmax": max_results,
            "sort": "relevance",
            "retmode": "json",
        }, timeout=TIMEOUT)
        search_resp.raise_for_status()
        id_list = search_resp.json().get("esearchresult", {}).get("idlist", [])
        if not id_list:
            return []

        # Step 2: fetch details
        fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        fetch_resp = requests.get(fetch_url, params={
            "db": "pubmed",
            "id": ",".join(id_list),
            "retmode": "xml",
        }, timeout=TIMEOUT)
        fetch_resp.raise_for_status()

        root = ET.fromstring(fetch_resp.text)
        papers = []
        for article in root.findall(".//PubmedArticle"):
            try:
                medline = article.find("MedlineCitation")
                art = medline.find("Article")
                pmid = medline.findtext("PMID", "")
                title = art.findtext("ArticleTitle", "")

                # Authors
                author_list = art.find("AuthorList")
                authors = []
                if author_list is not None:
                    for auth in author_list.findall("Author")[:3]:
                        last = auth.findtext("LastName", "")
                        initials = auth.findtext("Initials", "")
                        if last:
                            authors.append(f"{last} {initials}".strip())
                    if len(author_list.findall("Author")) > 3:
                        authors.append("et al.")

                # Journal + year
                journal_info = art.find("Journal")
                journal = journal_info.findtext("Title", "") if journal_info is not None else ""
                year_el = journal_info.find(".//Year") if journal_info is not None else None
                year = year_el.text if year_el is not None else ""

                # Abstract
                abstract_el = art.find("Abstract")
                abstract = ""
                if abstract_el is not None:
                    abstract_texts = abstract_el.findall("AbstractText")
                    abstract = " ".join(t.text or "" for t in abstract_texts)[:500]

                papers.append({
                    "pmid": pmid,
                    "title": title,
                    "authors": authors,
                    "journal": journal,
                    "year": year,
                    "abstract": abstract,
                    "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                })
            except Exception:
                continue

        logger.info("PubMed search '%s': %d papers found", query[:50], len(papers))
        return papers

    except Exception as exc:
        logger.warning("PubMed search failed for '%s': %s", query[:50], exc)
        return []


# ---------------------------------------------------------------------------
# ChEMBL REST API
# ---------------------------------------------------------------------------

def search_chembl_target(uniprot_id: str) -> dict:
    """Fetch ChEMBL target profile for a UniProt accession."""
    try:
        base = "https://www.ebi.ac.uk/chembl/api/data"

        # Get target by UniProt ID
        resp = requests.get(
            f"{base}/target.json",
            params={"target_components__accession": uniprot_id, "limit": 1},
            timeout=TIMEOUT,
        )
        resp.raise_for_status()
        results = resp.json().get("targets", [])
        if not results:
            return {}

        target = results[0]
        target_chembl_id = target.get("target_chembl_id", "")

        # Get bioactivity count
        n_bioactivities = 0
        try:
            act_resp = requests.get(
                f"{base}/activity.json",
                params={"target_chembl_id": target_chembl_id, "limit": 0},
                timeout=TIMEOUT,
            )
            act_resp.raise_for_status()
            n_bioactivities = act_resp.json().get("page_meta", {}).get("total_count", 0)
        except Exception:
            pass

        # Get approved drugs
        approved_drugs = []
        try:
            mech_resp = requests.get(
                f"{base}/mechanism.json",
                params={"target_chembl_id": target_chembl_id, "limit": 20},
                timeout=TIMEOUT,
            )
            mech_resp.raise_for_status()
            mechs = mech_resp.json().get("mechanisms", [])
            seen = set()
            for m in mechs:
                drug_name = m.get("molecule_name") or m.get("molecule_chembl_id", "")
                if drug_name and drug_name not in seen:
                    seen.add(drug_name)
                    approved_drugs.append({
                        "name": drug_name,
                        "mechanism": m.get("mechanism_of_action", ""),
                        "action_type": m.get("action_type", ""),
                    })
        except Exception:
            pass

        result = {
            "target_chembl_id": target_chembl_id,
            "target_type": target.get("target_type", ""),
            "organism": target.get("organism", ""),
            "pref_name": target.get("pref_name", ""),
            "n_bioactivities": n_bioactivities,
            "approved_drugs": approved_drugs,
            "n_approved_drugs": len(approved_drugs),
        }
        logger.info("ChEMBL target for %s: %s (%d bioactivities, %d drugs)",
                     uniprot_id, target_chembl_id, n_bioactivities, len(approved_drugs))
        return result

    except Exception as exc:
        logger.warning("ChEMBL target search failed for %s: %s", uniprot_id, exc)
        return {}


# ---------------------------------------------------------------------------
# UniProt REST API
# ---------------------------------------------------------------------------

def fetch_uniprot_info(uniprot_id: str) -> dict:
    """Fetch comprehensive protein information from UniProt."""
    try:
        resp = requests.get(
            f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json",
            timeout=TIMEOUT,
        )
        resp.raise_for_status()
        data = resp.json()

        protein_name = (
            data.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
            or (data.get("proteinDescription", {}).get("submittedName", [{}])[0].get("fullName", {}).get("value", ""))
        )
        gene_name = ""
        genes = data.get("genes", [])
        if genes:
            gene_name = genes[0].get("geneName", {}).get("value", "")

        organism = data.get("organism", {}).get("scientificName", "")
        seq_length = data.get("sequence", {}).get("length", 0)

        # Function
        function_text = ""
        for comment in data.get("comments", []):
            if comment.get("commentType") == "FUNCTION":
                texts = comment.get("texts", [])
                if texts:
                    function_text = texts[0].get("value", "")[:500]
                break

        # Subcellular locations
        locations = []
        for comment in data.get("comments", []):
            if comment.get("commentType") == "SUBCELLULAR LOCATION":
                for loc in comment.get("subcellularLocations", []):
                    val = loc.get("location", {}).get("value", "")
                    if val:
                        locations.append(val)

        # Diseases
        diseases = []
        for comment in data.get("comments", []):
            if comment.get("commentType") == "DISEASE":
                disease = comment.get("disease", {})
                name = disease.get("diseaseId", "") or disease.get("description", "")
                if name:
                    diseases.append({
                        "name": name,
                        "acronym": disease.get("acronym", ""),
                    })

        # Protein families
        families = []
        for comment in data.get("comments", []):
            if comment.get("commentType") == "SIMILARITY":
                texts = comment.get("texts", [])
                if texts:
                    families.append(texts[0].get("value", ""))

        result = {
            "protein_name": protein_name,
            "gene_name": gene_name,
            "organism": organism,
            "function": function_text,
            "subcellular_locations": locations,
            "diseases": diseases,
            "families": families,
            "sequence_length": seq_length,
        }
        logger.info("UniProt info for %s: %s (%s)", uniprot_id, gene_name, protein_name[:40])
        return result

    except Exception as exc:
        logger.warning("UniProt fetch failed for %s: %s", uniprot_id, exc)
        return {}


# ---------------------------------------------------------------------------
# ClinicalTrials.gov v2 API
# ---------------------------------------------------------------------------

def search_clinical_trials(gene_name: str, max_results: int = 5) -> list[dict]:
    """Search ClinicalTrials.gov for trials related to a gene/target."""
    try:
        resp = requests.get(
            "https://clinicaltrials.gov/api/v2/studies",
            params={
                "query.term": gene_name,
                "pageSize": max_results,
                "sort": "LastUpdatePostDate:desc",
                "fields": "NCTId,BriefTitle,OverallStatus,Phase,Condition",
                "format": "json",
            },
            timeout=TIMEOUT,
        )
        resp.raise_for_status()
        data = resp.json()

        trials = []
        for study in data.get("studies", []):
            proto = study.get("protocolSection", {})
            id_mod = proto.get("identificationModule", {})
            status_mod = proto.get("statusModule", {})
            design_mod = proto.get("designModule", {})
            cond_mod = proto.get("conditionsModule", {})

            nct_id = id_mod.get("nctId", "")
            title = id_mod.get("briefTitle", "")
            status = status_mod.get("overallStatus", "")
            phases = design_mod.get("phases", [])
            phase = phases[0] if phases else "N/A"
            conditions = cond_mod.get("conditions", [])

            trials.append({
                "nct_id": nct_id,
                "title": title,
                "phase": phase,
                "status": status,
                "conditions": conditions[:3],
                "url": f"https://clinicaltrials.gov/study/{nct_id}",
            })

        logger.info("ClinicalTrials search '%s': %d trials found", gene_name, len(trials))
        return trials

    except Exception as exc:
        logger.warning("ClinicalTrials search failed for '%s': %s", gene_name, exc)
        return []


# ---------------------------------------------------------------------------
# DGIdb API
# ---------------------------------------------------------------------------

def fetch_drug_interactions(gene_name: str) -> list[dict]:
    """Fetch drug-gene interactions from DGIdb."""
    try:
        resp = requests.get(
            "https://dgidb.org/api/v2/interactions.json",
            params={"genes": gene_name},
            timeout=TIMEOUT,
        )
        resp.raise_for_status()
        data = resp.json()

        interactions = []
        for matched_term in data.get("matchedTerms", []):
            for interaction in matched_term.get("interactions", [])[:10]:
                drug_name = interaction.get("drugName", "")
                int_types = interaction.get("interactionTypes", [])
                sources = interaction.get("sources", [])
                if drug_name:
                    interactions.append({
                        "drug_name": drug_name,
                        "interaction_types": int_types,
                        "sources": sources[:3],
                        "score": interaction.get("score", None),
                    })

        logger.info("DGIdb interactions for '%s': %d found", gene_name, len(interactions))
        return interactions

    except Exception as exc:
        logger.warning("DGIdb fetch failed for '%s': %s", gene_name, exc)
        return []


# ---------------------------------------------------------------------------
# Research Package Builder
# ---------------------------------------------------------------------------

def build_research_package(
    agent_type: str,
    uniprot_id: Optional[str] = None,
    gene_name: Optional[str] = None,
    target_name: Optional[str] = None,
) -> dict:
    """Assemble a research package tailored to each agent type.

    Parameters
    ----------
    agent_type : str
        One of: "target", "run_analysis", "candidate", "optimization"
    uniprot_id : str, optional
        UniProt accession ID (e.g. P00533)
    gene_name : str, optional
        Gene symbol (e.g. EGFR)
    target_name : str, optional
        Protein name (e.g. Epidermal growth factor receptor)
    """
    package = {"agent_type": agent_type, "research_available": False}
    search_term = gene_name or target_name or uniprot_id or ""

    if not search_term:
        return package

    if agent_type == "target":
        # Full research: UniProt + ChEMBL + PubMed + trials + interactions
        if uniprot_id:
            package["uniprot_info"] = fetch_uniprot_info(uniprot_id)
            package["chembl_profile"] = search_chembl_target(uniprot_id)
            # Use gene name from UniProt if not provided
            if not gene_name and package["uniprot_info"].get("gene_name"):
                gene_name = package["uniprot_info"]["gene_name"]
                search_term = gene_name

        package["pubmed_papers"] = search_pubmed(
            f"{search_term} target validation drug discovery", max_results=5
        )
        if gene_name:
            package["clinical_trials"] = search_clinical_trials(gene_name, max_results=5)
            package["drug_interactions"] = fetch_drug_interactions(gene_name)

    elif agent_type == "run_analysis":
        if uniprot_id:
            package["chembl_profile"] = search_chembl_target(uniprot_id)
        package["pubmed_papers"] = search_pubmed(
            f"{search_term} virtual screening hit rate benchmark", max_results=5
        )

    elif agent_type == "candidate":
        package["pubmed_papers"] = search_pubmed(
            f"{search_term} SAR structure-activity inhibitor", max_results=5
        )
        if gene_name:
            package["clinical_trials"] = search_clinical_trials(gene_name, max_results=3)
        if uniprot_id:
            package["chembl_profile"] = search_chembl_target(uniprot_id)

    elif agent_type == "optimization":
        package["pubmed_papers"] = search_pubmed(
            f"{search_term} lead optimization medicinal chemistry SAR", max_results=5
        )
        if uniprot_id:
            package["chembl_profile"] = search_chembl_target(uniprot_id)

    # Summary stats
    package["research_available"] = True
    package["research_summary"] = {
        "pubmed_papers_found": len(package.get("pubmed_papers", [])),
        "clinical_trials_found": len(package.get("clinical_trials", [])),
        "has_chembl_profile": bool(package.get("chembl_profile")),
        "has_uniprot_info": bool(package.get("uniprot_info")),
        "drug_interactions_found": len(package.get("drug_interactions", [])),
    }

    return package
