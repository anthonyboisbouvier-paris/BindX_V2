#!/usr/bin/env python3
"""
PharmacoDB -- Phase 3 Ingestion Pipeline
=========================================
Adds 10 additional public data sources to PharmacoDB:

 1. DGIdb      -- Drug-Gene Interaction Database
 2. SIDER      -- Side Effect Resource (MedDRA side effects)
 3. OpenFDA    -- FDA Adverse Event Reporting System
 4. ClinicalTrials.gov -- Clinical trial metadata
 5. UniChem    -- Cross-database compound identifiers
 6. ChEBI      -- Chemical Entities of Biological Interest (ontology)
 7. STRING     -- Protein-Protein Interactions
 8. Reactome   -- Biological Pathway Data
 9. DisGeNET   -- Disease-Gene Associations
10. KEGG       -- Pathway & Drug information

Requires: psycopg2, schema_phase3.sql already applied.
DB: localhost:5433, pharmaco, postgres/pharmaco_secret
"""

import gzip
import io
import json
import math
import os
import re
import sys
import time
import traceback
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

import psycopg2
import psycopg2.extras

# ============================================================
# CONFIG
# ============================================================

DB_CONFIG: Dict[str, Any] = {
    "host": "localhost",
    "port": 5433,
    "dbname": "pharmaco",
    "user": "postgres",
    "password": "pharmaco_secret",
}

# Minimum delay between API calls (seconds)
MIN_DELAY: float = 0.35

# Global retry config
MAX_RETRIES: int = 3
REQUEST_TIMEOUT: int = 30

# Batch insert page size
BATCH_PAGE_SIZE: int = 500

# ============================================================
# LOGGING
# ============================================================

def log(msg: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


# ============================================================
# DATABASE HELPERS
# ============================================================

def get_conn() -> psycopg2.extensions.connection:
    return psycopg2.connect(**DB_CONFIG)


def start_log(conn: psycopg2.extensions.connection, source: str, step: str) -> int:
    with conn.cursor() as cur:
        cur.execute("""
            INSERT INTO pharmaco_db.ingestion_log (source, step, status)
            VALUES (%s, %s, 'running') RETURNING id
        """, (source, step))
        log_id: int = cur.fetchone()[0]
    conn.commit()
    return log_id


def update_log(
    conn: psycopg2.extensions.connection,
    log_id: int,
    status: str,
    rows: Optional[int] = None,
    error: Optional[str] = None,
) -> None:
    with conn.cursor() as cur:
        cur.execute("""
            UPDATE pharmaco_db.ingestion_log
            SET status=%s, rows_inserted=%s, error_message=%s, completed_at=NOW()
            WHERE id=%s
        """, (status, rows, error, log_id))
    conn.commit()


# ============================================================
# HTTP HELPERS
# ============================================================

def api_get(
    url: str,
    retries: int = MAX_RETRIES,
    timeout: int = REQUEST_TIMEOUT,
    headers: Optional[Dict[str, str]] = None,
) -> Optional[Any]:
    """GET request returning parsed JSON, or None on failure."""
    hdrs: Dict[str, str] = {
        "Accept": "application/json",
        "User-Agent": "PharmacoDB-Ingestion/3.0 (academic research)",
    }
    if headers:
        hdrs.update(headers)
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=hdrs)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8", errors="replace"))
        except urllib.error.HTTPError as e:
            if e.code == 429:
                # Rate limited -- back off
                wait = 10 * (attempt + 1)
                log(f"    429 rate-limited, waiting {wait}s...")
                time.sleep(wait)
            elif e.code == 404:
                return None
            elif attempt < retries - 1:
                time.sleep(3 * (attempt + 1))
            else:
                return None
        except Exception:
            if attempt < retries - 1:
                time.sleep(3 * (attempt + 1))
            else:
                return None
    return None


def api_post(
    url: str,
    data: Dict[str, Any],
    retries: int = MAX_RETRIES,
    timeout: int = REQUEST_TIMEOUT,
) -> Optional[Any]:
    """POST JSON request returning parsed JSON response."""
    encoded = json.dumps(data).encode("utf-8")
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url,
                data=encoded,
                headers={
                    "Content-Type": "application/json",
                    "Accept": "application/json",
                    "User-Agent": "PharmacoDB-Ingestion/3.0",
                },
            )
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8", errors="replace"))
        except urllib.error.HTTPError as e:
            if e.code == 429:
                time.sleep(10 * (attempt + 1))
            elif e.code == 404:
                return None
            elif attempt < retries - 1:
                time.sleep(3 * (attempt + 1))
            else:
                return None
        except Exception:
            if attempt < retries - 1:
                time.sleep(3 * (attempt + 1))
            else:
                return None
    return None


def fetch_text(
    url: str,
    retries: int = MAX_RETRIES,
    timeout: int = 60,
) -> Optional[str]:
    """Fetch raw text content from a URL."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={
                "User-Agent": "PharmacoDB-Ingestion/3.0",
            })
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except Exception:
            if attempt < retries - 1:
                time.sleep(5 * (attempt + 1))
            else:
                return None
    return None


def fetch_gzip(url: str, retries: int = MAX_RETRIES, timeout: int = 120) -> Optional[bytes]:
    """Fetch gzipped content and return decompressed bytes."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={
                "User-Agent": "PharmacoDB-Ingestion/3.0",
            })
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                raw = resp.read()
                return gzip.decompress(raw)
        except Exception:
            if attempt < retries - 1:
                time.sleep(5 * (attempt + 1))
            else:
                return None
    return None


# ============================================================
# LOOKUP MAP BUILDERS (built once, shared across functions)
# ============================================================

def build_target_gene_map(conn: psycopg2.extensions.connection) -> Dict[str, int]:
    """gene_name (uppercase) -> target_id"""
    m: Dict[str, int] = {}
    with conn.cursor() as cur:
        cur.execute("""
            SELECT id, gene_name FROM pharmaco_db.targets
            WHERE gene_name IS NOT NULL
        """)
        for tid, gene in cur.fetchall():
            m[gene.upper()] = tid
    return m


def build_target_uniprot_map(conn: psycopg2.extensions.connection) -> Dict[str, int]:
    """uniprot_id -> target_id"""
    m: Dict[str, int] = {}
    with conn.cursor() as cur:
        cur.execute("""
            SELECT id, uniprot_id FROM pharmaco_db.targets
            WHERE uniprot_id IS NOT NULL
        """)
        for tid, uid in cur.fetchall():
            m[uid] = tid
    return m


def build_target_name_map(conn: psycopg2.extensions.connection) -> Dict[str, int]:
    """pref_name (uppercase) -> target_id"""
    m: Dict[str, int] = {}
    with conn.cursor() as cur:
        cur.execute("SELECT id, pref_name FROM pharmaco_db.targets WHERE pref_name IS NOT NULL")
        for tid, name in cur.fetchall():
            m[name.upper()] = tid
    return m


def build_compound_name_map(conn: psycopg2.extensions.connection) -> Dict[str, int]:
    """pref_name (uppercase) -> compound_id"""
    m: Dict[str, int] = {}
    with conn.cursor() as cur:
        cur.execute("""
            SELECT id, pref_name FROM pharmaco_db.compounds
            WHERE pref_name IS NOT NULL
        """)
        for cid, name in cur.fetchall():
            m[name.upper()] = cid
    return m


def build_compound_chembl_map(conn: psycopg2.extensions.connection) -> Dict[str, int]:
    """chembl_id -> compound_id"""
    m: Dict[str, int] = {}
    with conn.cursor() as cur:
        cur.execute("SELECT id, chembl_id FROM pharmaco_db.compounds WHERE chembl_id IS NOT NULL")
        for cid, chembl in cur.fetchall():
            m[chembl] = cid
    return m


def build_compound_pubchem_map(conn: psycopg2.extensions.connection) -> Dict[int, int]:
    """pubchem_cid -> compound_id"""
    m: Dict[int, int] = {}
    with conn.cursor() as cur:
        cur.execute("SELECT id, pubchem_cid FROM pharmaco_db.compounds WHERE pubchem_cid IS NOT NULL")
        for cid, pcid in cur.fetchall():
            m[pcid] = cid
    return m


def build_approved_drugs(conn: psycopg2.extensions.connection) -> List[Tuple[int, str]]:
    """Return (compound_id, pref_name) for approved drugs (max_phase >= 4)."""
    with conn.cursor() as cur:
        cur.execute("""
            SELECT id, pref_name FROM pharmaco_db.compounds
            WHERE max_phase >= 4 AND pref_name IS NOT NULL
            ORDER BY id
        """)
        return cur.fetchall()


def build_targets_with_uniprot(conn: psycopg2.extensions.connection) -> List[Tuple[int, str, str]]:
    """Return (target_id, uniprot_id, gene_name) for human targets."""
    with conn.cursor() as cur:
        cur.execute("""
            SELECT id, uniprot_id, gene_name FROM pharmaco_db.targets
            WHERE uniprot_id IS NOT NULL
            AND organism = 'Homo sapiens'
            ORDER BY id
        """)
        return cur.fetchall()


def build_targets_with_gene(conn: psycopg2.extensions.connection) -> List[Tuple[int, str]]:
    """Return (target_id, gene_name) for human targets with gene_name."""
    with conn.cursor() as cur:
        cur.execute("""
            SELECT id, gene_name FROM pharmaco_db.targets
            WHERE gene_name IS NOT NULL
            AND organism = 'Homo sapiens'
            ORDER BY id
        """)
        return cur.fetchall()


# ============================================================
# 1. DGIdb -- Drug-Gene Interaction Database
# ============================================================

def ingest_dgidb(conn: psycopg2.extensions.connection) -> int:
    """
    Ingest drug-gene interactions from DGIdb v5 GraphQL API.
    Strategy: query per-gene for our targets' gene_names (batches of 25).
    API: POST https://dgidb.org/api/graphql
    """
    log("=" * 60)
    log("1. DGIdb -- Drug-Gene Interactions (GraphQL v5)")
    log("=" * 60)
    log_id = start_log(conn, "dgidb", "interactions")

    target_gene_map = build_target_gene_map(conn)
    compound_name_map = build_compound_name_map(conn)
    compound_chembl_map = build_compound_chembl_map(conn)

    gene_list = list(target_gene_map.keys())
    log(f"  Querying DGIdb for {len(gene_list)} genes")

    total = 0
    errors = 0

    # DGIdb GraphQL supports batch gene queries
    batch_size = 25
    for batch_start in range(0, len(gene_list), batch_size):
        batch_genes = gene_list[batch_start : batch_start + batch_size]
        genes_json = json.dumps(batch_genes)

        if batch_start % 500 == 0:
            log(f"    Progress: {batch_start}/{len(gene_list)} genes ({total} interactions)")

        # DGIdb v5 GraphQL query
        graphql_query = {
            "query": f"""{{
                genes(names: {genes_json}) {{
                    nodes {{
                        name
                        interactions {{
                            drug {{ name conceptId }}
                            interactionScore
                            interactionTypes {{ type directionality }}
                            sources {{ fullName }}
                            publications {{ pmid }}
                        }}
                    }}
                }}
            }}"""
        }

        data = api_post("https://dgidb.org/api/graphql", graphql_query, timeout=60)
        if not data or "errors" in data:
            errors += 1
            time.sleep(MIN_DELAY)
            continue

        gene_nodes = data.get("data", {}).get("genes", {}).get("nodes", [])
        insert_batch: List[Tuple] = []

        for gene_node in gene_nodes:
            gene_name = (gene_node.get("name") or "").upper()
            target_id = target_gene_map.get(gene_name)

            for interaction in gene_node.get("interactions", []):
                drug_info = interaction.get("drug", {})
                drug_name = drug_info.get("name")
                drug_concept = drug_info.get("conceptId", "")

                # Extract ChEMBL ID from conceptId if available
                drug_chembl = None
                if drug_concept and "chembl" in drug_concept.lower():
                    drug_chembl = drug_concept

                int_types = interaction.get("interactionTypes", [])
                int_type = int_types[0].get("type") if int_types else None

                sources_list = [s.get("fullName") for s in interaction.get("sources", []) if s.get("fullName")]
                pmids_list = [str(p.get("pmid")) for p in interaction.get("publications", []) if p.get("pmid")]

                score = interaction.get("interactionScore")
                if score is not None:
                    try:
                        score = float(score)
                    except (ValueError, TypeError):
                        score = None

                # Try to match compound
                compound_id = None
                if drug_chembl:
                    compound_id = compound_chembl_map.get(drug_chembl)
                if compound_id is None and drug_name:
                    compound_id = compound_name_map.get(drug_name.upper())

                insert_batch.append((
                    target_id,
                    gene_name,
                    drug_name,
                    drug_chembl,
                    compound_id,
                    int_type,
                    score,
                    sources_list if sources_list else None,
                    pmids_list if pmids_list else None,
                ))

        if insert_batch:
            with conn.cursor() as cur:
                psycopg2.extras.execute_values(cur, """
                    INSERT INTO pharmaco_db.dgidb_interactions
                    (target_id, gene_name, drug_name, drug_chembl_id, compound_id,
                     interaction_type, interaction_score, sources, pmids)
                    VALUES %s
                    ON CONFLICT (gene_name, drug_name, interaction_type)
                    WHERE drug_name IS NOT NULL
                    DO NOTHING
                """, insert_batch, page_size=BATCH_PAGE_SIZE)
            conn.commit()
            total += len(insert_batch)

        time.sleep(MIN_DELAY)

    log(f"  DONE: {total} DGIdb interactions inserted ({errors} batch errors)")
    update_log(conn, log_id, "completed", total)
    return total


# ============================================================
# 2. SIDER -- Side Effect Resource
# ============================================================

def ingest_sider(conn: psycopg2.extensions.connection) -> int:
    """
    Ingest side effects from SIDER (MedDRA).
    Download: http://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz
    STITCH IDs map to PubChem CIDs: CIDm + CIDs prefixes, strip prefix and
    leading zeros to get PubChem CID.
    """
    log("=" * 60)
    log("2. SIDER -- Side Effects (MedDRA)")
    log("=" * 60)
    log_id = start_log(conn, "sider", "side_effects")

    compound_pubchem_map = build_compound_pubchem_map(conn)
    log(f"  Compound PubChem map: {len(compound_pubchem_map)} entries")

    # Try to download the SIDER TSV
    url = "http://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz"
    log(f"  Downloading SIDER data from {url}")
    raw_bytes = fetch_gzip(url, timeout=180)

    if not raw_bytes:
        # Fallback: try alternative URL
        log("  Primary SIDER URL failed, trying alternative...")
        url2 = "http://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz"
        raw_bytes = fetch_gzip(url2, timeout=180)

    if not raw_bytes:
        log("  ERROR: Could not download SIDER data")
        update_log(conn, log_id, "failed", 0, "Download failed")
        return 0

    text = raw_bytes.decode("utf-8", errors="replace")
    lines = text.strip().split("\n")
    log(f"  Downloaded {len(lines)} side effect records")

    total = 0
    batch: List[Tuple] = []

    for line in lines:
        parts = line.split("\t")
        if len(parts) < 6:
            continue

        stitch_id = parts[0]  # e.g., CID100000085
        # SIDER STITCH IDs: CIDsXXXXXXXXX (stereo) or CIDmXXXXXXXXX (merged)
        # Remove prefix and leading zeros to get PubChem CID
        cid_str = re.sub(r"^CID[sm]0*", "", stitch_id)
        pubchem_cid = None
        try:
            pubchem_cid = int(cid_str)
        except (ValueError, TypeError):
            pass

        compound_id = compound_pubchem_map.get(pubchem_cid) if pubchem_cid else None

        meddra_concept_type = parts[2] if len(parts) > 2 else None  # PT or LLT
        meddra_id = parts[4] if len(parts) > 4 else None
        side_effect_name = parts[5] if len(parts) > 5 else parts[3] if len(parts) > 3 else None

        if not side_effect_name:
            continue

        batch.append((
            compound_id,
            stitch_id,
            side_effect_name,
            meddra_concept_type,
            meddra_id,
            None,  # frequency not in this file
            "sider",
        ))

        if len(batch) >= 5000:
            with conn.cursor() as cur:
                psycopg2.extras.execute_values(cur, """
                    INSERT INTO pharmaco_db.side_effects
                    (compound_id, stitch_id, side_effect_name, meddra_concept_type,
                     meddra_id, frequency, source)
                    VALUES %s
                    ON CONFLICT (stitch_id, meddra_id, meddra_concept_type)
                    WHERE stitch_id IS NOT NULL AND meddra_id IS NOT NULL
                    DO NOTHING
                """, batch, page_size=BATCH_PAGE_SIZE)
            conn.commit()
            total += len(batch)
            batch = []

            if total % 50000 == 0:
                log(f"    Inserted {total:,} side effects so far...")

    # Flush remaining
    if batch:
        with conn.cursor() as cur:
            psycopg2.extras.execute_values(cur, """
                INSERT INTO pharmaco_db.side_effects
                (compound_id, stitch_id, side_effect_name, meddra_concept_type,
                 meddra_id, frequency, source)
                VALUES %s
                ON CONFLICT (stitch_id, meddra_id, meddra_concept_type)
                WHERE stitch_id IS NOT NULL AND meddra_id IS NOT NULL
                DO NOTHING
            """, batch, page_size=BATCH_PAGE_SIZE)
        conn.commit()
        total += len(batch)

    log(f"  DONE: {total:,} SIDER side effects inserted")
    update_log(conn, log_id, "completed", total)
    return total


# ============================================================
# 3. OpenFDA -- Adverse Events
# ============================================================

def ingest_openfda(conn: psycopg2.extensions.connection) -> int:
    """
    Ingest adverse event data from OpenFDA for approved drugs.
    API: https://api.fda.gov/drug/event.json?search=...&count=...
    Rate limit: 240 req/min (without key).
    """
    log("=" * 60)
    log("3. OpenFDA -- Adverse Events")
    log("=" * 60)
    log_id = start_log(conn, "openfda", "adverse_events")

    approved_drugs = build_approved_drugs(conn)
    log(f"  Found {len(approved_drugs)} approved drugs to query")

    total = 0
    errors = 0
    queried = 0

    for compound_id, drug_name in approved_drugs:
        queried += 1
        if queried % 50 == 0:
            log(f"    Progress: {queried}/{len(approved_drugs)} drugs ({total} events)")

        # URL-encode the drug name
        encoded_name = urllib.parse.quote(drug_name)
        url = (
            f"https://api.fda.gov/drug/event.json"
            f"?search=patient.drug.openfda.generic_name:\"{encoded_name}\""
            f"&count=patient.reaction.reactionmeddrapt.exact"
            f"&limit=50"
        )

        data = api_get(url, timeout=15)
        if not data or "results" not in data:
            errors += 1
            time.sleep(MIN_DELAY)
            continue

        batch: List[Tuple] = []
        for result in data["results"][:50]:  # Top 50 reactions
            reaction = result.get("term")
            count = result.get("count")
            if reaction and count:
                batch.append((
                    compound_id,
                    drug_name,
                    reaction,
                    count,
                    "openfda",
                ))

        if batch:
            with conn.cursor() as cur:
                psycopg2.extras.execute_values(cur, """
                    INSERT INTO pharmaco_db.adverse_events
                    (compound_id, drug_name, reaction, count, source)
                    VALUES %s
                    ON CONFLICT (drug_name, reaction)
                    WHERE drug_name IS NOT NULL
                    DO NOTHING
                """, batch, page_size=BATCH_PAGE_SIZE)
            conn.commit()
            total += len(batch)

        # OpenFDA rate limit: 240/min = 4/sec, so 0.3s is safe
        time.sleep(MIN_DELAY)

    log(f"  DONE: {total} adverse events inserted ({errors} drugs not found in FDA)")
    update_log(conn, log_id, "completed", total)
    return total


# ============================================================
# 4. ClinicalTrials.gov -- Clinical Trials
# ============================================================

def ingest_clinical_trials(conn: psycopg2.extensions.connection) -> int:
    """
    Ingest clinical trial metadata from ClinicalTrials.gov v2 API.
    API: https://clinicaltrials.gov/api/v2/studies?query.term=DRUG&pageSize=20
    """
    log("=" * 60)
    log("4. ClinicalTrials.gov -- Clinical Trials")
    log("=" * 60)
    log_id = start_log(conn, "clinicaltrials", "trials")

    approved_drugs = build_approved_drugs(conn)
    log(f"  Found {len(approved_drugs)} approved drugs to query")

    total = 0
    errors = 0
    queried = 0

    for compound_id, drug_name in approved_drugs:
        queried += 1
        if queried % 50 == 0:
            log(f"    Progress: {queried}/{len(approved_drugs)} drugs ({total} trials)")

        encoded = urllib.parse.quote(drug_name)
        url = (
            f"https://clinicaltrials.gov/api/v2/studies"
            f"?query.term={encoded}"
            f"&pageSize=20"
            f"&fields=NCTId,BriefTitle,OverallStatus,Phase,Condition,InterventionName,"
            f"StartDate,PrimaryCompletionDate,EnrollmentCount"
        )

        data = api_get(url, timeout=20)
        if not data:
            errors += 1
            time.sleep(MIN_DELAY)
            continue

        studies = data.get("studies", [])
        batch: List[Tuple] = []

        for study in studies:
            proto = study.get("protocolSection", {})
            ident = proto.get("identificationModule", {})
            status_mod = proto.get("statusModule", {})
            design = proto.get("designModule", {})
            conditions_mod = proto.get("conditionsModule", {})
            arms_mod = proto.get("armsInterventionsModule", {})

            nct_id = ident.get("nctId")
            if not nct_id:
                continue

            title = ident.get("briefTitle")
            status = status_mod.get("overallStatus")

            # Phase can be in designModule.phases (list)
            phases = design.get("phases", [])
            phase = phases[0] if phases else None

            # Conditions
            conditions = conditions_mod.get("conditions", [])

            # Interventions
            interventions_raw = arms_mod.get("interventions", [])
            intervention_names = [i.get("name", "") for i in interventions_raw]

            start_date = status_mod.get("startDateStruct", {}).get("date")
            completion_date = status_mod.get("primaryCompletionDateStruct", {}).get("date")

            enrollment_info = design.get("enrollmentInfo", {})
            enrollment = enrollment_info.get("count") if enrollment_info else None

            batch.append((
                compound_id,
                nct_id,
                title,
                status,
                phase,
                conditions if conditions else None,
                intervention_names if intervention_names else None,
                start_date,
                completion_date,
                enrollment,
                "clinicaltrials.gov",
            ))

        if batch:
            with conn.cursor() as cur:
                psycopg2.extras.execute_values(cur, """
                    INSERT INTO pharmaco_db.clinical_trials
                    (compound_id, nct_id, title, status, phase, conditions,
                     interventions, start_date, completion_date, enrollment, source)
                    VALUES %s
                    ON CONFLICT (nct_id) DO NOTHING
                """, batch, page_size=BATCH_PAGE_SIZE)
            conn.commit()
            total += len(batch)

        time.sleep(MIN_DELAY)

    log(f"  DONE: {total} clinical trials inserted ({errors} query failures)")
    update_log(conn, log_id, "completed", total)
    return total


# ============================================================
# 5. UniChem -- Cross-Database Compound Identifiers
# ============================================================

def ingest_unichem(conn: psycopg2.extensions.connection) -> int:
    """
    Ingest cross-references from UniChem for our compounds.
    API: POST https://www.ebi.ac.uk/unichem/api/v1/compounds
         with body {"compound": "<InChIKey>", "type": "inchikey"}
    Returns mappings to PubChem, DrugBank, ZINC, BindingDB, PDB, ChEBI, etc.
    """
    log("=" * 60)
    log("5. UniChem -- Cross-Database Identifiers")
    log("=" * 60)
    log_id = start_log(conn, "unichem", "cross_references")

    # Get compounds with InChI keys (required for UniChem v1 API)
    compounds_with_inchikey: List[Tuple[int, str, str]] = []
    with conn.cursor() as cur:
        cur.execute("""
            SELECT DISTINCT c.id, c.inchi_key, c.chembl_id
            FROM pharmaco_db.compounds c
            WHERE c.inchi_key IS NOT NULL
            AND (c.max_phase >= 1 OR c.id IN (
                SELECT DISTINCT compound_id FROM pharmaco_db.bioactivities LIMIT 50000
            ))
            ORDER BY c.id
            LIMIT 10000
        """)
        compounds_with_inchikey = cur.fetchall()

    log(f"  Querying UniChem for {len(compounds_with_inchikey)} compounds (InChIKey-based)")

    # Known source name mapping (from UniChem /api/v1/sources)
    # We identify sources by their baseIdUrl or name
    known_dbs: Dict[str, Tuple[str, str]] = {
        # baseIdUrl pattern -> (db_name, url_template)
        "chembldb": ("ChEMBL", "https://www.ebi.ac.uk/chembl/compound_report_card/{}"),
        "drugbank": ("DrugBank", "https://go.drugbank.com/drugs/{}"),
        "rcsb": ("PDB", "https://www.rcsb.org/ligand/{}"),
        "chebi": ("ChEBI", "https://www.ebi.ac.uk/chebi/searchId.do?chebiId={}"),
        "zinc": ("ZINC", "https://zinc15.docking.org/substances/{}"),
        "pubchem": ("PubChem", "https://pubchem.ncbi.nlm.nih.gov/compound/{}"),
        "bindingdb": ("BindingDB", "https://www.bindingdb.org/bind/chemsearch/marvin/MolStructure.jsp?monomerid={}"),
        "hmdb": ("HMDB", "https://hmdb.ca/metabolites/{}"),
        "pharmgkb": ("PharmGKB", "https://www.pharmgkb.org/chemical/{}"),
        "selleckchem": ("SelleckChem", "https://www.selleckchem.com/products/{}.html"),
    }

    total = 0
    errors = 0

    for i, (compound_id, inchi_key, chembl_id) in enumerate(compounds_with_inchikey):
        if i % 200 == 0 and i > 0:
            log(f"    Progress: {i}/{len(compounds_with_inchikey)} ({total} xrefs)")

        payload = {"compound": inchi_key, "type": "inchikey"}
        data = api_post(
            "https://www.ebi.ac.uk/unichem/api/v1/compounds",
            payload,
            timeout=15,
        )
        if not data or data.get("response") == "Not found":
            errors += 1
            time.sleep(MIN_DELAY)
            continue

        compounds_list = data.get("compounds", [])
        if not compounds_list:
            time.sleep(MIN_DELAY)
            continue

        batch: List[Tuple] = []
        # Each compound match has a "sources" list with cross-references
        for compound_match in compounds_list:
            sources = compound_match.get("sources", [])
            for source in sources:
                src_url = (source.get("url") or "").lower()
                src_compound_id = source.get("compoundId", "")
                src_name_raw = source.get("name", "")

                if not src_compound_id:
                    continue

                # Identify the database
                db_name = None
                xref_url = ""
                for key, (name, url_tpl) in known_dbs.items():
                    if key in src_url or key in src_name_raw.lower():
                        db_name = name
                        xref_url = url_tpl.format(src_compound_id)
                        break

                if not db_name:
                    continue

                # Skip self-references (ChEMBL -> ChEMBL)
                if db_name == "ChEMBL" and src_compound_id == chembl_id:
                    continue

                batch.append((
                    compound_id,
                    None,  # target_id
                    db_name,
                    src_compound_id,
                    xref_url,
                ))

        if batch:
            with conn.cursor() as cur:
                psycopg2.extras.execute_values(cur, """
                    INSERT INTO pharmaco_db.cross_references
                    (compound_id, target_id, db_name, db_id, url)
                    VALUES %s
                    ON CONFLICT (compound_id, target_id, db_name, db_id) DO NOTHING
                """, batch, page_size=BATCH_PAGE_SIZE)
            conn.commit()
            total += len(batch)

        time.sleep(MIN_DELAY)

    log(f"  DONE: {total} UniChem cross-references inserted ({errors} lookup failures)")
    update_log(conn, log_id, "completed", total)
    return total


# ============================================================
# 6. ChEBI -- Chemical Ontology Classification
# ============================================================

def ingest_chebi(conn: psycopg2.extensions.connection) -> int:
    """
    Enrich compounds with ChEBI ontology data.
    We use the ChEBI web services API to fetch roles/classifications.
    API: https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:ID (returns XML)
    Better: use ChEBI OLS API for ontology:
    https://www.ebi.ac.uk/ols4/api/ontologies/chebi/terms?obo_id=CHEBI:ID
    Strategy: look up ChEBI IDs from our cross_references table (populated by UniChem).
    """
    log("=" * 60)
    log("6. ChEBI -- Chemical Ontology Classification")
    log("=" * 60)
    log_id = start_log(conn, "chebi", "ontology")

    # Get compounds that have ChEBI cross-references
    chebi_compounds: List[Tuple[int, str]] = []
    with conn.cursor() as cur:
        cur.execute("""
            SELECT compound_id, db_id FROM pharmaco_db.cross_references
            WHERE db_name = 'ChEBI'
            AND compound_id IS NOT NULL
            ORDER BY compound_id
            LIMIT 5000
        """)
        chebi_compounds = cur.fetchall()

    if not chebi_compounds:
        log("  No ChEBI cross-references found (run UniChem first)")
        # Try getting ChEBI IDs from ChEMBL compounds directly
        with conn.cursor() as cur:
            cur.execute("""
                SELECT id, chembl_id FROM pharmaco_db.compounds
                WHERE max_phase >= 4
                LIMIT 500
            """)
            fallback = cur.fetchall()

        if not fallback:
            update_log(conn, log_id, "completed", 0)
            return 0

        log(f"  Fallback: fetching ChEBI IDs for {len(fallback)} approved drugs via ChEMBL API")
        chebi_compounds = []
        for cid, chembl_id in fallback:
            url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
            data = api_get(url, timeout=10)
            if data:
                xrefs = data.get("cross_references", [])
                for xr in xrefs:
                    if xr.get("xref_src") == "ChEBI":
                        chebi_id = xr.get("xref_id")
                        if chebi_id:
                            chebi_compounds.append((cid, chebi_id))
                            break
            time.sleep(MIN_DELAY)

    log(f"  Processing {len(chebi_compounds)} compounds with ChEBI IDs")
    total = 0

    for i, (compound_id, chebi_id) in enumerate(chebi_compounds):
        if i % 100 == 0 and i > 0:
            log(f"    Progress: {i}/{len(chebi_compounds)} ({total} annotations)")

        # Clean up ChEBI ID format
        chebi_clean = chebi_id.replace("CHEBI:", "")
        obo_id = f"CHEBI:{chebi_clean}"

        # Use OLS4 API to get ontology info
        url = (
            f"https://www.ebi.ac.uk/ols4/api/ontologies/chebi/terms"
            f"?obo_id={urllib.parse.quote(obo_id)}"
        )
        data = api_get(url, timeout=15)
        if not data:
            time.sleep(MIN_DELAY)
            continue

        # Parse embedded terms
        terms = data.get("_embedded", {}).get("terms", [])
        if not terms:
            time.sleep(MIN_DELAY)
            continue

        term = terms[0]
        label = term.get("label", "")
        description_list = term.get("description", [])
        description = description_list[0] if description_list else ""

        # Add cross-reference with ontology info
        batch: List[Tuple] = []
        batch.append((
            compound_id,
            None,
            "ChEBI_ontology",
            obo_id,
            f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={obo_id}",
        ))

        # Also try to get the "has role" classification
        # This is available in the annotations
        annotations = term.get("annotation", {})
        if annotations:
            # Store functional classification info
            for key in ["has_role", "is_a"]:
                roles = annotations.get(key, [])
                for role in roles[:5]:
                    batch.append((
                        compound_id,
                        None,
                        f"ChEBI_{key}",
                        str(role),
                        "",
                    ))

        if batch:
            with conn.cursor() as cur:
                psycopg2.extras.execute_values(cur, """
                    INSERT INTO pharmaco_db.cross_references
                    (compound_id, target_id, db_name, db_id, url)
                    VALUES %s
                    ON CONFLICT (compound_id, target_id, db_name, db_id) DO NOTHING
                """, batch, page_size=BATCH_PAGE_SIZE)
            conn.commit()
            total += len(batch)

        time.sleep(MIN_DELAY)

    log(f"  DONE: {total} ChEBI ontology annotations inserted")
    update_log(conn, log_id, "completed", total)
    return total


# ============================================================
# 7. STRING -- Protein-Protein Interactions
# ============================================================

def ingest_string(conn: psycopg2.extensions.connection) -> int:
    """
    Ingest protein-protein interactions from STRING.
    API: https://string-db.org/api/json/interaction_partners
    Strategy: batch query using newline-separated gene identifiers.
    STRING supports batch queries with multiple identifiers.
    """
    log("=" * 60)
    log("7. STRING -- Protein-Protein Interactions")
    log("=" * 60)
    log_id = start_log(conn, "string", "protein_interactions")

    target_gene_map = build_target_gene_map(conn)
    targets_with_gene = build_targets_with_gene(conn)
    log(f"  Found {len(targets_with_gene)} human targets with gene names")

    total = 0
    errors = 0
    processed = 0

    # STRING supports batching, but we query one at a time to get partner details
    for target_id, gene_name in targets_with_gene:
        processed += 1
        if processed % 100 == 0:
            log(f"    Progress: {processed}/{len(targets_with_gene)} ({total} interactions)")

        url = (
            f"https://string-db.org/api/json/interaction_partners"
            f"?identifiers={urllib.parse.quote(gene_name)}"
            f"&species=9606"
            f"&limit=30"
            f"&required_score=700"  # Only high-confidence (>=0.7)
        )

        data = api_get(url, timeout=20)
        if not data or not isinstance(data, list):
            errors += 1
            time.sleep(MIN_DELAY)
            continue

        batch: List[Tuple] = []
        for interaction in data:
            # STRING returns preferredName fields
            gene_a = interaction.get("preferredName_A", "")
            gene_b = interaction.get("preferredName_B", "")

            # Normalize: always store alphabetically smaller gene first
            if gene_a.upper() > gene_b.upper():
                gene_a, gene_b = gene_b, gene_a

            target_id_1 = target_gene_map.get(gene_a.upper())
            target_id_2 = target_gene_map.get(gene_b.upper())

            combined = interaction.get("score")
            # STRING scores are 0-1 float, convert to 0-1000 integer
            combined_int = int(combined * 1000) if combined else None

            escore = interaction.get("escore")
            dscore = interaction.get("dscore")
            tscore = interaction.get("tscore")
            ascore = interaction.get("ascore")

            experimental_int = int(float(escore) * 1000) if escore else None
            database_int = int(float(dscore) * 1000) if dscore else None
            textmining_int = int(float(tscore) * 1000) if tscore else None
            coexpression_int = int(float(ascore) * 1000) if ascore else None

            batch.append((
                target_id_1,
                target_id_2,
                gene_a.upper(),
                gene_b.upper(),
                combined_int,
                experimental_int,
                database_int,
                textmining_int,
                coexpression_int,
                "string",
            ))

        if batch:
            with conn.cursor() as cur:
                psycopg2.extras.execute_values(cur, """
                    INSERT INTO pharmaco_db.protein_interactions
                    (target_id_1, target_id_2, gene_1, gene_2,
                     combined_score, experimental_score, database_score,
                     textmining_score, coexpression_score, source)
                    VALUES %s
                    ON CONFLICT (gene_1, gene_2)
                    WHERE gene_1 < gene_2
                    DO NOTHING
                """, batch, page_size=BATCH_PAGE_SIZE)
            conn.commit()
            total += len(batch)

        # STRING rate limit: be respectful
        time.sleep(MIN_DELAY)

    log(f"  DONE: {total} STRING interactions inserted ({errors} query failures)")
    update_log(conn, log_id, "completed", total)
    return total


# ============================================================
# 8. Reactome -- Biological Pathways
# ============================================================

def ingest_reactome(conn: psycopg2.extensions.connection) -> int:
    """
    Ingest pathway data from Reactome for our targets.
    API: https://reactome.org/ContentService/data/mapping/UniProt/{ID}/pathways?species=9606
    Returns pathways a protein participates in.
    """
    log("=" * 60)
    log("8. Reactome -- Biological Pathways")
    log("=" * 60)
    log_id = start_log(conn, "reactome", "pathways")

    targets_with_uniprot = build_targets_with_uniprot(conn)
    log(f"  Found {len(targets_with_uniprot)} targets with UniProt IDs")

    total = 0
    errors = 0
    processed = 0

    for target_id, uniprot_id, gene_name in targets_with_uniprot:
        processed += 1
        if processed % 100 == 0:
            log(f"    Progress: {processed}/{len(targets_with_uniprot)} ({total} pathways)")

        # Reactome mapping endpoint (verified working)
        url = (
            f"https://reactome.org/ContentService/data/mapping"
            f"/UniProt/{uniprot_id}/pathways?species=9606"
        )

        data = api_get(url, timeout=15, headers={"Accept": "application/json"})
        if not data or not isinstance(data, list):
            errors += 1
            time.sleep(MIN_DELAY)
            continue

        batch: List[Tuple] = []
        for pathway in data:
            pathway_id = pathway.get("stId", "")
            pathway_name = pathway.get("displayName", "")
            species = pathway.get("speciesName", "Homo sapiens")

            # Reactome pathways have a schemaClass (Pathway, TopLevelPathway, etc.)
            category = pathway.get("schemaClass", "")

            batch.append((
                target_id,
                uniprot_id,
                pathway_id,
                pathway_name,
                category,
                species,
                "reactome",
            ))

        if batch:
            with conn.cursor() as cur:
                psycopg2.extras.execute_values(cur, """
                    INSERT INTO pharmaco_db.pathways
                    (target_id, uniprot_id, pathway_id, pathway_name,
                     pathway_category, species, source)
                    VALUES %s
                    ON CONFLICT (target_id, pathway_id, source)
                    WHERE target_id IS NOT NULL
                    DO NOTHING
                """, batch, page_size=BATCH_PAGE_SIZE)
            conn.commit()
            total += len(batch)

        time.sleep(MIN_DELAY)

    log(f"  DONE: {total} Reactome pathways inserted ({errors} query failures)")
    update_log(conn, log_id, "completed", total)
    return total


# ============================================================
# 9. DisGeNET -- Disease-Gene Associations
# ============================================================

def ingest_disgenet(conn: psycopg2.extensions.connection) -> int:
    """
    Ingest disease-gene associations from DisGeNET.
    Uses the free public data via the REST API (limited access without key).
    API: https://www.disgenet.org/api/gda/gene/{NCBI_GENE_ID}
    Alternative: use the public summary endpoint or curated TSV files.

    Since the API requires registration, we use the CURATED gene-disease
    associations available from the DisGeNET downloads page.
    Fallback: query via gene symbol using the public API endpoint.
    """
    log("=" * 60)
    log("9. DisGeNET -- Disease-Gene Associations")
    log("=" * 60)
    log_id = start_log(conn, "disgenet", "disease_associations")

    target_gene_map = build_target_gene_map(conn)

    # Get all our targets with NCBI gene IDs
    targets_with_ncbi: List[Tuple[int, str, int]] = []
    with conn.cursor() as cur:
        cur.execute("""
            SELECT id, gene_name, ncbi_gene_id FROM pharmaco_db.targets
            WHERE ncbi_gene_id IS NOT NULL
            AND gene_name IS NOT NULL
            AND organism = 'Homo sapiens'
            ORDER BY id
        """)
        targets_with_ncbi = cur.fetchall()

    log(f"  Found {len(targets_with_ncbi)} targets with NCBI gene IDs")

    # Try the DisGeNET public API first (no key required for basic queries)
    total = 0
    errors = 0
    processed = 0

    for target_id, gene_name, ncbi_gene_id in targets_with_ncbi:
        processed += 1
        if processed % 100 == 0:
            log(f"    Progress: {processed}/{len(targets_with_ncbi)} ({total} associations)")

        # DisGeNET API -- try gene search
        # The public API may require an API key; try without first
        url = (
            f"https://www.disgenet.org/api/gda/gene/{ncbi_gene_id}"
            f"?source=CURATED&format=json&limit=25"
        )

        data = api_get(url, timeout=15, headers={
            "Accept": "application/json",
        })

        if data is None:
            # API likely requires key. Try alternative: query the summary endpoint
            # Fallback: use the v7 OpenAPI which has less restrictions
            url2 = (
                f"https://www.disgenet.org/api/gda/summary"
                f"?gene={ncbi_gene_id}&source=CURATED&format=json"
            )
            data = api_get(url2, timeout=15)

        if not data or not isinstance(data, list):
            errors += 1
            time.sleep(MIN_DELAY)
            continue

        batch: List[Tuple] = []
        for assoc in data[:25]:  # Top 25 per gene
            disease_id = assoc.get("diseaseid", assoc.get("diseaseId", ""))
            disease_name = assoc.get("disease_name", assoc.get("diseaseName", ""))
            disease_type = assoc.get("disease_type", assoc.get("diseaseType", ""))
            score = assoc.get("score")
            ei = assoc.get("ei", assoc.get("evidenceIndex"))
            n_pmids = assoc.get("nofPmids", assoc.get("pmid_count", 0))
            n_snps = assoc.get("nofSnps", assoc.get("snp_count", 0))

            if score is not None:
                try:
                    score = float(score)
                except (ValueError, TypeError):
                    score = None

            if ei is not None:
                try:
                    ei = float(ei)
                except (ValueError, TypeError):
                    ei = None

            batch.append((
                target_id,
                gene_name,
                ncbi_gene_id,
                disease_id,
                disease_name,
                disease_type,
                score,
                ei,
                n_pmids,
                n_snps,
                "disgenet",
            ))

        if batch:
            with conn.cursor() as cur:
                psycopg2.extras.execute_values(cur, """
                    INSERT INTO pharmaco_db.disgenet_associations
                    (target_id, gene_name, gene_ncbi_id, disease_id,
                     disease_name, disease_type, association_score,
                     ei_score, num_pmids, num_snps, source)
                    VALUES %s
                    ON CONFLICT (gene_name, disease_id) DO NOTHING
                """, batch, page_size=BATCH_PAGE_SIZE)
            conn.commit()
            total += len(batch)

        time.sleep(MIN_DELAY)

    # If API failed for most targets, note it
    if errors > len(targets_with_ncbi) * 0.9:
        log("  WARNING: DisGeNET API may require an API key for bulk access.")
        log("  Register at https://www.disgenet.org/api/#/Authorization to get a key.")
        log("  Set DISGENET_API_KEY environment variable and re-run.")

    log(f"  DONE: {total} DisGeNET associations inserted ({errors} query failures)")
    update_log(conn, log_id, "completed", total)
    return total


# ============================================================
# 10. KEGG -- Pathway & Drug Info
# ============================================================

def ingest_kegg(conn: psycopg2.extensions.connection) -> int:
    """
    Ingest drug and pathway info from KEGG.
    APIs:
      - https://rest.kegg.jp/find/drug/DRUG_NAME  (search by name)
      - https://rest.kegg.jp/get/DRUG_ID  (get drug details)
      - https://rest.kegg.jp/link/pathway/hsa:GENE_ID  (gene to pathways)
      - https://rest.kegg.jp/get/PATHWAY_ID  (pathway details)

    KEGG returns flat-text responses, not JSON.
    """
    log("=" * 60)
    log("10. KEGG -- Drug & Pathway Info")
    log("=" * 60)
    log_id = start_log(conn, "kegg", "drugs_and_pathways")

    total = 0
    drug_total = 0
    pathway_total = 0

    # ---- Part A: KEGG Drug Info for approved drugs ----
    log("  Part A: KEGG Drug Info")
    approved_drugs = build_approved_drugs(conn)
    log(f"    Searching KEGG for {len(approved_drugs)} approved drugs")

    for i, (compound_id, drug_name) in enumerate(approved_drugs):
        if i % 50 == 0 and i > 0:
            log(f"      Progress: {i}/{len(approved_drugs)} drugs ({drug_total} KEGG entries)")

        # Search KEGG for the drug name
        search_url = f"https://rest.kegg.jp/find/drug/{urllib.parse.quote(drug_name)}"
        text = fetch_text(search_url, timeout=15)
        if not text or not text.strip():
            time.sleep(MIN_DELAY)
            continue

        # Parse: each line is "DRUG_ID\tDRUG_NAME; ..."
        kegg_drug_id = None
        for line in text.strip().split("\n"):
            parts = line.split("\t", 1)
            if len(parts) >= 1:
                kegg_drug_id = parts[0].strip()
                break

        if not kegg_drug_id:
            time.sleep(MIN_DELAY)
            continue

        # Get drug details
        detail_url = f"https://rest.kegg.jp/get/{kegg_drug_id}"
        detail_text = fetch_text(detail_url, timeout=15)
        if not detail_text:
            time.sleep(MIN_DELAY)
            continue

        # Parse KEGG flat-file format
        kegg_name = _kegg_extract_field(detail_text, "NAME")
        kegg_formula = _kegg_extract_field(detail_text, "FORMULA")
        categories = _kegg_extract_list(detail_text, "PRODUCT")
        pathway_ids_raw = _kegg_extract_list(detail_text, "PATHWAY")
        pathway_names_raw = _kegg_extract_list(detail_text, "PATHWAY")

        # Extract pathway IDs more carefully
        pathway_ids = []
        pathway_names = []
        for pw_line in pathway_ids_raw:
            # Format: "map12345  Pathway name"
            pw_match = re.match(r"(map\d+|hsa\d+)\s+(.*)", pw_line.strip())
            if pw_match:
                pathway_ids.append(pw_match.group(1))
                pathway_names.append(pw_match.group(2))

        # Therapeutic categories
        therapeutic = _kegg_extract_list(detail_text, "TARGET")

        with conn.cursor() as cur:
            cur.execute("""
                INSERT INTO pharmaco_db.kegg_drug_info
                (compound_id, kegg_drug_id, kegg_name, kegg_formula,
                 therapeutic_category, pathway_ids, pathway_names)
                VALUES (%s, %s, %s, %s, %s, %s, %s)
                ON CONFLICT (kegg_drug_id) WHERE kegg_drug_id IS NOT NULL
                DO NOTHING
            """, (
                compound_id,
                kegg_drug_id,
                kegg_name,
                kegg_formula,
                therapeutic if therapeutic else None,
                pathway_ids if pathway_ids else None,
                pathway_names if pathway_names else None,
            ))
        conn.commit()
        drug_total += 1

        time.sleep(MIN_DELAY)

    # ---- Part B: KEGG Pathways for our gene targets ----
    log("  Part B: KEGG Pathways for targets")
    targets_with_ncbi: List[Tuple[int, str, int]] = []
    with conn.cursor() as cur:
        cur.execute("""
            SELECT id, gene_name, ncbi_gene_id FROM pharmaco_db.targets
            WHERE ncbi_gene_id IS NOT NULL
            AND organism = 'Homo sapiens'
            ORDER BY id
        """)
        targets_with_ncbi = cur.fetchall()

    log(f"    Querying KEGG pathways for {len(targets_with_ncbi)} targets")

    for i, (target_id, gene_name, ncbi_gene_id) in enumerate(targets_with_ncbi):
        if i % 100 == 0 and i > 0:
            log(f"      Progress: {i}/{len(targets_with_ncbi)} ({pathway_total} pathways)")

        # KEGG gene-to-pathway mapping
        link_url = f"https://rest.kegg.jp/link/pathway/hsa:{ncbi_gene_id}"
        text = fetch_text(link_url, timeout=15)
        if not text or not text.strip():
            time.sleep(MIN_DELAY)
            continue

        pathway_ids_for_target = []
        for line in text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) >= 2:
                pw_id = parts[1].strip()
                if pw_id.startswith("path:"):
                    pw_id = pw_id[5:]
                pathway_ids_for_target.append(pw_id)

        # Fetch pathway names in batch
        batch: List[Tuple] = []
        for pw_id in pathway_ids_for_target[:20]:  # Limit to 20 pathways per target
            # Get pathway name
            pw_url = f"https://rest.kegg.jp/get/{pw_id}"
            pw_text = fetch_text(pw_url, timeout=10)
            pw_name = ""
            if pw_text:
                pw_name = _kegg_extract_field(pw_text, "NAME")
                if pw_name:
                    # Remove " - Homo sapiens (human)" suffix
                    pw_name = re.sub(r"\s*-\s*Homo sapiens.*$", "", pw_name)

            batch.append((
                target_id,
                None,  # uniprot_id (we link via target_id)
                pw_id,
                pw_name,
                "kegg_pathway",
                "Homo sapiens",
                "kegg",
            ))
            time.sleep(0.15)  # KEGG is sensitive to rate

        if batch:
            with conn.cursor() as cur:
                psycopg2.extras.execute_values(cur, """
                    INSERT INTO pharmaco_db.pathways
                    (target_id, uniprot_id, pathway_id, pathway_name,
                     pathway_category, species, source)
                    VALUES %s
                    ON CONFLICT (target_id, pathway_id, source)
                    WHERE target_id IS NOT NULL
                    DO NOTHING
                """, batch, page_size=BATCH_PAGE_SIZE)
            conn.commit()
            pathway_total += len(batch)

        time.sleep(MIN_DELAY)

    total = drug_total + pathway_total
    log(f"  DONE: {drug_total} KEGG drugs + {pathway_total} KEGG pathways = {total} total")
    update_log(conn, log_id, "completed", total)
    return total


def _kegg_extract_field(text: str, field: str) -> str:
    """Extract a single field value from KEGG flat-file text."""
    pattern = re.compile(rf"^{field}\s+(.*?)$", re.MULTILINE)
    match = pattern.search(text)
    if match:
        value = match.group(1).strip()
        # Remove trailing semicolons
        return value.rstrip(";").strip()
    return ""


def _kegg_extract_list(text: str, field: str) -> List[str]:
    """Extract a multi-line field from KEGG flat-file text.
    KEGG uses 12-space indentation for continuation lines."""
    result: List[str] = []
    in_field = False
    for line in text.split("\n"):
        if line.startswith(field):
            in_field = True
            # Value is after the field name + spaces
            value = line[len(field):].strip()
            if value:
                result.append(value)
        elif in_field:
            if line.startswith("            "):  # 12 spaces = continuation
                result.append(line.strip())
            elif line.startswith(" "):
                # Could still be continuation with different spacing
                stripped = line.strip()
                if stripped and not any(
                    line.startswith(f) for f in [
                        "ENTRY", "NAME", "FORMULA", "EXACT_MASS", "MOL_WEIGHT",
                        "REMARK", "ACTIVITY", "TARGET", "PATHWAY", "INTERACTION",
                        "PRODUCT", "SOURCE", "COMPONENT", "DBLINKS", "ATOM",
                        "BOND", "BRACKET", "COMMENT", "BRITE", "EFFICACY",
                        "DISEASE", "STR_MAP", "OTHER_MAP", "METABOLISM",
                    ]
                ):
                    result.append(stripped)
                else:
                    in_field = False
            else:
                in_field = False
    return result


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    start_time = time.time()
    log("=" * 70)
    log("PharmacoDB -- PHASE 3: Additional Public Databases")
    log(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log("=" * 70)

    conn = get_conn()

    # Apply schema first
    log("\n>>> Applying Phase 3 schema...")
    schema_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "schema_phase3.sql")
    if os.path.exists(schema_path):
        with open(schema_path, "r") as f:
            schema_sql = f.read()
        with conn.cursor() as cur:
            cur.execute(schema_sql)
        conn.commit()
        log("  Schema applied successfully")
    else:
        log(f"  WARNING: Schema file not found at {schema_path}")
        log("  Tables must already exist in the database")

    # Get current DB stats
    with conn.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.compounds")
        n_compounds = cur.fetchone()[0]
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.targets")
        n_targets = cur.fetchone()[0]
    log(f"\n  Current DB: {n_compounds:,} compounds, {n_targets:,} targets")

    results: Dict[str, int] = {}

    # Run all 10 ingestion sources. Each source is isolated so failures
    # in one do not prevent the others from running.
    sources = [
        ("DGIdb", ingest_dgidb),
        ("SIDER", ingest_sider),
        ("OpenFDA", ingest_openfda),
        ("ClinicalTrials.gov", ingest_clinical_trials),
        ("UniChem", ingest_unichem),
        ("ChEBI", ingest_chebi),
        ("STRING", ingest_string),
        ("Reactome", ingest_reactome),
        ("DisGeNET", ingest_disgenet),
        ("KEGG", ingest_kegg),
    ]

    for source_name, source_fn in sources:
        log(f"\n>>> STARTING: {source_name}")
        try:
            count = source_fn(conn)
            results[source_name] = count
        except Exception as e:
            log(f"  FATAL ERROR in {source_name}: {e}")
            traceback.print_exc()
            results[source_name] = -1
            # Rollback any pending transaction
            try:
                conn.rollback()
            except Exception:
                pass
            # Reconnect if needed
            try:
                conn.cursor().execute("SELECT 1")
            except Exception:
                log(f"  Reconnecting to database...")
                conn = get_conn()

    # ---- FINAL STATISTICS ----
    log("\n" + "=" * 70)
    log("  PHASE 3 INGESTION -- FINAL REPORT")
    log("=" * 70)
    grand_total = 0
    for source_name, count in results.items():
        status = f"{count:>10,}" if count >= 0 else "    FAILED"
        log(f"    {source_name:<25}: {status}")
        if count > 0:
            grand_total += count

    log(f"\n    {'GRAND TOTAL':<25}: {grand_total:>10,}")

    # Print new table counts
    log("\n  New table row counts:")
    new_tables = [
        ("dgidb_interactions", "pharmaco_db.dgidb_interactions"),
        ("side_effects", "pharmaco_db.side_effects"),
        ("clinical_trials", "pharmaco_db.clinical_trials"),
        ("adverse_events", "pharmaco_db.adverse_events"),
        ("protein_interactions", "pharmaco_db.protein_interactions"),
        ("pathways", "pharmaco_db.pathways"),
        ("kegg_drug_info", "pharmaco_db.kegg_drug_info"),
        ("disgenet_associations", "pharmaco_db.disgenet_associations"),
        ("cross_references", "pharmaco_db.cross_references"),
    ]
    for name, table in new_tables:
        try:
            with conn.cursor() as cur:
                cur.execute(f"SELECT COUNT(*) FROM {table}")
                cnt = cur.fetchone()[0]
                log(f"    {name:<25}: {cnt:>10,}")
        except Exception:
            conn.rollback()
            log(f"    {name:<25}:  (error)")

    # DB size
    try:
        with conn.cursor() as cur:
            cur.execute("SELECT pg_size_pretty(pg_database_size('pharmaco'))")
            db_size = cur.fetchone()[0]
            log(f"\n  Database size: {db_size}")
    except Exception:
        conn.rollback()

    conn.close()

    elapsed = time.time() - start_time
    hours = int(elapsed // 3600)
    minutes = int((elapsed % 3600) // 60)
    seconds = int(elapsed % 60)
    log(f"\n  Phase 3 total time: {hours}h {minutes}m {seconds}s")
    log("  PHASE 3 COMPLETE")
    log("=" * 70)


if __name__ == "__main__":
    main()
