#!/usr/bin/env python3
"""
MegaPharmaDB -- Phase 3 Consolidation
======================================
Fills the remaining empty tables in PharmacoDB (7 of 28):

1. KEGG Drug Info           -- KEGG REST API
2. OpenFDA Adverse Events   -- OpenFDA drug/event API
3. ClinicalTrials.gov       -- CT.gov API v2
4. Open Targets Diseases    -- GraphQL API
5. DisGeNET Associations    -- Open data TSV (curated GDAs)
6. Tox21 Toxicology         -- PubChem PUG REST (Tox21 assays)
7. Compound Structures      -- Continue RDKit fingerprint computation

Skipped: pharmacokinetics (no bulk public source)

Usage:
    python3 mega_consolidate_p3.py          # Run all
    python3 mega_consolidate_p3.py kegg     # Run one step
    python3 mega_consolidate_p3.py structures --limit 50000  # Batch limit

DB: localhost:5433, pharmaco, postgres (no password)
"""

from __future__ import annotations

import gzip
import io
import json
import csv
import os
import sys
import time
import traceback
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime
from pathlib import Path
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
    "password": "postgres",
}

DOWNLOAD_DIR: Path = Path("/tmp/megapharma_p3_downloads")
BATCH: int = 5000
USER_AGENT: str = "MegaPharmaDB/3.0 (drug-discovery-research)"


# ============================================================
# HELPERS
# ============================================================

def log(msg: str) -> None:
    """Print timestamped log message."""
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def get_conn() -> psycopg2.extensions.connection:
    """Create a new database connection."""
    return psycopg2.connect(**DB_CONFIG)


def fetch_url(
    url: str,
    data: Optional[bytes] = None,
    headers: Optional[Dict[str, str]] = None,
    timeout: int = 30,
    retries: int = 3,
) -> bytes:
    """Fetch URL with retries and exponential backoff."""
    hdrs: Dict[str, str] = {"User-Agent": USER_AGENT}
    if headers:
        hdrs.update(headers)
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, data=data, headers=hdrs)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as e:
            if attempt == retries - 1:
                raise
            wait: float = 2 ** attempt + 0.5
            log(f"    Retry {attempt + 1}/{retries} after {wait:.0f}s: {e}")
            time.sleep(wait)
    return b""  # unreachable


def fetch_json(
    url: str,
    data: Optional[Dict[str, Any]] = None,
    headers: Optional[Dict[str, str]] = None,
    timeout: int = 30,
    retries: int = 3,
) -> Any:
    """Fetch JSON from URL. If data is provided, POST as JSON."""
    hdrs: Dict[str, str] = {"Content-Type": "application/json"}
    if headers:
        hdrs.update(headers)
    raw_data: Optional[bytes] = json.dumps(data).encode() if data else None
    raw: bytes = fetch_url(url, data=raw_data, headers=hdrs, timeout=timeout, retries=retries)
    return json.loads(raw)


def download_file(url: str, dest: Path, desc: str = "") -> Optional[Path]:
    """Download a file to disk with caching."""
    if dest.exists() and dest.stat().st_size > 100:
        log(f"  Already downloaded: {dest.name} ({dest.stat().st_size / 1024 / 1024:.1f}MB)")
        return dest
    log(f"  Downloading {desc or url}...")
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        raw: bytes = fetch_url(url, timeout=300, retries=2)
        with open(dest, "wb") as f:
            f.write(raw)
        log(f"  Downloaded {len(raw) / 1024 / 1024:.1f} MB")
        return dest
    except Exception as e:
        log(f"  DOWNLOAD FAILED: {e}")
        if dest.exists():
            dest.unlink()
        return None


def table_count(conn: psycopg2.extensions.connection, table: str) -> int:
    """Get row count for a table."""
    with conn.cursor() as cur:
        cur.execute(f"SELECT COUNT(*) FROM pharmaco_db.{table}")
        return cur.fetchone()[0]


# ============================================================
# 1. KEGG DRUG INFO
# ============================================================

def fill_kegg_drug_info() -> int:
    """
    Fill kegg_drug_info table using KEGG REST API.

    Strategy:
    1. Download KEGG drug list (all drug IDs + names)
    2. Match KEGG drug names to our compounds.pref_name
    3. For matched drugs, fetch detailed info (formula, categories, pathways)
    """
    log("=" * 60)
    log("1. KEGG Drug Info")
    log("=" * 60)

    conn: psycopg2.extensions.connection = get_conn()

    # Check existing rows
    existing: int = table_count(conn, "kegg_drug_info")
    log(f"  Existing rows: {existing:,}")

    # Build compound name -> id mapping (approved drugs only)
    with conn.cursor() as cur:
        cur.execute("""
            SELECT UPPER(pref_name), id FROM pharmaco_db.compounds
            WHERE pref_name IS NOT NULL
        """)
        name_to_cid: Dict[str, int] = dict(cur.fetchall())
    log(f"  Compounds with names: {len(name_to_cid):,}")

    # Step 1: Get KEGG drug list
    log("  Fetching KEGG drug list...")
    try:
        raw: bytes = fetch_url("https://rest.kegg.jp/list/drug", timeout=60)
        kegg_list: str = raw.decode("utf-8")
    except Exception as e:
        log(f"  SKIP: KEGG list API error: {e}")
        conn.close()
        return 0

    # Parse drug list: "dr:D00001\tLepirudin (INN); ..."
    kegg_drugs: List[Tuple[str, str, List[str]]] = []  # (kegg_id, primary_name, all_names)
    for line in kegg_list.strip().split("\n"):
        parts: List[str] = line.split("\t", 1)
        if len(parts) < 2:
            continue
        kegg_id: str = parts[0].replace("dr:", "")
        raw_names: str = parts[1]
        # Names are semicolon-separated, some have annotations in parens
        all_names: List[str] = [n.strip() for n in raw_names.split(";")]
        primary_name: str = all_names[0] if all_names else ""
        kegg_drugs.append((kegg_id, primary_name, all_names))

    log(f"  KEGG drugs in list: {len(kegg_drugs):,}")

    # Step 2: Match to our compounds by name
    matched: List[Tuple[int, str, str]] = []  # (compound_id, kegg_id, kegg_name)
    unmatched_kegg: List[Tuple[str, str]] = []

    for kegg_id, primary_name, all_names in kegg_drugs:
        cid: Optional[int] = None
        matched_name: str = primary_name
        for name in all_names:
            # Strip annotations like "(INN)", "(JAN)", "(USAN)"
            clean: str = name.strip()
            for suffix in ["(INN)", "(JAN)", "(USAN)", "(BAN)", "(NF)", "(USP)"]:
                clean = clean.replace(suffix, "").strip()
            cid = name_to_cid.get(clean.upper())
            if cid:
                matched_name = clean
                break
            # Try without any parenthetical
            idx: int = clean.find("(")
            if idx > 0:
                cid = name_to_cid.get(clean[:idx].strip().upper())
                if cid:
                    matched_name = clean[:idx].strip()
                    break
        if cid:
            matched.append((cid, kegg_id, matched_name))
        else:
            unmatched_kegg.append((kegg_id, primary_name))

    log(f"  Matched to compounds: {len(matched):,}/{len(kegg_drugs):,}")

    # Step 3: Fetch detailed info for matched drugs (formula, categories, pathways)
    # Rate limit: KEGG allows ~10 req/s, we do 5/s to be safe
    log(f"  Fetching details for {len(matched):,} matched drugs...")
    insert_rows: List[Tuple] = []
    fetched: int = 0
    errors: int = 0

    for i, (cid, kegg_id, kegg_name) in enumerate(matched):
        formula: Optional[str] = None
        categories: List[str] = []
        pathway_ids: List[str] = []
        pathway_names: List[str] = []

        try:
            detail_url: str = f"https://rest.kegg.jp/get/{kegg_id}"
            raw_detail: bytes = fetch_url(detail_url, timeout=15, retries=2)
            detail_text: str = raw_detail.decode("utf-8")

            # Parse KEGG flat-file format
            current_field: str = ""
            for line in detail_text.split("\n"):
                if line.startswith("FORMULA"):
                    formula = line.split(None, 1)[1].strip() if len(line.split(None, 1)) > 1 else None
                elif line.startswith("PRODUCT"):
                    current_field = "PRODUCT"
                elif line.startswith("PATHWAY"):
                    current_field = "PATHWAY"
                    # First pathway line: "PATHWAY     map00010  Glycolysis..."
                    rest: str = line[12:].strip()
                    if rest:
                        parts_pw = rest.split(None, 1)
                        if parts_pw:
                            pathway_ids.append(parts_pw[0])
                            if len(parts_pw) > 1:
                                pathway_names.append(parts_pw[1])
                elif line.startswith("TARGET") or line.startswith("METABOLISM"):
                    current_field = ""
                elif line.startswith("CATEGORY"):
                    current_field = "CATEGORY"
                elif line.startswith("  ") and current_field == "PATHWAY":
                    rest = line.strip()
                    parts_pw = rest.split(None, 1)
                    if parts_pw:
                        pathway_ids.append(parts_pw[0])
                        if len(parts_pw) > 1:
                            pathway_names.append(parts_pw[1])
                elif line.startswith("  ") and current_field == "CATEGORY":
                    cat: str = line.strip()
                    if cat:
                        categories.append(cat)
                elif line.startswith("///"):
                    break
                elif not line.startswith(" "):
                    current_field = ""

            fetched += 1
        except Exception:
            errors += 1
            # Still insert basic info without details

        insert_rows.append((
            cid,
            kegg_id,
            kegg_name,
            formula,
            categories if categories else None,
            pathway_ids if pathway_ids else None,
            pathway_names if pathway_names else None,
        ))

        # Rate limit: 5 requests per second
        if (i + 1) % 5 == 0:
            time.sleep(1.0)

        if (i + 1) % 200 == 0:
            log(f"    Progress: {i + 1}/{len(matched)} (fetched={fetched}, errors={errors})")

    # Also add unmatched KEGG drugs with compound_id=NULL (for future mapping)
    # Only top-level info, no detail fetch
    for kegg_id, kegg_name in unmatched_kegg:
        insert_rows.append((None, kegg_id, kegg_name, None, None, None, None))

    # Bulk insert
    if insert_rows:
        with conn.cursor() as cur:
            # Clear and re-insert for idempotency
            cur.execute("TRUNCATE pharmaco_db.kegg_drug_info")
            psycopg2.extras.execute_values(
                cur,
                """
                INSERT INTO pharmaco_db.kegg_drug_info
                    (compound_id, kegg_drug_id, kegg_name, kegg_formula,
                     therapeutic_category, pathway_ids, pathway_names)
                VALUES %s
                ON CONFLICT DO NOTHING
                """,
                insert_rows,
                page_size=2000,
            )
        conn.commit()

    total_inserted: int = table_count(conn, "kegg_drug_info")
    conn.close()
    log(f"  KEGG: {total_inserted:,} rows ({len(matched):,} matched to compounds)")
    return total_inserted


# ============================================================
# 2. OPENFDA ADVERSE EVENTS
# ============================================================

def fill_adverse_events() -> int:
    """
    Fill adverse_events table using OpenFDA drug/event API.

    Strategy:
    1. Get approved drugs from compounds (is_drug=true, pref_name IS NOT NULL)
    2. For each, query OpenFDA for top reactions
    3. Store drug_name + reaction + count

    Rate limit: 240 req/min without API key = 4 req/s
    """
    log("=" * 60)
    log("2. OpenFDA Adverse Events")
    log("=" * 60)

    conn: psycopg2.extensions.connection = get_conn()
    existing: int = table_count(conn, "adverse_events")
    log(f"  Existing rows: {existing:,}")

    # Get approved drugs with names
    with conn.cursor() as cur:
        cur.execute("""
            SELECT id, pref_name FROM pharmaco_db.compounds
            WHERE is_drug = true AND pref_name IS NOT NULL
            AND max_phase >= 3
            ORDER BY max_phase DESC, id
        """)
        drugs: List[Tuple[int, str]] = cur.fetchall()
    log(f"  Approved drugs to query: {len(drugs):,}")

    # Check which drugs we already have data for
    already_done: set = set()
    if existing > 0:
        with conn.cursor() as cur:
            cur.execute("SELECT DISTINCT drug_name FROM pharmaco_db.adverse_events WHERE drug_name IS NOT NULL")
            already_done = {row[0].upper() for row in cur.fetchall()}
        log(f"  Already processed: {len(already_done):,} drugs")

    all_rows: List[Tuple] = []
    queried: int = 0
    errors: int = 0
    skipped: int = 0
    api_base: str = "https://api.fda.gov/drug/event.json"

    for i, (cid, drug_name) in enumerate(drugs):
        if drug_name.upper() in already_done:
            skipped += 1
            continue

        # URL-encode the drug name for the search
        encoded_name: str = urllib.parse.quote(drug_name)
        url: str = (
            f"{api_base}?"
            f"search=patient.drug.openfda.generic_name:\"{encoded_name}\""
            f"&count=patient.reaction.reactionmeddrapt.exact"
            f"&limit=25"
        )

        try:
            result: Any = fetch_json(url, timeout=15, retries=2)
            results_list: List[Dict] = result.get("results", [])

            for entry in results_list:
                reaction: str = entry.get("term", "")
                count: int = entry.get("count", 0)
                if reaction and count > 0:
                    all_rows.append((cid, drug_name, reaction, count, "openfda"))

            queried += 1
        except urllib.error.HTTPError as e:
            if e.code == 404:
                # Drug not found in OpenFDA -- expected for many drugs
                pass
            elif e.code == 429:
                # Rate limited -- wait and retry
                log(f"    Rate limited at drug {i}, waiting 60s...")
                time.sleep(60)
                errors += 1
            else:
                errors += 1
        except Exception:
            errors += 1

        # Rate limit: ~3 req/s to stay safe
        time.sleep(0.35)

        # Batch insert every 500 drugs to avoid memory issues
        if (queried + skipped) % 500 == 0 and all_rows:
            with conn.cursor() as cur:
                psycopg2.extras.execute_values(
                    cur,
                    """
                    INSERT INTO pharmaco_db.adverse_events
                        (compound_id, drug_name, reaction, count, source)
                    VALUES %s
                    ON CONFLICT DO NOTHING
                    """,
                    all_rows,
                    page_size=2000,
                )
            conn.commit()
            log(f"    Progress: {i + 1}/{len(drugs)} ({queried} queried, {len(all_rows)} AEs, {errors} errors)")
            all_rows = []

        if errors > 50:
            log(f"  Too many errors ({errors}), stopping early")
            break

    # Final insert
    if all_rows:
        with conn.cursor() as cur:
            psycopg2.extras.execute_values(
                cur,
                """
                INSERT INTO pharmaco_db.adverse_events
                    (compound_id, drug_name, reaction, count, source)
                VALUES %s
                ON CONFLICT DO NOTHING
                """,
                all_rows,
                page_size=2000,
            )
        conn.commit()

    total_inserted: int = table_count(conn, "adverse_events")
    conn.close()
    log(f"  OpenFDA: {total_inserted:,} adverse events from {queried:,} drugs")
    return total_inserted


# ============================================================
# 3. CLINICALTRIALS.GOV
# ============================================================

def fill_clinical_trials() -> int:
    """
    Fill clinical_trials table using ClinicalTrials.gov API v2.

    Strategy:
    1. Get compounds with max_phase >= 1 AND pref_name IS NOT NULL
    2. Query CT.gov API v2 by drug name
    3. Parse and insert trial metadata

    API v2: https://clinicaltrials.gov/api/v2/studies
    Rate: no explicit limit but be respectful (~2 req/s)
    """
    log("=" * 60)
    log("3. ClinicalTrials.gov")
    log("=" * 60)

    conn: psycopg2.extensions.connection = get_conn()
    existing: int = table_count(conn, "clinical_trials")
    log(f"  Existing rows: {existing:,}")

    # Get clinical compounds
    with conn.cursor() as cur:
        cur.execute("""
            SELECT id, pref_name, max_phase FROM pharmaco_db.compounds
            WHERE max_phase >= 1 AND pref_name IS NOT NULL
            ORDER BY max_phase DESC, id
        """)
        compounds: List[Tuple[int, str, int]] = cur.fetchall()
    log(f"  Clinical compounds to query: {len(compounds):,}")

    # Track already-fetched NCT IDs for dedup
    existing_ncts: set = set()
    if existing > 0:
        with conn.cursor() as cur:
            cur.execute("SELECT nct_id FROM pharmaco_db.clinical_trials")
            existing_ncts = {row[0] for row in cur.fetchall()}
        log(f"  Existing NCT IDs: {len(existing_ncts):,}")

    # Also track which compounds we already queried (by name)
    already_queried: set = set()
    if existing > 0:
        with conn.cursor() as cur:
            cur.execute("""
                SELECT DISTINCT c.pref_name FROM pharmaco_db.clinical_trials ct
                JOIN pharmaco_db.compounds c ON c.id = ct.compound_id
                WHERE c.pref_name IS NOT NULL
            """)
            already_queried = {row[0].upper() for row in cur.fetchall()}

    api_base: str = "https://clinicaltrials.gov/api/v2/studies"
    all_rows: List[Tuple] = []
    queried: int = 0
    errors: int = 0
    total_trials: int = 0

    for i, (cid, drug_name, max_phase) in enumerate(compounds):
        if drug_name.upper() in already_queried:
            continue

        # Query CT.gov API v2
        params: Dict[str, str] = {
            "query.intr": drug_name,
            "pageSize": "50",
            "fields": "NCTId,BriefTitle,OverallStatus,Phase,Condition,InterventionName,StartDate,PrimaryCompletionDate,EnrollmentCount",
        }
        url: str = api_base + "?" + urllib.parse.urlencode(params)

        try:
            result: Any = fetch_json(url, timeout=20, retries=2)
            studies: List[Dict] = result.get("studies", [])

            for study in studies:
                proto: Dict = study.get("protocolSection", {})
                ident: Dict = proto.get("identificationModule", {})
                status_mod: Dict = proto.get("statusModule", {})
                design: Dict = proto.get("designModule", {})
                conditions_mod: Dict = proto.get("conditionsModule", {})
                interventions_mod: Dict = proto.get("armsInterventionsModule", {})
                enrollment_info: Dict = design.get("enrollmentInfo", {})

                nct_id: str = ident.get("nctId", "")
                if not nct_id or nct_id in existing_ncts:
                    continue

                title: str = ident.get("briefTitle", "")
                status: str = status_mod.get("overallStatus", "")

                # Phase
                phases: List[str] = design.get("phases", [])
                phase: Optional[str] = phases[0] if phases else None

                # Conditions
                conditions: List[str] = conditions_mod.get("conditions", [])

                # Interventions
                interventions: List[str] = []
                for arm in interventions_mod.get("interventions", []):
                    name: str = arm.get("name", "")
                    if name:
                        interventions.append(name)

                # Dates
                start_date: Optional[str] = None
                sd_struct: Optional[Dict] = status_mod.get("startDateStruct")
                if sd_struct:
                    start_date = sd_struct.get("date")

                completion_date: Optional[str] = None
                cd_struct: Optional[Dict] = status_mod.get("primaryCompletionDateStruct")
                if cd_struct:
                    completion_date = cd_struct.get("date")

                # Enrollment
                enrollment: Optional[int] = None
                if enrollment_info:
                    try:
                        enrollment = int(enrollment_info.get("count", 0))
                    except (ValueError, TypeError):
                        pass

                all_rows.append((
                    cid,
                    nct_id,
                    title[:1000] if title else None,
                    status,
                    phase,
                    conditions if conditions else None,
                    interventions if interventions else None,
                    start_date,
                    completion_date,
                    enrollment,
                    "clinicaltrials.gov",
                ))
                existing_ncts.add(nct_id)
                total_trials += 1

            queried += 1
        except urllib.error.HTTPError as e:
            if e.code == 429:
                log(f"    Rate limited at compound {i}, waiting 30s...")
                time.sleep(30)
            errors += 1
        except Exception:
            errors += 1

        # Rate limit: ~2 req/s
        time.sleep(0.5)

        # Batch insert every 200 compounds
        if queried % 200 == 0 and queried > 0 and all_rows:
            with conn.cursor() as cur:
                psycopg2.extras.execute_values(
                    cur,
                    """
                    INSERT INTO pharmaco_db.clinical_trials
                        (compound_id, nct_id, title, status, phase,
                         conditions, interventions, start_date, completion_date,
                         enrollment, source)
                    VALUES %s
                    ON CONFLICT (nct_id) DO NOTHING
                    """,
                    all_rows,
                    page_size=2000,
                )
            conn.commit()
            log(f"    Progress: {i + 1}/{len(compounds)} ({queried} queried, {total_trials} trials, {errors} errors)")
            all_rows = []

        if errors > 50:
            log(f"  Too many errors ({errors}), stopping early")
            break

    # Final insert
    if all_rows:
        with conn.cursor() as cur:
            psycopg2.extras.execute_values(
                cur,
                """
                INSERT INTO pharmaco_db.clinical_trials
                    (compound_id, nct_id, title, status, phase,
                     conditions, interventions, start_date, completion_date,
                     enrollment, source)
                VALUES %s
                ON CONFLICT (nct_id) DO NOTHING
                """,
                all_rows,
                page_size=2000,
            )
        conn.commit()

    total_inserted: int = table_count(conn, "clinical_trials")
    conn.close()
    log(f"  ClinicalTrials: {total_inserted:,} trials from {queried:,} drugs")
    return total_inserted


# ============================================================
# 4. OPEN TARGETS DISEASE-TARGET ASSOCIATIONS
# ============================================================

def fill_disease_target_associations() -> int:
    """
    Fill disease_target_associations using Open Targets GraphQL API.

    Strategy:
    1. Get human targets with gene_name from our targets table
    2. Map gene_name -> Ensembl ID via Open Targets search
    3. Fetch top disease associations per target
    """
    log("=" * 60)
    log("4. Open Targets Disease-Target Associations")
    log("=" * 60)

    conn: psycopg2.extensions.connection = get_conn()
    existing: int = table_count(conn, "disease_target_associations")
    log(f"  Existing rows: {existing:,}")

    # Get human targets
    with conn.cursor() as cur:
        cur.execute("""
            SELECT gene_name, id FROM pharmaco_db.targets
            WHERE gene_name IS NOT NULL AND organism = 'Homo sapiens'
            ORDER BY id
        """)
        gene_map: Dict[str, int] = dict(cur.fetchall())
    log(f"  Human targets with gene names: {len(gene_map):,}")

    # Track already-done targets
    already_done: set = set()
    if existing > 0:
        with conn.cursor() as cur:
            cur.execute("SELECT DISTINCT target_id FROM pharmaco_db.disease_target_associations")
            already_done = {row[0] for row in cur.fetchall()}
        log(f"  Targets already processed: {len(already_done):,}")

    api_url: str = "https://api.platform.opentargets.org/api/v4/graphql"
    all_rows: List[Tuple] = []
    queried: int = 0
    errors: int = 0
    consecutive_errors: int = 0
    genes: List[Tuple[str, int]] = list(gene_map.items())

    for i, (gene_name, tid) in enumerate(genes):
        if tid in already_done:
            continue

        # Step 1: Search for the Ensembl ID
        search_query: Dict[str, Any] = {
            "query": """
            query($q: String!) {
              search(queryString: $q, entityNames: ["target"], page: {size: 1, index: 0}) {
                hits { id entity }
              }
            }
            """,
            "variables": {"q": gene_name},
        }

        try:
            search_result: Any = fetch_json(api_url, data=search_query, timeout=15, retries=2)
            hits: List[Dict] = search_result.get("data", {}).get("search", {}).get("hits", [])
            if not hits or hits[0].get("entity") != "target":
                queried += 1
                consecutive_errors = 0
                time.sleep(0.2)
                continue

            ensg_id: str = hits[0]["id"]

            # Step 2: Get disease associations
            assoc_query: Dict[str, Any] = {
                "query": """
                query($ensgId: String!) {
                  target(ensemblId: $ensgId) {
                    id
                    approvedSymbol
                    associatedDiseases(page: {size: 50, index: 0}) {
                      rows {
                        disease {
                          id
                          name
                          therapeuticAreas { id name }
                        }
                        score
                        datatypeScores { componentId score }
                      }
                    }
                  }
                }
                """,
                "variables": {"ensgId": ensg_id},
            }

            assoc_result: Any = fetch_json(api_url, data=assoc_query, timeout=15, retries=2)
            target_data: Dict = assoc_result.get("data", {}).get("target", {})
            if not target_data:
                queried += 1
                consecutive_errors = 0
                time.sleep(0.2)
                continue

            rows: List[Dict] = target_data.get("associatedDiseases", {}).get("rows", [])
            for row in rows:
                disease: Dict = row.get("disease", {})
                disease_id: str = disease.get("id", "")
                disease_name: str = disease.get("name", "")
                areas: List[Dict] = disease.get("therapeuticAreas", [])
                therapeutic_area: Optional[str] = areas[0].get("name") if areas else None
                overall_score: Optional[float] = row.get("score")

                # Parse datatype scores
                dt_scores: Dict[str, float] = {
                    d["componentId"]: d["score"]
                    for d in row.get("datatypeScores", [])
                }

                all_rows.append((
                    tid,
                    disease_id,
                    disease_name,
                    therapeutic_area,
                    overall_score,
                    dt_scores.get("ot_genetics_portal"),
                    dt_scores.get("cancer_gene_census") or dt_scores.get("intogen"),
                    dt_scores.get("chembl"),
                    dt_scores.get("europepmc"),
                    dt_scores.get("expression_atlas"),
                    dt_scores.get("phenodigm"),
                    "open_targets",
                ))

            queried += 1
            consecutive_errors = 0

        except Exception as e:
            errors += 1
            consecutive_errors += 1
            if consecutive_errors >= 10:
                log(f"  10 consecutive errors, pausing 30s... Last error: {e}")
                time.sleep(30)
                consecutive_errors = 0
            if errors > 100:
                log(f"  Too many errors ({errors}), stopping")
                break

        # Rate limit: ~3 req/s (2 requests per gene)
        time.sleep(0.3)

        # Batch insert every 200 genes
        if queried % 200 == 0 and queried > 0 and all_rows:
            with conn.cursor() as cur:
                psycopg2.extras.execute_values(
                    cur,
                    """
                    INSERT INTO pharmaco_db.disease_target_associations
                        (target_id, disease_id, disease_name, therapeutic_area,
                         overall_score, genetic_score, somatic_score,
                         known_drug_score, literature_score,
                         rna_expression_score, animal_model_score, source)
                    VALUES %s
                    ON CONFLICT DO NOTHING
                    """,
                    all_rows,
                    page_size=2000,
                )
            conn.commit()
            log(f"    Progress: {i + 1}/{len(genes)} ({queried} queried, {len(all_rows)} associations, {errors} errors)")
            all_rows = []

    # Final insert
    if all_rows:
        with conn.cursor() as cur:
            psycopg2.extras.execute_values(
                cur,
                """
                INSERT INTO pharmaco_db.disease_target_associations
                    (target_id, disease_id, disease_name, therapeutic_area,
                     overall_score, genetic_score, somatic_score,
                     known_drug_score, literature_score,
                     rna_expression_score, animal_model_score, source)
                VALUES %s
                ON CONFLICT DO NOTHING
                """,
                all_rows,
                page_size=2000,
            )
        conn.commit()

    total_inserted: int = table_count(conn, "disease_target_associations")
    conn.close()
    log(f"  Open Targets: {total_inserted:,} disease-target associations")
    return total_inserted


# ============================================================
# 5. DISGENET ASSOCIATIONS
# ============================================================

def fill_disgenet_associations() -> int:
    """
    Fill disgenet_associations from DisGeNET curated GDA TSV.

    DisGeNET provides curated gene-disease associations with scores.
    The curated set requires no API key (open data).

    Sources tried (in order):
    1. DisGeNET curated GDA download (TSV.gz)
    2. If download fails, try the "all" set
    """
    log("=" * 60)
    log("5. DisGeNET Associations")
    log("=" * 60)

    conn: psycopg2.extensions.connection = get_conn()
    existing: int = table_count(conn, "disgenet_associations")
    log(f"  Existing rows: {existing:,}")

    if existing > 0:
        log("  Already populated, skipping")
        conn.close()
        return existing

    # Download curated GDAs
    urls: List[Tuple[str, str]] = [
        (
            "https://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_gene_disease_associations.tsv.gz",
            "DisGeNET curated GDAs",
        ),
        (
            "https://www.disgenet.org/static/disgenet_ap1/files/downloads/all_gene_disease_associations.tsv.gz",
            "DisGeNET all GDAs",
        ),
    ]

    dest: Path = DOWNLOAD_DIR / "disgenet_curated_gda.tsv.gz"
    downloaded: Optional[Path] = None
    for url, desc in urls:
        downloaded = download_file(url, dest, desc)
        if downloaded:
            break

    if not downloaded:
        log("  SKIP: DisGeNET download failed (may need API key)")
        log("  NOTE: DisGeNET now requires free registration for bulk downloads.")
        log("  Trying alternative: mining disease associations from ChEMBL indications...")
        _fill_disgenet_from_chembl(conn)
        total: int = table_count(conn, "disgenet_associations")
        conn.close()
        return total

    # Build gene name -> target_id map
    with conn.cursor() as cur:
        cur.execute("SELECT gene_name, id FROM pharmaco_db.targets WHERE gene_name IS NOT NULL")
        gene_map: Dict[str, int] = dict(cur.fetchall())
    log(f"  Gene map: {len(gene_map):,} targets")

    # Parse TSV
    batch: List[Tuple] = []
    seen: set = set()
    try:
        with gzip.open(downloaded, "rt") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                gene: str = row.get("geneSymbol", "")
                tid: Optional[int] = gene_map.get(gene)
                if not tid:
                    continue

                disease_id: str = row.get("diseaseId", "")
                key: Tuple[str, str] = (gene, disease_id)
                if key in seen:
                    continue
                seen.add(key)

                disease_name: str = row.get("diseaseName", "")
                disease_type: str = row.get("diseaseType", "")
                score_str: str = row.get("score", "")
                ei_str: str = row.get("EI", "")
                n_pmids_str: str = row.get("NofPmids", "")
                n_snps_str: str = row.get("NofSnps", "")
                ncbi_id_str: str = row.get("geneId", "")

                batch.append((
                    tid,
                    gene,
                    int(ncbi_id_str) if ncbi_id_str and ncbi_id_str.isdigit() else None,
                    disease_id,
                    disease_name,
                    disease_type,
                    float(score_str) if score_str else None,
                    float(ei_str) if ei_str else None,
                    int(n_pmids_str) if n_pmids_str and n_pmids_str.isdigit() else None,
                    int(n_snps_str) if n_snps_str and n_snps_str.isdigit() else None,
                    "disgenet",
                ))
    except Exception as e:
        log(f"  DisGeNET parse error: {e}")
        traceback.print_exc()

    if batch:
        with conn.cursor() as cur:
            cur.execute("TRUNCATE pharmaco_db.disgenet_associations")
            psycopg2.extras.execute_values(
                cur,
                """
                INSERT INTO pharmaco_db.disgenet_associations
                    (target_id, gene_name, gene_ncbi_id, disease_id, disease_name,
                     disease_type, association_score, ei_score, num_pmids,
                     num_snps, source)
                VALUES %s
                ON CONFLICT DO NOTHING
                """,
                batch,
                page_size=5000,
            )
        conn.commit()

    total_inserted: int = table_count(conn, "disgenet_associations")
    conn.close()
    log(f"  DisGeNET: {total_inserted:,} gene-disease associations")
    return total_inserted


def _fill_disgenet_from_chembl(conn: psycopg2.extensions.connection) -> None:
    """
    Fallback: mine disease-gene associations from ChEMBL drug_indications.
    Maps target -> drug -> indication -> disease.
    """
    log("  Mining disease-gene associations from ChEMBL drug_indications...")
    with conn.cursor() as cur:
        cur.execute("""
            SELECT DISTINCT
                t.id AS target_id,
                t.gene_name,
                t.ncbi_gene_id,
                di.efo_id AS disease_id,
                di.mesh_heading AS disease_name,
                'disease' AS disease_type,
                CASE WHEN c.max_phase >= 4 THEN 0.9
                     WHEN c.max_phase >= 3 THEN 0.7
                     WHEN c.max_phase >= 2 THEN 0.5
                     ELSE 0.3 END AS score
            FROM pharmaco_db.drug_indications di
            JOIN pharmaco_db.compounds c ON c.id = di.compound_id
            JOIN pharmaco_db.bioactivities ba ON ba.compound_id = c.id
            JOIN pharmaco_db.targets t ON t.id = ba.target_id
            WHERE di.efo_id IS NOT NULL
                AND t.gene_name IS NOT NULL
                AND t.organism = 'Homo sapiens'
                AND ba.pchembl_value >= 6.0
        """)
        rows: List[Tuple] = cur.fetchall()

    if rows:
        batch: List[Tuple] = []
        seen: set = set()
        for row in rows:
            key: Tuple = (row[1], row[3])  # gene_name, disease_id
            if key in seen:
                continue
            seen.add(key)
            batch.append((
                row[0],  # target_id
                row[1],  # gene_name
                row[2],  # ncbi_gene_id
                row[3],  # disease_id
                row[4],  # disease_name
                row[5],  # disease_type
                row[6],  # score
                None,    # ei_score
                None,    # num_pmids
                None,    # num_snps
                "chembl_indications",
            ))

        with conn.cursor() as cur:
            psycopg2.extras.execute_values(
                cur,
                """
                INSERT INTO pharmaco_db.disgenet_associations
                    (target_id, gene_name, gene_ncbi_id, disease_id, disease_name,
                     disease_type, association_score, ei_score, num_pmids,
                     num_snps, source)
                VALUES %s
                ON CONFLICT DO NOTHING
                """,
                batch,
                page_size=5000,
            )
        conn.commit()
        log(f"  ChEMBL indications fallback: {len(batch):,} gene-disease pairs")
    else:
        log("  No drug_indications data available for fallback")


# ============================================================
# 6. TOX21 TOXICOLOGY DATA
# ============================================================

def fill_toxicology_data() -> int:
    """
    Fill toxicology_data table using PubChem PUG REST API for Tox21 assays.

    Strategy:
    1. For each Tox21 assay AID, fetch active+inactive compounds
    2. Match PubChem SIDs/CIDs to our compounds via pubchem_cid
    3. Insert activity outcomes

    Tox21 assay AIDs and their endpoints:
    - 720516: NR-AR (Androgen Receptor)
    - 743053: NR-ER (Estrogen Receptor alpha)
    - 720637: NR-AhR (Aryl hydrocarbon Receptor)
    - 720719: NR-ATAD5 (Genotoxicity)
    - 743122: SR-HSE (Heat Shock Element)
    - 743228: SR-MMP (Mitochondrial Membrane Potential)
    - 720725: SR-p53 (p53 pathway)
    - 743219: SR-ARE (Antioxidant Response Element / Nrf2)
    - 720552: NR-PPARg (PPARgamma)
    - 651741: NR-VDR (Vitamin D Receptor, confirmatory)
    """
    log("=" * 60)
    log("6. Tox21 Toxicology Data (PubChem)")
    log("=" * 60)

    conn: psycopg2.extensions.connection = get_conn()
    existing: int = table_count(conn, "toxicology_data")
    log(f"  Existing rows: {existing:,}")

    if existing > 0:
        log("  Already populated, skipping")
        conn.close()
        return existing

    # Build pubchem_cid -> compound_id map
    with conn.cursor() as cur:
        cur.execute("SELECT pubchem_cid, id FROM pharmaco_db.compounds WHERE pubchem_cid IS NOT NULL")
        pubchem_map: Dict[int, int] = dict(cur.fetchall())
    log(f"  PubChem CID map: {len(pubchem_map):,} compounds")

    if len(pubchem_map) == 0:
        log("  SKIP: No PubChem CIDs mapped. Run mega_consolidate_p2.py first.")
        conn.close()
        return 0

    # Tox21 assay definitions
    tox21_assays: List[Tuple[int, str, str]] = [
        (720516, "Tox21 NR-AR", "NR-AR"),
        (743053, "Tox21 NR-ER", "NR-ER"),
        (720637, "Tox21 NR-AhR", "NR-AhR"),
        (720719, "Tox21 NR-ATAD5", "NR-ATAD5"),
        (743122, "Tox21 SR-HSE", "SR-HSE"),
        (743228, "Tox21 SR-MMP", "SR-MMP"),
        (720725, "Tox21 SR-p53", "SR-p53"),
        (743219, "Tox21 SR-ARE", "SR-ARE"),
        (720552, "Tox21 NR-PPARg", "NR-PPARg"),
        (651741, "Tox21 NR-VDR", "NR-VDR"),
    ]

    all_rows: List[Tuple] = []
    total_matched: int = 0

    for aid, assay_name, endpoint in tox21_assays:
        log(f"  Processing AID {aid} ({endpoint})...")

        # Fetch active CIDs
        active_cids: set = set()
        inactive_cids: set = set()

        for outcome, outcome_label in [("active", "Active"), ("inactive", "Inactive")]:
            cid_set: set = active_cids if outcome == "active" else inactive_cids
            url: str = (
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{aid}"
                f"/cids/JSON?cids_type=standardized&outcome={outcome_label}"
            )
            try:
                result: Any = fetch_json(url, timeout=60, retries=3)
                cid_list: List[int] = (
                    result.get("InformationList", {})
                    .get("Information", [{}])[0]
                    .get("CID", [])
                )
                cid_set.update(cid_list)
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    log(f"    No {outcome} CIDs for AID {aid}")
                else:
                    log(f"    HTTP error for AID {aid} ({outcome}): {e.code}")
            except Exception as e:
                log(f"    Error fetching AID {aid} ({outcome}): {e}")

            time.sleep(0.5)  # PubChem rate limit

        log(f"    AID {aid}: {len(active_cids):,} active, {len(inactive_cids):,} inactive")

        # Match to our compounds
        matched_this_assay: int = 0
        for pcid in active_cids:
            cid: Optional[int] = pubchem_map.get(pcid)
            if cid:
                all_rows.append((cid, assay_name, endpoint, "active", None, None, "tox21"))
                matched_this_assay += 1

        for pcid in inactive_cids:
            cid = pubchem_map.get(pcid)
            if cid:
                all_rows.append((cid, assay_name, endpoint, "inactive", None, None, "tox21"))
                matched_this_assay += 1

        total_matched += matched_this_assay
        log(f"    Matched to our compounds: {matched_this_assay:,}")

        time.sleep(1.0)  # Be nice to PubChem

    # Bulk insert
    if all_rows:
        with conn.cursor() as cur:
            cur.execute("TRUNCATE pharmaco_db.toxicology_data")
            psycopg2.extras.execute_values(
                cur,
                """
                INSERT INTO pharmaco_db.toxicology_data
                    (compound_id, assay_name, assay_endpoint, activity_outcome,
                     ac50_um, efficacy, source)
                VALUES %s
                """,
                all_rows,
                page_size=5000,
            )
        conn.commit()

    total_inserted: int = table_count(conn, "toxicology_data")
    conn.close()
    log(f"  Tox21: {total_inserted:,} toxicology records ({total_matched:,} matched)")
    return total_inserted


# ============================================================
# 7. COMPOUND STRUCTURES (CONTINUATION)
# ============================================================

def fill_compound_structures(max_rows: int = 0) -> int:
    """
    Continue computing compound structures + fingerprints.

    Uses RDKit to compute: morgan_fp, maccs_keys, rdkit_fp,
    scaffold, num_atoms, num_bonds, num_rings, etc.

    Args:
        max_rows: Maximum rows to process (0 = unlimited).
    """
    log("=" * 60)
    log("7. Compound Structures (RDKit continuation)")
    log("=" * 60)

    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, AllChem, MACCSkeys
        from rdkit.Chem.Scaffolds import MurckoScaffold
        import base64
    except ImportError:
        log("  SKIP: RDKit not available (pip install rdkit)")
        return 0

    conn: psycopg2.extensions.connection = get_conn()
    existing: int = table_count(conn, "compound_structures")
    log(f"  Existing structures: {existing:,}")

    # Check actual schema columns
    with conn.cursor() as cur:
        cur.execute("""
            SELECT column_name FROM information_schema.columns
            WHERE table_schema = 'pharmaco_db' AND table_name = 'compound_structures'
            ORDER BY ordinal_position
        """)
        columns: List[str] = [row[0] for row in cur.fetchall()]
    log(f"  Table columns: {', '.join(columns)}")

    # Count remaining
    with conn.cursor() as cur:
        cur.execute("""
            SELECT COUNT(*) FROM pharmaco_db.compounds c
            LEFT JOIN pharmaco_db.compound_structures cs ON cs.compound_id = c.id
            WHERE cs.compound_id IS NULL AND c.canonical_smiles IS NOT NULL
        """)
        remaining: int = cur.fetchone()[0]
    log(f"  Remaining to process: {remaining:,}")

    if remaining == 0:
        log("  All compounds already processed")
        conn.close()
        return existing

    effective_limit: int = min(remaining, max_rows) if max_rows > 0 else remaining
    log(f"  Processing up to: {effective_limit:,}")

    # Use server-side cursor for large reads
    read_conn: psycopg2.extensions.connection = get_conn()
    read_cur = read_conn.cursor("smiles_cursor")
    read_cur.itersize = 5000

    limit_clause: str = f"LIMIT {effective_limit}" if max_rows > 0 else ""
    read_cur.execute(f"""
        SELECT c.id, c.canonical_smiles
        FROM pharmaco_db.compounds c
        LEFT JOIN pharmaco_db.compound_structures cs ON cs.compound_id = c.id
        WHERE cs.compound_id IS NULL AND c.canonical_smiles IS NOT NULL
        ORDER BY c.id
        {limit_clause}
    """)

    batch: List[Tuple] = []
    total: int = 0
    rdkit_errors: int = 0
    start_time: float = time.time()

    for row in read_cur:
        cid: int = row[0]
        smiles: str = row[1]

        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                rdkit_errors += 1
                continue

            # Fingerprints (stored as base64 bitstrings)
            morgan = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            morgan_b64: str = base64.b64encode(morgan.ToBitString().encode()).decode()

            maccs = MACCSkeys.GenMACCSKeys(mol)
            maccs_b64: str = base64.b64encode(maccs.ToBitString().encode()).decode()

            rdkit_fp = Chem.RDKFingerprint(mol)
            rdkit_b64: str = base64.b64encode(rdkit_fp.ToBitString().encode()).decode()

            # Descriptors
            num_atoms: int = mol.GetNumAtoms()
            num_bonds: int = mol.GetNumBonds()
            ring_info = mol.GetRingInfo()
            num_rings: int = ring_info.NumRings()
            num_het: int = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() != 6 and a.GetAtomicNum() != 1)
            frac_csp3: float = Descriptors.FractionCSP3(mol)
            tpsa: float = Descriptors.TPSA(mol)
            mr: float = Descriptors.MolMR(mol)
            n_val: int = Descriptors.NumValenceElectrons(mol)
            n_rad: int = Descriptors.NumRadicalElectrons(mol)
            charge: int = Chem.GetFormalCharge(mol)

            # Murcko scaffold
            scaffold_smi: Optional[str] = None
            try:
                scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                scaffold_smi = Chem.MolToSmiles(scaffold)
            except Exception:
                pass

            batch.append((
                cid,
                smiles,
                morgan_b64,
                maccs_b64,
                rdkit_b64,
                num_atoms,
                num_bonds,
                num_rings,
                num_het,
                round(frac_csp3, 4) if frac_csp3 is not None else None,
                round(tpsa, 2) if tpsa is not None else None,
                round(mr, 2) if mr is not None else None,
                n_val,
                n_rad,
                charge,
                scaffold_smi,
            ))

        except Exception:
            rdkit_errors += 1
            continue

        if len(batch) >= BATCH:
            with conn.cursor() as wr:
                psycopg2.extras.execute_values(
                    wr,
                    """
                    INSERT INTO pharmaco_db.compound_structures
                        (compound_id, canonical_smiles, morgan_fp_2048, maccs_fp, rdkit_fp,
                         num_atoms, num_bonds, num_rings, num_heteroatoms,
                         fraction_csp3, tpsa, molar_refractivity,
                         num_valence_electrons, num_radical_electrons, formal_charge,
                         murcko_scaffold)
                    VALUES %s
                    ON CONFLICT (compound_id) DO NOTHING
                    """,
                    batch,
                    page_size=2000,
                )
            conn.commit()
            total += len(batch)
            batch = []

            if total % 10000 == 0:
                elapsed: float = time.time() - start_time
                rate: float = total / elapsed if elapsed > 0 else 0
                eta: float = (effective_limit - total) / rate if rate > 0 else 0
                log(
                    f"    Progress: {total:,}/{effective_limit:,} "
                    f"({100 * total / effective_limit:.1f}%) "
                    f"rate={rate:.0f}/s "
                    f"errors={rdkit_errors:,} "
                    f"ETA={eta / 60:.1f}min"
                )

    # Final batch
    if batch:
        with conn.cursor() as wr:
            psycopg2.extras.execute_values(
                wr,
                """
                INSERT INTO pharmaco_db.compound_structures
                    (compound_id, canonical_smiles, morgan_fp_2048, maccs_fp, rdkit_fp,
                     num_atoms, num_bonds, num_rings, num_heteroatoms,
                     fraction_csp3, tpsa, molar_refractivity,
                     num_valence_electrons, num_radical_electrons, formal_charge,
                     murcko_scaffold)
                VALUES %s
                ON CONFLICT (compound_id) DO NOTHING
                """,
                batch,
                page_size=2000,
            )
        conn.commit()
        total += len(batch)

    read_cur.close()
    read_conn.close()

    new_total: int = table_count(conn, "compound_structures")
    conn.close()

    elapsed_total: float = time.time() - start_time
    log(
        f"  Compound structures: {total:,} new "
        f"({new_total:,} total, {rdkit_errors:,} errors) "
        f"in {elapsed_total / 60:.1f} min"
    )
    return new_total


# ============================================================
# FINAL STATS
# ============================================================

def print_final_stats() -> None:
    """Print comprehensive database statistics."""
    log("=" * 60)
    log("FINAL DATABASE STATISTICS")
    log("=" * 60)

    conn: psycopg2.extensions.connection = get_conn()
    tables: List[str] = [
        "compounds", "targets", "assays", "bioactivities",
        "drug_mechanisms", "drug_indications", "cross_references",
        "compound_tags", "target_tags", "activity_tags",
        "compound_structures", "atc_classification", "binding_sites",
        "target_interactions", "selectivity_profiles",
        "protein_interactions", "pathways",
        "dgidb_interactions", "gtop_interactions",
        "side_effects", "disgenet_associations",
        "disease_target_associations", "kegg_drug_info",
        "toxicology_data", "pharmacokinetics",
        "clinical_trials", "adverse_events",
    ]

    total_rows: int = 0
    empty_tables: List[str] = []

    with conn.cursor() as cur:
        for tbl in tables:
            try:
                cur.execute(f"SELECT COUNT(*) FROM pharmaco_db.{tbl}")
                cnt: int = cur.fetchone()[0]
                total_rows += cnt
                marker: str = ""
                if cnt == 0:
                    marker = " <<< EMPTY"
                    empty_tables.append(tbl)
                log(f"  {tbl:35s} {cnt:>12,}{marker}")
            except Exception:
                conn.rollback()
                log(f"  {tbl:35s} {'N/A':>12s}")

        # Summary stats
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.compounds WHERE pubchem_cid IS NOT NULL")
        pcid_count: int = cur.fetchone()[0]
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.compounds")
        total_compounds: int = cur.fetchone()[0]
        cur.execute("SELECT pg_size_pretty(pg_database_size('pharmaco'))")
        db_size: str = cur.fetchone()[0]

    log(f"\n  Total rows: {total_rows:,}")
    log(f"  PubChem CID coverage: {pcid_count:,}/{total_compounds:,} ({100 * pcid_count / max(total_compounds, 1):.1f}%)")
    log(f"  Database size: {db_size}")
    log(f"  Empty tables: {len(empty_tables)} ({', '.join(empty_tables) if empty_tables else 'none'})")
    conn.close()


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    start: float = time.time()
    log("=" * 60)
    log("MegaPharmaDB -- Phase 3 Consolidation")
    log(f"  Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    log("=" * 60)

    DOWNLOAD_DIR.mkdir(parents=True, exist_ok=True)

    # Parse arguments
    step: str = sys.argv[1] if len(sys.argv) > 1 else "all"
    max_structures: int = 0  # 0 = unlimited

    # Parse --limit flag for compound_structures
    for i, arg in enumerate(sys.argv):
        if arg == "--limit" and i + 1 < len(sys.argv):
            try:
                max_structures = int(sys.argv[i + 1])
            except ValueError:
                pass

    # Verify DB connection
    try:
        test_conn: psycopg2.extensions.connection = get_conn()
        with test_conn.cursor() as cur:
            cur.execute("SELECT 1 FROM pharmaco_db.compounds LIMIT 1")
        test_conn.close()
        log("  Database connection OK")
    except Exception as e:
        log(f"  ERROR: Cannot connect to database: {e}")
        log("  Make sure PostgreSQL is running on localhost:5433")
        sys.exit(1)

    steps: Dict[str, Any] = {
        "kegg": fill_kegg_drug_info,
        "openfda": fill_adverse_events,
        "clinical": fill_clinical_trials,
        "opentargets": fill_disease_target_associations,
        "disgenet": fill_disgenet_associations,
        "tox21": fill_toxicology_data,
        "structures": lambda: fill_compound_structures(max_rows=max_structures),
    }

    if step == "all":
        for name, func in steps.items():
            log(f"\n{'=' * 60}")
            try:
                func()
            except Exception as e:
                log(f"  ERROR in {name}: {e}")
                traceback.print_exc()
                log(f"  Continuing to next step...")
    elif step in steps:
        try:
            steps[step]()
        except Exception as e:
            log(f"  ERROR: {e}")
            traceback.print_exc()
    elif step == "stats":
        pass  # Just print stats
    else:
        log(f"  Unknown step: {step}")
        log(f"  Available: {', '.join(steps.keys())}, all, stats")
        sys.exit(1)

    print_final_stats()

    elapsed: float = time.time() - start
    log(f"\nTotal time: {elapsed / 60:.1f} minutes")


if __name__ == "__main__":
    main()
