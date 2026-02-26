#!/usr/bin/env python3
"""
Quick fixes for the remaining empty tables that mega_consolidate_p3.py struggled with.
Runs: OpenFDA (fixed URL), ClinicalTrials (fixed API), DisGeNET (ChEMBL fallback),
      Open Targets (GraphQL), Tox21 (fixed PubChem API)
"""
import json
import time
import urllib.request
import urllib.parse
import urllib.error
import traceback
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

import psycopg2
import psycopg2.extras

DB = {"host": "localhost", "port": 5433, "dbname": "pharmaco", "user": "postgres", "password": "postgres"}

def log(msg: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

def get_conn():
    return psycopg2.connect(**DB)

def fetch_json(url: str, retries: int = 3, delay: float = 2.0) -> Any:
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "PharmacoDB/1.0"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                return json.loads(resp.read().decode())
        except urllib.error.HTTPError as e:
            if e.code == 429:
                time.sleep(30)
                continue
            if attempt < retries - 1:
                time.sleep(delay)
                continue
            raise
        except Exception:
            if attempt < retries - 1:
                time.sleep(delay)
                continue
            raise

# ============================================================
# 1. OpenFDA - Fixed URL encoding
# ============================================================
def fill_openfda():
    log("=" * 60)
    log("OpenFDA Adverse Events")
    log("=" * 60)
    conn = get_conn()

    with conn.cursor() as cur:
        cur.execute("""
            SELECT id, pref_name FROM pharmaco_db.compounds
            WHERE is_drug = true AND pref_name IS NOT NULL AND max_phase >= 4
            ORDER BY id
        """)
        drugs = cur.fetchall()
    log(f"  Approved drugs: {len(drugs):,}")

    # Check already done
    with conn.cursor() as cur:
        cur.execute("SELECT DISTINCT compound_id FROM pharmaco_db.adverse_events")
        done = {r[0] for r in cur.fetchall()}
    log(f"  Already done: {len(done):,}")

    batch = []
    queried = 0
    found = 0
    errors = 0

    for i, (cid, name) in enumerate(drugs):
        if cid in done:
            continue

        # OpenFDA needs exact lower-case match with quotes
        safe_name = name.replace('"', '\\"').lower()
        url = (
            f"https://api.fda.gov/drug/event.json?"
            f"search=patient.drug.openfda.generic_name.exact:%22{urllib.parse.quote(safe_name)}%22"
            f"&count=patient.reaction.reactionmeddrapt.exact&limit=25"
        )

        try:
            data = fetch_json(url, retries=1, delay=1)
            results = data.get("results", [])
            for entry in results:
                reaction = entry.get("term", "")
                count = entry.get("count", 0)
                if reaction and count > 0:
                    batch.append((cid, name, reaction, count, "openfda"))
            if results:
                found += 1
            queried += 1
        except urllib.error.HTTPError as e:
            if e.code == 404:
                queried += 1  # Drug not found is normal
            else:
                errors += 1
        except Exception:
            errors += 1

        time.sleep(0.3)

        if (i + 1) % 200 == 0:
            if batch:
                with conn.cursor() as cur:
                    psycopg2.extras.execute_values(cur,
                        """INSERT INTO pharmaco_db.adverse_events
                           (compound_id, drug_name, reaction, count, source)
                           VALUES %s ON CONFLICT DO NOTHING""",
                        batch, page_size=2000)
                conn.commit()
                batch = []
            log(f"  Progress: {i+1}/{len(drugs)} queried={queried} found={found} AEs_total={table_count(conn)} errors={errors}")

        if errors > 100:
            log("  Too many errors, stopping")
            break

    if batch:
        with conn.cursor() as cur:
            psycopg2.extras.execute_values(cur,
                """INSERT INTO pharmaco_db.adverse_events
                   (compound_id, drug_name, reaction, count, source)
                   VALUES %s ON CONFLICT DO NOTHING""",
                batch, page_size=2000)
        conn.commit()

    total = table_count(conn)
    conn.close()
    log(f"  Done: {total:,} adverse event rows")
    return total

def table_count(conn, table="adverse_events"):
    with conn.cursor() as cur:
        cur.execute(f"SELECT COUNT(*) FROM pharmaco_db.{table}")
        return cur.fetchone()[0]

# ============================================================
# 2. ClinicalTrials.gov - API v2
# ============================================================
def fill_clinical_trials():
    log("=" * 60)
    log("ClinicalTrials.gov")
    log("=" * 60)
    conn = get_conn()

    with conn.cursor() as cur:
        cur.execute("""
            SELECT id, pref_name FROM pharmaco_db.compounds
            WHERE max_phase >= 2 AND pref_name IS NOT NULL
            ORDER BY max_phase DESC, id
        """)
        drugs = cur.fetchall()
    log(f"  Drugs to query: {len(drugs):,}")

    # Already done
    with conn.cursor() as cur:
        cur.execute("SELECT DISTINCT compound_id FROM pharmaco_db.clinical_trials")
        done = {r[0] for r in cur.fetchall()}
    log(f"  Already done: {len(done):,}")

    batch = []
    nct_seen = set()
    queried = 0
    found = 0
    errors = 0

    for i, (cid, name) in enumerate(drugs):
        if cid in done:
            continue

        # ClinicalTrials.gov API v2
        params = urllib.parse.urlencode({
            "query.intr": name,
            "pageSize": 20,
            "format": "json",
            "fields": "NCTId,BriefTitle,OverallStatus,Phase,Condition,InterventionName,StartDate,CompletionDate,EnrollmentCount",
        })
        url = f"https://clinicaltrials.gov/api/v2/studies?{params}"

        try:
            data = fetch_json(url, retries=2, delay=2)
            studies = data.get("studies", [])
            for study in studies:
                proto = study.get("protocolSection", {})
                ident = proto.get("identificationModule", {})
                status_mod = proto.get("statusModule", {})
                design = proto.get("designModule", {})
                conds = proto.get("conditionsModule", {})
                intv = proto.get("armsInterventionsModule", {})
                enroll = design.get("enrollmentInfo", {})

                nct_id = ident.get("nctId", "")
                if not nct_id or nct_id in nct_seen:
                    continue
                nct_seen.add(nct_id)

                conditions = conds.get("conditions", [])
                interventions_list = []
                for arm in intv.get("interventions", []):
                    interventions_list.append(arm.get("name", ""))

                batch.append((
                    cid, nct_id,
                    ident.get("briefTitle", "")[:500],
                    status_mod.get("overallStatus", ""),
                    ",".join(design.get("phases", [])),
                    conditions if conditions else None,
                    interventions_list if interventions_list else None,
                    status_mod.get("startDateStruct", {}).get("date", ""),
                    status_mod.get("completionDateStruct", {}).get("date", ""),
                    enroll.get("count"),
                    "clinicaltrials.gov",
                ))

            if studies:
                found += 1
            queried += 1
        except urllib.error.HTTPError as e:
            if e.code == 404:
                queried += 1
            else:
                errors += 1
        except Exception:
            errors += 1

        time.sleep(0.5)

        if (i + 1) % 200 == 0:
            if batch:
                with conn.cursor() as cur:
                    psycopg2.extras.execute_values(cur,
                        """INSERT INTO pharmaco_db.clinical_trials
                           (compound_id, nct_id, title, status, phase, conditions,
                            interventions, start_date, completion_date, enrollment, source)
                           VALUES %s ON CONFLICT (nct_id) DO NOTHING""",
                        batch, page_size=2000)
                conn.commit()
                batch = []
            with conn.cursor() as cur:
                cur.execute("SELECT COUNT(*) FROM pharmaco_db.clinical_trials")
                ct_count = cur.fetchone()[0]
            log(f"  Progress: {i+1}/{len(drugs)} queried={queried} found={found} trials={ct_count} errors={errors}")

        if errors > 100:
            log("  Too many errors, stopping")
            break

    if batch:
        with conn.cursor() as cur:
            psycopg2.extras.execute_values(cur,
                """INSERT INTO pharmaco_db.clinical_trials
                   (compound_id, nct_id, title, status, phase, conditions,
                    interventions, start_date, completion_date, enrollment, source)
                   VALUES %s ON CONFLICT (nct_id) DO NOTHING""",
                batch, page_size=2000)
        conn.commit()

    with conn.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.clinical_trials")
        total = cur.fetchone()[0]
    conn.close()
    log(f"  Done: {total:,} clinical trials")
    return total

# ============================================================
# 3. DisGeNET fallback from ChEMBL
# ============================================================
def fill_disgenet_from_chembl():
    log("=" * 60)
    log("DisGeNET (ChEMBL indications fallback)")
    log("=" * 60)
    conn = get_conn()

    with conn.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.disgenet_associations")
        if cur.fetchone()[0] > 0:
            log("  Already populated, skipping")
            conn.close()
            return

    log("  Mining disease-gene associations from drug_indications + bioactivities...")
    with conn.cursor() as cur:
        cur.execute("""
            INSERT INTO pharmaco_db.disgenet_associations
                (target_id, gene_name, gene_ncbi_id, disease_id, disease_name,
                 disease_type, association_score, source)
            SELECT DISTINCT ON (t.gene_name, di.efo_id)
                t.id,
                t.gene_name,
                t.ncbi_gene_id::integer,
                COALESCE(di.efo_id, di.mesh_id),
                COALESCE(di.efo_term, di.mesh_heading),
                'disease',
                CASE WHEN c.max_phase >= 4 THEN 0.9
                     WHEN c.max_phase >= 3 THEN 0.7
                     WHEN c.max_phase >= 2 THEN 0.5
                     ELSE 0.3 END,
                'chembl_indications'
            FROM pharmaco_db.drug_indications di
            JOIN pharmaco_db.compounds c ON c.id = di.compound_id
            JOIN pharmaco_db.drug_mechanisms dm ON dm.compound_id = c.id
            JOIN pharmaco_db.targets t ON t.id = dm.target_id
            WHERE (di.efo_id IS NOT NULL OR di.mesh_id IS NOT NULL)
                AND t.gene_name IS NOT NULL
                AND t.organism = 'Homo sapiens'
            ON CONFLICT (gene_name, disease_id) DO NOTHING
        """)
        conn.commit()

    # Also add from bioactivities (high confidence active compounds)
    with conn.cursor() as cur:
        cur.execute("""
            INSERT INTO pharmaco_db.disgenet_associations
                (target_id, gene_name, gene_ncbi_id, disease_id, disease_name,
                 disease_type, association_score, source)
            SELECT DISTINCT ON (t.gene_name, di.efo_id)
                t.id,
                t.gene_name,
                t.ncbi_gene_id::integer,
                COALESCE(di.efo_id, di.mesh_id),
                COALESCE(di.efo_term, di.mesh_heading),
                'disease',
                CASE WHEN ba.pchembl_value >= 8 THEN 0.6
                     WHEN ba.pchembl_value >= 7 THEN 0.4
                     ELSE 0.2 END,
                'chembl_bioactivity'
            FROM pharmaco_db.drug_indications di
            JOIN pharmaco_db.compounds c ON c.id = di.compound_id
            JOIN pharmaco_db.bioactivities ba ON ba.compound_id = c.id AND ba.pchembl_value >= 6.0
            JOIN pharmaco_db.targets t ON t.id = ba.target_id
            WHERE (di.efo_id IS NOT NULL OR di.mesh_id IS NOT NULL)
                AND t.gene_name IS NOT NULL
                AND t.organism = 'Homo sapiens'
            ON CONFLICT (gene_name, disease_id) DO NOTHING
        """)
        conn.commit()

    with conn.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.disgenet_associations")
        total = cur.fetchone()[0]
    conn.close()
    log(f"  Done: {total:,} disease-gene associations from ChEMBL")
    return total

# ============================================================
# 4. Open Targets via GraphQL
# ============================================================
def fill_open_targets():
    log("=" * 60)
    log("Open Targets Disease Associations")
    log("=" * 60)
    conn = get_conn()

    with conn.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.disease_target_associations")
        if cur.fetchone()[0] > 0:
            log("  Already populated")
            conn.close()
            return

    # Get gene names
    with conn.cursor() as cur:
        cur.execute("""
            SELECT id, gene_name FROM pharmaco_db.targets
            WHERE gene_name IS NOT NULL AND organism = 'Homo sapiens'
            ORDER BY id
        """)
        genes = cur.fetchall()
    log(f"  Targets to query: {len(genes):,}")

    api_url = "https://api.platform.opentargets.org/api/v4/graphql"
    batch = []
    queried = 0
    found = 0
    errors = 0
    consecutive_errors = 0

    for i, (tid, gene) in enumerate(genes):
        # Step 1: Get Ensembl ID from gene symbol
        search_query = json.dumps({
            "query": """{
                search(queryString: "%s", entityNames: ["target"], page: {index: 0, size: 1}) {
                    hits { id name }
                }
            }""" % gene.replace('"', '\\"'),
        }).encode()

        ensembl_id = None
        try:
            req = urllib.request.Request(api_url, data=search_query,
                headers={"Content-Type": "application/json", "User-Agent": "PharmacoDB/1.0"})
            with urllib.request.urlopen(req, timeout=15) as resp:
                data = json.loads(resp.read().decode())
            hits = data.get("data", {}).get("search", {}).get("hits", [])
            if hits:
                ensembl_id = hits[0].get("id")
        except Exception:
            errors += 1
            consecutive_errors += 1
            if consecutive_errors > 10:
                log(f"  10 consecutive errors, pausing 30s...")
                time.sleep(30)
                consecutive_errors = 0
            time.sleep(0.5)
            continue

        if not ensembl_id:
            queried += 1
            time.sleep(0.3)
            continue

        # Step 2: Get disease associations
        assoc_query = json.dumps({
            "query": """{
                target(ensemblId: "%s") {
                    associatedDiseases(page: {index: 0, size: 50}) {
                        rows {
                            disease { id name therapeuticAreas { id name } }
                            score
                            datasourceScores { componentId score }
                        }
                    }
                }
            }""" % ensembl_id,
        }).encode()

        try:
            req = urllib.request.Request(api_url, data=assoc_query,
                headers={"Content-Type": "application/json", "User-Agent": "PharmacoDB/1.0"})
            with urllib.request.urlopen(req, timeout=15) as resp:
                data = json.loads(resp.read().decode())

            target_data = data.get("data", {}).get("target", {})
            if not target_data:
                queried += 1
                time.sleep(0.3)
                continue

            rows = target_data.get("associatedDiseases", {}).get("rows", [])
            for row in rows:
                disease = row.get("disease", {})
                disease_id = disease.get("id", "")
                disease_name = disease.get("name", "")

                ta_list = disease.get("therapeuticAreas", [])
                ta = ta_list[0].get("name", "") if ta_list else ""

                overall = row.get("score", 0)

                # Parse datasource scores
                scores = {s.get("componentId", ""): s.get("score", 0)
                         for s in row.get("datasourceScores", [])}

                batch.append((
                    tid, disease_id, disease_name, ta, overall,
                    scores.get("ot_genetics_portal", None),
                    scores.get("cancer_gene_census", None),
                    scores.get("chembl", None),
                    scores.get("europepmc", None),
                    scores.get("expression_atlas", None),
                    scores.get("impc", None),
                    "open_targets",
                ))

            if rows:
                found += 1
            queried += 1
            consecutive_errors = 0
        except Exception:
            errors += 1
            consecutive_errors += 1

        time.sleep(0.35)

        if (i + 1) % 200 == 0:
            if batch:
                with conn.cursor() as cur:
                    psycopg2.extras.execute_values(cur,
                        """INSERT INTO pharmaco_db.disease_target_associations
                           (target_id, disease_id, disease_name, therapeutic_area,
                            overall_score, genetic_score, somatic_score,
                            known_drug_score, literature_score, rna_expression_score,
                            animal_model_score, source)
                           VALUES %s ON CONFLICT DO NOTHING""",
                        batch, page_size=2000)
                conn.commit()
                batch = []
            with conn.cursor() as cur:
                cur.execute("SELECT COUNT(*) FROM pharmaco_db.disease_target_associations")
                ot_count = cur.fetchone()[0]
            log(f"  Progress: {i+1}/{len(genes)} queried={queried} found={found} assocs={ot_count} errors={errors}")

        if errors > 200:
            log("  Too many errors, stopping")
            break

    if batch:
        with conn.cursor() as cur:
            psycopg2.extras.execute_values(cur,
                """INSERT INTO pharmaco_db.disease_target_associations
                   (target_id, disease_id, disease_name, therapeutic_area,
                    overall_score, genetic_score, somatic_score,
                    known_drug_score, literature_score, rna_expression_score,
                    animal_model_score, source)
                   VALUES %s ON CONFLICT DO NOTHING""",
                batch, page_size=2000)
        conn.commit()

    with conn.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.disease_target_associations")
        total = cur.fetchone()[0]
    conn.close()
    log(f"  Done: {total:,} disease-target associations")
    return total

# ============================================================
# 5. Tox21 from PubChem - Using BioAssay Concise CSV
# ============================================================
def fill_tox21():
    log("=" * 60)
    log("Tox21 Toxicology Data")
    log("=" * 60)
    conn = get_conn()

    with conn.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.toxicology_data")
        if cur.fetchone()[0] > 0:
            log("  Already populated")
            conn.close()
            return

    # Build pubchem_cid -> compound_id map
    with conn.cursor() as cur:
        cur.execute("SELECT pubchem_cid::bigint, id FROM pharmaco_db.compounds WHERE pubchem_cid IS NOT NULL")
        cid_map = {int(r[0]): r[1] for r in cur.fetchall()}
    log(f"  PubChem CID map: {len(cid_map):,} compounds")

    # Tox21 assays - use the summary endpoint
    assays = [
        (720516, "NR-AR", "Androgen Receptor"),
        (743053, "NR-ER", "Estrogen Receptor"),
        (720637, "NR-AhR", "Aryl Hydrocarbon Receptor"),
        (720719, "SR-ATAD5", "DNA Damage (ATAD5)"),
        (743122, "SR-HSE", "Heat Shock Response"),
        (743228, "SR-MMP", "Mitochondrial Membrane Potential"),
        (720725, "SR-p53", "p53 Tumor Suppressor"),
        (743219, "SR-ARE", "Antioxidant Response (Nrf2)"),
        (720552, "NR-PPARg", "PPARgamma"),
        (651741, "NR-VDR", "Vitamin D Receptor"),
    ]

    batch = []
    total_matched = 0

    for aid, assay_name, endpoint in assays:
        log(f"  Processing AID {aid} ({assay_name})...")

        for outcome in ["active", "inactive"]:
            # Use the concise data download
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{aid}/cids/JSON?cids_type={outcome}"
            try:
                data = fetch_json(url, retries=3, delay=3)
                info = data.get("InformationList", {}).get("Information", [{}])
                cids = []
                for entry in info:
                    cids.extend(entry.get("CID", []))

                matched = 0
                for cid in cids:
                    compound_id = cid_map.get(cid)
                    if compound_id:
                        batch.append((
                            compound_id, assay_name, endpoint,
                            outcome, None, None, "tox21"
                        ))
                        matched += 1

                log(f"    {outcome}: {len(cids):,} CIDs, {matched:,} matched")
                total_matched += matched
            except Exception as e:
                log(f"    Error fetching {outcome} for AID {aid}: {e}")

            time.sleep(1)

    if batch:
        with conn.cursor() as cur:
            psycopg2.extras.execute_values(cur,
                """INSERT INTO pharmaco_db.toxicology_data
                   (compound_id, assay_name, assay_endpoint, activity_outcome,
                    ac50_um, efficacy, source)
                   VALUES %s ON CONFLICT DO NOTHING""",
                batch, page_size=5000)
        conn.commit()

    with conn.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.toxicology_data")
        total = cur.fetchone()[0]
    conn.close()
    log(f"  Done: {total:,} toxicology records ({total_matched:,} matched)")
    return total


# ============================================================
# MAIN
# ============================================================
import sys

def main():
    log("=" * 60)
    log("PharmacoDB — Fill Remaining Tables")
    log("=" * 60)

    steps = {
        "disgenet": fill_disgenet_from_chembl,
        "openfda": fill_openfda,
        "clinical": fill_clinical_trials,
        "opentargets": fill_open_targets,
        "tox21": fill_tox21,
    }

    step = sys.argv[1] if len(sys.argv) > 1 else "all"

    if step == "all":
        for name, func in steps.items():
            try:
                func()
            except Exception as e:
                log(f"  ERROR in {name}: {e}")
                traceback.print_exc()
    elif step in steps:
        try:
            steps[step]()
        except Exception as e:
            log(f"  ERROR: {e}")
            traceback.print_exc()
    else:
        log(f"  Unknown: {step}. Available: {', '.join(steps.keys())}, all")

    # Print final stats
    conn = get_conn()
    for tbl in ["adverse_events", "clinical_trials", "disgenet_associations",
                "disease_target_associations", "toxicology_data"]:
        with conn.cursor() as cur:
            cur.execute(f"SELECT COUNT(*) FROM pharmaco_db.{tbl}")
            log(f"  {tbl}: {cur.fetchone()[0]:,}")
    conn.close()

if __name__ == "__main__":
    main()
