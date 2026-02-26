#!/usr/bin/env python3
"""
PharmacoDB — Additional Data Sources Ingestion
1. Guide to Pharmacology (IUPHAR/BPS) — Curated pharmacology
2. Open Targets — Target-disease associations
3. ChEMBL ATC Classification — Therapeutic classification
4. BindingDB — Additional binding data
5. ChEMBL Target Components — Protein domain info
"""

import json
import time
import csv
import io
import urllib.request
import urllib.error
import psycopg2
import psycopg2.extras
from datetime import datetime

DB_CONFIG = {
    "host": "localhost",
    "port": 5433,
    "dbname": "pharmaco",
    "user": "postgres",
    "password": "pharmaco_secret"
}

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"

def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

def api_get(url, retries=3, timeout=30):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={
                "Accept": "application/json",
                "User-Agent": "PharmacoDB-Ingestion/1.0"
            })
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode())
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(3 * (attempt + 1))
            else:
                return None
    return None

def fetch_text(url, retries=3, timeout=60):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={
                "User-Agent": "PharmacoDB-Ingestion/1.0"
            })
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read().decode('utf-8', errors='replace')
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(5 * (attempt + 1))
            else:
                log(f"    Failed to fetch {url}: {e}")
                return None
    return None

def get_conn():
    return psycopg2.connect(**DB_CONFIG)

def start_log(conn, source, step):
    with conn.cursor() as cur:
        cur.execute("""
            INSERT INTO pharmaco_db.ingestion_log (source, step, status)
            VALUES (%s, %s, 'running') RETURNING id
        """, (source, step))
        log_id = cur.fetchone()[0]
    conn.commit()
    return log_id

def update_log(conn, log_id, status, rows=None, error=None):
    with conn.cursor() as cur:
        cur.execute("""
            UPDATE pharmaco_db.ingestion_log
            SET status=%s, rows_inserted=%s, error_message=%s, completed_at=NOW()
            WHERE id=%s
        """, (status, rows, error, log_id))
    conn.commit()

# ============================================================
# 1. GUIDE TO PHARMACOLOGY (IUPHAR/BPS)
# ============================================================
def ingest_guide_to_pharmacology(conn):
    log("=" * 60)
    log("Ingesting Guide to Pharmacology (IUPHAR/BPS)")
    log("=" * 60)
    log_id = start_log(conn, 'guide_to_pharmacology', 'interactions')

    # GtoPdb has a REST API
    base = "https://www.guidetopharmacology.org/services"

    # Build lookup maps
    compound_map = {}
    with conn.cursor() as cur:
        cur.execute("SELECT pref_name, id FROM pharmaco_db.compounds WHERE pref_name IS NOT NULL")
        for name, cid in cur.fetchall():
            compound_map[name.upper()] = cid

    target_map = {}
    with conn.cursor() as cur:
        cur.execute("SELECT gene_name, id FROM pharmaco_db.targets WHERE gene_name IS NOT NULL")
        for gene, tid in cur.fetchall():
            target_map[gene.upper()] = tid
        cur.execute("SELECT pref_name, id FROM pharmaco_db.targets")
        for name, tid in cur.fetchall():
            if name:
                target_map[name.upper()] = tid

    log(f"  Compound map: {len(compound_map)}, Target map: {len(target_map)}")

    # Fetch ligands
    log("  Fetching GtoPdb ligands...")
    ligands_data = api_get(f"{base}/ligands?type=Synthetic+organic", timeout=120)
    if not ligands_data:
        log("  WARN: Could not fetch GtoPdb ligands")
        update_log(conn, log_id, 'failed', 0, "Could not fetch ligands")
        return 0

    log(f"  Got {len(ligands_data)} ligands")
    total = 0

    for i, lig in enumerate(ligands_data[:5000]):  # Process top 5000
        if i % 500 == 0 and i > 0:
            log(f"    Progress: {i}/{min(len(ligands_data), 5000)} ({total} interactions)")

        lig_id = lig.get("ligandId")
        lig_name = (lig.get("name") or "").upper()

        # Match to our compound
        compound_id = compound_map.get(lig_name)

        # Fetch interactions for this ligand
        interactions = api_get(f"{base}/ligands/{lig_id}/interactions")
        if not interactions:
            time.sleep(0.15)
            continue

        batch = []
        for inter in interactions:
            target_name = (inter.get("targetName") or "").upper()
            target_id = target_map.get(target_name)

            affinity_str = inter.get("affinity")
            affinity_value = None
            if affinity_str:
                try:
                    affinity_value = float(affinity_str)
                except (ValueError, TypeError):
                    pass

            batch.append((
                compound_id,
                target_id,
                lig_id,
                inter.get("targetId"),
                inter.get("type"),
                inter.get("affinityParameter"),
                affinity_value,
                inter.get("affinityHigh"),
                inter.get("affinityLow"),
                inter.get("species"),
                inter.get("endogenous", False),
                inter.get("primaryTarget", False),
            ))

        if batch:
            with conn.cursor() as cur:
                psycopg2.extras.execute_values(cur, """
                    INSERT INTO pharmaco_db.gtop_interactions
                    (compound_id, target_id, ligand_id, target_gtop_id,
                     interaction_type, affinity_type, affinity_value,
                     affinity_high, affinity_low, species, endogenous, primary_target)
                    VALUES %s
                """, batch, page_size=200)
            conn.commit()
            total += len(batch)

        time.sleep(0.2)

    log(f"  DONE: {total} GtoPdb interactions inserted")
    update_log(conn, log_id, 'completed', total)
    return total

# ============================================================
# 2. OPEN TARGETS — Target-Disease Associations
# ============================================================
def ingest_open_targets(conn):
    log("=" * 60)
    log("Ingesting Open Targets disease-target associations")
    log("=" * 60)
    log_id = start_log(conn, 'open_targets', 'disease_associations')

    # Open Targets Platform API
    OT_API = "https://api.platform.opentargets.org/api/v4/graphql"

    # Get all our targets with Ensembl Gene IDs (via UniProt → gene name)
    target_genes = {}
    with conn.cursor() as cur:
        cur.execute("""
            SELECT id, gene_name FROM pharmaco_db.targets
            WHERE gene_name IS NOT NULL
            ORDER BY id
        """)
        for tid, gene in cur.fetchall():
            target_genes[gene] = tid

    log(f"  Found {len(target_genes)} targets with gene names")
    total = 0
    processed = 0

    for gene_name, target_id in target_genes.items():
        processed += 1
        if processed % 200 == 0:
            log(f"    Progress: {processed}/{len(target_genes)} ({total} associations)")

        # GraphQL query for target associations
        query = {
            "query": """
            query TargetAssociations($ensgId: String!) {
                target(ensemblId: $ensgId) {
                    associatedDiseases(page: {size: 25, index: 0}) {
                        rows {
                            disease { id name therapeuticAreas { id name } }
                            score
                            datatypeScores {
                                componentId: id
                                score
                            }
                        }
                    }
                }
            }
            """,
            "variables": {"ensgId": gene_name}  # This won't work with gene symbol
        }

        # Use the search endpoint instead
        search_url = f"https://api.platform.opentargets.org/api/v4/graphql"
        search_query = {
            "query": """
            query SearchTarget($queryString: String!) {
                search(queryString: $queryString, entityNames: ["target"], page: {size: 1, index: 0}) {
                    hits { id }
                }
            }
            """,
            "variables": {"queryString": gene_name}
        }

        try:
            req = urllib.request.Request(search_url,
                data=json.dumps(search_query).encode(),
                headers={"Content-Type": "application/json", "User-Agent": "PharmacoDB/1.0"})
            with urllib.request.urlopen(req, timeout=15) as resp:
                search_data = json.loads(resp.read().decode())

            hits = search_data.get("data", {}).get("search", {}).get("hits", [])
            if not hits:
                time.sleep(0.15)
                continue

            ensembl_id = hits[0].get("id")
            if not ensembl_id:
                continue

            # Now fetch associations
            assoc_query = {
                "query": """
                query TargetAssoc($ensemblId: String!) {
                    target(ensemblId: $ensemblId) {
                        associatedDiseases(page: {size: 25, index: 0}) {
                            rows {
                                disease { id name therapeuticAreas { id name } }
                                score
                                datatypeScores { id score }
                            }
                        }
                    }
                }
                """,
                "variables": {"ensemblId": ensembl_id}
            }

            req2 = urllib.request.Request(search_url,
                data=json.dumps(assoc_query).encode(),
                headers={"Content-Type": "application/json", "User-Agent": "PharmacoDB/1.0"})
            with urllib.request.urlopen(req2, timeout=15) as resp2:
                assoc_data = json.loads(resp2.read().decode())

            rows = (assoc_data.get("data", {}).get("target", {})
                    .get("associatedDiseases", {}).get("rows", []))

            batch = []
            for row in rows:
                disease = row.get("disease", {})
                therapeutic_areas = disease.get("therapeuticAreas", [])
                ta_name = therapeutic_areas[0].get("name") if therapeutic_areas else None

                # Extract datatype scores
                dt_scores = {ds.get("id"): ds.get("score") for ds in row.get("datatypeScores", [])}

                batch.append((
                    target_id,
                    disease.get("id"),
                    disease.get("name"),
                    ta_name,
                    row.get("score"),
                    dt_scores.get("ot_genetics_portal"),
                    dt_scores.get("cancer_gene_census"),
                    dt_scores.get("chembl"),
                    dt_scores.get("europepmc"),
                    dt_scores.get("expression_atlas"),
                    dt_scores.get("impc"),
                ))

            if batch:
                with conn.cursor() as cur:
                    psycopg2.extras.execute_values(cur, """
                        INSERT INTO pharmaco_db.disease_target_associations
                        (target_id, disease_id, disease_name, therapeutic_area,
                         overall_score, genetic_score, somatic_score, known_drug_score,
                         literature_score, rna_expression_score, animal_model_score)
                        VALUES %s
                    """, batch, page_size=200)
                conn.commit()
                total += len(batch)

        except Exception as e:
            pass  # Skip silently

        time.sleep(0.2)  # Rate limit

    log(f"  DONE: {total} disease-target associations inserted")
    update_log(conn, log_id, 'completed', total)
    return total

# ============================================================
# 3. CHEMBL ATC CLASSIFICATION
# ============================================================
def ingest_atc_classification(conn):
    log("=" * 60)
    log("Ingesting ChEMBL ATC Classification")
    log("=" * 60)
    log_id = start_log(conn, 'chembl', 'atc_classification')

    compound_map = {}
    with conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.compounds")
        compound_map = dict(cur.fetchall())

    total = 0
    offset = 0

    while True:
        url = f"{CHEMBL_API}/atc_class.json?limit=500&offset={offset}"
        data = api_get(url)
        if not data or not data.get("atc"):
            break

        atc_list = data["atc"]
        log(f"  Fetched {len(atc_list)} ATC records (offset={offset})")

        # For each ATC, find associated molecules
        for atc in atc_list:
            atc_code = atc.get("level5")
            if not atc_code:
                continue

            # Fetch molecules with this ATC code
            mol_url = f"{CHEMBL_API}/molecule.json?atc_classifications={atc_code}&limit=100"
            mol_data = api_get(mol_url)
            if not mol_data or not mol_data.get("molecules"):
                time.sleep(0.2)
                continue

            batch = []
            for mol in mol_data["molecules"]:
                mol_chembl = mol.get("molecule_chembl_id")
                compound_id = compound_map.get(mol_chembl)
                if not compound_id:
                    continue

                batch.append((
                    compound_id, atc_code,
                    atc.get("level1"), atc.get("level1_description"),
                    atc.get("level2"), atc.get("level2_description"),
                    atc.get("level3"), atc.get("level3_description"),
                    atc.get("level4"), atc.get("level4_description"),
                    atc.get("level5"), atc.get("level5_description"),
                    atc.get("who_name"),
                ))

            if batch:
                with conn.cursor() as cur:
                    psycopg2.extras.execute_values(cur, """
                        INSERT INTO pharmaco_db.atc_classification
                        (compound_id, atc_code, level1, level1_desc, level2, level2_desc,
                         level3, level3_desc, level4, level4_desc, level5, level5_desc, who_name)
                        VALUES %s
                    """, batch, page_size=200)
                conn.commit()
                total += len(batch)

            time.sleep(0.2)

        if not data.get("page_meta", {}).get("next"):
            break
        offset += 500
        time.sleep(0.3)

    log(f"  DONE: {total} ATC classifications inserted")
    update_log(conn, log_id, 'completed', total)
    return total

# ============================================================
# 4. BINDINGDB — Additional binding data
# ============================================================
def ingest_bindingdb(conn):
    log("=" * 60)
    log("Ingesting BindingDB data (via TSV API)")
    log("=" * 60)
    log_id = start_log(conn, 'bindingdb', 'binding_data')

    # BindingDB provides data for specific targets via API
    # We'll query for our top targets
    with conn.cursor() as cur:
        cur.execute("""
            SELECT id, uniprot_id, gene_name FROM pharmaco_db.targets
            WHERE uniprot_id IS NOT NULL
            AND gene_name IS NOT NULL
            ORDER BY num_approved_drugs DESC NULLS LAST
            LIMIT 200
        """)
        targets = cur.fetchall()

    compound_smiles_map = {}
    with conn.cursor() as cur:
        cur.execute("SELECT canonical_smiles, id FROM pharmaco_db.compounds LIMIT 500000")
        compound_smiles_map = dict(cur.fetchall())

    log(f"  Querying BindingDB for {len(targets)} top targets")
    total = 0

    for i, (target_id, uniprot_id, gene_name) in enumerate(targets):
        if i % 20 == 0 and i > 0:
            log(f"    Progress: {i}/{len(targets)} ({total} records)")

        # BindingDB API for UniProt ID
        url = f"https://bindingdb.org/axis2/services/BDBService/getLigandsByUniprot?uniprot={uniprot_id}&response=application/json"
        data = api_get(url, timeout=30)

        if not data:
            time.sleep(0.5)
            continue

        # Parse BindingDB response
        affinities = data.get("getLigandsByUniprotResponse", {}).get("affinities", [])
        if not affinities:
            time.sleep(0.3)
            continue

        if isinstance(affinities, dict):
            affinities = [affinities]

        batch = []
        for aff in affinities[:100]:  # Limit per target
            smiles = aff.get("ligand_smiles")
            if not smiles:
                continue

            compound_id = compound_smiles_map.get(smiles)

            # Parse activity values
            ki = _parse_float(aff.get("ki"))
            ic50 = _parse_float(aff.get("ic50"))
            kd = _parse_float(aff.get("kd"))
            ec50 = _parse_float(aff.get("ec50"))

            # Add the best available measurement
            for atype, value in [("Ki", ki), ("IC50", ic50), ("Kd", kd), ("EC50", ec50)]:
                if value is not None and value > 0:
                    pchembl = _to_pchembl(value)
                    batch.append((
                        compound_id,
                        target_id,
                        None,  # assay_id
                        atype,
                        "=",
                        value,
                        "nM",
                        pchembl,
                        _classify_pchembl(pchembl),
                        "bindingdb",
                        aff.get("monomerid"),
                        None,  # data_validity
                        False,
                    ))

        if batch:
            with conn.cursor() as cur:
                psycopg2.extras.execute_values(cur, """
                    INSERT INTO pharmaco_db.bioactivities
                    (compound_id, target_id, assay_id, activity_type, relation,
                     value, units, pchembl_value, activity_class, source, source_id,
                     data_validity, potential_duplicate)
                    VALUES %s
                """, batch, page_size=200)
            conn.commit()
            total += len(batch)

        time.sleep(0.5)  # BindingDB rate limit

    log(f"  DONE: {total} BindingDB records inserted")
    update_log(conn, log_id, 'completed', total)
    return total

def _parse_float(v):
    if v is None:
        return None
    try:
        s = str(v).replace(">", "").replace("<", "").replace("~", "").strip()
        return float(s)
    except (ValueError, TypeError):
        return None

def _to_pchembl(value_nm):
    """Convert nM value to pChEMBL (-log10(M))."""
    if value_nm and value_nm > 0:
        import math
        return -math.log10(value_nm * 1e-9)
    return None

def _classify_pchembl(pchembl):
    if pchembl is None:
        return None
    if pchembl >= 6.0:
        return "active"
    elif pchembl >= 5.0:
        return "intermediate"
    return "inactive"

# ============================================================
# MAIN
# ============================================================
def main():
    log("=" * 60)
    log("PharmacoDB — Extra Sources Ingestion")
    log("=" * 60)

    conn = get_conn()
    try:
        ingest_atc_classification(conn)
        ingest_guide_to_pharmacology(conn)
        ingest_open_targets(conn)
        ingest_bindingdb(conn)
    except Exception as e:
        log(f"FATAL ERROR: {e}")
        import traceback
        traceback.print_exc()
    finally:
        conn.close()

if __name__ == "__main__":
    main()

