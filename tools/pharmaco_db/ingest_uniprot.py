#!/usr/bin/env python3
"""
PharmacoDB — UniProt Target Enrichment Pipeline
Enriches target records with UniProt protein data:
- Function, subcellular location, pathways, GO terms
- Gene names, PDB cross-references
- NCBI Gene IDs
"""

import json
import time
import sys
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

UNIPROT_API = "https://rest.uniprot.org/uniprotkb"
MAX_RETRIES = 3

def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

def api_get(url, retries=MAX_RETRIES):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={
                "Accept": "application/json",
                "User-Agent": "PharmacoDB-Ingestion/1.0"
            })
            with urllib.request.urlopen(req, timeout=30) as resp:
                return json.loads(resp.read().decode())
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(3 * (attempt + 1))
            else:
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

def enrich_targets(conn):
    log("=" * 60)
    log("UniProt Target Enrichment")
    log("=" * 60)
    log_id = start_log(conn, 'uniprot', 'target_enrichment')

    # Get all targets with UniProt IDs
    with conn.cursor() as cur:
        cur.execute("""
            SELECT id, uniprot_id, chembl_id FROM pharmaco_db.targets
            WHERE uniprot_id IS NOT NULL
            AND uniprot_function IS NULL
            ORDER BY id
        """)
        targets = cur.fetchall()

    log(f"  Found {len(targets)} targets to enrich")
    enriched = 0
    errors = 0

    for i, (target_id, uniprot_id, chembl_id) in enumerate(targets):
        if i % 100 == 0 and i > 0:
            log(f"  Progress: {i}/{len(targets)} ({enriched} enriched, {errors} errors)")

        url = f"{UNIPROT_API}/{uniprot_id}.json"
        data = api_get(url)
        if not data:
            errors += 1
            continue

        # Extract function
        function_text = None
        comments = data.get("comments", [])
        for c in comments:
            if c.get("commentType") == "FUNCTION":
                texts = c.get("texts", [])
                if texts:
                    function_text = texts[0].get("value", "")[:2000]
                    break

        # Extract subcellular location
        subcellular = None
        for c in comments:
            if c.get("commentType") == "SUBCELLULAR LOCATION":
                locs = c.get("subcellularLocations", [])
                parts = []
                for loc in locs:
                    location = loc.get("location", {})
                    parts.append(location.get("value", ""))
                subcellular = "; ".join(parts)[:500] if parts else None
                break

        # Extract gene name
        genes = data.get("genes", [])
        gene_name = None
        if genes:
            gene_name = genes[0].get("geneName", {}).get("value")

        # Extract GO terms
        go_mf = []
        go_bp = []
        go_cc = []
        xrefs = data.get("uniProtKBCrossReferences", [])
        for xref in xrefs:
            if xref.get("database") == "GO":
                go_id = xref.get("id", "")
                props = xref.get("properties", [])
                for p in props:
                    if p.get("key") == "GoTerm":
                        term = p.get("value", "")
                        if term.startswith("F:"):
                            go_mf.append(term[2:])
                        elif term.startswith("P:"):
                            go_bp.append(term[2:])
                        elif term.startswith("C:"):
                            go_cc.append(term[2:])

        # Extract PDB IDs
        pdb_ids = []
        for xref in xrefs:
            if xref.get("database") == "PDB":
                pdb_ids.append(xref.get("id"))

        # Extract KEGG pathway
        kegg_pathways = []
        for xref in xrefs:
            if xref.get("database") == "KEGG":
                kegg_pathways.append(xref.get("id"))

        # Extract Reactome pathway
        reactome_pathways = []
        for xref in xrefs:
            if xref.get("database") == "Reactome":
                reactome_pathways.append(xref.get("id"))

        # Extract NCBI Gene ID
        ncbi_gene_id = None
        for xref in xrefs:
            if xref.get("database") == "GeneID":
                try:
                    ncbi_gene_id = int(xref.get("id"))
                except (ValueError, TypeError):
                    pass

        # Extract protein family from keywords
        protein_family = None
        keywords = data.get("keywords", [])
        family_keywords = {"Kinase", "Transferase", "Hydrolase", "Receptor",
                          "Ion channel", "Transporter", "Protease", "Oxidoreductase",
                          "G-protein coupled receptor"}
        for kw in keywords:
            if kw.get("name") in family_keywords:
                protein_family = kw["name"]
                break

        # Update target record
        with conn.cursor() as cur:
            cur.execute("""
                UPDATE pharmaco_db.targets SET
                    gene_name = COALESCE(%s, gene_name),
                    uniprot_function = %s,
                    uniprot_subcellular = %s,
                    go_molecular_function = %s,
                    go_biological_process = %s,
                    go_cellular_component = %s,
                    pdb_ids = %s,
                    pathway_kegg = %s,
                    pathway_reactome = %s,
                    ncbi_gene_id = %s,
                    protein_family = COALESCE(protein_family, %s),
                    updated_at = NOW()
                WHERE id = %s
            """, (
                gene_name,
                function_text,
                subcellular,
                go_mf if go_mf else None,
                go_bp if go_bp else None,
                go_cc if go_cc else None,
                pdb_ids if pdb_ids else None,
                kegg_pathways if kegg_pathways else None,
                reactome_pathways if reactome_pathways else None,
                ncbi_gene_id,
                protein_family,
                target_id
            ))
        conn.commit()
        enriched += 1

        # Also add cross-references
        xref_batch = []
        if pdb_ids:
            for pdb_id in pdb_ids[:10]:  # Limit to 10 PDB entries
                xref_batch.append((None, target_id, "PDB", pdb_id,
                                  f"https://www.rcsb.org/structure/{pdb_id}"))
        if ncbi_gene_id:
            xref_batch.append((None, target_id, "NCBI_Gene", str(ncbi_gene_id),
                              f"https://www.ncbi.nlm.nih.gov/gene/{ncbi_gene_id}"))

        if xref_batch:
            with conn.cursor() as cur:
                psycopg2.extras.execute_values(cur, """
                    INSERT INTO pharmaco_db.cross_references
                    (compound_id, target_id, db_name, db_id, url)
                    VALUES %s
                    ON CONFLICT DO NOTHING
                """, xref_batch, page_size=100)
            conn.commit()

        # Rate limit: ~10 requests/second for UniProt
        time.sleep(0.1)

    log(f"  DONE: {enriched} targets enriched, {errors} errors")
    update_log(conn, log_id, 'completed', enriched)
    return enriched

def main():
    log("=" * 60)
    log("PharmacoDB — UniProt Enrichment Pipeline")
    log("=" * 60)

    conn = get_conn()
    try:
        enrich_targets(conn)
    except Exception as e:
        log(f"FATAL ERROR: {e}")
        import traceback
        traceback.print_exc()
    finally:
        conn.close()

if __name__ == "__main__":
    main()
