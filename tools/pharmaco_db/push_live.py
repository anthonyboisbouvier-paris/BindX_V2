#!/usr/bin/env python3
"""
PharmacoDB — Live Push to Supabase
Pushes data from local PostgreSQL to Supabase in parallel batches.
Can be run while ingestion is still in progress.
"""

import psycopg2
import psycopg2.extras
import time
import sys
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

LOCAL_DB = {
    "host": "localhost", "port": 5433,
    "dbname": "pharmaco", "user": "postgres", "password": "postgres"
}

SUPABASE_DB = {
    "host": "aws-0-us-west-2.pooler.supabase.com", "port": 6543,
    "dbname": "postgres",
    "user": "postgres.webkntghfzscrnuixfba",
    "password": "n1qsKxvAWJV954vR"
}

BATCH_SIZE = 500
PARALLEL_CONNECTIONS = 3  # Supabase pooler limit is typically 15-20

def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

def get_local_conn():
    return psycopg2.connect(**LOCAL_DB)

def get_supa_conn():
    return psycopg2.connect(**SUPABASE_DB, connect_timeout=10)

def count_rows(conn, table):
    with conn.cursor() as cur:
        cur.execute(f"SELECT COUNT(*) FROM pharmaco_db.{table}")
        return cur.fetchone()[0]

# ============================================================
# Generic batch push
# ============================================================

def push_table(table, columns, conflict_col=None, order_col="id"):
    """Push a table from local to Supabase in batches."""
    log(f"\n>>> Pushing {table}...")

    local = get_local_conn()
    local_count = count_rows(local, table)

    supa = get_supa_conn()
    supa_count = count_rows(supa, table)

    log(f"  Local: {local_count:,} rows | Supabase: {supa_count:,} rows")

    if local_count == 0:
        log(f"  No data to push, skipping.")
        local.close()
        supa.close()
        return 0

    if supa_count >= local_count:
        log(f"  Supabase already up to date, skipping.")
        local.close()
        supa.close()
        return 0

    col_list = ", ".join(columns)
    placeholders = ", ".join(["%s"] * len(columns))

    conflict_clause = ""
    if conflict_col:
        # Upsert: update all non-conflict columns
        update_cols = [c for c in columns if c != conflict_col]
        if update_cols:
            update_set = ", ".join([f"{c} = EXCLUDED.{c}" for c in update_cols])
            conflict_clause = f"ON CONFLICT ({conflict_col}) DO UPDATE SET {update_set}"
        else:
            conflict_clause = f"ON CONFLICT ({conflict_col}) DO NOTHING"

    insert_sql = f"""
        INSERT INTO pharmaco_db.{table} ({col_list})
        VALUES %s
        {conflict_clause}
    """

    # Read all from local, push in batches
    total_pushed = 0
    offset = 0

    while True:
        with local.cursor() as cur:
            cur.execute(f"SELECT {col_list} FROM pharmaco_db.{table} ORDER BY {order_col} LIMIT {BATCH_SIZE} OFFSET {offset}")
            rows = cur.fetchall()

        if not rows:
            break

        try:
            with supa.cursor() as cur:
                psycopg2.extras.execute_values(cur, insert_sql, rows, page_size=200)
            supa.commit()
            total_pushed += len(rows)

            if total_pushed % 5000 == 0:
                log(f"  Pushed {total_pushed:,} / {local_count:,} ({100*total_pushed//local_count}%)")
        except Exception as e:
            supa.rollback()
            log(f"  ERROR at offset {offset}: {e}")
            # Try to reconnect
            try:
                supa.close()
            except:
                pass
            supa = get_supa_conn()

        offset += BATCH_SIZE

    log(f"  DONE: {total_pushed:,} rows pushed to {table}")
    local.close()
    supa.close()
    return total_pushed

# ============================================================
# Table definitions
# ============================================================

TABLES = [
    {
        "table": "targets",
        "columns": [
            "id", "chembl_id", "uniprot_id", "pref_name", "gene_name", "organism", "tax_id",
            "target_type", "protein_family", "protein_class_l1", "protein_class_l2", "protein_class_l3",
            "uniprot_function", "uniprot_subcellular",
            "go_molecular_function", "go_biological_process", "go_cellular_component",
            "pathway_kegg", "pathway_reactome", "pdb_ids",
            "ncbi_gene_id", "ncbi_gene_summary",
            "is_druggable", "num_approved_drugs",
            "sequence_length", "mass", "ec_number",
            "signal_peptide", "transmembrane_regions", "active_site", "binding_sites",
            "disulfide_bonds", "glycosylation_sites", "phosphorylation_sites",
            "domains", "interpro_ids", "pfam_ids",
            "tissue_expression", "tissue_specificity",
            "disease_associations", "involvement_in_disease",
            "keywords", "protein_existence",
            "cross_refs_omim", "cross_refs_orphanet", "cross_refs_pharmgkb",
            "cross_refs_reactome", "cross_refs_string", "cross_refs_intact",
            "isoforms", "alternative_names", "chromosome", "gene_location",
        ],
        "conflict": "chembl_id",
    },
    {
        "table": "compounds",
        "columns": [
            "id", "canonical_smiles", "inchi", "inchi_key", "chembl_id", "pubchem_cid",
            "drugbank_id", "bindingdb_id", "pref_name",
            "molecular_weight", "alogp", "hba", "hbd", "psa", "rtb",
            "num_ro5_violations", "aromatic_rings", "heavy_atoms", "qed_weighted",
            "molecular_formula", "is_drug", "max_phase", "is_natural_product",
        ],
        "conflict": "chembl_id",
    },
    {
        "table": "assays",
        "columns": [
            "id", "chembl_id", "description", "assay_type", "assay_category",
            "src_id", "src_description", "journal", "year", "doi",
            "target_id", "target_chembl_id",
            "assay_organism", "assay_cell_type", "confidence_score",
        ],
        "conflict": "chembl_id",
    },
    {
        "table": "bioactivities",
        "columns": [
            "id", "compound_id", "target_id", "assay_id",
            "activity_type", "relation", "value", "units", "pchembl_value",
            "activity_class", "source", "source_id",
            "data_validity", "potential_duplicate",
        ],
        "conflict": None,  # No unique constraint, use id
    },
    {
        "table": "drug_mechanisms",
        "columns": [
            "id", "compound_id", "target_id", "mechanism_of_action",
            "action_type", "direct_interaction", "molecular_mechanism", "selectivity_comment",
        ],
        "conflict": None,
    },
    {
        "table": "drug_indications",
        "columns": [
            "id", "compound_id", "mesh_id", "mesh_heading",
            "efo_id", "efo_term", "max_phase", "indication_refs",
        ],
        "conflict": None,
    },
    {
        "table": "cross_references",
        "columns": [
            "id", "compound_id", "target_id", "db_name", "db_id", "url",
        ],
        "conflict": None,
    },
]

def main():
    start = time.time()
    log("=" * 60)
    log("PharmacoDB — LIVE PUSH TO SUPABASE")
    log(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log("=" * 60)

    # Test Supabase connection
    try:
        supa = get_supa_conn()
        supa.close()
        log("Supabase connection OK")
    except Exception as e:
        log(f"FATAL: Cannot connect to Supabase: {e}")
        sys.exit(1)

    total = 0
    for tdef in TABLES:
        try:
            n = push_table(
                tdef["table"],
                tdef["columns"],
                conflict_col=tdef.get("conflict"),
            )
            total += n
        except Exception as e:
            log(f"ERROR pushing {tdef['table']}: {e}")
            import traceback
            traceback.print_exc()

    elapsed = time.time() - start
    log(f"\n{'='*60}")
    log(f"PUSH COMPLETE: {total:,} total rows pushed in {int(elapsed)}s")

    # Show final Supabase stats
    try:
        supa = get_supa_conn()
        for tdef in TABLES:
            n = count_rows(supa, tdef["table"])
            log(f"  {tdef['table']:>20}: {n:>12,}")
        supa.close()
    except Exception as e:
        log(f"Could not fetch final stats: {e}")

    log("=" * 60)

if __name__ == "__main__":
    main()
