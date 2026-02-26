#!/usr/bin/env python3
"""
PharmacoDB — ChEMBL Bulk Ingestion Pipeline
Downloads compounds, targets, assays, and bioactivities from ChEMBL API.
Focuses on high-quality measured data for SAR/ML.
"""

import json
import time
import sys
import os
import urllib.request
import urllib.error
import psycopg2
import psycopg2.extras
from datetime import datetime

# ============================================================
# Configuration
# ============================================================
CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
DB_CONFIG = {
    "host": "localhost",
    "port": 5433,
    "dbname": "pharmaco",
    "user": "postgres",
    "password": "pharmaco_secret"
}
BATCH_SIZE = 500        # ChEMBL API max per page
INSERT_BATCH = 1000     # DB insert batch size
MAX_RETRIES = 3
RETRY_DELAY = 5

# Target families to prioritize (druggable genome)
PRIORITY_FAMILIES = [
    "Kinase", "GPCR", "Ion channel", "Nuclear receptor",
    "Protease", "Phosphodiesterase", "Transporter",
    "Epigenetic regulator", "Cytochrome P450"
]

# ============================================================
# Helpers
# ============================================================
def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

def api_get(url, retries=MAX_RETRIES):
    """Fetch JSON from ChEMBL API with retry logic."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={
                "Accept": "application/json",
                "User-Agent": "PharmacoDB-Ingestion/1.0"
            })
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.loads(resp.read().decode())
        except (urllib.error.URLError, urllib.error.HTTPError, Exception) as e:
            if attempt < retries - 1:
                wait = RETRY_DELAY * (attempt + 1)
                log(f"  Retry {attempt+1}/{retries} after {wait}s: {e}")
                time.sleep(wait)
            else:
                log(f"  FAILED after {retries} retries: {e}")
                return None
    return None

def get_conn():
    return psycopg2.connect(**DB_CONFIG)

def update_log(conn, log_id, status, rows=None, error=None):
    with conn.cursor() as cur:
        if rows is not None:
            cur.execute("""
                UPDATE pharmaco_db.ingestion_log
                SET status=%s, rows_inserted=%s, error_message=%s, completed_at=NOW()
                WHERE id=%s
            """, (status, rows, error, log_id))
        else:
            cur.execute("""
                UPDATE pharmaco_db.ingestion_log
                SET status=%s, error_message=%s, completed_at=NOW()
                WHERE id=%s
            """, (status, error, log_id))
    conn.commit()

def start_log(conn, source, step):
    with conn.cursor() as cur:
        cur.execute("""
            INSERT INTO pharmaco_db.ingestion_log (source, step, status)
            VALUES (%s, %s, 'running') RETURNING id
        """, (source, step))
        log_id = cur.fetchone()[0]
    conn.commit()
    return log_id

# ============================================================
# 1. TARGETS — All human single-protein targets
# ============================================================
def ingest_targets(conn):
    log("=" * 60)
    log("PHASE 1: Ingesting ChEMBL Targets")
    log("=" * 60)
    log_id = start_log(conn, 'chembl', 'targets')

    total_inserted = 0
    offset = 0
    batch = []

    while True:
        url = (f"{CHEMBL_API}/target.json?"
               f"target_type=SINGLE%20PROTEIN"
               f"&target_organism=Homo%20sapiens"
               f"&limit={BATCH_SIZE}&offset={offset}")
        data = api_get(url)
        if not data or not data.get("targets"):
            break

        targets = data["targets"]
        log(f"  Fetched {len(targets)} targets (offset={offset})")

        for t in targets:
            # Extract target components for UniProt
            components = t.get("target_components", [])
            uniprot_id = None
            if components:
                for comp in components:
                    accessions = comp.get("target_component_xrefs", [])
                    for xref in accessions:
                        if xref.get("xref_src_db") == "UniProt":
                            uniprot_id = xref.get("xref_id")
                            break
                    if uniprot_id:
                        break

            batch.append({
                "chembl_id": t.get("target_chembl_id"),
                "uniprot_id": uniprot_id,
                "pref_name": t.get("pref_name", "Unknown"),
                "target_type": t.get("target_type"),
                "organism": t.get("organism", "Homo sapiens"),
                "tax_id": t.get("tax_id"),
            })

        if len(batch) >= INSERT_BATCH:
            n = _insert_targets(conn, batch)
            total_inserted += n
            batch = []
            log(f"  Total targets: {total_inserted}")

        if not data.get("page_meta", {}).get("next"):
            break
        offset += BATCH_SIZE
        time.sleep(0.3)

    if batch:
        n = _insert_targets(conn, batch)
        total_inserted += n

    log(f"  DONE: {total_inserted} targets inserted")
    update_log(conn, log_id, 'completed', total_inserted)
    return total_inserted

def _insert_targets(conn, batch):
    if not batch:
        return 0
    with conn.cursor() as cur:
        psycopg2.extras.execute_values(cur, """
            INSERT INTO pharmaco_db.targets (chembl_id, uniprot_id, pref_name, target_type, organism, tax_id)
            VALUES %s
            ON CONFLICT (chembl_id) DO UPDATE SET
                uniprot_id = EXCLUDED.uniprot_id,
                pref_name = EXCLUDED.pref_name,
                updated_at = NOW()
        """, [(
            t["chembl_id"], t["uniprot_id"], t["pref_name"],
            t["target_type"], t["organism"], t["tax_id"]
        ) for t in batch], page_size=500)
    conn.commit()
    return len(batch)

# ============================================================
# 2. COMPOUNDS — Drug-like molecules with properties
# ============================================================
def ingest_compounds(conn):
    log("=" * 60)
    log("PHASE 2: Ingesting ChEMBL Compounds")
    log("=" * 60)
    log_id = start_log(conn, 'chembl', 'compounds')

    total_inserted = 0
    offset = 0
    batch = []

    while True:
        # Fetch molecules with structures
        url = (f"{CHEMBL_API}/molecule.json?"
               f"molecule_properties__mw_freebase__gte=100"
               f"&molecule_properties__mw_freebase__lte=900"
               f"&molecule_structures__canonical_smiles__isnull=false"
               f"&limit={BATCH_SIZE}&offset={offset}")
        data = api_get(url)
        if not data or not data.get("molecules"):
            break

        molecules = data["molecules"]
        log(f"  Fetched {len(molecules)} compounds (offset={offset}, total so far: {total_inserted})")

        for m in molecules:
            structs = m.get("molecule_structures") or {}
            props = m.get("molecule_properties") or {}
            smiles = structs.get("canonical_smiles")
            if not smiles:
                continue

            batch.append({
                "chembl_id": m.get("molecule_chembl_id"),
                "canonical_smiles": smiles,
                "inchi": structs.get("standard_inchi"),
                "inchi_key": structs.get("standard_inchi_key"),
                "pref_name": m.get("pref_name"),
                "molecular_weight": _float(props.get("mw_freebase")),
                "alogp": _float(props.get("alogp")),
                "hba": _int(props.get("hba")),
                "hbd": _int(props.get("hbd")),
                "psa": _float(props.get("psa")),
                "rtb": _int(props.get("rtb")),
                "num_ro5_violations": _int(props.get("num_ro5_violations")),
                "aromatic_rings": _int(props.get("aromatic_rings")),
                "heavy_atoms": _int(props.get("heavy_atom_count")),
                "qed_weighted": _float(props.get("qed_weighted")),
                "molecular_formula": props.get("full_molformula"),
                "max_phase": _int(m.get("max_phase")),
                "is_natural_product": m.get("natural_product") == 1 if m.get("natural_product") is not None else None,
            })

        if len(batch) >= INSERT_BATCH:
            n = _insert_compounds(conn, batch)
            total_inserted += n
            batch = []
            if total_inserted % 10000 == 0:
                log(f"  Total compounds: {total_inserted}")

        if not data.get("page_meta", {}).get("next"):
            break
        offset += BATCH_SIZE
        time.sleep(0.2)

    if batch:
        n = _insert_compounds(conn, batch)
        total_inserted += n

    log(f"  DONE: {total_inserted} compounds inserted")
    update_log(conn, log_id, 'completed', total_inserted)
    return total_inserted

def _insert_compounds(conn, batch):
    if not batch:
        return 0
    with conn.cursor() as cur:
        psycopg2.extras.execute_values(cur, """
            INSERT INTO pharmaco_db.compounds (
                chembl_id, canonical_smiles, inchi, inchi_key, pref_name,
                molecular_weight, alogp, hba, hbd, psa, rtb,
                num_ro5_violations, aromatic_rings, heavy_atoms, qed_weighted,
                molecular_formula, max_phase, is_natural_product
            ) VALUES %s
            ON CONFLICT (chembl_id) DO UPDATE SET
                canonical_smiles = EXCLUDED.canonical_smiles,
                molecular_weight = EXCLUDED.molecular_weight,
                alogp = EXCLUDED.alogp,
                updated_at = NOW()
        """, [(
            c["chembl_id"], c["canonical_smiles"], c["inchi"], c["inchi_key"],
            c["pref_name"], c["molecular_weight"], c["alogp"], c["hba"],
            c["hbd"], c["psa"], c["rtb"], c["num_ro5_violations"],
            c["aromatic_rings"], c["heavy_atoms"], c["qed_weighted"],
            c["molecular_formula"], c["max_phase"], c["is_natural_product"]
        ) for c in batch], page_size=500)
    conn.commit()
    return len(batch)

def _float(v):
    try:
        return float(v) if v is not None else None
    except (ValueError, TypeError):
        return None

def _int(v):
    try:
        return int(v) if v is not None else None
    except (ValueError, TypeError):
        return None

# ============================================================
# 3. ASSAYS — Experimental protocols
# ============================================================
def ingest_assays(conn):
    log("=" * 60)
    log("PHASE 3: Ingesting ChEMBL Assays")
    log("=" * 60)
    log_id = start_log(conn, 'chembl', 'assays')

    # Build target chembl_id -> internal id mapping
    target_map = {}
    with conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.targets")
        target_map = dict(cur.fetchall())
    log(f"  Target map: {len(target_map)} entries")

    total_inserted = 0
    offset = 0
    batch = []

    # Focus on binding and functional assays with high confidence
    for assay_type in ["B", "F"]:
        offset = 0
        while True:
            url = (f"{CHEMBL_API}/assay.json?"
                   f"assay_type={assay_type}"
                   f"&confidence_score__gte=7"
                   f"&assay_organism=Homo%20sapiens"
                   f"&limit={BATCH_SIZE}&offset={offset}")
            data = api_get(url)
            if not data or not data.get("assays"):
                break

            assays = data["assays"]
            log(f"  Fetched {len(assays)} {assay_type}-type assays (offset={offset})")

            for a in assays:
                target_chembl = a.get("target_chembl_id")
                batch.append({
                    "chembl_id": a.get("assay_chembl_id"),
                    "description": (a.get("description") or "")[:500],
                    "assay_type": a.get("assay_type"),
                    "assay_category": a.get("assay_category"),
                    "target_id": target_map.get(target_chembl),
                    "target_chembl_id": target_chembl,
                    "assay_organism": a.get("assay_organism"),
                    "assay_cell_type": a.get("assay_cell_type"),
                    "confidence_score": _int(a.get("confidence_score")),
                })

            if len(batch) >= INSERT_BATCH:
                n = _insert_assays(conn, batch)
                total_inserted += n
                batch = []
                if total_inserted % 10000 == 0:
                    log(f"  Total assays: {total_inserted}")

            if not data.get("page_meta", {}).get("next"):
                break
            offset += BATCH_SIZE
            time.sleep(0.3)

    if batch:
        n = _insert_assays(conn, batch)
        total_inserted += n

    log(f"  DONE: {total_inserted} assays inserted")
    update_log(conn, log_id, 'completed', total_inserted)
    return total_inserted

def _insert_assays(conn, batch):
    if not batch:
        return 0
    with conn.cursor() as cur:
        psycopg2.extras.execute_values(cur, """
            INSERT INTO pharmaco_db.assays (
                chembl_id, description, assay_type, assay_category,
                target_id, target_chembl_id,
                assay_organism, assay_cell_type, confidence_score
            ) VALUES %s
            ON CONFLICT (chembl_id) DO NOTHING
        """, [(
            a["chembl_id"], a["description"], a["assay_type"], a["assay_category"],
            a["target_id"], a["target_chembl_id"],
            a["assay_organism"], a["assay_cell_type"], a["confidence_score"]
        ) for a in batch], page_size=500)
    conn.commit()
    return len(batch)

# ============================================================
# 4. BIOACTIVITIES — The core data (IC50, Ki, EC50, Kd)
# ============================================================
def ingest_bioactivities(conn):
    log("=" * 60)
    log("PHASE 4: Ingesting ChEMBL Bioactivities")
    log("=" * 60)
    log_id = start_log(conn, 'chembl', 'bioactivities')

    # Build lookup maps
    compound_map = {}
    with conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.compounds")
        compound_map = dict(cur.fetchall())
    log(f"  Compound map: {len(compound_map)} entries")

    target_map = {}
    with conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.targets")
        target_map = dict(cur.fetchall())
    log(f"  Target map: {len(target_map)} entries")

    assay_map = {}
    with conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.assays")
        assay_map = dict(cur.fetchall())
    log(f"  Assay map: {len(assay_map)} entries")

    total_inserted = 0
    total_skipped = 0
    offset = 0
    batch = []

    # Fetch activities with pchembl_value (standardized, high quality)
    while True:
        url = (f"{CHEMBL_API}/activity.json?"
               f"pchembl_value__isnull=false"
               f"&target_organism=Homo%20sapiens"
               f"&limit={BATCH_SIZE}&offset={offset}")
        data = api_get(url)
        if not data or not data.get("activities"):
            break

        activities = data["activities"]
        log(f"  Fetched {len(activities)} activities (offset={offset}, inserted={total_inserted}, skipped={total_skipped})")

        for act in activities:
            mol_chembl = act.get("molecule_chembl_id")
            target_chembl = act.get("target_chembl_id")
            assay_chembl = act.get("assay_chembl_id")

            compound_id = compound_map.get(mol_chembl)
            if not compound_id:
                total_skipped += 1
                continue

            target_id = target_map.get(target_chembl)
            assay_id = assay_map.get(assay_chembl)

            # Standardize value to nM
            value = _float(act.get("value"))
            units = act.get("units")

            batch.append({
                "compound_id": compound_id,
                "target_id": target_id,
                "assay_id": assay_id,
                "activity_type": act.get("standard_type") or act.get("type"),
                "relation": act.get("standard_relation") or act.get("relation"),
                "value": _float(act.get("standard_value")),
                "units": act.get("standard_units") or units,
                "pchembl_value": _float(act.get("pchembl_value")),
                "activity_class": _classify_activity(_float(act.get("pchembl_value"))),
                "source": "chembl",
                "source_id": act.get("activity_id"),
                "data_validity": act.get("data_validity_comment"),
                "potential_duplicate": act.get("potential_duplicate", False),
            })

        if len(batch) >= INSERT_BATCH:
            n = _insert_activities(conn, batch)
            total_inserted += n
            batch = []
            if total_inserted % 50000 == 0:
                log(f"  >>> Total activities: {total_inserted}")

        if not data.get("page_meta", {}).get("next"):
            break
        offset += BATCH_SIZE
        time.sleep(0.2)

    if batch:
        n = _insert_activities(conn, batch)
        total_inserted += n

    log(f"  DONE: {total_inserted} activities inserted, {total_skipped} skipped")
    update_log(conn, log_id, 'completed', total_inserted)
    return total_inserted

def _classify_activity(pchembl):
    if pchembl is None:
        return None
    if pchembl >= 6.0:  # < 1 uM
        return "active"
    elif pchembl >= 5.0:  # 1-10 uM
        return "intermediate"
    else:
        return "inactive"

def _insert_activities(conn, batch):
    if not batch:
        return 0
    with conn.cursor() as cur:
        psycopg2.extras.execute_values(cur, """
            INSERT INTO pharmaco_db.bioactivities (
                compound_id, target_id, assay_id,
                activity_type, relation, value, units, pchembl_value,
                activity_class, source, source_id,
                data_validity, potential_duplicate
            ) VALUES %s
        """, [(
            a["compound_id"], a["target_id"], a["assay_id"],
            a["activity_type"], a["relation"], a["value"], a["units"],
            a["pchembl_value"], a["activity_class"], a["source"],
            str(a["source_id"]) if a["source_id"] else None,
            a["data_validity"], a["potential_duplicate"]
        ) for a in batch], page_size=500)
    conn.commit()
    return len(batch)

# ============================================================
# 5. DRUG MECHANISMS — Mechanism of action for approved drugs
# ============================================================
def ingest_mechanisms(conn):
    log("=" * 60)
    log("PHASE 5: Ingesting Drug Mechanisms")
    log("=" * 60)
    log_id = start_log(conn, 'chembl', 'mechanisms')

    compound_map = {}
    with conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.compounds")
        compound_map = dict(cur.fetchall())

    target_map = {}
    with conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.targets")
        target_map = dict(cur.fetchall())

    total = 0
    offset = 0

    while True:
        url = f"{CHEMBL_API}/mechanism.json?limit={BATCH_SIZE}&offset={offset}"
        data = api_get(url)
        if not data or not data.get("mechanisms"):
            break

        mechs = data["mechanisms"]
        log(f"  Fetched {len(mechs)} mechanisms (offset={offset})")

        batch = []
        for m in mechs:
            mol_chembl = m.get("molecule_chembl_id")
            target_chembl = m.get("target_chembl_id")
            compound_id = compound_map.get(mol_chembl)
            if not compound_id:
                continue

            batch.append((
                compound_id,
                target_map.get(target_chembl),
                m.get("mechanism_of_action"),
                m.get("action_type"),
                m.get("direct_interaction"),
                m.get("molecular_mechanism"),
                m.get("selectivity_comment"),
            ))

        if batch:
            with conn.cursor() as cur:
                psycopg2.extras.execute_values(cur, """
                    INSERT INTO pharmaco_db.drug_mechanisms
                    (compound_id, target_id, mechanism_of_action, action_type,
                     direct_interaction, molecular_mechanism, selectivity_comment)
                    VALUES %s
                """, batch, page_size=500)
            conn.commit()
            total += len(batch)

        if not data.get("page_meta", {}).get("next"):
            break
        offset += BATCH_SIZE
        time.sleep(0.3)

    log(f"  DONE: {total} mechanisms inserted")
    update_log(conn, log_id, 'completed', total)
    return total

# ============================================================
# 6. DRUG INDICATIONS
# ============================================================
def ingest_indications(conn):
    log("=" * 60)
    log("PHASE 6: Ingesting Drug Indications")
    log("=" * 60)
    log_id = start_log(conn, 'chembl', 'indications')

    compound_map = {}
    with conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.compounds")
        compound_map = dict(cur.fetchall())

    total = 0
    offset = 0

    while True:
        url = f"{CHEMBL_API}/drug_indication.json?limit={BATCH_SIZE}&offset={offset}"
        data = api_get(url)
        if not data or not data.get("drug_indications"):
            break

        indications = data["drug_indications"]
        log(f"  Fetched {len(indications)} indications (offset={offset})")

        batch = []
        for ind in indications:
            mol_chembl = ind.get("molecule_chembl_id")
            compound_id = compound_map.get(mol_chembl)
            if not compound_id:
                continue

            batch.append((
                compound_id,
                ind.get("mesh_id"),
                ind.get("mesh_heading"),
                ind.get("efo_id"),
                ind.get("efo_term"),
                _int(ind.get("max_phase_for_ind")),
                ind.get("indication_refs"),
            ))

        if batch:
            with conn.cursor() as cur:
                psycopg2.extras.execute_values(cur, """
                    INSERT INTO pharmaco_db.drug_indications
                    (compound_id, mesh_id, mesh_heading, efo_id, efo_term,
                     max_phase, indication_refs)
                    VALUES %s
                """, batch, page_size=500)
            conn.commit()
            total += len(batch)

        if not data.get("page_meta", {}).get("next"):
            break
        offset += BATCH_SIZE
        time.sleep(0.3)

    log(f"  DONE: {total} indications inserted")
    update_log(conn, log_id, 'completed', total)
    return total

# ============================================================
# MAIN
# ============================================================
def main():
    log("=" * 60)
    log("PharmacoDB — ChEMBL Ingestion Pipeline")
    log("=" * 60)

    conn = get_conn()

    try:
        # Phase 1: Targets first (needed for FK resolution)
        n_targets = ingest_targets(conn)

        # Phase 2: Compounds
        n_compounds = ingest_compounds(conn)

        # Phase 3: Assays
        n_assays = ingest_assays(conn)

        # Phase 4: Bioactivities (the big one)
        n_activities = ingest_bioactivities(conn)

        # Phase 5: Drug mechanisms
        n_mechanisms = ingest_mechanisms(conn)

        # Phase 6: Drug indications
        n_indications = ingest_indications(conn)

        # Summary
        log("=" * 60)
        log("CHEMBL INGESTION COMPLETE")
        log(f"  Targets:       {n_targets:>10,}")
        log(f"  Compounds:     {n_compounds:>10,}")
        log(f"  Assays:        {n_assays:>10,}")
        log(f"  Activities:    {n_activities:>10,}")
        log(f"  Mechanisms:    {n_mechanisms:>10,}")
        log(f"  Indications:   {n_indications:>10,}")
        log("=" * 60)

    except Exception as e:
        log(f"FATAL ERROR: {e}")
        import traceback
        traceback.print_exc()
    finally:
        conn.close()

if __name__ == "__main__":
    main()
