#!/usr/bin/env python3
"""
PharmacoDB -- Comprehensive UniProt Target Enrichment Pipeline (V2)

Replaces the basic ingest_uniprot.py with full-spectrum protein annotation:

  - Batch API: 50 accessions per request via /uniprotkb/accessions
  - Extracts 40+ fields from each UniProt JSON entry:
      * Sequence length & mass
      * EC number, protein existence level
      * Signal peptide, transmembrane, active site, binding sites
      * Disulfide bonds, glycosylation, phosphorylation counts
      * Domains, InterPro, Pfam identifiers
      * GO terms (MF/BP/CC), KEGG, Reactome pathways
      * PDB cross-references
      * Disease associations, tissue specificity, keywords
      * OMIM, Orphanet, PharmGKB, STRING, IntAct cross-refs
      * Isoforms, alternative names, gene synonyms
      * NCBI Gene ID, chromosome (from Ensembl)
  - Retry logic with exponential backoff (3 attempts)
  - Rate limiting: 0.5s between batch requests
  - Progress logging every 100 targets
  - Graceful handling of missing/null fields
  - Cross-references table insertion (PDB, NCBI Gene, STRING, IntAct)

Prerequisites:
  - Run schema_uniprot_extra.sql FIRST to add new columns
  - Local PostgreSQL: localhost:5433, db=pharmaco, user=postgres

Usage:
  python3 ingest_uniprot_full.py [--limit N] [--force] [--batch-size N]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
import urllib.error
import urllib.request
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

import psycopg2
import psycopg2.extras

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DB_CONFIG: Dict[str, Any] = {
    "host": "localhost",
    "port": 5433,
    "dbname": "pharmaco",
    "user": "postgres",
    "password": "pharmaco_secret",
}

UNIPROT_BATCH_API: str = "https://rest.uniprot.org/uniprotkb/accessions"
UNIPROT_SINGLE_API: str = "https://rest.uniprot.org/uniprotkb"

MAX_RETRIES: int = 3
BATCH_SIZE: int = 80  # UniProt allows up to 100, stay under
RATE_LIMIT_SECONDS: float = 0.3
PARALLEL_WORKERS: int = 2  # concurrent batch fetches (UniProt rate limits aggressively)
REQUEST_TIMEOUT: int = 60


# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------

def log(msg: str) -> None:
    """Print timestamped log message."""
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def log_error(msg: str) -> None:
    """Print timestamped error message."""
    print(f"[{datetime.now().strftime('%H:%M:%S')}] ERROR: {msg}", flush=True)


# ---------------------------------------------------------------------------
# HTTP helpers with retry + exponential backoff
# ---------------------------------------------------------------------------

def _build_opener() -> urllib.request.OpenerDirector:
    """Build a URL opener that follows redirects (for secondary accessions)."""
    return urllib.request.build_opener(urllib.request.HTTPRedirectHandler)


def api_get_json(
    url: str,
    retries: int = MAX_RETRIES,
    follow_redirects: bool = True,
) -> Optional[Dict]:
    """
    Fetch JSON from URL with retry logic.
    Follows 303 redirects (common for secondary UniProt accessions).
    Returns parsed JSON dict or None on failure after all retries.
    """
    opener = _build_opener() if follow_redirects else None

    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={
                "Accept": "application/json",
                "User-Agent": "PharmacoDB-Ingestion/2.0 (BindX; contact@bindx.io)",
            })
            if opener:
                resp = opener.open(req, timeout=REQUEST_TIMEOUT)
            else:
                resp = urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT)
            raw = resp.read().decode("utf-8")
            resp.close()
            return json.loads(raw)
        except urllib.error.HTTPError as e:
            if e.code == 400:
                # Bad request -- likely invalid accession(s), do not retry
                log_error(f"HTTP 400 for {url[:120]}... (bad request, skipping)")
                return None
            if e.code == 429:
                # Rate limited -- back off more aggressively
                wait = 5 * (attempt + 1)
                log(f"  Rate limited (429), waiting {wait}s...")
                time.sleep(wait)
                continue
            if e.code == 404:
                # Not found -- accession deleted/merged
                return None
            log_error(f"HTTP {e.code} for {url[:120]}... (attempt {attempt + 1}/{retries})")
        except Exception as e:
            log_error(f"Request failed: {e} (attempt {attempt + 1}/{retries})")

        if attempt < retries - 1:
            wait = 2 ** attempt  # 1s, 2s, 4s
            time.sleep(wait)

    return None


def fetch_single_with_redirect(accession: str) -> Optional[Dict]:
    """
    Fetch a single UniProt entry, following 303 redirects.
    Used for secondary/obsolete accessions that the batch API ignores.
    """
    url = f"{UNIPROT_SINGLE_API}/{accession}.json"
    entry = api_get_json(url, follow_redirects=True)
    if entry:
        # Tag the entry with the original accession so we can match it back
        entry["_queried_accession"] = accession
    return entry


def fetch_batch(accessions: List[str]) -> List[Dict]:
    """
    Fetch a batch of UniProt entries using the accessions endpoint.
    The batch API only returns primary accessions, so secondary/obsolete
    accessions are fetched individually with redirect following.
    Returns list of parsed entry dicts (tagged with _queried_accession).
    """
    if not accessions:
        return []

    # Step 1: Try batch API for all accessions
    acc_str = ",".join(accessions)
    url = f"{UNIPROT_BATCH_API}?accessions={acc_str}"

    batch_results: List[Dict] = []
    data = api_get_json(url, follow_redirects=False)
    if data and "results" in data:
        batch_results = data["results"]

    # Index returned entries by primary + secondary accessions
    found_accessions: set = set()
    for entry in batch_results:
        primary = entry.get("primaryAccession", "")
        found_accessions.add(primary)
        for sec in entry.get("secondaryAccessions", []):
            found_accessions.add(sec)
        # Tag with queried accession (primary in batch case)
        entry["_queried_accession"] = primary

    # Step 2: Find which accessions were NOT returned by the batch
    missing = [acc for acc in accessions if acc not in found_accessions]

    if missing:
        # Fetch missing ones concurrently (secondary accessions with 303 redirects)
        from concurrent.futures import ThreadPoolExecutor, as_completed
        with ThreadPoolExecutor(max_workers=min(len(missing), 8)) as pool:
            futures = {pool.submit(fetch_single_with_redirect, acc): acc for acc in missing}
            for future in as_completed(futures):
                entry = future.result()
                if entry:
                    batch_results.append(entry)

    return batch_results


# ---------------------------------------------------------------------------
# JSON parsing: extract all fields from a UniProt entry
# ---------------------------------------------------------------------------

def safe_str(value: Any, max_len: int = 4000) -> Optional[str]:
    """Safely convert value to a trimmed string."""
    if value is None:
        return None
    s = str(value).strip()
    return s[:max_len] if s else None


def parse_entry(entry: Dict) -> Dict[str, Any]:
    """
    Parse a single UniProt JSON entry into a flat dict of database fields.
    All keys map to pharmaco_db.targets columns.
    """
    result: Dict[str, Any] = {}

    accession = entry.get("primaryAccession", "")
    result["_accession"] = accession

    # ---- Sequence metadata ----
    seq = entry.get("sequence", {})
    result["sequence_length"] = seq.get("length")
    result["mass"] = seq.get("molWeight")

    # ---- Protein existence ----
    result["protein_existence"] = safe_str(entry.get("proteinExistence"))

    # ---- Protein description & EC number ----
    desc = entry.get("proteinDescription", {})
    rec_name = desc.get("recommendedName", {})

    ec_numbers = rec_name.get("ecNumbers", [])
    if ec_numbers:
        result["ec_number"] = ec_numbers[0].get("value")
    else:
        result["ec_number"] = None

    # Alternative names
    alt_names: List[str] = []
    for alt in desc.get("alternativeNames", []):
        full = alt.get("fullName", {}).get("value")
        if full:
            alt_names.append(full)
        short = alt.get("shortNames", [])
        for s in short:
            v = s.get("value")
            if v:
                alt_names.append(v)
    # Also check submittedName
    for sub in desc.get("submissionNames", []):
        full = sub.get("fullName", {}).get("value")
        if full:
            alt_names.append(full)
    result["alternative_names"] = alt_names if alt_names else None

    # ---- Gene info ----
    genes = entry.get("genes", [])
    gene_name = None
    if genes:
        gene_name = genes[0].get("geneName", {}).get("value")
    result["gene_name"] = gene_name

    # ---- Comments extraction ----
    comments = entry.get("comments", [])
    result.update(_parse_comments(comments))

    # ---- Features extraction ----
    features = entry.get("features", [])
    result.update(_parse_features(features))

    # ---- Keywords ----
    keywords = entry.get("keywords", [])
    kw_names = [kw.get("name") for kw in keywords if kw.get("name")]
    result["keywords"] = kw_names if kw_names else None

    # ---- Protein family from keywords ----
    result["protein_family"] = _extract_protein_family(keywords)

    # ---- Cross-references ----
    xrefs = entry.get("uniProtKBCrossReferences", [])
    result.update(_parse_xrefs(xrefs))

    return result


def _parse_comments(comments: List[Dict]) -> Dict[str, Any]:
    """Extract structured data from UniProt comments."""
    result: Dict[str, Any] = {}

    # -- FUNCTION --
    function_parts: List[str] = []
    for c in comments:
        if c.get("commentType") == "FUNCTION":
            for t in c.get("texts", []):
                val = t.get("value", "").strip()
                if val:
                    function_parts.append(val)
    result["uniprot_function"] = " ".join(function_parts)[:4000] if function_parts else None

    # -- SUBCELLULAR LOCATION --
    subcellular_parts: List[str] = []
    for c in comments:
        if c.get("commentType") == "SUBCELLULAR LOCATION":
            for loc in c.get("subcellularLocations", []):
                location = loc.get("location", {})
                val = location.get("value", "").strip()
                if val:
                    subcellular_parts.append(val)
                # Also capture topology
                topology = loc.get("topology", {})
                tval = topology.get("value", "").strip()
                if tval and tval not in subcellular_parts:
                    subcellular_parts.append(tval)
    result["uniprot_subcellular"] = "; ".join(subcellular_parts)[:2000] if subcellular_parts else None

    # -- TISSUE SPECIFICITY --
    tissue_texts: List[str] = []
    for c in comments:
        if c.get("commentType") == "TISSUE SPECIFICITY":
            for t in c.get("texts", []):
                val = t.get("value", "").strip()
                if val:
                    tissue_texts.append(val)
    result["tissue_specificity"] = " ".join(tissue_texts)[:2000] if tissue_texts else None

    # Parse tissue names from the specificity text for the array column
    result["tissue_expression"] = _extract_tissue_names(result["tissue_specificity"])

    # -- DISEASE --
    disease_names: List[str] = []
    disease_descriptions: List[str] = []
    for c in comments:
        if c.get("commentType") == "DISEASE":
            disease = c.get("disease", {})
            name = disease.get("diseaseId", "").strip()
            if name:
                disease_names.append(name)
            desc = disease.get("description", "").strip()
            if desc:
                disease_descriptions.append(f"{name}: {desc[:300]}" if name else desc[:300])
    result["disease_associations"] = disease_names if disease_names else None
    result["involvement_in_disease"] = "; ".join(disease_descriptions)[:4000] if disease_descriptions else None

    # -- ALTERNATIVE PRODUCTS (isoforms) --
    isoform_count = 0
    for c in comments:
        if c.get("commentType") == "ALTERNATIVE PRODUCTS":
            isoforms = c.get("isoforms", [])
            isoform_count = len(isoforms)
    result["isoforms"] = isoform_count if isoform_count > 0 else None

    return result


def _extract_tissue_names(tissue_text: Optional[str]) -> Optional[List[str]]:
    """
    Extract tissue names from UniProt tissue specificity free text.
    Returns a list of recognized tissue names found in the text.
    """
    if not tissue_text:
        return None

    # Common tissue/organ names to look for in the text
    known_tissues = [
        "brain", "liver", "kidney", "heart", "lung", "spleen", "pancreas",
        "intestine", "colon", "stomach", "skin", "muscle", "bone", "bone marrow",
        "thymus", "lymph node", "blood", "plasma", "serum", "placenta",
        "testis", "ovary", "uterus", "prostate", "breast", "mammary",
        "adrenal", "thyroid", "pituitary", "hypothalamus", "cerebellum",
        "hippocampus", "cortex", "retina", "cornea", "lens", "cochlea",
        "trachea", "esophagus", "bladder", "gallbladder", "salivary gland",
        "adipose", "cartilage", "tendon", "aorta", "vein", "artery",
        "peripheral nerve", "spinal cord", "dorsal root ganglion",
        "fetal brain", "fetal liver", "fetal kidney",
        "skeletal muscle", "smooth muscle", "cardiac muscle",
        "small intestine", "large intestine", "duodenum", "ileum", "jejunum",
        "appendix", "rectum", "tonsil", "lymphocyte", "monocyte",
        "neutrophil", "platelet", "erythrocyte", "fibroblast",
        "epithelium", "endothelium", "leukocyte",
    ]

    text_lower = tissue_text.lower()
    found: List[str] = []
    for tissue in known_tissues:
        if tissue in text_lower:
            found.append(tissue)

    return found if found else None


def _parse_features(features: List[Dict]) -> Dict[str, Any]:
    """Extract structured data from UniProt features."""
    result: Dict[str, Any] = {}

    # -- Signal peptide --
    signal_peptides = [f for f in features if f.get("type") == "Signal"]
    if signal_peptides:
        sp = signal_peptides[0]
        loc = sp.get("location", {})
        start = loc.get("start", {}).get("value", "?")
        end = loc.get("end", {}).get("value", "?")
        result["signal_peptide"] = f"{start}-{end}"
    else:
        result["signal_peptide"] = None

    # -- Transmembrane regions --
    tm_regions: List[str] = []
    for f in features:
        if f.get("type") == "Transmembrane":
            loc = f.get("location", {})
            start = loc.get("start", {}).get("value", "?")
            end = loc.get("end", {}).get("value", "?")
            desc = f.get("description", "")
            label = f"{start}-{end}"
            if desc:
                label += f" ({desc})"
            tm_regions.append(label)
    result["transmembrane_regions"] = tm_regions if tm_regions else None

    # -- Active site --
    active_sites: List[str] = []
    for f in features:
        if f.get("type") == "Active site":
            loc = f.get("location", {})
            pos = loc.get("start", {}).get("value", "?")
            desc = f.get("description", "")
            label = f"pos {pos}"
            if desc:
                label += f" ({desc})"
            active_sites.append(label)
    result["active_site"] = "; ".join(active_sites)[:1000] if active_sites else None

    # -- Binding sites --
    binding_sites: List[str] = []
    for f in features:
        if f.get("type") == "Binding site":
            loc = f.get("location", {})
            start = loc.get("start", {}).get("value", "?")
            end = loc.get("end", {}).get("value", "?")
            desc = f.get("description", "")
            # Get ligand from cross-refs
            ligands = []
            for xref in f.get("featureCrossReferences", []):
                if xref.get("database") == "ChEBI":
                    ligands.append(xref.get("id", ""))
            label = f"{start}-{end}"
            if desc:
                label += f" ({desc})"
            if ligands:
                label += f" [ChEBI:{','.join(ligands)}]"
            binding_sites.append(label)
    result["binding_sites"] = binding_sites if binding_sites else None

    # -- Disulfide bonds (count) --
    disulfides = [f for f in features if f.get("type") == "Disulfide bond"]
    result["disulfide_bonds"] = len(disulfides) if disulfides else None

    # -- Glycosylation sites (count) --
    glyco = [f for f in features if f.get("type") == "Glycosylation"]
    result["glycosylation_sites"] = len(glyco) if glyco else None

    # -- Phosphorylation sites (count) --
    # Phosphorylation is under "Modified residue" with "Phospho" in description
    mod_residues = [f for f in features if f.get("type") == "Modified residue"]
    phospho_count = sum(
        1 for f in mod_residues
        if "phospho" in f.get("description", "").lower()
    )
    result["phosphorylation_sites"] = phospho_count if phospho_count > 0 else None

    # -- Domains --
    domains: List[str] = []
    for f in features:
        if f.get("type") == "Domain":
            desc = f.get("description", "").strip()
            loc = f.get("location", {})
            start = loc.get("start", {}).get("value", "?")
            end = loc.get("end", {}).get("value", "?")
            if desc:
                domains.append(f"{desc} ({start}-{end})")
            else:
                domains.append(f"domain ({start}-{end})")
    result["domains"] = domains if domains else None

    return result


def _parse_xrefs(xrefs: List[Dict]) -> Dict[str, Any]:
    """Extract structured data from UniProt cross-references."""
    result: Dict[str, Any] = {}

    go_mf: List[str] = []
    go_bp: List[str] = []
    go_cc: List[str] = []
    pdb_ids: List[str] = []
    kegg_pathways: List[str] = []
    reactome_pathways: List[str] = []
    pfam_ids: List[str] = []
    interpro_ids: List[str] = []
    omim_ids: List[str] = []
    orphanet_ids: List[str] = []
    pharmgkb_ids: List[str] = []
    intact_ids: List[str] = []
    string_id: Optional[str] = None
    ncbi_gene_id: Optional[int] = None
    chromosome: Optional[str] = None

    for xref in xrefs:
        db = xref.get("database", "")
        xid = xref.get("id", "")
        props = {p.get("key"): p.get("value") for p in xref.get("properties", [])}

        if db == "GO":
            term = props.get("GoTerm", "")
            if term.startswith("F:"):
                go_mf.append(term[2:])
            elif term.startswith("P:"):
                go_bp.append(term[2:])
            elif term.startswith("C:"):
                go_cc.append(term[2:])

        elif db == "PDB":
            pdb_ids.append(xid)

        elif db == "KEGG":
            kegg_pathways.append(xid)

        elif db == "Reactome":
            pathway_name = props.get("PathwayName", "")
            if pathway_name:
                reactome_pathways.append(f"{xid}:{pathway_name}")
            else:
                reactome_pathways.append(xid)

        elif db == "Pfam":
            entry_name = props.get("EntryName", "")
            if entry_name:
                pfam_ids.append(f"{xid}:{entry_name}")
            else:
                pfam_ids.append(xid)

        elif db == "InterPro":
            entry_name = props.get("EntryName", "")
            if entry_name:
                interpro_ids.append(f"{xid}:{entry_name}")
            else:
                interpro_ids.append(xid)

        elif db == "MIM":
            omim_ids.append(xid)

        elif db == "Orphanet":
            disease = props.get("Disease", "")
            if disease:
                orphanet_ids.append(f"{xid}:{disease}")
            else:
                orphanet_ids.append(xid)

        elif db == "PharmGKB":
            pharmgkb_ids.append(xid)

        elif db == "STRING":
            string_id = xid

        elif db == "IntAct":
            intact_ids.append(xid)

        elif db == "GeneID":
            try:
                ncbi_gene_id = int(xid)
            except (ValueError, TypeError):
                pass

        elif db == "Ensembl":
            # Try to extract chromosome from Ensembl gene ID pattern
            # Ensembl IDs don't directly encode chromosome, but we store
            # the gene ID for reference
            pass

        elif db == "HGNC":
            # Could use HGNC for chromosome info in future
            pass

    result["go_molecular_function"] = go_mf if go_mf else None
    result["go_biological_process"] = go_bp if go_bp else None
    result["go_cellular_component"] = go_cc if go_cc else None
    result["pdb_ids"] = pdb_ids if pdb_ids else None
    result["pathway_kegg"] = kegg_pathways if kegg_pathways else None
    result["pathway_reactome"] = [r for r in reactome_pathways] if reactome_pathways else None
    result["pfam_ids"] = pfam_ids if pfam_ids else None
    result["interpro_ids"] = interpro_ids if interpro_ids else None
    result["cross_refs_omim"] = omim_ids if omim_ids else None
    result["cross_refs_orphanet"] = orphanet_ids if orphanet_ids else None
    result["cross_refs_pharmgkb"] = pharmgkb_ids if pharmgkb_ids else None
    result["cross_refs_string"] = string_id
    result["cross_refs_intact"] = intact_ids if intact_ids else None
    result["cross_refs_reactome"] = [r for r in reactome_pathways] if reactome_pathways else None
    result["ncbi_gene_id"] = ncbi_gene_id
    result["chromosome"] = chromosome
    result["gene_location"] = None  # Would need Ensembl API for precise location

    return result


def _extract_protein_family(keywords: List[Dict]) -> Optional[str]:
    """Determine protein family from UniProt keywords."""
    family_keywords = {
        "Kinase": "Kinase",
        "Tyrosine-protein kinase": "Kinase",
        "Serine/threonine-protein kinase": "Kinase",
        "Transferase": "Transferase",
        "Hydrolase": "Hydrolase",
        "Receptor": "Receptor",
        "G-protein coupled receptor": "GPCR",
        "Ion channel": "Ion channel",
        "Transporter": "Transporter",
        "Protease": "Protease",
        "Oxidoreductase": "Oxidoreductase",
        "Ligase": "Ligase",
        "Lyase": "Lyase",
        "Isomerase": "Isomerase",
        "Phosphatase": "Phosphatase",
        "Metalloprotease": "Protease",
        "Serine protease": "Protease",
        "Cysteine protease": "Protease",
        "Aspartyl protease": "Protease",
        "Nuclear receptor": "Nuclear receptor",
        "Transcription factor": "Transcription factor",
        "Chaperone": "Chaperone",
        "Deubiquitinase": "Deubiquitinase",
        "Ubiquitin ligase": "Ubiquitin ligase",
    }

    # Higher priority families first
    priority = [
        "GPCR", "Kinase", "Ion channel", "Nuclear receptor",
        "Protease", "Phosphatase", "Transporter",
        "Receptor", "Transferase", "Hydrolase", "Oxidoreductase",
        "Ligase", "Lyase", "Isomerase",
        "Transcription factor", "Chaperone",
        "Deubiquitinase", "Ubiquitin ligase",
    ]

    kw_names = {kw.get("name") for kw in keywords if kw.get("name")}
    matched_families = set()

    for kw_name in kw_names:
        if kw_name in family_keywords:
            matched_families.add(family_keywords[kw_name])

    if not matched_families:
        return None

    # Return highest priority match
    for fam in priority:
        if fam in matched_families:
            return fam

    return matched_families.pop()


# ---------------------------------------------------------------------------
# Database helpers
# ---------------------------------------------------------------------------

def get_conn() -> psycopg2.extensions.connection:
    """Get a database connection."""
    return psycopg2.connect(**DB_CONFIG)


def start_log(conn: psycopg2.extensions.connection, source: str, step: str) -> int:
    """Create an ingestion log entry, return log ID."""
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
    """Update an ingestion log entry."""
    with conn.cursor() as cur:
        cur.execute("""
            UPDATE pharmaco_db.ingestion_log
            SET status=%s, rows_inserted=%s, error_message=%s, completed_at=NOW()
            WHERE id=%s
        """, (status, rows, error, log_id))
    conn.commit()


def get_targets_to_enrich(
    conn: psycopg2.extensions.connection,
    force: bool = False,
    limit: Optional[int] = None,
) -> List[Tuple[int, str, Optional[str]]]:
    """
    Fetch targets that need UniProt enrichment.
    Returns list of (target_id, uniprot_id, chembl_id).
    """
    where_clause = "WHERE t.uniprot_id IS NOT NULL"
    if not force:
        # Only enrich targets that haven't been processed yet.
        # We use uniprot_enriched_at (not uniprot_function) because many
        # valid proteins lack a FUNCTION comment (unreviewed TrEMBL entries).
        where_clause += " AND t.uniprot_enriched_at IS NULL"

    query = f"""
        SELECT t.id, t.uniprot_id, t.chembl_id
        FROM pharmaco_db.targets t
        {where_clause}
        ORDER BY t.id
    """
    if limit:
        query += f" LIMIT {limit}"

    with conn.cursor() as cur:
        cur.execute(query)
        return cur.fetchall()


# ---------------------------------------------------------------------------
# Main UPDATE query builder
# ---------------------------------------------------------------------------

# Columns that we update in the targets table
# Maps: column_name -> (sql_type_hint, is_array)
UPDATE_COLUMNS: Dict[str, Tuple[str, bool]] = {
    "gene_name":              ("text",     False),
    "uniprot_function":       ("text",     False),
    "uniprot_subcellular":    ("text",     False),
    "go_molecular_function":  ("text[]",   True),
    "go_biological_process":  ("text[]",   True),
    "go_cellular_component":  ("text[]",   True),
    "pdb_ids":                ("text[]",   True),
    "pathway_kegg":           ("text[]",   True),
    "pathway_reactome":       ("text[]",   True),
    "ncbi_gene_id":           ("integer",  False),
    "protein_family":         ("text",     False),
    "sequence_length":        ("integer",  False),
    "mass":                   ("integer",  False),
    "ec_number":              ("text",     False),
    "signal_peptide":         ("text",     False),
    "transmembrane_regions":  ("text[]",   True),
    "active_site":            ("text",     False),
    "binding_sites":          ("text[]",   True),
    "disulfide_bonds":        ("integer",  False),
    "glycosylation_sites":    ("integer",  False),
    "phosphorylation_sites":  ("integer",  False),
    "domains":                ("text[]",   True),
    "interpro_ids":           ("text[]",   True),
    "pfam_ids":               ("text[]",   True),
    "tissue_expression":      ("text[]",   True),
    "tissue_specificity":     ("text",     False),
    "disease_associations":   ("text[]",   True),
    "involvement_in_disease": ("text",     False),
    "keywords":               ("text[]",   True),
    "protein_existence":      ("text",     False),
    "cross_refs_omim":        ("text[]",   True),
    "cross_refs_orphanet":    ("text[]",   True),
    "cross_refs_pharmgkb":    ("text[]",   True),
    "cross_refs_reactome":    ("text[]",   True),
    "cross_refs_string":      ("text",     False),
    "cross_refs_intact":      ("text[]",   True),
    "isoforms":               ("integer",  False),
    "alternative_names":      ("text[]",   True),
    "chromosome":             ("text",     False),
    "gene_location":          ("text",     False),
}


def build_update_query() -> str:
    """Build the parameterized UPDATE query for all enrichment columns."""
    set_clauses: List[str] = []
    for col, (_, _) in UPDATE_COLUMNS.items():
        # Use COALESCE for columns that may already have data (gene_name, protein_family)
        if col in ("gene_name", "protein_family"):
            set_clauses.append(f"{col} = COALESCE(%s, {col})")
        else:
            set_clauses.append(f"{col} = %s")

    set_clauses.append("updated_at = NOW()")
    set_clauses.append("uniprot_enriched_at = NOW()")

    return f"""
        UPDATE pharmaco_db.targets SET
            {', '.join(set_clauses)}
        WHERE uniprot_id = %s
    """


def build_update_params(parsed: Dict[str, Any], accession: str) -> Tuple:
    """Build the parameter tuple for the UPDATE query."""
    params = []
    for col, (_, _) in UPDATE_COLUMNS.items():
        value = parsed.get(col)
        params.append(value)
    params.append(accession)
    return tuple(params)


# ---------------------------------------------------------------------------
# Cross-references insertion
# ---------------------------------------------------------------------------

def insert_cross_references(
    conn: psycopg2.extensions.connection,
    target_id: int,
    parsed: Dict[str, Any],
) -> None:
    """Insert cross-references into the cross_references table."""
    xref_batch: List[Tuple] = []

    # PDB structures
    pdb_ids = parsed.get("pdb_ids")
    if pdb_ids:
        for pdb_id in pdb_ids[:20]:  # Limit to 20
            xref_batch.append((
                None, target_id, "PDB", pdb_id,
                f"https://www.rcsb.org/structure/{pdb_id}",
            ))

    # NCBI Gene
    ncbi_gene_id = parsed.get("ncbi_gene_id")
    if ncbi_gene_id:
        xref_batch.append((
            None, target_id, "NCBI_Gene", str(ncbi_gene_id),
            f"https://www.ncbi.nlm.nih.gov/gene/{ncbi_gene_id}",
        ))

    # STRING
    string_id = parsed.get("cross_refs_string")
    if string_id:
        xref_batch.append((
            None, target_id, "STRING", string_id,
            f"https://string-db.org/network/{string_id}",
        ))

    # IntAct
    intact_ids = parsed.get("cross_refs_intact")
    if intact_ids:
        for intact_id in intact_ids[:5]:
            xref_batch.append((
                None, target_id, "IntAct", intact_id,
                f"https://www.ebi.ac.uk/intact/query/{intact_id}",
            ))

    # OMIM
    omim_ids = parsed.get("cross_refs_omim")
    if omim_ids:
        for omim_id in omim_ids[:10]:
            xref_batch.append((
                None, target_id, "OMIM", omim_id,
                f"https://omim.org/entry/{omim_id}",
            ))

    # Orphanet
    orphanet_ids = parsed.get("cross_refs_orphanet")
    if orphanet_ids:
        for oid in orphanet_ids[:10]:
            raw_id = oid.split(":")[0] if ":" in oid else oid
            xref_batch.append((
                None, target_id, "Orphanet", raw_id,
                f"https://www.orpha.net/consor/cgi-bin/OC_Exp.php?lng=en&Expert={raw_id}",
            ))

    if xref_batch:
        with conn.cursor() as cur:
            psycopg2.extras.execute_values(cur, """
                INSERT INTO pharmaco_db.cross_references
                (compound_id, target_id, db_name, db_id, url)
                VALUES %s
                ON CONFLICT DO NOTHING
            """, xref_batch, page_size=100)


# ---------------------------------------------------------------------------
# Main enrichment loop
# ---------------------------------------------------------------------------

def enrich_targets(
    conn: psycopg2.extensions.connection,
    force: bool = False,
    limit: Optional[int] = None,
    batch_size: int = BATCH_SIZE,
) -> int:
    """
    Main enrichment function.
    Fetches UniProt data in batches and updates the targets table.
    Returns number of successfully enriched targets.
    """
    log("=" * 70)
    log("UniProt Comprehensive Target Enrichment (V2)")
    log("=" * 70)

    log_id = start_log(conn, "uniprot_full", "comprehensive_enrichment")

    # Get targets to enrich
    targets = get_targets_to_enrich(conn, force=force, limit=limit)
    total = len(targets)
    log(f"  Found {total} targets to enrich (force={force}, limit={limit})")

    if total == 0:
        log("  Nothing to do.")
        update_log(conn, log_id, "completed", 0)
        return 0

    # Build target lookup: uniprot_id -> (target_id, chembl_id)
    target_lookup: Dict[str, Tuple[int, Optional[str]]] = {}
    all_accessions: List[str] = []
    for target_id, uniprot_id, chembl_id in targets:
        target_lookup[uniprot_id] = (target_id, chembl_id)
        all_accessions.append(uniprot_id)

    # Pre-build the UPDATE query
    update_query = build_update_query()

    enriched = 0
    errors = 0
    skipped = 0
    batch_num = 0
    start_time = time.time()

    # Pre-fetch batches concurrently for speed
    from concurrent.futures import ThreadPoolExecutor, as_completed
    import threading

    # Build list of batch accession groups
    batch_groups = []
    for i in range(0, total, batch_size):
        batch_groups.append((i, all_accessions[i : i + batch_size]))

    # Prefetch queue: fetch PARALLEL_WORKERS batches ahead
    prefetch_cache: Dict[int, List[Dict]] = {}
    cache_lock = threading.Lock()

    def prefetch_batch(idx_accs):
        idx, accs = idx_accs
        entries = fetch_batch(accs)
        with cache_lock:
            prefetch_cache[idx] = entries or []
        return idx

    # Process with prefetch pool
    fetch_pool = ThreadPoolExecutor(max_workers=PARALLEL_WORKERS)
    # Submit first N batches
    pending_futures = {}
    submit_ptr = 0
    for _ in range(min(PARALLEL_WORKERS * 2, len(batch_groups))):
        if submit_ptr < len(batch_groups):
            f = fetch_pool.submit(prefetch_batch, batch_groups[submit_ptr])
            pending_futures[submit_ptr] = f
            submit_ptr += 1

    for batch_idx, (i, batch_accessions) in enumerate(batch_groups):
        batch_num += 1

        # Progress logging
        if i > 0 and i % 200 == 0:
            elapsed = time.time() - start_time
            rate = enriched / elapsed if elapsed > 0 else 0
            eta_seconds = (total - i) / rate if rate > 0 else 0
            eta_min = int(eta_seconds // 60)
            eta_sec = int(eta_seconds % 60)
            log(
                f"  Progress: {i}/{total} "
                f"({enriched} enriched, {errors} errors, {skipped} skipped) "
                f"[{rate:.1f}/s, ETA {eta_min}m{eta_sec:02d}s]"
            )

        # Wait for this batch's prefetch to complete
        if batch_idx in pending_futures:
            pending_futures[batch_idx].result()
            del pending_futures[batch_idx]

        # Submit next batch for prefetch (with slight delay)
        if submit_ptr < len(batch_groups):
            time.sleep(RATE_LIMIT_SECONDS)
            f = fetch_pool.submit(prefetch_batch, batch_groups[submit_ptr])
            pending_futures[submit_ptr] = f
            submit_ptr += 1

        # Get entries from prefetch cache
        with cache_lock:
            entries = prefetch_cache.pop(batch_idx, [])

        if not entries:
            errors += len(batch_accessions)
            log_error(f"  Batch {batch_num} returned 0 entries for {len(batch_accessions)} accessions")
            continue

        # Index entries by accession for lookup
        entry_map: Dict[str, Dict] = {}
        for entry in entries:
            acc = entry.get("primaryAccession", "")
            if acc:
                entry_map[acc] = entry
            for sec_acc in entry.get("secondaryAccessions", []):
                if sec_acc not in entry_map:
                    entry_map[sec_acc] = entry
            queried = entry.get("_queried_accession", "")
            if queried and queried not in entry_map:
                entry_map[queried] = entry

        # Process each target in the batch
        for accession in batch_accessions:
            entry = entry_map.get(accession)
            if not entry:
                skipped += 1
                continue

            try:
                parsed = parse_entry(entry)
                params = build_update_params(parsed, accession)

                with conn.cursor() as cur:
                    cur.execute(update_query, params)

                target_id, _ = target_lookup[accession]
                insert_cross_references(conn, target_id, parsed)

                conn.commit()
                enriched += 1

            except Exception as e:
                conn.rollback()
                errors += 1
                log_error(f"  Failed to process {accession}: {e}")
                if errors <= 5:
                    traceback.print_exc()

    fetch_pool.shutdown(wait=True)

    # Final stats
    elapsed = time.time() - start_time
    elapsed_min = int(elapsed // 60)
    elapsed_sec = int(elapsed % 60)

    log("=" * 70)
    log(f"  COMPLETED in {elapsed_min}m {elapsed_sec}s")
    log(f"  Total targets:   {total}")
    log(f"  Enriched:        {enriched}")
    log(f"  Errors:          {errors}")
    log(f"  Skipped (404):   {skipped}")
    log(f"  Rate:            {enriched / elapsed:.1f} targets/sec" if elapsed > 0 else "  Rate: N/A")
    log("=" * 70)

    update_log(conn, log_id, "completed", enriched)
    return enriched


# ---------------------------------------------------------------------------
# Post-enrichment statistics
# ---------------------------------------------------------------------------

def print_enrichment_stats(conn: psycopg2.extensions.connection) -> None:
    """Print summary statistics after enrichment."""
    log("")
    log("=" * 70)
    log("POST-ENRICHMENT STATISTICS")
    log("=" * 70)

    queries = [
        ("Total targets", "SELECT COUNT(*) FROM pharmaco_db.targets"),
        ("With UniProt ID", "SELECT COUNT(*) FROM pharmaco_db.targets WHERE uniprot_id IS NOT NULL"),
        ("With function", "SELECT COUNT(*) FROM pharmaco_db.targets WHERE uniprot_function IS NOT NULL"),
        ("With subcellular", "SELECT COUNT(*) FROM pharmaco_db.targets WHERE uniprot_subcellular IS NOT NULL"),
        ("With GO terms (MF)", "SELECT COUNT(*) FROM pharmaco_db.targets WHERE go_molecular_function IS NOT NULL"),
        ("With GO terms (BP)", "SELECT COUNT(*) FROM pharmaco_db.targets WHERE go_biological_process IS NOT NULL"),
        ("With PDB structures", "SELECT COUNT(*) FROM pharmaco_db.targets WHERE pdb_ids IS NOT NULL"),
        ("With KEGG pathways", "SELECT COUNT(*) FROM pharmaco_db.targets WHERE pathway_kegg IS NOT NULL"),
        ("With Reactome pathways", "SELECT COUNT(*) FROM pharmaco_db.targets WHERE pathway_reactome IS NOT NULL"),
        ("With EC number", "SELECT COUNT(*) FROM pharmaco_db.targets WHERE ec_number IS NOT NULL"),
        ("With domains", "SELECT COUNT(*) FROM pharmaco_db.targets WHERE domains IS NOT NULL"),
        ("With disease assoc.", "SELECT COUNT(*) FROM pharmaco_db.targets WHERE disease_associations IS NOT NULL"),
        ("With tissue specificity", "SELECT COUNT(*) FROM pharmaco_db.targets WHERE tissue_specificity IS NOT NULL"),
        ("With keywords", "SELECT COUNT(*) FROM pharmaco_db.targets WHERE keywords IS NOT NULL"),
        ("With Pfam IDs", "SELECT COUNT(*) FROM pharmaco_db.targets WHERE pfam_ids IS NOT NULL"),
        ("With InterPro IDs", "SELECT COUNT(*) FROM pharmaco_db.targets WHERE interpro_ids IS NOT NULL"),
        ("With OMIM refs", "SELECT COUNT(*) FROM pharmaco_db.targets WHERE cross_refs_omim IS NOT NULL"),
        ("With STRING ref", "SELECT COUNT(*) FROM pharmaco_db.targets WHERE cross_refs_string IS NOT NULL"),
        ("With signal peptide", "SELECT COUNT(*) FROM pharmaco_db.targets WHERE signal_peptide IS NOT NULL"),
        ("With transmembrane", "SELECT COUNT(*) FROM pharmaco_db.targets WHERE transmembrane_regions IS NOT NULL"),
        ("With active site", "SELECT COUNT(*) FROM pharmaco_db.targets WHERE active_site IS NOT NULL"),
        ("With binding sites", "SELECT COUNT(*) FROM pharmaco_db.targets WHERE binding_sites IS NOT NULL"),
        ("With isoforms", "SELECT COUNT(*) FROM pharmaco_db.targets WHERE isoforms IS NOT NULL AND isoforms > 0"),
    ]

    with conn.cursor() as cur:
        for label, query in queries:
            cur.execute(query)
            count = cur.fetchone()[0]
            log(f"  {label:<30} {count:>7,}")

    # Top protein families
    log("")
    log("  Top protein families:")
    with conn.cursor() as cur:
        cur.execute("""
            SELECT protein_family, COUNT(*) as cnt
            FROM pharmaco_db.targets
            WHERE protein_family IS NOT NULL
            GROUP BY protein_family
            ORDER BY cnt DESC
            LIMIT 15
        """)
        for family, count in cur.fetchall():
            log(f"    {family:<30} {count:>5,}")

    # Cross-reference counts
    log("")
    log("  Cross-references by database:")
    with conn.cursor() as cur:
        cur.execute("""
            SELECT db_name, COUNT(*) as cnt
            FROM pharmaco_db.cross_references
            WHERE target_id IS NOT NULL
            GROUP BY db_name
            ORDER BY cnt DESC
        """)
        for db_name, count in cur.fetchall():
            log(f"    {db_name:<20} {count:>8,}")

    log("=" * 70)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    """Main entry point with CLI argument parsing."""
    parser = argparse.ArgumentParser(
        description="PharmacoDB -- Comprehensive UniProt Target Enrichment",
    )
    parser.add_argument(
        "--limit", type=int, default=None,
        help="Limit number of targets to enrich (for testing)",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-enrich ALL targets (even those already processed via uniprot_enriched_at)",
    )
    parser.add_argument(
        "--batch-size", type=int, default=BATCH_SIZE,
        help=f"Number of accessions per batch request (default: {BATCH_SIZE})",
    )
    parser.add_argument(
        "--stats-only", action="store_true",
        help="Only print enrichment statistics, do not fetch data",
    )
    args = parser.parse_args()

    log("=" * 70)
    log("PharmacoDB -- Comprehensive UniProt Enrichment Pipeline (V2)")
    log(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log("=" * 70)

    conn = get_conn()
    try:
        if args.stats_only:
            print_enrichment_stats(conn)
            return

        enriched = enrich_targets(
            conn,
            force=args.force,
            limit=args.limit,
            batch_size=args.batch_size,
        )

        if enriched > 0:
            print_enrichment_stats(conn)

        log("")
        log("DONE.")

    except KeyboardInterrupt:
        log("\nInterrupted by user.")
        sys.exit(1)
    except Exception as e:
        log_error(f"FATAL: {e}")
        traceback.print_exc()
        sys.exit(1)
    finally:
        conn.close()


if __name__ == "__main__":
    main()
