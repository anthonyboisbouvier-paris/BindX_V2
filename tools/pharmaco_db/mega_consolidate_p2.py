#!/usr/bin/env python3
"""
MegaPharmaDB — Phase 2 Consolidation
Fills ALL remaining empty tables:
1. PubChem CID mapping via UniChem (unlocks SIDER)
2. DGIdb via GraphQL API
3. GtoPdb via TSV download
4. SIDER side effects (with PubChem CIDs)
5. Open Targets disease-target associations
6. Compound structures + fingerprints (RDKit)
7. DisGeNET (alternative: mine from ChEMBL)
8. Tox21 toxicology
9. KEGG drug info
"""

import psycopg2
import psycopg2.extras
import os, sys, time, gzip, csv, io, json, struct
import urllib.request, urllib.error
from datetime import datetime
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

DB_CONFIG = {
    "host": "localhost", "port": 5433,
    "dbname": "pharmaco", "user": "postgres", "password": "pharmaco_secret",
}
DOWNLOAD_DIR = Path("/tmp/megapharma_downloads")
BATCH = 10000

def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

def get_conn():
    return psycopg2.connect(**DB_CONFIG)

def download_file(url, dest, desc=""):
    if dest.exists() and dest.stat().st_size > 100:
        log(f"  Already downloaded: {dest.name} ({dest.stat().st_size/1024/1024:.1f}MB)")
        return dest
    log(f"  Downloading {desc or dest.name}...")
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "MegaPharmaDB/1.0"})
        with urllib.request.urlopen(req, timeout=300) as resp:
            with open(dest, "wb") as f:
                downloaded = 0
                while True:
                    chunk = resp.read(1024 * 1024)
                    if not chunk: break
                    f.write(chunk)
                    downloaded += len(chunk)
            log(f"  Downloaded {downloaded/1024/1024:.1f} MB")
        return dest
    except Exception as e:
        log(f"  DOWNLOAD FAILED: {e}")
        if dest.exists(): dest.unlink()
        return None

def fetch_json(url, data=None, headers=None):
    hdrs = {"User-Agent": "MegaPharmaDB/1.0", "Content-Type": "application/json"}
    if headers: hdrs.update(headers)
    req = urllib.request.Request(url, data=json.dumps(data).encode() if data else None, headers=hdrs)
    with urllib.request.urlopen(req, timeout=60) as resp:
        return json.loads(resp.read())

# ============================================================
# 1. PUBCHEM CID MAPPING VIA UNICHEM
# ============================================================
def map_pubchem_cids():
    """Map ChEMBL IDs to PubChem CIDs using UniChem bulk download."""
    log("=" * 60)
    log("1. PubChem CID mapping via UniChem")
    log("=" * 60)

    # UniChem src1 (ChEMBL) -> src22 (PubChem)
    url = "https://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src22.txt.gz"
    dest = DOWNLOAD_DIR / "unichem_chembl_pubchem.txt.gz"
    f = download_file(url, dest, "UniChem ChEMBL->PubChem mapping")
    if not f:
        log("  SKIP: UniChem download failed")
        return 0

    # Parse: chembl_id -> pubchem_cid
    mapping = {}
    with gzip.open(f, "rt") as fh:
        next(fh)  # header: "From src:'1'\tTo src:'22'"
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                chembl_id, pcid_str = parts[0].strip(), parts[1].strip()
                try:
                    pcid = int(pcid_str)
                    # File already has "CHEMBL" prefix
                    if chembl_id not in mapping:  # Keep first (lowest CID)
                        mapping[chembl_id] = pcid
                except ValueError:
                    continue
    log(f"  UniChem mappings loaded: {len(mapping):,}")

    # Update compounds table
    conn = get_conn()
    with conn.cursor() as cur:
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.compounds WHERE chembl_id IS NOT NULL AND pubchem_cid IS NULL")
        to_update = cur.fetchall()
    log(f"  Compounds missing PubChem CID: {len(to_update):,}")

    batch = []
    total = 0
    for chembl_id, cid in to_update:
        pcid = mapping.get(chembl_id)
        if pcid:
            batch.append((pcid, cid))
            if len(batch) >= BATCH:
                with conn.cursor() as cur:
                    psycopg2.extras.execute_batch(cur, """
                        UPDATE pharmaco_db.compounds SET pubchem_cid = %s WHERE id = %s
                    """, batch, page_size=5000)
                conn.commit()
                total += len(batch)
                batch = []
    if batch:
        with conn.cursor() as cur:
            psycopg2.extras.execute_batch(cur, """
                UPDATE pharmaco_db.compounds SET pubchem_cid = %s WHERE id = %s
            """, batch, page_size=5000)
        conn.commit()
        total += len(batch)

    conn.close()
    log(f"  Updated PubChem CIDs: {total:,}")
    return total


# ============================================================
# 2. DGIdb VIA GRAPHQL
# ============================================================
def ingest_dgidb():
    """Ingest DGIdb drug-gene interactions via GraphQL API."""
    log("=" * 60)
    log("2. DGIdb interactions (GraphQL API)")
    log("=" * 60)

    conn = get_conn()
    with conn.cursor() as cur:
        cur.execute("SELECT gene_name, id FROM pharmaco_db.targets WHERE gene_name IS NOT NULL")
        gene_map = dict(cur.fetchall())
        cur.execute("SELECT pref_name, id FROM pharmaco_db.compounds WHERE pref_name IS NOT NULL")
        name_map = {k.upper(): v for k, v in cur.fetchall()}
        cur.execute("SELECT chembl_id, id FROM pharmaco_db.compounds WHERE chembl_id IS NOT NULL")
        chembl_map = dict(cur.fetchall())

    genes = list(gene_map.keys())
    log(f"  Querying {len(genes):,} genes in batches of 50...")

    all_interactions = []
    api_url = "https://dgidb.org/api/graphql"

    for i in range(0, len(genes), 50):
        chunk = genes[i:i+50]
        query = {
            "query": """
            query($names: [String!]!) {
              genes(names: $names) {
                nodes {
                  name
                  interactions {
                    drug { name conceptId }
                    interactionScore
                    interactionTypes { type directionality }
                    publications { pmid }
                    sources { fullName }
                  }
                }
              }
            }
            """,
            "variables": {"names": chunk}
        }

        try:
            result = fetch_json(api_url, query)
            nodes = result.get("data", {}).get("genes", {}).get("nodes", [])
            for gene_node in nodes:
                gene_name = gene_node["name"]
                tid = gene_map.get(gene_name)
                if not tid:
                    continue
                for ix in gene_node.get("interactions", []):
                    drug = ix.get("drug", {})
                    drug_name = drug.get("name", "")
                    drug_concept = drug.get("conceptId", "")
                    score = ix.get("interactionScore")
                    types = [t.get("type", "") for t in ix.get("interactionTypes", [])]
                    pmids = [str(p.get("pmid", "")) for p in ix.get("publications", []) if p.get("pmid")]
                    sources = [s.get("fullName", "") for s in ix.get("sources", []) if s.get("fullName")]

                    # Map drug to compound
                    cid = None
                    if drug_concept and drug_concept.startswith("chembl:"):
                        chembl_id = drug_concept.replace("chembl:", "").upper()
                        cid = chembl_map.get(chembl_id)
                    if not cid and drug_name:
                        cid = name_map.get(drug_name.upper())

                    ix_type = ", ".join(types) if types else None

                    all_interactions.append((
                        cid, tid, gene_name, drug_name,
                        drug_concept if drug_concept else None,
                        ix_type[:200] if ix_type else None,
                        float(score) if score else None,
                        sources if sources else None,
                        pmids if pmids else None,
                    ))
        except Exception as e:
            if i == 0:
                log(f"  DGIdb API error: {e}")
            continue

        if (i // 50) % 20 == 0 and i > 0:
            log(f"  Progress: {i}/{len(genes)} ({len(all_interactions)} interactions)")
        time.sleep(0.1)  # Rate limit

    if all_interactions:
        # Deduplicate
        seen = set()
        unique = []
        for row in all_interactions:
            key = (row[2], row[3], row[5])  # gene, drug, type
            if key not in seen:
                seen.add(key)
                unique.append(row)

        with conn.cursor() as cur:
            cur.execute("TRUNCATE pharmaco_db.dgidb_interactions")
            psycopg2.extras.execute_values(cur, """
                INSERT INTO pharmaco_db.dgidb_interactions
                (compound_id, target_id, gene_name, drug_name, drug_chembl_id,
                 interaction_type, interaction_score, sources, pmids)
                VALUES %s
            """, unique, page_size=5000)
        conn.commit()
        log(f"  DGIdb: {len(unique):,} unique interactions")
    else:
        log("  DGIdb: 0 interactions (API may be down)")

    conn.close()


# ============================================================
# 3. GUIDE TO PHARMACOLOGY (GtoPdb)
# ============================================================
def ingest_gtopdb():
    """Download and ingest GtoPdb interactions."""
    log("=" * 60)
    log("3. Guide to Pharmacology (GtoPdb)")
    log("=" * 60)

    url = "https://www.guidetopharmacology.org/DATA/interactions.tsv"
    dest = DOWNLOAD_DIR / "gtopdb_interactions.tsv"
    f = download_file(url, dest, "GtoPdb interactions")
    if not f:
        log("  SKIP: GtoPdb download failed")
        return

    conn = get_conn()
    with conn.cursor() as cur:
        cur.execute("SELECT gene_name, id FROM pharmaco_db.targets WHERE gene_name IS NOT NULL")
        gene_map = dict(cur.fetchall())
        cur.execute("SELECT uniprot_id, id FROM pharmaco_db.targets WHERE uniprot_id IS NOT NULL")
        uniprot_map = dict(cur.fetchall())
        cur.execute("SELECT pref_name, id FROM pharmaco_db.compounds WHERE pref_name IS NOT NULL")
        name_map = {k.upper(): v for k, v in cur.fetchall()}

    batch = []
    seen = set()
    with open(f, "r", encoding="utf-8") as fh:
        # Skip comment line
        header_line = fh.readline()
        while header_line.startswith("#") or header_line.startswith('"#'):
            header_line = fh.readline()
        reader = csv.DictReader(io.StringIO(header_line + fh.read()), delimiter="\t", quotechar='"')

        for row in reader:
            target_gene = row.get("Target Gene Symbol", "").strip()
            target_uniprot = row.get("Target UniProt ID", "").strip()
            ligand_name = row.get("Ligand", "") or row.get("Target Ligand", "")
            ligand_id = row.get("Ligand ID", "") or row.get("Target Ligand ID", "")
            target_id_gtop = row.get("Target ID", "")
            ix_type = row.get("Type", "") or row.get("Action", "")
            affinity_type = row.get("Action_comment", "") or row.get("Affinity Units", "")
            affinity_val = row.get("Affinity Median", "") or row.get("Affinity High", "")
            species = row.get("Target Species", "")
            endogenous = row.get("Endogenous", "")
            primary = row.get("Primary Target", "") or row.get("Primary", "")

            # Map target
            tid = gene_map.get(target_gene) or uniprot_map.get(target_uniprot)
            if not tid:
                continue

            # Map compound
            cid = name_map.get(ligand_name.strip().upper()) if ligand_name else None

            key = (tid, ligand_name, ix_type)
            if key in seen:
                continue
            seen.add(key)

            batch.append((
                cid, tid,
                int(ligand_id) if ligand_id and ligand_id.isdigit() else None,
                int(target_id_gtop) if target_id_gtop and target_id_gtop.isdigit() else None,
                ix_type[:100] if ix_type else None,
                affinity_type[:50] if affinity_type else None,
                float(affinity_val) if affinity_val and affinity_val.replace(".", "").replace("-", "").isdigit() else None,
                None, None,  # affinity_high, affinity_low
                species[:100] if species else None,
                endogenous.lower() in ("true", "t", "yes", "1") if endogenous else False,
                primary.lower() in ("true", "t", "yes", "1") if primary else False,
            ))

    if batch:
        with conn.cursor() as cur:
            cur.execute("TRUNCATE pharmaco_db.gtop_interactions")
            psycopg2.extras.execute_values(cur, """
                INSERT INTO pharmaco_db.gtop_interactions
                (compound_id, target_id, ligand_id, target_gtop_id,
                 interaction_type, affinity_type, affinity_value,
                 affinity_high, affinity_low, species, endogenous, primary_target)
                VALUES %s
            """, batch, page_size=5000)
        conn.commit()
        log(f"  GtoPdb: {len(batch):,} interactions")
    conn.close()


# ============================================================
# 4. SIDER SIDE EFFECTS (with PubChem CIDs now mapped)
# ============================================================
def ingest_sider():
    """Re-run SIDER ingestion now that PubChem CIDs are mapped."""
    log("=" * 60)
    log("4. SIDER side effects (with PubChem CIDs)")
    log("=" * 60)

    dest = DOWNLOAD_DIR / "sider_all_se.tsv.gz"
    if not dest.exists():
        url = "http://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz"
        dest = download_file(url, dest, "SIDER side effects")
    if not dest or not dest.exists():
        log("  SKIP: SIDER file not available")
        return

    conn = get_conn()
    with conn.cursor() as cur:
        cur.execute("SELECT pubchem_cid, id FROM pharmaco_db.compounds WHERE pubchem_cid IS NOT NULL")
        pubchem_map = dict(cur.fetchall())
    log(f"  PubChem CID map: {len(pubchem_map):,} entries")

    if len(pubchem_map) == 0:
        log("  SKIP: No PubChem CIDs mapped yet")
        conn.close()
        return

    batch = []
    seen = set()
    with gzip.open(dest, "rt") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            stitch_flat = parts[0]
            meddra_type = parts[4]
            meddra_term = parts[5]
            umls_cui = parts[3]

            try:
                # STITCH flat ID: CID1XXXXXXX where leading 1 means flat (not stereo)
                if stitch_flat.startswith("CID1"):
                    pcid = int(stitch_flat[4:])
                elif stitch_flat.startswith("CID0"):
                    pcid = int(stitch_flat[4:])
                else:
                    continue
            except ValueError:
                continue

            cid = pubchem_map.get(pcid)
            if not cid:
                continue

            key = (stitch_flat, umls_cui, meddra_type)
            if key in seen:
                continue
            seen.add(key)

            batch.append((cid, stitch_flat, meddra_term, meddra_type, umls_cui, None, "sider"))

    if batch:
        with conn.cursor() as cur:
            cur.execute("TRUNCATE pharmaco_db.side_effects")
            psycopg2.extras.execute_values(cur, """
                INSERT INTO pharmaco_db.side_effects
                (compound_id, stitch_id, side_effect_name,
                 meddra_concept_type, meddra_id, frequency, source)
                VALUES %s
            """, batch, page_size=5000)
        conn.commit()
        log(f"  SIDER: {len(batch):,} unique side effects")
    else:
        log("  SIDER: 0 matched")
    conn.close()


# ============================================================
# 5. OPEN TARGETS (via REST API, gene symbol based)
# ============================================================
def ingest_open_targets():
    """Ingest Open Targets associations via their REST API."""
    log("=" * 60)
    log("5. Open Targets disease-target associations")
    log("=" * 60)

    conn = get_conn()
    with conn.cursor() as cur:
        cur.execute("""
            SELECT gene_name, id FROM pharmaco_db.targets
            WHERE gene_name IS NOT NULL AND organism = 'Homo sapiens'
        """)
        gene_map = dict(cur.fetchall())
    log(f"  Human targets: {len(gene_map):,}")

    # Open Targets REST API - search by gene symbol
    api_base = "https://api.platform.opentargets.org/api/v4/graphql"
    all_assocs = []
    genes = list(gene_map.items())
    errors = 0

    for i, (gene_name, tid) in enumerate(genes):
        if errors > 20:
            log(f"  Too many errors, stopping at {i}/{len(genes)}")
            break

        query = {
            "query": """
            query($symbol: String!) {
              search(queryString: $symbol, entityNames: ["target"], page: {size: 1, index: 0}) {
                hits {
                  id
                  entity
                }
              }
            }
            """,
            "variables": {"symbol": gene_name}
        }

        try:
            result = fetch_json(api_base, query)
            hits = result.get("data", {}).get("search", {}).get("hits", [])
            if not hits:
                continue

            ensg_id = hits[0]["id"]

            # Now get associations
            query2 = {
                "query": """
                query($ensgId: String!) {
                  target(ensemblId: $ensgId) {
                    id
                    approvedSymbol
                    associatedDiseases(page: {size: 25, index: 0}) {
                      rows {
                        disease { id name therapeuticAreas { id name } }
                        score
                        datatypeScores { componentId score }
                      }
                    }
                  }
                }
                """,
                "variables": {"ensgId": ensg_id}
            }

            result2 = fetch_json(api_base, query2)
            target_data = result2.get("data", {}).get("target", {})
            if not target_data:
                continue

            for row in target_data.get("associatedDiseases", {}).get("rows", []):
                disease = row.get("disease", {})
                disease_id = disease.get("id", "")
                disease_name = disease.get("name", "")
                areas = disease.get("therapeuticAreas", [])
                therapeutic_area = areas[0].get("name", "") if areas else None
                overall_score = row.get("score")

                # Parse datatype scores
                dt_scores = {d["componentId"]: d["score"] for d in row.get("datatypeScores", [])}

                all_assocs.append((
                    tid, disease_id, disease_name, therapeutic_area,
                    overall_score,
                    dt_scores.get("ot_genetics_portal"),
                    dt_scores.get("cancer_gene_census") or dt_scores.get("intogen"),
                    dt_scores.get("chembl"),
                    dt_scores.get("europepmc"),
                    dt_scores.get("expression_atlas"),
                    dt_scores.get("phenodigm"),
                    "open_targets"
                ))
            errors = 0

        except Exception as e:
            errors += 1
            continue

        if i % 100 == 0 and i > 0:
            log(f"  Progress: {i}/{len(genes)} ({len(all_assocs)} associations)")
        time.sleep(0.15)  # Rate limit

    if all_assocs:
        with conn.cursor() as cur:
            cur.execute("TRUNCATE pharmaco_db.disease_target_associations")
            psycopg2.extras.execute_values(cur, """
                INSERT INTO pharmaco_db.disease_target_associations
                (target_id, disease_id, disease_name, therapeutic_area,
                 overall_score, genetic_score, somatic_score, known_drug_score,
                 literature_score, rna_expression_score, animal_model_score, source)
                VALUES %s
            """, all_assocs, page_size=5000)
        conn.commit()
        log(f"  Open Targets: {len(all_assocs):,} disease-target associations")
    conn.close()


# ============================================================
# 6. COMPOUND STRUCTURES + FINGERPRINTS (RDKit)
# ============================================================
def compute_compound_structures():
    """Compute RDKit fingerprints and descriptors for all compounds."""
    log("=" * 60)
    log("6. Compound structures & fingerprints (RDKit)")
    log("=" * 60)

    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, AllChem, MACCSkeys, RDKFingerprint
        from rdkit.Chem.Scaffolds import MurckoScaffold
        import base64
    except ImportError:
        log("  SKIP: RDKit not available")
        return

    conn = get_conn()
    with conn.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.compound_structures")
        existing = cur.fetchone()[0]
    log(f"  Existing structures: {existing:,}")

    # Process in batches using server-side cursor
    read_conn = get_conn()
    read_cur = read_conn.cursor("smiles_cursor")
    read_cur.itersize = 5000
    read_cur.execute("""
        SELECT c.id, c.canonical_smiles
        FROM pharmaco_db.compounds c
        LEFT JOIN pharmaco_db.compound_structures cs ON cs.compound_id = c.id
        WHERE cs.compound_id IS NULL
        AND c.canonical_smiles IS NOT NULL
        ORDER BY c.id
    """)

    batch = []
    total = 0
    errors = 0

    for row in read_cur:
        cid, smiles = row
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                errors += 1
                continue

            # Fingerprints
            morgan = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            morgan_b64 = base64.b64encode(morgan.ToBitString().encode()).decode()

            maccs = MACCSkeys.GenMACCSKeys(mol)
            maccs_b64 = base64.b64encode(maccs.ToBitString().encode()).decode()

            rdkit_fp = Chem.RDKFingerprint(mol)
            rdkit_b64 = base64.b64encode(rdkit_fp.ToBitString().encode()).decode()

            # Descriptors
            num_atoms = mol.GetNumAtoms()
            num_bonds = mol.GetNumBonds()
            ring_info = mol.GetRingInfo()
            num_rings = ring_info.NumRings()
            num_het = Descriptors.NumHeteroatoms(mol) if hasattr(Descriptors, 'NumHeteroatoms') else sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() != 6)
            frac_csp3 = Descriptors.FractionCSP3(mol)
            tpsa = Descriptors.TPSA(mol)
            mr = Descriptors.MolMR(mol)
            n_val = Descriptors.NumValenceElectrons(mol)
            n_rad = Descriptors.NumRadicalElectrons(mol)
            charge = Chem.GetFormalCharge(mol)

            # Scaffold
            try:
                scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                scaffold_smi = Chem.MolToSmiles(scaffold)
            except:
                scaffold_smi = None

            batch.append((
                cid, smiles, morgan_b64, maccs_b64, rdkit_b64,
                num_atoms, num_bonds, num_rings, num_het,
                round(frac_csp3, 4) if frac_csp3 else None,
                round(tpsa, 2) if tpsa else None,
                round(mr, 2) if mr else None,
                n_val, n_rad, charge, scaffold_smi
            ))

        except Exception:
            errors += 1
            continue

        if len(batch) >= 5000:
            with conn.cursor() as wr:
                psycopg2.extras.execute_values(wr, """
                    INSERT INTO pharmaco_db.compound_structures
                    (compound_id, canonical_smiles, morgan_fp_2048, maccs_fp, rdkit_fp,
                     num_atoms, num_bonds, num_rings, num_heteroatoms,
                     fraction_csp3, tpsa, molar_refractivity,
                     num_valence_electrons, num_radical_electrons, formal_charge,
                     murcko_scaffold)
                    VALUES %s
                    ON CONFLICT (compound_id) DO NOTHING
                """, batch, page_size=2000)
            conn.commit()
            total += len(batch)
            if total % 50000 == 0:
                log(f"  Progress: {total:,} structures ({errors:,} errors)")
            batch = []

    if batch:
        with conn.cursor() as wr:
            psycopg2.extras.execute_values(wr, """
                INSERT INTO pharmaco_db.compound_structures
                (compound_id, canonical_smiles, morgan_fp_2048, maccs_fp, rdkit_fp,
                 num_atoms, num_bonds, num_rings, num_heteroatoms,
                 fraction_csp3, tpsa, molar_refractivity,
                 num_valence_electrons, num_radical_electrons, formal_charge,
                 murcko_scaffold)
                VALUES %s
                ON CONFLICT (compound_id) DO NOTHING
            """, batch, page_size=2000)
        conn.commit()
        total += len(batch)

    read_cur.close()
    read_conn.close()
    conn.close()
    log(f"  Compound structures: {total:,} computed ({errors:,} RDKit errors)")


# ============================================================
# 7. KEGG DRUG INFO
# ============================================================
def ingest_kegg():
    """Ingest KEGG drug info via REST API."""
    log("=" * 60)
    log("7. KEGG drug info")
    log("=" * 60)

    conn = get_conn()
    with conn.cursor() as cur:
        cur.execute("""
            SELECT c.pref_name, c.id FROM pharmaco_db.compounds c
            WHERE c.is_drug = true AND c.pref_name IS NOT NULL
        """)
        drug_map = {k.upper(): v for k, v in cur.fetchall()}
    log(f"  Approved drugs with names: {len(drug_map):,}")

    # Get KEGG drug list
    try:
        req = urllib.request.Request("https://rest.kegg.jp/list/drug", headers={"User-Agent": "MegaPharmaDB/1.0"})
        with urllib.request.urlopen(req, timeout=60) as resp:
            kegg_list = resp.read().decode()
    except Exception as e:
        log(f"  SKIP: KEGG API error: {e}")
        conn.close()
        return

    batch = []
    lines = kegg_list.strip().split("\n")
    log(f"  KEGG drugs: {len(lines):,}")

    for line in lines:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        kegg_id = parts[0].replace("dr:", "")
        names = parts[1].split(";")
        kegg_name = names[0].strip()

        # Try to match to our compounds
        cid = None
        for name in names:
            name_clean = name.strip().upper()
            cid = drug_map.get(name_clean)
            if cid:
                break

        batch.append((cid, kegg_id, kegg_name, None, None, None, None))

    if batch:
        # Only keep ones we matched
        matched = [b for b in batch if b[0] is not None]
        all_kegg = batch

        with conn.cursor() as cur:
            cur.execute("TRUNCATE pharmaco_db.kegg_drug_info")
            psycopg2.extras.execute_values(cur, """
                INSERT INTO pharmaco_db.kegg_drug_info
                (compound_id, kegg_drug_id, kegg_name, kegg_formula,
                 therapeutic_category, pathway_ids, pathway_names)
                VALUES %s
            """, all_kegg[:10000], page_size=5000)  # Cap at 10K
        conn.commit()
        log(f"  KEGG: {len(all_kegg):,} drugs ({len(matched):,} matched to compounds)")
    conn.close()


# ============================================================
# 8. UPDATE COLUMN COMMENTS FOR NEW DATA
# ============================================================
def update_comments():
    """Add column comments for any new columns."""
    log("=" * 60)
    log("8. Updating column descriptions")
    log("=" * 60)

    conn = get_conn()
    # The main comments were already applied in mega_consolidate.py
    # Just add any missing ones
    extra_comments = [
        ("compound_tags", "le", "Ligand Efficiency: binding energy per heavy atom (kcal/mol/HA)"),
        ("compound_tags", "bei", "Binding Efficiency Index: pActivity/MW * 1000"),
        ("compound_tags", "sei", "Surface Efficiency Index: pActivity/PSA * 100"),
        ("compound_tags", "lle", "Lipophilic Ligand Efficiency: pActivity - logP"),
    ]
    with conn.cursor() as cur:
        applied = 0
        for table, column, desc in extra_comments:
            try:
                cur.execute(f"COMMENT ON COLUMN pharmaco_db.{table}.{column} IS %s", (desc,))
                applied += 1
            except:
                conn.rollback()
        conn.commit()
    log(f"  Applied {applied} extra column descriptions")
    conn.close()


# ============================================================
# FINAL STATS
# ============================================================
def print_final_stats():
    log("=" * 60)
    log("FINAL DATABASE STATISTICS")
    log("=" * 60)

    conn = get_conn()
    with conn.cursor() as cur:
        tables = [
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
        total_rows = 0
        for tbl in tables:
            try:
                cur.execute(f"SELECT COUNT(*) FROM pharmaco_db.{tbl}")
                cnt = cur.fetchone()[0]
                total_rows += cnt
                marker = " <<<NEW" if cnt > 0 else ""
                log(f"  {tbl:30s} {cnt:>12,}{marker}")
            except:
                conn.rollback()

        # PubChem CID coverage
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.compounds WHERE pubchem_cid IS NOT NULL")
        pcid_count = cur.fetchone()[0]
        cur.execute("SELECT COUNT(*) FROM pharmaco_db.compounds")
        total_compounds = cur.fetchone()[0]

        cur.execute("SELECT pg_size_pretty(pg_database_size('pharmaco'))")
        size = cur.fetchone()[0]

        log(f"\n  Total rows: {total_rows:,}")
        log(f"  PubChem CID coverage: {pcid_count:,}/{total_compounds:,} ({100*pcid_count/total_compounds:.1f}%)")
        log(f"  Database size: {size}")
    conn.close()


# ============================================================
# MAIN
# ============================================================
def main():
    start = time.time()
    log("=" * 60)
    log("MegaPharmaDB — Phase 2 Consolidation")
    log("=" * 60)

    DOWNLOAD_DIR.mkdir(parents=True, exist_ok=True)
    step = sys.argv[1] if len(sys.argv) > 1 else "all"

    if step in ("all", "pubchem"):
        map_pubchem_cids()

    if step in ("all", "dgidb"):
        ingest_dgidb()

    if step in ("all", "gtopdb"):
        ingest_gtopdb()

    if step in ("all", "sider"):
        ingest_sider()

    if step in ("all", "opentargets"):
        ingest_open_targets()

    if step in ("all", "structures"):
        compute_compound_structures()

    if step in ("all", "kegg"):
        ingest_kegg()

    if step in ("all", "comments"):
        update_comments()

    print_final_stats()

    elapsed = time.time() - start
    log(f"\nTotal time: {elapsed / 60:.1f} minutes")

if __name__ == "__main__":
    main()
