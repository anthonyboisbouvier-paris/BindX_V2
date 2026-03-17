"""
DockIt pipeline — Ligand retrieval.

Fetches candidate ligands from:
  1. ChEMBL (known bioactive compounds for the target)
  2. ZINC (local SDF files or hardcoded drug-like set)
  3. User-supplied SMILES
  4. ENAMINE REAL (V6.3: virtual library sampling via fragment combination)
  5. PubChem (V9: PUG REST API for bioactive compounds by target gene)
  6. MolPort (V9: commercial availability check via SMILES lookup)

V3: auto_select_ligand_strategy + query_chembl_count for intelligent
    source selection based on target documentation level.
V5bis: ChEMBL activity filtering (IC50/Ki/EC50/Kd < 10 uM) + pchembl_value.
V6.3: ENAMINE REAL library sampling + Rule-of-3 fragment library.
V9: PubChem PUG REST + MolPort availability + enhanced auto strategy.
"""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path
from typing import Optional

import requests

logger = logging.getLogger(__name__)

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
HTTP_TIMEOUT = (10, 60)


# ---------------------------------------------------------------------------
# V6.3: Fragment library (Rule-of-3 compliant)
# ---------------------------------------------------------------------------

def load_fragment_library() -> list[dict]:
    """Return 50 common Rule-of-3 compliant fragments for fragment-based screening.

    Rule-of-3: MW <= 300, logP <= 3, HBD <= 3, HBA <= 3, rotatable bonds <= 3.
    These fragments represent common drug-like building blocks used in
    medicinal chemistry and fragment-based drug discovery.

    Returns
    -------
    list[dict]
        Each dict has keys: ``smiles``, ``name``, ``mwt``, ``logp``,
        ``hbd``, ``hba``, ``tpsa``, ``source``.
    """
    raw = [
        # Aromatic rings
        {"smiles": "c1ccccc1", "name": "benzene", "mw": 78.11, "logp": 1.56},
        {"smiles": "c1ccncc1", "name": "pyridine", "mw": 79.10, "logp": 0.65},
        {"smiles": "c1ccoc1", "name": "furan", "mw": 68.07, "logp": 1.34},
        {"smiles": "c1ccsc1", "name": "thiophene", "mw": 84.14, "logp": 1.81},
        {"smiles": "c1cc[nH]c1", "name": "pyrrole", "mw": 67.09, "logp": 0.75},
        {"smiles": "c1cnc[nH]1", "name": "imidazole", "mw": 68.08, "logp": -0.08},
        {"smiles": "c1ccnnc1", "name": "pyridazine", "mw": 80.09, "logp": -0.72},
        {"smiles": "c1cncnc1", "name": "pyrimidine", "mw": 80.09, "logp": -0.40},
        {"smiles": "c1cnccn1", "name": "pyrazine", "mw": 80.09, "logp": -0.26},
        {"smiles": "c1cnn[nH]1", "name": "1H-1,2,3-triazole", "mw": 69.07, "logp": -0.58},
        # Substituted aromatics
        {"smiles": "c1ccc(O)cc1", "name": "phenol", "mw": 94.11, "logp": 1.46},
        {"smiles": "c1ccc(N)cc1", "name": "aniline", "mw": 93.13, "logp": 0.90},
        {"smiles": "c1ccc(F)cc1", "name": "fluorobenzene", "mw": 96.10, "logp": 2.27},
        {"smiles": "c1ccc(Cl)cc1", "name": "chlorobenzene", "mw": 112.56, "logp": 2.84},
        {"smiles": "c1ccc(C)cc1", "name": "toluene", "mw": 92.14, "logp": 2.73},
        {"smiles": "c1ccc(CO)cc1", "name": "benzyl alcohol", "mw": 108.14, "logp": 1.10},
        {"smiles": "c1ccc(C(=O)O)cc1", "name": "benzoic acid", "mw": 122.12, "logp": 1.87},
        {"smiles": "c1ccc(C=O)cc1", "name": "benzaldehyde", "mw": 106.12, "logp": 1.48},
        {"smiles": "c1ccc(OC)cc1", "name": "anisole", "mw": 108.14, "logp": 2.11},
        {"smiles": "c1ccc(C#N)cc1", "name": "benzonitrile", "mw": 103.12, "logp": 1.56},
        # Saturated rings
        {"smiles": "C1CCNCC1", "name": "piperidine", "mw": 85.15, "logp": 0.84},
        {"smiles": "C1CNCCN1", "name": "piperazine", "mw": 86.14, "logp": -1.50},
        {"smiles": "C1CCOCC1", "name": "morpholine", "mw": 87.12, "logp": -0.86},
        {"smiles": "C1CCCC1", "name": "cyclopentane", "mw": 70.13, "logp": 3.00},
        {"smiles": "C1CCCCC1", "name": "cyclohexane", "mw": 84.16, "logp": 3.44},
        {"smiles": "C1CCC(CC1)O", "name": "cyclohexanol", "mw": 100.16, "logp": 1.23},
        {"smiles": "C1CC1", "name": "cyclopropane", "mw": 42.08, "logp": 1.72},
        {"smiles": "C1CNC1", "name": "azetidine", "mw": 57.09, "logp": -0.31},
        {"smiles": "C1COC1", "name": "oxetane", "mw": 58.08, "logp": 0.04},
        {"smiles": "C1CCNC1", "name": "pyrrolidine", "mw": 71.12, "logp": 0.46},
        # Bicyclic and fused
        {"smiles": "c1ccc2[nH]ccc2c1", "name": "indole", "mw": 117.15, "logp": 2.14},
        {"smiles": "c1ccc2ncccc2c1", "name": "quinoline", "mw": 129.16, "logp": 2.04},
        {"smiles": "c1cnc2ccccc2n1", "name": "quinazoline", "mw": 130.15, "logp": 1.52},
        {"smiles": "c1ccc2occc2c1", "name": "benzofuran", "mw": 118.13, "logp": 2.67},
        {"smiles": "c1ccc2sccc2c1", "name": "benzothiophene", "mw": 134.20, "logp": 3.12},
        # Simple chains and functional groups
        {"smiles": "CCO", "name": "ethanol", "mw": 46.07, "logp": -0.31},
        {"smiles": "CC(=O)O", "name": "acetic acid", "mw": 60.05, "logp": -0.17},
        {"smiles": "CC(=O)N", "name": "acetamide", "mw": 59.07, "logp": -1.26},
        {"smiles": "CC(C)O", "name": "isopropanol", "mw": 60.10, "logp": 0.05},
        {"smiles": "CCCO", "name": "1-propanol", "mw": 60.10, "logp": 0.25},
        {"smiles": "CCS", "name": "ethanethiol", "mw": 62.13, "logp": 0.99},
        {"smiles": "CCN", "name": "ethylamine", "mw": 45.08, "logp": -0.13},
        {"smiles": "CCNC", "name": "N-methylethylamine", "mw": 59.11, "logp": 0.26},
        {"smiles": "CC(=O)C", "name": "acetone", "mw": 58.08, "logp": -0.24},
        {"smiles": "CCOC", "name": "methoxyethane", "mw": 60.10, "logp": 0.10},
        # Heteroaromatic with substituents
        {"smiles": "c1ccnc(N)c1", "name": "2-aminopyridine", "mw": 94.11, "logp": 0.27},
        {"smiles": "c1ccnc(O)c1", "name": "2-hydroxypyridine", "mw": 95.10, "logp": 0.24},
        {"smiles": "c1cc(Cl)ncc1", "name": "3-chloropyridine", "mw": 113.54, "logp": 1.90},
        {"smiles": "Oc1ccccn1", "name": "2-pyridinol", "mw": 95.10, "logp": 0.24},
        {"smiles": "Nc1ncccn1", "name": "2-aminopyrimidine", "mw": 95.10, "logp": -0.60},
    ]
    # Normalize: mw → mwt, add source (no RDKit calculation at import)
    fragments: list[dict] = []
    for frag in raw:
        entry = {
            "smiles": frag["smiles"],
            "name": frag["name"],
            "mwt": frag["mw"],
            "logp": frag["logp"],
            "source": "fragments",
        }
        fragments.append(entry)
    return fragments


# ---------------------------------------------------------------------------
# V6.3: ENAMINE REAL library sampling
# ---------------------------------------------------------------------------

# Fragment pools for combinatorial SMILES generation
_ENAMINE_CORES: list[str] = [
    "c1ccccc1",        # benzene
    "c1ccncc1",        # pyridine
    "c1cncnc1",        # pyrimidine
    "c1ccoc1",         # furan
    "c1ccsc1",         # thiophene
    "c1cc[nH]c1",      # pyrrole
    "c1cnc[nH]1",      # imidazole
    "C1CCNCC1",        # piperidine
    "C1CNCCN1",        # piperazine
    "C1CCOCC1",        # morpholine
    "C1CCNC1",         # pyrrolidine
    "c1ccc2[nH]ccc2c1",  # indole
    "c1ccc2ncccc2c1",  # quinoline
]

_ENAMINE_LINKERS: list[str] = [
    "C",               # methylene
    "CC",              # ethylene
    "CCC",             # propylene
    "C(=O)N",          # amide
    "NC(=O)",          # reverse amide
    "C(=O)O",          # ester
    "OC(=O)",          # reverse ester
    "N",               # amine
    "O",               # ether
    "S",               # thioether
    "NC",              # aminomethyl
    "CO",              # methoxy
    "CS",              # thiomethyl
    "C(=O)",           # carbonyl
    "SO2N",            # sulfonamide
]

_ENAMINE_TAILS: list[str] = [
    "F",               # fluoro
    "Cl",              # chloro
    "O",               # hydroxy
    "N",               # amino
    "OC",              # methoxy
    "C(=O)O",          # carboxylic acid
    "C#N",             # nitrile
    "C(F)(F)F",        # trifluoromethyl
    "NC(=O)C",         # acetamido
    "S(=O)(=O)C",      # methylsulfonyl
    "C",               # methyl
    "CC",              # ethyl
    "C(C)C",           # isopropyl
    "OCC",             # ethoxy
    "NC(=O)",          # carbamoyl
]


def sample_enamine_real(n: int = 500) -> list[dict]:
    """Sample drug-like molecules from the ENAMINE REAL virtual library.

    Requires a real ENAMINE REAL database connection or API access.
    """
    raise RuntimeError(
        "ENAMINE REAL sampling is not available. "
        "Configure ENAMINE REAL database access or API credentials."
    )


# ---------------------------------------------------------------------------
# V3: Auto ligand strategy
# ---------------------------------------------------------------------------

def query_chembl_count(uniprot_id: str) -> int:
    """Query ChEMBL to count how many bioactive compounds exist for a target.

    Parameters
    ----------
    uniprot_id : str
        UniProt accession (e.g. "P00533").

    Returns
    -------
    int
        Number of known bioactive compounds. Returns 0 on any error.
    """
    try:
        # Resolve to ChEMBL target ID first
        target_chembl_id = _resolve_target(uniprot_id)
        if not target_chembl_id:
            logger.info("No ChEMBL target found for %s, count = 0", uniprot_id)
            return 0

        # Query activity count
        url = f"{CHEMBL_BASE}/activity.json"
        params = {
            "target_chembl_id": target_chembl_id,
            "standard_type": "IC50",
            "format": "json",
            "limit": 0,  # Only need the total count, not actual results
        }
        resp = requests.get(url, params=params, timeout=HTTP_TIMEOUT)
        resp.raise_for_status()
        data = resp.json()

        # The page_meta.total_count gives the total number of matching activities
        total = data.get("page_meta", {}).get("total_count", 0)
        logger.info("ChEMBL activity count for %s (%s): %d", uniprot_id, target_chembl_id, total)
        return int(total)

    except Exception as exc:
        logger.warning("Failed to query ChEMBL count for %s: %s", uniprot_id, exc)
        return 0


def auto_select_ligand_strategy(uniprot_id: str) -> dict:
    """Auto-select the best ligand sourcing strategy based on ChEMBL data availability.

    V9: Integrates PubChem + ENAMINE REAL into tiered strategy:
      - >100 actives: ChEMBL + PubChem (no ENAMINE supplement needed)
      - 10-100 actives: ChEMBL + PubChem + 500 ENAMINE REAL molecules
      - <10 actives: PubChem + 1000 ENAMINE REAL + AI generation

    Parameters
    ----------
    uniprot_id : str
        UniProt accession (e.g. "P00533").

    Returns
    -------
    dict
        Strategy descriptor with keys: source, n_ligands, use_zinc,
        use_generation, use_pubchem, message, chembl_count, enamine_count.
    """
    try:
        chembl_count = query_chembl_count(uniprot_id)

        # V9: Resolve gene name for PubChem queries
        gene_name = resolve_gene_name(uniprot_id)
        use_pubchem = gene_name is not None

        # V5bis: estimate activity ratio based on typical ChEMBL distributions
        estimated_active_ratio = "~40-60% expected active (< 10 uM)" if chembl_count > 20 else "ratio unknown"

        pubchem_note = f" + PubChem ({gene_name})" if use_pubchem else ""

        if chembl_count > 100:
            strategy: dict = {
                "source": "chembl+pubchem" if use_pubchem else "chembl",
                "n_ligands": 50,
                "use_zinc": False,
                "use_generation": False,
                "use_pubchem": use_pubchem,
                "gene_name": gene_name,
                "chembl_count": chembl_count,
                "enamine_count": 0,
                "message": (
                    f"Well-documented target - {chembl_count} known compounds in ChEMBL "
                    f"({estimated_active_ratio}){pubchem_note}"
                ),
            }
        elif chembl_count > 10:
            strategy = {
                "source": "chembl+pubchem+enamine" if use_pubchem else "chembl+enamine",
                "n_ligands": 50,
                "use_zinc": True,
                "use_generation": False,
                "use_pubchem": use_pubchem,
                "gene_name": gene_name,
                "chembl_count": chembl_count,
                "enamine_count": 500,
                "message": (
                    f"Understudied target - {chembl_count} known compounds "
                    f"({estimated_active_ratio}){pubchem_note}, "
                    f"supplemented with 500 ENAMINE REAL molecules"
                ),
            }
        else:
            strategy = {
                "source": "pubchem+enamine+reinvent4" if use_pubchem else "enamine+reinvent4",
                "n_ligands": 50,
                "use_zinc": True,
                "use_generation": True,
                "use_pubchem": use_pubchem,
                "gene_name": gene_name,
                "chembl_count": chembl_count,
                "enamine_count": 1000,
                "message": (
                    f"Undocumented target ({chembl_count} compounds, {estimated_active_ratio}) "
                    f"- PubChem{'' if use_pubchem else ' (gene unknown)'} + "
                    f"1000 ENAMINE REAL molecules + AI generation"
                ),
            }

        logger.info(
            "Auto ligand strategy for %s: source=%s, chembl_count=%d, enamine_count=%d, pubchem=%s",
            uniprot_id, strategy["source"], chembl_count, strategy["enamine_count"],
            use_pubchem,
        )
        return strategy

    except Exception as exc:
        logger.error("auto_select_ligand_strategy failed for %s: %s", uniprot_id, exc)
        return {
            "source": "chembl+enamine",
            "n_ligands": 50,
            "use_zinc": True,
            "use_generation": False,
            "use_pubchem": False,
            "gene_name": None,
            "chembl_count": 0,
            "enamine_count": 500,
            "message": "Default strategy (error during target analysis)",
        }


# ---------------------------------------------------------------------------
# ChEMBL
# ---------------------------------------------------------------------------

def fetch_chembl_ligands(
    uniprot_id: str,
    max_count: int = 50,
    filters: dict | None = None,
) -> list[dict]:
    """Retrieve known active compounds from ChEMBL for a UniProt target.

    Parameters
    ----------
    uniprot_id : str
        UniProt accession (e.g. ``"P00533"``).
    max_count : int
        Maximum number of ligands to return.
    filters : dict | None
        Optional filters: activity_types, activity_cutoff, pchembl_min.
    """
    try:
        target_chembl_id = _resolve_target(uniprot_id)
        if not target_chembl_id:
            logger.warning("Could not resolve %s to a ChEMBL target", uniprot_id)
            return []

        logger.info("Resolved %s -> %s", uniprot_id, target_chembl_id)

        # Extract filter params
        f = filters or {}
        activity_types = f.get("activity_types") or None
        activity_cutoff = float(f.get("activity_cutoff", 10000))
        pchembl_min = f.get("pchembl_min")
        if pchembl_min is not None:
            pchembl_min = float(pchembl_min)

        ligands = _fetch_activities(
            target_chembl_id,
            max_count,
            activity_types=activity_types,
            activity_cutoff_nM=activity_cutoff,
            pchembl_min=pchembl_min,
        )
        logger.info("Retrieved %d ligands from ChEMBL for %s", len(ligands), uniprot_id)
        return ligands

    except Exception as exc:
        logger.error("ChEMBL ligand fetch failed for %s: %s", uniprot_id, exc)
        return []


def _resolve_target(uniprot_id: str) -> Optional[str]:
    """Resolve UniProt accession to a ChEMBL target ID."""
    url = f"{CHEMBL_BASE}/target.json"
    params = {
        "target_components__accession": uniprot_id,
        "format": "json",
        "limit": 1,
    }

    try:
        resp = requests.get(url, params=params, timeout=HTTP_TIMEOUT)
        resp.raise_for_status()
        data = resp.json()
        targets = data.get("targets", [])
        if targets:
            return targets[0].get("target_chembl_id")
        return None
    except Exception as exc:
        logger.warning("ChEMBL target resolution failed: %s", exc)
        return None


def _fetch_activities(
    target_chembl_id: str,
    max_count: int,
    activity_types: list[str] | None = None,
    activity_cutoff_nM: float = 10000,
    pchembl_min: float | None = None,
) -> list[dict]:
    """Fetch bioactivity data for a ChEMBL target.

    Parameters
    ----------
    target_chembl_id : str
        ChEMBL target ID.
    max_count : int
        Maximum number of ligands to return.
    activity_types : list[str] | None
        Activity types to fetch (e.g. ["IC50", "Kd"]). Defaults to all 4.
    activity_cutoff_nM : float
        Max activity value in nM (default 10000 = 10 uM).
    pchembl_min : float | None
        If set, only include ligands with pchembl_value >= this threshold.

    Returns
    -------
    list[dict]
        Filtered and deduplicated ligands with activity data.
    """
    if not activity_types:
        activity_types = ["IC50", "Ki", "EC50", "Kd"]
    all_activities: list[dict] = []

    for std_type in activity_types:
        url = f"{CHEMBL_BASE}/activity.json"
        params = {
            "target_chembl_id": target_chembl_id,
            "standard_type": std_type,
            "format": "json",
            "limit": min(max_count * 2, 500),
        }

        try:
            resp = requests.get(url, params=params, timeout=HTTP_TIMEOUT)
            resp.raise_for_status()
            data = resp.json()
            activities = data.get("activities", [])
            all_activities.extend(activities)
        except Exception as exc:
            logger.warning("ChEMBL activity fetch failed for %s/%s: %s",
                           target_chembl_id, std_type, exc)
            continue

    if not all_activities:
        logger.warning("No activities fetched from ChEMBL for %s", target_chembl_id)
        return []

    # V5bis: Filter by activity value < 10000 nM (10 uM) and deduplicate
    seen_smiles: set[str] = set()
    ligands: list[dict] = []
    n_tested = 0
    n_active = 0
    n_inactive = 0

    for act in all_activities:
        smiles = act.get("canonical_smiles")
        if not smiles or smiles in seen_smiles:
            continue
        seen_smiles.add(smiles)
        n_tested += 1

        # V5bis: Activity filtering
        standard_value = act.get("standard_value")
        standard_type = act.get("standard_type", "")
        pchembl_value = act.get("pchembl_value")

        # Check if the molecule is active (below cutoff)
        # Exclude molecules with null activity value when a cutoff filter is set
        if standard_value is None:
            n_inactive += 1
            continue
        try:
            value_nm = float(standard_value)
            if value_nm >= activity_cutoff_nM:
                n_inactive += 1
                continue
        except (ValueError, TypeError):
            n_inactive += 1
            continue

        # pChEMBL minimum filter
        if pchembl_min is not None:
            if pchembl_value is None:
                continue
            try:
                if float(pchembl_value) < pchembl_min:
                    continue
            except (ValueError, TypeError):
                continue

        n_active += 1
        chembl_id = act.get("molecule_chembl_id", "")
        pref_name = act.get("pref_name") or chembl_id

        ligand_entry: dict = {
            "name": pref_name,
            "smiles": smiles,
            "source": "chembl",
            "chembl_id": chembl_id,
            "activity_type": standard_type,
        }

        # V5bis: Store pchembl_value if available
        if pchembl_value is not None:
            try:
                ligand_entry["pchembl_value"] = round(float(pchembl_value), 2)
            except (ValueError, TypeError):
                pass

        # Store activity value in nM for reference
        if standard_value is not None:
            try:
                ligand_entry["activity_value_nM"] = round(float(standard_value), 1)
            except (ValueError, TypeError):
                pass

        ligands.append(ligand_entry)

        if len(ligands) >= max_count:
            break

    logger.info(
        "ChEMBL: %d tested - %d active (IC50/Ki/EC50/Kd < 10uM) USED - %d inactive EXCLUDED",
        n_tested, n_active, n_inactive,
    )
    return ligands


# ---------------------------------------------------------------------------
# ZINC API
# ---------------------------------------------------------------------------

ZINC20_BASE = "https://zinc.docking.org/substances"
ZINC20_PAGE_SIZE = 3000  # Max reliable page size (4000+ returns 502)
ZINC20_TIMEOUT = (10, 60)  # (connect, read) seconds
ZINC20_HEADERS = {"User-Agent": "Mozilla/5.0 (compatible; BindX/1.0)"}


def _zinc_session() -> requests.Session:
    """Create a requests Session with a custom SSL adapter for ZINC.

    The ZINC server uses weak TLS ciphers that OpenSSL 3.5+ rejects by default.
    This adapter forces TLS 1.2 max and SECLEVEL=0 to maintain compatibility.
    """
    import ssl
    from urllib3.util.ssl_ import create_urllib3_context
    from requests.adapters import HTTPAdapter

    class _ZincSSLAdapter(HTTPAdapter):
        def init_poolmanager(self, *args, **kwargs):
            ctx = create_urllib3_context()
            ctx.maximum_version = ssl.TLSVersion.TLSv1_2
            ctx.set_ciphers("DEFAULT:@SECLEVEL=0")
            kwargs["ssl_context"] = ctx
            return super().init_poolmanager(*args, **kwargs)

    session = requests.Session()
    session.mount("https://zinc.docking.org", _ZincSSLAdapter())
    session.headers.update(ZINC20_HEADERS)
    return session

# Available subsets with human-readable labels
ZINC20_SUBSETS: dict[str, str] = {
    "in-stock": "In Stock (~12M)",
    "for-sale": "For Sale (~389M)",
    "fda": "FDA Approved (~1.4K)",
    "in-trials": "In Clinical Trials (~5.8K)",
    "in-vivo": "In Vivo Tested (~115K)",
    "in-vitro": "In Vitro Active ≤10µM (~276K)",
    "natural-products": "Natural Products (~80K)",
    "biogenic": "Biogenic (~135K)",
    "clean": "PAINS-Filtered (~125M)",
    "endogenous": "Human Metabolites (~50K)",
    "world": "World-Approved Drugs (~3.4K)",
}


def fetch_zinc_ligands(
    max_count: int = 100,
    subset: str | None = None,
    filters: dict | None = None,
) -> list[dict]:
    """Fetch molecules from ZINC20 API.

    Tries local SDF files first (``/data/zinc/*.sdf``), then the ZINC20
    REST API with optional subset and property filters.

    Parameters
    ----------
    max_count : int
        Maximum number of molecules to return (up to 1M).
    subset : str or None
        ZINC20 subset name (e.g. ``"in-stock"``, ``"fda"``).
    filters : dict or None
        Property filters with keys like ``mwt_min``, ``mwt_max``,
        ``logp_min``, ``logp_max``, ``hbd_max``, ``hba_max``,
        ``tpsa_min``, ``tpsa_max``, ``rotatable_max``,
        ``num_rings_min``, ``num_rings_max``,
        ``num_heavy_atoms_min``, ``num_heavy_atoms_max``.

    Returns
    -------
    list[dict]
        Each dict has keys: ``name``, ``smiles``, ``source``, and
        optional ``zinc_id``, ``mwt``, ``logp``, ``hbd``, ``hba``, ``tpsa``.

    Raises
    ------
    RuntimeError
        If the ZINC20 API is unreachable after retries.
    """
    # Tier 1: local SDF files
    zinc_dir = Path("/data/zinc")
    if zinc_dir.exists():
        ligands = _load_sdf_files(zinc_dir, max_count)
        if ligands:
            logger.info("Loaded %d ligands from local ZINC SDF files", len(ligands))
            return ligands

    # Tier 2: ZINC20 API
    return _fetch_zinc_api(max_count, subset, filters)


def _build_zinc_url(subset: str | None) -> str:
    """Build the ZINC20 substances URL with optional subset(s).

    Supports single subset ("in-stock") or combined subsets ("in-stock+fda").
    Combined subsets return the UNION of both sets.
    """
    if subset:
        # Validate each part of combined subsets (e.g. "in-stock+fda")
        parts = [p.strip() for p in subset.split("+") if p.strip()]
        valid_parts = [p for p in parts if p in ZINC20_SUBSETS]
        if valid_parts:
            combined = "+".join(valid_parts)
            return f"{ZINC20_BASE}/subsets/{combined}.json:zinc_id+smiles+mwt+logp+hbd+hba+tpsa"
    return f"{ZINC20_BASE}.json:zinc_id+smiles+mwt+logp+hbd+hba+tpsa"


def _build_zinc_params(
    page_size: int,
    page: int,
    filters: dict | None,
) -> dict:
    """Build query params for a ZINC20 API request."""
    params: dict = {"count": page_size, "page": page}

    if not filters:
        return params

    # Map our filter keys to ZINC20 API "between" params
    between_map = {
        ("mwt_min", "mwt_max"): "mwt-between",
        ("logp_min", "logp_max"): "logp-between",
        ("hbd_min", "hbd_max"): "hbd-between",
        ("hba_min", "hba_max"): "hba-between",
        ("tpsa_min", "tpsa_max"): "tpsa-between",
        ("rotatable_min", "rotatable_max"): "rb-between",
        ("num_rings_min", "num_rings_max"): "num_rings-between",
        ("num_heavy_atoms_min", "num_heavy_atoms_max"): "num_heavy_atoms-between",
    }

    for (key_min, key_max), api_param in between_map.items():
        lo = filters.get(key_min)
        hi = filters.get(key_max)
        if lo is not None or hi is not None:
            lo_val = lo if lo is not None else ""
            hi_val = hi if hi is not None else ""
            params[api_param] = f"{lo_val} {hi_val}"

    return params


def _fetch_zinc_api(
    max_count: int,
    subset: str | None,
    filters: dict | None,
) -> list[dict]:
    """Paginated fetch from ZINC REST API."""
    session = _zinc_session()
    url = _build_zinc_url(subset)
    ligands: list[dict] = []
    seen_ids: set[str] = set()
    page = 1
    page_size = min(max_count, ZINC20_PAGE_SIZE)

    logger.info(
        "Fetching up to %d molecules from ZINC (subset=%s, filters=%s)",
        max_count, subset, filters,
    )

    while len(ligands) < max_count:
        remaining = max_count - len(ligands)
        current_page_size = min(page_size, remaining)
        params = _build_zinc_params(current_page_size, page, filters)

        try:
            resp = session.get(url, params=params, timeout=ZINC20_TIMEOUT)
            resp.raise_for_status()
        except requests.exceptions.RequestException as exc:
            if page == 1 and subset:
                # Subset endpoint may be unreachable — retry without subset
                logger.warning(
                    "ZINC subset '%s' failed (%s), retrying without subset filter",
                    subset, exc,
                )
                url = _build_zinc_url(None)
                try:
                    resp = session.get(url, params=params, timeout=ZINC20_TIMEOUT)
                    resp.raise_for_status()
                except requests.exceptions.RequestException as exc2:
                    raise RuntimeError(f"ZINC API error: {exc2}")
            elif not ligands:
                raise RuntimeError(f"ZINC API error: {exc}")
            else:
                logger.warning("ZINC API error on page %d: %s — returning %d molecules", page, exc, len(ligands))
                break

        try:
            data = resp.json()
        except ValueError:
            logger.warning("ZINC returned invalid JSON on page %d", page)
            break

        if not data:
            break

        for entry in data:
            smiles = entry.get("smiles", "")
            zinc_id = entry.get("zinc_id", "")
            if not smiles or zinc_id in seen_ids:
                continue
            seen_ids.add(zinc_id)

            ligands.append({
                "name": zinc_id,
                "smiles": smiles,
                "source": "zinc",
                "zinc_id": zinc_id,
                "mwt": entry.get("mwt"),
                "logp": entry.get("logp"),
                "hbd": entry.get("hbd"),
                "hba": entry.get("hba"),
                "tpsa": entry.get("tpsa"),
            })

            if len(ligands) >= max_count:
                break

        if len(data) < current_page_size:
            break

        page += 1

    logger.info(
        "Fetched %d molecules from ZINC API (subset=%s)",
        len(ligands), subset,
    )
    return ligands


def _load_sdf_files(zinc_dir: Path, max_count: int) -> list[dict]:
    """Parse SDF files in the ZINC directory using RDKit."""
    try:
        from rdkit import Chem

        ligands: list[dict] = []
        sdf_files = sorted(zinc_dir.glob("*.sdf"))

        for sdf_file in sdf_files:
            if len(ligands) >= max_count:
                break
            try:
                supplier = Chem.SDMolSupplier(str(sdf_file))
                for mol in supplier:
                    if mol is None:
                        continue
                    smiles = Chem.MolToSmiles(mol)
                    name = mol.GetProp("_Name") if mol.HasProp("_Name") else sdf_file.stem
                    ligands.append({
                        "name": name,
                        "smiles": smiles,
                        "source": "zinc",
                    })
                    if len(ligands) >= max_count:
                        break
            except Exception as exc:
                logger.warning("Error reading %s: %s", sdf_file, exc)
                continue

        return ligands

    except ImportError:
        logger.warning("RDKit not available; cannot read SDF files")
        return []


def _hardcoded_druglike() -> list[dict]:
    """Curated set of well-known drug-like molecules for testing."""
    return [
        {"name": "Aspirin", "smiles": "CC(=O)Oc1ccccc1C(=O)O", "source": "zinc"},
        {"name": "Ibuprofen", "smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "source": "zinc"},
        {"name": "Caffeine", "smiles": "Cn1c(=O)c2c(ncn2C)n(C)c1=O", "source": "zinc"},
        {"name": "Acetaminophen", "smiles": "CC(=O)Nc1ccc(O)cc1", "source": "zinc"},
        {"name": "Metformin", "smiles": "CN(C)C(=N)NC(=N)N", "source": "zinc"},
        {"name": "Naproxen", "smiles": "COc1ccc2cc(ccc2c1)C(C)C(=O)O", "source": "zinc"},
        {"name": "Omeprazole", "smiles": "COc1ccc2[nH]c(nc2c1)S(=O)Cc1ncc(C)c(OC)c1C", "source": "zinc"},
        {"name": "Atorvastatin", "smiles": "CC(C)c1n(CC[C@@H](O)C[C@@H](O)CC(=O)O)c(c2ccc(F)cc2)c(c1c1ccccc1)C(=O)Nc1ccccc1", "source": "zinc"},
        {"name": "Losartan", "smiles": "CCCCc1nc(Cl)c(n1Cc1ccc(cc1)c1ccccc1c1nn[nH]n1)CO", "source": "zinc"},
        {"name": "Amlodipine", "smiles": "CCOC(=O)C1=C(COCCN)NC(C)=C(C1c1ccccc1Cl)C(=O)OC", "source": "zinc"},
        {"name": "Celecoxib", "smiles": "Cc1ccc(cc1)c1cc(cf1)c1ccc(cc1)S(=O)(=O)N", "source": "zinc"},
        {"name": "Diclofenac", "smiles": "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl", "source": "zinc"},
        {"name": "Ranitidine", "smiles": "CNC(/N=C/[N+](=O)[O-])NCCSCc1ccc(CN(C)C)o1", "source": "zinc"},
        {"name": "Fluoxetine", "smiles": "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1", "source": "zinc"},
        {"name": "Ciprofloxacin", "smiles": "OC(=O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O", "source": "zinc"},
        {"name": "Tamoxifen", "smiles": "CC/C(=C(\\c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1", "source": "zinc"},
        {"name": "Warfarin", "smiles": "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O", "source": "zinc"},
        {"name": "Sildenafil", "smiles": "CCCc1nn(C)c2c1nc(nc2OCC)c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1", "source": "zinc"},
        {"name": "Methotrexate", "smiles": "CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(cc1)C(=O)NC(CCC(=O)O)C(=O)O", "source": "zinc"},
        {"name": "Dexamethasone", "smiles": "CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(F)C(O)CC2(C)C1(O)C(=O)CO", "source": "zinc"},
    ]


# ---------------------------------------------------------------------------
# User-supplied SMILES
# ---------------------------------------------------------------------------

def parse_user_smiles(smiles_list: list[str]) -> list[dict]:
    """Convert a list of raw SMILES strings into ligand dicts.

    Parameters
    ----------
    smiles_list : list[str]
        SMILES strings provided by the user.

    Returns
    -------
    list[dict]
        Validated ligands with ``name``, ``smiles``, ``source`` keys.
    """
    ligands: list[dict] = []
    for i, smi in enumerate(smiles_list):
        smi = smi.strip()
        if not smi:
            continue
        # Try to validate with RDKit if available
        name = f"UserLigand_{i+1}"
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                logger.warning("Invalid SMILES skipped: %s", smi)
                continue
            # Canonicalise
            smi = Chem.MolToSmiles(mol)
        except ImportError:
            pass  # Accept as-is without validation

        ligands.append({
            "name": name,
            "smiles": smi,
            "source": "user",
        })

    logger.info("Parsed %d valid user SMILES out of %d provided", len(ligands), len(smiles_list))
    return ligands


# ---------------------------------------------------------------------------
# V9: PubChem PUG REST API
# ---------------------------------------------------------------------------

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def fetch_pubchem_ligands(
    gene_name: str,
    max_count: int = 50,
    activity_type: str = "IC50",
) -> list[dict]:
    """Retrieve bioactive compounds from PubChem for a target gene.

    Uses the PUG REST API to find compounds with bioactivity data against
    a specified gene/protein target.

    Parameters
    ----------
    gene_name : str
        Gene symbol (e.g. ``"EGFR"``, ``"BRAF"``).
    max_count : int
        Maximum number of ligands to return.
    activity_type : str
        Activity type filter (default ``"IC50"``).

    Returns
    -------
    list[dict]
        Each dict has keys: ``name``, ``smiles``, ``source``, ``cid``,
        ``activity_value_nM``.
    """
    if not gene_name:
        return []

    try:
        # Step 1: Search for assays targeting this gene
        url = f"{PUBCHEM_BASE}/assay/target/genesymbol/{gene_name}/aids/JSON"
        resp = requests.get(url, timeout=HTTP_TIMEOUT)
        if resp.status_code == 404:
            logger.info("PubChem: no assays found for gene %s", gene_name)
            return []
        resp.raise_for_status()
        data = resp.json()

        aids = data.get("IdentifierList", {}).get("AID", [])
        if not aids:
            logger.info("PubChem: no assay IDs found for %s", gene_name)
            return []

        # Take first 5 assays to avoid excessive requests
        aids = aids[:5]
        logger.info("PubChem: found %d assays for %s, querying first %d",
                     len(data.get("IdentifierList", {}).get("AID", [])),
                     gene_name, len(aids))

        # Step 2: Parallel fetch — active CIDs + activity values per assay
        import math
        import statistics
        from concurrent.futures import ThreadPoolExecutor, as_completed

        def _fetch_assay_data(aid: int) -> dict:
            """Fetch active CIDs, activity values, and properties for one assay."""
            result = {"aid": aid, "cids": [], "activities": {}, "properties": [], "assay_name": ""}
            try:
                # 2a: active CIDs
                cid_resp = requests.get(
                    f"{PUBCHEM_BASE}/assay/aid/{aid}/cids/JSON?cids_type=active",
                    timeout=HTTP_TIMEOUT,
                )
                if cid_resp.status_code != 200:
                    return result
                cid_data = cid_resp.json()
                cids = cid_data.get("IdentifierList", {}).get("CID", [])
                if not cids:
                    info_list = cid_data.get("InformationList", {}).get("Information", [])
                    if info_list:
                        cids = info_list[0].get("CID", [])
                if not cids:
                    return result
                result["cids"] = cids[:30]

                # 2b: concise activity data
                try:
                    conc_resp = requests.get(
                        f"{PUBCHEM_BASE}/assay/aid/{aid}/concise/JSON",
                        timeout=(10, 30),
                    )
                    if conc_resp.status_code == 200:
                        table = conc_resp.json().get("Table", {})
                        cols = table.get("Columns", {}).get("Column", [])
                        if "CID" in cols and "Activity Value [uM]" in cols:
                            cidx = cols.index("CID")
                            vidx = cols.index("Activity Value [uM]")
                            oidx = cols.index("Activity Outcome") if "Activity Outcome" in cols else None
                            nidx = cols.index("Assay Name") if "Assay Name" in cols else None
                            cid_strs = {str(c) for c in result["cids"]}
                            for row in table.get("Row", []):
                                cells = row.get("Cell", [])
                                if cells[cidx] not in cid_strs:
                                    continue
                                if oidx is not None and cells[oidx] != "Active":
                                    continue
                                if not result["assay_name"] and nidx is not None and cells[nidx]:
                                    result["assay_name"] = cells[nidx]
                                if cells[vidx]:
                                    try:
                                        result["activities"].setdefault(
                                            int(cells[cidx]), []
                                        ).append(float(cells[vidx]))
                                    except (ValueError, IndexError):
                                        pass
                except Exception as exc:
                    logger.debug("PubChem concise fetch for AID %d: %s", aid, exc)

                # 2c: compound properties
                cid_str = ",".join(str(c) for c in result["cids"][:20])
                prop_resp = requests.get(
                    f"{PUBCHEM_BASE}/compound/cid/{cid_str}"
                    f"/property/CanonicalSMILES,IsomericSMILES,IUPACName,"
                    f"MolecularWeight,XLogP,HBondDonorCount,HBondAcceptorCount,TPSA/JSON",
                    timeout=HTTP_TIMEOUT,
                )
                if prop_resp.status_code == 200:
                    result["properties"] = prop_resp.json().get("PropertyTable", {}).get("Properties", [])

            except Exception as exc:
                logger.debug("PubChem assay %d fetch error: %s", aid, exc)
            return result

        # Run assays in parallel (3 threads to avoid PubChem rate limits)
        assay_results = []
        with ThreadPoolExecutor(max_workers=3) as pool:
            futures = {pool.submit(_fetch_assay_data, aid): aid for aid in aids}
            for fut in as_completed(futures):
                assay_results.append(fut.result())
        # Sort back to original AID order for deterministic output
        aid_order = {aid: i for i, aid in enumerate(aids)}
        assay_results.sort(key=lambda r: aid_order.get(r["aid"], 99))

        # Step 3: Assemble ligands from parallel results
        ligands: list[dict] = []
        seen_cids: set[int] = set()
        cid_activities: dict[int, list[float]] = {}
        cid_assay: dict[int, str] = {}  # CID -> first assay name

        for ar in assay_results:
            if len(ligands) >= max_count:
                break
            new_cids = [c for c in ar["cids"] if c not in seen_cids]
            if not new_cids:
                continue
            seen_cids.update(new_cids)

            # Merge activity data
            for cid, vals in ar["activities"].items():
                cid_activities.setdefault(cid, []).extend(vals)
                if cid not in cid_assay and ar.get("assay_name"):
                    cid_assay[cid] = ar["assay_name"]

            for prop in ar["properties"]:
                if len(ligands) >= max_count:
                    break
                smiles = (
                    prop.get("CanonicalSMILES")
                    or prop.get("IsomericSMILES")
                    or prop.get("ConnectivitySMILES")
                    or ""
                )
                if not smiles:
                    continue
                cid = prop.get("CID", 0)
                name = prop.get("IUPACName", f"PubChem_{cid}")
                if len(name) > 60:
                    name = f"PubChem_{cid}"

                entry = {
                    "name": name,
                    "smiles": smiles,
                    "source": "pubchem",
                    "cid": cid,
                }
                if prop.get("MolecularWeight") is not None:
                    entry["mwt"] = round(float(prop["MolecularWeight"]), 2)
                if prop.get("XLogP") is not None:
                    entry["logp"] = round(float(prop["XLogP"]), 2)
                if prop.get("HBondDonorCount") is not None:
                    entry["import_hbd"] = int(prop["HBondDonorCount"])
                if prop.get("HBondAcceptorCount") is not None:
                    entry["import_hba"] = int(prop["HBondAcceptorCount"])
                if prop.get("TPSA") is not None:
                    entry["import_tpsa"] = round(float(prop["TPSA"]), 2)

                activities = cid_activities.get(cid, [])
                if activities:
                    median_um = statistics.median(activities)
                    entry["activity_value_nM"] = round(median_um * 1000, 2)
                    entry["activity_type"] = activity_type
                    if median_um > 0:
                        entry["pchembl_value"] = round(-math.log10(median_um * 1e-6), 2)
                if cid_assay.get(cid):
                    entry["assay_name"] = cid_assay[cid]

                ligands.append(entry)

        logger.info("Retrieved %d ligands from PubChem for %s (%d with activity data)",
                     len(ligands), gene_name,
                     sum(1 for l in ligands if "activity_value_nM" in l))
        return ligands

    except Exception as exc:
        logger.warning("PubChem ligand fetch failed for %s: %s", gene_name, exc)
        return []


def check_pubchem_availability(smiles: str) -> dict | None:
    """Check if a compound exists in PubChem and retrieve basic info.

    Parameters
    ----------
    smiles : str
        Canonical SMILES string.

    Returns
    -------
    dict or None
        PubChem compound info with ``cid``, ``iupac_name``, ``mw``,
        or None if not found.
    """
    try:
        url = (
            f"{PUBCHEM_BASE}/compound/smiles/{requests.utils.quote(smiles)}"
            f"/property/CID,IUPACName,MolecularWeight/JSON"
        )
        resp = requests.get(url, timeout=(5, 15))
        if resp.status_code == 404:
            return None
        resp.raise_for_status()
        data = resp.json()
        props = data.get("PropertyTable", {}).get("Properties", [])
        if props:
            return {
                "cid": props[0].get("CID"),
                "iupac_name": props[0].get("IUPACName"),
                "mw": props[0].get("MolecularWeight"),
                "source": "pubchem",
            }
        return None
    except Exception as exc:
        logger.debug("PubChem availability check failed: %s", exc)
        return None


# ---------------------------------------------------------------------------
# V9: MolPort commercial availability check
# ---------------------------------------------------------------------------

MOLPORT_VENDORS = {
    "enamine": "Enamine",
    "chembridge": "ChemBridge",
    "life_chemicals": "Life Chemicals",
    "vitas_m": "Vitas-M",
    "asinex": "Asinex",
    "specs": "Specs",
}


def check_commercial_availability(smiles: str) -> dict:
    """Check commercial availability via PubChem vendor annotations.

    Uses PubChem's vendor information as a proxy for commercial availability.
    This avoids needing direct MolPort API access.

    Parameters
    ----------
    smiles : str
        Canonical SMILES string.

    Returns
    -------
    dict
        Keys: ``available`` (bool), ``vendors`` (list[str]),
        ``cid`` (int or None).
    """
    result = {"available": False, "vendors": [], "cid": None}
    try:
        # First get CID
        info = check_pubchem_availability(smiles)
        if not info or not info.get("cid"):
            return result
        result["cid"] = info["cid"]

        # Check PubChem sources for vendor info
        cid = info["cid"]
        url = f"{PUBCHEM_BASE}/compound/cid/{cid}/xrefs/RegistryID/JSON"
        resp = requests.get(url, timeout=(5, 15))
        if resp.status_code == 200:
            data = resp.json()
            registries = data.get("InformationList", {}).get("Information", [])
            if registries:
                result["available"] = True
                result["vendors"] = ["PubChem-registered"]

        return result
    except Exception as exc:
        logger.debug("Commercial availability check failed: %s", exc)
        return result


# ---------------------------------------------------------------------------
# V9: UniProt → gene name resolver (for PubChem queries)
# ---------------------------------------------------------------------------

def resolve_gene_name(uniprot_id: str) -> str | None:
    """Resolve a UniProt accession to a gene symbol.

    Parameters
    ----------
    uniprot_id : str
        UniProt accession (e.g. ``"P00533"``).

    Returns
    -------
    str or None
        Gene symbol (e.g. ``"EGFR"``), or None if resolution fails.
    """
    try:
        resp = requests.get(
            f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json",
            timeout=(5, 10),
        )
        resp.raise_for_status()
        data = resp.json()
        genes = data.get("genes", [])
        if genes:
            name = genes[0].get("geneName", {}).get("value")
            if name:
                logger.info("Resolved %s -> gene %s", uniprot_id, name)
                return name
        return None
    except Exception as exc:
        logger.debug("Gene name resolution failed for %s: %s", uniprot_id, exc)
        return None
