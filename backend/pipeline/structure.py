"""
DockIt pipeline -- Protein structure retrieval.

Fetches a 3D structure for a given UniProt accession or raw sequence:
  1. Try RCSB PDB (experimental, best resolution with co-crystal ligand) -- returns (path, "pdb_experimental").
  2. Try AlphaFold DB (fast, pre-computed) -- returns (path, "alphafold").
  3. Fallback to ESMFold API (slower, on-demand folding) -- returns (path, "esmfold").
  4. UniProt PDB cross-ref -- returns (path, "pdb_experimental").
  5. Last-resort mock for offline development -- returns (path, "mock").

V3: fetch_structure now returns (Path, str) tuple with source name.
    New fetch_structure_from_sequence for direct sequence input.
V5bis: query_rcsb_pdb added as primary lookup (RCSB Search API v2).
       fetch_structure now returns (Path, str, dict|None) with optional pdb_info.
V6: predict_disorder() added for IDR prediction (IUPred3 + mock fallback).
"""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path
from typing import Optional

import requests

logger = logging.getLogger(__name__)

ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
ESMFOLD_API = "https://api.esm.metaprotein.ai/fold"

# Timeout for HTTP calls (connect, read) in seconds
HTTP_TIMEOUT = (10, 120)


# ---------------------------------------------------------------------------
# RCSB PDB lookup (V5bis)
# ---------------------------------------------------------------------------

def query_rcsb_pdb(uniprot_id: str) -> Optional[dict]:
    """Query RCSB PDB for experimental structures by UniProt accession.

    Backward-compatible wrapper — returns the single best structure.
    """
    results = query_rcsb_pdb_multi(uniprot_id)
    return results[0] if results else None


def _rcsb_search(uniprot_id: str, extra_nodes: list = None, rows: int = 10) -> list[str]:
    """Run an RCSB Search API v2 query and return PDB IDs."""
    nodes = [
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": uniprot_id,
            },
        },
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entry_info.resolution_combined",
                "operator": "less",
                "value": 3.5,
            },
        },
    ]
    if extra_nodes:
        nodes.extend(extra_nodes)

    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": nodes,
        },
        "return_type": "entry",
        "request_options": {
            "results_content_type": ["experimental"],
            "sort": [
                {"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}
            ],
            "paginate": {"start": 0, "rows": rows},
        },
    }

    resp = requests.post(
        "https://search.rcsb.org/rcsbsearch/v2/query",
        json=query,
        timeout=15,
    )
    if resp.status_code != 200:
        return []
    return [e.get("identifier", "") for e in resp.json().get("result_set", []) if e.get("identifier")]


# Common solvent/ion IDs to exclude from "ligand" detection
_SOLVENT_IDS = frozenset({
    "HOH", "SO4", "PO4", "GOL", "EDO", "ACE", "NAG", "MAN", "GAL",
    "PEG", "DMS", "BME", "CL", "NA", "MG", "ZN", "CA", "K", "MN",
    "FE", "CU", "CO", "NI", "IOD", "BR", "FMT", "ACT", "TRS", "MPD",
    "PG4", "EPE", "MES", "CIT", "TAR", "SUC", "MLI", "NH4",
})


def _fetch_pdb_detail(pdb_id: str) -> Optional[dict]:
    """Fetch details + ligand info for a single PDB entry."""
    try:
        detail_resp = requests.get(
            f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}",
            timeout=10,
        )
        if detail_resp.status_code != 200:
            return None
        detail = detail_resp.json()
    except Exception:
        return None

    resolution_raw = detail.get("rcsb_entry_info", {}).get("resolution_combined", [99.0])
    resolution: float = (
        resolution_raw[0]
        if isinstance(resolution_raw, list) and resolution_raw
        else float(resolution_raw) if resolution_raw else 99.0
    )
    method: str = (
        detail.get("exptl", [{}])[0].get("method", "UNKNOWN")
        if detail.get("exptl")
        else "UNKNOWN"
    )

    # Check multiple nonpolymer entities for drug-like ligands
    ligand_id: Optional[str] = None
    ligand_name: Optional[str] = None
    n_nonpoly = detail.get("rcsb_entry_info", {}).get("nonpolymer_entity_count", 0) or 0
    for eidx in range(1, min(n_nonpoly + 1, 6)):
        try:
            lig_resp = requests.get(
                f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id}/{eidx}",
                timeout=5,
            )
            if lig_resp.status_code != 200:
                continue
            lig_data = lig_resp.json()
            comp_id = lig_data.get("pdbx_entity_nonpoly", {}).get("comp_id")
            if comp_id and comp_id.upper() not in _SOLVENT_IDS:
                ligand_id = comp_id
                ligand_name = lig_data.get("pdbx_entity_nonpoly", {}).get("name")
                break
        except Exception:
            continue

    score: float = -resolution
    if ligand_id:
        score += 15.0
    if "X-RAY" in method.upper():
        score += 5.0

    return {
        "pdb_id": pdb_id,
        "resolution": round(resolution, 2),
        "method": method,
        "ligand_id": ligand_id,
        "ligand_name": ligand_name,
        "download_url": f"https://files.rcsb.org/download/{pdb_id}.pdb",
        "score": round(score, 2),
    }


def query_rcsb_pdb_multi(uniprot_id: str, max_results: int = 5) -> list[dict]:
    """Query RCSB PDB for experimental structures by UniProt accession.

    Runs two searches: best resolution + structures with bound ligands.
    Returns up to ``max_results`` unique structures sorted by score.

    Parameters
    ----------
    uniprot_id : str
        UniProt accession (e.g. ``"P00533"``).
    max_results : int
        Maximum number of structures to return.

    Returns
    -------
    list[dict]
        Each dict has keys: ``pdb_id``, ``resolution``, ``method``, ``ligand_id``,
        ``ligand_name``, ``download_url``, ``score``.
    """
    try:
        logger.info("Querying RCSB PDB for %s (multi-structure) ...", uniprot_id)

        # Search 1: best resolution (any)
        best_res_ids = _rcsb_search(uniprot_id, rows=8)

        # Search 2: structures with bound small molecules (holo)
        holo_filter = {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entry_info.nonpolymer_entity_count",
                "operator": "greater",
                "value": 0,
            },
        }
        holo_ids = _rcsb_search(uniprot_id, extra_nodes=[holo_filter], rows=8)

        # Merge unique PDB IDs (best resolution first, then holo)
        seen = set()
        merged_ids: list[str] = []
        for pid in best_res_ids + holo_ids:
            if pid not in seen:
                seen.add(pid)
                merged_ids.append(pid)

        if not merged_ids:
            logger.info("RCSB PDB: no results for %s", uniprot_id)
            return []

        # Fetch details for all candidates
        scored: list[dict] = []
        for pdb_id in merged_ids:
            info = _fetch_pdb_detail(pdb_id)
            if info:
                scored.append(info)

        # Sort by score descending, take top N
        scored.sort(key=lambda x: x["score"], reverse=True)
        top = scored[:max_results]

        for s in top:
            logger.info(
                "RCSB PDB: %s (%.2f A, %s, ligand=%s, score=%.1f)",
                s["pdb_id"], s["resolution"], s["method"],
                s.get("ligand_id") or "none", s["score"],
            )
        return top

    except Exception as exc:
        logger.warning("RCSB PDB query failed for %s: %s", uniprot_id, exc)
        return []


def _download_pdb_from_rcsb(pdb_info: dict, pdb_path: Path) -> bool:
    """Download PDB file from RCSB given pdb_info dict.

    Parameters
    ----------
    pdb_info : dict
        Must contain ``download_url`` key.
    pdb_path : Path
        Destination file path.

    Returns
    -------
    bool
        True if download was successful and file contains ATOM records.
    """
    url = pdb_info.get("download_url", "")
    if not url:
        return False
    try:
        resp = requests.get(url, timeout=HTTP_TIMEOUT)
        resp.raise_for_status()
        content = resp.text
        if "ATOM" in content and len(content) > 100:
            pdb_path.write_text(content)
            logger.info(
                "Downloaded PDB %s from RCSB (%d bytes)",
                pdb_info.get("pdb_id", "?"),
                len(content),
            )
            return True
        logger.warning(
            "RCSB PDB file for %s has no ATOM records or is too small",
            pdb_info.get("pdb_id", "?"),
        )
        return False
    except Exception as exc:
        logger.warning(
            "Failed to download PDB from RCSB (%s): %s",
            pdb_info.get("pdb_id", "?"),
            exc,
        )
        return False


# ---------------------------------------------------------------------------
# Disorder prediction (V6)
# ---------------------------------------------------------------------------

def predict_disorder(sequence: str) -> dict:
    """Predict intrinsically disordered regions (IDRs) in a protein sequence.

    Uses IUPred3 API with a mock fallback based on amino acid composition
    heuristics (Dunker et al. disorder-promoting residues).

    An IDR is defined as a consecutive stretch of >= 20 residues with
    disorder score > 0.5.

    Parameters
    ----------
    sequence : str
        Raw amino acid sequence (single-letter code, e.g. ``"MTEYKLVVV..."``).

    Returns
    -------
    dict
        Keys:
        - ``disorder_scores``: list[float] -- per-residue disorder scores [0, 1].
        - ``idr_regions``: list[tuple[int, int]] -- (start, end) pairs of IDRs.
        - ``fraction_disordered``: float -- fraction of residues in IDRs.
        - ``method``: ``"iupred3"`` or ``"mock"``.
    """
    if not sequence:
        return {
            "disorder_scores": [],
            "idr_regions": [],
            "fraction_disordered": 0.0,
            "method": "mock",
        }

    # ------------------------------------------------------------------
    # Try IUPred3 API first
    # ------------------------------------------------------------------
    try:
        logger.info("Querying IUPred3 API for disorder prediction (%d residues) ...", len(sequence))
        resp = requests.post(
            "https://iupred3.elte.hu/iupred3/",
            data={"sequence": sequence, "iupred_type": "long"},
            timeout=(5, 30),
        )
        if resp.status_code == 200:
            # IUPred3 returns a text/plain response with tab-separated columns:
            # position  residue  iupred_score
            scores: list[float] = []
            for line in resp.text.strip().splitlines():
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        scores.append(round(float(parts[2]), 3))
                    except ValueError:
                        continue

            if len(scores) == len(sequence):
                idr_regions = _find_idr_regions(scores)
                fraction = sum(e - s for s, e in idr_regions) / max(len(sequence), 1)
                logger.info(
                    "IUPred3 prediction complete: %d IDRs, %.1f%% disordered",
                    len(idr_regions), fraction * 100,
                )
                return {
                    "disorder_scores": scores,
                    "idr_regions": idr_regions,
                    "fraction_disordered": round(fraction, 3),
                    "method": "iupred3",
                }
            else:
                logger.debug(
                    "IUPred3 returned %d scores for %d residues, falling back to mock",
                    len(scores), len(sequence),
                )
    except Exception as exc:
        logger.debug("IUPred3 API unavailable: %s", exc)

    # ------------------------------------------------------------------
    # Mock fallback: composition-based heuristic
    # ------------------------------------------------------------------
    # Disorder-promoting amino acids (from Dunker et al.):
    #   G, P, S, Q, E, K, A
    # Order-promoting amino acids:
    #   W, Y, F, I, L, V, C, N
    logger.info("Using mock disorder prediction (composition heuristic) for %d residues", len(sequence))
    disorder_promoting = set("GPSQEKA")
    scores = []
    window = 21  # sliding window size

    for i in range(len(sequence)):
        start = max(0, i - window // 2)
        end = min(len(sequence), i + window // 2 + 1)
        region = sequence[start:end]
        dp_fraction = sum(1 for aa in region if aa in disorder_promoting) / len(region)
        # Map to 0-1 score: pure disorder-promoting -> ~0.8, pure order -> ~0.1
        score = 0.1 + 0.7 * dp_fraction
        scores.append(round(score, 3))

    idr_regions = _find_idr_regions(scores)
    fraction = sum(e - s for s, e in idr_regions) / max(len(sequence), 1)

    logger.info(
        "Mock disorder prediction complete: %d IDRs, %.1f%% disordered",
        len(idr_regions), fraction * 100,
    )

    return {
        "disorder_scores": scores,
        "idr_regions": idr_regions,
        "fraction_disordered": round(fraction, 3),
        "method": "mock",
    }


def _find_idr_regions(scores: list[float], min_length: int = 20, threshold: float = 0.5) -> list[tuple[int, int]]:
    """Identify IDR regions from per-residue disorder scores.

    An IDR is a consecutive stretch of >= *min_length* residues with
    disorder score > *threshold*.

    Parameters
    ----------
    scores : list[float]
        Per-residue disorder scores.
    min_length : int
        Minimum number of consecutive residues to qualify as an IDR.
    threshold : float
        Score above which a residue is considered disordered.

    Returns
    -------
    list[tuple[int, int]]
        List of (start, end) pairs (0-indexed, end exclusive).
    """
    idr_regions: list[tuple[int, int]] = []
    in_idr = False
    start_pos = 0

    for i, s in enumerate(scores):
        if s > threshold and not in_idr:
            in_idr = True
            start_pos = i
        elif s <= threshold and in_idr:
            if i - start_pos >= min_length:
                idr_regions.append((start_pos, i))
            in_idr = False

    # Handle IDR that extends to the end of the sequence
    if in_idr and len(scores) - start_pos >= min_length:
        idr_regions.append((start_pos, len(scores)))

    return idr_regions


# ---------------------------------------------------------------------------
# Main entry points
# ---------------------------------------------------------------------------

def fetch_structure(uniprot_id: str, work_dir: Path) -> tuple[Path, str]:
    """Obtain a PDB file for *uniprot_id* and save it under *work_dir*.

    The function also stores ``pdb_info`` as a module-level attribute on the
    returned path object (``pdb_path._pdb_info``) for callers that need
    experimental metadata. Use ``get_pdb_info(pdb_path)`` to retrieve it safely.

    Strategy order:
      1. RCSB PDB (experimental, best resolution) -- V5bis
      2. AlphaFold DB (pre-computed)
      3. ESMFold (on-demand folding)
      4. UniProt PDB cross-ref (legacy fallback)
      5. Mock PDB (offline development)

    Parameters
    ----------
    uniprot_id : str
        UniProt accession (e.g. ``"P00533"``).
    work_dir : Path
        Directory where the PDB file will be written.

    Returns
    -------
    tuple[Path, str]
        (absolute path to the .pdb file, source name).
        Source is one of: "pdb_experimental", "alphafold", "esmfold", "mock".

    Raises
    ------
    RuntimeError
        If all retrieval strategies fail and mock is not possible.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = work_dir / f"{uniprot_id}.pdb"

    # If already downloaded in a previous run, reuse it
    if pdb_path.exists() and pdb_path.stat().st_size > 100:
        logger.info("PDB already cached: %s", pdb_path)
        # Try to detect source from cached pdb_info file
        cached_source = _load_cached_source(work_dir, uniprot_id)
        return pdb_path, cached_source or "alphafold"

    # ----- Strategy 1 (V5bis): RCSB PDB Search API -----
    pdb_info = _try_rcsb_pdb(uniprot_id, pdb_path, work_dir)
    if pdb_info is not None:
        _save_pdb_info(work_dir, uniprot_id, pdb_info, "pdb_experimental")
        logger.info(
            "Structure saved from RCSB PDB %s: %s (%d bytes)",
            pdb_info.get("pdb_id", "?"),
            pdb_path,
            pdb_path.stat().st_size,
        )
        return pdb_path, "pdb_experimental"

    # ----- Strategy 2: AlphaFold DB -----
    pdb_content = _try_alphafold(uniprot_id)
    if pdb_content is not None:
        pdb_path.write_text(pdb_content)
        _save_pdb_info(work_dir, uniprot_id, None, "alphafold")
        logger.info(
            "Structure saved from AlphaFold: %s (%d bytes)",
            pdb_path,
            pdb_path.stat().st_size,
        )
        return pdb_path, "alphafold"

    # ----- Strategy 3: ESMFold API -----
    pdb_content = _try_esmfold(uniprot_id)
    if pdb_content is not None:
        pdb_path.write_text(pdb_content)
        _save_pdb_info(work_dir, uniprot_id, None, "esmfold")
        logger.info(
            "Structure saved from ESMFold: %s (%d bytes)",
            pdb_path,
            pdb_path.stat().st_size,
        )
        return pdb_path, "esmfold"

    # ----- Strategy 4: UniProt PDB cross-ref -----
    pdb_content = _try_uniprot_pdb(uniprot_id)
    if pdb_content is not None:
        pdb_path.write_text(pdb_content)
        _save_pdb_info(work_dir, uniprot_id, None, "pdb_experimental")
        logger.info(
            "Structure saved from PDB: %s (%d bytes)",
            pdb_path,
            pdb_path.stat().st_size,
        )
        return pdb_path, "pdb_experimental"

    # ----- Strategy 5: Mock PDB (offline fallback) -----
    logger.warning(
        "All structure sources failed for %s, generating mock PDB", uniprot_id
    )
    mock_content = _generate_mock_pdb(uniprot_id)
    pdb_path.write_text(mock_content)
    _save_pdb_info(work_dir, uniprot_id, None, "mock")
    logger.info(
        "Mock structure saved: %s (%d bytes)", pdb_path, pdb_path.stat().st_size
    )
    return pdb_path, "mock"


def fetch_structure_from_sequence(sequence: str, work_dir: Path) -> tuple[Path, str]:
    """Obtain a PDB file by folding a raw amino acid sequence with ESMFold.

    Parameters
    ----------
    sequence : str
        Raw amino acid sequence (e.g. "MTEYKLVVV...").
    work_dir : Path
        Directory where the PDB file will be written.

    Returns
    -------
    tuple[Path, str]
        (absolute path to the .pdb file, source name).
        Source is "esmfold" or "mock".
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    # Use a hash of the sequence as filename
    seq_hash = hashlib.sha256(sequence.encode()).hexdigest()[:12]
    pdb_path = work_dir / f"seq_{seq_hash}.pdb"

    # If already folded in a previous run, reuse it
    if pdb_path.exists() and pdb_path.stat().st_size > 100:
        logger.info("Folded structure already cached: %s", pdb_path)
        return pdb_path, "esmfold"

    # ----- Strategy 1: ESMFold API with direct sequence -----
    try:
        logger.info("Folding sequence with ESMFold (%d residues) ...", len(sequence))
        resp = requests.post(
            "https://api.esmatlas.com/foldSequence/v1/pdb/",
            data=sequence,
            headers={"Content-Type": "text/plain"},
            timeout=(10, 300),
        )
        resp.raise_for_status()
        content = resp.text
        if "ATOM" in content:
            pdb_path.write_text(content)
            logger.info(
                "ESMFold folded structure saved: %s (%d bytes)",
                pdb_path,
                pdb_path.stat().st_size,
            )
            return pdb_path, "esmfold"
        else:
            logger.warning(
                "ESMFold response does not contain ATOM records for sequence"
            )
    except Exception as exc:
        logger.warning("ESMFold failed for direct sequence: %s", exc)

    # ----- Strategy 2: Mock PDB (offline fallback) -----
    logger.warning("ESMFold failed for sequence, generating mock PDB")
    mock_content = _generate_mock_pdb(f"SEQ_{seq_hash}")
    pdb_path.write_text(mock_content)
    logger.info(
        "Mock structure saved for sequence: %s (%d bytes)",
        pdb_path,
        pdb_path.stat().st_size,
    )
    return pdb_path, "mock"


def get_pdb_info(work_dir: Path, uniprot_id: str) -> Optional[dict]:
    """Retrieve stored PDB experimental metadata for a given job.

    Parameters
    ----------
    work_dir : Path
        Job working directory.
    uniprot_id : str
        UniProt accession.

    Returns
    -------
    dict or None
        The pdb_info dict (pdb_id, resolution, method, ligand_id, etc.)
        or None if not available.
    """
    import json

    info_path = work_dir / f"{uniprot_id}_pdb_info.json"
    if info_path.exists():
        try:
            data = json.loads(info_path.read_text())
            return data.get("pdb_info")
        except Exception:
            return None
    return None


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _try_rcsb_pdb(
    uniprot_id: str, pdb_path: Path, work_dir: Path
) -> Optional[dict]:
    """Try RCSB PDB search and download; return pdb_info on success or None."""
    try:
        pdb_info = query_rcsb_pdb(uniprot_id)
        if pdb_info is None:
            return None
        if _download_pdb_from_rcsb(pdb_info, pdb_path):
            return pdb_info
        return None
    except Exception as exc:
        logger.warning("RCSB PDB strategy failed for %s: %s", uniprot_id, exc)
        return None


def _save_pdb_info(
    work_dir: Path,
    uniprot_id: str,
    pdb_info: Optional[dict],
    source: str,
) -> None:
    """Persist pdb_info and source to a JSON sidecar file for later retrieval."""
    import json

    info_path = work_dir / f"{uniprot_id}_pdb_info.json"
    try:
        data = {"source": source, "pdb_info": pdb_info}
        info_path.write_text(json.dumps(data, indent=2))
    except Exception as exc:
        logger.debug("Failed to save pdb_info: %s", exc)


def _load_cached_source(work_dir: Path, uniprot_id: str) -> Optional[str]:
    """Load the cached structure source from the sidecar JSON file."""
    import json

    info_path = work_dir / f"{uniprot_id}_pdb_info.json"
    if info_path.exists():
        try:
            data = json.loads(info_path.read_text())
            return data.get("source")
        except Exception:
            pass
    return None


def _try_alphafold(uniprot_id: str) -> Optional[str]:
    """Fetch PDB from the AlphaFold Protein Structure Database."""
    url = ALPHAFOLD_API.format(uniprot_id=uniprot_id)
    try:
        logger.info("Trying AlphaFold DB for %s ...", uniprot_id)
        resp = requests.get(url, timeout=HTTP_TIMEOUT)
        resp.raise_for_status()
        data = resp.json()
        # The API returns a list; take the first entry
        if isinstance(data, list) and len(data) > 0:
            entry = data[0]
        else:
            entry = data

        pdb_url = entry.get("pdbUrl")
        if not pdb_url:
            logger.warning("AlphaFold response missing pdbUrl for %s", uniprot_id)
            return None

        logger.info("Downloading PDB from %s", pdb_url)
        pdb_resp = requests.get(pdb_url, timeout=HTTP_TIMEOUT)
        pdb_resp.raise_for_status()
        content = pdb_resp.text
        if len(content) < 100:
            logger.warning(
                "AlphaFold PDB suspiciously small (%d bytes)", len(content)
            )
            return None
        return content

    except requests.RequestException as exc:
        logger.warning("AlphaFold DB failed for %s: %s", uniprot_id, exc)
        return None
    except (KeyError, ValueError, IndexError) as exc:
        logger.warning("AlphaFold response parse error: %s", exc)
        return None


def _try_esmfold(uniprot_id: str) -> Optional[str]:
    """Fetch a folded structure from ESMFold API (requires a FASTA sequence)."""
    try:
        # First get the sequence from UniProt
        sequence = _fetch_uniprot_sequence(uniprot_id)
        if not sequence:
            logger.warning("Cannot run ESMFold without sequence for %s", uniprot_id)
            return None

        logger.info(
            "Trying ESMFold for %s (%d residues) ...", uniprot_id, len(sequence)
        )

        # ESMFold via the public API
        resp = requests.post(
            "https://api.esmatlas.com/foldSequence/v1/pdb/",
            data=sequence,
            headers={"Content-Type": "text/plain"},
            timeout=(10, 300),
        )
        resp.raise_for_status()
        content = resp.text
        if "ATOM" in content:
            logger.info("ESMFold returned valid PDB for %s", uniprot_id)
            return content
        logger.warning("ESMFold response does not contain ATOM records")
        return None

    except Exception as exc:
        logger.warning("ESMFold failed for %s: %s", uniprot_id, exc)
        return None


def _try_uniprot_pdb(uniprot_id: str) -> Optional[str]:
    """Try to find an experimental PDB via UniProt cross-references."""
    try:
        logger.info("Trying UniProt PDB cross-ref for %s ...", uniprot_id)
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        resp = requests.get(url, timeout=HTTP_TIMEOUT)
        resp.raise_for_status()
        data = resp.json()

        # Look for PDB cross-references
        for xref in data.get("uniProtKBCrossReferences", []):
            if xref.get("database") == "PDB":
                pdb_id = xref.get("id", "").lower()
                if pdb_id:
                    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
                    pdb_resp = requests.get(pdb_url, timeout=HTTP_TIMEOUT)
                    if pdb_resp.status_code == 200 and "ATOM" in pdb_resp.text:
                        logger.info(
                            "Got experimental PDB %s for %s", pdb_id, uniprot_id
                        )
                        return pdb_resp.text
        return None
    except Exception as exc:
        logger.warning(
            "UniProt PDB cross-ref failed for %s: %s", uniprot_id, exc
        )
        return None


def _fetch_uniprot_sequence(uniprot_id: str) -> Optional[str]:
    """Retrieve the amino-acid sequence from UniProt."""
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
        resp = requests.get(url, timeout=HTTP_TIMEOUT)
        resp.raise_for_status()
        lines = resp.text.strip().split("\n")
        # Skip header line(s) starting with '>'
        seq = "".join(line.strip() for line in lines if not line.startswith(">"))
        return seq if len(seq) > 10 else None
    except Exception as exc:
        logger.warning(
            "UniProt sequence fetch failed for %s: %s", uniprot_id, exc
        )
        return None


def _generate_mock_pdb(identifier: str) -> str:
    """Generate a minimal mock PDB file for offline development.

    Parameters
    ----------
    identifier : str
        A label to embed in the PDB HEADER.

    Returns
    -------
    str
        A valid-ish PDB string with a few ATOM records.
    """
    # Deterministic coordinates based on identifier hash
    h = int(hashlib.md5(identifier.encode()).hexdigest()[:8], 16)
    x_base = (h % 100) / 2.0
    y_base = ((h >> 8) % 100) / 2.0
    z_base = ((h >> 16) % 100) / 2.0

    lines = [
        f"HEADER    MOCK STRUCTURE FOR {identifier}",
        f"REMARK   1 Generated by DockIt (mock mode) for development/testing.",
        f"REMARK   2 Not a real structure.",
    ]
    # Generate a small alpha-helix-like set of atoms
    residues = [
        "ALA", "GLY", "VAL", "LEU", "ILE",
        "PRO", "PHE", "TRP", "MET", "SER",
        "THR", "CYS", "TYR", "ASN", "GLN",
        "ASP", "GLU", "LYS", "ARG", "HIS",
    ]
    atom_idx = 1
    for i in range(30):
        res = residues[i % len(residues)]
        res_num = i + 1
        # Helical coordinates
        import math

        angle = i * 100.0 * math.pi / 180.0
        rise = i * 1.5
        x = x_base + 5.0 * math.cos(angle)
        y = y_base + 5.0 * math.sin(angle)
        z = z_base + rise
        lines.append(
            f"ATOM  {atom_idx:5d}  CA  {res} A{res_num:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
        )
        atom_idx += 1
    lines.append("END")
    return "\n".join(lines) + "\n"
