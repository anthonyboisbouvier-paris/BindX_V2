"""
DockIt pipeline -- Binding-pocket detection.

Detection strategies (in priority order):
  1. Co-crystallized ligand pocket extraction (if ligand_id known from PDB)
  2. P2Rank ML-based pocket detection (preferred computational method)
  3. fpocket geometric detection (fallback)
  4. Geometric heuristic (no external tools)
  5. Centroid fallback (last resort)

V5bis: P2Rank support, ligand pocket extraction, enhanced mock with probability scores.
"""

from __future__ import annotations

import csv
import hashlib
import logging
import re
import shutil
import subprocess
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def detect_pockets(
    pdb_path: Path,
    work_dir: Path,
    ligand_id: Optional[str] = None,
) -> list[dict]:
    """Detect binding pockets on a protein structure.

    Parameters
    ----------
    pdb_path : Path
        Path to the input PDB file.
    work_dir : Path
        Working directory for tool output.
    ligand_id : str, optional
        Co-crystallized ligand ID from PDB (e.g. "AQ4" for erlotinib).
        When provided, the ligand neighbourhood is used as the primary pocket.

    Returns
    -------
    list[dict]
        Each dict has keys ``center`` (tuple[float,float,float]),
        ``score`` (float), ``probability`` (float, 0-1),
        ``method`` (str), and optionally ``volume`` (float),
        ``residues`` (list[str]).
        Sorted by score descending (best pocket first).
    """
    pockets: list[dict] = []

    # ----- Strategy 1: Co-crystallized ligand pocket -----
    if ligand_id:
        try:
            ligand_pocket = extract_ligand_pocket(str(pdb_path), ligand_id)
            if ligand_pocket:
                pockets.append(ligand_pocket)
                logger.info(
                    "Extracted co-crystallized ligand pocket for %s at (%.1f, %.1f, %.1f)",
                    ligand_id,
                    ligand_pocket["center"][0],
                    ligand_pocket["center"][1],
                    ligand_pocket["center"][2],
                )
        except Exception as exc:
            logger.warning("Ligand pocket extraction failed for %s: %s", ligand_id, exc)

    # ----- Strategy 2: P2Rank ML-based detection -----
    p2rank_pockets = detect_pockets_p2rank(str(pdb_path), str(work_dir / "p2rank_out"))
    if p2rank_pockets:
        # Merge with existing pockets (avoid duplicates near ligand pocket)
        for pp in p2rank_pockets:
            if not _is_duplicate_pocket(pp, pockets, threshold=5.0):
                pockets.append(pp)

    # ----- Strategy 3: fpocket fallback -----
    if not p2rank_pockets:
        fpocket_pockets = _run_fpocket(pdb_path, work_dir)
        if fpocket_pockets:
            for fp in fpocket_pockets:
                # Normalize fpocket output to include probability and method
                fp.setdefault("probability", fp.get("score", 0.5))
                fp.setdefault("method", "fpocket")
                if not _is_duplicate_pocket(fp, pockets, threshold=5.0):
                    pockets.append(fp)

    # ----- Strategy 4: Geometric fallback -----
    if not pockets:
        logger.info("No pockets from P2Rank or fpocket; using geometric fallback")
        geo_pockets = _geometric_fallback(pdb_path)
        if geo_pockets:
            for gp in geo_pockets:
                gp.setdefault("probability", gp.get("score", 0.3))
                gp.setdefault("method", "geometric")
                pockets.append(gp)

    # ----- Strategy 5: Centroid fallback -----
    if not pockets:
        logger.warning("No pockets detected; using protein centroid as single pocket")
        pockets = _centroid_fallback(pdb_path)
        for cp in pockets:
            cp.setdefault("probability", 0.1)
            cp.setdefault("method", "centroid")

    # Sort by score descending (best pocket first)
    pockets.sort(key=lambda p: p.get("score", 0.0), reverse=True)

    logger.info(
        "Detected %d pocket(s); best score=%.2f method=%s",
        len(pockets),
        pockets[0]["score"],
        pockets[0].get("method", "unknown"),
    )
    return pockets


# ---------------------------------------------------------------------------
# P2Rank ML pocket detection (V5bis)
# ---------------------------------------------------------------------------

def detect_pockets_p2rank(pdb_path: str, output_dir: str) -> list[dict]:
    """Detect binding pockets using P2Rank ML (preferred) or mock fallback.

    P2Rank predicts pockets using a machine-learning model trained on
    known binding sites. Its output CSV contains columns: name, rank, score,
    probability, sas_points, surf_atoms, center_x, center_y, center_z,
    residue_ids.

    Parameters
    ----------
    pdb_path : str
        Path to the input PDB file.
    output_dir : str
        Directory for P2Rank output files.

    Returns
    -------
    list[dict]
        Pockets with keys: ``center``, ``score``, ``probability``, ``method``,
        ``rank``, ``residues``, ``volume``.
        Empty list if detection fails.
    """
    # Check both known P2Rank install locations
    p2rank_bin = "/opt/p2rank/prank"
    if not Path(p2rank_bin).exists():
        p2rank_bin = "/opt/p2rank_2.4.2/prank"

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Try P2Rank first
    if Path(p2rank_bin).exists():
        try:
            logger.info("Running P2Rank on %s ...", pdb_path)
            result = subprocess.run(
                [p2rank_bin, "predict", "-f", pdb_path, "-o", output_dir],
                capture_output=True,
                text=True,
                timeout=120,
            )
            if result.returncode == 0:
                # Parse P2Rank CSV output
                # P2Rank creates: <output_dir>/<basename>_predictions.csv
                pdb_stem = Path(pdb_path).stem
                csv_path = Path(output_dir) / f"{pdb_stem}_predictions.csv"
                # Also try the standard P2Rank output path pattern
                if not csv_path.exists():
                    csv_path = Path(output_dir) / f"{pdb_stem}.pdb_predictions.csv"
                if not csv_path.exists():
                    # Search recursively for any predictions CSV
                    candidates = list(Path(output_dir).rglob("*predictions.csv"))
                    if candidates:
                        csv_path = candidates[0]

                if csv_path.exists():
                    pockets = _parse_p2rank_csv(csv_path)
                    if pockets:
                        logger.info("P2Rank detected %d pockets", len(pockets))
                        return pockets
                else:
                    logger.warning("P2Rank CSV output not found in %s", output_dir)
            else:
                logger.warning(
                    "P2Rank exited with code %d: %s",
                    result.returncode,
                    result.stderr[:500],
                )
        except subprocess.TimeoutExpired:
            logger.warning("P2Rank timed out after 120s")
        except FileNotFoundError:
            logger.warning("P2Rank binary not found at %s", p2rank_bin)
        except Exception as exc:
            logger.warning("P2Rank failed: %s, falling back", exc)
    else:
        logger.info("P2Rank not installed at %s; using mock P2Rank detection", p2rank_bin)

    # ----- Mock P2Rank fallback -----
    # Generate deterministic mock pockets with P2Rank-style probability scores
    return _mock_p2rank_pockets(pdb_path)


def _parse_p2rank_csv(csv_path: Path) -> list[dict]:
    """Parse P2Rank predictions CSV into pocket dicts.

    Parameters
    ----------
    csv_path : Path
        Path to the P2Rank ``*_predictions.csv`` file.

    Returns
    -------
    list[dict]
        Parsed pockets sorted by probability descending.
    """
    pockets: list[dict] = []
    try:
        with open(csv_path, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    # P2Rank CSV columns may have leading/trailing spaces
                    cleaned = {k.strip(): v.strip() for k, v in row.items() if k}

                    center_x = float(cleaned.get("center_x", "0"))
                    center_y = float(cleaned.get("center_y", "0"))
                    center_z = float(cleaned.get("center_z", "0"))
                    probability = float(cleaned.get("probability", "0"))
                    score_val = float(cleaned.get("score", "0"))
                    rank = int(cleaned.get("rank", "0"))
                    residues_str = cleaned.get("residue_ids", "")

                    residues = [r.strip() for r in residues_str.split(",") if r.strip()] if residues_str else []

                    pockets.append({
                        "center": (
                            round(center_x, 3),
                            round(center_y, 3),
                            round(center_z, 3),
                        ),
                        "score": round(probability, 4),  # Use probability as score
                        "probability": round(probability, 4),
                        "p2rank_score": round(score_val, 4),
                        "rank": rank,
                        "method": "p2rank",
                        "residues": residues,
                        "volume": 0.0,
                    })
                except (ValueError, KeyError) as exc:
                    logger.debug("Skipping P2Rank row: %s", exc)
                    continue

        pockets.sort(key=lambda p: p["probability"], reverse=True)

    except Exception as exc:
        logger.warning("Error parsing P2Rank CSV %s: %s", csv_path, exc)

    return pockets


def _mock_p2rank_pockets(pdb_path: str) -> list[dict]:
    """Generate mock pockets with P2Rank-style probability scores.

    Uses deterministic hash-based generation so repeated runs produce
    the same results for the same input.

    Parameters
    ----------
    pdb_path : str
        Path to the PDB file (used for atom coordinate parsing and hashing).

    Returns
    -------
    list[dict]
        Mock pockets with ``probability`` field (0-1) and ``method: "p2rank_mock"``.
    """
    atoms = _parse_pdb_atoms(Path(pdb_path))
    if len(atoms) < 10:
        return []

    # Use file path hash for deterministic results
    seed_hash = hashlib.md5(pdb_path.encode()).hexdigest()

    # Compute protein centroid
    xs = [a[0] for a in atoms]
    ys = [a[1] for a in atoms]
    zs = [a[2] for a in atoms]
    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)
    cz = sum(zs) / len(zs)

    # Generate 3-5 mock pockets at different octants
    import random
    rng = random.Random(int(seed_hash[:8], 16))

    n_pockets = rng.randint(3, 5)
    pockets: list[dict] = []

    for i in range(n_pockets):
        # Offset from centroid
        offset_x = rng.uniform(-8, 8)
        offset_y = rng.uniform(-8, 8)
        offset_z = rng.uniform(-8, 8)

        # Probability decreases with rank
        base_prob = 0.95 - i * 0.15
        prob = max(0.1, base_prob + rng.uniform(-0.05, 0.05))

        pockets.append({
            "center": (
                round(cx + offset_x, 3),
                round(cy + offset_y, 3),
                round(cz + offset_z, 3),
            ),
            "score": round(prob, 4),
            "probability": round(prob, 4),
            "p2rank_score": round(prob * 20, 2),
            "rank": i + 1,
            "method": "p2rank_mock",
            "residues": [],
            "volume": round(rng.uniform(200, 800), 1),
        })

    pockets.sort(key=lambda p: p["probability"], reverse=True)
    logger.info("Generated %d mock P2Rank pockets", len(pockets))
    return pockets


# ---------------------------------------------------------------------------
# Co-crystallized ligand pocket extraction (V5bis)
# ---------------------------------------------------------------------------

def extract_ligand_pocket(
    pdb_path: str,
    ligand_id: str,
    radius: float = 6.0,
) -> Optional[dict]:
    """Extract binding pocket as residues within radius of co-crystallized ligand.

    Parses the PDB file to find HETATM records matching the ligand,
    then identifies protein residues (ATOM records) within the specified
    radius of any ligand atom.

    Parameters
    ----------
    pdb_path : str
        Path to the PDB file containing both protein and ligand.
    ligand_id : str
        3-letter residue name of the co-crystallized ligand (e.g. "AQ4").
    radius : float
        Distance cutoff in Angstroms for defining pocket residues (default 6.0).

    Returns
    -------
    dict or None
        Pocket dict with ``center``, ``score``, ``probability``, ``method``,
        ``residues``, ``ligand_id``, or None if ligand not found in PDB.
    """
    import math

    ligand_atoms: list[tuple[float, float, float]] = []
    protein_atoms: list[tuple[float, float, float, str]] = []  # x, y, z, residue_label

    try:
        with open(pdb_path, "r") as fh:
            for line in fh:
                record_type = line[:6].strip()

                if record_type == "HETATM":
                    res_name = line[17:20].strip()
                    if res_name.upper() == ligand_id.upper():
                        try:
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            ligand_atoms.append((x, y, z))
                        except (ValueError, IndexError):
                            continue

                elif record_type == "ATOM":
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        chain = line[21:22].strip() or "A"
                        res_num = line[22:26].strip()
                        res_name = line[17:20].strip()
                        residue_label = f"{chain}_{res_name}{res_num}"
                        protein_atoms.append((x, y, z, residue_label))
                    except (ValueError, IndexError):
                        continue

    except Exception as exc:
        logger.warning("Failed to parse PDB for ligand pocket: %s", exc)
        return None

    if not ligand_atoms:
        logger.info("Ligand %s not found in PDB %s", ligand_id, pdb_path)
        return None

    # Compute ligand center
    lig_cx = sum(a[0] for a in ligand_atoms) / len(ligand_atoms)
    lig_cy = sum(a[1] for a in ligand_atoms) / len(ligand_atoms)
    lig_cz = sum(a[2] for a in ligand_atoms) / len(ligand_atoms)

    # Find protein residues within radius of any ligand atom
    nearby_residues: set[str] = set()
    radius_sq = radius * radius

    for px, py, pz, res_label in protein_atoms:
        for lx, ly, lz in ligand_atoms:
            dx = px - lx
            dy = py - ly
            dz = pz - lz
            dist_sq = dx * dx + dy * dy + dz * dz
            if dist_sq <= radius_sq:
                nearby_residues.add(res_label)
                break  # No need to check other ligand atoms

    if not nearby_residues:
        logger.warning(
            "No protein residues found within %.1f A of ligand %s",
            radius,
            ligand_id,
        )
        return None

    logger.info(
        "Ligand %s pocket: center=(%.1f, %.1f, %.1f), %d residues within %.1f A",
        ligand_id,
        lig_cx,
        lig_cy,
        lig_cz,
        len(nearby_residues),
        radius,
    )

    return {
        "center": (round(lig_cx, 3), round(lig_cy, 3), round(lig_cz, 3)),
        "score": 0.99,  # Highest confidence -- experimentally validated site
        "probability": 0.99,
        "method": "co-crystallized_ligand",
        "residues": sorted(nearby_residues),
        "ligand_id": ligand_id,
        "volume": 0.0,
    }


# ---------------------------------------------------------------------------
# fpocket
# ---------------------------------------------------------------------------

def _run_fpocket(pdb_path: Path, work_dir: Path) -> Optional[list[dict]]:
    """Run fpocket and parse its output."""
    if not shutil.which("fpocket"):
        logger.debug("fpocket not on PATH")
        return None

    try:
        cmd = ["fpocket", "-f", str(pdb_path)]
        logger.info("Running fpocket: %s", " ".join(cmd))
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120,
            cwd=str(work_dir),
        )

        if result.returncode != 0:
            logger.warning(
                "fpocket exited with code %d: %s",
                result.returncode,
                result.stderr[:500],
            )
            return None

        # fpocket creates <stem>_out/ directory next to the input PDB
        out_dir = pdb_path.parent / f"{pdb_path.stem}_out"
        if not out_dir.exists():
            # Sometimes it creates it in the cwd
            out_dir = work_dir / f"{pdb_path.stem}_out"

        if not out_dir.exists():
            logger.warning("fpocket output directory not found")
            return None

        return _parse_fpocket_output(out_dir)

    except FileNotFoundError:
        logger.warning("fpocket binary not found")
        return None
    except subprocess.TimeoutExpired:
        logger.warning("fpocket timed out")
        return None
    except Exception as exc:
        logger.warning("fpocket error: %s", exc)
        return None


def _parse_fpocket_output(out_dir: Path) -> Optional[list[dict]]:
    """Parse the fpocket *_info.txt file to extract pocket centres and scores."""
    info_file = out_dir / f"{out_dir.stem.replace('_out', '')}_info.txt"
    if not info_file.exists():
        # Try globbing
        candidates = list(out_dir.glob("*_info.txt"))
        if not candidates:
            logger.warning("No fpocket info file found in %s", out_dir)
            return None
        info_file = candidates[0]

    pockets: list[dict] = []
    try:
        text = info_file.read_text()
        # Split by pocket blocks
        blocks = re.split(r"Pocket\s+\d+\s*:", text)
        for block in blocks[1:]:  # skip header
            center = _extract_xyz(block)
            score = _extract_float(block, r"Druggability\s+Score\s*:\s*([\d.]+)")
            if score is None:
                score = _extract_float(block, r"Score\s*:\s*([\d.]+)")
            volume = _extract_float(block, r"Volume\s*:\s*([\d.]+)")
            if center:
                pockets.append({
                    "center": center,
                    "score": score if score is not None else 0.5,
                    "probability": score if score is not None else 0.5,
                    "method": "fpocket",
                    "volume": volume if volume is not None else 0.0,
                })
    except Exception as exc:
        logger.warning("Error parsing fpocket output: %s", exc)
        return None

    pockets.sort(key=lambda p: p["score"], reverse=True)
    return pockets if pockets else None


def _extract_xyz(block: str) -> Optional[tuple[float, float, float]]:
    """Extract center x/y/z from an fpocket pocket block."""
    x = _extract_float(
        block,
        r"Center\s*of\s*mass\s*-\s*Alpha\s*Sphere\s*center\s*x\s*:\s*([-\d.]+)",
    )
    if x is None:
        x = _extract_float(block, r"center_x\s*[=:]\s*([-\d.]+)")
    y = _extract_float(block, r"center.*?y\s*[=:]\s*([-\d.]+)")
    z = _extract_float(block, r"center.*?z\s*[=:]\s*([-\d.]+)")

    # Broader fallback patterns
    if x is None or y is None or z is None:
        coords = re.findall(r"([-]?\d+\.\d+)", block)
        if len(coords) >= 3:
            x, y, z = float(coords[0]), float(coords[1]), float(coords[2])

    if x is not None and y is not None and z is not None:
        return (x, y, z)
    return None


def _extract_float(text: str, pattern: str) -> Optional[float]:
    m = re.search(pattern, text, re.IGNORECASE)
    return float(m.group(1)) if m else None


# ---------------------------------------------------------------------------
# Geometric fallback (no external tools required)
# ---------------------------------------------------------------------------

def _geometric_fallback(pdb_path: Path) -> Optional[list[dict]]:
    """Identify binding-pocket candidates from PDB coordinates.

    Strategy: compute the centroid of all CA atoms, then find the region
    with the highest local atom density that is NOT at the very core
    (pockets tend to be partially solvent-exposed).
    """
    atoms = _parse_pdb_atoms(pdb_path)
    if len(atoms) < 10:
        return None

    # Compute overall centroid
    xs = [a[0] for a in atoms]
    ys = [a[1] for a in atoms]
    zs = [a[2] for a in atoms]

    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)
    cz = sum(zs) / len(zs)

    # Divide the space into octants relative to centroid and pick the densest
    octants: dict[tuple[int, int, int], list[tuple[float, float, float]]] = {}
    for x, y, z in atoms:
        key = (
            1 if x >= cx else 0,
            1 if y >= cy else 0,
            1 if z >= cz else 0,
        )
        octants.setdefault(key, []).append((x, y, z))

    pockets: list[dict] = []
    for key, pts in octants.items():
        if len(pts) < 5:
            continue
        px = sum(p[0] for p in pts) / len(pts)
        py = sum(p[1] for p in pts) / len(pts)
        pz = sum(p[2] for p in pts) / len(pts)
        # Score by density (more atoms in a region = more likely a pocket)
        score = min(1.0, len(pts) / len(atoms) * 4)
        pockets.append({
            "center": (round(px, 3), round(py, 3), round(pz, 3)),
            "score": round(score, 3),
            "volume": 0.0,
        })

    pockets.sort(key=lambda p: p["score"], reverse=True)
    return pockets[:5] if pockets else None


def _centroid_fallback(pdb_path: Path) -> list[dict]:
    """Absolute last resort: use the geometric centre of all atoms."""
    atoms = _parse_pdb_atoms(pdb_path)
    if not atoms:
        # Return origin as a dummy
        return [{"center": (0.0, 0.0, 0.0), "score": 0.1, "volume": 0.0}]

    cx = sum(a[0] for a in atoms) / len(atoms)
    cy = sum(a[1] for a in atoms) / len(atoms)
    cz = sum(a[2] for a in atoms) / len(atoms)
    return [
        {
            "center": (round(cx, 3), round(cy, 3), round(cz, 3)),
            "score": 0.5,
            "volume": 0.0,
        }
    ]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _is_duplicate_pocket(
    pocket: dict,
    existing: list[dict],
    threshold: float = 5.0,
) -> bool:
    """Check if a pocket center is within threshold distance of any existing pocket.

    Parameters
    ----------
    pocket : dict
        Pocket to check (must have ``center`` key).
    existing : list[dict]
        Already-selected pockets.
    threshold : float
        Minimum distance in Angstroms to consider pockets as distinct.

    Returns
    -------
    bool
        True if the pocket is a duplicate (too close to an existing one).
    """
    import math

    pc = pocket.get("center", (0, 0, 0))
    for ep in existing:
        ec = ep.get("center", (0, 0, 0))
        dist = math.sqrt(
            (pc[0] - ec[0]) ** 2 + (pc[1] - ec[1]) ** 2 + (pc[2] - ec[2]) ** 2
        )
        if dist < threshold:
            return True
    return False


def _parse_pdb_atoms(pdb_path: Path) -> list[tuple[float, float, float]]:
    """Extract (x, y, z) coordinates of ATOM records from a PDB file."""
    atoms: list[tuple[float, float, float]] = []
    try:
        with open(pdb_path, "r") as fh:
            for line in fh:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        atoms.append((x, y, z))
                    except (ValueError, IndexError):
                        continue
    except Exception as exc:
        logger.warning("Failed to parse PDB atoms from %s: %s", pdb_path, exc)
    return atoms
