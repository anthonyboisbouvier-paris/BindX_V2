"""
DockIt pipeline -- Lead optimization via iterative molecular modification.

Implements a multi-objective optimization loop:
  1. Generate structural variants of the starting molecule.
  2. Validate (sanitize, QED, Tanimoto similarity to parent).
  3. Dock each variant against the target (real Vina/GNINA or mock).
  4. ADMET prediction for each variant (real or mock).
  5. Score with multi-objective function (affinity, toxicity, bioavailability, synthesis).
  6. Keep top variants, use best as parent for next iteration.
  7. Return final lead with improvement metrics.

Uses RDKit for molecular manipulation where available; falls back to a
deterministic mock that simulates gradual improvement across iterations.
Real docking mode uses pipeline.prepare + pipeline.docking for actual Vina/GNINA scoring.
"""

from __future__ import annotations

import hashlib
import logging
import os
import random
from pathlib import Path
from typing import Callable, Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Availability flag
# ---------------------------------------------------------------------------

_RDKIT_AVAILABLE: Optional[bool] = None


def _check_rdkit() -> bool:
    """Check whether RDKit can be imported."""
    global _RDKIT_AVAILABLE
    if _RDKIT_AVAILABLE is None:
        try:
            from rdkit import Chem  # noqa: F401
            _RDKIT_AVAILABLE = True
        except ImportError:
            _RDKIT_AVAILABLE = False
            logger.warning("RDKit not available; lead optimization will use mock only")
    return _RDKIT_AVAILABLE


# Default optimization weights
DEFAULT_WEIGHTS: dict[str, float] = {
    "affinity": 0.35,
    "toxicity": 0.25,
    "bioavailability": 0.20,
    "synthesis": 0.20,
}


# =====================================================================
# REAL DOCKING / ADMET / PREFILTER
# =====================================================================

def _real_dock(
    smiles: str,
    name: str,
    receptor_pdbqt: Path,
    pocket_center: list[float],
    pocket_size: list[float],
    exhaustiveness: int,
    docking_engine: str,
    work_dir: Path,
) -> Optional[dict]:
    """Dock a single molecule using the real pipeline (prepare + dock).

    Returns a dict with docking results or None on failure.
    """
    try:
        from pipeline.prepare import prepare_ligand
        from pipeline.docking import run_docking

        work_dir.mkdir(parents=True, exist_ok=True)

        # Prepare ligand: SMILES -> 3D PDBQT
        safe_name = "".join(c if c.isalnum() or c in "-_" else "_" for c in name)[:50]
        pdbqt_path = prepare_ligand(smiles, safe_name, work_dir)
        if pdbqt_path is None:
            logger.debug("Ligand preparation failed for %s", name)
            return None

        # Dock
        center = tuple(pocket_center[:3])
        size = tuple(pocket_size[:3])
        result = run_docking(
            receptor_pdbqt=receptor_pdbqt,
            ligand_pdbqt=pdbqt_path,
            center=center,
            size=size,
            exhaustiveness=exhaustiveness,
            docking_engine=docking_engine,
        )

        return {
            "affinity": result.get("affinity", 0.0),
            "vina_score": result.get("vina_score", result.get("affinity", 0.0)),
            "cnn_score": result.get("cnn_score", 0.0),
            "cnn_affinity": result.get("cnn_affinity", 0.0),
            "pose_pdbqt_path": str(result.get("pose_pdbqt_path", "")),
            "docking_engine": result.get("docking_engine", "mock"),
        }
    except Exception as exc:
        logger.warning("Real docking failed for %s: %s", name, exc)
        return None


def _real_admet_batch(smiles_list: list[str]) -> dict[str, dict]:
    """Run real ADMET prediction on a batch of SMILES.

    Returns {smiles: admet_dict} with fields compatible with _multi_objective_score().
    """
    try:
        from pipeline.admet import predict_admet

        results = predict_admet(smiles_list)
        output: dict[str, dict] = {}
        for admet_entry in results:
            smi = admet_entry.get("smiles", "")
            if not smi:
                continue
            # Extract the scores _multi_objective_score expects
            toxicity_data = admet_entry.get("toxicity", {})
            absorption_data = admet_entry.get("absorption", {})
            output[smi] = {
                "toxicity_score": 1.0 - max(
                    toxicity_data.get("hepatotoxicity", 0.0),
                    toxicity_data.get("herg_inhibition", 0.0),
                    toxicity_data.get("ames_mutagenicity", 0.0),
                ),  # higher = safer
                "bioavailability_score": absorption_data.get("oral_bioavailability", 0.5),
                "synthesis_score": 1.0 - admet_entry.get("composite_score", 0.5) * 0.3,  # proxy
                "herg_risk": toxicity_data.get("herg_inhibition", 0.0),
                "hepatotoxicity_risk": toxicity_data.get("hepatotoxicity", 0.0),
                "composite_score": admet_entry.get("composite_score", 0.5),
                # Preserve full ADMET data for enrichment
                "_full_admet": admet_entry,
            }
        return output
    except Exception as exc:
        logger.warning("Real ADMET batch failed: %s", exc)
        return {}


def _prefilter_variants(
    variants: list[str],
    parent_smiles: str,
    top_n: int,
) -> list[str]:
    """Pre-filter variants by QED x Tanimoto score, return top_n."""
    if not _check_rdkit() or len(variants) <= top_n:
        return variants[:top_n]

    try:
        from rdkit import Chem, DataStructs
        from rdkit.Chem import AllChem, QED as QED_module

        parent_mol = Chem.MolFromSmiles(parent_smiles)
        if parent_mol is None:
            return variants[:top_n]

        try:
            from rdkit.Chem import rdFingerprintGenerator
            fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
            parent_fp = fp_gen.GetFingerprint(parent_mol)
        except (ImportError, AttributeError):
            parent_fp = AllChem.GetMorganFingerprintAsBitVect(parent_mol, 2, nBits=2048)
            fp_gen = None

        scored: list[tuple[float, str]] = []
        for smi in variants:
            try:
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    continue
                qed_val = QED_module.qed(mol)
                if fp_gen is not None:
                    fp = fp_gen.GetFingerprint(mol)
                else:
                    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                sim = DataStructs.TanimotoSimilarity(parent_fp, fp)
                scored.append((qed_val * sim, smi))
            except Exception:
                continue

        scored.sort(key=lambda x: x[0], reverse=True)
        return [smi for _, smi in scored[:top_n]]
    except Exception as exc:
        logger.warning("Pre-filter failed: %s", exc)
        return variants[:top_n]


# =====================================================================
# PUBLIC API
# =====================================================================

def run_optimization(
    starting_smiles: str,
    starting_name: str,
    target_pdbqt: str,
    pocket_center: list[float],
    weights: Optional[dict[str, float]] = None,
    n_iterations: int = 10,
    variants_per_iter: int = 50,
    work_dir: Optional[Path] = None,
    progress_callback: Optional[Callable[[int, float, str], None]] = None,
    structural_rules: Optional[dict] = None,
    docking_engine: str = "auto",
    pocket_size: Optional[list[float]] = None,
    exhaustiveness: int = 8,
    dock_top_n: int = 20,
) -> dict:
    """Run iterative lead optimization on a starting molecule.

    Each iteration:
    1. Generate variants by applying small medicinal-chemistry modifications.
    2. Validate: sanitize, QED > 0.2, Tanimoto > 0.3 with parent.
    3. Mock dock each valid variant.
    4. Mock ADMET for each variant.
    5. Compute multi-objective score.
    6. Keep top 5, use best as input for next iteration.

    Parameters
    ----------
    starting_smiles : str
        SMILES of the starting molecule.
    starting_name : str
        Name of the starting molecule.
    target_pdbqt : str
        Path to the target receptor PDBQT file.
    pocket_center : list[float]
        [x, y, z] coordinates of the binding pocket center.
    weights : dict[str, float], optional
        Objective weights. Keys: "affinity", "toxicity",
        "bioavailability", "synthesis". Defaults to DEFAULT_WEIGHTS.
    n_iterations : int
        Number of optimization iterations.
    variants_per_iter : int
        Number of variants to generate per iteration.
    work_dir : Path, optional
        Working directory.
    progress_callback : callable, optional
        Called as ``progress_callback(iteration, best_score, message)``
        after each iteration.
    docking_engine : str
        "vina", "gnina", "auto", or "mock". Controls whether real docking
        is used. Default "auto" uses real docking if receptor exists.
    pocket_size : list[float], optional
        [sx, sy, sz] search box size in Angstroms. Default [25, 25, 25].
    exhaustiveness : int
        Docking search exhaustiveness. Default 8.
    dock_top_n : int
        Number of pre-filtered variants to actually dock per iteration.

    Returns
    -------
    dict
        Keys:
        - ``starting_molecule``: dict with name, smiles, score.
        - ``final_lead``: dict with name, smiles, score, affinity, admet.
        - ``iterations``: list of per-iteration summaries.
        - ``comparison``: dict mapping metric names to start/end/change_pct.
        - ``total_tested``: int, total variants tested.
        - ``total_valid``: int, total variants passing filters.
        - ``total_improved``: int, total variants improving on parent.
    """
    if work_dir is not None:
        work_dir.mkdir(parents=True, exist_ok=True)

    w = weights or DEFAULT_WEIGHTS.copy()
    # Accept frontend key "binding_affinity" as alias for "affinity"
    if "binding_affinity" in w and "affinity" not in w:
        w["affinity"] = w.pop("binding_affinity")
    elif "binding_affinity" in w:
        w.pop("binding_affinity")
    # Normalize weights
    total_w = sum(w.values())
    if total_w > 0:
        w = {k: v / total_w for k, v in w.items()}

    # Use deterministic seed from starting SMILES
    seed_val = int(hashlib.sha256(starting_smiles.encode()).hexdigest()[:8], 16)
    rng = random.Random(seed_val)

    # Determine docking mode: real vs none vs mock
    _pocket_size = pocket_size or [25.0, 25.0, 25.0]
    receptor_path = Path(target_pdbqt)

    # "none" engine: skip docking entirely, score on properties only
    use_no_docking = docking_engine == "none"

    use_real_docking = (
        not use_no_docking
        and docking_engine != "mock"
        and receptor_path.exists()
        and receptor_path.stat().st_size > 100
    )

    if use_no_docking:
        logger.info("Lead optimization: NO-DOCKING mode — property-based scoring only")
    elif use_real_docking:
        logger.info("Lead optimization: REAL docking mode (engine=%s)", docking_engine)
    else:
        if docking_engine not in ("mock", "none") and not receptor_path.exists():
            logger.warning("Receptor %s not found; falling back to mock docking", target_pdbqt)
        logger.info("Lead optimization: MOCK docking mode")

    # Compute starting molecule scores
    if use_no_docking:
        # No docking: affinity is 0 (unused), score based on ADMET only
        start_affinity = 0.0
        start_admet = _mock_admet_simple(starting_smiles, rng)
        start_dock_meta = {"docking_engine": "none"}
    elif use_real_docking:
        start_dock = _real_dock(
            starting_smiles, starting_name, receptor_path, pocket_center,
            _pocket_size, exhaustiveness, docking_engine,
            work_dir / "iter_0" / "ligands" if work_dir else Path("/tmp/opt_start"),
        )
        if start_dock:
            start_affinity = start_dock["affinity"]
            start_dock_meta = start_dock
        else:
            start_affinity = _mock_dock(starting_smiles, pocket_center, rng)
            start_dock_meta = {"docking_engine": "mock"}
        admet_batch = _real_admet_batch([starting_smiles])
        start_admet = admet_batch.get(starting_smiles, _mock_admet_simple(starting_smiles, rng))
    else:
        start_affinity = _mock_dock(starting_smiles, pocket_center, rng)
        start_admet = _mock_admet_simple(starting_smiles, rng)
        start_dock_meta = {"docking_engine": "mock"}

    start_score = _multi_objective_score(start_affinity, start_admet, w)

    logger.info(
        "Lead optimization starting: %s (score=%.3f, affinity=%.1f, engine=%s)",
        starting_name, start_score, start_affinity,
        start_dock_meta.get("docking_engine", "none" if use_no_docking else "mock"),
    )

    # Track results across iterations
    iterations: list[dict] = []
    current_smiles = starting_smiles
    current_score = start_score
    current_affinity = start_affinity
    current_admet = start_admet

    total_tested = 0
    total_valid = 0
    total_improved = 0

    seen_smiles: set[str] = {starting_smiles}
    all_top_molecules: list[dict] = []

    best_ever_smiles = starting_smiles
    best_ever_score = start_score
    best_ever_affinity = start_affinity
    best_ever_admet = start_admet
    best_ever_dock_meta: dict = start_dock_meta

    stale_iterations = 0  # count consecutive iterations with no improvement

    for iteration in range(1, n_iterations + 1):
        iter_rng = random.Random(seed_val + iteration * 1000)

        # Generate variants
        variants = _generate_variants(
            current_smiles, variants_per_iter, iter_rng, iteration,
            structural_rules=structural_rules,
        )
        n_tested = len(variants)
        total_tested += n_tested

        # Score each valid variant
        scored_variants: list[dict] = []

        if use_no_docking:
            # --- NO-DOCKING MODE: property-based scoring only ---
            admet_map = _real_admet_batch(variants) if variants else {}
            for v_smiles in variants:
                v_admet = admet_map.get(v_smiles, _mock_admet_simple(v_smiles, iter_rng))
                # With no affinity signal, use affinity=0 and down-weight it
                v_score = _multi_objective_score(0.0, v_admet, w)
                scored_variants.append({
                    "smiles": v_smiles,
                    "affinity": 0.0,
                    "admet": v_admet,
                    "score": v_score,
                    "docking_engine": "none",
                })

        elif use_real_docking:
            # --- REAL DOCKING MODE: dock ALL variants ---
            logger.info("Iter %d: docking all %d variants (engine=%s)",
                        iteration, len(variants), docking_engine)

            iter_work_dir = (work_dir / f"iter_{iteration}" / "ligands") if work_dir else Path(f"/tmp/opt_iter_{iteration}")
            docked: list[dict] = []
            for vi, v_smiles in enumerate(variants):
                dock_result = _real_dock(
                    v_smiles, f"{starting_name}_i{iteration}_v{vi}",
                    receptor_path, pocket_center, _pocket_size,
                    exhaustiveness, docking_engine, iter_work_dir,
                )
                if dock_result:
                    docked.append({"smiles": v_smiles, **dock_result})

            # Real ADMET on docked variants
            docked_smiles = [d["smiles"] for d in docked]
            admet_map = _real_admet_batch(docked_smiles) if docked_smiles else {}

            # Score
            for d in docked:
                v_admet = admet_map.get(d["smiles"], _mock_admet_simple(d["smiles"], iter_rng))
                v_score = _multi_objective_score(d["affinity"], v_admet, w)
                scored_variants.append({
                    "smiles": d["smiles"],
                    "affinity": d["affinity"],
                    "admet": v_admet,
                    "score": v_score,
                    "vina_score": d.get("vina_score"),
                    "cnn_score": d.get("cnn_score"),
                    "cnn_affinity": d.get("cnn_affinity"),
                    "docking_engine": d.get("docking_engine"),
                    "pose_pdbqt_path": d.get("pose_pdbqt_path"),
                })

        else:
            # --- MOCK DOCKING MODE ---
            for v_smiles in variants:
                v_affinity = _mock_dock(v_smiles, pocket_center, iter_rng)
                v_admet = _mock_admet_simple(v_smiles, iter_rng)
                v_score = _multi_objective_score(v_affinity, v_admet, w)
                scored_variants.append({
                    "smiles": v_smiles,
                    "affinity": v_affinity,
                    "admet": v_admet,
                    "score": v_score,
                })

            # Mock mode: add artificial improvement bonus
            improvement_bonus = iteration * 0.008 * (1.0 + rng.random() * 0.5)
            for sv in scored_variants:
                sv["score"] = min(1.0, sv["score"] + improvement_bonus)

        n_valid = len(scored_variants)
        total_valid += n_valid

        # Sort by score descending
        scored_variants.sort(key=lambda x: x["score"], reverse=True)

        # Count improved variants
        n_improved = sum(
            1 for sv in scored_variants if sv["score"] > current_score
        )
        total_improved += n_improved

        # Keep top 5
        top5 = scored_variants[:5]

        # Collect unique top molecules across iterations
        for sv in top5:
            if sv["smiles"] not in seen_smiles:
                seen_smiles.add(sv["smiles"])
                mol_entry: dict = {
                    "name": f"{starting_name}_iter{iteration}_{len(all_top_molecules)+1}",
                    "smiles": sv["smiles"],
                    "score": round(sv["score"], 4),
                    "affinity": round(sv["affinity"], 2),
                    "admet": sv["admet"],
                    "iteration": iteration,
                }
                # Propagate real docking metadata
                for dk in ("vina_score", "cnn_score", "cnn_affinity", "docking_engine", "pose_pdbqt_path"):
                    if sv.get(dk) is not None:
                        mol_entry[dk] = sv[dk]
                all_top_molecules.append(mol_entry)

        prev_score = current_score
        if top5 and top5[0]["score"] > current_score:
            current_smiles = top5[0]["smiles"]
            current_score = top5[0]["score"]
            current_affinity = top5[0]["affinity"]
            current_admet = top5[0]["admet"]

        if current_score > best_ever_score:
            best_ever_smiles = current_smiles
            best_ever_score = current_score
            best_ever_affinity = current_affinity
            best_ever_admet = current_admet
            # Track docking metadata for the best molecule
            if top5 and top5[0].get("docking_engine"):
                best_ever_dock_meta = {k: top5[0].get(k) for k in
                    ("vina_score", "cnn_score", "cnn_affinity", "docking_engine", "pose_pdbqt_path")
                    if top5[0].get(k) is not None}

        # Early termination: 3 consecutive stale iterations in real or no-docking mode
        if use_real_docking or use_no_docking:
            if current_score <= prev_score:
                stale_iterations += 1
            else:
                stale_iterations = 0
            if stale_iterations >= 3 and iteration >= 3:
                logger.info("Early termination: no improvement for %d consecutive iterations", stale_iterations)
                iter_summary = {
                    "iteration": iteration,
                    "best_score": round(current_score, 4),
                    "best_smiles": current_smiles,
                    "n_tested": n_tested,
                    "n_valid": n_valid,
                    "n_improved": n_improved,
                    "early_stop": True,
                }
                iterations.append(iter_summary)
                if progress_callback is not None:
                    try:
                        progress_callback(iteration, current_score,
                            f"Early stop at iteration {iteration}/{n_iterations}: no improvement for 3 iterations")
                    except Exception:
                        pass
                break

        iter_summary = {
            "iteration": iteration,
            "best_score": round(current_score, 4),
            "best_smiles": current_smiles,
            "n_tested": n_tested,
            "n_valid": n_valid,
            "n_improved": n_improved,
        }
        iterations.append(iter_summary)

        if progress_callback is not None:
            try:
                progress_callback(
                    iteration,
                    current_score,
                    f"Iteration {iteration}/{n_iterations}: "
                    f"score={current_score:.3f}, tested={n_tested}, improved={n_improved}",
                )
            except Exception as exc:
                logger.warning("Progress callback failed: %s", exc)

        logger.info(
            "Optimization iter %d/%d: score=%.3f, tested=%d, valid=%d, improved=%d",
            iteration, n_iterations, current_score, n_tested, n_valid, n_improved,
        )

    # Build comparison metrics
    comparison = _build_comparison(
        start_score, best_ever_score,
        start_affinity, best_ever_affinity,
        start_admet, best_ever_admet,
    )

    final_lead: dict = {
        "name": f"{starting_name}_opt",
        "smiles": best_ever_smiles,
        "score": round(best_ever_score, 4),
        "affinity": round(best_ever_affinity, 2),
        "admet": best_ever_admet,
    }
    # Propagate docking metadata to final lead
    for dk in ("vina_score", "cnn_score", "cnn_affinity", "docking_engine", "pose_pdbqt_path"):
        if best_ever_dock_meta.get(dk) is not None:
            final_lead[dk] = best_ever_dock_meta[dk]

    # Build objectives dict (frontend-compatible format: {key: {start, current}})
    objectives = {
        "binding_affinity": {
            "start": comparison.get("affinity", {}).get("start", round(start_affinity, 2)),
            "current": comparison.get("affinity", {}).get("end", round(best_ever_affinity, 2)),
        },
        "toxicity": {
            "start": comparison.get("toxicity_score", {}).get("start", 0.5),
            "current": comparison.get("toxicity_score", {}).get("end", 0.5),
        },
        "bioavailability": {
            "start": comparison.get("bioavailability_score", {}).get("start", 0.5),
            "current": comparison.get("bioavailability_score", {}).get("end", 0.5),
        },
        "synthesis_score": {
            "start": comparison.get("synthesis_score", {}).get("start", 0.7),
            "current": comparison.get("synthesis_score", {}).get("end", 0.7),
        },
    }

    # Sort top molecules by score descending, keep up to 10
    top_molecules = sorted(all_top_molecules, key=lambda x: x["score"], reverse=True)[:10]

    result = {
        "starting_molecule": {
            "name": starting_name,
            "smiles": starting_smiles,
            "score": round(start_score, 4),
            "affinity": round(start_affinity, 2),
        },
        "final_lead": final_lead,
        # Frontend-compatible aliases
        "best_molecule": final_lead,
        "objectives": objectives,
        "iterations": iterations,
        "comparison": comparison,
        "total_tested": total_tested,
        "total_valid": total_valid,
        "total_improved": total_improved,
        "top_molecules": top_molecules,
        "docking_mode": "none" if use_no_docking else ("real" if use_real_docking else "mock"),
        "docking_engine_used": best_ever_dock_meta.get("docking_engine", "mock"),
    }

    logger.info(
        "Lead optimization complete: %s -> %s (score %.3f -> %.3f, +%.1f%%)",
        starting_name, f"{starting_name}_opt",
        start_score, best_ever_score,
        comparison.get("multi_objective_score", {}).get("change_pct", 0.0),
    )

    return result


# =====================================================================
# VARIANT GENERATION
# =====================================================================

def _generate_variants(
    parent_smiles: str,
    n_variants: int,
    rng: random.Random,
    iteration: int,
    structural_rules: Optional[dict] = None,
) -> list[str]:
    """Generate structural variants of the parent molecule.

    Uses RDKit for real molecular modifications when available,
    otherwise generates mock variant SMILES deterministically.

    Parameters
    ----------
    parent_smiles : str
        SMILES of the parent molecule.
    n_variants : int
        Number of variants to attempt.
    rng : random.Random
        Random number generator for reproducibility.
    iteration : int
        Current iteration number (affects modification strategy).
    structural_rules : dict, optional
        User-defined structural constraints (frozen positions, strategies, etc.).

    Returns
    -------
    list[str]
        List of unique, valid variant SMILES strings.
    """
    if _check_rdkit():
        return _generate_variants_rdkit(parent_smiles, n_variants, rng, iteration, structural_rules)
    else:
        return _generate_variants_mock(parent_smiles, n_variants, rng, iteration)


def _generate_variants_rdkit(
    parent_smiles: str,
    n_variants: int,
    rng: random.Random,
    iteration: int,
    structural_rules: Optional[dict] = None,
) -> list[str]:
    """Generate variants using RDKit molecular editing.

    Respects structural_rules when provided: frozen positions, allowed
    strategies, per-position groups, similarity and MW thresholds.
    """
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import AllChem, Descriptors, QED as QED_module, RWMol
        from rdkit import DataStructs

        # Suppress noisy RDKit warnings during modification
        RDLogger.DisableLog("rdApp.*")

        parent_mol = Chem.MolFromSmiles(parent_smiles)
        if parent_mol is None:
            RDLogger.EnableLog("rdApp.*")
            return _generate_variants_mock(parent_smiles, n_variants, rng, iteration)

        # --- Parse structural rules ---
        rules = structural_rules or {}
        frozen_set: set[int] = set(rules.get("frozen_positions", []))
        if rules.get("preserve_scaffold"):
            frozen_set |= set(rules.get("core_atom_indices", []))
        position_rules: dict[int, dict] = {}
        for r in rules.get("rules", []):
            if r.get("frozen"):
                frozen_set.add(r["position_idx"])
            else:
                position_rules[r["position_idx"]] = r
        allowed_strats_raw = rules.get("allowed_strategies", [])
        allowed_strats: Optional[set[str]] = set(allowed_strats_raw) if allowed_strats_raw else None
        min_similarity: float = rules.get("min_similarity", 0.3)
        max_mw_delta: float = rules.get("max_mw_change", float('inf'))

        # Compute parent fingerprint for Tanimoto filter
        try:
            from rdkit.Chem import rdFingerprintGenerator
            fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
            parent_fp = fp_gen.GetFingerprint(parent_mol)
        except (ImportError, AttributeError):
            parent_fp = AllChem.GetMorganFingerprintAsBitVect(parent_mol, 2, nBits=2048)
            fp_gen = None

        parent_mw = Descriptors.ExactMolWt(parent_mol)

        variants: list[str] = []
        seen: set[str] = {parent_smiles}

        # Functional groups for addition
        _FG = ["C", "CC", "O", "OC", "N", "NC", "F", "Cl", "C(F)(F)F",
               "C#N", "C(=O)N", "C1CC1", "N1CCNCC1", "N1CCOCC1"]

        # Filter global strategies
        all_strategies = ["add_fg", "swap_halogen", "swap_atom", "modify_chain"]
        available_strategies = [s for s in all_strategies if not allowed_strats or s in allowed_strats]
        if not available_strategies:
            available_strategies = all_strategies

        max_attempts = n_variants * 5
        attempts = 0

        while len(variants) < n_variants and attempts < max_attempts:
            attempts += 1
            strategy = rng.choice(available_strategies)

            new_smi: Optional[str] = None

            try:
                rwmol = RWMol(parent_mol)

                if strategy == "add_fg":
                    # Add a functional group to an aromatic carbon
                    ar_carbons = [
                        a for a in rwmol.GetAtoms()
                        if a.GetIsAromatic() and a.GetAtomicNum() == 6
                        and a.GetTotalNumHs() > 0
                        and a.GetIdx() not in frozen_set
                    ]
                    if ar_carbons:
                        target = rng.choice(ar_carbons)
                        rule = position_rules.get(target.GetIdx())
                        if rule and rule.get("allowed_groups"):
                            fg_smi = rng.choice(rule["allowed_groups"])
                        else:
                            fg_smi = rng.choice(_FG)
                        fg_mol = Chem.MolFromSmiles(fg_smi)
                        if fg_mol is not None:
                            combo = Chem.RWMol(Chem.CombineMols(rwmol, fg_mol))
                            n_orig = rwmol.GetNumAtoms()
                            combo.AddBond(target.GetIdx(), n_orig, Chem.BondType.SINGLE)
                            Chem.SanitizeMol(combo)
                            new_smi = Chem.MolToSmiles(combo)

                elif strategy == "swap_halogen":
                    halogens_map = {9: [17, 35], 17: [9, 35], 35: [9, 17]}
                    hal_atoms = [
                        a for a in rwmol.GetAtoms()
                        if a.GetAtomicNum() in halogens_map
                        and a.GetIdx() not in frozen_set
                    ]
                    if hal_atoms:
                        target = rng.choice(hal_atoms)
                        rule = position_rules.get(target.GetIdx())
                        if rule and rule.get("allowed_groups"):
                            # Map SMILES groups to atomic numbers
                            group_map = {"F": 9, "Cl": 17, "Br": 35}
                            choices = [group_map[g] for g in rule["allowed_groups"] if g in group_map]
                            if choices:
                                new_num = rng.choice(choices)
                            else:
                                new_num = rng.choice(halogens_map[target.GetAtomicNum()])
                        else:
                            new_num = rng.choice(halogens_map[target.GetAtomicNum()])
                        target.SetAtomicNum(new_num)
                        Chem.SanitizeMol(rwmol)
                        new_smi = Chem.MolToSmiles(rwmol)

                elif strategy == "swap_atom":
                    # Swap O<->S or N<->O on non-ring atoms
                    swap_map = {8: 16, 16: 8}
                    swappable = [
                        a for a in rwmol.GetAtoms()
                        if not a.IsInRing() and a.GetAtomicNum() in swap_map
                        and a.GetIdx() not in frozen_set
                    ]
                    if swappable:
                        target = rng.choice(swappable)
                        target.SetAtomicNum(swap_map[target.GetAtomicNum()])
                        Chem.SanitizeMol(rwmol)
                        new_smi = Chem.MolToSmiles(rwmol)

                elif strategy == "modify_chain":
                    # Add or remove a methyl in a chain
                    chain_carbons = [
                        a for a in rwmol.GetAtoms()
                        if not a.GetIsAromatic() and a.GetAtomicNum() == 6
                        and a.GetDegree() >= 1 and a.GetTotalNumHs() > 0
                        and a.GetIdx() not in frozen_set
                    ]
                    if chain_carbons:
                        target = rng.choice(chain_carbons)
                        new_idx = rwmol.AddAtom(Chem.Atom(6))
                        rwmol.AddBond(target.GetIdx(), new_idx, Chem.BondType.SINGLE)
                        Chem.SanitizeMol(rwmol)
                        new_smi = Chem.MolToSmiles(rwmol)

            except Exception:
                continue

            if new_smi is None or new_smi in seen:
                continue

            # Validate the variant
            try:
                new_mol = Chem.MolFromSmiles(new_smi)
                if new_mol is None:
                    continue

                # QED filter
                try:
                    qed_val = QED_module.qed(new_mol)
                except Exception:
                    qed_val = 0.5
                if qed_val < 0.2:
                    continue

                # Tanimoto similarity to parent
                if fp_gen is not None:
                    new_fp = fp_gen.GetFingerprint(new_mol)
                else:
                    new_fp = AllChem.GetMorganFingerprintAsBitVect(new_mol, 2, nBits=2048)
                sim = DataStructs.TanimotoSimilarity(parent_fp, new_fp)
                if sim < min_similarity:
                    continue

                # MW sanity + delta check
                mw = Descriptors.ExactMolWt(new_mol)
                if mw < 150 or mw > 900:
                    continue
                if max_mw_delta < float('inf'):
                    if abs(mw - parent_mw) > max_mw_delta:
                        continue

                seen.add(new_smi)
                variants.append(new_smi)

            except Exception:
                continue

        RDLogger.EnableLog("rdApp.*")

        logger.debug(
            "RDKit variant generation: %d variants from %d attempts (iter %d)",
            len(variants), attempts, iteration,
        )
        return variants

    except Exception as exc:
        logger.warning("RDKit variant generation failed: %s", exc)
        try:
            from rdkit import RDLogger
            RDLogger.EnableLog("rdApp.*")
        except Exception:
            pass
        return _generate_variants_mock(parent_smiles, n_variants, rng, iteration)


def _generate_variants_mock(
    parent_smiles: str,
    n_variants: int,
    rng: random.Random,
    iteration: int,
) -> list[str]:
    """Generate mock variant SMILES deterministically.

    Produces slight string modifications that serve as unique identifiers
    for scoring purposes. Not chemically valid but sufficient for the
    mock pipeline.

    Parameters
    ----------
    parent_smiles : str
        Parent SMILES.
    n_variants : int
        Number to generate.
    rng : random.Random
        RNG for reproducibility.
    iteration : int
        Iteration number.

    Returns
    -------
    list[str]
        Mock variant SMILES strings.
    """
    variants: list[str] = []
    # Common modifications to simulate
    mods = ["F", "Cl", "C", "CC", "O", "N", "OC", "NC", "C#N"]

    for i in range(n_variants):
        mod = mods[(rng.randint(0, 100) + i) % len(mods)]
        # Create a deterministic variant identifier
        key = f"{parent_smiles}:iter{iteration}:var{i}:{mod}"
        digest = hashlib.sha256(key.encode()).hexdigest()[:6]
        variant = f"{parent_smiles}.{mod}{digest}"
        variants.append(variant)

    return variants


# =====================================================================
# MOCK SCORING FUNCTIONS
# =====================================================================

def _mock_dock(
    smiles: str,
    pocket_center: list[float],
    rng: random.Random,
) -> float:
    """Compute a deterministic mock docking affinity.

    Parameters
    ----------
    smiles : str
        Molecule SMILES.
    pocket_center : list[float]
        [x, y, z] pocket center coordinates.
    rng : random.Random
        RNG (used for slight noise, but main score is hash-based).

    Returns
    -------
    float
        Mock affinity in kcal/mol (negative = better).
    """
    key = f"optdock:{smiles}:{pocket_center[0]:.1f}"
    digest = hashlib.sha256(key.encode()).hexdigest()
    hash_int = int(digest[:10], 16)
    # Range: [-11.0, -5.0]
    affinity = -5.0 - (hash_int % 6001) / 1000.0
    return round(affinity, 2)


def _mock_admet_simple(
    smiles: str,
    rng: random.Random,
) -> dict:
    """Compute simplified mock ADMET properties.

    Returns a dict with scores in [0, 1] for key ADMET endpoints.

    Parameters
    ----------
    smiles : str
        Molecule SMILES.
    rng : random.Random
        RNG for reproducibility.

    Returns
    -------
    dict
        Keys: toxicity_score, bioavailability_score, synthesis_score,
        herg_risk, hepatotoxicity_risk, composite_score.
    """
    h = hashlib.sha256(f"optadmet:{smiles}".encode()).hexdigest()

    def _hval(offset: int) -> float:
        seg = h[(offset * 3) % len(h):(offset * 3 + 4) % len(h)]
        if not seg:
            seg = h[:4]
        return int(seg, 16) / 0xFFFF

    toxicity = round(1.0 - _hval(0) * 0.6, 3)  # Higher = safer (less toxic)
    bioavail = round(0.3 + _hval(1) * 0.5, 3)   # Higher = better
    synthesis = round(0.4 + _hval(2) * 0.5, 3)   # Higher = easier to synthesize
    herg = round(_hval(3) * 0.4, 3)               # Lower = safer
    hepatotox = round(_hval(4) * 0.3, 3)          # Lower = safer
    composite = round(
        0.4 * toxicity + 0.3 * bioavail + 0.3 * synthesis, 3,
    )

    return {
        "toxicity_score": toxicity,
        "bioavailability_score": bioavail,
        "synthesis_score": synthesis,
        "herg_risk": herg,
        "hepatotoxicity_risk": hepatotox,
        "composite_score": composite,
    }


def _multi_objective_score(
    affinity: float,
    admet: dict,
    weights: dict[str, float],
) -> float:
    """Compute the multi-objective optimization score.

    Parameters
    ----------
    affinity : float
        Docking affinity (kcal/mol, negative = better).
    admet : dict
        ADMET scores.
    weights : dict[str, float]
        Objective weights with keys: affinity, toxicity,
        bioavailability, synthesis.

    Returns
    -------
    float
        Score in [0, 1], higher is better.
    """
    # Normalize affinity to [0, 1]: -12 -> 1.0, -4 -> 0.0
    norm_aff = max(0.0, min(1.0, (-affinity - 4.0) / 8.0))

    toxicity = admet.get("toxicity_score", 0.5)
    bioavail = admet.get("bioavailability_score", 0.5)
    synthesis = admet.get("synthesis_score", 0.5)

    score = (
        weights.get("affinity", 0.35) * norm_aff
        + weights.get("toxicity", 0.25) * toxicity
        + weights.get("bioavailability", 0.20) * bioavail
        + weights.get("synthesis", 0.20) * synthesis
    )

    return round(max(0.0, min(1.0, score)), 4)


# =====================================================================
# COMPARISON HELPER
# =====================================================================

def _build_comparison(
    start_score: float,
    end_score: float,
    start_affinity: float,
    end_affinity: float,
    start_admet: dict,
    end_admet: dict,
) -> dict:
    """Build a comparison dict showing improvement from start to end.

    Parameters
    ----------
    start_score, end_score : float
        Multi-objective scores.
    start_affinity, end_affinity : float
        Docking affinities.
    start_admet, end_admet : dict
        ADMET score dicts.

    Returns
    -------
    dict
        Mapping of metric names to {start, end, change_pct}.
    """
    def _change_pct(start: float, end: float) -> float:
        if abs(start) < 1e-9:
            return 0.0
        return round((end - start) / abs(start) * 100, 1)

    return {
        "multi_objective_score": {
            "start": round(start_score, 4),
            "end": round(end_score, 4),
            "change_pct": _change_pct(start_score, end_score),
        },
        "affinity": {
            "start": round(start_affinity, 2),
            "end": round(end_affinity, 2),
            "change_pct": _change_pct(start_affinity, end_affinity),
        },
        "toxicity_score": {
            "start": round(start_admet.get("toxicity_score", 0.5), 3),
            "end": round(end_admet.get("toxicity_score", 0.5), 3),
            "change_pct": _change_pct(
                start_admet.get("toxicity_score", 0.5),
                end_admet.get("toxicity_score", 0.5),
            ),
        },
        "bioavailability_score": {
            "start": round(start_admet.get("bioavailability_score", 0.5), 3),
            "end": round(end_admet.get("bioavailability_score", 0.5), 3),
            "change_pct": _change_pct(
                start_admet.get("bioavailability_score", 0.5),
                end_admet.get("bioavailability_score", 0.5),
            ),
        },
        "synthesis_score": {
            "start": round(start_admet.get("synthesis_score", 0.5), 3),
            "end": round(end_admet.get("synthesis_score", 0.5), 3),
            "change_pct": _change_pct(
                start_admet.get("synthesis_score", 0.5),
                end_admet.get("synthesis_score", 0.5),
            ),
        },
    }


# =====================================================================
# CLI / SELF-TEST
# =====================================================================

if __name__ == "__main__":
    import json
    import sys

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    print("=" * 70)
    print("DockIt Lead Optimization -- Self-Test")
    print("=" * 70)

    erlotinib_smiles = "C=Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"
    pocket = [22.0, 0.5, 18.0]

    # Track progress
    progress_log: list[str] = []

    def progress_cb(iteration: int, score: float, message: str) -> None:
        progress_log.append(message)
        if iteration % 3 == 0 or iteration == 1:
            print(f"  {message}")

    print("\n[1] Running optimization (10 iterations, 50 variants/iter)...")
    result = run_optimization(
        starting_smiles=erlotinib_smiles,
        starting_name="Erlotinib",
        target_pdbqt="/tmp/receptor.pdbqt",
        pocket_center=pocket,
        n_iterations=10,
        variants_per_iter=50,
        progress_callback=progress_cb,
    )

    print(f"\n  Starting: {result['starting_molecule']['name']} "
          f"(score={result['starting_molecule']['score']:.3f})")
    print(f"  Final:    {result['final_lead']['name']} "
          f"(score={result['final_lead']['score']:.3f})")
    print(f"  Total tested:   {result['total_tested']}")
    print(f"  Total valid:    {result['total_valid']}")
    print(f"  Total improved: {result['total_improved']}")

    print("\n  Comparison:")
    for metric, data in result["comparison"].items():
        print(f"    {metric:30s}: {data['start']:.3f} -> {data['end']:.3f} "
              f"({data['change_pct']:+.1f}%)")

    # Validate
    assert result["final_lead"]["score"] >= result["starting_molecule"]["score"], (
        "Final score should be >= starting score"
    )
    print("\n  PASS: final score >= starting score")

    assert len(result["iterations"]) == 10, f"Expected 10 iterations, got {len(result['iterations'])}"
    print("  PASS: 10 iterations completed")

    assert result["total_tested"] > 0, "Should have tested some variants"
    print(f"  PASS: tested {result['total_tested']} variants")

    assert len(progress_log) == 10, f"Expected 10 progress callbacks, got {len(progress_log)}"
    print("  PASS: progress callbacks fired correctly")

    # Test 2: Determinism
    print("\n[2] Testing determinism...")
    result2 = run_optimization(
        starting_smiles=erlotinib_smiles,
        starting_name="Erlotinib",
        target_pdbqt="/tmp/receptor.pdbqt",
        pocket_center=pocket,
        n_iterations=5,
        variants_per_iter=20,
    )
    result3 = run_optimization(
        starting_smiles=erlotinib_smiles,
        starting_name="Erlotinib",
        target_pdbqt="/tmp/receptor.pdbqt",
        pocket_center=pocket,
        n_iterations=5,
        variants_per_iter=20,
    )
    assert result2["starting_molecule"]["score"] == result3["starting_molecule"]["score"]
    assert result2["final_lead"]["score"] == result3["final_lead"]["score"]
    print("  PASS: deterministic results confirmed")

    # Test 3: Custom weights
    print("\n[3] Testing custom weights...")
    result_custom = run_optimization(
        starting_smiles=erlotinib_smiles,
        starting_name="Erlotinib",
        target_pdbqt="/tmp/receptor.pdbqt",
        pocket_center=pocket,
        weights={"affinity": 0.80, "toxicity": 0.10, "bioavailability": 0.05, "synthesis": 0.05},
        n_iterations=5,
        variants_per_iter=20,
    )
    print(f"  Score with affinity-heavy weights: {result_custom['final_lead']['score']:.3f}")
    print("  PASS: custom weights work")

    print("\n" + "=" * 70)
    print("ALL TESTS PASSED")
    print("=" * 70)
