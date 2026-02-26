"""
DockIt pipeline -- Post-docking enrichment pipeline.

Extracts the multi-step enrichment logic from tasks.py into a single
reusable function.  Each step is wrapped in try/except so that partial
failures never abort the entire pipeline.

Steps (in order):
  1. ADMET predictions
  2. Composite scoring (v2 with ADMET, v1 without)
  3. Consensus detail enrichment (z-scores, agreement)
  4. Hard cutoffs (eliminatory filters)
  5. Butina clustering
  6. Pareto multi-objective ranking
  7. Combined off-target screening (top 5)
  8. Specialized hERG screening (top 20)
  9. Retrosynthesis planning (top 5, optional)
 10. Confidence scoring (all molecules)

The module also provides ``generate_3d_conformer`` for producing 3D SDF
files from SMILES, useful for optimization molecules that never went
through a docking engine.
"""

from __future__ import annotations

import logging
import tempfile
from collections import Counter
from pathlib import Path
from typing import Callable, Optional

logger = logging.getLogger(__name__)

# Type alias for the progress callback: (percent: int, message: str) -> None
ProgressCallback = Callable[[int, str], None]


# ---------------------------------------------------------------------------
# 3D conformer generation
# ---------------------------------------------------------------------------

def generate_3d_conformer(smiles: str, output_path: Path) -> Optional[Path]:
    """Generate a 3D SDF file from a SMILES string using RDKit.

    Uses the ETKDGv3 algorithm for initial embedding followed by MMFF
    force-field optimization.  The conformer is placed at the RDKit default
    coordinates (not translated to any pocket).  Only real docking poses
    are valid for 3D visualization in the binding site.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.
    output_path : Path
        Destination path for the SDF file.

    Returns
    -------
    Path or None
        The output path on success, ``None`` on failure.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning("RDKit could not parse SMILES: %.60s", smiles)
            return None

        mol = Chem.AddHs(mol)

        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        embed_status = AllChem.EmbedMolecule(mol, params)
        if embed_status != 0:
            # Fallback: try with random coordinates
            logger.debug("ETKDGv3 embedding failed, retrying with random coords")
            embed_status = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            if embed_status != 0:
                logger.warning("RDKit embedding failed for: %.60s", smiles)
                return None

        # Optimize geometry with MMFF; ignore return code (may not converge)
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        except Exception:
            logger.debug("MMFF optimization failed, keeping unoptimized conformer")

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        writer = Chem.SDWriter(str(output_path))
        writer.write(mol)
        writer.close()

        if output_path.exists() and output_path.stat().st_size > 0:
            logger.info("Generated 3D conformer: %s", output_path)
            return output_path

        logger.warning("SDF file is empty or missing: %s", output_path)
        return None

    except ImportError:
        logger.warning("RDKit not available -- cannot generate 3D conformer")
        return None
    except Exception as exc:
        logger.error("3D conformer generation failed for %.60s: %s", smiles, exc)
        return None


# ---------------------------------------------------------------------------
# Main enrichment pipeline
# ---------------------------------------------------------------------------

def enrich_results(
    molecules: list[dict],
    work_dir: Optional[Path] = None,
    structure_source: str = "cached",
    pocket_method: str = "unknown",
    disorder_info: Optional[dict] = None,
    enable_retrosynthesis: bool = False,
    progress_callback: Optional[ProgressCallback] = None,
    # V12: granular analysis flags
    enable_admet: bool = True,
    enable_selectivity: bool = True,
    enable_herg: bool = True,
    enable_safety: bool = True,
) -> tuple[list[dict], list[dict], dict]:
    """Run the full post-docking enrichment pipeline on molecule dicts.

    Each molecule dict must have at minimum ``smiles`` and ``name``.
    Optional fields such as ``affinity``, ``composite_score``, ``source``,
    ``vina_score``, ``cnn_score``, etc. will be used when present.

    The function never raises -- every step is wrapped in try/except.
    Partial failures are logged and recorded in the returned summary.

    Parameters
    ----------
    molecules : list[dict]
        Flat list of molecule result dicts (from docking, generation, etc.).
    work_dir : Path or None
        Working directory for intermediate files.  A temporary directory is
        created when ``None``.
    structure_source : str
        How the protein structure was obtained (e.g. ``"alphafold"``,
        ``"esmfold"``, ``"cached"``).  Passed to confidence scoring.
    pocket_method : str
        Pocket detection method (e.g. ``"fpocket"``).
    disorder_info : dict or None
        Disorder analysis result.  A ``fraction_disordered > 0.3`` applies
        a confidence penalty.
    enable_retrosynthesis : bool
        Whether to run retrosynthesis planning on the top candidates.
    progress_callback : callable or None
        Optional ``(percent: int, message: str) -> None`` callback invoked
        at the start of each pipeline step.

    Returns
    -------
    tuple[list[dict], list[dict], dict]
        ``(passed_molecules, eliminated_molecules, enrichment_summary)``

        * *passed_molecules* -- molecules surviving hard cutoffs, enriched
          with all downstream annotations, sorted by composite score.
        * *eliminated_molecules* -- molecules removed by hard cutoffs, each
          carrying an ``elimination_reason`` field.
        * *enrichment_summary* -- dict tracking which steps succeeded, counts,
          and sub-summaries (cutoff reasons, hERG risk, Pareto front size...).
    """

    # -- Housekeeping -------------------------------------------------------
    if work_dir is None:
        work_dir = Path(tempfile.mkdtemp(prefix="dockit_enrich_"))
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    summary: dict = {
        "steps_completed": [],
        "steps_failed": [],
        "n_input": len(molecules),
    }

    # Defensive copy -- we mutate dicts in place but never remove/add keys
    # from the caller's list unexpectedly.
    scored: list[dict] = list(molecules)
    eliminated: list[dict] = []

    def _progress(pct: int, msg: str) -> None:
        if progress_callback is not None:
            try:
                progress_callback(pct, msg)
            except Exception:
                pass

    # -----------------------------------------------------------------------
    # Step 1: ADMET predictions
    # -----------------------------------------------------------------------
    admet_results: list[dict] = []
    if enable_admet:
        _progress(0, "Predicting ADMET properties")
        try:
            from pipeline.admet import predict_admet

            unique_smiles = list({
                m.get("smiles", "") for m in scored if m.get("smiles")
            })
            if unique_smiles:
                admet_results = predict_admet(unique_smiles)
                logger.info("ADMET predictions for %d unique SMILES", len(admet_results))
            summary["n_admet"] = len(admet_results)
            summary["steps_completed"].append("admet")
        except Exception as exc:
            logger.warning("ADMET prediction failed: %s", exc)
            summary["steps_failed"].append("admet")
    else:
        logger.info("ADMET disabled -- skipping")

    _progress(15, f"ADMET complete for {len(admet_results)} molecules")

    # -----------------------------------------------------------------------
    # Step 2: Composite scoring
    # -----------------------------------------------------------------------
    _progress(20, "Computing composite scores")
    try:
        if admet_results:
            from pipeline.scoring import score_results_v2
            scored = score_results_v2(scored, admet_results)
        else:
            from pipeline.scoring import score_results
            scored = score_results(scored)

        summary["n_scored"] = len(scored)
        summary["steps_completed"].append("scoring")
        logger.info("Scoring complete: %d molecules", len(scored))
    except Exception as exc:
        logger.warning("Scoring failed: %s", exc)
        summary["steps_failed"].append("scoring")

    _progress(30, "Scoring complete")

    # -----------------------------------------------------------------------
    # Step 3: Consensus detail enrichment (z-scores, agreement)
    # -----------------------------------------------------------------------
    try:
        from pipeline.scoring import enrich_consensus_detail

        if scored:
            enrich_consensus_detail(scored)
        summary["steps_completed"].append("consensus_detail")
    except Exception as exc:
        logger.warning("Consensus detail enrichment failed: %s", exc)
        summary["steps_failed"].append("consensus_detail")

    _progress(35, "Consensus detail enrichment done")

    # -----------------------------------------------------------------------
    # Step 4: Hard cutoffs
    # -----------------------------------------------------------------------
    _progress(40, "Applying hard cutoffs")
    try:
        from pipeline.scoring import apply_hard_cutoffs

        scored, eliminated = apply_hard_cutoffs(scored)

        reason_counts = dict(Counter(
            mol.get("elimination_reason", "unknown") for mol in eliminated
        ))
        summary["hard_cutoffs"] = {
            "passed": len(scored),
            "eliminated": len(eliminated),
            "reasons": reason_counts,
        }
        summary["steps_completed"].append("hard_cutoffs")
        logger.info(
            "Hard cutoffs: %d passed, %d eliminated",
            len(scored), len(eliminated),
        )
    except Exception as exc:
        logger.warning("Hard cutoffs failed: %s", exc)
        summary["steps_failed"].append("hard_cutoffs")

    _progress(45, f"{len(scored)} molecules passed cutoffs")

    # -----------------------------------------------------------------------
    # Step 5: Butina clustering
    # -----------------------------------------------------------------------
    _progress(50, "Clustering molecules")
    try:
        from pipeline.scoring import cluster_results

        if scored:
            cluster_results(scored)
            n_clusters = len(set(m.get("cluster_id", 0) for m in scored))
        else:
            n_clusters = 0

        summary["n_clusters"] = n_clusters
        summary["steps_completed"].append("clustering")
        logger.info("Butina clustering: %d families", n_clusters)
    except Exception as exc:
        logger.warning("Butina clustering failed: %s", exc)
        summary["steps_failed"].append("clustering")

    _progress(55, "Clustering complete")

    # -----------------------------------------------------------------------
    # Step 6: Pareto multi-objective ranking
    # -----------------------------------------------------------------------
    _progress(58, "Pareto multi-objective ranking")
    try:
        from pipeline.scoring import pareto_ranking

        if scored:
            pareto_ranking(scored)
            front_size = sum(1 for m in scored if m.get("pareto_front", False))
        else:
            front_size = 0

        summary["pareto_front_size"] = front_size
        summary["steps_completed"].append("pareto_ranking")
        logger.info("Pareto ranking: %d on front out of %d", front_size, len(scored))
    except Exception as exc:
        logger.warning("Pareto ranking failed: %s", exc)
        summary["steps_failed"].append("pareto_ranking")

    _progress(62, "Pareto ranking complete")

    # -----------------------------------------------------------------------
    # Step 7: Combined off-target screening (top 5)
    # -----------------------------------------------------------------------
    if enable_selectivity:
        _progress(65, "Combined off-target screening (top 5)")
        try:
            from pipeline.off_target import combined_off_target_screening

            top5 = sorted(
                scored,
                key=lambda x: x.get("composite_score", 0),
                reverse=True,
            )[:5]

            combined_ot_summary: list[dict] = []
            for mol in top5:
                smi = mol.get("smiles", "")
                if not smi:
                    continue
                combined_result = combined_off_target_screening(smi, work_dir)
                mol["combined_off_target"] = combined_result
                combined_ot_summary.append({
                    "name": mol.get("name", "unknown"),
                    "combined_selectivity": combined_result.get("combined_selectivity", 0.0),
                    "tier1_hits": len(combined_result.get("tier1_hits", [])),
                    "tier2_safe": combined_result.get("tier2_safe_count", 0),
                })

            # Merge back into full scored list by SMILES
            combined_by_smiles = {
                m.get("smiles", ""): m.get("combined_off_target")
                for m in top5
                if m.get("combined_off_target")
            }
            for mol in scored:
                smi = mol.get("smiles", "")
                if smi in combined_by_smiles:
                    mol["combined_off_target"] = combined_by_smiles[smi]

            summary["combined_off_target_summary"] = combined_ot_summary
            summary["steps_completed"].append("combined_off_target")
            logger.info("Combined off-target: screened %d molecules", len(top5))
        except Exception as exc:
            logger.warning("Combined off-target screening failed: %s", exc)
            summary["steps_failed"].append("combined_off_target")
    else:
        logger.info("Off-target selectivity disabled -- skipping")

    _progress(72, "Off-target screening complete")

    # -----------------------------------------------------------------------
    # Step 8: Specialized hERG screening (top 20)
    # -----------------------------------------------------------------------
    if enable_herg:
        _progress(75, "Specialized hERG screening (top 20)")
        try:
            from pipeline.admet import predict_herg_specialized

            top20 = sorted(
                scored,
                key=lambda x: x.get("composite_score", 0),
                reverse=True,
            )[:20]

            herg_summary_counts: dict[str, int] = {"LOW": 0, "MODERATE": 0, "HIGH": 0}
            for mol in top20:
                smi = mol.get("smiles", "")
                if not smi:
                    continue
                herg_result = predict_herg_specialized(smi)
                mol["herg_specialized"] = herg_result
                herg_summary_counts[herg_result.get("risk_level", "LOW")] += 1

            # Merge back by SMILES
            herg_by_smiles = {
                m.get("smiles", ""): m.get("herg_specialized")
                for m in top20
                if m.get("herg_specialized")
            }
            for mol in scored:
                smi = mol.get("smiles", "")
                if smi in herg_by_smiles:
                    mol["herg_specialized"] = herg_by_smiles[smi]

            summary["herg_specialized_summary"] = herg_summary_counts
            summary["steps_completed"].append("herg_specialized")
            logger.info(
                "Specialized hERG: %d screened (LOW=%d, MOD=%d, HIGH=%d)",
                len(top20),
                herg_summary_counts["LOW"],
                herg_summary_counts["MODERATE"],
                herg_summary_counts["HIGH"],
            )
        except Exception as exc:
            logger.warning("Specialized hERG screening failed: %s", exc)
            summary["steps_failed"].append("herg_specialized")
    else:
        logger.info("hERG screening disabled -- skipping")

    _progress(80, "hERG screening complete")

    # -----------------------------------------------------------------------
    # Step 9: Retrosynthesis (top 5, optional)
    # -----------------------------------------------------------------------
    if enable_retrosynthesis:
        _progress(82, "Planning retrosynthesis for top molecules")
        try:
            from pipeline.retrosynthesis import plan_synthesis

            top_for_synthesis = sorted(
                scored,
                key=lambda x: x.get("composite_score", 0),
                reverse=True,
            )[:5]

            for i, mol in enumerate(top_for_synthesis):
                smiles = mol.get("smiles", "")
                name = mol.get("name", f"mol_{i}")
                if not smiles:
                    continue
                _progress(82 + i * 2, f"Retrosynthesis: {name}")
                try:
                    route = plan_synthesis(smiles, max_depth=6, timeout_sec=120)
                    mol["synthesis_route"] = route
                except Exception as synth_exc:
                    logger.warning("Retrosynthesis failed for %s: %s", name, synth_exc)
                    mol["synthesis_route"] = None

            # Merge back by SMILES
            synth_by_smiles = {
                m.get("smiles", ""): m.get("synthesis_route")
                for m in top_for_synthesis
                if m.get("synthesis_route") is not None
            }
            for mol in scored:
                smi = mol.get("smiles", "")
                if smi in synth_by_smiles and "synthesis_route" not in mol:
                    mol["synthesis_route"] = synth_by_smiles[smi]

            summary["n_retrosynthesis"] = len([
                m for m in top_for_synthesis if m.get("synthesis_route")
            ])
            summary["steps_completed"].append("retrosynthesis")
        except Exception as exc:
            logger.warning("Retrosynthesis failed: %s", exc)
            summary["steps_failed"].append("retrosynthesis")
    else:
        logger.info("Retrosynthesis disabled -- skipping")

    _progress(92, "Retrosynthesis complete" if enable_retrosynthesis else "Skipping retrosynthesis")

    # -----------------------------------------------------------------------
    # Step 10: Confidence scoring (all molecules)
    # -----------------------------------------------------------------------
    if enable_safety:
        _progress(94, "Computing confidence scores")
        try:
            from pipeline.confidence import calculate_confidence

            for mol in scored:
                mol.setdefault("pocket_method", pocket_method)
                mol["confidence"] = calculate_confidence(
                    mol, structure_source, disorder_info=disorder_info,
                )

            # Also score eliminated molecules (useful for reporting)
            for mol in eliminated:
                mol.setdefault("pocket_method", pocket_method)
                mol["confidence"] = calculate_confidence(
                    mol, structure_source, disorder_info=disorder_info,
                )

            n_confidence = len(scored) + len(eliminated)
            summary["steps_completed"].append("confidence")
            summary["disorder_penalty_applied"] = bool(
                disorder_info and disorder_info.get("fraction_disordered", 0) > 0.3
            )
            logger.info("Confidence scores computed for %d molecules", n_confidence)
        except Exception as exc:
            logger.warning("Confidence scoring failed: %s", exc)
            summary["steps_failed"].append("confidence")
    else:
        logger.info("Safety/confidence scoring disabled -- skipping")

    _progress(100, "Enrichment complete")

    # -----------------------------------------------------------------------
    # Final sort: scored molecules by composite_score descending
    # -----------------------------------------------------------------------
    scored.sort(key=lambda x: x.get("composite_score", 0), reverse=True)

    summary["n_passed"] = len(scored)
    summary["n_eliminated"] = len(eliminated)

    logger.info(
        "Enrichment pipeline complete: %d passed, %d eliminated, "
        "steps OK=%s, steps FAILED=%s",
        len(scored),
        len(eliminated),
        summary["steps_completed"],
        summary["steps_failed"],
    )

    return scored, eliminated, summary
