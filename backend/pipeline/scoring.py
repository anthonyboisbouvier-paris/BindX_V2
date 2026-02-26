"""
DockIt pipeline — Molecular property computation and composite scoring.

Uses RDKit to compute drug-likeness descriptors and generate 2D structure
images (SVG). Falls back gracefully when RDKit is not available.

V2: adds compute_composite_score_v2 integrating ADMET + novelty.
V5bis: adds apply_hard_cutoffs, cluster_results, PAINS detection, SA score.
V6.1: adds enrich_consensus_detail with z-scores, agreement, per-method ranks.
"""

from __future__ import annotations

import logging
import math
from typing import Optional

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# NaN-safe helpers
# ---------------------------------------------------------------------------

def _safe_float(val) -> float:
    """Convert to float, treating NaN and None as 0.0."""
    try:
        v = float(val if val is not None else 0.0)
        return v if not math.isnan(v) else 0.0
    except (TypeError, ValueError):
        return 0.0


def _safe_zscore(val: float, mu: float, sd: float) -> float:
    """Compute z-score safely, returning 0.0 on any edge case."""
    if sd <= 0:
        return 0.0
    try:
        z = (val - mu) / sd
        return 0.0 if (math.isnan(z) or math.isinf(z)) else round(z, 3)
    except (TypeError, ValueError, ZeroDivisionError):
        return 0.0


# ---------------------------------------------------------------------------
# PAINS SMARTS patterns (Pan-Assay INterference compoundS)
# Simplified set of common PAINS substructure alerts
# ---------------------------------------------------------------------------

_PAINS_SMARTS: list[str] = [
    "[#6]1:[#6]:[#6](:[#6]:[#6]:[#6]:1)-[#7]=[#7]-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2",  # azo dye
    "[#6]-[#16](=[#8])(-[#6])-[#7]",                                                     # sulfonamide
    "c1cc([OH])c([OH])cc1",                                                                  # catechol (free OH only)
    "[#6](=[#8])-[#6]=[#6]-[#6](=[#8])",                                                 # quinone-like
    "[CH]=[CH][C](=O)",                                                                    # Michael acceptor (enone)
    "c1ccc2c(c1)[nH]c1ccccc12",                                                           # carbazole
    "[#7+]=[#6]-[#16]",                                                                    # isothiourea
    "c1ccnc(SSc2ncccc2)c1",                                                                # 2-amino pyridine disulfide
    "[N;R1]=[N;R1]",                                                                       # azo in ring
    "c1ccc(-c2cc(=O)c3ccccc3o2)cc1",                                                       # flavone (promiscuous)
]


# ---------------------------------------------------------------------------
# Molecular properties
# ---------------------------------------------------------------------------

def compute_properties(smiles: str) -> dict:
    """Compute physicochemical and drug-likeness properties for a molecule.

    Parameters
    ----------
    smiles : str
        Canonical SMILES string.

    Returns
    -------
    dict
        Keys: ``logP``, ``MW``, ``tpsa``, ``qed``, ``hbd``, ``hba``,
        ``rotatable_bonds``. Values may be ``None`` if computation fails.
    """
    props: dict = {
        "logP": None,
        "MW": None,
        "tpsa": None,
        "qed": None,
        "hbd": None,
        "hba": None,
        "rotatable_bonds": None,
    }

    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, QED, rdMolDescriptors

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning("RDKit cannot parse SMILES for properties: %s", smiles[:80])
            return props

        props["logP"] = round(Descriptors.MolLogP(mol), 2)
        props["MW"] = round(Descriptors.ExactMolWt(mol), 2)
        props["tpsa"] = round(Descriptors.TPSA(mol), 2)
        props["hbd"] = rdMolDescriptors.CalcNumHBD(mol)
        props["hba"] = rdMolDescriptors.CalcNumHBA(mol)
        props["rotatable_bonds"] = rdMolDescriptors.CalcNumRotatableBonds(mol)

        try:
            props["qed"] = round(QED.qed(mol), 3)
        except Exception:
            props["qed"] = 0.5  # Safe default

    except ImportError:
        logger.warning("RDKit not available; returning empty properties")
    except Exception as exc:
        logger.warning("Property computation error for %s: %s", smiles[:60], exc)

    return props


# ---------------------------------------------------------------------------
# PAINS detection
# ---------------------------------------------------------------------------

def compute_pains_alert(smiles: str) -> bool:
    """Check if a molecule matches any PAINS substructure pattern.

    Parameters
    ----------
    smiles : str
        Canonical SMILES string.

    Returns
    -------
    bool
        True if one or more PAINS patterns match, False otherwise.
    """
    try:
        from rdkit import Chem

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        for smarts in _PAINS_SMARTS:
            try:
                pattern = Chem.MolFromSmarts(smarts)
                if pattern is not None and mol.HasSubstructMatch(pattern):
                    return True
            except Exception:
                continue
        return False

    except ImportError:
        return False
    except Exception as exc:
        logger.debug("PAINS check failed for %s: %s", smiles[:60], exc)
        return False


# ---------------------------------------------------------------------------
# Synthetic Accessibility (SA) score
# ---------------------------------------------------------------------------

def compute_sa_score(smiles: str) -> Optional[float]:
    """Compute the Synthetic Accessibility score for a molecule.

    Uses RDKit's SA_Score contribution model (Ertl & Schuffenhauer, 2009).
    Score range is 1 (easy to synthesize) to 10 (difficult).

    Parameters
    ----------
    smiles : str
        Canonical SMILES string.

    Returns
    -------
    float or None
        SA score in [1, 10], or None if computation fails.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # Try to use the dedicated SA_Score module if available
        try:
            from rdkit.Chem import RDConfig
            import os
            import sys
            sa_module_path = os.path.join(RDConfig.RDContribDir, "SA_Score")
            if sa_module_path not in sys.path:
                sys.path.insert(0, sa_module_path)
            import sascorer  # type: ignore[import-untyped]
            return round(sascorer.calculateScore(mol), 2)
        except (ImportError, Exception):
            pass

        # Fallback: heuristic SA score based on molecular complexity
        n_atoms = mol.GetNumHeavyAtoms()
        n_rings = rdMolDescriptors.CalcNumRings(mol)
        n_stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
        n_arom_rings = rdMolDescriptors.CalcNumAromaticRings(mol)

        # Heuristic: complexity-based estimate
        # Baseline score of 3 (moderate complexity)
        score = 3.0

        # Penalize large molecules
        if n_atoms > 35:
            score += (n_atoms - 35) * 0.05
        elif n_atoms < 15:
            score -= 0.5  # Small molecules are easier

        # Penalize many rings (macrocycles, complex scaffolds)
        if n_rings > 4:
            score += (n_rings - 4) * 0.4

        # Penalize stereocenters
        score += n_stereo * 0.3

        # Penalize non-aromatic rings (harder to build)
        non_arom_rings = n_rings - n_arom_rings
        if non_arom_rings > 2:
            score += (non_arom_rings - 2) * 0.3

        # Reward common drug-like features
        if 1 <= n_arom_rings <= 3:
            score -= 0.3

        return round(max(1.0, min(10.0, score)), 2)

    except ImportError:
        return None
    except Exception as exc:
        logger.debug("SA score computation failed for %s: %s", smiles[:60], exc)
        return None


# ---------------------------------------------------------------------------
# Composite score
# ---------------------------------------------------------------------------

def compute_composite_score(
    affinity: float,
    qed: Optional[float] = None,
    logp: Optional[float] = None,
) -> float:
    """Compute a weighted composite docking score.

    Higher is better.  Docking affinity is the dominant term so that
    known target-specific binders rank above generic druglike molecules.

    Formula
    -------
    ``0.65 * norm_affinity + 0.20 * qed + 0.15 * logp_penalty``

    - ``norm_affinity``: ``min(-14, affinity) / -14`` clamped to [0, 1].
    - ``logp_penalty``: Gaussian-like centered at 2.5 (sigma 2.5).

    Parameters
    ----------
    affinity : float
        Binding affinity in kcal/mol (negative = better).
    qed : float or None
        Quantitative Estimate of Drug-likeness [0, 1].
    logp : float or None
        Octanol-water partition coefficient.

    Returns
    -------
    float
        Composite score in approximately [0, 1].
    """
    import math

    # Normalise affinity: cap at -14 kcal/mol, map to [0, 1]
    capped = min(-14.0, affinity) if affinity < 0 else affinity
    norm_affinity = max(0.0, min(1.0, capped / -14.0))

    # QED (default 0.5 if unknown)
    qed_val = qed if qed is not None else 0.5

    # Continuous logP penalty: Gaussian centered at 2.5, sigma=2.5
    if logp is not None:
        logp_penalty = math.exp(-0.5 * ((logp - 2.5) / 2.5) ** 2)
    else:
        logp_penalty = 0.5

    score = 0.65 * norm_affinity + 0.20 * qed_val + 0.15 * logp_penalty
    return round(score, 4)


# ---------------------------------------------------------------------------
# V2 Composite score (Vina + ADMET + drug-likeness + novelty)
# ---------------------------------------------------------------------------

def compute_composite_score_v2(
    affinity: float,
    admet_score: Optional[float] = None,
    qed: Optional[float] = None,
    novelty: Optional[float] = None,
) -> float:
    """Compute the V2 weighted composite score.

    Formula: vina * 0.55 + admet * 0.20 + drug_likeness * 0.15 + novelty * 0.10

    Docking affinity is the dominant term.

    Parameters
    ----------
    affinity : float
        Binding affinity in kcal/mol (negative = better).
    admet_score : float or None
        ADMET composite score [0, 1].
    qed : float or None
        Drug-likeness QED [0, 1].
    novelty : float or None
        Novelty score [0, 1] (1 = very novel).

    Returns
    -------
    float
        Composite V2 score in [0, 1].
    """
    capped = min(-14.0, affinity) if affinity < 0 else affinity
    norm_affinity = max(0.0, min(1.0, capped / -14.0))

    admet_val = admet_score if admet_score is not None else 0.5
    qed_val = qed if qed is not None else 0.5
    novelty_val = novelty if novelty is not None else 0.0

    score = 0.55 * norm_affinity + 0.20 * admet_val + 0.15 * qed_val + 0.10 * novelty_val
    return round(score, 4)


# ---------------------------------------------------------------------------
# 2D SVG generation
# ---------------------------------------------------------------------------

def generate_2d_svg(smiles: str, width: int = 300, height: int = 200) -> Optional[str]:
    """Render a 2D structure image as an SVG string.

    Parameters
    ----------
    smiles : str
        Canonical SMILES.
    width : int
        Image width in pixels.
    height : int
        Image height in pixels.

    Returns
    -------
    str or None
        SVG markup string, or ``None`` if rendering fails.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return svg

    except ImportError:
        logger.debug("RDKit not available for SVG generation")
        return None
    except Exception as exc:
        logger.warning("SVG generation error for %s: %s", smiles[:60], exc)
        return None


# ---------------------------------------------------------------------------
# Score a full results list
# ---------------------------------------------------------------------------

def score_results(docking_results: list[dict]) -> list[dict]:
    """Enrich docking results with properties and composite scores.

    Modifies the dicts in-place and returns the list sorted by
    composite_score descending (best first).

    Parameters
    ----------
    docking_results : list[dict]
        Output from ``dock_all_ligands``.

    Returns
    -------
    list[dict]
        Same list, enriched and re-sorted.
    """
    for entry in docking_results:
        smiles = entry.get("smiles", "")
        affinity = entry.get("affinity", 0.0)

        # Compute properties
        props = compute_properties(smiles)
        entry.update(props)

        # V5bis: PAINS alert and SA score
        entry["pains_alert"] = compute_pains_alert(smiles)
        entry["sa_score"] = compute_sa_score(smiles)

        # Composite score
        entry["composite_score"] = compute_composite_score(
            affinity,
            qed=props.get("qed"),
            logp=props.get("logP"),
        )

        # 2D SVG
        entry["svg_2d"] = generate_2d_svg(smiles)

    # Sort by composite score descending
    docking_results.sort(key=lambda r: r.get("composite_score", 0.0), reverse=True)
    return docking_results


def score_results_v2(
    docking_results: list[dict],
    admet_results: Optional[list[dict]] = None,
) -> list[dict]:
    """V2 scoring: enrich results with properties, ADMET, and V2 composite.

    Parameters
    ----------
    docking_results : list[dict]
        Output from docking (Vina or DiffDock).
    admet_results : list[dict], optional
        ADMET predictions keyed by SMILES.

    Returns
    -------
    list[dict]
        Enriched and sorted results.
    """
    # Build ADMET lookup by SMILES
    admet_by_smiles: dict = {}
    if admet_results:
        for a in admet_results:
            smi = a.get("smiles", "")
            if smi:
                admet_by_smiles[smi] = a

    for entry in docking_results:
        smiles = entry.get("smiles", "")
        affinity = entry.get("affinity", 0.0)

        # Compute basic properties
        props = compute_properties(smiles)
        entry.update(props)

        # V5bis: PAINS alert and SA score
        entry["pains_alert"] = compute_pains_alert(smiles)
        entry["sa_score"] = compute_sa_score(smiles)

        # 2D SVG
        entry["svg_2d"] = generate_2d_svg(smiles)

        # ADMET data
        admet_data = admet_by_smiles.get(smiles)
        if admet_data:
            entry["admet"] = admet_data
            admet_composite = admet_data.get("composite_score", 0.5)
        else:
            admet_composite = None

        novelty = entry.get("novelty_score")

        # V2 composite score
        entry["composite_score"] = compute_composite_score_v2(
            affinity,
            admet_score=admet_composite,
            qed=props.get("qed"),
            novelty=novelty,
        )

    docking_results.sort(key=lambda r: r.get("composite_score", 0.0), reverse=True)
    return docking_results


# ---------------------------------------------------------------------------
# V5bis: Hard cutoffs (eliminatory, non-compensable)
# ---------------------------------------------------------------------------

def apply_hard_cutoffs(molecules: list[dict]) -> tuple[list[dict], list[dict]]:
    """Apply eliminatory hard cutoffs BEFORE composite scoring.

    These cutoffs are non-compensable -- a molecule that triggers any one
    of them is removed from the ranking regardless of its other properties.

    Parameters
    ----------
    molecules : list[dict]
        Scored molecules (must already have properties computed).

    Returns
    -------
    tuple[list[dict], list[dict]]
        ``(passed, eliminated)`` where each eliminated molecule includes
        ``eliminated: True`` and ``elimination_reason: str``.

    Hard cutoffs
    ------------
    - hERG pIC50 predicted > 6.0 (cardiac risk)
    - Lipinski violations > 1 (oral bioavailability)
    - QED < 0.25 (drug-likeness)
    - SA score > 6.0 (not synthesizable)
    - PAINS alert (likely false positive)
    - CNN score < 0.2 (incorrect pose, if available)
    """
    passed: list[dict] = []
    eliminated: list[dict] = []

    for mol in molecules:
        reason: Optional[str] = None

        # hERG check from ADMET data
        admet = mol.get("admet", {})
        if isinstance(admet, dict):
            toxicity = admet.get("toxicity")
            if isinstance(toxicity, dict):
                herg = toxicity.get("herg_inhibition")
            else:
                herg = admet.get("herg_inhibition")
            if herg is not None and herg > 0.7:  # pIC50 > 6 proxy (0.7 threshold for heuristic ADMET)
                reason = "hERG inhibition risk (cardiac)"

        # Lipinski violations
        if reason is None:
            violations = 0
            mw = mol.get("mw") or mol.get("MW") or 0
            logp = mol.get("logp") or mol.get("logP") or 0
            hbd = mol.get("hbd") or 0
            hba = mol.get("hba") or 0
            if mw > 500:
                violations += 1
            if logp > 5:
                violations += 1
            if hbd > 5:
                violations += 1
            if hba > 10:
                violations += 1
            if violations > 2:
                reason = f"Lipinski violations ({violations})"

        # QED
        if reason is None:
            qed = mol.get("qed")
            if qed is not None and qed < 0.25:
                reason = f"QED too low ({qed:.2f})"

        # SA score
        if reason is None:
            sa = mol.get("sa_score")
            if sa is not None and sa > 6.0:
                reason = f"SA score too high ({sa:.1f})"

        # PAINS
        if reason is None:
            if mol.get("pains_alert"):
                reason = "PAINS alert (likely false positive)"

        # CNN score (GNINA) — only apply when CNN score is reliable.
        # GNINA with PDBQT input may produce 0.0 CNN scores (clamped from
        # negative values), which indicates an input format issue, not a bad pose.
        # Only eliminate if CNN was genuinely computed and is low.
        if reason is None:
            cnn = mol.get("cnn_score")
            engine = mol.get("docking_engine", "")
            if cnn is not None and engine == "gnina" and cnn > 0.0 and cnn < 0.2:
                reason = f"CNN pose score too low ({cnn:.2f})"

        if reason:
            mol["eliminated"] = True
            mol["elimination_reason"] = reason
            eliminated.append(mol)
        else:
            mol["eliminated"] = False
            passed.append(mol)

    logger.info(
        "Hard cutoffs: %d/%d passed. %d eliminated",
        len(passed), len(molecules), len(eliminated),
    )
    return passed, eliminated


# ---------------------------------------------------------------------------
# V6.1: Consensus detail enrichment (z-scores, agreement, per-method ranks)
# ---------------------------------------------------------------------------

def enrich_consensus_detail(molecules: list[dict]) -> list[dict]:
    """Add detailed consensus scoring metadata to each molecule.

    For the batch of molecules, computes z-scores of vina_score, cnn_score,
    and cnn_affinity, determines how many scoring methods agree that a
    molecule is in the top 30%, and records per-method ranks.

    Modifies molecules in-place and returns the list.

    Parameters
    ----------
    molecules : list[dict]
        Scored molecules. Expected keys: ``vina_score``, ``cnn_score``,
        ``cnn_affinity``. Missing values default to 0.0.

    Returns
    -------
    list[dict]
        Same list with ``consensus_detail`` dict added to each molecule::

            {
                "z_vina": float,
                "z_cnn_score": float,
                "z_cnn_affinity": float,
                "agreement": "3/3" | "2/3" | "1/3" | "0/3",
                "per_method_ranks": {
                    "vina": int,
                    "cnn_score": int,
                    "cnn_affinity": int,
                }
            }
    """
    n = len(molecules)
    if n == 0:
        return molecules

    # -- Collect raw scores ------------------------------------------------
    vina_vals: list[float] = []
    cnn_score_vals: list[float] = []
    cnn_aff_vals: list[float] = []

    for mol in molecules:
        vina_vals.append(_safe_float(mol.get("vina_score")))
        cnn_score_vals.append(_safe_float(mol.get("cnn_score")))
        cnn_aff_vals.append(_safe_float(mol.get("cnn_affinity")))

    # -- Compute mean and std for z-scores ---------------------------------
    def _mean_std(vals: list[float]) -> tuple[float, float]:
        if not vals:
            return 0.0, 1.0
        mu = sum(vals) / len(vals)
        variance = sum((v - mu) ** 2 for v in vals) / len(vals)
        sigma = math.sqrt(variance) if variance > 0 else 1.0
        return mu, sigma

    vina_mu, vina_sd = _mean_std(vina_vals)
    cnn_mu, cnn_sd = _mean_std(cnn_score_vals)
    aff_mu, aff_sd = _mean_std(cnn_aff_vals)

    # -- Compute per-method ranks ------------------------------------------
    # vina_score: lower is better -> sort ascending
    vina_order = sorted(range(n), key=lambda i: vina_vals[i])
    # cnn_score: higher is better -> sort descending
    cnn_order = sorted(range(n), key=lambda i: -cnn_score_vals[i])
    # cnn_affinity: higher is better -> sort descending
    aff_order = sorted(range(n), key=lambda i: -cnn_aff_vals[i])

    vina_rank: dict[int, int] = {idx: rank + 1 for rank, idx in enumerate(vina_order)}
    cnn_rank: dict[int, int] = {idx: rank + 1 for rank, idx in enumerate(cnn_order)}
    aff_rank: dict[int, int] = {idx: rank + 1 for rank, idx in enumerate(aff_order)}

    # -- Top 30% threshold -------------------------------------------------
    top_30_cutoff = max(1, int(math.ceil(n * 0.3)))

    # -- Enrich each molecule ----------------------------------------------
    for i, mol in enumerate(molecules):
        # Z-scores: for vina, more negative is better so we negate
        # (z_vina > 0 means better than average)
        z_vina = -_safe_zscore(vina_vals[i], vina_mu, vina_sd)
        z_cnn_score = _safe_zscore(cnn_score_vals[i], cnn_mu, cnn_sd)
        z_cnn_aff = _safe_zscore(cnn_aff_vals[i], aff_mu, aff_sd)

        # Agreement: how many methods rank this in top 30%
        r_vina = vina_rank[i]
        r_cnn = cnn_rank[i]
        r_aff = aff_rank[i]
        in_top = sum(1 for r in [r_vina, r_cnn, r_aff] if r <= top_30_cutoff)

        mol["consensus_detail"] = {
            "z_vina": z_vina,
            "z_cnn_score": z_cnn_score,
            "z_cnn_affinity": z_cnn_aff,
            "agreement": f"{in_top}/3",
            "per_method_ranks": {
                "vina": r_vina,
                "cnn_score": r_cnn,
                "cnn_affinity": r_aff,
            },
        }

    logger.info(
        "Consensus detail enrichment: %d molecules, top-30%% cutoff rank=%d",
        n, top_30_cutoff,
    )
    return molecules


# ---------------------------------------------------------------------------
# V5bis: Butina clustering on Morgan fingerprints
# ---------------------------------------------------------------------------

def cluster_results(molecules: list[dict], cutoff: float = 0.4) -> list[dict]:
    """Cluster results using the Butina algorithm on Morgan fingerprints.

    Assigns ``cluster_id`` and ``is_representative`` to each molecule.
    The representative of each cluster is the molecule with the highest
    ``composite_score``.

    Parameters
    ----------
    molecules : list[dict]
        Scored molecules with ``smiles`` and ``composite_score``.
    cutoff : float
        Tanimoto distance cutoff for clustering (default 0.4).

    Returns
    -------
    list[dict]
        Same list with added ``cluster_id``, ``cluster_size``, and
        ``is_representative`` fields.
    """
    if not molecules:
        return molecules

    try:
        from rdkit import Chem, DataStructs
        from rdkit.Chem import AllChem
        from rdkit.ML.Cluster import Butina

        # Generate fingerprints
        fps: list = []
        valid_indices: list[int] = []
        for i, mol in enumerate(molecules):
            smiles = mol.get("smiles", "")
            rdmol = Chem.MolFromSmiles(smiles) if smiles else None
            if rdmol:
                fp = AllChem.GetMorganFingerprintAsBitVect(rdmol, 2, nBits=2048)
                fps.append(fp)
                valid_indices.append(i)
            else:
                mol["cluster_id"] = -1
                mol["cluster_size"] = 1
                mol["is_representative"] = False

        if len(fps) < 2:
            for mol in molecules:
                mol["cluster_id"] = 0
                mol["cluster_size"] = 1
                mol["is_representative"] = True
            return molecules

        # Compute distance matrix (lower triangle)
        n = len(fps)
        dists: list[float] = []
        for i in range(1, n):
            for j in range(i):
                dist = 1.0 - DataStructs.TanimotoSimilarity(fps[i], fps[j])
                dists.append(dist)

        # Butina clustering
        clusters = Butina.ClusterData(dists, n, cutoff, isDistData=True)

        # Assign cluster IDs
        for cluster_idx, members in enumerate(clusters):
            # Find representative (highest composite_score)
            best_member = max(
                members,
                key=lambda m: molecules[valid_indices[m]].get("composite_score", 0),
            )
            for member in members:
                mol_idx = valid_indices[member]
                molecules[mol_idx]["cluster_id"] = cluster_idx
                molecules[mol_idx]["cluster_size"] = len(members)
                molecules[mol_idx]["is_representative"] = (member == best_member)

        n_clusters = len(clusters)
        logger.info(
            "Butina clustering: %d molecules -> %d chemical families (cutoff=%.2f)",
            len(molecules), n_clusters, cutoff,
        )
        return molecules

    except ImportError:
        logger.warning("RDKit not available for clustering, skipping")
        for mol in molecules:
            mol["cluster_id"] = 0
            mol["cluster_size"] = 1
            mol["is_representative"] = True
        return molecules
    except Exception as exc:
        logger.warning("Clustering failed: %s, assigning all to cluster 0", exc)
        for mol in molecules:
            mol["cluster_id"] = 0
            mol["cluster_size"] = 1
            mol["is_representative"] = True
        return molecules


# ---------------------------------------------------------------------------
# V6.2: Pareto multi-objective ranking
# ---------------------------------------------------------------------------

def _extract_pareto_objectives(mol: dict, n_molecules: int) -> dict[str, float]:
    """Extract and normalize four objective scores for Pareto ranking.

    All objectives are normalized to [0, 1] where higher is better.

    Parameters
    ----------
    mol : dict
        Scored molecule dictionary.
    n_molecules : int
        Total number of molecules (used to normalize consensus_rank).

    Returns
    -------
    dict[str, float]
        Four objective scores: affinity, safety, bioavailability, synthesis.
    """
    # --- Affinity ---
    # Prefer consensus_rank if available (lower rank = better).
    # Normalize: best rank (1) -> 1.0, worst rank (n) -> 0.0.
    consensus_rank = mol.get("consensus_rank")
    if consensus_rank is not None and n_molecules > 1:
        affinity = 1.0 - (float(consensus_rank) - 1.0) / (n_molecules - 1.0)
        affinity = max(0.0, min(1.0, affinity))
    else:
        # Fallback to vina_score: more negative is better.
        # Min-max normalize assuming range [-12, 0].
        vina = _safe_float(mol.get("vina_score")) or _safe_float(mol.get("affinity"))
        affinity = max(0.0, min(1.0, vina / -12.0)) if vina < 0 else 0.0

    # --- Safety ---
    # Use ADMET composite_score if available (already 0-1).
    # Combine with off_target selectivity_score / 10 if available.
    safety_components: list[float] = []
    admet = mol.get("admet")
    if isinstance(admet, dict):
        admet_composite = admet.get("composite_score")
        if admet_composite is not None:
            safety_components.append(float(admet_composite))
    off_target = mol.get("off_target_results")
    if isinstance(off_target, dict):
        selectivity = off_target.get("selectivity_score")
        if selectivity is not None:
            safety_components.append(float(selectivity) / 10.0)
    if safety_components:
        safety = sum(safety_components) / len(safety_components)
        safety = max(0.0, min(1.0, safety))
    else:
        safety = 0.5

    # --- Bioavailability ---
    # Use QED value (already 0-1).
    qed = mol.get("qed")
    if qed is not None:
        bioavailability = max(0.0, min(1.0, float(qed)))
    else:
        bioavailability = 0.5

    # --- Synthesis ---
    # Use 1 - (sa_score / 10). SA score range is [1, 10]: 1 = easy, 10 = hard.
    sa_score = mol.get("sa_score")
    if sa_score is not None:
        synthesis = max(0.0, min(1.0, 1.0 - float(sa_score) / 10.0))
    else:
        synthesis = 0.5

    return {
        "affinity": round(affinity, 4),
        "safety": round(safety, 4),
        "bioavailability": round(bioavailability, 4),
        "synthesis": round(synthesis, 4),
    }


def _dominates(a: dict[str, float], b: dict[str, float]) -> bool:
    """Return True if objective vector *a* Pareto-dominates *b*.

    A dominates B iff A >= B on all objectives AND A > B on at least one.

    Parameters
    ----------
    a, b : dict[str, float]
        Objective dictionaries with keys: affinity, safety, bioavailability, synthesis.

    Returns
    -------
    bool
    """
    keys = ("affinity", "safety", "bioavailability", "synthesis")
    at_least_one_strict = False
    for k in keys:
        va = a[k]
        vb = b[k]
        if va < vb:
            return False
        if va > vb:
            at_least_one_strict = True
    return at_least_one_strict


def pareto_ranking(molecules: list[dict]) -> list[dict]:
    """Assign Pareto ranks to molecules based on four objectives.

    Iteratively peels Pareto fronts:
      - Non-dominated set = rank 0 (the Pareto front)
      - Remove front, find next non-dominated set = rank 1
      - Continue until all molecules are ranked

    Each molecule receives:
      - ``pareto_rank`` (int): 0 = front, 1 = second layer, etc.
      - ``pareto_front`` (bool): True if pareto_rank == 0.
      - ``pareto_objectives`` (dict): the four normalized objectives.

    Returns molecules sorted by pareto_rank ascending, then by
    composite_score descending within each rank.

    Parameters
    ----------
    molecules : list[dict]
        Scored molecules (must already have properties, ADMET, etc.).

    Returns
    -------
    list[dict]
        Same list, enriched with Pareto data and re-sorted.
    """
    if not molecules:
        return molecules

    n = len(molecules)

    # Step 1: Extract objectives for every molecule
    objectives: list[dict[str, float]] = []
    for mol in molecules:
        obj = _extract_pareto_objectives(mol, n)
        objectives.append(obj)
        mol["pareto_objectives"] = obj

    # Step 2: Iterative Pareto front peeling
    remaining_indices: set[int] = set(range(n))
    rank = 0

    while remaining_indices:
        # Find the non-dominated set among remaining molecules
        non_dominated: list[int] = []
        remaining_list = list(remaining_indices)

        for i in remaining_list:
            dominated = False
            for j in remaining_list:
                if i == j:
                    continue
                if _dominates(objectives[j], objectives[i]):
                    dominated = True
                    break
            if not dominated:
                non_dominated.append(i)

        # Assign rank to this front
        for idx in non_dominated:
            molecules[idx]["pareto_rank"] = rank
            molecules[idx]["pareto_front"] = (rank == 0)
            remaining_indices.discard(idx)

        rank += 1

        # Safety: prevent infinite loops if something goes wrong
        if rank > n:
            logger.warning("Pareto ranking exceeded molecule count, breaking")
            for idx in remaining_indices:
                molecules[idx]["pareto_rank"] = rank
                molecules[idx]["pareto_front"] = False
            break

    # Step 3: Sort by pareto_rank ascending, then composite_score descending
    molecules.sort(
        key=lambda m: (
            m.get("pareto_rank", 999),
            -m.get("composite_score", 0.0),
        )
    )

    front_size = sum(1 for m in molecules if m.get("pareto_front", False))
    logger.info(
        "Pareto ranking: %d molecules, %d on Pareto front, %d ranks total",
        n, front_size, rank,
    )

    return molecules
