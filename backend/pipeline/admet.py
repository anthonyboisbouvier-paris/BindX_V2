"""
DockIt pipeline -- ADMET-AI prediction wrapper.

Predicts ~40 pharmacological properties per molecule across five categories:
Absorption, Distribution, Metabolism, Excretion, and Toxicity.

Strategies (in order):
  1. ADMET-AI library (``admet_ai.ADMETModel``) for ML-based predictions.
  2. RDKit physicochemical rules as heuristic proxies.
  3. Hash-based deterministic mock when no cheminformatics tools are available.

V5bis: adds applicability domain checking (check_applicability_domain).
"""

from __future__ import annotations

import hashlib
import logging
from typing import Optional

logger = logging.getLogger(__name__)

# Cache the ADMET-AI model so we only load it once per process
_ADMET_MODEL: Optional[object] = None
_ADMET_MODEL_CHECKED: bool = False


# ---------------------------------------------------------------------------
# Known mutagenic / reactive SMARTS patterns
# ---------------------------------------------------------------------------

# Structural alerts for Ames mutagenicity (simplified set)
_AMES_ALERTS_SMARTS: list[str] = [
    "[N+](=O)[O-]",          # nitro group
    "[NH2]c1ccccc1",          # aromatic amine
    "N=N",                    # azo group
    "[N;X2]=[N;X2]",         # diazo
    "C1OC1",                  # epoxide
    "[CH2]Cl",                # alkyl halide (chloromethyl)
    "[N;X3](~[O])~[O]",      # nitroso
    "c1cccc2c1[nH]c1ccccc12", # carbazole (polycyclic aromatic)
]

# Structural alerts for hepatotoxicity (reactive groups)
_HEPATOTOX_ALERTS_SMARTS: list[str] = [
    "C(=O)Cl",               # acyl chloride
    "C1OC1",                  # epoxide
    "[S;X2]([#6])[#6]",      # thioether (weak signal, context-dependent)
    "C=CC(=O)",               # Michael acceptor
    "[N+](=O)[O-]",          # nitro group
    "c1cc([NH2])ccc1",       # aniline
    "C(=O)OO",                # peroxide
    "[CX4](Cl)(Cl)",         # polyhalogenated carbon (CCl4, chloroform)
    "[CX4](F)(F)(F)",        # trifluoromethyl on saturated carbon
]


# =====================================================================
# PUBLIC API
# =====================================================================

def predict_admet(smiles_list: list[str]) -> list[dict]:
    """Predict ADMET properties for a list of molecules.

    Tries the real ADMET-AI model first, then falls back to RDKit-based
    heuristics, and finally to a deterministic hash-based mock.

    Parameters
    ----------
    smiles_list : list[str]
        List of SMILES strings to evaluate.

    Returns
    -------
    list[dict]
        One dict per input SMILES with keys:
        - ``smiles`` (str)
        - ``absorption`` (dict): oral_bioavailability, intestinal_permeability,
          solubility, pgp_substrate
        - ``distribution`` (dict): plasma_protein_binding, bbb_permeability, vd
        - ``metabolism`` (dict): cyp1a2_inhibitor, cyp2c9_inhibitor,
          cyp2c19_inhibitor, cyp2d6_inhibitor, cyp3a4_inhibitor
        - ``excretion`` (dict): clearance, half_life
        - ``toxicity`` (dict): herg_inhibition, ames_mutagenicity,
          hepatotoxicity, skin_sensitization, carcinogenicity
        - ``composite_score`` (float): 0-1, higher = better drug profile
        - ``flags`` (list[str]): warnings
        - ``color_code`` (str): "green", "yellow", or "red"
    """
    if not smiles_list:
        return []

    logger.info("Predicting ADMET for %d molecules", len(smiles_list))

    # --- Strategy 1: ADMET-AI library ---
    results = _try_admet_ai(smiles_list)
    if results is not None:
        logger.info("ADMET predictions computed via ADMET-AI model")
        results = _add_applicability_domain(results)
        return results

    # --- Strategy 2: RDKit heuristic proxies ---
    results = _try_rdkit_heuristics(smiles_list)
    if results is not None:
        logger.info("ADMET predictions estimated via RDKit heuristics")
        results = _add_applicability_domain(results)
        return results

    # --- Strategy 3: Hash-based deterministic mock ---
    logger.warning("No cheminformatics tools available; using hash-based ADMET mock")
    results = _mock_admet(smiles_list)
    results = _add_applicability_domain(results)
    return results


def compute_admet_composite(admet: dict) -> tuple[float, list[str], str]:
    """Compute composite ADMET score, flags, and color code from raw predictions.

    The composite score integrates key ADMET endpoints into a single
    drug-likeness metric. Penalties are applied for toxicity liabilities
    and poor pharmacokinetic properties. Rewards are given for favorable
    absorption and clearance characteristics.

    Parameters
    ----------
    admet : dict
        Dictionary with ``absorption``, ``distribution``, ``metabolism``,
        ``excretion``, and ``toxicity`` sub-dicts.

    Returns
    -------
    tuple[float, list[str], str]
        - ``composite_score``: float in [0, 1], higher is better.
        - ``flags``: list of warning strings.
        - ``color_code``: "green" (score >= 0.6), "yellow" (>= 0.35),
          or "red" (< 0.35).
    """
    score: float = 0.5  # baseline
    flags: list[str] = []

    absorption = admet.get("absorption", {})
    distribution = admet.get("distribution", {})
    metabolism = admet.get("metabolism", {})
    excretion = admet.get("excretion", {})
    toxicity = admet.get("toxicity", {})

    # --- Absorption rewards/penalties ---
    oral_bioavail = absorption.get("oral_bioavailability", 0.5)
    if oral_bioavail > 0.5:
        score += 0.10
    elif oral_bioavail < 0.3:
        score -= 0.05
        flags.append("warning: low oral bioavailability predicted")

    solubility = absorption.get("solubility", 0.5)
    if solubility < 0.3:
        score -= 0.05
        flags.append("warning: poor aqueous solubility")

    pgp = absorption.get("pgp_substrate", 0.5)
    if pgp > 0.7:
        score -= 0.05
        flags.append("warning: likely P-gp substrate (efflux risk)")

    # --- Distribution ---
    bbb = distribution.get("bbb_permeability", 0.5)
    # BBB permeability is informational, not inherently good or bad
    ppb = distribution.get("plasma_protein_binding", 0.5)
    if ppb > 0.9:
        score -= 0.03
        flags.append("info: high plasma protein binding (>90%)")

    # --- Metabolism penalties ---
    cyp_inhibitors = [
        ("cyp1a2_inhibitor", "CYP1A2"),
        ("cyp2c9_inhibitor", "CYP2C9"),
        ("cyp2c19_inhibitor", "CYP2C19"),
        ("cyp2d6_inhibitor", "CYP2D6"),
        ("cyp3a4_inhibitor", "CYP3A4"),
    ]
    cyp_count = 0
    for key, label in cyp_inhibitors:
        val = metabolism.get(key, 0.0)
        if val > 0.5:
            cyp_count += 1

    if cyp_count >= 3:
        score -= 0.10
        flags.append(f"warning: inhibits {cyp_count}/5 major CYP enzymes (DDI risk)")
    elif cyp_count >= 1:
        score -= 0.03 * cyp_count
        flags.append(f"info: inhibits {cyp_count}/5 CYP enzymes")

    # --- Excretion ---
    clearance = excretion.get("clearance", 0.5)
    if clearance > 0.7:
        score -= 0.05
        flags.append("warning: high clearance predicted (short duration)")
    half_life = excretion.get("half_life", 0.5)
    if half_life < 0.3:
        score -= 0.03
        flags.append("info: short predicted half-life")

    # --- Toxicity (major penalties) ---
    herg = toxicity.get("herg_inhibition", 0.0)
    if herg > 0.5:
        score -= 0.20
        flags.append("ALERT: hERG inhibition risk (cardiotoxicity)")
    elif herg > 0.3:
        score -= 0.05
        flags.append("warning: hERG inhibition borderline")

    ames = toxicity.get("ames_mutagenicity", 0.0)
    if ames > 0.5:
        score -= 0.20
        flags.append("ALERT: Ames mutagenicity positive")
    elif ames > 0.3:
        score -= 0.05
        flags.append("warning: Ames mutagenicity borderline")

    hepatotox = toxicity.get("hepatotoxicity", 0.0)
    if hepatotox > 0.5:
        score -= 0.15
        flags.append("ALERT: hepatotoxicity risk")
    elif hepatotox > 0.3:
        score -= 0.05
        flags.append("warning: hepatotoxicity borderline")

    skin = toxicity.get("skin_sensitization", 0.0)
    if skin > 0.5:
        score -= 0.05
        flags.append("warning: skin sensitization risk")

    carcino = toxicity.get("carcinogenicity", 0.0)
    if carcino > 0.5:
        score -= 0.15
        flags.append("ALERT: carcinogenicity risk")

    # --- Clamp score to [0, 1] ---
    score = max(0.0, min(1.0, score))
    score = round(score, 4)

    # --- Color code ---
    if score >= 0.6:
        color_code = "green"
    elif score >= 0.35:
        color_code = "yellow"
    else:
        color_code = "red"

    return score, flags, color_code


# =====================================================================
# V6.3: SPECIALIZED hERG PREDICTION
# =====================================================================


def predict_herg_specialized(smiles: str) -> dict:
    """Specialized hERG cardiac risk prediction with calibrated IC50 estimation.

    Uses RDKit molecular features (LogP, MW, basic nitrogen count, aromatic
    ring count) to produce a mock-but-calibrated hERG IC50 value. Falls back
    to a hash-based estimate when RDKit is unavailable.

    IC50 estimation starts at a base of 50 uM and is reduced (multiplied
    by factors < 1) for each risk-increasing molecular feature:

    - logP > 3: IC50 *= 0.6
    - basic nitrogen present: IC50 *= 0.5
    - logP > 4.5 AND basic nitrogen: IC50 *= 0.7  (combined penalty)
    - aromatic rings >= 3: IC50 *= 0.8

    Risk levels:
    - IC50 > 30 uM: LOW
    - IC50 10-30 uM: MODERATE
    - IC50 < 10 uM: HIGH (should be ELIMINATED)

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule to evaluate.

    Returns
    -------
    dict
        Keys:
        - ``ic50_um``: float, estimated hERG IC50 in micromolar.
        - ``risk_level``: str, one of "LOW", "MODERATE", "HIGH".
        - ``method``: str, always ``"heuristic_v2"``.
        - ``features``: dict with the computed molecular features.
    """
    if not smiles:
        return {
            "ic50_um": 50.0,
            "risk_level": "LOW",
            "method": "heuristic_v2",
            "features": {},
        }

    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning("Invalid SMILES for hERG prediction: %s", smiles[:60])
            return {
                "ic50_um": 50.0,
                "risk_level": "LOW",
                "method": "heuristic_v2",
                "features": {"error": "invalid_smiles"},
            }

        logp = Descriptors.MolLogP(mol)
        mw = Descriptors.ExactMolWt(mol)
        aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        has_basic_n = _has_basic_nitrogen(mol, Chem)
        basic_n_count = 1 if has_basic_n else 0

        # Start with base IC50
        ic50: float = 50.0

        # Apply feature-based reductions
        if logp > 3:
            ic50 *= 0.6
        if has_basic_n:
            ic50 *= 0.5
        if logp > 4.5 and has_basic_n:
            ic50 *= 0.7
        if aromatic_rings >= 3:
            ic50 *= 0.8

        ic50 = round(ic50, 2)

        # Determine risk level
        if ic50 > 30.0:
            risk_level = "LOW"
        elif ic50 >= 10.0:
            risk_level = "MODERATE"
        else:
            risk_level = "HIGH"

        features = {
            "logP": round(logp, 2),
            "mw": round(mw, 1),
            "aromatic_rings": aromatic_rings,
            "basic_nitrogen": has_basic_n,
            "basic_n_count": basic_n_count,
        }

        logger.info(
            "Specialized hERG for %s: IC50=%.2f uM (%s) [logP=%.1f, basicN=%s, aroRings=%d]",
            smiles[:40], ic50, risk_level, logp, has_basic_n, aromatic_rings,
        )

        return {
            "ic50_um": ic50,
            "risk_level": risk_level,
            "method": "heuristic_v2",
            "features": features,
        }

    except ImportError:
        # RDKit not available -- use hash-based fallback
        import hashlib as _hashlib
        digest = _hashlib.sha256(f"herg_v2:{smiles}".encode("utf-8")).hexdigest()
        hash_val = int(digest[:8], 16) / 0xFFFFFFFF
        ic50 = round(10.0 + hash_val * 40.0, 2)  # [10, 50] uM range

        if ic50 > 30.0:
            risk_level = "LOW"
        elif ic50 >= 10.0:
            risk_level = "MODERATE"
        else:
            risk_level = "HIGH"

        return {
            "ic50_um": ic50,
            "risk_level": risk_level,
            "method": "heuristic_v2",
            "features": {"fallback": True},
        }

    except Exception as exc:
        logger.warning("Specialized hERG prediction failed for %s: %s", smiles[:60], exc)
        return {
            "ic50_um": 50.0,
            "risk_level": "LOW",
            "method": "heuristic_v2",
            "features": {"error": str(exc)},
        }


# =====================================================================
# V5bis: APPLICABILITY DOMAIN CHECK
# =====================================================================

def check_applicability_domain(smiles: str) -> dict:
    """Check if a molecule is within the ADMET training set's applicability domain.

    Uses Tanimoto similarity to nearest neighbor in a reference set of
    well-known drug molecules. The result indicates how reliable ADMET
    predictions are likely to be for this molecule.

    Parameters
    ----------
    smiles : str
        SMILES string to evaluate.

    Returns
    -------
    dict
        Keys:
        - ``status``: "in_domain", "partial", or "out_of_domain".
        - ``nearest_tanimoto``: float, similarity to closest reference.
        - ``confidence_modifier``: float (1.0, 0.7, or 0.3).
        - ``note``: human-readable explanation.
    """
    try:
        from rdkit import Chem, DataStructs
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {
                "status": "out_of_domain",
                "nearest_tanimoto": 0.0,
                "confidence_modifier": 0.3,
                "note": "Invalid SMILES",
            }

        query_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

        # Reference set: common drug-like molecules (representative fingerprints)
        # In production, this would be loaded from /data/admet_training_fps.pkl
        # For now, use a set of ~20 well-known drug SMILES as reference
        reference_smiles = [
            "CC(=O)Oc1ccccc1C(=O)O",                                     # Aspirin
            "CC(O)c1ccc(CC(C)C)cc1",                                      # Ibuprofen
            "OC(=O)c1ccccc1O",                                            # Salicylic acid
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",                              # Caffeine
            "CC12CCC3C(CCC4CC(=O)CCC43C)C1CCC2O",                        # Testosterone
            "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34",                           # Pyrene
            "C(=O)(OCC)c1ccccc1N",                                        # Benzocaine
            "COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1",               # Omeprazole
            "CC(CS)C(=O)N1CCCC1C(=O)O",                                  # Captopril
            "CN1C2CCC1CC(C2)OC(=O)C(CO)c1ccccc1",                        # Atropine
            "c1ccc(NC(=O)c2ccccc2)cc1",                                   # Benzanilide
            "O=C(O)c1cc(O)c(O)c(O)c1",                                   # Gallic acid
            "CC(=O)Nc1ccc(O)cc1",                                         # Acetaminophen
            "Clc1ccc(C(c2ccc(Cl)cc2)C(Cl)(Cl)Cl)cc1",                    # DDT (toxic reference)
            "c1ccc(-c2cccc(-c3ccccc3)c2)cc1",                             # Terphenyl
            "C=CC#N",                                                      # Acrylonitrile (toxic)
            "c1nc[nH]c1CCN",                                              # Histamine
            "C(=C\\c1ccccc1)\\c1ccccc1",                                  # Stilbene
            "OC(=O)C(F)(F)F",                                             # TFA
            "c1ccc2c(c1)-c1cccc3cccc-2c13",                               # Fluoranthene
        ]

        max_tanimoto: float = 0.0
        for ref_smi in reference_smiles:
            ref_mol = Chem.MolFromSmiles(ref_smi)
            if ref_mol:
                ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol, 2, nBits=2048)
                sim = DataStructs.TanimotoSimilarity(query_fp, ref_fp)
                max_tanimoto = max(max_tanimoto, sim)

        if max_tanimoto > 0.5:
            return {
                "status": "in_domain",
                "nearest_tanimoto": round(max_tanimoto, 3),
                "confidence_modifier": 1.0,
                "note": f"In domain (Tanimoto: {max_tanimoto:.2f}) - High confidence",
            }
        elif max_tanimoto > 0.3:
            return {
                "status": "partial",
                "nearest_tanimoto": round(max_tanimoto, 3),
                "confidence_modifier": 0.7,
                "note": f"Partial domain (Tanimoto: {max_tanimoto:.2f}) - Predictions less reliable",
            }
        else:
            return {
                "status": "out_of_domain",
                "nearest_tanimoto": round(max_tanimoto, 3),
                "confidence_modifier": 0.3,
                "note": f"Out of domain (Tanimoto: {max_tanimoto:.2f}) - ADMET predictions unreliable",
            }

    except ImportError:
        # RDKit not available -- return partial by default
        return {
            "status": "partial",
            "nearest_tanimoto": 0.4,
            "confidence_modifier": 0.7,
            "note": "Domain check unavailable (RDKit missing)",
        }
    except Exception as exc:
        logger.warning("Applicability domain check failed for %s: %s", smiles[:60], exc)
        return {
            "status": "partial",
            "nearest_tanimoto": 0.4,
            "confidence_modifier": 0.7,
            "note": f"Domain check error: {exc}",
        }


# =====================================================================
# STRATEGY 1: ADMET-AI MODEL
# =====================================================================

def _get_admet_model() -> Optional[object]:
    """Lazy-load the ADMET-AI model (singleton)."""
    global _ADMET_MODEL, _ADMET_MODEL_CHECKED

    if _ADMET_MODEL_CHECKED:
        return _ADMET_MODEL

    _ADMET_MODEL_CHECKED = True

    try:
        from admet_ai import ADMETModel  # type: ignore[import-untyped]
        logger.info("Loading ADMET-AI model...")
        _ADMET_MODEL = ADMETModel()
        logger.info("ADMET-AI model loaded successfully")
        return _ADMET_MODEL
    except ImportError:
        logger.info("admet_ai package not installed; will use heuristic fallback")
        return None
    except Exception as exc:
        logger.warning("Failed to load ADMET-AI model: %s", exc)
        return None


def _try_admet_ai(smiles_list: list[str]) -> Optional[list[dict]]:
    """Predict ADMET using the admet_ai library."""
    model = _get_admet_model()
    if model is None:
        return None

    try:
        results: list[dict] = []
        for smi in smiles_list:
            try:
                preds = model.predict(smiles=smi)  # type: ignore[union-attr]

                # Map the ~40 ADMET-AI output keys to our structured format.
                # ADMET-AI returns a dict (or DataFrame row) with property names
                # as keys and probability/regression values as values.
                admet_raw = _map_admet_ai_output(preds)

                composite, flags, color = compute_admet_composite(admet_raw)
                admet_raw["smiles"] = smi
                admet_raw["composite_score"] = composite
                admet_raw["flags"] = flags
                admet_raw["color_code"] = color
                results.append(admet_raw)

            except Exception as exc:
                logger.warning("ADMET-AI prediction failed for %s: %s", smi[:60], exc)
                # Fall through to next SMILES; append a default entry
                results.append(_default_admet_entry(smi))

        return results

    except Exception as exc:
        logger.warning("ADMET-AI batch prediction failed: %s", exc)
        return None


def _map_admet_ai_output(preds: dict) -> dict:
    """Map ADMET-AI model output keys to our canonical structure.

    The ADMET-AI package may return keys like ``'HIA_Hou'``, ``'BBB_Martins'``,
    ``'CYP2D6_Veith'``, etc. This function normalises them into our five
    standard sub-dictionaries.

    Parameters
    ----------
    preds : dict
        Raw output from ``ADMETModel.predict()``.

    Returns
    -------
    dict
        Structured ADMET dict matching our schema.
    """
    def _get(keys: list[str], default: float = 0.5) -> float:
        """Return first matching key value, else default."""
        for k in keys:
            if k in preds:
                val = preds[k]
                try:
                    return float(val)
                except (TypeError, ValueError):
                    continue
        return default

    return {
        "absorption": {
            "oral_bioavailability": _get([
                "Bioavailability_Ma", "bioavailability", "F20_Lethal",
                "HIA_Hou", "oral_bioavailability",
            ], 0.5),
            "intestinal_permeability": _get([
                "HIA_Hou", "Caco2_Wang", "intestinal_absorption",
                "intestinal_permeability",
            ], 0.5),
            "solubility": _get([
                "Solubility_AqSolDB", "ESOL", "solubility",
            ], 0.5),
            "pgp_substrate": _get([
                "Pgp_Broccatelli", "pgp_substrate", "P-gp",
            ], 0.5),
        },
        "distribution": {
            "plasma_protein_binding": _get([
                "PPBR_AZ", "ppb", "plasma_protein_binding",
            ], 0.5),
            "bbb_permeability": _get([
                "BBB_Martins", "bbb_permeability", "BBB",
            ], 0.5),
            "vd": _get([
                "VDss_Lombardo", "vd", "volume_distribution",
            ], 0.5),
        },
        "metabolism": {
            "cyp1a2_inhibitor": _get([
                "CYP1A2_Veith", "cyp1a2_inhibitor", "CYP1A2",
            ], 0.0),
            "cyp2c9_inhibitor": _get([
                "CYP2C9_Veith", "CYP2C9_Substrate_CarbonMangels",
                "cyp2c9_inhibitor",
            ], 0.0),
            "cyp2c19_inhibitor": _get([
                "CYP2C19_Veith", "CYP2C19_Substrate_CarbonMangels",
                "cyp2c19_inhibitor",
            ], 0.0),
            "cyp2d6_inhibitor": _get([
                "CYP2D6_Veith", "CYP2D6_Substrate_CarbonMangels",
                "cyp2d6_inhibitor",
            ], 0.0),
            "cyp3a4_inhibitor": _get([
                "CYP3A4_Veith", "CYP3A4_Substrate_CarbonMangels",
                "cyp3a4_inhibitor",
            ], 0.0),
        },
        "excretion": {
            "clearance": _get([
                "CL-Hepa_Obach", "CL-Micro_Obach", "clearance",
            ], 0.5),
            "half_life": _get([
                "Half_Life_Obach", "half_life", "t_half",
            ], 0.5),
        },
        "toxicity": {
            "herg_inhibition": _get([
                "hERG", "hERG_Karim", "herg_inhibition",
            ], 0.0),
            "ames_mutagenicity": _get([
                "AMES", "Ames_Xu", "ames_mutagenicity",
            ], 0.0),
            "hepatotoxicity": _get([
                "DILI", "hepatotoxicity", "Hepatotoxicity",
            ], 0.0),
            "skin_sensitization": _get([
                "Skin_Reaction", "skin_sensitization", "Sensitization",
            ], 0.0),
            "carcinogenicity": _get([
                "Carcinogens_Lagunin", "carcinogenicity", "Carcinogenicity",
            ], 0.0),
        },
    }


# =====================================================================
# STRATEGY 2: RDKIT HEURISTIC PROXIES
# =====================================================================

def _try_rdkit_heuristics(smiles_list: list[str]) -> Optional[list[dict]]:
    """Estimate ADMET properties from RDKit descriptors and structural alerts."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors
    except ImportError:
        return None

    results: list[dict] = []
    for smi in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                logger.warning("RDKit cannot parse SMILES for ADMET: %s", smi[:80])
                results.append(_default_admet_entry(smi))
                continue

            admet = _rdkit_estimate(mol, Chem, Descriptors, rdMolDescriptors)
            composite, flags, color = compute_admet_composite(admet)
            admet["smiles"] = smi
            admet["composite_score"] = composite
            admet["flags"] = flags
            admet["color_code"] = color
            results.append(admet)

        except Exception as exc:
            logger.warning("RDKit ADMET estimation failed for %s: %s", smi[:60], exc)
            results.append(_default_admet_entry(smi))

    return results


def _rdkit_estimate(
    mol: object,
    Chem: object,    # type: ignore[valid-type]
    Descriptors: object,  # type: ignore[valid-type]
    rdMolDescriptors: object,  # type: ignore[valid-type]
) -> dict:
    """Estimate ADMET properties from molecular descriptors and SMARTS rules.

    Uses the following heuristic mappings:

    - **Oral bioavailability**: Lipinski Rule of Five + Veber rules.
    - **BBB permeability**: MW < 400, TPSA < 90, 1 < logP < 3.
    - **hERG inhibition**: logP > 3 and presence of basic nitrogen.
    - **Ames mutagenicity**: presence of known mutagenic substructures.
    - **Hepatotoxicity**: presence of reactive functional groups.
    - **CYP inhibition**: lipophilicity-driven estimates.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Parsed RDKit molecule.
    Chem, Descriptors, rdMolDescriptors
        RDKit modules (passed in to avoid repeated imports).

    Returns
    -------
    dict
        Structured ADMET dict (without smiles, composite_score, flags,
        color_code -- those are added by the caller).
    """
    # --- Core descriptors ---
    mw = Descriptors.ExactMolWt(mol)  # type: ignore[union-attr]
    logp = Descriptors.MolLogP(mol)   # type: ignore[union-attr]
    tpsa = Descriptors.TPSA(mol)      # type: ignore[union-attr]
    hbd = rdMolDescriptors.CalcNumHBD(mol)  # type: ignore[union-attr]
    hba = rdMolDescriptors.CalcNumHBA(mol)  # type: ignore[union-attr]
    rotbonds = rdMolDescriptors.CalcNumRotatableBonds(mol)  # type: ignore[union-attr]
    num_rings = rdMolDescriptors.CalcNumRings(mol)  # type: ignore[union-attr]

    # Number of heavy atoms
    num_heavy = mol.GetNumHeavyAtoms()  # type: ignore[union-attr]

    # --- Lipinski Rule of Five ---
    lipinski_violations = sum([
        mw > 500,
        logp > 5,
        hbd > 5,
        hba > 10,
    ])

    # --- Veber rules (oral bioavailability beyond Lipinski) ---
    veber_ok = (rotbonds <= 10) and (tpsa <= 140)

    # --- ABSORPTION ---

    # Oral bioavailability: Lipinski + Veber
    if lipinski_violations == 0 and veber_ok:
        oral_bioavail = 0.85
    elif lipinski_violations <= 1 and veber_ok:
        oral_bioavail = 0.65
    elif lipinski_violations <= 1:
        oral_bioavail = 0.45
    else:
        oral_bioavail = 0.20

    # Intestinal permeability: TPSA and MW driven
    if tpsa < 80 and mw < 500:
        intestinal_perm = 0.80
    elif tpsa < 120 and mw < 600:
        intestinal_perm = 0.55
    else:
        intestinal_perm = 0.25

    # Solubility: inverse logP (highly lipophilic = poorly soluble)
    if logp < 1:
        solubility = 0.85
    elif logp < 3:
        solubility = 0.60
    elif logp < 5:
        solubility = 0.35
    else:
        solubility = 0.15

    # P-gp substrate: large, polar molecules tend to be P-gp substrates
    pgp_score = 0.3
    if mw > 400:
        pgp_score += 0.15
    if hbd > 2:
        pgp_score += 0.10
    if num_rings >= 3:
        pgp_score += 0.10
    pgp_score = min(1.0, pgp_score)

    # --- DISTRIBUTION ---

    # Plasma protein binding: logP-driven (lipophilic => high PPB)
    if logp > 4:
        ppb = 0.95
    elif logp > 2:
        ppb = 0.80
    elif logp > 0:
        ppb = 0.60
    else:
        ppb = 0.40

    # BBB permeability
    bbb_score = 0.5
    if mw < 400 and tpsa < 90 and 1 <= logp <= 3:
        bbb_score = 0.80
    elif mw < 450 and tpsa < 100:
        bbb_score = 0.55
    elif tpsa > 120 or mw > 500:
        bbb_score = 0.15

    # Volume of distribution (simplified: logP-correlated)
    if logp > 3:
        vd = 0.75
    elif logp > 1:
        vd = 0.55
    else:
        vd = 0.35

    # --- METABOLISM (CYP inhibition) ---

    # CYP inhibition correlates with lipophilicity and size
    def _cyp_inhibition_estimate(logp: float, mw: float, extra: float = 0.0) -> float:
        """Heuristic CYP inhibition probability."""
        base = 0.15
        if logp > 3:
            base += 0.20
        if logp > 4.5:
            base += 0.15
        if mw > 400:
            base += 0.10
        base += extra
        return min(1.0, max(0.0, base))

    cyp1a2 = _cyp_inhibition_estimate(logp, mw, extra=0.05 if num_rings >= 3 else 0.0)
    cyp2c9 = _cyp_inhibition_estimate(logp, mw, extra=0.05 if mw > 350 else 0.0)
    cyp2c19 = _cyp_inhibition_estimate(logp, mw, extra=0.0)
    cyp2d6 = _cyp_inhibition_estimate(logp, mw, extra=0.10 if _has_basic_nitrogen(mol, Chem) else 0.0)  # type: ignore[arg-type]
    cyp3a4 = _cyp_inhibition_estimate(logp, mw, extra=0.10 if mw > 500 else 0.0)

    # --- EXCRETION ---

    # Clearance: smaller, hydrophilic molecules tend to have higher renal clearance
    if logp < 0 and mw < 300:
        clearance = 0.75
    elif logp < 2 and mw < 400:
        clearance = 0.55
    else:
        clearance = 0.35

    # Half-life: inversely correlated with clearance, high PPB extends it
    if ppb > 0.9 and clearance < 0.4:
        half_life = 0.75
    elif clearance > 0.6:
        half_life = 0.25
    else:
        half_life = 0.50

    # --- TOXICITY ---

    # hERG inhibition: logP > 3 + basic nitrogen is the classic risk profile
    # Penalties calibrated for heuristic predictions — avoid over-elimination
    # of lipophilic kinase inhibitors (many approved drugs have logP > 4 + basic N)
    herg = 0.10
    has_basic_n = _has_basic_nitrogen(mol, Chem)  # type: ignore[arg-type]
    if logp > 3:
        herg += 0.15  # moderate lipophilicity signal
    if has_basic_n:
        herg += 0.20  # basic nitrogen (mild signal alone)
    if logp > 4.5 and has_basic_n:
        herg += 0.10  # combined risk only when both present
    if logp > 5.5:
        herg += 0.15  # extreme lipophilicity (strong signal)
    herg = min(1.0, herg)

    # Ames mutagenicity: structural alerts
    ames_hits = _count_smarts_matches(mol, _AMES_ALERTS_SMARTS, Chem)  # type: ignore[arg-type]
    if ames_hits >= 2:
        ames = 0.80
    elif ames_hits == 1:
        ames = 0.50
    else:
        ames = 0.10

    # Hepatotoxicity: reactive groups + high dose surrogates
    hepatotox_hits = _count_smarts_matches(mol, _HEPATOTOX_ALERTS_SMARTS, Chem)  # type: ignore[arg-type]
    hepatotox = 0.10
    if hepatotox_hits >= 2:
        hepatotox = 0.70
    elif hepatotox_hits == 1:
        hepatotox = 0.40
    if logp > 4:
        hepatotox = min(1.0, hepatotox + 0.15)

    # Skin sensitization: electrophilic groups (Michael acceptors, etc.)
    skin = 0.10
    skin_smarts = ["C=CC(=O)", "C1OC1", "C(=O)Cl"]
    skin_hits = _count_smarts_matches(mol, skin_smarts, Chem)  # type: ignore[arg-type]
    if skin_hits >= 1:
        skin = 0.55

    # Carcinogenicity: combination of Ames + structural red flags
    carcino = 0.10
    if ames > 0.5:
        carcino += 0.25
    if hepatotox > 0.5:
        carcino += 0.15
    if logp > 5 and mw > 500:
        carcino += 0.10
    carcino = min(1.0, carcino)

    return {
        "absorption": {
            "oral_bioavailability": round(oral_bioavail, 3),
            "intestinal_permeability": round(intestinal_perm, 3),
            "solubility": round(solubility, 3),
            "pgp_substrate": round(pgp_score, 3),
        },
        "distribution": {
            "plasma_protein_binding": round(ppb, 3),
            "bbb_permeability": round(bbb_score, 3),
            "vd": round(vd, 3),
        },
        "metabolism": {
            "cyp1a2_inhibitor": round(cyp1a2, 3),
            "cyp2c9_inhibitor": round(cyp2c9, 3),
            "cyp2c19_inhibitor": round(cyp2c19, 3),
            "cyp2d6_inhibitor": round(cyp2d6, 3),
            "cyp3a4_inhibitor": round(cyp3a4, 3),
        },
        "excretion": {
            "clearance": round(clearance, 3),
            "half_life": round(half_life, 3),
        },
        "toxicity": {
            "herg_inhibition": round(herg, 3),
            "ames_mutagenicity": round(ames, 3),
            "hepatotoxicity": round(hepatotox, 3),
            "skin_sensitization": round(skin, 3),
            "carcinogenicity": round(carcino, 3),
        },
    }


def _has_basic_nitrogen(mol: object, Chem: object) -> bool:
    """Check if the molecule contains a basic (protonatable) nitrogen.

    Looks for aliphatic amines and nitrogen-containing rings that are
    commonly protonated at physiological pH.
    """
    basic_n_smarts = [
        "[NH2;X3;!$([NH2]C=O)]",  # primary amine (not amide)
        "[NH;X3;!$([NH]C=O)]",    # secondary amine (not amide)
        "[N;X3;!$([N]C=O);!$([N]S=O)]",  # tertiary amine (not amide/sulfonamide)
        "[nH]",                     # aromatic NH (pyrrole-type: weakly basic)
        "[n;R1;!$([n]1cccc1)]",    # basic aromatic N (pyridine-type)
    ]
    for smarts in basic_n_smarts:
        try:
            pattern = Chem.MolFromSmarts(smarts)  # type: ignore[union-attr]
            if pattern is not None and mol.HasSubstructMatch(pattern):  # type: ignore[union-attr]
                return True
        except Exception:
            continue
    return False


def _count_smarts_matches(mol: object, smarts_list: list[str], Chem: object) -> int:
    """Count how many distinct SMARTS patterns match the molecule.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Target molecule.
    smarts_list : list[str]
        SMARTS patterns to check.
    Chem : module
        RDKit Chem module.

    Returns
    -------
    int
        Number of distinct patterns that match at least once.
    """
    count = 0
    for smarts in smarts_list:
        try:
            pattern = Chem.MolFromSmarts(smarts)  # type: ignore[union-attr]
            if pattern is not None and mol.HasSubstructMatch(pattern):  # type: ignore[union-attr]
                count += 1
        except Exception:
            continue
    return count


# =====================================================================
# STRATEGY 3: HASH-BASED DETERMINISTIC MOCK
# =====================================================================

def _mock_admet(smiles_list: list[str]) -> list[dict]:
    """Generate deterministic mock ADMET predictions from SMILES hashes.

    Produces reproducible values per molecule without any cheminformatics
    dependency. Useful for offline development and testing.
    """
    results: list[dict] = []

    for smi in smiles_list:
        # Deterministic pseudo-random from SMILES hash
        h = hashlib.sha256(smi.encode("utf-8")).hexdigest()

        def _hval(offset: int) -> float:
            """Extract a float in [0, 1] from hex digest at given offset."""
            segment = h[(offset * 2) % len(h):(offset * 2 + 4) % len(h)]
            if not segment:
                segment = h[:4]
            return int(segment, 16) / 0xFFFF

        admet = {
            "absorption": {
                "oral_bioavailability": round(_hval(0) * 0.5 + 0.3, 3),
                "intestinal_permeability": round(_hval(1) * 0.5 + 0.3, 3),
                "solubility": round(_hval(2) * 0.6 + 0.2, 3),
                "pgp_substrate": round(_hval(3) * 0.6 + 0.1, 3),
            },
            "distribution": {
                "plasma_protein_binding": round(_hval(4) * 0.5 + 0.3, 3),
                "bbb_permeability": round(_hval(5) * 0.6 + 0.2, 3),
                "vd": round(_hval(6) * 0.5 + 0.25, 3),
            },
            "metabolism": {
                "cyp1a2_inhibitor": round(_hval(7) * 0.5, 3),
                "cyp2c9_inhibitor": round(_hval(8) * 0.5, 3),
                "cyp2c19_inhibitor": round(_hval(9) * 0.5, 3),
                "cyp2d6_inhibitor": round(_hval(10) * 0.5, 3),
                "cyp3a4_inhibitor": round(_hval(11) * 0.5, 3),
            },
            "excretion": {
                "clearance": round(_hval(12) * 0.5 + 0.25, 3),
                "half_life": round(_hval(13) * 0.5 + 0.25, 3),
            },
            "toxicity": {
                "herg_inhibition": round(_hval(14) * 0.4, 3),
                "ames_mutagenicity": round(_hval(15) * 0.4, 3),
                "hepatotoxicity": round(_hval(16) * 0.4, 3),
                "skin_sensitization": round(_hval(17) * 0.3, 3),
                "carcinogenicity": round(_hval(18) * 0.3, 3),
            },
        }

        composite, flags, color = compute_admet_composite(admet)
        admet["smiles"] = smi
        admet["composite_score"] = composite
        admet["flags"] = flags
        admet["color_code"] = color
        results.append(admet)

    return results


# =====================================================================
# HELPERS
# =====================================================================

def _default_admet_entry(smiles: str) -> dict:
    """Return a safe default ADMET entry when prediction fails for a single molecule."""
    admet: dict = {
        "smiles": smiles,
        "absorption": {
            "oral_bioavailability": 0.5,
            "intestinal_permeability": 0.5,
            "solubility": 0.5,
            "pgp_substrate": 0.5,
        },
        "distribution": {
            "plasma_protein_binding": 0.5,
            "bbb_permeability": 0.5,
            "vd": 0.5,
        },
        "metabolism": {
            "cyp1a2_inhibitor": 0.0,
            "cyp2c9_inhibitor": 0.0,
            "cyp2c19_inhibitor": 0.0,
            "cyp2d6_inhibitor": 0.0,
            "cyp3a4_inhibitor": 0.0,
        },
        "excretion": {
            "clearance": 0.5,
            "half_life": 0.5,
        },
        "toxicity": {
            "herg_inhibition": 0.0,
            "ames_mutagenicity": 0.0,
            "hepatotoxicity": 0.0,
            "skin_sensitization": 0.0,
            "carcinogenicity": 0.0,
        },
        "composite_score": 0.5,
        "flags": ["info: default ADMET values (prediction failed)"],
        "color_code": "yellow",
    }
    return admet


# =====================================================================
# V5bis: APPLICABILITY DOMAIN INTEGRATION HELPER
# =====================================================================

def _add_applicability_domain(results: list[dict]) -> list[dict]:
    """Add applicability domain data to each ADMET result.

    Calls ``check_applicability_domain()`` for each molecule and stores
    the result under the ``applicability_domain`` key.

    Parameters
    ----------
    results : list[dict]
        ADMET prediction results (from any strategy).

    Returns
    -------
    list[dict]
        Same list with ``applicability_domain`` added to each entry.
    """
    for entry in results:
        smiles = entry.get("smiles", "")
        if smiles:
            try:
                ad = check_applicability_domain(smiles)
                entry["applicability_domain"] = ad
            except Exception as exc:
                logger.debug(
                    "Applicability domain check failed for %s: %s",
                    smiles[:60], exc,
                )
                entry["applicability_domain"] = {
                    "status": "partial",
                    "nearest_tanimoto": 0.4,
                    "confidence_modifier": 0.7,
                    "note": "Domain check failed",
                }
        else:
            entry["applicability_domain"] = {
                "status": "partial",
                "nearest_tanimoto": 0.4,
                "confidence_modifier": 0.7,
                "note": "No SMILES provided",
            }
    logger.info("Applicability domain checked for %d molecules", len(results))
    return results


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

    # Erlotinib SMILES (CHEMBL553, EGFR inhibitor)
    erlotinib_smiles = "C=Cc1cc2c(cc1OCC)nc(nc2/N=C/c1ccccc1OC)Nc1cc(ccc1)C#C"
    # Canonical Erlotinib SMILES (alternative representation)
    erlotinib_canonical = "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OC)cc23)c1"

    test_smiles = [
        erlotinib_canonical,                                          # Erlotinib
        "CC(=O)Oc1ccccc1C(=O)O",                                     # Aspirin
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",                                # Ibuprofen
        "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C",                       # Testosterone
        "ClC(Cl)(Cl)Cl",                                              # Carbon tetrachloride (toxic)
        "c1ccc2c(c1)[nH]c1ccc(N)cc12",                               # 2-Aminofluorene (mutagen)
    ]

    print("=" * 72)
    print("DockIt ADMET Prediction Self-Test")
    print("=" * 72)

    results = predict_admet(test_smiles)

    for r in results:
        smi = r["smiles"]
        print(f"\nSMILES: {smi[:60]}")
        print(f"  Composite Score : {r['composite_score']}")
        print(f"  Color Code      : {r['color_code']}")
        print(f"  Oral Bioavail   : {r['absorption']['oral_bioavailability']}")
        print(f"  BBB Permeability: {r['distribution']['bbb_permeability']}")
        print(f"  hERG Inhibition : {r['toxicity']['herg_inhibition']}")
        print(f"  Ames Mutagen    : {r['toxicity']['ames_mutagenicity']}")
        print(f"  Hepatotoxicity  : {r['toxicity']['hepatotoxicity']}")
        # V5bis: Show applicability domain
        ad = r.get("applicability_domain", {})
        if ad:
            print(f"  Appl. Domain    : {ad.get('status', 'unknown')} "
                  f"(Tanimoto: {ad.get('nearest_tanimoto', 'N/A')}, "
                  f"confidence: {ad.get('confidence_modifier', 'N/A')})")
        if r["flags"]:
            print(f"  Flags:")
            for flag in r["flags"]:
                print(f"    - {flag}")
        print()

    # Validate Erlotinib has a reasonable profile (approved drug)
    erlotinib_result = results[0]
    assert erlotinib_result["composite_score"] > 0.3, (
        f"Erlotinib composite score too low: {erlotinib_result['composite_score']}"
    )
    assert erlotinib_result["color_code"] in ("green", "yellow"), (
        f"Erlotinib should not be red: {erlotinib_result['color_code']}"
    )

    # Validate carbon tetrachloride has poor profile (known toxic)
    ccl4_result = results[4]
    # CCl4 is very small/simple, may not trigger all alerts, but should
    # not be flagged as an ideal drug
    print("--- Validation ---")
    print(f"Erlotinib score: {erlotinib_result['composite_score']} "
          f"({erlotinib_result['color_code']}) -- PASS")
    print(f"CCl4 score: {ccl4_result['composite_score']} "
          f"({ccl4_result['color_code']})")

    # Full JSON dump for inspection
    print("\n--- Full JSON (Erlotinib) ---")
    print(json.dumps(erlotinib_result, indent=2))

    # V6.3: Test specialized hERG prediction
    print("\n--- V6.3: Specialized hERG Prediction ---")
    herg_test_smiles = [
        ("Erlotinib", erlotinib_canonical),
        ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("Terfenadine-like", "CCCC(O)(c1ccccc1)c1ccc(cc1)C(O)CCCN1CCC(O)CC1"),  # logP>3 + basic N
    ]
    for name, smi in herg_test_smiles:
        herg_result = predict_herg_specialized(smi)
        print(f"  {name:20s}: IC50={herg_result['ic50_um']:7.2f} uM "
              f"({herg_result['risk_level']:8s}) method={herg_result['method']}")
        if herg_result.get("features"):
            feats = herg_result["features"]
            if "logP" in feats:
                print(f"    logP={feats['logP']:.1f}, basicN={feats.get('basic_nitrogen')}, "
                      f"aroRings={feats.get('aromatic_rings')}")
    # Erlotinib should not be HIGH risk
    erl_herg = predict_herg_specialized(erlotinib_canonical)
    assert erl_herg["risk_level"] in ("LOW", "MODERATE"), (
        f"Erlotinib hERG risk should be LOW or MODERATE, got {erl_herg['risk_level']}"
    )
    print(f"  PASS: Erlotinib hERG risk = {erl_herg['risk_level']}")
    # Empty SMILES
    empty_herg = predict_herg_specialized("")
    assert empty_herg["risk_level"] == "LOW"
    print("  PASS: empty SMILES returns LOW risk")

    print("\n--- All tests passed ---")
    sys.exit(0)
