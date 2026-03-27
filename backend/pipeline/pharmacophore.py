"""
BindX pipeline -- Reference ligand similarity scoring.

Scores library molecules against a single reference ligand using Gobbi 2D
pharmacophore fingerprints (Tanimoto similarity) and feature-type comparison.

Main entry point:
  - score_against_reference(smiles, ref_smiles) → dict  (called per molecule in a run)
"""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)

# Feature types tracked in pharmacophore models
FEATURE_TYPES = ["donor", "acceptor", "hydrophobic", "aromatic", "positive", "negative"]

FEATURE_SMARTS = {
    "donor": ["[NH]", "[NH2]", "[OH]", "[nH]"],
    "acceptor": ["[#7;!$([nH])]", "[O;!$([OH])]", "[F]"],
    "hydrophobic": ["[CH3]", "[CH2]", "[S;!$([S;H1])]"],
    "aromatic": ["a1aaaa1", "a1aaaaa1"],
    "positive": ["[NH3+]", "[NH2+]", "[nH+]", "[#7;+]"],
    "negative": ["[O-]", "[S-]", "C(=O)[O-]", "S(=O)(=O)[O-]"],
}

# Priority order for SVG coloring: most specific features first so they
# take visual precedence when an atom participates in multiple types.
_SVG_PRIORITY = ["positive", "negative", "donor", "acceptor", "aromatic", "hydrophobic"]


def _get_gobbi_fingerprint(mol):
    """Generate Gobbi 2D pharmacophore fingerprint for an RDKit mol object."""
    from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
    factory = Gobbi_Pharm2D.factory
    fp = Generate.Gen2DFingerprint(mol, factory)
    return fp


def _count_features(mol) -> dict:
    """Count pharmacophore features using RDKit's standard BaseFeatures.fdef.

    Uses MolChemicalFeatureFactory (industry standard) instead of custom SMARTS.
    Each feature = one functional group (not individual atoms).
    """
    from rdkit.Chem import ChemicalFeatures
    from rdkit import RDConfig
    import os

    fdef = os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
    factory = ChemicalFeatures.BuildFeatureFactory(fdef)
    feats = factory.GetFeaturesForMol(mol)

    _FAMILY_MAP = {
        "Donor": "donor",
        "Acceptor": "acceptor",
        "Hydrophobe": "hydrophobic",
        "LumpedHydrophobe": "hydrophobic",
        "Aromatic": "aromatic",
        "PosIonizable": "positive",
        "NegIonizable": "negative",
    }

    counts = {ft: 0 for ft in FEATURE_TYPES}
    for feat in feats:
        mapped = _FAMILY_MAP.get(feat.GetFamily())
        if mapped:
            counts[mapped] += 1
    return counts


def score_against_reference(smiles: str, ref_smiles: str) -> dict:
    """Score a molecule against a single reference ligand.

    Parameters
    ----------
    smiles : str
        SMILES of the molecule to score.
    ref_smiles : str
        SMILES of the reference ligand.

    Returns
    -------
    dict
        Keys: pharmacophore_fit (0-1), feature_counts, ref_feature_counts,
        matched_types (list), missing_types (list).
    """
    from rdkit import Chem, DataStructs
    from rdkit.DataStructs import ExplicitBitVect

    mol = Chem.MolFromSmiles(smiles)
    ref_mol = Chem.MolFromSmiles(ref_smiles)

    if mol is None or ref_mol is None:
        return {
            "pharmacophore_fit": 0.0,
            "feature_counts": {},
            "ref_feature_counts": {},
            "matched_types": [],
            "missing_types": list(FEATURE_TYPES),
        }

    try:
        mol_fp = _get_gobbi_fingerprint(mol)
        ref_fp = _get_gobbi_fingerprint(ref_mol)

        mol_on_bits = set(mol_fp.GetOnBits())
        ref_on_bits = set(ref_fp.GetOnBits())

        if not mol_on_bits and not ref_on_bits:
            pharmacophore_fit = 1.0
        elif not mol_on_bits or not ref_on_bits:
            pharmacophore_fit = 0.0
        else:
            max_bit = max(max(mol_on_bits), max(ref_on_bits)) + 1
            bv_mol = ExplicitBitVect(max_bit)
            for b in mol_on_bits:
                bv_mol.SetBit(b)
            bv_ref = ExplicitBitVect(max_bit)
            for b in ref_on_bits:
                bv_ref.SetBit(b)
            pharmacophore_fit = round(DataStructs.TanimotoSimilarity(bv_mol, bv_ref), 4)
    except Exception as e:
        logger.warning("Pharmacophore scoring failed for %s: %s", smiles, e)
        pharmacophore_fit = 0.0

    mol_features = _count_features(mol)
    ref_features = _count_features(ref_mol)

    # Feature types present in the reference
    ref_types = [ft for ft in FEATURE_TYPES if ref_features.get(ft, 0) > 0]
    matched_types = [ft for ft in ref_types if mol_features.get(ft, 0) > 0]
    missing_types = [ft for ft in ref_types if mol_features.get(ft, 0) == 0]

    return {
        "pharmacophore_fit": pharmacophore_fit,
        "feature_counts": mol_features,
        "ref_feature_counts": ref_features,
        "matched_types": matched_types,
        "missing_types": missing_types,
    }


# Feature type → highlight color (RGB tuples for RDKit MolDraw2DSVG)
_FEATURE_COLORS = {
    "donor":       (0.23, 0.51, 0.96),   # blue
    "acceptor":    (0.94, 0.27, 0.27),   # red
    "hydrophobic": (0.98, 0.75, 0.15),   # yellow
    "aromatic":    (0.66, 0.33, 0.97),   # violet
    "positive":    (0.06, 0.73, 0.51),   # green
    "negative":    (0.98, 0.47, 0.22),   # orange
}


def generate_feature_map_svg(smiles: str, width: int = 400, height: int = 250) -> str:
    """Generate a clean 2D structure SVG for a molecule using rdMolDraw2D."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.Draw import rdMolDraw2D

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    AllChem.Compute2DCoords(mol)

    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    opts.clearBackground = True
    opts.bondLineWidth = 1.5

    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()
