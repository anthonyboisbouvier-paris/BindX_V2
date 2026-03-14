"""
BindX pipeline -- Pharmacophore feature extraction.

Extracts 3D pharmacophoric features (H-bond donors/acceptors, hydrophobic
regions, aromatic rings, charge centers) from molecule conformers using RDKit.
Falls back to 2D feature counting when 3D conformer generation fails.
"""

from __future__ import annotations

import logging
from typing import Optional

logger = logging.getLogger(__name__)


def extract_pharmacophore(smiles: str) -> dict:
    """Extract pharmacophore features from a SMILES string.

    Parameters
    ----------
    smiles : str
        Input SMILES string.

    Returns
    -------
    dict
        Keys: smiles, features (list of dicts), feature_counts,
        fingerprint (list of int), n_features.
    """
    try:
        from rdkit import Chem
        return _extract_rdkit(smiles)
    except ImportError:
        logger.warning("RDKit not available; returning empty pharmacophore")
        return _fallback(smiles)


def compute_pharmacophore_similarity(smiles_list: list[str]) -> dict:
    """Compute pairwise pharmacophore similarity for a list of molecules.

    Returns
    -------
    dict
        Keys: similarity_matrix (list of lists), labels (list of str).
    """
    try:
        from rdkit import Chem
        return _similarity_rdkit(smiles_list)
    except ImportError:
        return {"similarity_matrix": [], "labels": smiles_list}


def _fallback(smiles: str) -> dict:
    return {
        "smiles": smiles,
        "features": [],
        "feature_counts": {},
        "fingerprint": [],
        "n_features": 0,
    }


def _extract_rdkit(smiles: str) -> dict:
    """Full pharmacophore extraction using RDKit."""
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdMolChemicalFeatures

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return _fallback(smiles)

    # Generate 3D conformer
    mol_3d = Chem.AddHs(mol)
    conf_ok = False
    try:
        result = AllChem.EmbedMolecule(mol_3d, AllChem.ETKDGv3())
        if result == 0:
            AllChem.MMFFOptimizeMolecule(mol_3d, maxIters=200)
            conf_ok = True
    except Exception:
        pass

    # Use built-in feature factory
    try:
        from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
        factory = Gobbi_Pharm2D.factory
    except Exception:
        factory = None

    features = []
    feature_counts = {
        "donor": 0,
        "acceptor": 0,
        "hydrophobic": 0,
        "aromatic": 0,
        "positive": 0,
        "negative": 0,
    }

    # Extract features using SMARTS patterns
    FEATURE_SMARTS = {
        "donor": ["[NH]", "[NH2]", "[OH]", "[nH]"],
        "acceptor": ["[#7;!$([nH])]", "[O;!$([OH])]", "[F]"],
        "hydrophobic": ["[CH3]", "[CH2]", "c1ccccc1", "[S;!$([S;H1])]"],
        "aromatic": ["a1aaaa1", "a1aaaaa1"],
        "positive": ["[NH3+]", "[NH2+]", "[NX3;H2]", "[nH+]", "[#7;+]"],
        "negative": ["[O-]", "[S-]", "C(=O)[O-]", "S(=O)(=O)[O-]"],
    }

    for feat_type, smarts_list in FEATURE_SMARTS.items():
        for smarts in smarts_list:
            pat = Chem.MolFromSmarts(smarts)
            if pat is None:
                continue
            matches = mol.GetSubstructMatches(pat)
            for match in matches:
                center = list(match)
                position = None
                if conf_ok and mol_3d.GetNumConformers() > 0:
                    conf = mol_3d.GetConformer()
                    # Get centroid of matched atoms
                    try:
                        xs, ys, zs = [], [], []
                        for idx in match:
                            pos = conf.GetAtomPosition(idx)
                            xs.append(pos.x)
                            ys.append(pos.y)
                            zs.append(pos.z)
                        position = [
                            round(sum(xs) / len(xs), 3),
                            round(sum(ys) / len(ys), 3),
                            round(sum(zs) / len(zs), 3),
                        ]
                    except Exception:
                        pass

                features.append({
                    "type": feat_type,
                    "atom_indices": list(match),
                    "position": position,
                    "smarts": smarts,
                })
                feature_counts[feat_type] = feature_counts.get(feat_type, 0) + 1

    # Generate pharmacophore fingerprint
    fingerprint = []
    if factory:
        try:
            fp = Generate.Gen2DFingerprint(mol, factory)
            fingerprint = list(fp.GetOnBits())
        except Exception:
            pass

    return {
        "smiles": smiles,
        "features": features,
        "feature_counts": feature_counts,
        "fingerprint": fingerprint,
        "n_features": len(features),
    }


def _similarity_rdkit(smiles_list: list[str]) -> dict:
    """Compute pairwise pharmacophore fingerprint similarity."""
    from rdkit import Chem, DataStructs
    from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate

    factory = Gobbi_Pharm2D.factory
    fps = []
    valid_smiles = []

    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        try:
            fp = Generate.Gen2DFingerprint(mol, factory)
            fps.append(fp)
            valid_smiles.append(smi)
        except Exception:
            continue

    n = len(fps)
    matrix = [[0.0] * n for _ in range(n)]

    for i in range(n):
        matrix[i][i] = 1.0
        for j in range(i + 1, n):
            try:
                sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            except Exception:
                sim = 0.0
            matrix[i][j] = round(sim, 4)
            matrix[j][i] = round(sim, 4)

    return {
        "similarity_matrix": matrix,
        "labels": valid_smiles,
    }
