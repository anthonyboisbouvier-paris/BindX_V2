"""
BindX V9 — Analytics router.

Endpoint: t-SNE chemical space projection for a phase's molecules.
"""

from __future__ import annotations

import logging
from typing import Dict, List
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from auth_v9 import require_v9_user
from database_v9 import get_v9_db
from models_v9 import MoleculeORM_V9
from routers.v9.deps import get_phase_owned

logger = logging.getLogger(__name__)

router = APIRouter()

# In-memory cache: phase_id -> {result, molecule_count}
_tsne_cache: Dict[str, dict] = {}


@router.get("/phases/{phase_id}/analytics/tsne")
async def get_tsne_projection(
    phase_id: UUID,
    user: dict = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Compute t-SNE 2D projection of molecules in a phase using Morgan fingerprints."""

    phase = await get_phase_owned(phase_id, user["sub"], db)

    # Check cache
    cache_key = str(phase_id)
    stmt = select(MoleculeORM_V9.id, MoleculeORM_V9.smiles).where(
        MoleculeORM_V9.phase_id == phase.id
    )
    result = await db.execute(stmt)
    rows = result.all()

    if not rows:
        return {"points": [], "status": "empty"}

    if len(rows) < 5:
        return {"points": [], "status": "too_few", "message": "Need at least 5 molecules for t-SNE"}

    # Check cache validity (same molecule count)
    if cache_key in _tsne_cache and _tsne_cache[cache_key]["molecule_count"] == len(rows):
        return {"points": _tsne_cache[cache_key]["result"], "status": "ok"}

    # Compute t-SNE
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        import numpy as np
        from sklearn.manifold import TSNE
    except ImportError as e:
        logger.warning("t-SNE dependencies unavailable: %s", e)
        return {"points": [], "status": "unavailable", "message": "sklearn or rdkit not available"}

    mol_ids = []
    fps = []
    for row in rows:
        mol = Chem.MolFromSmiles(row.smiles)
        if mol is None:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        fps.append(np.array(fp, dtype=np.float32))
        mol_ids.append(str(row.id))

    if len(fps) < 5:
        return {"points": [], "status": "too_few", "message": "Not enough valid SMILES for t-SNE"}

    try:
        fp_matrix = np.array(fps)
        perplexity = min(30, len(fps) - 1)
        tsne = TSNE(
            n_components=2,
            perplexity=perplexity,
            n_iter=1000,
            random_state=42,
            init="pca",
            learning_rate="auto",
        )
        coords = tsne.fit_transform(fp_matrix)

        points = [
            {"molecule_id": mol_ids[i], "x": float(coords[i, 0]), "y": float(coords[i, 1])}
            for i in range(len(mol_ids))
        ]

        # Cache result
        _tsne_cache[cache_key] = {"result": points, "molecule_count": len(rows)}

        return {"points": points, "status": "ok"}

    except Exception as exc:
        logger.error("t-SNE computation failed: %s", exc)
        return {"points": [], "status": "error", "message": str(exc)}
