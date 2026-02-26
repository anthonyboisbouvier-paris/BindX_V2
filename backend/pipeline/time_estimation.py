"""
DockIt pipeline — Time estimation for pipeline steps.

V9: Provides realistic time estimates for each pipeline step based on
input parameters (number of ligands, mode, options). Helps the user
understand expected runtime when configuring a screening run.

Estimates are based on empirical benchmarks from the DockIt pipeline:
- Structure retrieval: AlphaFold DB ~3s, ESMFold ~20-60s, PDB ~2s
- Pocket detection: fpocket ~2-10s depending on protein size
- Receptor preparation: Open Babel ~2-5s
- Ligand fetching: ChEMBL ~5-15s, PubChem ~5-20s, Enamine ~1-3s
- Docking (Vina): ~5-30s per ligand (CPU, no GPU)
- Docking (GNINA): ~10-45s per ligand
- ADMET predictions: ~0.5-1s per molecule
- Retrosynthesis: ~10-30s per molecule
- Off-target screening: ~2-5s per molecule
- Report generation: ~5-15s
- AI agents (OpenAI): ~3-10s per call
"""

from __future__ import annotations

import logging
from typing import Optional

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Per-step time constants (seconds) — empirical estimates, CPU-only
# ---------------------------------------------------------------------------

# Structure retrieval
TIME_ALPHAFOLD_FETCH = 3.0       # AlphaFold DB HTTP fetch
TIME_PDB_FETCH = 2.0             # RCSB PDB fetch
TIME_ESMFOLD_PREDICT = 40.0      # ESMFold via HuggingFace (varies with length)
TIME_ESMFOLD_PER_RESIDUE = 0.1   # Additional per-residue ESMFold cost

# Pocket detection
TIME_FPOCKET_BASE = 3.0          # fpocket base cost
TIME_FPOCKET_PER_RESIDUE = 0.005 # Additional per-residue

# Receptor preparation
TIME_PREPARE_RECEPTOR = 3.0      # Open Babel PDB -> PDBQT

# Ligand fetching
TIME_CHEMBL_QUERY = 10.0         # ChEMBL API query + parsing
TIME_PUBCHEM_QUERY = 12.0        # PubChem PUG REST
TIME_ZINC_LOCAL = 1.0            # Local SDF loading
TIME_ENAMINE_SAMPLE = 2.0        # Enamine REAL combinatorial generation
TIME_USER_SMILES_PARSE = 0.5     # User SMILES validation

# Docking per ligand
TIME_VINA_PER_LIGAND = 15.0      # AutoDock Vina (avg, CPU)
TIME_GNINA_PER_LIGAND = 25.0     # GNINA with CNN scoring
TIME_DIFFDOCK_PER_LIGAND = 30.0  # DiffDock AI docking

# Post-docking per molecule
TIME_SCORING_PER_MOL = 0.3       # RDKit properties + composite score
TIME_ADMET_PER_MOL = 0.8         # ADMET prediction per molecule
TIME_RETROSYNTHESIS_PER_MOL = 20.0  # Retrosynthesis planning
TIME_OFF_TARGET_PER_MOL = 3.0    # Off-target screening
TIME_HERG_PER_MOL = 1.0          # Specialized hERG screening
TIME_INTERACTION_PER_MOL = 2.0   # Interaction analysis

# Generation
TIME_GENERATION_BASE = 5.0       # Molecule generation setup
TIME_GENERATION_PER_MOL = 0.5    # Per generated molecule

# Report
TIME_REPORT_PDF = 8.0            # PDF report generation
TIME_REPORT_ZIP = 3.0            # ZIP archive creation

# Confidence + clustering
TIME_CONFIDENCE_PER_MOL = 0.1
TIME_CLUSTERING_BASE = 2.0
TIME_CLUSTERING_PER_MOL = 0.05

# AI Agents (OpenAI GPT-4o)
TIME_AGENT_CALL = 6.0            # Single agent call (avg)

# Disorder prediction
TIME_DISORDER_PREDICT = 2.0


def estimate_pipeline_time(
    n_ligands: int = 50,
    mode: str = "rapid",
    docking_engine: str = "auto",
    enable_generation: bool = False,
    n_generated_molecules: int = 100,
    enable_retrosynthesis: bool = False,
    enable_diffdock: bool = False,
    structure_source: str = "alphafold",
    use_pubchem: bool = False,
    use_enamine: bool = False,
    use_chembl: bool = True,
    include_agents: bool = False,
    protein_length: int = 400,
) -> dict:
    """Estimate total pipeline runtime and per-step breakdown.

    Parameters
    ----------
    n_ligands : int
        Number of ligands to dock.
    mode : str
        Pipeline mode: "rapid", "standard", or "deep".
    docking_engine : str
        Docking engine: "vina", "gnina", "auto".
    enable_generation : bool
        Whether AI molecule generation is enabled.
    n_generated_molecules : int
        Number of molecules to generate.
    enable_retrosynthesis : bool
        Whether retrosynthesis planning is enabled.
    enable_diffdock : bool
        Whether DiffDock AI docking is used.
    structure_source : str
        Expected structure source.
    use_pubchem : bool
        Whether PubChem will be queried.
    use_enamine : bool
        Whether Enamine REAL sampling is used.
    use_chembl : bool
        Whether ChEMBL will be queried.
    include_agents : bool
        Whether AI agents (OpenAI) will be called.
    protein_length : int
        Approximate protein sequence length (residues).

    Returns
    -------
    dict
        Keys: ``total_seconds``, ``total_human`` (formatted string),
        ``steps`` (list of step dicts with ``name``, ``seconds``, ``detail``),
        ``confidence`` ("high"|"medium"|"low").
    """
    is_standard = mode in ("standard", "advanced")
    is_deep = mode == "deep"
    steps: list[dict] = []

    # 1. Structure retrieval
    if structure_source in ("alphafold", "alphafold_db"):
        struct_time = TIME_ALPHAFOLD_FETCH
        struct_detail = "AlphaFold DB fetch"
    elif structure_source in ("pdb", "pdb_experimental", "pdb_holo"):
        struct_time = TIME_PDB_FETCH
        struct_detail = "RCSB PDB fetch"
    elif structure_source == "esmfold":
        struct_time = TIME_ESMFOLD_PREDICT + protein_length * TIME_ESMFOLD_PER_RESIDUE
        struct_detail = f"ESMFold prediction ({protein_length} residues)"
    else:
        struct_time = TIME_ALPHAFOLD_FETCH
        struct_detail = "Structure retrieval"
    steps.append({"name": "Structure Retrieval", "seconds": struct_time, "detail": struct_detail})

    # 1b. Disorder prediction
    steps.append({"name": "Disorder Prediction", "seconds": TIME_DISORDER_PREDICT, "detail": "IUPred-style disorder scan"})

    # 2. Pocket detection
    pocket_time = TIME_FPOCKET_BASE + protein_length * TIME_FPOCKET_PER_RESIDUE
    steps.append({"name": "Pocket Detection", "seconds": pocket_time, "detail": f"fpocket on {protein_length}-residue protein"})

    # 3. Receptor preparation
    steps.append({"name": "Receptor Preparation", "seconds": TIME_PREPARE_RECEPTOR, "detail": "PDB to PDBQT conversion"})

    # 4. Ligand fetching
    ligand_time = 0.0
    ligand_sources: list[str] = []
    if use_chembl:
        ligand_time += TIME_CHEMBL_QUERY
        ligand_sources.append("ChEMBL")
    if use_pubchem:
        ligand_time += TIME_PUBCHEM_QUERY
        ligand_sources.append("PubChem")
    if use_enamine:
        ligand_time += TIME_ENAMINE_SAMPLE
        ligand_sources.append("Enamine REAL")
    ligand_time += TIME_ZINC_LOCAL  # always loaded as fallback
    ligand_time += TIME_USER_SMILES_PARSE
    steps.append({
        "name": "Ligand Collection",
        "seconds": ligand_time,
        "detail": f"Querying {', '.join(ligand_sources) or 'databases'} for {n_ligands} molecules",
    })

    # 5. Docking
    if enable_diffdock and is_standard:
        dock_per = TIME_DIFFDOCK_PER_LIGAND
        dock_engine_name = "DiffDock"
    elif docking_engine == "gnina":
        dock_per = TIME_GNINA_PER_LIGAND
        dock_engine_name = "GNINA"
    else:
        dock_per = TIME_VINA_PER_LIGAND
        dock_engine_name = "AutoDock Vina"

    if is_deep:
        # Deep mode: 5 passes
        dock_time = dock_per * n_ligands
        steps.append({
            "name": "Massive Docking (5 passes)",
            "seconds": dock_time,
            "detail": f"{n_ligands} ligands x {dock_per:.0f}s/ligand ({dock_engine_name})",
        })
    else:
        dock_time = dock_per * n_ligands
        steps.append({
            "name": "Molecular Docking",
            "seconds": dock_time,
            "detail": f"{n_ligands} ligands x {dock_per:.0f}s/ligand ({dock_engine_name})",
        })

    # 6. Generation (standard/deep mode)
    if (is_standard or is_deep) and enable_generation:
        gen_time = TIME_GENERATION_BASE + n_generated_molecules * TIME_GENERATION_PER_MOL
        gen_dock_time = dock_per * min(n_generated_molecules, 30)  # top N re-docked
        steps.append({
            "name": "AI Molecule Generation",
            "seconds": gen_time + gen_dock_time,
            "detail": f"Generate {n_generated_molecules} + dock top 30",
        })

    # 7. ADMET (standard/deep)
    total_mols = n_ligands + (n_generated_molecules if enable_generation else 0)
    if is_standard or is_deep:
        admet_time = TIME_ADMET_PER_MOL * total_mols
        steps.append({
            "name": "ADMET Predictions",
            "seconds": admet_time,
            "detail": f"{total_mols} molecules",
        })

    # 8. Scoring + clustering
    score_time = TIME_SCORING_PER_MOL * total_mols + TIME_CLUSTERING_BASE + TIME_CLUSTERING_PER_MOL * total_mols
    steps.append({
        "name": "Scoring & Clustering",
        "seconds": score_time,
        "detail": f"RDKit properties + Butina clustering for {total_mols} molecules",
    })

    # 9. Off-target (standard/deep) — top 5
    if is_standard or is_deep:
        ot_time = TIME_OFF_TARGET_PER_MOL * 5 + TIME_HERG_PER_MOL * 20
        steps.append({
            "name": "Safety Screening",
            "seconds": ot_time,
            "detail": "Off-target (top 5) + hERG (top 20)",
        })

    # 10. Retrosynthesis (standard/deep) — top 5
    if (is_standard or is_deep) and enable_retrosynthesis:
        retro_time = TIME_RETROSYNTHESIS_PER_MOL * 5
        steps.append({
            "name": "Retrosynthesis Planning",
            "seconds": retro_time,
            "detail": "Top 5 candidates",
        })

    # 11. Confidence + interactions
    conf_time = TIME_CONFIDENCE_PER_MOL * total_mols + TIME_INTERACTION_PER_MOL * 10
    steps.append({
        "name": "Confidence & Interactions",
        "seconds": conf_time,
        "detail": f"Confidence for {total_mols} + interactions for top 10",
    })

    # 12. Report
    steps.append({
        "name": "Report Generation",
        "seconds": TIME_REPORT_PDF + TIME_REPORT_ZIP,
        "detail": "PDF report + ZIP archive",
    })

    # 13. AI Agents (optional)
    if include_agents:
        agent_time = TIME_AGENT_CALL * 2  # target + run analysis
        steps.append({
            "name": "AI Agent Analysis",
            "seconds": agent_time,
            "detail": "GPT-4o scientific advisory",
        })

    # Compute totals
    total_seconds = sum(s["seconds"] for s in steps)

    # Round step seconds
    for s in steps:
        s["seconds"] = round(s["seconds"], 1)

    # Determine confidence
    if n_ligands <= 100 and not is_deep:
        confidence = "high"
    elif n_ligands <= 500:
        confidence = "medium"
    else:
        confidence = "low"

    return {
        "total_seconds": round(total_seconds, 1),
        "total_human": _format_duration(total_seconds),
        "steps": steps,
        "confidence": confidence,
        "mode": mode,
        "n_ligands": n_ligands,
    }


def _format_duration(seconds: float) -> str:
    """Format seconds into a human-readable duration string."""
    if seconds < 60:
        return f"{seconds:.0f} seconds"
    elif seconds < 3600:
        minutes = seconds / 60
        return f"{minutes:.0f} min"
    else:
        hours = seconds / 3600
        remaining_min = (seconds % 3600) / 60
        if remaining_min > 5:
            return f"{hours:.0f}h {remaining_min:.0f}min"
        return f"{hours:.1f} hours"


# ---------------------------------------------------------------------------
# Quick estimate for display in the UI run config form
# ---------------------------------------------------------------------------

def quick_estimate(
    n_ligands: int = 50,
    mode: str = "rapid",
    docking_engine: str = "auto",
) -> dict:
    """Simplified time estimate for the run configuration UI.

    Returns a compact estimate with min/max range.

    Parameters
    ----------
    n_ligands : int
        Number of ligands.
    mode : str
        Pipeline mode.
    docking_engine : str
        Docking engine.

    Returns
    -------
    dict
        Keys: ``min_seconds``, ``max_seconds``, ``estimate_human``,
        ``breakdown`` (list of ``{step, time}``).
    """
    # Base estimates per mode
    if mode == "rapid":
        overhead = 30  # structure + pockets + prepare + scoring + report
        if docking_engine == "gnina":
            per_ligand = 25
        else:
            per_ligand = 15
        min_s = overhead + per_ligand * n_ligands * 0.7
        max_s = overhead + per_ligand * n_ligands * 1.5

    elif mode in ("standard", "advanced"):
        overhead = 60  # structure + pockets + prepare + ADMET + offtarget + scoring + report
        if docking_engine == "gnina":
            per_ligand = 28
        else:
            per_ligand = 18
        min_s = overhead + per_ligand * n_ligands * 0.7
        max_s = overhead + per_ligand * n_ligands * 1.5

    elif mode == "deep":
        overhead = 120
        per_ligand = 20
        min_s = overhead + per_ligand * n_ligands * 0.7
        max_s = overhead + per_ligand * n_ligands * 1.5

    else:
        overhead = 30
        per_ligand = 15
        min_s = overhead + per_ligand * n_ligands * 0.7
        max_s = overhead + per_ligand * n_ligands * 1.5

    avg_s = (min_s + max_s) / 2

    breakdown = [
        {"step": "Setup", "time": _format_duration(overhead)},
        {"step": "Docking", "time": _format_duration(per_ligand * n_ligands)},
    ]
    if mode in ("standard", "advanced", "deep"):
        breakdown.append({"step": "Analysis", "time": _format_duration(overhead * 0.5)})

    return {
        "min_seconds": round(min_s),
        "max_seconds": round(max_s),
        "estimate_human": _format_duration(avg_s),
        "breakdown": breakdown,
    }
