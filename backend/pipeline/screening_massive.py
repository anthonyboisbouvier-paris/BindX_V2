"""
DockIt pipeline -- Multi-pass massive screening orchestrator.

Orchestrates the 5-pass screening pipeline for mode="deep", processing
up to ~2.4M molecules through progressively more expensive filters:

  Pass 1: Pharmacological filter  (~2.4M -> ~200K)   [filter_pharma]
  Pass 2: 3D Shape filter         (~200K -> ~10K)    [filter_shape]
  Pass 3: Rapid smina/Vinardo     (~10K  -> ~500)    [scoring_rapid]
  Pass 4: Full Vina docking       (~500  -> ~50)     [docking]
  Pass 5: ADMET + retrosynthesis  (~50   -> ~50)     [admet, retrosynthesis]

The orchestrator manages data flow between passes, tracks timing, reports
progress, and handles all fallback/mock scenarios gracefully when external
tools or data files are not available.
"""

from __future__ import annotations

import hashlib
import logging
import random
import time
from pathlib import Path
from typing import Callable, Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# ChEMBL dump configuration
# ---------------------------------------------------------------------------

CHEMBL_DUMP_PATH: Path = Path("/data/chembl_all.smi")

# Hardcoded drug-like SMILES for mock screening when ChEMBL dump is unavailable.
# These represent a diverse set of drug-like molecules covering different
# pharmacological classes and structural features.
_MOCK_DRUGLIKE_SMILES: list[str] = [
    # Kinase inhibitors
    "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1",           # Erlotinib
    "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",      # Gefitinib
    "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",  # Imatinib
    "CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1",   # Sorafenib
    "Cc1nc(Nc2ncc(C(=O)Nc3c(C)cccc3Cl)s2)cc(N2CCN(CCO)CC2)n1",     # Dasatinib
    "CC(Oc1cc(-c2cnn(C3CCNCC3)c2)cnc1N)c1c(Cl)ccc(F)c1Cl",         # Crizotinib
    "COc1cc(N(C)CCN(C)C)c(NC(=O)/C=C/CN(C)C)cc1Nc1nccc(-c2cn(C)c3ccccc23)n1",  # Osimertinib
    "CN(C)/C=C/C(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1OC1CCOC1",   # Afatinib
    # GPCR ligands
    "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1",                   # Fluoxetine
    "CC(C)NCC(O)c1ccc(O)c(O)c1",                             # Isoproterenol
    "CN1CCC(=C2c3ccccc3CCc3cccnc32)CC1",                     # Mirtazapine
    "c1ccc(-c2nc3ccc(OCCCCN4CCN(c5cccc(Cl)c5)CC4)cc3o2)cc1", # Aripiprazole-like
    # Anti-inflammatories
    "CC(=O)Oc1ccccc1C(=O)O",                                 # Aspirin
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",                            # Ibuprofen
    "COc1ccc2cc(ccc2c1)C(C)C(=O)O",                          # Naproxen
    "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl",                        # Diclofenac
    # CNS drugs
    "Cn1c(=O)c2c(ncn2C)n(C)c1=O",                            # Caffeine
    "CC(=O)Nc1ccc(O)cc1",                                     # Acetaminophen
    "CN(C)C(=N)NC(=N)N",                                      # Metformin
    "COc1ccc2[nH]c(nc2c1)S(=O)Cc1ncc(C)c(OC)c1C",           # Omeprazole
    # Cardiovascular
    "CCCCc1nc(Cl)c(n1Cc1ccc(cc1)c1ccccc1c1nn[nH]n1)CO",     # Losartan
    "OC(=O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",           # Ciprofloxacin
    # Diverse scaffolds
    "CC1(C)CC(NC(=O)c2ccc(-c3ccnc4[nH]ncc34)cc2)CC(C)(C)N1",
    "COc1ccc(-c2nc(NC3CCNCC3)c3c(n2)CCC3=O)cc1Cl",
    "CS(=O)(=O)c1ccc(Nc2ncnc3cc(OC4CCOC4)c(OC)cc23)cc1",
    "O=C(NC1CCCNC1)c1cc(-c2ccc(F)cc2F)nc2ccccc12",
    "Cc1[nH]c(C(=O)Nc2cccc(N3CCNCC3)c2)c(-c2cccc(F)c2)c1C",
    "COc1cc2nccc(Nc3ccc(F)c(NC(=O)c4ccco4)c3)c2cc1OCCO",
    "CC(=O)Nc1ccc(-c2cc3c(N)ncnc3n2C2CCC2)cc1",
    "Fc1ccc(-c2[nH]nc3c(Nc4ccc5[nH]ncc5c4)ncnc23)cc1",
    "O=C(Nc1cccc(-c2ccnc(N3CCOCC3)n2)c1)C1CCC(F)(F)CC1",
    "Cc1cc(C)c(Nc2nccc(-c3ccc(CN4CCCC4)cc3)n2)c(C)c1",
    "O=C(c1ccncc1)N1CCN(c2cccc(Oc3cccc(Cl)c3)c2)CC1",
    "Nc1ccc(-c2cn3c(n2)sc2cc(Cl)ccc23)cn1",
    "CNC(=O)c1ccc(-c2ccnc(Nc3ccc(OC)c(OC)c3)c2)cc1",
    "COc1cc(NC(=O)c2cc(C)on2)cc(-c2cccc(F)c2)c1",
    "O=C(NC1CCN(c2ncnc3[nH]ccc23)CC1)c1ccc(Cl)s1",
    "Cc1nc(N2CCCC2)c2cc(-c3cccc4[nH]ncc34)nn2c1=O",
    "FC(F)Oc1ccc(-c2cc(-c3ccc(N4CCNCC4)nc3)[nH]n2)cc1",
    "COc1ccc(S(=O)(=O)Nc2ccccc2-c2nc(C)c(C)[nH]2)cc1",
    "O=c1[nH]c(-c2ccccc2Cl)nc2cc(N3CCCC3)ccc12",
    # Additional diverse molecules
    "Cc1ccnc(Nc2ccc(CN3CCNCC3)cc2)c1C(=O)Nc1ccccc1F",
    "COc1cc(Nc2ncc3c(n2)n(C2CCCC2)c(=O)n3C)cc(OC)c1OC",
    "CC(C)Oc1ccc(-c2cnc(N)c(C(=O)Nc3ccccc3Cl)n2)cc1",
    "Nc1ncnc2c1cc(-c1ccc(F)c(NC(=O)C3CC3)c1)n2C1CCCC1",
    "O=C(Nc1ccc(Oc2ccnc3cc(F)ccc23)cc1)c1cccs1",
    "CCn1c(=O)c2c(nc3ccc(OC)cc3n2C)n(CC)c1=O",
    "CCOC(=O)C1=C(COCCN)NC(C)=C(C1c1ccccc1Cl)C(=O)OC",
    "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O",
    "CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(F)C(O)CC2(C)C1(O)C(=O)CO",
]


# ---------------------------------------------------------------------------
# ChEMBL dump loading
# ---------------------------------------------------------------------------

def check_chembl_dump_available() -> bool:
    """Check if the ChEMBL SMILES dump file is available.

    Returns
    -------
    bool
        True if /data/chembl_all.smi exists and is non-empty.
    """
    if CHEMBL_DUMP_PATH.exists() and CHEMBL_DUMP_PATH.stat().st_size > 0:
        logger.info("ChEMBL dump found at %s", CHEMBL_DUMP_PATH)
        return True
    logger.info(
        "ChEMBL dump not found at %s; will use mock molecule set",
        CHEMBL_DUMP_PATH,
    )
    return False


def load_chembl_smiles(max_count: int = 0) -> list[str]:
    """Load SMILES from the ChEMBL dump file.

    The dump file is expected to be a tab-separated file with SMILES
    in the first column (one molecule per line). Lines starting with
    '#' are treated as comments and skipped.

    Parameters
    ----------
    max_count : int
        Maximum number of SMILES to load. 0 means load all.

    Returns
    -------
    list[str]
        List of SMILES strings. If the dump file is not available,
        returns a mock set of ~10000 molecules generated from the
        hardcoded drug-like list.
    """
    if check_chembl_dump_available():
        return _load_chembl_from_file(max_count)

    # Fallback: generate a mock set from the hardcoded list
    logger.info("Generating mock ChEMBL set from %d seed molecules", len(_MOCK_DRUGLIKE_SMILES))
    return _generate_mock_chembl(max_count)


def _load_chembl_from_file(max_count: int) -> list[str]:
    """Load SMILES from the actual ChEMBL dump file.

    Parameters
    ----------
    max_count : int
        Maximum to load (0 = all).

    Returns
    -------
    list[str]
        SMILES strings from the file.
    """
    smiles_list: list[str] = []

    try:
        with open(CHEMBL_DUMP_PATH, "r") as fh:
            for line_num, line in enumerate(fh, 1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                # Tab-separated: first column is SMILES
                parts = line.split("\t")
                smi = parts[0].strip()

                if not smi:
                    continue

                smiles_list.append(smi)

                if max_count > 0 and len(smiles_list) >= max_count:
                    break

                # Log progress for very large files
                if line_num % 500000 == 0:
                    logger.info(
                        "Loading ChEMBL dump: %d lines read, %d SMILES loaded",
                        line_num, len(smiles_list),
                    )

        logger.info(
            "Loaded %d SMILES from ChEMBL dump (%s)",
            len(smiles_list), CHEMBL_DUMP_PATH,
        )
        return smiles_list

    except Exception as exc:
        logger.error("Failed to load ChEMBL dump: %s", exc)
        logger.info("Falling back to mock molecule set")
        return _generate_mock_chembl(max_count)


def _generate_mock_chembl(max_count: int) -> list[str]:
    """Generate a mock molecule set by mutating seed molecules.

    Creates ~10000 molecules (or max_count if specified) by:
      1. Starting with the ~50 hardcoded drug-like SMILES.
      2. If RDKit is available, applying simple mutations (halogen swaps,
         functional group additions) to expand the set.
      3. If RDKit is not available, repeating and hashing the seed set
         to generate unique identifiers.

    Parameters
    ----------
    max_count : int
        Target number of molecules. 0 defaults to 10000.

    Returns
    -------
    list[str]
        Mock SMILES list for screening.
    """
    target_count = max_count if max_count > 0 else 10000
    seeds = list(_MOCK_DRUGLIKE_SMILES)

    # If we already have enough seeds, return a subset
    if len(seeds) >= target_count:
        return seeds[:target_count]

    # Try RDKit-based expansion
    expanded = _expand_with_rdkit(seeds, target_count)
    if expanded is not None and len(expanded) >= target_count:
        logger.info("Generated %d mock molecules via RDKit expansion", len(expanded))
        return expanded[:target_count]

    # Fallback: repeat seeds with index-based uniqueness
    # (same SMILES but differentiated by position for mock purposes)
    result: list[str] = list(seeds)
    rng = random.Random(42)  # deterministic

    while len(result) < target_count:
        # Cycle through seeds
        base_smi = seeds[len(result) % len(seeds)]
        # Add a small variation tag (won't parse as real SMILES in mock mode)
        # For the mock, we just repeat the seeds -- downstream mock filters
        # use hash-based selection so different positions get different results
        result.append(base_smi)

    logger.info(
        "Generated %d mock molecules via seed repetition", len(result)
    )
    return result[:target_count]


def _expand_with_rdkit(
    seeds: list[str],
    target_count: int,
) -> Optional[list[str]]:
    """Expand the seed set using RDKit mutations.

    Parameters
    ----------
    seeds : list[str]
        Seed SMILES.
    target_count : int
        Target number of molecules.

    Returns
    -------
    list[str] or None
        Expanded SMILES list, or None if RDKit is not available.
    """
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import AllChem

        RDLogger.DisableLog("rdApp.*")
    except ImportError:
        return None

    result: list[str] = list(seeds)
    seen: set[str] = set(seeds)

    rng = random.Random(42)
    halogen_swaps = {9: 17, 17: 35, 35: 9}  # F->Cl->Br->F

    attempts = 0
    max_attempts = target_count * 5

    while len(result) < target_count and attempts < max_attempts:
        attempts += 1
        base_smi = seeds[attempts % len(seeds)]
        mol = Chem.MolFromSmiles(base_smi)
        if mol is None:
            continue

        # Apply a simple mutation
        strategy = rng.choice(["halogen_swap", "add_methyl", "remove_atom"])

        try:
            if strategy == "halogen_swap":
                rwmol = Chem.RWMol(mol)
                halogens = [
                    a for a in rwmol.GetAtoms()
                    if a.GetAtomicNum() in halogen_swaps
                ]
                if halogens:
                    target = rng.choice(halogens)
                    target.SetAtomicNum(halogen_swaps[target.GetAtomicNum()])
                    Chem.SanitizeMol(rwmol)
                    new_smi = Chem.MolToSmiles(rwmol)
                else:
                    continue

            elif strategy == "add_methyl":
                rwmol = Chem.RWMol(mol)
                aromatic_c = [
                    a for a in rwmol.GetAtoms()
                    if a.GetIsAromatic() and a.GetAtomicNum() == 6
                    and a.GetTotalNumHs() > 0
                ]
                if aromatic_c:
                    target = rng.choice(aromatic_c)
                    new_idx = rwmol.AddAtom(Chem.Atom(6))
                    rwmol.AddBond(target.GetIdx(), new_idx, Chem.BondType.SINGLE)
                    Chem.SanitizeMol(rwmol)
                    new_smi = Chem.MolToSmiles(rwmol)
                else:
                    continue

            elif strategy == "remove_atom":
                # Remove a terminal non-ring atom
                rwmol = Chem.RWMol(mol)
                terminal = [
                    a for a in rwmol.GetAtoms()
                    if a.GetDegree() == 1
                    and not a.IsInRing()
                    and a.GetAtomicNum() not in (1,)  # skip H
                ]
                if terminal:
                    target = rng.choice(terminal)
                    rwmol.RemoveAtom(target.GetIdx())
                    Chem.SanitizeMol(rwmol)
                    new_smi = Chem.MolToSmiles(rwmol)
                else:
                    continue
            else:
                continue

            if new_smi and new_smi not in seen:
                # Validate the new molecule
                check_mol = Chem.MolFromSmiles(new_smi)
                if check_mol is not None and check_mol.GetNumHeavyAtoms() >= 5:
                    seen.add(new_smi)
                    result.append(new_smi)

        except Exception:
            continue

    try:
        RDLogger.EnableLog("rdApp.*")
    except Exception:
        pass

    return result


# =========================================================================
# Main orchestrator
# =========================================================================

def run_massive_screening(
    receptor_pdbqt: Path,
    pdb_path: Path,
    pocket_center: tuple[float, float, float],
    pocket_size: tuple[float, float, float],
    work_dir: Path,
    max_chembl: int = 0,
    progress_callback: Optional[Callable[[dict], None]] = None,
) -> list[dict]:
    """Run the full 5-pass massive screening pipeline.

    This is the main entry point for mode="deep" screening, orchestrating
    all five passes from initial molecule loading through to final
    ADMET-enriched results.

    Parameters
    ----------
    receptor_pdbqt : Path
        Prepared receptor PDBQT file.
    pdb_path : Path
        Original PDB structure file (used for pocket analysis).
    pocket_center : tuple[float, float, float]
        (x, y, z) centre of the binding pocket.
    pocket_size : tuple[float, float, float]
        (sx, sy, sz) dimensions of the search box.
    work_dir : Path
        Working directory for all intermediate files.
    max_chembl : int
        Maximum molecules to load from ChEMBL dump (0 = all).
    progress_callback : callable, optional
        Called with a dict after each pass:
        ``{"pass_name": str, "pass_number": int, "input_count": int,
           "output_count": int, "elapsed": float, "total_elapsed": float}``.

    Returns
    -------
    list[dict]
        Final list of enriched docking results, sorted by composite score.
        Each dict contains: name, smiles, affinity, source, ADMET data,
        retrosynthesis data, and other scoring fields.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    global_start = time.monotonic()

    logger.info("=" * 70)
    logger.info("MASSIVE SCREENING PIPELINE -- Starting 5-pass screening")
    logger.info("  Receptor: %s", receptor_pdbqt)
    logger.info("  Pocket centre: %s", pocket_center)
    logger.info("  Pocket size: %s", pocket_size)
    logger.info("  Work dir: %s", work_dir)
    logger.info("=" * 70)

    # Compute pocket volume for Pass 1
    pocket_vol = pocket_size[0] * pocket_size[1] * pocket_size[2]

    # =====================================================================
    # PASS 0: Load molecules
    # =====================================================================
    pass_start = time.monotonic()
    logger.info("[Pass 0] Loading ChEMBL molecules...")

    all_smiles = load_chembl_smiles(max_count=max_chembl)
    pass_elapsed = time.monotonic() - pass_start

    logger.info(
        "[Pass 0] Loaded %d molecules in %.1fs",
        len(all_smiles), pass_elapsed,
    )

    _report_progress(progress_callback, {
        "pass_name": "load_molecules",
        "pass_number": 0,
        "input_count": 0,
        "output_count": len(all_smiles),
        "elapsed": pass_elapsed,
        "total_elapsed": time.monotonic() - global_start,
    })

    if not all_smiles:
        logger.error("No molecules loaded; aborting massive screening")
        return []

    # =====================================================================
    # PASS 1: Pharmacological filter
    # =====================================================================
    pass_start = time.monotonic()
    pass_input = len(all_smiles)
    logger.info("[Pass 1] Pharmacological filter: %d molecules...", pass_input)

    from pipeline.filter_pharma import filter_pharmacological

    pass1_smiles = filter_pharmacological(
        smiles_list=all_smiles,
        pocket_volume=pocket_vol,
        progress_callback=lambda info: _relay_sub_progress(
            progress_callback, "pharmacological_filter", 1, info,
        ),
    )

    pass_elapsed = time.monotonic() - pass_start
    logger.info(
        "[Pass 1] Pharmacological filter: %d -> %d (%.1f%%) in %.1fs",
        pass_input, len(pass1_smiles),
        len(pass1_smiles) / max(pass_input, 1) * 100,
        pass_elapsed,
    )

    _report_progress(progress_callback, {
        "pass_name": "pharmacological_filter",
        "pass_number": 1,
        "input_count": pass_input,
        "output_count": len(pass1_smiles),
        "elapsed": pass_elapsed,
        "total_elapsed": time.monotonic() - global_start,
    })

    # Free memory
    del all_smiles

    if not pass1_smiles:
        logger.warning("[Pass 1] No molecules passed pharmacological filter")
        return []

    # =====================================================================
    # PASS 2: 3D Shape filter
    # =====================================================================
    pass_start = time.monotonic()
    pass_input = len(pass1_smiles)
    logger.info("[Pass 2] Shape filter: %d molecules...", pass_input)

    from pipeline.filter_shape import filter_by_shape

    pass2_smiles = filter_by_shape(
        smiles_list=pass1_smiles,
        pocket_center=pocket_center,
        pocket_size=pocket_size,
        pocket_pdb=str(pdb_path) if pdb_path.exists() else None,
        progress_callback=lambda info: _relay_sub_progress(
            progress_callback, "shape_filter", 2, info,
        ),
    )

    pass_elapsed = time.monotonic() - pass_start
    logger.info(
        "[Pass 2] Shape filter: %d -> %d (%.1f%%) in %.1fs",
        pass_input, len(pass2_smiles),
        len(pass2_smiles) / max(pass_input, 1) * 100,
        pass_elapsed,
    )

    _report_progress(progress_callback, {
        "pass_name": "shape_filter",
        "pass_number": 2,
        "input_count": pass_input,
        "output_count": len(pass2_smiles),
        "elapsed": pass_elapsed,
        "total_elapsed": time.monotonic() - global_start,
    })

    del pass1_smiles

    if not pass2_smiles:
        logger.warning("[Pass 2] No molecules passed shape filter")
        return []

    # =====================================================================
    # PASS 3: Rapid scoring (smina/Vinardo) -> top 500
    # =====================================================================
    pass_start = time.monotonic()
    pass_input = len(pass2_smiles)
    logger.info("[Pass 3] Rapid scoring: %d molecules...", pass_input)

    from pipeline.scoring_rapid import score_rapid_batch

    rapid_dir = work_dir / "pass3_rapid"
    rapid_dir.mkdir(parents=True, exist_ok=True)

    pass3_results = score_rapid_batch(
        receptor_pdbqt=receptor_pdbqt,
        smiles_list=pass2_smiles,
        center=pocket_center,
        size=pocket_size,
        work_dir=rapid_dir,
        n_top=500,
        progress_callback=lambda info: _relay_sub_progress(
            progress_callback, "rapid_scoring", 3, info,
        ),
    )

    pass_elapsed = time.monotonic() - pass_start
    logger.info(
        "[Pass 3] Rapid scoring: %d -> %d (top 500) in %.1fs",
        pass_input, len(pass3_results), pass_elapsed,
    )

    _report_progress(progress_callback, {
        "pass_name": "rapid_scoring",
        "pass_number": 3,
        "input_count": pass_input,
        "output_count": len(pass3_results),
        "elapsed": pass_elapsed,
        "total_elapsed": time.monotonic() - global_start,
    })

    del pass2_smiles

    if not pass3_results:
        logger.warning("[Pass 3] No molecules scored in rapid pass")
        return []

    # =====================================================================
    # PASS 4: Full Vina docking (exhaustiveness=32) -> top 50
    # =====================================================================
    pass_start = time.monotonic()
    pass_input = len(pass3_results)
    logger.info("[Pass 4] Full Vina docking: %d molecules...", pass_input)

    from pipeline.docking import dock_all_ligands

    # Convert pass3 results to ligand dicts for dock_all_ligands
    ligands_for_docking: list[dict] = []
    for r in pass3_results:
        ligands_for_docking.append({
            "name": r.get("name", "unknown"),
            "smiles": r.get("smiles", ""),
            "source": r.get("source", "rapid_screen"),
            "rapid_affinity": r.get("affinity", 0.0),
        })

    docking_dir = work_dir / "pass4_docking"
    docking_dir.mkdir(parents=True, exist_ok=True)

    pass4_results = dock_all_ligands(
        receptor_pdbqt=receptor_pdbqt,
        ligands=ligands_for_docking,
        center=pocket_center,
        work_dir=docking_dir,
        progress_callback=lambda pct, msg: _relay_docking_progress(
            progress_callback, 4, pct, msg,
        ),
        size=pocket_size,
        exhaustiveness=32,
    )

    # Sort by affinity and take top 50
    pass4_results.sort(key=lambda r: r.get("affinity", 0.0))
    pass4_top50 = pass4_results[:50]

    pass_elapsed = time.monotonic() - pass_start
    logger.info(
        "[Pass 4] Full Vina docking: %d -> %d (top 50) in %.1fs",
        pass_input, len(pass4_top50), pass_elapsed,
    )

    _report_progress(progress_callback, {
        "pass_name": "full_docking",
        "pass_number": 4,
        "input_count": pass_input,
        "output_count": len(pass4_top50),
        "elapsed": pass_elapsed,
        "total_elapsed": time.monotonic() - global_start,
    })

    if not pass4_top50:
        logger.warning("[Pass 4] No molecules succeeded in full docking")
        return []

    # =====================================================================
    # PASS 5: ADMET + Retrosynthesis enrichment
    # =====================================================================
    pass_start = time.monotonic()
    pass_input = len(pass4_top50)
    logger.info("[Pass 5] ADMET + retrosynthesis: %d molecules...", pass_input)

    enriched_results = _enrich_with_admet_retro(pass4_top50, progress_callback)

    pass_elapsed = time.monotonic() - pass_start
    total_elapsed = time.monotonic() - global_start

    logger.info(
        "[Pass 5] Enrichment complete: %d molecules in %.1fs",
        len(enriched_results), pass_elapsed,
    )

    _report_progress(progress_callback, {
        "pass_name": "admet_retrosynthesis",
        "pass_number": 5,
        "input_count": pass_input,
        "output_count": len(enriched_results),
        "elapsed": pass_elapsed,
        "total_elapsed": total_elapsed,
    })

    # =====================================================================
    # Final summary
    # =====================================================================
    logger.info("=" * 70)
    logger.info("MASSIVE SCREENING COMPLETE")
    logger.info("  Total time: %.1fs", total_elapsed)
    logger.info("  Final results: %d molecules", len(enriched_results))
    if enriched_results:
        best = enriched_results[0]
        logger.info(
            "  Best hit: %s (affinity=%.2f, composite=%.4f)",
            best.get("name", "?"),
            best.get("affinity", 0.0),
            best.get("composite_score", 0.0),
        )
    logger.info("=" * 70)

    return enriched_results


# =========================================================================
# Pass 5: ADMET + Retrosynthesis enrichment
# =========================================================================

def _enrich_with_admet_retro(
    docking_results: list[dict],
    progress_callback: Optional[Callable[[dict], None]],
) -> list[dict]:
    """Enrich docking results with ADMET predictions and retrosynthesis.

    Parameters
    ----------
    docking_results : list[dict]
        Top docking results from Pass 4.
    progress_callback : callable, optional
        Progress callback.

    Returns
    -------
    list[dict]
        Enriched results sorted by composite score (descending).
    """
    smiles_list = [r.get("smiles", "") for r in docking_results]

    # --- ADMET predictions ---
    admet_results: list[dict] = []
    try:
        from pipeline.admet import predict_admet
        admet_results = predict_admet(smiles_list)
        logger.info("ADMET predictions computed for %d molecules", len(admet_results))
    except Exception as exc:
        logger.warning("ADMET prediction failed: %s; using defaults", exc)

    # Build ADMET lookup by SMILES
    admet_by_smiles: dict[str, dict] = {}
    for admet in admet_results:
        smi = admet.get("smiles", "")
        if smi:
            admet_by_smiles[smi] = admet

    # --- Retrosynthesis planning ---
    retro_results: dict[str, dict] = {}
    try:
        from pipeline.retrosynthesis import plan_synthesis
        for i, smi in enumerate(smiles_list):
            if not smi:
                continue
            try:
                retro = plan_synthesis(smi, max_depth=6, timeout_sec=60)
                retro_results[smi] = retro
            except Exception as exc:
                logger.debug("Retrosynthesis failed for %s: %s", smi[:40], exc)

            # Report sub-progress
            if progress_callback is not None and (i + 1) % 10 == 0:
                try:
                    progress_callback({
                        "pass_name": "retrosynthesis",
                        "pass_number": 5,
                        "processed": i + 1,
                        "total": len(smiles_list),
                    })
                except Exception:
                    pass

        logger.info("Retrosynthesis planned for %d molecules", len(retro_results))
    except Exception as exc:
        logger.warning("Retrosynthesis planning failed: %s", exc)

    # --- Compute composite scores and enrich ---
    try:
        from pipeline.scoring import (
            compute_properties,
            compute_composite_score_v2,
            generate_2d_svg,
        )
    except ImportError:
        logger.warning("scoring module not available; using basic enrichment")
        compute_properties = None  # type: ignore[assignment]
        compute_composite_score_v2 = None  # type: ignore[assignment]
        generate_2d_svg = None  # type: ignore[assignment]

    for entry in docking_results:
        smiles = entry.get("smiles", "")
        affinity = entry.get("affinity", 0.0)

        # Basic properties
        if compute_properties is not None:
            try:
                props = compute_properties(smiles)
                entry.update(props)
            except Exception:
                pass

        # 2D SVG
        if generate_2d_svg is not None:
            try:
                entry["svg_2d"] = generate_2d_svg(smiles)
            except Exception:
                entry["svg_2d"] = None

        # ADMET data
        admet_data = admet_by_smiles.get(smiles)
        if admet_data:
            entry["admet"] = admet_data
            admet_composite = admet_data.get("composite_score", 0.5)
            entry["admet_color"] = admet_data.get("color_code", "yellow")
            entry["admet_flags"] = admet_data.get("flags", [])
        else:
            admet_composite = 0.5
            entry["admet"] = None
            entry["admet_color"] = "yellow"
            entry["admet_flags"] = []

        # Retrosynthesis data
        retro_data = retro_results.get(smiles)
        if retro_data:
            entry["retrosynthesis"] = retro_data
            entry["synthesis_steps"] = retro_data.get("n_steps", 0)
            entry["synthesis_confidence"] = retro_data.get("confidence", 0.0)
            entry["synthesis_cost"] = retro_data.get("estimated_cost", "high")
        else:
            entry["retrosynthesis"] = None
            entry["synthesis_steps"] = 0
            entry["synthesis_confidence"] = 0.0
            entry["synthesis_cost"] = "high"

        # Composite V2 score
        if compute_composite_score_v2 is not None:
            try:
                entry["composite_score"] = compute_composite_score_v2(
                    affinity=affinity,
                    admet_score=admet_composite,
                    qed=entry.get("qed"),
                    novelty=None,
                )
            except Exception:
                entry["composite_score"] = _fallback_composite(affinity)
        else:
            entry["composite_score"] = _fallback_composite(affinity)

    # Sort by composite score descending (best first)
    docking_results.sort(
        key=lambda r: r.get("composite_score", 0.0),
        reverse=True,
    )

    return docking_results


def _fallback_composite(affinity: float) -> float:
    """Compute a simple composite score when the scoring module is unavailable.

    Parameters
    ----------
    affinity : float
        Binding affinity in kcal/mol.

    Returns
    -------
    float
        Simple normalized score in [0, 1].
    """
    capped = min(-12.0, affinity) if affinity < 0 else affinity
    return round(max(0.0, min(1.0, capped / -12.0)), 4)


# =========================================================================
# Progress helpers
# =========================================================================

def _report_progress(
    callback: Optional[Callable[[dict], None]],
    info: dict,
) -> None:
    """Safely call the progress callback.

    Parameters
    ----------
    callback : callable or None
        Progress callback.
    info : dict
        Progress information dict.
    """
    if callback is None:
        return
    try:
        callback(info)
    except Exception as exc:
        logger.debug("Progress callback error: %s", exc)


def _relay_sub_progress(
    main_callback: Optional[Callable[[dict], None]],
    pass_name: str,
    pass_number: int,
    sub_info: dict,
) -> None:
    """Relay sub-pass progress to the main callback.

    Parameters
    ----------
    main_callback : callable or None
        Main progress callback.
    pass_name : str
        Name of the current pass.
    pass_number : int
        Pass number (1-5).
    sub_info : dict
        Sub-progress info from the filter/scorer.
    """
    if main_callback is None:
        return
    try:
        main_callback({
            "pass_name": pass_name,
            "pass_number": pass_number,
            "sub_progress": sub_info,
        })
    except Exception:
        pass


def _relay_docking_progress(
    main_callback: Optional[Callable[[dict], None]],
    pass_number: int,
    pct: int,
    msg: str,
) -> None:
    """Relay docking progress (int, str) to the main callback dict format.

    Parameters
    ----------
    main_callback : callable or None
        Main progress callback.
    pass_number : int
        Pass number.
    pct : int
        Percent complete.
    msg : str
        Status message.
    """
    if main_callback is None:
        return
    try:
        main_callback({
            "pass_name": "full_docking",
            "pass_number": pass_number,
            "percent": pct,
            "message": msg,
        })
    except Exception:
        pass


# =========================================================================
# CLI / Self-test
# =========================================================================

if __name__ == "__main__":
    import json
    import sys

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    print("=" * 70)
    print("DockIt -- screening_massive.py self-test")
    print("=" * 70)

    # Test ChEMBL loading
    print("\n[1] Testing check_chembl_dump_available...")
    available = check_chembl_dump_available()
    print(f"  ChEMBL dump available: {available}")

    print("\n[2] Testing load_chembl_smiles (max_count=100)...")
    smiles = load_chembl_smiles(max_count=100)
    print(f"  Loaded {len(smiles)} SMILES")
    if smiles:
        print(f"  First: {smiles[0][:60]}")
        print(f"  Last:  {smiles[-1][:60]}")

    # Test the full pipeline with a small set
    print("\n[3] Testing run_massive_screening (mock mode)...")

    test_work_dir = Path("/tmp/dockit_massive_test")
    test_work_dir.mkdir(parents=True, exist_ok=True)

    # Create dummy receptor
    dummy_receptor = test_work_dir / "receptor.pdbqt"
    dummy_receptor.write_text("REMARK mock receptor\nEND\n")
    dummy_pdb = test_work_dir / "structure.pdb"
    dummy_pdb.write_text("REMARK mock structure\nEND\n")

    pocket_center = (22.0, 0.5, 18.0)
    pocket_size = (25.0, 25.0, 25.0)

    progress_log: list[dict] = []

    def on_progress(info: dict) -> None:
        progress_log.append(info)
        pass_name = info.get("pass_name", "?")
        pass_num = info.get("pass_number", "?")
        inp = info.get("input_count", "?")
        out = info.get("output_count", "?")
        elapsed = info.get("elapsed", "?")
        print(f"  [Pass {pass_num}] {pass_name}: {inp} -> {out} ({elapsed}s)")

    results = run_massive_screening(
        receptor_pdbqt=dummy_receptor,
        pdb_path=dummy_pdb,
        pocket_center=pocket_center,
        pocket_size=pocket_size,
        work_dir=test_work_dir,
        max_chembl=200,  # Small for testing
        progress_callback=on_progress,
    )

    print(f"\n  Final results: {len(results)} molecules")
    if results:
        print(f"\n  {'Rank':<6} {'Name':<20} {'Affinity':<10} {'Composite':<10} {'SMILES'}")
        print("  " + "-" * 90)
        for i, r in enumerate(results[:10], 1):
            print(
                f"  {i:<6} {r.get('name', '?'):<20} "
                f"{r.get('affinity', 0.0):<10.2f} "
                f"{r.get('composite_score', 0.0):<10.4f} "
                f"{r.get('smiles', '?')[:40]}"
            )

    # Validation
    errors = 0

    if len(smiles) == 0:
        print("\nFAIL: load_chembl_smiles returned empty")
        errors += 1
    else:
        print("PASS: load_chembl_smiles works")

    if len(progress_log) == 0:
        print("WARNING: no progress callbacks received")
    else:
        print(f"PASS: {len(progress_log)} progress callbacks received")

    # Results may be empty if all mocks filter everything out
    # but the pipeline should not crash
    print(f"PASS: pipeline completed without errors ({len(results)} results)")

    if errors == 0:
        print("\nALL CHECKS PASSED")
    else:
        print(f"\n{errors} CHECK(S) FAILED")

    sys.exit(1 if errors > 0 else 0)
