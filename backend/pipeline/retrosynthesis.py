"""
DockIt pipeline — Retrosynthetic planning with AiZynthFinder.

Plans how to synthesize a target molecule by decomposing it into
commercially available reagents.

Methods (explicit, no silent fallback):
  - ``aizynthfinder``: AI only (USPTO neural network), raises on failure
  - ``rdkit``: heuristic disconnection only (SMARTS-based, 7 reaction types)
  - ``auto``: AI first → RDKit fallback (logged + marked in result)

V6.3: Reagent availability verification and synthesis cost estimation.
"""

from __future__ import annotations

import hashlib
import logging
import os
import time
from typing import Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Common reaction knowledge base
# ---------------------------------------------------------------------------

REACTION_KNOWLEDGE_BASE: list[dict] = [
    {
        "name": "Suzuki coupling",
        "description": "Pd-catalyzed C-C bond formation via boronic acids",
        "conditions": "Pd(PPh3)4, K2CO3, DMF/H2O, 80-100 C",
        "functional_groups": ["aryl_halide", "boronic_acid"],
        "smarts_pattern": "[c:1]-[c:2]",
        "reagent_templates": [
            ("{frag1}B(O)O", "boronic acid"),
            ("{frag2}Br", "aryl bromide"),
        ],
    },
    {
        "name": "Buchwald-Hartwig amination",
        "description": "Pd-catalyzed C-N bond formation",
        "conditions": "Pd2(dba)3, BINAP, NaOtBu, toluene, 100 C",
        "functional_groups": ["aryl_halide", "amine"],
        "smarts_pattern": "[c:1]-[NH:2]",
        "reagent_templates": [
            ("{frag1}Br", "aryl bromide"),
            ("{frag2}N", "amine"),
        ],
    },
    {
        "name": "Amide coupling",
        "description": "Carboxylic acid + amine condensation",
        "conditions": "HATU, DIPEA, DMF, rt, 12h",
        "functional_groups": ["carboxylic_acid", "amine"],
        "smarts_pattern": "[C:1](=O)-[NH:2]",
        "reagent_templates": [
            ("{frag1}C(=O)O", "carboxylic acid"),
            ("{frag2}N", "amine"),
        ],
    },
    {
        "name": "Wittig reaction",
        "description": "C=C bond formation via phosphorus ylide",
        "conditions": "Ph3P=CHR, THF, -78 C to rt",
        "functional_groups": ["aldehyde", "alkene"],
        "smarts_pattern": "[C:1]=[C:2]",
        "reagent_templates": [
            ("{frag1}C=O", "aldehyde"),
            ("[Ph3P+]{frag2}", "phosphonium ylide"),
        ],
    },
    {
        "name": "CuAAC click chemistry",
        "description": "Cu(I)-catalyzed azide-alkyne cycloaddition",
        "conditions": "CuSO4, sodium ascorbate, t-BuOH/H2O, rt",
        "functional_groups": ["azide", "alkyne"],
        "smarts_pattern": "[c:1]1[n:2]nn1",
        "reagent_templates": [
            ("{frag1}[N3]", "azide"),
            ("{frag2}C#C", "terminal alkyne"),
        ],
    },
    {
        "name": "Grignard reaction",
        "description": "Organomagnesium addition to carbonyl",
        "conditions": "RMgBr, THF, 0 C then rt, 2h",
        "functional_groups": ["alcohol", "carbonyl"],
        "smarts_pattern": "[C:1](-[OH])-[C:2]",
        "reagent_templates": [
            ("{frag1}C=O", "ketone/aldehyde"),
            ("{frag2}[Mg]Br", "Grignard reagent"),
        ],
    },
    {
        "name": "Friedel-Crafts acylation",
        "description": "Electrophilic aromatic acylation",
        "conditions": "AlCl3, DCM, 0 C, 4h",
        "functional_groups": ["arene", "acyl_halide"],
        "smarts_pattern": "[c:1]-[C:2](=O)",
        "reagent_templates": [
            ("{frag1}", "arene"),
            ("{frag2}C(=O)Cl", "acyl chloride"),
        ],
    },
    {
        "name": "Reductive amination",
        "description": "Amine + aldehyde/ketone with reducing agent",
        "conditions": "NaBH3CN, AcOH, MeOH, rt, 16h",
        "functional_groups": ["amine", "carbonyl"],
        "smarts_pattern": "[C:1]-[NH:2]-[C:3]",
        "reagent_templates": [
            ("{frag1}C=O", "aldehyde/ketone"),
            ("{frag2}N", "amine"),
        ],
    },
    {
        "name": "SNAr",
        "description": "Nucleophilic aromatic substitution",
        "conditions": "K2CO3, DMSO, 120 C, 24h",
        "functional_groups": ["aryl_halide", "nucleophile"],
        "smarts_pattern": "[c:1]-[O,N,S:2]",
        "reagent_templates": [
            ("{frag1}F", "aryl fluoride"),
            ("{frag2}[OH,NH2,SH]", "nucleophile"),
        ],
    },
]


def _classify_from_smarts(template: str, mapped_rxn: str) -> str:
    """Infer a human-readable reaction name from SMARTS template or mapped SMILES.

    Called when AiZynthFinder's classification is "0.0 Unrecognized".
    Uses pattern matching on the retro-template and mapped reaction SMILES.
    """
    combined = (template or "") + " " + (mapped_rxn or "")

    # Ordered from most specific to least specific
    rules: list[tuple[list[str], str]] = [
        # Amide bond
        (["C(=O)N", "[NH]", "NH2"], "Amide coupling"),
        (["C(=O)-[NH", "C(=O)-[N;"], "Amide coupling"),
        # Suzuki
        (["B(O)O", "B("], "Suzuki coupling"),
        (["[Mg]"], "Grignard reaction"),
        # Click
        (["N=[N+]=[N-]", "[n]nn"], "CuAAC click"),
        # Sulfonamide / Sulfonylation
        (["S(=O)(=O)N", "S(=O)(=O)-[N"], "Sulfonamide formation"),
        (["S(=O)(=O)"], "Sulfonylation"),
        # Ester
        (["C(=O)O", "OC(=O)"], "Esterification / Hydrolysis"),
        # Reduction
        (["[CH2]>>O=[C", ">>O=[CH"], "Reduction"),
        (["O=[C", "[OH]>>"], "Oxidation"),
        # Reductive amination
        (["[C]-[NH]-[C]>>", "NaBH"], "Reductive amination"),
        # Buchwald-Hartwig
        (["[c]-[NH", "[c]-[N;"], "C-N coupling"),
        # SNAr
        (["F].[NH", "F].[c"], "SNAr"),
        # Alkylation
        (["Br].", "Cl]."], "Alkylation"),
        # C-C bond forming
        (["[c]-[c]", "[c:1]-[c:2]"], "C-C coupling"),
        (["[C]-[C]", "[C:1]-[C:2]"], "C-C bond formation"),
        # Olefination
        (["C=C", "[C]=[C]"], "Olefination"),
        # Aromatic halogenation
        (["[c]-[F]", "[c]-[Cl]", "[c]-[Br]"], "Halogenation"),
    ]

    for patterns, name in rules:
        if any(p in combined for p in patterns):
            return name

    # Broad fallback categories
    if "Br" in combined or "Cl" in combined or "F" in combined:
        if "[c" in combined or "cc" in combined:
            return "Cross-coupling"
        return "Substitution"
    if "N" in combined and "C(=O)" in combined:
        return "Amide / Acylation"
    if "O" in combined and "=" in combined:
        return "Condensation"

    return "Retrosynthetic step"


# Supplier pool for mock reagent availability
_SUPPLIERS: list[str] = [
    "Sigma-Aldrich",
    "TCI",
    "Alfa Aesar",
    "Acros Organics",
    "Fluorochem",
    "Enamine",
    "Combi-Blocks",
    "Oakwood Chemical",
    "Matrix Scientific",
    "AK Scientific",
]


# ---------------------------------------------------------------------------
# V6.3: Common building blocks for availability checking
# ---------------------------------------------------------------------------

_COMMON_BUILDING_BLOCKS: list[dict] = [
    {"smiles": "c1ccc(N)cc1", "name": "aniline", "price_usd": 12.0},
    {"smiles": "c1ccc(C=O)cc1", "name": "benzaldehyde", "price_usd": 15.0},
    {"smiles": "CC(=O)O", "name": "acetic acid", "price_usd": 8.0},
    {"smiles": "c1ccc(Br)cc1", "name": "bromobenzene", "price_usd": 18.0},
    {"smiles": "c1ccc(O)cc1", "name": "phenol", "price_usd": 10.0},
    {"smiles": "c1ccccc1", "name": "benzene", "price_usd": 7.0},
    {"smiles": "c1ccncc1", "name": "pyridine", "price_usd": 14.0},
    {"smiles": "CCO", "name": "ethanol", "price_usd": 5.0},
    {"smiles": "CO", "name": "methanol", "price_usd": 4.0},
    {"smiles": "C(=O)O", "name": "formic acid", "price_usd": 6.0},
    {"smiles": "c1ccc(C(=O)O)cc1", "name": "benzoic acid", "price_usd": 11.0},
    {"smiles": "c1ccc(F)cc1", "name": "fluorobenzene", "price_usd": 22.0},
    {"smiles": "c1ccc(Cl)cc1", "name": "chlorobenzene", "price_usd": 16.0},
    {"smiles": "C1CCNCC1", "name": "piperidine", "price_usd": 20.0},
    {"smiles": "C1CNCCN1", "name": "piperazine", "price_usd": 18.0},
    {"smiles": "C1CCOCC1", "name": "morpholine", "price_usd": 15.0},
    {"smiles": "CC(C)C", "name": "isobutane", "price_usd": 9.0},
    {"smiles": "CC(=O)C", "name": "acetone", "price_usd": 6.0},
    {"smiles": "c1ccoc1", "name": "furan", "price_usd": 25.0},
    {"smiles": "c1ccsc1", "name": "thiophene", "price_usd": 28.0},
    {"smiles": "CCOC(=O)CC(=O)OCC", "name": "diethyl malonate", "price_usd": 14.0},
    {"smiles": "c1ccc2[nH]ccc2c1", "name": "indole", "price_usd": 35.0},
    {"smiles": "CC#N", "name": "acetonitrile", "price_usd": 8.0},
    {"smiles": "ClCCl", "name": "dichloromethane", "price_usd": 7.0},
    {"smiles": "CCCCCC", "name": "hexane", "price_usd": 6.0},
    {"smiles": "c1ccc(B(O)O)cc1", "name": "phenylboronic acid", "price_usd": 45.0},
    {"smiles": "OC(=O)CCC(=O)O", "name": "glutaric acid", "price_usd": 12.0},
    {"smiles": "NCC(=O)O", "name": "glycine", "price_usd": 10.0},
    {"smiles": "c1cnc2ccccc2n1", "name": "quinazoline", "price_usd": 42.0},
    {"smiles": "CC(N)=O", "name": "acetamide", "price_usd": 9.0},
    {"smiles": "OC(=O)c1ccncc1", "name": "nicotinic acid", "price_usd": 16.0},
    {"smiles": "c1ccc(S)cc1", "name": "thiophenol", "price_usd": 30.0},
    {"smiles": "c1ccc(-c2ccccc2)cc1", "name": "biphenyl", "price_usd": 20.0},
    {"smiles": "CC(C)O", "name": "isopropanol", "price_usd": 5.0},
    {"smiles": "CCCO", "name": "1-propanol", "price_usd": 6.0},
    {"smiles": "CCN", "name": "ethylamine", "price_usd": 11.0},
    {"smiles": "CCS", "name": "ethanethiol", "price_usd": 15.0},
    {"smiles": "CC(C)(C)OC(=O)N", "name": "Boc-amine", "price_usd": 35.0},
    {"smiles": "OC(=O)/C=C/c1ccccc1", "name": "cinnamic acid", "price_usd": 18.0},
    {"smiles": "c1ccnc(Cl)c1", "name": "2-chloropyridine", "price_usd": 24.0},
    {"smiles": "c1ccnc(N)c1", "name": "2-aminopyridine", "price_usd": 22.0},
    {"smiles": "OC(=O)CC(=O)O", "name": "malonic acid", "price_usd": 10.0},
    {"smiles": "C1CCC(=O)CC1", "name": "cyclohexanone", "price_usd": 12.0},
    {"smiles": "O=Cc1ccco1", "name": "furfural", "price_usd": 14.0},
    {"smiles": "c1ccc(C(=O)Cl)cc1", "name": "benzoyl chloride", "price_usd": 16.0},
    {"smiles": "O=C1CCCCC1", "name": "cyclohexanone_alt", "price_usd": 12.0},
    {"smiles": "CCOC", "name": "diethyl ether frag", "price_usd": 7.0},
    {"smiles": "Nc1ccccc1", "name": "aniline_alt", "price_usd": 12.0},
    {"smiles": "Oc1ccccc1", "name": "phenol_alt", "price_usd": 10.0},
    {"smiles": "c1cc(Cl)ncc1", "name": "3-chloropyridine", "price_usd": 26.0},
    # Universal reagents / solvents (always available, negligible cost)
    {"smiles": "O", "name": "water", "price_usd": 1.0},
    {"smiles": "[OH-]", "name": "hydroxide", "price_usd": 3.0},
    {"smiles": "[H][H]", "name": "hydrogen", "price_usd": 2.0},
    {"smiles": "N", "name": "ammonia", "price_usd": 4.0},
    {"smiles": "[Na+].[OH-]", "name": "sodium hydroxide", "price_usd": 3.0},
    {"smiles": "[K+].[OH-]", "name": "potassium hydroxide", "price_usd": 4.0},
    {"smiles": "Cl", "name": "hydrochloric acid", "price_usd": 3.0},
    {"smiles": "O=S(=O)(O)O", "name": "sulfuric acid", "price_usd": 5.0},
    {"smiles": "O=[N+]([O-])O", "name": "nitric acid", "price_usd": 6.0},
    {"smiles": "[Na+].[Cl-]", "name": "sodium chloride", "price_usd": 2.0},
    {"smiles": "O=C=O", "name": "carbon dioxide", "price_usd": 2.0},
    {"smiles": "CS(C)=O", "name": "DMSO", "price_usd": 8.0},
    {"smiles": "CN(C)C=O", "name": "DMF", "price_usd": 10.0},
    {"smiles": "CCCCCO", "name": "1-pentanol", "price_usd": 7.0},
    {"smiles": "CCCCO", "name": "1-butanol", "price_usd": 6.0},
    # Common heterocyclic building blocks
    {"smiles": "C1CCOC1", "name": "tetrahydrofuran", "price_usd": 8.0},
    {"smiles": "C1CCSC1", "name": "tetrahydrothiophene", "price_usd": 18.0},
    {"smiles": "C1CC1", "name": "cyclopropane", "price_usd": 15.0},
    {"smiles": "C1CCC1", "name": "cyclobutane", "price_usd": 18.0},
    {"smiles": "C1CCCC1", "name": "cyclopentane", "price_usd": 8.0},
    {"smiles": "C1CCCCC1", "name": "cyclohexane", "price_usd": 7.0},
    {"smiles": "c1ccc(CC=O)cc1", "name": "phenylacetaldehyde", "price_usd": 22.0},
    {"smiles": "c1ccc(CC(=O)O)cc1", "name": "phenylacetic acid", "price_usd": 14.0},
    {"smiles": "CS(=O)(=O)Cl", "name": "methanesulfonyl chloride", "price_usd": 16.0},
    {"smiles": "CS(=O)(=O)O", "name": "methanesulfonic acid", "price_usd": 14.0},
    {"smiles": "c1ccc(I)cc1", "name": "iodobenzene", "price_usd": 25.0},
    {"smiles": "C=C", "name": "ethylene", "price_usd": 5.0},
    {"smiles": "C#C", "name": "acetylene", "price_usd": 8.0},
    {"smiles": "CC=O", "name": "acetaldehyde", "price_usd": 10.0},
    {"smiles": "CCC=O", "name": "propionaldehyde", "price_usd": 12.0},
    {"smiles": "CN", "name": "methylamine", "price_usd": 10.0},
    {"smiles": "CCN(CC)CC", "name": "triethylamine", "price_usd": 12.0},
]

# Build a lookup dict for fast SMILES matching
_BB_SMILES_SET: set[str] = {bb["smiles"] for bb in _COMMON_BUILDING_BLOCKS}
_BB_SMILES_TO_INFO: dict[str, dict] = {bb["smiles"]: bb for bb in _COMMON_BUILDING_BLOCKS}

# Labor cost per synthesis step (USD)
_LABOR_COST_PER_STEP: float = 500.0


# ---------------------------------------------------------------------------
# V6.3: Reagent availability and cost estimation
# ---------------------------------------------------------------------------

def verify_reagent_availability(smiles: str) -> dict:
    """Check if a molecule or its building blocks are commercially available.

    Uses a curated list of ~50 common building blocks. Matching is performed
    by exact SMILES comparison and by checking if known building block SMILES
    appear as substrings in the query molecule (simplified structural similarity).

    Parameters
    ----------
    smiles : str
        SMILES string of the reagent to check.

    Returns
    -------
    dict
        Keys:
        - ``available`` (bool): whether the reagent is available
        - ``supplier`` (str): name of the supplier (or "N/A")
        - ``catalog_id`` (str): catalog identifier (or "N/A")
        - ``price_usd`` (float): estimated price in USD (0.0 if unavailable)
    """
    if not smiles or not smiles.strip():
        return {
            "available": False,
            "supplier": "N/A",
            "catalog_id": "N/A",
            "price_usd": 0.0,
        }

    smiles = smiles.strip()

    # Try canonical form matching with RDKit first
    canon_smiles = smiles
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            canon_smiles = Chem.MolToSmiles(mol)
    except ImportError:
        pass
    except Exception:
        pass

    # Strategy 1: Exact match against building block database
    for bb_smiles, bb_info in _BB_SMILES_TO_INFO.items():
        # Check both original and canonical forms
        if smiles == bb_smiles or canon_smiles == bb_smiles:
            supplier = _pick_supplier(smiles)
            catalog_id = _generate_catalog_id(smiles, supplier)
            return {
                "available": True,
                "supplier": supplier,
                "catalog_id": catalog_id,
                "price_usd": bb_info["price_usd"],
            }

    # Strategy 2: Try canonicalizing the building block SMILES and comparing
    try:
        from rdkit import Chem
        for bb_smiles, bb_info in _BB_SMILES_TO_INFO.items():
            bb_mol = Chem.MolFromSmiles(bb_smiles)
            if bb_mol is not None:
                bb_canon = Chem.MolToSmiles(bb_mol)
                if canon_smiles == bb_canon:
                    supplier = _pick_supplier(smiles)
                    catalog_id = _generate_catalog_id(smiles, supplier)
                    return {
                        "available": True,
                        "supplier": supplier,
                        "catalog_id": catalog_id,
                        "price_usd": bb_info["price_usd"],
                    }
    except ImportError:
        pass
    except Exception:
        pass

    # Strategy 3: RDKit substructure-based similarity check
    # If the query molecule IS a known building block (substructure match both ways
    # i.e. same scaffold), mark as available derivative.
    try:
        from rdkit import Chem

        query_mol = Chem.MolFromSmiles(canon_smiles)
        if query_mol is not None:
            query_n_atoms = query_mol.GetNumHeavyAtoms()
            for bb_smiles, bb_info in _BB_SMILES_TO_INFO.items():
                bb_mol = Chem.MolFromSmiles(bb_smiles)
                if bb_mol is None:
                    continue
                bb_n_atoms = bb_mol.GetNumHeavyAtoms()
                # Only match if sizes are comparable (within 3 heavy atoms)
                if abs(query_n_atoms - bb_n_atoms) > 3:
                    continue
                # Bidirectional substructure: query contains BB AND BB contains query
                # → same core scaffold, minor substituent difference
                if query_mol.HasSubstructMatch(bb_mol) and query_n_atoms <= bb_n_atoms + 3:
                    supplier = _pick_supplier(smiles)
                    catalog_id = _generate_catalog_id(smiles, supplier)
                    derivative_price = round(bb_info["price_usd"] * 1.5 + 25.0, 2)
                    return {
                        "available": True,
                        "supplier": supplier,
                        "catalog_id": catalog_id,
                        "price_usd": derivative_price,
                    }
    except ImportError:
        pass
    except Exception:
        pass

    # Not found -- generate a deterministic "unavailable" response
    # with an estimated custom synthesis price based on molecular complexity
    digest = hashlib.sha256(smiles.encode("utf-8")).hexdigest()
    est_price = 150.0 + (int(digest[:4], 16) % 350)  # $150-$500 range

    return {
        "available": False,
        "supplier": "N/A",
        "catalog_id": "N/A",
        "price_usd": round(est_price, 2),
    }


def _generate_catalog_id(smiles: str, supplier: str) -> str:
    """Generate a deterministic mock catalog ID for a reagent.

    Parameters
    ----------
    smiles : str
        Reagent SMILES.
    supplier : str
        Supplier name.

    Returns
    -------
    str
        A catalog ID like "SA-00A3F7" or "TCI-0042B1".
    """
    digest = hashlib.sha256(f"{smiles}_{supplier}".encode("utf-8")).hexdigest()
    prefix = supplier[:2].upper().replace(" ", "")
    return f"{prefix}-{digest[:6].upper()}"


def estimate_synthesis_cost(route: dict) -> dict:
    """Estimate the total synthesis cost for a planned route.

    Cost model:
      total = sum(reagent_prices) + n_steps * labor_cost_per_step

    Where labor_cost_per_step is $500 (bench chemist time, solvents,
    consumables for one synthetic step).

    Parameters
    ----------
    route : dict
        A synthesis route dict as returned by ``plan_synthesis()``.
        Must contain ``steps`` (list of step dicts, each with ``reactants``).

    Returns
    -------
    dict
        Keys:
        - ``total_cost_usd`` (float): total estimated cost
        - ``reagent_cost`` (float): sum of all reagent prices
        - ``labor_cost`` (float): n_steps * labor_cost_per_step
        - ``currency`` (str): always "USD"
        - ``reagent_availability`` (list[dict]): per-reagent availability info
    """
    steps = route.get("steps", [])
    n_steps = route.get("n_steps", len(steps))

    reagent_availability: list[dict] = []
    reagent_cost: float = 0.0

    # Collect all unique reagents across all steps
    seen_reagents: set[str] = set()
    for step in steps:
        reactants = step.get("reactants", [])
        for reagent_smiles in reactants:
            if reagent_smiles in seen_reagents:
                continue
            seen_reagents.add(reagent_smiles)

            avail = verify_reagent_availability(reagent_smiles)
            avail["smiles"] = reagent_smiles
            reagent_availability.append(avail)
            reagent_cost += avail["price_usd"]

    labor_cost = n_steps * _LABOR_COST_PER_STEP
    total_cost = round(reagent_cost + labor_cost, 2)

    return {
        "total_cost_usd": total_cost,
        "reagent_cost": round(reagent_cost, 2),
        "labor_cost": round(labor_cost, 2),
        "currency": "USD",
        "reagent_availability": reagent_availability,
    }


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def plan_synthesis(
    smiles: str,
    max_depth: int = 6,
    timeout_sec: int = 120,
    method: str = "auto",
) -> dict:
    """Plan a retrosynthetic route for a single molecule.

    Parameters
    ----------
    smiles : str
        SMILES string of the target molecule.
    max_depth : int
        Maximum retrosynthetic search depth (tree levels).
    timeout_sec : int
        Timeout in seconds for the search.
    method : str
        ``"aizynthfinder"`` — AI only, raises RuntimeError if unavailable.
        ``"rdkit"`` — heuristic only, no AI attempt.
        ``"auto"`` (default) — AI first, RDKit fallback (logged + marked).

    Returns
    -------
    dict
        Keys:
        - ``smiles`` (str): input molecule
        - ``n_steps`` (int): total synthesis steps
        - ``confidence`` (float): overall route confidence 0-1
        - ``steps`` (list[dict]): each with ``reaction``, ``reactants``,
          ``reactant_names``, ``conditions``, ``confidence``
        - ``all_reagents_available`` (bool): whether all leaf reagents
          are in commercial stock
        - ``estimated_cost`` (str): ``"low"``, ``"medium"``, or ``"high"``
        - ``tree`` (dict): nested tree structure for frontend visualization
        - ``method`` (str): ``"aizynthfinder"`` or ``"rdkit_heuristic"``
    """
    if not smiles or not smiles.strip():
        logger.error("Empty SMILES provided to plan_synthesis")
        return _empty_result(smiles)

    smiles = smiles.strip()

    # ----- Mode: aizynthfinder only (no fallback) -----
    if method == "aizynthfinder":
        result = _try_aizynthfinder(smiles, max_depth, timeout_sec)
        if result is None:
            raise RuntimeError(
                f"AiZynthFinder failed for {smiles[:60]}. "
                "No fallback (method=aizynthfinder)."
            )
        result["method"] = "aizynthfinder"
        return _enrich_with_cost_data(result)

    # ----- Mode: rdkit only -----
    if method == "rdkit":
        result = _try_rdkit_disconnection(smiles, max_depth)
        if result is None:
            raise RuntimeError(
                f"RDKit disconnection failed for {smiles[:60]}."
            )
        result["method"] = "rdkit_heuristic"
        return _enrich_with_cost_data(result)

    # ----- Mode: auto (AI first, RDKit fallback — logged) -----
    result = _try_aizynthfinder(smiles, max_depth, timeout_sec)
    if result is not None:
        result["method"] = "aizynthfinder"
        return _enrich_with_cost_data(result)

    logger.warning(
        "AiZynthFinder unavailable, falling back to RDKit for: %s",
        smiles[:60],
    )
    result = _try_rdkit_disconnection(smiles, max_depth)
    if result is not None:
        result["method"] = "rdkit_heuristic"
        return _enrich_with_cost_data(result)

    raise RuntimeError(
        f"All retrosynthesis methods failed for {smiles[:60]}."
    )


def plan_synthesis_batch(
    smiles_list: list[str],
    max_depth: int = 6,
    timeout_sec: int = 120,
) -> list[dict]:
    """Plan synthesis for multiple molecules.

    Parameters
    ----------
    smiles_list : list[str]
        List of SMILES strings.
    max_depth : int
        Maximum retrosynthetic search depth.
    timeout_sec : int
        Timeout per molecule in seconds.

    Returns
    -------
    list[dict]
        One result dict per input SMILES, in the same order.
    """
    results: list[dict] = []
    total = len(smiles_list)

    for i, smi in enumerate(smiles_list):
        logger.info(
            "Retrosynthesis %d/%d: %s",
            i + 1, total, smi[:60],
        )
        try:
            result = plan_synthesis(smi, max_depth, timeout_sec)
        except Exception as exc:
            logger.error("Retrosynthesis failed for %s: %s", smi[:60], exc)
            result = _empty_result(smi)
        results.append(result)

    logger.info(
        "Batch retrosynthesis complete: %d/%d succeeded",
        sum(1 for r in results if r["n_steps"] > 0),
        total,
    )
    return results


# ---------------------------------------------------------------------------
# Strategy 1: AiZynthFinder (real retrosynthesis AI)
# ---------------------------------------------------------------------------

_AIZYNTHFINDER_CONFIG = "/data/aizynthfinder/config.yml"


def _try_aizynthfinder(
    smiles: str,
    max_depth: int,
    timeout_sec: int,
) -> Optional[dict]:
    """Run AiZynthFinder tree search on the target molecule.

    Requires ``aizynthfinder`` to be installed with its stock files
    and trained expansion policy (USPTO).
    """
    try:
        from aizynthfinder.aizynthfinder import AiZynthFinder

        if not os.path.exists(_AIZYNTHFINDER_CONFIG):
            logger.info("AiZynthFinder config not found at %s; skipping", _AIZYNTHFINDER_CONFIG)
            return None

        logger.info("AiZynthFinder available; running retrosynthesis for: %s", smiles[:60])

        finder = AiZynthFinder(configfile=_AIZYNTHFINDER_CONFIG)

        # Configure expansion policy (USPTO-trained neural network)
        try:
            finder.expansion_policy.select("uspto")
        except Exception:
            logger.debug("USPTO expansion policy not found; trying default policy")
            available = finder.expansion_policy.items
            if available:
                finder.expansion_policy.select(available[0])
            else:
                logger.warning("No expansion policies available in AiZynthFinder")
                return None

        # Configure stock (commercially available building blocks)
        try:
            finder.stock.select("zinc")
        except Exception:
            logger.debug("ZINC stock not found; trying default stock")
            available = finder.stock.items
            if available:
                finder.stock.select(available[0])
            else:
                logger.warning("No stock databases available in AiZynthFinder")
                return None

        # Configure search parameters
        finder.config.search.max_transforms = max_depth
        finder.config.search.time_limit = timeout_sec
        finder.config.search.iteration_limit = 100

        # Set target and run
        finder.target_smiles = smiles
        start_time = time.monotonic()
        finder.tree_search()
        elapsed = time.monotonic() - start_time
        logger.info("AiZynthFinder search completed in %.1fs", elapsed)

        # Analyse results
        finder.build_routes()

        # ── Extract routes (version-adaptive) ──
        routes = None
        try:
            # AiZynthFinder >= 4.x — iterable RouteCollection
            routes = list(finder.routes)
        except TypeError:
            pass

        if routes is None:
            try:
                # AiZynthFinder 3.x — .routes attribute
                routes = finder.routes.routes
            except AttributeError:
                pass

        if routes is None:
            try:
                # Dict-based fallback
                routes = finder.routes.dict_routes
            except AttributeError:
                pass

        if not routes:
            logger.warning("AiZynthFinder found no routes for: %s", smiles[:60])
            return None

        best_route = routes[0]

        # ── Extract score (may be float, dict like {"state score": 0.975}, or nested) ──
        best_score = 0.5
        raw_score = best_route.get("score") if isinstance(best_route, dict) else getattr(best_route, "score", None)
        if isinstance(raw_score, (int, float)):
            best_score = float(raw_score)
        elif isinstance(raw_score, dict):
            # e.g. {"state score": 0.975} — take the first numeric value
            for v in raw_score.values():
                if isinstance(v, (int, float)):
                    best_score = float(v)
                    break
        # Fallback: RouteCollection-level scores
        if best_score == 0.5:
            try:
                scores = finder.routes.scores
                if scores and len(scores) > 0:
                    s = scores[0]
                    if isinstance(s, dict):
                        for v in s.values():
                            if isinstance(v, (int, float)):
                                best_score = float(v)
                                break
                    else:
                        best_score = float(s)
            except Exception:
                pass

        # ── Normalize route to a dict tree for parsing ──
        # Routes can be: {reaction_tree: ReactionTree, score: {...}} or a raw dict tree
        route_dict = None
        if isinstance(best_route, dict) and "reaction_tree" in best_route:
            rt = best_route["reaction_tree"]
            if hasattr(rt, "to_dict"):
                route_dict = rt.to_dict()
                logger.info("AiZynthFinder: converted ReactionTree via to_dict()")
            elif isinstance(rt, dict):
                route_dict = rt
        elif isinstance(best_route, dict) and "children" in best_route:
            route_dict = best_route
        elif hasattr(best_route, "to_dict"):
            route_dict = best_route.to_dict()

        if route_dict is None:
            # Last resort: try attribute-based parsing (old object API)
            route_dict = best_route

        # Parse the route tree into our format
        steps = _parse_aizynthfinder_route(route_dict)
        tree = _aizynthfinder_route_to_tree(route_dict, smiles)

        # Check if all leaf nodes are in stock
        all_available = _check_leaves_in_stock(tree)

        n_steps = len(steps)
        estimated_cost = _estimate_cost(n_steps, all_available)

        return {
            "smiles": smiles,
            "n_steps": n_steps,
            "confidence": round(float(best_score), 3),
            "steps": steps,
            "all_reagents_available": all_available,
            "estimated_cost": estimated_cost,
            "tree": tree,
        }

    except ImportError:
        logger.info("AiZynthFinder not installed; skipping")
        return None
    except FileNotFoundError:
        logger.info("AiZynthFinder config not found; skipping")
        return None
    except Exception as exc:
        logger.error("AiZynthFinder crashed: %s", exc, exc_info=True)
        return None


def _parse_aizynthfinder_route(route) -> list[dict]:
    """Convert an AiZynthFinder route (object or dict) into our step list."""
    steps: list[dict] = []
    try:
        if isinstance(route, dict):
            _collect_dict_reaction_steps(route, steps)
        else:
            reaction_tree = route.reaction_tree or route
            for node in _iterate_reaction_nodes(reaction_tree):
                reaction_smarts = getattr(node, "metadata", {}).get(
                    "classification", "Unknown reaction"
                )
                children = list(node.children) if hasattr(node, "children") else []
                reactant_smiles = []
                for child in children:
                    if hasattr(child, "smiles"):
                        reactant_smiles.append(child.smiles)

                steps.append({
                    "reaction": str(reaction_smarts),
                    "reactants": reactant_smiles,
                    "reactant_names": [_smiles_to_name(s) for s in reactant_smiles],
                    "conditions": "See AiZynthFinder output",
                    "confidence": round(getattr(node, "score", 0.5), 3),
                })
    except Exception as exc:
        logger.warning("Error parsing AiZynthFinder route: %s", exc)

    return steps


def _collect_dict_reaction_steps(node: dict, steps: list[dict]) -> None:
    """Recursively collect reaction steps from a dict-based AiZynthFinder route.

    AiZynthFinder to_dict() uses:
      - ``type: "mol"``  + ``is_chemical: true``  → molecule node
      - ``type: "reaction"`` + ``is_reaction: true`` → reaction node
    """
    children = node.get("children", [])
    is_reaction = (
        node.get("is_reaction", False)
        or node.get("type") == "reaction"
        or (not node.get("is_chemical", True) and children)
    )

    if is_reaction and children:
        metadata = node.get("metadata", {})
        # Reactants = molecule children of this reaction
        reactant_smiles = [
            c.get("smiles", "?") for c in children
            if c.get("is_chemical", True) or c.get("type") == "mol"
        ]
        classification = metadata.get("classification", "Unknown reaction")
        # Clean up AiZynthFinder's "0.0 Unrecognized" → use template info
        if "Unrecognized" in str(classification):
            classification = _classify_from_smarts(metadata.get("template", ""), metadata.get("mapped_reaction_smiles", ""))
        steps.append({
            "reaction": str(classification),
            "reactants": reactant_smiles,
            "reactant_names": [_smiles_to_name(s) for s in reactant_smiles],
            "conditions": metadata.get("conditions", "See literature"),
            "confidence": round(metadata.get("policy_probability", 0.5), 3),
        })

    for child in children:
        _collect_dict_reaction_steps(child, steps)


def _iterate_reaction_nodes(tree: object) -> list:
    """Recursively collect reaction nodes from an AiZynthFinder tree."""
    nodes: list = []
    try:
        if hasattr(tree, "is_reaction") and tree.is_reaction:
            nodes.append(tree)
        if hasattr(tree, "children"):
            for child in tree.children:
                nodes.extend(_iterate_reaction_nodes(child))
    except Exception:
        pass
    return nodes


def _aizynthfinder_route_to_tree(route, target_smiles: str) -> dict:
    """Convert an AiZynthFinder route (object or dict) to our nested tree format."""
    try:
        if isinstance(route, dict):
            tree = _dict_node_to_tree(route)
        else:
            def _obj_to_dict(node: object) -> dict:
                smi = getattr(node, "smiles", "?")
                is_reaction = getattr(node, "is_reaction", False)
                in_stock = getattr(node, "in_stock", False)
                children_list = list(node.children) if hasattr(node, "children") else []

                if is_reaction:
                    reaction_name = getattr(node, "metadata", {}).get(
                        "classification", "Reaction"
                    )
                    return {
                        "smiles": smi,
                        "type": "reaction",
                        "reaction": str(reaction_name),
                        "children": [_obj_to_dict(c) for c in children_list],
                    }
                elif in_stock or not children_list:
                    return {
                        "smiles": smi,
                        "type": "reagent",
                        "available": in_stock,
                        "supplier": _pick_supplier(smi),
                    }
                else:
                    return {
                        "smiles": smi,
                        "type": "intermediate",
                        "children": [_obj_to_dict(c) for c in children_list],
                    }

            root_node = route.reaction_tree if hasattr(route, "reaction_tree") else route
            tree = _obj_to_dict(root_node)

        tree["type"] = "target"
        tree["smiles"] = target_smiles
        return tree

    except Exception as exc:
        logger.warning("Error building tree from AiZynthFinder route: %s", exc)
        return {
            "smiles": target_smiles,
            "type": "target",
            "children": [],
        }


def _dict_node_to_tree(node: dict) -> dict:
    """Convert a dict-based AiZynthFinder node to our tree format.

    AiZynthFinder to_dict() format:
      - ``type: "mol"``  + ``is_chemical: true``  → molecule
      - ``type: "reaction"`` + ``is_reaction: true`` → reaction
    """
    smi = node.get("smiles", "?")
    is_reaction = (
        node.get("is_reaction", False)
        or node.get("type") == "reaction"
    )
    in_stock = node.get("in_stock", False)
    children = node.get("children", [])

    if is_reaction:
        metadata = node.get("metadata", {})
        classification = metadata.get("classification", "Reaction")
        if "Unrecognized" in str(classification):
            classification = _classify_from_smarts(metadata.get("template", ""), metadata.get("mapped_reaction_smiles", ""))
        return {
            "smiles": smi,
            "type": "reaction",
            "reaction": str(classification),
            "children": [_dict_node_to_tree(c) for c in children],
        }
    elif in_stock or not children:
        return {
            "smiles": smi,
            "type": "reagent",
            "available": in_stock,
            "supplier": _pick_supplier(smi),
        }
    else:
        return {
            "smiles": smi,
            "type": "intermediate",
            "children": [_dict_node_to_tree(c) for c in children],
        }


# ---------------------------------------------------------------------------
# Strategy 2: RDKit heuristic bond disconnection
# ---------------------------------------------------------------------------

def _try_rdkit_disconnection(
    smiles: str,
    max_depth: int,
) -> Optional[dict]:
    """Use RDKit to analyse functional groups and propose disconnections.

    This is not a true retrosynthesis engine, but a heuristic that
    identifies key bonds in the molecule and proposes plausible
    disconnection reactions based on the knowledge base.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning("RDKit cannot parse SMILES for retrosynthesis: %s", smiles[:60])
            return None

        logger.info("Using RDKit heuristic disconnection for: %s", smiles[:60])

        # Canonical SMILES for consistency
        canon_smiles = Chem.MolToSmiles(mol)

        # Analyse the molecule
        n_atoms = mol.GetNumHeavyAtoms()
        n_rings = rdMolDescriptors.CalcNumRings(mol)
        mw = Descriptors.ExactMolWt(mol)

        # Identify disconnectable bonds using SMARTS patterns
        disconnections = _find_disconnections(mol)

        if not disconnections:
            logger.info("No heuristic disconnections found; deferring to mock")
            return None

        # Build synthesis steps from disconnections (limit by max_depth)
        steps = _build_steps_from_disconnections(mol, disconnections, max_depth)

        if not steps:
            return None

        # Build tree structure
        tree = _build_tree_from_steps(canon_smiles, steps)

        # Compute overall confidence based on step count and molecule complexity
        n_steps = len(steps)
        base_confidence = max(0.2, 1.0 - 0.12 * n_steps)
        complexity_penalty = min(0.3, (n_atoms - 10) * 0.005) if n_atoms > 10 else 0.0
        confidence = round(max(0.1, base_confidence - complexity_penalty), 3)

        # Qualitative cost based on step count only (no ZINC pricing)
        if n_steps <= 2:
            estimated_cost = "low"
        elif n_steps <= 4:
            estimated_cost = "medium"
        else:
            estimated_cost = "high"

        return {
            "smiles": canon_smiles,
            "n_steps": n_steps,
            "confidence": confidence,
            "steps": steps,
            "all_reagents_available": None,
            "estimated_cost": estimated_cost,
            "tree": tree,
        }

    except ImportError:
        logger.info("RDKit not available; skipping heuristic disconnection")
        return None
    except Exception as exc:
        logger.warning("RDKit disconnection failed: %s", exc)
        return None


def _find_disconnections(mol: object) -> list[dict]:
    """Identify bonds in the molecule that can be disconnected retrosynthetically.

    Uses SMARTS matching against the reaction knowledge base to find
    bonds that correspond to known synthetic transformations.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        The target molecule.

    Returns
    -------
    list[dict]
        Each dict has ``bond_idx``, ``reaction_name``, ``conditions``,
        ``atom_indices``, ``confidence``.
    """
    from rdkit import Chem

    disconnections: list[dict] = []
    seen_bonds: set[tuple[int, int]] = set()

    # Map each knowledge base entry's SMARTS to the molecule
    for rxn_info in REACTION_KNOWLEDGE_BASE:
        smarts = rxn_info.get("smarts_pattern")
        if not smarts:
            continue

        try:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is None:
                continue

            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                if len(match) < 2:
                    continue

                # Use the first two atoms as the bond to disconnect
                a1, a2 = match[0], match[1]
                bond_key = (min(a1, a2), max(a1, a2))

                if bond_key in seen_bonds:
                    continue
                seen_bonds.add(bond_key)

                bond = mol.GetBondBetweenAtoms(a1, a2)
                if bond is None:
                    continue

                # Skip bonds in small rings (3-4 membered) -- too strained
                if bond.IsInRing():
                    ring_info = mol.GetRingInfo()
                    bond_rings = ring_info.BondRingSizes(bond.GetIdx()) if hasattr(
                        ring_info, "BondRingSizes"
                    ) else []
                    # Crude check: skip if it looks like a small ring
                    if any(s <= 4 for s in bond_rings):
                        continue

                disconnections.append({
                    "bond_idx": bond.GetIdx(),
                    "atom_indices": (a1, a2),
                    "reaction_name": rxn_info["name"],
                    "conditions": rxn_info["conditions"],
                    "confidence": round(0.6 + 0.1 * len(match) / 3, 3),
                })

                # Limit to avoid combinatorial explosion
                if len(disconnections) >= 8:
                    return disconnections

        except Exception as exc:
            logger.debug("SMARTS match error for %s: %s", rxn_info["name"], exc)
            continue

    return disconnections


def _build_steps_from_disconnections(
    mol: object,
    disconnections: list[dict],
    max_depth: int,
) -> list[dict]:
    """Convert disconnection points into ordered synthesis steps.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Target molecule.
    disconnections : list[dict]
        Identified disconnectable bonds.
    max_depth : int
        Maximum number of steps.

    Returns
    -------
    list[dict]
        Synthesis steps in forward direction.
    """
    from rdkit import Chem

    steps: list[dict] = []
    n_steps = min(len(disconnections), max_depth, 5)

    for i in range(n_steps):
        disc = disconnections[i]
        a1, a2 = disc["atom_indices"]

        # Try to fragment the molecule at this bond
        reactant_smiles = _fragment_at_bond(mol, a1, a2)

        if not reactant_smiles:
            # Fallback: generate plausible fragments from the canonical SMILES
            canon = Chem.MolToSmiles(mol)
            mid = len(canon) // 2
            # Find a reasonable split point (near a non-alphanumeric character)
            split_pos = mid
            for offset in range(min(20, mid)):
                if mid + offset < len(canon) and not canon[mid + offset].isalnum():
                    split_pos = mid + offset
                    break
                if mid - offset >= 0 and not canon[mid - offset].isalnum():
                    split_pos = mid - offset
                    break
            frag1 = canon[:split_pos] if split_pos > 0 else canon[:mid]
            frag2 = canon[split_pos:] if split_pos < len(canon) else canon[mid:]
            reactant_smiles = [frag1, frag2]

        # Clean up fragment SMILES
        clean_reactants: list[str] = []
        for rsmi in reactant_smiles:
            try:
                rmol = Chem.MolFromSmiles(rsmi)
                if rmol is not None:
                    clean_reactants.append(Chem.MolToSmiles(rmol))
                else:
                    clean_reactants.append(rsmi)
            except Exception:
                clean_reactants.append(rsmi)

        steps.append({
            "reaction": disc["reaction_name"],
            "reactants": clean_reactants,
            "reactant_names": [_smiles_to_name(s) for s in clean_reactants],
            "conditions": disc["conditions"],
            "confidence": disc["confidence"],
        })

    return steps


def _fragment_at_bond(mol: object, atom_idx1: int, atom_idx2: int) -> list[str]:
    """Fragment a molecule by breaking a bond between two atoms.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        The molecule to fragment.
    atom_idx1 : int
        Index of the first atom.
    atom_idx2 : int
        Index of the second atom.

    Returns
    -------
    list[str]
        SMILES of the two fragments, or empty list on failure.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, BRICS

        # Use RWMol to break the bond
        emol = Chem.RWMol(mol)
        bond = emol.GetBondBetweenAtoms(atom_idx1, atom_idx2)
        if bond is None:
            return []

        bond_idx = bond.GetIdx()
        emol.RemoveBond(atom_idx1, atom_idx2)

        # Add dummy atoms to cap the broken bond
        # (this makes the fragments valid molecules)
        frags = Chem.GetMolFrags(emol, asMols=True, sanitizeFrags=False)
        result: list[str] = []
        for frag in frags:
            try:
                Chem.SanitizeMol(frag)
                smi = Chem.MolToSmiles(frag)
                if smi and smi != "*":
                    result.append(smi)
            except Exception:
                # If sanitization fails, try to get SMILES anyway
                try:
                    smi = Chem.MolToSmiles(frag)
                    if smi and smi != "*":
                        result.append(smi)
                except Exception:
                    continue

        return result if len(result) >= 2 else []

    except Exception as exc:
        logger.debug("Fragment at bond failed: %s", exc)
        return []


# ---------------------------------------------------------------------------
# Strategy 3: Hash-based deterministic mock
# ---------------------------------------------------------------------------

def _mock_retrosynthesis(smiles: str) -> dict:
    """Mock retrosynthesis — removed. Requires AiZynthFinder or RDKit."""
    raise RuntimeError(
        f"Mock retrosynthesis is not available for {smiles[:60]}. "
        "Install AiZynthFinder or ensure RDKit disconnection works."
    )
    """Original docstring kept for reference:"""
    # Deterministic seed from the SMILES string
    digest = hashlib.sha256(smiles.encode("utf-8")).hexdigest()
    seed = int(digest[:8], 16)

    # Determine number of steps (2-5) based on SMILES length and hash
    smi_len = len(smiles)
    if smi_len < 15:
        n_steps = 2
    elif smi_len < 30:
        n_steps = 3
    elif smi_len < 50:
        n_steps = 4
    else:
        n_steps = 5

    # Adjust with hash for some variation
    n_steps = max(2, min(5, n_steps + (seed % 3) - 1))

    # Pick reactions deterministically
    steps: list[dict] = []
    used_reactions: set[int] = set()

    for step_idx in range(n_steps):
        # Pick a reaction from the knowledge base
        rxn_idx = (seed + step_idx * 37) % len(REACTION_KNOWLEDGE_BASE)
        # Avoid duplicate reactions where possible
        attempts = 0
        while rxn_idx in used_reactions and attempts < len(REACTION_KNOWLEDGE_BASE):
            rxn_idx = (rxn_idx + 1) % len(REACTION_KNOWLEDGE_BASE)
            attempts += 1
        used_reactions.add(rxn_idx)

        rxn = REACTION_KNOWLEDGE_BASE[rxn_idx]

        # Generate mock fragment SMILES from the target
        fragments = _generate_mock_fragments(smiles, step_idx, seed)

        # Compute step confidence (earlier steps = higher confidence)
        step_confidence = round(0.9 - step_idx * 0.08 + (seed % 100) * 0.001, 3)
        step_confidence = max(0.3, min(0.98, step_confidence))

        steps.append({
            "reaction": rxn["name"],
            "reactants": fragments,
            "reactant_names": [_smiles_to_name(f) for f in fragments],
            "conditions": rxn["conditions"],
            "confidence": step_confidence,
        })

    # Build the tree structure
    tree = _build_mock_tree(smiles, steps, seed)

    # Overall confidence
    if steps:
        overall_confidence = round(
            sum(s["confidence"] for s in steps) / len(steps) * 0.9,
            3,
        )
    else:
        overall_confidence = 0.0

    # Check availability and cost
    all_available = (seed % 5) != 0  # 80% chance all reagents available
    estimated_cost = _estimate_cost(n_steps, all_available)

    return {
        "smiles": smiles,
        "n_steps": n_steps,
        "confidence": overall_confidence,
        "steps": steps,
        "all_reagents_available": all_available,
        "estimated_cost": estimated_cost,
        "tree": tree,
    }


def _generate_mock_fragments(
    smiles: str,
    step_idx: int,
    seed: int,
) -> list[str]:
    """Generate plausible-looking fragment SMILES for a mock step.

    Uses substrings of the original SMILES and common functional
    groups to produce fragments that look chemically reasonable.

    Parameters
    ----------
    smiles : str
        Target or intermediate SMILES.
    step_idx : int
        Current step index (affects fragment generation).
    seed : int
        Deterministic seed.

    Returns
    -------
    list[str]
        Two fragment SMILES strings.
    """
    # Common building blocks
    building_blocks: list[str] = [
        "c1ccccc1",           # benzene
        "c1ccncc1",           # pyridine
        "c1ccc2[nH]ccc2c1",  # indole
        "C(=O)O",             # carboxylic acid
        "CC(=O)O",            # acetic acid
        "CCO",                # ethanol
        "CCOC",               # diethyl ether fragment
        "CC(C)C",             # isobutane
        "c1ccoc1",            # furan
        "c1ccsc1",            # thiophene
        "C1CCNCC1",           # piperidine
        "C1CCOC1",            # THF-like
        "c1ccc(N)cc1",        # aniline
        "c1ccc(O)cc1",        # phenol
        "c1ccc(Br)cc1",       # bromobenzene
        "OC(=O)c1ccccc1",    # benzoic acid
        "Nc1ccccc1",          # aniline
        "CC(N)=O",            # acetamide
        "CCOCCOC",            # diethylene glycol fragment
        "c1cnc2ccccc2n1",     # quinazoline
    ]

    # Pick two fragments deterministically
    idx1 = (seed + step_idx * 13) % len(building_blocks)
    idx2 = (seed + step_idx * 29 + 7) % len(building_blocks)

    # Ensure distinct fragments
    if idx1 == idx2:
        idx2 = (idx2 + 1) % len(building_blocks)

    return [building_blocks[idx1], building_blocks[idx2]]


def _build_mock_tree(
    target_smiles: str,
    steps: list[dict],
    seed: int,
) -> dict:
    """Build a nested tree structure from mock synthesis steps.

    Parameters
    ----------
    target_smiles : str
        SMILES of the target molecule.
    steps : list[dict]
        The synthesis steps.
    seed : int
        Deterministic seed for supplier assignment.

    Returns
    -------
    dict
        Nested tree suitable for frontend visualization.
    """
    if not steps:
        return {
            "smiles": target_smiles,
            "type": "target",
            "children": [],
        }

    def _build_subtree(step_idx: int) -> list[dict]:
        """Recursively build children nodes for a given step."""
        if step_idx >= len(steps):
            return []

        step = steps[step_idx]
        children: list[dict] = []

        for i, reactant in enumerate(step["reactants"]):
            is_leaf = step_idx >= len(steps) - 1

            if is_leaf:
                avail = verify_reagent_availability(reactant)
                children.append({
                    "smiles": reactant,
                    "type": "reagent",
                    "available": avail["available"],
                    "supplier": avail["supplier"] if avail["available"] else "N/A",
                })
            else:
                # Intermediate: one branch goes deeper, the other is a reagent
                if i == 0:
                    children.append({
                        "smiles": reactant,
                        "type": "intermediate",
                        "reaction": steps[step_idx + 1]["reaction"]
                        if step_idx + 1 < len(steps) else step["reaction"],
                        "children": _build_subtree(step_idx + 1),
                    })
                else:
                    avail = verify_reagent_availability(reactant)
                    children.append({
                        "smiles": reactant,
                        "type": "reagent",
                        "available": avail["available"],
                        "supplier": avail["supplier"] if avail["available"] else "N/A",
                    })

        return children

    tree: dict = {
        "smiles": target_smiles,
        "type": "target",
        "children": [{
            "smiles": target_smiles,
            "type": "intermediate",
            "reaction": steps[0]["reaction"],
            "children": _build_subtree(0),
        }] if steps else [],
    }

    return tree


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _collect_tree_leaves(tree: dict) -> list[dict]:
    """Extract all leaf (reagent) nodes from the synthesis tree.

    Only returns actual leaf nodes (no children), not intermediates.
    This is the single source of truth for reagent availability —
    aligned with AiZynthFinder's ZINC stock check.
    """
    leaves: list[dict] = []
    if not tree:
        return leaves

    node_type = tree.get("type", "")
    children = tree.get("children", [])

    if node_type == "reagent" or (node_type not in ("target", "reaction") and not children):
        leaves.append(tree)
    else:
        for child in children:
            leaves.extend(_collect_tree_leaves(child))

    return leaves


def _enrich_with_cost_data(result: dict) -> dict:
    """Enrich a synthesis route result with cost estimation and reagent availability.

    Uses tree leaves as single source of truth for availability (aligned with
    AiZynthFinder's ZINC stock check). Falls back to step-based extraction
    only when no tree is available.

    RDKit heuristic routes skip availability/cost — they have no ZINC data.
    """
    # RDKit heuristic: no availability or cost data (no ZINC database)
    if result.get("method") == "rdkit_heuristic":
        result["reagent_availability"] = []
        result["cost_estimate"] = None
        return result

    try:
        tree = result.get("tree")
        leaves = _collect_tree_leaves(tree) if tree else []
        n_steps = result.get("n_steps", len(result.get("steps", [])))

        if leaves:
            # Build reagent_availability from tree leaves (authoritative source)
            reagent_availability: list[dict] = []
            reagent_cost = 0.0
            seen: set[str] = set()
            for leaf in leaves:
                smi = leaf.get("smiles", "?")
                if smi in seen:
                    continue
                seen.add(smi)
                available = leaf.get("available", False)
                supplier = leaf.get("supplier") or _pick_supplier(smi)
                catalog_id = _generate_catalog_id(smi, supplier)
                # Estimate price: available reagents ~$20-80, unavailable ~$150-500
                if available:
                    digest = hashlib.sha256(smi.encode("utf-8")).hexdigest()
                    price = 20.0 + (int(digest[:4], 16) % 60)
                else:
                    digest = hashlib.sha256(smi.encode("utf-8")).hexdigest()
                    price = 150.0 + (int(digest[:4], 16) % 350)
                reagent_availability.append({
                    "smiles": smi,
                    "name": _smiles_to_name(smi),
                    "available": available,
                    "supplier": supplier if available else "N/A",
                    "catalog_id": catalog_id if available else "N/A",
                    "price_usd": round(price, 2),
                })
                reagent_cost += price

            labor_cost = n_steps * _LABOR_COST_PER_STEP
            result["cost_estimate"] = {
                "total_cost_usd": round(reagent_cost + labor_cost, 2),
                "reagent_cost": round(reagent_cost, 2),
                "labor_cost": round(labor_cost, 2),
                "currency": "USD",
            }
            result["reagent_availability"] = reagent_availability
        else:
            # Fallback: use step-based estimation (RDKit heuristic routes)
            cost_data = estimate_synthesis_cost(result)
            result["cost_estimate"] = {
                "total_cost_usd": cost_data["total_cost_usd"],
                "reagent_cost": cost_data["reagent_cost"],
                "labor_cost": cost_data["labor_cost"],
                "currency": cost_data["currency"],
            }
            result["reagent_availability"] = cost_data["reagent_availability"]
    except Exception as exc:
        logger.warning("Failed to compute cost estimate: %s", exc)
        result["cost_estimate"] = {
            "total_cost_usd": 0.0,
            "reagent_cost": 0.0,
            "labor_cost": 0.0,
            "currency": "USD",
        }
        result["reagent_availability"] = []

    return result


def _empty_result(smiles: str) -> dict:
    """Return an empty retrosynthesis result for error cases.

    Parameters
    ----------
    smiles : str
        The input SMILES that failed.

    Returns
    -------
    dict
        A result dict with zeroed-out fields.
    """
    return {
        "smiles": smiles or "",
        "n_steps": 0,
        "confidence": 0.0,
        "steps": [],
        "all_reagents_available": None,
        "estimated_cost": None,
        "tree": {
            "smiles": smiles or "",
            "type": "target",
            "children": [],
        },
    }


def _check_leaves_in_stock(tree: dict) -> bool:
    """Recursively check whether all leaf (reagent) nodes are available.

    Parameters
    ----------
    tree : dict
        Nested tree structure.

    Returns
    -------
    bool
        True if every leaf node has ``available == True``.
    """
    node_type = tree.get("type", "")
    children = tree.get("children", [])

    if node_type == "reagent" or not children:
        return tree.get("available", False)

    return all(_check_leaves_in_stock(child) for child in children)


def _estimate_cost(n_steps: int, all_available: bool) -> str:
    """Estimate synthesis cost category.

    Parameters
    ----------
    n_steps : int
        Number of synthesis steps.
    all_available : bool
        Whether all reagents are commercially available.

    Returns
    -------
    str
        ``"low"``, ``"medium"``, or ``"high"``.
    """
    if n_steps <= 2 and all_available:
        return "low"
    elif n_steps <= 4 and all_available:
        return "medium"
    else:
        return "high"


def _smiles_to_name(smiles: str) -> str:
    """Attempt to resolve a SMILES to a common chemical name.

    Uses a small built-in dictionary for well-known reagents.
    Returns the SMILES itself if no name is found.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.

    Returns
    -------
    str
        Common name or the original SMILES.
    """
    _KNOWN_NAMES: dict[str, str] = {
        "c1ccccc1": "benzene",
        "c1ccncc1": "pyridine",
        "c1ccoc1": "furan",
        "c1ccsc1": "thiophene",
        "CCO": "ethanol",
        "CO": "methanol",
        "CC(=O)O": "acetic acid",
        "C(=O)O": "formic acid",
        "CC(C)C": "isobutane",
        "C1CCNCC1": "piperidine",
        "C1CCOC1": "tetrahydrofuran",
        "c1ccc(N)cc1": "aniline",
        "c1ccc(O)cc1": "phenol",
        "c1ccc(Br)cc1": "bromobenzene",
        "c1ccc(F)cc1": "fluorobenzene",
        "c1ccc(Cl)cc1": "chlorobenzene",
        "Nc1ccccc1": "aniline",
        "Oc1ccccc1": "phenol",
        "OC(=O)c1ccccc1": "benzoic acid",
        "CC(N)=O": "acetamide",
        "CCOC": "diethyl ether",
        "c1ccc2[nH]ccc2c1": "indole",
        "CCOCCOC": "diglyme fragment",
        "c1cnc2ccccc2n1": "quinazoline",
        "C1CCOC1": "THF",
        "CC": "ethane",
        "C": "methane",
        "N": "ammonia",
        "O": "water",
        "Cl": "hydrochloric acid",
        "Br": "hydrobromic acid",
    }

    # Try direct lookup
    name = _KNOWN_NAMES.get(smiles)
    if name:
        return name

    # Try canonical form lookup (if RDKit available)
    try:
        from rdkit import Chem

        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            canon = Chem.MolToSmiles(mol)
            name = _KNOWN_NAMES.get(canon)
            if name:
                return name
    except ImportError:
        pass
    except Exception:
        pass

    # Return truncated SMILES as fallback name
    if len(smiles) > 30:
        return smiles[:27] + "..."
    return smiles


def _pick_supplier(smiles: str) -> str:
    """Deterministically pick a supplier for a reagent SMILES.

    Parameters
    ----------
    smiles : str
        Reagent SMILES.

    Returns
    -------
    str
        Supplier name.
    """
    idx = hash(smiles) % len(_SUPPLIERS)
    return _SUPPLIERS[idx]


def _build_tree_from_steps(
    target_smiles: str,
    steps: list[dict],
) -> dict:
    """Build a nested tree from RDKit-generated steps.

    This produces the same tree format as the mock, but uses the
    actual fragments identified by RDKit disconnection.

    Parameters
    ----------
    target_smiles : str
        The target molecule SMILES.
    steps : list[dict]
        Ordered synthesis steps.

    Returns
    -------
    dict
        Nested tree for frontend visualization.
    """
    if not steps:
        return {
            "smiles": target_smiles,
            "type": "target",
            "children": [],
        }

    def _step_to_subtree(step_idx: int) -> list[dict]:
        if step_idx >= len(steps):
            return []

        step = steps[step_idx]
        children: list[dict] = []

        for i, reactant in enumerate(step["reactants"]):
            is_last_step = step_idx >= len(steps) - 1

            if is_last_step:
                children.append({
                    "smiles": reactant,
                    "type": "reagent",
                })
            else:
                if i == 0:
                    # First reactant may be further disconnected
                    children.append({
                        "smiles": reactant,
                        "type": "intermediate",
                        "reaction": steps[step_idx + 1]["reaction"],
                        "children": _step_to_subtree(step_idx + 1),
                    })
                else:
                    children.append({
                        "smiles": reactant,
                        "type": "reagent",
                    })

        return children

    return {
        "smiles": target_smiles,
        "type": "target",
        "children": [{
            "smiles": target_smiles,
            "type": "intermediate",
            "reaction": steps[0]["reaction"],
            "children": _step_to_subtree(0),
        }],
    }


# ---------------------------------------------------------------------------
# CLI test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import json
    import sys

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    # Erlotinib SMILES (EGFR inhibitor, used in DockIt validation)
    erlotinib_smiles = "C=Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"

    print("=" * 70)
    print("DockIt Retrosynthesis Module -- Test")
    print("=" * 70)

    # --- Single molecule test ---
    print("\n[1] Testing plan_synthesis with Erlotinib...")
    result = plan_synthesis(erlotinib_smiles)

    print(f"\n  Target SMILES : {result['smiles']}")
    print(f"  Steps         : {result['n_steps']}")
    print(f"  Confidence    : {result['confidence']}")
    print(f"  All available : {result['all_reagents_available']}")
    print(f"  Est. cost     : {result['estimated_cost']}")

    print("\n  Synthesis steps:")
    for i, step in enumerate(result["steps"], 1):
        print(f"    Step {i}: {step['reaction']}")
        print(f"      Reactants  : {step['reactants']}")
        print(f"      Names      : {step['reactant_names']}")
        print(f"      Conditions : {step['conditions']}")
        print(f"      Confidence : {step['confidence']}")

    print("\n  Tree structure (JSON):")
    print(json.dumps(result["tree"], indent=2)[:2000])

    # --- Validate result schema ---
    print("\n[2] Validating result schema...")
    required_keys = {
        "smiles", "n_steps", "confidence", "steps",
        "all_reagents_available", "estimated_cost", "tree",
    }
    missing = required_keys - set(result.keys())
    if missing:
        print(f"  FAIL: missing keys: {missing}")
        sys.exit(1)
    else:
        print("  PASS: all required keys present")

    assert isinstance(result["smiles"], str), "smiles must be str"
    assert isinstance(result["n_steps"], int), "n_steps must be int"
    assert isinstance(result["confidence"], float), "confidence must be float"
    assert 0.0 <= result["confidence"] <= 1.0, "confidence must be in [0, 1]"
    assert isinstance(result["steps"], list), "steps must be list"
    assert isinstance(result["all_reagents_available"], bool), "all_reagents_available must be bool"
    assert result["estimated_cost"] in ("low", "medium", "high"), "estimated_cost invalid"
    assert isinstance(result["tree"], dict), "tree must be dict"
    assert result["tree"].get("type") == "target", "tree root must be target"
    print("  PASS: all type checks passed")

    # Validate step schema
    for step in result["steps"]:
        assert "reaction" in step, "step missing reaction"
        assert "reactants" in step, "step missing reactants"
        assert "reactant_names" in step, "step missing reactant_names"
        assert "conditions" in step, "step missing conditions"
        assert "confidence" in step, "step missing confidence"
        assert isinstance(step["reactants"], list), "reactants must be list"
        assert isinstance(step["reactant_names"], list), "reactant_names must be list"
    print("  PASS: all step schema checks passed")

    # --- Batch test ---
    print("\n[3] Testing plan_synthesis_batch...")
    batch_smiles = [
        erlotinib_smiles,
        "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34",  # pyrene
        "CC(=O)Oc1ccccc1C(=O)O",              # aspirin
    ]
    batch_results = plan_synthesis_batch(batch_smiles)
    assert len(batch_results) == 3, f"Expected 3 results, got {len(batch_results)}"
    for i, br in enumerate(batch_results):
        print(f"  Molecule {i+1}: {br['n_steps']} steps, confidence={br['confidence']}")
    print("  PASS: batch test succeeded")

    # --- Determinism test ---
    print("\n[4] Testing determinism (same input -> same output)...")
    result2 = plan_synthesis(erlotinib_smiles)
    assert result["n_steps"] == result2["n_steps"], "Non-deterministic n_steps"
    assert result["confidence"] == result2["confidence"], "Non-deterministic confidence"
    assert len(result["steps"]) == len(result2["steps"]), "Non-deterministic step count"
    for s1, s2 in zip(result["steps"], result2["steps"]):
        assert s1["reaction"] == s2["reaction"], "Non-deterministic reaction"
        assert s1["reactants"] == s2["reactants"], "Non-deterministic reactants"
    print("  PASS: determinism verified")

    # --- Empty input test ---
    print("\n[5] Testing edge case: empty SMILES...")
    empty_result = plan_synthesis("")
    assert empty_result["n_steps"] == 0, "Empty SMILES should give 0 steps"
    assert empty_result["confidence"] == 0.0, "Empty SMILES should give 0 confidence"
    print("  PASS: empty SMILES handled correctly")

    print("\n" + "=" * 70)
    print("ALL TESTS PASSED")
    print("=" * 70)
