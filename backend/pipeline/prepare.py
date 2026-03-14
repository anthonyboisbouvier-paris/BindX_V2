"""
DockIt pipeline — Receptor and ligand preparation.

Receptor: 7-step state-of-the-art pipeline (P2→P9) with graceful fallbacks.
Ligand: SMILES -> PDBQT via RDKit + Meeko/obabel.

Each receptor step is optional — if a dependency is missing, it is skipped
with a warning. The pipeline remains robust even in degraded mode.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Receptor preparation — constants
# ---------------------------------------------------------------------------

# HETATM residues to REMOVE (waters, buffers, crystallisation artifacts)
_REMOVE_RESIDUES = frozenset({
    "HOH", "WAT", "DOD",           # waters
    "GOL", "PEG", "PG4", "P6G",    # polyethylene glycol / glycerol
    "SO4", "PO4", "ACT", "FMT",    # sulfate, phosphate, acetate, formate
    "DMS", "EDO", "BME", "TRS",    # DMSO, ethylene glycol, BME, tris
    "CIT", "MPD", "IOD", "BR",     # citrate, MPD, iodide, bromide
    "CL", "NA", "K",               # common counter-ions (not catalytic)
})

# HETATM residues to KEEP — catalytic cofactors
_KEEP_COFACTORS = frozenset({
    "NAD", "NAP", "NDP", "FAD", "FMN", "HEM", "HEC",
    "ATP", "ADP", "AMP", "GTP", "GDP", "SAM", "SAH",
    "COA", "ACO", "PLP", "TPP",
})

# HETATM metals to KEEP — catalytic metal ions
_KEEP_METALS = frozenset({
    "ZN", "MG", "FE", "FE2", "MN", "CA", "CO", "NI", "CU", "CU1",
    "MO", "SE", "W", "CD",
})


# ---------------------------------------------------------------------------
# P2 — Intelligent HETATM cleaning (replaces blind _strip_hetatm)
# ---------------------------------------------------------------------------

def _clean_pdb(
    pdb_path: Path,
    work_dir: Path,
    *,
    keep_metals: bool = True,
    keep_cofactors: bool = True,
) -> tuple[Path, dict]:
    """Selectively remove HETATM records, keeping metals & cofactors.

    Returns (cleaned_pdb_path, cleaning_report).
    """
    cleaned_path = work_dir / f"{pdb_path.stem}_clean.pdb"
    kept_atoms = 0
    removed = 0
    removed_residues: set[str] = set()
    kept_het: set[str] = set()

    with open(pdb_path, "r") as fin, open(cleaned_path, "w") as fout:
        for line in fin:
            record = line[:6].strip()

            if record == "HETATM":
                res_name = line[17:20].strip()

                # Always remove waters and buffers
                if res_name in _REMOVE_RESIDUES:
                    removed_residues.add(res_name)
                    removed += 1
                    continue

                # Keep cofactors if requested
                if keep_cofactors and res_name in _KEEP_COFACTORS:
                    fout.write(line)
                    kept_het.add(res_name)
                    continue

                # Keep metals if requested
                if keep_metals and res_name in _KEEP_METALS:
                    fout.write(line)
                    kept_het.add(res_name)
                    continue

                # Everything else (ligands, unknown) — remove
                removed_residues.add(res_name)
                removed += 1
                continue

            fout.write(line)
            if record == "ATOM":
                kept_atoms += 1

    report = {
        "atoms_kept": kept_atoms,
        "hetatm_removed": removed,
        "residues_removed": sorted(removed_residues),
        "het_kept": sorted(kept_het),
    }

    if removed > 0:
        logger.info(
            "P2 Cleaning: removed %d HETATM (%s), kept %d ATOM + %s from %s",
            removed, ", ".join(sorted(removed_residues)), kept_atoms,
            ", ".join(sorted(kept_het)) or "no HET",
            pdb_path.name,
        )

    return cleaned_path, report


# ---------------------------------------------------------------------------
# P4 — Resolve alternate conformations (altloc)
# ---------------------------------------------------------------------------

def _resolve_altloc(pdb_path: Path, work_dir: Path) -> tuple[Path, int]:
    """Keep only the highest-occupancy alternate conformation per residue.

    PDB column 17 (0-indexed: 16) is the altloc indicator.
    Returns (resolved_pdb_path, n_resolved).
    """
    output_path = work_dir / f"{pdb_path.stem}_altloc.pdb"
    n_resolved = 0

    # First pass: find best altloc per (chain, resSeq, iCode)
    best_alt: dict[tuple, str] = {}  # (chain, resSeq, iCode) -> best altloc char
    with open(pdb_path, "r") as fh:
        for line in fh:
            record = line[:6].strip()
            if record not in ("ATOM", "HETATM"):
                continue
            altloc = line[16] if len(line) > 16 else " "
            if altloc == " ":
                continue
            chain = line[21] if len(line) > 21 else " "
            res_seq = line[22:26].strip()
            icode = line[26] if len(line) > 26 else " "
            key = (chain, res_seq, icode)

            try:
                occ = float(line[54:60])
            except (ValueError, IndexError):
                occ = 1.0

            if key not in best_alt:
                best_alt[key] = (altloc, occ)
            elif occ > best_alt[key][1]:
                best_alt[key] = (altloc, occ)
            elif occ == best_alt[key][1] and altloc < best_alt[key][0]:
                # Equal occupancy — prefer A
                best_alt[key] = (altloc, occ)

    # Resolve: best altloc char per residue key
    best_chars = {k: v[0] for k, v in best_alt.items()}

    # Second pass: filter lines
    with open(pdb_path, "r") as fin, open(output_path, "w") as fout:
        for line in fin:
            record = line[:6].strip()
            if record not in ("ATOM", "HETATM"):
                fout.write(line)
                continue

            altloc = line[16] if len(line) > 16 else " "
            if altloc == " ":
                fout.write(line)
                continue

            chain = line[21] if len(line) > 21 else " "
            res_seq = line[22:26].strip()
            icode = line[26] if len(line) > 26 else " "
            key = (chain, res_seq, icode)

            if key in best_chars and altloc == best_chars[key]:
                # Write this line but clear the altloc indicator
                fout.write(line[:16] + " " + line[17:])
                n_resolved += 1
            # else: skip this alternate conformation

    if n_resolved > 0:
        logger.info("P4 Altloc: resolved %d residues (kept best occupancy)", len(best_chars))

    return output_path, n_resolved


# ---------------------------------------------------------------------------
# P3 + P5 — Non-standard residues + missing atoms (PDBFixer)
# ---------------------------------------------------------------------------

def _fix_with_pdbfixer(pdb_path: Path, work_dir: Path) -> tuple[Path, dict]:
    """Use PDBFixer for non-standard residues and missing atoms.

    Returns (fixed_pdb_path, fix_report).
    """
    report: dict = {"nonstandard_replaced": [], "missing_atoms_added": 0, "available": False}

    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except ImportError:
        logger.warning("P3/P5: PDBFixer not available — skipping non-standard residue fix and missing atoms")
        return pdb_path, report

    report["available"] = True

    try:
        fixer = PDBFixer(filename=str(pdb_path))

        # P3: Non-standard residues (MSE→MET, SEP→SER, etc.)
        fixer.findNonstandardResidues()
        if fixer.nonstandardResidues:
            conversions = []
            for residue, replacement_name in fixer.nonstandardResidues:
                conversions.append(f"{residue.name}→{replacement_name}")
            report["nonstandard_replaced"] = conversions
            logger.info("P3: replacing non-standard residues: %s", ", ".join(conversions))
            fixer.replaceNonstandardResidues()

        # P5: Missing atoms (side-chain atoms only — NOT missing residues/loops)
        fixer.findMissingResidues()
        # Clear missing residues to avoid rebuilding long loops (adds noise)
        fixer.missingResidues = {}

        fixer.findMissingAtoms()
        n_missing = sum(len(atoms) for atoms in fixer.missingAtoms.values())
        if n_missing > 0:
            logger.info("P5: adding %d missing heavy atoms", n_missing)
            fixer.addMissingAtoms()
            report["missing_atoms_added"] = n_missing

        # Save fixed structure
        output_path = work_dir / f"{pdb_path.stem}_fixed.pdb"
        with open(output_path, "w") as fh:
            PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)

        return output_path, report

    except Exception as exc:
        logger.warning("P3/P5: PDBFixer error — skipping: %s", exc)
        return pdb_path, report


# ---------------------------------------------------------------------------
# P6+P7 — Protonation analysis (PropKa) + add hydrogens (PDBFixer)
# ---------------------------------------------------------------------------

def _protonate_receptor(pdb_path: Path, work_dir: Path, ph: float = 7.4) -> tuple[Path, dict]:
    """Analyze pKa with PropKa, then add hydrogens at target pH.

    Returns (protonated_pdb_path, protonation_report).
    """
    report: dict = {"shifted_pka_residues": [], "ph": ph}

    # Step 1: PropKa pKa analysis (informational — warns about shifted residues)
    try:
        import propka.run as propka_run

        # propka.run.single() returns a MolecularContainer
        # write_pka=False avoids writing .pka file to disk
        mol_container = propka_run.single(str(pdb_path), write_pka=False)

        # Access pKa groups from the averaged conformation ('AVR')
        conformation = mol_container.conformations.get("AVR")
        if conformation is None and mol_container.conformation_names:
            conformation = mol_container.conformations[mol_container.conformation_names[0]]

        shifted = []
        if conformation:
            for group in conformation.get_titratable_groups():
                shift = abs(group.pka_value - group.model_pka)
                if shift > 1.0:
                    entry = {
                        "residue": f"{group.residue_type}{group.atom.res_num}",
                        "chain": getattr(group.atom, "chain_id", "?"),
                        "pka": round(group.pka_value, 1),
                        "model_pka": round(group.model_pka, 1),
                        "shift": round(shift, 1),
                    }
                    shifted.append(entry)
                    logger.warning(
                        "P6 pKa shift: %s (chain %s) pKa=%.1f (model=%.1f, shift=%.1f)",
                        entry["residue"], entry["chain"],
                        group.pka_value, group.model_pka, shift,
                    )

        report["shifted_pka_residues"] = shifted
        if shifted:
            logger.info("P6: %d residues with shifted pKa (>1 unit from model)", len(shifted))
        else:
            logger.info("P6: no significant pKa shifts detected")

    except ImportError:
        logger.warning("P6: PropKa not available — skipping pKa analysis")
    except Exception as exc:
        logger.warning("P6: PropKa error — skipping pKa analysis: %s", exc)

    # Step 2: Add hydrogens at target pH via PDBFixer
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile

        fixer = PDBFixer(filename=str(pdb_path))
        fixer.addMissingHydrogens(ph)

        output_path = work_dir / f"{pdb_path.stem}_h.pdb"
        with open(output_path, "w") as fh:
            PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)

        logger.info("P7: added hydrogens at pH %.1f", ph)
        report["hydrogens_added"] = True
        return output_path, report

    except ImportError:
        logger.warning("P7: PDBFixer not available — skipping hydrogen addition")
        report["hydrogens_added"] = False
        return pdb_path, report
    except Exception as exc:
        logger.warning("P7: hydrogen addition error — skipping: %s", exc)
        report["hydrogens_added"] = False
        return pdb_path, report


# ---------------------------------------------------------------------------
# P8 — Energy minimisation (OpenMM, conditional)
# ---------------------------------------------------------------------------

def _minimize_receptor(
    pdb_path: Path,
    work_dir: Path,
    max_iterations: int = 500,
) -> tuple[Path, dict]:
    """Minimise receptor with harmonic restraints on Cα atoms.

    Returns (minimized_pdb_path, minimization_report).
    """
    report: dict = {"minimization_converged": False, "method": "none"}

    try:
        from openmm.app import PDBFile, ForceField, Modeller, NoCutoff, HBonds
        from openmm import LangevinMiddleIntegrator, CustomExternalForce
        import openmm.unit as unit
    except ImportError:
        logger.warning("P8: OpenMM not available — skipping minimisation")
        return pdb_path, report

    try:
        pdb = PDBFile(str(pdb_path))
        # Use implicit solvent (GBn2) for vacuum minimisation — NOT explicit water tip3pfb
        forcefield = ForceField("amber14-all.xml", "implicit/gbn2.xml")
        modeller = Modeller(pdb.topology, pdb.positions)

        # Remove water if any got through
        modeller.deleteWater()

        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=NoCutoff,
            constraints=HBonds,
        )

        # Add harmonic restraints on Cα atoms to prevent cavity collapse
        # Use Cartesian distance formula (not periodicdistance — we use NoCutoff, no PBC)
        restraint = CustomExternalForce(
            "k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
        )
        restraint.addGlobalParameter(
            "k", 10.0 * unit.kilocalories_per_mole / unit.angstroms ** 2
        )
        restraint.addPerParticleParameter("x0")
        restraint.addPerParticleParameter("y0")
        restraint.addPerParticleParameter("z0")

        n_restrained = 0
        for atom in modeller.topology.atoms():
            if atom.name == "CA":
                pos = modeller.positions[atom.index]
                # Pass position Vec3 directly — OpenMM auto-strips units
                restraint.addParticle(
                    atom.index, pos,
                )
                n_restrained += 1

        system.addForce(restraint)

        integrator = LangevinMiddleIntegrator(
            300 * unit.kelvin,
            1.0 / unit.picoseconds,
            0.002 * unit.picoseconds,
        )

        from openmm import Platform
        # Use CPU platform (always available)
        platform = Platform.getPlatformByName("CPU")
        from openmm.app import Simulation
        simulation = Simulation(
            modeller.topology, system, integrator, platform
        )
        simulation.context.setPositions(modeller.positions)

        # Minimise
        simulation.minimizeEnergy(maxIterations=max_iterations)

        # Save minimised structure
        output_path = work_dir / f"{pdb_path.stem}_min.pdb"
        positions = simulation.context.getState(getPositions=True).getPositions()
        with open(output_path, "w") as fh:
            PDBFile.writeFile(modeller.topology, positions, fh, keepIds=True)

        logger.info(
            "P8: minimised with %d Cα restraints (%d iterations max)",
            n_restrained, max_iterations,
        )
        report["minimization_converged"] = True
        report["method"] = "amber14_restrained"
        report["ca_restraints"] = n_restrained
        return output_path, report

    except Exception as exc:
        logger.warning("P8: minimisation error — skipping: %s", exc)
        return pdb_path, report


# ---------------------------------------------------------------------------
# P9 — PDBQT conversion (Meeko > obabel > manual)
# ---------------------------------------------------------------------------

def _meeko_pdb_to_pdbqt(pdb_path: Path, pdbqt_path: Path) -> bool:
    """Use Meeko mk_prepare_receptor for PDB → PDBQT.

    Important: mk_prepare_receptor uses -o as a BASENAME (not full path),
    requires -p flag to write PDBQT, and outputs {basename}_rigid.pdbqt.
    """
    # Compute basename and expected output filename
    basename = str(pdbqt_path).replace(".pdbqt", "")
    rigid_output = Path(f"{basename}_rigid.pdbqt")

    def _try_meeko_cmd(cmd: list[str]) -> bool:
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            if result.returncode == 0 and rigid_output.exists() and rigid_output.stat().st_size > 100:
                # Rename {basename}_rigid.pdbqt → {pdbqt_path}
                rigid_output.rename(pdbqt_path)
                logger.info("P9: Meeko mk_prepare_receptor successful: %s", pdbqt_path)
                return True
            if result.returncode != 0:
                logger.warning("P9: Meeko mk_prepare_receptor failed (rc=%d): %s",
                               result.returncode, result.stderr[:300])
        except Exception as exc:
            logger.warning("P9: Meeko error: %s", exc)
        return False

    # Try CLI executable first
    exe = None
    if shutil.which("mk_prepare_receptor.py"):
        exe = "mk_prepare_receptor.py"
    elif shutil.which("mk_prepare_receptor"):
        exe = "mk_prepare_receptor"

    if exe:
        cmd = [exe, "-i", str(pdb_path), "-o", basename, "-p"]
        if _try_meeko_cmd(cmd):
            return True

    # Fallback: try as Python module
    cmd = ["python", "-m", "meeko.cli.mk_prepare_receptor",
           "-i", str(pdb_path), "-o", basename, "-p"]
    return _try_meeko_cmd(cmd)


def _obabel_pdb_to_pdbqt(pdb_path: Path, pdbqt_path: Path) -> bool:
    """Use Open Babel to convert PDB to PDBQT."""
    try:
        cmd = ["obabel", str(pdb_path), "-O", str(pdbqt_path), "-xr"]
        logger.info("Running: %s", " ".join(cmd))
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode == 0 and pdbqt_path.exists() and pdbqt_path.stat().st_size > 100:
            logger.info("P9: obabel conversion successful: %s", pdbqt_path)
            return True
        logger.warning("P9: obabel failed (rc=%d): %s", result.returncode, result.stderr[:300])
        return False
    except Exception as exc:
        logger.warning("P9: obabel error: %s", exc)
        return False


def _adfr_prepare_receptor(pdb_path: Path, pdbqt_path: Path) -> bool:
    """Use ADFRsuite prepare_receptor script."""
    try:
        cmd = ["prepare_receptor", "-r", str(pdb_path), "-o", str(pdbqt_path)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode == 0 and pdbqt_path.exists():
            return True
        return False
    except Exception:
        return False


def _manual_pdb_to_pdbqt(pdb_path: Path, pdbqt_path: Path) -> None:
    """Manual conversion: add ROOT/ENDROOT and copy ATOM lines.

    Simplified conversion without partial charges or proper atom types.
    Last-resort fallback when no conversion tool is available.
    """
    lines_out: list[str] = []
    lines_out.append("ROOT\n")

    with open(pdb_path, "r") as fh:
        for line in fh:
            if line.startswith("ATOM"):
                padded = line.rstrip("\n").ljust(78)
                element = padded[76:78].strip()
                if not element:
                    atom_name = padded[12:16].strip()
                    element = atom_name[0] if atom_name else "C"
                pdbqt_line = padded[:54] + "  0.00  0.00" + f"    +0.000 {element:>2s}\n"
                lines_out.append(pdbqt_line)

    lines_out.append("ENDROOT\n")
    pdbqt_path.write_text("".join(lines_out))
    logger.info("P9: manual PDBQT written: %s (%d atom lines)", pdbqt_path, len(lines_out) - 2)


# ---------------------------------------------------------------------------
# Main entry point — 7-step receptor preparation pipeline
# ---------------------------------------------------------------------------

def prepare_receptor(
    pdb_path: Path,
    work_dir: Path,
    *,
    structure_source: str = "pdb",
    keep_metals: bool = True,
    keep_cofactors: bool = True,
    minimize: bool | None = None,
    ph: float = 7.4,
) -> tuple[Path, dict]:
    """Prepare a receptor for docking via a 7-step state-of-the-art pipeline.

    Steps: P2 smart cleaning → P4 altloc → P3 non-standard residues →
    P5 missing atoms → P6+P7 protonation → P8 minimisation → P9 PDBQT.

    Each step is optional: if its dependency is missing, it is skipped
    with a warning log. The pipeline remains functional in degraded mode.

    Parameters
    ----------
    pdb_path : Path
        Input PDB file.
    work_dir : Path
        Output directory.
    structure_source : str
        "pdb", "alphafold", or "esmfold". Controls minimisation default.
    keep_metals : bool
        Keep catalytic metal ions (ZN, MG, FE, etc.).
    keep_cofactors : bool
        Keep cofactors (NAD, FAD, HEM, ATP, etc.).
    minimize : bool or None
        Force-enable/disable minimisation. None = auto (True for AF/ESMFold).
    ph : float
        Target pH for protonation (default 7.4).

    Returns
    -------
    tuple[Path, dict]
        (pdbqt_path, prep_report) — report contains details of each step.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    pdbqt_path = work_dir / f"{pdb_path.stem}_receptor.pdbqt"
    prepared_pdb_path = work_dir / f"{pdb_path.stem}_prepared.pdb"

    prep_report: dict = {
        "steps_completed": [],
        "steps_skipped": [],
        "nonstandard_residues": [],
        "missing_atoms_added": 0,
        "shifted_pka_residues": [],
        "minimization_converged": False,
        "conversion_method": "none",
        "structure_source": structure_source,
    }

    # Check cache
    if pdbqt_path.exists() and pdbqt_path.stat().st_size > 100:
        logger.info("Receptor PDBQT already exists: %s", pdbqt_path)
        prep_report["steps_completed"].append("cache_hit")
        return pdbqt_path, prep_report

    logger.info("Starting receptor preparation pipeline for %s (source=%s)",
                pdb_path.name, structure_source)

    current_pdb = pdb_path

    # --- P2: Intelligent HETATM cleaning ---
    current_pdb, clean_report = _clean_pdb(
        current_pdb, work_dir,
        keep_metals=keep_metals,
        keep_cofactors=keep_cofactors,
    )
    prep_report["steps_completed"].append("P2_cleaning")
    prep_report["cleaning"] = clean_report

    # --- P4: Resolve alternate conformations ---
    current_pdb, n_altloc = _resolve_altloc(current_pdb, work_dir)
    if n_altloc > 0:
        prep_report["steps_completed"].append("P4_altloc")
        prep_report["altloc_resolved"] = n_altloc
    else:
        prep_report["steps_skipped"].append(("P4_altloc", "no alternate conformations"))

    # --- P3 + P5: Non-standard residues + missing atoms (PDBFixer) ---
    current_pdb, fix_report = _fix_with_pdbfixer(current_pdb, work_dir)
    if fix_report.get("nonstandard_replaced") or fix_report.get("missing_atoms_added"):
        prep_report["steps_completed"].append("P3_nonstandard")
        prep_report["steps_completed"].append("P5_missing_atoms")
        prep_report["nonstandard_residues"] = fix_report.get("nonstandard_replaced", [])
        prep_report["missing_atoms_added"] = fix_report.get("missing_atoms_added", 0)
    else:
        pdbfixer_ok = fix_report.get("available", False)
        if not fix_report.get("nonstandard_replaced"):
            reason = "none found" if pdbfixer_ok else "PDBFixer not installed"
            prep_report["steps_skipped"].append(("P3_nonstandard", reason))
        if not fix_report.get("missing_atoms_added"):
            reason = "none found" if pdbfixer_ok else "PDBFixer not installed"
            prep_report["steps_skipped"].append(("P5_missing_atoms", reason))

    # --- P6 + P7: Protonation analysis + hydrogens ---
    current_pdb, prot_report = _protonate_receptor(current_pdb, work_dir, ph=ph)
    if prot_report.get("hydrogens_added"):
        prep_report["steps_completed"].append("P6_pka_analysis")
        prep_report["steps_completed"].append("P7_hydrogens")
    else:
        prep_report["steps_skipped"].append(("P6+P7_protonation", "PDBFixer/PropKa unavailable"))
    prep_report["shifted_pka_residues"] = prot_report.get("shifted_pka_residues", [])

    # --- P8: Minimisation (conditional) ---
    should_minimize = minimize
    if should_minimize is None:
        # Auto: minimise predicted structures, skip for experimental
        should_minimize = structure_source in ("alphafold", "esmfold")

    if should_minimize:
        current_pdb, min_report = _minimize_receptor(current_pdb, work_dir)
        prep_report["minimization_converged"] = min_report.get("minimization_converged", False)
        if min_report.get("minimization_converged"):
            prep_report["steps_completed"].append("P8_minimization")
        else:
            prep_report["steps_skipped"].append(("P8_minimization", min_report.get("method", "failed")))
    else:
        prep_report["steps_skipped"].append(("P8_minimization", f"disabled (source={structure_source})"))

    # Save the prepared PDB (useful for GPU docking which takes PDB directly)
    import shutil as _shutil
    _shutil.copy2(current_pdb, prepared_pdb_path)

    # --- P9: PDBQT conversion ---
    # Priority: Meeko > obabel > ADFRsuite > manual
    converted = False

    if not converted:
        converted = _meeko_pdb_to_pdbqt(current_pdb, pdbqt_path)
        if converted:
            prep_report["conversion_method"] = "meeko"

    if not converted and shutil.which("obabel"):
        converted = _obabel_pdb_to_pdbqt(current_pdb, pdbqt_path)
        if converted:
            prep_report["conversion_method"] = "obabel"

    if not converted and shutil.which("prepare_receptor"):
        converted = _adfr_prepare_receptor(current_pdb, pdbqt_path)
        if converted:
            prep_report["conversion_method"] = "adfrsuite"

    if not converted:
        logger.info("P9: no conversion tools available — using manual PDB→PDBQT")
        _manual_pdb_to_pdbqt(current_pdb, pdbqt_path)
        prep_report["conversion_method"] = "manual"

    prep_report["steps_completed"].append("P9_conversion")

    logger.info(
        "Receptor preparation complete: %d steps done, %d skipped, conversion=%s",
        len(prep_report["steps_completed"]),
        len(prep_report["steps_skipped"]),
        prep_report["conversion_method"],
    )

    return pdbqt_path, prep_report


# ---------------------------------------------------------------------------
# Ligand standardisation helpers (L1–L3, L8)
# ---------------------------------------------------------------------------

MAX_ROTATABLE_BONDS = 15


def _standardize_mol(mol):
    """L1 — Standardize: largest fragment, normalize, uncharge."""
    from rdkit.Chem.MolStandardize import rdMolStandardize

    mol = rdMolStandardize.FragmentParent(mol)
    normalizer = rdMolStandardize.Normalizer()
    mol = normalizer.normalize(mol)
    uncharger = rdMolStandardize.Uncharger()
    mol = uncharger.uncharge(mol)
    return mol


def _best_tautomer(mol):
    """L2 — Return the canonical (highest-scored) tautomer."""
    from rdkit.Chem.MolStandardize import rdMolStandardize

    enumerator = rdMolStandardize.TautomerEnumerator()
    return enumerator.Canonicalize(mol)


def _protonate_at_ph(smiles: str, ph: float = 7.4) -> str:
    """L3 — Return SMILES protonated at target pH via Dimorphite-DL."""
    from dimorphite_dl.protonate.run import Protonate

    protonator = Protonate(
        smiles_input=[smiles],
        ph_min=ph, ph_max=ph, precision=0.5,
        max_variants=1,
    )
    protonated = list(protonator)
    if not protonated:
        raise ValueError(f"Dimorphite-DL returned no protonated form for: {smiles[:80]}")
    return protonated[0]


# ---------------------------------------------------------------------------
# Ligand preparation
# ---------------------------------------------------------------------------

def prepare_ligand(smiles: str, name: str, work_dir: Path) -> tuple[Optional[Path], dict]:
    """Convert a SMILES string to a 3D PDBQT file for docking.

    Parameters
    ----------
    smiles : str
        Canonical SMILES of the ligand.
    name : str
        Ligand identifier (used for file naming).
    work_dir : Path
        Output directory.

    Returns
    -------
    tuple[Optional[Path], dict]
        (Path to PDBQT file or None, preparation metadata dict).
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    # Sanitize name for filenames
    safe_name = "".join(c if c.isalnum() or c in "-_" else "_" for c in name)[:60]
    pdbqt_path = work_dir / f"{safe_name}.pdbqt"

    if pdbqt_path.exists() and pdbqt_path.stat().st_size > 50:
        return pdbqt_path, {"input_smiles": smiles, "cached": True}

    # Step 1: Generate 3D SDF from SMILES using RDKit
    sdf_path, meta = _smiles_to_sdf(smiles, safe_name, work_dir)
    if sdf_path is None:
        logger.warning("Failed to generate SDF for %s (%s)", name, smiles)
        return None, meta

    # Step 2: Convert SDF to PDBQT
    # Try Meeko first
    success = _meeko_sdf_to_pdbqt(sdf_path, pdbqt_path)
    if success:
        return pdbqt_path, meta

    # Try Open Babel
    success = _obabel_sdf_to_pdbqt(sdf_path, pdbqt_path)
    if success:
        return pdbqt_path, meta

    # Manual fallback: SDF -> PDB-like -> PDBQT
    success = _manual_sdf_to_pdbqt(sdf_path, pdbqt_path, smiles)
    if success:
        return pdbqt_path, meta

    logger.warning("All PDBQT conversion methods failed for %s", name)
    meta["failure"] = "pdbqt_conversion"
    return None, meta


def _smiles_to_sdf(smiles: str, name: str, work_dir: Path) -> tuple[Optional[Path], dict]:
    """Generate a 3D SDF file from SMILES using RDKit.

    Pipeline: parse → L1 standardize → L2 tautomer → L3 protonate pH 7.4
              → L8 rotatable bonds filter → AddHs → ETKDGv3 → MMFF94 → SDF

    Returns
    -------
    tuple[Optional[Path], dict]
        (sdf_path or None, preparation metadata dict)
    """
    meta: dict = {"input_smiles": smiles, "standardized": False,
                  "tautomer_changed": False, "protonated": False,
                  "minimization": "none", "rotatable_bonds": None}
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors

        # --- Parse raw SMILES ---
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning("RDKit could not parse SMILES: %s", smiles)
            meta["failure"] = "parse"
            return None, meta

        input_canonical = Chem.MolToSmiles(mol)

        # --- L1: Standardize (fragment parent + normalize + uncharge) ---
        mol = _standardize_mol(mol)
        post_std = Chem.MolToSmiles(mol)
        meta["standardized"] = (post_std != input_canonical)

        # --- L2: Canonical tautomer ---
        pre_taut = post_std
        mol = _best_tautomer(mol)
        post_taut = Chem.MolToSmiles(mol)
        meta["tautomer_changed"] = (post_taut != pre_taut)

        # --- L3: Protonate at pH 7.4 (requires SMILES round-trip) ---
        clean_smiles = post_taut
        try:
            protonated_smiles = _protonate_at_ph(clean_smiles)
            mol = Chem.MolFromSmiles(protonated_smiles)
            if mol is None:
                logger.warning("Re-parse after protonation failed for %s", name)
                mol = Chem.MolFromSmiles(clean_smiles)
                if mol is None:
                    meta["failure"] = "protonation"
                    return None, meta
            else:
                meta["protonated"] = True
        except Exception:
            mol = Chem.MolFromSmiles(clean_smiles)
            if mol is None:
                meta["failure"] = "protonation"
                return None, meta

        prepared_smiles = Chem.MolToSmiles(mol)
        meta["prepared_smiles"] = prepared_smiles

        # --- L8: Reject molecules with too many rotatable bonds ---
        rot_bonds = Descriptors.NumRotatableBonds(mol)
        meta["rotatable_bonds"] = rot_bonds
        if rot_bonds > MAX_ROTATABLE_BONDS:
            logger.warning(
                "Skipping %s: %d rotatable bonds (max %d)",
                name, rot_bonds, MAX_ROTATABLE_BONDS,
            )
            meta["failure"] = "rotatable_bonds"
            return None, meta

        # --- L4: 3D conformer (ETKDGv3) ---
        mol = Chem.AddHs(mol)
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        if result == -1:
            result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            if result == -1:
                params = AllChem.ETKDGv3()
                params.useRandomCoords = True
                result = AllChem.EmbedMolecule(mol, params)
                if result == -1:
                    logger.warning("EmbedMolecule failed for %s", name)
                    meta["failure"] = "conformer"
                    return None, meta

        # --- L5: Minimize (MMFF94, fallback UFF) ---
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            meta["minimization"] = "MMFF94"
        except Exception:
            try:
                AllChem.UFFOptimizeMolecule(mol, maxIters=200)
                meta["minimization"] = "UFF"
            except Exception:
                meta["minimization"] = "unoptimized"

        sdf_path = work_dir / f"{name}.sdf"
        writer = Chem.SDWriter(str(sdf_path))
        writer.write(mol)
        writer.close()

        if sdf_path.exists() and sdf_path.stat().st_size > 10:
            return sdf_path, meta
        meta["failure"] = "sdf_write"
        return None, meta

    except ImportError:
        logger.warning("RDKit not available; cannot generate 3D coordinates from SMILES")
        meta["failure"] = "no_rdkit"
        return None, meta
    except Exception as exc:
        logger.warning("SDF generation error for %s: %s", name, exc)
        meta["failure"] = f"error:{exc}"
        return None, meta


def _meeko_sdf_to_pdbqt(sdf_path: Path, pdbqt_path: Path) -> bool:
    """Convert SDF to PDBQT using Meeko."""
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
        from rdkit import Chem

        mol = Chem.SDMolSupplier(str(sdf_path), removeHs=False)[0]
        if mol is None:
            return False

        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(mol)

        for setup in mol_setups:
            pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
            if is_ok:
                pdbqt_path.write_text(pdbqt_string)
                logger.info("Meeko conversion successful: %s", pdbqt_path)
                return True
            else:
                logger.warning("Meeko write error: %s", error_msg)

        return False

    except ImportError:
        logger.debug("Meeko not available")
        return False
    except Exception as exc:
        logger.warning("Meeko error: %s", exc)
        return False


def _obabel_sdf_to_pdbqt(sdf_path: Path, pdbqt_path: Path) -> bool:
    """Convert SDF to PDBQT using Open Babel."""
    if not shutil.which("obabel"):
        return False
    try:
        cmd = ["obabel", str(sdf_path), "-O", str(pdbqt_path), "--gen3d"]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if result.returncode == 0 and pdbqt_path.exists() and pdbqt_path.stat().st_size > 50:
            logger.info("obabel SDF->PDBQT successful: %s", pdbqt_path)
            return True
        return False
    except Exception as exc:
        logger.warning("obabel SDF->PDBQT error: %s", exc)
        return False


def _manual_sdf_to_pdbqt(sdf_path: Path, pdbqt_path: Path, smiles: str) -> bool:
    """Manual SDF to PDBQT conversion (minimal/approximate)."""
    try:
        from rdkit import Chem

        mol = Chem.SDMolSupplier(str(sdf_path), removeHs=False)[0]
        if mol is None:
            return False

        conf = mol.GetConformer()
        lines: list[str] = []
        lines.append("ROOT\n")

        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            element = atom.GetSymbol()
            atom_name = f"{element}{i+1}"[:4].ljust(4)
            # PDB-like ATOM line
            line = (
                f"ATOM  {i+1:5d} {atom_name} LIG A   1    "
                f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}"
                f"  1.00  0.00    +0.000 {element:>2s}\n"
            )
            lines.append(line)

        lines.append("ENDROOT\n")
        lines.append("TORSDOF 0\n")
        pdbqt_path.write_text("".join(lines))
        logger.info("Manual SDF->PDBQT written: %s", pdbqt_path)
        return True

    except Exception as exc:
        logger.warning("Manual SDF->PDBQT failed: %s", exc)
        return False
