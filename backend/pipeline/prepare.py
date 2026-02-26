"""
DockIt pipeline — Receptor and ligand preparation.

Converts PDB -> PDBQT (receptor) and SMILES -> PDBQT (ligand) using
Open Babel and Meeko, with fallbacks when those tools are unavailable.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Receptor preparation
# ---------------------------------------------------------------------------

def _strip_hetatm(pdb_path: Path, work_dir: Path) -> Path:
    """Strip HETATM records from a PDB file for receptor preparation.

    Removes waters (HOH), co-crystallized ligands, and crystallization
    artifacts that would occupy the binding pocket and block docking.
    Keeps only protein ATOM records, CONECT, TER, END, and header lines.

    Parameters
    ----------
    pdb_path : Path
        Original PDB file (possibly from AlphaFold or PDB).
    work_dir : Path
        Directory for the cleaned PDB.

    Returns
    -------
    Path
        Path to the cleaned PDB file (protein-only).
    """
    cleaned_path = work_dir / f"{pdb_path.stem}_clean.pdb"
    kept = 0
    removed = 0
    removed_residues: set[str] = set()

    with open(pdb_path, "r") as fin, open(cleaned_path, "w") as fout:
        for line in fin:
            record = line[:6].strip()
            if record == "HETATM":
                res_name = line[17:20].strip()
                removed_residues.add(res_name)
                removed += 1
                continue
            fout.write(line)
            if record == "ATOM":
                kept += 1

    if removed > 0:
        logger.info(
            "Stripped %d HETATM atoms (%s) from %s, kept %d ATOM records",
            removed, ", ".join(sorted(removed_residues)), pdb_path.name, kept,
        )
    else:
        logger.debug("No HETATM records found in %s", pdb_path.name)

    return cleaned_path


def prepare_receptor(pdb_path: Path, work_dir: Path) -> Path:
    """Convert a PDB file to PDBQT format for AutoDock Vina.

    The PDB is first cleaned by stripping all HETATM records (waters,
    co-crystallized ligands, crystallization artifacts) that would block
    the binding pocket during docking.

    Parameters
    ----------
    pdb_path : Path
        Input PDB file.
    work_dir : Path
        Output directory.

    Returns
    -------
    Path
        Path to the generated ``.pdbqt`` file.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    pdbqt_path = work_dir / f"{pdb_path.stem}_receptor.pdbqt"

    if pdbqt_path.exists() and pdbqt_path.stat().st_size > 100:
        logger.info("Receptor PDBQT already exists: %s", pdbqt_path)
        return pdbqt_path

    # Step 0: Strip HETATM records (waters, ligands, artifacts)
    cleaned_pdb = _strip_hetatm(pdb_path, work_dir)

    # Strategy 1: Open Babel
    if shutil.which("obabel"):
        success = _obabel_pdb_to_pdbqt(cleaned_pdb, pdbqt_path)
        if success:
            return pdbqt_path

    # Strategy 2: ADFRsuite prepare_receptor
    if shutil.which("prepare_receptor"):
        success = _adfr_prepare_receptor(cleaned_pdb, pdbqt_path)
        if success:
            return pdbqt_path

    # Strategy 3: Manual fallback (no external tools)
    logger.info("No conversion tools available; using manual PDB->PDBQT conversion")
    _manual_pdb_to_pdbqt(cleaned_pdb, pdbqt_path)
    return pdbqt_path


def _obabel_pdb_to_pdbqt(pdb_path: Path, pdbqt_path: Path) -> bool:
    """Use Open Babel to convert PDB to PDBQT."""
    try:
        cmd = ["obabel", str(pdb_path), "-O", str(pdbqt_path), "-xr"]
        logger.info("Running: %s", " ".join(cmd))
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode == 0 and pdbqt_path.exists() and pdbqt_path.stat().st_size > 100:
            logger.info("obabel conversion successful: %s", pdbqt_path)
            return True
        logger.warning("obabel failed (rc=%d): %s", result.returncode, result.stderr[:300])
        return False
    except Exception as exc:
        logger.warning("obabel error: %s", exc)
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

    This is a simplified conversion that does NOT add partial charges
    or atom types, but produces a file that Vina can at least attempt
    to read. For production use, a proper tool chain is required.
    """
    lines_out: list[str] = []
    lines_out.append("ROOT\n")

    with open(pdb_path, "r") as fh:
        for line in fh:
            if line.startswith("ATOM"):
                # Ensure line is at least 78 chars (PDB fixed-width format)
                padded = line.rstrip("\n").ljust(78)
                # Append autodock atom type (use element symbol from cols 76-78)
                element = padded[76:78].strip()
                if not element:
                    # Guess from atom name (cols 12-16)
                    atom_name = padded[12:16].strip()
                    element = atom_name[0] if atom_name else "C"
                # Build PDBQT line: add charge (0.000) and AD type
                pdbqt_line = padded[:54] + "  0.00  0.00" + f"    +0.000 {element:>2s}\n"
                lines_out.append(pdbqt_line)

    lines_out.append("ENDROOT\n")
    pdbqt_path.write_text("".join(lines_out))
    logger.info("Manual PDBQT written: %s (%d atom lines)", pdbqt_path, len(lines_out) - 2)


# ---------------------------------------------------------------------------
# Ligand preparation
# ---------------------------------------------------------------------------

def prepare_ligand(smiles: str, name: str, work_dir: Path) -> Optional[Path]:
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
    Optional[Path]
        Path to PDBQT file, or ``None`` if preparation fails.
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    # Sanitize name for filenames
    safe_name = "".join(c if c.isalnum() or c in "-_" else "_" for c in name)[:60]
    pdbqt_path = work_dir / f"{safe_name}.pdbqt"

    if pdbqt_path.exists() and pdbqt_path.stat().st_size > 50:
        return pdbqt_path

    # Step 1: Generate 3D SDF from SMILES using RDKit
    sdf_path = _smiles_to_sdf(smiles, safe_name, work_dir)
    if sdf_path is None:
        logger.warning("Failed to generate SDF for %s (%s)", name, smiles)
        return None

    # Step 2: Convert SDF to PDBQT
    # Try Meeko first
    success = _meeko_sdf_to_pdbqt(sdf_path, pdbqt_path)
    if success:
        return pdbqt_path

    # Try Open Babel
    success = _obabel_sdf_to_pdbqt(sdf_path, pdbqt_path)
    if success:
        return pdbqt_path

    # Manual fallback: SDF -> PDB-like -> PDBQT
    success = _manual_sdf_to_pdbqt(sdf_path, pdbqt_path, smiles)
    if success:
        return pdbqt_path

    logger.warning("All PDBQT conversion methods failed for %s", name)
    return None


def _smiles_to_sdf(smiles: str, name: str, work_dir: Path) -> Optional[Path]:
    """Generate a 3D SDF file from SMILES using RDKit."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning("RDKit could not parse SMILES: %s", smiles)
            return None

        mol = Chem.AddHs(mol)
        # Generate 3D coordinates
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        if result == -1:
            # Try with random coordinates
            result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            if result == -1:
                # Last attempt with less strict parameters
                params = AllChem.ETKDGv3()
                params.useRandomCoords = True
                result = AllChem.EmbedMolecule(mol, params)
                if result == -1:
                    logger.warning("EmbedMolecule failed for %s", name)
                    return None

        # Optimize geometry
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            try:
                AllChem.UFFOptimizeMolecule(mol, maxIters=200)
            except Exception:
                pass  # Use unoptimized coords

        sdf_path = work_dir / f"{name}.sdf"
        writer = Chem.SDWriter(str(sdf_path))
        writer.write(mol)
        writer.close()

        if sdf_path.exists() and sdf_path.stat().st_size > 10:
            return sdf_path
        return None

    except ImportError:
        logger.warning("RDKit not available; cannot generate 3D coordinates from SMILES")
        return None
    except Exception as exc:
        logger.warning("SDF generation error for %s: %s", name, exc)
        return None


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
