#!/usr/bin/env python3
"""
Benchmark: Validate the new 7-step receptor preparation pipeline.

Uses well-known PDB entries to test each step:
  1. 1M17 — EGFR kinase + erlotinib (co-crystallised, has waters/buffers)
  2. 3HTB — CDK2 + dinaciclib (has MSE selenomethionine, altloc, metals)
  3. 5WIU — BACE1 (has alternate conformations)

For each protein:
  - Download PDB from RCSB
  - Run the full preparation pipeline
  - Validate that each step produces expected outputs
  - Report step-by-step results

Usage:
    python benchmark_receptor_prep.py [--all]
"""
import json
import sys
import tempfile
import urllib.request
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline.prepare import prepare_receptor


# ---------------------------------------------------------------------------
# Test cases
# ---------------------------------------------------------------------------

TEST_CASES = [
    {
        "pdb_id": "1M17",
        "name": "EGFR + erlotinib",
        "description": "Classic kinase target. Has waters, SO4, ligand to remove.",
        "structure_source": "pdb",
        "expect": {
            "should_clean": True,          # Waters + SO4 + ligand
            "should_have_pdbqt": True,
            "minimization": "skipped",     # PDB structure → no minimization
        },
    },
    {
        "pdb_id": "3HTB",
        "name": "CDK2 + dinaciclib",
        "description": "Has MSE (selenomethionine) non-standard residues.",
        "structure_source": "pdb",
        "expect": {
            "should_clean": True,
            "should_have_pdbqt": True,
            "minimization": "skipped",
        },
    },
    {
        "pdb_id": "5WIU",
        "name": "BACE1",
        "description": "Structure with alternate conformations to test altloc resolution.",
        "structure_source": "pdb",
        "expect": {
            "should_clean": True,
            "should_have_pdbqt": True,
            "minimization": "skipped",
        },
    },
]


def download_pdb(pdb_id: str, output_dir: Path) -> Path:
    """Download PDB from RCSB."""
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    output_path = output_dir / f"{pdb_id.upper()}.pdb"

    if output_path.exists():
        print(f"  [cache] {pdb_id} already downloaded")
        return output_path

    print(f"  [download] {url}")
    urllib.request.urlretrieve(url, str(output_path))

    # Validate
    content = output_path.read_text()
    atom_count = content.count("\nATOM ")
    hetatm_count = content.count("\nHETATM")
    print(f"  [raw] {atom_count} ATOM, {hetatm_count} HETATM records")

    return output_path


def validate_pdb_file(path: Path, label: str) -> dict:
    """Count atoms/residues in a PDB file."""
    if not path.exists():
        return {"exists": False, "label": label}

    content = path.read_text()
    atoms = content.count("\nATOM ")
    hetatm = content.count("\nHETATM")
    lines = len(content.splitlines())
    size_kb = path.stat().st_size / 1024

    # Count unique residues
    residues = set()
    for line in content.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            res = line[17:20].strip()
            if res:
                residues.add(res)

    # Check for hydrogens
    h_count = sum(1 for line in content.splitlines()
                  if (line.startswith("ATOM") or line.startswith("HETATM"))
                  and line[76:78].strip() == "H")

    return {
        "exists": True,
        "label": label,
        "atoms": atoms,
        "hetatm": hetatm,
        "lines": lines,
        "size_kb": round(size_kb, 1),
        "unique_residues": sorted(residues),
        "hydrogen_count": h_count,
    }


def validate_pdbqt_file(path: Path) -> dict:
    """Validate PDBQT file format."""
    if not path.exists():
        return {"exists": False, "valid": False}

    content = path.read_text()
    lines = content.splitlines()
    atom_lines = [l for l in lines if l.startswith("ATOM") or l.startswith("HETATM")]

    # PDBQT has AutoDock atom types in columns 78-79
    has_ad_types = False
    for line in atom_lines[:10]:
        if len(line) >= 78:
            ad_type = line[77:79].strip()
            if ad_type:
                has_ad_types = True
                break

    return {
        "exists": True,
        "valid": len(atom_lines) > 0,
        "atom_count": len(atom_lines),
        "has_autodock_types": has_ad_types,
        "size_kb": round(path.stat().st_size / 1024, 1),
    }


def run_benchmark(case: dict, base_dir: Path) -> dict:
    """Run the preparation pipeline for one test case."""
    pdb_id = case["pdb_id"]
    print(f"\n{'='*60}")
    print(f"BENCHMARK: {pdb_id} — {case['name']}")
    print(f"  {case['description']}")
    print(f"{'='*60}")

    result = {
        "pdb_id": pdb_id,
        "name": case["name"],
        "status": "pending",
        "steps": {},
        "errors": [],
    }

    work_dir = base_dir / pdb_id
    work_dir.mkdir(parents=True, exist_ok=True)

    # 1. Download
    try:
        pdb_path = download_pdb(pdb_id, work_dir)
        raw_info = validate_pdb_file(pdb_path, "raw")
        result["raw_pdb"] = raw_info
        print(f"  [raw] {raw_info['atoms']} atoms, {raw_info['hetatm']} hetatm, "
              f"{len(raw_info['unique_residues'])} unique residue types")
    except Exception as e:
        result["status"] = "download_failed"
        result["errors"].append(str(e))
        print(f"  [ERROR] Download failed: {e}")
        return result

    # 2. Run preparation pipeline
    print(f"\n  --- Running 7-step preparation pipeline ---")
    try:
        pdbqt_path, prep_report = prepare_receptor(
            pdb_path, work_dir,
            structure_source=case["structure_source"],
            keep_metals=True,
            keep_cofactors=True,
            minimize=None,   # auto (should be False for PDB sources)
            ph=7.4,
        )
        result["prep_report"] = prep_report
        result["status"] = "completed"
    except Exception as e:
        result["status"] = "pipeline_failed"
        result["errors"].append(str(e))
        print(f"  [ERROR] Pipeline failed: {e}")
        import traceback
        traceback.print_exc()
        return result

    # 3. Report step-by-step results
    steps_done = prep_report.get("steps_completed", [])
    steps_skip = prep_report.get("steps_skipped", [])

    print(f"\n  --- Pipeline Results ---")
    print(f"  Steps completed ({len(steps_done)}):")
    for s in steps_done:
        print(f"    ✓ {s}")
    print(f"  Steps skipped ({len(steps_skip)}):")
    for s in steps_skip:
        if isinstance(s, (list, tuple)):
            print(f"    ○ {s[0]}: {s[1]}")
        else:
            print(f"    ○ {s}")

    # 4. Validate cleaning
    cleaning = prep_report.get("cleaning", {})
    if cleaning:
        print(f"\n  --- Cleaning Report ---")
        print(f"    HETATM removed: {cleaning.get('hetatm_removed', 0)}")
        print(f"    Residues removed: {cleaning.get('residues_removed', [])}")
        print(f"    Het kept (metals/cofactors): {cleaning.get('het_kept', [])}")

    # 5. Validate non-standard residues
    nonstd = prep_report.get("nonstandard_residues", [])
    if nonstd:
        print(f"\n  --- Non-standard Residues Replaced ---")
        for r in nonstd:
            print(f"    {r}")

    # 6. Validate protonation
    shifted = prep_report.get("shifted_pka_residues", [])
    if shifted:
        print(f"\n  --- Shifted pKa Residues ---")
        for r in shifted:
            if isinstance(r, dict):
                print(f"    {r['residue']} chain {r.get('chain','?')} "
                      f"pKa={r['pka']} (model={r.get('model_pka','?')}, shift={r.get('shift','?')})")
            else:
                print(f"    {r}")

    # 7. Validate output files
    # Check cleaned PDB
    clean_candidates = list(work_dir.glob("*_clean.pdb"))
    if clean_candidates:
        clean_info = validate_pdb_file(clean_candidates[0], "cleaned")
        result["cleaned_pdb"] = clean_info
        print(f"\n  --- Cleaned PDB ---")
        print(f"    {clean_info['atoms']} atoms, {clean_info['hetatm']} hetatm, "
              f"{clean_info['size_kb']} KB")

        # Verify waters were removed
        has_hoh = "HOH" in clean_info["unique_residues"]
        print(f"    Waters removed: {'NO (FAIL)' if has_hoh else 'YES (OK)'}")
        if has_hoh:
            result["errors"].append("Waters (HOH) still present after cleaning")

    # Check prepared PDB (with hydrogens)
    prepared_candidates = list(work_dir.glob("*_prepared.pdb"))
    if prepared_candidates:
        prep_info = validate_pdb_file(prepared_candidates[0], "prepared")
        result["prepared_pdb"] = prep_info
        print(f"\n  --- Prepared PDB (with H) ---")
        print(f"    {prep_info['atoms']} atoms, {prep_info['hydrogen_count']} hydrogens, "
              f"{prep_info['size_kb']} KB")

        # Compare atom counts (should have more atoms after adding H)
        if raw_info.get("atoms", 0) > 0:
            ratio = prep_info["atoms"] / raw_info["atoms"]
            print(f"    Atom ratio (prepared/raw): {ratio:.2f}x")
            if ratio < 0.9:
                result["errors"].append(f"Too many atoms lost: ratio={ratio:.2f}")

    # Check PDBQT
    pdbqt_info = validate_pdbqt_file(pdbqt_path)
    result["pdbqt"] = pdbqt_info
    print(f"\n  --- PDBQT Output ---")
    if pdbqt_info["exists"]:
        print(f"    Valid: {pdbqt_info['valid']}")
        print(f"    Atoms: {pdbqt_info['atom_count']}")
        print(f"    AutoDock types: {pdbqt_info['has_autodock_types']}")
        print(f"    Size: {pdbqt_info['size_kb']} KB")
        print(f"    Conversion method: {prep_report.get('conversion_method', 'unknown')}")
    else:
        print(f"    [FAIL] PDBQT file not found!")
        result["errors"].append("PDBQT file not generated")

    # 8. Overall validation
    print(f"\n  --- Validation ---")
    expect = case.get("expect", {})

    if expect.get("should_clean") and cleaning.get("hetatm_removed", 0) == 0:
        result["errors"].append("Expected cleaning to remove HETATM but none removed")
        print(f"    [WARN] Expected HETATM removal but got 0")

    if expect.get("should_have_pdbqt") and not pdbqt_info.get("valid"):
        result["errors"].append("Expected valid PDBQT but output is invalid")
        print(f"    [FAIL] PDBQT validation failed")

    if expect.get("minimization") == "skipped" and "P8_minimization" in steps_done:
        result["errors"].append("Minimization should have been skipped for PDB source")
        print(f"    [WARN] Minimization ran but should have been skipped")

    if not result["errors"]:
        print(f"    ✓ ALL CHECKS PASSED")
        result["status"] = "passed"
    else:
        print(f"    ✗ {len(result['errors'])} issue(s)")
        result["status"] = "issues"

    return result


def main():
    print("=" * 60)
    print("BindX V9 — Receptor Preparation Pipeline Benchmark")
    print("=" * 60)

    # Create persistent benchmark directory
    bench_dir = Path(__file__).parent / "receptor_prep_results"
    bench_dir.mkdir(exist_ok=True)

    cases = TEST_CASES
    if "--all" not in sys.argv and len(sys.argv) > 1:
        # Allow filtering by PDB ID
        target = sys.argv[1].upper()
        cases = [c for c in cases if c["pdb_id"] == target]
        if not cases:
            print(f"No test case found for PDB ID: {target}")
            sys.exit(1)

    results = []
    for case in cases:
        result = run_benchmark(case, bench_dir)
        results.append(result)

    # Summary
    print(f"\n{'='*60}")
    print("BENCHMARK SUMMARY")
    print(f"{'='*60}")

    passed = sum(1 for r in results if r["status"] == "passed")
    issues = sum(1 for r in results if r["status"] == "issues")
    failed = sum(1 for r in results if r["status"] in ("pipeline_failed", "download_failed"))

    for r in results:
        icon = "✓" if r["status"] == "passed" else "⚠" if r["status"] == "issues" else "✗"
        steps_done = len(r.get("prep_report", {}).get("steps_completed", []))
        method = r.get("prep_report", {}).get("conversion_method", "?")
        print(f"  {icon} {r['pdb_id']} ({r['name']}) — {r['status']} — "
              f"{steps_done} steps, conv={method}")
        if r["errors"]:
            for e in r["errors"]:
                print(f"      → {e}")

    print(f"\nTotal: {passed} passed, {issues} with issues, {failed} failed")

    # Save detailed results
    report_path = bench_dir / "benchmark_report.json"
    # Serialize, handle tuples
    def _serialize(obj):
        if isinstance(obj, set):
            return sorted(obj)
        if isinstance(obj, tuple):
            return list(obj)
        if isinstance(obj, Path):
            return str(obj)
        raise TypeError(f"Cannot serialize {type(obj)}")

    with open(report_path, "w") as f:
        json.dump(results, f, indent=2, default=_serialize)
    print(f"\nDetailed report: {report_path}")

    sys.exit(0 if failed == 0 else 1)


if __name__ == "__main__":
    main()
