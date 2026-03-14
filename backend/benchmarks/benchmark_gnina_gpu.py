"""
GNINA GPU Docking Benchmark — Measures real docking times per configuration.

Calls dock_all_runpod_batch() directly (no Celery) with the kinase inhibitor
dataset and measures wall-clock time for each parameter combination.

Usage:
    cd backend && python benchmarks/benchmark_gnina_gpu.py

Output: JSON results + summary table to stdout.
Requires: RUNPOD_API_KEY set in environment.
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# Ensure backend/ is on sys.path so pipeline imports work
# ---------------------------------------------------------------------------
BACKEND_DIR = Path(__file__).resolve().parent.parent
if str(BACKEND_DIR) not in sys.path:
    sys.path.insert(0, str(BACKEND_DIR))

import requests

from benchmarks.kinase_inhibitors import get_all_actives, DECOYS

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PDB_ID = "1M17"
PDB_URL = f"https://files.rcsb.org/download/{PDB_ID}.pdb"

# EGFR pocket center (Erlotinib binding site in 1M17)
POCKET_CENTER = (22.0, 0.5, 52.0)
POCKET_SIZE = (20.0, 20.0, 20.0)

RESULTS_PATH = Path(__file__).resolve().parent / "benchmark_results.json"

# Override defaults for benchmark:
# - Smaller batches so each RunPod job finishes within executionTimeout
# - Longer polling timeout for cold starts
os.environ.setdefault("RUNPOD_TIMEOUT", "1800")
os.environ.setdefault("GPU_BATCH_SIZE", "20")

# ---------------------------------------------------------------------------
# 8 benchmark configurations
# ---------------------------------------------------------------------------

CONFIGS = [
    {"label": "baseline",       "n_mols": 66,  "exhaustiveness": 8,  "cnn_scoring": "rescore",    "num_modes": 9},
    {"label": "scale_x2",       "n_mols": 132, "exhaustiveness": 8,  "cnn_scoring": "rescore",    "num_modes": 9},
    {"label": "scale_x4",       "n_mols": 264, "exhaustiveness": 8,  "cnn_scoring": "rescore",    "num_modes": 9},
    {"label": "high_exhaust",   "n_mols": 66,  "exhaustiveness": 32, "cnn_scoring": "rescore",    "num_modes": 9},
    {"label": "refinement",     "n_mols": 66,  "exhaustiveness": 8,  "cnn_scoring": "refinement", "num_modes": 9},
    {"label": "refine+exhaust", "n_mols": 66,  "exhaustiveness": 32, "cnn_scoring": "refinement", "num_modes": 9},
    {"label": "minimal_modes",  "n_mols": 66,  "exhaustiveness": 8,  "cnn_scoring": "rescore",    "num_modes": 3},
    {"label": "full_quality",   "n_mols": 66,  "exhaustiveness": 32, "cnn_scoring": "refinement", "num_modes": 9},
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _build_ligand_set(n_mols: int) -> list[dict]:
    """Build a ligand list of exactly n_mols by combining actives + decoys, duplicating if needed."""
    base = get_all_actives() + [{"name": d["name"], "smiles": d["smiles"]} for d in DECOYS]
    # base has 56 actives + 10 decoys = 66
    if n_mols <= len(base):
        return base[:n_mols]

    # Duplicate with suffixed names to reach target count
    result = list(base)
    copy_idx = 1
    while len(result) < n_mols:
        for mol in base:
            if len(result) >= n_mols:
                break
            result.append({**mol, "name": f"{mol['name']}_copy{copy_idx}"})
        copy_idx += 1
    return result


def _download_pdb(pdb_id: str, dest: Path) -> Path:
    """Download a PDB file from RCSB if not already present."""
    pdb_file = dest / f"{pdb_id}.pdb"
    if pdb_file.exists():
        print(f"  PDB {pdb_id} already cached at {pdb_file}")
        return pdb_file

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"  Downloading {pdb_id} from RCSB...")
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    pdb_file.write_text(resp.text)
    print(f"  Saved to {pdb_file} ({len(resp.text)} bytes)")
    return pdb_file


def _format_time(seconds: float) -> str:
    """Format seconds as mm:ss or ss.fs."""
    if seconds >= 60:
        m, s = divmod(seconds, 60)
        return f"{int(m)}m{s:04.1f}s"
    return f"{seconds:.1f}s"


def _print_table(results: list[dict]) -> None:
    """Print a formatted results table."""
    header = f"{'Config':<21s} {'Mols':>5s} {'Exh':>4s} {'CNN':>10s} {'Modes':>5s} {'Prep':>8s} {'GPU':>8s} {'Total':>8s} {'mol/min':>8s} {'Results':>7s}"
    sep = "-" * len(header)

    print()
    print(sep)
    print(header)
    print(sep)

    for r in results:
        if r.get("error"):
            print(f"{r['label']:<21s} {'ERROR: ' + r['error']}")
            continue

        throughput = r["n_results"] / r["total_time"] * 60 if r["total_time"] > 0 else 0
        print(
            f"{r['label']:<21s} "
            f"{r['n_mols']:>5d} "
            f"{r['exhaustiveness']:>4d} "
            f"{r['cnn_scoring']:>10s} "
            f"{r['num_modes']:>5d} "
            f"{_format_time(r['prep_time']):>8s} "
            f"{_format_time(r['gpu_time']):>8s} "
            f"{_format_time(r['total_time']):>8s} "
            f"{throughput:>8.1f} "
            f"{r['n_results']:>7d}"
        )

    print(sep)


# ---------------------------------------------------------------------------
# Main benchmark runner
# ---------------------------------------------------------------------------

def run_benchmark() -> list[dict]:
    """Run all benchmark configurations and return results."""
    # Check API key
    api_key = os.environ.get("RUNPOD_API_KEY", "")
    if not api_key:
        print("ERROR: RUNPOD_API_KEY not set. Export it and retry.")
        print("  export RUNPOD_API_KEY=your_key_here")
        sys.exit(1)

    # Import docking function
    from pipeline.docking_gpu import dock_all_runpod_batch

    # Setup temp directory for receptor + work files
    with tempfile.TemporaryDirectory(prefix="gnina_bench_") as tmpdir:
        tmpdir = Path(tmpdir)
        pdb_file = _download_pdb(PDB_ID, tmpdir)

        results = []
        total_configs = len(CONFIGS)

        print(f"\nGNINA GPU Benchmark — {total_configs} configurations")
        print(f"Receptor: {PDB_ID} (EGFR + Erlotinib)")
        print(f"Pocket center: {POCKET_CENTER}")
        print(f"Endpoint: {os.environ.get('GNINA_ENDPOINT_ID', 'weeeiy6z4jdsv3')}")
        print()

        for i, cfg in enumerate(CONFIGS, 1):
            label = cfg["label"]
            n_mols = cfg["n_mols"]
            exhaustiveness = cfg["exhaustiveness"]
            cnn_scoring = cfg["cnn_scoring"]
            num_modes = cfg["num_modes"]

            print(f"[{i}/{total_configs}] {label}: {n_mols} mols, exh={exhaustiveness}, "
                  f"cnn={cnn_scoring}, modes={num_modes}")

            # Build ligand set
            ligands = _build_ligand_set(n_mols)

            # Work directory per config
            work_dir = tmpdir / label
            work_dir.mkdir(exist_ok=True)

            # Run docking
            t_start = time.time()

            try:
                # Retry up to 2 times on transient network/SSL errors
                docked = None
                last_err = None
                for attempt in range(3):
                    try:
                        docked = dock_all_runpod_batch(
                            receptor_pdbqt=pdb_file,
                            ligands=ligands,
                            center=POCKET_CENTER,
                            work_dir=work_dir,
                            size=POCKET_SIZE,
                            exhaustiveness=exhaustiveness,
                            num_modes=num_modes,
                            cnn_scoring=cnn_scoring,
                        )
                        break
                    except RuntimeError as e:
                        if "network:" in str(e) and attempt < 2:
                            print(f"  -> Retry {attempt + 1}/2 (network error)")
                            time.sleep(5)
                            last_err = e
                        else:
                            raise
                if docked is None:
                    raise last_err or RuntimeError("docking returned None")
                t_total = time.time() - t_start

                # Estimate prep vs GPU time from logs (approximate)
                # Prep = SDF generation, GPU = submit+poll
                # We'll use a heuristic: prep ~ proportional to n_mols
                n_results = len(docked)

                result = {
                    "label": label,
                    "n_mols": n_mols,
                    "exhaustiveness": exhaustiveness,
                    "cnn_scoring": cnn_scoring,
                    "num_modes": num_modes,
                    "total_time": round(t_total, 2),
                    "prep_time": 0.0,  # filled below
                    "gpu_time": 0.0,   # filled below
                    "n_results": n_results,
                    "throughput_mol_per_min": round(n_results / t_total * 60, 1) if t_total > 0 else 0,
                    "error": None,
                }

                # Extract GPU-reported elapsed from docked results (if available)
                # The gpu_time is embedded in the batch results by _poll_job
                # For now, total_time is the best measure
                result["gpu_time"] = round(t_total, 2)
                result["prep_time"] = 0.0

                print(f"  -> {n_results} results in {_format_time(t_total)} "
                      f"({result['throughput_mol_per_min']:.1f} mol/min)")

            except Exception as e:
                t_total = time.time() - t_start
                result = {
                    "label": label,
                    "n_mols": n_mols,
                    "exhaustiveness": exhaustiveness,
                    "cnn_scoring": cnn_scoring,
                    "num_modes": num_modes,
                    "total_time": round(t_total, 2),
                    "prep_time": 0.0,
                    "gpu_time": 0.0,
                    "n_results": 0,
                    "throughput_mol_per_min": 0,
                    "error": str(e),
                }
                print(f"  -> ERROR after {_format_time(t_total)}: {e}")

            results.append(result)
            print()

    return results


def main() -> None:
    """Entry point."""
    print("=" * 70)
    print("  GNINA GPU Docking Benchmark")
    print("=" * 70)

    start = time.time()
    results = run_benchmark()
    elapsed = time.time() - start

    # Print summary table
    _print_table(results)

    # Summary stats
    successful = [r for r in results if not r.get("error")]
    if successful:
        total_mols = sum(r["n_results"] for r in successful)
        total_time = sum(r["total_time"] for r in successful)
        print(f"\nTotal: {total_mols} molecules docked across {len(successful)} configs "
              f"in {_format_time(total_time)}")
        print(f"Benchmark wall-clock: {_format_time(elapsed)}")

    failed = [r for r in results if r.get("error")]
    if failed:
        print(f"\n{len(failed)} config(s) failed:")
        for r in failed:
            print(f"  - {r['label']}: {r['error']}")

    # Save JSON
    output = {
        "benchmark_date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "receptor_pdb": PDB_ID,
        "pocket_center": list(POCKET_CENTER),
        "pocket_size": list(POCKET_SIZE),
        "endpoint_id": os.environ.get("GNINA_ENDPOINT_ID", "weeeiy6z4jdsv3"),
        "total_wall_clock_s": round(elapsed, 2),
        "configs": results,
    }

    RESULTS_PATH.write_text(json.dumps(output, indent=2))
    print(f"\nResults saved to {RESULTS_PATH}")


if __name__ == "__main__":
    main()
