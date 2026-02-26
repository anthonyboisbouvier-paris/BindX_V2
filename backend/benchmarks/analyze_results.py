#!/usr/bin/env python3
"""
DockIt Benchmark Analysis Script

Analyzes benchmark results across 5 kinase targets:
- Enrichment factors (EF1%, EF10%, avg rank ratio)
- Spearman rank correlation between docking scores and experimental IC50
- Score range validation
- GPU vs CPU consistency
- Publication-ready summary tables
"""

import json
import subprocess
import sys
from pathlib import Path
from typing import Optional

import numpy as np
from scipy import stats

try:
    from rdkit import Chem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))
from benchmarks.kinase_inhibitors import TARGETS, DECOYS


def _canonicalize(smiles: str) -> Optional[str]:
    """Return RDKit canonical SMILES, or None if invalid."""
    if not HAS_RDKIT or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol else None


def _build_lookup(compounds: list[dict], is_active: bool) -> dict[str, dict]:
    """Build canonical SMILES -> compound info lookup."""
    lookup = {}
    for c in compounds:
        canon = _canonicalize(c["smiles"])
        if canon:
            lookup[canon] = {
                "name": c["name"],
                "ic50_nM": c.get("ic50_nM") if is_active else None,
                "is_active": is_active,
            }
    return lookup


def fetch_results(job_id: str) -> list[dict]:
    """Fetch results from DockIt API."""
    r = subprocess.run(
        ["curl", "-s", f"http://localhost:8000/api/jobs/{job_id}/results"],
        capture_output=True, text=True,
    )
    data = json.loads(r.stdout)
    return data.get("results", [])


def identify_molecule(mol: dict, active_lookup: dict, decoy_lookup: dict) -> tuple[str, Optional[float], bool, bool]:
    """Identify a molecule as active/decoy using canonical SMILES matching.

    Returns (name, ic50_nM, is_active, is_decoy).
    """
    smiles = mol.get("smiles", "")
    canon = _canonicalize(smiles)

    if canon and canon in active_lookup:
        info = active_lookup[canon]
        return info["name"], info["ic50_nM"], True, False

    if canon and canon in decoy_lookup:
        info = decoy_lookup[canon]
        return info["name"], None, False, True

    return mol.get("name", "Unknown"), None, False, False


def analyze_target(target_key: str, job_id: str, verbose: bool = True) -> dict:
    """Full analysis for one target."""
    target_data = TARGETS[target_key]
    actives = target_data["actives"]
    decoys_list = DECOYS

    # Build canonical SMILES lookups
    active_lookup = _build_lookup(actives, is_active=True)
    decoy_lookup = _build_lookup(decoys_list, is_active=False)

    mols = fetch_results(job_id)
    if not mols:
        print(f"  ERROR: No results for {target_key}")
        return {}

    # Identify molecules using canonical SMILES matching
    for m in mols:
        name, ic50, is_act, is_dec = identify_molecule(m, active_lookup, decoy_lookup)
        m["_identified"] = name
        m["_ic50"] = ic50
        m["_is_active"] = is_act
        m["_is_decoy"] = is_dec

    # Engine stats
    engines = {}
    for m in mols:
        eng = m.get("docking_engine", "unknown")
        engines[eng] = engines.get(eng, 0) + 1

    # Composite ranking (as returned by API)
    active_ranks_comp = [i for i, m in enumerate(mols, 1) if m["_is_active"]]
    decoy_ranks_comp = [i for i, m in enumerate(mols, 1) if m["_is_decoy"]]

    # Vina ranking (sort by Vina score, best = most negative)
    docked = [m for m in mols if m.get("docking_engine", "") in ("gnina", "gnina_gpu")]
    docked.sort(key=lambda m: float(m.get("vina_score") or 0))
    active_ranks_vina = [i for i, m in enumerate(docked, 1) if m["_is_active"]]
    decoy_ranks_vina = [i for i, m in enumerate(docked, 1) if m["_is_decoy"]]

    # CNN affinity ranking
    docked_cnn = sorted(docked, key=lambda m: -float(m.get("cnn_affinity") or 0))
    active_ranks_cnn = [i for i, m in enumerate(docked_cnn, 1) if m["_is_active"]]
    decoy_ranks_cnn = [i for i, m in enumerate(docked_cnn, 1) if m["_is_decoy"]]

    n = len(docked)

    # Enrichment calculations
    def enrichment(active_ranks, decoy_ranks, n_total):
        if not active_ranks or not decoy_ranks:
            return 0, 0, 0
        avg_act = np.mean(active_ranks)
        avg_dec = np.mean(decoy_ranks)
        ef_ratio = avg_dec / avg_act if avg_act > 0 else 0
        # EF10%
        top_n = max(1, n_total // 10)
        actives_in_top = sum(1 for r in active_ranks if r <= top_n)
        n_actives = len(active_ranks)
        ef10 = (actives_in_top / top_n) / (n_actives / n_total) if n_actives > 0 and top_n > 0 else 0
        return ef_ratio, ef10, avg_act

    ef_vina, ef10_vina, avg_act_vina = enrichment(active_ranks_vina, decoy_ranks_vina, n)
    ef_cnn, ef10_cnn, avg_act_cnn = enrichment(active_ranks_cnn, decoy_ranks_cnn, n)

    # Spearman correlation: docking score vs IC50
    spearman_vina = None
    spearman_cnn = None
    active_docked = [m for m in docked if m["_is_active"] and m["_ic50"] is not None]

    if len(active_docked) >= 4:
        vina_scores = [float(m.get("vina_score") or 0) for m in active_docked]
        cnn_affs = [float(m.get("cnn_affinity") or 0) for m in active_docked]
        ic50s = [m["_ic50"] for m in active_docked]
        log_ic50s = [np.log10(ic + 0.01) for ic in ic50s]

        # Vina: more negative = better binding, lower IC50 = better
        # Correlation between Vina (neg) and log(IC50) should be positive
        if len(set(vina_scores)) > 1 and len(set(log_ic50s)) > 1:
            rho_v, p_v = stats.spearmanr(vina_scores, log_ic50s)
            spearman_vina = {"rho": round(rho_v, 3), "p": round(p_v, 4), "n": len(active_docked)}

        # CNN affinity: higher = better binding, lower IC50 = better
        # Correlation between CNN_aff and log(IC50) should be negative
        if len(set(cnn_affs)) > 1 and len(set(log_ic50s)) > 1:
            rho_c, p_c = stats.spearmanr(cnn_affs, log_ic50s)
            spearman_cnn = {"rho": round(rho_c, 3), "p": round(p_c, 4), "n": len(active_docked)}

    # Score ranges
    vina_scores_all = [float(m.get("vina_score") or 0) for m in docked if m.get("vina_score")]
    cnn_scores_all = [float(m.get("cnn_score") or 0) for m in docked if m.get("cnn_score")]
    cnn_aff_all = [float(m.get("cnn_affinity") or 0) for m in docked if m.get("cnn_affinity")]

    active_vina = [float(m.get("vina_score") or 0) for m in docked if m["_is_active"]]
    decoy_vina = [float(m.get("vina_score") or 0) for m in docked if m["_is_decoy"]]

    result = {
        "target": target_key,
        "uniprot": target_data["uniprot"],
        "n_total": len(mols),
        "n_docked": n,
        "engines": engines,
        "n_actives_found": len(active_ranks_vina),
        "n_decoys_found": len(decoy_ranks_vina),
        "ef_vina": round(ef_vina, 2),
        "ef10_vina": round(ef10_vina, 1),
        "ef_cnn": round(ef_cnn, 2),
        "ef10_cnn": round(ef10_cnn, 1),
        "avg_active_vina": round(np.mean(active_vina), 2) if active_vina else None,
        "avg_decoy_vina": round(np.mean(decoy_vina), 2) if decoy_vina else None,
        "vina_range": [round(min(vina_scores_all), 2), round(max(vina_scores_all), 2)] if vina_scores_all else None,
        "cnn_range": [round(min(cnn_scores_all), 3), round(max(cnn_scores_all), 3)] if cnn_scores_all else None,
        "cnn_aff_range": [round(min(cnn_aff_all), 2), round(max(cnn_aff_all), 2)] if cnn_aff_all else None,
        "spearman_vina": spearman_vina,
        "spearman_cnn": spearman_cnn,
    }

    if verbose:
        print(f"\n{'='*80}")
        print(f"  {target_key} ({target_data['uniprot']})")
        print(f"{'='*80}")
        print(f"  Molecules: {len(mols)} total, {n} docked ({engines})")
        print(f"  Actives found: {len(active_ranks_vina)}/{len(actives)}")
        print(f"  Decoys found:  {len(decoy_ranks_vina)}/{len(decoys_list)}")

        # Detailed per-molecule table
        print(f"\n  {'Rank':>4} {'Name':>25} {'Vina':>8} {'CNN':>6} {'CNN_aff':>8} {'IC50':>10} {'Engine':>10} {'Type':>4}")
        print(f"  {'-'*80}")
        for i, m in enumerate(docked, 1):
            name = m["_identified"][:25]
            vina = m.get("vina_score", "")
            cnn = m.get("cnn_score", "")
            cnn_aff = m.get("cnn_affinity", "")
            ic50 = f"{m['_ic50']:.1f}" if m["_ic50"] else ""
            eng = m.get("docking_engine", "?")
            marker = "ACT" if m["_is_active"] else ("DEC" if m["_is_decoy"] else "")
            print(f"  {i:>4}. {name:>25} {vina:>8} {cnn:>6} {cnn_aff:>8} {ic50:>10} {eng:>10} [{marker}]")

        print(f"\n  Enrichment (Vina):      {ef_vina:.2f}x  EF10%={ef10_vina:.1f}  {'PASS' if ef_vina > 1 else 'FAIL'}")
        print(f"  Enrichment (CNN aff):   {ef_cnn:.2f}x  EF10%={ef10_cnn:.1f}  {'PASS' if ef_cnn > 1 else 'FAIL'}")

        if active_vina and decoy_vina:
            print(f"  Active avg Vina:  {np.mean(active_vina):.2f} kcal/mol")
            print(f"  Decoy avg Vina:   {np.mean(decoy_vina):.2f} kcal/mol")

        if spearman_vina:
            print(f"  Spearman (Vina vs IC50):     rho={spearman_vina['rho']:.3f}, p={spearman_vina['p']:.4f}, n={spearman_vina['n']}")
        if spearman_cnn:
            print(f"  Spearman (CNN_aff vs IC50):  rho={spearman_cnn['rho']:.3f}, p={spearman_cnn['p']:.4f}, n={spearman_cnn['n']}")

        if vina_scores_all:
            print(f"  Vina range: [{min(vina_scores_all):.2f}, {max(vina_scores_all):.2f}]")

    return result


def gpu_vs_cpu_analysis(gpu_job_id: str, cpu_job_id: str):
    """Compare GPU and CPU docking results for the same target."""
    gpu_mols = fetch_results(gpu_job_id)
    cpu_mols = fetch_results(cpu_job_id)

    actives = TARGETS["EGFR"]["actives"]
    active_lookup = _build_lookup(actives, is_active=True)
    decoy_lookup = _build_lookup(DECOYS, is_active=False)

    # Build lookup by name
    gpu_by_smi = {}
    for m in gpu_mols:
        name, _, _, _ = identify_molecule(m, active_lookup, decoy_lookup)
        gpu_by_smi[name] = m

    cpu_by_smi = {}
    for m in cpu_mols:
        name, _, _, _ = identify_molecule(m, active_lookup, decoy_lookup)
        cpu_by_smi[name] = m

    common = set(gpu_by_smi.keys()) & set(cpu_by_smi.keys())

    print(f"\n{'='*80}")
    print(f"  GPU vs CPU CONSISTENCY (EGFR, n={len(common)} shared molecules)")
    print(f"{'='*80}")
    print(f"  {'Drug':>25} {'GPU Vina':>10} {'CPU Vina':>10} {'Diff':>8} {'GPU CNN':>8} {'CPU CNN':>8}")
    print(f"  {'-'*75}")

    vina_diffs = []
    cnn_diffs = []
    gpu_ranks = []
    cpu_ranks = []

    # Sort both by Vina for rank comparison
    gpu_sorted = sorted(common, key=lambda n: float(gpu_by_smi[n].get("vina_score") or 0))
    cpu_sorted = sorted(common, key=lambda n: float(cpu_by_smi[n].get("vina_score") or 0))

    for name in sorted(common):
        gm = gpu_by_smi[name]
        cm = cpu_by_smi[name]
        gv = float(gm.get("vina_score") or 0)
        cv = float(cm.get("vina_score") or 0)
        gc = float(gm.get("cnn_affinity") or 0)
        cc = float(cm.get("cnn_affinity") or 0)

        diff = abs(gv - cv)
        vina_diffs.append(diff)
        if gc and cc:
            cnn_diffs.append(abs(gc - cc))

        gr = gpu_sorted.index(name) + 1
        cr = cpu_sorted.index(name) + 1
        gpu_ranks.append(gr)
        cpu_ranks.append(cr)

        print(f"  {name:>25} {gv:>10.2f} {cv:>10.2f} {diff:>8.2f} {gc:>8.2f} {cc:>8.2f}")

    if vina_diffs:
        print(f"\n  Vina score diff: mean={np.mean(vina_diffs):.2f}, max={max(vina_diffs):.2f} kcal/mol")
    if cnn_diffs:
        print(f"  CNN aff diff:    mean={np.mean(cnn_diffs):.2f}, max={max(cnn_diffs):.2f} pK")

    if len(gpu_ranks) >= 4:
        rho, p = stats.spearmanr(gpu_ranks, cpu_ranks)
        print(f"  Rank correlation: Spearman rho={rho:.3f}, p={p:.4f}")

    # Top-5 overlap
    top5_gpu = set(gpu_sorted[:5])
    top5_cpu = set(cpu_sorted[:5])
    overlap5 = len(top5_gpu & top5_cpu)
    top10_gpu = set(gpu_sorted[:10])
    top10_cpu = set(cpu_sorted[:10])
    overlap10 = len(top10_gpu & top10_cpu)
    print(f"  Top-5 overlap: {overlap5}/5 ({100*overlap5/5:.0f}%)")
    print(f"  Top-10 overlap: {overlap10}/10 ({100*overlap10/10:.0f}%)")


def summary_table(results: list[dict]):
    """Print publication-ready summary table."""
    print(f"\n{'='*100}")
    print(f"  PUBLICATION SUMMARY TABLE — DockIt 5-Target Benchmark")
    print(f"{'='*100}")

    headers = ["Metric"] + [r["target"] for r in results]
    col_w = 14

    print(f"  {'Metric':<35}", end="")
    for r in results:
        print(f" {r['target']:>{col_w}}", end="")
    print()
    print(f"  {'-'*35}", end="")
    for _ in results:
        print(f" {'-'*col_w}", end="")
    print()

    rows = [
        ("UniProt", lambda r: r["uniprot"]),
        ("N molecules", lambda r: str(r["n_total"])),
        ("N docked (GNINA)", lambda r: str(r["n_docked"])),
        ("GNINA success %", lambda r: f"{100*r['n_docked']/r['n_total']:.0f}%"),
        ("Actives identified", lambda r: str(r["n_actives_found"])),
        ("Vina range (kcal/mol)", lambda r: f"[{r['vina_range'][0]}, {r['vina_range'][1]}]" if r.get("vina_range") else "N/A"),
        ("CNN score range", lambda r: f"[{r['cnn_range'][0]}, {r['cnn_range'][1]}]" if r.get("cnn_range") else "N/A"),
        ("CNN aff range (pK)", lambda r: f"[{r['cnn_aff_range'][0]}, {r['cnn_aff_range'][1]}]" if r.get("cnn_aff_range") else "N/A"),
        ("Active avg Vina", lambda r: f"{r['avg_active_vina']}" if r.get("avg_active_vina") else "N/A"),
        ("Decoy avg Vina", lambda r: f"{r['avg_decoy_vina']}" if r.get("avg_decoy_vina") else "N/A"),
        ("EF (Vina)", lambda r: f"{r['ef_vina']:.2f}x"),
        ("EF10% (Vina)", lambda r: f"{r['ef10_vina']:.1f}"),
        ("EF (CNN aff)", lambda r: f"{r['ef_cnn']:.2f}x"),
        ("EF10% (CNN aff)", lambda r: f"{r['ef10_cnn']:.1f}"),
        ("Spearman rho (Vina/IC50)", lambda r: f"{r['spearman_vina']['rho']:.3f} (p={r['spearman_vina']['p']:.3f})" if r.get("spearman_vina") else "N/A"),
        ("Spearman rho (CNN/IC50)", lambda r: f"{r['spearman_cnn']['rho']:.3f} (p={r['spearman_cnn']['p']:.3f})" if r.get("spearman_cnn") else "N/A"),
    ]

    for label, fn in rows:
        print(f"  {label:<35}", end="")
        for r in results:
            val = fn(r)
            print(f" {val:>{col_w}}", end="")
        print()


if __name__ == "__main__":
    # Job IDs from command line or files
    job_ids = {}
    targets_to_run = ["EGFR", "CDK2", "BRAF_V600E", "JAK2", "KRAS_G12C"]
    file_map = {"EGFR": "egfr", "CDK2": "cdk2", "BRAF_V600E": "braf_v600e", "JAK2": "jak2", "KRAS_G12C": "kras_g12c"}

    for target_key in targets_to_run:
        fpath = Path(f"/tmp/bench_final_{file_map[target_key]}.txt")
        if fpath.exists():
            job_ids[target_key] = fpath.read_text().strip()

    if not job_ids:
        print("No job IDs found in /tmp/bench_final_*.txt")
        sys.exit(1)

    print(f"Analyzing {len(job_ids)} targets...")

    results = []
    for target_key, job_id in job_ids.items():
        r = analyze_target(target_key, job_id)
        if r:
            results.append(r)

    if results:
        summary_table(results)

    # GPU vs CPU if available
    gpu_file = Path("/tmp/bench_final_egfr.txt")
    cpu_file = Path("/tmp/bench_final_egfr_cpu.txt")
    if gpu_file.exists() and cpu_file.exists():
        gpu_id = gpu_file.read_text().strip()
        cpu_id = cpu_file.read_text().strip()
        gpu_vs_cpu_analysis(gpu_id, cpu_id)

    # Save results as JSON
    output = Path("/tmp/benchmark_results.json")
    with open(output, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to {output}")
