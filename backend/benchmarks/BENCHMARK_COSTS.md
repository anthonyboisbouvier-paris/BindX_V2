# BindX — Benchmark Costs & Performance Tracking

## GNINA GPU (RunPod Serverless) — 2026-03-10

**Endpoint**: `weeeiy6z4jdsv3` | **GPU**: A40 / L40S | **Rate**: $0.00016/s

### Throughput by configuration

| Config | Batch | Mols | Time (s) | Throughput (mol/min) | Cost ($) | Rejects (%) | Notes |
|--------|-------|------|----------|---------------------|----------|-------------|-------|
| rescore, exh=8 | 20 | 20 | 38 | 31.6 | 0.006 | 0% | Baseline |
| rescore, exh=16 | 20 | 20 | 44 | 27.3 | 0.007 | 0% | +16% time |
| rescore, exh=32 | 20 | 20 | 52 | 23.1 | 0.008 | 0% | +37% time |
| refinement, exh=8 | 20 | 20 | 102 | 11.8 | 0.016 | 17% | x2.7 slower, positive Vina rejects |
| refinement, exh=16 | 20 | 20 | 125 | 9.6 | 0.020 | 17% | x3.3 slower |
| refinement, exh=32 | 20 | 20 | TIMEOUT | - | - | - | >600s, never completes |
| rescore, exh=8 | 100 | 100 | 105 | 57.1 | 0.017 | 0% | Sub-linear scaling |
| rescore, exh=8 | 250 | 250 | 230 | 65.2 | 0.037 | 0% | Optimal batch size |

### Key findings

- **Rescore** is the production default — fast, no rejects
- **Refinement** is ~2.7x slower and rejects ~17% of molecules (positive Vina scores)
- **Refinement + exh>=32** systematically timeouts (>600s for 20 mols)
- **Scaling** is sub-linear: throughput increases with batch size up to ~250
- **Cost model**: `time_seconds * $0.00016` (RunPod A40 serverless rate)

### Guards implemented

- `refinement` mode: batch_size capped at 10, polling timeout doubled
- `refinement + exh>=32`: exhaustiveness capped at 16 with warning log
- Partial results saved batch-by-batch (timeout doesn't lose all results)

---

*Updated: 2026-03-10. Add new rows after each benchmark run.*
