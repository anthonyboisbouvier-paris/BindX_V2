# Generation Features — Implementation Specs

## Status Overview

| Feature | Frontend UI | Backend Pipeline | Backend Task | Status |
|---------|------------|-----------------|-------------|--------|
| Scaffold Hopping (batch) | OK | OK (`generation.py`) | OK (`_run_generation`) | **WORKING** |
| Fragment Growing | OK | OK (`generation.py`) | OK (`_run_generation`) | **WORKING** |
| De Novo SMILES | OK | OK (`generation.py`) | OK (`_run_generation`) | **WORKING** |
| Per-Molecule R-group control | OK (ScaffoldAnalyzer) | OK (`scaffold_analysis.py`) | OK (`_run_lead_optimization`) | **WORKING** |
| Lead Optimization (iterative) | OK (RunCreator) | OK (`lead_optimization.py`) | OK (`_run_lead_optimization`) | **WORKING** |
| Post-generation analyses | OK (toggles) | OK (chaining) | OK (`_run_generation`) | **WORKING** |
| Pharmacophore mapping | OK | OK (`pharmacophore.py`) | OK (`_run_calculation`) | **WORKING** |

---

## SPEC 1: Fragment Growing

### Description
Extend molecules at reactive positions (BRICS bonds, terminal groups) by adding fragments from a curated library.

### Backend Implementation
**File**: `backend/pipeline/generation.py`

1. Add `generate_fragment_growing()` function:
   - Input: `smiles` (parent molecule), `fragment_library` (default BRICS fragments), `max_new_atoms` (default 15)
   - Process:
     1. Detect BRICS bonds on parent molecule (use `BRICS.BRICSDecompose`)
     2. For each BRICS site, enumerate possible fragment attachments
     3. Filter by: MW delta < max_mw_change, QED > min_qed, Lipinski compliance
     4. Score remaining candidates by estimated binding (shape similarity to seed)
     5. Return top N candidates
   - Output: list of `{smiles, name, parent_smiles, fragment_added, estimated_affinity, qed, novelty_score}`

2. Reference implementation: V1 `lead_optimization.py:668-881` (`_generate_variants_rdkit`)
   - Strategy `add_fg`: adds functional groups at aromatic C-H positions
   - Strategy `modify_chain`: extends/modifies aliphatic chains
   - 15 functional groups: C, CC, O, OC, N, NC, F, Cl, CF3, CN, etc.

### Task Integration
**File**: `backend/tasks_v9.py` — `_run_generation()`

- When `config.method == "fragment_growing"`:
  - Call `generate_fragment_growing()` instead of `generate_molecules()`
  - Same result storage pattern as scaffold hopping

### Frontend
Already implemented in RunCreator — just remove `implemented: false` flag.

---

## SPEC 2: De Novo SMILES Generation

### Description
Generate completely novel molecular structures using REINVENT4 RL-based generation, constrained by target pocket shape and pharmacophoric features.

### Backend Implementation
**File**: `backend/pipeline/generation.py`

The `generate_molecules()` function already supports REINVENT4 as primary backend:
1. Check if REINVENT4 is available (Docker/RunPod)
2. If yes: run full RL generation with receptor constraints
3. If no: fall back to RDKit-based scaffold mutation

For de novo mode specifically:
- Do NOT pass `seed_smiles` to `generate_molecules()` — this triggers unconstrained generation
- Use larger `n_molecules` (200+) since de novo has lower hit rate
- Apply stricter post-filtering: QED > 0.5, Lipinski pass, PAINS clean

### Task Integration
**File**: `backend/tasks_v9.py` — `_run_generation()`

- When `config.method == "de_novo"`:
  - Call `generate_molecules()` with `seed_smiles=None` (unconstrained)
  - Higher `n_molecules` (200), lower `n_top` (50)

### Frontend
Already implemented — just remove `implemented: false` flag.

---

## SPEC 3: Per-Molecule R-group Control

### Description
Allow users to define structural modification rules per molecule position: freeze specific positions, select allowed strategies (add_fg, swap_halogen, swap_atom, modify_chain), whitelist specific R-groups.

### Current State
- **Frontend**: ScaffoldAnalyzer component exists and works (lazy-loaded in RunCreator)
  - Calls `/api/agent/analyze-scaffold` (V8 legacy endpoint)
  - Returns position cards with freeze/strategy/R-group controls
  - Outputs `scaffold_rules` object to RunCreator config

- **Backend**: `pipeline/scaffold_analysis.py` exists (395 lines, V1 port)
  - `analyze_scaffold(smiles)` returns Murcko scaffold, BRICS bonds, R-group positions

- **Missing**: Connection between RunCreator rules → generation task

### Implementation Plan

1. **Backend endpoint** (new): `POST /api/v9/scaffold-analysis`
   ```python
   @router.post("/scaffold-analysis")
   async def scaffold_analysis(body: dict, user_id = Depends(require_v9_user)):
       from pipeline.scaffold_analysis import analyze_scaffold
       result = analyze_scaffold(body["smiles"])
       return result
   ```

2. **Frontend**: Update ScaffoldAnalyzer to call V9 endpoint instead of V8 legacy

3. **Task integration**: `_run_generation()` reads `config.scaffold_rules`:
   ```python
   structural_rules = config.get("scaffold_rules")
   if structural_rules:
       from pipeline.lead_optimization import run_optimization
       result = run_optimization(
           starting_smiles=seed_smiles[0],
           starting_name=parent_mols[0]["name"],
           target_pdbqt=str(pdb_path),
           pocket_center=pocket_center,
           structural_rules=structural_rules,
           n_iterations=config.get("iterations", 3),
           variants_per_iter=config.get("variants_per_iteration", 10),
       )
   ```

4. **Rules format** (from ScaffoldAnalyzer):
   ```json
   {
     "rules": [
       {"position_idx": 3, "strategy": "add_fg", "allowed_groups": ["F", "Cl"], "frozen": false},
       {"position_idx": 7, "strategy": "any", "allowed_groups": [], "frozen": true}
     ],
     "frozen_positions": [7],
     "allowed_strategies": ["add_fg", "swap_halogen"],
     "preserve_scaffold": true,
     "min_similarity": 0.3,
     "max_mw_change": 100,
     "core_atom_indices": [0, 1, 2, 3, 4, 5]
   }
   ```

---

## SPEC 4: Lead Optimization (Multi-iteration)

### Description
Iterative multi-objective optimization: generate variants → dock → score → filter → repeat. Uses 4 optimization weights (binding_affinity, toxicity, bioavailability, synthesis_ease).

### Current State
- **Backend**: `pipeline/lead_optimization.py` exists (1049 lines, V1 port)
  - `run_optimization()` — full multi-iteration loop with progress callback
  - 4 variant generation strategies (RDKit-based)
  - Multi-objective ranking with configurable weights
- **Frontend**: NO UI (V1 had `OptimizationView.jsx` — 605 lines)

### Implementation Plan

1. **Frontend component**: Create `OptimizationView.jsx` (port from V1)
   - Setup phase: weight sliders (4 objectives), iteration count, variants per iter
   - Running phase: progress bar, score evolution chart, ETA
   - Complete phase: before/after comparison, score chart, download results

2. **Backend task**: New handler in `_run_generation()` when mode is "optimization":
   ```python
   if config.get("mode") == "optimization":
       from pipeline.lead_optimization import run_optimization
       result = run_optimization(
           starting_smiles=seed_smiles[0],
           starting_name=parent_mols[0]["name"],
           target_pdbqt=str(pdb_path),
           pocket_center=pocket_center,
           weights=config.get("weights", {"binding_affinity": 0.35, "toxicity": 0.25, "bioavailability": 0.20, "synthesis_ease": 0.20}),
           n_iterations=config.get("iterations", 5),
           variants_per_iter=config.get("variants_per_iteration", 50),
           structural_rules=config.get("scaffold_rules"),
           progress_callback=lambda pct, score, step: _db_sync(_update_run(run_id, progress=int(pct), current_step=step)),
       )
   ```

3. **New run type** or subtype: `type="generation"` with `config.mode="optimization"`

---

## SPEC 5: Post-Generation Analyses

### Description
Auto-run docking, ADMET, and/or scoring on generated molecules immediately after generation completes.

### Implementation Plan

1. **Backend task**: After `_run_generation()` completes, check config flags:
   ```python
   if config.get("include_docking"):
       # Create and execute a calculation run for docking on new molecules
   if config.get("include_admet"):
       # Create and execute ADMET calculation
   if config.get("include_scoring"):
       # Create and execute scoring
   ```

2. **Approach**: Chain Celery tasks or execute sequentially within the same task
   - Option A: Sequential execution within `_run_generation` (simpler, single progress bar)
   - Option B: Create child runs (more granular, separate progress tracking)

3. **Frontend**: Already has UI toggles — just connect them by passing config to the run

---

## SPEC 6: Pharmacophore Mapping

### Description
Map 3D pharmacophoric features (H-bond donors/acceptors, hydrophobic regions, aromatic rings, charge centers) and compare across molecule series.

### Implementation Plan

1. **Backend module**: `pipeline/pharmacophore.py` (new)
   - Use RDKit `Chem.Pharm2D` or `rdkit.Chem.Features`
   - Extract pharmacophore features from 3D conformer
   - Compute pharmacophore fingerprints for clustering
   - Output: feature map (type, position, radius), pairwise similarity matrix

2. **Frontend**: Visualization in ProteinViewer (3D spheres for features)

3. **Task integration**: New calculation subtype `pharmacophore` in enrichment

---

## Priority Order

1. **Per-Molecule R-group control** (highest value, most code already exists)
2. **Fragment Growing** (reuse V1 variant generation code)
3. **Post-Generation Analyses** (chain existing pipelines)
4. **Lead Optimization UI** (port V1 OptimizationView)
5. **De Novo SMILES** (REINVENT4 dependency)
6. **Pharmacophore** (new module required)
