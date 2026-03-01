# BindX V9 — E2E Test Protocol

> **Version**: 1.0 (2026-03-01)
> **Objectif**: Protocole de test end-to-end reproductible, utilisable par Claude Code pour valider chaque livraison.
> **Principe**: Simuler exactement le parcours d'un utilisateur via l'API, en respectant toutes les contraintes du frontend.

---

## Table des matieres

1. [Principes du protocole](#1-principes)
2. [Pre-requis](#2-pre-requis)
3. [Scenario complet](#3-scenario-complet)
4. [Regles metier a respecter](#4-regles-metier)
5. [Checklist par page](#5-checklist-par-page)
6. [Commandes de test](#6-commandes-de-test)
7. [Maintenance](#7-maintenance)

---

## 1. Principes

### Ce que ce protocole garantit

- Chaque test suit **exactement** le parcours utilisateur du frontend
- Aucun raccourci API qui contournerait les contraintes UI
- Validation a chaque etape : l'etat backend correspond a ce que le frontend afficherait
- Screenshots via navigateur headless pour validation visuelle (quand disponible)

### Regle d'or

> **Ne jamais faire via l'API ce que l'utilisateur ne pourrait pas faire via le frontend.**

Exemples de violations interdites :
- Lancer un run sans avoir configure la cible
- Passer de Phase A a Phase C sans creer Phase B
- Lancer un run calcul sans selectionner de molecules
- Creer un run sur une phase frozen
- Lancer un 2e run pendant qu'un run est en cours

---

## 2. Pre-requis

### Comptes de test

```
Email:    anthony.boisbouvier@gmail.com
Password: test1234
```

### Infrastructure

```bash
# Containers Docker running
docker-compose up --build -d
# Verifier : backend healthy, frontend healthy, redis healthy, celery running
docker-compose ps
```

### Obtenir un token

```bash
TOKEN=$(curl -s 'https://webkntghfzscrnuixfba.supabase.co/auth/v1/token?grant_type=password' \
  -H 'apikey: sb_publishable_Zy18OaVJ4k_DgaafyiYaZA_9jB0ATmY' \
  -H 'Content-Type: application/json' \
  -d '{"email":"anthony.boisbouvier@gmail.com","password":"test1234"}' \
  | python3 -c "import sys,json; print(json.load(sys.stdin)['access_token'])")
```

### Base URL

```
API directe:  http://localhost:8000
Via nginx:    http://localhost:3000  (comme le frontend)
```

> Toujours tester via **port 3000** (nginx proxy) pour reproduire les conditions reelles.

---

## 3. Scenario complet

### Phase 0 : Nettoyage (optionnel)

Si un projet de test existe deja, le supprimer ou en creer un nouveau.

---

### Etape 1 : Creer un projet

**Action frontend** : Page `/welcome` → "New Project" → nom + description

```bash
PROJECT=$(curl -s -X POST http://localhost:3000/api/v9/projects \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"name":"E2E Test Protocol","description":"Automated validation"}' \
  | python3 -c "import sys,json; d=json.load(sys.stdin); print(d['id'])")
echo "Project: $PROJECT"
```

**Validations** :
- [ ] Projet cree avec `status: "active"`
- [ ] Une campagne par defaut est creee automatiquement
- [ ] `target_preview` est `null` (pas encore configure)

```bash
curl -s http://localhost:3000/api/v9/projects/$PROJECT \
  -H "Authorization: Bearer $TOKEN" | python3 -c "
import sys,json
p = json.load(sys.stdin)
assert p['status'] == 'draft', f'Expected draft, got {p[\"status\"]}'
assert len(p['campaigns']) >= 1, 'No default campaign'
assert p.get('target_preview') is None, 'Target should be None'
print('PASS: Project created correctly')
print(f'  Campaign ID: {p[\"campaigns\"][0][\"id\"]}')
"
```

---

### Etape 2a : Resoudre la structure 3D (VRAI pipeline)

**Action frontend** : Page `/project/:id/target-setup` → UniProt ID → Preview

**Regle** : Le frontend bloque la creation de campagne tant que `target_preview` est null.

> **IMPORTANT** : Cette etape appelle les VRAIES APIs externes (RCSB PDB, UniProt, ChEMBL, P2Rank).
> Duree : 30 a 120 secondes.

```bash
# Appel reel : UniProt → RCSB PDB → P2Rank pocket detection
PREVIEW=$(curl -s -X POST http://localhost:3000/api/preview-target \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"uniprot_id": "P00533"}' \
  --max-time 180)
```

**Validations** :
- [ ] Structure resolue (`structure.source != 'none'`)
- [ ] PDB ID present (`structure.pdb_id`)
- [ ] Resolution disponible (`structure.resolution`)
- [ ] Au moins 1 pocket detectee
- [ ] Nom de proteine present (`protein_name`)

### Etape 2b : Valider la pocket pour le docking

**Regle** : La pocket selectionnee DOIT avoir des coordonnees `center: [x, y, z]` valides.
Sans ces coordonnees, GNINA ne peut pas definir sa boite de recherche.

**Validations** :
- [ ] `pockets[0].center` est un tableau de 3 floats
- [ ] `pockets[0].method` indique la methode de detection (co-crystallized_ligand, p2rank, etc.)
- [ ] `pockets[0].probability` > 0

### Etape 2c : Sauvegarder la cible dans le projet

**Action frontend** : TargetSetup → "Validate Target"

```bash
# Construire le payload comme le frontend le fait
# Inclure selected_pocket_index = 0 (premiere pocket)
TARGET_DATA=$(echo "$PREVIEW" | python3 -c "...")  # voir e2e_test.sh
curl -s -X PUT http://localhost:3000/api/v9/projects/$PROJECT \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d "$TARGET_DATA"
```

**Validations** :
- [ ] `target_preview` sauvegarde non null
- [ ] `target_name` present
- [ ] `target_pdb_id` present
- [ ] `pockets_detected` sauvegarde avec coordonnees `center`
- [ ] `target_preview.selected_pocket_index` = 0
- [ ] Frontend afficherait maintenant le bouton "Create Campaign"

---

### Etape 3 : Recuperer la campagne par defaut + creer Phase A

**Action frontend** : Page `/project/:id` → la campagne par defaut existe → Clic "Phase A"

```bash
# Recuperer la campagne
CAMPAIGN=$(curl -s http://localhost:3000/api/v9/projects/$PROJECT \
  -H "Authorization: Bearer $TOKEN" | python3 -c "
import sys,json
p = json.load(sys.stdin)
print(p['campaigns'][0]['id'])
")

# Creer Phase A (hit_discovery)
PHASE_A=$(curl -s -X POST http://localhost:3000/api/v9/campaigns/$CAMPAIGN/phases \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"type": "hit_discovery"}' \
  | python3 -c "import sys,json; print(json.load(sys.stdin)['id'])")
echo "Phase A: $PHASE_A"
```

**Validations** :
- [ ] Phase creee avec `type: "hit_discovery"`, `status: "active"`
- [ ] Pas de molecules (molecule_count = 0)
- [ ] Pas de runs

---

### Etape 4 : Run Import (SMILES)

**Action frontend** : PhaseDashboard → "New Run" → Import → SMILES list

**Regles** :
- Import ne necessite PAS de selection de molecules
- Il faut fournir `smiles_list` OU un fichier OU `source: phase_selection`

```bash
# 8 molecules de test (SMILES valides)
IMPORT_RUN=$(curl -s -X POST "http://localhost:3000/api/v9/phases/$PHASE_A/runs" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "type": "import",
    "config": {
      "smiles_list": [
        "c1ccccc1",
        "CC(=O)Oc1ccccc1C(=O)O",
        "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
        "CC12CCC3C(CCC4CC(=O)CCC43C)C1CCC2O",
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
        "CC(=O)NC1=CC=C(O)C=C1",
        "C1CCCCC1"
      ]
    }
  }' | python3 -c "import sys,json; d=json.load(sys.stdin); print(d['id'])")
echo "Import run: $IMPORT_RUN"
```

**Attendre la completion** (Celery traite le run) :

```bash
for i in $(seq 1 30); do
  STATUS=$(curl -s "http://localhost:3000/api/v9/runs/$IMPORT_RUN" \
    -H "Authorization: Bearer $TOKEN" | python3 -c "import sys,json; print(json.load(sys.stdin)['status'])")
  echo "Poll $i: $STATUS"
  [ "$STATUS" = "completed" ] && break
  [ "$STATUS" = "failed" ] && echo "FAIL: Import failed" && break
  sleep 2
done
```

**Validations** :
- [ ] Run status = `completed`
- [ ] Run a des logs (>= 2 entries)
- [ ] Molecules creees dans la phase

```bash
curl -s "http://localhost:3000/api/v9/phases/$PHASE_A/molecules?limit=10" \
  -H "Authorization: Bearer $TOKEN" | python3 -c "
import sys,json
d = json.load(sys.stdin)
assert d['total'] > 0, f'No molecules imported (total={d[\"total\"]})'
print(f'PASS: {d[\"total\"]} molecules imported')
for m in d['molecules'][:3]:
    print(f'  {m[\"name\"]}: {m[\"smiles\"][:40]}')
"
```

---

### Etape 5 : Verifier qu'on ne peut PAS lancer un 2e run pendant le 1er

**Regle frontend** : Si un run est `created` ou `running`, le bouton "New Run" est desactive.

> Note : Cette verification doit etre faite pendant que le run tourne.
> Si le run est deja termine, cette etape est informative seulement.

```bash
# Le backend n'empeche pas 2 runs simultanes (CDC dit "1 run par user V9")
# C'est le FRONTEND qui bloque via :
#   const activeRun = runs.find(r => r.status === 'created' || r.status === 'running')
#   if (activeRun) { toast('A run is already in progress...'); return }
echo "INFO: Concurrent run check is frontend-only (no backend enforcement)"
```

---

### Etape 6 : Selectionner des molecules + Run Calculation

**Action frontend** : Dashboard → cocher des molecules → "New Run" → Calculation

**Regles** :
- Run calculation NECESSITE des molecules selectionnees (`input_molecule_ids` non vide)
- Run calculation NECESSITE au moins 1 `calculation_type`
- Types valides : docking, admet, scoring, enrichment, clustering, off_target, confidence, retrosynthesis, safety

```bash
# Recuperer les IDs des 4 premieres molecules (simule la selection dans le tableau)
MOL_IDS=$(curl -s "http://localhost:3000/api/v9/phases/$PHASE_A/molecules?limit=4" \
  -H "Authorization: Bearer $TOKEN" | python3 -c "
import sys,json
d = json.load(sys.stdin)
ids = [m['id'] for m in d['molecules'][:4]]
print(json.dumps(ids))
")
echo "Selected molecules: $MOL_IDS"

# Lancer un run calcul (admet + scoring — pas docking car GPU indisponible)
CALC_RUN=$(curl -s -X POST "http://localhost:3000/api/v9/phases/$PHASE_A/runs" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d "{
    \"type\": \"calculation\",
    \"calculation_types\": [\"admet\", \"scoring\"],
    \"input_molecule_ids\": $MOL_IDS
  }" | python3 -c "import sys,json; d=json.load(sys.stdin); print(d['id'])")
echo "Calculation run: $CALC_RUN"
```

**Attendre completion** :

```bash
for i in $(seq 1 60); do
  STATUS=$(curl -s "http://localhost:3000/api/v9/runs/$CALC_RUN" \
    -H "Authorization: Bearer $TOKEN" | python3 -c "import sys,json; print(json.load(sys.stdin)['status'])")
  echo "Poll $i: $STATUS"
  [ "$STATUS" = "completed" ] && break
  [ "$STATUS" = "failed" ] && echo "FAIL: Calculation failed" && break
  sleep 2
done
```

**Validations** :
- [ ] Run status = `completed`
- [ ] Run a `calculation_types: ["admet", "scoring"]`
- [ ] Run a des logs
- [ ] Les 4 molecules selectionnees ont maintenant des `properties`
- [ ] Properties contiennent les cles ADMET (logP, MW, solubility, etc.)
- [ ] `composite_score` calcule par le scoring

```bash
curl -s "http://localhost:3000/api/v9/phases/$PHASE_A/molecules?limit=20" \
  -H "Authorization: Bearer $TOKEN" | python3 -c "
import sys,json

ALIASES = {'hbd':'HBD','hba':'HBA','qed':'QED','tpsa':'TPSA',
  'bbb_permeability':'BBB','herg_inhibition':'hERG','color_code':'safety_color_code'}
SKIP = {'flags','note','status','smiles','confidence_modifier','nearest_tanimoto','docking_status'}

def deep_flatten(obj, out=None):
    if out is None: out = {}
    for k,v in obj.items():
        if k in SKIP: continue
        if isinstance(v, dict): deep_flatten(v, out)
        else: out[ALIASES.get(k, k)] = v
    return out

d = json.load(sys.stdin)
with_props = [m for m in d['molecules'] if m.get('properties')]
print(f'Molecules with properties: {len(with_props)}/{d[\"total\"]}')
assert len(with_props) >= 4, f'Expected >=4 molecules with properties, got {len(with_props)}'

# Verify key ADMET columns exist
m = with_props[0]
flat = deep_flatten(m['properties'])
expected = ['logP', 'MW', 'HBD', 'HBA', 'QED', 'composite_score']
missing = [k for k in expected if k not in flat]
assert not missing, f'Missing keys: {missing}'
print(f'PASS: Properties contain {sorted(flat.keys())}')
"
```

---

### Etape 7 : Verifier les runs (logs, config, types)

**Ce que le frontend affiche** : RunHistory avec type label, config summary, logs

```bash
curl -s "http://localhost:3000/api/v9/phases/$PHASE_A/runs" \
  -H "Authorization: Bearer $TOKEN" | python3 -c "
import sys,json
runs = json.load(sys.stdin)
print(f'Total runs: {len(runs)}')
for r in runs:
    logs = r.get('logs', [])
    ct = r.get('calculation_types')
    cfg = r.get('config', {})
    print(f'  {r[\"type\"]:12s} | status={r[\"status\"]:10s} | logs={len(logs)} | calc_types={ct}')
    assert 'id' in r, 'Missing id'
    assert 'status' in r, 'Missing status'
    assert 'created_at' in r, 'Missing created_at'
    if r['type'] == 'calculation':
        assert ct is not None, 'calculation run missing calculation_types'
        assert len(logs) > 0, 'calculation run has no logs'
    if r['type'] == 'import':
        assert len(logs) > 0, 'import run has no logs'
print('PASS: All runs have correct structure')
"
```

---

### Etape 8 : Bookmarker des molecules

**Action frontend** : Dashboard → cocher molecules → "Bookmark Selected"

**Regles** :
- Necessite au moins 1 molecule selectionnee
- Phase ne doit PAS etre frozen

```bash
# Bookmark les 4 molecules qui ont des properties
curl -s -X POST "http://localhost:3000/api/v9/phases/$PHASE_A/molecules/bookmark-batch" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d "{\"molecule_ids\": $MOL_IDS, \"bookmarked\": true}" | python3 -c "
import sys,json
d = json.load(sys.stdin)
print(f'PASS: Bookmarked {d.get(\"updated\", \"?\")} molecules')
"
```

**Validation** :

```bash
curl -s "http://localhost:3000/api/v9/phases/$PHASE_A/molecules/stats" \
  -H "Authorization: Bearer $TOKEN" | python3 -c "
import sys,json
s = json.load(sys.stdin)
assert s['bookmarked'] >= 4, f'Expected >=4 bookmarked, got {s[\"bookmarked\"]}'
print(f'PASS: Stats — total={s[\"total\"]}, bookmarked={s[\"bookmarked\"]}')
"
```

---

### Etape 9 : Freeze Phase A

**Action frontend** : PhaseHeader → Freeze → FreezeDialog confirmation

**Regles** :
- Apres freeze : pas de nouveau run, pas de bookmark, pas d'annotation
- La selection bookmarkee est verrouillee pour la phase suivante

```bash
curl -s -X POST "http://localhost:3000/api/v9/phases/$PHASE_A/freeze" \
  -H "Authorization: Bearer $TOKEN" | python3 -c "
import sys,json
d = json.load(sys.stdin)
assert d['status'] == 'frozen', f'Expected frozen, got {d[\"status\"]}'
assert d.get('bookmarked_molecule_count', 0) >= 4, 'Expected >=4 bookmarked molecules'
print(f'PASS: Phase A frozen with {d[\"bookmarked_molecule_count\"]} bookmarked molecules')
"
```

---

### Etape 10 : Verifier que la phase frozen bloque les actions

**Regles frontend** :
- "New Run" desactive
- Bookmark desactive
- Annotations desactivees

```bash
# Tenter de creer un run sur phase frozen → doit echouer (400)
HTTP_CODE=$(curl -s -o /dev/null -w "%{http_code}" -X POST "http://localhost:3000/api/v9/phases/$PHASE_A/runs" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"type": "import", "config": {"smiles_list": ["CCCC"]}}')
echo "Create run on frozen phase: HTTP $HTTP_CODE"
[ "$HTTP_CODE" = "400" ] && echo "PASS: Frozen phase blocks run creation" || echo "FAIL: Expected 400"
```

---

### Etape 11 : Creer Phase B (hit_to_lead)

**Action frontend** : ProjectHome → Clic sur slot "Phase B" dans le funnel

**Regle** : On doit creer Phase B AVANT de pouvoir y envoyer des molecules

```bash
PHASE_B=$(curl -s -X POST "http://localhost:3000/api/v9/campaigns/$CAMPAIGN/phases" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"type": "hit_to_lead"}' \
  | python3 -c "import sys,json; print(json.load(sys.stdin)['id'])")
echo "Phase B: $PHASE_B"
```

**Validation** :
- [ ] Phase B creee avec status `active`
- [ ] 0 molecules

---

### Etape 12 : Envoyer les bookmarks de Phase A vers Phase B

**Action frontend** : PhaseDashboard (Phase A) → "Send to Next Phase" → import internal

**Regles** :
- Source = `phase_selection`
- `source_phase_id` = Phase A ID
- Phase A doit etre frozen
- Phase B doit exister

```bash
NEXT_RUN=$(curl -s -X POST "http://localhost:3000/api/v9/phases/$PHASE_B/runs" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d "{
    \"type\": \"import\",
    \"config\": {
      \"source\": \"phase_selection\",
      \"source_phase_id\": \"$PHASE_A\"
    }
  }" | python3 -c "import sys,json; d=json.load(sys.stdin); print(d['id'])")
echo "Next phase import run: $NEXT_RUN"
```

**Attendre completion + validation** :

```bash
for i in $(seq 1 30); do
  STATUS=$(curl -s "http://localhost:3000/api/v9/runs/$NEXT_RUN" \
    -H "Authorization: Bearer $TOKEN" | python3 -c "import sys,json; print(json.load(sys.stdin)['status'])")
  [ "$STATUS" = "completed" ] && break
  [ "$STATUS" = "failed" ] && echo "FAIL" && break
  sleep 2
done

curl -s "http://localhost:3000/api/v9/phases/$PHASE_B/molecules?limit=10" \
  -H "Authorization: Bearer $TOKEN" | python3 -c "
import sys,json
d = json.load(sys.stdin)
assert d['total'] >= 4, f'Expected >=4 molecules in Phase B, got {d[\"total\"]}'
print(f'PASS: Phase B has {d[\"total\"]} molecules (from Phase A bookmarks)')
"
```

---

### Etape 13 : Run calcul sur Phase B

**Action frontend** : PhaseDashboard (Phase B) → selectionner molecules → New Run → Calculation

```bash
# Selectionner toutes les molecules de Phase B
MOL_IDS_B=$(curl -s "http://localhost:3000/api/v9/phases/$PHASE_B/molecules?limit=10" \
  -H "Authorization: Bearer $TOKEN" | python3 -c "
import sys,json
d = json.load(sys.stdin)
ids = [m['id'] for m in d['molecules']]
print(json.dumps(ids))
")

# Run ADMET + safety
CALC_B=$(curl -s -X POST "http://localhost:3000/api/v9/phases/$PHASE_B/runs" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d "{
    \"type\": \"calculation\",
    \"calculation_types\": [\"admet\", \"safety\", \"scoring\"],
    \"input_molecule_ids\": $MOL_IDS_B
  }" | python3 -c "import sys,json; d=json.load(sys.stdin); print(d['id'])")
echo "Phase B calc run: $CALC_B"

# Attendre
for i in $(seq 1 60); do
  STATUS=$(curl -s "http://localhost:3000/api/v9/runs/$CALC_B" \
    -H "Authorization: Bearer $TOKEN" | python3 -c "import sys,json; print(json.load(sys.stdin)['status'])")
  [ "$STATUS" = "completed" ] && break
  [ "$STATUS" = "failed" ] && break
  sleep 2
done
echo "Phase B calc status: $STATUS"
```

---

### Etape 14 : Bookmark + Freeze Phase B + Creer Phase C

```bash
# Bookmark les 2 premieres molecules de Phase B
MOL_2=$(curl -s "http://localhost:3000/api/v9/phases/$PHASE_B/molecules?limit=2" \
  -H "Authorization: Bearer $TOKEN" | python3 -c "
import sys,json; print(json.dumps([m['id'] for m in json.load(sys.stdin)['molecules'][:2]]))")

curl -s -X POST "http://localhost:3000/api/v9/phases/$PHASE_B/molecules/bookmark-batch" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d "{\"molecule_ids\": $MOL_2, \"bookmarked\": true}" > /dev/null

# Freeze Phase B
curl -s -X POST "http://localhost:3000/api/v9/phases/$PHASE_B/freeze" \
  -H "Authorization: Bearer $TOKEN" | python3 -c "
import sys,json; d=json.load(sys.stdin)
print(f'PASS: Phase B frozen, {d.get(\"bookmarked_molecule_count\",0)} bookmarked')
"

# Creer Phase C
PHASE_C=$(curl -s -X POST "http://localhost:3000/api/v9/campaigns/$CAMPAIGN/phases" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"type": "lead_optimization"}' \
  | python3 -c "import sys,json; print(json.load(sys.stdin)['id'])")
echo "Phase C: $PHASE_C"

# Importer les bookmarks de Phase B
IMPORT_C=$(curl -s -X POST "http://localhost:3000/api/v9/phases/$PHASE_C/runs" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d "{
    \"type\": \"import\",
    \"config\": {\"source\": \"phase_selection\", \"source_phase_id\": \"$PHASE_B\"}
  }" | python3 -c "import sys,json; d=json.load(sys.stdin); print(d['id'])")

# Attendre
for i in $(seq 1 30); do
  STATUS=$(curl -s "http://localhost:3000/api/v9/runs/$IMPORT_C" \
    -H "Authorization: Bearer $TOKEN" | python3 -c "import sys,json; print(json.load(sys.stdin)['status'])")
  [ "$STATUS" = "completed" ] && break
  [ "$STATUS" = "failed" ] && break
  sleep 2
done
echo "Phase C import: $STATUS"
```

---

### Etape 15 : Validation finale — etat complet

```bash
curl -s "http://localhost:3000/api/v9/projects/$PROJECT" \
  -H "Authorization: Bearer $TOKEN" | python3 -c "
import sys,json
p = json.load(sys.stdin)
print(f'Project: {p[\"name\"]}')
print(f'Target: {p.get(\"target_name\",\"NOT SET\")} ({p.get(\"target_pdb_id\",\"?\")})')
print(f'Campaigns: {len(p[\"campaigns\"])}')
for c in p['campaigns']:
    print(f'  Campaign: {c[\"name\"]}')
    for ph in c['phases']:
        print(f'    Phase {ph[\"type\"]}: status={ph[\"status\"]}')
"
```

---

## 4. Regles metier a respecter

### Regles bloquantes (le frontend les applique)

| # | Regle | Controle | Consequence si viole |
|---|-------|----------|---------------------|
| R1 | Target doit etre configure avant de creer une campagne | `project.target_preview != null` | Bouton "Create Campaign" cache |
| R2 | Phase frozen → pas de nouveau run | `phase.status != 'frozen'` | Bouton "New Run" desactive |
| R3 | Phase frozen → pas de bookmark/annotation | `isFrozen` check | Callbacks a undefined |
| R4 | Run calcul necessite selection de molecules | `input_molecule_ids.length > 0` | Bouton "Launch" desactive |
| R5 | Run calcul necessite >= 1 calculation_type | `calculation_types.length > 0` | Bouton "Launch" desactive |
| R6 | 1 seul run actif a la fois (par user) | `activeRun != null` | Toast warning, bloque |
| R7 | Import necessite source de donnees | `smiles_list` ou `source=phase_selection` | Bouton "Launch" desactive |
| R8 | Phase B doit exister avant "Send to next phase" | `phases[myIdx+1]` | Toast warning |
| R9 | 1 phase par type par campagne | Backend 409 | Erreur |
| R10 | Bookmark necessite selection | `selectedCount > 0` | Bouton desactive |
| R11 | Unfreeze avec downstream runs → warning | `hasDownstreamRuns` check | Dialog rouge |

### Regles fonctionnelles (logique metier)

| # | Regle | Detail |
|---|-------|--------|
| F1 | 1 molecule = 1 ligne | Deduplique par canonical_smiles |
| F2 | Run calcul ajoute des colonnes | Pas des lignes (sauf import/generation) |
| F3 | Toutes les colonnes avec donnees visibles par defaut | `detectAvailableColumns()` |
| F4 | Flatten recursif des properties | `deepFlatten()` avec aliasing |
| F5 | Color scale en percentile | Pas de seuils fixes |
| F6 | Progression sequentielle : A → B → C | Pas de saut de phase |
| F7 | Freeze verrouille la selection bookmarkee | Input pour phase suivante |
| F8 | Run logs visibles dans RunHistory | `run.logs[]` avec level/message/timestamp |

---

## 5. Checklist par page

### ProjectListPage (`/welcome`)

- [ ] Liste des projets visible
- [ ] KPI par projet (total molecules, runs, phases)
- [ ] Mini 3D protein viewer avec couleur par projet
- [ ] Nom du projet au-dessus de la 3D
- [ ] Bouton "New Project"
- [ ] Polices >= 12px (text-xs minimum)

### ProjectHome (`/project/:id`)

- [ ] Target Setup card si `target_preview == null`
- [ ] Funnel A → B → C visible
- [ ] Stats par phase (molecule count, bookmarked, runs)
- [ ] Bouton "Create Campaign" visible SEULEMENT si target configure
- [ ] Navigation vers PhaseDashboard au clic sur une phase

### PhaseDashboard (`/project/:id/phase/:phaseId`)

- [ ] Breadcrumb correct (Project > Campaign > Phase)
- [ ] PhaseHeader avec stats (total, bookmarked, filtered)
- [ ] Bouton "New Run" actif (si phase non frozen et pas de run en cours)
- [ ] MoleculeTable avec colonnes detectees automatiquement
- [ ] Colonnes numeriques avec color scale
- [ ] Colonnes booleennes avec checkmark/X
- [ ] SelectionToolbar visible quand selection > 0
- [ ] RunHistory avec logs, type label, config summary
- [ ] Terminal logs agrege
- [ ] Freeze/Unfreeze dialog fonctionnel
- [ ] "Send to Next Phase" visible si bookmarks > 0 et phase non frozen
- [ ] MoleculeDetailPanel (clic sur molecule)
- [ ] Filtres type-aware (number: min/max, boolean: Yes/No/All)
- [ ] Polices >= 12px partout

### RunCreator (modal)

- [ ] Step 1 : choix type (Import / Calculation / Generation)
- [ ] Step 2 : config type-aware
- [ ] Step 3 : confirmation avec estimation
- [ ] Bouton "Launch" desactive si contraintes non satisfaites
- [ ] Message si molecules non selectionnees pour calculation/generation

---

## 6. Commandes de test

### Script complet automatise

Tout le scenario ci-dessus peut etre execute d'un bloc :

```bash
bash docs/e2e_test.sh
```

> A creer : `docs/e2e_test.sh` avec toutes les etapes ci-dessus en script executable.

### Test rapide (smoke test)

```bash
# Verification minimale : API up, auth OK, un projet existe
curl -sf http://localhost:3000/api/v9/health && echo "API: OK" || echo "API: FAIL"
TOKEN=$(curl -s 'https://webkntghfzscrnuixfba.supabase.co/auth/v1/token?grant_type=password' \
  -H 'apikey: sb_publishable_Zy18OaVJ4k_DgaafyiYaZA_9jB0ATmY' \
  -H 'Content-Type: application/json' \
  -d '{"email":"anthony.boisbouvier@gmail.com","password":"test1234"}' \
  | python3 -c "import sys,json; print(json.load(sys.stdin)['access_token'])")
curl -sf http://localhost:3000/api/v9/projects -H "Authorization: Bearer $TOKEN" | python3 -c "
import sys,json; print(f'Projects: {len(json.load(sys.stdin))}')
" && echo "Auth: OK" || echo "Auth: FAIL"
```

---

## 7. Maintenance

### Quand mettre a jour ce document

- Ajout d'une nouvelle page/composant → ajouter a la checklist par page
- Nouvelle regle metier → ajouter aux tables R# ou F#
- Nouveau type de run → ajouter au scenario
- Changement de schema API → mettre a jour les commandes curl
- Nouveau calcul disponible → ajouter a la liste des calculation_types

### Convention de versioning

```
Version X.Y :
  X = changement majeur du scenario (nouvelles etapes)
  Y = ajustements mineurs (nouvelles validations, corrections)
```

### Fichiers lies

- `docs/CDC_V9_FINAL.md` — Source of truth specs
- `docs/e2e_test.sh` — Script executable du scenario (a creer)
- `tests/` — Tests unitaires backend
