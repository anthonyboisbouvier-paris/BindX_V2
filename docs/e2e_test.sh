#!/bin/bash
# BindX V9 — E2E Test Script
# Usage: bash docs/e2e_test.sh
# Follows E2E_TEST_PROTOCOL.md step by step

set -e

# Colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

pass() { echo -e "${GREEN}PASS${NC}: $1"; }
fail() { echo -e "${RED}FAIL${NC}: $1"; FAILURES=$((FAILURES+1)); }
info() { echo -e "${CYAN}INFO${NC}: $1"; }
step() { echo -e "\n${YELLOW}=== STEP $1 ===${NC}"; }

FAILURES=0
BASE_URL="http://localhost:3000"
SUPABASE_URL="https://webkntghfzscrnuixfba.supabase.co"
ANON_KEY="sb_publishable_Zy18OaVJ4k_DgaafyiYaZA_9jB0ATmY"

# ---------------------------------------------------------------------------
# PRE-REQUIS
# ---------------------------------------------------------------------------
step "0: Pre-requisites"

# Check Docker containers
info "Checking Docker containers..."
docker-compose ps --format json 2>/dev/null | python3 -c "
import sys,json
lines = sys.stdin.read().strip().split('\n')
for line in lines:
    c = json.loads(line)
    name = c.get('Name','?')
    state = c.get('State','?')
    if 'healthy' not in c.get('Health','') and state != 'running':
        print(f'  WARNING: {name} is {state}')
    else:
        print(f'  OK: {name}')
" 2>/dev/null || info "Could not parse container status"

# Check API health
HTTP=$(curl -s -o /dev/null -w "%{http_code}" "$BASE_URL/api/v9/health")
[ "$HTTP" = "200" ] && pass "API health check" || fail "API not responding (HTTP $HTTP)"

# Get auth token
info "Authenticating..."
TOKEN=$(curl -s "$SUPABASE_URL/auth/v1/token?grant_type=password" \
  -H "apikey: $ANON_KEY" \
  -H "Content-Type: application/json" \
  -d '{"email":"anthony.boisbouvier@gmail.com","password":"test1234"}' \
  | python3 -c "import sys,json; print(json.load(sys.stdin)['access_token'])" 2>/dev/null)

[ ${#TOKEN} -gt 100 ] && pass "Authentication (token ${#TOKEN} chars)" || { fail "Authentication failed"; exit 1; }

AUTH="-H \"Authorization: Bearer $TOKEN\""

# Helper function for API calls
api() {
  local method="$1" path="$2" data="$3"
  if [ -n "$data" ]; then
    curl -s -X "$method" "$BASE_URL$path" \
      -H "Authorization: Bearer $TOKEN" \
      -H "Content-Type: application/json" \
      -d "$data"
  else
    curl -s -X "$method" "$BASE_URL$path" \
      -H "Authorization: Bearer $TOKEN"
  fi
}

api_code() {
  local method="$1" path="$2" data="$3"
  if [ -n "$data" ]; then
    curl -s -o /dev/null -w "%{http_code}" -X "$method" "$BASE_URL$path" \
      -H "Authorization: Bearer $TOKEN" \
      -H "Content-Type: application/json" \
      -d "$data"
  else
    curl -s -o /dev/null -w "%{http_code}" -X "$method" "$BASE_URL$path" \
      -H "Authorization: Bearer $TOKEN"
  fi
}

wait_run() {
  local run_id="$1" max_wait="${2:-60}"
  for i in $(seq 1 $max_wait); do
    local status=$(api GET "/api/v9/runs/$run_id" | python3 -c "import sys,json; print(json.load(sys.stdin)['status'])" 2>/dev/null)
    if [ "$status" = "completed" ]; then echo "completed"; return 0; fi
    if [ "$status" = "failed" ]; then echo "failed"; return 1; fi
    sleep 2
  done
  echo "timeout"
  return 1
}

# ---------------------------------------------------------------------------
# STEP 1: Create project
# ---------------------------------------------------------------------------
step "1: Create project"

PROJECT=$(api POST "/api/v9/projects" '{"name":"E2E Test '$(date +%H%M%S)'","description":"Automated E2E test"}' \
  | python3 -c "import sys,json; print(json.load(sys.stdin)['id'])" 2>/dev/null)
[ -n "$PROJECT" ] && pass "Project created: $PROJECT" || { fail "Project creation"; exit 1; }

# Validate project state
api GET "/api/v9/projects/$PROJECT" | python3 -c "
import sys,json
p = json.load(sys.stdin)
errors = []
if p['status'] != 'active': errors.append(f'status={p[\"status\"]} (expected active)')
if len(p['campaigns']) < 1: errors.append('no default campaign')
if p.get('target_preview') is not None: errors.append('target_preview should be None')
if errors:
    print('FAIL: ' + ', '.join(errors))
    sys.exit(1)
print(f'PASS: Project state OK (campaign={p[\"campaigns\"][0][\"id\"]})')
" && pass "Project state valid" || fail "Project state invalid"

CAMPAIGN=$(api GET "/api/v9/projects/$PROJECT" | python3 -c "import sys,json; print(json.load(sys.stdin)['campaigns'][0]['id'])")

# ---------------------------------------------------------------------------
# STEP 2a: Resolve target structure (REAL pipeline)
# ---------------------------------------------------------------------------
step "2a: Resolve target structure (UniProt → PDB → Pockets)"

# Rule R1: Target must be configured before creating campaigns
# This calls the REAL pipeline: UniProt → RCSB PDB → P2Rank pocket detection
info "Calling /api/preview-target with UniProt P00533 (EGFR)..."
info "This calls external APIs (RCSB, UniProt, ChEMBL, P2Rank) — may take 30-120s"

PREVIEW=$(curl -s -X POST "$BASE_URL/api/preview-target" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"uniprot_id": "P00533"}' \
  --max-time 180)

# Validate structure was resolved
echo "$PREVIEW" | python3 -c "
import sys,json
try:
    d = json.load(sys.stdin)
except:
    print('FAIL: Could not parse preview-target response')
    sys.exit(1)

errors = []

# Must have a resolved structure
structure = d.get('structure', {})
if not structure:
    errors.append('No structure resolved')
elif structure.get('source') == 'none':
    errors.append(f'Structure source is none: {structure.get(\"explanation\",\"?\")}')
else:
    print(f'  Structure: {structure.get(\"pdb_id\",\"?\")} ({structure.get(\"source\",\"?\")}, {structure.get(\"resolution\",\"?\")}A)')

# Must have pockets with center coordinates
pockets = d.get('pockets', [])
if not pockets:
    errors.append('No pockets detected')
else:
    print(f'  Pockets: {len(pockets)} detected')
    for i, p in enumerate(pockets[:3]):
        center = p.get('center', [])
        method = p.get('method', '?')
        prob = p.get('probability', 0)
        print(f'    Pocket {i+1}: method={method}, probability={prob:.2f}, center={center}')
        if not center or len(center) != 3:
            errors.append(f'Pocket {i+1} has no valid center: {center}')

# Must have protein name
protein_name = d.get('protein_name', '')
if not protein_name:
    errors.append('No protein name')
else:
    print(f'  Protein: {protein_name}')

if errors:
    for e in errors: print(f'  ERROR: {e}')
    sys.exit(1)
" && pass "Target preview resolved" || fail "Target preview"

# Extract key data for project update
TARGET_DATA=$(echo "$PREVIEW" | python3 -c "
import sys,json
d = json.load(sys.stdin)
s = d.get('structure', {})
pockets = d.get('pockets', [])
# Build the project update payload
update = {
    'target_input_type': 'uniprot',
    'target_input_value': 'P00533',
    'target_name': d.get('protein_name', 'EGFR'),
    'target_pdb_id': s.get('pdb_id'),
    'structure_source': s.get('source', 'experimental'),
    'structure_resolution': s.get('resolution'),
    'structure_method': s.get('method'),
    'cocrystal_ligand': s.get('ligand_id'),
    'target_preview': {
        'uniprot': d.get('uniprot', {'id': 'P00533'}),
        'structure': s,
        'structures': d.get('structures', []),
        'pockets': pockets,
        'selected_pocket_index': 0,  # Select first (best) pocket
        'chembl': d.get('chembl_info', {})
    },
    'pockets_detected': pockets,
    'chembl_actives_count': d.get('chembl_info', {}).get('n_actives', 0)
}
print(json.dumps(update))
")

# ---------------------------------------------------------------------------
# STEP 2b: Validate pocket has docking-ready center coordinates
# ---------------------------------------------------------------------------
step "2b: Validate pocket center for docking"

echo "$PREVIEW" | python3 -c "
import sys,json
d = json.load(sys.stdin)
pockets = d.get('pockets', [])
if not pockets:
    print('FAIL: No pockets to select')
    sys.exit(1)

# Pocket 1 (best) must have valid center
p = pockets[0]
center = p.get('center', [])
if len(center) != 3:
    print(f'FAIL: Pocket 1 center invalid: {center}')
    sys.exit(1)
x, y, z = center
if not all(isinstance(c, (int, float)) for c in [x, y, z]):
    print(f'FAIL: Pocket 1 center not numeric: {center}')
    sys.exit(1)

method = p.get('method', '?')
prob = p.get('probability', 0)
print(f'Selected pocket: method={method}, prob={prob:.2f}, center=[{x:.1f}, {y:.1f}, {z:.1f}]')
print(f'  → Docking box: center=({x:.1f}, {y:.1f}, {z:.1f}), size=22x22x22 A')
" && pass "Pocket center valid for docking" || fail "Pocket center validation"

# ---------------------------------------------------------------------------
# STEP 2c: Save target to project (like frontend 'Validate Target' click)
# ---------------------------------------------------------------------------
step "2c: Save target to project"

api PUT "/api/v9/projects/$PROJECT" "$TARGET_DATA" | python3 -c "
import sys,json
p = json.load(sys.stdin)
errors = []
if not p.get('target_preview'):
    errors.append('target_preview not saved')
if not p.get('target_name'):
    errors.append('target_name not saved')
if not p.get('target_pdb_id'):
    errors.append('target_pdb_id not saved')
if not p.get('pockets_detected'):
    errors.append('pockets_detected not saved')
else:
    # Verify selected pocket has center
    pockets = p['pockets_detected']
    if pockets and 'center' in pockets[0]:
        print(f'  Saved: {p[\"target_name\"]} / PDB {p[\"target_pdb_id\"]}')
        print(f'  Pockets: {len(pockets)}, selected pocket center: {pockets[0][\"center\"]}')
    else:
        errors.append('Saved pockets missing center data')

if errors:
    for e in errors: print(f'  ERROR: {e}')
    sys.exit(1)
" && pass "Target saved to project" || fail "Target save"

# ---------------------------------------------------------------------------
# STEP 3: Create Phase A (hit_discovery)
# ---------------------------------------------------------------------------
step "3: Create Phase A"

PHASE_A=$(api POST "/api/v9/campaigns/$CAMPAIGN/phases" '{"type":"hit_discovery"}' \
  | python3 -c "import sys,json; print(json.load(sys.stdin)['id'])" 2>/dev/null)
[ -n "$PHASE_A" ] && pass "Phase A created: $PHASE_A" || { fail "Phase A creation"; exit 1; }

# Rule R9: Cannot create duplicate phase type
HTTP=$(api_code POST "/api/v9/campaigns/$CAMPAIGN/phases" '{"type":"hit_discovery"}')
[ "$HTTP" = "409" ] && pass "Rule R9: Duplicate phase type blocked (409)" || fail "Rule R9: Expected 409, got $HTTP"

# ---------------------------------------------------------------------------
# STEP 4: Run Import (SMILES)
# ---------------------------------------------------------------------------
step "4: Run Import"

IMPORT_RUN=$(api POST "/api/v9/phases/$PHASE_A/runs" '{
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
}' | python3 -c "import sys,json; print(json.load(sys.stdin)['id'])" 2>/dev/null)
[ -n "$IMPORT_RUN" ] && pass "Import run created: $IMPORT_RUN" || { fail "Import run creation"; exit 1; }

info "Waiting for import to complete..."
IMPORT_STATUS=$(wait_run "$IMPORT_RUN" 30)
[ "$IMPORT_STATUS" = "completed" ] && pass "Import completed" || fail "Import $IMPORT_STATUS"

# Validate molecules imported
api GET "/api/v9/phases/$PHASE_A/molecules?limit=1" | python3 -c "
import sys,json
d = json.load(sys.stdin)
total = d['total']
if total >= 6:
    print(f'PASS: {total} molecules imported')
else:
    print(f'FAIL: Only {total} molecules')
    sys.exit(1)
" || fail "Molecule count check"

# Validate run has logs
api GET "/api/v9/runs/$IMPORT_RUN" | python3 -c "
import sys,json
r = json.load(sys.stdin)
logs = r.get('logs', [])
if len(logs) >= 2:
    print(f'PASS: Import run has {len(logs)} logs')
else:
    print(f'FAIL: Only {len(logs)} logs')
    sys.exit(1)
" || fail "Import run logs"

# ---------------------------------------------------------------------------
# STEP 5: Rule R7 — Import requires data source
# ---------------------------------------------------------------------------
step "5: Validate import constraints"

HTTP=$(api_code POST "/api/v9/phases/$PHASE_A/runs" '{"type":"import","config":{}}')
[ "$HTTP" = "422" ] && pass "Rule R7: Empty import blocked (422)" || fail "Rule R7: Expected 422, got $HTTP"

# ---------------------------------------------------------------------------
# STEP 6: Run Calculation (requires molecule selection)
# ---------------------------------------------------------------------------
step "6: Run Calculation"

# Get molecule IDs for selection
MOL_IDS=$(api GET "/api/v9/phases/$PHASE_A/molecules?limit=4" | python3 -c "
import sys,json
d = json.load(sys.stdin)
ids = [m['id'] for m in d['molecules'][:4]]
print(json.dumps(ids))
")
info "Selected molecules: $MOL_IDS"

# Rule R4: Calculation requires molecules
HTTP=$(api_code POST "/api/v9/phases/$PHASE_A/runs" '{"type":"calculation","calculation_types":["admet"]}')
[ "$HTTP" = "422" ] && pass "Rule R4: Calc without molecules blocked (422)" || fail "Rule R4: Expected 422, got $HTTP"

# Rule R5: Calculation requires calculation_types
HTTP=$(api_code POST "/api/v9/phases/$PHASE_A/runs" "{\"type\":\"calculation\",\"input_molecule_ids\":$MOL_IDS}")
[ "$HTTP" = "422" ] && pass "Rule R5: Calc without types blocked (422)" || fail "Rule R5: Expected 422, got $HTTP"

# Valid calculation run — includes docking to test pocket/structure integration
CALC_RUN=$(api POST "/api/v9/phases/$PHASE_A/runs" "{
  \"type\": \"calculation\",
  \"calculation_types\": [\"docking\", \"admet\", \"scoring\"],
  \"input_molecule_ids\": $MOL_IDS
}" | python3 -c "import sys,json; print(json.load(sys.stdin)['id'])" 2>/dev/null)
[ -n "$CALC_RUN" ] && pass "Calculation run created (docking+admet+scoring): $CALC_RUN" || { fail "Calc run creation"; exit 1; }

info "Waiting for calculation to complete..."
CALC_STATUS=$(wait_run "$CALC_RUN" 120)
[ "$CALC_STATUS" = "completed" ] && pass "Calculation completed" || fail "Calculation $CALC_STATUS"

# Validate properties on molecules
api GET "/api/v9/phases/$PHASE_A/molecules?limit=20" | python3 -c "
import sys,json
d = json.load(sys.stdin)
with_props = [m for m in d['molecules'] if m.get('properties')]
if len(with_props) >= 4:
    print(f'PASS: {len(with_props)} molecules have properties')
else:
    print(f'FAIL: Only {len(with_props)} molecules have properties')
    sys.exit(1)
" || fail "Molecule properties"

# ---------------------------------------------------------------------------
# STEP 6b: Validate docking status
# ---------------------------------------------------------------------------
step "6b: Validate docking integration"

# Docking without GPU produces 'pending_gpu' status — this is expected
# The key validation is that the docking_status property EXISTS (structure+pocket were accessible)
api GET "/api/v9/phases/$PHASE_A/molecules?limit=20" | python3 -c "
import sys,json
d = json.load(sys.stdin)
with_props = [m for m in d['molecules'] if m.get('properties')]
docking_ok = 0
for m in with_props:
    props = m['properties']
    # Check docking_status exists (proves structure+pocket chain worked)
    ds = props.get('docking_status')
    if ds:
        docking_ok += 1
        status = ds.get('status', '?') if isinstance(ds, dict) else ds
        if docking_ok == 1:
            print(f'  Docking status: {status}')
    # Check ADMET properties exist
    admet = props.get('admet') or props.get('physicochemical')
    if admet and docking_ok == 1:
        if isinstance(admet, dict):
            top_keys = list(admet.keys())[:5]
            print(f'  ADMET keys: {top_keys}')

if docking_ok >= 4:
    print(f'PASS: {docking_ok} molecules have docking_status')
else:
    print(f'FAIL: Only {docking_ok} molecules have docking_status (expected >=4)')
    sys.exit(1)
" && pass "Docking integration validated" || fail "Docking integration"

# Validate run has calculation_types including docking
api GET "/api/v9/runs/$CALC_RUN" | python3 -c "
import sys,json
r = json.load(sys.stdin)
ct = r.get('calculation_types', [])
if 'docking' not in ct:
    print(f'FAIL: calculation_types missing docking: {ct}')
    sys.exit(1)
if 'admet' not in ct:
    print(f'FAIL: calculation_types missing admet: {ct}')
    sys.exit(1)
print(f'PASS: Run has calculation_types: {ct}')
" && pass "Run calculation_types correct" || fail "Run calculation_types"

# ---------------------------------------------------------------------------
# STEP 7: Validate run structure (logs, types, config)
# ---------------------------------------------------------------------------
step "7: Validate run structure"

api GET "/api/v9/phases/$PHASE_A/runs" | python3 -c "
import sys,json
runs = json.load(sys.stdin)
errors = []
for r in runs:
    if 'logs' not in r: errors.append(f'{r[\"type\"]} run missing logs field')
    if r['type'] == 'calculation':
        if not r.get('calculation_types'): errors.append('calc run missing calculation_types')
    if 'created_at' not in r: errors.append(f'{r[\"type\"]} missing created_at')
if errors:
    for e in errors: print(f'FAIL: {e}')
    sys.exit(1)
print(f'PASS: All {len(runs)} runs have correct structure')
" || fail "Run structure"

# ---------------------------------------------------------------------------
# STEP 8: Bookmark molecules
# ---------------------------------------------------------------------------
step "8: Bookmark molecules"

api POST "/api/v9/phases/$PHASE_A/molecules/bookmark-batch" \
  "{\"molecule_ids\": $MOL_IDS, \"bookmarked\": true}" | python3 -c "
import sys,json
d = json.load(sys.stdin)
print(f'Bookmarked: {d.get(\"updated\",\"?\")}')
" > /dev/null
pass "Bookmark request sent"

# Validate stats
api GET "/api/v9/phases/$PHASE_A/molecules/stats" | python3 -c "
import sys,json
s = json.load(sys.stdin)
if s['bookmarked'] >= 4:
    print(f'PASS: {s[\"bookmarked\"]} bookmarked / {s[\"total\"]} total')
else:
    print(f'FAIL: Only {s[\"bookmarked\"]} bookmarked')
    sys.exit(1)
" || fail "Bookmark stats"

# ---------------------------------------------------------------------------
# STEP 9: Freeze Phase A
# ---------------------------------------------------------------------------
step "9: Freeze Phase A"

api POST "/api/v9/phases/$PHASE_A/freeze" | python3 -c "
import sys,json
d = json.load(sys.stdin)
if d['status'] == 'frozen':
    print(f'PASS: Frozen with {d.get(\"bookmarked_molecule_count\",0)} bookmarked')
else:
    print(f'FAIL: Status is {d[\"status\"]}')
    sys.exit(1)
" || fail "Freeze Phase A"

# ---------------------------------------------------------------------------
# STEP 10: Verify frozen phase blocks actions
# ---------------------------------------------------------------------------
step "10: Verify frozen constraints"

# Rule R2: Cannot create run on frozen phase
HTTP=$(api_code POST "/api/v9/phases/$PHASE_A/runs" '{"type":"import","config":{"smiles_list":["CCCC"]}}')
[ "$HTTP" = "400" ] && pass "Rule R2: Run on frozen phase blocked (400)" || fail "Rule R2: Expected 400, got $HTTP"

# Rule R2: Cannot import to frozen phase
HTTP=$(api_code POST "/api/v9/phases/$PHASE_A/runs/import-database" '{"databases":["chembl"],"uniprot_id":"P00533"}')
[ "$HTTP" = "400" ] && pass "Rule R2: Import to frozen phase blocked (400)" || fail "Rule R2 import: Expected 400, got $HTTP"

# ---------------------------------------------------------------------------
# STEP 11: Create Phase B (BEFORE sending molecules)
# ---------------------------------------------------------------------------
step "11: Create Phase B"

# Rule F6: Sequential progression A → B → C
PHASE_B=$(api POST "/api/v9/campaigns/$CAMPAIGN/phases" '{"type":"hit_to_lead"}' \
  | python3 -c "import sys,json; print(json.load(sys.stdin)['id'])" 2>/dev/null)
[ -n "$PHASE_B" ] && pass "Phase B created: $PHASE_B" || { fail "Phase B creation"; exit 1; }

# ---------------------------------------------------------------------------
# STEP 12: Send bookmarks from Phase A to Phase B
# ---------------------------------------------------------------------------
step "12: Phase A → Phase B transition"

NEXT_RUN=$(api POST "/api/v9/phases/$PHASE_B/runs" "{
  \"type\": \"import\",
  \"config\": {\"source\": \"phase_selection\", \"source_phase_id\": \"$PHASE_A\"}
}" | python3 -c "import sys,json; print(json.load(sys.stdin)['id'])" 2>/dev/null)
[ -n "$NEXT_RUN" ] && pass "Phase transition import created" || { fail "Phase transition"; exit 1; }

info "Waiting for phase transition import..."
NEXT_STATUS=$(wait_run "$NEXT_RUN" 30)
[ "$NEXT_STATUS" = "completed" ] && pass "Phase transition completed" || fail "Phase transition $NEXT_STATUS"

# Validate Phase B has molecules from Phase A bookmarks
api GET "/api/v9/phases/$PHASE_B/molecules?limit=1" | python3 -c "
import sys,json
d = json.load(sys.stdin)
if d['total'] >= 4:
    print(f'PASS: Phase B has {d[\"total\"]} molecules')
else:
    print(f'FAIL: Phase B has only {d[\"total\"]} molecules')
    sys.exit(1)
" || fail "Phase B molecule count"

# ---------------------------------------------------------------------------
# STEP 13: Run calculation on Phase B
# ---------------------------------------------------------------------------
step "13: Phase B calculation"

MOL_IDS_B=$(api GET "/api/v9/phases/$PHASE_B/molecules?limit=10" | python3 -c "
import sys,json; print(json.dumps([m['id'] for m in json.load(sys.stdin)['molecules']]))")

CALC_B=$(api POST "/api/v9/phases/$PHASE_B/runs" "{
  \"type\": \"calculation\",
  \"calculation_types\": [\"admet\", \"safety\", \"scoring\"],
  \"input_molecule_ids\": $MOL_IDS_B
}" | python3 -c "import sys,json; print(json.load(sys.stdin)['id'])" 2>/dev/null)
[ -n "$CALC_B" ] && pass "Phase B calc run created" || { fail "Phase B calc"; exit 1; }

info "Waiting for Phase B calculation..."
CALC_B_STATUS=$(wait_run "$CALC_B" 60)
[ "$CALC_B_STATUS" = "completed" ] && pass "Phase B calculation completed" || fail "Phase B calc $CALC_B_STATUS"

# ---------------------------------------------------------------------------
# STEP 14: Phase B → Phase C (full progression)
# ---------------------------------------------------------------------------
step "14: Phase B → C progression"

# Bookmark 2 molecules in Phase B
MOL_2=$(api GET "/api/v9/phases/$PHASE_B/molecules?limit=2" | python3 -c "
import sys,json; print(json.dumps([m['id'] for m in json.load(sys.stdin)['molecules'][:2]]))")
api POST "/api/v9/phases/$PHASE_B/molecules/bookmark-batch" "{\"molecule_ids\":$MOL_2,\"bookmarked\":true}" > /dev/null
pass "Bookmarked 2 molecules in Phase B"

# Freeze Phase B
api POST "/api/v9/phases/$PHASE_B/freeze" | python3 -c "
import sys,json; d=json.load(sys.stdin)
print(f'Frozen with {d.get(\"bookmarked_molecule_count\",0)} bookmarked')
" > /dev/null
pass "Phase B frozen"

# Create Phase C
PHASE_C=$(api POST "/api/v9/campaigns/$CAMPAIGN/phases" '{"type":"lead_optimization"}' \
  | python3 -c "import sys,json; print(json.load(sys.stdin)['id'])" 2>/dev/null)
[ -n "$PHASE_C" ] && pass "Phase C created: $PHASE_C" || { fail "Phase C creation"; exit 1; }

# Import from Phase B
IMPORT_C=$(api POST "/api/v9/phases/$PHASE_C/runs" "{
  \"type\":\"import\",
  \"config\":{\"source\":\"phase_selection\",\"source_phase_id\":\"$PHASE_B\"}
}" | python3 -c "import sys,json; print(json.load(sys.stdin)['id'])" 2>/dev/null)

info "Waiting for Phase C import..."
IMPORT_C_STATUS=$(wait_run "$IMPORT_C" 30)
[ "$IMPORT_C_STATUS" = "completed" ] && pass "Phase C import completed" || fail "Phase C import $IMPORT_C_STATUS"

api GET "/api/v9/phases/$PHASE_C/molecules?limit=1" | python3 -c "
import sys,json
d = json.load(sys.stdin)
if d['total'] >= 2:
    print(f'PASS: Phase C has {d[\"total\"]} molecules')
else:
    print(f'FAIL: Phase C has only {d[\"total\"]} molecules')
    sys.exit(1)
" || fail "Phase C molecule count"

# ---------------------------------------------------------------------------
# STEP 15: Final validation — full project state
# ---------------------------------------------------------------------------
step "15: Final validation"

api GET "/api/v9/projects/$PROJECT" | python3 -c "
import sys,json
p = json.load(sys.stdin)
print(f'Project: {p[\"name\"]}')
print(f'Target:  {p.get(\"target_name\",\"NOT SET\")} ({p.get(\"target_pdb_id\",\"?\")})')
phases = []
for c in p['campaigns']:
    for ph in c['phases']:
        phases.append(ph)
        print(f'  Phase {ph[\"type\"]:20s} status={ph[\"status\"]}')

errors = []
if len(phases) != 3: errors.append(f'Expected 3 phases, got {len(phases)}')
types = [ph['type'] for ph in phases]
if 'hit_discovery' not in types: errors.append('Missing hit_discovery')
if 'hit_to_lead' not in types: errors.append('Missing hit_to_lead')
if 'lead_optimization' not in types: errors.append('Missing lead_optimization')

hd = next((ph for ph in phases if ph['type'] == 'hit_discovery'), None)
htl = next((ph for ph in phases if ph['type'] == 'hit_to_lead'), None)
if hd and hd['status'] != 'frozen': errors.append(f'Phase A should be frozen, is {hd[\"status\"]}')
if htl and htl['status'] != 'frozen': errors.append(f'Phase B should be frozen, is {htl[\"status\"]}')

if errors:
    for e in errors: print(f'  ERROR: {e}')
    sys.exit(1)
print('All phases present and in correct state')
"
[ $? -eq 0 ] && pass "Final state validation" || fail "Final state validation"

# ---------------------------------------------------------------------------
# SUMMARY
# ---------------------------------------------------------------------------
echo ""
echo "==========================================="
if [ $FAILURES -eq 0 ]; then
  echo -e "${GREEN}ALL TESTS PASSED${NC}"
else
  echo -e "${RED}$FAILURES TEST(S) FAILED${NC}"
fi
echo "==========================================="
echo "Project ID: $PROJECT"
echo "Phase A:    $PHASE_A"
echo "Phase B:    $PHASE_B"
echo "Phase C:    $PHASE_C"
echo "==========================================="

exit $FAILURES
