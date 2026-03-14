#!/usr/bin/env python3
"""
BindX V2 — 100K Molecule Stress Test
=====================================
Tests scalability of FastAPI backend + Supabase PostgreSQL with 100K molecules.

Steps:
  1. Authenticate → get JWT token
  2. Create test project + campaign + phase via API
  3. Bulk insert 100K molecules directly into PostgreSQL
  4. Run performance tests (list, sort, stats, concurrent)
  5. Cleanup
  6. Print full report
"""

import json
import random
import statistics
import string
import time
import uuid
import urllib.request
import urllib.parse
import urllib.error
import asyncio
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone

import psycopg2
import psycopg2.extras

# ============================================================================
# Configuration
# ============================================================================

BACKEND_URL = "http://localhost:8000"

# Supabase Auth
SUPABASE_URL = "https://webkntghfzscrnuixfba.supabase.co"
SUPABASE_ANON_KEY = "sb_publishable_Zy18OaVJ4k_DgaafyiYaZA_9jB0ATmY"
AUTH_EMAIL = "anthony.boisbouvier@gmail.com"
AUTH_PASSWORD = "test1234"

# Direct PostgreSQL connection (bypass asyncpg driver prefix)
PG_DSN = (
    "postgresql://postgres.webkntghfzscrnuixfba:n1qsKxvAWJV954vR"
    "@aws-0-us-west-2.pooler.supabase.com:5432/postgres"
)

TOTAL_MOLECULES = 100_000
BATCH_SIZE = 1000
PASS_THRESHOLD_MS = 2000  # 2 seconds = PASS

# SMILES building blocks for realistic molecule generation
SMILES_CORES = [
    "c1ccccc1", "c1ccncc1", "C1CCCCC1", "c1ccoc1",
    "c1ccsc1", "C1CCNCC1", "c1ccc2ccccc2c1", "c1cnc2ccccc2n1",
    "c1ccc(cc1)", "C1CC1", "C1CCC1", "c1cc[nH]c1",
    "c1cccc2cccccc12", "c1cncc1", "C1CCOCC1",
]

SMILES_SUBSTITUENTS = [
    "C", "CC", "CCC", "CCCC", "OC", "NC", "FC", "ClC",
    "BrC", "CC(=O)", "C(=O)O", "C(=O)N", "S(=O)(=O)",
    "c1ccccc1", "CN", "CO", "CF",
]


# ============================================================================
# Utility helpers
# ============================================================================

def log(msg: str, level: str = "INFO"):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] [{level}] {msg}")


def sep(title: str = ""):
    line = "=" * 70
    if title:
        print(f"\n{line}")
        print(f"  {title}")
        print(f"{line}")
    else:
        print(line)


def generate_smiles(i: int) -> str:
    """Generate a unique SMILES-like string. Simple patterns, no RDKit needed."""
    core = random.choice(SMILES_CORES)
    n_subs = random.randint(1, 3)
    subs = [random.choice(SMILES_SUBSTITUENTS) for _ in range(n_subs)]
    # Tag with index to guarantee uniqueness (canonical_smiles must be unique per phase)
    tag = f"[{i:06d}]"  # non-standard atom map-like tag for uniqueness
    return f"{core}{''.join(subs)}.{tag}" if subs else f"{core}{tag}"


def generate_molecule_batch(
    phase_id: str,
    start_idx: int,
    batch_size: int,
) -> list[tuple]:
    """Generate a batch of molecule rows as tuples for psycopg2 executemany."""
    rows = []
    now = datetime.now(timezone.utc)
    for i in range(start_idx, start_idx + batch_size):
        mol_id = str(uuid.uuid4())
        smiles = generate_smiles(i)
        name = f"MOL_{i:06d}"
        bookmarked = random.random() < 0.05  # 5% bookmarked
        ai_generated = random.random() < 0.1  # 10% AI-generated
        generation_level = random.randint(0, 2)
        rows.append((
            mol_id,
            phase_id,
            smiles,
            smiles,  # canonical_smiles = smiles (simplified)
            name,
            bookmarked,
            ai_generated,
            generation_level,
            now,
        ))
    return rows


def generate_property_batch(molecule_ids: list[str]) -> list[tuple]:
    """Generate ADMET-like properties for a batch of molecules."""
    props = []
    now = datetime.now(timezone.utc)
    prop_defs = [
        ("docking_score", lambda: round(random.uniform(-12.0, -2.0), 3)),
        ("cnn_score", lambda: round(random.uniform(0.0, 1.0), 4)),
        ("cnn_affinity", lambda: round(random.uniform(0.5, 12.0), 3)),
        ("composite_score", lambda: round(random.uniform(0.0, 1.0), 4)),
        ("QED", lambda: round(random.uniform(0.0, 1.0), 4)),
        ("logP", lambda: round(random.uniform(-2.0, 6.0), 3)),
        ("MW", lambda: round(random.uniform(150.0, 750.0), 2)),
        ("TPSA", lambda: round(random.uniform(0.0, 150.0), 2)),
        ("HBD", lambda: random.randint(0, 5)),
        ("HBA", lambda: random.randint(0, 10)),
        ("safety_color_code", lambda: random.choice(["green", "yellow", "red"])),
    ]
    for mol_id in molecule_ids:
        for prop_name, val_fn in prop_defs:
            prop_id = str(uuid.uuid4())
            val = val_fn()
            prop_value = json.dumps({"value": val})
            props.append((prop_id, mol_id, prop_name, prop_value, now))
    return props


# ============================================================================
# HTTP helpers (pure stdlib)
# ============================================================================

def http_post(url: str, data: dict, headers: dict = None) -> tuple[int, dict]:
    body = json.dumps(data).encode("utf-8")
    h = {"Content-Type": "application/json"}
    if headers:
        h.update(headers)
    req = urllib.request.Request(url, data=body, headers=h, method="POST")
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            return resp.status, json.loads(resp.read())
    except urllib.error.HTTPError as e:
        return e.code, json.loads(e.read())


def http_get(url: str, headers: dict = None) -> tuple[int, dict, float]:
    """Returns (status, body, elapsed_ms)."""
    h = {}
    if headers:
        h.update(headers)
    req = urllib.request.Request(url, headers=h, method="GET")
    t0 = time.perf_counter()
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            body = json.loads(resp.read())
            elapsed = (time.perf_counter() - t0) * 1000
            return resp.status, body, elapsed
    except urllib.error.HTTPError as e:
        elapsed = (time.perf_counter() - t0) * 1000
        return e.code, json.loads(e.read() or b"{}"), elapsed


def timed_get(url: str, headers: dict) -> float:
    """Return response time in ms for a GET request."""
    _, _, elapsed = http_get(url, headers)
    return elapsed


# ============================================================================
# Step 1: Authentication
# ============================================================================

def authenticate() -> str:
    sep("STEP 1: Authentication")
    url = f"{SUPABASE_URL}/auth/v1/token?grant_type=password"
    headers = {
        "apikey": SUPABASE_ANON_KEY,
        "Content-Type": "application/json",
    }
    log(f"Authenticating as {AUTH_EMAIL}...")
    t0 = time.perf_counter()
    status, body = http_post(url, {"email": AUTH_EMAIL, "password": AUTH_PASSWORD}, headers)
    elapsed = (time.perf_counter() - t0) * 1000

    if status != 200 or "access_token" not in body:
        log(f"Auth failed: {status} — {body}", "ERROR")
        raise RuntimeError("Authentication failed")

    token = body["access_token"]
    log(f"Auth OK in {elapsed:.0f}ms — token: {token[:40]}...")
    return token


# ============================================================================
# Step 2: Create test project + campaign + phase
# ============================================================================

def create_test_hierarchy(token: str) -> tuple[str, str, str]:
    sep("STEP 2: Create project / campaign / phase")
    auth_headers = {"Authorization": f"Bearer {token}"}

    # Create project
    log("Creating project 'Stress Test 100K'...")
    status, body = http_post(
        f"{BACKEND_URL}/api/v9/projects",
        {"name": "Stress Test 100K", "target_pdb_id": "TEST", "target_input_value": ""},
        auth_headers,
    )
    if status not in (200, 201):
        log(f"Project creation failed: {status} — {body}", "ERROR")
        raise RuntimeError("Project creation failed")
    project_id = body["id"]
    log(f"Project created: {project_id}")

    # Create campaign
    log("Creating campaign...")
    status, body = http_post(
        f"{BACKEND_URL}/api/v9/projects/{project_id}/campaigns",
        {"name": "Stress Test Campaign"},
        auth_headers,
    )
    if status not in (200, 201):
        log(f"Campaign creation failed: {status} — {body}", "ERROR")
        raise RuntimeError("Campaign creation failed")
    campaign_id = body["id"]
    log(f"Campaign created: {campaign_id}")

    # Get auto-created phases (campaigns create a hit_discovery phase automatically?)
    # Actually, per the schema, phases are created separately
    # Create a phase
    log("Creating hit_discovery phase...")
    status, body = http_post(
        f"{BACKEND_URL}/api/v9/campaigns/{campaign_id}/phases",
        {"type": "hit_discovery"},
        auth_headers,
    )
    if status not in (200, 201):
        log(f"Phase creation failed: {status} — {body}", "ERROR")
        raise RuntimeError("Phase creation failed")
    phase_id = body["id"]
    log(f"Phase created: {phase_id}")

    return project_id, campaign_id, phase_id


# ============================================================================
# Step 3: Bulk insert 100K molecules via direct PostgreSQL
# ============================================================================

def bulk_insert_molecules(phase_id: str) -> dict:
    sep("STEP 3: Bulk Insert 100K Molecules")

    # Strip asyncpg driver prefix for psycopg2
    dsn = PG_DSN

    log(f"Connecting to PostgreSQL...")
    conn = psycopg2.connect(dsn)
    conn.autocommit = False
    cur = conn.cursor()

    log(f"Starting bulk insert of {TOTAL_MOLECULES:,} molecules in batches of {BATCH_SIZE:,}...")
    total_start = time.perf_counter()

    molecule_ids_all = []
    batch_times = []

    n_batches = TOTAL_MOLECULES // BATCH_SIZE
    molecules_inserted = 0

    mol_insert_sql = """
        INSERT INTO molecules
            (id, phase_id, smiles, canonical_smiles, name, bookmarked, ai_generated, generation_level, created_at)
        VALUES %s
        ON CONFLICT (phase_id, canonical_smiles) DO NOTHING
    """

    for batch_num in range(n_batches):
        start_idx = batch_num * BATCH_SIZE
        batch_start = time.perf_counter()

        rows = generate_molecule_batch(phase_id, start_idx, BATCH_SIZE)
        mol_ids_in_batch = [r[0] for r in rows]

        psycopg2.extras.execute_values(cur, mol_insert_sql, rows, page_size=500)
        conn.commit()

        inserted_this_batch = cur.rowcount
        molecules_inserted += inserted_this_batch
        molecule_ids_all.extend(mol_ids_in_batch)

        batch_elapsed = (time.perf_counter() - batch_start) * 1000
        batch_times.append(batch_elapsed)

        if (batch_num + 1) % 10 == 0:
            pct = ((batch_num + 1) / n_batches) * 100
            avg_batch = statistics.mean(batch_times[-10:])
            log(
                f"  Batch {batch_num+1}/{n_batches} ({pct:.0f}%) | "
                f"inserted {molecules_inserted:,} | "
                f"batch avg {avg_batch:.0f}ms"
            )

    total_mol_elapsed = time.perf_counter() - total_start
    log(f"Molecules inserted: {molecules_inserted:,} in {total_mol_elapsed:.1f}s")
    log(f"Average batch time: {statistics.mean(batch_times):.0f}ms")

    # Now insert properties for all molecules
    log(f"\nInserting properties ({len(molecule_ids_all):,} molecules x 11 props = {len(molecule_ids_all)*11:,} rows)...")
    prop_insert_sql = """
        INSERT INTO molecule_properties
            (id, molecule_id, property_name, property_value, created_at)
        VALUES %s
        ON CONFLICT (molecule_id, property_name, run_id) DO NOTHING
    """

    # Note: run_id is null for these test properties, constraint is (molecule_id, property_name, run_id)
    # We need to handle null run_id in the unique constraint
    prop_insert_sql_no_run = """
        INSERT INTO molecule_properties
            (id, molecule_id, property_name, property_value, created_at)
        VALUES %s
        ON CONFLICT DO NOTHING
    """

    prop_start = time.perf_counter()
    prop_batch_size = 500  # 500 molecules x 11 = 5500 rows per batch
    prop_batches_inserted = 0

    for i in range(0, len(molecule_ids_all), prop_batch_size):
        batch_ids = molecule_ids_all[i:i + prop_batch_size]
        prop_rows = generate_property_batch(batch_ids)
        psycopg2.extras.execute_values(cur, prop_insert_sql_no_run, prop_rows, page_size=2000)
        conn.commit()
        prop_batches_inserted += 1

        if prop_batches_inserted % 20 == 0:
            pct = (i / len(molecule_ids_all)) * 100
            log(f"  Properties batch {prop_batches_inserted} ({pct:.0f}%)")

    prop_elapsed = time.perf_counter() - prop_start
    log(f"Properties inserted in {prop_elapsed:.1f}s")

    # Verify counts
    cur.execute("SELECT COUNT(*) FROM molecules WHERE phase_id = %s", (phase_id,))
    mol_count = cur.fetchone()[0]
    cur.execute(
        "SELECT COUNT(*) FROM molecule_properties WHERE molecule_id IN "
        "(SELECT id FROM molecules WHERE phase_id = %s)",
        (phase_id,),
    )
    prop_count = cur.fetchone()[0]
    log(f"Verification: {mol_count:,} molecules, {prop_count:,} properties in DB")

    cur.close()
    conn.close()

    return {
        "molecules_inserted": mol_count,
        "properties_inserted": prop_count,
        "mol_insert_time_s": total_mol_elapsed,
        "prop_insert_time_s": prop_elapsed,
        "total_insert_time_s": total_mol_elapsed + prop_elapsed,
        "avg_batch_time_ms": statistics.mean(batch_times),
        "mol_throughput_per_s": mol_count / total_mol_elapsed,
    }


# ============================================================================
# Step 4: Performance tests via HTTP API
# ============================================================================

def run_perf_test(
    name: str,
    url: str,
    headers: dict,
    n_runs: int = 3,
) -> dict:
    """Run a GET endpoint N times, return stats."""
    times = []
    statuses = []
    first_body = None

    for i in range(n_runs):
        t0 = time.perf_counter()
        req = urllib.request.Request(url, headers=headers, method="GET")
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                body = json.loads(resp.read())
                status = resp.status
        except urllib.error.HTTPError as e:
            body = {}
            status = e.code
        elapsed = (time.perf_counter() - t0) * 1000
        times.append(elapsed)
        statuses.append(status)
        if i == 0:
            first_body = body

    avg = statistics.mean(times)
    mn = min(times)
    mx = max(times)
    # p95 — with only 3 runs, p95 = max; with more runs use sorted
    p95 = sorted(times)[int(len(times) * 0.95)] if len(times) >= 20 else mx

    passed = avg < PASS_THRESHOLD_MS
    status_ok = all(s in (200, 201) for s in statuses)

    result = {
        "name": name,
        "url": url,
        "runs": n_runs,
        "times_ms": [round(t, 1) for t in times],
        "avg_ms": round(avg, 1),
        "min_ms": round(mn, 1),
        "max_ms": round(mx, 1),
        "p95_ms": round(p95, 1),
        "all_200": status_ok,
        "passed": passed and status_ok,
        "first_body_preview": first_body,
    }
    status_str = "PASS" if (passed and status_ok) else "FAIL"
    log(
        f"  [{status_str}] {name}: avg={avg:.0f}ms min={mn:.0f}ms max={mx:.0f}ms p95={p95:.0f}ms"
    )
    return result


def run_concurrent_test(url: str, headers: dict, n_concurrent: int = 10) -> dict:
    """Fire N concurrent GET requests and measure throughput."""
    sep(f"Concurrent Test ({n_concurrent} simultaneous requests)")
    log(f"URL: {url}")

    times = []
    errors = 0

    def make_request(_):
        t0 = time.perf_counter()
        req = urllib.request.Request(url, headers=headers, method="GET")
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                resp.read()  # consume body
                return (time.perf_counter() - t0) * 1000, resp.status
        except Exception as e:
            return (time.perf_counter() - t0) * 1000, 0

    wall_start = time.perf_counter()
    with ThreadPoolExecutor(max_workers=n_concurrent) as pool:
        futures = [pool.submit(make_request, i) for i in range(n_concurrent)]
        results = [f.result() for f in as_completed(futures)]
    wall_elapsed = time.perf_counter() - wall_start

    for elapsed_ms, status in results:
        times.append(elapsed_ms)
        if status not in (200, 201):
            errors += 1

    avg = statistics.mean(times)
    mn = min(times)
    mx = max(times)
    p95 = sorted(times)[int(len(times) * 0.95)]
    throughput = n_concurrent / wall_elapsed

    passed = avg < PASS_THRESHOLD_MS and errors == 0
    status_str = "PASS" if passed else "FAIL"

    log(f"  [{status_str}] avg={avg:.0f}ms min={mn:.0f}ms max={mx:.0f}ms p95={p95:.0f}ms")
    log(f"  Wall time: {wall_elapsed*1000:.0f}ms | Throughput: {throughput:.1f} req/s | Errors: {errors}")

    return {
        "name": f"Concurrent ({n_concurrent} requests)",
        "concurrent": n_concurrent,
        "times_ms": [round(t, 1) for t in times],
        "avg_ms": round(avg, 1),
        "min_ms": round(mn, 1),
        "max_ms": round(mx, 1),
        "p95_ms": round(p95, 1),
        "wall_elapsed_ms": round(wall_elapsed * 1000, 1),
        "throughput_rps": round(throughput, 2),
        "errors": errors,
        "passed": passed,
    }


def run_cursor_pagination_test(
    phase_id: str,
    token: str,
    n_pages: int = 5,
) -> dict:
    """Follow cursor pagination for N pages, measure per-page timing."""
    sep("Cursor Pagination Test (5 pages deep)")
    headers = {"Authorization": f"Bearer {token}"}
    base_url = f"{BACKEND_URL}/api/v9/phases/{phase_id}/molecules"

    times = []
    cursor = None

    for page in range(1, n_pages + 1):
        if cursor:
            url = f"{base_url}?limit=50&sort_by=created_at&sort_dir=desc&cursor={urllib.parse.quote(cursor)}"
        else:
            url = f"{base_url}?limit=50&sort_by=created_at&sort_dir=desc"

        t0 = time.perf_counter()
        req = urllib.request.Request(url, headers=headers, method="GET")
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                body = json.loads(resp.read())
                status = resp.status
        except urllib.error.HTTPError as e:
            body = {}
            status = e.code
        elapsed = (time.perf_counter() - t0) * 1000
        times.append(elapsed)

        cursor = body.get("next_cursor")
        mol_count = len(body.get("molecules", []))
        log(f"  Page {page}: {elapsed:.0f}ms | {mol_count} molecules | has_more={body.get('has_more')} | status={status}")

        if not cursor:
            log(f"  No more cursor after page {page}")
            break

    avg = statistics.mean(times)
    passed = avg < PASS_THRESHOLD_MS

    return {
        "name": "Cursor pagination (5 pages)",
        "pages": len(times),
        "times_ms": [round(t, 1) for t in times],
        "avg_ms": round(avg, 1),
        "max_ms": round(max(times), 1),
        "passed": passed,
    }


def get_system_info() -> dict:
    """Collect basic system resource info."""
    info = {}
    try:
        import subprocess
        # Memory
        result = subprocess.run(["free", "-m"], capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            lines = result.stdout.strip().split("\n")
            if len(lines) >= 2:
                parts = lines[1].split()
                info["mem_total_mb"] = int(parts[1]) if len(parts) > 1 else "N/A"
                info["mem_used_mb"] = int(parts[2]) if len(parts) > 2 else "N/A"
                info["mem_free_mb"] = int(parts[3]) if len(parts) > 3 else "N/A"
        # CPU
        result = subprocess.run(["nproc"], capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            info["cpu_cores"] = int(result.stdout.strip())
        # Load
        result = subprocess.run(["uptime"], capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            info["uptime"] = result.stdout.strip()
    except Exception:
        pass
    return info


# ============================================================================
# Step 5: Cleanup
# ============================================================================

def cleanup(project_id: str, token: str):
    sep("STEP 5: Cleanup")
    log(f"Deleting test project {project_id} (and all cascaded data)...")

    dsn = PG_DSN
    try:
        conn = psycopg2.connect(dsn)
        conn.autocommit = False
        cur = conn.cursor()

        # Count before delete
        cur.execute(
            "SELECT COUNT(*) FROM molecules m "
            "JOIN phases ph ON m.phase_id = ph.id "
            "JOIN campaigns c ON ph.campaign_id = c.id "
            "WHERE c.project_id = %s",
            (project_id,)
        )
        mol_count = cur.fetchone()[0]
        log(f"Deleting {mol_count:,} molecules (cascades from project delete)...")

        cur.execute("DELETE FROM projects WHERE id = %s", (project_id,))
        conn.commit()
        log(f"Project {project_id} deleted successfully.")
        cur.close()
        conn.close()
    except Exception as e:
        log(f"Cleanup error (non-critical): {e}", "WARN")


# ============================================================================
# Step 6: Report
# ============================================================================

def print_report(
    insert_stats: dict,
    perf_results: list[dict],
    concurrent_result: dict,
    cursor_result: dict,
    system_info: dict,
):
    sep("FINAL PERFORMANCE REPORT")

    print(f"""
SYSTEM INFO
-----------
CPU Cores     : {system_info.get('cpu_cores', 'N/A')}
Memory Total  : {system_info.get('mem_total_mb', 'N/A')} MB
Memory Used   : {system_info.get('mem_used_mb', 'N/A')} MB
Uptime        : {system_info.get('uptime', 'N/A')}

INSERT PERFORMANCE (100K Molecules)
------------------------------------
Molecules inserted      : {insert_stats['molecules_inserted']:,}
Properties inserted     : {insert_stats['properties_inserted']:,}
Molecule insert time    : {insert_stats['mol_insert_time_s']:.1f}s
Property insert time    : {insert_stats['prop_insert_time_s']:.1f}s
Total insert time       : {insert_stats['total_insert_time_s']:.1f}s
Avg batch time          : {insert_stats['avg_batch_time_ms']:.0f}ms / {BATCH_SIZE} molecules
Molecule throughput     : {insert_stats['mol_throughput_per_s']:.0f} molecules/s
""")

    print("API ENDPOINT PERFORMANCE")
    print("-" * 70)
    print(f"{'Test':<40} {'Avg':>8} {'Min':>8} {'Max':>8} {'P95':>8} {'Result':>8}")
    print("-" * 70)

    all_results = perf_results + [cursor_result, concurrent_result]
    pass_count = 0
    fail_count = 0

    for r in all_results:
        status = "PASS" if r.get("passed") else "FAIL"
        if r.get("passed"):
            pass_count += 1
        else:
            fail_count += 1
        avg = r.get("avg_ms", 0)
        mn = r.get("min_ms", 0)
        mx = r.get("max_ms", 0)
        p95 = r.get("p95_ms", mx)
        print(f"  {r['name']:<38} {avg:>7.0f}ms {mn:>7.0f}ms {mx:>7.0f}ms {p95:>7.0f}ms  [{status}]")

    print("-" * 70)
    print(f"\nCONCURRENCY TEST")
    print(f"  Throughput       : {concurrent_result.get('throughput_rps', 0):.1f} req/s")
    print(f"  Wall time        : {concurrent_result.get('wall_elapsed_ms', 0):.0f}ms")
    print(f"  Errors           : {concurrent_result.get('errors', 0)}")

    print(f"""
SUMMARY
-------
Tests passed  : {pass_count}/{pass_count + fail_count}
Tests failed  : {fail_count}/{pass_count + fail_count}
Pass threshold: {PASS_THRESHOLD_MS}ms average response time
""")

    if fail_count == 0:
        print("RESULT: ALL TESTS PASSED — Platform handles 100K molecules at production scale.")
    else:
        print("RESULT: SOME TESTS FAILED — See recommendations below.")

    # Recommendations
    print("\nRECOMMENDATIONS")
    print("-" * 40)
    for r in all_results:
        if not r.get("passed"):
            name = r["name"]
            avg = r.get("avg_ms", 0)
            if "docking_score" in name or "composite_score" in name:
                print(f"  - [{name}] avg={avg:.0f}ms: JSONB lateral join sort is slow at 100K scale.")
                print(f"    → Consider materialized column or GIN index on property values.")
            elif "concurrent" in name.lower():
                print(f"  - [{name}] avg={avg:.0f}ms: High concurrency saturation detected.")
                print(f"    → Increase DB connection pool or add read replica.")
            elif "stats" in name.lower():
                print(f"  - [{name}] avg={avg:.0f}ms: COUNT aggregation slow.")
                print(f"    → Consider a cached counter table updated via triggers.")
            else:
                print(f"  - [{name}] avg={avg:.0f}ms: Exceeds {PASS_THRESHOLD_MS}ms threshold.")
                print(f"    → Check index usage with EXPLAIN ANALYZE.")

    if fail_count == 0:
        print("  No bottlenecks detected. Indexes are working effectively.")

    sep()


# ============================================================================
# Main
# ============================================================================

def main():
    sep("BindX V2 — 100K Molecule Stress Test")
    log(f"Target backend: {BACKEND_URL}")
    log(f"Molecules to insert: {TOTAL_MOLECULES:,}")
    log(f"Pass threshold: {PASS_THRESHOLD_MS}ms")

    system_info = get_system_info()
    log(f"System: {system_info.get('cpu_cores', '?')} cores, {system_info.get('mem_total_mb', '?')} MB RAM")

    # Step 1: Auth
    token = authenticate()
    auth_headers = {"Authorization": f"Bearer {token}"}

    # Step 2: Create hierarchy
    project_id, campaign_id, phase_id = create_test_hierarchy(token)
    log(f"\nTest hierarchy:")
    log(f"  project_id  = {project_id}")
    log(f"  campaign_id = {campaign_id}")
    log(f"  phase_id    = {phase_id}")

    # Step 3: Bulk insert
    insert_stats = bulk_insert_molecules(phase_id)

    # Step 4: Performance tests
    sep("STEP 4: API Performance Tests")
    log(f"Pass threshold: avg < {PASS_THRESHOLD_MS}ms")

    perf_results = []

    # Test 1: List page 1 (created_at desc)
    log("\n[Test 1] List molecules page 1 (created_at desc, limit=50)...")
    url1 = f"{BACKEND_URL}/api/v9/phases/{phase_id}/molecules?limit=50&sort_by=created_at&sort_dir=desc"
    perf_results.append(run_perf_test("List page 1 (created_at)", url1, auth_headers, n_runs=5))

    # Test 2: Cursor pagination (5 pages)
    cursor_result = run_cursor_pagination_test(phase_id, token, n_pages=5)

    # Test 3: Sort by docking_score (JSONB lateral join)
    log("\n[Test 3] Sort by docking_score (JSONB sort, 100K rows)...")
    url3 = f"{BACKEND_URL}/api/v9/phases/{phase_id}/molecules?limit=50&sort_by=docking_score&sort_dir=asc"
    perf_results.append(run_perf_test("Sort by docking_score (JSONB)", url3, auth_headers, n_runs=3))

    # Test 4: Sort by composite_score
    log("\n[Test 4] Sort by composite_score (JSONB sort, 100K rows)...")
    url4 = f"{BACKEND_URL}/api/v9/phases/{phase_id}/molecules?limit=50&sort_by=composite_score&sort_dir=desc"
    perf_results.append(run_perf_test("Sort by composite_score (JSONB)", url4, auth_headers, n_runs=3))

    # Test 5: Sort by QED
    log("\n[Test 5] Sort by QED (JSONB sort)...")
    url5 = f"{BACKEND_URL}/api/v9/phases/{phase_id}/molecules?limit=50&sort_by=QED&sort_dir=desc"
    perf_results.append(run_perf_test("Sort by QED (JSONB)", url5, auth_headers, n_runs=3))

    # Test 6: Stats endpoint
    log("\n[Test 6] Stats endpoint (COUNT on 100K rows)...")
    url6 = f"{BACKEND_URL}/api/v9/phases/{phase_id}/molecules/stats"
    perf_results.append(run_perf_test("Stats (COUNT, bookmarked, ai_gen)", url6, auth_headers, n_runs=5))

    # Test 7: Phase detail (includes molecule_count)
    log("\n[Test 7] Phase detail (includes molecule_count)...")
    url7 = f"{BACKEND_URL}/api/v9/phases/{phase_id}"
    perf_results.append(run_perf_test("Phase detail (molecule_count)", url7, auth_headers, n_runs=5))

    # Test 8: Search query
    log("\n[Test 8] Search filter (ILIKE on 100K molecules)...")
    url8 = f"{BACKEND_URL}/api/v9/phases/{phase_id}/molecules?limit=50&search=MOL_05"
    perf_results.append(run_perf_test("Search (ILIKE filter)", url8, auth_headers, n_runs=3))

    # Test 9: Bookmarked only filter
    log("\n[Test 9] Bookmarked-only filter...")
    url9 = f"{BACKEND_URL}/api/v9/phases/{phase_id}/molecules?limit=50&bookmarked_only=true"
    perf_results.append(run_perf_test("Bookmarked-only filter", url9, auth_headers, n_runs=3))

    # Test 10: Concurrent requests
    url_concurrent = f"{BACKEND_URL}/api/v9/phases/{phase_id}/molecules?limit=50&sort_by=created_at&sort_dir=desc"
    concurrent_result = run_concurrent_test(url_concurrent, auth_headers, n_concurrent=10)

    # Extra: 20 concurrent
    log("\n[Test 11] 20 simultaneous requests...")
    concurrent_result_20 = run_concurrent_test(url_concurrent, auth_headers, n_concurrent=20)
    # We'll use 10-concurrent as the main metric but log 20 too

    # Step 5: Cleanup
    cleanup(project_id, token)

    # Step 6: Report
    print_report(insert_stats, perf_results, concurrent_result, cursor_result, system_info)

    # Also print raw JSON for reference
    print("\nRAW RESULTS (JSON)")
    print("-" * 40)
    raw = {
        "insert_stats": insert_stats,
        "perf_results": perf_results,
        "cursor_pagination": cursor_result,
        "concurrent_10": concurrent_result,
        "concurrent_20": {
            "avg_ms": concurrent_result_20.get("avg_ms"),
            "max_ms": concurrent_result_20.get("max_ms"),
            "throughput_rps": concurrent_result_20.get("throughput_rps"),
            "errors": concurrent_result_20.get("errors"),
        },
    }
    print(json.dumps(raw, indent=2, default=str))


if __name__ == "__main__":
    main()
