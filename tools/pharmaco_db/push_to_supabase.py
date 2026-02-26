#!/usr/bin/env python3
"""
PharmacoDB — Push to Supabase
Exports local PostgreSQL data and pushes to Supabase.
Requires: SUPABASE_DB_URL environment variable with connection string.

Usage:
  export SUPABASE_DB_URL="postgresql://postgres:PASSWORD@db.pijmvnlqljtnwhcltqwn.supabase.co:5432/postgres"
  python3 push_to_supabase.py
"""

import os
import sys
import subprocess
from datetime import datetime

LOCAL_DB = "postgresql://postgres:pharmaco_secret@localhost:5433/pharmaco"
SUPABASE_DB_URL = os.environ.get("SUPABASE_DB_URL")

def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)

def main():
    if not SUPABASE_DB_URL:
        print("ERROR: Set SUPABASE_DB_URL environment variable")
        print("  export SUPABASE_DB_URL='postgresql://postgres:PASSWORD@db.pijmvnlqljtnwhcltqwn.supabase.co:5432/postgres'")
        sys.exit(1)

    log("PharmacoDB — Push to Supabase")
    log("=" * 50)

    # Step 1: Dump the pharmaco_db schema from local
    log("Step 1: Dumping local database...")
    dump_file = "/tmp/pharmaco_db_dump.sql"
    result = subprocess.run([
        "docker", "exec", "pharmaco_db",
        "pg_dump", "-U", "postgres",
        "-d", "pharmaco",
        "-n", "pharmaco_db",
        "--no-owner", "--no-privileges",
        "--clean", "--if-exists"
    ], capture_output=True, text=True)

    if result.returncode != 0:
        log(f"ERROR: pg_dump failed: {result.stderr}")
        sys.exit(1)

    with open(dump_file, "w") as f:
        f.write(result.stdout)

    size_mb = os.path.getsize(dump_file) / 1024 / 1024
    log(f"  Dump size: {size_mb:.1f} MB")

    # Step 2: Push to Supabase
    log("Step 2: Pushing to Supabase...")
    result = subprocess.run(
        ["psql", SUPABASE_DB_URL, "-f", dump_file],
        capture_output=True, text=True,
        timeout=3600
    )

    if result.returncode != 0:
        log(f"WARNING: Some errors during push (may be OK):")
        # Only show first 20 error lines
        errors = [l for l in result.stderr.split('\n') if 'ERROR' in l]
        for e in errors[:20]:
            log(f"  {e}")
    else:
        log("  Push completed successfully!")

    log("DONE")

if __name__ == "__main__":
    main()
