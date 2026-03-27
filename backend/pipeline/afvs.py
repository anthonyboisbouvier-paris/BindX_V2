"""
BindX V9 — AFVS Pipeline Module (BDX-41)

AdaptiveFlow Virtual Screening controller.
Orchestrates AFVS via SSH on an EC2 controller instance,
submits AWS Batch jobs, polls completion, collects top results,
and imports hits into BindX Phase A dashboard.

Environment variables:
  AFVS_EC2_HOST        — EC2 controller IP/hostname
  AFVS_EC2_USER        — SSH user (default: ubuntu)
  AFVS_SSH_KEY_PATH    — Path to SSH private key
  AFVS_S3_JOB_BUCKET   — S3 bucket for AFVS job data + outputs
  AFVS_S3_DATA_BUCKET  — S3 bucket for Enamine REAL Space library
"""

from __future__ import annotations

import csv
import io
import json
import logging
import os
import time
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Environment / Config
# ---------------------------------------------------------------------------

AFVS_EC2_HOST = os.environ.get("AFVS_EC2_HOST", "")
AFVS_EC2_USER = os.environ.get("AFVS_EC2_USER", "ubuntu")
AFVS_SSH_KEY_PATH = os.environ.get("AFVS_SSH_KEY_PATH", "/app/secrets/afvs_ec2_key.pem")
AFVS_S3_JOB_BUCKET = os.environ.get("AFVS_S3_JOB_BUCKET", "bindx-afvs-outputs")
AFVS_S3_DATA_BUCKET = os.environ.get("AFVS_S3_DATA_BUCKET", "")

# Fixed paths on EC2 controller
AFVS_TOOLS_DIR = f"/home/{os.environ.get('AFVS_EC2_USER', 'ubuntu')}/AFVS/tools"
AFVS_INPUT_DIR = f"/home/{os.environ.get('AFVS_EC2_USER', 'ubuntu')}/AFVS/input-files"
AFVS_WORKFLOW_DIR = f"/home/{os.environ.get('AFVS_EC2_USER', 'ubuntu')}/AFVS/workflow"
AFVS_VENV = f"/home/{os.environ.get('AFVS_EC2_USER', 'ubuntu')}/afvs_env/bin/activate"

# ---------------------------------------------------------------------------
# Strategy presets — paper-validated (Gorgulla 2023, VirtualFlow 2.0)
# reps=1, exhaust_pre=1, exhaust_pri=8 for all. Differentiation = coverage.
# ---------------------------------------------------------------------------

STRATEGY_PRESETS = {
    "lean": {
        "exhaustiveness_prescreen": 1,
        "exhaustiveness_primary": 8,
        "reps_per_tranche": 1,
        "tranche_pct": 0.05,
        "top_n_results": 5_000,
        "budget_cap_usd": 250,
    },
    "balanced": {
        "exhaustiveness_prescreen": 1,
        "exhaustiveness_primary": 8,
        "reps_per_tranche": 1,
        "tranche_pct": 0.2,
        "top_n_results": 10_000,
        "budget_cap_usd": 700,
    },
    "aggressive": {
        "exhaustiveness_prescreen": 1,
        "exhaustiveness_primary": 8,
        "reps_per_tranche": 1,
        "tranche_pct": 1.0,
        "top_n_results": 50_000,
        "budget_cap_usd": 3_500,
    },
}

# Frontend docking scenario → AFVS docking method mapping
DOCKING_METHODS = {
    "smina_vinardo": "smina_rigid",
    "qvina2": "qvina02",
}


def estimate_cost(tranche_pct: float, reps: int) -> Optional[float]:
    """Indicative cost estimate. Returns None until calibrated."""
    return None


# ---------------------------------------------------------------------------
# AFVSController — SSH wrapper for EC2 controller
# ---------------------------------------------------------------------------

class AFVSController:
    """SSH-based controller for AdaptiveFlow on EC2.

    All AFVS operations are performed via SSH to the EC2 controller
    instance, which has AdaptiveFlow installed and IAM role access to
    AWS Batch and S3.
    """

    def __init__(self):
        self._client = None

    # --- SSH connection management ---

    def connect(self):
        """Establish SSH connection to EC2 controller with keepalive."""
        import paramiko

        if self._client is not None:
            try:
                self._client.exec_command("echo ok", timeout=5)
                return
            except Exception:
                self._close()

        if not AFVS_EC2_HOST:
            raise RuntimeError("AFVS_EC2_HOST not configured — set it in .env")

        key_path = Path(AFVS_SSH_KEY_PATH)
        if not key_path.exists():
            raise RuntimeError(f"SSH key not found: {AFVS_SSH_KEY_PATH}")

        client = paramiko.SSHClient()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        client.connect(
            hostname=AFVS_EC2_HOST,
            username=AFVS_EC2_USER,
            key_filename=str(key_path),
            timeout=30,
        )
        transport = client.get_transport()
        if transport:
            transport.set_keepalive(60)
        self._client = client
        logger.info("SSH connected to %s@%s", AFVS_EC2_USER, AFVS_EC2_HOST)

    def _close(self):
        if self._client:
            try:
                self._client.close()
            except Exception:
                pass
            self._client = None

    def run_command(self, cmd: str, timeout: int = 300) -> tuple[str, str, int]:
        """Execute a command via SSH. Auto-reconnects up to 3 times."""
        for attempt in range(3):
            try:
                self.connect()
                stdin, stdout, stderr = self._client.exec_command(cmd, timeout=timeout)
                exit_code = stdout.channel.recv_exit_status()
                return stdout.read().decode(), stderr.read().decode(), exit_code
            except Exception as e:
                logger.warning("SSH command attempt %d failed: %s", attempt + 1, e)
                self._close()
                if attempt == 2:
                    raise RuntimeError(f"SSH command failed after 3 attempts: {e}")
                time.sleep(5)

    def upload_file(self, local_path: str, remote_path: str):
        """Upload a file via SFTP."""
        self.connect()
        sftp = self._client.open_sftp()
        try:
            sftp.put(local_path, remote_path)
            logger.info("Uploaded %s → %s", local_path, remote_path)
        finally:
            sftp.close()

    # --- AFVS job lifecycle ---

    def cleanup_previous(self):
        """Remove previous workflow directory to avoid conflicts."""
        self.run_command(f"rm -rf {AFVS_WORKFLOW_DIR}")
        logger.info("Cleaned up previous AFVS workflow")

    def setup_job(
        self,
        run_id: str,
        receptor_pdb_path: str,
        docking_box: dict,
        config: dict,
    ) -> tuple[str, str, int]:
        """Set up a complete AFVS job on EC2.

        1. Upload receptor to scenario folder
        2. Write docking config.txt
        3. Write all.ctrl (AdaptiveFlow master config)
        4. Run afvs_prepare_folders.py
        5. Run afvs_prepare_workunits.py

        Returns (job_name, scenario_name, num_workunits).
        """
        short_id = run_id[:8].replace("-", "")
        job_name = f"bx{short_id}"
        scenario_name = f"bx_{short_id}"
        method = DOCKING_METHODS.get(
            config.get("docking_scenario", "smina_vinardo"), "smina_rigid"
        )
        exhaust = config.get("exhaustiveness_prescreen", 1)
        reps = config.get("reps_per_tranche", 1)

        # Validate data bucket
        if not AFVS_S3_DATA_BUCKET:
            raise RuntimeError(
                "AFVS_S3_DATA_BUCKET not configured. "
                "Request REAL Space access at virtual-flow.org and set the bucket name."
            )

        # 1. Create scenario folder + upload receptor
        scenario_dir = f"{AFVS_INPUT_DIR}/{scenario_name}"
        self.run_command(f"mkdir -p {scenario_dir}")
        self.upload_file(receptor_pdb_path, f"{scenario_dir}/receptor.pdb")
        logger.info("Uploaded receptor to %s", scenario_dir)

        # 2. Write docking config.txt (Vina-family format)
        box = docking_box
        config_txt = (
            f"receptor = receptor.pdb\n"
            f"center_x = {box['center_x']}\n"
            f"center_y = {box['center_y']}\n"
            f"center_z = {box['center_z']}\n"
            f"size_x = {box.get('size_x', 22)}\n"
            f"size_y = {box.get('size_y', 22)}\n"
            f"size_z = {box.get('size_z', 22)}\n"
            f"exhaustiveness = {exhaust}\n"
            f"num_modes = 1\n"
        )
        self._write_remote_file(f"{scenario_dir}/config.txt", config_txt)

        # 3. Write all.ctrl (AdaptiveFlow master config)
        all_ctrl = self._build_all_ctrl(job_name, scenario_name, method, reps)
        self._write_remote_file(f"{AFVS_TOOLS_DIR}/templates/all.ctrl", all_ctrl)
        logger.info("Wrote all.ctrl for job %s", job_name)

        # 4. Prepare folders
        stdout, stderr, code = self.run_command(
            f"cd {AFVS_TOOLS_DIR} && source {AFVS_VENV} && "
            f"python3 ./afvs_prepare_folders.py 2>&1",
            timeout=300,
        )
        if code != 0:
            raise RuntimeError(f"afvs_prepare_folders failed: {(stderr or stdout)[:500]}")

        # 5. Prepare workunits (packages config + receptor into tar.gz, uploads to S3)
        stdout, stderr, code = self.run_command(
            f"cd {AFVS_TOOLS_DIR} && source {AFVS_VENV} && "
            f"python3 ./afvs_prepare_workunits.py 2>&1",
            timeout=600,
        )
        if code != 0:
            raise RuntimeError(f"afvs_prepare_workunits failed: {(stderr or stdout)[:500]}")

        # Get workunit count from workflow/status.json
        num_workunits = self._get_workunit_count()
        logger.info("Prepared %d workunits for job %s", num_workunits, job_name)

        return job_name, scenario_name, num_workunits

    def submit_jobs(self, num_workunits: int):
        """Submit all workunits to AWS Batch."""
        if num_workunits <= 0:
            raise RuntimeError("No workunits to submit")

        stdout, stderr, code = self.run_command(
            f"cd {AFVS_TOOLS_DIR} && source {AFVS_VENV} && "
            f"python3 ./afvs_submit_jobs.py 1 {num_workunits} 2>&1",
            timeout=600,
        )
        if code != 0:
            raise RuntimeError(f"Job submission failed: {(stderr or stdout)[:500]}")
        logger.info("Submitted %d workunits to AWS Batch", num_workunits)

    def get_job_status(self) -> dict:
        """Query AWS Batch for job completion status.

        Returns dict with keys: total, succeeded, failed, running, pending, done.
        """
        status_script = r'''
import json, boto3, sys
try:
    with open("../workflow/status.json") as f:
        status = json.load(f)
except FileNotFoundError:
    print(json.dumps({"error": "no workflow", "done": False}))
    sys.exit(0)

batch = boto3.client("batch", region_name="us-east-1")
total = 0; succeeded = 0; failed = 0; running = 0; pending = 0
for wu_id, wu in status.get("workunits", {}).items():
    job_id = wu.get("batch_job_id")
    if not job_id:
        pending += 1
        total += 1
        continue
    try:
        resp = batch.describe_jobs(jobs=[job_id])
        for job in resp.get("jobs", []):
            s = job["status"]
            if s == "SUCCEEDED":
                succeeded += 1
            elif s == "FAILED":
                failed += 1
            elif s in ("RUNNING", "STARTING"):
                running += 1
            else:
                pending += 1
            total += 1
    except Exception:
        total += 1
        pending += 1

done = total > 0 and (succeeded + failed == total)
pct = round(succeeded / total * 100) if total > 0 else 0
print(json.dumps({
    "total": total, "succeeded": succeeded, "failed": failed,
    "running": running, "pending": pending, "done": done, "pct": pct
}))
'''
        self._write_remote_file("/tmp/bindx_status.py", status_script)
        stdout, stderr, code = self.run_command(
            f"cd {AFVS_TOOLS_DIR} && source {AFVS_VENV} && "
            f"python3 /tmp/bindx_status.py",
            timeout=60,
        )
        try:
            last_line = stdout.strip().split("\n")[-1]
            return json.loads(last_line)
        except (json.JSONDecodeError, IndexError):
            logger.warning("Failed to parse status output: %s", stdout[:200])
            return {"error": (stderr or stdout)[:200], "done": False}

    def collect_top_results(
        self, job_name: str, scenario_name: str, top_n: int
    ) -> str:
        """Collect top N results from CSV.gz files on S3.

        Downloads all result CSV.gz files, sorts by score, takes top N,
        and writes a single CSV to S3.

        Returns S3 path of the results CSV.
        """
        output_key = f"AF/AFVS/jobs/{job_name}/bindx_top_results.csv"

        collect_script = f'''
import boto3, csv, gzip, io, json, sys

bucket = "{AFVS_S3_JOB_BUCKET}"
top_n = {top_n}
output_key = "{output_key}"

with open("../workflow/config.json") as f:
    config = json.load(f)

prefix = config.get("object_store_job_prefix", "AF/AFVS/jobs")
job_name = config.get("job_name", "{job_name}")
scenario = "{scenario_name}"
search_prefix = f"{{prefix}}/{{job_name}}/{{scenario}}/csv/"

s3 = boto3.client("s3")
all_scores = []
files_read = 0

paginator = s3.get_paginator("list_objects_v2")
for page in paginator.paginate(Bucket=bucket, Prefix=search_prefix):
    for obj in page.get("Contents", []):
        key = obj["Key"]
        if not key.endswith(".csv.gz"):
            continue
        try:
            response = s3.get_object(Bucket=bucket, Key=key)
            with gzip.open(io.BytesIO(response["Body"].read()), "rt") as gz:
                reader = csv.DictReader(gz)
                for row in reader:
                    smi = row.get("smi", "")
                    score_raw = row.get("score_min") or row.get("score_1") or row.get("score", "")
                    if smi and score_raw:
                        try:
                            all_scores.append((float(score_raw), smi))
                        except ValueError:
                            pass
            files_read += 1
        except Exception as e:
            print(f"Warning: {{key}}: {{e}}", file=sys.stderr)

all_scores.sort()
top = all_scores[:top_n]

output = io.StringIO()
writer = csv.writer(output)
writer.writerow(["smiles", "docking_score"])
for score, smi in top:
    writer.writerow([smi, score])

s3.put_object(Bucket=bucket, Key=output_key, Body=output.getvalue().encode())
print(f"{{len(top)}} results from {{files_read}} files -> s3://{{bucket}}/{{output_key}}")
'''
        self._write_remote_file("/tmp/bindx_collect.py", collect_script)
        stdout, stderr, code = self.run_command(
            f"cd {AFVS_TOOLS_DIR} && source {AFVS_VENV} && "
            f"python3 /tmp/bindx_collect.py 2>&1",
            timeout=1800,
        )
        if code != 0:
            raise RuntimeError(f"Result collection failed: {(stderr or stdout)[:500]}")

        logger.info("Collected top results: %s", stdout.strip())
        return f"s3://{AFVS_S3_JOB_BUCKET}/{output_key}"

    def cancel_all_jobs(self):
        """Cancel all AWS Batch jobs for the current workflow."""
        self.run_command(
            f"cd {AFVS_TOOLS_DIR} && source {AFVS_VENV} && "
            f"python3 ./afvs_killall_awsbatch.py 2>/dev/null || true",
            timeout=120,
        )
        logger.info("Cancelled all AFVS Batch jobs")

    def cleanup_job(self, scenario_name: str):
        """Remove job-specific files from EC2."""
        self.run_command(f"rm -rf {AFVS_INPUT_DIR}/{scenario_name}")
        self.run_command(f"rm -rf {AFVS_WORKFLOW_DIR}")
        self.run_command("rm -f /tmp/bindx_status.py /tmp/bindx_collect.py")

    # --- Internal helpers ---

    def _write_remote_file(self, remote_path: str, content: str):
        """Write content to a file on EC2 via heredoc."""
        # Use a unique delimiter to avoid conflicts with content
        self.run_command(
            f"cat > {remote_path} << 'BINDX_HEREDOC_EOF'\n{content}\nBINDX_HEREDOC_EOF"
        )

    def _get_workunit_count(self) -> int:
        """Read workunit count from workflow/status.json."""
        stdout, _, code = self.run_command(
            f"python3 -c \""
            f"import json; "
            f"d=json.load(open('{AFVS_WORKFLOW_DIR}/status.json')); "
            f"print(len(d.get('workunits',{{}})))"
            f"\""
        )
        try:
            return int(stdout.strip())
        except (ValueError, AttributeError):
            return 0

    def _build_all_ctrl(
        self, job_name: str, scenario_name: str, method: str, reps: int
    ) -> str:
        """Generate AdaptiveFlow all.ctrl content (real format)."""
        return f"""# BindX AFVS — auto-generated for job {job_name}
# Do not edit manually — overwritten on each BindX AFVS run

# --- Job ---
job_name={job_name}

# --- Batch system ---
batchsystem=awsbatch
threads_per_docking=1
threads_to_use=16
program_timeout=90

# --- AWS Batch ---
aws_batch_prefix=af
aws_batch_number_of_queues=1
aws_batch_jobdef=af-jobdef-afvs
aws_batch_array_job_size=200
aws_ecr_repository_name=af-afvs-ecr
aws_region=us-east-1
aws_batch_subjob_vcpus=8
aws_batch_subjob_memory=15000
aws_batch_subjob_timeout=20800

# --- Storage ---
data_storage_mode=s3
job_storage_mode=s3
data_collection_addressing_mode=hash
data_collection_identifier=Enamine_REAL_Space_2022q12-sparse
job_addressing_mode=metatranche
object_store_job_bucket={AFVS_S3_JOB_BUCKET}
object_store_job_prefix=AF/AFVS/jobs
object_store_data_bucket={AFVS_S3_DATA_BUCKET}
object_store_data_collection_prefix=Enamine_REAL_Space_2022q12
collection_folder=/home/{AFVS_EC2_USER}/collections

# --- Output ---
summary_formats=csv.gz
print_attrs_in_summary=smi,heavy_atom_count

# --- Workflow ---
collection_list_type=standard
collection_list_file=templates/todo.all
dockings_per_subjob=1000
ligand_library_format=pdbqt
dynamic_tranche_filtering=0
tempdir_default=/dev/shm

# --- Docking scenario ---
docking_scenario_names={scenario_name}
docking_scenario_methods={method}
docking_scenario_replicas=1
docking_scenario_batchsizes=1
docking_scenario_basefolder=../input-files

# --- Prescreen ---
prescreen_mode=1
prescreen_ligands_per_tranche={reps}
run_atom_check=1
"""


# ---------------------------------------------------------------------------
# Result parsing
# ---------------------------------------------------------------------------

def parse_afvs_results_csv(csv_content: str) -> list[dict]:
    """Parse AFVS top_results CSV into molecule dicts.

    Expected columns: smiles, docking_score
    Returns list of {smiles, canonical_smiles, name, docking_score}.
    """
    try:
        from rdkit import Chem
        has_rdkit = True
    except ImportError:
        has_rdkit = False

    results = []
    reader = csv.DictReader(io.StringIO(csv_content))

    for row in reader:
        smiles = (row.get("smiles") or row.get("SMILES") or row.get("smi") or "").strip()
        if not smiles:
            continue

        if has_rdkit:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            canonical = Chem.MolToSmiles(mol)
        else:
            canonical = smiles

        name = row.get("name") or row.get("zinc_id") or row.get("compound_id") or None
        score_raw = (
            row.get("docking_score") or row.get("score_min")
            or row.get("score") or row.get("affinity")
        )

        try:
            score = float(score_raw) if score_raw else None
        except (ValueError, TypeError):
            score = None

        results.append({
            "smiles": smiles,
            "canonical_smiles": canonical,
            "name": name,
            "docking_score": score,
        })

    return results


def download_s3_csv(s3_path: str) -> str:
    """Download a CSV file from S3 and return its content as string.

    Uses boto3. Requires AWS credentials (env vars or IAM role).
    Falls back to SSH download via EC2 if boto3 fails locally.
    """
    import boto3

    parts = s3_path.replace("s3://", "").split("/", 1)
    bucket = parts[0]
    key = parts[1] if len(parts) > 1 else ""

    s3 = boto3.client("s3")
    response = s3.get_object(Bucket=bucket, Key=key)
    return response["Body"].read().decode("utf-8")


def download_s3_csv_via_ec2(controller: AFVSController, s3_path: str) -> str:
    """Download CSV from S3 via the EC2 controller (uses its IAM role)."""
    stdout, stderr, code = controller.run_command(
        f"aws s3 cp {s3_path} - 2>/dev/null",
        timeout=300,
    )
    if code != 0:
        raise RuntimeError(f"S3 download failed: {stderr[:300]}")
    return stdout
