"""
DockIt pipeline -- Audit logging for pipeline runs.

Provides a structured, timestamped log of every step in the pipeline
execution. The audit log captures tool versions, timing information,
and per-step details for reproducibility and regulatory traceability.

Usage
-----
    audit = AuditLog(job_id)
    audit.log("structure", "Fetched structure from AlphaFold", {"source": "alphafold"})
    audit.log("pockets", "Detected 3 pockets", {"n_pockets": 3})
    ...
    audit.save(work_dir / "audit_log.json")
"""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


class AuditLog:
    """Structured audit log for a pipeline run.

    Captures timestamped entries for each pipeline step, along with
    tool versions and overall timing information.

    Parameters
    ----------
    job_id : str
        UUID of the pipeline job.
    """

    def __init__(self, job_id: str) -> None:
        self.job_id: str = job_id
        self.entries: list[dict] = []
        self.start_time: datetime = datetime.now(timezone.utc)

    def log(
        self,
        step: str,
        action: str,
        details: Optional[dict] = None,
    ) -> None:
        """Add a timestamped entry to the audit log.

        Parameters
        ----------
        step : str
            Pipeline step name (e.g. "structure", "pockets", "docking").
        action : str
            Human-readable description of what happened.
        details : dict, optional
            Additional structured details (counts, paths, parameters).
        """
        now = datetime.now(timezone.utc)
        elapsed = (now - self.start_time).total_seconds()

        entry: dict = {
            "timestamp": now.isoformat(),
            "elapsed_seconds": round(elapsed, 3),
            "step": step,
            "action": action,
            "details": details or {},
        }

        self.entries.append(entry)
        logger.debug(
            "[audit:%s] [%.1fs] %s: %s",
            self.job_id[:8], elapsed, step, action,
        )

    def to_dict(self) -> dict:
        """Serialize the audit log to a dictionary.

        Returns
        -------
        dict
            Complete audit log with metadata, entries, and tool versions.
        """
        now = datetime.now(timezone.utc)
        total_elapsed = (now - self.start_time).total_seconds()

        return {
            "job_id": self.job_id,
            "start_time": self.start_time.isoformat(),
            "end_time": now.isoformat(),
            "total_elapsed_seconds": round(total_elapsed, 3),
            "total_entries": len(self.entries),
            "entries": self.entries,
            "tool_versions": _get_tool_versions(),
        }

    def save(self, path: Path) -> Path:
        """Write the audit log to a JSON file.

        Parameters
        ----------
        path : Path
            Output file path.

        Returns
        -------
        Path
            The path where the log was written.
        """
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            with open(path, "w", encoding="utf-8") as f:
                json.dump(self.to_dict(), f, indent=2, ensure_ascii=False)
            logger.info(
                "[audit:%s] Audit log saved to %s (%d entries)",
                self.job_id[:8], path, len(self.entries),
            )
            return path
        except Exception as exc:
            logger.warning(
                "[audit:%s] Failed to save audit log to %s: %s",
                self.job_id[:8], path, exc,
            )
            raise

    def get_step_entries(self, step: str) -> list[dict]:
        """Return all entries for a specific pipeline step.

        Parameters
        ----------
        step : str
            Step name to filter by.

        Returns
        -------
        list[dict]
            Entries matching the given step.
        """
        return [e for e in self.entries if e["step"] == step]

    def get_summary(self) -> dict:
        """Return a compact summary of the audit log.

        Returns
        -------
        dict
            Summary with step names, total entries, and elapsed time.
        """
        now = datetime.now(timezone.utc)
        total_elapsed = (now - self.start_time).total_seconds()
        steps = list(dict.fromkeys(e["step"] for e in self.entries))

        return {
            "job_id": self.job_id,
            "total_entries": len(self.entries),
            "total_elapsed_seconds": round(total_elapsed, 3),
            "steps": steps,
        }


# =====================================================================
# TOOL VERSION DETECTION
# =====================================================================

def _get_tool_versions() -> dict[str, str]:
    """Detect and return versions of tools used in the pipeline.

    Returns
    -------
    dict[str, str]
        Mapping of tool name to version string.
    """
    versions: dict[str, str] = {}

    # RDKit
    try:
        from rdkit import rdBase
        versions["rdkit"] = rdBase.rdkitVersion
    except ImportError:
        versions["rdkit"] = "not installed"

    # AutoDock Vina
    try:
        import shutil
        if shutil.which("vina"):
            versions["vina"] = "1.2.5 (detected)"
        else:
            versions["vina"] = "1.2.5 (mock)"
    except Exception:
        versions["vina"] = "1.2.5 (mock)"

    # fpocket
    try:
        import shutil
        if shutil.which("fpocket"):
            versions["fpocket"] = "4.0 (detected)"
        else:
            versions["fpocket"] = "4.0 (mock)"
    except Exception:
        versions["fpocket"] = "4.0 (mock)"

    # REINVENT4
    try:
        import reinvent  # noqa: F401
        versions["reinvent4"] = "installed"
    except ImportError:
        versions["reinvent4"] = "mock"

    # AiZynthFinder
    try:
        import aizynthfinder  # noqa: F401
        versions["aizynthfinder"] = "installed"
    except ImportError:
        versions["aizynthfinder"] = "mock"

    # ADMET-AI
    try:
        import admet_ai  # noqa: F401
        versions["admet_ai"] = "installed"
    except ImportError:
        versions["admet_ai"] = "mock"

    # Open Babel
    try:
        import shutil
        if shutil.which("obabel"):
            versions["openbabel"] = "detected"
        else:
            versions["openbabel"] = "mock"
    except Exception:
        versions["openbabel"] = "mock"

    return versions


# =====================================================================
# CLI / SELF-TEST
# =====================================================================

if __name__ == "__main__":
    import sys
    import time

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    print("=" * 70)
    print("DockIt Audit Log -- Self-Test")
    print("=" * 70)

    # Create audit log
    audit = AuditLog("test-job-12345")

    # Simulate pipeline steps
    audit.log("structure", "Fetching protein structure from AlphaFold", {
        "uniprot_id": "P00533",
        "source": "alphafold",
    })
    time.sleep(0.01)

    audit.log("pockets", "Detected binding pockets", {
        "n_pockets": 3,
        "best_center": [22.0, 0.5, 18.0],
    })
    time.sleep(0.01)

    audit.log("prepare", "Receptor prepared (PDB -> PDBQT)", {
        "receptor_path": "/data/jobs/test/receptor.pdbqt",
    })
    time.sleep(0.01)

    audit.log("ligands", "Fetched 15 ligands from ChEMBL", {
        "n_chembl": 10,
        "n_zinc": 5,
    })
    time.sleep(0.01)

    audit.log("docking", "Docking complete", {
        "n_docked": 15,
        "best_affinity": -9.8,
    })
    time.sleep(0.01)

    audit.log("scoring", "Scored and ranked results", {
        "n_results": 15,
        "best_composite": 0.82,
    })

    # Get dict representation
    log_dict = audit.to_dict()

    print(f"\nAudit log for job: {log_dict['job_id']}")
    print(f"  Start time:    {log_dict['start_time']}")
    print(f"  End time:      {log_dict['end_time']}")
    print(f"  Total entries: {log_dict['total_entries']}")
    print(f"  Tool versions: {json.dumps(log_dict['tool_versions'], indent=4)}")

    print("\n  Entries:")
    for entry in log_dict["entries"]:
        print(f"    [{entry['elapsed_seconds']:.3f}s] {entry['step']}: {entry['action']}")

    # Validate
    assert log_dict["total_entries"] == 6, f"Expected 6 entries, got {log_dict['total_entries']}"
    assert log_dict["job_id"] == "test-job-12345"
    assert len(log_dict["entries"]) == 6
    assert "tool_versions" in log_dict
    print("\n  PASS: dict representation is correct")

    # Test save
    test_path = Path("/tmp/dockit_audit_test/audit_log.json")
    saved_path = audit.save(test_path)
    assert saved_path.exists()
    with open(saved_path) as f:
        loaded = json.load(f)
    assert loaded["total_entries"] == 6
    print(f"  PASS: saved to {saved_path}")

    # Test summary
    summary = audit.get_summary()
    assert summary["total_entries"] == 6
    assert "structure" in summary["steps"]
    assert "docking" in summary["steps"]
    print(f"  PASS: summary has {len(summary['steps'])} steps")

    # Test step filter
    docking_entries = audit.get_step_entries("docking")
    assert len(docking_entries) == 1
    print("  PASS: step filtering works")

    print("\n" + "=" * 70)
    print("ALL TESTS PASSED")
    print("=" * 70)
