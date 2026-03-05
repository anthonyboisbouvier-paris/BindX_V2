"""
BindX V9 — Health check router.

GET /api/v9/health     — verifies DB connection and presence of all 9 V9 tables.
GET /api/v9/health/gpu — checks RunPod GPU docking availability.
"""

from __future__ import annotations

import logging
import os

from fastapi import APIRouter
from sqlalchemy import text

from database_v9 import engine_v9

logger = logging.getLogger(__name__)

router = APIRouter()

V9_TABLES = [
    "projects",
    "campaigns",
    "phases",
    "runs",
    "molecules",
    "molecule_properties",
    "calculation_cache",
    "artifacts",
    "audit_log",
]


@router.get("/health")
async def health_check_v9():
    """Check V9 database connectivity and table presence."""
    if engine_v9 is None:
        return {
            "status": "error",
            "database": {
                "connected": False,
                "error": "DATABASE_V9_URL not configured",
            },
        }

    try:
        async with engine_v9.connect() as conn:
            # Check which V9 tables exist
            result = await conn.execute(
                text(
                    "SELECT table_name FROM information_schema.tables "
                    "WHERE table_schema = 'public' AND table_name = ANY(:tables)"
                ),
                {"tables": V9_TABLES},
            )
            found = sorted([row[0] for row in result.fetchall()])

        return {
            "status": "ok",
            "database": {
                "connected": True,
                "tables_found": found,
                "tables_expected": len(V9_TABLES),
                "tables_missing": sorted(set(V9_TABLES) - set(found)),
            },
        }
    except Exception as e:
        logger.error("V9 health check failed: %s", e)
        return {
            "status": "error",
            "database": {
                "connected": False,
                "error": str(e),
            },
        }


@router.get("/health/gpu")
async def gpu_health_check():
    """Check if GPU docking via RunPod is available."""
    api_key = os.environ.get("RUNPOD_API_KEY", "")
    endpoint_id = os.environ.get("GNINA_ENDPOINT_ID", "")
    has_key = bool(api_key and api_key != "your-runpod-api-key-here")
    has_endpoint = bool(endpoint_id)

    return {
        "gpu_available": has_key and has_endpoint,
        "endpoint_id": endpoint_id if has_endpoint else None,
        "engines": {
            "gnina_gpu": has_key and has_endpoint,
            "gnina_cpu": True,
            "vina": True,
            "diffdock": False,
        },
    }
