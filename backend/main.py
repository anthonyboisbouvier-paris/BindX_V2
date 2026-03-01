"""
BindX — FastAPI application entry point.

V9 architecture: all business logic in routers/.
Legacy V8 endpoints (preview-target, detect-pockets, scaffold, agents) kept
in routers/v8_legacy.py until migrated.
"""

from __future__ import annotations

import logging
import os
from contextlib import asynccontextmanager
from pathlib import Path
from typing import AsyncGenerator

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from database import init_db

# V9 routers
from routers.v9 import health as v9_health_router
from routers.v9 import projects as v9_projects_router
from routers.v9 import campaigns as v9_campaigns_router
from routers.v9 import phases as v9_phases_router
from routers.v9 import runs as v9_runs_router
from routers.v9 import molecules as v9_molecules_router

# Legacy V8 endpoints still used by V9 frontend
from routers.v8_legacy import router as v8_legacy_router, shutdown_pools as v8_shutdown_pools

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
)
logger = logging.getLogger("dockit.api")

DATA_DIR = Path(os.environ.get("DOCKIT_DATA_DIR", "/data")) / "jobs"


# ---------------------------------------------------------------------------
# Lifespan
# ---------------------------------------------------------------------------

@asynccontextmanager
async def lifespan(app: FastAPI) -> AsyncGenerator[None, None]:
    """Application lifespan: initialise DB on startup, clean up on shutdown."""
    logger.info("BindX API starting up")
    init_db()
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    yield
    logger.info("BindX API shutting down")
    v8_shutdown_pools()


# ---------------------------------------------------------------------------
# FastAPI app
# ---------------------------------------------------------------------------

app = FastAPI(
    title="BindX API",
    description="Drug discovery platform — V9",
    version="9.0.0",
    lifespan=lifespan,
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# V9 routers
app.include_router(v9_health_router.router, prefix="/api/v9", tags=["v9"])
app.include_router(v9_projects_router.router, prefix="/api/v9", tags=["v9-projects"])
app.include_router(v9_campaigns_router.router, prefix="/api/v9", tags=["v9-campaigns"])
app.include_router(v9_phases_router.router, prefix="/api/v9", tags=["v9-phases"])
app.include_router(v9_runs_router.router, prefix="/api/v9", tags=["v9-runs"])
app.include_router(v9_molecules_router.router, prefix="/api/v9", tags=["v9-molecules"])

# Legacy V8 endpoints: /api/preview-target, /api/detect-pockets, etc.
app.include_router(v8_legacy_router, prefix="/api", tags=["v8-legacy"])


# ---------------------------------------------------------------------------
# Health check
# ---------------------------------------------------------------------------

@app.get("/api/health")
async def health_check() -> dict:
    """Health-check endpoint."""
    return {"status": "ok", "service": "bindx-api", "version": "9.0.0"}
