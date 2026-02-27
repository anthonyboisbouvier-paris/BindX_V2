"""
BindX V9 — Async database connection (PostgreSQL via Supabase).

Provides async engine, session factory, and FastAPI dependency.
Separate from V8 database.py (SQLite/sync) — no shared state.
"""

from __future__ import annotations

import logging
import os
from typing import AsyncGenerator

from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker, create_async_engine
from sqlalchemy.orm import DeclarativeBase

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Connection URL
# ---------------------------------------------------------------------------
# Format: postgresql+asyncpg://user:pass@host:port/dbname
# For Supabase direct connection (not pooler), add ?prepared_statement_cache_size=0
# if using pgbouncer/supavisor in transaction mode.

DATABASE_V9_URL = os.environ.get(
    "DATABASE_V9_URL",
    "",
)

# ---------------------------------------------------------------------------
# Engine & Session
# ---------------------------------------------------------------------------

engine_v9 = create_async_engine(
    DATABASE_V9_URL,
    pool_pre_ping=True,
    pool_size=5,
    pool_recycle=300,
    echo=False,
) if DATABASE_V9_URL else None

AsyncSessionV9 = async_sessionmaker(
    bind=engine_v9,
    class_=AsyncSession,
    expire_on_commit=False,
) if engine_v9 else None


# ---------------------------------------------------------------------------
# Base
# ---------------------------------------------------------------------------

class BaseV9(DeclarativeBase):
    """Declarative base for V9 models. Separate from V8 Base."""
    pass


# ---------------------------------------------------------------------------
# FastAPI dependency
# ---------------------------------------------------------------------------

async def get_v9_db() -> AsyncGenerator[AsyncSession, None]:
    """Yield an async session, commit on success, rollback on error."""
    if AsyncSessionV9 is None:
        raise RuntimeError(
            "V9 database not configured. Set DATABASE_V9_URL environment variable."
        )
    async with AsyncSessionV9() as session:
        try:
            yield session
            await session.commit()
        except Exception:
            await session.rollback()
            raise
