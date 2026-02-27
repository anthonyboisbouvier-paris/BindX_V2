"""
BindX V9 — Supabase JWT Authentication.

Verifies Supabase-issued JWTs (HS256) and extracts user_id.
Uses PyJWT already installed for V8.
"""

from __future__ import annotations

import logging
import os
from typing import Optional
from uuid import UUID

import jwt
from fastapi import Depends, HTTPException, Request

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SUPABASE_JWT_SECRET: str = os.environ.get("SUPABASE_JWT_SECRET", "")


# ---------------------------------------------------------------------------
# JWT helpers
# ---------------------------------------------------------------------------

def _decode_supabase_jwt(token: str) -> dict:
    """Decode and validate a Supabase JWT.

    Raises jwt.InvalidTokenError on any failure.
    """
    return jwt.decode(
        token,
        SUPABASE_JWT_SECRET,
        algorithms=["HS256"],
        audience="authenticated",
    )


def get_v9_user_id(request: Request) -> Optional[str]:
    """Extract user_id from Supabase JWT if present. Returns None if no/invalid token."""
    auth_header = request.headers.get("Authorization", "")
    if not auth_header.startswith("Bearer "):
        return None
    token = auth_header[7:]
    try:
        payload = _decode_supabase_jwt(token)
        return payload.get("sub")
    except jwt.InvalidTokenError:
        return None


# ---------------------------------------------------------------------------
# FastAPI dependency
# ---------------------------------------------------------------------------

async def require_v9_user(request: Request) -> str:
    """FastAPI dependency that requires a valid Supabase JWT.

    Returns the user_id (UUID string from the 'sub' claim).
    Raises 401 if token is missing or invalid.
    """
    auth_header = request.headers.get("Authorization", "")
    if not auth_header.startswith("Bearer "):
        raise HTTPException(status_code=401, detail="Missing authorization header")

    token = auth_header[7:]
    try:
        payload = _decode_supabase_jwt(token)
    except jwt.ExpiredSignatureError:
        raise HTTPException(status_code=401, detail="Token expired")
    except jwt.InvalidTokenError as e:
        raise HTTPException(status_code=401, detail=f"Invalid token: {e}")

    user_id = payload.get("sub")
    if not user_id:
        raise HTTPException(status_code=401, detail="Token missing 'sub' claim")

    return user_id
