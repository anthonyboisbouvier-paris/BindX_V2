"""
BindX V9 — Supabase JWT Authentication.

Verifies Supabase-issued JWTs (ES256) via JWKS endpoint.
"""

from __future__ import annotations

import logging
import os
from typing import Optional

import jwt
from jwt import PyJWKClient
from fastapi import Depends, HTTPException, Request

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SUPABASE_URL: str = os.environ.get("SUPABASE_URL", "")
SUPABASE_JWT_SECRET: str = os.environ.get("SUPABASE_JWT_SECRET", "")

# JWKS client — fetches and caches public keys from Supabase
_jwks_url = f"{SUPABASE_URL}/auth/v1/.well-known/jwks.json" if SUPABASE_URL else ""
_jwks_client: Optional[PyJWKClient] = None

if _jwks_url:
    _jwks_client = PyJWKClient(_jwks_url, cache_keys=True)
    logger.info("V9 auth: JWKS client configured for %s", _jwks_url)


# ---------------------------------------------------------------------------
# JWT helpers
# ---------------------------------------------------------------------------

def _decode_supabase_jwt(token: str) -> dict:
    """Decode and validate a Supabase JWT.

    Tries ES256 (JWKS) first, falls back to HS256 (shared secret).
    Raises jwt.InvalidTokenError on any failure.
    """
    # Try ES256 via JWKS (current Supabase default)
    if _jwks_client:
        try:
            signing_key = _jwks_client.get_signing_key_from_jwt(token)
            return jwt.decode(
                token,
                signing_key.key,
                algorithms=["ES256"],
                audience="authenticated",
            )
        except (jwt.InvalidTokenError, Exception) as e:
            logger.debug("V9 auth: ES256 decode failed (%s), trying HS256 fallback", e)

    # Fallback to HS256 (legacy or if JWKS not configured)
    if SUPABASE_JWT_SECRET:
        return jwt.decode(
            token,
            SUPABASE_JWT_SECRET,
            algorithms=["HS256"],
            audience="authenticated",
        )

    raise jwt.InvalidTokenError("No valid verification method available")


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
        logger.warning("V9 auth: invalid token: %s", e)
        raise HTTPException(status_code=401, detail=f"Invalid token: {e}")

    user_id = payload.get("sub")
    if not user_id:
        raise HTTPException(status_code=401, detail="Token missing 'sub' claim")

    return user_id
