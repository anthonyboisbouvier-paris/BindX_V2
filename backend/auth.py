"""
DockIt -- Authentication utilities.

Provides JWT token creation/verification and password hashing.
Auth is OPTIONAL -- all endpoints work without a token.

V7: Phase C (Projects + Accounts).
"""

from __future__ import annotations

import logging
import os
from datetime import datetime, timedelta
from typing import Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# JWT configuration
# ---------------------------------------------------------------------------

JWT_SECRET: str = os.environ.get("DOCKIT_JWT_SECRET", "dockit-dev-secret-change-in-prod")
JWT_ALGORITHM: str = "HS256"
JWT_EXPIRE_HOURS: int = 72


# ---------------------------------------------------------------------------
# Password hashing
# ---------------------------------------------------------------------------

def hash_password(password: str) -> str:
    """Hash a password using bcrypt.

    Falls back to SHA-256 if bcrypt is not installed (dev only, not
    production-safe).

    Parameters
    ----------
    password : str
        The plaintext password.

    Returns
    -------
    str
        The hashed password string.
    """
    try:
        import bcrypt
        return bcrypt.hashpw(password.encode("utf-8"), bcrypt.gensalt()).decode("utf-8")
    except ImportError:
        # Fallback: simple SHA256 (dev only, not production-safe)
        import hashlib
        logger.warning("bcrypt not installed, using SHA256 fallback (NOT production-safe)")
        return "sha256:" + hashlib.sha256(password.encode("utf-8")).hexdigest()


def verify_password(password: str, password_hash: str) -> bool:
    """Verify a password against its hash.

    Supports both bcrypt hashes and the SHA-256 fallback format.

    Parameters
    ----------
    password : str
        The plaintext password to check.
    password_hash : str
        The stored hash (bcrypt or ``sha256:...`` format).

    Returns
    -------
    bool
        ``True`` if the password matches.
    """
    try:
        import bcrypt
        if password_hash.startswith("sha256:"):
            import hashlib
            return password_hash == "sha256:" + hashlib.sha256(password.encode("utf-8")).hexdigest()
        return bcrypt.checkpw(password.encode("utf-8"), password_hash.encode("utf-8"))
    except ImportError:
        import hashlib
        if password_hash.startswith("sha256:"):
            return password_hash == "sha256:" + hashlib.sha256(password.encode("utf-8")).hexdigest()
        return False


# ---------------------------------------------------------------------------
# JWT tokens
# ---------------------------------------------------------------------------

def create_token(user_id: str) -> str:
    """Create a JWT token for a user.

    Parameters
    ----------
    user_id : str
        The user UUID to embed in the token ``sub`` claim.

    Returns
    -------
    str
        An encoded JWT string.

    Raises
    ------
    ImportError
        If neither ``PyJWT`` nor ``python-jose`` is installed.
    """
    try:
        import jwt
    except ImportError:
        try:
            import jose.jwt as jwt  # type: ignore[no-redef]
        except ImportError:
            logger.error("No JWT library available (install PyJWT or python-jose)")
            raise ImportError("Install PyJWT: pip install PyJWT")

    payload = {
        "sub": user_id,
        "exp": datetime.utcnow() + timedelta(hours=JWT_EXPIRE_HOURS),
        "iat": datetime.utcnow(),
    }
    return jwt.encode(payload, JWT_SECRET, algorithm=JWT_ALGORITHM)


def decode_token(token: str) -> Optional[str]:
    """Decode a JWT token and return the ``user_id``, or ``None`` if invalid.

    Parameters
    ----------
    token : str
        The raw JWT string (without ``Bearer `` prefix).

    Returns
    -------
    str or None
        The user ID from the ``sub`` claim, or ``None`` on any error
        (expired, tampered, missing library, etc.).
    """
    try:
        import jwt
    except ImportError:
        try:
            import jose.jwt as jwt  # type: ignore[no-redef]
        except ImportError:
            return None

    try:
        payload = jwt.decode(token, JWT_SECRET, algorithms=[JWT_ALGORITHM])
        return payload.get("sub")
    except Exception:
        return None
