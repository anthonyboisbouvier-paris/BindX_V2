"""
DockIt — Celery application configuration.

Connects to Redis as both broker and result backend.
"""

from __future__ import annotations

import os
import sys

from celery import Celery

# Ensure /app is in the Python path for worker processes
app_dir = os.path.dirname(os.path.abspath(__file__))
if app_dir not in sys.path:
    sys.path.insert(0, app_dir)

# Allow overriding the Redis URL via env vars (useful for local dev vs Docker)
REDIS_BROKER = os.environ.get("CELERY_BROKER_URL", "redis://redis:6379/0")
REDIS_BACKEND = os.environ.get("CELERY_RESULT_BACKEND", "redis://redis:6379/1")

celery_app = Celery(
    "dockit",
    broker=REDIS_BROKER,
    backend=REDIS_BACKEND,
    include=["tasks", "tasks_v9"],
)

celery_app.conf.update(
    task_serializer="json",
    result_serializer="json",
    accept_content=["json"],
    task_track_started=True,
    task_acks_late=True,
    worker_prefetch_multiplier=1,
    # Longer visibility timeout for slow docking jobs
    broker_transport_options={"visibility_timeout": 3600},
)
