"""Shared portable paths for local PL-Seq analysis programs."""

from __future__ import annotations

import os
from pathlib import Path


def repo_root() -> Path:
    """Return the repository root, optionally overridden by ``PLSEQ_ROOT``."""
    override = os.environ.get("PLSEQ_ROOT")
    if override:
        return Path(override).expanduser().resolve()
    return Path(__file__).resolve().parents[2]

