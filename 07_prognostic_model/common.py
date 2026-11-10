"""Input fingerprinting used by release/build_release.py."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path


def sha256(path: str | Path, max_bytes: int | None = None) -> str:
    """Hash a file, optionally only its first `max_bytes` bytes (plus the file size)."""
    h = hashlib.sha256()
    read = 0
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(1 << 22)
            if not chunk:
                break
            h.update(chunk)
            read += len(chunk)
            if max_bytes is not None and read >= max_bytes:
                break
    tag = "sha256" if max_bytes is None else f"sha256-first{max_bytes}B"
    return f"{tag}:{h.hexdigest()}:size={os.path.getsize(path)}"


def input_fingerprint(path: str | Path) -> dict:
    """Full hash for small files, head-hash + size for anything over 1 GB."""
    size = os.path.getsize(path)
    limit = None if size < (1 << 30) else (1 << 26)
    return {"path": str(path), "size_bytes": size, "hash": sha256(path, limit)}
