#!/usr/bin/env python3
"""Verify the small Level 3 programs and interval inputs by SHA-256."""

from __future__ import annotations

import hashlib
from pathlib import Path


def main() -> int:
    source = Path(__file__).resolve().parents[1] / "data" / "source"
    failures = []
    checked = 0
    for line in (source / "SHA256SUMS").read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        expected, name = line.split(maxsplit=1)
        path = source / name
        if not path.is_file():
            failures.append(f"missing: {name}")
            continue
        observed = hashlib.sha256(path.read_bytes()).hexdigest()
        checked += 1
        if observed != expected:
            failures.append(f"checksum mismatch: {name}")
    if failures:
        print("Level 3 source checksum verification failed:")
        for failure in failures:
            print(f"  - {failure}")
        return 1
    print(f"Level 3 source checksum verification passed ({checked} files).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
