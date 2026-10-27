#!/usr/bin/env python3
"""Verify deterministic PNG outputs against the reviewed publication figures."""

from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path


def main() -> int:
    root = Path(subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip())
    manifest = root / "config/figure_checksums.sha256"
    failures: list[str] = []
    checked = 0

    for line in manifest.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        expected, relative = line.split(maxsplit=1)
        path = root / relative
        if not path.is_file():
            failures.append(f"missing: {relative}")
            continue
        observed = hashlib.sha256(path.read_bytes()).hexdigest()
        checked += 1
        if observed != expected:
            failures.append(f"checksum mismatch: {relative}\n  expected {expected}\n  observed {observed}")

    if failures:
        print("Figure checksum verification failed:")
        for failure in failures:
            print(f"  - {failure}")
        return 1
    print(f"Figure checksum verification passed ({checked} reviewed PNGs).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

