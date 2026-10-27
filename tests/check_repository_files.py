#!/usr/bin/env python3
"""Reject sequencing files, large binaries, and unsafe paths before commit."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


FORBIDDEN_SUFFIXES = (
    ".bam",
    ".bai",
    ".cram",
    ".crai",
    ".sam",
    ".fastq",
    ".fastq.gz",
    ".fq",
    ".fq.gz",
    ".fa",
    ".fa.gz",
    ".fasta",
    ".fasta.gz",
    ".fna",
    ".sif",
    ".img",
    ".simg",
)
FORBIDDEN_PARTS = {
    ".quarto",
    "__pycache__",
    "node_modules",
    "scratch",
    "work",
}
IGNORED_WORKING_PARTS = {".git", ".quarto", "__pycache__", "node_modules"}
MAX_BYTES = 25 * 1024 * 1024


def staged_paths(root: Path) -> list[Path]:
    result = subprocess.run(
        ["git", "diff", "--cached", "--name-only", "--diff-filter=ACMR"],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    )
    return [root / line for line in result.stdout.splitlines() if line]


def working_paths(root: Path) -> list[Path]:
    return [
        path
        for path in root.rglob("*")
        if path.is_file() and not IGNORED_WORKING_PARTS.intersection(path.relative_to(root).parts)
    ]


def violations(paths: list[Path], root: Path) -> list[str]:
    problems: list[str] = []
    for path in paths:
        relative = path.relative_to(root)
        lowered = relative.as_posix().lower()
        if any(lowered.endswith(suffix) for suffix in FORBIDDEN_SUFFIXES):
            problems.append(f"forbidden file type: {relative}")
        if FORBIDDEN_PARTS.intersection(relative.parts):
            problems.append(f"forbidden generated path: {relative}")
        if path.exists() and path.stat().st_size > MAX_BYTES:
            mib = path.stat().st_size / (1024 * 1024)
            problems.append(f"file exceeds 25 MiB ({mib:.1f} MiB): {relative}")
    return problems


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--staged", action="store_true", help="inspect only files staged in Git")
    args = parser.parse_args()

    root = Path(subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip())
    paths = staged_paths(root) if args.staged else working_paths(root)
    problems = violations(paths, root)
    if problems:
        print("Repository safety check failed:", file=sys.stderr)
        for problem in problems:
            print(f"  - {problem}", file=sys.stderr)
        return 1
    print(f"Repository safety check passed ({len(paths)} files inspected).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
