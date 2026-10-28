#!/usr/bin/env python3
"""Create a submission-only manuscript with all main figures at the end."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def extract_figure_boxes(text: str) -> tuple[str, list[str]]:
    lines = text.splitlines()
    body: list[str] = []
    figures: list[str] = []
    current: list[str] | None = None

    for line in lines:
        if current is None and line.strip() == "::: {.figure-box}":
            current = [line]
            continue

        if current is not None:
            current.append(line)
            if line.strip() == ":::":
                figures.append("\n".join(current))
                current = None
            continue

        body.append(line)

    if current is not None:
        raise ValueError("Unclosed {.figure-box} block in manuscript source")
    if not figures:
        raise ValueError("No {.figure-box} blocks found in manuscript source")

    body_text = "\n".join(body).rstrip()
    # Repository-relative downloads do not survive manuscript upload. Omit the
    # TSV call-to-action and repository README note from the submission copy,
    # while leaving both intact in the online reader source.
    body_text = re.sub(
        r" ?\[Download TSV\]\(\.\./supplementary_tables/[^)]+\)\.",
        "",
        body_text,
    )
    body_text = re.sub(
        r"\nTable-specific notes and shared identifier definitions[^\n]+\n",
        "\n",
        body_text,
    )
    # Preserve the visible labels of other local links (for example,
    # Supplementary Figure 1) without retaining unusable local targets.
    body_text = re.sub(r"\[([^]]+)\]\(\.\./[^)]+\)", r"\1", body_text)
    figure_text = "\n\n\\clearpage\n\n".join(figures)
    submission = (
        f"{body_text}\n\n"
        "\\clearpage\n\n"
        f"{figure_text}\n"
    )
    return submission, figures


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    submission, figures = extract_figure_boxes(args.source.read_text())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(submission)
    print(f"Wrote {args.output} with {len(figures)} deferred main figures")


if __name__ == "__main__":
    main()
