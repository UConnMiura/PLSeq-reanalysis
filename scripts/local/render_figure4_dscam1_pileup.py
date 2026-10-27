#!/usr/bin/env python3
"""Render a compact Dscam1 long/short 3'UTR read-pileup pilot.

Inputs are deliberately compact: mean fractional-coverage bedGraphs and
deterministically sampled BED12 reads prepared on the HPC.  Replicates are
equally weighted in the coverage tracks.  Read pileups are illustrative,
whereas the adjacent labels report all primary informative reads.
"""

from pathlib import Path
import csv
import re

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle


ROOT = Path("/Users/pedromiura/Documents/codex/PMID37679364_reanalysis")
DATA = ROOT / "data/figure4_pileup/Dscam1"
OUT = ROOT / "figures/figure4_pileup"
SUMMARY = DATA / "Dscam1.pileup_summary.tsv"

CHROM = "chr2R"
VIEW_START, VIEW_END = 7_317_850, 7_322_780
DISTAL = (7_317_924, 7_319_545)
COMMON = (7_319_553, 7_320_764)
CE = (7_321_583, 7_321_649)
FLANK = (7_322_559, 7_322_718)

COL_LONG = "#5B4B9A"
COL_SHORT = "#D18F2E"
COL_CE = "#C83E4D"
COL_SKIP = "#9A9A9A"
COL_OTHER = "#C8C8C8"


def read_bedgraph(path):
    rows = []
    with path.open() as fh:
        for line in fh:
            chrom, start, end, value = line.rstrip().split("\t")
            if chrom == CHROM and int(end) > VIEW_START and int(start) < VIEW_END:
                rows.append((int(start), int(end), float(value)))
    return rows


def read_bed12(path):
    reads = []
    with path.open() as fh:
        for line in fh:
            f = line.rstrip().split("\t")
            start = int(f[1])
            sizes = [int(x) for x in f[10].rstrip(",").split(",")]
            offsets = [int(x) for x in f[11].rstrip(",").split(",")]
            blocks = [(start + off, start + off + size) for size, off in zip(sizes, offsets)]
            category = f[3].split("|")[2]
            reads.append((category, blocks))
    return reads


def summary_stats():
    rows = []
    with SUMMARY.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows.append(row)
    answer = {}
    for fraction in ("long", "short"):
        subset = [r for r in rows if r["fraction"] == fraction]
        answer[fraction] = {
            "included": sum(int(r["included"]) for r in subset),
            "skipped": sum(int(r["skipped"]) for r in subset),
            "informative": sum(int(r["informative"]) for r in subset),
            "psi": [float(r["PSI"]) for r in subset],
        }
    return answer


def draw_coverage(ax, fraction, color):
    rows = read_bedgraph(DATA / f"Dscam1.{fraction}.mean_fractional_coverage.bedgraph")
    for start, end, value in rows:
        ax.fill_between([start, end], [0, 0], [value, value], color=color, linewidth=0)
    ax.set_ylim(0, 1.04)
    ax.set_yticks([0, 0.5, 1])
    ax.set_ylabel("Fractional\ncoverage", fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="both", labelsize=7, length=2)
    ax.set_xticklabels([])


def merge_nearby_blocks(blocks, max_gap=12):
    """Visually suppress tiny alignment gaps caused by ONT indels."""
    clipped = [(max(s, VIEW_START), min(e, VIEW_END)) for s, e in blocks if e > VIEW_START and s < VIEW_END]
    if not clipped:
        return []
    merged = [list(clipped[0])]
    for start, end in clipped[1:]:
        if start - merged[-1][1] <= max_gap:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [tuple(x) for x in merged]


def draw_reads(ax, fraction):
    all_reads = []
    for rep in range(1, 5):
        path = DATA / f"Dscam1.{fraction}.rep{rep}.sampled.bed"
        all_reads.extend((rep, category, blocks) for category, blocks in read_bed12(path))
    # Included first makes the two splice outcomes immediately legible.
    order = {"included": 0, "skipped": 1, "unresolved": 2, "boundary_fail": 3}
    all_reads.sort(key=lambda x: (order.get(x[1], 9), x[0], x[2][0][0]))
    height = 0.70
    for y, (rep, category, blocks) in enumerate(all_reads):
        blocks = merge_nearby_blocks(blocks)
        if not blocks:
            continue
        color = COL_SKIP
        ax.plot([blocks[0][0], blocks[-1][1]], [y, y], color=color, lw=0.25, alpha=0.55, zorder=1)
        for start, end in blocks:
            ax.add_patch(Rectangle((start, y - height / 2), end - start, height,
                                   facecolor=color, edgecolor="none", alpha=0.85, zorder=2))
    ax.set_ylim(-1, max(len(all_reads), 1))
    ax.invert_yaxis()
    ax.set_yticks([])
    ax.set_ylabel("Sampled\nreads", fontsize=8)
    ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)


def draw_model(ax):
    y = 0.5
    ax.plot([DISTAL[0], FLANK[1]], [y, y], color="#6A6A6A", lw=1.1, zorder=1)
    # Remove the visually distracting 8-nt annotation seam between the two UTR
    # segments. It is not meant to represent an intron in this schematic.
    utr_join = COMMON[0]
    for start, end, color, label in [
        (DISTAL[0], utr_join, "#D9D9D9", "distal 3′UTR"),
        (*COMMON, "#BEBEBE", "common 3′UTR"),
        (*CE, COL_CE, "CE23"),
        (*FLANK, "#BEBEBE", "flanking exon"),
    ]:
        ax.add_patch(Rectangle((start, y - 0.15), end - start, 0.30,
                               facecolor=color, edgecolor="#555555", linewidth=0.5))
        ax.text((start + end) / 2, y - 0.25, label, ha="center", va="top", fontsize=7)
    # Dscam1 is on the minus strand. Directional arrowheads belong only on
    # intronic line segments, never on the exon/UTR rectangles.
    introns = [(COMMON[1], CE[0]), (CE[1], FLANK[0])]
    for start, end in introns:
        spacing = 430
        for x in range(start + 180, end - 80, spacing):
            arrow = FancyArrowPatch(
                (x + 55, y), (x - 55, y), arrowstyle="-|>",
                mutation_scale=8.5, linewidth=1.25, color="#666666", zorder=3
            )
            ax.add_patch(arrow)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_xlabel(f"{CHROM} genomic coordinate (dm6)", fontsize=8)
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.tick_params(axis="x", labelsize=7, length=2)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    stats = summary_stats()
    fig = plt.figure(figsize=(11.2, 7.8))
    gs = fig.add_gridspec(5, 1, height_ratios=[1.0, 2.2, 1.0, 2.2, 0.8], hspace=0.34)
    axes = [fig.add_subplot(gs[i, 0]) for i in range(5)]
    for ax in axes:
        ax.set_xlim(VIEW_START, VIEW_END)

    draw_coverage(axes[0], "long", COL_LONG)
    draw_reads(axes[1], "long")
    draw_coverage(axes[2], "short", COL_SHORT)
    draw_reads(axes[3], "short")
    draw_model(axes[4])
    for ax, fraction, title, color in [
        (axes[0], "long", "Long 3′UTR molecules", COL_LONG),
        (axes[2], "short", "Short 3′UTR molecules", COL_SHORT),
    ]:
        s = stats[fraction]
        mean_psi = sum(s["psi"]) / len(s["psi"])
        ax.text(0.006, 1.15, title, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=10, weight="bold", color=color, clip_on=False)
        ax.text(0.994, 1.15,
                f"primary informative reads: {s['informative']:,}  |  mean PSI: {mean_psi:.3f}",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
                color="#333333", clip_on=False)

    fig.suptitle("Dscam1 cassette exon 23 differs sharply between alternative 3′UTR isoforms",
                 x=0.08, y=0.988, ha="left", fontsize=13, weight="bold")
    fig.text(0.08, 0.956,
             "Mean fractional coverage gives each of four embryo replicates equal weight; read rows are deterministic samples.",
             ha="left", fontsize=8.5, color="#444444")
    fig.text(0.08, 0.012,
             "Reads are shown in gray; CE23 inclusion is visible from the aligned blocks. Small alignment gaps (≤12 nt) are merged only for display.",
             ha="left", fontsize=7.5, color="#555555")
    fig.subplots_adjust(left=0.08, right=0.985, top=0.895, bottom=0.075)

    for suffix in ("png", "pdf"):
        path = OUT / f"Dscam1_long_short_read_pileup.{suffix}"
        fig.savefig(path, dpi=300 if suffix == "png" else None, bbox_inches="tight")
        print(path)
    plt.close(fig)


if __name__ == "__main__":
    main()
