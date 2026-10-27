#!/usr/bin/env python3
"""Dscam1 pileup pilot with representative dm6 transcript annotations."""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import render_figure4_dscam1_pileup as base

ROOT = base.ROOT
BED12 = ROOT / "reference/ref/Fly_dm6.bed"
OUT = ROOT / "figures/figure4_pileup"

MODELS = [
    ("FBtr0334801", "Long 3′UTR · CE23 skipped"),
    ("FBtr0089016", "Short 3′UTR · CE23 included"),
]


def load_models():
    wanted = {x[0] for x in MODELS}
    found = {}
    with BED12.open() as fh:
        for line in fh:
            f = line.rstrip().split("\t")
            if f[3] not in wanted:
                continue
            start = int(f[1])
            sizes = [int(x) for x in f[10].rstrip(",").split(",")]
            offsets = [int(x) for x in f[11].rstrip(",").split(",")]
            found[f[3]] = [(start + o, start + o + n) for n, o in zip(sizes, offsets)]
    return found


def draw_direction_arrows(ax, intron_start, intron_end, y):
    if intron_end - intron_start < 100:
        return
    spacing = 430
    for x in range(intron_start + spacing // 2, intron_end, spacing):
        # Dscam1 is on the minus strand, so arrows point left.
        ax.add_patch(FancyArrowPatch(
            (x + 55, y), (x - 55, y), arrowstyle="-|>", mutation_scale=8.5,
            linewidth=1.25, color="#666666", zorder=2
        ))


def draw_transcripts(ax):
    found = load_models()
    for row, (transcript, category) in enumerate(MODELS):
        y = 1.35 - row * 0.72
        blocks = [(max(s, base.VIEW_START), min(e, base.VIEW_END))
                  for s, e in found[transcript]
                  if e > base.VIEW_START and s < base.VIEW_END]
        ax.plot([blocks[0][0], blocks[-1][1]], [y, y], color="#666666", lw=1.1, zorder=1)
        for (left_start, left_end), (right_start, right_end) in zip(blocks, blocks[1:]):
            draw_direction_arrows(ax, left_end, right_start, y)
        for start, end in blocks:
            overlap_ce = max(start, base.CE[0]) < min(end, base.CE[1])
            color = base.COL_CE if overlap_ce else "#D0D0D0"
            ax.add_patch(Rectangle((start, y - 0.14), end - start, 0.28,
                                   facecolor=color, edgecolor="#555555", linewidth=0.55, zorder=3))
        ax.text(base.VIEW_START, y + 0.23, f"{category}  ({transcript})",
                ha="left", va="bottom", fontsize=7.2, color="#333333")
    ax.set_ylim(0.15, 1.75)
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_xlabel(f"{base.CHROM} genomic coordinate (dm6)", fontsize=8)
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.tick_params(axis="x", labelsize=7, length=2)


def main():
    stats = base.summary_stats()
    fig = plt.figure(figsize=(11.2, 8.15))
    gs = fig.add_gridspec(5, 1, height_ratios=[1.0, 2.2, 1.0, 2.2, 1.25], hspace=0.34)
    axes = [fig.add_subplot(gs[i, 0]) for i in range(5)]
    for ax in axes:
        ax.set_xlim(base.VIEW_START, base.VIEW_END)

    base.draw_coverage(axes[0], "long", base.COL_LONG)
    base.draw_reads(axes[1], "long")
    base.draw_coverage(axes[2], "short", base.COL_SHORT)
    base.draw_reads(axes[3], "short")
    draw_transcripts(axes[4])

    for ax, fraction, title, color in [
        (axes[0], "long", "Long 3′UTR molecules", base.COL_LONG),
        (axes[2], "short", "Short 3′UTR molecules", base.COL_SHORT),
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
             "Observed PL-Seq molecules with representative dm6 transcript annotations below.",
             ha="left", fontsize=8.5, color="#444444")
    fig.text(0.08, 0.010,
             "Transcript models are selected representatives from Fly_dm6.bed; red marks annotated CE23.",
             ha="left", fontsize=7.5, color="#555555")
    fig.subplots_adjust(left=0.08, right=0.985, top=0.895, bottom=0.070)

    for suffix in ("png", "pdf"):
        path = OUT / f"Dscam1_representative_dm6_transcripts.{suffix}"
        fig.savefig(path, dpi=300 if suffix == "png" else None, bbox_inches="tight")
        print(path)
    plt.close(fig)


if __name__ == "__main__":
    main()
