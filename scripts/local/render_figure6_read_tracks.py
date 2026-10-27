#!/usr/bin/env python3
"""Render control-versus-elav5 PL-Seq read tracks for Khc-73 and Dscam1."""

from pathlib import Path
import csv
import math
import tomllib

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
from project_paths import repo_root

ROOT = repo_root()
DATA = ROOT / "data/pileup_inputs/figure6"
CONFIG = ROOT / "config/figure3e_pileup_catalog_manifest.tsv"
ATLAS = ROOT / "config/figure3e_locus_atlas_manifest.tsv"
FLY_BED = ROOT / "data/source/Fly_dm6.bed"
STYLE = ROOT / "config/nature_communications.toml"
OUT = ROOT / "figures/main"

GENES = [("Khc-73", "15"), ("Dscam1", "23")]
CONDITIONS = [("control", "Control"), ("elav5", r"$elav^{5}$")]


def load_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def configure_style():
    with STYLE.open("rb") as handle:
        style = tomllib.load(handle)
    colors = style["colors"]
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": style["typography"]["body_pt"],
        "axes.linewidth": style["lines"]["axis_pt"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })
    return style


def read_bedgraph(path, chrom, start, end):
    rows = []
    with path.open() as handle:
        for line in handle:
            c, a, b, value = line.rstrip().split("\t")
            a, b = int(a), int(b)
            if c == chrom and b > start and a < end:
                rows.append((a, b, float(value)))
    return rows


def read_bed12(path):
    reads = []
    with path.open() as handle:
        for line in handle:
            f = line.rstrip().split("\t")
            if len(f) < 12:
                continue
            start = int(f[1])
            sizes = [int(x) for x in f[10].rstrip(",").split(",")]
            offsets = [int(x) for x in f[11].rstrip(",").split(",")]
            reads.append([(start + o, start + o + n) for n, o in zip(sizes, offsets)])
    return reads


def merge_blocks(blocks, view_start, view_end, max_gap=12):
    clipped = [(max(a, view_start), min(b, view_end)) for a, b in blocks
               if b > view_start and a < view_end]
    if not clipped:
        return []
    merged = [list(clipped[0])]
    for start, end in clipped[1:]:
        if start - merged[-1][1] <= max_gap:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [tuple(x) for x in merged]


def fraction_stats(path, fraction):
    rows = [r for r in load_tsv(path) if r["fraction"] == fraction]
    psi = [float(r["PSI"]) for r in rows if r["PSI"]]
    return {
        "informative": sum(int(r["informative"]) for r in rows),
        "mean_psi": sum(psi) / len(psi) if psi else math.nan,
        "replicates": len(psi),
    }


def load_models():
    selected, wanted = {}, set()
    for row in load_tsv(ATLAS):
        if row["gene"] in {g for g, _ in GENES} and row["status"] == "SELECTED" and row["transcript"] != "NA":
            selected.setdefault(row["gene"], []).append((row["category"], row["transcript"]))
            wanted.add(row["transcript"])
    structures = {}
    with FLY_BED.open() as handle:
        for line in handle:
            f = line.rstrip().split("\t")
            if f[3] not in wanted:
                continue
            start = int(f[1])
            sizes = [int(x) for x in f[10].rstrip(",").split(",")]
            offsets = [int(x) for x in f[11].rstrip(",").split(",")]
            structures[f[3]] = [(start + o, start + o + n) for n, o in zip(sizes, offsets)]
    return selected, structures


MODELS, MODEL_BLOCKS = load_models()


def draw_coverage(ax, rows, color, show_ylabel):
    for start, end, value in rows:
        ax.fill_between([start, end], [0, 0], [value, value], color=color, linewidth=0)
    ax.set_ylim(0, 1.03)
    ax.set_yticks([0, 1])
    ax.set_ylabel("Fractional\ncoverage" if show_ylabel else "", fontsize=5.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.tick_params(axis="y", labelsize=5, length=1.8)


def draw_reads(ax, paths, start, end, show_ylabel):
    reads = []
    for path in paths:
        reads.extend(read_bed12(path))
    if len(reads) > 72:
        stride = len(reads) / 72
        reads = [reads[math.floor(i * stride)] for i in range(72)]
    reads.sort(key=lambda x: (x[0][0], x[-1][1]))
    for y, blocks in enumerate(reads):
        blocks = merge_blocks(blocks, start, end)
        if not blocks:
            continue
        ax.plot([blocks[0][0], blocks[-1][1]], [y, y], color="#999999",
                lw=0.45, alpha=0.55, zorder=1)
        for a, b in blocks:
            ax.add_patch(Rectangle((a, y - 0.32), b - a, 0.64,
                                   facecolor="#999999", edgecolor="none", alpha=0.88, zorder=2))
    ax.set_ylim(-1, max(len(reads), 1))
    ax.invert_yaxis()
    ax.set_yticks([])
    ax.set_ylabel("Sampled\nreads" if show_ylabel else "", fontsize=5.5)
    ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)


def draw_models(ax, cfg, gene, ce_number):
    start, end = int(cfg["view_start"]), int(cfg["view_end"])
    ce = (int(cfg["ce_start"]), int(cfg["ce_end"]))
    models = MODELS.get(gene, [])
    for idx, (category, transcript) in enumerate(models):
        y = len(models) - idx - 0.5
        blocks = [(max(a, start), min(b, end)) for a, b in MODEL_BLOCKS[transcript]
                  if b > start and a < end]
        if not blocks:
            continue
        # The intron line is behind the exon blocks.
        ax.plot([blocks[0][0], blocks[-1][1]], [y, y], color="#777777", lw=1.15, zorder=1)
        for (a, b), (c, d) in zip(blocks, blocks[1:]):
            if c - b < 90:
                continue
            x = (b + c) / 2
            dx = min(60, (c - b) / 5)
            p0, p1 = ((x - dx, y), (x + dx, y)) if cfg["strand"] == "+" else ((x + dx, y), (x - dx, y))
            ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=7,
                                         linewidth=1.25, color="#666666", zorder=2))
        for a, b in blocks:
            color = "#C83E4D" if max(a, ce[0]) < min(b, ce[1]) else "#D0D0D0"
            ax.add_patch(Rectangle((a, y - 0.16), b - a, 0.32,
                                   facecolor=color, edgecolor="#555555", linewidth=0.45, zorder=3))
        ax.text(start, y + 0.23, category.replace(" | ", " · "), fontsize=5,
                ha="left", va="bottom")
    ax.set_ylim(-0.15, max(len(models), 1) + 0.08)
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.tick_params(axis="x", labelsize=5, length=2)
    ax.set_xlabel(f"{cfg['chrom']} coordinate (dm6)", fontsize=5.5)


def draw_condition(fig, slot, cfg, gene, ce_number, condition, colors, show_ylabel):
    folder = DATA / gene / condition
    start, end, chrom = int(cfg["view_start"]), int(cfg["view_end"]), cfg["chrom"]
    summary = folder / f"{gene}.pileup_summary.tsv"
    grid = slot.subgridspec(5, 1, height_ratios=[0.72, 1.25, 0.72, 1.25, 1.25], hspace=0.38)
    axes = [fig.add_subplot(grid[i, 0]) for i in range(5)]
    for ax in axes:
        ax.set_xlim(start, end)

    for idx, fraction in enumerate(("long", "short")):
        cov_ax, read_ax = axes[idx * 2], axes[idx * 2 + 1]
        color = colors["long_utr"] if fraction == "long" else colors["short_utr"]
        draw_coverage(cov_ax, read_bedgraph(folder / f"{gene}.{fraction}.mean_fractional_coverage.bedgraph",
                                            chrom, start, end), color, show_ylabel)
        draw_reads(read_ax, sorted(folder.glob(f"{gene}.{fraction}.rep*.sampled.bed")),
                   start, end, show_ylabel)
        stats = fraction_stats(summary, fraction)
        label = "Long 3′UTR" if fraction == "long" else "Short 3′UTR"
        cov_ax.text(0.002, 1.22, label, transform=cov_ax.transAxes, ha="left", va="bottom",
                    fontsize=6, fontweight="bold", color=color, clip_on=False)
        cov_ax.text(0.998, 1.22,
                    f"PSI = {stats['mean_psi']:.3f}   informative reads = {stats['informative']:,}",
                    transform=cov_ax.transAxes, ha="right", va="bottom", fontsize=5.2,
                    color="#333333", clip_on=False)

    draw_models(axes[4], cfg, gene, ce_number)
    return axes


def build():
    style = configure_style()
    colors = style["colors"]
    configs = {row["gene"]: row for row in load_tsv(CONFIG)}
    fig = plt.figure(figsize=(183 / 25.4, 235 / 25.4))
    outer = fig.add_gridspec(2, 2, left=0.09, right=0.98, top=0.925, bottom=0.055,
                             hspace=0.34, wspace=0.14)

    for row, (gene, ce_number) in enumerate(GENES):
        row_axes = []
        for col, (condition, _) in enumerate(CONDITIONS):
            row_axes.append(draw_condition(fig, outer[row, col], configs[gene], gene, ce_number,
                                           condition, colors, show_ylabel=(col == 0)))
        # Anchor both headings to their first coverage track with the same
        # vertical offset so A and B have identical title-to-plot spacing.
        heading_y = row_axes[0][0].get_position().y1 + 0.032
        fig.text(0.018, heading_y, ("A" if row == 0 else "B"),
                 fontsize=11, fontweight="bold", ha="left", va="bottom")
        fig.text(0.052, heading_y, f"{gene} CE{ce_number}",
                 fontsize=8.5, fontweight="bold", ha="left", va="bottom")

        if row == 0:
            column_centers = [
                (axes[0].get_position().x0 + axes[0].get_position().x1) / 2
                for axes in row_axes
            ]

    fig.text(column_centers[0], 0.972, "Control", ha="center", va="top",
             fontsize=9, fontweight="bold")
    fig.text(column_centers[1], 0.972, r"$\mathbf{\mathit{elav}^{5}}$", ha="center", va="top",
             fontsize=9, fontweight="bold")
    fig.text(0.09, 0.012,
             "Coverage is the replicate-weighted mean. Gray rows are deterministic samples of aligned molecules; "
             "PSI and informative counts use all classified primary alignments.",
             fontsize=5.2, color="#555555", ha="left", va="bottom")

    OUT.mkdir(parents=True, exist_ok=True)
    pdf = OUT / "Figure5_control_elav5_read_tracks.pdf"
    png = OUT / "Figure5_control_elav5_read_tracks.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=300)
    plt.close(fig)
    print(pdf)
    print(png)


if __name__ == "__main__":
    build()
