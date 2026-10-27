#!/usr/bin/env python3
"""Build a complete Nature Communications-style Figure 3E pilot in code."""

from pathlib import Path
import csv
import math
import tomllib

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

ROOT = Path("/Users/pedromiura/Documents/codex/PMID37679364_reanalysis")
STYLE_PATH = ROOT / "figure_system/nature_communications.toml"
CONFIG = ROOT / "metadata/figure3e_pileup_catalog_manifest.tsv"
STATS = ROOT / "hpc_results/figure3e/figure3e_FASTQ_statistics.tsv"
COUNTS = ROOT / "hpc_results/figure3e/figure3e_reconstructed_counts.tsv"
DATA = ROOT / "data/figure3e_pileup_catalog"
ATLAS = ROOT / "figures/figure3e_loci/atlas/figure3e_locus_atlas_manifest.tsv"
FLY_BED = ROOT / "reference/ref/Fly_dm6.bed"
OUT = ROOT / "figures/nature_communications_pilot"


def load_tsv(path):
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_style():
    with STYLE_PATH.open("rb") as fh:
        return tomllib.load(fh)


def mm(value):
    return value / 25.4


def configure(style):
    colors = style["colors"]
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": style["typography"]["body_pt"],
        "axes.labelsize": style["typography"]["body_pt"],
        "xtick.labelsize": style["typography"]["small_pt"],
        "ytick.labelsize": style["typography"]["small_pt"],
        "axes.linewidth": style["lines"]["axis_pt"],
        "lines.linewidth": style["lines"]["data_pt"],
        "text.color": colors["text"],
        "axes.labelcolor": colors["text"],
        "xtick.color": colors["text"],
        "ytick.color": colors["text"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    })


def panel_label(fig, x, y, label, style):
    fig.text(x, y, label, fontsize=style["typography"]["panel_label_pt"],
             fontweight="bold", ha="left", va="top")


def plot_overview(ax, cfg_rows, stats_rows, style):
    stats = {r["gene"]: r for r in stats_rows}
    genes = [r["gene"] for r in cfg_rows]
    gene_labels = [f"{g} CE{stats[g]['CE']}" for g in genes]
    long = [float(stats[g]["mean_long"]) for g in genes]
    short = [float(stats[g]["mean_short"]) for g in genes]
    qvals = [float(stats[g]["BH_q"]) for g in genes]
    y = list(range(len(genes)))[::-1]
    c = style["colors"]
    for yi, lo, sh, q in zip(y, long, short, qvals):
        ax.plot([sh, lo], [yi, yi], color=c["nonsignificant"], lw=1.0, zorder=1)
        ax.scatter(sh, yi, s=9, facecolor=c["short_utr"], edgecolor="white", linewidth=0.35, zorder=3)
        ax.scatter(lo, yi, s=9, facecolor=c["long_utr"], edgecolor="white", linewidth=0.35, zorder=3)
        significant = q < 0.05
        ax.text(1.085, yi, f"{q:.2g}", ha="center", va="center", fontsize=5,
                color="#9E2F3A" if significant else "#666666",
                fontweight="bold" if significant else "normal",
                bbox=(dict(boxstyle="round,pad=0.16", facecolor="#F7E4E6",
                           edgecolor="#DDA9AF", linewidth=0.55)
                      if significant else None))
    ax.set_yticks(y, gene_labels)
    ax.set_xlim(-0.02, 1.13)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xlabel("Cassette-exon PSI")
    ax.set_title("Splicing differs between long- and short-3′UTR molecules", loc="left",
                 fontsize=7, fontweight="bold", pad=3)
    ax.text(1.085, len(genes) - 0.15, "BH q", ha="center", va="bottom",
            fontsize=5, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=2, width=1)
    ax.scatter([], [], s=10, color=c["long_utr"], label="Long 3′UTR")
    ax.scatter([], [], s=10, color=c["short_utr"], label="Short 3′UTR")
    ax.legend(frameon=False, ncol=2, loc="lower right", fontsize=5,
              handletextpad=0.4, columnspacing=0.8, borderaxespad=0.2)


def read_bedgraph(path):
    rows = []
    with path.open() as fh:
        for line in fh:
            chrom, start, end, value = line.rstrip().split("\t")
            rows.append((chrom, int(start), int(end), float(value)))
    return rows


def read_bed12(path):
    ans = []
    with path.open() as fh:
        for line in fh:
            f = line.rstrip().split("\t")
            if len(f) < 12:
                continue
            start = int(f[1])
            sizes = [int(x) for x in f[10].rstrip(",").split(",")]
            offsets = [int(x) for x in f[11].rstrip(",").split(",")]
            ans.append([(start + o, start + o + n) for n, o in zip(sizes, offsets)])
    return ans


def fraction_summary(path, fraction):
    rows = [r for r in load_tsv(path) if r["fraction"] == fraction]
    psi = [float(r["PSI"]) for r in rows if r["PSI"]]
    return {
        "informative": sum(int(r["informative"]) for r in rows),
        "mean_psi": sum(psi) / len(psi) if psi else math.nan,
    }


def load_models():
    selected = {}
    wanted = set()
    for row in load_tsv(ATLAS):
        if row["status"] == "SELECTED":
            selected.setdefault(row["gene"], []).append((row["category"], row["transcript"]))
            wanted.add(row["transcript"])
    blocks = {}
    with FLY_BED.open() as fh:
        for line in fh:
            f = line.rstrip().split("\t")
            if f[3] not in wanted:
                continue
            start = int(f[1]); sizes = [int(x) for x in f[10].rstrip(",").split(",")]
            offsets = [int(x) for x in f[11].rstrip(",").split(",")]
            blocks[f[3]] = [(start + o, start + o + n) for n, o in zip(sizes, offsets)]
    return selected, blocks


MODELS, MODEL_BLOCKS = load_models()


def plot_coverage(ax, rows, chrom, start, end, color):
    for c, a, b, v in rows:
        if c == chrom and b > start and a < end:
            ax.fill_between([a, b], [0, 0], [v, v], color=color, linewidth=0)
    ax.set_ylim(0, 1.03); ax.set_yticks([0, 1]); ax.set_ylabel("Coverage", fontsize=5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.tick_params(length=1.8, width=1)


def plot_reads(ax, paths, start, end, color):
    reads = []
    for path in paths:
        reads.extend(read_bed12(path))
    # Preserve deterministic sampling but cap display density for final-size clarity.
    if len(reads) > 64:
        stride = len(reads) / 64
        reads = [reads[math.floor(i * stride)] for i in range(64)]
    for y, blocks in enumerate(reads):
        blocks = [(max(a, start), min(b, end)) for a, b in blocks if b > start and a < end]
        if not blocks: continue
        ax.plot([blocks[0][0], blocks[-1][1]], [y, y], color=color, lw=1.0, alpha=0.42)
        for a, b in blocks:
            ax.add_patch(Rectangle((a, y - 0.32), b - a, 0.64,
                                   facecolor=color, edgecolor="none", alpha=0.85))
    ax.set_ylim(-1, max(len(reads), 1)); ax.invert_yaxis(); ax.set_yticks([])
    ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)


def plot_models(ax, cfg, gene, ce_number, style):
    start, end = int(cfg["view_start"]), int(cfg["view_end"])
    ce = (int(cfg["ce_start"]), int(cfg["ce_end"]))
    models = MODELS[gene]
    for idx, (category, tx) in enumerate(models):
        y = len(models) - idx - 0.5
        blocks = [(max(a, start), min(b, end)) for a, b in MODEL_BLOCKS[tx] if b > start and a < end]
        if not blocks: continue
        ax.plot([blocks[0][0], blocks[-1][1]], [y, y], color="#666666", lw=1.0,
                zorder=1)
        for (a, b), (c, d) in zip(blocks, blocks[1:]):
            if c - b > 100:
                x = (b + c) / 2; dx = min(60, (c - b) / 5)
                p0, p1 = ((x-dx,y),(x+dx,y)) if cfg["strand"] == "+" else ((x+dx,y),(x-dx,y))
                ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=6,
                                             linewidth=1.0, color="#666666", zorder=2))
        for a, b in blocks:
            color = style["colors"]["significant"] if max(a, ce[0]) < min(b, ce[1]) else style["colors"]["annotation"]
            ax.add_patch(Rectangle((a, y - 0.16), b - a, 0.32, facecolor=color,
                                   edgecolor="#555555", linewidth=0.5, zorder=3))
        ax.text(start, y + 0.24, category.replace(" | ", " · "), fontsize=5,
                ha="left", va="bottom")
    ax.set_ylim(-0.15, max(len(models), 1) + 0.1); ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.tick_params(axis="x", labelsize=5, length=2, width=1)
    ax.set_xlabel(f"{cfg['chrom']} coordinate (dm6)", fontsize=5)


def plot_locus(fig, slot, cfg, stat, style):
    gene = cfg["gene"]; folder = DATA / gene
    start, end, chrom = int(cfg["view_start"]), int(cfg["view_end"]), cfg["chrom"]
    gs = slot.subgridspec(5, 1, height_ratios=[0.75, 1.2, 0.75, 1.2, 1.15], hspace=0.26)
    axes = [fig.add_subplot(gs[i, 0]) for i in range(5)]
    for ax in axes: ax.set_xlim(start, end)
    plot_coverage(axes[0], read_bedgraph(folder / f"{gene}.long.mean_fractional_coverage.bedgraph"),
                  chrom, start, end, style["colors"]["long_utr"])
    plot_reads(axes[1], sorted(folder.glob(f"{gene}.long.rep*.sampled.bed")), start, end,
               style["colors"]["reads"])
    plot_coverage(axes[2], read_bedgraph(folder / f"{gene}.short.mean_fractional_coverage.bedgraph"),
                  chrom, start, end, style["colors"]["short_utr"])
    plot_reads(axes[3], sorted(folder.glob(f"{gene}.short.rep*.sampled.bed")), start, end,
               style["colors"]["reads"])
    plot_models(axes[4], cfg, gene, stat["CE"], style)
    summary_path = folder / f"{gene}.pileup_summary.tsv"
    fraction_values = {
        "long": fraction_summary(summary_path, "long"),
        "short": fraction_summary(summary_path, "short"),
    }
    for ax, fraction, label, color in (
        (axes[0], "long", "Long 3′UTR", style["colors"]["long_utr"]),
        (axes[2], "short", "Short 3′UTR", style["colors"]["short_utr"]),
    ):
        values = fraction_values[fraction]
        ax.text(0.003, 1.13, label, transform=ax.transAxes,
                color=color, fontsize=6, fontweight="bold",
                ha="left", va="bottom", clip_on=False)
        ax.text(0.997, 1.13,
                f"PSI = {values['mean_psi']:.3f}   informative reads = {values['informative']:,}",
                transform=ax.transAxes, color="#333333", fontsize=5.2,
                ha="right", va="bottom", clip_on=False)
    axes[0].set_title(f"{gene} CE{stat['CE']}  |  ΔPSI = {float(stat['delta_long_minus_short']):+.3f}; "
                      f"BH q = {float(stat['BH_q']):.2g}", loc="left", fontsize=7,
                      fontweight="bold", pad=12)


def build():
    style = read_style(); configure(style)
    configs = load_tsv(CONFIG); cfg = {r["gene"]: r for r in configs}
    stats_rows = load_tsv(STATS); stats = {r["gene"]: r for r in stats_rows}
    # Figure 2 already contains the complete 28-gene statistical summary.
    # Figure 3 therefore contains only the complementary read-level evidence.
    fig = plt.figure(figsize=(mm(style["canvas"]["width_mm"]), mm(230.0)))
    outer = fig.add_gridspec(2, 2, hspace=0.17, wspace=0.16,
                             left=0.10, right=0.965, top=0.965, bottom=0.055)
    plot_locus(fig, outer[0, 0], cfg["Dys"], stats["Dys"], style)
    plot_locus(fig, outer[0, 1], cfg["Dscam1"], stats["Dscam1"], style)
    plot_locus(fig, outer[1, 0], cfg["Khc-73"], stats["Khc-73"], style)
    plot_locus(fig, outer[1, 1], cfg["Eip63E"], stats["Eip63E"], style)

    panel_label(fig, 0.018, 0.985, "A", style)
    panel_label(fig, 0.505, 0.985, "B", style)
    panel_label(fig, 0.018, 0.505, "C", style)
    panel_label(fig, 0.505, 0.505, "D", style)
    OUT.mkdir(parents=True, exist_ok=True)
    pdf = OUT / "Figure3E_Nature_Communications_pilot.pdf"
    png = OUT / "Figure3E_Nature_Communications_pilot.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=style["export"]["preview_dpi"])
    plt.close(fig)
    print(pdf); print(png)


if __name__ == "__main__":
    build()
