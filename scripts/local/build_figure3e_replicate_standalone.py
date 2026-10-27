#!/usr/bin/env python3
"""Standalone replicate-resolved Figure 3E in the current visual style."""

from pathlib import Path
from collections import defaultdict
import csv
import math
import tomllib

import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT = Path("/Users/pedromiura/Documents/codex/PMID37679364_reanalysis")
STYLE = ROOT / "figure_system/nature_communications.toml"
CONFIG = ROOT / "metadata/figure3e_pileup_catalog_manifest.tsv"
COUNTS = ROOT / "hpc_results/figure3e/figure3e_reconstructed_counts.tsv"
STATS = ROOT / "hpc_results/figure3e/figure3e_FASTQ_statistics.tsv"
OUT = ROOT / "figures/figure3e_replicate_standalone"


def load_tsv(path):
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def number(x):
    try:
        value = float(x)
        return value if math.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def main():
    with STYLE.open("rb") as fh:
        style = tomllib.load(fh)
    colors = style["colors"]
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 6,
        "axes.labelsize": 7,
        "axes.linewidth": 1.0,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })

    cfg_rows = load_tsv(CONFIG)
    cfg = {r["gene"]: r for r in cfg_rows}
    stats = {r["gene"]: r for r in load_tsv(STATS)}
    genes = [r["gene"] for r in cfg_rows]

    values = defaultdict(dict)
    for r in load_tsv(COUNTS):
        gene = r["gene"]
        if gene not in cfg or r["probe_set"] != cfg[gene]["probe_set"]:
            continue
        if r["isoform"] not in {"Long", "Short"}:
            continue
        psi = number(r["PSI"])
        if psi is not None:
            values[(gene, r["replicate"])][r["isoform"]] = psi

    width_mm, height_mm = 180, 205
    fig, axes = plt.subplots(1, 2, figsize=(width_mm / 25.4, height_mm / 25.4),
                             gridspec_kw={"wspace": 0.40})
    halves = [genes[:14], genes[14:]]

    for ax, subset in zip(axes, halves):
        ypos = list(range(len(subset)))[::-1]
        for y, gene in zip(ypos, subset):
            pairs = []
            for (g, rep), row in values.items():
                if g == gene and "Long" in row and "Short" in row:
                    pairs.append((int(rep), row["Long"], row["Short"]))
            pairs.sort()
            n = len(pairs)
            offsets = ([0.0] if n == 1 else
                       [(-0.21 + i * 0.42 / (n - 1)) for i in range(n)])
            for offset, (rep, long_psi, short_psi) in zip(offsets, pairs):
                yy = y + offset
                ax.plot([short_psi, long_psi], [yy, yy], color="#B2B2B2",
                        lw=1.0, alpha=0.78, zorder=1)
                ax.scatter(long_psi, yy, s=9, color=colors["long_utr"],
                           edgecolor="white", linewidth=0.3, zorder=3)
                ax.scatter(short_psi, yy, s=9, color=colors["short_utr"],
                           edgecolor="white", linewidth=0.3, zorder=3)
            q = float(stats[gene]["BH_q"])
            significant = q < 0.05
            ax.text(1.09, y, f"{q:.2g}", ha="center", va="center", fontsize=5.5,
                    color="#222222" if significant else "#666666",
                    fontweight="bold" if significant else "normal",
                    bbox=(dict(boxstyle="round,pad=0.16", facecolor="#F7E4E6",
                               edgecolor="#DDA9AF", linewidth=0.55)
                          if significant else None))

        labels = [f"{g} CE{stats[g]['CE']}" for g in subset]
        ax.set_yticks(ypos, labels)
        ax.set_xlim(-0.03, 1.15)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
        ax.set_xlabel("Cassette-exon PSI")
        ax.text(1.09, len(subset) - 0.12, "BH q", ha="center", va="bottom",
                fontsize=6, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(length=2.5, width=1)

    axes[0].scatter([], [], s=16, color=colors["long_utr"], label="Long 3′UTR")
    axes[0].scatter([], [], s=16, color=colors["short_utr"], label="Short 3′UTR")
    fig.legend(frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, 0.935),
               fontsize=6, handletextpad=0.4, columnspacing=1.0)
    fig.suptitle("Replicate-resolved splicing differences between long- and short-3′UTR molecules",
                 fontsize=9, fontweight="bold", y=0.977)
    fig.text(0.5, 0.947,
             "Each thin line joins long- and short-3′UTR PSI measurements from the same biological replicate.",
             ha="center", va="top", fontsize=6, color="#555555")
    fig.subplots_adjust(left=0.15, right=0.965, top=0.89, bottom=0.08)

    OUT.mkdir(parents=True, exist_ok=True)
    pdf = OUT / "Figure3E_replicate_resolved_standalone.pdf"
    png = OUT / "Figure3E_replicate_resolved_standalone.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=300)
    plt.close(fig)
    print(pdf); print(png)


if __name__ == "__main__":
    main()
