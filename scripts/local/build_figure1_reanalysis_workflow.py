#!/usr/bin/env python3
"""Create manuscript Figure 1: reproducible reanalysis workflow and timeline."""

from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from project_paths import repo_root

ROOT = repo_root()
OUT = ROOT / "figures/main"

INK = "#242424"
MUTED = "#676767"
LIGHT = "#E8E8E8"
PURPLE = "#5B4B9A"
ORANGE = "#D69127"
RED = "#C94A57"
BLUE = "#3D718E"
GREEN = "#4E8062"
PALE = "#F7F7F7"


def setup():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "text.color": INK,
        "axes.labelcolor": INK,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })


def rounded(ax, x, y, w, h, face, edge=LIGHT, lw=0.8, radius=0.018):
    patch = FancyBboxPatch((x, y), w, h,
                           boxstyle=f"round,pad=0.008,rounding_size={radius}",
                           facecolor=face, edgecolor=edge, linewidth=lw)
    ax.add_patch(patch)
    return patch


def arrow(ax, x0, y0, x1, y1, color="#9A9A9A", lw=1.2):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                                 mutation_scale=9, linewidth=lw, color=color,
                                 shrinkA=2, shrinkB=2))


def panel_label(fig, x, y, label):
    fig.text(x, y, label, fontsize=15, fontweight="bold", ha="left", va="top")


def workflow(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    stages = [
        ("Published study", "Paper, supplements,\nSRA accessions and BEDs", "#F1F1F1"),
        ("Recover", "Map samples; audit files;\nreconstruct missing context", "#EAF1F5"),
        ("Recompute", "FASTQ -> dm6 BAM;\nlong/short molecule filters", "#EEEAF7"),
        ("Reanalyze", "Replicate-level effects;\nBH-adjusted q values", "#F8ECEE"),
        ("Reuse", "Read-level figures, Quarto,\nGitHub, Zenodo and preprint", "#EAF3ED"),
    ]
    left, gap = 0.015, 0.018
    width = (0.97 - gap * 4) / 5
    for i, (title, text, color) in enumerate(stages):
        x = left + i * (width + gap)
        rounded(ax, x, 0.27, width, 0.55, color, edge="#D3D3D3")
        ax.text(x + width / 2, 0.68, title, ha="center", va="center",
                fontsize=8.5, fontweight="bold")
        ax.text(x + width / 2, 0.47, text, ha="center", va="center",
                fontsize=6.8, color="#3F3F3F", linespacing=1.35)
        if i < 4:
            arrow(ax, x + width + 0.003, 0.545,
                  x + width + gap - 0.003, 0.545)
    ax.text(0.015, 0.94, "From archived publication to reusable analysis",
            fontsize=10.2, fontweight="bold", ha="left", va="top")
    ax.text(0.015, 0.10,
            "Human review at every stage: biological questions, validation, interpretation and final presentation",
            fontsize=6.8, color=MUTED, ha="left", va="center")


def improvements(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.01, 0.97, "What the reanalysis adds", fontsize=10.2,
            fontweight="bold", ha="left", va="top")
    rows = [
        ("Evidence", "Summary PSI", "Replicate-resolved PSI + observed molecules", PURPLE),
        ("Statistics", "Nominal significance / stars", "Effect sizes + BH-adjusted q values", RED),
        ("Interpretation", "Pooled changes", "Genotype-by-3'UTR interaction", BLUE),
        ("Presentation", "Manually assembled figures", "Code-generated main and locus-catalog figures", ORANGE),
        ("Reproducibility", "Static publication outputs", "Apptainer, Slurm, manifests, checksums and Quarto", GREEN),
    ]
    ax.text(0.20, 0.84, "Published representation", fontsize=6.8,
            fontweight="bold", color=MUTED, ha="center")
    ax.text(0.72, 0.84, "Reusable reanalysis product", fontsize=6.8,
            fontweight="bold", color=MUTED, ha="center")
    for i, (category, before, after, color) in enumerate(rows):
        y = 0.69 - i * 0.135
        ax.text(0.01, y, category, fontsize=7.1, fontweight="bold",
                ha="left", va="center", color=color)
        rounded(ax, 0.11, y - 0.045, 0.31, 0.09, "#F3F3F3", edge="#DDDDDD", radius=0.012)
        rounded(ax, 0.54, y - 0.045, 0.44, 0.09, "white", edge=color, lw=0.9, radius=0.012)
        ax.text(0.265, y, before, fontsize=6.3, color="#555555", ha="center", va="center")
        ax.text(0.76, y, after, fontsize=6.3, color=INK, ha="center", va="center")
        arrow(ax, 0.435, y, 0.525, y, color=color, lw=1.0)


def timeline(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.01, 0.96, "A four-day reanalysis and dissemination sprint",
            fontsize=10.2, fontweight="bold", ha="left", va="top")
    ax.text(0.99, 0.96, "24-27 October 2026", fontsize=7.1,
            color=MUTED, ha="right", va="top")

    x0, x1, yline = 0.06, 0.94, 0.54
    ax.plot([x0, x1], [yline, yline], color="#B8B8B8", lw=1.3)
    days = [
        (x0, "Day 1", "Paper + supplement audit\nPhase 1/2 tables\nFASTQ reconstruction pilot", BLUE),
        (x0 + (x1-x0)/3, "Day 2", "HPC alignment + filtering\nReplicate-aware statistics\nFigures 2-4 redesign", PURPLE),
        (x0 + 2*(x1-x0)/3, "Day 3", "Figure 5 + locus catalog\nPublication-style figures\nManuscript draft", GREEN),
        (x1, "Day 4", "GitHub repository published\nZenodo release archived\nbioRxiv submission", RED),
    ]
    for i, (x, day, text, color) in enumerate(days):
        ax.scatter([x], [yline], s=52, color=color, edgecolor="white", linewidth=1.0, zorder=3)
        ha = "left" if i == 0 else ("right" if i == len(days) - 1 else "center")
        tx = x if i not in (0, len(days) - 1) else x + (0.005 if i == 0 else -0.005)
        ax.text(tx, yline + 0.16, day, fontsize=7.4, fontweight="bold", ha=ha, color=color)
        ax.text(tx, yline - 0.13, text, fontsize=6.3, ha=ha, va="top",
                color="#444444", linespacing=1.35)

    # Timing qualification is retained in the manuscript legend rather than
    # rendered as miniature text inside the figure.


def build():
    setup(); OUT.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(7.2, 8.15))
    grid = fig.add_gridspec(3, 1, height_ratios=[0.92, 1.28, 0.92], hspace=0.19,
                           left=0.065, right=0.975, top=0.94, bottom=0.055)
    axes = [fig.add_subplot(grid[i, 0]) for i in range(3)]
    workflow(axes[0]); improvements(axes[1]); timeline(axes[2])
    panel_label(fig, 0.018, 0.946, "A")
    panel_label(fig, 0.018, 0.653, "B")
    panel_label(fig, 0.018, 0.300, "C")
    fig.suptitle("Reanalysis converts an archived PL-Seq study into a reusable research resource",
                 x=0.065, y=0.985, ha="left", fontsize=12.5, fontweight="bold")
    pdf = OUT / "Figure1_reanalysis_workflow_and_timeline.pdf"
    png = OUT / "Figure1_reanalysis_workflow_and_timeline.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=240)
    plt.close(fig)
    print(pdf); print(png)


if __name__ == "__main__":
    build()
