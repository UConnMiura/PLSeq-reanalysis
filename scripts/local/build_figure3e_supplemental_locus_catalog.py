#!/usr/bin/env python3
"""Build the complete 28-event Figure 3E locus catalog in the final panel style."""

from pathlib import Path
import csv
import math
import os
from collections import defaultdict
from functools import lru_cache

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyArrowPatch, Rectangle
from PIL import Image, ImageDraw

import build_nature_communications_figure3e_pilot as base

# Every locus is rendered once in an assembled page and once as a standalone
# panel. Cache parsed tracks so the second rendering does not reread BED files.
base.read_bedgraph = lru_cache(maxsize=None)(base.read_bedgraph)
base.read_bed12 = lru_cache(maxsize=None)(base.read_bed12)

# Curated annotation exception. FBtr0301328 only touches the deposited distal
# interval by two boundary bases, which is not evidence for a long 3'UTR.
# FBtr0343762 is a single-exon overlapping model that does not provide a
# credible CASK cassette-exon/APA structure. The supported display is therefore
# the one conventional short-UTR CASK transcript model.
base.MODELS["CASK"] = [("Short | CE included", "FBtr0301328")]


ROOT = base.ROOT
OUT = ROOT / "figures/figure3e_supplemental_locus_catalog"
PAGE_PNG = OUT / "page_png"
INDIVIDUAL_PDF = OUT / "individual_pdf"
INDIVIDUAL_PNG = OUT / "individual_png"


def load_inputs():
    configs = base.load_tsv(base.CONFIG)
    stats = {row["gene"]: row for row in base.load_tsv(base.STATS)}
    return configs, stats


def significance_badge(ax, q):
    significant = q < 0.05
    ax.text(
        0.995, 1.29, f"BH q = {q:.2g}", transform=ax.transAxes,
        ha="right", va="bottom", fontsize=6.2,
        color="#222222" if significant else "#666666",
        fontweight="bold" if significant else "normal",
        bbox=(dict(boxstyle="round,pad=0.20", facecolor="#F7E4E6",
                   edgecolor="#DDA9AF", linewidth=0.6)
              if significant else None),
        clip_on=False,
    )


def plot_models(ax, cfg, gene, style):
    start, end = int(cfg["view_start"]), int(cfg["view_end"])
    ce = (int(cfg["ce_start"]), int(cfg["ce_end"]))
    models = base.MODELS.get(gene, [])
    for idx, (category, tx) in enumerate(models):
        y = len(models) - idx - 0.5
        blocks = [(max(a, start), min(b, end)) for a, b in base.MODEL_BLOCKS.get(tx, [])
                  if b > start and a < end]
        if not blocks:
            continue
        # Intron first; exon blocks are drawn above it so the line never crosses them.
        ax.plot([blocks[0][0], blocks[-1][1]], [y, y], color="#666666",
                lw=1.15, zorder=1)
        for (_, left_end), (right_start, _) in zip(blocks, blocks[1:]):
            gap = right_start - left_end
            if gap <= 100:
                continue
            spacing = max(300, (end - start) // 10)
            for x in range(left_end + spacing // 2, right_start, spacing):
                dx = min(60, max(20, gap // 12))
                p0, p1 = (((x - dx, y), (x + dx, y)) if cfg["strand"] == "+"
                          else ((x + dx, y), (x - dx, y)))
                ax.add_patch(FancyArrowPatch(
                    p0, p1, arrowstyle="-|>", mutation_scale=6.5,
                    linewidth=1.15, color="#666666", zorder=2,
                ))
        for a, b in blocks:
            is_ce = max(a, ce[0]) < min(b, ce[1])
            ax.add_patch(Rectangle(
                (a, y - 0.16), b - a, 0.32,
                facecolor=(style["colors"]["significant"] if is_ce
                           else style["colors"]["annotation"]),
                edgecolor="#555555", linewidth=0.5, zorder=3,
            ))
        ax.text(start, y + 0.24, category.replace(" | ", " · "),
                fontsize=5.4, ha="left", va="bottom", color="#333333")
    ax.set_ylim(-0.15, max(len(models), 1) + 0.1)
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.tick_params(axis="x", labelsize=5.2, length=2, width=0.8)
    ax.set_xlabel(f"{cfg['chrom']} coordinate (dm6)", fontsize=5.4)


def plot_locus(fig, slot, cfg, stat, style):
    gene = cfg["gene"]
    folder = base.DATA / gene
    start, end, chrom = int(cfg["view_start"]), int(cfg["view_end"]), cfg["chrom"]
    # Deliberate vertical whitespace separates labels from observed reads.
    gs = slot.subgridspec(5, 1, height_ratios=[0.72, 1.10, 0.72, 1.10, 1.18], hspace=0.34)
    axes = [fig.add_subplot(gs[i, 0]) for i in range(5)]
    for ax in axes:
        ax.set_xlim(start, end)

    base.plot_coverage(
        axes[0], base.read_bedgraph(folder / f"{gene}.long.mean_fractional_coverage.bedgraph"),
        chrom, start, end, style["colors"]["long_utr"],
    )
    base.plot_reads(axes[1], sorted(folder.glob(f"{gene}.long.rep*.sampled.bed")),
                    start, end, style["colors"]["reads"])
    base.plot_coverage(
        axes[2], base.read_bedgraph(folder / f"{gene}.short.mean_fractional_coverage.bedgraph"),
        chrom, start, end, style["colors"]["short_utr"],
    )
    base.plot_reads(axes[3], sorted(folder.glob(f"{gene}.short.rep*.sampled.bed")),
                    start, end, style["colors"]["reads"])
    plot_models(axes[4], cfg, gene, style)

    axes[0].text(0.003, 1.12, "Long 3′UTR", transform=axes[0].transAxes,
                 color=style["colors"]["long_utr"], fontsize=6.2,
                 fontweight="bold", ha="left", va="bottom", clip_on=False)
    axes[2].text(0.003, 1.12, "Short 3′UTR", transform=axes[2].transAxes,
                 color=style["colors"]["short_utr"], fontsize=6.2,
                 fontweight="bold", ha="left", va="bottom", clip_on=False)

    delta, q = float(stat["delta_long_minus_short"]), float(stat["BH_q"])
    axes[0].set_title(f"{gene} CE{stat['CE']}   |   ΔPSI = {delta:+.3f}",
                      loc="left", fontsize=7.3, fontweight="bold", pad=15)
    significance_badge(axes[0], q)
    return axes


def add_panel_border(fig, slot):
    box = slot.get_position(fig)
    pad_x, pad_y = 0.008, 0.010
    fig.add_artist(Rectangle(
        (box.x0 - pad_x, box.y0 - pad_y), box.width + 2 * pad_x, box.height + 2 * pad_y,
        transform=fig.transFigure, fill=False, edgecolor="#D8D8D8",
        linewidth=0.65, zorder=20, clip_on=False,
    ))


def render_page(rows, stats, style, page_number, pdf):
    fig = plt.figure(figsize=(7.2, 9.4))
    grid = fig.add_gridspec(2, 2, left=0.075, right=0.975, top=0.955,
                           bottom=0.055, wspace=0.20, hspace=0.20)
    for i, cfg in enumerate(rows):
        slot = grid[i // 2, i % 2]
        plot_locus(fig, slot, cfg, stats[cfg["gene"]], style)
    fig.text(0.975, 0.985, f"page {page_number} of 7",
             fontsize=5.2, color="#666666", ha="right", va="top")
    fig.text(0.975, 0.018,
             "Coverage is the replicate-mean fractional coverage; gray tracks show sampled observed molecules. "
             "Transcript models are representative dm6 annotations.",
             fontsize=5.2, color="#555555", ha="right", va="bottom")
    pdf.savefig(fig)
    path = PAGE_PNG / f"Figure3E_supplemental_catalog_page_{page_number:02d}.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def render_individual(cfg, stat, style):
    fig = plt.figure(figsize=(7.2, 4.65))
    grid = fig.add_gridspec(1, 1, left=0.10, right=0.975, top=0.92, bottom=0.12)
    slot = grid[0, 0]
    plot_locus(fig, slot, cfg, stat, style)
    stem = f"{int(cfg['order']):02d}_{cfg['gene']}_CE{stat['CE']}"
    png, pdf = INDIVIDUAL_PNG / f"{stem}.png", INDIVIDUAL_PDF / f"{stem}.pdf"
    fig.savefig(png, dpi=220)
    fig.savefig(pdf)
    plt.close(fig)
    return png


def contact_sheet(page_paths):
    images = [Image.open(path).convert("RGB") for path in page_paths]
    width = 850
    resized = [im.resize((width, round(im.height * width / im.width)), Image.Resampling.LANCZOS)
               for im in images]
    gap, header = 24, 60
    canvas = Image.new("RGB", (width * 2 + gap, header + sum(
        max(resized[i].height, resized[i + 1].height if i + 1 < len(resized) else 0) + gap
        for i in range(0, len(resized), 2))), "white")
    draw = ImageDraw.Draw(canvas)
    draw.text((18, 18), "Figure 3E supplemental locus catalog · all 28 events", fill="#222222")
    y = header
    for i in range(0, len(resized), 2):
        canvas.paste(resized[i], (0, y))
        if i + 1 < len(resized):
            canvas.paste(resized[i + 1], (width + gap, y))
        y += max(resized[i].height, resized[i + 1].height if i + 1 < len(resized) else 0) + gap
    out = OUT / "Figure3E_supplemental_locus_catalog_contact_sheet.jpg"
    canvas.save(out, quality=91, optimize=True)
    return out


def write_manifest(configs, stats):
    out = OUT / "Figure3E_supplemental_locus_catalog_manifest.tsv"
    with out.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["order", "gene", "cassette_exon", "probe_set", "delta_PSI", "BH_q",
                         "significant_BH_0.05", "page", "panel"])
        for i, cfg in enumerate(configs):
            stat = stats[cfg["gene"]]
            writer.writerow([cfg["order"], cfg["gene"], f"CE{stat['CE']}", cfg["probe_set"],
                             stat["delta_long_minus_short"], stat["BH_q"],
                             str(float(stat["BH_q"]) < 0.05), i // 4 + 1, "abcd"[i % 4]])
    return out


def main():
    style = base.read_style()
    base.configure(style)
    for path in (OUT, PAGE_PNG, INDIVIDUAL_PDF, INDIVIDUAL_PNG):
        path.mkdir(parents=True, exist_ok=True)
    configs, stats = load_inputs()
    multipage = OUT / "Figure3E_supplemental_locus_catalog_all_28.pdf"
    expected_pages = [PAGE_PNG / f"Figure3E_supplemental_catalog_page_{i:02d}.png"
                      for i in range(1, math.ceil(len(configs) / 4) + 1)]
    force = os.environ.get("PLSEQ_FORCE_REBUILD", "0") == "1"
    if not force and multipage.exists() and all(path.exists() for path in expected_pages):
        page_paths = expected_pages
    else:
        page_paths = []
        with PdfPages(multipage) as pdf:
            for start in range(0, len(configs), 4):
                page_paths.append(render_page(configs[start:start + 4], stats, style,
                                              start // 4 + 1, pdf))
    # The assembled seven-page PDF is the archival supplemental figure. Optional
    # standalone exports can be requested for later layout reuse.
    if os.environ.get("PLSEQ_EXPORT_INDIVIDUAL", "0") == "1":
        for cfg in configs:
            render_individual(cfg, stats[cfg["gene"]], style)
    contact = contact_sheet(page_paths)
    manifest = write_manifest(configs, stats)
    print(multipage)
    print(contact)
    print(manifest)


if __name__ == "__main__":
    main()
