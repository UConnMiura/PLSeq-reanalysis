#!/usr/bin/env python3
"""Render the complete 28-event PL-Seq long/short read-pileup catalog."""

from pathlib import Path
from collections import defaultdict
import csv
import math
import os

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
from PIL import Image, ImageDraw

ROOT = Path("/Users/pedromiura/Documents/codex/PMID37679364_reanalysis")
DATA = ROOT / "data/figure3e_pileup_catalog"
CONFIG = ROOT / "metadata/figure3e_pileup_catalog_manifest.tsv"
STATS = ROOT / "hpc_results/figure3e/figure3e_FASTQ_statistics.tsv"
OUT = ROOT / "figures/figure3e_pileup_catalog"
ATLAS_MODELS = ROOT / "figures/figure3e_loci/atlas/figure3e_locus_atlas_manifest.tsv"
FLY_BED = ROOT / "reference/ref/Fly_dm6.bed"
USE_TRANSCRIPT_MODELS = os.environ.get("PLSEQ_TRANSCRIPT_MODELS", "0") == "1"
if USE_TRANSCRIPT_MODELS:
    OUT = ROOT / "figures/figure3e_pileup_catalog_representative_transcripts"

COL_LONG = "#5B4B9A"
COL_SHORT = "#D18F2E"
COL_CE = "#C83E4D"
COL_READ = "#999999"


def load_tsv(path):
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fnum(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def read_bedgraph(path, chrom, start, end):
    rows = []
    with path.open() as fh:
        for line in fh:
            c, a, b, value = line.rstrip().split("\t")
            a, b = int(a), int(b)
            if c == chrom and b > start and a < end:
                rows.append((a, b, float(value)))
    return rows


def read_bed12(path):
    reads = []
    with path.open() as fh:
        for line in fh:
            f = line.rstrip().split("\t")
            if len(f) < 12:
                continue
            start = int(f[1])
            sizes = [int(x) for x in f[10].rstrip(",").split(",")]
            offsets = [int(x) for x in f[11].rstrip(",").split(",")]
            reads.append([(start + o, start + o + n) for n, o in zip(sizes, offsets)])
    return reads


def merge_blocks(blocks, view_start, view_end, max_gap=12):
    blocks = [(max(s, view_start), min(e, view_end)) for s, e in blocks
              if e > view_start and s < view_end]
    if not blocks:
        return []
    ans = [list(blocks[0])]
    for start, end in blocks[1:]:
        if start - ans[-1][1] <= max_gap:
            ans[-1][1] = max(ans[-1][1], end)
        else:
            ans.append([start, end])
    return [tuple(x) for x in ans]


def fraction_summary(rows, fraction):
    x = [r for r in rows if r["fraction"] == fraction]
    psi = [fnum(r["PSI"]) for r in x]
    psi = [v for v in psi if math.isfinite(v)]
    return {
        "informative": sum(int(r["informative"]) for r in x),
        "psi": sum(psi) / len(psi) if psi else math.nan,
        "n_psi": len(psi),
        "n_reps": len(x),
        "mismatch": sum(r["count_match"].lower() != "true" for r in x),
    }


def draw_coverage(ax, rows, color):
    for start, end, value in rows:
        ax.fill_between([start, end], [0, 0], [value, value], color=color, linewidth=0)
    ax.set_ylim(0, 1.04)
    ax.set_yticks([0, 0.5, 1])
    ax.set_ylabel("Fractional\ncoverage", fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="both", labelsize=7, length=2)
    ax.tick_params(axis="x", labelbottom=False)


def draw_reads(ax, paths, view_start, view_end):
    reads = []
    for path in paths:
        reads.extend(read_bed12(path))
    reads.sort(key=lambda b: (b[0][0] if b else view_start, b[-1][1] if b else view_start))
    for y, blocks in enumerate(reads):
        blocks = merge_blocks(blocks, view_start, view_end)
        if not blocks:
            continue
        ax.plot([blocks[0][0], blocks[-1][1]], [y, y], color=COL_READ,
                lw=0.24, alpha=0.50, zorder=1)
        for start, end in blocks:
            ax.add_patch(Rectangle((start, y - 0.34), end - start, 0.68,
                                   facecolor=COL_READ, edgecolor="none", alpha=0.86, zorder=2))
    ax.set_ylim(-1, max(len(reads), 1))
    ax.invert_yaxis()
    ax.set_yticks([])
    ax.set_ylabel("Sampled\nreads", fontsize=8)
    ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)


def union_intervals(intervals):
    intervals = sorted(intervals)
    ans = []
    for start, end in intervals:
        if not ans or start > ans[-1][1]:
            ans.append([start, end])
        else:
            ans[-1][1] = max(ans[-1][1], end)
    return ans


def load_representative_models():
    selected = defaultdict(list)
    wanted = set()
    for row in load_tsv(ATLAS_MODELS):
        if row["status"] == "SELECTED" and row["transcript"] != "NA":
            selected[row["gene"]].append((row["category"], row["transcript"]))
            wanted.add(row["transcript"])
    structures = {}
    with FLY_BED.open() as fh:
        for line in fh:
            f = line.rstrip().split("\t")
            if f[3] not in wanted:
                continue
            start = int(f[1])
            sizes = [int(x) for x in f[10].rstrip(",").split(",")]
            offsets = [int(x) for x in f[11].rstrip(",").split(",")]
            structures[f[3]] = [(start + o, start + o + n) for n, o in zip(sizes, offsets)]
    return selected, structures


REPRESENTATIVES, TRANSCRIPT_BLOCKS = load_representative_models()


def draw_transcript_models(ax, cfg, published_ce):
    gene = cfg["gene"]
    view_start, view_end = int(cfg["view_start"]), int(cfg["view_end"])
    models = REPRESENTATIVES[gene]
    n = max(len(models), 1)
    top = n - 0.35
    for row, (category, transcript) in enumerate(models):
        y = top - row
        blocks = [(max(s, view_start), min(e, view_end))
                  for s, e in TRANSCRIPT_BLOCKS[transcript]
                  if e > view_start and s < view_end]
        if not blocks:
            continue
        ax.plot([blocks[0][0], blocks[-1][1]], [y, y], color="#666666", lw=1.1, zorder=1)
        for (_, left_end), (right_start, _) in zip(blocks, blocks[1:]):
            if right_start - left_end < 80:
                continue
            spacing = max(320, (view_end - view_start) // 12)
            for x in range(left_end + spacing // 2, right_start, spacing):
                dx = min(55, max(18, (right_start - left_end) // 12))
                p0, p1 = ((x - dx, y), (x + dx, y)) if cfg["strand"] == "+" else ((x + dx, y), (x - dx, y))
                ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=8.5,
                                             linewidth=1.25, color="#666666", zorder=2))
        ce_start, ce_end = int(cfg["ce_start"]), int(cfg["ce_end"])
        for start, end in blocks:
            overlap_ce = max(start, ce_start) < min(end, ce_end)
            color = COL_CE if overlap_ce else "#D0D0D0"
            ax.add_patch(Rectangle((start, y - 0.18), end - start, 0.36,
                                   facecolor=color, edgecolor="#555555", linewidth=0.55, zorder=3))
        ax.text(view_start, y + 0.27, f"{category.replace(' | ', ' · ')}  ({transcript})",
                ha="left", va="bottom", fontsize=6.8, color="#333333")
    ax.set_ylim(-0.1, max(n, 1) + 0.15)
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_xlabel(f"{cfg['chrom']} genomic coordinate (dm6)", fontsize=8)
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.tick_params(axis="x", labelsize=7, length=2)


def draw_model(ax, cfg, published_ce):
    y = 0.52
    strand = cfg["strand"]
    view_start, view_end = int(cfg["view_start"]), int(cfg["view_end"])
    intervals = {
        "distal 3′UTR": [int(cfg["distal_start"]), int(cfg["distal_end"])],
        "common 3′UTR": [int(cfg["common_start"]), int(cfg["common_end"])],
        f"CE{published_ce}": [int(cfg["ce_start"]), int(cfg["ce_end"])],
        "upstream exon": [int(cfg["up_start"]), int(cfg["up_end"])],
        "downstream exon": [int(cfg["down_start"]), int(cfg["down_end"])],
    }
    # Make nearly adjacent UTR annotations touch in the schematic.
    d, c = intervals["distal 3′UTR"], intervals["common 3′UTR"]
    gap = max(d[0], c[0]) - min(d[1], c[1])
    if 0 < gap <= 25:
        join = (min(d[1], c[1]) + max(d[0], c[0])) // 2
        if d[0] < c[0]: d[1], c[0] = join, join
        else: c[1], d[0] = join, join

    # Avoid drawing duplicate flanking rectangles when they are the common UTR.
    shown = []
    for label in ("distal 3′UTR", "common 3′UTR", f"CE{published_ce}",
                  "upstream exon", "downstream exon"):
        start, end = intervals[label]
        # Treat near-identical annotations as the same displayed block.  Small
        # boundary differences commonly arise after making adjacent UTR blocks
        # touch and should not create a second rectangle or overlapping label.
        duplicate = False
        for a, b, _ in shown:
            overlap = max(0, min(end, b) - max(start, a))
            shorter = max(1, min(end - start, b - a))
            if overlap / shorter >= 0.95:
                duplicate = True
                break
        if duplicate:
            continue
        color = COL_CE if label.startswith("CE") else ("#BEBEBE" if label == "common 3′UTR" else "#D9D9D9")
        short_label = "flanking exon" if label.endswith("exon") else label
        shown.append((start, end, short_label))
        ax.add_patch(Rectangle((start, y - 0.15), end - start, 0.30,
                               facecolor=color, edgecolor="#555555", linewidth=0.5, zorder=3))
        ax.text((start + end) / 2, y - 0.24, short_label,
                ha="center", va="top", fontsize=6.5)

    occupied = union_intervals([(a, b) for a, b, _ in shown])
    ax.plot([min(a for a, _, _ in shown), max(b for _, b, _ in shown)], [y, y],
            color="#666666", lw=1.1, zorder=1)
    # Arrowheads are restricted to intronic line segments.
    for (_, left_end), (right_start, _) in zip(occupied, occupied[1:]):
        if right_start - left_end < 80:
            continue
        spacing = max(320, (view_end - view_start) // 12)
        for x in range(left_end + spacing // 2, right_start, spacing):
            dx = min(55, max(18, (right_start - left_end) // 12))
            p0, p1 = ((x - dx, y), (x + dx, y)) if strand == "+" else ((x + dx, y), (x - dx, y))
            ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=8.5,
                                         linewidth=1.25, color="#666666", zorder=2))
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_xlabel(f"{cfg['chrom']} genomic coordinate (dm6)", fontsize=8)
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.tick_params(axis="x", labelsize=7, length=2)


def render_gene(cfg, stat):
    gene = cfg["gene"]
    folder = DATA / gene
    summary_rows = load_tsv(folder / f"{gene}.pileup_summary.tsv")
    long = fraction_summary(summary_rows, "long")
    short = fraction_summary(summary_rows, "short")
    view_start, view_end = int(cfg["view_start"]), int(cfg["view_end"])

    fig = plt.figure(figsize=(11.2, 8.8 if USE_TRANSCRIPT_MODELS else 7.8))
    bottom_height = 2.15 if USE_TRANSCRIPT_MODELS else 0.8
    gs = fig.add_gridspec(5, 1, height_ratios=[1.0, 2.2, 1.0, 2.2, bottom_height], hspace=0.34)
    axes = [fig.add_subplot(gs[i, 0]) for i in range(5)]
    for ax in axes: ax.set_xlim(view_start, view_end)

    draw_coverage(axes[0], read_bedgraph(folder / f"{gene}.long.mean_fractional_coverage.bedgraph",
                                         cfg["chrom"], view_start, view_end), COL_LONG)
    draw_reads(axes[1], sorted(folder.glob(f"{gene}.long.rep*.sampled.bed")), view_start, view_end)
    draw_coverage(axes[2], read_bedgraph(folder / f"{gene}.short.mean_fractional_coverage.bedgraph",
                                         cfg["chrom"], view_start, view_end), COL_SHORT)
    draw_reads(axes[3], sorted(folder.glob(f"{gene}.short.rep*.sampled.bed")), view_start, view_end)
    if USE_TRANSCRIPT_MODELS:
        draw_transcript_models(axes[4], cfg, stat["CE"])
    else:
        draw_model(axes[4], cfg, stat["CE"])

    for ax, label, color, values in [
        (axes[0], "Long 3′UTR molecules", COL_LONG, long),
        (axes[2], "Short 3′UTR molecules", COL_SHORT, short),
    ]:
        psi = "NA" if not math.isfinite(values["psi"]) else f"{values['psi']:.3f}"
        rep_note = "" if values["n_psi"] == values["n_reps"] else f" ({values['n_psi']}/{values['n_reps']} defined)"
        ax.text(0.006, 1.15, label, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=10, weight="bold", color=color, clip_on=False)
        ax.text(0.994, 1.15,
                f"primary informative reads: {values['informative']:,}  |  mean PSI: {psi}{rep_note}",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
                color="#333333", clip_on=False)

    delta, q = fnum(stat["delta_long_minus_short"]), fnum(stat["BH_q"])
    fig.suptitle(f"{gene} cassette exon {stat['CE']}: long- versus short-3′UTR molecules",
                 x=0.08, y=0.988, ha="left", fontsize=13, weight="bold")
    fig.text(0.08, 0.956,
             f"{cfg['probe_set']} embryo capture; ΔPSI (long − short) = {delta:+.3f}; BH q = {q:.3g}. "
             "Replicates are equally weighted in coverage.",
             ha="left", fontsize=8.5, color="#444444")
    footer = ("Representative transcript models are real Fly_dm6.bed annotations; red marks the cassette exon."
              if USE_TRANSCRIPT_MODELS else
              "Reads are shown in gray; exon inclusion is visible from aligned blocks. Small alignment gaps (≤12 nt) are merged only for display.")
    fig.text(0.08, 0.012, footer,
             ha="left", fontsize=7.5, color="#555555")
    fig.subplots_adjust(left=0.08, right=0.985, top=0.895, bottom=0.075)

    png = OUT / "individual_png" / f"{int(cfg['order']):02d}_{gene}.png"
    pdf = OUT / "individual_pdf" / f"{int(cfg['order']):02d}_{gene}.pdf"
    png.parent.mkdir(parents=True, exist_ok=True); pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=180, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, {
        "order": cfg["order"], "gene": gene, "CE": stat["CE"], "probe_set": cfg["probe_set"],
        "n_replicates": long["n_reps"], "long_informative": long["informative"],
        "short_informative": short["informative"], "mean_long_PSI": long["psi"],
        "mean_short_PSI": short["psi"], "delta_long_minus_short": delta, "BH_q": q,
        "undefined_fraction_replicates": (long["n_reps"] - long["n_psi"]) + (short["n_reps"] - short["n_psi"]),
        "legacy_count_mismatches": long["mismatch"] + short["mismatch"],
        "suggested_tier": "",
    }


def make_contact_sheet(pngs):
    thumb_w = 900
    images = []
    for path in pngs:
        img = Image.open(path).convert("RGB")
        h = round(img.height * thumb_w / img.width)
        images.append((path, img.resize((thumb_w, h), Image.Resampling.LANCZOS)))
    cols = 2
    cell_h = max(i.height for _, i in images) + 52
    rows = math.ceil(len(images) / cols)
    sheet = Image.new("RGB", (thumb_w * cols, cell_h * rows), "white")
    draw = ImageDraw.Draw(sheet)
    for idx, (path, img) in enumerate(images):
        x, y = (idx % cols) * thumb_w, (idx // cols) * cell_h
        draw.text((x + 12, y + 8), path.stem.replace("_", "  ", 1), fill="#222222")
        sheet.paste(img, (x, y + 40))
    out = OUT / "figure3e_read_pileup_catalog_contact_sheet.jpg"
    sheet.save(out, quality=90, optimize=True)
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    configs = load_tsv(CONFIG)
    stats = {r["gene"]: r for r in load_tsv(STATS)}
    pngs, audit = [], []
    for cfg in configs:
        print(f"Rendering {cfg['order']}/28: {cfg['gene']}", flush=True)
        png, row = render_gene(cfg, stats[cfg["gene"]])
        pngs.append(png); audit.append(row)

    for row in audit:
        depth = row["long_informative"] + row["short_informative"]
        q = row["BH_q"]
        effect = abs(row["delta_long_minus_short"])
        if q < 0.05 and effect >= 0.4 and depth >= 1000 and row["undefined_fraction_replicates"] == 0:
            row["suggested_tier"] = "MAIN_CANDIDATE"
        elif depth < 500 or row["undefined_fraction_replicates"] > 0:
            row["suggested_tier"] = "LOW_DEPTH_OR_PARTIAL"
        else:
            row["suggested_tier"] = "SUPPLEMENT"
    audit.sort(key=lambda r: (r["suggested_tier"] != "MAIN_CANDIDATE", -abs(r["delta_long_minus_short"])))
    audit_path = OUT / "figure3e_read_pileup_catalog_audit.tsv"
    with audit_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(audit[0]), delimiter="\t")
        w.writeheader(); w.writerows(audit)
    contact = make_contact_sheet(pngs)
    print(audit_path)
    print(contact)


if __name__ == "__main__":
    main()
