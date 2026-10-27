#!/usr/bin/env python3
"""Compose four quantitative annotation variants of the Figure 3e locus atlas."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


TEAL = "#178C8C"
ORANGE = "#E07A3F"
DARK = "#263238"
MID = "#667085"
LIGHT = "#D7DCE2"


def get_font(size: int, bold: bool = False):
    paths = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Helvetica.ttc",
    ]
    for path in paths:
        try:
            return ImageFont.truetype(path, size)
        except OSError:
            pass
    return ImageFont.load_default()


def q_text(q: float) -> str:
    return f"q={q:.1e}" if q < 0.001 else f"q={q:.3f}"


def read_inputs(atlas_dir: Path, stats_path: Path):
    with stats_path.open() as handle:
        stats = {row["gene"]: row for row in csv.DictReader(handle, delimiter="\t")}
    with (atlas_dir / "figure3e_locus_atlas_manifest.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    metadata = {}
    for row in rows:
        order = int(row["order"])
        metadata.setdefault(order, {
            "gene": row["gene"], "CE": row["published_CE"], "strand": row["strand"]
        })
    panels = []
    for order in sorted(metadata):
        gene = metadata[order]["gene"]
        path = next((atlas_dir / "panels").glob(f"{order:02d}_{gene}/panel.png"))
        panels.append((metadata[order], stats[gene], Image.open(path).convert("RGB")))
    return panels


def title(draw, x0, y0, width, meta, stat, style):
    draw.text((x0 + 14, y0 + 8), f"{meta['gene']}  CE{meta['CE']}",
              font=get_font(23, True), fill=DARK)
    if style == "numeric":
        delta = float(stat["delta_long_minus_short"])
        q = float(stat["BH_q"])
        draw.text((x0 + width - 14, y0 + 10),
                  f"ΔPSI={delta:+.3f}   {q_text(q)}",
                  font=get_font(19, True), fill=TEAL if delta >= 0 else ORANGE,
                  anchor="ra")
    else:
        draw.text((x0 + width - 14, y0 + 11), f"{meta['strand']} strand",
                  font=get_font(17), fill=MID, anchor="ra")


def dumbbell(draw, box, stat):
    x0, y0, x1, y1 = box
    left, right = x0 + 48, x1 - 155
    cy = (y0 + y1) // 2
    long = float(stat["mean_long"]); short = float(stat["mean_short"])
    xl = left + long * (right - left); xs = left + short * (right - left)
    draw.line((left, cy, right, cy), fill="#C7CDD4", width=2)
    draw.line((xs, cy, xl, cy), fill="#9CA6B0", width=4)
    draw.ellipse((xl-8, cy-8, xl+8, cy+8), fill=TEAL)
    draw.polygon([(xs, cy-10), (xs-9, cy+8), (xs+9, cy+8)], fill=ORANGE)
    draw.text((left, cy+13), "0", font=get_font(14), fill=MID, anchor="ma")
    draw.text((right, cy+13), "1", font=get_font(14), fill=MID, anchor="ma")
    draw.text((x1-12, cy-2), q_text(float(stat["BH_q"])),
              font=get_font(18, True), fill=DARK, anchor="ra")


def delta_bar(draw, box, stat):
    x0, y0, x1, y1 = box
    left, right = x0 + 48, x1 - 155
    cy = (y0 + y1) // 2
    center = (left + right) / 2
    delta = float(stat["delta_long_minus_short"])
    endpoint = center + delta * (right-left)/2
    draw.line((left, cy, right, cy), fill="#C7CDD4", width=2)
    draw.line((center, cy-10, center, cy+10), fill="#7B8794", width=2)
    draw.rectangle((min(center, endpoint), cy-7, max(center, endpoint), cy+7),
                   fill=TEAL if delta >= 0 else ORANGE)
    draw.text((left, cy+13), "−1", font=get_font(14), fill=MID, anchor="ma")
    draw.text((right, cy+13), "+1", font=get_font(14), fill=MID, anchor="ma")
    draw.text((x1-12, cy-4), f"Δ={delta:+.2f}  {q_text(float(stat['BH_q']))}",
              font=get_font(17, True), fill=DARK, anchor="ra")


def scatter(draw, box, stat):
    x0, y0, x1, y1 = box
    size = min(92, y1-y0-8)
    left = x0 + 52; top = y0 + 4; bottom = top + size; right = left + size
    draw.rectangle((left, top, right, bottom), outline="#AEB6BF", width=2)
    draw.line((left, bottom, right, top), fill="#B9C0C8", width=2)
    short = float(stat["mean_short"]); long = float(stat["mean_long"])
    px = left + short * size; py = bottom - long * size
    draw.ellipse((px-7, py-7, px+7, py+7), fill=TEAL, outline="white", width=2)
    draw.text((left-7, bottom), "L", font=get_font(15, True), fill=TEAL, anchor="rs")
    draw.text((right, bottom+3), "S", font=get_font(15, True), fill=ORANGE, anchor="ma")
    delta = float(stat["delta_long_minus_short"])
    draw.text((left+size+24, top+18), f"ΔPSI {delta:+.2f}",
              font=get_font(18, True), fill=TEAL if delta >= 0 else ORANGE)
    draw.text((left+size+24, top+48), q_text(float(stat["BH_q"])),
              font=get_font(18), fill=DARK)


def compose(panels, destination: Path, style: str):
    panel_w = max(im.width for _, _, im in panels)
    raw_h = max(im.height for _, _, im in panels)
    title_h = 50
    metric_h = 0 if style == "numeric" else (68 if style != "scatter" else 108)
    cell_h = title_h + raw_h + metric_h
    header_h, footer_h = 145, 60
    canvas = Image.new("RGB", (panel_w*4, header_h+cell_h*7+footer_h), "white")
    draw = ImageDraw.Draw(canvas)
    names = {
        "numeric": "Numeric effect and significance",
        "dumbbell": "Long–short PSI dumbbells",
        "delta": "Signed ΔPSI bars",
        "scatter": "Mean long-versus-short PSI insets",
    }
    draw.text((canvas.width//2, 24), "Figure 3e locus atlas",
              font=get_font(46, True), fill=DARK, anchor="ma")
    draw.text((canvas.width//2, 83), names[style],
              font=get_font(29), fill=MID, anchor="ma")
    for i, (meta, stat, im) in enumerate(panels):
        row, col = divmod(i, 4)
        x0, y0 = col*panel_w, header_h+row*cell_h
        draw.rectangle((x0+4, y0+2, x0+panel_w-5, y0+cell_h-4),
                       outline=LIGHT, width=2)
        title(draw, x0, y0, panel_w, meta, stat, style)
        canvas.paste(im, (x0+(panel_w-im.width)//2, y0+title_h))
        metric_box = (x0+8, y0+title_h+raw_h, x0+panel_w-8, y0+cell_h-6)
        if style == "dumbbell": dumbbell(draw, metric_box, stat)
        elif style == "delta": delta_bar(draw, metric_box, stat)
        elif style == "scatter": scatter(draw, metric_box, stat)
    draw.text((canvas.width//2, canvas.height-28),
              "ΔPSI = mean(long 3'UTR PSI) − mean(short 3'UTR PSI); q values are BH-adjusted.",
              font=get_font(22), fill=MID, anchor="mm")
    canvas.save(destination, dpi=(200, 200))
    canvas.save(destination.with_suffix(".pdf"), resolution=200)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas-dir", type=Path, required=True)
    ap.add_argument("--statistics", type=Path, required=True)
    args = ap.parse_args()
    panels = read_inputs(args.atlas_dir, args.statistics)
    out = args.atlas_dir / "metric_variants"
    out.mkdir(exist_ok=True)
    for style in ("numeric", "dumbbell", "delta", "scatter"):
        compose(panels, out / f"figure3e_atlas_{style}.png", style)
    print(f"Wrote four variants to {out}")


if __name__ == "__main__":
    main()
