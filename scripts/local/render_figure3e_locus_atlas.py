#!/usr/bin/env python3
"""Render a 28-gene Figure 3e locus atlas using pyGenomeTracks."""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


LONG_RGB = "23,140,140"
SHORT_RGB = "224,122,63"
CATEGORY_ORDER = (
    "Long | CE included",
    "Long | CE skipped",
    "Short | CE included",
    "Short | CE skipped",
)
PUBLISHED_CE = (33, 15, 9, 8, 18, 17, 11, 12, 9, 10, 15, 9, 13, 2,
                10, 4, 12, 15, 14, 12, 6, 6, 10, 15, 4, 4, 16, 23)


def read_bed(path: Path) -> list[list[str]]:
    rows = []
    with path.open() as handle:
        for line in handle:
            if line.strip() and not line.startswith(("#", "track", "browser")):
                # Some supplied BEDs mix tabs and runs of spaces.
                rows.append(line.split())
    return rows


def blocks(row: list[str]) -> list[tuple[int, int]]:
    start = int(row[1])
    sizes = [int(x) for x in row[10].rstrip(",").split(",")]
    offsets = [int(x) for x in row[11].rstrip(",").split(",")]
    return [(start + offset, start + offset + size)
            for size, offset in zip(sizes, offsets)]


def overlaps(intervals: list[tuple[int, int]], start: int, end: int) -> bool:
    return any(a < end and b > start for a, b in intervals)


def font(size: int, bold: bool = False):
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Helvetica.ttc",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size)
        except OSError:
            pass
    return ImageFont.load_default()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--pygenometracks", type=Path, required=True)
    args = ap.parse_args()

    ref = args.reference.resolve()
    out = args.output.resolve()
    panels = out / "panels"
    panels.mkdir(parents=True, exist_ok=True)

    ce_rows = read_bed(ref / "SE_PSI.figure3e.suppdata9.bed")
    up = {row[3]: row for row in read_bed(ref / "up.figure3e.reconstructed.bed")}
    down = {row[3]: row for row in read_bed(ref / "down.figure3e.reconstructed.bed")}
    common = {row[3]: row for row in read_bed(ref / "UTR.bed")}
    distal = {row[3]: row for row in read_bed(ref / "dis.bed")}

    by_chrom: dict[str, list[list[str]]] = defaultdict(list)
    for row in read_bed(ref / "Fly_dm6.bed"):
        if len(row) >= 12:
            by_chrom[row[0]].append(row)

    manifest = []
    panel_paths = []
    panel_meta = []
    env = os.environ.copy()
    env["MPLCONFIGDIR"] = "/private/tmp/plseq-atlas-mpl"
    env["PATH"] = f"{args.pygenometracks.parent}:{env.get('PATH', '')}"

    for order, ce in enumerate(ce_rows, 1):
        target = ce[3]
        gene = target.rsplit("_exon", 1)[0]
        bed_exon_label = target.rsplit("_exon", 1)[1]
        ce_number = str(PUBLISHED_CE[order - 1])
        if target not in up or target not in down or gene not in common or gene not in distal:
            raise RuntimeError(f"Missing feature interval for {target}")

        chrom = ce[0]
        ce_start, ce_end = int(ce[1]), int(ce[2])
        com_start, com_end = int(common[gene][1]), int(common[gene][2])
        dis_start, dis_end = int(distal[gene][1]), int(distal[gene][2])

        candidates = []
        for tx in by_chrom[chrom]:
            tx_blocks = blocks(tx)
            if not overlaps(tx_blocks, com_start, com_end):
                continue
            is_long = overlaps(tx_blocks, dis_start, dis_end)
            includes_ce = overlaps(tx_blocks, ce_start, ce_end)
            category = ("Long" if is_long else "Short") + " | CE " + (
                "included" if includes_ce else "skipped")
            candidates.append((category, tx))

        grouped: dict[str, list[list[str]]] = defaultdict(list)
        for category, tx in candidates:
            grouped[category].append(tx)

        # Choose the model with the smallest full transcript span in each category;
        # this generally favors a compact representative while preserving the locus.
        selected = {}
        for category in CATEGORY_ORDER:
            rows = grouped.get(category, [])
            if rows:
                selected[category] = min(rows, key=lambda x: (int(x[2]) - int(x[1]), x[3]))

        gene_dir = panels / f"{order:02d}_{gene}"
        gene_dir.mkdir(exist_ok=True)
        tx_bed = gene_dir / "transcripts.bed"
        with tx_bed.open("w") as handle:
            for category in CATEGORY_ORDER:
                if category not in selected:
                    continue
                tx = selected[category].copy()
                transcript_id = tx[3]
                tx[3] = category
                tx[8] = LONG_RGB if category.startswith("Long") else SHORT_RGB
                handle.write("\t".join(tx) + "\n")
                manifest.append({
                    "order": order, "gene": gene, "published_CE": ce_number,
                    "bed_target": target, "bed_exon_label": bed_exon_label,
                    "strand": tx[5], "category": category,
                    "transcript": transcript_id,
                    "available_models_in_category": len(grouped[category]),
                    "status": "SELECTED",
                })
            for category in CATEGORY_ORDER:
                if category not in selected:
                    manifest.append({
                        "order": order, "gene": gene, "published_CE": ce_number,
                        "bed_target": target, "bed_exon_label": bed_exon_label,
                        "strand": candidates[0][1][5] if candidates else "NA",
                        "category": category, "transcript": "NA",
                        "available_models_in_category": 0,
                        "status": "NO_ANNOTATED_MODEL",
                    })

        def write_feature(filename: str, row: list[str], label: str) -> Path:
            path = gene_dir / filename
            fields = row[:6] if len(row) >= 6 else row[:4] + ["0", "."]
            fields[3], fields[4] = label, "0"
            path.write_text("\t".join(fields) + "\n")
            return path

        ce_bed = write_feature("ce.bed", ce, f"CE{ce_number}")
        com_bed = write_feature("common.bed", common[gene], "Common")
        dis_bed = write_feature("distal.bed", distal[gene], "Distal")

        ini = gene_dir / "tracks.ini"
        strand = candidates[0][1][5] if candidates else "."
        ini.write_text(f"""[isoforms]
file = {tx_bed}
file_type = bed
title =
height = 3.4
display = stacked
labels = true
style = UCSC
color = bed_rgb
border_color = bed_rgb
color_utr = bed_rgb
color_backbone = bed_rgb
fontsize = 8
line_width = 1
arrow_interval = 10
color_arrow = bed_rgb
all_labels_inside = true

[cassette exon]
file = {ce_bed}
file_type = bed
title = CE
height = 0.35
color = #DDE2E8
border_color = #B8C0CC
display = collapsed
labels = false

[common UTR]
file = {com_bed}
file_type = bed
title = Common
height = 0.35
color = #DDE2E8
border_color = #B8C0CC
display = collapsed
labels = false

[distal UTR]
file = {dis_bed}
file_type = bed
title = Distal
height = 0.35
color = #DDE2E8
border_color = #B8C0CC
display = collapsed
labels = false

""")

        all_intervals = [ce, up[target], down[target], common[gene], distal[gene]]
        start = min(int(row[1]) for row in all_intervals)
        end = max(int(row[2]) for row in all_intervals)
        padding = max(100, round((end - start) * 0.04))
        region = f"{chrom}:{max(0, start-padding)}-{end+padding}"
        panel = gene_dir / "panel.png"
        subprocess.run(
            [str(args.pygenometracks), "--tracks", str(ini), "--region", region,
             "--width", "14", "--height", "6.2", "--dpi", "130",
             "--trackLabelFraction", "0.18", "--outFileName", str(panel)],
            check=True, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        panel_paths.append(panel)
        panel_meta.append((gene, ce_number, chrom, strand, (end - start) / 1000))
        print(f"{order:02d}/28 {gene}: {len(selected)} combinations")

    with (out / "figure3e_locus_atlas_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)

    opened = [Image.open(path).convert("RGB") for path in panel_paths]
    cell_w = max(image.width for image in opened)
    raw_h = max(image.height for image in opened)
    cell_title_h = 48
    cell_h = raw_h + cell_title_h
    header_h = 150
    footer_h = 80
    canvas = Image.new("RGB", (cell_w * 4, header_h + cell_h * 7 + footer_h), "white")
    draw = ImageDraw.Draw(canvas)
    draw.text((canvas.width // 2, 26),
              "Figure 3e loci: cassette exons and alternative 3'UTRs",
              fill="#263238", font=font(48, bold=True), anchor="ma")
    draw.text((canvas.width // 2, 88),
              "Teal: long 3'UTR   Orange: short 3'UTR   Light gray: CE, common and distal reference intervals",
              fill="#52606D", font=font(28), anchor="ma")
    cell_title_font = font(24, bold=True)
    cell_subtitle_font = font(18)
    for index, image in enumerate(opened):
        row, col = divmod(index, 4)
        cell_left = col * cell_w
        cell_top = header_h + row * cell_h
        draw.rectangle(
            (cell_left + 4, cell_top + 2, cell_left + cell_w - 5,
             cell_top + cell_h - 4),
            outline="#D7DCE2", width=2,
        )
        x = col * cell_w + (cell_w - image.width) // 2
        y0 = header_h + row * cell_h
        y = y0 + cell_title_h
        gene, ce_number, chrom, strand, span_kb = panel_meta[index]
        draw.text((col * cell_w + 16, y0 + 5), f"{gene}  CE{ce_number}",
                  fill="#263238", font=cell_title_font)
        draw.text((col * cell_w + cell_w - 16, y0 + 9),
                  f"{chrom}  {strand} strand  |  {span_kb:.1f} kb window",
                  fill="#667085", font=cell_subtitle_font, anchor="ra")
        canvas.paste(image, (x, y))
    draw.text((canvas.width // 2, canvas.height - 35),
              "Each panel has an independent genomic scale; arrows indicate transcriptional direction.",
              fill="#52606D", font=font(25), anchor="mm")
    canvas.save(out / "figure3e_28_gene_locus_atlas.png", dpi=(200, 200))
    canvas.save(out / "figure3e_28_gene_locus_atlas.pdf", resolution=200)

    with (out / "README.txt").open("w") as handle:
        handle.write(
            "Figure 3e 28-gene locus atlas\n"
            "================================\n\n"
            "Each mini-panel was rendered independently by pyGenomeTracks.\n"
            "One real FlyBase BED12 model was selected for each available long/short 3'UTR\n"
            "and cassette-exon included/skipped combination. No missing combination was invented.\n"
            "Panels retain genomic left-to-right orientation; transcript arrows show strand.\n"
            "Because every locus uses an independent coordinate range, apparent feature widths\n"
            "must be interpreted only within a panel. See the manifest for transcript provenance.\n"
        )
    print(f"Wrote atlas: {out / 'figure3e_28_gene_locus_atlas.png'}")


if __name__ == "__main__":
    main()
