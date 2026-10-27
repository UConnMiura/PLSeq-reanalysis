#!/usr/bin/env python3
"""Render the Dys CE33/alternative-3'UTR locus with pyGenomeTracks."""

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path


TRANSCRIPTS = {
    "FBtr0334562": ("Long 3'UTR | CE33 included", "23,140,140"),
    "FBtr0110917": ("Long 3'UTR | CE33 skipped", "23,140,140"),
    "FBtr0110919": ("Short 3'UTR | CE33 included", "224,122,63"),
    "FBtr0083764": ("Short 3'UTR | CE33 skipped", "224,122,63"),
}


def matching_line(path: Path, names: set[str]) -> list[list[str]]:
    rows = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4 and fields[3] in names:
                rows.append(fields)
    return rows


def feature(path: Path, name: str) -> list[str]:
    rows = matching_line(path, {name})
    if len(rows) != 1:
        raise RuntimeError(f"Expected exactly one {name!r} row in {path}; found {len(rows)}")
    return rows[0]


def write_bed6(path: Path, row: list[str], label: str, rgb: str) -> None:
    fields = row[:6] if len(row) >= 6 else row[:4] + ["0", "."]
    fields[3] = label
    fields[4] = "0"
    path.write_text("\t".join(fields) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--pygenometracks", type=Path, required=True)
    args = parser.parse_args()

    ref = args.reference.resolve()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)

    ce = feature(ref / "SE_PSI.figure3e.suppdata9.bed", "Dys_exon33")
    upstream = feature(ref / "up.figure3e.reconstructed.bed", "Dys_exon33")
    downstream = feature(ref / "down.figure3e.reconstructed.bed", "Dys_exon33")
    common = feature(ref / "UTR.bed", "Dys")
    distal = feature(ref / "dis.bed", "Dys")

    transcript_rows = matching_line(ref / "Fly_dm6.bed", set(TRANSCRIPTS))
    observed = {row[3] for row in transcript_rows}
    missing = set(TRANSCRIPTS) - observed
    if missing:
        raise RuntimeError(f"Missing expected Dys transcripts: {sorted(missing)}")

    transcript_bed = out / "Dys_four_isoform_models.bed"
    by_id = {row[3]: row for row in transcript_rows}
    with transcript_bed.open("w") as handle:
        for transcript_id in TRANSCRIPTS:
            label, rgb = TRANSCRIPTS[transcript_id]
            row = by_id[transcript_id].copy()
            row[3] = f"{label} | {transcript_id}"
            row[8] = rgb
            handle.write("\t".join(row) + "\n")

    write_bed6(out / "Dys_CE33.bed", ce, "CE33", "123,44,191")
    write_bed6(out / "Dys_upstream_filter.bed", upstream, "Upstream exon", "100,116,139")
    write_bed6(out / "Dys_downstream_common_UTR.bed", downstream,
               "Downstream/common 3'UTR", "148,163,184")
    write_bed6(out / "Dys_common_UTR.bed", common, "Common 3'UTR", "148,163,184")
    write_bed6(out / "Dys_distal_UTR.bed", distal, "Long-specific distal 3'UTR", "23,140,140")

    ini = out / "Dys_CE33_3UTR_tracks.ini"
    ini.write_text(f"""[Dys isoforms]
file = {transcript_bed}
file_type = bed
title = Dys transcript isoforms
height = 5.2
display = stacked
labels = true
style = UCSC
color = bed_rgb
border_color = bed_rgb
color_utr = bed_rgb
color_backbone = bed_rgb
fontsize = 9
line_width = 1
arrow_interval = 12
color_arrow = bed_rgb
all_labels_inside = true

[spacer]
height = 0.25

[cassette exon]
file = {out / 'Dys_CE33.bed'}
file_type = bed
title = Cassette exon
height = 0.55
color = #7B2CBF
border_color = #7B2CBF
display = collapsed
labels = true

[upstream exon]
file = {out / 'Dys_upstream_filter.bed'}
file_type = bed
title = Upstream filter
height = 0.55
color = #64748B
border_color = #64748B
display = collapsed
labels = false

[common UTR]
file = {out / 'Dys_common_UTR.bed'}
file_type = bed
title = Common 3'UTR filter
height = 0.55
color = #94A3B8
border_color = #94A3B8
display = collapsed
labels = false

[distal UTR]
file = {out / 'Dys_distal_UTR.bed'}
file_type = bed
title = Long-specific distal 3'UTR
height = 0.55
color = #178C8C
border_color = #178C8C
display = collapsed
labels = false

[x-axis]
where = bottom
""")

    region_start = min(int(upstream[1]), int(ce[1])) - 800
    region_end = max(int(distal[2]), int(common[2])) + 350
    region = f"{ce[0]}:{region_start}-{region_end}"
    png = out / "Dys_CE33_alternative_3UTR.png"
    pdf = out / "Dys_CE33_alternative_3UTR.pdf"
    env = os.environ.copy()
    env["MPLCONFIGDIR"] = "/private/tmp/plseq-pygenometracks-mpl"
    env["PATH"] = f"{args.pygenometracks.parent}:{env.get('PATH', '')}"
    for destination, dpi in ((png, "220"), (pdf, "220")):
        subprocess.run(
            [str(args.pygenometracks), "--tracks", str(ini), "--region", region,
             "--width", "30", "--height", "12", "--dpi", dpi,
             "--trackLabelFraction", "0.18", "--outFileName", str(destination)],
            check=True, env=env,
        )

    (out / "README.txt").write_text(
        "Dys CE33 and alternative 3'UTR pyGenomeTracks pilot\n"
        "===================================================\n\n"
        f"Region: {region}\n"
        "Transcript models are real Fly_dm6.bed FlyBase BED12 records.\n"
        "Teal transcripts use the long 3'UTR; orange transcripts use the short 3'UTR.\n"
        "Purple marks the Supplementary Data 9 cassette-exon interval.\n"
        "Gray and teal interval tracks show the filtering regions used for isoform assignment.\n"
        "The selected four transcripts represent CE inclusion/skipping crossed with long/short 3'UTRs.\n"
    )
    print(f"Rendered {png}")


if __name__ == "__main__":
    main()
