#!/usr/bin/env python3
"""Prepare replicate-aware PL-Seq read pileups from filtered BAMs.

The script preserves CIGAR-derived alignment blocks, classifies cassette-exon
usage using the thresholds in the published ``calculate_PSI.py`` workflow,
checks the resulting informative-read counts against its PSI CSV outputs, and
writes compact BED12 and fractional-coverage files for local rendering.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from collections import defaultdict
from pathlib import Path

import pysam


def overlap_len(blocks, start, end):
    return sum(max(0, min(b, end) - max(a, start)) for a, b in blocks)


def classify(blocks, ce_start, ce_end):
    """Faithful operational equivalent of calculate_PSI.py's thresholds."""
    left = overlap_len(blocks, ce_start - 50, ce_start) / 50.0
    right = overlap_len(blocks, ce_end, ce_end + 50) / 50.0
    if left >= 0.1 or right >= 0.1:
        return "boundary_fail"
    rate = overlap_len(blocks, ce_start, ce_end) / float(ce_end - ce_start)
    if rate > 0.8:
        return "included"
    if rate < 0.1:
        return "skipped"
    return "unresolved"


def bed12(read, name, color):
    blocks = read.get_blocks()
    chrom_start = blocks[0][0]
    chrom_end = blocks[-1][1]
    sizes = ",".join(str(e - s) for s, e in blocks) + ","
    starts = ",".join(str(s - chrom_start) for s, _ in blocks) + ","
    strand = "-" if read.is_reverse else "+"
    return [
        read.reference_name, chrom_start, chrom_end, name, 0, strand,
        chrom_start, chrom_end, color, len(blocks), sizes, starts,
    ]


def stable_rank(read_name, seed):
    return hashlib.sha256(f"{seed}|{read_name}".encode()).hexdigest()


def allocate_sample(records, limit, seed):
    """Proportional deterministic sample, retaining a real minority class."""
    if len(records) <= limit:
        return list(records)
    by_status = defaultdict(list)
    for rec in records:
        by_status[rec[0]].append(rec)
    useful = {k: v for k, v in by_status.items() if k in {"included", "skipped"}}
    chosen = []
    remaining = limit
    allocations = {}
    total = sum(len(v) for v in useful.values())
    for status, vals in useful.items():
        allocations[status] = max(1, round(limit * len(vals) / total))
    while sum(allocations.values()) > limit:
        key = max(allocations, key=lambda k: allocations[k])
        if allocations[key] > 1:
            allocations[key] -= 1
        else:
            break
    for status, n in allocations.items():
        vals = sorted(useful[status], key=lambda r: stable_rank(r[1].query_name, seed))
        chosen.extend(vals[:n])
        remaining -= min(n, len(vals))
    if remaining > 0:
        leftovers = [r for r in records if r not in chosen]
        leftovers.sort(key=lambda r: stable_rank(r[1].query_name, seed))
        chosen.extend(leftovers[:remaining])
    return chosen[:limit]


def read_expected(path, event):
    with path.open() as handle:
        for row in csv.reader(handle):
            if row and row[0] == event:
                # calculate_PSI.py writes an empty fourth field when the
                # denominator is zero.  That is a legitimate undefined PSI,
                # not a malformed record; retain the 0/0 counts and represent
                # the PSI as missing in the audit table.
                psi = float(row[3]) if row[3].strip() else None
                return int(row[1]), int(row[2]), psi
    raise ValueError(f"{event} absent from {path}")


def write_bedgraph(path, chrom, start, values):
    with path.open("w") as out:
        if not values:
            return
        run_start = 0
        current = values[0]
        for i in range(1, len(values) + 1):
            changed = i == len(values) or not math.isclose(values[i], current, abs_tol=1e-12)
            if changed:
                out.write(f"{chrom}\t{start + run_start}\t{start + i}\t{current:.8f}\n")
                if i < len(values):
                    run_start, current = i, values[i]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--bam-root", type=Path, required=True)
    ap.add_argument("--psi-root", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--event", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--region-start", type=int, required=True)
    ap.add_argument("--region-end", type=int, required=True)
    ap.add_argument("--ce-start", type=int, required=True)
    ap.add_argument("--ce-end", type=int, required=True)
    ap.add_argument("--reads-per-replicate", type=int, default=20)
    ap.add_argument("--seed", type=int, default=37679364)
    ap.add_argument(
        "--allow-count-mismatch", action="store_true",
        help="Write outputs and return success when the primary-only classifier differs slightly from legacy PSI counts."
    )
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    with args.manifest.open() as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))

    summary = []
    coverage_by_fraction = defaultdict(list)
    colors = {"included": "228,87,86", "skipped": "115,128,140",
              "unresolved": "190,190,190", "boundary_fail": "215,215,215"}

    for row in manifest:
        run, replicate = row["run_accession"], row["replicate"]
        for fraction in ("long", "short"):
            # Prefer a compact gene-specific primary BAM when one exists.  For
            # catalog generation, fall back to the indexed full fraction BAM;
            # the regional fetch and flag checks below produce the same input
            # without creating 224 redundant BAM files.
            bam_path = args.bam_root / run / f"{args.gene}.{fraction}.primary.bam"
            if not bam_path.exists():
                bam_path = args.bam_root / run / f"{fraction}.bam"
            if not bam_path.exists():
                raise FileNotFoundError(f"No BAM found for {run} {args.gene} {fraction}")
            records = []
            coverage = [0] * (args.region_end - args.region_start)
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                for read in bam.fetch(args.chrom, args.region_start, args.region_end):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    blocks = read.get_blocks()
                    status = classify(blocks, args.ce_start, args.ce_end)
                    records.append((status, read))
                    for start, end in blocks:
                        lo = max(start, args.region_start) - args.region_start
                        hi = min(end, args.region_end) - args.region_start
                        for i in range(max(0, lo), min(len(coverage), hi)):
                            coverage[i] += 1
            denom = len(records)
            coverage_by_fraction[fraction].append(
                [v / denom if denom else 0.0 for v in coverage])

            counts = defaultdict(int)
            for status, _ in records:
                counts[status] += 1
            expected_i, expected_n, expected_psi = read_expected(
                args.psi_root / run / f"{fraction}.psi.csv", args.event)
            observed_i = counts["included"]
            observed_n = counts["included"] + counts["skipped"]
            summary.append({
                "run": run, "replicate": replicate, "gene": args.gene,
                "fraction": fraction, "primary": len(records),
                "included": observed_i, "skipped": counts["skipped"],
                "unresolved": counts["unresolved"],
                "boundary_fail": counts["boundary_fail"],
                "informative": observed_n,
                "PSI": observed_i / observed_n if observed_n else "",
                "expected_included": expected_i, "expected_informative": expected_n,
                "expected_PSI": expected_psi,
                "count_match": observed_i == expected_i and observed_n == expected_n,
            })

            sample = allocate_sample(records, args.reads_per_replicate,
                                     args.seed + int(replicate))
            bed_path = args.output / f"{args.gene}.{fraction}.rep{replicate}.sampled.bed"
            with bed_path.open("w") as out:
                for status, read in sorted(sample, key=lambda r: (r[0], r[1].reference_start)):
                    name = f"{run}|rep{replicate}|{status}|{read.query_name}"
                    out.write("\t".join(map(str, bed12(read, name, colors[status]))) + "\n")

    for fraction, arrays in coverage_by_fraction.items():
        mean = [sum(vals) / len(vals) for vals in zip(*arrays)]
        write_bedgraph(args.output / f"{args.gene}.{fraction}.mean_fractional_coverage.bedgraph",
                       args.chrom, args.region_start, mean)

    summary_path = args.output / f"{args.gene}.pileup_summary.tsv"
    with summary_path.open("w", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=list(summary[0]))
        writer.writeheader(); writer.writerows(summary)

    failures = [r for r in summary if not r["count_match"]]
    print(f"Wrote {args.gene} pileup inputs to {args.output}")
    if failures:
        print(f"WARNING: {len(failures)} replicate/fraction count checks differ; inspect summary")
        if not args.allow_count_mismatch:
            raise SystemExit(2)
        print("Continuing because --allow-count-mismatch was specified")
        return
    print("All included/informative counts match calculate_PSI.py outputs")


if __name__ == "__main__":
    main()
