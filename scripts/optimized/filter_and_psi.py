#!/usr/bin/env python3
"""Faster, operationally equivalent PL-Seq filtering and PSI reconstruction.

The published scripts remain the reference implementation.  This program keeps
their filtering order and threshold arithmetic while removing repeated BAM
fetches, temporary CSV files, quadratic list deletion, and unused failure BAMs.
It requires pysam and is intended to run in the archived PL-Seq Apptainer image.
"""

import argparse
import csv
import hashlib
import os
import re
import shutil
import sys
import time
from collections import OrderedDict

import pysam


def read_bed(path):
    rows = []
    with open(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 4:
                raise ValueError("{}:{} is not BED4".format(path, line_number))
            rows.append((fields[0], int(fields[1]), int(fields[2]), fields[3]))
    if not rows:
        raise ValueError("BED file is empty: {}".format(path))
    return rows


def legacy_covered_bases(blocks, start, end):
    """Reproduce exon_covering.py's block-index and boundary arithmetic."""
    if not blocks:
        return 0

    if blocks[-1][0] < start:
        pre_index = len(blocks) - 1
        pre_delta = 0
    else:
        pre_index = None
        pre_delta = 0
        for index, block in enumerate(blocks):
            if block[0] >= start:
                if index > 0 and blocks[index - 1][1] > start:
                    pre_delta = blocks[index - 1][1] - start
                pre_index = index
                break

    if blocks[-1][1] < end:
        post_index = len(blocks)
        post_delta = 0
    else:
        post_index = None
        post_delta = 0
        for index, block in enumerate(blocks):
            if block[1] >= end:
                if block[0] < end:
                    post_delta = end - block[0]
                post_index = index
                break

    if pre_index is None or post_index is None:
        return 0
    middle = sum(block_end - block_start
                 for block_start, block_end in blocks[pre_index:post_index])
    return pre_delta + post_delta + middle


def sort_and_index(unsorted_bam, output_bam, threads):
    pysam.sort("-@", str(threads), "-o", output_bam, unsorted_bam)
    pysam.index("-@", str(threads), output_bam)
    os.remove(unsorted_bam)


def filter_exon_coverage(input_bam, bed_path, output_bam, work_dir,
                         label, threshold, threads):
    """Keep each BED-overlapping record whose aligned coverage meets threshold."""
    unsorted_bam = os.path.join(work_dir, label + ".unsorted.bam")
    counts = OrderedDict()
    with pysam.AlignmentFile(input_bam, "rb") as source:
        with pysam.AlignmentFile(unsorted_bam, "wb", template=source) as target:
            for chrom, start, end, name in read_bed(bed_path):
                total = 0
                passed = 0
                for read in source.fetch(chrom, start, end):
                    total += 1
                    covered = legacy_covered_bases(read.get_blocks(), start, end)
                    if float(covered) / float(end - start) >= threshold:
                        target.write(read)
                        passed += 1
                counts[name] = (passed, total)
                print("{} valid lines:  {}".format(name, passed))
                print("{} total lines:  {}".format(name, total))
    sort_and_index(unsorted_bam, output_bam, threads)
    return counts


def interval_lookup(rows):
    lookup = {}
    for chrom, start, end, unused_name in rows:
        lookup.setdefault(chrom, []).append((start, end))
    for chrom in lookup:
        lookup[chrom].sort()
    return lookup


def overlaps_any_reference_span(read, lookup):
    """Match the reference-span overlap used by `samtools view -L`."""
    if read.is_unmapped or read.reference_end is None:
        return False
    chrom = read.reference_name
    for start, end in lookup.get(chrom, []):
        if start >= read.reference_end:
            break
        if read.reference_start < end and read.reference_end > start:
            return True
    return False


def split_distal_overlap(input_bam, distal_bed, overlap_bam,
                         short_candidate_bam, threads):
    lookup = interval_lookup(read_bed(distal_bed))
    overlap_count = 0
    candidate_count = 0
    with pysam.AlignmentFile(input_bam, "rb") as source:
        with pysam.AlignmentFile(overlap_bam, "wb", template=source) as overlap_out:
            with pysam.AlignmentFile(short_candidate_bam, "wb", template=source) as short_out:
                for read in source.fetch(until_eof=True):
                    if overlaps_any_reference_span(read, lookup):
                        overlap_out.write(read)
                        overlap_count += 1
                    else:
                        short_out.write(read)
                        candidate_count += 1
    pysam.index("-@", str(threads), short_candidate_bam)
    return overlap_count, candidate_count


def cigar_terminal_lengths(cigar_string):
    values = re.findall(r"\d+", cigar_string or "")
    if not values:
        raise ValueError("alignment has no parseable CIGAR: {!r}".format(cigar_string))
    return int(values[0]), int(values[-1])


def passes_legacy_polya(read):
    """Reproduce active second_round_filtering() in polyA_filtering.py."""
    sequence = read.query_sequence
    if sequence is None:
        return False
    first_length, last_length = cigar_terminal_lengths(read.cigarstring)
    pre = sequence[first_length - 10:first_length + 8]
    post = sequence[-last_length - 8:-last_length + 10]
    pre_condition = pre[:10].count("T") >= 8 and pre[-8:].count("T") < 8
    post_condition = post[-10:].count("A") >= 8 and post[:8].count("A") < 8
    return ((pre_condition or post_condition) and
            ((pre[:10].count("T") >= 8) != (post[-10:].count("A") >= 8)))


def filter_polya(input_bam, output_bam, threads):
    total = 0
    passed = 0
    with pysam.AlignmentFile(input_bam, "rb") as source:
        with pysam.AlignmentFile(output_bam, "wb", template=source) as target:
            for read in source.fetch(until_eof=True):
                total += 1
                if passes_legacy_polya(read):
                    target.write(read)
                    passed += 1
    pysam.index("-@", str(threads), output_bam)
    print("Total lines: {}, valid lines: {}".format(total, passed))
    return total, passed


def legacy_overlap_indices(blocks, start, end):
    """Return indices/deltas calculated by calculate_PSI.py for one read."""
    if not blocks:
        return 0, 0, 0, 0

    if blocks[-1][0] < start:
        pre_index = len(blocks) - 1
        pre_delta = 0
    else:
        pre_index = None
        pre_delta = 0
        for index, block in enumerate(blocks):
            if block[0] >= start:
                if index > 0 and blocks[index - 1][1] > start:
                    pre_delta = blocks[index - 1][1] - start
                pre_index = index
                break

    if blocks[-1][1] < end:
        post_index = len(blocks)
        post_delta = 0
    else:
        post_index = None
        post_delta = 0
        for index, block in enumerate(blocks):
            if block[1] >= end:
                if block[0] < end:
                    post_delta = end - block[0]
                post_index = index
                break

    if pre_index is None or post_index is None:
        return 0, 0, 0, 0
    return pre_index, post_index, pre_delta, post_delta


def legacy_psi_components(block_lists, start, end):
    """Reproduce calculate_PSI.py, including its shared-list behavior.

    In the archived program, ``calculate_difference_of_block`` appends the same
    growing ``diff_row`` object for every read and never clears it.  Every row
    therefore refers to the final flattened block-length list.  This is not the
    biologically natural calculation, but preserving it is necessary for exact
    reconstruction of the deposited PSI values.  A corrected metric should be
    exposed separately rather than silently changing the legacy result.
    """
    indices = [legacy_overlap_indices(blocks, start, end)
               for blocks in block_lists]
    cumulative_lengths = []
    for blocks in block_lists:
        cumulative_lengths.extend(block_end - block_start
                                  for block_start, block_end in blocks)
    components = []
    for pre_index, post_index, pre_delta, post_delta in indices:
        middle = sum(cumulative_lengths[pre_index:post_index])
        components.append((pre_delta, post_delta, middle))
    return components


def psi_for_event(source, chrom, start, end, threshold):
    fetched_blocks = []
    for read in source.fetch(chrom, start - 50, end + 50):
        fetched_blocks.append(read.get_blocks())

    left_components = legacy_psi_components(
        fetched_blocks, start - 50, start - 1)
    right_components = legacy_psi_components(
        fetched_blocks, end + 1, end + 50)
    blocks_to_test = []
    for blocks, left, right in zip(
            fetched_blocks, left_components, right_components):
        # Preserve the original scripts' operator precedence exactly.
        left_rate = left[0] + left[1] + left[2] / 50.0
        right_rate = right[0] + right[1] + right[2] / 50.0
        if left_rate < 0.1 and right_rate < 0.1:
            blocks_to_test.append(blocks)

    included = 0
    skipped = 0
    event_components = legacy_psi_components(blocks_to_test, start, end)
    for component in event_components:
        rate = component[0] + component[1] + component[2] / float(end - start)
        if rate > threshold:
            included += 1
        elif rate < 0.1:
            skipped += 1
    psi = float(included) / float(included + skipped) if included + skipped else None
    return included, skipped, psi


def calculate_psi(input_bam, cassette_bed, output_csv, threshold):
    with pysam.AlignmentFile(input_bam, "rb") as source:
        with open(output_csv, "w", newline="") as handle:
            writer = csv.writer(handle)
            for chrom, start, end, name in read_bed(cassette_bed):
                included, skipped, psi = psi_for_event(
                    source, chrom, start, end, threshold)
                writer.writerow([name, included, included + skipped, psi])


def count_records(path):
    count = 0
    with pysam.AlignmentFile(path, "rb") as source:
        for unused_read in source.fetch(until_eof=True):
            count += 1
    return count


def write_sha256(paths, output_path):
    with open(output_path, "w") as output:
        for path in paths:
            digest = hashlib.sha256()
            with open(path, "rb") as handle:
                for block in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(block)
            output.write("{}  {}\n".format(digest.hexdigest(), path))


def require_file(path):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise FileNotFoundError("missing or empty required file: {}".format(path))


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-bam", required=True)
    parser.add_argument("--upstream-bed", required=True)
    parser.add_argument("--downstream-bed", required=True)
    parser.add_argument("--universal-utr-bed", required=True)
    parser.add_argument("--extended-utr-bed", required=True)
    parser.add_argument("--cassette-exon-bed", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--work-dir")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--coverage-threshold", type=float, default=0.8)
    parser.add_argument("--keep-intermediates", action="store_true")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    started = time.time()
    required = [args.input_bam, args.input_bam + ".bai", args.upstream_bed,
                args.downstream_bed, args.universal_utr_bed,
                args.extended_utr_bed, args.cassette_exon_bed]
    for path in required:
        require_file(path)

    output_dir = os.path.abspath(args.output_dir)
    work_dir = os.path.abspath(args.work_dir or os.path.join(output_dir, "work"))
    bam_dir = os.path.join(output_dir, "bam")
    psi_dir = os.path.join(output_dir, "psi")
    os.makedirs(work_dir, exist_ok=True)
    os.makedirs(bam_dir, exist_ok=True)
    os.makedirs(psi_dir, exist_ok=True)

    pysam.quickcheck(args.input_bam)
    stage_paths = OrderedDict()
    stage_paths["raw"] = args.input_bam

    up_bam = os.path.join(work_dir, "01_up_pass.bam")
    print("Stage 1: upstream constitutive-exon filter")
    filter_exon_coverage(args.input_bam, args.upstream_bed, up_bam,
                         work_dir, "01_up", args.coverage_threshold, args.threads)
    stage_paths["upstream_pass"] = up_bam

    down_bam = os.path.join(work_dir, "02_down_pass.bam")
    print("Stage 2: downstream constitutive-exon filter")
    filter_exon_coverage(up_bam, args.downstream_bed, down_bam,
                         work_dir, "02_down", args.coverage_threshold, args.threads)
    stage_paths["downstream_pass"] = down_bam

    utr_bam = os.path.join(work_dir, "03_utr_pass.bam")
    print("Stage 3: universal 3'UTR filter")
    filter_exon_coverage(down_bam, args.universal_utr_bed, utr_bam,
                         work_dir, "03_utr", args.coverage_threshold, args.threads)
    stage_paths["universal_utr_pass"] = utr_bam

    long_bam = os.path.join(bam_dir, "long.bam")
    print("Stage 4a: long-3'UTR filter")
    filter_exon_coverage(utr_bam, args.extended_utr_bed, long_bam,
                         work_dir, "04_long", args.coverage_threshold, args.threads)

    overlap_bam = os.path.join(work_dir, "04_distal_overlap.bam")
    short_candidates = os.path.join(work_dir, "04_short_candidates.bam")
    print("Stage 4b: construct the true non-distal short-candidate complement")
    split_distal_overlap(utr_bam, args.extended_utr_bed, overlap_bam,
                         short_candidates, args.threads)

    short_bam = os.path.join(bam_dir, "short.bam")
    print("Stage 4c: poly(A) filter short candidates")
    filter_polya(short_candidates, short_bam, args.threads)
    pysam.quickcheck(long_bam, short_bam)

    print("Stage 5: calculate all, long, and short PSI")
    all_csv = os.path.join(psi_dir, "all.psi.csv")
    long_csv = os.path.join(psi_dir, "long.psi.csv")
    short_csv = os.path.join(psi_dir, "short.psi.csv")
    calculate_psi(down_bam, args.cassette_exon_bed, all_csv,
                  args.coverage_threshold)
    calculate_psi(long_bam, args.cassette_exon_bed, long_csv,
                  args.coverage_threshold)
    calculate_psi(short_bam, args.cassette_exon_bed, short_csv,
                  args.coverage_threshold)

    stage_paths["distal_overlap"] = overlap_bam
    stage_paths["short_candidates"] = short_candidates
    stage_paths["long_final"] = long_bam
    stage_paths["short_final"] = short_bam
    count_path = os.path.join(psi_dir, "filter_stage_counts.tsv")
    with open(count_path, "w") as handle:
        handle.write("stage\talignment_records\n")
        for stage, path in stage_paths.items():
            handle.write("{}\t{}\n".format(stage, count_records(path)))

    checksum_path = os.path.join(psi_dir, "SHA256SUMS")
    write_sha256([all_csv, long_csv, short_csv, count_path], checksum_path)
    with open(os.path.join(psi_dir, "PIPELINE_COMPLETE"), "w") as handle:
        handle.write("optimized PL-Seq filtering completed\n")

    elapsed = time.time() - started
    with open(os.path.join(psi_dir, "runtime.tsv"), "w") as handle:
        handle.write("implementation\telapsed_seconds\tpysam_version\n")
        handle.write("optimized\t{:.3f}\t{}\n".format(elapsed, pysam.__version__))

    if not args.keep_intermediates:
        # Retain stage-count evidence and final products, but remove large transient BAMs.
        for path in [up_bam, down_bam, utr_bam, overlap_bam, short_candidates]:
            for candidate in [path, path + ".bai"]:
                if os.path.exists(candidate):
                    os.remove(candidate)
        try:
            os.rmdir(work_dir)
        except OSError:
            pass

    print("Completed in {:.1f} seconds".format(elapsed))
    print("Output: {}".format(output_dir))


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        raise
