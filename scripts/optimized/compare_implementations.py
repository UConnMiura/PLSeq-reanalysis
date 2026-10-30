#!/usr/bin/env python3
"""Strictly compare legacy and optimized PL-Seq BAM and PSI products."""

import argparse
import collections
import csv
import os
import sys

import pysam


def bam_records(path):
    records = collections.Counter()
    with pysam.AlignmentFile(path, "rb") as handle:
        for read in handle.fetch(until_eof=True):
            records[read.to_string()] += 1
    return records


def compare_bam(left, right):
    left_records = bam_records(left)
    right_records = bam_records(right)
    if left_records == right_records:
        print("PASS BAM: {} records: {}".format(os.path.basename(left),
                                                sum(left_records.values())))
        return True
    missing = left_records - right_records
    extra = right_records - left_records
    print("FAIL BAM: {} missing={} extra={}".format(
        os.path.basename(left), sum(missing.values()), sum(extra.values())))
    for label, difference in [("missing", missing), ("extra", extra)]:
        for record, count in list(difference.items())[:5]:
            print("  {} x{}: {}".format(label, count, record))
    return False


def read_psi(path):
    rows = {}
    with open(path, newline="") as handle:
        for row in csv.reader(handle):
            if row:
                rows[row[0]] = row[1:]
    return rows


def compare_psi(left, right):
    left_rows = read_psi(left)
    right_rows = read_psi(right)
    if left_rows == right_rows:
        print("PASS PSI: {} events: {}".format(os.path.basename(left), len(left_rows)))
        return True
    print("FAIL PSI: {}".format(os.path.basename(left)))
    for event in sorted(set(left_rows) | set(right_rows)):
        if left_rows.get(event) != right_rows.get(event):
            print("  {} legacy={} optimized={}".format(
                event, left_rows.get(event), right_rows.get(event)))
    return False


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-bam-dir", required=True)
    parser.add_argument("--legacy-psi-dir", required=True)
    parser.add_argument("--optimized-dir", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    passed = True
    for name in ["long.bam", "short.bam"]:
        passed &= compare_bam(
            os.path.join(args.legacy_bam_dir, name),
            os.path.join(args.optimized_dir, "bam", name))
    for name in ["all.psi.csv", "long.psi.csv", "short.psi.csv"]:
        passed &= compare_psi(
            os.path.join(args.legacy_psi_dir, name),
            os.path.join(args.optimized_dir, "psi", name))
    if not passed:
        print("Implementations are NOT equivalent", file=sys.stderr)
        return 1
    print("All final BAM records and PSI rows are exactly equivalent")
    return 0


if __name__ == "__main__":
    sys.exit(main())
