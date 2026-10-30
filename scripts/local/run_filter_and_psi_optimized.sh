#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_filter_and_psi_optimized.sh --input-bam BAM --beds-dir DIR --output-dir DIR [options]

Options:
  --python PATH       Python with pysam installed (default: python3)
  --threads N         Sorting/indexing threads (default: 4)
  --work-dir DIR      Intermediate working directory
  --keep-intermediates

The BED directory must contain upstream_exons.bed, downstream_exons.bed,
universal_utr.bed, extended_utr.bed, and cassette_exons.bed.
EOF
}

PYTHON=python3
THREADS=4
INPUT_BAM=
BEDS_DIR=
OUTPUT_DIR=
WORK_DIR=
KEEP=

while (($#)); do
    case "$1" in
        --input-bam) INPUT_BAM=$2; shift 2 ;;
        --beds-dir) BEDS_DIR=$2; shift 2 ;;
        --output-dir) OUTPUT_DIR=$2; shift 2 ;;
        --python) PYTHON=$2; shift 2 ;;
        --threads) THREADS=$2; shift 2 ;;
        --work-dir) WORK_DIR=$2; shift 2 ;;
        --keep-intermediates) KEEP=--keep-intermediates; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n "$INPUT_BAM" && -n "$BEDS_DIR" && -n "$OUTPUT_DIR" ]] || {
    usage >&2
    exit 2
}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROGRAM="$SCRIPT_DIR/../optimized/filter_and_psi.py"

command=(
    "$PYTHON" "$PROGRAM"
    --input-bam "$INPUT_BAM"
    --upstream-bed "$BEDS_DIR/upstream_exons.bed"
    --downstream-bed "$BEDS_DIR/downstream_exons.bed"
    --universal-utr-bed "$BEDS_DIR/universal_utr.bed"
    --extended-utr-bed "$BEDS_DIR/extended_utr.bed"
    --cassette-exon-bed "$BEDS_DIR/cassette_exons.bed"
    --output-dir "$OUTPUT_DIR"
    --threads "$THREADS"
)

[[ -z "$WORK_DIR" ]] || command+=(--work-dir "$WORK_DIR")
[[ -z "$KEEP" ]] || command+=("$KEEP")

"${command[@]}"
