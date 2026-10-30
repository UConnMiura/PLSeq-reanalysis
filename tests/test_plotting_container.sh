#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 /path/to/plseq-plotting.sif /path/to/output" >&2
    exit 2
fi

image=$(readlink -f "$1")
output=$(mkdir -p "$2" && readlink -f "$2")
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

apptainer test "$image"
apptainer exec \
    --bind "$repo:$repo" \
    --bind "$output:$output" \
    "$image" \
    python3 "$repo/tests/smoke_test_plotting_container.py" --output "$output"

apptainer exec \
    --bind "$output:$output" \
    "$image" \
    pdfinfo "$output/plseq_plotting_smoke.pdf" \
    | grep -q '^Pages:[[:space:]]*1$'
test -s "$output/plseq_plotting_smoke.png"
test -s "$output/plseq_plotting_smoke.pdf"
echo "Plotting container smoke test passed"
