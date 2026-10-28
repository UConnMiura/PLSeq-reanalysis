#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$SCRIPT_DIR"

BUILD_DIR="tmp/pdfs/biorxiv"
OUTPUT_DIR="output/pdf"
SOURCE="$BUILD_DIR/manuscript_biorxiv.qmd"
OUTPUT="$OUTPUT_DIR/PLSeq_reanalysis_biorxiv_submission.pdf"

mkdir -p "$BUILD_DIR" "$OUTPUT_DIR"

python3 prepare_biorxiv_manuscript.py manuscript.qmd "$SOURCE"

pandoc "$SOURCE" \
  --from=markdown-implicit_figures \
  --standalone \
  --citeproc \
  --resource-path=".:$SCRIPT_DIR" \
  --csl=nature.csl \
  --lua-filter=figure-box.lua \
  --include-in-header=reader-preamble.tex \
  --pdf-engine=xelatex \
  -V mainfont='Helvetica' \
  -V colorlinks=true \
  -V linkcolor=PaperBlue \
  -V urlcolor=PaperBlue \
  -V citecolor=PaperBlue \
  -V papersize=letter \
  -V fontsize=10pt \
  -o "$OUTPUT"

printf 'Wrote %s\n' "$OUTPUT"
