#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$SCRIPT_DIR"

mkdir -p output/pdf
MAINFONT="${PLSEQ_MAINFONT:-Libertinus Sans}"

pandoc manuscript.qmd \
  --from=markdown-implicit_figures \
  --standalone \
  --citeproc \
  --csl=nature.csl \
  --lua-filter=figure-box.lua \
  --include-in-header=reader-preamble.tex \
  --pdf-engine=xelatex \
  -V mainfont="$MAINFONT" \
  -V colorlinks=true \
  -V linkcolor=PaperBlue \
  -V urlcolor=PaperBlue \
  -V citecolor=PaperBlue \
  -V papersize=letter \
  -V fontsize=10pt \
  -o output/pdf/PLSeq_reanalysis_reader_manuscript.pdf

printf 'Wrote %s\n' "$SCRIPT_DIR/output/pdf/PLSeq_reanalysis_reader_manuscript.pdf"
