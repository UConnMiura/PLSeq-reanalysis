#!/usr/bin/env bash
#SBATCH --job-name=p37679364_fastpsi
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --time=04:00:00

set -euo pipefail

PROJECT="${PLSEQ_HPC_PROJECT:?Set PLSEQ_HPC_PROJECT to the HPC project directory}"
CONTAINER="${PLSEQ_APPTAINER_IMAGE:-$PROJECT/software/pullalong-core.sif}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

if ! command -v apptainer >/dev/null 2>&1 && command -v module >/dev/null 2>&1; then
    module load "${PLSEQ_APPTAINER_MODULE:-apptainer}"
fi
command -v apptainer >/dev/null 2>&1 || {
    echo "ERROR: apptainer is unavailable" >&2
    exit 2
}

manifest="$PROJECT/metadata/run_manifest.tsv"
line=$((SLURM_ARRAY_TASK_ID + 2))
run=$(sed -n "${line}p" "$manifest" | cut -f1)
[[ -n "$run" ]] || { echo "ERROR: no manifest row for task $SLURM_ARRAY_TASK_ID" >&2; exit 2; }

PROGRAM="$PROJECT/software/PLSeq-reanalysis/filter_and_psi.py"
REF="$PROJECT/references/beds"
INPUT="$PROJECT/bam/${run}.bam"
OUTPUT="$PROJECT/optimized/${run}"
SCRATCH_BASE="${SLURM_TMPDIR:-$PROJECT/work/optimized/${run}}"
WORK="$SCRATCH_BASE/filter_work"

mkdir -p "$OUTPUT" "$WORK" "$PROJECT/logs"

echo "Run: $run"
echo "Input: $INPUT"
echo "Output: $OUTPUT"
echo "Work: $WORK"
echo "Started: $(date --iso-8601=seconds)"

apptainer exec \
  --bind "$PROJECT:$PROJECT" \
  --bind "$SCRATCH_BASE:$SCRATCH_BASE" \
  "$CONTAINER" \
  python3 "$PROGRAM" \
    --input-bam "$INPUT" \
    --upstream-bed "$REF/upstream_exons.bed" \
    --downstream-bed "$REF/downstream_exons.bed" \
    --universal-utr-bed "$REF/universal_utr.bed" \
    --extended-utr-bed "$REF/extended_utr.bed" \
    --cassette-exon-bed "$REF/cassette_exons.bed" \
    --output-dir "$OUTPUT" \
    --work-dir "$WORK" \
    --threads "$THREADS"

echo "Completed: $(date --iso-8601=seconds)"
