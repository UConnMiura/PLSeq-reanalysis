#!/usr/bin/env bash
#SBATCH --job-name=p37679364_psi
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --time=12:00:00

set -euo pipefail

PROJECT="${PLSEQ_HPC_PROJECT:?Set PLSEQ_HPC_PROJECT to the HPC project directory}"
CONTAINER="${PLSEQ_APPTAINER_IMAGE:-$PROJECT/software/pullalong-core.sif}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

if ! command -v apptainer >/dev/null 2>&1 && command -v module >/dev/null 2>&1; then
    module load "${PLSEQ_APPTAINER_MODULE:-apptainer}"
fi
command -v apptainer >/dev/null 2>&1 || {
    echo "ERROR: apptainer is unavailable; load it or set PLSEQ_APPTAINER_MODULE" >&2
    exit 2
}

manifest="$PROJECT/metadata/run_manifest.tsv"
line=$((SLURM_ARRAY_TASK_ID + 2))
run=$(sed -n "${line}p" "$manifest" | cut -f1)
[[ -n "$run" ]] || { echo "ERROR: no manifest row for task $SLURM_ARRAY_TASK_ID" >&2; exit 2; }

REF="$PROJECT/references"
INPUT="$PROJECT/bam/${run}.bam"
WORK="$PROJECT/work/filter/${run}"
BAMOUT="$PROJECT/filtered_bam/${run}"
PSIOUT="$PROJECT/psi/${run}"

mkdir -p "$WORK" "$BAMOUT" "$PSIOUT" "$PROJECT/logs"
cd "$WORK"  # The published filtering programs use fixed temporary filenames.

for required in \
    "$INPUT" \
    "$INPUT.bai" \
    "$CONTAINER" \
    "$REF/beds/upstream_exons.bed" \
    "$REF/beds/downstream_exons.bed" \
    "$REF/beds/universal_utr.bed" \
    "$REF/beds/extended_utr.bed" \
    "$REF/beds/cassette_exons.bed" \
    "$PROJECT/software/Pull-a-long/exon_covering.py" \
    "$PROJECT/software/Pull-a-long/polyA_filtering.py" \
    "$PROJECT/software/Pull-a-long/calculate_PSI.py"
do
    [[ -s "$required" ]] || { echo "ERROR: missing or empty: $required" >&2; exit 3; }
done

run_in_image() {
    apptainer exec --bind "$PROJECT:$PROJECT" "$CONTAINER" "$@"
}

repo="$PROJECT/software/Pull-a-long"
beds="$REF/beds"

echo "Run: $run"
echo "Started: $(date --iso-8601=seconds)"
run_in_image samtools quickcheck -v "$INPUT"

echo "Stage 1: upstream constitutive-exon filter"
run_in_image python3 "$repo/exon_covering.py" \
  -i "$INPUT" -b "$beds/upstream_exons.bed" -p 0.8 -t "$THREADS" \
  -o "$WORK/01_up_pass.bam" -r "$WORK/01_up_fail.bam"

echo "Stage 2: downstream constitutive-exon filter"
run_in_image python3 "$repo/exon_covering.py" \
  -i "$WORK/01_up_pass.bam" -b "$beds/downstream_exons.bed" -p 0.8 -t "$THREADS" \
  -o "$WORK/02_down_pass.bam" -r "$WORK/02_down_fail.bam"

echo "Stage 3: universal 3'UTR filter"
run_in_image python3 "$repo/exon_covering.py" \
  -i "$WORK/02_down_pass.bam" -b "$beds/universal_utr.bed" -p 0.8 -t "$THREADS" \
  -o "$WORK/03_utr_pass.bam" -r "$WORK/03_utr_fail.bam"

echo "Stage 4a: long-3'UTR filter"
run_in_image python3 "$repo/exon_covering.py" \
  -i "$WORK/03_utr_pass.bam" -b "$beds/extended_utr.bed" -p 0.8 -t "$THREADS" \
  -o "$BAMOUT/long.bam" -r "$WORK/04_distal_partial_fail.bam"

echo "Stage 4b: construct the true non-distal short-candidate complement"
run_in_image samtools view -h -b "$WORK/03_utr_pass.bam" \
  -L "$beds/extended_utr.bed" \
  -U "$WORK/04_short_candidates.bam" \
  -o "$WORK/04_distal_overlap.bam"
run_in_image samtools index -@ "$THREADS" "$WORK/04_short_candidates.bam"

echo "Stage 4c: poly(A) filter short candidates"
run_in_image python3 "$repo/polyA_filtering.py" \
  -i "$WORK/04_short_candidates.bam" -t "$THREADS" -o "$BAMOUT/short.bam"

echo "Stage 5: validate final BAMs"
run_in_image samtools quickcheck -v "$BAMOUT/long.bam" "$BAMOUT/short.bam"

echo "Stage 6: calculate all, long, and short PSI"
run_in_image python3 "$repo/calculate_PSI.py" \
  -i "$WORK/02_down_pass.bam" -b "$beds/cassette_exons.bed" -t "$THREADS" \
  -o "$PSIOUT/all.psi.csv"
run_in_image python3 "$repo/calculate_PSI.py" \
  -i "$BAMOUT/long.bam" -b "$beds/cassette_exons.bed" -t "$THREADS" \
  -o "$PSIOUT/long.psi.csv"
run_in_image python3 "$repo/calculate_PSI.py" \
  -i "$BAMOUT/short.bam" -b "$beds/cassette_exons.bed" -t "$THREADS" \
  -o "$PSIOUT/short.psi.csv"

echo "Stage 7: write QC counts"
COUNT_FILE="$PSIOUT/filter_stage_counts.tsv"
printf 'stage\talignment_records\n' > "$COUNT_FILE"
count_bam() {
    local stage=$1 bam=$2
    printf '%s\t%s\n' "$stage" "$(run_in_image samtools view -c "$bam")" >> "$COUNT_FILE"
}
count_bam raw                "$INPUT"
count_bam upstream_pass      "$WORK/01_up_pass.bam"
count_bam downstream_pass    "$WORK/02_down_pass.bam"
count_bam universal_utr_pass "$WORK/03_utr_pass.bam"
count_bam distal_overlap     "$WORK/04_distal_overlap.bam"
count_bam short_candidates   "$WORK/04_short_candidates.bam"
count_bam long_final         "$BAMOUT/long.bam"
count_bam short_final        "$BAMOUT/short.bam"

sha256sum "$PSIOUT/all.psi.csv" "$PSIOUT/long.psi.csv" \
  "$PSIOUT/short.psi.csv" "$COUNT_FILE" > "$PSIOUT/SHA256SUMS"
touch "$PSIOUT/PIPELINE_COMPLETE"
echo "Completed: $(date --iso-8601=seconds)"
