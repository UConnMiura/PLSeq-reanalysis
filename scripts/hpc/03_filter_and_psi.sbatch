#!/usr/bin/env bash
#SBATCH --job-name=p37679364_psi
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miura@uchc.edu
#SBATCH --output=logs/plseq_psi.%A_%a.out
#SBATCH --error=logs/plseq_psi.%A_%a.err

set -euo pipefail
PROJECT="${PLSEQ_HPC_PROJECT:?Set PLSEQ_HPC_PROJECT to the Mantis project directory}"
CONTAINER="${PLSEQ_APPTAINER_IMAGE:-$PROJECT/software/pullalong-core.sif}"

module purge
module load apptainer/1.4.4-gcc-11.4.0-znzbvwm
apptainer --version

manifest="$PROJECT/metadata/run_manifest.tsv"
line=$((SLURM_ARRAY_TASK_ID + 2))
run=$(sed -n "${line}p" "$manifest" | cut -f1)
[[ -n "$run" ]]

scratch="${SLURM_TMPDIR:-/tmp}/${run}_filter"
mkdir -p "$scratch" "$PROJECT/psi" "$PROJECT/filtered_bam" "$PROJECT/logs"
cd "$scratch"  # Required: authors' scripts use fixed temporary filenames.

repo="$PROJECT/software/Pull-a-long"
beds="$PROJECT/references/beds"
threads="$SLURM_CPUS_PER_TASK"
input="$PROJECT/bam/${run}.bam"

# The BED package must contain one reviewed BED4 record set per filtering stage.
apptainer exec "$CONTAINER" python "$repo/exon_covering.py" -i "$input" -b "$beds/upstream_exons.bed" -o up.bam -r up_fail.bam -t "$threads"
apptainer exec "$CONTAINER" python "$repo/exon_covering.py" -i up.bam -b "$beds/downstream_exons.bed" -o down.bam -r down_fail.bam -t "$threads"
apptainer exec "$CONTAINER" python "$repo/exon_covering.py" -i down.bam -b "$beds/universal_utr.bed" -o utr.bam -r utr_fail.bam -t "$threads"
apptainer exec "$CONTAINER" python "$repo/exon_covering.py" -i utr.bam -b "$beds/extended_utr.bed" -o long.bam -r non_long.bam -t "$threads"
apptainer exec "$CONTAINER" python "$repo/polyA_filtering.py" -i non_long.bam -o short.bam -t "$threads"
apptainer exec "$CONTAINER" python "$repo/calculate_PSI.py" -i long.bam -b "$beds/cassette_exons.bed" -o "$PROJECT/psi/${run}.long.csv" -t "$threads"
apptainer exec "$CONTAINER" python "$repo/calculate_PSI.py" -i short.bam -b "$beds/cassette_exons.bed" -o "$PROJECT/psi/${run}.short.csv" -t "$threads"

apptainer exec "$CONTAINER" samtools quickcheck -v long.bam short.bam
mv long.bam "$PROJECT/filtered_bam/${run}.long.bam"
mv long.bam.bai "$PROJECT/filtered_bam/${run}.long.bam.bai"
mv short.bam "$PROJECT/filtered_bam/${run}.short.bam"
mv short.bam.bai "$PROJECT/filtered_bam/${run}.short.bam.bai"
