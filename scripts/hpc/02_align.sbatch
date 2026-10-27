#!/usr/bin/env bash
#SBATCH --job-name=p37679364_align
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miura@uchc.edu
#SBATCH --output=logs/plseq_align.%A_%a.out
#SBATCH --error=logs/plseq_align.%A_%a.err

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

mkdir -p "$PROJECT/bam" "$PROJECT/qc" "$PROJECT/logs"
scratch="${SLURM_TMPDIR:-/tmp}/${run}"
mkdir -p "$scratch"

apptainer exec "$CONTAINER" minimap2 -ax splice -B 3 -O 3,20 \
  --junc-bed "$PROJECT/references/Fly_dm6.bed" \
  -t "$SLURM_CPUS_PER_TASK" \
  "$PROJECT/references/dm6.fa" "$PROJECT/fastq/${run}_1.fastq.gz" \
  | apptainer exec "$CONTAINER" samtools sort -@ "$SLURM_CPUS_PER_TASK" \
      -m 3G -T "$scratch/${run}.sort" -o "$scratch/${run}.bam"

apptainer exec "$CONTAINER" samtools index -@ "$SLURM_CPUS_PER_TASK" "$scratch/${run}.bam"
apptainer exec "$CONTAINER" samtools quickcheck -v "$scratch/${run}.bam"
apptainer exec "$CONTAINER" samtools flagstat "$scratch/${run}.bam" > "$PROJECT/qc/${run}.flagstat.txt"
apptainer exec "$CONTAINER" samtools stats "$scratch/${run}.bam" > "$PROJECT/qc/${run}.samtools_stats.txt"
mv "$scratch/${run}.bam" "$scratch/${run}.bam.bai" "$PROJECT/bam/"
