#!/usr/bin/env bash
#SBATCH --job-name=p37679364_align
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

set -euo pipefail
PROJECT="${PLSEQ_HPC_PROJECT:?Set PLSEQ_HPC_PROJECT to the Mantis project directory}"
CONTAINER="${PLSEQ_APPTAINER_IMAGE:-$PROJECT/software/pullalong-core.sif}"

if ! command -v apptainer >/dev/null 2>&1 && command -v module >/dev/null 2>&1; then
    module load "${PLSEQ_APPTAINER_MODULE:-apptainer}"
fi
command -v apptainer >/dev/null 2>&1 || {
    echo "ERROR: apptainer is unavailable; load it or set PLSEQ_APPTAINER_MODULE" >&2
    exit 2
}
apptainer --version

manifest="$PROJECT/metadata/run_manifest.tsv"
line=$((SLURM_ARRAY_TASK_ID + 2))
run=$(sed -n "${line}p" "$manifest" | cut -f1)
[[ -n "$run" ]]

mkdir -p "$PROJECT/bam" "$PROJECT/qc" "$PROJECT/logs"
for required in "$CONTAINER" "$PROJECT/references/Fly_dm6.bed" \
                "$PROJECT/references/dm6.fa" "$PROJECT/fastq/${run}_1.fastq.gz"; do
    [[ -s "$required" ]] || { echo "ERROR: missing or empty: $required" >&2; exit 3; }
done

scratch_root="${SLURM_TMPDIR:-$PROJECT/work/scratch}"
scratch="$scratch_root/${SLURM_JOB_ID:-manual}_${run}"
mkdir -p "$scratch"

apptainer exec --bind "$PROJECT:$PROJECT" "$CONTAINER" minimap2 -ax splice -B 3 -O 3,20 \
  --junc-bed "$PROJECT/references/Fly_dm6.bed" \
  -t "$SLURM_CPUS_PER_TASK" \
  "$PROJECT/references/dm6.fa" "$PROJECT/fastq/${run}_1.fastq.gz" \
  | apptainer exec --bind "$PROJECT:$PROJECT" "$CONTAINER" samtools sort -@ "$SLURM_CPUS_PER_TASK" \
      -m 3G -T "$scratch/${run}.sort" -o "$scratch/${run}.bam"

apptainer exec --bind "$PROJECT:$PROJECT" "$CONTAINER" samtools index -@ "$SLURM_CPUS_PER_TASK" "$scratch/${run}.bam"
apptainer exec --bind "$PROJECT:$PROJECT" "$CONTAINER" samtools quickcheck -v "$scratch/${run}.bam"
apptainer exec --bind "$PROJECT:$PROJECT" "$CONTAINER" samtools flagstat "$scratch/${run}.bam" > "$PROJECT/qc/${run}.flagstat.txt"
apptainer exec --bind "$PROJECT:$PROJECT" "$CONTAINER" samtools stats "$scratch/${run}.bam" > "$PROJECT/qc/${run}.samtools_stats.txt"
mv "$scratch/${run}.bam" "$scratch/${run}.bam.bai" "$PROJECT/bam/"
