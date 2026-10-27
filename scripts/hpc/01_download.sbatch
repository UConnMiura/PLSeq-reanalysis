#!/usr/bin/env bash
#SBATCH --job-name=p37679364_dl
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=miura@uchc.edu
#SBATCH --output=/home/FCAM/pmiura/eofiles/plseq_download.%A_%a.out
#SBATCH --error=/home/FCAM/pmiura/eofiles/plseq_download.%A_%a.err

set -euo pipefail
: "${PROJECT:?Set PROJECT to the Mantis project directory}"

manifest="$PROJECT/metadata/run_manifest.tsv"
line=$((SLURM_ARRAY_TASK_ID + 2))
IFS=$'\t' read -r run url bytes md5 < <(sed -n "${line}p" "$manifest")
[[ -n "$run" && -n "$url" && -n "$md5" ]]

mkdir -p "$PROJECT/fastq" "$PROJECT/logs"
target="$PROJECT/fastq/${run}_1.fastq.gz"
curl -L --fail --retry 5 --retry-delay 10 "https://${url}" -o "$target.part"
printf '%s  %s\n' "$md5" "$target.part" | md5sum -c -
mv "$target.part" "$target"
gzip -t "$target"
