#!/usr/bin/env bash
#SBATCH --job-name=p37679364_dl
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --time=24:00:00

set -euo pipefail
PROJECT="${PLSEQ_HPC_PROJECT:?Set PLSEQ_HPC_PROJECT to the Mantis project directory}"

manifest="$PROJECT/metadata/run_manifest.tsv"
line=$((SLURM_ARRAY_TASK_ID + 2))
IFS=$'\t' read -r run url bytes md5 < <(sed -n "${line}p" "$manifest")
[[ -n "$run" && -n "$url" && -n "$md5" ]]

mkdir -p "$PROJECT/fastq" "$PROJECT/logs"
target="$PROJECT/fastq/${run}_1.fastq.gz"
case "$url" in
  http://*|https://*) download_url="$url" ;;
  *) download_url="https://${url}" ;;
esac

if [[ -s "$target" ]]; then
    printf '%s  %s\n' "$md5" "$target" | md5sum -c -
    echo "Already downloaded and validated: $target"
    exit 0
fi

curl -L --fail --retry 5 --retry-delay 10 "$download_url" -o "$target.part"
printf '%s  %s\n' "$md5" "$target.part" | md5sum -c -
if [[ -n "${bytes:-}" ]]; then
    observed_bytes=$(stat -c %s "$target.part")
    [[ "$observed_bytes" == "$bytes" ]] || {
        echo "ERROR: expected $bytes bytes, observed $observed_bytes" >&2
        exit 4
    }
fi
mv "$target.part" "$target"
