#!/bin/bash
#SBATCH --job-name=plseq-cat
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --array=0-27%4

set -euo pipefail

PROJECT="${PLSEQ_HPC_PROJECT:?Set PLSEQ_HPC_PROJECT to the Mantis project directory}"
IMAGE="${PLSEQ_APPTAINER_IMAGE:-$PROJECT/software/pullalong-core.sif}"
CATALOG="$PROJECT/metadata/figure3e_pileup_catalog_manifest.tsv"
SAMPLE_MAP="$PROJECT/metadata/analyzed_samples.tsv"
PREP="$PROJECT/scripts/prepare_figure4_pileup.py"
OUT_ROOT="$PROJECT/pileup_inputs/figure3e_catalog"

if ! command -v apptainer >/dev/null 2>&1 && command -v module >/dev/null 2>&1; then
    module load "${PLSEQ_APPTAINER_MODULE:-apptainer}"
fi
command -v apptainer >/dev/null 2>&1 || {
    echo "ERROR: apptainer is unavailable; load it or set PLSEQ_APPTAINER_MODULE" >&2
    exit 2
}

mkdir -p "$PROJECT/logs" "$OUT_ROOT"

IFS=$'\t' read -r array_index order gene event probe_set strand chrom \
  view_start view_end ce_start ce_end up_start up_end down_start down_end \
  common_start common_end distal_start distal_end selection_depth \
  < <(awk -F '\t' -v i="$SLURM_ARRAY_TASK_ID" 'NR>1 && $1==i {print; exit}' "$CATALOG")

if [[ -z "${gene:-}" ]]; then
    echo "No catalog row for array index $SLURM_ARRAY_TASK_ID" >&2
    exit 1
fi

gene_out="$OUT_ROOT/$gene"
mkdir -p "$gene_out"
sample_manifest="$gene_out/sample_manifest.tsv"

printf 'run_accession\treplicate\n' > "$sample_manifest"
awk -F '\t' -v probe="$probe_set" '
  BEGIN {OFS="\t"}
  NR>1 && $8=="control" && $10==probe {print $1,$11}
' "$SAMPLE_MAP" | sort -t $'\t' -k2,2n >> "$sample_manifest"

echo "Preparing $gene ($event), $probe_set, $(($(wc -l < "$sample_manifest") - 1)) replicates"

apptainer exec --bind "$PROJECT:$PROJECT" "$IMAGE" \
  python3 "$PREP" \
    --manifest "$sample_manifest" \
    --bam-root "$PROJECT/filtered_bam" \
    --psi-root "$PROJECT/psi" \
    --output "$gene_out" \
    --gene "$gene" \
    --event "$event" \
    --chrom "$chrom" \
    --region-start "$view_start" \
    --region-end "$view_end" \
    --ce-start "$ce_start" \
    --ce-end "$ce_end" \
    --reads-per-replicate 20 \
    --allow-count-mismatch

echo "Completed $gene at $(date --iso-8601=seconds)"
