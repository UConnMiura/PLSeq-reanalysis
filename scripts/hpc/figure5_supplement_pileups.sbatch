#!/bin/bash
#SBATCH --job-name=fig5-supp
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=0-27%4
#SBATCH --output=/projects/Enterica/Pedro/paper/plseq/logs/fig5_supp_%A_%a.out
#SBATCH --error=/projects/Enterica/Pedro/paper/plseq/logs/fig5_supp_%A_%a.err

set -euo pipefail

PROJECT=/projects/Enterica/Pedro/paper/plseq
IMAGE="$PROJECT/software/pullalong-core.sif"
CATALOG="$PROJECT/metadata/figure3e_pileup_catalog_manifest.tsv"
PREP="$PROJECT/scripts/prepare_figure4_pileup.py"
CONTROL_MANIFEST="$PROJECT/metadata/figure6_control/manifest.tsv"
MUTANT_MANIFEST="$PROJECT/metadata/figure6_elav5/filter_manifest.tsv"
OUT_ROOT="$PROJECT/pileup_inputs/figure5_supplement"

module load apptainer/1.4.4-gcc-11.4.0-znzbvwm

mkdir -p "$PROJECT/logs" "$OUT_ROOT"

IFS=$'\t' read -r array_index order gene event probe_set strand chrom \
  view_start view_end ce_start ce_end up_start up_end down_start down_end \
  common_start common_end distal_start distal_end selection_depth \
  < <(awk -F '\t' -v i="$SLURM_ARRAY_TASK_ID" \
      'NR>1 && $1==i {print; exit}' "$CATALOG")

if [[ -z "${gene:-}" ]]; then
    echo "ERROR: no catalog row for array index $SLURM_ARRAY_TASK_ID" >&2
    exit 2
fi

# These two events are displayed in main Figure 5 and are not duplicated in
# the supplemental atlas.
if [[ "$gene" == "Khc-73" || "$gene" == "Dscam1" ]]; then
    echo "Skipping $gene: displayed in main Figure 5"
    exit 0
fi

prepare_condition () {
    condition=$1
    manifest=$2
    bam_root=$3
    psi_root=$4
    output="$OUT_ROOT/$gene/$condition"

    mkdir -p "$output"
    echo "Preparing $gene ($event): $condition"

    apptainer exec \
      --bind "$PROJECT:$PROJECT" \
      "$IMAGE" \
      python3 "$PREP" \
        --manifest "$manifest" \
        --bam-root "$bam_root" \
        --psi-root "$psi_root" \
        --output "$output" \
        --gene "$gene" \
        --event "$event" \
        --chrom "$chrom" \
        --region-start "$view_start" \
        --region-end "$view_end" \
        --ce-start "$ce_start" \
        --ce-end "$ce_end" \
        --reads-per-replicate 20 \
        --seed 37679364 \
        --allow-count-mismatch
}

prepare_condition \
  control \
  "$CONTROL_MANIFEST" \
  "$PROJECT/filtered_bam" \
  "$PROJECT/psi"

prepare_condition \
  elav5 \
  "$MUTANT_MANIFEST" \
  "$PROJECT/filtered_bam/figure6_elav5" \
  "$PROJECT/psi/figure6_elav5"

echo "Completed $gene at $(date --iso-8601=seconds)"
