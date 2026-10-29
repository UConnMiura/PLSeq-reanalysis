#!/usr/bin/env bash
set -euo pipefail

stage=${1:-}
case "$stage" in
    download) script=01_download.sbatch; stem=download ;;
    align) script=02_align.sbatch; stem=align ;;
    filter) script=03_filter_and_psi.sbatch; stem=filter ;;
    figure3e-pileups) script=figure3e_pileup_catalog.sbatch; stem=figure3e_pileups ;;
    figure5-supplement) script=figure5_supplement_pileups.sbatch; stem=figure5_supplement ;;
    *) echo "Usage: $0 {download|align|filter|figure3e-pileups|figure5-supplement}" >&2; exit 2 ;;
esac

PROJECT="${PLSEQ_HPC_PROJECT:?Set PLSEQ_HPC_PROJECT to the prepared HPC project directory}"
manifest="$PROJECT/metadata/run_manifest.tsv"
[[ -s "$manifest" ]] || { echo "ERROR: missing $manifest; run 00_setup_project.sh" >&2; exit 3; }
mkdir -p "$PROJECT/logs"

limit="${PLSEQ_ARRAY_LIMIT:-4}"
case "$stage" in
    figure3e-pileups|figure5-supplement) array="0-27%${limit}" ;;
    *)
        tasks=$(($(wc -l < "$manifest") - 1))
        ((tasks > 0)) || { echo "ERROR: manifest has no data rows" >&2; exit 3; }
        array="0-$((tasks - 1))%${limit}"
        ;;
esac

args=(--parsable --array="$array"
      --output="$PROJECT/logs/plseq_${stem}.%A_%a.out"
      --error="$PROJECT/logs/plseq_${stem}.%A_%a.err")
[[ -n "${PLSEQ_PARTITION:-}" ]] && args+=(--partition="$PLSEQ_PARTITION")
[[ -n "${PLSEQ_QOS:-}" ]] && args+=(--qos="$PLSEQ_QOS")
if [[ -n "${PLSEQ_MAIL_USER:-}" ]]; then
    args+=(--mail-type=END,FAIL --mail-user="$PLSEQ_MAIL_USER")
fi
if [[ -n "${PLSEQ_SBATCH_ARGS:-}" ]]; then
    read -r -a extra <<< "$PLSEQ_SBATCH_ARGS"
    args+=("${extra[@]}")
fi

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
job_id=$(sbatch "${args[@]}" "$SCRIPT_DIR/$script")
echo "Submitted $stage array ($array): $job_id"
