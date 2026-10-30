#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
tmp=$(mktemp -d "${TMPDIR:-/tmp}/plseq-hpc-test.XXXXXX")
trap 'rm -rf "$tmp"' EXIT

printf '>chrTest\nACGT\n' > "$tmp/dm6.fa"
bash "$ROOT/scripts/hpc/00_setup_project.sh" \
  --project "$tmp/project" --dm6-fasta "$tmp/dm6.fa" >/dev/null

[[ $(awk 'END {print NR-1}' "$tmp/project/metadata/run_manifest.tsv") == 17 ]]
[[ $(awk 'END {print NR-1}' "$tmp/project/metadata/figure6_control/manifest.tsv") == 4 ]]
[[ $(awk 'END {print NR-1}' "$tmp/project/metadata/figure6_elav5/manifest.tsv") == 4 ]]

for required in \
  references/dm6.fa \
  references/Fly_dm6.bed \
  references/beds/upstream_exons.bed \
  references/beds/downstream_exons.bed \
  references/beds/universal_utr.bed \
  references/beds/extended_utr.bed \
  references/beds/cassette_exons.bed \
  software/Pull-a-long/exon_covering.py \
  software/Pull-a-long/polyA_filtering.py \
  software/Pull-a-long/calculate_PSI.py \
  software/PLSeq-reanalysis/filter_and_psi.py \
  scripts/prepare_figure4_pileup.py
do
  [[ -s "$tmp/project/$required" ]] || {
    echo "ERROR: setup did not create $required" >&2
    exit 1
  }
done

if grep -R -n -E 'miura@uchc|/projects/Enterica/Pedro|#SBATCH --partition=general|#SBATCH --qos=general' \
  "$ROOT/scripts/hpc" "$ROOT/config/paths.example.yml"; then
  echo "ERROR: personal or site-specific setting remains in a portable template" >&2
  exit 1
fi

mkdir -p "$tmp/bin"
printf '%s\n' \
  '#!/usr/bin/env bash' \
  'printf "%s\n" "$@" > "$MOCK_SBATCH_OUT"' \
  'printf "123456\n"' \
  > "$tmp/bin/sbatch"
chmod +x "$tmp/bin/sbatch"

submitted_job=$(
  MOCK_SBATCH_OUT="$tmp/sbatch.args" \
  PATH="$tmp/bin:$PATH" \
  PLSEQ_HPC_PROJECT="$tmp/project" \
  PLSEQ_ARRAY_LIMIT=3 \
  PLSEQ_PARTITION=test-partition \
  PLSEQ_QOS=test-qos \
  bash "$ROOT/scripts/hpc/submit_level3.sh" download 2>/dev/null
)

[[ "$submitted_job" == 123456 ]]

grep -Fx -- '--array=0-16%3' "$tmp/sbatch.args" >/dev/null
grep -Fx -- '--partition=test-partition' "$tmp/sbatch.args" >/dev/null
grep -Fx -- '--qos=test-qos' "$tmp/sbatch.args" >/dev/null
grep -F -- "$tmp/project/logs/plseq_download.%A_%a.out" "$tmp/sbatch.args" >/dev/null

MOCK_SBATCH_OUT="$tmp/sbatch-fast.args" \
PATH="$tmp/bin:$PATH" \
PLSEQ_HPC_PROJECT="$tmp/project" \
bash "$ROOT/scripts/hpc/submit_level3.sh" filter-optimized >/dev/null

grep -F -- '03_filter_and_psi_optimized.sbatch' "$tmp/sbatch-fast.args" >/dev/null

echo "HPC portability check passed."
