#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  00_setup_project.sh --project DIR --dm6-fasta FILE [--download-container]

Creates the Level 3 HPC directory layout from a repository clone. The dm6
FASTA must be supplied explicitly. With --download-container, the exact SIF
archived in Zenodo record 23769044 is downloaded and checksum-verified.
EOF
}

PROJECT="${PLSEQ_HPC_PROJECT:-}"
DM6_FASTA=""
DOWNLOAD_CONTAINER=0
while (($#)); do
    case "$1" in
        --project) PROJECT=$2; shift 2 ;;
        --dm6-fasta) DM6_FASTA=$2; shift 2 ;;
        --download-container) DOWNLOAD_CONTAINER=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n "$PROJECT" ]] || { echo "ERROR: provide --project or PLSEQ_HPC_PROJECT" >&2; exit 2; }
[[ -s "$DM6_FASTA" ]] || { echo "ERROR: provide an existing dm6 FASTA with --dm6-fasta" >&2; exit 2; }

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/../.." && pwd)
PROJECT=$(mkdir -p "$PROJECT" && cd "$PROJECT" && pwd)

python3 "$REPO_ROOT/tests/check_source_checksums.py"

mkdir -p \
  "$PROJECT"/{bam,fastq,filtered_bam,logs,metadata,qc,psi,work} \
  "$PROJECT/references/beds" \
  "$PROJECT/software/Pull-a-long" \
  "$PROJECT/software/PLSeq-reanalysis" \
  "$PROJECT/scripts" \
  "$PROJECT/metadata/figure6_control" \
  "$PROJECT/metadata/figure6_elav5"

install -m 0644 "$DM6_FASTA" "$PROJECT/references/dm6.fa"
install -m 0644 "$REPO_ROOT/data/source/Fly_dm6.bed" "$PROJECT/references/Fly_dm6.bed"
install -m 0644 "$REPO_ROOT/data/source/up.figure3e.reconstructed.bed" "$PROJECT/references/beds/upstream_exons.bed"
install -m 0644 "$REPO_ROOT/data/source/down.figure3e.reconstructed.bed" "$PROJECT/references/beds/downstream_exons.bed"
install -m 0644 "$REPO_ROOT/data/source/UTR.figure3e.bed" "$PROJECT/references/beds/universal_utr.bed"
install -m 0644 "$REPO_ROOT/data/source/dis.figure3e.bed" "$PROJECT/references/beds/extended_utr.bed"
install -m 0644 "$REPO_ROOT/data/source/SE_PSI.figure3e.suppdata9.bed" "$PROJECT/references/beds/cassette_exons.bed"

for program in exon_covering.py polyA_filtering.py calculate_PSI.py; do
    install -m 0644 "$REPO_ROOT/data/source/$program" "$PROJECT/software/Pull-a-long/$program"
done
install -m 0755 "$REPO_ROOT/scripts/optimized/filter_and_psi.py" \
  "$PROJECT/software/PLSeq-reanalysis/filter_and_psi.py"
install -m 0755 "$REPO_ROOT/scripts/local/prepare_figure4_pileup.py" "$PROJECT/scripts/prepare_figure4_pileup.py"
install -m 0644 "$REPO_ROOT/config/figure3e_pileup_catalog_manifest.tsv" \
  "$PROJECT/metadata/figure3e_pileup_catalog_manifest.tsv"
install -m 0644 "$REPO_ROOT/supplementary_tables/Supplementary_Table_1_samples.tsv" \
  "$PROJECT/metadata/analyzed_samples.tsv"

# The Level 3 workflow defaults to the 17 embryo runs actually reanalyzed, not
# all 49 accessions deposited with the source study.
awk -F '\t' 'BEGIN {OFS="\t"; print "run","url","bytes","md5"}
  NR>1 {print $1,$15,$17,$16}' \
  "$REPO_ROOT/supplementary_tables/Supplementary_Table_1_samples.tsv" \
  > "$PROJECT/metadata/run_manifest.tsv"

awk -F '\t' 'BEGIN {OFS="\t"; print "run_accession","replicate"}
  NR>1 && $8=="control" && $10=="p3" {print $1,$11}' \
  "$REPO_ROOT/supplementary_tables/Supplementary_Table_1_samples.tsv" \
  > "$PROJECT/metadata/figure6_control/manifest.tsv"
awk -F '\t' 'BEGIN {OFS="\t"; print "run_accession","replicate"}
  NR>1 && $8=="elav5" && $10=="p3" {print $1,$11}' \
  "$REPO_ROOT/supplementary_tables/Supplementary_Table_1_samples.tsv" \
  > "$PROJECT/metadata/figure6_elav5/manifest.tsv"

if ((DOWNLOAD_CONTAINER)); then
    command -v curl >/dev/null 2>&1 || { echo "ERROR: curl is required" >&2; exit 3; }
    image="$PROJECT/software/pullalong-core.sif"
    curl -L --fail --retry 5 \
      'https://zenodo.org/records/23769044/files/pullalong-core.sif?download=1' \
      -o "$image.part"
    printf '%s  %s\n' \
      '6037ffe045c08e7435e8a997b890d791f6f3554ab01cef16b0b895761c199beb' \
      "$image.part" | sha256sum -c -
    mv "$image.part" "$image"
fi

cat <<EOF
Level 3 project prepared at: $PROJECT
Run manifest: $PROJECT/metadata/run_manifest.tsv ($(($(wc -l < "$PROJECT/metadata/run_manifest.tsv") - 1)) runs)

Next:
  export PLSEQ_HPC_PROJECT="$PROJECT"
  export PLSEQ_APPTAINER_IMAGE="$PROJECT/software/pullalong-core.sif"
  bash "$REPO_ROOT/scripts/hpc/submit_level3.sh" download
EOF
