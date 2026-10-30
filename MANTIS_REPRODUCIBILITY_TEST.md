# Clean Mantis reproducibility test

## Purpose

Run a clean, one-sample test of the PL-Seq reconstruction workflow on the UConn School of Medicine HPC cluster. This test must start from the public GitHub repository, obtain the archived Apptainer image and public dm6 reference independently, and avoid using files from the original analysis directory.

The pilot sample is **SRR16318553**. It tests the three computational levels in sequence:

1. Download the archived FASTQ file.
2. Align it to dm6 with the exact archived software environment.
3. Reconstruct the long- and short-3'UTR molecule classes and calculate PSI.

Do not modify or delete anything under:

```text
/projects/Enterica/Pedro/paper/plseq
```

Use the exact Apptainer image throughout. Do not substitute host-installed minimap2, SAMtools, pysam, or other analysis software.

## Expected public inputs

- GitHub repository: <https://github.com/UConnMiura/PLSeq-reanalysis>
- Repository commit tested when this runbook was written: `bb8724a7ed0709c0b9021074c19dd4fb792ec687`
- Zenodo container record: <https://doi.org/10.5281/zenodo.23769044>
- Expected container SHA-256: `6037ffe045c08e7435e8a997b890d791f6f3554ab01cef16b0b895761c199beb`
- UCSC dm6 FASTA: <https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz>
- Expected dm6 archive MD5: `4a5777324403eff92c3650c30af30120`

## Rules for the executing agent

- Run each numbered section in order.
- Stop immediately if a checksum, required-file check, BAM validation, or Slurm job fails.
- Report the exact command and complete error output when stopping.
- Do not silently replace a missing dependency or reference with a file from the original project.
- Do not submit the alignment until the download job has completed successfully.
- Do not submit filtering/PSI reconstruction until the alignment job has completed successfully.
- Record all Slurm job IDs.
- A queued or running job is not a completed test. Use `sacct` to verify `COMPLETED` with exit code `0:0`.

## 1. Define an isolated test location and clone the repository

Run on Mantis:

```bash
TEST_ROOT=/projects/Enterica/Pedro/paper/plseq_reproducibility_test
REPO="$TEST_ROOT/PLSeq-reanalysis"
PROJECT="$TEST_ROOT/run"

mkdir -p "$TEST_ROOT"
cd "$TEST_ROOT"

git clone https://github.com/UConnMiura/PLSeq-reanalysis.git
cd "$REPO"

git rev-parse HEAD
git status --short
```

Checkpoint:

- Record the returned commit hash.
- The expected commit is `bb8724a7ed0709c0b9021074c19dd4fb792ec687`, or a documented newer commit.
- `git status --short` should be empty.
- If the test directory already contains a clone or run directory, stop and report it rather than overwriting it.

## 2. Download and verify the dm6 reference genome

The dm6 FASTA is deliberately not stored in GitHub or Zenodo. It is a public reference dependency obtained directly from UCSC and verified by checksum.

```bash
mkdir -p "$TEST_ROOT/reference_download"
cd "$TEST_ROOT/reference_download"

curl -L \
  https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz \
  -o dm6.fa.gz

printf '4a5777324403eff92c3650c30af30120  dm6.fa.gz\n' \
  | md5sum -c -

gzip -dc dm6.fa.gz > dm6.fa

test -s dm6.fa
head -1 dm6.fa
```

Checkpoint: `md5sum` must report `dm6.fa.gz: OK`, and the first FASTA line should begin with `>`.

## 3. Create the portable project and obtain the archived container

```bash
cd "$REPO"

export PLSEQ_HPC_PROJECT="$PROJECT"

bash scripts/hpc/00_setup_project.sh \
  --project "$PROJECT" \
  --dm6-fasta "$TEST_ROOT/reference_download/dm6.fa" \
  --download-container

export PLSEQ_APPTAINER_IMAGE="$PROJECT/software/pullalong-core.sif"
```

Checkpoint:

- Setup should create a 17-run analyzed-sample manifest.
- Confirm the main directories and files exist:

```bash
test -s "$PROJECT/metadata/run_manifest.tsv"
test -s "$PROJECT/references/dm6.fa"
test -s "$PLSEQ_APPTAINER_IMAGE"

wc -l "$PROJECT/metadata/run_manifest.tsv"
find "$PROJECT" -maxdepth 2 -type d | sort
```

The manifest should contain 18 lines: one header plus 17 analyzed runs.

## 4. Verify the exact Apptainer image

```bash
cd "$PROJECT/software"

printf '%s  %s\n' \
  '6037ffe045c08e7435e8a997b890d791f6f3554ab01cef16b0b895761c199beb' \
  'pullalong-core.sif' \
  | sha256sum -c -
```

Checkpoint: the command must report `pullalong-core.sif: OK`.

Load Apptainer and inspect the key versions:

```bash
module load apptainer/1.4.4-gcc-11.4.0-znzbvwm

apptainer exec "$PLSEQ_APPTAINER_IMAGE" \
  minimap2 --version

apptainer exec "$PLSEQ_APPTAINER_IMAGE" \
  samtools --version | head -2

apptainer exec "$PLSEQ_APPTAINER_IMAGE" \
  python3 -c 'import pysam; print("pysam", pysam.__version__)'
```

Expected versions:

```text
minimap2 2.17-r941
samtools 1.6
Using htslib 1.6
pysam 0.15.4
```

## 5. Restrict the pilot to SRR16318553

Preserve the complete public manifest, then make a one-sample working manifest:

```bash
cd "$PROJECT/metadata"

cp run_manifest.tsv run_manifest.full_17.tsv

{
  head -1 run_manifest.full_17.tsv
  awk -F '\t' '$1 == "SRR16318553"' run_manifest.full_17.tsv
} > run_manifest.tsv

cat run_manifest.tsv
```

Checkpoint:

```bash
test "$(wc -l < run_manifest.tsv)" -eq 2
grep -q '^SRR16318553[[:space:]]' run_manifest.tsv
```

## 6. Configure this Mantis session

Run these exports again after any logout or new shell:

```bash
cd "$REPO"

export PLSEQ_HPC_PROJECT="$PROJECT"
export PLSEQ_APPTAINER_IMAGE="$PROJECT/software/pullalong-core.sif"
export PLSEQ_APPTAINER_MODULE='apptainer/1.4.4-gcc-11.4.0-znzbvwm'
export PLSEQ_PARTITION=general
export PLSEQ_QOS=general
export PLSEQ_ARRAY_LIMIT=1
```

Optional verification:

```bash
printf 'Repository: %s\nProject: %s\nContainer: %s\n' \
  "$REPO" "$PLSEQ_HPC_PROJECT" "$PLSEQ_APPTAINER_IMAGE"
```

## 7. Submit and verify the download stage

```bash
cd "$REPO"
download_job=$(bash scripts/hpc/submit_level3.sh download)
echo "Download job: $download_job"
```

The wrapper prints its explanatory message to standard error and returns only
the numeric Slurm job ID on standard output, so command substitution stores a
value that can be passed directly to `squeue` and `sacct`. Monitor it with:

```bash
squeue -j "$download_job"
sacct -j "$download_job" \
  --format=JobID,State,Elapsed,MaxRSS,ExitCode
```

Do not continue until the relevant batch/array task reports `COMPLETED` and `0:0`.

Then verify:

```bash
ls -lh "$PROJECT/fastq/SRR16318553_1.fastq.gz"
awk -F '\t' '$1 == "SRR16318553" {print $4 "  ../fastq/" $1 "_1.fastq.gz"}' \
  "$PROJECT/metadata/run_manifest.tsv" \
  | (cd "$PROJECT/metadata" && md5sum -c -)
```

The deposited MD5 is authoritative. Do not also run `gzip -t`, which performs
a second complete read/decompression pass over the archive without adding a
meaningful integrity guarantee after an exact MD5 match.

If the job fails, inspect the corresponding files under `$PROJECT/logs/` and stop.

## 8. Submit and verify the alignment stage

Only after the download checkpoint passes:

```bash
cd "$REPO"
align_job=$(bash scripts/hpc/submit_level3.sh align)
echo "Alignment job: $align_job"
```

Monitor:

```bash
squeue -j "$align_job"
sacct -j "$align_job" \
  --format=JobID,State,Elapsed,MaxRSS,ExitCode
```

Do not continue until the relevant task reports `COMPLETED` and `0:0`.

Validate the BAM with the archived container:

```bash
apptainer exec \
  --bind "$PROJECT:$PROJECT" \
  "$PLSEQ_APPTAINER_IMAGE" \
  samtools quickcheck -v \
  "$PROJECT/bam/SRR16318553.bam"

test -s "$PROJECT/bam/SRR16318553.bam.bai"
```

`samtools quickcheck` should print no error.

## 9. Submit and verify filtering plus PSI reconstruction

Only after the alignment checkpoint passes:

```bash
cd "$REPO"
filter_job=$(bash scripts/hpc/submit_level3.sh filter)
echo "Filtering job: $filter_job"
```

Monitor:

```bash
squeue -j "$filter_job"
sacct -j "$filter_job" \
  --format=JobID,State,Elapsed,MaxRSS,ExitCode
```

Do not declare success until the task reports `COMPLETED` and `0:0`.

## 10. Verify final products

All of the following must exist and be non-empty:

```bash
for required in \
  "$PROJECT/bam/SRR16318553.bam" \
  "$PROJECT/bam/SRR16318553.bam.bai" \
  "$PROJECT/filtered_bam/SRR16318553/long.bam" \
  "$PROJECT/filtered_bam/SRR16318553/short.bam" \
  "$PROJECT/psi/SRR16318553/all.psi.csv" \
  "$PROJECT/psi/SRR16318553/long.psi.csv" \
  "$PROJECT/psi/SRR16318553/short.psi.csv" \
  "$PROJECT/psi/SRR16318553/filter_stage_counts.tsv" \
  "$PROJECT/psi/SRR16318553/SHA256SUMS" \
  "$PROJECT/psi/SRR16318553/PIPELINE_COMPLETE"
do
  test -s "$required" || {
    echo "ERROR: missing or empty output: $required" >&2
    exit 1
  }
done
```

Validate the final BAMs:

```bash
apptainer exec \
  --bind "$PROJECT:$PROJECT" \
  "$PLSEQ_APPTAINER_IMAGE" \
  samtools quickcheck -v \
  "$PROJECT/filtered_bam/SRR16318553/long.bam" \
  "$PROJECT/filtered_bam/SRR16318553/short.bam"
```

Inspect the reconstructed results:

```bash
echo '===== stage counts ====='
column -t -s $'\t' \
  "$PROJECT/psi/SRR16318553/filter_stage_counts.tsv"

echo '===== all molecules ====='
head -10 "$PROJECT/psi/SRR16318553/all.psi.csv"

echo '===== long 3-prime UTR molecules ====='
head -10 "$PROJECT/psi/SRR16318553/long.psi.csv"

echo '===== short 3-prime UTR molecules ====='
head -10 "$PROJECT/psi/SRR16318553/short.psi.csv"

echo '===== completion marker ====='
cat "$PROJECT/psi/SRR16318553/PIPELINE_COMPLETE"
```

If `SHA256SUMS` contains paths relative to its directory, verify them as follows:

```bash
cd "$PROJECT/psi/SRR16318553"
sha256sum -c SHA256SUMS
```

## 11. Final report

Report:

- Git commit tested.
- dm6 checksum result.
- SIF checksum result.
- Reported minimap2, SAMtools/HTSlib, and pysam versions.
- Download, alignment, and filtering Slurm job IDs and final `sacct` states.
- Whether every expected output exists.
- Whether all three BAMs pass `samtools quickcheck`.
- The contents of `filter_stage_counts.tsv`.
- The SRR16318553 rows from `all.psi.csv`, `long.psi.csv`, and `short.psi.csv` for the events of interest.
- Any warnings, deviations, or manual interventions.

Success means that a clean public clone, the checksum-verified dm6 reference, and the checksum-verified archived SIF reproduce the expected BAM and PSI product structure without reading from the original analysis directory.

## Scope of the reproducibility claim

This test establishes **execution reproducibility with the exact archived container**. The original Apptainer definition file was not retained, so rebuilding a byte-identical SIF from a recipe is not currently claimed. The SIF itself, its checksum, software-version record, source data, workflow scripts, and public reference-genome checksum together provide the executable environment needed for this analysis.
