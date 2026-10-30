# Reproducible reanalysis of targeted long-read PL-Seq data

[![Repository safety](https://github.com/UConnMiura/PLSeq-reanalysis/actions/workflows/repository-safety.yml/badge.svg)](https://github.com/UConnMiura/PLSeq-reanalysis/actions/workflows/repository-safety.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23712455.svg)](https://doi.org/10.5281/zenodo.23712455)

**[Read the formatted manuscript](https://uconnmiura.github.io/PLSeq-reanalysis/)** ·
**[Browse the main figures](https://github.com/UConnMiura/PLSeq-reanalysis/tree/main/figures/main)** ·
**[View the Zenodo archive](https://doi.org/10.5281/zenodo.23712455)**

This repository contains the data products, analysis code, publication figures,
and manuscript source for a reproducible reanalysis of the Pull-a-Long-Seq
study associated with [PMID 37679364](https://pubmed.ncbi.nlm.nih.gov/37679364/).

Pull-a-Long-Seq, or PL-Seq, uses oligonucleotide probes to enrich selected
full-length cDNAs before Oxford Nanopore sequencing. This targeted strategy
provides increased read depth for genes of interest while retaining connections
between alternative exons and alternative 3′ ends on individual sequenced
molecules.

The present work revisits the deposited PL-Seq data to improve statistical
inference, read-level visualization, computational reproducibility, and reuse.
It is intended primarily as a reproducible research record and biological
reanalysis rather than as a general-purpose software package.

## What this reanalysis adds

- Reconstruction of selected PL-Seq results from deposited FASTQ files.
- Explicit separation of reads assigned to long and short 3′UTR classes.
- Replicate-resolved presentation of cassette-exon inclusion.
- Benjamini–Hochberg correction across each defined family of tested events.
- Reporting of adjusted q values instead of nominal significance symbols alone.
- Read-level locus plots that preserve individual exon structures.
- Programmatic transcript annotations derived from dm6 reference models.
- Reframing of the Elav analysis as a genotype-by-3′UTR interaction.
- Fully scripted figure assembly without manual Illustrator composition.
- A Quarto manuscript connecting analyses, figures, legends, and references.

The reanalysis does not replace the original study. It provides a more
inspectable and reusable representation of the deposited data and highlights
additional biological patterns revealed by replicate-aware statistics and
read-level visualization.

## Main outputs

### Figure 1: Reanalysis workflow

Figure 1 summarizes the progression from the published paper and supplementary
tables through FASTQ reconstruction, containerized alignment, statistical
reanalysis, read-level visualization, and manuscript generation.

### Figure 2: Replicate-resolved exon–3′UTR relationships

Figure 2 shows long- and short-3′UTR PSI values for the 28 cassette-exon events
examined in 16–18-hour embryos. Individual biological replicates are displayed
rather than only group summaries. Benjamini–Hochberg-adjusted q values identify
events that remain significant after correction across the tested event family.

### Figure 3: Read-level examples

Figure 3 combines the quantitative analysis with representative read-level
views for selected genes. The panels show mean fractional coverage across
replicates, sampled primary alignments, cassette-exon inclusion or skipping,
common and distal 3′UTR regions, representative dm6 transcript models when
available, and the number of informative primary alignments.

### Figure 4: Elav-dependent exon–3′UTR coordination

Figure 4 reanalyzes the control-versus-*elav* mutant experiment. The principal
estimand is the genotype-by-3′UTR interaction:

```text
(PSI_long − PSI_short)_elav5 − (PSI_long − PSI_short)_control
```

This asks whether loss of Elav changes the relationship between cassette-exon
inclusion and 3′UTR identity, rather than testing only pooled PSI or overall
distal polyadenylation-site usage. The figure also distinguishes interaction
effects from marginal changes in PSI and distal polyadenylation-site usage.

### Figure 5: Read-level validation in control and *elav* mutant samples

Figure 5 shows control and *elav* mutant reads for *Khc-73* and *Dscam1*,
separated into long- and short-3′UTR classes. These panels provide a direct view
of the read structures underlying the interaction effects shown in Figure 4.

### Supplemental locus atlas

The supplemental atlas extends the read-level presentation across all 28
cassette-exon events, providing a systematic visual catalog rather than relying
only on selected examples.

## Three levels of reproducibility

### Level 1: Inspect the results

No sequencing files or specialized environment are required. Compact tables in
`data/derived/` contain reconstructed counts, PSI values, statistical results,
and audit summaries. Final figures are in `figures/`, and the manuscript source
is in `manuscript/`.

### Level 2: Rebuild the figures and manuscript locally

Compact sampled-read BED files and fractional-coverage tracks permit the
principal figures to be rebuilt without downloading FASTQ files or storing BAM
files locally.

```bash
conda env create -f environment/environment.yml
conda activate plseq-reanalysis
bash scripts/install_git_hooks.sh

make check
make figures-core
make manuscript-pdf
```

The reviewed PNG outputs have recorded SHA-256 checksums. Automated tests
verify that regenerated figures match the reviewed versions. `poppler` is part
of the environment because the preflight step calls `pdfinfo`; XeLaTeX must be
available separately for PDF creation. `make manuscript-pdf` builds only the
PDF, whereas `make manuscript` intentionally regenerates the core figures first.

The reviewed figures were created with Arial/Helvetica metrics. On systems
without either font, Matplotlib may fall back to DejaVu Sans and produce a
visually different raster, causing the strict checksum test to fail. The
tracked publication figures remain the reviewed reference outputs while a
fully bundled open-font revision is developed. Figure-generation commands
overwrite the corresponding tracked files, so run them in a clean branch or
worktree when comparing environments.

For systems with Apptainer, the repository also includes a rebuildable plotting
container definition and an executable smoke test:

```bash
apptainer build --fakeroot plseq-plotting.sif \
  environment/container/plseq-plotting.def
bash tests/test_plotting_container.sh \
  plseq-plotting.sif /tmp/plseq-plotting-test
```

This plotting image is separate from the historical `pullalong-core.sif` used
for Level 3 alignment and filtering.

### Level 3: Reconstruct the analysis from FASTQ files

The complete reconstruction defaults to the 17 embryo FASTQ runs actually used
in the reanalysis. `config/accessions.tsv` is the broader 49-run source-study
inventory and is not used as the default Level 3 run manifest. The Slurm
workflow in `scripts/hpc/` performs:

1. FASTQ download and checksum validation.
2. Spliced alignment to dm6 using minimap2.
3. BAM sorting, indexing, and quality control.
4. Filtering through required flanking exons.
5. Classification by common and distal 3′UTR coverage.
6. Poly(A)-based filtering of short-3′UTR candidates.
7. Long- and short-3′UTR PSI calculation.
8. Preparation of compact plotting inputs.

The exact AMD64 Apptainer image is archived in Zenodo v3 at
[doi:10.5281/zenodo.23769044](https://doi.org/10.5281/zenodo.23769044). It is
not stored in GitHub. The original definition file was not retained; the SIF,
checksum, embedded labels, and version report are the authoritative record of
the executed environment.

### Level 3 setup

Clone the repository on the cluster and obtain the UCSC dm6 FASTA:

```bash
curl -L https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz \
  -o dm6.fa.gz
printf '4a5777324403eff92c3650c30af30120  dm6.fa.gz\n' | md5sum -c -
gzip -dc dm6.fa.gz > dm6.fa
```

Create the complete project layout and download the verified container:

```bash
export PLSEQ_HPC_PROJECT=/path/to/hpc/project

bash scripts/hpc/00_setup_project.sh \
  --project "$PLSEQ_HPC_PROJECT" \
  --dm6-fasta "$PWD/dm6.fa" \
  --download-container

export PLSEQ_APPTAINER_IMAGE="$PLSEQ_HPC_PROJECT/software/pullalong-core.sif"
```

The setup script creates the expected reference, BED, filtering-program,
metadata, output, and log directories. It maps the small reviewed source files
from `data/source/` into their runtime locations and builds the 17-run manifest
from Supplementary Table 1.

Submit each independently restartable stage after reviewing its predecessor:

```bash
bash scripts/hpc/submit_level3.sh download
bash scripts/hpc/submit_level3.sh align
bash scripts/hpc/submit_level3.sh filter
```

The wrapper creates `logs/`, derives the array range from the manifest, and
defaults to four simultaneous tasks. Optional site settings are passed through
environment variables rather than embedded personal or cluster-specific values:

```bash
export PLSEQ_ARRAY_LIMIT=4
export PLSEQ_PARTITION=general       # omit when the cluster default is correct
export PLSEQ_QOS=general             # omit when no QOS is required
export PLSEQ_APPTAINER_MODULE='apptainer/1.4.4-gcc-11.4.0-znzbvwm'
export PLSEQ_MAIL_USER='your.name@example.edu'  # optional; no email by default
```

Additional `sbatch` options may be supplied through `PLSEQ_SBATCH_ARGS`.
Because the three stages are submitted separately, successfully downloaded or
aligned samples can be processed without waiting indefinitely for one failed
array task; rerunning a stage reuses FASTQs that already match their deposited
MD5 checksum. The workflow does not additionally run `gzip -t`: both operations
read the complete archive, and the authoritative deposited MD5 is sufficient.

### Validated optimized filtering

The published three-program filtering route remains the reference workflow and
is selected with `filter`. The validated optimized implementation preserves
the filtering order, thresholds, and legacy PSI arithmetic while eliminating
repeated interval fetches, temporary CSV files, quadratic list deletion, and
unused failure BAMs. It writes to `optimized/<run>/`, leaving reference outputs
untouched:

```bash
bash scripts/hpc/submit_level3.sh filter-optimized
```

The wrapper runs one sample per Slurm array task. Set the maximum number of
samples processed concurrently before submission; four is the default and is a
conservative starting point on Mantis:

```bash
export PLSEQ_ARRAY_LIMIT=4
filter_job=$(bash scripts/hpc/submit_level3.sh filter-optimized)
sacct -j "$filter_job" --format=JobID,State,Elapsed,MaxRSS,ExitCode
```

This is sample-level parallelism: each sample remains independent and receives
the four CPUs requested by the batch script. Increase the array limit only when
the number of concurrent jobs and aggregate filesystem load are appropriate
for the cluster allocation.

The same program can run locally when Python and pysam are available:

```bash
bash scripts/local/run_filter_and_psi_optimized.sh \
  --input-bam /path/to/input.bam \
  --beds-dir /path/to/project/references/beds \
  --output-dir /path/to/optimized-output
```

Before treating optimized products as interchangeable, compare them with the
legacy final BAMs and PSI tables using
`scripts/optimized/compare_implementations.py`. This comparison uses complete
SAM alignment records, including multiplicity, rather than filenames or BAM
file checksums. Complete local and archived-container Mantis validations,
including exact output comparisons and resource use, are recorded in
[`docs/optimized_filter_validation.md`](docs/optimized_filter_validation.md).

After filtering is complete, the same wrapper can submit the optional compact
pileup-input stages with consistent site settings and absolute log paths:

```bash
bash scripts/hpc/submit_level3.sh figure3e-pileups
bash scripts/hpc/submit_level3.sh figure5-supplement
```

## Repository organization

```text
PLSeq-reanalysis/
├── config/             # Accessions, event manifests, style and checksums
├── data/
│   ├── source/         # Compact BED annotations and filtering programs
│   ├── derived/        # Counts, PSI values, statistics and audits
│   └── ...             # Compact read-level plotting inputs
├── docs/               # Provenance, AI use and release validation
├── environment/        # Local environment and container documentation
├── figures/
│   ├── main/
│   └── supplementary/
├── manuscript/         # Quarto source, bibliography and PDF support
├── scripts/
│   ├── hpc/            # Apptainer and SLURM processing workflow
│   └── local/          # Statistics and figure generation
└── tests/               # Repository and figure-integrity checks
```

## Data intentionally excluded from GitHub

The repository does not contain FASTQ, BAM, BAI, CRAM, CRAI, SAM, dm6 genome
FASTA, or Apptainer SIF files. It also excludes temporary filtering products,
logs, and Quarto build caches.

Original sequencing reads remain available through their public accessions.
Genome references and large derived files can be downloaded or regenerated
using the documented workflow. A pre-commit safety check and GitHub Actions
workflow reject prohibited sequencing files, genome FASTA files, container
images, and unexpectedly large files.

## Software environments

Local statistical analysis and figure generation use
`environment/environment.yml`.

FASTQ processing and BAM filtering were performed using Apptainer. The verified
analysis environment contained:

- minimap2 2.17-r941
- samtools 1.6
- pysam 0.15.4
- featureCounts 2.0.6

The HPC scripts require an explicitly defined project and container:

```bash
export PLSEQ_HPC_PROJECT=/path/to/hpc/project
export PLSEQ_APPTAINER_IMAGE="$PLSEQ_HPC_PROJECT/software/pullalong-core.sif"
```

Container checksum, inspection labels, and recorded versions are in
`environment/container/`. Site-specific Slurm partition, QOS, module, email,
and array-concurrency settings are deliberately not hardcoded into the batch
scripts.

## Validation

The repository has been tested from a fresh clone outside the development
directory. Four principal PNG figures regenerated byte-for-byte from the
included compact inputs. Reviewed checksums are stored in
`config/figure_checksums.sha256`.

Validation can be repeated with:

```bash
make check
```

GitHub Actions performs the same repository-safety, checksum, Python-syntax,
and shell-syntax checks after each push.

## Use of AI-assisted research software

ChatGPT Codex assisted with workflow planning, code drafting and revision, HPC
troubleshooting, statistical reframing, programmatic figure design, manuscript
organization, and reproducibility documentation. Paper2Agent-derived resources
were used to make the original paper and supplementary information accessible
for structured review.

All biological interpretations, statistical choices, figure designs, and
manuscript text remain subject to author review. AI-generated code is retained
as version-controlled source and evaluated against deposited tables,
reconstructed counts, and explicit validation checks. Additional details are
provided in `docs/AI-use.md`.

## Citation

Versioned GitHub releases are archived by Zenodo. For citations intended to
resolve to the latest archived version, use the concept DOI
[doi:10.5281/zenodo.23712455](https://doi.org/10.5281/zenodo.23712455).
Citation metadata is also provided in `CITATION.cff`.

## License

The analysis software and repository documentation are released under the MIT
License. Third-party publications, sequencing data, genome annotations, and
software retain their original licenses and terms of use.
