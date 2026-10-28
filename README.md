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
make manuscript
```

The reviewed PNG outputs have recorded SHA-256 checksums. Automated tests
verify that regenerated figures match the reviewed versions.

### Level 3: Reconstruct the analysis from FASTQ files

The complete reconstruction begins with the public sequencing accessions in
`config/accessions.tsv`. The SLURM workflow in `scripts/hpc/` performs:

1. FASTQ download and checksum validation.
2. Spliced alignment to dm6 using minimap2.
3. BAM sorting, indexing, and quality control.
4. Filtering through required flanking exons.
5. Classification by common and distal 3′UTR coverage.
6. Poly(A)-based filtering of short-3′UTR candidates.
7. Long- and short-3′UTR PSI calculation.
8. Preparation of compact plotting inputs.

The HPC processing steps were run through a pinned Apptainer environment. The
binary container is not stored in GitHub.

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
