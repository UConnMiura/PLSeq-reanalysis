# Reproducible reanalysis of targeted long-read PL-Seq data

This repository contains the code, compact analysis inputs, statistical outputs,
figure sources, and Quarto manuscript for a reproducible reanalysis of the
PL-Seq study associated with PMID 37679364.

The project re-examines coordination between cassette-exon inclusion and
alternative 3′UTR usage in targeted Oxford Nanopore cDNA sequencing. It adds
replicate-resolved visualization, Benjamini–Hochberg multiple-testing control,
read-level locus displays, and a reproducible manuscript and figure build.

## Repository scope

The repository intentionally does **not** contain FASTQ, BAM, CRAM, reference
FASTA, or Apptainer SIF files. Original reads remain available from public
sequence archives under the accessions listed in `config/`. Large derived files
can be regenerated on an HPC system using the Apptainer-based workflow.

The principal layers are:

1. `scripts/hpc/`: FASTQ download, minimap2 alignment, molecule filtering, and
   PSI reconstruction on a SLURM cluster.
2. `data/source/`: compact event definitions and analysis helper programs.
3. `data/derived/`: compact reconstructed counts and statistical results.
4. `scripts/local/`: statistical analyses and deterministic figure generation.
5. `manuscript/`: the Quarto manuscript, references, and PDF build support.

## Safety against accidental large-file commits

The `.gitignore` blocks sequencing data, alignments, genome FASTA files,
container images, caches, and build products. In addition,
`tests/check_repository_files.py` inspects both the working tree and Git index
and fails if a forbidden file type or unexpectedly large file is present.

Run it before every commit:

```bash
python3 tests/check_repository_files.py
```

The repository will also install the same check as a local pre-commit hook via:

```bash
bash scripts/install_git_hooks.sh
```

## Current reproducibility status

This is the clean-room repository scaffold. Compact evidence tables, current
figure outputs, manuscript sources, and the scripts used during the analysis
have been collected. The next development step is to remove remaining
machine-specific absolute paths and connect the stages through one documented
configuration-driven build.

## Data and software availability

The original FASTQ files are not redistributed here. Accession mappings and
download provenance will be maintained in `config/accessions.tsv`. The
Apptainer definition and its SHA-256 checksum will be deposited in
`environment/`; the binary image may be archived separately with a permanent
release if needed.

## Citation

Citation metadata will be finalized in `CITATION.cff` when the manuscript title
and public repository URL are fixed. A versioned GitHub release should be
archived through Zenodo to obtain a permanent DOI.

