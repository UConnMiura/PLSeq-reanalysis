# Repository inventory and provenance boundary

This document records the boundary between files retained in Git and files that
must be obtained or regenerated elsewhere.

## Included in Git

- Compact event-defining BED files for the 28 cassette-exon/3′UTR pairs.
- The dm6 BED12 transcript annotation used for representative gene models.
- Reconstructed count, PSI, statistical, and audit tables.
- Python and R analysis programs.
- SLURM templates for the Mantis HPC workflow.
- Main and supplemental publication figures.
- Quarto manuscript source and bibliography.
- Repository safety and validation tests.

## Excluded from Git

- Original and reconstructed FASTQ files.
- BAM, BAI, CRAM, CRAI, and SAM files.
- dm6 genome FASTA files and generated indexes.
- Apptainer SIF images.
- Temporary filtering products, caches, and logs.
- Generated Quarto and LaTeX build directories.

## Current external working locations

The paths below document the environment used during development; they are not
required paths for future users.

- Mantis project: `/projects/Enterica/Pedro/paper/plseq`
- Exact analysis container: archived as `pullalong-core.sif` in Zenodo record
  <https://doi.org/10.5281/zenodo.23769044>; the binary itself remains excluded
  from GitHub.
- Original local working analysis: `PMID37679364_reanalysis`

Portable workflow configuration will use `config/paths.yml`, created locally
from `config/paths.example.yml`. The local configuration is ignored by Git.

## Publication release plan

1. Tag a tested GitHub release.
2. Archive that release with Zenodo and add the resulting DOI to the manuscript.
3. Cite original sequencing reads by public run accession.
4. Deposit large derived files separately only if they add value beyond the
   documented regeneration workflow.
