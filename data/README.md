# Data organization

## `source/`

Small files that define the analyzed cassette exons, constitutive flanking
exons, common 3′UTR regions, and distal 3′UTR regions. These files are retained
because they define the biological events and are required to reproduce the
filtering decisions.

The dm6 genome FASTA is deliberately excluded. Its source, version, checksum,
and retrieval command will be recorded in the data-provenance document.

## `derived/`

Compact count tables, PSI tables, statistical results, and audit summaries
derived from the accessioned reads. These permit independent checking and
figure regeneration without downloading and aligning the complete FASTQ set.

## Files stored elsewhere

| File class | Location | Reason |
|---|---|---|
| Original FASTQ | SRA/ENA accessions | Already publicly archived |
| Reconstructed BAM/BAI | Regenerated on HPC; optional archival deposit | Too large for ordinary Git |
| dm6 FASTA and indexes | Authoritative genome source | Publicly downloadable |
| Apptainer SIF | Optional release/Zenodo asset | Binary and comparatively large |
| Temporary filtered BAMs | Not archived by default | Regenerable intermediate |

