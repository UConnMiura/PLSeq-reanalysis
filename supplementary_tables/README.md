# Supplementary tables

These tab-separated files provide the machine-readable source data and
analysis results accompanying the PL-Seq reanalysis manuscript. Each file has
one header row, one observation per row, and uses `NA` for unavailable values.
Genomic coordinates are zero-based and half-open unless explicitly stated
otherwise.

The tables are rebuilt with:

```bash
node scripts/build_supplementary_tables.mjs
```

## Tables

### Supplementary Table 1: samples

`Supplementary_Table_1_samples.tsv` contains one row for each of the 17 public
FASTQ runs analyzed in this study. It provides SRA, BioSample, and BioProject
identifiers; experimental group, genotype, probe set, and biological replicate;
and FASTQ location, checksum, and file size.

### Supplementary Table 2: event definitions

`Supplementary_Table_2_event_definitions.tsv` contains one row for each of 28
cassette-exon events. It reports the published cassette-exon label and the dm6
coordinates used for the cassette exon, flanking exons, common 3'UTR, distal
3'UTR, and plotted locus window. `event_id` combines the gene with the
published CE label. `bed_target_id` retains the identifier used by the interval
and counting workflow because historical BED identifiers do not always use the
published exon-numbering convention.

### Supplementary Table 3: PL-Seq counts

`Supplementary_Table_3_plseq_counts.tsv` contains one row per run, event, and
molecule class. The `analysis_context` field distinguishes the all-gene p1/p2
embryo reconstruction from the p3 read-level reconstruction of *Khc-73* and
*Dscam1*. `inclusion_reads`, `skipping_reads`, `informative_reads`, and `psi`
report classifier-derived values. The `expected_*` fields retain the
corresponding `calculate_PSI.py` outputs used for reconciliation.

### Supplementary Table 4: embryo long-versus-short statistics

`Supplementary_Table_4_embryo_long_short_statistics.tsv` contains one row for
each of 28 events used in Figures 2 and 3. It reports replicate count, long- and
short-3'UTR mean PSI and SEM, the paired effect and 95% confidence interval,
nominal P value, and BH-adjusted q value across the 28-event family.

### Supplementary Table 5: *elav* interaction statistics

`Supplementary_Table_5_elav_interaction.tsv` contains one row for each of 22
events with estimable control and *elav* conditional PSI values. The interaction
is defined as `(long PSI - short PSI)_elav5 - (long PSI - short PSI)_control`.
The table reports the four conditional means, interaction estimate and 95%
confidence interval, nominal P value, and BH-adjusted q value.

### Supplementary Table 6: *elav* marginal statistics

`Supplementary_Table_6_elav_marginal_statistics.tsv` contains two rows per
event: one for QAPA-derived change in distal polyadenylation usage and one for
change in pooled cassette-exon PSI. BH correction is performed separately for
the two 22-test metric families.

### Supplementary Table 7: reconstruction audit

`Supplementary_Table_7_reconstruction_audit.tsv` compares deposited and
FASTQ-reconstructed inclusion counts, informative counts, and PSI values for
the p1/p2 embryo analysis. Differences, minimum-depth status, and the final
reconciliation status are reported explicitly.

## Shared identifiers and missing values

- `run_accession` joins observations to Supplementary Table 1.
- `event_id` and `gene` join event-level records to Supplementary Table 2.
- `biological_replicate` is numbered within each experimental group.
- `molecule_class` is `all`, `long_3utr`, or `short_3utr`.
- `NA` denotes a value that is unavailable or not applicable; it is not zero.
- Boolean fields use `TRUE` and `FALSE`.

## Provenance

The tables are generated from the versioned manifests and derived results in
`config/` and `data/derived/`. Original sequence reads remain in SRA/ENA and
are not included in this repository.
