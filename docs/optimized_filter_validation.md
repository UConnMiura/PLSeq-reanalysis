# Optimized filtering validation

Status: experimental, local validation only  
Validation date: 2026-10-30

The optimized filtering implementation is additive. The archived
`exon_covering.py`, `polyA_filtering.py`, and `calculate_PSI.py` programs remain
the publication-reference workflow.

## What was changed

The optimized implementation preserves the original stage order and numerical
rules but removes several sources of computational overhead:

- each BED-overlapping alignment is evaluated once per filtering stage rather
  than fetched three times;
- block coordinates are handled in memory instead of being written to and read
  from temporary CSV files;
- failed reads are not repeatedly removed from Python lists;
- unused failure BAMs are not sorted and indexed;
- poly(A) candidates are streamed rather than retained as three full lists;
- intermediate work can be placed on node-local scratch on Slurm systems.

The implementation deliberately reproduces a shared-list behavior in the
archived `calculate_PSI.py`. Although that arithmetic is not the natural
per-read formulation, retaining it is required to reconstruct the published
PSI values exactly. A scientifically corrected PSI implementation, if added,
must be labeled separately and must not silently replace the reconstruction.

## Local equivalence test

The non-empty Khc-73 demonstration BAM distributed with the Pull-a-Long-Seq
software was processed with the same 28-event BED package by both the archived
and optimized implementations. Both paths used the same local Python/pysam and
SAMtools environment for this comparison.

| Stage | Legacy records | Optimized records |
|---|---:|---:|
| Raw | 5,297 | 5,297 |
| Upstream-exon pass | 443 | 443 |
| Downstream-exon pass | 413 | 413 |
| Universal-UTR pass | 397 | 397 |
| Distal overlap | 261 | 261 |
| Short candidates | 136 | 136 |
| Long final | 127 | 127 |
| Short final | 99 | 99 |

The validation then compared complete SAM record strings, including duplicate
multiplicity, for the final long and short BAMs. It also compared the complete
CSV fields for all 28 events in each of the all-, long-, and short-molecule PSI
tables. All comparisons were exact.

Observed wall-clock time on this small fixture was approximately 19 seconds for
the archived workflow and 0.8 seconds for the optimized workflow. This is a
useful development benchmark, not a claim about full-dataset speedup.

An additional local SRR16474418 pilot produced exact final BAM and PSI
equivalence, but its current 28-event filtering path yielded no final long or
short molecules and is therefore a weaker test than the non-empty demo.

## Required next validation

Before the optimized implementation is used as a publication-equivalent
replacement, it should be tested on Mantis with the archived Apptainer SIF and
at least one complete embryo BAM. Required comparisons are:

1. every stage count;
2. complete long-BAM alignment records and multiplicity;
3. complete short-BAM alignment records and multiplicity;
4. all-, long-, and short-molecule PSI tables;
5. runtime and peak memory from Slurm accounting.

The optimized Slurm route writes to `optimized/<run>/` and therefore cannot
overwrite the legacy outputs under `filtered_bam/` and `psi/`.
