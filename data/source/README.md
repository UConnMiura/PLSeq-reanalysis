# Level 3 source inputs

This directory contains the small programs and interval files used by the
FASTQ-level reconstruction. `scripts/hpc/00_setup_project.sh` copies them to
the runtime names expected by the Slurm workflow.

| Repository file | Runtime role |
|---|---|
| `Fly_dm6.bed` | dm6 BED12 junction annotation supplied to minimap2 |
| `up.figure3e.reconstructed.bed` | required upstream flanking exons |
| `down.figure3e.reconstructed.bed` | required downstream flanking exons |
| `UTR.figure3e.bed` | common/universal 3'UTR intervals |
| `dis.figure3e.bed` | distal/extended 3'UTR intervals |
| `SE_PSI.figure3e.suppdata9.bed` | cassette-exon intervals reconstructed from Supplementary Data 9 |
| `exon_covering.py` | flanking-exon and UTR coverage filter |
| `polyA_filtering.py` | short-3'UTR poly(A) filter |
| `calculate_PSI.py` | cassette-exon inclusion and PSI calculation |

The three Python programs are the copies used during this reanalysis. They were
obtained with the PL-Seq analysis materials, but no independent upstream
release tag or commit identifier was available. Their identities are therefore
defined by the repository commit and `SHA256SUMS`, rather than by an inferred
software version.

The upstream/downstream flanking intervals were reconstructed from the retained
dm6 transcript annotation and reviewed for the 28 Figure 3E events. The
cassette-exon intervals were derived from Supplementary Data 9. The common and
distal UTR intervals are the retained PL-Seq event intervals used in the
reconstruction. These are targeted event definitions, not a general dm6 gene
annotation package.

Verify all files from this directory with:

```bash
cd data/source
sha256sum -c SHA256SUMS
```
