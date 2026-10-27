# Release validation

## Clean-clone reconstruction test

On 2026-10-27, commit `42fefdd` was cloned into a newly created temporary
directory outside the development checkout. The `figures-core` target was run
with the pinned local publication Python environment.

The following outputs were regenerated successfully:

- Figure 1 workflow and timeline.
- Figure 2 replicate-resolved long- versus short-3′UTR PSI.
- Figure 3 read-level locus views.
- Figure 5 control-versus-*elav5* read-level views.
- The Figure 3 journal-art preflight report.

All four PNG files were byte-for-byte identical to the reviewed files in the
repository. The expected SHA-256 values are recorded in
`config/figure_checksums.sha256` and tested by
`tests/verify_figure_checksums.py`.

The regenerated PDF files differed at the binary level because Matplotlib
embeds creation timestamps in PDF metadata. This does not indicate a graphical
difference; the PNG equality confirms deterministic rendered content for the
reviewed outputs.

The clean clone also passed the repository safety audit and contained no FASTQ,
BAM, genome FASTA, or Apptainer SIF files.

