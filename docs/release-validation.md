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

## Public archive

GitHub release `v0.1.1` was archived by Zenodo on 2026-10-27.

- Version DOI: `10.5281/zenodo.23712456`
- Concept DOI: `10.5281/zenodo.23712455`
- Record: <https://zenodo.org/records/23712456>

Zenodo version v3 added the exact 155 MB AMD64 Apptainer image used for the
FASTQ-level reconstruction, together with its SHA-256 checksum, embedded image
labels, command-derived software versions, and a container-specific README.

- Container archive DOI: `10.5281/zenodo.23769044`
- Concept DOI: `10.5281/zenodo.23712455`
- SIF SHA-256: `6037ffe045c08e7435e8a997b890d791f6f3554ab01cef16b0b895761c199beb`
- Record: <https://zenodo.org/records/23769044>

The original image definition file was not retained. Verification therefore
establishes identity with the executed binary environment, not independent
reconstruction of that binary from a definition recipe.
