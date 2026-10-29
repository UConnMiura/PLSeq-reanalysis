# Archived Apptainer image

The exact AMD64 Apptainer image used for FASTQ alignment and BAM/PSI
reconstruction is archived in Zenodo v3:

- Record: <https://doi.org/10.5281/zenodo.23769044>
- File: `pullalong-core.sif`
- SHA-256: `6037ffe045c08e7435e8a997b890d791f6f3554ab01cef16b0b895761c199beb`

Download and verify it with:

```bash
curl -L 'https://zenodo.org/records/23769044/files/pullalong-core.sif?download=1' \
  -o pullalong-core.sif
sha256sum -c pullalong-core.sif.sha256
```

The original definition file was not retained. The archived SIF, checksum,
embedded labels, and command-derived version inventory are therefore the
authoritative record of the executed environment. The binary is reproducibly
identifiable and executable, but it is not independently rebuildable from a
definition file contained in this repository.

Third-party programs bundled in the image retain their respective licenses.
