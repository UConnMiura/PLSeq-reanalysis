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

## Rebuildable plotting image

`plseq-plotting.def` defines a separate plotting and document environment. It
contains Python 3.12, Matplotlib, NumPy, SciPy, Pillow, pandas, pysam,
pyGenomeTracks, R with ggplot2/dplyr/tidyr, Pandoc, and Poppler. It does not
replace `pullalong-core.sif`; alignment and PSI reconstruction remain tied to
the archived historical image described above.

Build the plotting image on an AMD64 Linux system with Apptainer:

```bash
apptainer build --fakeroot plseq-plotting.sif \
  environment/container/plseq-plotting.def
```

The definition includes an embedded `%test` section. The repository also
provides an artifact-level smoke test that generates and validates both PNG
and PDF output:

```bash
bash tests/test_plotting_container.sh \
  plseq-plotting.sif \
  /tmp/plseq-plotting-test
```

The image binary is intentionally excluded by `.gitignore`; Git tracks the
definition and tests, not the 1.2-GB SIF. The validated Mantis build has the
following SHA-256 checksum, also recorded in `plseq-plotting.sif.sha256`:

```text
6b7ad4435e3e7e1b62bf63c921f4fb6b0f2adfac101bf12a4bf1edb5264aeb4b
```

Record the checksum of a new release build with:

```bash
sha256sum plseq-plotting.sif > plseq-plotting.sif.sha256
```
