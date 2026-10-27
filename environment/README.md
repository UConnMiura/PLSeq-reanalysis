# Software environments

## Local analysis and figure environment

`environment.yml` records the principal Python and R packages used for the
publication figures. Create it with:

```bash
conda env create -f environment/environment.yml
conda activate plseq-reanalysis
```

The manuscript reader PDF additionally requires XeLaTeX. Quarto may be used for
multi-format manuscript output; the reader PDF build itself invokes Pandoc and
XeLaTeX through `manuscript/build_reader_pdf.sh`.

## HPC processing environment

FASTQ alignment and BAM filtering were executed inside the pinned Apptainer
image `pullalong-core.sif`. The image used minimap2 2.17-r941, samtools 1.6,
pysam 0.15.4, and featureCounts 2.0.6. HPC scripts load Apptainer 1.4.4 and
require the image path through `PLSEQ_APPTAINER_IMAGE`, or default to
`$PLSEQ_HPC_PROJECT/software/pullalong-core.sif`.

The original image definition was not recovered into the local publication
workspace and must not be reconstructed from memory. Before public release,
either the exact definition should be recovered or the verified SIF should be
archived separately and its SHA-256 checksum recorded here.

