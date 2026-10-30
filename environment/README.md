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
XeLaTeX through `manuscript/build_reader_pdf.sh`. The default PDF font is
Libertinus Sans and can be overridden with `PLSEQ_MAINFONT`. Poppler is included
because the figure preflight check calls `pdfinfo`.

## HPC processing environment

FASTQ alignment and BAM filtering were executed inside the pinned Apptainer
image `pullalong-core.sif`. The image used minimap2 2.17-r941, samtools 1.6,
pysam 0.15.4, and featureCounts 2.0.6. HPC scripts load Apptainer 1.4.4 and
require the image path through `PLSEQ_APPTAINER_IMAGE`, or default to
`$PLSEQ_HPC_PROJECT/software/pullalong-core.sif`.

The exact image is archived at <https://doi.org/10.5281/zenodo.23769044>.
Checksum, inspection, and version records are tracked in `environment/container/`.
The original image definition was not retained, so the archived SIF is the
authoritative executed environment and must not be described as rebuildable
from a definition file.

Figure generation can alternatively use the rebuildable
`environment/container/plseq-plotting.def` image. This newer image supplies
Matplotlib, pyGenomeTracks, R plotting packages, Pandoc, and Poppler while
keeping the historical alignment/filtering container unchanged. Build and
test instructions are in `environment/container/README.md`.
