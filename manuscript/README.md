# Manuscript source

This directory contains the reproducible Quarto manuscript draft.

The rendered HTML manuscript is published automatically at
<https://uconnmiura.github.io/PLSeq-reanalysis/>. GitHub Actions rebuilds the
site whenever manuscript sources or publication figures change on `main`.

The draft currently defines five main figures:

1. workflow, analytical improvements, and project timing;
2. replicate-resolved long- versus short-3′UTR PSI across all 28 genes;
3. selected read-level locus views;
4. Elav-dependent genotype-by-3′UTR analysis;
5. read-level control-versus-*elav* views for *Khc-73* and *Dscam1*.

The manuscript uses the pinned `nature.csl` style, which produces the numbered
Nature reference format requested by Nature Communications.

Render from this directory with:

```bash
quarto render
```

To build the polished reader/preprint PDF with a journal-style title page and
each figure plus its legend enclosed in one light border, run:

```bash
bash build_reader_pdf.sh
```

Items intentionally left for completion before submission include the
repository URL, archived release DOI, complete accession table, exact software
version table, and final supplementary numbering.
