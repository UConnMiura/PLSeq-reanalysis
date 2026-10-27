options(stringsAsFactors = FALSE, scipen = 999)

root <- "/Users/pedromiura/Documents/codex/PMID37679364_reanalysis"
src <- file.path(
  "/Users/pedromiura/Documents/codex/PMID37679364_agent/dist",
  "pmid-37679364-agent/skill/pmid-37679364-paper/assets/supp_table"
)
out <- file.path(root, "figures", "reproduced")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

paper_col <- c(Long = "#79C5BB", Short = "#F7AA78",
               Control = "#71C69E", Mutant = "#F4B3C2",
               elav = "#E6BD72")

read_sheet <- function(book, sheet) {
  read.csv(file.path(src, sprintf("41467-2023-41207-moesm%d-esm-sheet-%d.csv",
                                  book, sheet)), check.names = FALSE)
}

sem <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x) / sqrt(length(x))
}

sig_label <- function(p) {
  ifelse(!is.finite(p), "n.s.",
         ifelse(p < 0.005, "***", ifelse(p < 0.01, "**",
         ifelse(p < 0.05, "*", "n.s."))))
}

open_devices <- function(stem, width, height) {
  pdf(file.path(out, paste0(stem, ".pdf")), width = width, height = height,
      useDingbats = FALSE)
}

save_png <- function(stem, width, height, draw) {
  png(file.path(out, paste0(stem, ".png")), width = width * 300,
      height = height * 300, res = 300)
  draw()
  dev.off()
}

grouped_bars <- function(summary, raw, group_col, group_levels, colors,
                         ylim = c(0, 1.08), ylab = "PSI", main = "",
                         label_col = "label", p_col = "p") {
  labels <- unique(summary[[label_col]])
  summary[[label_col]] <- factor(summary[[label_col]], levels = labels)
  summary[[group_col]] <- factor(summary[[group_col]], levels = group_levels)
  mat <- matrix(NA_real_, nrow = length(group_levels), ncol = length(labels),
                dimnames = list(group_levels, labels))
  err <- mat
  for (i in seq_len(nrow(summary))) {
    mat[as.character(summary[[group_col]][i]), as.character(summary[[label_col]][i])] <- summary$mean[i]
    err[as.character(summary[[group_col]][i]), as.character(summary[[label_col]][i])] <- summary$sem[i]
  }
  bp <- barplot(mat, beside = TRUE, col = colors, border = NA, ylim = ylim,
                ylab = ylab, main = main, las = 2, cex.names = 0.62,
                axis.lty = 0)
  arrows(bp, mat, bp, mat + err, angle = 90, code = 2, length = 0.025)
  for (j in seq_along(labels)) {
    for (i in seq_along(group_levels)) {
      z <- raw[raw[[label_col]] == labels[j] & raw[[group_col]] == group_levels[i], ]
      if (!nrow(z)) next
      points(rep(bp[i, j], nrow(z)) + seq(-0.05, 0.05, length.out = nrow(z)),
             z$value, pch = 21, bg = "white", col = "#555555", cex = 0.48)
    }
    p <- unique(summary[[p_col]][summary[[label_col]] == labels[j]])
    has_both <- sum(is.finite(mat[, j])) == length(group_levels)
    lab <- if (has_both) sig_label(p[1]) else "N/A"
    top <- suppressWarnings(max(mat[, j] + err[, j], na.rm = TRUE))
    if (!is.finite(top)) top <- 0
    text(mean(bp[, j]), min(ylim[2] - 0.015, top + 0.055), lab, cex = 0.62)
  }
  legend("topright", legend = group_levels, fill = colors, bty = "n", cex = 0.75)
}

# Figure 3e-f: rebuild from inclusion/total counts, not the workbook's summary cells.
figure3_data <- function(sheet, panel) {
  d <- read_sheet(13, sheet)
  ce_map <- read_sheet(12, 1)
  names(ce_map)[1:5] <- c("Chromosome", "CE_start", "CE_end", "Gene", "CE_number")
  ce_key <- paste(ce_map$Gene, ce_map$CE_start, ce_map$CE_end, sep = "|")
  d$PSI[d$Total_counts < 2] <- NA_real_
  d$key <- paste(d$Gene, d$CE_start, d$CE_end, sep = "|")
  keys <- unique(d$key)
  rows <- list()
  raw <- list()
  for (k in keys) {
    z <- d[d$key == k, ]
    gene <- z$Gene[1]
    ce <- ce_map$CE_number[match(k, ce_key)]
    if (!length(ce) || is.na(ce)) ce <- z$CE_start[1]
    w <- reshape(z[, c("Run", "Isoform", "PSI")], idvar = "Run",
                 timevar = "Isoform", direction = "wide")
    ok <- complete.cases(w$PSI.Long, w$PSI.Short)
    p <- if (sum(ok) >= 2) t.test(w$PSI.Long[ok], w$PSI.Short[ok], paired = TRUE)$p.value else NA_real_
    for (iso in c("Long", "Short")) {
      vals <- z$PSI[z$Isoform == iso]
      rows[[length(rows) + 1]] <- data.frame(
        label = paste0(gene, "\n", ce), Gene = gene, CE = ce, Isoform = iso,
        n = sum(is.finite(vals)), mean = mean(vals, na.rm = TRUE), sem = sem(vals), p = p
      )
      raw[[length(raw) + 1]] <- data.frame(
        label = paste0(gene, "\n", ce), Isoform = iso, value = vals
      )
    }
  }
  list(summary = do.call(rbind, rows), raw = do.call(rbind, raw), panel = panel)
}

f3e <- figure3_data(1, "16–18 hr embryos")
f3f <- figure3_data(2, "Adult heads")

draw_f3 <- function() {
  par(mfrow = c(2, 1), mar = c(7.7, 4.1, 2.2, 1), oma = c(0, 0, 1, 0))
  grouped_bars(f3e$summary, f3e$raw, "Isoform", c("Long", "Short"),
               paper_col[c("Long", "Short")], main = f3e$panel,
               ylab = "CE PSI")
  grouped_bars(f3f$summary, f3f$raw, "Isoform", c("Long", "Short"),
               paper_col[c("Long", "Short")], main = f3f$panel,
               ylab = "CE PSI")
  mtext("Figure 3e–f reproduction — Source Data, count-derived PSI", outer = TRUE,
        cex = 1.05, font = 2)
}
open_devices("figure3_ef_reproduction", 13, 8.5); draw_f3(); dev.off()
save_png("figure3_ef_reproduction", 13, 8.5, draw_f3)

# Figure 5a-b: PL-Seq replicate values in Source Data sheets 4 and 5.
read_fig5_psi <- function(sheet, gene, mutant_label) {
  d <- read_sheet(15, sheet)
  d <- d[, 1:4]
  names(d) <- c("Gene", "value", "Isoform", "Sample")
  d <- d[is.finite(d$value) & nzchar(d$Sample), ]
  d$Sample2 <- ifelse(d$Sample == "Control", "Control", "Mutant")
  d$label <- d$Isoform
  rows <- do.call(rbind, lapply(c("Short", "Long"), function(iso) {
    do.call(rbind, lapply(c("Control", "Mutant"), function(grp) {
      x <- d$value[d$Isoform == iso & d$Sample2 == grp]
      data.frame(label = iso, Sample2 = grp, n = length(x),
                 mean = if (length(x)) mean(x) else NA_real_, sem = sem(x), p = NA_real_)
    }))
  }))
  # The published comparison is control-short versus mutant-short.
  x <- d$value[d$Isoform == "Short" & d$Sample2 == "Control"]
  y <- d$value[d$Isoform == "Short" & d$Sample2 == "Mutant"]
  rows$p[rows$label == "Short"] <- t.test(x, y)$p.value
  list(summary = rows, raw = d[, c("label", "Sample2", "value")],
       title = paste0(gene, " (", mutant_label, ")"))
}

f5a <- read_fig5_psi(4, "Dscam1 exon 19", "long-3UTR deletion")
f5b <- read_fig5_psi(5, "Khc-73 exon 15", "long-3UTR deletion")

read_fig5_qpcr <- function(sheet, title) {
  d <- read_sheet(15, sheet)[, 1:2]
  names(d) <- c("Sample", "value")
  d <- d[is.finite(d$value) & nzchar(d$Sample), ]
  d$Sample <- gsub("△", "delta-", d$Sample, fixed = TRUE)
  d$Sample <- gsub(" ", "", d$Sample, fixed = TRUE)
  lev <- unique(d$Sample)
  p <- t.test(d$value[d$Sample == lev[1]], d$value[d$Sample == lev[2]])$p.value
  s <- do.call(rbind, lapply(lev, function(g) {
    x <- d$value[d$Sample == g]
    data.frame(label = "ext/uni", Sample2 = g, n = length(x), mean = mean(x), sem = sem(x), p = p)
  }))
  d$label <- "ext/uni"; d$Sample2 <- d$Sample
  list(summary = s, raw = d[, c("label", "Sample2", "value")], title = title, levels = lev)
}

f5c <- read_fig5_qpcr(6, "Dscam1 exon 19 deletion")
f5d <- read_fig5_qpcr(7, "Khc-73 exon 15 deletion")

draw_f5 <- function() {
  par(mfrow = c(2, 2), mar = c(5, 4.2, 2.7, 1), oma = c(0, 0, 1.2, 0))
  grouped_bars(f5a$summary, f5a$raw, "Sample2", c("Control", "Mutant"),
               paper_col[c("Control", "Mutant")], main = f5a$title, ylab = "PSI")
  grouped_bars(f5b$summary, f5b$raw, "Sample2", c("Control", "Mutant"),
               paper_col[c("Control", "Mutant")], main = f5b$title, ylab = "PSI")
  grouped_bars(f5c$summary, f5c$raw, "Sample2", f5c$levels,
               c("white", "#E8E8E8"), ylim = c(0, 1.5), main = f5c$title, ylab = "ext/uni")
  grouped_bars(f5d$summary, f5d$raw, "Sample2", f5d$levels,
               c("white", "#E8E8E8"), ylim = c(0, 1.5), main = f5d$title, ylab = "ext/uni")
  mtext("Figure 5 quantitative-panel reproduction — Source Data", outer = TRUE,
        cex = 1.05, font = 2)
}
open_devices("figure5_quantitative_reproduction", 10, 8); draw_f5(); dev.off()
save_png("figure5_quantitative_reproduction", 10, 8, draw_f5)

# Figure 6c: raw values occupy the left-hand columns of Source Data sheets 8–12.
read_fig6_metric <- function(sheet, metric, value_name, isoform_col = NULL) {
  d <- read_sheet(15, sheet)
  if (sheet == 8) {
    raw <- d[, 1:3]; names(raw) <- c("Gene", "Sample", "value")
  } else if (sheet == 12) {
    raw <- d[, 1:3]; names(raw) <- c("Gene", "value", "Sample")
  } else {
    raw <- d[, 1:4]; names(raw) <- c("Gene", "Isoform", "Sample", "value")
  }
  raw <- raw[is.finite(raw$value) & nzchar(raw$Gene) & nzchar(raw$Sample), ]
  raw$Sample2 <- ifelse(raw$Sample == "Control", "Control", "elav")
  raw$label <- raw$Gene
  genes <- unique(raw$Gene)
  rows <- do.call(rbind, lapply(genes, function(g) {
    x <- raw$value[raw$Gene == g & raw$Sample2 == "Control"]
    y <- raw$value[raw$Gene == g & raw$Sample2 == "elav"]
    p <- if (length(x) >= 2 && length(y) >= 2) t.test(x, y)$p.value else NA_real_
    rbind(data.frame(label = g, Sample2 = "Control", n = length(x), mean = mean(x), sem = sem(x), p = p),
          data.frame(label = g, Sample2 = "elav", n = length(y), mean = mean(y), sem = sem(y), p = p))
  }))
  list(summary = rows, raw = raw[, c("label", "Sample2", "value")],
       title = metric, ylab = value_name)
}

f6 <- list(
  read_fig6_metric(8, "dPAU", "dPAU"),
  read_fig6_metric(9, expression(PSI[All]), "PSI"),
  read_fig6_metric(11, expression(PSI[Short]), "PSI"),
  read_fig6_metric(10, expression(PSI[Long]), "PSI"),
  read_fig6_metric(12, expression(abs(PSI[Long]-PSI[Short])), "Absolute PSI difference")
)

draw_f6 <- function() {
  par(mfrow = c(5, 1), mar = c(5.8, 4, 1.8, 1), oma = c(0, 0, 1.2, 0))
  for (z in f6) {
    grouped_bars(z$summary, z$raw, "Sample2", c("Control", "elav"),
                 c(paper_col["Control"], paper_col["elav"]),
                 main = z$title, ylab = z$ylab)
  }
  mtext("Figure 6c reproduction — Source Data replicate values", outer = TRUE,
        cex = 1.05, font = 2)
}
open_devices("figure6c_reproduction", 13, 15); draw_f6(); dev.off()
save_png("figure6c_reproduction", 13, 15, draw_f6)

writeLines(c(
  "Reproduced quantitative panels:",
  "- Figure 3e-f: PSI recalculated from Supplementary Data 10 counts; observations with total reads <2 set to missing; paired two-sided t-tests.",
  "- Figure 5a-d: plotted replicate values from Source Data sheets 4-7; two-sided t-tests.",
  "- Figure 6c: plotted replicate values from Source Data sheets 8-12; two-sided Welch t-tests in this script.",
  "", 
  "Not reconstructed here:",
  "- Figure 3a-d enrichment, gene-body coverage, and read-length panels.",
  "- Figure 5 coverage tracks and gene/CRISPR schematics.",
  "- Figure 6a-b overlap/table and Figure 6d-e coverage tracks.",
  "", 
  "These outputs reproduce statistical plots from deposited processed/source data. They are not yet independent FASTQ-to-figure reproductions."
), file.path(out, "README.txt"))

cat("Wrote Figure 3, 5, and 6 reproductions to", out, "\n")
