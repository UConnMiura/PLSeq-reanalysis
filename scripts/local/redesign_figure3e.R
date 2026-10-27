options(stringsAsFactors = FALSE, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
full_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", full_args, value = TRUE)[1])
default_root <- normalizePath(file.path(dirname(script_file), "../.."), mustWork = TRUE)
project <- if (length(args)) args[1] else Sys.getenv("PLSEQ_ROOT", default_root)
input_dir <- file.path(project, "data", "derived")
output_dir <- file.path(project, "figures", "supplementary", "figure3e_redesign")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

counts <- read.delim(file.path(input_dir, "figure3e_reconstructed_counts.tsv"),
                     na.strings = "NA", check.names = FALSE)
stats <- read.delim(file.path(input_dir, "figure3e_FASTQ_statistics.tsv"),
                    na.strings = "NA", check.names = FALSE)

gene_order <- stats$gene
counts <- counts[counts$isoform %in% c("Long", "Short") &
                   counts$gene %in% gene_order, ]
counts$gene <- factor(counts$gene, levels = gene_order)

long_col <- "#178C8C"
short_col <- "#E07A3F"
grid_col <- "#E5E7EB"
text_col <- "#263238"

format_q <- function(x) {
  ifelse(is.na(x), "q = NA",
         ifelse(x < 0.001, paste0("q = ", formatC(x, format = "e", digits = 1)),
                paste0("q = ", formatC(x, format = "f", digits = 3))))
}
stats$q_label <- format_q(stats$BH_q)

draw_faithful <- function() {
  old <- par(mar = c(8.6, 4.5, 4.4, 1.1), mgp = c(2.6, 0.75, 0),
             las = 1, family = "sans", fg = text_col, col.axis = text_col,
             col.lab = text_col, col.main = text_col)
  on.exit(par(old))

  means <- rbind(Long = stats$mean_long, Short = stats$mean_short)
  sems <- rbind(Long = stats$sem_long, Short = stats$sem_short)
  labels <- paste0(stats$gene, "\nCE ", stats$CE)
  top <- pmin(1.02, apply(means + ifelse(is.finite(sems), sems, 0), 2, max,
                          na.rm = TRUE) + 0.055)

  bp <- barplot(means, beside = TRUE, col = c(long_col, short_col), border = NA,
                ylim = c(0, 1.17), names.arg = labels, las = 2,
                cex.names = 0.67, ylab = "Cassette-exon PSI", axes = FALSE,
                main = "Figure 3e — faithful presentation")
  abline(h = seq(0, 1, 0.25), col = grid_col, lwd = 0.8)
  abline(h = 0, col = text_col, lwd = 1.2)
  axis(2, at = seq(0, 1, 0.25), labels = sprintf("%.2f", seq(0, 1, 0.25)))

  ok <- is.finite(means) & is.finite(sems)
  arrows(bp[ok], means[ok] - sems[ok], bp[ok], means[ok] + sems[ok],
         angle = 90, code = 3, length = 0.018, lwd = 0.75, col = text_col)

  for (j in seq_len(nrow(stats))) {
    for (i in seq_along(c("Long", "Short"))) {
      iso <- c("Long", "Short")[i]
      z <- counts[counts$gene == stats$gene[j] & counts$isoform == iso &
                    is.finite(counts$PSI_depth2), ]
      if (nrow(z)) {
        jitter <- if (nrow(z) == 1) 0 else seq(-0.045, 0.045, length.out = nrow(z))
        points(rep(bp[i, j], nrow(z)) + jitter, z$PSI_depth2,
               pch = 21, bg = "white", col = text_col, lwd = 0.6, cex = 0.48)
      }
    }
    text(mean(bp[, j]), top[j], stats$q_label[j], srt = 90,
         adj = c(0, 0.5), cex = 0.52, col = text_col)
  }

  legend("topright", inset = c(0, -0.02), xpd = NA,
         legend = c("Long 3'UTR", "Short 3'UTR"),
         fill = c(long_col, short_col), border = NA, bty = "n",
         horiz = TRUE, cex = 0.82)
  mtext("Bars: mean ± SEM; points: biological replicates; q: BH-adjusted paired t-test",
        side = 3, line = 0.45, cex = 0.78, col = "#52606D")
}

draw_paired_panel <- function(panel_stats, show_legend = FALSE) {
  n <- nrow(panel_stats)
  plot(NA, xlim = c(0, 1.27), ylim = c(0.35, n + 0.65), axes = FALSE,
       xlab = "Cassette-exon PSI", ylab = "",
       xaxs = "i", yaxs = "i")
  abline(v = seq(0, 1, 0.25), col = grid_col, lwd = 0.8)
  axis(1, at = seq(0, 1, 0.25), labels = sprintf("%.2f", seq(0, 1, 0.25)))
  axis(2, at = seq_len(n),
       labels = paste0(rev(panel_stats$gene), "  CE ", rev(panel_stats$CE)),
       las = 1, tick = FALSE, cex.axis = 0.77)
  abline(v = 0, col = text_col, lwd = 1.1)

  panel_rev <- panel_stats[n:1, , drop = FALSE]
  for (row in seq_len(n)) {
    g <- panel_rev$gene[row]
    z <- counts[counts$gene == g & is.finite(counts$PSI_depth2), ]
    wide <- reshape(z[, c("run", "replicate", "isoform", "PSI_depth2")],
                    idvar = c("run", "replicate"), timevar = "isoform",
                    direction = "wide")
    wide <- wide[complete.cases(wide$PSI_depth2.Long, wide$PSI_depth2.Short), ]
    if (nrow(wide)) {
      offsets <- if (nrow(wide) == 1) 0 else seq(-0.20, 0.20, length.out = nrow(wide))
      yy <- row + offsets
      segments(wide$PSI_depth2.Short, yy, wide$PSI_depth2.Long, yy,
               col = "#AAB2B8", lwd = 0.9)
      points(wide$PSI_depth2.Long, yy, pch = 16, col = long_col, cex = 0.72)
      points(wide$PSI_depth2.Short, yy, pch = 17, col = short_col, cex = 0.72)
    }
    text(1.035, row, panel_rev$q_label[row], adj = 0,
         cex = 0.68, col = text_col)
  }
  text(1.035, n + 0.55, "BH-adjusted", adj = 0, cex = 0.66,
       font = 2, col = "#52606D")
  box(bty = "l", col = text_col)
  if (show_legend) {
    legend("top", inset = c(0, -0.085), xpd = NA, horiz = TRUE, bty = "n",
           legend = c("Long 3'UTR", "Short 3'UTR"),
           pch = c(16, 17),
           col = c(long_col, short_col),
           pt.cex = 0.85, cex = 0.80)
  }
}

draw_paired <- function() {
  old <- par(mfrow = c(1, 2), mar = c(4.0, 7.3, 4.5, 1.2),
             oma = c(1.2, 0, 2.8, 0), mgp = c(2.35, 0.72, 0),
             family = "sans", fg = text_col, col.axis = text_col,
             col.lab = text_col, col.main = text_col)
  on.exit(par(old))
  draw_paired_panel(stats[1:14, ], show_legend = TRUE)
  draw_paired_panel(stats[15:28, ])
  mtext("Figure 3e — paired-replicate presentation", side = 3, outer = TRUE,
        line = 1.35, cex = 1.22, font = 2, col = text_col)
  mtext("Each line joins long- and short-3'UTR PSI from the same biological replicate",
        side = 3, outer = TRUE, line = 0.25, cex = 0.78, col = "#52606D")
  mtext("Genes retain the published Figure 3e order; PSI values with <2 informative reads are omitted",
        side = 1, outer = TRUE, line = 0.15, cex = 0.72, col = "#52606D")
}

pdf(file.path(output_dir, "figure3e_faithful_numeric_q.pdf"), width = 14.5,
    height = 6.2, useDingbats = FALSE)
draw_faithful()
dev.off()
png(file.path(output_dir, "figure3e_faithful_numeric_q.png"), width = 4350,
    height = 1860, res = 300)
draw_faithful()
dev.off()

pdf(file.path(output_dir, "figure3e_paired_numeric_q.pdf"), width = 13.2,
    height = 8.4, useDingbats = FALSE)
draw_paired()
dev.off()
png(file.path(output_dir, "figure3e_paired_numeric_q.png"), width = 3960,
    height = 2520, res = 300)
draw_paired()
dev.off()

write.table(stats[, c("gene", "CE", "n_pairs", "mean_long", "sem_long",
                      "mean_short", "sem_short", "delta_long_minus_short",
                      "nominal_p", "BH_q")],
            file.path(output_dir, "figure3e_plot_statistics.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")

writeLines(c(
  "Figure 3e presentation redesigns",
  "=================================",
  "",
  "Both figures use the local FASTQ reconstruction and published gene order.",
  "Numerical q-values are Benjamini-Hochberg adjusted across all 28 genes.",
  "The faithful version retains grouped means, SEM bars, and replicate points.",
  "The paired version connects long- and short-3'UTR PSI within each replicate.",
  "PSI measurements with fewer than two informative reads remain excluded.",
  "Known event-coordinate limitations in the underlying reconstruction are unchanged."
), file.path(output_dir, "README.txt"))

cat("Wrote Figure 3e redesigns to:", output_dir, "\n")
