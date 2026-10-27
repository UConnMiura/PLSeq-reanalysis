options(stringsAsFactors = FALSE, scipen = 999)

root <- "/Users/pedromiura/Documents/codex/PMID37679364_reanalysis"
src <- file.path(
  "/Users/pedromiura/Documents/codex/PMID37679364_agent/dist",
  "pmid-37679364-agent/skill/pmid-37679364-paper/assets/supp_table"
)
out_fig <- file.path(root, "figures", "reproduced")
out_tab <- file.path(root, "results")
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)

read_sheet <- function(book, sheet) {
  read.csv(file.path(src, sprintf("41467-2023-41207-moesm%d-esm-sheet-%d.csv",
                                  book, sheet)), check.names = FALSE)
}

sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))

q_symbol <- function(q) {
  ifelse(!is.finite(q), "N/A",
         ifelse(q < 0.001, "***",
         ifelse(q < 0.01, "**",
         ifelse(q < 0.05, "*", "n.s."))))
}

d <- read_sheet(13, 1)
d$PSI[d$Total_counts < 2] <- NA_real_
d$key <- paste(d$Gene, d$CE_start, d$CE_end, sep = "|")

ce <- read_sheet(12, 1)
names(ce)[1:5] <- c("Chromosome", "CE_start", "CE_end", "Gene", "CE_number")
ce$key <- paste(ce$Gene, ce$CE_start, ce$CE_end, sep = "|")

stats <- list()
for (k in unique(d$key)) {
  z <- d[d$key == k, ]
  w <- reshape(z[, c("Run", "Isoform", "PSI")], idvar = "Run",
               timevar = "Isoform", direction = "wide")
  w <- w[complete.cases(w$PSI.Long, w$PSI.Short), ]
  tt <- if (nrow(w) >= 2) t.test(w$PSI.Long, w$PSI.Short, paired = TRUE) else NULL
  ce_num <- ce$CE_number[match(k, ce$key)]
  # These two assayed events use coordinates distinct from the compact
  # Supplementary Data 9 lookup row; retain the exon labels printed in Fig. 3e.
  if (z$Gene[1] == "CASK") ce_num <- 17
  if (z$Gene[1] == "gish") ce_num <- 14
  stats[[length(stats) + 1]] <- data.frame(
    key = k,
    Gene = z$Gene[1],
    CE = ce_num,
    n_pairs = nrow(w),
    mean_long = mean(w$PSI.Long),
    sem_long = sem(w$PSI.Long),
    mean_short = mean(w$PSI.Short),
    sem_short = sem(w$PSI.Short),
    delta_long_minus_short = mean(w$PSI.Long - w$PSI.Short),
    ci95_low = if (is.null(tt)) NA_real_ else tt$conf.int[1],
    ci95_high = if (is.null(tt)) NA_real_ else tt$conf.int[2],
    nominal_p = if (is.null(tt)) NA_real_ else tt$p.value
  )
}

stats <- do.call(rbind, stats)
stats$BH_q <- p.adjust(stats$nominal_p, method = "BH")
stats$BH_symbol <- q_symbol(stats$BH_q)
stats$nominal_significant <- stats$nominal_p < 0.05
stats$BH_significant <- stats$BH_q < 0.05
# Exact custom order and displayed cassette-exon labels from published Fig. 3e.
published_genes <- c(
  "Dys", "Khc-73", "pod1", "X11L", "spen", "CASK", "Mbs", "mtd",
  "nuf", "CG8671", "tutl", "Camta", "par-1", "Evi5", "shi", "CdGAPr",
  "Ptp10D", "Nc73EF", "gish", "Vav", "stai", "CG15312", "Crag", "lap",
  "Eip63E", "Calx", "Lar", "Dscam1"
)
published_ce <- c(
  33, 15, 9, 8, 18, 17, 11, 12, 9, 10, 15, 9, 13, 2,
  10, 4, 12, 15, 14, 12, 6, 6, 10, 15, 4, 4, 16, 23
)
stats <- stats[match(published_genes, stats$Gene), ]
stopifnot(identical(stats$Gene, published_genes))
stats$CE <- published_ce
write.csv(stats, file.path(out_tab, "figure3e_embryo_bh_results.csv"), row.names = FALSE)

draw_plot <- function() {
  old <- par(mar = c(7.8, 4.4, 4.2, 6.2), mgp = c(2.6, 0.8, 0), las = 1,
             xpd = FALSE)
  on.exit(par(old))

  means <- rbind(Long = stats$mean_long, Short = stats$mean_short)
  errors <- rbind(Long = stats$sem_long, Short = stats$sem_short)
  labels <- paste0(stats$Gene, "\nCE ", stats$CE)
  cols <- c(Long = "#79C5BB", Short = "#F7AA78")

  bp <- barplot(means, beside = TRUE, col = cols, border = NA,
                ylim = c(0, 1.12), names.arg = labels, las = 2,
                cex.names = 0.60, ylab = "Cassette-exon PSI",
                main = "Figure 3e: 16-18 hr embryos", axes = FALSE)
  axis(2, las = 1)
  axis(1, at = colMeans(bp), labels = FALSE, tick = FALSE, lwd = 1.2)
  box(bty = "l", lwd = 1.2)
  arrows(bp, means, bp, means + errors, angle = 90, code = 2,
         length = 0.025, lwd = 0.8)

  # Show the actual replicate measurements used in each paired test.
  for (j in seq_len(nrow(stats))) {
    z <- d[d$key == stats$key[j], ]
    for (i in seq_along(c("Long", "Short"))) {
      iso <- c("Long", "Short")[i]
      values <- z$PSI[z$Isoform == iso & is.finite(z$PSI)]
      points(rep(bp[i, j], length(values)) + seq(-0.045, 0.045, length.out = length(values)),
             values, pch = 21, bg = "white", col = "#555555", cex = 0.45)
    }
    text(mean(bp[, j]), 1.055, stats$BH_symbol[j], cex = 0.62)
  }

  par(xpd = NA)
  legend(x = max(bp) + 1.0, y = 1.08, xjust = 0, yjust = 1,
         legend = c("Long 3'UTR", "Short 3'UTR"),
         fill = cols, border = NA, bty = "n", horiz = FALSE, cex = 0.78)
  par(xpd = FALSE)
  mtext("Two-sided paired t-tests; BH correction across 28 embryo comparisons",
        side = 3, line = 0.45, cex = 0.76)
  mtext("* q<0.05, ** q<0.01, *** q<0.001; bars are mean +/- SEM",
        side = 1, line = 6.8, cex = 0.72)
}

pdf(file.path(out_fig, "figure3e_embryo_BH.pdf"), width = 13, height = 5.8,
    useDingbats = FALSE)
draw_plot()
dev.off()

png(file.path(out_fig, "figure3e_embryo_BH.png"), width = 3900, height = 1740,
    res = 300)
draw_plot()
dev.off()

cat(sprintf("Figure 3e: %d/28 nominal p<0.05; %d/28 BH q<0.05\n",
            sum(stats$nominal_significant, na.rm = TRUE),
            sum(stats$BH_significant, na.rm = TRUE)))
