options(stringsAsFactors = FALSE, scipen = 999)

full_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", full_args, value = TRUE)[1])
default_root <- normalizePath(file.path(dirname(script_file), "../.."), mustWork = TRUE)
root <- Sys.getenv("PLSEQ_ROOT", default_root)
src <- Sys.getenv("PLSEQ_SUPP_TABLE_DIR")
if (!nzchar(src)) stop("Set PLSEQ_SUPP_TABLE_DIR to the extracted supplementary-table directory")
res <- file.path(root, "data", "derived")
fig <- file.path(root, "figures", "supplementary", "phase1")
dir.create(res, recursive = TRUE, showWarnings = FALSE)
dir.create(fig, recursive = TRUE, showWarnings = FALSE)

read_sheet <- function(n, sheet) {
  read.csv(file.path(src, sprintf("41467-2023-41207-moesm%d-esm-sheet-%d.csv", n, sheet)), check.names = FALSE)
}

safe_t <- function(x, y, paired = FALSE, equal.var = FALSE) {
  ok <- if (paired) complete.cases(x, y) else c(rep(TRUE, length(x)), rep(TRUE, length(y)))
  if (paired) { x <- x[ok]; y <- y[ok] } else { x <- x[is.finite(x)]; y <- y[is.finite(y)] }
  if (length(x) < 2 || length(y) < 2) return(c(p = NA, estimate = NA, lower = NA, upper = NA, effect = NA))
  tt <- try(t.test(x, y, paired = paired, var.equal = equal.var, conf.level = .95), silent = TRUE)
  if (inherits(tt, "try-error")) return(c(p = NA, estimate = mean(x)-mean(y), lower = NA, upper = NA, effect = NA))
  d <- mean(x) - mean(y)
  effect <- if (paired) d / sd(x-y) else d / sqrt(((length(x)-1)*var(x)+(length(y)-1)*var(y))/(length(x)+length(y)-2))
  c(p = tt$p.value, estimate = d, lower = tt$conf.int[1], upper = tt$conf.int[2], effect = effect)
}

audit_psi <- function(d, panel) {
  calc <- with(d, ifelse(Total_counts > 0, Inclusion_counts / Total_counts, NA_real_))
  data.frame(panel = panel, rows = nrow(d), max_abs_error = max(abs(calc-d$PSI), na.rm=TRUE),
             mismatches_gt_1e_8 = sum(abs(calc-d$PSI) > 1e-8, na.rm=TRUE))
}

paired_isoform <- function(d, panel) {
  d <- d[d$Isoform %in% c("Short", "Long"), ]
  w <- reshape(d[, c("Run","Gene","Isoform","PSI")], idvar=c("Run","Gene"), timevar="Isoform", direction="wide")
  genes <- unique(w$Gene)
  out <- do.call(rbind, lapply(genes, function(g) {
    z <- w[w$Gene == g, ]
    z <- z[complete.cases(z$PSI.Short, z$PSI.Long), ]
    st <- safe_t(z$PSI.Short, z$PSI.Long, paired=TRUE)
    data.frame(panel=panel, Gene=g, n=sum(complete.cases(z$PSI.Short,z$PSI.Long)),
      mean_short=mean(z$PSI.Short,na.rm=TRUE), mean_long=mean(z$PSI.Long,na.rm=TRUE),
      delta_short_minus_long=unname(st["estimate"]), ci95_low=unname(st["lower"]), ci95_high=unname(st["upper"]),
      cohen_dz=unname(st["effect"]), p_value=unname(st["p"]))
  }))
  out$bh_fdr <- p.adjust(out$p_value, method="BH")
  out$raw_sig_0.05 <- out$p_value < .05
  out$bh_sig_0.05 <- out$bh_fdr < .05
  out[order(out$p_value), ]
}

two_group_isoform <- function(d, panel, group1="Control", group2="elav") {
  out <- list(); k <- 1
  for (iso in unique(d$Isoform)) for (g in unique(d$Gene)) {
    z <- d[d$Gene==g & d$Isoform==iso,]
    x <- z$PSI[z$Sample==group1]; y <- z$PSI[z$Sample==group2]
    x <- x[is.finite(x)]; y <- y[is.finite(y)]
    if (length(x) < 2 || length(y) < 2) next
    st <- safe_t(x,y,paired=FALSE,equal.var=FALSE)
    out[[k]] <- data.frame(panel=panel, Gene=g, Isoform=iso, n_control=length(x), n_mutant=length(y),
      mean_control=mean(x), mean_mutant=mean(y), delta_control_minus_mutant=unname(st["estimate"]),
      ci95_low=unname(st["lower"]), ci95_high=unname(st["upper"]), cohen_d=unname(st["effect"]), p_value=unname(st["p"]))
    k <- k+1
  }
  ans <- do.call(rbind,out); ans$bh_fdr <- ave(ans$p_value, ans$Isoform, FUN=function(x)p.adjust(x,"BH"))
  ans$raw_sig_0.05 <- ans$p_value < .05; ans$bh_sig_0.05 <- ans$bh_fdr < .05
  ans[order(ans$Isoform,ans$p_value),]
}

emb <- read_sheet(13,1); head_d <- read_sheet(13,2); deletion <- read_sheet(13,3); elav <- read_sheet(13,4); endov <- read_sheet(13,5)
psi_audit <- rbind(audit_psi(emb,"Figure 3 embryo"), audit_psi(head_d,"Figure 3 head"), audit_psi(deletion,"Figure 5 deletion"), audit_psi(elav,"Figure 6 elav"), audit_psi(endov,"Supplement Endov"))
write.csv(psi_audit,file.path(res,"phase1_psi_arithmetic_audit.csv"),row.names=FALSE)

# The methods require at least two reads for a PSI estimate. Retain raw rows for the
# arithmetic audit above, then mark lower-coverage observations missing for inference.
for (nm in c("emb","head_d","deletion","elav","endov")) {
  x <- get(nm); x$PSI[x$Total_counts < 2] <- NA_real_; assign(nm,x)
}

fig3 <- rbind(paired_isoform(emb,"Embryo"),paired_isoform(head_d,"Adult head"))
write.csv(fig3,file.path(res,"phase1_figure3_recomputed.csv"),row.names=FALSE)

# Figure 6: Welch two-sample tests reproduce the Source Data convention.
fig6 <- two_group_isoform(elav,"Figure 6")
write.csv(fig6,file.path(res,"phase1_figure6_recomputed.csv"),row.names=FALSE)

# Recompute the within-run absolute difference between long- and short-isoform PSI, then compare genotypes.
w6 <- reshape(elav[elav$Isoform %in% c("Long","Short"),c("Run","Gene","Sample","Isoform","PSI")], idvar=c("Run","Gene","Sample"), timevar="Isoform", direction="wide")
w6$abs_delta <- abs(w6$PSI.Short-w6$PSI.Long)
abs_rows <- list(); k <- 1
for (g in unique(w6$Gene)) {
  z <- w6[w6$Gene==g,]; x <- z$abs_delta[z$Sample=="Control"]; y <- z$abs_delta[z$Sample=="elav"]
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if(length(x)<2 || length(y)<2) next
  st <- safe_t(x,y,equal.var=FALSE)
  abs_rows[[k]] <- data.frame(Gene=g,n_control=length(x),n_mutant=length(y),mean_control=mean(x),mean_mutant=mean(y),
    delta_control_minus_mutant=unname(st["estimate"]),ci95_low=unname(st["lower"]),ci95_high=unname(st["upper"]),cohen_d=unname(st["effect"]),p_value=unname(st["p"])); k<-k+1
}
fig6_abs <- do.call(rbind,abs_rows); fig6_abs$bh_fdr <- p.adjust(fig6_abs$p_value,"BH"); fig6_abs$raw_sig_0.05 <- fig6_abs$p_value<.05; fig6_abs$bh_sig_0.05 <- fig6_abs$bh_fdr<.05
write.csv(fig6_abs[order(fig6_abs$p_value),],file.path(res,"phase1_figure6_abs_delta_recomputed.csv"),row.names=FALSE)

# Figure 5 and Endov are descriptive audits because the deletion removes the assayed long isoform.
desc <- function(d,panel) aggregate(PSI~Sample+Gene+Isoform,d,function(x)c(n=length(x),mean=mean(x),sem=sd(x)/sqrt(length(x))))
f5 <- desc(deletion,"Figure 5"); f5out <- data.frame(Sample=f5$Sample,Gene=f5$Gene,Isoform=f5$Isoform,n=f5$PSI[,"n"],mean=f5$PSI[,"mean"],sem=f5$PSI[,"sem"])
write.csv(f5out,file.path(res,"phase1_figure5_descriptive.csv"),row.names=FALSE)
e1 <- aggregate(PSI~Isoform,endov,function(x)c(n=length(x),mean=mean(x),sem=sd(x)/sqrt(length(x))))
eout <- data.frame(Isoform=e1$Isoform,n=e1$PSI[,"n"],mean=e1$PSI[,"mean"],sem=e1$PSI[,"sem"])
write.csv(eout,file.path(res,"phase1_endov_descriptive.csv"),row.names=FALSE)

# Phase 2 threshold sensitivity from supplied QAPA/rMATS result tables.
sets <- data.frame(label=c("Embryo 2-4h vs 16-18h","Ovary vs head","CS vs elav,fne","CS vs elav","mES DIV-8 vs DIV7"), moesm=c(4,5,7,9,10))
qapa_grid <- list(); rmats_grid <- list(); overlap <- list(); coverage <- list(); iq<-1; ir<-1; io<-1; ic<-1
for(i in seq_len(nrow(sets))) {
  q <- read_sheet(sets$moesm[i],1); r <- read_sheet(sets$moesm[i],2)
  q$abs_log2_fc <- abs(log2(q$Fold_change)); r$abs_dpsi <- abs(r$IncLevelDifference)
  q_pub <- is.finite(q$P_adj)&q$P_adj<.05&is.finite(q$Fold_change)&(q$Fold_change>2|q$Fold_change<.5)
  r_pub <- is.finite(r$FDR)&r$FDR<.05&is.finite(r$abs_dpsi)&r$abs_dpsi>.2
  qgene <- if("Gene_Name" %in% names(q)) q$Gene_Name else q$Gene
  rgene <- r$geneSymbol
  overlap[[io]] <- data.frame(contrast=sets$label[i],qapa_events=sum(q_pub),qapa_genes=length(unique(qgene[q_pub])),rmats_events=sum(r_pub),rmats_genes=length(unique(rgene[r_pub])),overlap_genes=length(intersect(unique(qgene[q_pub]),unique(rgene[r_pub])))); io<-io+1
  coverage[[ic]] <- data.frame(contrast=sets$label[i],qapa_rows=nrow(q),qapa_all_rows_meet_published=sum(q_pub,na.rm=TRUE)==nrow(q),qapa_min_padj=min(q$P_adj,na.rm=TRUE),qapa_min_abs_log2_fc=min(q$abs_log2_fc,na.rm=TRUE),rmats_rows=nrow(r),rmats_all_rows_meet_published=sum(r_pub,na.rm=TRUE)==nrow(r),rmats_min_fdr=min(r$FDR,na.rm=TRUE),rmats_min_abs_dpsi=min(r$abs_dpsi,na.rm=TRUE)); ic<-ic+1
  for(fdr in c(.01,.05)) for(fc in c(1.5,2,3)) {
    keep <- is.finite(q$P_adj)&q$P_adj<fdr&is.finite(q$Fold_change)&(q$Fold_change>fc|q$Fold_change<1/fc)
    qapa_grid[[iq]] <- data.frame(contrast=sets$label[i],fdr=fdr,fold_cutoff=fc,events=sum(keep),genes=length(unique(qgene[keep]))); iq<-iq+1
  }
  for(fdr in c(.01,.05)) for(ds in c(.1,.2,.3,.4)) {
    keep <- is.finite(r$FDR)&r$FDR<fdr&is.finite(r$abs_dpsi)&r$abs_dpsi>ds
    rmats_grid[[ir]] <- data.frame(contrast=sets$label[i],fdr=fdr,abs_dpsi_cutoff=ds,events=sum(keep),genes=length(unique(rgene[keep]))); ir<-ir+1
  }
}
write.csv(do.call(rbind,qapa_grid),file.path(res,"phase2_qapa_sensitivity.csv"),row.names=FALSE)
write.csv(do.call(rbind,rmats_grid),file.path(res,"phase2_rmats_sensitivity.csv"),row.names=FALSE)
write.csv(do.call(rbind,overlap),file.path(res,"phase2_apa_as_overlap.csv"),row.names=FALSE)
write.csv(do.call(rbind,coverage),file.path(res,"phase2_table_coverage_audit.csv"),row.names=FALSE)

# Figures
png(file.path(fig,"figure3_effects_and_multiplicity.png"),width=1800,height=900,res=180)
par(mfrow=c(1,2),mar=c(5,10,3,2))
for(p in c("Embryo","Adult head")) {
  z<-fig3[fig3$panel==p,]; z<-z[order(z$delta_short_minus_long),]; cols<-ifelse(z$bh_sig_0.05,"#0072B2",ifelse(z$raw_sig_0.05,"#E69F00","#999999"))
  plot(z$delta_short_minus_long,seq_len(nrow(z)),xlim=range(c(z$ci95_low,z$ci95_high),na.rm=T),pch=19,col=cols,yaxt="n",xlab="Mean PSI(short) - PSI(long), 95% CI",ylab="",main=p); axis(2,at=seq_len(nrow(z)),labels=z$Gene,las=2,cex.axis=.62); abline(v=0,lty=2); segments(z$ci95_low,seq_len(nrow(z)),z$ci95_high,seq_len(nrow(z)),col=cols)
}
dev.off()

qg<-do.call(rbind,qapa_grid); rg<-do.call(rbind,rmats_grid)
png(file.path(fig,"phase2_threshold_sensitivity.png"),width=1800,height=900,res=180)
par(mfrow=c(1,2),mar=c(8,5,3,1))
z<-qg[qg$fdr==.05,]; mat<-xtabs(events~contrast+fold_cutoff,z); matplot(t(mat),type="b",pch=19,lty=1,xaxt="n",xlab="Fold-change cutoff",ylab="QAPA events",main="QAPA sensitivity (FDR < 0.05)"); axis(1,at=seq_len(ncol(mat)),labels=colnames(mat)); legend("topright",legend=rownames(mat),col=seq_len(nrow(mat)),lty=1,pch=19,cex=.65)
z<-rg[rg$fdr==.05,]; mat<-xtabs(events~contrast+abs_dpsi_cutoff,z); matplot(t(mat),type="b",pch=19,lty=1,xaxt="n",xlab="Absolute delta PSI cutoff",ylab="rMATS events",main="rMATS sensitivity (FDR < 0.05)"); axis(1,at=seq_len(ncol(mat)),labels=colnames(mat)); legend("topright",legend=rownames(mat),col=seq_len(nrow(mat)),lty=1,pch=19,cex=.65)
dev.off()

summary <- data.frame(
  metric=c("PSI arithmetic mismatches","Figure 3 raw p<0.05","Figure 3 BH FDR<0.05","Figure 6 isoform tests raw p<0.05","Figure 6 isoform tests BH FDR<0.05","Figure 6 abs-delta raw p<0.05","Figure 6 abs-delta BH FDR<0.05"),
  value=c(sum(psi_audit$mismatches_gt_1e_8),sum(fig3$raw_sig_0.05,na.rm=T),sum(fig3$bh_sig_0.05,na.rm=T),sum(fig6$raw_sig_0.05,na.rm=T),sum(fig6$bh_sig_0.05,na.rm=T),sum(fig6_abs$raw_sig_0.05,na.rm=T),sum(fig6_abs$bh_sig_0.05,na.rm=T)))
write.csv(summary,file.path(res,"headline_summary.csv"),row.names=FALSE)
cat("Analysis complete\n")
