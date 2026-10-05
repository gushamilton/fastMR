#!/usr/bin/env Rscript
# Single-instrument (Wald ratio) agreement: Rscript wald_check.R RESDIR_REP SIZE   (env TRAITS = the store set)
# The showcase designs have no nsnp = 1 pairs, so this builds them: each exposure keeps only its strongest
# TwoSampleMR instrument (lowest exposure p among mr_keep rows), giving one-SNP pairs. Then compares
#   TwoSampleMR::mr() (mr_wald_ratio)  vs  fastMR::fast_mr() on the same harmonised rows (exactness), and
#   vs  fastMR::fast_mr_compressed() on the .cpr stores (the C arm's path; differs only by store quantisation).
a <- commandArgs(TRUE); R <- a[1]; size <- a[2]
src_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(src_dir, "common.R")); setup_libs(c("tsmr", "fastmr", "cs"))
suppressMessages({library(data.table); library(TwoSampleMR); library(fastMR); library(CompreSSoR)})
ref <- if (dir.exists(file.path(R, paste0(size, "_A")))) "A" else "A8"
dat <- as.data.table(readRDS(file.path(R, paste0(size, "_", ref), "harmonised.rds")))[mr_keep == TRUE]
top <- dat[order(pval.exposure), .(SNP = SNP[1]), by = id.exposure]
d1 <- as.data.frame(dat[top, on = .(id.exposure, SNP)])
meth <- c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode")
tsmr <- as.data.table(TwoSampleMR::mr(d1))
fm <- as.data.table(fastMR::fast_mr(d1, methods = meth, nboot = 1000, seed = 1))
ex <- unique(top$id.exposure); oc <- unique(d1$id.outcome)
ef <- setNames(file.path(TRAITS, "cpr", paste0(ex, ".cpr")), ex); of <- setNames(file.path(TRAITS, "cpr", paste0(oc, ".cpr")), oc)
inst <- setNames(lapply(ex, function(e) toupper(top[id.exposure == e, SNP])), ex)
fc <- as.data.table(fastMR::fast_mr_compressed(ef, of, inst, methods = meth, nboot = 1000, seed = 1, threads = 8, io_threads = 8))
k <- function(m) { m <- copy(m); m[, pair := paste(id.exposure, id.outcome, sep = "|")]; m }
w <- k(tsmr)[method == "Wald ratio", .(pair, b, se, pval)]
cmp <- function(x, lab) {
  x <- k(x); iv <- x[method == "Inverse variance weighted" & nsnp == 1, .(pair, b, se, pval)]
  j <- merge(w, iv, by = "pair", suffixes = c(".tsmr", ".x"))
  other <- x[nsnp == 1 & method != "Inverse variance weighted"]
  data.table(arm = lab, size = size, tsmr_wald_rows = nrow(w), arm_nsnp1_ivw_rows = nrow(iv), matched = nrow(j),
             max_abs_db = max(abs(j$b.tsmr - j$b.x)), max_abs_dse = max(abs(j$se.tsmr - j$se.x)), max_abs_dp = max(abs(j$pval.tsmr - j$pval.x)),
             med_db_se = median(abs(j$b.tsmr - j$b.x) / j$se.tsmr), max_db_se = max(abs(j$b.tsmr - j$b.x) / j$se.tsmr),
             other_methods_rows = nrow(other), other_methods_all_na = all(is.na(other$b)),
             tsmr_other_method_rows = nrow(k(tsmr)[method != "Wald ratio"]))
}
out <- rbind(cmp(fm, "fast_mr (same harmonised rows)"), cmp(fc, "fast_mr_compressed (.cpr stores)"))
print(out); fwrite(out, file.path(R, sprintf("agreement_%s_wald.csv", size)))
