#!/usr/bin/env Rscript
# Agreement across arms for one replicate dir: Rscript agree.R RESDIR_REP SIZE
a <- commandArgs(TRUE); R <- a[1]; size <- a[2]
src_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(src_dir, "common.R")); setup_libs(c("tsmr", "fastmr", "cs"))
suppressMessages(library(data.table))
arms <- c("A", "A8", "B", "C1", "C8"); L <- list()
for (x in arms) { f <- file.path(R, paste0(size, "_", x), "result.rds"); if (file.exists(f)) L[[x]] <- readRDS(f) }
jac <- function(a, b) length(intersect(a, b)) / max(1, length(union(a, b)))
out <- list()
ref <- if (!is.null(L[["A"]])) "A" else "A8"; stopifnot(!is.null(L[[ref]]))
for (x in setdiff(names(L), ref)) {
  ex <- names(L[[ref]]$instruments)
  out$instr[[x]] <- data.table(arm = x, exposure = ex, n_ref = lengths(L[[ref]]$instruments[ex]),
                               n_arm = lengths(L[[x]]$instruments[ex]),
                               jaccard = mapply(function(e) jac(L[[ref]]$instruments[[e]], L[[x]]$instruments[[e]]), ex))
  pk <- intersect(names(L[[ref]]$harm), names(L[[x]]$harm))
  out$harm[[x]] <- data.table(arm = x, pair = pk, n_ref = lengths(L[[ref]]$harm[pk]), n_arm = lengths(L[[x]]$harm[pk]),
                              jaccard = mapply(function(p) jac(L[[ref]]$harm[[p]], L[[x]]$harm[[p]]), pk))
}
# estimates: A vs others in SE units, per method
key <- function(m) {
  m <- as.data.table(m); m[, k := paste(id.exposure, id.outcome, method, sep = "|")]; m[, .(k, method, b, se, pval, nsnp)] }
ma <- key(L[[ref]]$mr)
for (x in setdiff(names(L), ref)) {
  mx <- key(L[[x]]$mr); j <- merge(ma, mx, by = c("k", "method"), suffixes = c(".A", ".x"))
  j[, `:=`(db_se = abs(b.A - b.x) / se.A, dse_rel = abs(se.A - se.x) / se.A, same_nsnp = nsnp.A == nsnp.x)]
  out$est[[x]] <- j[, .(arm = x, pairs = .N, frac_same_nsnp = mean(same_nsnp), med_db_se = median(db_se, na.rm = TRUE),
                        max_db_se = max(db_se, na.rm = TRUE), med_dse_rel = median(dse_rel, na.rm = TRUE), max_dse_rel = max(dse_rel, na.rm = TRUE)), by = method]
}
# exactness: fast_mr() on arm A's own harmonised data frame vs TSMR (same input, different code); IVW/Egger deterministic
if (requireNamespace("fastMR", quietly = TRUE)) {
  dat <- readRDS(file.path(R, paste0(size, "_", ref), "harmonised.rds"))
  fm <- as.data.table(fastMR::fast_mr(dat, methods = c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode"), nboot = 1000, seed = 1))
  mk <- key(fm); j <- merge(ma, mk, by = c("k", "method"), suffixes = c(".tsmr", ".fastmr"))
  out$exact <- j[, .(pairs = .N, max_abs_db = max(abs(b.tsmr - b.fastmr), na.rm = TRUE), max_abs_dse = max(abs(se.tsmr - se.fastmr), na.rm = TRUE),
                     max_abs_dp = max(abs(pval.tsmr - pval.fastmr), na.rm = TRUE), med_abs_dse = median(abs(se.tsmr - se.fastmr), na.rm = TRUE)), by = method]
}
for (n in names(out)) { cat("\n==", n, "==\n"); print(if (is.data.table(out[[n]])) out[[n]] else rbindlist(out[[n]])) }
saveRDS(out, file.path(R, paste0("agreement_", size, ".rds")))
for (n in names(out)) fwrite(if (is.data.table(out[[n]])) out[[n]] else rbindlist(out[[n]]), file.path(R, sprintf("agreement_%s_%s.csv", size, n)))
