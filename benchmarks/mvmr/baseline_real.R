# Baseline for the real benchmark: fast_mr_multivariable() looped over the
# 2,940 outcomes for a subset of exposures (M1 designs), extrapolated.
# Usage: Rscript baseline_real.R <lib> <stores_dir> <inputs.rds> <out.csv> [n_exposures]
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[[1]], "/user/work/fh6520/r_packages", .libPaths()))
suppressPackageStartupMessages(library(fastMR))
stores <- args[[2]]; inp <- readRDS(args[[3]]); out_csv <- args[[4]]
n_e <- if (length(args) >= 5) as.integer(args[[5]]) else 5L
set.seed(1)
ex <- sample(names(inp$m1), n_e)
keys <- unique(unlist(inp$m1[ex]))
lab <- inp$excl$id.outcome
paths <- file.path(stores, lab)
codecs <- fastMR:::fastmr_compressed_validate_stores(paths, 8L)
dat <- fastMR:::fastmr_io_map(paths, rep(list(keys), length(paths)), c("beta", "standard_error"), 8L,
                              codecs = unname(codecs[paths]))
B <- S <- matrix(NA_real_, length(keys), length(lab), dimnames = list(keys, lab))
for (k in seq_along(dat)) { h <- match(keys, dat[[k]]$variant_key); B[, k] <- dat[[k]]$beta[h]; S[, k] <- dat[[k]]$standard_error[h] }
Braw <- B
f <- fastMR:::fastmr_mvmr_key_fields(keys); pos <- as.numeric(f[, 2])
for (k in seq_along(lab)) {
  hit <- f[, 1] == as.character(inp$excl$chromosome[k]) & pos >= inp$excl$start[k] & pos <= inp$excl$end[k]
  B[hit, k] <- NA
}
plt <- inp$covs$PLT
t0 <- proc.time()[["elapsed"]]
n_fit <- 0
for (e in ex) {
  s <- inp$m1[[e]]
  X <- cbind(exposure = Braw[s, e], PLT = plt$beta[match(s, plt$SNP)])
  for (k in seq_along(lab)) {
    ok <- is.finite(B[s, k]) & is.finite(S[s, k])
    r <- fast_mr_multivariable(X[ok, , drop = FALSE], B[s, k][ok], S[s, k][ok])
    n_fit <- n_fit + 1
  }
}
t <- proc.time()[["elapsed"]] - t0
cpu <- sub(".*: *", "", grep("model name", readLines("/proc/cpuinfo"), value = TRUE)[1])
out <- data.frame(case = "fast_mr_multivariable_loop", exposures = n_e, outcomes = length(lab),
                  fits = n_fit, seconds = t, per_fit_ms = 1000 * t / n_fit,
                  extrapolated_2610_seconds = t / n_e * length(inp$m1), cpu = cpu)
print(out)
write.csv(out, out_csv, row.names = FALSE)
