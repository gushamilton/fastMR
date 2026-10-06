#!/usr/bin/env Rscript
# Multi-trait batch extract (read_sumstats_batch, 20 stores, K keys, T threads) per build, fresh process:
# elapsed, time inside CompreSSoR's shared-identity verification, and the same call with that
# verification stubbed out (lib7+ only; for attribution, not a supported mode).
# Rscript prof_batch_extract.R LIB CPR_DIR K T [noverify|share]   (share = options(CompreSSoR.batch_share_panels = TRUE),
# the 0.7.0 default that CompreSSoR #55 made opt-in). RPROF=1 adds an Rprof summary.
a <- commandArgs(TRUE); .libPaths(c(a[1], .libPaths())); suppressMessages({library(CompreSSoR); library(data.table)})
K <- as.integer(a[3]); T <- as.integer(a[4]); nover <- isTRUE(a[5] == "noverify")
if (isTRUE(a[5] == "share")) options(CompreSSoR.batch_share_panels = TRUE)
ids <- c(sprintf("exp%02d", 1:10), sprintf("out%02d", 1:10)); P <- file.path(a[2], paste0(ids, ".cpr"))
keys <- fread(sprintf("/user/work/fh6520/showcase/storage/multi/lists/keys_%d.tsv", K), colClasses = list(character = c("chrom", "ref", "alt")))
v <- paste(keys$chrom, keys$pos, keys$ref, keys$alt, sep = ":")
cols <- c("chromosome", "base_pair_location", "other_allele", "effect_allele", "beta", "standard_error", "effect_allele_frequency", "p_value")   # CPR_COLS of the multi suite
vt <- 0
if (exists("pcodec_verify_identity_files", asNamespace("CompreSSoR"))) {
  real <- CompreSSoR:::pcodec_verify_identity_files
  f <- if (nover) function(stores, threads = 1L) invisible(TRUE) else function(stores, threads = 1L) {
    t0 <- proc.time()[["elapsed"]]; on.exit(vt <<- vt + proc.time()[["elapsed"]] - t0); real(stores, threads) }
  assignInNamespace("pcodec_verify_identity_files", f, "CompreSSoR")
}
prof <- nzchar(Sys.getenv("RPROF"))
if (prof) Rprof(pf <- tempfile(), interval = 0.01)
t <- system.time(r <- read_sumstats_batch(P, variants = v, columns = cols, threads = T))[["elapsed"]]
if (prof) { Rprof(NULL); s <- summaryRprof(pf); print(head(s$by.total[order(-s$by.total$total.time), ], 35)); print(head(s$by.self, 15)) }
cat(sprintf("CompreSSoR=%s K=%d T=%d mode=%s elapsed=%.3f verify=%.3f rows=%d\n", packageVersion("CompreSSoR"), K, T,
            if (is.na(a[5])) "default" else a[5], t, vt, sum(vapply(r, nrow, 0L))))
