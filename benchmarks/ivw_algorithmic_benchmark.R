args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]]) else getwd()
.libPaths(c(file.path(root, ".local", "Rlib"), .libPaths()))
suppressPackageStartupMessages(library(fastMR))
d <- read.delim(file.path(root, "inst", "extdata", "il6_crp_primary_100.tsv"),
                check.names = FALSE, stringsAsFactors = FALSE)
n <- 50L
g <- list(
  exposure_beta = matrix(rep(d[["beta.exposure"]], n), nrow = n, byrow = TRUE),
  outcome_beta = matrix(rep(d[["beta.outcome"]], n), nrow = n, byrow = TRUE),
  exposure_se = matrix(rep(d[["se.exposure"]], n), nrow = n, byrow = TRUE),
  outcome_se = matrix(rep(d[["se.outcome"]], n), nrow = n, byrow = TRUE)
)
for (threads in c(1L, 5L, 10L)) {
  invisible(fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se,
                         g$outcome_se, methods = "ivw", nboot = 0, threads = threads))
}
rows <- lapply(c(1L, 5L, 10L), function(threads) {
  elapsed <- replicate(5L, {
    system.time(fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se,
                             g$outcome_se, methods = "ivw", nboot = 0,
                             threads = threads))[["elapsed"]]
  })
  data.frame(threads = threads, pairs = 2500L, repeats = length(elapsed),
             median_fastMR_seconds = median(elapsed), min_fastMR_seconds = min(elapsed),
             scalar_pre_blas_seconds = 1.313,
             speedup_vs_scalar = 1.313 / median(elapsed))
})
result <- do.call(rbind, rows)
dir.create(file.path(root, "outputs"), showWarnings = FALSE, recursive = TRUE)
write.csv(result, file.path(root, "outputs", "ivw_algorithmic_benchmark.csv"), row.names = FALSE)
lines <- c(
  "# IVW algorithmic optimisation benchmark", "",
  "IL6 fixture: 82 SNPs; 50×50 grid; nboot=0; five warm repeats per thread count.",
  "The scalar baseline is the measured pre-BLAS fastMR path on the same Mac mini and command (1.313 s).",
  "The new path batches IVW numerator/denominator cross-products through BLAS and flattens tidy grid output once.",
  "", "| threads | median fastMR s | min fastMR s | pre-BLAS scalar s | speedup |",
  "|---:|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(result))) lines <- c(lines, sprintf(
  "| %d | %.3f | %.3f | %.3f | %.2fx |", result$threads[i],
  result$median_fastMR_seconds[i], result$min_fastMR_seconds[i],
  result$scalar_pre_blas_seconds[i], result$speedup_vs_scalar[i]))
writeLines(lines, file.path(root, "outputs", "ivw_algorithmic_benchmark.md"))
print(result, row.names = FALSE)
