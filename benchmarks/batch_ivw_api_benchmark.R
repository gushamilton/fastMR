#!/usr/bin/env Rscript

# Small, generic smoke benchmark for the public masked and sparse APIs.
# Usage: Rscript benchmarks/batch_ivw_api_benchmark.R [seed]

seed <- if (length(commandArgs(trailingOnly = TRUE))) {
  as.integer(commandArgs(trailingOnly = TRUE)[[1L]])
} else 20260807L
set.seed(seed)

exposures <- 128L
outcomes <- 128L
snps <- 800L
exposure_beta <- matrix(rnorm(exposures * snps, sd = 0.05), exposures, snps)
outcome_beta <- matrix(rnorm(outcomes * snps, sd = 0.05), outcomes, snps)
outcome_se <- matrix(runif(outcomes * snps, 0.02, 0.08), outcomes, snps)
exposure_present <- matrix(runif(exposures * snps) < 0.04, exposures, snps)
outcome_present <- matrix(runif(outcomes * snps) < 0.95, outcomes, snps)

row_ptr <- integer(exposures + 1L)
col_index <- integer()
stored_beta <- numeric()
for (i in seq_len(exposures)) {
  snp <- which(exposure_present[i, ])
  col_index <- c(col_index, snp - 1L)
  stored_beta <- c(stored_beta, exposure_beta[i, snp])
  row_ptr[i + 1L] <- length(col_index)
}

timed <- function(label, expression) {
  started <- proc.time()[["elapsed"]]
  result <- force(expression)
  elapsed <- proc.time()[["elapsed"]] - started
  cat(sprintf("%s: %.4f s (%d x %d output cells)\n",
              label, elapsed, exposures, outcomes))
  invisible(result)
}

masked <- timed("masked", fastMR::fast_mr_masked_ivw(
  exposure_beta, outcome_beta, outcome_se,
  exposure_present, outcome_present, max_memory_mb = 512
))
sparse <- timed("sparse", fastMR::fast_mr_sparse_ivw(
  row_ptr, col_index, stored_beta,
  outcome_beta, outcome_se, outcome_present, threads = 2L,
  max_memory_mb = 512
))
cat(sprintf("maximum beta difference on comparable results: %.3e\n",
            max(abs(masked$beta - sparse$beta), na.rm = TRUE)))
