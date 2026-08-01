args <- commandArgs(trailingOnly = TRUE)
root <- args[[1L]]
fast_lib <- file.path(root, ".local", "Rlib")
tsmr_lib <- "/Users/fergushamilton/projects/twosamplemr-fast/.local/Rlib"
.libPaths(c(fast_lib, tsmr_lib, .libPaths()))
suppressPackageStartupMessages({
  library(fastMR)
  library(TwoSampleMR)
})

fixture <- file.path(root, "inst", "extdata", "il6_crp_primary_100.tsv")
d <- read.delim(fixture, check.names = FALSE, stringsAsFactors = FALSE)
methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median",
             "mr_simple_mode", "mr_weighted_mode")
params <- modifyList(TwoSampleMR::default_parameters(), list(nboot = 100L))
seed <- 20260801L
nboot <- 100L
n <- 50L

make_grid <- function(d, n = 50L) {
  exp_base <- matrix(rep(d$beta.exposure, n), nrow = n, byrow = TRUE)
  out_base <- matrix(rep(d$beta.outcome, n), nrow = n, byrow = TRUE)
  exp_scale <- 1 + (seq_len(n) - (n + 1) / 2) * 0.001
  out_scale <- 1 + (seq_len(n) - (n + 1) / 2) * 0.001
  list(
    exposure_beta = exp_base * exp_scale,
    outcome_beta = out_base * out_scale,
    exposure_se = matrix(rep(d$se.exposure, n), nrow = n, byrow = TRUE),
    outcome_se = matrix(rep(d$se.outcome, n), nrow = n, byrow = TRUE)
  )
}

g <- make_grid(d, n)
pairs <- n * n

run_tsmr <- function() {
  results <- vector("list", pairs)
  index <- 0L
  for (exposure in seq_len(n)) {
    for (outcome in seq_len(n)) {
      index <- index + 1L
      dat <- data.frame(
        SNP = d$SNP,
        beta.exposure = g$exposure_beta[exposure, ],
        beta.outcome = g$outcome_beta[outcome, ],
        se.exposure = g$exposure_se[exposure, ],
        se.outcome = g$outcome_se[outcome, ],
        id.exposure = paste0("E", exposure),
        id.outcome = paste0("O", outcome),
        exposure = paste0("E", exposure),
        outcome = paste0("O", outcome),
        mr_keep = TRUE,
        stringsAsFactors = FALSE
      )
      results[[index]] <- suppressMessages(suppressWarnings(
        TwoSampleMR::mr(dat, method_list = methods, parameters = params)
      ))
    }
  }
  do.call(rbind, results)
}

fast_methods <- c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode")
invisible(fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se,
                       g$outcome_se, methods = fast_methods, nboot = 0L,
                       seed = seed, threads = 10L))
gc()

started <- proc.time()[["elapsed"]]
fast <- fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se,
                     g$outcome_se, methods = fast_methods, nboot = nboot,
                     seed = seed, threads = 10L)
fast_elapsed <- proc.time()[["elapsed"]] - started

started <- proc.time()[["elapsed"]]
tsmr <- run_tsmr()
tsmr_elapsed <- proc.time()[["elapsed"]] - started

fast_ivw <- fast[fast[["method_code"]] == "ivw", ]
tsmr_ivw <- tsmr[tsmr[["method"]] == "Inverse variance weighted", ]
correctness <- max(abs(fast_ivw[["b"]] - tsmr_ivw[["b"]]), na.rm = TRUE)
fast[["pair"]] <- (fast[["exposure_index"]] - 1L) * n + fast[["outcome_index"]]
tsmr[["pair"]] <- rep(seq_len(pairs), each = length(methods))
method_map <- c("Inverse variance weighted" = "ivw", "MR Egger" = "egger",
                "Weighted median" = "weighted_median", "Simple mode" = "simple_mode",
                "Weighted mode" = "weighted_mode")
tsmr[["method_code"]] <- unname(method_map[tsmr[["method"]]])
parity <- merge(fast[, c("pair", "method_code", "b", "se", "pval")],
                tsmr[, c("pair", "method_code", "b", "se", "pval")],
                by = c("pair", "method_code"), suffixes = c("_fast", "_tsmr"), sort = TRUE)
parity[["db"]] <- parity[["b_fast"]] - parity[["b_tsmr"]]
parity[["dse"]] <- parity[["se_fast"]] - parity[["se_tsmr"]]
parity[["dp"]] <- parity[["pval_fast"]] - parity[["pval_tsmr"]]
parity[["relative_se_difference"]] <- abs(parity[["dse"]]) / pmax(abs(parity[["se_tsmr"]]), 1e-15)
parity_summary <- do.call(rbind, lapply(split(parity, parity[["method_code"]]), function(x) {
  data.frame(method_code = x[["method_code"]][[1L]], pairs = nrow(x),
             max_abs_beta_difference = max(abs(x[["db"]]), na.rm = TRUE),
             max_abs_se_difference = max(abs(x[["dse"]]), na.rm = TRUE),
             median_abs_se_difference = median(abs(x[["dse"]]), na.rm = TRUE),
             max_relative_se_difference = max(x[["relative_se_difference"]], na.rm = TRUE),
             max_abs_pval_difference = max(abs(x[["dp"]]), na.rm = TRUE))
}))
print(parity_summary, row.names = FALSE, digits = 6)
write.csv(parity_summary, file.path(root, "outputs", "native_tsmr_grid_parity.csv"), row.names = FALSE)
saveRDS(list(fast = fast, TwoSampleMR = tsmr, parity = parity),
        file.path(root, "outputs", "native_tsmr_grid_results.rds"))
result <- data.frame(
  workload = "IL6 82-row fixture; 50x50 grid; five methods; nboot=100",
  pairs = pairs,
  nboot = nboot,
  seed = seed,
  fastMR_seconds = fast_elapsed,
  TwoSampleMR_seconds = tsmr_elapsed,
  speedup = tsmr_elapsed / fast_elapsed,
  fastMR_pairs_per_second = pairs / fast_elapsed,
  TwoSampleMR_pairs_per_second = pairs / tsmr_elapsed,
  max_ivw_beta_difference = correctness
)
print(result, row.names = FALSE)
write.csv(result, file.path(root, "outputs", "native_tsmr_grid_benchmark.csv"), row.names = FALSE)
writeLines(c(
  "# Native TwoSampleMR grid benchmark", "",
  sprintf("Fixture: %d rows; grid: %d pairs; methods: %d; nboot: %d; seed: %d.", nrow(d), pairs, length(methods), nboot, seed),
  "", "| implementation | wall seconds | pairs/s |", "|---|---:|---:|",
  sprintf("| fastMR exact R/C++ (threads=10) | %.6f | %.3f |", fast_elapsed, pairs / fast_elapsed),
  sprintf("| TwoSampleMR native R workflow | %.6f | %.3f |", tsmr_elapsed, pairs / tsmr_elapsed),
  "", sprintf("Speedup: %.3fx.", tsmr_elapsed / fast_elapsed),
  sprintf("Maximum absolute IVW beta difference: %.3e.", correctness)
), file.path(root, "outputs", "native_tsmr_grid_benchmark.md"))
