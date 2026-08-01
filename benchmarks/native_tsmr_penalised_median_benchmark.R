args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]]) else getwd()
.libPaths(c(file.path(root, ".local", "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(fastMR); library(TwoSampleMR) })

d <- read.delim(file.path(root, "inst", "extdata", "il6_crp_primary_100.tsv"),
                check.names = FALSE, stringsAsFactors = FALSE)
scenarios <- data.frame(
  scenario = c("balanced_50x50", "one_exposure_250_outcomes",
               "one_outcome_250_exposures", "ten_exposures_100_outcomes",
               "hundred_exposures_ten_outcomes"),
  exposures = c(50L, 1L, 250L, 10L, 100L),
  outcomes = c(50L, 250L, 1L, 100L, 10L),
  stringsAsFactors = FALSE
)
nboot <- 100L
seed <- 20260904L
threads <- 5L

make_grid <- function(n_exposure, n_outcome) {
  list(
    exposure_beta = matrix(rep(d$beta.exposure, n_exposure), nrow = n_exposure,
                           byrow = TRUE) *
      (1 + (seq_len(n_exposure) - (n_exposure + 1) / 2) * 0.001),
    outcome_beta = matrix(rep(d$beta.outcome, n_outcome), nrow = n_outcome,
                          byrow = TRUE) *
      (1 + (seq_len(n_outcome) - (n_outcome + 1) / 2) * 0.001),
    exposure_se = matrix(rep(d$se.exposure, n_exposure), nrow = n_exposure,
                         byrow = TRUE),
    outcome_se = matrix(rep(d$se.outcome, n_outcome), nrow = n_outcome,
                        byrow = TRUE)
  )
}

run_native <- function(g, n_exposure, n_outcome, params) {
  out <- vector("list", n_exposure * n_outcome)
  index <- 0L
  for (exposure in seq_len(n_exposure)) {
    for (outcome in seq_len(n_outcome)) {
      index <- index + 1L
      z <- data.frame(
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
      out[[index]] <- suppressMessages(suppressWarnings(
        TwoSampleMR::mr(z, method_list = "mr_penalised_weighted_median",
                        parameters = params)))
    }
  }
  do.call(rbind, out)
}

params <- modifyList(TwoSampleMR::default_parameters(), list(nboot = nboot))
rows <- vector("list", nrow(scenarios))
for (i in seq_len(nrow(scenarios))) {
  s <- scenarios[i, ]
  g <- make_grid(s$exposures, s$outcomes)
  pairs <- s$exposures * s$outcomes
  invisible(fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se,
                         g$outcome_se, methods = "penalised_weighted_median",
                         nboot = 0L, seed = seed, threads = threads))
  started <- proc.time()[["elapsed"]]
  fast <- fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se,
                       g$outcome_se, methods = "penalised_weighted_median",
                       nboot = nboot, seed = seed, threads = threads)
  fast_seconds <- proc.time()[["elapsed"]] - started
  set.seed(seed)
  started <- proc.time()[["elapsed"]]
  native <- run_native(g, s$exposures, s$outcomes, params)
  native_seconds <- proc.time()[["elapsed"]] - started
  fast <- fast[order(fast$exposure_index, fast$outcome_index), ]
  native <- native[order(seq_len(nrow(native))), ]
  rows[[i]] <- data.frame(
    scenario = s$scenario, exposures = s$exposures, outcomes = s$outcomes,
    pairs = pairs, nboot = nboot, threads = threads,
    fastMR_seconds = fast_seconds, TwoSampleMR_seconds = native_seconds,
    speedup = native_seconds / fast_seconds,
    fastMR_pairs_per_second = pairs / fast_seconds,
    TwoSampleMR_pairs_per_second = pairs / native_seconds,
    max_abs_beta_difference = max(abs(fast$b - native$b), na.rm = TRUE),
    median_abs_se_difference = median(abs(fast$se - native$se), na.rm = TRUE),
    max_abs_se_difference = max(abs(fast$se - native$se), na.rm = TRUE),
    max_abs_pval_difference = max(abs(fast$pval - native$pval), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  print(rows[[i]], row.names = FALSE, digits = 6)
}

result <- do.call(rbind, rows)
dir.create(file.path(root, "outputs"), showWarnings = FALSE, recursive = TRUE)
write.csv(result, file.path(root, "outputs", "native_tsmr_penalised_median_benchmark.csv"),
          row.names = FALSE)
lines <- c(
  "# Native TwoSampleMR penalised weighted median benchmark", "",
  sprintf("IL6 fixture: %d SNPs; nboot=%d; fastMR threads=%d; seed=%d.",
          nrow(d), nboot, threads, seed),
  "The native comparator is the installed TwoSampleMR R workflow.", "",
  "| scenario | pairs | fastMR s | TwoSampleMR s | speedup | max beta delta | median SE delta |",
  "|---|---:|---:|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(result))) {
  x <- result[i, ]
  lines <- c(lines, sprintf(
    "| %s | %d | %.3f | %.3f | %.2fx | %.3e | %.3e |",
    x$scenario, x$pairs, x$fastMR_seconds, x$TwoSampleMR_seconds,
    x$speedup, x$max_abs_beta_difference, x$median_abs_se_difference))
}
writeLines(lines, file.path(root, "outputs", "native_tsmr_penalised_median_benchmark.md"))
