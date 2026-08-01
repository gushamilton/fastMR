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
method_specs <- data.frame(
  code = c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode"),
  tsmr = c("mr_ivw", "mr_egger_regression", "mr_weighted_median",
            "mr_simple_mode", "mr_weighted_mode"),
  display = c("Inverse variance weighted", "MR Egger", "Weighted median",
              "Simple mode", "Weighted mode"),
  stringsAsFactors = FALSE
)
scenarios <- data.frame(
  scenario = c("balanced_50x50", "one_exposure_250_outcomes",
               "one_outcome_250_exposures", "ten_exposures_100_outcomes",
               "hundred_exposures_ten_outcomes"),
  exposures = c(50L, 1L, 250L, 10L, 100L),
  outcomes = c(50L, 250L, 1L, 100L, 10L),
  stringsAsFactors = FALSE
)
nboot <- 100L
seed <- 20260802L
threads <- 5L

make_grid <- function(d, n_exposure, n_outcome) {
  exp_base <- matrix(rep(d$beta.exposure, n_exposure), nrow = n_exposure, byrow = TRUE)
  out_base <- matrix(rep(d$beta.outcome, n_outcome), nrow = n_outcome, byrow = TRUE)
  exp_scale <- 1 + (seq_len(n_exposure) - (n_exposure + 1) / 2) * 0.001
  out_scale <- 1 + (seq_len(n_outcome) - (n_outcome + 1) / 2) * 0.001
  list(
    exposure_beta = exp_base * exp_scale,
    outcome_beta = out_base * out_scale,
    exposure_se = matrix(rep(d$se.exposure, n_exposure), nrow = n_exposure, byrow = TRUE),
    outcome_se = matrix(rep(d$se.outcome, n_outcome), nrow = n_outcome, byrow = TRUE)
  )
}

run_tsmr <- function(g, n_exposure, n_outcome, tsmr_method, params) {
  pairs <- n_exposure * n_outcome
  results <- vector("list", pairs)
  index <- 0L
  for (exposure in seq_len(n_exposure)) {
    for (outcome in seq_len(n_outcome)) {
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
        TwoSampleMR::mr(dat, method_list = tsmr_method, parameters = params)
      ))
    }
  }
  out <- do.call(rbind, results)
  out$pair <- seq_len(pairs)
  out
}

params <- modifyList(TwoSampleMR::default_parameters(), list(nboot = nboot))
rows <- vector("list", nrow(scenarios) * nrow(method_specs))
row_index <- 0L

for (scenario_index in seq_len(nrow(scenarios))) {
  scenario <- scenarios[scenario_index, ]
  g <- make_grid(d, scenario$exposures, scenario$outcomes)
  pairs <- scenario$exposures * scenario$outcomes

  for (method_index in seq_len(nrow(method_specs))) {
    spec <- method_specs[method_index, ]
    row_index <- row_index + 1L
    cat(sprintf("[%02d/%02d] %s / %s (%d pairs)\n", row_index, length(rows),
                scenario$scenario, spec$code, pairs))

    invisible(fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se,
                           g$outcome_se, methods = spec$code, nboot = 0L,
                           seed = seed, threads = threads))
    fast_started <- proc.time()[["elapsed"]]
    fast <- fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se,
                         g$outcome_se, methods = spec$code, nboot = nboot,
                         seed = seed, threads = threads)
    fast_seconds <- proc.time()[["elapsed"]] - fast_started

    set.seed(seed)
    native_started <- proc.time()[["elapsed"]]
    native <- run_tsmr(g, scenario$exposures, scenario$outcomes,
                       spec$tsmr, params)
    native_seconds <- proc.time()[["elapsed"]] - native_started

    fast <- fast[order(fast$exposure_index, fast$outcome_index), ]
    native <- native[order(native$pair), ]
    beta_delta <- fast$b - native$b
    se_delta <- fast$se - native$se
    pval_delta <- fast$pval - native$pval
    rows[[row_index]] <- data.frame(
      scenario = scenario$scenario,
      exposures = scenario$exposures,
      outcomes = scenario$outcomes,
      pairs = pairs,
      method_code = spec$code,
      method = spec$display,
      nboot = nboot,
      threads = threads,
      fastMR_seconds = fast_seconds,
      TwoSampleMR_seconds = native_seconds,
      speedup = native_seconds / fast_seconds,
      fastMR_pairs_per_second = pairs / fast_seconds,
      TwoSampleMR_pairs_per_second = pairs / native_seconds,
      max_abs_beta_difference = max(abs(beta_delta), na.rm = TRUE),
      median_abs_se_difference = median(abs(se_delta), na.rm = TRUE),
      max_abs_se_difference = max(abs(se_delta), na.rm = TRUE),
      max_relative_se_difference = max(abs(se_delta) / pmax(abs(native$se), 1e-15), na.rm = TRUE),
      max_abs_pval_difference = max(abs(pval_delta), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    print(rows[[row_index]], row.names = FALSE, digits = 6)
  }
}

result <- do.call(rbind, rows)
write.csv(result, file.path(root, "outputs", "native_tsmr_method_benchmark.csv"), row.names = FALSE)
markdown <- c(
  "# Native TwoSampleMR method benchmarks", "",
  sprintf("IL6 fixture: %d SNPs; nboot=%d; fastMR threads=%d; seed=%d.",
          nrow(d), nboot, threads, seed),
  "The native comparator is the installed TwoSampleMR R workflow, timed one method at a time.",
  "Bootstrap SE differences are Monte Carlo differences when comparing a shared fastMR grid layout with native per-pair streams.",
  "", "| scenario | method | pairs | fastMR s | TwoSampleMR s | speedup | max beta delta | median SE delta |", "|---|---|---:|---:|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(result))) {
  x <- result[i, ]
  markdown <- c(markdown, sprintf(
    "| %s | %s | %d | %.3f | %.3f | %.2fx | %.3e | %.3e |",
    x$scenario, x$method_code, x$pairs, x$fastMR_seconds,
    x$TwoSampleMR_seconds, x$speedup, x$max_abs_beta_difference,
    x$median_abs_se_difference))
}
writeLines(markdown, file.path(root, "outputs", "native_tsmr_method_benchmark.md"))
print(result, row.names = FALSE, digits = 6)
