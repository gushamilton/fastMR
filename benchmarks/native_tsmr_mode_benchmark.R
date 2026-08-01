args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]]) else getwd()
.libPaths(c(file.path(root, ".local", "Rlib"),
            "/Users/fergushamilton/projects/twosamplemr-fast/.local/Rlib", .libPaths()))
suppressPackageStartupMessages({ library(fastMR); library(TwoSampleMR) })
d <- read.delim(file.path(root, "inst", "extdata", "il6_crp_primary_100.tsv"),
                check.names = FALSE, stringsAsFactors = FALSE)
specs <- data.frame(
  method = c("simple_mode", "weighted_mode"),
  native = c("mr_simple_mode", "mr_weighted_mode"), stringsAsFactors = FALSE
)
scenarios <- data.frame(
  scenario = c("balanced_50x50", "one_exposure_250_outcomes",
               "one_outcome_250_exposures", "ten_exposures_100_outcomes",
               "hundred_exposures_ten_outcomes"),
  exposures = c(50L, 1L, 250L, 10L, 100L),
  outcomes = c(50L, 250L, 1L, 100L, 10L), stringsAsFactors = FALSE
)
nboot <- 100L; seed <- 20260818L; threads <- 5L
make_grid <- function(ne, no) list(
  exposure_beta = matrix(rep(d[["beta.exposure"]], ne), nrow = ne, byrow = TRUE) *
    (1 + (seq_len(ne) - (ne + 1) / 2) * 0.001),
  outcome_beta = matrix(rep(d[["beta.outcome"]], no), nrow = no, byrow = TRUE) *
    (1 + (seq_len(no) - (no + 1) / 2) * 0.001),
  exposure_se = matrix(rep(d[["se.exposure"]], ne), nrow = ne, byrow = TRUE),
  outcome_se = matrix(rep(d[["se.outcome"]], no), nrow = no, byrow = TRUE)
)
run_native <- function(g, ne, no, method, params) {
  out <- vector("list", ne * no); k <- 0L
  for (i in seq_len(ne)) for (j in seq_len(no)) {
    k <- k + 1L
    z <- data.frame(SNP = d[["SNP"]], beta.exposure = g$exposure_beta[i, ],
                    beta.outcome = g$outcome_beta[j, ], se.exposure = g$exposure_se[i, ],
                    se.outcome = g$outcome_se[j, ], id.exposure = "E", id.outcome = "O",
                    exposure = "E", outcome = "O", mr_keep = TRUE)
    out[[k]] <- suppressMessages(suppressWarnings(
      TwoSampleMR::mr(z, method_list = method, parameters = params)))
  }
  do.call(rbind, out)
}
params <- modifyList(TwoSampleMR::default_parameters(), list(nboot = nboot))
rows <- list(); k <- 0L
for (scenario_index in seq_len(nrow(scenarios))) {
  s <- scenarios[scenario_index, ]; g <- make_grid(s$exposures, s$outcomes)
  pairs <- s$exposures * s$outcomes
  for (method_index in seq_len(nrow(specs))) {
    k <- k + 1L; spec <- specs[method_index, ]
    invisible(fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se,
                           g$outcome_se, methods = spec$method, nboot = nboot,
                           seed = seed, threads = threads))
    started <- proc.time()[["elapsed"]]
    fast <- fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se,
                         g$outcome_se, methods = spec$method, nboot = nboot,
                         seed = seed, threads = threads)
    fast_seconds <- proc.time()[["elapsed"]] - started
    set.seed(seed); started <- proc.time()[["elapsed"]]
    native <- run_native(g, s$exposures, s$outcomes, spec$native, params)
    native_seconds <- proc.time()[["elapsed"]] - started
    rows[[k]] <- data.frame(
      scenario = s$scenario, method = spec$method, pairs = pairs,
      fastMR_seconds = fast_seconds, TwoSampleMR_seconds = native_seconds,
      speedup = native_seconds / fast_seconds,
      max_abs_beta_difference = max(abs(fast$b - native$b), na.rm = TRUE),
      max_abs_se_difference = max(abs(fast$se - native$se), na.rm = TRUE),
      max_abs_pval_difference = max(abs(fast$pval - native$pval), na.rm = TRUE)
    )
    print(rows[[k]], row.names = FALSE)
  }
}
result <- do.call(rbind, rows)
dir.create(file.path(root, "outputs"), showWarnings = FALSE, recursive = TRUE)
write.csv(result, file.path(root, "outputs", "native_tsmr_mode_benchmark.csv"), row.names = FALSE)
lines <- c("# Native TwoSampleMR exact mode benchmark", "",
  sprintf("IL6 fixture: %d SNPs; nboot=%d; fastMR threads=%d; seed=%d.",
          nrow(d), nboot, threads, seed), "",
  "| scenario | method | pairs | fastMR s | TwoSampleMR s | speedup | max beta delta | max SE delta |",
  "|---|---|---:|---:|---:|---:|---:|---:|")
for (i in seq_len(nrow(result))) lines <- c(lines, sprintf(
  "| %s | %s | %d | %.3f | %.3f | %.2fx | %.3e | %.3e |",
  result$scenario[i], result$method[i], result$pairs[i], result$fastMR_seconds[i],
  result$TwoSampleMR_seconds[i], result$speedup[i],
  result$max_abs_beta_difference[i], result$max_abs_se_difference[i]))
writeLines(lines, file.path(root, "outputs", "native_tsmr_mode_benchmark.md"))
