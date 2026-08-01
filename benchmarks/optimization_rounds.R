args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]]) else getwd()
.libPaths(c(file.path(root, ".local", "Rlib"), .libPaths()))
suppressPackageStartupMessages(library(fastMR))

d <- read.delim(file.path(root, "inst/extdata", "il6_crp_primary_100.tsv"),
                check.names = FALSE, stringsAsFactors = FALSE)
n <- 50L
grid <- list(
  exposure_beta = matrix(rep(d$beta.exposure, n), nrow = n, byrow = TRUE),
  outcome_beta = matrix(rep(d$beta.outcome, n), nrow = n, byrow = TRUE),
  exposure_se = matrix(rep(d$se.exposure, n), nrow = n, byrow = TRUE),
  outcome_se = matrix(rep(d$se.outcome, n), nrow = n, byrow = TRUE)
)
methods <- c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode")
timed <- function(repeats = 3L) {
  invisible(fast_mr_grid(grid$exposure_beta, grid$outcome_beta, grid$exposure_se,
                         grid$outcome_se, methods = methods, nboot = 100L,
                         seed = 20260803, threads = 10L))
  values <- vapply(seq_len(repeats), function(i) {
    gc(FALSE)
    system.time(fast_mr_grid(grid$exposure_beta, grid$outcome_beta,
      grid$exposure_se, grid$outcome_se, methods = methods, nboot = 100L,
      seed = 20260803, threads = 10L))["elapsed"]
  }, numeric(1))
  median(values)
}

history <- data.frame(
  stage = c("baseline", "round_1", "round_2", "round_3", "round_4", "round_5", "final_validation"),
  change = c(
    "Pre-cycle O2 build, nested mixed-grid results, original mode scan",
    "Compiler -O3 (rejected: R CMD check portability warning)",
    "Static fallback thread chunks (rejected)",
    "Flat mixed-grid Result storage",
    "Linearized exact-mode interpolation positions",
    "Fused simple/weighted exact-mode max scans",
    "Final retained implementation"
  ),
  seconds = c(1.544, 1.523, 1.728, 1.538, 1.534, 1.498, NA_real_),
  accepted = c(TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
)
history$delta_from_baseline_pct <- 100 * (history$seconds / history$seconds[[1L]] - 1)
history$seconds[[nrow(history)]] <- timed()
history$delta_from_baseline_pct[[nrow(history)]] <-
  100 * (history$seconds[[nrow(history)]] / history$seconds[[1L]] - 1)
dir.create(file.path(root, "outputs"), showWarnings = FALSE, recursive = TRUE)
write.csv(history, file.path(root, "outputs", "optimization_rounds.csv"), row.names = FALSE)
writeLines(c(
  "# Five-round fastMR optimization cycle", "",
  "Workload: IL6/CRP fixture, 50x50 exposure/outcome grid, five methods, nboot=100, ten requested threads.",
  "The rejected round is retained to show that the slower scheduler was not kept.", "",
  "| stage | change | seconds | delta vs baseline | accepted |",
  "|---|---|---:|---:|:---:|",
  vapply(seq_len(nrow(history)), function(i) sprintf(
    "| %s | %s | %.3f | %.2f%% | %s |", history$stage[i], history$change[i],
    history$seconds[i], history$delta_from_baseline_pct[i], history$accepted[i]), character(1)),
  "", sprintf("Final validation median: %.3f seconds (%.2f%% versus baseline).",
               history$seconds[[nrow(history)]], history$delta_from_baseline_pct[[nrow(history)]])
), file.path(root, "outputs", "optimization_rounds.md"))
print(history, row.names = FALSE)
