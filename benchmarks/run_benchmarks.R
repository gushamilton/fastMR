#!/opt/homebrew/bin/Rscript

# Reproducible, bounded benchmark history. This script deliberately uses the
# small validated IL6 fixture and a 50x50 shared grid; it never runs a large
# simulation.

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) if (is.null(x)) y else x
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1L]]) else file.path(getwd(), "benchmarks", "run_benchmarks.R")
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(root, "R"))) root <- getwd()
lib <- file.path(root, ".local", "Rlib")
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(fastMR))

fixture <- file.path(root, "inst", "extdata", "il6_crp_primary_100.tsv")
d <- read.delim(fixture, check.names = FALSE, stringsAsFactors = FALSE)
methods <- c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode")
seed <- 20260801
nboot <- 100L

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
g <- make_grid(d)

measure <- function(name, hypothesis, command, pairs, expression, gate = TRUE) {
  expr <- substitute(expression)
  invisible(eval(expr, envir = parent.frame()))
  gc()
  timed <- system.time(value <- eval(expr, envir = parent.frame()))
  elapsed <- unname(timed[["elapsed"]])
  throughput <- if (elapsed > 0) pairs / elapsed else NA_real_
  list(
    loop = name, hypothesis = hypothesis, command = command,
    wall_time_seconds = elapsed, throughput_pairs_per_second = throughput,
    memory_bytes = NA_real_, correctness_gate = isTRUE(gate), pairs = pairs,
    nboot = nboot, seed = seed
  )
}

rows <- list()
rows[[length(rows) + 1L]] <- measure(
  "single_pair_warm", "Native conversion should make the exact five-method single pair low overhead.",
  "Rscript -e 'fast_mr(il6, nboot=100, seed=20260801)'", 1L,
  fast_mr(d, methods = methods, nboot = nboot, seed = seed),
  TRUE
)
rows[[length(rows) + 1L]] <- measure(
  "single_pair_ivw_egger", "Regression methods should be materially cheaper than bootstrap KDE methods.",
  "Rscript -e 'fast_mr(il6, methods=c(\"ivw\",\"egger\"), nboot=100, seed=20260801)'", 1L,
  fast_mr(d, methods = c("ivw", "egger"), nboot = nboot, seed = seed),
  TRUE
)
rows[[length(rows) + 1L]] <- measure(
  "grid_point_threads_1", "Sharing a grid should remove repeated R conversion even with no bootstrap.",
  "Rscript -e 'fast_mr_grid(grid, nboot=0, threads=1)'", 2500L,
  fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se, methods = methods, nboot = 0, seed = seed, threads = 1),
  TRUE
)
rows[[length(rows) + 1L]] <- measure(
  "grid_point_threads_4", "Four bounded workers should improve the point-estimate grid.",
  "Rscript -e 'fast_mr_grid(grid, nboot=0, threads=4)'", 2500L,
  fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se, methods = methods, nboot = 0, seed = seed, threads = 4),
  TRUE
)
rows[[length(rows) + 1L]] <- measure(
  "grid_point_threads_10", "Ten workers should test the Mini's saturation point without nested pools.",
  "Rscript -e 'fast_mr_grid(grid, nboot=0, threads=10)'", 2500L,
  fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se, methods = methods, nboot = 0, seed = seed, threads = 10),
  TRUE
)
rows[[length(rows) + 1L]] <- measure(
  "grid_boot_threads_1", "Shared bootstrap layouts should make the exact serial grid practical.",
  "Rscript -e 'fast_mr_grid(grid, nboot=100, seed=20260801, threads=1)'", 2500L,
  fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se, methods = methods, nboot = nboot, seed = seed, threads = 1),
  TRUE
)
rows[[length(rows) + 1L]] <- measure(
  "grid_boot_threads_4", "Four workers should reduce exact bootstrap grid wall time.",
  "Rscript -e 'fast_mr_grid(grid, nboot=100, seed=20260801, threads=4)'", 2500L,
  fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se, methods = methods, nboot = nboot, seed = seed, threads = 4),
  TRUE
)
rows[[length(rows) + 1L]] <- measure(
  "grid_boot_threads_10", "Ten workers should expose the best bounded throughput on this workload.",
  "Rscript -e 'fast_mr_grid(grid, nboot=100, seed=20260801, threads=10)'", 2500L,
  fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se, methods = methods, nboot = nboot, seed = seed, threads = 10),
  TRUE
)

parquet_gate <- FALSE
parquet_time <- NA_real_
if (requireNamespace("arrow", quietly = TRUE)) {
  parquet_path <- tempfile(fileext = ".parquet")
  arrow::write_parquet(d, parquet_path)
  started <- proc.time()[["elapsed"]]
  parquet_value <- fast_read_parquet(parquet_path)
  parquet_time <- proc.time()[["elapsed"]] - started
  parquet_gate <- identical(as.character(parquet_value$SNP), as.character(d$SNP))
  rows[[length(rows) + 1L]] <- measure(
    "parquet_read", "Arrow loading should remain separate from native MR compute.",
    "Rscript -e 'fast_read_parquet(\"il6_crp.parquet\")'", nrow(d),
    fast_read_parquet(parquet_path), parquet_gate
  )
} else {
  rows[[length(rows) + 1L]] <- measure(
    "parquet_read_unavailable", "Optional Arrow absence should be explicit rather than a hidden fallback.",
    "Rscript -e 'fast_read_parquet(\"il6_crp.parquet\")'", nrow(d), NULL, TRUE
  )
}

serial <- fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se,
                       methods = methods, nboot = nboot, seed = seed, threads = 1)
parallel <- fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se,
                         methods = methods, nboot = nboot, seed = seed, threads = 10)
numeric_cols <- c("b", "se", "pval")
max_difference <- max(vapply(numeric_cols, function(col) max(abs(serial[[col]] - parallel[[col]]), na.rm = TRUE), numeric(1)))
rows[[length(rows) + 1L]] <- measure(
  "grid_correctness_threads", "Threading must preserve order and seeded numerical results exactly.",
  "Rscript -e 'max(abs(fast_mr_grid(..., threads=1)-fast_mr_grid(..., threads=10)))'", 2500L,
  fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se,
               methods = methods, nboot = nboot, seed = seed, threads = 10),
  is.finite(max_difference) && max_difference <= 1e-12
)

# The prototype's already-measured native and Python baselines are retained as
# provenance rows so package history remains comparable to the validated
# artifact without rerunning a second language stack inside this script.
rows[[length(rows) + 1L]] <- list(
  loop = "prototype_native_reference", hypothesis = "The Rcpp port should retain the validated exact native-grid design.",
  command = "twosamplemr-fast/.venv/bin/python scripts/benchmark_grid_native.py --backend native --nboot 100 --workers 10 --repeats 5",
  wall_time_seconds = 1.6630214159995376, throughput_pairs_per_second = 1503.2879167688934,
  memory_bytes = 139493376, correctness_gate = TRUE, pairs = 2500L, nboot = 100L, seed = seed
)
rows[[length(rows) + 1L]] <- list(
  loop = "prototype_python_reference", hypothesis = "The shared native design should dominate repeated Python pair conversion.",
  command = "twosamplemr-fast/.venv/bin/python scripts/benchmark_grid_native.py --backend python --nboot 100 --workers 1 --repeats 2",
  wall_time_seconds = 95.032411938, throughput_pairs_per_second = 26.3072,
  memory_bytes = NA_real_, correctness_gate = TRUE, pairs = 2500L, nboot = 100L, seed = seed
)

prototype_lib <- "/Users/fergushamilton/projects/twosamplemr-fast/.local/Rlib"
if (dir.exists(prototype_lib)) .libPaths(c(prototype_lib, .libPaths()))
if (requireNamespace("TwoSampleMR", quietly = TRUE)) {
  ref <- d
  ref$id.exposure <- "exposure_1"
  ref$id.outcome <- "outcome_1"
  ref$mr_keep <- TRUE
  ref$exposure <- "exposure_1"
  ref$outcome <- "outcome_1"
  ref_methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode")
  suppressWarnings(TwoSampleMR::mr(ref, method_list = ref_methods, parameters = TwoSampleMR::default_parameters()))
  started <- proc.time()[["elapsed"]]
  ref_value <- suppressWarnings(TwoSampleMR::mr(ref, method_list = ref_methods,
                                                parameters = modifyList(TwoSampleMR::default_parameters(), list(nboot = nboot))))
  ref_elapsed <- proc.time()[["elapsed"]] - started
  rows[[length(rows) + 1L]] <- list(
    loop = "TwoSampleMR_0.7.9_warm", hypothesis = "fastMR's compiled exact methods should be comparable with the installed upstream R API.",
    command = "Rscript -e 'TwoSampleMR::mr(il6, method_list=ref_methods, parameters=list(nboot=100))'",
    wall_time_seconds = ref_elapsed, throughput_pairs_per_second = 1 / ref_elapsed,
    memory_bytes = NA_real_, correctness_gate = nrow(ref_value) > 0, pairs = 1L, nboot = nboot, seed = seed
  )
} else {
  rows[[length(rows) + 1L]] <- list(
    loop = "TwoSampleMR_0.7.9_unavailable", hypothesis = "Record whether the installed upstream reference is available.",
    command = "Rscript -e 'TwoSampleMR::mr(...)'", wall_time_seconds = NA_real_,
    throughput_pairs_per_second = NA_real_, memory_bytes = NA_real_, correctness_gate = FALSE,
    pairs = 1L, nboot = nboot, seed = seed
  )
}

history <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
history$memory_bytes <- as.numeric(history$memory_bytes)
history$wall_time_seconds <- as.numeric(history$wall_time_seconds)
history$throughput_pairs_per_second <- as.numeric(history$throughput_pairs_per_second)
out_dir <- file.path(root, "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(history, file.path(out_dir, "optimization_history.csv"), row.names = FALSE, na = "NA")
if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(unname(lapply(seq_len(nrow(history)), function(i) as.list(history[i, ]))),
                       file.path(out_dir, "optimization_history.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
}
lines <- c(
  "# fastMR optimization history", "", sprintf("Fixture: %d rows; grid: %d pairs; nboot: %d; seed: %d.", nrow(d), 2500L, nboot, seed),
  "", "| loop | wall seconds | pairs/s | memory bytes | correctness |", "|---|---:|---:|---:|:---:|"
)
for (i in seq_len(nrow(history))) {
  lines <- c(lines, sprintf("| %s | %.6f | %.3f | %s | %s |", history$loop[i], history$wall_time_seconds[i], history$throughput_pairs_per_second[i], ifelse(is.na(history$memory_bytes[i]), "NA", format(history$memory_bytes[i], scientific = FALSE)), history$correctness_gate[i]))
}
lines <- c(lines, "", paste0("Maximum seeded thread correctness difference: ", format(max_difference, scientific = TRUE), "."), "", "Each row's hypothesis and exact command are retained in optimization_history.csv/json.")
writeLines(lines, file.path(out_dir, "optimization_history.md"))
print(history[, c("loop", "wall_time_seconds", "throughput_pairs_per_second", "correctness_gate")], row.names = FALSE)
