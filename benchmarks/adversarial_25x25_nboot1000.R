args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("usage: adversarial_25x25_nboot1000.R ROOT OUTPUT_DIR [FAST_LIBRARY]", call. = FALSE)
}

root <- normalizePath(args[[1L]], mustWork = TRUE)
out_dir <- args[[2L]]
fast_lib <- if (length(args) >= 3L) normalizePath(args[[3L]], mustWork = TRUE) else
  file.path(root, ".local-adversarial", "Rlib")
tsmr_lib <- "/Users/fergushamilton/projects/twosamplemr-fast/.local/Rlib"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(fast_lib, tsmr_lib, .libPaths()))

suppressPackageStartupMessages({
  library(fastMR)
  library(TwoSampleMR)
})

seed <- 20260802L
n_exposure <- 25L
n_outcome <- 25L
n_snp <- 82L
nboot <- 1000L
fast_threads <- as.integer(Sys.getenv("FASTMR_THREADS", "10"))
methods <- c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode")
native_methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median",
                    "mr_simple_mode", "mr_weighted_mode")
method_map <- c("Inverse variance weighted" = "ivw", "MR Egger" = "egger",
                "Weighted median" = "weighted_median", "Simple mode" = "simple_mode",
                "Weighted mode" = "weighted_mode")

timed <- function(expr) {
  started <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, seconds = proc.time()[["elapsed"]] - started)
}

# Deliberately different from the package's IL6/CRP fixture: heavy-tailed
# effects, heteroskedastic standard errors, row/column scale differences, and
# mixed signs. The effect magnitudes are kept away from exact zero so that
# bootstrap instability is a property of the method, not a malformed input.
make_simulation <- function() {
  set.seed(seed)
  col_scale <- exp(seq(log(0.75), log(1.25), length.out = n_snp))
  exp_scale <- seq(0.82, 1.18, length.out = n_exposure)
  out_scale <- seq(0.86, 1.14, length.out = n_outcome)
  exposure_beta <- matrix(rt(n_exposure * n_snp, df = 5), nrow = n_exposure)
  outcome_beta <- matrix(rt(n_outcome * n_snp, df = 6), nrow = n_outcome)
  exposure_beta <- sweep(exposure_beta, 1L, exp_scale, "*")
  outcome_beta <- sweep(outcome_beta, 1L, out_scale, "*")
  exposure_beta <- sweep(exposure_beta, 2L, col_scale, "*") * 0.025
  outcome_beta <- sweep(outcome_beta, 2L, rev(col_scale), "*") * 0.030
  exposure_beta <- exposure_beta + matrix(rnorm(n_exposure * n_snp, 0, 0.001),
                                           nrow = n_exposure)
  outcome_beta <- outcome_beta + matrix(rnorm(n_outcome * n_snp, 0, 0.0015),
                                         nrow = n_outcome)
  exposure_se <- matrix(exp(rnorm(n_exposure * n_snp, log(0.012), 0.30)),
                        nrow = n_exposure)
  outcome_se <- matrix(exp(rnorm(n_outcome * n_snp, log(0.018), 0.35)),
                       nrow = n_outcome)
  list(exposure_beta = exposure_beta, outcome_beta = outcome_beta,
       exposure_se = exposure_se, outcome_se = outcome_se)
}

sim <- make_simulation()
n_pairs <- n_exposure * n_outcome
expected_rows <- n_pairs * length(methods)

# Warm up compilation/registration before timing.
invisible(fast_mr_grid(sim$exposure_beta[1:2, , drop = FALSE],
                       sim$outcome_beta[1:2, , drop = FALSE],
                       sim$exposure_se[1:2, , drop = FALSE],
                       sim$outcome_se[1:2, , drop = FALSE],
                       methods = methods, nboot = 0L, seed = seed, threads = 1L))

fast_parallel_timed <- timed(fast_mr_grid(
  sim$exposure_beta, sim$outcome_beta, sim$exposure_se, sim$outcome_se,
  methods = methods, nboot = nboot, seed = seed, threads = fast_threads
))
fast_parallel <- fast_parallel_timed$value

fast_serial_timed <- timed(fast_mr_grid(
  sim$exposure_beta, sim$outcome_beta, sim$exposure_se, sim$outcome_se,
  methods = methods, nboot = nboot, seed = seed, threads = 1L
))
fast_serial <- fast_serial_timed$value

if (nrow(fast_parallel) != expected_rows || nrow(fast_serial) != expected_rows) {
  stop(sprintf("fastMR returned %d/%d rows; expected %d", nrow(fast_parallel),
              nrow(fast_serial), expected_rows), call. = FALSE)
}

run_native_grid <- function(g) {
  params <- modifyList(TwoSampleMR::default_parameters(), list(nboot = nboot))
  snps <- data.frame(SNP = paste0("adv", seq_len(ncol(g$exposure_beta))),
                    stringsAsFactors = FALSE)
  results <- vector("list", n_pairs)
  k <- 0L
  for (i in seq_len(nrow(g$exposure_beta))) {
    for (j in seq_len(nrow(g$outcome_beta))) {
      k <- k + 1L
      dat <- cbind(snps,
        beta.exposure = g$exposure_beta[i, ],
        beta.outcome = g$outcome_beta[j, ],
        se.exposure = g$exposure_se[i, ],
        se.outcome = g$outcome_se[j, ],
        id.exposure = paste0("E", i), id.outcome = paste0("O", j),
        exposure = paste0("E", i), outcome = paste0("O", j),
        mr_keep = TRUE)
      results[[k]] <- suppressWarnings(suppressMessages(TwoSampleMR::mr(
        dat, method_list = native_methods, parameters = params
      )))
    }
  }
  out <- do.call(rbind, results)
  out$pair <- rep(seq_len(n_pairs), each = length(native_methods))
  out
}

set.seed(seed)
native_timed <- timed(run_native_grid(sim))
native <- native_timed$value
if (nrow(native) != expected_rows) {
  stop(sprintf("TwoSampleMR returned %d rows; expected %d", nrow(native), expected_rows),
       call. = FALSE)
}

add_pair <- function(x) {
  x$pair <- (as.integer(x$exposure_index) - 1L) * n_outcome +
    as.integer(x$outcome_index)
  x$method_code <- as.character(x$method_code)
  x
}
fast_parallel <- add_pair(fast_parallel)
fast_serial <- add_pair(fast_serial)
native$method_code <- unname(method_map[as.character(native$method)])

join_fields <- c("pair", "method_code")
fast_compare <- merge(
  fast_parallel[, c(join_fields, "b", "se", "pval")],
  fast_serial[, c(join_fields, "b", "se", "pval")],
  by = join_fields, suffixes = c("_parallel", "_serial"), sort = FALSE
)
fast_compare$abs_beta_delta <- abs(fast_compare$b_parallel - fast_compare$b_serial)
fast_compare$abs_se_delta <- abs(fast_compare$se_parallel - fast_compare$se_serial)
fast_compare$abs_pval_delta <- abs(fast_compare$pval_parallel - fast_compare$pval_serial)

native_compare <- merge(
  fast_parallel[, c(join_fields, "b", "se", "pval")],
  native[, c(join_fields, "b", "se", "pval")],
  by = join_fields, suffixes = c("_fast", "_native"), sort = FALSE
)
native_compare$abs_beta_delta <- abs(native_compare$b_fast - native_compare$b_native)
native_compare$abs_se_delta <- abs(native_compare$se_fast - native_compare$se_native)
native_compare$abs_pval_delta <- abs(native_compare$pval_fast - native_compare$pval_native)
native_compare$relative_se_delta <- native_compare$abs_se_delta /
  pmax(abs(native_compare$se_native), 1e-12)

finite_max <- function(x) if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_
point_rows <- native_compare[native_compare$method_code %in% methods &
                               is.finite(native_compare$b_fast) &
                               is.finite(native_compare$b_native), , drop = FALSE]
point_max_beta_delta <- finite_max(point_rows$abs_beta_delta)
if (!is.finite(point_max_beta_delta) || point_max_beta_delta > 1e-6) {
  stop(sprintf("point-estimate parity failed: max beta delta %.6g", point_max_beta_delta),
       call. = FALSE)
}

fast_max_delta <- max(c(fast_compare$abs_beta_delta,
                        fast_compare$abs_se_delta,
                        fast_compare$abs_pval_delta), na.rm = TRUE)
if (!is.finite(fast_max_delta) || fast_max_delta > 1e-10) {
  stop(sprintf("seeded thread reproducibility failed: max delta %.6g", fast_max_delta),
       call. = FALSE)
}

method_summary <- do.call(rbind, lapply(methods, function(method) {
  x <- native_compare[native_compare$method_code == method, , drop = FALSE]
  data.frame(
    method = method,
    pairs = nrow(x),
    max_abs_beta_delta = finite_max(x$abs_beta_delta),
    median_abs_beta_delta = median(x$abs_beta_delta, na.rm = TRUE),
    max_abs_se_delta = finite_max(x$abs_se_delta),
    p95_abs_se_delta = unname(quantile(x$abs_se_delta, 0.95, na.rm = TRUE)),
    median_relative_se_delta = median(x$relative_se_delta, na.rm = TRUE),
    p95_relative_se_delta = unname(quantile(x$relative_se_delta, 0.95, na.rm = TRUE)),
    median_abs_pval_delta = median(x$abs_pval_delta, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

# Small-input and malformed-input checks are kept separate from the timed grid.
edge_cases <- list()
record_edge <- function(name, passed, detail) {
  edge_cases[[length(edge_cases) + 1L]] <<- data.frame(
    check = name, passed = isTRUE(passed), detail = detail,
    stringsAsFactors = FALSE
  )
}

for (snp_count in c(1L, 2L, 3L, 7L)) {
  x <- sim$exposure_beta[1, seq_len(snp_count), drop = FALSE]
  y <- sim$outcome_beta[1, seq_len(snp_count), drop = FALSE]
  sx <- sim$exposure_se[1, seq_len(snp_count), drop = FALSE]
  sy <- sim$outcome_se[1, seq_len(snp_count), drop = FALSE]
  result <- tryCatch(fast_mr_grid(x, y, sx, sy, methods = methods,
                                  nboot = 0L, seed = seed, threads = 1L),
                     error = function(e) e)
  passed <- !inherits(result, "error") && nrow(result) == length(methods)
  record_edge(paste0("small_input_", snp_count, "_snps"), passed,
              if (inherits(result, "error")) conditionMessage(result) else
                paste(nrow(result), "rows"))
}

bad <- sim$exposure_beta[1:2, , drop = FALSE]
bad[1, 1] <- NA_real_
bad_result <- tryCatch(fast_mr_grid(bad, sim$outcome_beta[1:2, , drop = FALSE],
                                    sim$exposure_se[1:2, , drop = FALSE],
                                    sim$outcome_se[1:2, , drop = FALSE],
                                    methods = "ivw", nboot = 0L),
                       error = function(e) e)
record_edge("non_finite_input_rejected", inherits(bad_result, "error"),
            if (inherits(bad_result, "error")) conditionMessage(bad_result) else
              "unexpectedly accepted NA input")

mismatch_result <- tryCatch(fast_mr_grid(sim$exposure_beta[1:2, , drop = FALSE],
                                         sim$outcome_beta[1:2, -1, drop = FALSE],
                                         sim$exposure_se[1:2, , drop = FALSE],
                                         sim$outcome_se[1:2, -1, drop = FALSE],
                                         methods = "ivw", nboot = 0L),
                            error = function(e) e)
record_edge("mismatched_snp_count_rejected", inherits(mismatch_result, "error"),
            if (inherits(mismatch_result, "error")) conditionMessage(mismatch_result) else
              "unexpectedly accepted mismatched matrices")

harm_exposure <- data.frame(
  SNP = paste0("h", 1:4), beta.exposure = c(.10, .12, .08, .15),
  se.exposure = rep(.02, 4), effect_allele.exposure = c("A", "C", "A", "A"),
  other_allele.exposure = c("G", "T", "C", "G"), eaf.exposure = c(.20, .30, .45, .49),
  id.exposure = "E", stringsAsFactors = FALSE
)
harm_outcome <- data.frame(
  SNP = paste0("h", 1:4), beta.outcome = c(.02, -.03, .01, .04),
  se.outcome = rep(.01, 4), effect_allele.outcome = c("G", "G", "T", "A"),
  other_allele.outcome = c("A", "C", "G", "C"), eaf.outcome = c(.80, .70, .55, .51),
  id.outcome = "O", stringsAsFactors = FALSE
)
harm_result <- tryCatch(fast_harmonise_data(harm_exposure, harm_outcome, action = 2),
                        error = function(e) e)
harm_passed <- !inherits(harm_result, "error") &&
  nrow(harm_result) == 4L && any(harm_result$mr_keep %in% TRUE) &&
  any(harm_result$beta.outcome < 0)
record_edge("flipped_alleles_harmonised", harm_passed,
            if (inherits(harm_result, "error")) conditionMessage(harm_result) else
              paste("rows", nrow(harm_result), "kept", sum(harm_result$mr_keep, na.rm = TRUE)))

edge_table <- do.call(rbind, edge_cases)
if (!all(edge_table$passed)) {
  stop("one or more adversarial edge checks failed", call. = FALSE)
}

summary_row <- data.frame(
  exposures = n_exposure, outcomes = n_outcome, snps = n_snp,
  pairs = n_pairs, methods = length(methods), nboot = nboot,
  fast_threads = fast_threads, expected_rows = expected_rows,
  fast_rows = nrow(fast_parallel), native_rows = nrow(native),
  fast_parallel_seconds = fast_parallel_timed$seconds,
  fast_serial_seconds = fast_serial_timed$seconds,
  native_seconds = native_timed$seconds,
  speedup_vs_native = native_timed$seconds / fast_parallel_timed$seconds,
  max_thread_repro_delta = fast_max_delta,
  max_native_point_beta_delta = point_max_beta_delta,
  edge_checks = nrow(edge_table), edge_checks_passed = sum(edge_table$passed),
  stringsAsFactors = FALSE
)

write.csv(summary_row, file.path(out_dir, "summary.csv"), row.names = FALSE)
write.csv(method_summary, file.path(out_dir, "method_summary.csv"), row.names = FALSE)
write.csv(fast_compare, file.path(out_dir, "thread_reproducibility.csv"), row.names = FALSE)
write.csv(native_compare, file.path(out_dir, "native_parity.csv"), row.names = FALSE)
write.csv(edge_table, file.path(out_dir, "edge_checks.csv"), row.names = FALSE)

report <- c(
  "# Adversarial 25 x 25 fastMR validation",
  "",
  sprintf("Simulation: %d exposures x %d outcomes, %d SNPs, five main methods, nboot=%d.",
          n_exposure, n_outcome, n_snp, nboot),
  sprintf("Rows: %d fastMR and %d native TwoSampleMR; expected %d.",
          nrow(fast_parallel), nrow(native), expected_rows),
  sprintf("Timings: fastMR %d threads %.3fs; fastMR 1 thread %.3fs; native TwoSampleMR %.3fs.",
          fast_threads, fast_parallel_timed$seconds, fast_serial_timed$seconds,
          native_timed$seconds),
  sprintf("Speedup: %.2fx versus native TwoSampleMR.", summary_row$speedup_vs_native),
  sprintf("Seeded thread reproducibility maximum delta: %.3e.", fast_max_delta),
  sprintf("Maximum native point-estimate beta delta: %.3e.", point_max_beta_delta),
  sprintf("Edge checks: %d/%d passed.", sum(edge_table$passed), nrow(edge_table)),
  "",
  "Bootstrap SE and p-value differences are not treated as point-estimate failures: fastMR and TwoSampleMR use independent bootstrap implementations/streams. The method_summary.csv file reports their per-method distributions.",
  "",
  "The input includes deliberately varied signs, scales, heavy-tailed effects, and heteroskedastic standard errors. Separate edge checks cover 1/2/3/7-SNP inputs, rejected non-finite and mismatched matrices, and flipped alleles through fast_harmonise_data(action = 2)."
)
writeLines(report, file.path(out_dir, "README.md"))
cat(paste(report, collapse = "\n"), "\n")

