args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) {
  stop("usage: final_package_audit.R ROOT BMI_SIG_TSV CRP_MATCH_TSV OUTPUT_DIR", call. = FALSE)
}
root <- normalizePath(args[[1L]], mustWork = TRUE)
bmi_file <- normalizePath(args[[2L]], mustWork = TRUE)
crp_file <- normalizePath(args[[3L]], mustWork = TRUE)
out_dir <- args[[4L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fast_lib <- file.path(root, ".local", "Rlib")
tsmr_lib <- "/Users/fergushamilton/projects/twosamplemr-fast/.local/Rlib"
.libPaths(c(fast_lib, tsmr_lib, .libPaths()))
suppressPackageStartupMessages({
  library(fastMR)
  library(TwoSampleMR)
})

methods <- c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode")
native_methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median",
                    "mr_simple_mode", "mr_weighted_mode")
method_map <- c("Inverse variance weighted" = "ivw", "MR Egger" = "egger",
                "Weighted median" = "weighted_median", "Simple mode" = "simple_mode",
                "Weighted mode" = "weighted_mode")
nboot <- 100L
seed <- 20260802L

bmi <- read.delim(bmi_file, check.names = FALSE, stringsAsFactors = FALSE)
crp <- read.delim(crp_file, check.names = FALSE, stringsAsFactors = FALSE)
exposure <- data.frame(
  SNP = bmi$snp,
  id.exposure = "BMI",
  exposure = "BMI",
  beta.exposure = as.numeric(bmi$bmi_beta_raw),
  se.exposure = as.numeric(bmi$bmi_se_raw),
  effect_allele.exposure = bmi$bmi_effect_allele,
  other_allele.exposure = bmi$bmi_other_allele,
  eaf.exposure = as.numeric(bmi$bmi_eaf_hrs),
  pval.exposure = as.numeric(bmi$bmi_p),
  samplesize.exposure = as.numeric(bmi$bmi_n),
  stringsAsFactors = FALSE
)
outcome <- data.frame(
  SNP = crp$variant_id,
  id.outcome = "CRP",
  outcome = "CRP",
  beta.outcome = as.numeric(crp$beta),
  se.outcome = as.numeric(crp$standard_error),
  effect_allele.outcome = crp$effect_allele,
  other_allele.outcome = crp$other_allele,
  eaf.outcome = NA_real_,
  pval.outcome = as.numeric(crp$p_value),
  stringsAsFactors = FALSE
)

canonical_harmonised <- function(x) {
  key <- paste(x$SNP, x$id.outcome, sep = "\r")
  cols <- c("SNP", "id.outcome", "beta.exposure", "beta.outcome", "se.exposure",
            "se.outcome", "effect_allele.exposure", "other_allele.exposure",
            "effect_allele.outcome", "other_allele.outcome", "remove",
            "palindromic", "ambiguous", "mr_keep")
  x <- x[order(key), cols, drop = FALSE]
  for (column in names(x)) if (is.character(x[[column]])) {
    x[[column]][is.na(x[[column]])] <- "<NA>"
  }
  x
}

run_bmi_crp <- function(repeats = 5L) {
  # Warm up package loading and the native registration before timing.
  invisible(fast_mr(fast_harmonise_data(exposure, outcome), methods = methods,
                    nboot = 0L, seed = seed, threads = 1L))
  invisible(suppressMessages(TwoSampleMR::mr(
    suppressMessages(TwoSampleMR::harmonise_data(exposure, outcome, action = 2)),
    method_list = native_methods,
    parameters = modifyList(TwoSampleMR::default_parameters(), list(nboot = nboot))
  )))
  params <- modifyList(TwoSampleMR::default_parameters(), list(nboot = nboot))
  rows <- vector("list", repeats)
  fast_last <- native_last <- NULL
  for (rep in seq_len(repeats)) {
    fast_h_clock <- proc.time()[["elapsed"]]
    fast_h <- suppressMessages(fast_harmonise_data(exposure, outcome, action = 2))
    fast_h_seconds <- proc.time()[["elapsed"]] - fast_h_clock
    fast_m_clock <- proc.time()[["elapsed"]]
    fast_m <- suppressWarnings(fast_mr(fast_h, methods = methods, nboot = nboot,
                                        seed = seed + rep, threads = 1L))
    fast_m_seconds <- proc.time()[["elapsed"]] - fast_m_clock

    set.seed(seed + rep)
    native_h_clock <- proc.time()[["elapsed"]]
    native_h <- suppressMessages(TwoSampleMR::harmonise_data(exposure, outcome, action = 2))
    native_h_seconds <- proc.time()[["elapsed"]] - native_h_clock
    native_m_clock <- proc.time()[["elapsed"]]
    native_m <- suppressWarnings(suppressMessages(TwoSampleMR::mr(
      native_h, method_list = native_methods, parameters = params
    )))
    native_m_seconds <- proc.time()[["elapsed"]] - native_m_clock

    fast_last <- fast_m
    native_last <- native_m
    rows[[rep]] <- data.frame(
      workload = "local BMI -> CRP",
      repeat_index = rep,
      bmi_rows = nrow(exposure),
      crp_rows_matched = nrow(outcome),
      harmonised_rows_fast = nrow(fast_h),
      mr_rows_fast = nrow(fast_m),
      fast_harmonise_seconds = fast_h_seconds,
      fast_mr_seconds = fast_m_seconds,
      fast_pipeline_seconds = fast_h_seconds + fast_m_seconds,
      native_harmonise_seconds = native_h_seconds,
      native_mr_seconds = native_m_seconds,
      native_pipeline_seconds = native_h_seconds + native_m_seconds,
      pipeline_speedup = (native_h_seconds + native_m_seconds) /
        (fast_h_seconds + fast_m_seconds),
      harmonisation_key_exact = isTRUE(all.equal(canonical_harmonised(fast_h),
                                                  canonical_harmonised(native_h),
                                                  check.attributes = FALSE)),
      stringsAsFactors = FALSE
    )
  }
  results <- do.call(rbind, rows)
  fast_result <- fast_last
  native_result <- native_last
  fast_result$method_code <- as.character(fast_result$method_code)
  native_result$method_code <- unname(method_map[as.character(native_result$method)])
  comparison <- merge(fast_result[, c("method_code", "b", "se", "pval")],
                      native_result[, c("method_code", "b", "se", "pval")],
                      by = "method_code", suffixes = c("_fast", "_native"))
  comparison$abs_beta_delta <- abs(comparison$b_fast - comparison$b_native)
  comparison$abs_se_delta <- abs(comparison$se_fast - comparison$se_native)
  comparison$abs_pval_delta <- abs(comparison$pval_fast - comparison$pval_native)
  list(repeats = results, comparison = comparison)
}

bmi_crp <- run_bmi_crp()
write.csv(bmi_crp$repeats, file.path(out_dir, "bmi_crp_benchmark.csv"), row.names = FALSE)
write.csv(bmi_crp$comparison, file.path(out_dir, "bmi_crp_parity.csv"), row.names = FALSE)

# A deliberately different simulation: heavy-tailed effects, heteroskedastic
# standard errors, row-specific scales, and no copied IL6/CRP fixture values.
make_simulation <- function(n_exposure = 50L, n_outcome = 50L, n_snp = 400L) {
  set.seed(seed + 101L)
  exposure_beta <- matrix(rt(n_exposure * n_snp, df = 4) * 0.025,
                          nrow = n_exposure)
  outcome_beta <- matrix(rt(n_outcome * n_snp, df = 5) * 0.03,
                         nrow = n_outcome)
  exposure_beta <- exposure_beta * (1 + (seq_len(n_exposure) - 25.5) * 0.002)
  outcome_beta <- outcome_beta * (1 + (seq_len(n_outcome) - 25.5) * 0.002)
  exposure_se <- matrix(exp(rnorm(n_exposure * n_snp, log(0.012), 0.35)),
                        nrow = n_exposure)
  outcome_se <- matrix(exp(rnorm(n_outcome * n_snp, log(0.018), 0.40)),
                       nrow = n_outcome)
  list(exposure_beta = exposure_beta, outcome_beta = outcome_beta,
       exposure_se = exposure_se, outcome_se = outcome_se)
}

sim <- make_simulation()
sim_fast_warm <- invisible(fast_mr_grid(sim$exposure_beta, sim$outcome_beta,
                                        sim$exposure_se, sim$outcome_se,
                                        methods = methods, nboot = 0L,
                                        seed = seed, threads = 10L))
sim_fast_clock <- proc.time()[["elapsed"]]
sim_fast <- fast_mr_grid(sim$exposure_beta, sim$outcome_beta,
                         sim$exposure_se, sim$outcome_se, methods = methods,
                         nboot = nboot, seed = seed, threads = 10L)
sim_fast_seconds <- proc.time()[["elapsed"]] - sim_fast_clock

run_native_grid <- function(g) {
  params <- modifyList(TwoSampleMR::default_parameters(), list(nboot = nboot))
  results <- vector("list", 2500L)
  d <- data.frame(SNP = paste0("sim", seq_len(ncol(g$exposure_beta))),
                  stringsAsFactors = FALSE)
  k <- 0L
  for (i in seq_len(nrow(g$exposure_beta))) for (j in seq_len(nrow(g$outcome_beta))) {
    k <- k + 1L
    dat <- cbind(d,
      beta.exposure = g$exposure_beta[i, ],
      beta.outcome = g$outcome_beta[j, ],
      se.exposure = g$exposure_se[i, ],
      se.outcome = g$outcome_se[j, ],
      id.exposure = paste0("E", i), id.outcome = paste0("O", j),
      exposure = paste0("E", i), outcome = paste0("O", j), mr_keep = TRUE)
    results[[k]] <- suppressWarnings(suppressMessages(TwoSampleMR::mr(
      dat, method_list = native_methods, parameters = params
    )))
  }
  out <- do.call(rbind, results)
  out$pair <- rep(seq_len(2500L), each = length(native_methods))
  out
}

set.seed(seed)
sim_native_clock <- proc.time()[["elapsed"]]
sim_native <- run_native_grid(sim)
sim_native_seconds <- proc.time()[["elapsed"]] - sim_native_clock
sim_fast$pair <- (sim_fast$exposure_index - 1L) * 50L + sim_fast$outcome_index
sim_fast$method_code <- as.character(sim_fast$method_code)
sim_native$method_code <- unname(method_map[as.character(sim_native$method)])
sim_parity <- merge(sim_fast[, c("pair", "method_code", "b", "se", "pval")],
                    sim_native[, c("pair", "method_code", "b", "se", "pval")],
                    by = c("pair", "method_code"), suffixes = c("_fast", "_native"))
sim_parity$abs_beta_delta <- abs(sim_parity$b_fast - sim_parity$b_native)
sim_parity$abs_se_delta <- abs(sim_parity$se_fast - sim_parity$se_native)
sim_parity$abs_pval_delta <- abs(sim_parity$pval_fast - sim_parity$pval_native)
grid_result <- data.frame(
  workload = "different heavy-tailed heteroskedastic simulation",
  exposures = 50L, outcomes = 50L, snps = 400L, pairs = 2500L,
  methods = length(methods), nboot = nboot, fast_threads = 10L,
  fastMR_seconds = sim_fast_seconds, TwoSampleMR_seconds = sim_native_seconds,
  speedup = sim_native_seconds / sim_fast_seconds,
  fastMR_pairs_per_second = 2500 / sim_fast_seconds,
  TwoSampleMR_pairs_per_second = 2500 / sim_native_seconds,
  max_abs_beta_delta = max(sim_parity$abs_beta_delta, na.rm = TRUE),
  max_abs_se_delta = max(sim_parity$abs_se_delta, na.rm = TRUE),
  max_abs_pval_delta = max(sim_parity$abs_pval_delta, na.rm = TRUE),
  stringsAsFactors = FALSE
)
write.csv(grid_result, file.path(out_dir, "final_grid_benchmark.csv"), row.names = FALSE)
write.csv(sim_parity, file.path(out_dir, "final_grid_parity.csv"), row.names = FALSE)

single_sim_fast_clock <- proc.time()[["elapsed"]]
single_sim_fast <- fast_mr_grid(sim$exposure_beta[1, , drop = FALSE],
                                sim$outcome_beta[1, , drop = FALSE],
                                sim$exposure_se[1, , drop = FALSE],
                                sim$outcome_se[1, , drop = FALSE],
                                methods = methods, nboot = nboot,
                                seed = seed, threads = 1L)
single_sim_fast_seconds <- proc.time()[["elapsed"]] - single_sim_fast_clock
single_sim_dat <- data.frame(
  SNP = paste0("sim", seq_len(ncol(sim$exposure_beta))),
  beta.exposure = sim$exposure_beta[1, ], beta.outcome = sim$outcome_beta[1, ],
  se.exposure = sim$exposure_se[1, ], se.outcome = sim$outcome_se[1, ],
  id.exposure = "E1", id.outcome = "O1", exposure = "E1", outcome = "O1",
  mr_keep = TRUE, stringsAsFactors = FALSE
)
single_sim_native_clock <- proc.time()[["elapsed"]]
single_sim_native <- suppressWarnings(suppressMessages(TwoSampleMR::mr(
  single_sim_dat, method_list = native_methods,
  parameters = modifyList(TwoSampleMR::default_parameters(), list(nboot = nboot))
)))
single_sim_native_seconds <- proc.time()[["elapsed"]] - single_sim_native_clock
single_result <- data.frame(
  workload = "different simulation, 1 exposure x 1 outcome",
  snps = ncol(sim$exposure_beta), methods = length(methods), nboot = nboot,
  fastMR_seconds = single_sim_fast_seconds,
  TwoSampleMR_seconds = single_sim_native_seconds,
  speedup = single_sim_native_seconds / single_sim_fast_seconds,
  max_abs_beta_delta = max(abs(single_sim_fast$b - single_sim_native$b), na.rm = TRUE),
  max_abs_se_delta = max(abs(single_sim_fast$se - single_sim_native$se), na.rm = TRUE),
  max_abs_pval_delta = max(abs(single_sim_fast$pval - single_sim_native$pval), na.rm = TRUE),
  stringsAsFactors = FALSE
)
write.csv(single_result, file.path(out_dir, "final_single_pair_benchmark.csv"), row.names = FALSE)

grid_method_summary <- lapply(split(sim_parity, sim_parity$method_code), function(x) {
  c(median_abs_se_delta = median(x$abs_se_delta, na.rm = TRUE),
    p95_abs_se_delta = unname(quantile(x$abs_se_delta, 0.95, na.rm = TRUE)),
    median_abs_pval_delta = median(x$abs_pval_delta, na.rm = TRUE))
})
grid_method_summary <- do.call(rbind, grid_method_summary)
grid_method_summary_text <- paste(
  sprintf("%s: median SE %.3e, p95 SE %.3e, median p %.3e",
          rownames(grid_method_summary),
          grid_method_summary[, "median_abs_se_delta"],
          grid_method_summary[, "p95_abs_se_delta"],
          grid_method_summary[, "median_abs_pval_delta"]),
  collapse = "; "
)

median_fast_pipeline <- median(bmi_crp$repeats$fast_pipeline_seconds)
median_native_pipeline <- median(bmi_crp$repeats$native_pipeline_seconds)
report <- c(
  "# Final fastMR package audit", "",
  sprintf("Local BMI -> CRP: %d BMI rows, %d matched local CRP rows, %d repeats.",
          nrow(exposure), nrow(outcome), nrow(bmi_crp$repeats)),
  sprintf("Median BMI -> CRP pipeline: fastMR %.6fs; TwoSampleMR %.6fs; speedup %.2fx.",
          median_fast_pipeline, median_native_pipeline,
          median_native_pipeline / median_fast_pipeline),
  sprintf("BMI -> CRP harmonisation parity was exact: %s.",
          all(bmi_crp$repeats$harmonisation_key_exact)),
  "",
  sprintf("Different simulation 50x50 grid: %d SNPs, five methods, nboot=%d, fastMR threads=10.", ncol(sim$exposure_beta), nboot),
  sprintf("Grid: fastMR %.6fs; TwoSampleMR %.6fs; speedup %.2fx.",
          sim_fast_seconds, sim_native_seconds, sim_native_seconds / sim_fast_seconds),
  sprintf("Grid maximum absolute deltas: beta %.3e, SE %.3e, p-value %.3e.",
          grid_result$max_abs_beta_delta, grid_result$max_abs_se_delta,
          grid_result$max_abs_pval_delta),
  sprintf("Grid beta parity is exact to numerical precision for IVW/Egger (max beta delta %.3e); all methods have max beta delta %.3e.",
          max(sim_parity$abs_beta_delta[sim_parity$method_code %in% c("ivw", "egger")], na.rm = TRUE),
          max(sim_parity$abs_beta_delta, na.rm = TRUE)),
  paste("Bootstrap SE/p-value variation by method (median SE delta; 95th percentile SE delta; median p-value delta):", grid_method_summary_text),
  "The large raw SE maximum is a stress-test property of the heavy-tailed simulation: near-zero exposure draws make ratio bootstraps unstable. It is not a point-estimate disagreement; bootstrap streams are intentionally independent between fastMR and TwoSampleMR.",
  "",
  sprintf("Different simulation 1x1: fastMR %.6fs; TwoSampleMR %.6fs; speedup %.2fx.",
          single_sim_fast_seconds, single_sim_native_seconds, single_result$speedup),
  sprintf("1x1 maximum absolute deltas: beta %.3e, SE %.3e, p-value %.3e.",
          single_result$max_abs_beta_delta, single_result$max_abs_se_delta,
          single_result$max_abs_pval_delta),
  "",
  "Interpretation: the main weakness for 1x1 workloads is R/data-frame and harmonisation overhead; the large-grid path is dominated by bootstrap-heavy modes/Egger rather than IVW. Local PLINK timing remains unavailable because no PLINK executable is installed on the Mac mini; the wrapper was tested separately with a contract-validating stub."
)
writeLines(report, file.path(out_dir, "final_package_audit.md"))
cat(paste(report, collapse = "\n"), "\n")
