args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]]) else getwd()
.libPaths(c(file.path(root, ".local", "Rlib"),
            "/Users/fergushamilton/projects/twosamplemr-fast/.local/Rlib",
            .libPaths()))

suppressPackageStartupMessages({
  library(fastMR)
  library(TwoSampleMR)
})

# This is deliberately independent of the IL6/CRP fixture. The outcome is a
# causal signal plus pleiotropic noise, with heterogeneous standard errors and
# a mix of positive and negative instrument effects.
seed <- 20260824L
nboot <- 100L
threads <- 5L
fast_repeats <- 5L
native_repeats <- 3L

method_specs <- data.frame(
  method = c(
    "ivw", "ivw_fe", "ivw_mre", "egger", "egger_bootstrap", "uwr",
    "sign", "simple_median", "weighted_median", "penalised_weighted_median",
    "simple_mode", "weighted_mode", "wald_ratio"
  ),
  native = c(
    "mr_ivw", "mr_ivw_fe", "mr_ivw_mre", "mr_egger_regression",
    "mr_egger_regression_bootstrap", "mr_uwr", "mr_sign", "mr_simple_median",
    "mr_weighted_median", "mr_penalised_weighted_median", "mr_simple_mode",
    "mr_weighted_mode", "mr_wald_ratio"
  ),
  stringsAsFactors = FALSE
)

make_simulation <- function(n_snp) {
  set.seed(seed + n_snp)
  beta_exposure <- rnorm(n_snp, mean = 0, sd = 0.045)
  beta_exposure[abs(beta_exposure) < 0.008] <- 0.008
  se_exposure <- runif(n_snp, 0.008, 0.025)
  se_outcome <- runif(n_snp, 0.010, 0.030)
  causal_effect <- 0.28
  pleiotropy <- rt(n_snp, df = 5) * 0.008
  beta_outcome <- causal_effect * beta_exposure + pleiotropy +
    rnorm(n_snp, 0, se_outcome)
  list(
    SNP = paste0("sim", seq_len(n_snp)),
    beta.exposure = beta_exposure,
    beta.outcome = beta_outcome,
    se.exposure = se_exposure,
    se.outcome = se_outcome
  )
}

as_grid <- function(d) {
  list(
    exposure_beta = matrix(d$beta.exposure, nrow = 1L),
    outcome_beta = matrix(d$beta.outcome, nrow = 1L),
    exposure_se = matrix(d$se.exposure, nrow = 1L),
    outcome_se = matrix(d$se.outcome, nrow = 1L)
  )
}

as_tidy <- function(d) {
  data.frame(
    SNP = d$SNP,
    beta.exposure = d$beta.exposure,
    beta.outcome = d$beta.outcome,
    se.exposure = d$se.exposure,
    se.outcome = d$se.outcome,
    id.exposure = "sim_exposure",
    id.outcome = "sim_outcome",
    exposure = "sim_exposure",
    outcome = "sim_outcome",
    mr_keep = TRUE,
    stringsAsFactors = FALSE
  )
}

timed_fast <- function(g, method) {
  invisible(fast_mr_grid(g$exposure_beta, g$outcome_beta,
                         g$exposure_se, g$outcome_se, methods = method,
                         nboot = nboot, seed = seed, threads = threads))
  elapsed <- numeric(fast_repeats)
  result <- NULL
  for (i in seq_len(fast_repeats)) {
    started <- proc.time()[["elapsed"]]
    result <- fast_mr_grid(g$exposure_beta, g$outcome_beta,
                           g$exposure_se, g$outcome_se, methods = method,
                           nboot = nboot, seed = seed, threads = threads)
    elapsed[i] <- proc.time()[["elapsed"]] - started
  }
  list(result = result, seconds = median(elapsed))
}

timed_native <- function(dat, method) {
  params <- modifyList(TwoSampleMR::default_parameters(), list(nboot = nboot))
  elapsed <- numeric(native_repeats)
  result <- NULL
  for (i in seq_len(native_repeats)) {
    set.seed(seed)
    started <- proc.time()[["elapsed"]]
    result <- suppressMessages(suppressWarnings(
      TwoSampleMR::mr(dat, method_list = method, parameters = params)))
    elapsed[i] <- proc.time()[["elapsed"]] - started
  }
  list(result = result, seconds = median(elapsed))
}

one_result <- function(d, spec, scenario) {
  g <- as_grid(d)
  dat <- as_tidy(d)
  fast <- timed_fast(g, spec$method)
  native <- timed_native(dat, spec$native)
  if (nrow(fast$result) != 1L || nrow(native$result) != 1L) {
    stop(sprintf("Expected one result for %s; got fastMR=%d native=%d",
                spec$method, nrow(fast$result), nrow(native$result)))
  }
  data.frame(
    scenario = scenario,
    method_code = spec$method,
    n_snps = length(d$SNP),
    nboot = nboot,
    threads = threads,
    fastMR_seconds = fast$seconds,
    TwoSampleMR_seconds = native$seconds,
    speedup = native$seconds / fast$seconds,
    fastMR_beta = fast$result$b,
    TwoSampleMR_beta = native$result$b,
    abs_beta_difference = abs(fast$result$b - native$result$b),
    abs_se_difference = abs(fast$result$se - native$result$se),
    abs_pval_difference = abs(fast$result$pval - native$result$pval),
    stringsAsFactors = FALSE
  )
}

rows <- vector("list", nrow(method_specs))
multi <- make_simulation(400L)
for (i in seq_len(nrow(method_specs) - 1L)) {
  cat(sprintf("[%02d/%02d] %s / 400 SNP simulation\n", i,
              nrow(method_specs), method_specs$method[i]))
  rows[[i]] <- one_result(multi, method_specs[i, ], "causal_400_snp")
  print(rows[[i]], row.names = FALSE, digits = 6)
}

cat(sprintf("[%02d/%02d] wald_ratio / single-SNP simulation\n",
            nrow(method_specs), nrow(method_specs)))
rows[[nrow(method_specs)]] <- one_result(make_simulation(1L),
                                         method_specs[nrow(method_specs), ],
                                         "causal_1_snp")
print(rows[[nrow(method_specs)]], row.names = FALSE, digits = 6)

result <- do.call(rbind, rows)
dir.create(file.path(root, "outputs"), showWarnings = FALSE, recursive = TRUE)
write.csv(result, file.path(root, "outputs", "simulated_method_benchmark.csv"),
          row.names = FALSE)

lines <- c(
  "# Independent simulated method benchmark", "",
  sprintf("Causal simulation seed=%d; nboot=%d; fastMR threads=%d; fastMR repeats=%d; native repeats=%d.",
          seed, nboot, threads, fast_repeats, native_repeats),
  "The 400-SNP simulation uses a causal effect, heavy-tailed pleiotropic noise,",
  "heterogeneous standard errors, and mixed-sign instrument effects. Wald ratio",
  "is tested separately on a one-SNP simulation, as required by the method.", "",
  "| scenario | method | SNPs | fastMR s | TwoSampleMR s | speedup | abs beta delta | abs SE delta | abs p-value delta |",
  "|---|---|---:|---:|---:|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(result))) {
  x <- result[i, ]
  lines <- c(lines, sprintf(
    "| %s | %s | %d | %.6f | %.6f | %.2fx | %.3e | %.3e | %.3e |",
    x$scenario, x$method_code, x$n_snps, x$fastMR_seconds,
    x$TwoSampleMR_seconds, x$speedup, x$abs_beta_difference,
    x$abs_se_difference, x$abs_pval_difference))
}
writeLines(lines, file.path(root, "outputs", "simulated_method_benchmark.md"))
