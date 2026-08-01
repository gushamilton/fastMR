args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]]) else getwd()
.libPaths(c(file.path(root, ".local", "Rlib"), .libPaths()))
suppressPackageStartupMessages(library(fastMR))

same_numeric <- function(a, b, tolerance = 0) {
  both_na <- is.na(a) & is.na(b)
  both_finite <- is.finite(a) & is.finite(b)
  all(both_na | (both_finite & abs(a - b) <= tolerance))
}

make_case <- function(case, n_exposure, n_outcome, n_snp) {
  x <- matrix(rnorm(n_exposure * n_snp, sd = 0.08), nrow = n_exposure)
  y <- matrix(rnorm(n_outcome * n_snp, sd = 0.08), nrow = n_outcome)
  sx <- matrix(runif(n_exposure * n_snp, 0.005, 0.04), nrow = n_exposure)
  sy <- matrix(runif(n_outcome * n_snp, 0.005, 0.04), nrow = n_outcome)
  if (case %% 7L == 0L) x[1L, seq_len(min(4L, n_snp))] <- 0
  if (case %% 7L == 1L) x[1L, seq_len(min(4L, n_snp))] <-
    c(1e-12, -1e-12, 5e-12, -5e-12)[seq_len(min(4L, n_snp))]
  if (case %% 7L == 2L) x[1L, ] <- -abs(x[1L, ])
  if (case %% 7L == 3L) y <- matrix(rep(x[1L, ] * 0.25, n_outcome),
                                      nrow = n_outcome, byrow = TRUE)
  if (case %% 7L == 4L) sy <- sy * 100
  if (case %% 7L == 5L) sx <- sx * 100
  if (case %% 7L == 6L) x[1L, ] <- rep(c(-0.1, 0.1, 0.2), length.out = n_snp)
  list(exposure_beta = x, outcome_beta = y,
       exposure_se = sx, outcome_se = sy)
}

methods <- c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode")
rows <- vector("list", 150L)
index <- 0L
for (repeat_index in 1:3) {
  set.seed(20261010L + repeat_index)
  for (case in seq_len(50L)) {
    n_exposure <- sample(2:6, 1L)
    n_outcome <- sample(2:6, 1L)
    n_snp <- if (case %% 10L == 0L) 3L else sample(4:25, 1L)
    g <- make_case(case, n_exposure, n_outcome, n_snp)
    serial <- fast_mr_grid(g$exposure_beta, g$outcome_beta,
                           g$exposure_se, g$outcome_se, methods = methods,
                           nboot = 7L, seed = 20261020L + repeat_index,
                           threads = 1L)
    parallel <- fast_mr_grid(g$exposure_beta, g$outcome_beta,
                             g$exposure_se, g$outcome_se, methods = methods,
                             nboot = 7L, seed = 20261020L + repeat_index,
                             threads = 5L)
    numeric_columns <- c("b", "se", "pval", "Q", "Q_df", "Q_pval", "sigma",
                         "intercept", "intercept_se", "intercept_pval")
    deltas <- vapply(numeric_columns, function(name) {
      a <- serial[[name]]; b <- parallel[[name]]
      ok <- is.finite(a) & is.finite(b)
      if (any(ok)) max(abs(a[ok] - b[ok])) else 0
    }, numeric(1))
    index <- index + 1L
    rows[[index]] <- data.frame(
      repeat_index = repeat_index, case = case,
      exposures = n_exposure, outcomes = n_outcome, snps = n_snp,
      max_delta = max(deltas), exact = all(deltas == 0),
      finite_serial = all(vapply(numeric_columns, function(name)
        all(is.finite(serial[[name]]) | is.na(serial[[name]])), logical(1))),
      stringsAsFactors = FALSE
    )
  }
}
result <- do.call(rbind, rows)
stopifnot(all(result$exact), all(result$finite_serial), max(result$max_delta) == 0)
dir.create(file.path(root, "outputs"), showWarnings = FALSE, recursive = TRUE)
write.csv(result, file.path(root, "outputs", "adversarial_mixed_grid.csv"), row.names = FALSE)
writeLines(c(
  "# Mixed-method grid adversarial thread gate", "",
  sprintf("%d randomized/edge panels, repeated %d times; nboot=7; serial versus five threads.",
          nrow(result), max(result$repeat_index)),
  sprintf("Maximum serial-versus-five-thread delta: %.3e; failures: %d.",
          max(result$max_delta), sum(!result$exact | !result$finite_serial))),
  file.path(root, "outputs", "adversarial_mixed_grid.md"))
print(data.frame(panels = nrow(result), max_delta = max(result$max_delta),
                 failures = sum(!result$exact | !result$finite_serial)), row.names = FALSE)
