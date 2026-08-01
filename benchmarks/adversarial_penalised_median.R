args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]]) else getwd()
.libPaths(c(file.path(root, ".local", "Rlib"), .libPaths()))
suppressPackageStartupMessages(library(fastMR))

same_numeric <- function(a, b, tolerance = 0) {
  both_na <- is.na(a) & is.na(b)
  both_finite <- is.finite(a) & is.finite(b)
  all(both_na | (both_finite & abs(a - b) <= tolerance))
}

make_case <- function(case, n_snp) {
  x <- matrix(rnorm(4L * n_snp, sd = 0.08), nrow = 4L)
  y <- matrix(rnorm(5L * n_snp, sd = 0.08), nrow = 5L)
  sx <- matrix(runif(4L * n_snp, 0.005, 0.04), nrow = 4L)
  sy <- matrix(runif(5L * n_snp, 0.005, 0.04), nrow = 5L)
  if (case %% 8L == 0L) x[1L, seq_len(min(4L, n_snp))] <- 0
  if (case %% 8L == 1L) x[2L, seq_len(min(4L, n_snp))] <- c(1e-12, -1e-12, 5e-12, -5e-12)[seq_len(min(4L, n_snp))]
  if (case %% 8L == 2L) x[3L, ] <- -abs(x[3L, ])
  if (case %% 8L == 3L) y <- x[1L, ] * 0.25 + matrix(0, nrow = 5L, ncol = n_snp, byrow = TRUE)
  if (case %% 8L == 4L) sy <- sy * 100
  if (case %% 8L == 5L) sx <- sx * 100
  if (case %% 8L == 6L) x[4L, ] <- 0
  if (case %% 8L == 7L) {
    x[1L, ] <- rep(c(-0.1, 0.1, 0.2), length.out = n_snp)
    y <- matrix(rep(x[1L, ] * 0.4, 5L), nrow = 5L, byrow = TRUE)
  }
  list(exposure_beta = x, outcome_beta = y,
       exposure_se = sx, outcome_se = sy)
}

rows <- list()
index <- 0L
for (repeat_index in 1:3) {
  set.seed(20260905L + repeat_index)
  for (case in seq_len(100L)) {
    n_snp <- if (case %% 10L == 0L) 3L else sample(4:20, 1L)
    g <- make_case(case, n_snp)
    serial <- fast_mr_grid(g$exposure_beta, g$outcome_beta,
                           g$exposure_se, g$outcome_se,
                           methods = "penalised_weighted_median",
                           nboot = 25L, seed = 20260910L + repeat_index,
                           threads = 1L)
    parallel <- fast_mr_grid(g$exposure_beta, g$outcome_beta,
                             g$exposure_se, g$outcome_se,
                             methods = "penalised_weighted_median",
                             nboot = 25L, seed = 20260910L + repeat_index,
                             threads = 5L)
    index <- index + 1L
    deltas <- c(abs(serial$b - parallel$b), abs(serial$se - parallel$se),
                abs(serial$pval - parallel$pval))
    finite_deltas <- deltas[is.finite(deltas)]
    rows[[index]] <- data.frame(
      repeat_index = repeat_index, case = case, snps = n_snp,
      max_delta = if (length(finite_deltas)) max(finite_deltas) else 0,
      exact = same_numeric(serial$b, parallel$b) &&
        same_numeric(serial$se, parallel$se) &&
        same_numeric(serial$pval, parallel$pval),
      finite_serial = all(is.finite(serial$b) | is.na(serial$b)) &&
        all(is.finite(serial$se) | is.na(serial$se)) &&
        all(is.finite(serial$pval) | is.na(serial$pval)),
      stringsAsFactors = FALSE
    )
  }
}
result <- do.call(rbind, rows)
stopifnot(all(result$exact), all(result$finite_serial), max(result$max_delta) == 0)
dir.create(file.path(root, "outputs"), showWarnings = FALSE, recursive = TRUE)
write.csv(result, file.path(root, "outputs", "adversarial_penalised_median.csv"), row.names = FALSE)
writeLines(c(
  "# Penalised weighted median adversarial thread gate", "",
  sprintf("%d randomized/edge panels, repeated %d times; nboot=25; serial versus five threads.",
          nrow(result), max(result$repeat_index)),
  sprintf("Maximum serial-versus-five-thread delta: %.3e; failures: %d.",
          max(result$max_delta), sum(!result$exact | !result$finite_serial))),
  file.path(root, "outputs", "adversarial_penalised_median.md"))
print(data.frame(panels = nrow(result), max_delta = max(result$max_delta),
                 failures = sum(!result$exact | !result$finite_serial)), row.names = FALSE)
