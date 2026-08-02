args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]]) else getwd()
.libPaths(c(file.path(root, ".local", "Rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(fastMR)
  library(TwoSampleMR)
})

n <- 800L
p <- 4L
s <- seq_len(n)
X <- cbind(
  exposure_1 = sin(s / 17) + s / n / 4,
  exposure_2 = cos(s / 23) - s / n / 5,
  exposure_3 = sin(s / 31 + 0.4),
  exposure_4 = cos(s / 13 - 0.2) + s / n / 6
)
y <- 0.35 * X[, 1] - 0.2 * X[, 2] + 0.15 * X[, 3] +
  0.1 * X[, 4] + sin(s / 7) / 40
sy <- 0.04 + (s %% 11) / 2000
P <- vapply(seq_len(p), function(j) {
  ifelse(s %% (j + 2L) == 0L, 1e-9, 1e-4)
}, numeric(n))
colnames(P) <- colnames(X)
mvdat <- list(
  exposure_beta = X,
  exposure_pval = P,
  exposure_se = matrix(0.02, nrow = n, ncol = p,
                       dimnames = dimnames(X)),
  outcome_beta = y,
  outcome_pval = rep(1e-5, n),
  outcome_se = sy,
  expname = data.frame(id.exposure = colnames(X),
                       exposure = colnames(X), stringsAsFactors = FALSE),
  outname = data.frame(id.outcome = "outcome_1", outcome = "Outcome",
                       stringsAsFactors = FALSE)
)

timed <- function(fun, repeats) {
  invisible(fun())
  elapsed <- numeric(repeats)
  result <- NULL
  for (i in seq_len(repeats)) {
    started <- proc.time()[["elapsed"]]
    result <- fun()
    elapsed[i] <- proc.time()[["elapsed"]] - started
  }
  list(result = result, seconds = median(elapsed))
}

compare <- function(component, fast_fun, native_fun, key) {
  fast <- timed(fast_fun, 20L)
  native <- timed(native_fun, 5L)
  a <- fast$result
  b <- native$result$result
  a <- a[match(key, a$id.exposure), , drop = FALSE]
  b <- b[match(key, b$id.exposure), , drop = FALSE]
  finite_max <- function(x) if (any(is.finite(x))) max(abs(x), na.rm = TRUE) else NA_real_
  data.frame(
    component = component,
    fastMR_seconds = fast$seconds,
    TwoSampleMR_seconds = native$seconds,
    speedup = native$seconds / fast$seconds,
    fast_rows = nrow(a),
    native_rows = nrow(b),
    row_count_match = nrow(a) == nrow(b),
    key_match = identical(as.character(a$id.exposure), as.character(b$id.exposure)),
    nsnp_match = identical(as.numeric(a$nsnp), as.numeric(b$nsnp)),
    max_abs_beta_delta = finite_max(a$b - b$b),
    max_abs_se_delta = finite_max(a$se - b$se),
    max_abs_pval_delta = finite_max(a$pval - b$pval),
    stringsAsFactors = FALSE
  )
}

key <- colnames(X)
native_mv_ivw <- function() suppressMessages(suppressWarnings(
  TwoSampleMR::mv_ivw(mvdat, pval_threshold = 5e-8)))
native_shared <- function() suppressMessages(suppressWarnings(
  TwoSampleMR::mv_multiple(mvdat, pval_threshold = 5e-8,
                           instrument_specific = FALSE, plots = FALSE)))
native_specific <- function() suppressMessages(suppressWarnings(
  TwoSampleMR::mv_multiple(mvdat, pval_threshold = 5e-8,
                           instrument_specific = TRUE, plots = FALSE)))
native_intercept <- function() suppressMessages(suppressWarnings(
  TwoSampleMR::mv_multiple(mvdat, pval_threshold = 5e-8,
                           intercept = TRUE, instrument_specific = FALSE,
                           plots = FALSE)))

rows <- list(
  compare("mv_ivw", function() fast_mr_multivariable_ivw(
    X, y, sy, P, pval_threshold = 5e-8), native_mv_ivw, key),
  compare("mv_multiple_shared", function() fast_mr_multivariable(
    X, y, sy, P, pval_threshold = 5e-8), native_shared, key),
  compare("mv_multiple_specific", function() fast_mr_multivariable(
    X, y, sy, P, pval_threshold = 5e-8, instrument_specific = TRUE),
    native_specific, key),
  compare("mv_multiple_intercept", function() fast_mr_multivariable(
    X, y, sy, P, pval_threshold = 5e-8, intercept = TRUE),
    native_intercept, key)
)
result <- do.call(rbind, rows)
dir.create(file.path(root, "outputs"), showWarnings = FALSE, recursive = TRUE)
write.csv(result, file.path(root, "outputs", "multivariable_native_parity.csv"),
          row.names = FALSE)
lines <- c(
  "# fastMR multivariable MR versus native TwoSampleMR", "",
  sprintf("TwoSampleMR 0.7.9; deterministic %d-SNP × %d-exposure simulation; timings are medians of 20 fast and 5 native calls.", n, p), "",
  "| component | fastMR s | TwoSampleMR s | speedup | rows | key/nsnp match | max beta delta | max SE delta | max p delta |",
  "|---|---:|---:|---:|---:|---|---:|---:|---:|"
)
for (i in seq_len(nrow(result))) {
  x <- result[i, ]
  lines <- c(lines, sprintf(
    "| %s | %.6f | %.6f | %.2fx | %d/%d | %s/%s | %.3e | %.3e | %.3e |",
    x$component, x$fastMR_seconds, x$TwoSampleMR_seconds, x$speedup,
    x$fast_rows, x$native_rows, x$key_match, x$nsnp_match,
    x$max_abs_beta_delta, x$max_abs_se_delta, x$max_abs_pval_delta))
}
writeLines(lines, file.path(root, "outputs", "multivariable_native_parity.md"))
print(result, row.names = FALSE, digits = 6)
