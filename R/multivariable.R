#' Basic multivariable IVW MR
#'
#' This is the low-dependency equivalent of the basic `mv_ivw` path in
#' TwoSampleMR. It fits a weighted regression through the origin of the
#' outcome effects on all exposure effects and reports the coefficient for
#' each exposure. If `exposure_pval` is supplied, each reported coefficient
#' uses the instruments passing that exposure's threshold.
#'
#' @param exposure_beta Numeric matrix, SNPs by exposures.
#' @param outcome_beta Numeric outcome-effect vector.
#' @param outcome_se Numeric outcome-standard-error vector.
#' @param exposure_pval Optional matrix of exposure p-values.
#' @param pval_threshold P-value threshold for per-exposure instrument sets.
#' @return A tidy data frame with one row per exposure.
#' @export
fast_mr_multivariable_ivw <- function(exposure_beta, outcome_beta, outcome_se,
                                      exposure_pval = NULL, pval_threshold = 5e-8) {
  X <- as.matrix(exposure_beta)
  storage.mode(X) <- "double"
  y <- fastmr_numeric(outcome_beta, "outcome_beta")
  sy <- fastmr_numeric(outcome_se, "outcome_se")
  if (!is.matrix(X) || nrow(X) != length(y) || length(y) != length(sy)) {
    stop("exposure_beta must be a SNP-by-exposure matrix matching outcome vectors", call. = FALSE)
  }
  if (length(pval_threshold) != 1L || !is.finite(pval_threshold) || pval_threshold <= 0) {
    stop("pval_threshold must be positive and finite", call. = FALSE)
  }
  P <- NULL
  if (!is.null(exposure_pval)) {
    P <- as.matrix(exposure_pval)
    storage.mode(P) <- "double"
    if (!all(dim(P) == dim(X))) stop("exposure_pval must match exposure_beta", call. = FALSE)
  }
  p <- ncol(X)
  result <- vector("list", p)
  for (j in seq_len(p)) {
    keep <- is.finite(y) & is.finite(sy) & sy > 0 & apply(X, 1L, function(row) all(is.finite(row)))
    if (!is.null(P)) keep <- keep & is.finite(P[, j]) & P[, j] < pval_threshold
    n <- sum(keep)
    beta <- rep(NA_real_, p)
    se <- rep(NA_real_, p)
    if (n > p) {
      Xj <- X[keep, , drop = FALSE]
      yj <- y[keep]
      w <- 1 / sy[keep]^2
      xtwx <- crossprod(Xj, Xj * w)
      inv <- tryCatch(solve(xtwx), error = function(e) NULL)
      if (!is.null(inv)) {
        beta <- drop(inv %*% crossprod(Xj, w * yj))
        residual <- yj - drop(Xj %*% beta)
        sigma <- sqrt(sum(w * residual^2) / (n - p))
        se <- sqrt(diag(inv)) * sigma
      }
    }
    stat <- beta[j] / se[j]
    result[[j]] <- data.frame(
      id.exposure = if (!is.null(colnames(X))) colnames(X)[j] else as.character(j),
      method = "Multivariable IVW",
      method_code = "mv_ivw",
      nsnp = n,
      b = beta[j],
      se = se[j],
      pval = 2 * stats::pnorm(abs(stat), lower.tail = FALSE),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, result)
}
