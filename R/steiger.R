fastmr_r_from_pn <- function(p, n) {
  p <- fastmr_numeric(p, "p-value")
  n <- fastmr_numeric(n, "sample size")
  if (length(n) == 1L && length(p) > 1L) n <- rep(n, length(p))
  if (length(p) != length(n)) stop("p-values and sample sizes must have equal length", call. = FALSE)
  f <- suppressWarnings(stats::qf(p, 1, n - 1, lower.tail = FALSE))
  r2 <- f / (n - 2 + f)
  bad <- !is.finite(f)
  if (any(bad)) r2[bad] <- NA_real_
  sqrt(r2)
}

fastmr_steiger_rtest_p <- function(r_exp, r_out, n_exp, n_out) {
  if (!is.finite(r_exp) || !is.finite(r_out) ||
      !is.finite(n_exp) || !is.finite(n_out) || n_exp <= 3 || n_out <= 3) {
    return(NA_real_)
  }
  z <- (0.5 * log((1 + r_exp) / (1 - r_exp)) -
          0.5 * log((1 + r_out) / (1 - r_out))) /
    sqrt(1 / (n_exp - 3) + 1 / (n_out - 3))
  2 * stats::pnorm(abs(z), lower.tail = FALSE)
}

#' Calculate the Steiger directionality test from SNP correlations
#'
#' Missing correlations are approximated from p-values and sample sizes using
#' the same quantitative-trait approximation as native TwoSampleMR. The
#' implementation intentionally omits the native plotting object.
#'
#' @param p_exp Exposure p-values.
#' @param p_out Outcome p-values.
#' @param n_exp Exposure sample sizes.
#' @param n_out Outcome sample sizes.
#' @param r_exp Optional SNP-exposure correlations.
#' @param r_out Optional SNP-outcome correlations.
#' @param r_xxo Exposure reliability/correlation correction, between 0 and 1.
#' @param r_yyo Outcome reliability/correlation correction, between 0 and 1.
#' @return A list with native-compatible Steiger R-squared and direction fields.
#' @export
fast_mr_steiger <- function(p_exp, p_out, n_exp, n_out,
                            r_exp = NA_real_, r_out = NA_real_,
                            r_xxo = 1, r_yyo = 1) {
  p_exp <- fastmr_numeric(p_exp, "p_exp")
  p_out <- fastmr_numeric(p_out, "p_out")
  n_exp <- fastmr_numeric(n_exp, "n_exp")
  n_out <- fastmr_numeric(n_out, "n_out")
  r_exp <- fastmr_numeric(r_exp, "r_exp")
  r_out <- fastmr_numeric(r_out, "r_out")
  if (length(n_exp) == 1L && length(p_exp) > 1L) n_exp <- rep(n_exp, length(p_exp))
  if (length(n_out) == 1L && length(p_out) > 1L) n_out <- rep(n_out, length(p_out))
  n <- max(length(p_exp), length(p_out), length(n_exp), length(n_out),
           length(r_exp), length(r_out))
  recycle <- function(x) if (length(x) == 1L) rep(x, n) else x
  p_exp <- recycle(p_exp); p_out <- recycle(p_out)
  n_exp <- recycle(n_exp); n_out <- recycle(n_out)
  r_exp <- recycle(r_exp); r_out <- recycle(r_out)
  if (any(lengths(list(p_exp, p_out, n_exp, n_out, r_exp, r_out)) != n)) {
    stop("Steiger inputs must have equal lengths or length one", call. = FALSE)
  }
  r_exp <- abs(r_exp)
  r_out <- abs(r_out)
  missing_exp <- is.na(r_exp) & !is.na(p_exp) & !is.na(n_exp)
  missing_out <- is.na(r_out) & !is.na(p_out) & !is.na(n_out)
  if (any(missing_exp)) r_exp[missing_exp] <- fastmr_r_from_pn(p_exp[missing_exp], n_exp[missing_exp])
  if (any(missing_out)) r_out[missing_out] <- fastmr_r_from_pn(p_out[missing_out], n_out[missing_out])
  keep <- !is.na(r_exp) | !is.na(r_out)
  total_exp <- sqrt(sum(r_exp[keep]^2, na.rm = TRUE))
  total_out <- sqrt(sum(r_out[keep]^2, na.rm = TRUE))
  if (length(r_xxo) != 1L || !is.finite(r_xxo) || r_xxo < 0 || r_xxo > 1) {
    stop("r_xxo must be one finite value between 0 and 1", call. = FALSE)
  }
  if (length(r_yyo) != 1L || !is.finite(r_yyo) || r_yyo < 0 || r_yyo > 1) {
    stop("r_yyo must be one finite value between 0 and 1", call. = FALSE)
  }
  adjusted_exp <- sqrt(total_exp^2 / r_xxo^2)
  adjusted_out <- sqrt(total_out^2 / r_yyo^2)
  n_exp_mean <- mean(n_exp, na.rm = TRUE)
  n_out_mean <- mean(n_out, na.rm = TRUE)
  test <- fastmr_steiger_rtest_p(total_exp, total_out, n_exp_mean, n_out_mean)
  test_adjusted <- fastmr_steiger_rtest_p(adjusted_exp, adjusted_out,
                                           n_exp_mean, n_out_mean)
  a <- max(total_exp, total_out)
  b <- min(total_exp, total_out)
  vz <- a * log(a) - b * log(b) + a * b * (log(b) - log(a))
  vz0 <- -2 * b - b * log(a) - a * b * log(a) + 2 * a * b
  vz1 <- abs(vz - vz0)
  list(
    r2_exp = total_exp^2,
    r2_out = total_out^2,
    r2_exp_adj = adjusted_exp^2,
    r2_out_adj = adjusted_out^2,
    correct_causal_direction = total_exp > total_out,
    steiger_test = test,
    correct_causal_direction_adj = adjusted_exp > adjusted_out,
    steiger_test_adj = test_adjusted,
    vz = vz,
    vz0 = vz0,
    vz1 = vz1,
    sensitivity_ratio = vz1 / vz0,
    sensitivity_plot = NULL
  )
}

#' Run the tidy Steiger directionality test
#'
#' @param data A harmonised TwoSampleMR-style data frame. Supply either
#'   `r.exposure`/`r.outcome` or p-value and sample-size columns for both traits.
#' @return One tidy directionality row per exposure/outcome pair, or `NULL`
#'   when neither correlation nor p-value/sample-size inputs are available.
#' @export
fast_mr_directionality_test <- function(data) {
  if (!is.data.frame(data)) stop("data must be a data.frame", call. = FALSE)
  has_r <- all(c("r.exposure", "r.outcome") %in% names(data))
  has_pn <- all(c("pval.exposure", "pval.outcome",
                  "samplesize.exposure", "samplesize.outcome") %in% names(data))
  if (!has_r && !has_pn) {
    message("r.exposure and/or r.outcome not present.")
    message("Cannot calculate approximate SNP correlations without p-values and sample sizes.")
    return(NULL)
  }
  groups <- fastmr_diagnostic_groups(data)
  rows <- lapply(groups, function(group) {
    x <- group$data
    p_exp <- if ("pval.exposure" %in% names(x)) x$pval.exposure else rep(NA_real_, nrow(x))
    p_out <- if ("pval.outcome" %in% names(x)) x$pval.outcome else rep(NA_real_, nrow(x))
    n_exp <- if ("samplesize.exposure" %in% names(x)) x$samplesize.exposure else rep(NA_real_, nrow(x))
    n_out <- if ("samplesize.outcome" %in% names(x)) x$samplesize.outcome else rep(NA_real_, nrow(x))
    r_exp <- if ("r.exposure" %in% names(x)) x$r.exposure else rep(NA_real_, nrow(x))
    r_out <- if ("r.outcome" %in% names(x)) x$r.outcome else rep(NA_real_, nrow(x))
    result <- fast_mr_steiger(p_exp, p_out, n_exp, n_out, r_exp, r_out)
    data.frame(
      id.exposure = group$id.exposure,
      id.outcome = group$id.outcome,
      exposure = group$exposure,
      outcome = group$outcome,
      snp_r2.exposure = result$r2_exp,
      snp_r2.outcome = result$r2_out,
      correct_causal_direction = result$correct_causal_direction,
      steiger_pval = result$steiger_test,
      stringsAsFactors = FALSE
    )
  })
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}
