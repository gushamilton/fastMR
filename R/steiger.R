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

fastmr_steiger_rtest_p_vector <- function(r_exp, r_out, n_exp, n_out) {
  n <- max(length(r_exp), length(r_out), length(n_exp), length(n_out))
  recycle <- function(x) if (length(x) == 1L) rep(x, n) else x
  r_exp <- recycle(r_exp)
  r_out <- recycle(r_out)
  n_exp <- recycle(n_exp)
  n_out <- recycle(n_out)
  p <- rep(NA_real_, n)
  valid <- is.finite(r_exp) & is.finite(r_out) &
    is.finite(n_exp) & is.finite(n_out) & n_exp > 3 & n_out > 3
  if (any(valid)) {
    z <- suppressWarnings((0.5 * log((1 + r_exp[valid]) / (1 - r_exp[valid])) -
      0.5 * log((1 + r_out[valid]) / (1 - r_out[valid]))) /
      sqrt(1 / (n_exp[valid] - 3) + 1 / (n_out[valid] - 3)))
    p[valid] <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  }
  p
}

fastmr_steiger_r2_from_bsen <- function(beta, se, n) {
  f <- (beta / se)^2
  f / (n - 2 + f)
}

fastmr_steiger_effective_n <- function(ncase, ncontrol) {
  2 / (1 / ncase + 1 / ncontrol)
}

fastmr_steiger_population_af <- function(af, prop, odds_ratio, prevalence) {
  eps <- 1e-15
  a <- odds_ratio - 1
  b <- (af + prop) * (1 - odds_ratio) - 1
  c_value <- odds_ratio * af * prop
  z <- numeric(length(odds_ratio))
  linear <- abs(a) < eps
  z[linear] <- -c_value[linear] / b[linear]
  quadratic <- !linear
  if (any(quadratic)) {
    discriminant <- pmax(0, b[quadratic]^2 -
      4 * a[quadratic] * c_value[quadratic])
    sqrt_discriminant <- sqrt(discriminant)
    two_a <- 2 * a[quadratic]
    z_pos <- (-b[quadratic] + sqrt_discriminant) / two_a
    z_neg <- (-b[quadratic] - sqrt_discriminant) / two_a
    af_q <- af[quadratic]
    prop_q <- prop[quadratic]
    tolerance <- -1e-7
    valid_pos <- z_pos >= tolerance & (prop_q - z_pos) >= tolerance &
      (af_q - z_pos) >= tolerance &
      (1 + z_pos - af_q - prop_q) >= tolerance
    z[quadratic] <- ifelse(valid_pos, z_pos, z_neg)
  }
  af_controls <- (af - z) / (1 - prop)
  af_cases <- z / prop
  af_controls * (1 - prevalence) + af_cases * prevalence
}

fastmr_steiger_r2_from_lor <- function(beta, eaf, ncase, ncontrol, prevalence) {
  proportion <- ncase / (ncase + ncontrol)
  population_af <- fastmr_steiger_population_af(
    eaf, proportion, exp(beta), prevalence)
  genetic_variance <- beta^2 * population_af * (1 - population_af)
  residual_variance <- pi^2 / 3
  genetic_variance / (genetic_variance + residual_variance)
}

fastmr_steiger_unique <- function(x) {
  length(unique(x)) == 1L
}

fastmr_steiger_add_rsq_one <- function(data, what) {
  units_name <- paste0("units.", what)
  rsq_name <- paste0("rsq.", what)
  effective_n_name <- paste0("effective_n.", what)
  if (!units_name %in% names(data)) data[[units_name]] <- NA_character_
  if (rsq_name %in% names(data)) return(data)

  p_name <- paste0("pval.", what)
  beta_name <- paste0("beta.", what)
  se_name <- paste0("se.", what)
  eaf_name <- paste0("eaf.", what)
  sample_name <- paste0("samplesize.", what)
  data[[rsq_name]] <- rep(NA_real_, nrow(data))
  if (p_name %in% names(data)) {
    p <- suppressWarnings(as.numeric(data[[p_name]]))
    p[!is.na(p) & p < 9.99999999999999e-301] <-
      9.99999999999999e-301
    data[[p_name]] <- p
  } else {
    p <- rep(NA_real_, nrow(data))
  }
  beta <- if (beta_name %in% names(data))
    suppressWarnings(as.numeric(data[[beta_name]])) else rep(NA_real_, nrow(data))
  se <- if (se_name %in% names(data))
    suppressWarnings(as.numeric(data[[se_name]])) else rep(NA_real_, nrow(data))
  eaf <- if (eaf_name %in% names(data))
    suppressWarnings(as.numeric(data[[eaf_name]])) else rep(NA_real_, nrow(data))
  samplesize <- if (sample_name %in% names(data))
    suppressWarnings(as.numeric(data[[sample_name]])) else rep(NA_real_, nrow(data))
  units <- as.character(data[[units_name]])

  if (length(units) && !is.na(units[[1L]]) && units[[1L]] == "log odds") {
    prevalence_name <- paste0("prevalence.", what)
    if (!prevalence_name %in% names(data)) {
      data[[prevalence_name]] <- rep(0.1, nrow(data))
      warning(paste0("Assuming ", what,
                     " prevalence of 0.1. Alternatively, add prevalence.", what,
                     " column and re-run."), call. = FALSE)
    }
    prevalence <- suppressWarnings(as.numeric(data[[prevalence_name]]))
    ncase_name <- paste0("ncase.", what)
    ncontrol_name <- paste0("ncontrol.", what)
    ncase <- if (ncase_name %in% names(data))
      suppressWarnings(as.numeric(data[[ncase_name]])) else rep(NA_real_, nrow(data))
    ncontrol <- if (ncontrol_name %in% names(data))
      suppressWarnings(as.numeric(data[[ncontrol_name]])) else rep(NA_real_, nrow(data))
    valid <- is.finite(beta) & is.finite(eaf) & is.finite(ncase) &
      is.finite(ncontrol) & is.finite(prevalence)
    if (any(valid)) {
      data[[rsq_name]][valid] <- fastmr_steiger_r2_from_lor(
        beta[valid], eaf[valid], ncase[valid], ncontrol[valid], prevalence[valid])
      effective_n <- if (effective_n_name %in% names(data))
        suppressWarnings(as.numeric(data[[effective_n_name]])) else rep(NA_real_, nrow(data))
      effective_n[valid] <- fastmr_steiger_effective_n(ncase[valid], ncontrol[valid])
      data[[effective_n_name]] <- effective_n
    }
    return(data)
  }

  is_sd <- length(units) && all(!is.na(units) & grepl("SD", units)) &&
    all(!is.na(eaf))
  if (is_sd) {
    data[[rsq_name]] <- 2 * beta^2 * eaf * (1 - eaf)
    data[[effective_n_name]] <- samplesize
    return(data)
  }

  valid <- !is.na(p) & !is.na(samplesize)
  if (any(valid)) {
    data[[rsq_name]][valid] <- fastmr_steiger_r2_from_bsen(
      beta[valid], se[valid], samplesize[valid])
    data[[effective_n_name]] <- samplesize
  }
  data
}

#' Add per-SNP Steiger directionality flags and p-values
#'
#' This is the dependency-light local equivalent of
#' [TwoSampleMR::steiger_filtering()]. It supports supplied `rsq.*` columns,
#' standard-error/sample-size approximation, SD-scaled quantitative traits,
#' and log-odds traits with allele frequencies and case/control counts.
#'
#' @param data A TwoSampleMR-style harmonised data frame.
#' @return The input rows with `rsq.exposure`, `rsq.outcome`, effective sample
#'   sizes, `steiger_dir`, and `steiger_pval` added.
#' @export
fast_mr_steiger_filtering <- function(data) {
  if (!is.data.frame(data)) stop("data must be a data.frame", call. = FALSE)
  groups <- fastmr_diagnostic_groups(data)
  rows <- lapply(groups, function(group) {
    x <- group$data
    if (!"units.exposure" %in% names(x)) x$units.exposure <- NA_character_
    if (!"units.outcome" %in% names(x)) x$units.outcome <- NA_character_
    if (!fastmr_steiger_unique(x$exposure) || !fastmr_steiger_unique(x$outcome) ||
        !fastmr_steiger_unique(x$units.exposure) ||
        !fastmr_steiger_unique(x$units.outcome)) {
      stop("each exposure/outcome pair must have unique labels and units",
           call. = FALSE)
    }
    x <- fastmr_steiger_add_rsq_one(x, "exposure")
    x <- fastmr_steiger_add_rsq_one(x, "outcome")
    if (!"effective_n.exposure" %in% names(x)) {
      x$effective_n.exposure <- NA_real_
    }
    if (!"effective_n.outcome" %in% names(x)) {
      x$effective_n.outcome <- NA_real_
    }
    x$steiger_dir <- x$rsq.exposure > x$rsq.outcome
    x$steiger_pval <- fastmr_steiger_rtest_p_vector(
      sqrt(x$rsq.exposure), sqrt(x$rsq.outcome),
      x$effective_n.exposure, x$effective_n.outcome)
    x
  })
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
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
