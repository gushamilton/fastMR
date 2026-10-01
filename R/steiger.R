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

fastmr_steiger_recycle <- function(values) {
  lengths <- vapply(values, length, integer(1))
  size <- if (length(lengths)) max(lengths) else 0L
  if (size == 0L) return(lapply(values, function(value) numeric()))
  invalid <- lengths != 0L & lengths != 1L & lengths != size
  if (any(invalid)) {
    stop("Steiger R-squared inputs must have equal lengths or length one",
         call. = FALSE)
  }
  lapply(values, function(value) {
    if (length(value) == 0L) return(rep(NA_real_, size))
    if (length(value) == 1L && size > 1L) return(rep(value, size))
    value
  })
}

fastmr_steiger_set_reason <- function(reason, rows, value) {
  rows <- rows & reason == "ok"
  if (length(value) == length(reason)) {
    reason[rows] <- value[rows]
  } else {
    reason[rows] <- value
  }
  reason
}

#' Calculate per-variant Steiger R-squared values
#'
#' This composable primitive calculates per-variant variance explained using
#' one explicit trait model. Inputs of length one are recycled to the longest
#' input; all other non-empty inputs must have that length. Rows with invalid
#' required inputs return `NA` rather than stopping the whole calculation, with
#' a stable code in `reason` explaining why. A finite calculated value outside
#' the R-squared range is retained for diagnosis with `reason = "invalid_rsq"`.
#'
#' `continuous_bsen` calculates `F / (n - 2 + F)`, where
#' `F = (beta / se)^2`. `standardized` calculates
#' `2 * beta^2 * eaf * (1 - eaf)`. `binary_lor` uses the same log-odds
#' approximation as [fast_mr_steiger_filtering()] and returns the harmonic
#' effective sample size `2 / (1 / ncase + 1 / ncontrol)`.
#'
#' Missing required values use `missing_<input>` reason codes. Non-finite or
#' out-of-domain values use `invalid_<input>` codes. Sample size is not needed
#' for a valid `standardized` R-squared value, so `effective_n` is `NA` when
#' `n` is omitted while `valid` remains `TRUE`.
#'
#' @param beta Variant effect estimates. For `binary_lor`, these must be
#'   log-odds effects.
#' @param se Standard errors. Required by `continuous_bsen` and ignored by the
#'   other models.
#' @param n Sample sizes. Required by `continuous_bsen`; optional for
#'   `standardized`, where it is returned as `effective_n` when supplied.
#' @param eaf Effect-allele frequencies. Required by `standardized` and
#'   `binary_lor`.
#' @param model The variance model: `continuous_bsen` uses beta, standard error,
#'   and sample size; `standardized` uses beta and effect-allele frequency; and
#'   `binary_lor` uses log-odds beta, effect-allele frequency, prevalence, and
#'   case/control counts.
#' @param prevalence Population outcome prevalence. Required by `binary_lor`;
#'   there is deliberately no default.
#' @param ncase,ncontrol Case and control counts. Required by `binary_lor`.
#' @return A data frame with `rsq`, `effective_n`, `valid`, and `reason` for
#'   each input row. `reason` is `"ok"` for valid R-squared estimates.
#' @export
fast_mr_steiger_r2 <- function(beta, se = NULL, n = NULL, eaf = NULL,
                               model = c("continuous_bsen", "standardized",
                                         "binary_lor"),
                               prevalence = NULL, ncase = NULL,
                               ncontrol = NULL) {
  model <- match.arg(model)
  raw_values <- list(
    beta = beta, se = se, n = n, eaf = eaf, prevalence = prevalence,
    ncase = ncase, ncontrol = ncontrol)
  used <- switch(
    model,
    continuous_bsen = c("beta", "se", "n"),
    standardized = c("beta", "n", "eaf"),
    binary_lor = c("beta", "eaf", "prevalence", "ncase", "ncontrol")
  )
  values <- lapply(names(raw_values), function(name) {
    if (!name %in% used) return(numeric())
    fastmr_numeric(raw_values[[name]], name)
  })
  names(values) <- names(raw_values)
  values <- fastmr_steiger_recycle(values)
  size <- length(values$beta)
  rsq <- rep(NA_real_, size)
  effective_n <- rep(NA_real_, size)
  reason <- rep("ok", size)

  reason <- fastmr_steiger_set_reason(
    reason, is.na(values$beta), "missing_beta")
  reason <- fastmr_steiger_set_reason(
    reason, !is.finite(values$beta), "invalid_beta")

  if (model == "continuous_bsen") {
    reason <- fastmr_steiger_set_reason(
      reason, is.na(values$se), "missing_se")
    reason <- fastmr_steiger_set_reason(
      reason, !is.finite(values$se) | values$se <= 0, "invalid_se")
    reason <- fastmr_steiger_set_reason(
      reason, is.na(values$n), "missing_n")
    reason <- fastmr_steiger_set_reason(
      reason, !is.finite(values$n) | values$n <= 2, "invalid_n")
    valid_inputs <- reason == "ok"
    if (any(valid_inputs)) {
      rsq[valid_inputs] <- fastmr_steiger_r2_from_bsen(
        values$beta[valid_inputs], values$se[valid_inputs],
        values$n[valid_inputs])
      effective_n[valid_inputs] <- values$n[valid_inputs]
    }
  } else if (model == "standardized") {
    reason <- fastmr_steiger_set_reason(
      reason, is.na(values$eaf), "missing_eaf")
    reason <- fastmr_steiger_set_reason(
      reason, !is.finite(values$eaf) | values$eaf <= 0 | values$eaf >= 1,
      "invalid_eaf")
    valid_inputs <- reason == "ok"
    if (any(valid_inputs)) {
      rsq[valid_inputs] <- 2 * values$beta[valid_inputs]^2 *
        values$eaf[valid_inputs] * (1 - values$eaf[valid_inputs])
    }
    valid_n <- !is.na(values$n) & is.finite(values$n) & values$n > 0
    effective_n[valid_n] <- values$n[valid_n]
  } else {
    reason <- fastmr_steiger_set_reason(
      reason, is.na(values$eaf), "missing_eaf")
    reason <- fastmr_steiger_set_reason(
      reason, !is.finite(values$eaf) | values$eaf <= 0 | values$eaf >= 1,
      "invalid_eaf")
    reason <- fastmr_steiger_set_reason(
      reason, is.na(values$prevalence), "missing_prevalence")
    reason <- fastmr_steiger_set_reason(
      reason, !is.finite(values$prevalence) | values$prevalence <= 0 |
        values$prevalence >= 1, "invalid_prevalence")
    reason <- fastmr_steiger_set_reason(
      reason, is.na(values$ncase), "missing_ncase")
    reason <- fastmr_steiger_set_reason(
      reason, !is.finite(values$ncase) | values$ncase <= 0, "invalid_ncase")
    reason <- fastmr_steiger_set_reason(
      reason, is.na(values$ncontrol), "missing_ncontrol")
    reason <- fastmr_steiger_set_reason(
      reason, !is.finite(values$ncontrol) | values$ncontrol <= 0,
      "invalid_ncontrol")
    valid_inputs <- reason == "ok"
    if (any(valid_inputs)) {
      rsq[valid_inputs] <- fastmr_steiger_r2_from_lor(
        values$beta[valid_inputs], values$eaf[valid_inputs],
        values$ncase[valid_inputs], values$ncontrol[valid_inputs],
        values$prevalence[valid_inputs])
      effective_n[valid_inputs] <- fastmr_steiger_effective_n(
        values$ncase[valid_inputs], values$ncontrol[valid_inputs])
    }
  }

  invalid_result <- reason == "ok" &
    (!is.finite(rsq) | rsq < 0 | rsq > 1)
  reason[invalid_result] <- "invalid_rsq"
  data.frame(
    rsq = rsq,
    effective_n = effective_n,
    valid = reason == "ok",
    reason = reason,
    stringsAsFactors = FALSE
  )
}

fastmr_steiger_unique <- function(x) {
  length(unique(x)) == 1L
}

fastmr_steiger_add_rsq_one <- function(data, what) {
  units_name <- paste0("units.", what)
  rsq_name <- paste0("rsq.", what)
  effective_n_name <- paste0("effective_n.", what)
  if (!units_name %in% names(data)) data[[units_name]] <- NA_character_
  valid_name <- paste0("rsq_valid.", what)
  reason_name <- paste0("rsq_reason.", what)
  if (rsq_name %in% names(data)) {
    rsq <- suppressWarnings(as.numeric(data[[rsq_name]]))
    valid <- is.finite(rsq) & rsq >= 0 & rsq <= 1
    reason <- ifelse(is.na(rsq), "missing_rsq",
                     ifelse(valid, "ok", "invalid_rsq"))
    data[[rsq_name]] <- rsq
    data[[valid_name]] <- valid
    data[[reason_name]] <- reason
    return(data)
  }

  p_name <- paste0("pval.", what)
  beta_name <- paste0("beta.", what)
  se_name <- paste0("se.", what)
  eaf_name <- paste0("eaf.", what)
  sample_name <- paste0("samplesize.", what)
  if (p_name %in% names(data)) {
    p <- suppressWarnings(as.numeric(data[[p_name]]))
    p[!is.na(p) & p < 9.99999999999999e-301] <-
      9.99999999999999e-301
    data[[p_name]] <- p
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
    prevalence <- if (prevalence_name %in% names(data))
      suppressWarnings(as.numeric(data[[prevalence_name]])) else NULL
    ncase_name <- paste0("ncase.", what)
    ncontrol_name <- paste0("ncontrol.", what)
    ncase <- if (ncase_name %in% names(data))
      suppressWarnings(as.numeric(data[[ncase_name]])) else rep(NA_real_, nrow(data))
    ncontrol <- if (ncontrol_name %in% names(data))
      suppressWarnings(as.numeric(data[[ncontrol_name]])) else rep(NA_real_, nrow(data))
    result <- fast_mr_steiger_r2(
      beta, eaf = eaf, model = "binary_lor", prevalence = prevalence,
      ncase = ncase, ncontrol = ncontrol)
    data[[rsq_name]] <- result$rsq
    data[[effective_n_name]] <- result$effective_n
    data[[valid_name]] <- result$valid
    data[[reason_name]] <- result$reason
    return(data)
  }

  is_sd <- length(units) && all(!is.na(units) & grepl("SD", units))
  if (is_sd) {
    result <- fast_mr_steiger_r2(
      beta, n = samplesize, eaf = eaf, model = "standardized")
    data[[rsq_name]] <- result$rsq
    data[[effective_n_name]] <- result$effective_n
    data[[valid_name]] <- result$valid
    data[[reason_name]] <- result$reason
    return(data)
  }

  result <- fast_mr_steiger_r2(
    beta, se = se, n = samplesize, model = "continuous_bsen")
  data[[rsq_name]] <- result$rsq
  data[[effective_n_name]] <- result$effective_n
  data[[valid_name]] <- result$valid
  data[[reason_name]] <- result$reason
  data
}

#' Add per-SNP Steiger directionality flags and p-values
#'
#' This is the dependency-light local equivalent of
#' [TwoSampleMR::steiger_filtering()]. It supports supplied `rsq.*` columns,
#' standard-error/sample-size approximation, SD-scaled quantitative traits,
#' and log-odds traits with allele frequencies and case/control counts.
#'
#' Binary traits require explicit `prevalence.*`, `ncase.*`, and `ncontrol.*`
#' columns; prevalence is never assumed. Quantitative beta/SE/sample-size
#' inputs do not require p-values.
#'
#' @param data A TwoSampleMR-style harmonised data frame.
#' @return The input rows with `rsq.exposure`, `rsq.outcome`, effective sample
#'   sizes, `steiger_dir`, and `steiger_pval` added. Per-trait `rsq_valid.*` and
#'   `rsq_reason.*` columns describe R-squared validity; `steiger_pval_valid`
#'   and `steiger_pval_reason` explicitly describe p-value availability.
#' @export
fast_mr_steiger_filtering <- function(data) {
  if (!is.data.frame(data)) stop("data must be a data.frame", call. = FALSE)
  fastmr_prepare_vectors(data)
  g <- fastmr_diagnostic_group_index(data)
  if (!length(g$rows)) return(data.frame())
  n <- nrow(data)
  gid <- rep.int(seq_along(g$rows), lengths(g$rows))
  x <- data[unlist(g$rows, use.names = FALSE), , drop = FALSE]
  if (!"units.exposure" %in% names(x)) x$units.exposure <- NA_character_
  if (!"units.outcome" %in% names(x)) x$units.outcome <- NA_character_
  first <- cumsum(c(1L, lengths(g$rows)))[seq_along(g$rows)]
  for (column in c("exposure", "outcome", "units.exposure", "units.outcome")) {
    value <- x[[column]]
    if (is.null(value)) {
      stop("each exposure/outcome pair must have unique labels and units",
           call. = FALSE)
    }
    code <- match(value, unique(value))
    if (any(code != code[first][gid])) {
      stop("each exposure/outcome pair must have unique labels and units",
           call. = FALSE)
    }
  }
  for (what in c("exposure", "outcome")) {
    units <- as.character(x[[paste0("units.", what)]])
    binary <- !is.na(units) & units == "log odds"
    standardized <- !binary & !is.na(units) & grepl("SD", units)
    classes <- list(binary, standardized, !binary & !standardized)
    original <- names(x)
    # Compute every model from the pre-update data so a column created for one
    # model is not mistaken for a supplied rsq column by the next.
    base <- x
    for (rows in classes) {
      if (!any(rows)) next
      part <- fastmr_steiger_add_rsq_one(base[rows, , drop = FALSE], what)
      touched <- names(part)[!names(part) %in% original |
        names(part) %in% paste0(c("rsq.", "pval.", "effective_n.", "rsq_valid.", "rsq_reason."), what)]
      for (column in touched) {
        value <- part[[column]]
        if (is.null(x[[column]]) || !identical(typeof(x[[column]]), typeof(value)) ||
            !is.null(attributes(x[[column]]))) {
          x[[column]] <- value[rep(NA_integer_, n)]
        }
        x[[column]][rows] <- value
      }
    }
  }
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
  p_reason <- rep("ok", n)
  p_reason <- fastmr_steiger_set_reason(
    p_reason, !x$rsq_valid.exposure,
    paste0("exposure_", x$rsq_reason.exposure))
  p_reason <- fastmr_steiger_set_reason(
    p_reason, !x$rsq_valid.outcome,
    paste0("outcome_", x$rsq_reason.outcome))
  p_reason <- fastmr_steiger_set_reason(
    p_reason, is.na(x$effective_n.exposure), "missing_exposure_n")
  p_reason <- fastmr_steiger_set_reason(
    p_reason, !is.finite(x$effective_n.exposure) |
      x$effective_n.exposure <= 3, "invalid_exposure_n")
  p_reason <- fastmr_steiger_set_reason(
    p_reason, is.na(x$effective_n.outcome), "missing_outcome_n")
  p_reason <- fastmr_steiger_set_reason(
    p_reason, !is.finite(x$effective_n.outcome) |
      x$effective_n.outcome <= 3, "invalid_outcome_n")
  p_reason <- fastmr_steiger_set_reason(
    p_reason, !is.finite(x$steiger_pval), "invalid_steiger_pval")
  x$steiger_pval_valid <- p_reason == "ok"
  x$steiger_pval_reason <- p_reason
  x
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
