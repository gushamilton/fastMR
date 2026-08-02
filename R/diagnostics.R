fastmr_diagnostic_groups <- function(data) {
  if (!is.data.frame(data)) stop("data must be a data.frame", call. = FALSE)
  # Validate the same required columns and numeric conversions as fast_mr.
  fastmr_prepare_vectors(data)
  n <- nrow(data)
  id.exp <- if ("id.exposure" %in% names(data)) as.character(data$id.exposure) else rep("", n)
  id.out <- if ("id.outcome" %in% names(data)) as.character(data$id.outcome) else rep("", n)
  id.exp[is.na(id.exp)] <- ""
  id.out[is.na(id.out)] <- ""
  pairs <- unique(data.frame(id.exposure = id.exp, id.outcome = id.out,
                             stringsAsFactors = FALSE))
  groups <- vector("list", nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    index <- which(id.exp == pairs$id.exposure[[i]] & id.out == pairs$id.outcome[[i]])
    label.exp <- if ("exposure" %in% names(data)) as.character(data$exposure[index[[1L]]]) else pairs$id.exposure[[i]]
    label.out <- if ("outcome" %in% names(data)) as.character(data$outcome[index[[1L]]]) else pairs$id.outcome[[i]]
    if (is.na(label.exp)) label.exp <- pairs$id.exposure[[i]]
    if (is.na(label.out)) label.out <- pairs$id.outcome[[i]]
    groups[[i]] <- list(
      data = data[index, , drop = FALSE],
      id.exposure = pairs$id.exposure[[i]],
      id.outcome = pairs$id.outcome[[i]],
      exposure = label.exp,
      outcome = label.out
    )
  }
  groups
}

fastmr_diagnostic_keep <- function(data) {
  n <- nrow(data)
  keep <- if ("mr_keep" %in% names(data)) {
    !is.na(data$mr_keep) & as.logical(data$mr_keep)
  } else {
    rep(TRUE, n)
  }
  prepared <- fastmr_prepare_vectors(data)
  valid <- is.finite(prepared$beta.exposure) & is.finite(prepared$beta.outcome) &
    is.finite(prepared$se.exposure) & is.finite(prepared$se.outcome) &
    prepared$se.exposure > 0 & prepared$se.outcome > 0
  if (any(keep & !valid)) {
    stop("kept rows must have finite beta values and positive standard errors", call. = FALSE)
  }
  snp <- as.character(data$SNP)
  snp[is.na(snp)] <- ""
  if (any(keep & !nzchar(snp))) stop("kept rows must have non-empty SNP identifiers", call. = FALSE)
  keep
}

fastmr_diagnostic_result <- function(group, result, method) {
  registry <- fastmr_method_registry()
  data.frame(
    id.exposure = group$id.exposure,
    id.outcome = group$id.outcome,
    outcome = group$outcome,
    exposure = group$exposure,
    method = registry$method[match(method, registry$code)],
    Q = fastmr_scalar(result, "Q"),
    Q_df = fastmr_scalar(result, "Q_df"),
    Q_pval = fastmr_scalar(result, "Q_pval"),
    stringsAsFactors = FALSE
  )
}

#' Calculate TwoSampleMR-compatible heterogeneity statistics
#'
#' @param data A harmonised TwoSampleMR-style data frame.
#' @param methods Heterogeneity-capable fastMR methods. The default matches
#'   native `TwoSampleMR::mr_heterogeneity()` (`ivw` and `egger`).
#' @param threads Maximum native worker count passed to [fast_mr()].
#' @return A tidy data frame with `Q`, `Q_df`, and `Q_pval` per method and pair.
#' @export
fast_mr_heterogeneity <- function(data, methods = c("ivw", "egger"), threads = 1) {
  methods <- fastmr_normalize_methods(methods)
  supported <- c("ivw", "ivw_fe", "ivw_mre", "egger", "uwr")
  unsupported <- setdiff(methods, supported)
  if (length(unsupported)) {
    stop("heterogeneity is not defined for method(s): ",
         paste(unsupported, collapse = ", "), call. = FALSE)
  }
  groups <- fastmr_diagnostic_groups(data)
  rows <- vector("list", length(groups) * length(methods))
  k <- 0L
  for (group in groups) {
    result <- fast_mr(group$data, methods = methods, nboot = 0,
                      threads = threads)
    for (method in methods) {
      k <- k + 1L
      rows[[k]] <- fastmr_diagnostic_result(
        group, result[result$method_code == method, , drop = FALSE], method)
    }
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

#' Calculate the MR-Egger intercept pleiotropy test
#'
#' @param data A harmonised TwoSampleMR-style data frame.
#' @param threads Maximum native worker count passed to [fast_mr()].
#' @return A tidy data frame compatible with
#'   `TwoSampleMR::mr_pleiotropy_test()`.
#' @export
fast_mr_pleiotropy_test <- function(data, threads = 1) {
  groups <- fastmr_diagnostic_groups(data)
  rows <- lapply(groups, function(group) {
    result <- fast_mr(group$data, methods = "egger", nboot = 0,
                      threads = threads)
    data.frame(
      id.exposure = group$id.exposure,
      id.outcome = group$id.outcome,
      outcome = group$outcome,
      exposure = group$exposure,
      egger_intercept = fastmr_scalar(result, "intercept"),
      se = fastmr_scalar(result, "intercept_se"),
      pval = fastmr_scalar(result, "intercept_pval"),
      stringsAsFactors = FALSE
    )
  })
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

fastmr_diagnostic_sample_size <- function(data) {
  candidates <- c("samplesize.outcome", "samplesize", "sample_size")
  present <- candidates[candidates %in% names(data)]
  if (!length(present)) return(NA_real_)
  candidate <- present[[1L]]
  value <- suppressWarnings(as.numeric(data[[candidate]][[1L]]))
  if (length(value) && is.finite(value)) value else NA_real_
}

fastmr_wald_rows <- function(group, rows) {
  prepared <- fastmr_prepare_vectors(rows)
  x <- prepared$beta.exposure
  y <- prepared$beta.outcome
  sy <- prepared$se.outcome
  beta <- y / x
  # Match both TwoSampleMR::mr_wald_ratio() and fastMR's native Wald path:
  # the standard error treats the exposure estimate as fixed.
  se <- sy / abs(x)
  p <- rep(NA_real_, length(beta))
  valid <- is.finite(beta) & is.finite(se) & se > 0
  p[valid] <- 2 * stats::pnorm(abs(beta[valid] / se[valid]), lower.tail = FALSE)
  data.frame(
    exposure = rep(group$exposure, length(beta)),
    outcome = rep(group$outcome, length(beta)),
    id.exposure = rep(group$id.exposure, length(beta)),
    id.outcome = rep(group$id.outcome, length(beta)),
    samplesize = rep(fastmr_diagnostic_sample_size(rows), length(beta)),
    SNP = as.character(rows$SNP),
    b = beta,
    se = se,
    p = p,
    stringsAsFactors = FALSE
  )
}

#' Calculate single-SNP MR estimates and aggregate estimates
#'
#' @param data A harmonised TwoSampleMR-style data frame.
#' @param single_method The single-SNP method; defaults to `wald_ratio`.
#' @param all_method Methods used for the aggregate `All - ...` rows.
#' @param threads Maximum native worker count passed to [fast_mr()].
#' @return A tidy data frame compatible with `TwoSampleMR::mr_singlesnp()`.
#' @export
fast_mr_singlesnp <- function(data, single_method = "wald_ratio",
                              all_method = c("ivw", "egger"), threads = 1) {
  single_method <- fastmr_normalize_methods(single_method)
  all_method <- fastmr_normalize_methods(all_method)
  if (length(single_method) != 1L || single_method != "wald_ratio") {
    stop("single_method must be the wald_ratio method", call. = FALSE)
  }
  groups <- fastmr_diagnostic_groups(data)
  rows <- list()
  k <- 0L
  registry <- fastmr_method_registry()
  for (group in groups) {
    keep <- fastmr_diagnostic_keep(group$data)
    snp <- as.character(group$data$SNP)
    snp[is.na(snp)] <- ""
    selected <- which(keep & !duplicated(snp))
    if (length(selected)) {
      single_rows <- fastmr_wald_rows(group, group$data[selected, , drop = FALSE])
      for (index in seq_len(nrow(single_rows))) {
      k <- k + 1L
        rows[[k]] <- single_rows[index, , drop = FALSE]
      }
    }
    aggregate <- fast_mr(group$data, methods = all_method, nboot = 0,
                         threads = threads)
    for (method in all_method) {
      k <- k + 1L
      result <- aggregate[aggregate$method_code == method, , drop = FALSE]
      rows[[k]] <- data.frame(
        exposure = group$exposure,
        outcome = group$outcome,
        id.exposure = group$id.exposure,
        id.outcome = group$id.outcome,
        samplesize = fastmr_diagnostic_sample_size(group$data),
        SNP = paste("All -", registry$method[match(method, registry$code)]),
        b = fastmr_scalar(result, "b"),
        se = fastmr_scalar(result, "se"),
        p = fastmr_scalar(result, "pval"),
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

fastmr_leaveoneout_regression <- function(group, method) {
  keep <- fastmr_diagnostic_keep(group$data)
  snp <- as.character(group$data$SNP)
  snp[is.na(snp)] <- ""
  selected <- which(keep & !duplicated(snp))
  rows <- group$data[selected, , drop = FALSE]
  prepared <- fastmr_prepare_vectors(rows)
  n <- nrow(rows)
  if (n == 0L) return(data.frame())
  x <- prepared$beta.exposure
  y <- prepared$beta.outcome
  sx <- prepared$se.exposure
  sy <- prepared$se.outcome
  if (method == "uwr") {
    w <- rep(1, n)
  } else {
    w <- 1 / (sy * sy)
  }
  sum_w <- sum(w)
  sum_wx <- sum(w * x)
  sum_wxx <- sum(w * x * x)
  sum_wy <- sum(w * y)
  sum_wxy <- sum(w * x * y)
  sum_wyy <- sum(w * y * y)
  denominator <- sum_wxx - w * x * x
  numerator <- sum_wxy - w * x * y
  y_sum <- sum_wyy - w * y * y
  beta <- rep(NA_real_, n)
  se <- rep(NA_real_, n)
  p <- rep(NA_real_, n)
  valid <- n > 2L & is.finite(denominator) & denominator > 0
  beta[valid] <- numerator[valid] / denominator[valid]
  rss <- y_sum - numerator * numerator / denominator
  df <- n - 2L
  sigma <- sqrt(pmax(0, rss / df))
  base_se <- sqrt(1 / denominator)
  residual_se <- base_se * sigma
  if (method == "ivw_fe") {
    se[valid] <- base_se[valid]
  } else if (method == "ivw_mre") {
    se[valid] <- residual_se[valid]
  } else {
    correction <- pmin(1, sigma)
    se[valid] <- residual_se[valid] / correction[valid]
  }
  valid_p <- valid & is.finite(beta) & is.finite(se) & se > 0
  p[valid_p] <- 2 * stats::pnorm(abs(beta[valid_p] / se[valid_p]), lower.tail = FALSE)
  data.frame(
    exposure = rep(group$exposure, n),
    outcome = rep(group$outcome, n),
    id.exposure = rep(group$id.exposure, n),
    id.outcome = rep(group$id.outcome, n),
    samplesize = rep(fastmr_diagnostic_sample_size(rows), n),
    SNP = snp[selected],
    b = beta,
    se = se,
    p = p,
    stringsAsFactors = FALSE
  )
}

#' Calculate leave-one-SNP-out MR estimates
#'
#' @param data A harmonised TwoSampleMR-style data frame.
#' @param method A leave-one-out-capable method, normally `ivw` or `egger`.
#' @param threads Maximum native worker count passed to [fast_mr()].
#' @return A tidy data frame compatible with `TwoSampleMR::mr_leaveoneout()`.
#' @export
fast_mr_leaveoneout <- function(data, method = "ivw", threads = 1) {
  method <- fastmr_normalize_methods(method)
  if (length(method) != 1L || !method %in% c("ivw", "ivw_fe", "ivw_mre", "egger", "uwr")) {
    stop("method must be one heterogeneity-capable regression method", call. = FALSE)
  }
  groups <- fastmr_diagnostic_groups(data)
  rows <- list()
  k <- 0L
  for (group in groups) {
    keep <- fastmr_diagnostic_keep(group$data)
    snp <- as.character(group$data$SNP)
    snp[is.na(snp)] <- ""
    selected <- which(keep & !duplicated(snp))
    if (method %in% c("ivw", "ivw_fe", "ivw_mre", "uwr")) {
      leave_rows <- fastmr_leaveoneout_regression(group, method)
      for (index in seq_len(nrow(leave_rows))) {
      k <- k + 1L
        rows[[k]] <- leave_rows[index, , drop = FALSE]
      }
    } else {
      for (index in selected) {
        remaining <- group$data
        remaining <- remaining[!(as.character(remaining$SNP) == snp[[index]]), , drop = FALSE]
        result <- fast_mr(remaining, methods = method, nboot = 0,
                          threads = threads)
        k <- k + 1L
        rows[[k]] <- data.frame(
          exposure = group$exposure,
          outcome = group$outcome,
          id.exposure = group$id.exposure,
          id.outcome = group$id.outcome,
          samplesize = fastmr_diagnostic_sample_size(remaining),
          SNP = snp[[index]],
          b = fastmr_scalar(result, "b"),
          se = fastmr_scalar(result, "se"),
          p = fastmr_scalar(result, "pval"),
          stringsAsFactors = FALSE
        )
      }
    }
    result <- fast_mr(group$data, methods = method, nboot = 0,
                      threads = threads)
    k <- k + 1L
    rows[[k]] <- data.frame(
      exposure = group$exposure,
      outcome = group$outcome,
      id.exposure = group$id.exposure,
      id.outcome = group$id.outcome,
      samplesize = fastmr_diagnostic_sample_size(group$data),
      SNP = "All",
      b = fastmr_scalar(result, "b"),
      se = fastmr_scalar(result, "se"),
      p = fastmr_scalar(result, "pval"),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}
