# Batched multivariable IVW: one design (SNPs x exposures) against many
# outcomes, or many designs against one shared outcome panel.  The native
# kernel is fastmr_mvmr_batch_native() in src/mvmr.cpp.

fastmr_mvmr_se_models <- c("multiplicative", "multiplicative_floored", "fixed")

fastmr_mvmr_se_code <- function(se_model) {
  match(se_model, fastmr_mvmr_se_models) - 1L
}

fastmr_mvmr_check_flag <- function(value, argument) {
  if (length(value) != 1L || !is.logical(value) || is.na(value)) {
    stop(argument, " must be TRUE or FALSE", call. = FALSE)
  }
  value
}

# Per-panel-row and per-outcome factors of the outcome standard errors,
# se[i, k] ~= row_se[i] * scale[k], from an additive fit on the log scale
# over finite positive cells, with each outcome's largest absolute log
# deviation (NA for an outcome without finite cells).
fastmr_mvmr_se_factors <- function(outcome_se) {
  log_se <- log(outcome_se)
  log_se[!is.finite(log_se)] <- NA_real_
  row_log <- rowMeans(log_se, na.rm = TRUE)
  centred <- log_se - row_log
  scale_log <- colMeans(centred, na.rm = TRUE)
  deviation <- abs(sweep(centred, 2L, scale_log))
  max_dev <- suppressWarnings(apply(deviation, 2L, max, na.rm = TRUE))
  max_dev[!is.finite(max_dev)] <- NA_real_
  list(row_se = exp(row_log), scale = exp(scale_log), deviation = max_dev)
}

# Correlation input -> p x p x E array.  `cor` is NULL (zero off-diagonal),
# one p x p matrix for every design, or a list with one matrix per design.
fastmr_mvmr_cor_array <- function(cor, terms, design_count) {
  p <- length(terms)
  one <- function(m, label) {
    if (length(m) == 1L && p == 1L) m <- matrix(m, 1L, 1L)
    m <- as.matrix(m)
    storage.mode(m) <- "double"
    if (!all(dim(m) == p)) {
      stop(label, " must be a ", p, " x ", p, " correlation matrix", call. = FALSE)
    }
    if (!is.null(rownames(m)) && !is.null(colnames(m)) &&
        all(terms %in% rownames(m)) && all(terms %in% colnames(m))) {
      m <- m[terms, terms, drop = FALSE]
    }
    if (any(!is.finite(m)) || any(abs(m) > 1 + 1e-12) ||
        any(abs(diag(m) - 1) > 1e-12) || max(abs(m - t(m))) > 1e-12) {
      stop(label, " must be a symmetric correlation matrix with unit diagonal",
           call. = FALSE)
    }
    m
  }
  if (is.null(cor)) {
    return(array(rep(diag(p), design_count), c(p, p, design_count)))
  }
  if (is.list(cor) && !is.data.frame(cor)) {
    if (length(cor) != design_count) {
      stop("a list exposure_cor needs one matrix per design", call. = FALSE)
    }
    return(array(unlist(lapply(seq_along(cor), function(i) {
      one(cor[[i]], "exposure_cor")
    })), c(p, p, design_count)))
  }
  m <- one(cor, "exposure_cor")
  array(rep(m, design_count), c(p, p, design_count))
}

# Sanderson-Windmeijer conditional F statistic for each column of one design,
# as MVMR::strength_mvmr(): delta from the unweighted regression of column j
# on the other columns through the origin, then
# Q_j = sum_i (x_ij - x_i,-j' delta)^2 / var_i(delta) with
# var_i(delta) = Sigma_i[j, j] - 2 delta' Sigma_i[-j, j] + delta' Sigma_i[-j, -j] delta,
# Sigma_i = diag(se_i) R diag(se_i), and F_j = Q_j / (L - (p - 1)).
# With p = 1 it is the mean squared z statistic.
fastmr_mvmr_conditional_f <- function(beta, se, cor) {
  beta <- as.matrix(beta)
  se <- as.matrix(se)
  p <- ncol(beta)
  L <- nrow(beta)
  out <- rep(NA_real_, p)
  if (L <= p - 1L) return(out)
  if (p == 1L) return(sum((beta[, 1L] / se[, 1L])^2) / L)
  for (j in seq_len(p)) {
    others <- setdiff(seq_len(p), j)
    Z <- beta[, others, drop = FALSE]
    delta <- tryCatch(qr.solve(crossprod(Z), crossprod(Z, beta[, j])),
                      error = function(e) NULL)
    if (is.null(delta)) next
    delta <- drop(delta)
    residual <- beta[, j] - drop(Z %*% delta)
    sj <- se[, j]
    so <- se[, others, drop = FALSE]
    variance <- sj^2
    for (a in seq_along(others)) {
      variance <- variance - 2 * delta[a] * cor[others[a], j] * so[, a] * sj
      for (b in seq_along(others)) {
        variance <- variance + delta[a] * delta[b] * cor[others[a], others[b]] *
          so[, a] * so[, b]
      }
    }
    out[j] <- sum(residual^2 / variance) / (L - (p - 1L))
  }
  out
}

fastmr_mvmr_weak_warning <- function(weak, total, threshold, exposures = NULL) {
  if (!weak) return(invisible(NULL))
  shown <- if (length(exposures)) {
    paste0(" (", paste(utils::head(exposures, 10L), collapse = ", "),
           if (length(exposures) > 10L) ", ..." else "", ")")
  } else ""
  warning(
    weak, " of ", total, " design(s) have a conditional F statistic below ",
    threshold, " for at least one exposure", shown, ". The instruments are weak ",
    "conditional on the other exposures, so multivariable estimates are biased ",
    "toward the confounded (observational) association, and in one-sample or ",
    "overlapping-sample designs inflated rather than diluted. Consider ",
    "method = \"residualised\", stronger or exposure-specific instruments, or a ",
    "weak-instrument-robust estimator.", call. = FALSE
  )
}

# Validated native call on a stacked batch.  `rows` are 1-based panel rows.
fastmr_mvmr_engine <- function(row_ptr, rows, design, outcome_beta, outcome_se,
                               se_model, threads, design_se = NULL,
                               cor_array = NULL, weights = "exact",
                               shared_tolerance = 1e-8, return_vcov = FALSE) {
  shared <- NULL
  if (identical(weights, "shared")) {
    shared <- fastmr_mvmr_se_factors(outcome_se)
    # The shared path only checks outcome betas, so an outcome with a beta
    # whose standard error is invalid stays on the exact path.
    invalid_se <- colSums(is.finite(outcome_beta) &
                            !(is.finite(outcome_se) & outcome_se > 0)) > 0
    shared$use <- !is.na(shared$deviation) & shared$deviation <= shared_tolerance &
      is.finite(shared$scale) & !invalid_se
    row_se <- shared$row_se
    row_se[!is.finite(row_se)] <- NA_real_
  }
  native <- fastmr_mvmr_batch_native(
    as.integer(row_ptr), as.integer(rows) - 1L, design, outcome_beta, outcome_se,
    fastmr_mvmr_se_code(se_model), as.integer(threads),
    design_se = design_se,
    design_cor = if (is.null(design_se)) NULL else as.numeric(cor_array),
    shared_row_se = if (is.null(shared)) NULL else row_se,
    shared_outcome_scale = if (is.null(shared)) NULL else shared$scale,
    shared_outcome = if (is.null(shared)) NULL else shared$use,
    return_vcov = return_vcov
  )
  native$shared_deviation <- if (is.null(shared)) NULL else shared$deviation
  native
}

fastmr_mvmr_pval <- function(b, se) {
  out <- 2 * stats::pnorm(abs(b / se), lower.tail = FALSE)
  out[!is.finite(b) | !is.finite(se) | se == 0] <- NA_real_
  out
}

fastmr_mvmr_validate_common <- function(se_model, weights, shared_tolerance,
                                        threads, weak_f) {
  controls <- fastmr_validate_controls(0L, NULL, threads)
  if (length(shared_tolerance) != 1L || !is.numeric(shared_tolerance) ||
      is.na(shared_tolerance) || shared_tolerance < 0) {
    stop("shared_tolerance must be one non-negative number", call. = FALSE)
  }
  if (length(weak_f) != 1L || !is.numeric(weak_f) || is.na(weak_f)) {
    stop("weak_f must be one number", call. = FALSE)
  }
  controls$threads
}

fastmr_mvmr_outcome_panel <- function(outcome_beta, outcome_se, snp_count) {
  if (is.null(dim(outcome_beta))) outcome_beta <- matrix(outcome_beta, ncol = 1L)
  if (is.null(dim(outcome_se))) outcome_se <- matrix(outcome_se, ncol = 1L)
  outcome_beta <- fastmr_matrix_numeric(outcome_beta, "outcome_beta")
  outcome_se <- fastmr_matrix_numeric(outcome_se, "outcome_se")
  if (!identical(dim(outcome_beta), dim(outcome_se))) {
    stop("outcome_beta and outcome_se must have the same dimensions", call. = FALSE)
  }
  if (!is.null(snp_count) && nrow(outcome_beta) != snp_count) {
    stop("outcome matrices must have one row per SNP of exposure_beta", call. = FALSE)
  }
  if (!ncol(outcome_beta)) stop("at least one outcome is required", call. = FALSE)
  list(beta = outcome_beta, se = outcome_se)
}

#' Multivariable IVW for one design against many outcomes
#'
#' Fits the multivariable inverse-variance-weighted model of TwoSampleMR's
#' `mv_ivw()` / `mv_multiple()` (weighted regression through the origin of
#' the outcome effects on the exposure effects, weights `1 / se_y^2`) for
#' every column of an outcome matrix in one native call. With outcome-specific
#' weights each outcome has its own `p x p` normal equations; the kernel
#' accumulates them in one pass over the SNPs, solves them by Cholesky and
#' takes a second pass for the residual sum of squares, in parallel over
#' outcomes.
#'
#' SNPs whose outcome beta or standard error is missing or invalid for an
#' outcome are dropped for that outcome only, so `nsnp` can differ by outcome
#' (set a cell to `NA` to mask a SNP for one outcome, e.g. its cis window).
#' The exposure effects must be finite.
#'
#' **Standard errors.** `se_model = "multiplicative"` (default) is the
#' TwoSampleMR `mv_ivw()` / `mv_multiple()` and [fast_mr_multivariable()]
#' convention: the `lm()` standard error, i.e. the fixed-effect standard error
#' times the residual standard error `sigma = sqrt(Q / (nsnp - p))`, *not*
#' floored at one. `"multiplicative_floored"` multiplies by `max(1, sigma)`,
#' the univariable TwoSampleMR `mr_ivw()` convention (and the one used by
#' [fast_mr()]'s IVW). `"fixed"` uses the fixed-effect standard error.
#' P-values are two-sided normal, as in TwoSampleMR.
#'
#' **Degenerate fits.** With fewer valid SNPs than exposures, or a design
#' whose weighted cross-product has a Cholesky pivot below `1e-12` of its
#' diagonal (collinear columns), every estimate is `NA`. With exactly `p`
#' SNPs the estimates are reported (the fit is exact, as `lm()`), `Q`,
#' `sigma` and `Q_A` are `NA`, and the standard error is `NA` for
#' `"multiplicative"` and the fixed-effect value for the other two models.
#' ([fast_mr_multivariable()] returns `NA` estimates when `nsnp <= p`.)
#'
#' **Shared weights.** `weights = "shared"` is an optional fast path for
#' outcome standard errors that factor as `se[i, k] = s_i * c_k` (proportional
#' across SNPs). The factors are fitted on the log scale; an outcome uses the
#' shared path when every log standard error is within `shared_tolerance` of
#' its fitted value and all of its SNPs are present, and the exact path
#' otherwise. Under exact proportionality the two paths agree to rounding;
#' with a larger tolerance the shared path is an approximation (weights within
#' a factor of about `exp(2 * shared_tolerance)`), so leave it `"exact"`
#' unless speed matters more than the last digits.
#'
#' **Diagnostics.** When `exposure_se` is supplied, the result includes the
#' Sanderson-Windmeijer conditional F statistic of each exposure (as
#' `MVMR::strength_mvmr()`, computed on all rows of the design) and Sanderson's
#' `Q_A` heterogeneity statistic for each outcome (as `MVMR::pleiotropy_mvmr()`
#' but with the fitted IVW coefficients: residuals divided by
#' `se_y^2 + b' Sigma_i b`, `nsnp - p` degrees of freedom). Both need the
#' covariance of each SNP's exposure effects, `Sigma_i = diag(se_i) R
#' diag(se_i)`, where `R` (`exposure_cor`) is the phenotypic correlation of the
#' exposures scaled by their sample overlap. If it is absent it is taken as the
#' identity (independent exposure samples) with a warning; this overstates the
#' conditional F when the exposure GWAS share participants. A warning is given
#' when any conditional F is below `weak_f`.
#'
#' @param exposure_beta Numeric SNP x exposure matrix (all finite). Column
#'   names name the exposures.
#' @param outcome_beta,outcome_se Numeric SNP x outcome matrices (or vectors
#'   for one outcome).
#' @param exposure_se Optional SNP x exposure standard errors for the
#'   diagnostics.
#' @param exposure_cor Optional exposure correlation matrix (see Details).
#' @param se_model `"multiplicative"`, `"multiplicative_floored"` or `"fixed"`.
#' @param weights `"exact"` (default) or `"shared"`.
#' @param shared_tolerance Largest absolute log-scale deviation from
#'   proportional standard errors accepted by the shared path.
#' @param threads Native worker count.
#' @param weak_f Conditional F threshold for the weak-instrument warning
#'   (`-Inf` silences it).
#' @return A list of class `fastmr_mvmr` with `b`, `se` and `pval` (exposure x
#'   outcome matrices), `nsnp`, `Q`, `Q_df`, `Q_pval`, `sigma`, `Q_A`,
#'   `Q_A_pval` (one value per outcome), `conditional_F` (one per exposure, or
#'   `NULL`), `shared` (logical, outcomes fitted on the shared path) and
#'   `se_model`.
#' @seealso [fast_mvmr_ivw_batch()], [fast_mvmr_compressed()]
#' @export
fast_mvmr_ivw <- function(exposure_beta, outcome_beta, outcome_se,
                          exposure_se = NULL, exposure_cor = NULL,
                          se_model = c("multiplicative", "multiplicative_floored", "fixed"),
                          weights = c("exact", "shared"), shared_tolerance = 1e-8,
                          threads = 1L, weak_f = 10) {
  se_model <- match.arg(se_model)
  weights <- match.arg(weights)
  if (is.null(dim(exposure_beta))) exposure_beta <- matrix(exposure_beta, ncol = 1L)
  X <- fastmr_matrix_numeric(exposure_beta, "exposure_beta")
  if (is.null(colnames(X))) colnames(X) <- paste0("exposure", seq_len(ncol(X)))
  batch <- fast_mvmr_ivw_batch(
    list(list(rows = seq_len(nrow(X)), beta = X, se = exposure_se)),
    outcome_beta, outcome_se, exposure_cor = exposure_cor, se_model = se_model,
    weights = weights, shared_tolerance = shared_tolerance, threads = threads,
    weak_f = weak_f, panel_rows = nrow(X)
  )
  one <- function(a) {
    m <- t(matrix(a[1L, , ], dim(a)[2L], dim(a)[3L]))
    dimnames(m) <- list(dimnames(a)[[3L]], dimnames(a)[[2L]])
    m
  }
  vec <- function(m) stats::setNames(m[1L, ], colnames(m))
  structure(list(
    b = one(batch$b), se = one(batch$se), pval = one(batch$pval),
    nsnp = vec(batch$nsnp), Q = vec(batch$Q), Q_df = vec(batch$Q_df),
    Q_pval = vec(batch$Q_pval), sigma = vec(batch$sigma),
    Q_A = vec(batch$Q_A), Q_A_pval = vec(batch$Q_A_pval),
    conditional_F = if (is.null(batch$conditional_F)) NULL else
      stats::setNames(batch$conditional_F[1L, ], colnames(X)),
    shared = vec(batch$shared), se_model = se_model
  ), class = "fastmr_mvmr")
}

#' Batched multivariable IVW over many designs and a shared outcome panel
#'
#' The batched form of [fast_mvmr_ivw()]: `designs` holds one design per
#' primary exposure (its instrument rows of a shared outcome panel and the
#' SNP x exposure effect matrix at those rows), and every design is fitted
#' against every outcome column of the panel in one native call, in parallel
#' over (design, outcome block) pairs. All designs must have the same
#' exposure columns (e.g. the primary exposure followed by the same covariate
#' traits); their row sets may differ. Conventions are those of
#' [fast_mvmr_ivw()].
#'
#' @param designs A list of designs, each a list with `rows` (panel row
#'   indices, 1-based, or row names of `outcome_beta`), `beta` (an
#'   `length(rows)` x p matrix of exposure effects) and optionally `se` (same
#'   shape) and `cor` (a p x p exposure correlation overriding
#'   `exposure_cor`).
#' @param outcome_beta,outcome_se SNP x outcome panel matrices.
#' @param exposure_cor Optional p x p correlation used by every design without
#'   its own `cor`.
#' @param se_model,weights,shared_tolerance,threads,weak_f As in
#'   [fast_mvmr_ivw()].
#' @param return_vcov If `TRUE`, also return the coefficient covariance
#'   matrices (`vcov`, design x outcome x p x p).
#' @param panel_rows Internal; expected panel row count.
#' @return A list with `b`, `se`, `pval` (design x outcome x exposure arrays),
#'   `nsnp`, `Q`, `Q_df`, `Q_pval`, `sigma`, `Q_A`, `Q_A_pval`, `shared`
#'   (design x outcome matrices), `conditional_F` (design x exposure, or
#'   `NULL`), `shared_deviation` (per outcome, with `weights = "shared"`) and
#'   `se_model`.
#' @export
fast_mvmr_ivw_batch <- function(designs, outcome_beta, outcome_se,
                                exposure_cor = NULL,
                                se_model = c("multiplicative", "multiplicative_floored", "fixed"),
                                weights = c("exact", "shared"),
                                shared_tolerance = 1e-8, threads = 1L, weak_f = 10,
                                return_vcov = FALSE, panel_rows = NULL) {
  se_model <- match.arg(se_model)
  weights <- match.arg(weights)
  threads <- fastmr_mvmr_validate_common(se_model, weights, shared_tolerance,
                                         threads, weak_f)
  return_vcov <- fastmr_mvmr_check_flag(return_vcov, "return_vcov")
  panel <- fastmr_mvmr_outcome_panel(outcome_beta, outcome_se, panel_rows)
  if (!is.list(designs) || !length(designs)) {
    stop("designs must be a non-empty list", call. = FALSE)
  }
  design_names <- names(designs)
  if (is.null(design_names)) design_names <- as.character(seq_along(designs))
  terms <- NULL
  has_se <- NULL
  blocks <- vector("list", length(designs))
  for (e in seq_along(designs)) {
    d <- designs[[e]]
    if (!is.list(d) || is.null(d$rows) || is.null(d$beta)) {
      stop("each design must be a list with rows and beta", call. = FALSE)
    }
    beta <- d$beta
    if (is.null(dim(beta))) beta <- matrix(beta, ncol = 1L)
    beta <- fastmr_matrix_numeric(beta, "design beta")
    rows <- d$rows
    if (is.character(rows)) {
      matched <- match(rows, rownames(panel$beta))
      if (anyNA(matched)) {
        stop("design ", design_names[e], " names rows absent from outcome_beta",
             call. = FALSE)
      }
      rows <- matched
    }
    if (!is.numeric(rows) || anyNA(rows) || any(rows < 1) ||
        any(rows > nrow(panel$beta)) || any(rows != floor(rows))) {
      stop("design ", design_names[e], " rows must address outcome_beta", call. = FALSE)
    }
    if (nrow(beta) != length(rows)) {
      stop("design ", design_names[e], " beta must have one row per design row",
           call. = FALSE)
    }
    if (anyDuplicated(rows)) {
      stop("design ", design_names[e], " has duplicated rows", call. = FALSE)
    }
    if (any(!is.finite(beta))) {
      stop("design ", design_names[e], " exposure effects must be finite", call. = FALSE)
    }
    cn <- colnames(beta)
    if (is.null(cn)) cn <- paste0("exposure", seq_len(ncol(beta)))
    if (is.null(terms)) terms <- cn
    if (!identical(cn, terms)) {
      stop("every design must have the same exposure columns", call. = FALSE)
    }
    this_se <- !is.null(d$se)
    if (is.null(has_se)) has_se <- this_se
    if (!identical(has_se, this_se)) {
      stop("supply se for every design or for none", call. = FALSE)
    }
    se <- NULL
    if (this_se) {
      se <- d$se
      if (is.null(dim(se))) se <- matrix(se, ncol = 1L)
      se <- fastmr_matrix_numeric(se, "design se")
      if (!identical(dim(se), dim(beta)) || any(!is.finite(se) | se <= 0)) {
        stop("design ", design_names[e],
             " se must match beta and be finite and positive", call. = FALSE)
      }
    }
    blocks[[e]] <- list(rows = as.integer(rows), beta = beta, se = se, cor = d$cor)
  }
  p <- length(terms)
  design_count <- length(blocks)
  counts <- vapply(blocks, function(b) length(b$rows), integer(1))
  row_ptr <- c(0L, cumsum(counts))
  rows <- unlist(lapply(blocks, `[[`, "rows"), use.names = FALSE)
  design <- do.call(rbind, lapply(blocks, `[[`, "beta"))
  if (is.null(design)) design <- matrix(numeric(), 0L, p)
  design_se <- NULL
  cor_array <- NULL
  conditional_F <- NULL
  if (isTRUE(has_se)) {
    design_se <- do.call(rbind, lapply(blocks, `[[`, "se"))
    own <- vapply(blocks, function(b) !is.null(b$cor), logical(1))
    if (is.null(exposure_cor) && !all(own) && p > 1L) {
      warning("exposure_cor not supplied: assuming uncorrelated exposure effect ",
              "estimates (no sample overlap between exposure GWAS). The ",
              "conditional F and Q_A ignore that covariance and the conditional ",
              "F is overstated if the samples overlap.", call. = FALSE)
    }
    base <- fastmr_mvmr_cor_array(exposure_cor, terms, 1L)[, , 1L]
    cor_list <- lapply(blocks, function(b) if (is.null(b$cor)) base else b$cor)
    cor_array <- fastmr_mvmr_cor_array(cor_list, terms, design_count)
    conditional_F <- t(vapply(seq_len(design_count), function(e) {
      b <- blocks[[e]]
      if (!length(b$rows)) return(rep(NA_real_, p))
      fastmr_mvmr_conditional_f(b$beta, b$se, cor_array[, , e])
    }, numeric(p)))
    if (p == 1L) conditional_F <- matrix(conditional_F, ncol = 1L)
    dimnames(conditional_F) <- list(design_names, terms)
  }
  native <- fastmr_mvmr_engine(
    row_ptr, rows, design, panel$beta, panel$se, se_model, threads,
    design_se = design_se, cor_array = cor_array, weights = weights,
    shared_tolerance = shared_tolerance, return_vcov = return_vcov
  )
  K <- ncol(panel$beta)
  outcome_names <- colnames(panel$beta)
  if (is.null(outcome_names)) outcome_names <- as.character(seq_len(K))
  mat <- function(x) matrix(x, design_count, K, dimnames = list(design_names, outcome_names))
  arr <- function(x) array(x, c(design_count, K, p),
                           dimnames = list(design_names, outcome_names, terms))
  nsnp <- mat(native$nsnp)
  Q <- mat(native$Q)
  Q_df <- nsnp - p
  Q_df[!is.finite(Q)] <- NA_real_
  QA <- if (isTRUE(has_se)) mat(native$Q_A) else NULL
  chisq <- function(q) {
    if (is.null(q)) return(NULL)
    out <- stats::pchisq(q, Q_df, lower.tail = FALSE)
    out[!is.finite(q)] <- NA_real_
    out
  }
  b <- arr(native$beta)
  se <- arr(native$se)
  if (!is.null(conditional_F) && is.finite(weak_f)) {
    weak <- rowSums(conditional_F < weak_f, na.rm = TRUE) > 0
    fastmr_mvmr_weak_warning(sum(weak), design_count, weak_f,
                             if (design_count > 1L) design_names[weak] else NULL)
  }
  out <- list(
    b = b, se = se, pval = fastmr_mvmr_pval(b, se),
    nsnp = nsnp, Q = Q, Q_df = Q_df, Q_pval = chisq(Q), sigma = mat(native$sigma),
    Q_A = QA, Q_A_pval = chisq(QA),
    shared = mat(native$shared == 1L),
    conditional_F = conditional_F,
    shared_deviation = if (is.null(native$shared_deviation)) NULL else
      stats::setNames(native$shared_deviation, outcome_names),
    se_model = se_model
  )
  if (return_vcov) {
    out$vcov <- array(native$vcov, c(design_count, K, p, p),
                      dimnames = list(design_names, outcome_names, terms, terms))
  }
  out
}

#' @export
print.fastmr_mvmr <- function(x, ...) {
  cat("fastMR multivariable IVW:", nrow(x$b), "exposure(s) x", ncol(x$b),
      "outcome(s); se_model =", x$se_model, "\n")
  invisible(x)
}
