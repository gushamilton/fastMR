#' Mendelian randomization with correlated (LD) instruments
#'
#' Summary-data MR methods that account for linkage disequilibrium between
#' the instruments, intended for cis-MR at a single locus (for example a
#' drug-target gene region). The methods follow Burgess, Dudbridge & Thompson
#' (2016) for generalised-least-squares (GLS) IVW and MR-Egger, and Burgess et
#' al. (2017) for principal-components (PC) IVW.
#'
#' With `omega = diag(sy) %*% ld %*% diag(sy)` the outcome-effect covariance
#' matrix:
#'
#' * `"ivw_gls"`: `b = (bx' W bx)^-1 bx' W by` with `W = omega^-1`;
#'   fixed-effect `se = (bx' W bx)^-1/2`.
#' * `"egger_gls"`: GLS regression of `sign(bx) * by` on `|bx|` with an
#'   intercept, `ld` sign-flipped by `sign(bx)`, as in
#'   `MendelianRandomization::mr_egger(correl = TRUE)`.
#' * `"pc_ivw"`: the instruments are replaced by the leading `K` principal
#'   components of `Phi = (bx/sy)(bx/sy)' * ld`, retaining the smallest `K`
#'   whose cumulative variance share reaches `pc_threshold`; IVW is then
#'   applied to the projected associations with covariance `W' omega W`.
#'   This avoids inverting a near-singular LD matrix.
#'
#' Random-effects (`model = "random"`, the default) multiplies the fixed-effect
#' standard error by `max(1, phi)`, where `phi = sqrt(Q/df)` and `Q` is the
#' LD-aware residual heterogeneity statistic (df = `n-1` for IVW, `n-2` for
#' Egger, `K-1` for PC-IVW), the same floor-at-one convention as TwoSampleMR
#' and MendelianRandomization. `model = "fixed"` applies no inflation.
#'
#' @section LD matrix:
#' `ld` must be the **signed** correlation matrix of the instruments (not
#' `r^2`), with rows and columns in the same order as the effect vectors and
#' **aligned to the exposure effect allele** of each SNP: `ld[i, j]` is the
#' correlation between the effect-allele dosages of SNPs `i` and `j` in the
#' alleles in which `bx` and `by` are expressed. Misaligned signs give wrong
#' answers silently, so align alleles before computing or flipping the matrix.
#' If `snp` or the names of `bx` are supplied they must match `dimnames(ld)`.
#'
#' @section Near-singular LD:
#' For `ivw_gls` and `egger_gls`, if the smallest eigenvalue of `ld` is below
#' `ld_tol`, `ld` is shrunk towards the identity,
#' `(1 - lambda) * ld + lambda * I`, with the smallest `lambda` that lifts the
#' minimum eigenvalue to `ld_tol`; the value used is returned in the
#' `ld_shrinkage` column (0 when untouched). Use `ld_tol = 0` or
#' `ld_action = "error"` to refuse instead. `pc_ivw` does not need the inverse
#' of `ld`, only of the retained-component covariance.
#'
#' @param bx,by Numeric SNP-exposure and SNP-outcome effects (same alleles).
#' @param sx,sy Their standard errors. `sx` is only used for the GLS F-statistic.
#' @param ld Signed SNP correlation matrix, exposure-allele aligned.
#' @param method One or more of `"ivw_gls"`, `"egger_gls"`, `"pc_ivw"`.
#' @param model `"random"` (default, multiplicative, floored at 1) or `"fixed"`.
#' @param pc_threshold Cumulative variance share to retain for `pc_ivw`.
#' @param n_pc Optional fixed number of principal components for `pc_ivw`.
#' @param pc_center If `TRUE`, replicate the component selection of
#'   `MendelianRandomization::mr_pcgmm` (`prcomp(Phi)`: columns of `Phi`
#'   centred, and variance shares computed from squared eigenvalues, which
#'   retains fewer components). Default `FALSE`: exact eigen-decomposition of
#'   the positive semi-definite `Phi`, with shares proportional to its
#'   eigenvalues (more conservative).
#' @param ld_tol Minimum acceptable eigenvalue of `ld` for GLS methods.
#' @param ld_action `"shrink"` (default) or `"error"` when `ld` is below `ld_tol`.
#' @param distribution `"normal"` (default) or `"t"` for p-values.
#' @param snp Optional SNP identifiers to check against `dimnames(ld)`.
#' @return A data frame with one row per method: `method`, `model`, `nsnp`, `b`,
#'   `se`, `se_fixed`, `pval`, `phi`, `Q`, `Q_df`, `Q_pval`, `F`, `npc`,
#'   `pc_var_explained`, `intercept`, `intercept_se`, `intercept_pval`,
#'   `ld_shrinkage`, `ld_min_eigen`.
#' @references
#' Burgess S, Dudbridge F, Thompson SG (2016). Combining information on
#' multiple instrumental variables in Mendelian randomization. Stat Med
#' 35:1880-1906.
#'
#' Burgess S, Zuber V, Valdes-Marquez E, Sun BB, Hopewell JC (2017). Mendelian
#' randomization with fine-mapped genetic data. Genet Epidemiol 41:714-725.
#' @export
fast_mr_correlated <- function(bx, by, sx, sy, ld,
                               method = c("ivw_gls", "egger_gls", "pc_ivw"),
                               model = c("random", "fixed"),
                               pc_threshold = 0.99, n_pc = NULL, pc_center = FALSE,
                               ld_tol = 1e-8, ld_action = c("shrink", "error"),
                               distribution = c("normal", "t"), snp = NULL) {
  method <- match.arg(method, c("ivw_gls", "egger_gls", "pc_ivw"), several.ok = TRUE)
  model <- match.arg(model)
  ld_action <- match.arg(ld_action)
  distribution <- match.arg(distribution)
  nm <- names(bx)
  bx <- fastmr_numeric(bx, "bx"); by <- fastmr_numeric(by, "by")
  sx <- fastmr_numeric(sx, "sx"); sy <- fastmr_numeric(sy, "sy")
  n <- length(bx)
  if (length(by) != n || length(sx) != n || length(sy) != n) {
    stop("bx, by, sx and sy must have the same length", call. = FALSE)
  }
  if (anyNA(c(bx, by, sx, sy)) || any(!is.finite(c(bx, by, sx, sy))) || any(sx <= 0) || any(sy <= 0)) {
    stop("bx, by, sx and sy must be finite with positive standard errors", call. = FALSE)
  }
  if (!is.numeric(pc_threshold) || length(pc_threshold) != 1L || pc_threshold <= 0 || pc_threshold > 1) {
    stop("pc_threshold must be in (0, 1]", call. = FALSE)
  }
  ld <- as.matrix(ld)
  if (!is.numeric(ld) || nrow(ld) != n || ncol(ld) != n) {
    stop(sprintf("ld must be a %d x %d numeric matrix matching the %d SNPs", n, n, n), call. = FALSE)
  }
  ids <- if (!is.null(snp)) as.character(snp) else nm
  if (!is.null(ids) && length(ids) != n) stop("snp must have one entry per SNP", call. = FALSE)
  if (!is.null(ids) && !is.null(rownames(ld))) {
    if (!identical(rownames(ld), ids) || !identical(colnames(ld), ids)) {
      stop("dimnames(ld) do not match the SNP identifiers; reorder/subset ld to the effect vectors", call. = FALSE)
    }
  }
  if (anyNA(ld) || max(abs(ld - t(ld))) > 1e-6) stop("ld must be finite and symmetric", call. = FALSE)
  if (max(abs(diag(ld) - 1)) > 1e-6 || max(abs(ld)) > 1 + 1e-6) {
    stop("ld must be a correlation matrix (unit diagonal, |r| <= 1), not r^2 or covariance", call. = FALSE)
  }
  ld <- (ld + t(ld)) / 2
  rows <- lapply(method, function(m) {
    switch(m,
      ivw_gls = fastmr_corr_gls(bx, by, sx, sy, ld, FALSE, model, ld_tol, ld_action, distribution),
      egger_gls = fastmr_corr_gls(bx, by, sx, sy, ld, TRUE, model, ld_tol, ld_action, distribution),
      pc_ivw = fastmr_corr_pc(bx, by, sx, sy, ld, model, pc_threshold, n_pc, pc_center, distribution))
  })
  do.call(rbind, rows)
}

fastmr_corr_blank <- function(method, model, n) {
  data.frame(method = method, model = model, nsnp = n, b = NA_real_, se = NA_real_,
             se_fixed = NA_real_, pval = NA_real_, phi = NA_real_, Q = NA_real_,
             Q_df = NA_real_, Q_pval = NA_real_, F = NA_real_, npc = NA_integer_,
             pc_var_explained = NA_real_, intercept = NA_real_, intercept_se = NA_real_,
             intercept_pval = NA_real_, ld_shrinkage = 0, ld_min_eigen = NA_real_,
             stringsAsFactors = FALSE)
}

fastmr_corr_pval <- function(z, df, distribution) {
  if (distribution == "t" && df > 0) 2 * stats::pt(-abs(z), df) else 2 * stats::pnorm(-abs(z))
}

# Smallest-shrinkage ridge to the identity lifting the minimum eigenvalue to tol.
fastmr_corr_condition <- function(ld, tol, action) {
  ev_min <- min(eigen(ld, symmetric = TRUE, only.values = TRUE)$values)
  lambda <- 0
  if (ev_min < tol) {
    if (action == "error") {
      stop(sprintf("ld is not positive definite at tolerance %g (min eigenvalue %.3g); use ld_action = \"shrink\", prune SNPs, or method = \"pc_ivw\"", tol, ev_min), call. = FALSE)
    }
    lambda <- (tol - ev_min) / (1 - ev_min)
    ld <- (1 - lambda) * ld + lambda * diag(nrow(ld))
  }
  list(ld = ld, lambda = lambda, min_eigen = ev_min)
}

fastmr_corr_gls <- function(bx, by, sx, sy, ld, egger, model, ld_tol, ld_action, distribution) {
  n <- length(bx)
  method <- if (egger) "egger_gls" else "ivw_gls"
  out <- fastmr_corr_blank(method, model, n)
  p <- if (egger) 2L else 1L
  if (n < p) return(out)
  cond <- fastmr_corr_condition(ld, ld_tol, ld_action)
  out$ld_shrinkage <- cond$lambda
  out$ld_min_eigen <- cond$min_eigen
  R <- cond$ld
  if (egger) {
    s <- ifelse(bx < 0, -1, 1)
    by <- by * s; bx <- abs(bx); R <- R * (s %o% s)
  }
  omega <- (sy %o% sy) * R
  oinv <- chol2inv(chol(omega))
  X <- if (egger) cbind(1, bx) else matrix(bx, ncol = 1L)
  xtw <- crossprod(X, oinv)
  cov_fixed <- solve(xtw %*% X)
  theta <- drop(cov_fixed %*% (xtw %*% by))
  res <- by - drop(X %*% theta)
  df <- n - p
  Q <- drop(crossprod(res, oinv %*% res))
  phi <- if (df > 0) sqrt(Q / df) else NA_real_
  infl <- if (model == "random" && is.finite(phi)) max(1, phi) else 1
  j <- p
  se_fixed <- sqrt(cov_fixed[j, j])
  se <- se_fixed * infl
  out$b <- theta[j]; out$se <- se; out$se_fixed <- se_fixed
  out$pval <- fastmr_corr_pval(theta[j] / se, df, distribution)
  out$phi <- phi
  if (df > 0) { out$Q <- Q; out$Q_df <- df; out$Q_pval <- stats::pchisq(Q, df, lower.tail = FALSE) }
  # GLS first-stage F as in mr_ivw(correl = TRUE)
  out$F <- sum((backsolve(chol(R), bx / sx, transpose = TRUE))^2) / n
  if (egger) {
    out$intercept <- theta[1]
    out$intercept_se <- sqrt(cov_fixed[1, 1]) * infl
    out$intercept_pval <- fastmr_corr_pval(theta[1] / out$intercept_se, df, distribution)
  }
  out
}

fastmr_corr_pc <- function(bx, by, sx, sy, ld, model, threshold, n_pc, center, distribution) {
  n <- length(bx)
  out <- fastmr_corr_blank("pc_ivw", model, n)
  phi_mat <- ((bx / sy) %o% (bx / sy)) * ld
  if (center) phi_mat <- scale(phi_mat, center = TRUE, scale = FALSE)
  e <- eigen(if (center) crossprod(phi_mat) else phi_mat, symmetric = TRUE)
  vals <- e$values
  vals <- pmax(vals, 0)
  share <- cumsum(vals) / sum(vals)
  K <- if (is.null(n_pc)) which(share >= threshold - 1e-12)[1] else as.integer(n_pc)
  if (is.na(K) || K < 1L || K > n) stop("invalid number of principal components", call. = FALSE)
  W <- e$vectors[, seq_len(K), drop = FALSE]
  bx0 <- drop(crossprod(W, bx)); by0 <- drop(crossprod(W, by))
  omega <- (sy %o% sy) * ld
  om0 <- crossprod(W, omega %*% W)
  om0 <- (om0 + t(om0)) / 2
  oinv <- tryCatch(chol2inv(chol(om0)), error = function(e) {
    stop("covariance of the retained principal components is singular; reduce pc_threshold or n_pc", call. = FALSE)
  })
  den <- drop(crossprod(bx0, oinv %*% bx0))
  theta <- drop(crossprod(bx0, oinv %*% by0)) / den
  res <- by0 - theta * bx0
  df <- K - 1L
  Q <- drop(crossprod(res, oinv %*% res))
  phi <- if (df > 0) sqrt(Q / df) else NA_real_
  infl <- if (model == "random" && is.finite(phi)) max(1, phi) else 1
  se_fixed <- sqrt(1 / den)
  se <- se_fixed * infl
  out$b <- theta; out$se <- se; out$se_fixed <- se_fixed
  out$pval <- fastmr_corr_pval(theta / se, df, distribution)
  out$phi <- phi
  if (df > 0) { out$Q <- Q; out$Q_df <- df; out$Q_pval <- stats::pchisq(Q, df, lower.tail = FALSE) }
  out$npc <- K
  out$pc_var_explained <- share[K]
  out$ld_min_eigen <- min(eigen(ld, symmetric = TRUE, only.values = TRUE)$values)
  out$F <- drop(crossprod(bx0, solve(crossprod(W, ((sx %o% sx) * ld) %*% W), bx0))) / K
  out
}

#' Effective number of independent signals in an LD matrix
#'
#' @param ld Signed SNP correlation matrix.
#' @param threshold Cumulative eigenvalue share for the principal-component count.
#' @return A named numeric vector: `nsnp`, `n_pc` (components of `ld` needed to
#'   reach `threshold`), `galwey` (Galwey 2009, `(sum sqrt(l))^2 / sum(l)` over
#'   positive eigenvalues) and `li_ji` (Li & Ji 2005).
#' @export
fast_mr_ld_neff <- function(ld, threshold = 0.99) {
  ld <- as.matrix(ld)
  if (nrow(ld) != ncol(ld)) stop("ld must be square", call. = FALSE)
  ev <- pmax(eigen((ld + t(ld)) / 2, symmetric = TRUE, only.values = TRUE)$values, 0)
  li_ji <- sum(as.numeric(ev >= 1) + (ev - floor(ev)))
  c(nsnp = nrow(ld),
    n_pc = which(cumsum(ev) / sum(ev) >= threshold - 1e-12)[1],
    galwey = sum(sqrt(ev))^2 / sum(ev),
    li_ji = li_ji)
}
