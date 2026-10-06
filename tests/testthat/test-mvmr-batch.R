mvmr_fixture <- function(n = 60L, p = 3L, K = 25L, seed = 11L, na = 80L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p, 0, 0.1), n, p,
              dimnames = list(NULL, c("exposure", if (p > 1L) paste0("cov", seq_len(p - 1L)))))
  B <- X %*% matrix(rnorm(p * K, 0, 0.5), p, K) + matrix(rnorm(n * K, 0, 0.02), n, K)
  S <- matrix(runif(n * K, 0.01, 0.03), n, K)
  if (na) B[sample(length(B), na)] <- NA
  colnames(B) <- colnames(S) <- paste0("out", seq_len(K))
  Sx <- matrix(runif(n * p, 0.004, 0.012), n, p, dimnames = dimnames(X))
  list(X = X, B = B, S = S, Sx = Sx)
}

lm_reference <- function(X, y, s) {
  ok <- is.finite(y) & is.finite(s) & s > 0
  fit <- summary(stats::lm(y[ok] ~ 0 + X[ok, , drop = FALSE], weights = 1 / s[ok]^2))
  list(b = fit$coefficients[, 1], se = fit$coefficients[, 2], n = sum(ok),
       Q = sum(fit$residuals^2), sigma = fit$sigma)
}

rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300))

test_that("fast_mvmr_ivw matches lm(), fast_mr_multivariable() and mv_ivw for every outcome", {
  f <- mvmr_fixture()
  r <- fast_mvmr_ivw(f$X, f$B, f$S)
  expect_s3_class(r, "fastmr_mvmr")
  expect_identical(dim(r$b), c(3L, 25L))
  for (k in seq_len(ncol(f$B))) {
    ref <- lm_reference(f$X, f$B[, k], f$S[, k])
    expect_lt(rel(r$b[, k], ref$b), 1e-12)
    expect_lt(rel(r$se[, k], ref$se), 1e-12)
    expect_equal(r$nsnp[[k]], ref$n)
    expect_lt(rel(r$sigma[[k]], ref$sigma), 1e-12)
    expect_equal(r$Q_df[[k]], ref$n - 3)
    ok <- is.finite(f$B[, k])
    legacy <- fast_mr_multivariable(f$X[ok, ], f$B[ok, k], f$S[ok, k])
    expect_lt(rel(r$b[, k], legacy$b), 1e-12)
    expect_lt(rel(r$se[, k], legacy$se), 1e-12)
    expect_lt(rel(r$pval[, k], legacy$pval), 1e-10)
  }
  expect_equal(r$Q_pval, stats::pchisq(r$Q, r$Q_df, lower.tail = FALSE))
})

test_that("fast_mvmr_ivw matches TwoSampleMR mv_multiple() and mv_ivw() when installed", {
  skip_if_not(requireNamespace("TwoSampleMR", quietly = TRUE), "TwoSampleMR not installed")
  skip_if_not_installed("ggplot2")
  mv_multiple <- getExportedValue("TwoSampleMR", "mv_multiple")
  mv_ivw <- getExportedValue("TwoSampleMR", "mv_ivw")
  f <- mvmr_fixture(na = 0L, K = 6L)
  r <- fast_mvmr_ivw(f$X, f$B, f$S)
  P <- matrix(1e-10, nrow(f$X), ncol(f$X), dimnames = dimnames(f$X))
  for (k in seq_len(ncol(f$B))) {
    mvdat <- list(
      exposure_beta = f$X, exposure_pval = P, exposure_se = f$Sx,
      outcome_beta = f$B[, k], outcome_se = f$S[, k],
      expname = data.frame(id.exposure = colnames(f$X), exposure = colnames(f$X)),
      outname = data.frame(id.outcome = "o", outcome = "o")
    )
    a <- mv_multiple(mvdat)$result
    a <- a[match(colnames(f$X), a$id.exposure), ]
    b <- suppressWarnings(mv_ivw(mvdat))$result
    b <- b[match(colnames(f$X), b$id.exposure), ]
    expect_lt(rel(r$b[, k], a$b), 1e-12)
    expect_lt(rel(r$se[, k], a$se), 1e-12)
    expect_lt(rel(r$b[, k], b$b), 1e-12)
    expect_lt(rel(r$se[, k], b$se), 1e-12)
  }
})

test_that("standard-error models: floored and fixed", {
  f <- mvmr_fixture(na = 0L)
  m <- fast_mvmr_ivw(f$X, f$B, f$S)
  fl <- fast_mvmr_ivw(f$X, f$B, f$S, se_model = "multiplicative_floored")
  fx <- fast_mvmr_ivw(f$X, f$B, f$S, se_model = "fixed")
  expect_equal(fl$b, m$b)
  sigma <- matrix(m$sigma, nrow(m$se), ncol(m$se), byrow = TRUE)
  expect_equal(fx$se, m$se / sigma, tolerance = 1e-12)
  expect_equal(fl$se, m$se / sigma * pmax(1, sigma), tolerance = 1e-12)
})

test_that("the batch equals per-design fits and does not depend on threads", {
  f <- mvmr_fixture(n = 80L, K = 40L)
  designs <- list(
    a = list(rows = 1:30, beta = f$X[1:30, ], se = f$Sx[1:30, ]),
    b = list(rows = c(70:41, 5), beta = f$X[c(70:41, 5), ], se = f$Sx[c(70:41, 5), ]),
    c = list(rows = rownames(f$B)[1:2], beta = f$X[1:2, ], se = f$Sx[1:2, ])
  )
  rownames(f$B) <- rownames(f$S) <- paste0("snp", seq_len(nrow(f$B)))
  designs$c$rows <- rownames(f$B)[1:2]
  one <- suppressWarnings(fast_mvmr_ivw_batch(designs, f$B, f$S, threads = 1))
  four <- suppressWarnings(fast_mvmr_ivw_batch(designs, f$B, f$S, threads = 4))
  expect_identical(one, four)
  for (e in 1:2) {
    d <- designs[[e]]
    single <- suppressWarnings(fast_mvmr_ivw(d$beta, f$B[d$rows, ], f$S[d$rows, ],
                                             exposure_se = d$se))
    expect_equal(t(one$b[e, , ]), single$b, ignore_attr = TRUE)
    expect_equal(one$Q_A[e, ], single$Q_A, ignore_attr = TRUE)
    expect_equal(one$conditional_F[e, ], single$conditional_F, ignore_attr = TRUE)
  }
  # Design c has fewer SNPs than exposures: everything NA.
  expect_true(all(is.na(one$b["c", , ])))
  expect_true(all(one$nsnp["c", ] <= 2))
})

test_that("degenerate designs: collinear, n <= p, single SNP, empty outcome", {
  f <- mvmr_fixture(na = 0L, K = 3L)
  X <- f$X
  X[, 3] <- X[, 1] * 2 - X[, 2]
  r <- fast_mvmr_ivw(X, f$B, f$S)
  expect_true(all(is.na(r$b)))
  expect_true(all(is.na(r$se)))
  # Exactly p SNPs: estimates exact (as lm), Q/sigma undefined.
  rows <- 1:3
  exact <- fast_mvmr_ivw(f$X[rows, ], f$B[rows, ], f$S[rows, ])
  expect_equal(exact$b[, 1], drop(solve(f$X[rows, ], f$B[rows, 1])), tolerance = 1e-10)
  expect_true(all(is.na(exact$se)))
  expect_true(all(is.na(exact$Q)))
  fixed <- fast_mvmr_ivw(f$X[rows, ], f$B[rows, ], f$S[rows, ], se_model = "fixed")
  expect_true(all(is.finite(fixed$se)))
  # Fewer SNPs than exposures.
  short <- fast_mvmr_ivw(f$X[1:2, ], f$B[1:2, ], f$S[1:2, ])
  expect_true(all(is.na(short$b)))
  # One exposure, one SNP: the Wald ratio with its first-order se (floored).
  wald <- fast_mvmr_ivw(f$X[1, 1], f$B[1, 1], f$S[1, 1], se_model = "multiplicative_floored")
  expect_equal(wald$b[[1]], unname(f$B[1, 1] / f$X[1, 1]))
  expect_equal(wald$se[[1]], unname(f$S[1, 1] / abs(f$X[1, 1])))
  ivw <- fast_mr(data.frame(SNP = "s", beta.exposure = f$X[1, 1], beta.outcome = f$B[1, 1],
                            se.exposure = 0.01, se.outcome = f$S[1, 1],
                            id.exposure = "e", id.outcome = "o"), methods = "ivw")
  expect_equal(wald$b[[1]], ivw$b)
  expect_equal(wald$se[[1]], ivw$se)
  # An outcome with no valid SNP.
  B <- f$B
  B[, 2] <- NA
  S <- f$S
  S[1:10, 3] <- -1
  r <- fast_mvmr_ivw(f$X, B, S)
  expect_equal(r$nsnp[[2]], 0)
  expect_true(all(is.na(r$b[, 2])))
  expect_equal(r$nsnp[[3]], nrow(f$X) - 10)
  expect_error(fast_mvmr_ivw(replace(f$X, 1, NA), f$B, f$S), "finite")
})

test_that("p = 1 multiplicative-floored fits equal fast_mr() IVW", {
  f <- mvmr_fixture(na = 0L, K = 4L, p = 1L)
  r <- fast_mvmr_ivw(f$X, f$B, f$S, se_model = "multiplicative_floored")
  for (k in 1:4) {
    d <- data.frame(SNP = paste0("s", seq_len(nrow(f$X))), beta.exposure = f$X[, 1],
                    beta.outcome = f$B[, k], se.exposure = f$Sx[, 1], se.outcome = f$S[, k],
                    id.exposure = "e", id.outcome = "o")
    ivw <- fast_mr(d, methods = "ivw")
    expect_lt(rel(r$b[1, k], ivw$b), 1e-12)
    expect_lt(rel(r$se[1, k], ivw$se), 1e-12)
    expect_lt(rel(r$Q[[k]], ivw$Q), 1e-10)
  }
})

test_that("shared-weight path equals the exact path for proportional SEs", {
  f <- mvmr_fixture(na = 0L)
  S <- outer(runif(nrow(f$X), 0.01, 0.03), runif(ncol(f$B), 0.5, 2))
  dimnames(S) <- dimnames(f$B)
  B <- f$B
  B[3, 5] <- NA  # outcome 5 falls back to the exact path
  exact <- fast_mvmr_ivw(f$X, B, S, exposure_se = f$Sx, exposure_cor = diag(3))
  shared <- fast_mvmr_ivw(f$X, B, S, exposure_se = f$Sx, exposure_cor = diag(3),
                          weights = "shared")
  expect_false(shared$shared[[5]])
  expect_true(all(shared$shared[-5]))
  expect_false(any(exact$shared))
  for (name in c("b", "se", "Q", "sigma", "Q_A")) {
    expect_lt(rel(shared[[name]], exact[[name]]), 1e-12)
  }
  # Non-proportional SEs: the shared path is declined at the default tolerance.
  declined <- fast_mvmr_ivw(f$X, f$B, f$S, weights = "shared")
  expect_false(any(declined$shared))
  expect_equal(declined$b, fast_mvmr_ivw(f$X, f$B, f$S)$b)
  approx <- fast_mvmr_ivw(f$X, f$B, f$S, weights = "shared", shared_tolerance = Inf)
  expect_true(all(approx$shared))
})

# Transcription of MVMR::strength_mvmr() / pleiotropy_mvmr() (gencov as a
# list of per-SNP covariance matrices), with the L - (p - 1) divisor.
mvmr_reference <- function(X, Sx, y, sy, R) {
  L <- nrow(X); p <- ncol(X)
  Sig <- lapply(seq_len(L), function(i) diag(Sx[i, ]) %*% R %*% diag(Sx[i, ]))
  Fs <- numeric(p)
  for (j in seq_len(p)) {
    delta <- stats::lm.fit(X[, -j, drop = FALSE], X[, j])$coefficients
    d <- numeric(p); d[j] <- -1; d[-j] <- delta
    v <- vapply(Sig, function(S) drop(t(d) %*% S %*% d), numeric(1))
    Fs[j] <- sum((X[, j] - X[, -j, drop = FALSE] %*% delta)^2 / v) / (L - (p - 1))
  }
  A <- stats::lm(y ~ 0 + X, weights = 1 / sy^2)$coefficients
  vA <- sy^2 + vapply(Sig, function(S) drop(t(A) %*% S %*% A), numeric(1))
  list(F = Fs, QA = sum((y - X %*% A)^2 / vA))
}

test_that("conditional F and Q_A match the MVMR package formulas", {
  f <- mvmr_fixture(na = 0L, K = 3L)
  R <- matrix(c(1, 0.3, 0.1, 0.3, 1, -0.2, 0.1, -0.2, 1), 3,
              dimnames = list(colnames(f$X), colnames(f$X)))
  r <- fast_mvmr_ivw(f$X, f$B, f$S, exposure_se = f$Sx, exposure_cor = R)
  for (k in 1:3) {
    ref <- mvmr_reference(f$X, f$Sx, f$B[, k], f$S[, k], R)
    expect_lt(rel(r$conditional_F, ref$F), 1e-10)
    expect_lt(rel(r$Q_A[[k]], ref$QA), 1e-10)
  }
  expect_equal(r$Q_A_pval, stats::pchisq(r$Q_A, r$Q_df, lower.tail = FALSE))
  # Correlation given in another order is matched by name.
  shuffled <- R[c(3, 1, 2), c(3, 1, 2)]
  expect_equal(fast_mvmr_ivw(f$X, f$B, f$S, exposure_se = f$Sx,
                             exposure_cor = shuffled)$conditional_F,
               r$conditional_F)
})

test_that("missing correlation and weak conditional instruments warn", {
  f <- mvmr_fixture(na = 0L, K = 2L)
  expect_warning(fast_mvmr_ivw(f$X, f$B, f$S, exposure_se = f$Sx), "exposure_cor not supplied")
  # Nearly collinear exposures: weak conditional instruments.
  X <- f$X
  X[, 2] <- X[, 1] + rnorm(nrow(X), 0, 0.002)
  expect_warning(
    fast_mvmr_ivw(X, f$B, f$S, exposure_se = f$Sx, exposure_cor = diag(3)),
    "conditional F statistic below 10.*biased toward the confounded"
  )
  expect_silent(fast_mvmr_ivw(X, f$B, f$S, exposure_se = f$Sx, exposure_cor = diag(3),
                              weak_f = -Inf))
  expect_silent(fast_mvmr_ivw(f$X, f$B, f$S))
  expect_error(fast_mvmr_ivw(f$X, f$B, f$S, exposure_se = f$Sx,
                             exposure_cor = matrix(2, 3, 3)), "correlation")
})
