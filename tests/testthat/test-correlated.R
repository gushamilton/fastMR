sim_correlated <- function(n = 12, rho = 0.9, block = 4, seed = 1, theta = 0.3, dep = TRUE) {
  set.seed(seed)
  R <- diag(n)
  idx <- split(seq_len(n), ceiling(seq_len(n) / block))
  for (g in idx) R[g, g] <- rho^abs(outer(seq_along(g), seq_along(g), "-"))
  # random sign flips mimic exposure-allele alignment of the LD matrix
  s <- sample(c(-1, 1), n, TRUE)
  R <- R * (s %o% s)
  sx <- runif(n, 0.01, 0.03)
  sy <- runif(n, 0.01, 0.05)
  bx <- abs(rnorm(n, 0.1, 0.05)) * s
  Sy <- (sy %o% sy) * R
  by <- theta * bx + drop(t(chol(Sy)) %*% rnorm(n))
  list(bx = bx, by = by, sx = sx, sy = sy, ld = R)
}

mr_input <- function(d) {
  MendelianRandomization::mr_input(bx = d$bx, bxse = d$sx, by = d$by, byse = d$sy,
                                   correlation = d$ld)
}

test_that("ivw_gls matches mr_ivw(correl = TRUE) fixed and random", {
  skip_if_not_installed("MendelianRandomization")
  for (seed in 1:5) {
    d <- sim_correlated(seed = seed, rho = c(0.3, 0.6, 0.9, 0.9, 0.5)[seed])
    obj <- mr_input(d)
    for (m in c("random", "fixed")) {
      ref <- MendelianRandomization::mr_ivw(obj, model = m)
      got <- fast_mr_correlated(d$bx, d$by, d$sx, d$sy, d$ld, "ivw_gls", model = m)
      expect_equal(got$b, ref@Estimate, tolerance = 1e-8)
      expect_equal(got$se, ref@StdError, tolerance = 1e-8)
      expect_equal(got$Q, ref@Heter.Stat[1], tolerance = 1e-8)
      expect_equal(got$Q_pval, ref@Heter.Stat[2], tolerance = 1e-8)
      expect_equal(got$F, ref@Fstat, tolerance = 1e-8)
    }
  }
})

test_that("egger_gls matches mr_egger(correl = TRUE)", {
  skip_if_not_installed("MendelianRandomization")
  for (seed in 1:5) {
    d <- sim_correlated(seed = seed + 10, rho = 0.8)
    ref <- MendelianRandomization::mr_egger(mr_input(d))
    got <- fast_mr_correlated(d$bx, d$by, d$sx, d$sy, d$ld, "egger_gls")
    expect_equal(got$b, ref@Estimate, tolerance = 1e-8)
    expect_equal(got$se, ref@StdError.Est, tolerance = 1e-8)
    expect_equal(got$intercept, ref@Intercept, tolerance = 1e-8)
    expect_equal(got$intercept_se, ref@StdError.Int, tolerance = 1e-8)
    expect_equal(got$pval, ref@Pvalue.Est, tolerance = 1e-8)
    expect_equal(got$Q, ref@Heter.Stat[1], tolerance = 1e-8)
  }
})

test_that("pc_ivw selects the same PCs as mr_pcgmm and matches the published formulae", {
  skip_if_not_installed("MendelianRandomization")
  for (seed in 1:4) {
    # mr_pcgmm uses |bx|, so orient every SNP to a positive exposure effect
    d <- sim_correlated(seed = seed + 20, rho = 0.9, n = 12)
    s <- sign(d$bx); d$bx <- abs(d$bx); d$by <- d$by * s; d$ld <- d$ld * (s %o% s)
    ref <- suppressMessages(MendelianRandomization::mr_pcgmm(mr_input(d), nx = 1e4, ny = 1e4,
                                                              robust = FALSE))
    got <- fast_mr_correlated(d$bx, d$by, d$sx, d$sy, d$ld, "pc_ivw", pc_center = TRUE)
    expect_equal(got$npc, ref@PCs)
    # independent transcription of Burgess et al. (2017)
    Phi <- ((d$bx / d$sy) %o% (d$bx / d$sy)) * d$ld
    pc <- eigen(Phi, symmetric = TRUE)
    K <- which(cumsum(pc$values) / sum(pc$values) >= 0.99)[1]
    W <- pc$vectors[, 1:K, drop = FALSE]
    Om <- t(W) %*% ((d$sy %o% d$sy) * d$ld) %*% W
    bx0 <- drop(t(W) %*% d$bx); by0 <- drop(t(W) %*% d$by)
    Oi <- solve(Om)
    th <- drop(t(bx0) %*% Oi %*% by0 / (t(bx0) %*% Oi %*% bx0))
    fx <- sqrt(1 / drop(t(bx0) %*% Oi %*% bx0))
    r <- by0 - th * bx0
    se <- fx * max(1, sqrt(drop(t(r) %*% Oi %*% r) / (K - 1)))
    got0 <- fast_mr_correlated(d$bx, d$by, d$sx, d$sy, d$ld, "pc_ivw")
    expect_equal(got0$npc, K)
    expect_equal(got0$b, th, tolerance = 1e-8)
    expect_equal(got0$se, se, tolerance = 1e-8)
  }
})

test_that("pc_ivw with all components equals GLS IVW", {
  d <- sim_correlated(seed = 3, rho = 0.7)
  a <- fast_mr_correlated(d$bx, d$by, d$sx, d$sy, d$ld, "ivw_gls")
  b <- fast_mr_correlated(d$bx, d$by, d$sx, d$sy, d$ld, "pc_ivw", n_pc = length(d$bx))
  expect_equal(b$b, a$b, tolerance = 1e-8)
  expect_equal(b$se, a$se, tolerance = 1e-8)
  expect_equal(b$Q, a$Q, tolerance = 1e-8)
})

test_that("identity LD reduces to standard IVW and Egger", {
  d <- sim_correlated(seed = 4, rho = 0)
  d$ld <- diag(length(d$bx))
  dat <- data.frame(SNP = paste0("rs", seq_along(d$bx)), beta.exposure = d$bx,
                    beta.outcome = d$by, se.exposure = d$sx, se.outcome = d$sy,
                    id.exposure = "x", id.outcome = "y")
  std <- fast_mr(dat, methods = c("ivw", "egger"), nboot = 0)
  got <- fast_mr_correlated(d$bx, d$by, d$sx, d$sy, d$ld, c("ivw_gls", "egger_gls"))
  iv <- std[std$method == grep("Inverse variance weighted", std$method, value = TRUE)[1], ]
  eg <- std[grepl("Egger", std$method) & !grepl("intercept", std$method), ][1, ]
  expect_equal(got$b[1], iv$b[1], tolerance = 1e-8)
  expect_equal(got$se[1], iv$se[1], tolerance = 1e-8)
  expect_equal(got$b[2], eg$b[1], tolerance = 1e-8)
  expect_equal(got$se[2], eg$se[1], tolerance = 1e-8)
})

test_that("inputs are validated", {
  d <- sim_correlated(seed = 5)
  f <- function(...) fast_mr_correlated(d$bx, d$by, d$sx, d$sy, d$ld, ...)
  expect_error(fast_mr_correlated(d$bx, d$by, d$sx, d$sy, d$ld[-1, -1]), "ld must be a")
  expect_error(fast_mr_correlated(d$bx, d$by[-1], d$sx, d$sy, d$ld), "same length")
  expect_error(fast_mr_correlated(d$bx, d$by, d$sx, d$sy, 2 * d$ld), "correlation matrix")
  nms <- paste0("rs", seq_along(d$bx))
  l <- d$ld; dimnames(l) <- list(nms, nms)
  expect_error(fast_mr_correlated(d$bx, d$by, d$sx, d$sy, l, snp = rev(nms)), "do not match")
  expect_silent(fast_mr_correlated(d$bx, d$by, d$sx, d$sy, l, snp = nms))
  expect_error(f(pc_threshold = 2), "pc_threshold")
})

test_that("near-singular LD is shrunk and recorded, or refused", {
  n <- 6; d <- sim_correlated(n = n, seed = 9, rho = 0.9, block = 6)
  # SNP 2 duplicates SNP 1 -> exactly singular LD
  Rs <- d$ld; Rs[2, ] <- Rs[1, ]; Rs[, 2] <- Rs[, 1]; Rs[2, 2] <- 1; Rs[1, 2] <- Rs[2, 1] <- 1
  expect_error(fast_mr_correlated(d$bx, d$by, d$sx, d$sy, Rs, "ivw_gls", ld_action = "error"),
               "not positive definite")
  got <- fast_mr_correlated(d$bx, d$by, d$sx, d$sy, Rs, c("ivw_gls", "pc_ivw"))
  expect_gt(got$ld_shrinkage[1], 0)
  expect_equal(got$ld_shrinkage[2], 0)
  expect_true(all(is.finite(got$b)))
})

test_that("fast_mr_ld_neff behaves on simple matrices", {
  expect_equal(unname(fast_mr_ld_neff(diag(5))["galwey"]), 5)
  ones <- matrix(0.999, 5, 5); diag(ones) <- 1
  expect_lt(unname(fast_mr_ld_neff(ones)["galwey"]), 1.2)
})
