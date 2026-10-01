test_that("Steiger directionality returns a tidy per-pair result", {
  d <- diagnostic_fixture()
  d$r.exposure <- seq_len(nrow(d)) / 1000
  d$r.outcome <- seq_len(nrow(d)) / 1200
  d$samplesize.exposure <- 10000
  d$samplesize.outcome <- 12000
  result <- fast_mr_directionality_test(d)
  expect_equal(names(result), c("id.exposure", "id.outcome", "exposure", "outcome",
                                "snp_r2.exposure", "snp_r2.outcome",
                                "correct_causal_direction", "steiger_pval"))
  expect_true(result$correct_causal_direction)
  expect_true(is.finite(result$steiger_pval))
  expect_equal(result$snp_r2.exposure, sum(d$r.exposure^2), tolerance = 1e-15)
})

test_that("Steiger can approximate correlations from p-values and sample sizes", {
  d <- diagnostic_fixture()
  d$pval.exposure <- rep(1e-8, nrow(d))
  d$pval.outcome <- rep(1e-6, nrow(d))
  d$samplesize.exposure <- rep(10000, nrow(d))
  d$samplesize.outcome <- rep(12000, nrow(d))
  result <- fast_mr_directionality_test(d)
  expect_true(is.finite(result$snp_r2.exposure))
  expect_true(is.finite(result$snp_r2.outcome))
  expect_true(is.finite(result$steiger_pval))
})

test_that("Steiger reports missing input clearly", {
  d <- diagnostic_fixture()
  expect_null(fast_mr_directionality_test(d))
  expect_error(fast_mr_steiger(0.1, 0.2, 100, 100, r_xxo = 1.1), "r_xxo")
})

test_that("vectorized Steiger R-squared supports continuous beta, SE, and N", {
  result <- fast_mr_steiger_r2(
    beta = c(0.1, 0.2), se = 0.05, n = 1000,
    model = "continuous_bsen"
  )
  f_statistic <- (c(0.1, 0.2) / 0.05)^2
  expect_equal(result$rsq, f_statistic / (998 + f_statistic),
               tolerance = 1e-15)
  expect_equal(result$effective_n, rep(1000, 2))
  expect_true(all(result$valid))
  expect_equal(result$reason, rep("ok", 2))
})

test_that("vectorized Steiger R-squared supports standardized effects", {
  result <- fast_mr_steiger_r2(
    beta = c(0.1, -0.2), eaf = 0.25, model = "standardized"
  )
  expect_equal(result$rsq,
               2 * c(0.1, -0.2)^2 * 0.25 * 0.75,
               tolerance = 1e-15)
  expect_true(all(result$valid))
  expect_true(all(is.na(result$effective_n)))

  with_n <- fast_mr_steiger_r2(
    beta = 0.1, se = 1:3, eaf = c(0.25, 0.4), n = 10000,
    model = "standardized"
  )
  expect_equal(with_n$effective_n, rep(10000, 2))
})

test_that("vectorized Steiger R-squared requires binary prevalence and counts", {
  beta <- c(log(1.1), log(0.9))
  result <- fast_mr_steiger_r2(
    beta = beta, eaf = 0.3, prevalence = c(0.1, 0.2),
    ncase = 2000, ncontrol = 8000, model = "binary_lor"
  )
  expect_true(all(is.finite(result$rsq)))
  expect_true(all(result$rsq >= 0 & result$rsq <= 1))
  expect_equal(result$effective_n, rep(3200, 2))
  expect_true(all(result$valid))

  missing_prevalence <- fast_mr_steiger_r2(
    beta = beta, eaf = 0.3, ncase = 2000, ncontrol = 8000,
    model = "binary_lor"
  )
  expect_true(all(is.na(missing_prevalence$rsq)))
  expect_false(any(missing_prevalence$valid))
  expect_equal(missing_prevalence$reason,
               rep("missing_prevalence", length(beta)))
})

test_that("vectorized Steiger R-squared reports invalid rows deterministically", {
  result <- fast_mr_steiger_r2(
    beta = c(0.1, NA, 0.3, 0.4),
    se = c(0.05, 0.05, 0, 0.05),
    n = c(1000, 1000, 1000, 2),
    model = "continuous_bsen"
  )
  expect_equal(result$valid, c(TRUE, FALSE, FALSE, FALSE))
  expect_equal(result$reason,
               c("ok", "missing_beta", "invalid_se", "invalid_n"))
  expect_true(all(is.na(result$rsq[-1])))
  expect_error(
    fast_mr_steiger_r2(c(0.1, 0.2), se = c(0.01, 0.02, 0.03), n = 1000),
    "equal lengths or length one"
  )
})

test_that("Steiger filtering adds per-SNP quantitative-trait diagnostics", {
  d <- diagnostic_fixture()
  d$units.exposure <- "SD"
  d$units.outcome <- "SD"
  d$eaf.exposure <- 0.3
  d$eaf.outcome <- 0.4
  d$samplesize.exposure <- 10000
  d$samplesize.outcome <- 12000
  result <- fast_mr_steiger_filtering(d)
  expect_equal(nrow(result), nrow(d))
  expect_equal(result$rsq.exposure,
               2 * d$beta.exposure^2 * d$eaf.exposure * (1 - d$eaf.exposure),
               tolerance = 1e-15)
  expect_equal(result$rsq.outcome,
               2 * d$beta.outcome^2 * d$eaf.outcome * (1 - d$eaf.outcome),
               tolerance = 1e-15)
  expect_equal(result$effective_n.exposure, rep(10000, nrow(d)))
  expect_equal(result$effective_n.outcome, rep(12000, nrow(d)))
  # Native psych::r.test also returns NaN when an exaggerated summary-statistic
  # R-squared exceeds one; the useful finite rows must still be populated.
  expect_true(any(is.finite(result$steiger_pval)))
})

test_that("Steiger filtering supports log-odds metadata and supplied R-squared", {
  d <- diagnostic_fixture()[1:3, , drop = FALSE]
  d$units.exposure <- "log odds"
  d$units.outcome <- "log odds"
  d$eaf.exposure <- 0.3
  d$eaf.outcome <- 0.4
  d$ncase.exposure <- 2000
  d$ncontrol.exposure <- 8000
  d$ncase.outcome <- 3000
  d$ncontrol.outcome <- 9000
  d$prevalence.exposure <- 0.2
  d$prevalence.outcome <- 0.25
  result <- fast_mr_steiger_filtering(d)
  expect_equal(result$effective_n.exposure, rep(3200, 3))
  expect_equal(result$effective_n.outcome, rep(4500, 3))
  expect_true(all(is.finite(result$rsq.exposure)))
  expect_true(all(is.finite(result$rsq.outcome)))

  supplied <- diagnostic_fixture()[1:3, , drop = FALSE]
  supplied$rsq.exposure <- rep(0.02, 3)
  supplied$rsq.outcome <- rep(0.01, 3)
  supplied$effective_n.exposure <- rep(10000, 3)
  supplied$effective_n.outcome <- rep(12000, 3)
  supplied_result <- fast_mr_steiger_filtering(supplied)
  expect_equal(supplied_result$rsq.exposure, supplied$rsq.exposure)
  expect_true(all(supplied_result$steiger_dir))
})

test_that("Steiger filtering delegates to the explicit models", {
  continuous <- diagnostic_fixture()[1:2, , drop = FALSE]
  continuous$samplesize.exposure <- 10000
  continuous$samplesize.outcome <- 12000
  result <- fast_mr_steiger_filtering(continuous)
  expect_equal(
    result$rsq.exposure,
    fast_mr_steiger_r2(
      continuous$beta.exposure, continuous$se.exposure, 10000,
      model = "continuous_bsen"
    )$rsq
  )
  expect_true(all(result$rsq_valid.exposure))
  expect_true(all(result$steiger_pval_valid))

  standardized <- diagnostic_fixture()[1:2, , drop = FALSE]
  standardized$units.exposure <- "SD"
  standardized$units.outcome <- "SD"
  standardized$eaf.exposure <- 0.3
  standardized$eaf.outcome <- 0.4
  without_n <- fast_mr_steiger_filtering(standardized)
  expect_true(all(without_n$rsq_valid.exposure))
  expect_false(any(without_n$steiger_pval_valid))
  expect_equal(without_n$steiger_pval_reason,
               rep("missing_exposure_n", 2))
})

test_that("Steiger filtering never assumes binary prevalence", {
  binary <- diagnostic_fixture()[1:2, , drop = FALSE]
  binary$units.exposure <- "log odds"
  binary$units.outcome <- "log odds"
  binary$eaf.exposure <- 0.3
  binary$eaf.outcome <- 0.4
  binary$ncase.exposure <- 2000
  binary$ncontrol.exposure <- 8000
  binary$ncase.outcome <- 3000
  binary$ncontrol.outcome <- 9000

  result <- expect_no_warning(fast_mr_steiger_filtering(binary))
  expect_false("prevalence.exposure" %in% names(result))
  expect_true(all(is.na(result$rsq.exposure)))
  expect_equal(result$rsq_reason.exposure,
               rep("missing_prevalence", 2))
  expect_false(any(result$steiger_pval_valid))
})
