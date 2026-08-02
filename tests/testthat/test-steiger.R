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
