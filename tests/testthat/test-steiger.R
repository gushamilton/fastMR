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
