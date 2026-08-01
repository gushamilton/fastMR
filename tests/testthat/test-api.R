test_that("fast_mr returns tidy TwoSampleMR-style results", {
  d <- il6_fixture()
  result <- fast_mr(d, nboot = 0, seed = 20260801)
  expect_s3_class(result, "data.frame")
  expect_equal(result$method_code, c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode"))
  expect_equal(result$nsnp, rep(82, 5))
  expect_named(result, c(
    "id.exposure", "id.outcome", "method", "method_code", "nsnp", "b", "se", "pval",
    "Q", "Q_df", "Q_pval", "sigma", "intercept", "intercept_se", "intercept_pval",
    "ratio_se_mean", "bootstrap", "phi", "flipped", "se_exposure_mean"
  ))
  expect_equal(result$b[result$method_code == "ivw"], 0.0012105753, tolerance = 1e-7)
  expect_equal(result$b[result$method_code == "egger"], 0.0034024061, tolerance = 1e-7)
})

test_that("dispatcher accepts common TwoSampleMR method spellings", {
  d <- il6_fixture()[1:4, ]
  result <- fast_mr(d, methods = c("mr_ivw", "mr_wald_ratio"), nboot = 0)
  expect_equal(result$method_code, c("ivw", "wald_ratio"))
  expect_equal(result$nsnp[[2]], 4)
  expect_true(is.na(result$b[[2]]))
  expect_true(all(c("code", "method", "description") %in% names(fastmr_method_registry())))
})

test_that("grouping preserves tidy pair metadata", {
  d <- il6_fixture()[1:8, ]
  d$id.exposure <- rep(c("A", "B"), each = 4)
  d$id.outcome <- rep(c("Y", "Z"), 4)
  result <- fast_mr(d, methods = "ivw", nboot = 0)
  expect_equal(nrow(result), 4)
  expect_equal(paste(result$id.exposure, result$id.outcome), c("A Y", "A Z", "B Y", "B Z"))
})
