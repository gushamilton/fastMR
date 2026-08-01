test_that("compact mixed grid preserves diagnostics and method ordering", {
  g <- grid_fixture(2L, 2L)
  result <- fast_mr_grid(
    g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se,
    methods = c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode"),
    nboot = 3L, seed = 20260911, threads = 5L
  )
  expect_equal(nrow(result), 20L)
  expect_equal(result$method_code[1:5],
               c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode"))
  expect_true(all(is.finite(result$Q[result$method_code == "ivw"])))
  expect_true(all(is.finite(result$intercept[result$method_code == "egger"])))
  expect_true(all(is.finite(result$bootstrap[result$method_code == "weighted_median"])))
})
