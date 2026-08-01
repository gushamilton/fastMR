test_that("penalised weighted median is registered and finite", {
  d <- il6_fixture()[1:12, ]
  result <- fast_mr(d, methods = c("mr_penalised_weighted_median"),
                    nboot = 25, seed = 20260901)
  expect_equal(result$method, "Penalised weighted median")
  expect_true(all(is.finite(result$b)))
  expect_true(all(is.finite(result$se)))
  expect_true(all(is.finite(result$pval)))
  expect_true("penalised_weighted_median" %in% fastmr_method_registry()$code)
})

test_that("penalised weighted median grid is deterministic across threads", {
  g <- grid_fixture(3L, 4L)
  serial <- fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se,
                         g$outcome_se, methods = "penalised_weighted_median",
                         nboot = 25, seed = 20260902, threads = 1)
  parallel <- fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se,
                           g$outcome_se, methods = "penalised_weighted_median",
                           nboot = 25, seed = 20260902, threads = 5)
  expect_equal(serial$b, parallel$b, tolerance = 0)
  expect_equal(serial$se, parallel$se, tolerance = 0)
  expect_equal(serial$pval, parallel$pval, tolerance = 0)
})

test_that("penalty multiplier changes penalised weighting", {
  d <- il6_fixture()[1:12, ]
  default <- fast_mr(d, methods = "penalised_weighted_median",
                     nboot = 0, seed = 20260903)
  strong <- fast_mr(d, methods = "penalised_weighted_median",
                    nboot = 0, seed = 20260903, penk = 2)
  expect_true(is.finite(default$b))
  expect_true(is.finite(strong$b))
  expect_false(isTRUE(all.equal(default$b, strong$b)))
})
