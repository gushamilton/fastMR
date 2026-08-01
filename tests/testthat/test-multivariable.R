test_that("basic multivariable IVW returns one row per exposure", {
  X <- cbind(x1 = c(.1, .2, .3, .4, .5), x2 = c(.2, .1, .4, .3, .6))
  y <- 0.5 * X[, 1] - 0.25 * X[, 2]
  result <- fast_mr_multivariable_ivw(X, y, rep(.1, 5))
  expect_equal(nrow(result), 2)
  expect_equal(result$method_code, rep("mv_ivw", 2))
  expect_equal(result$b, c(.5, -.25), tolerance = 1e-10)
})
