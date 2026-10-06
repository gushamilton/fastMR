test_that("basic multivariable IVW returns one row per exposure", {
  X <- cbind(x1 = c(.1, .2, .3, .4, .5), x2 = c(.2, .1, .4, .3, .6))
  y <- 0.5 * X[, 1] - 0.25 * X[, 2]
  result <- fast_mr_multivariable_ivw(X, y, rep(.1, 5))
  expect_equal(nrow(result), 2)
  expect_equal(result$method_code, rep("mv_ivw", 2))
  expect_equal(result$b, c(.5, -.25), tolerance = 1e-10)
})

test_that("multivariable path supports shared and exposure-specific instruments", {
  X <- cbind(x1 = c(.1, .2, .3, .4, .5, .6),
             x2 = c(.2, .1, .4, .3, .6, .5))
  y <- 0.5 * X[, 1] - 0.25 * X[, 2] + c(.01, -.01, .02, -.02, .01, -.01)
  sy <- rep(.1, nrow(X))
  P <- matrix(c(1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9,
                1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-6), ncol = 2)
  shared <- fast_mr_multivariable(X, y, sy, P, pval_threshold = 5e-8)
  specific <- fast_mr_multivariable(X, y, sy, P, pval_threshold = 5e-8,
                                    instrument_specific = TRUE)
  intercept <- fast_mr_multivariable(X, y, sy, P, intercept = TRUE)
  expect_equal(shared$nsnp, c(6, 5))
  expect_equal(specific$nsnp, c(6, 5))
  expect_equal(intercept$nsnp, c(6, 5))
  expect_true(all(is.finite(shared$b)))
  expect_true(all(is.finite(specific$b)))
  expect_true(all(is.finite(intercept$b)))
  expect_equal(shared$method_code, rep("mv_multiple", 2))
})

test_that("collinear exposures warn and blank column names get labels", {
  set.seed(11)
  x1 <- rnorm(30)
  X <- cbind(a = x1, x1 * 2)
  y <- 0.3 * x1 + rnorm(30, sd = 0.05)
  sy <- rep(0.05, 30)
  expect_warning(ivw <- fast_mr_multivariable_ivw(X, y, sy), "singular")
  expect_equal(ivw$id.exposure, c("a", "2"))
  expect_true(all(is.na(ivw$b)))
  expect_warning(mv <- fast_mr_multivariable(X, y, sy), "singular")
  expect_equal(mv$id.exposure, c("a", "2"))
  ok <- cbind(a = x1, b = rnorm(30))
  expect_no_warning(fast_mr_multivariable_ivw(ok, y, sy))
  expect_no_warning(fast_mr_multivariable(ok, y, sy))
})
