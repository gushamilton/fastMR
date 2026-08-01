test_that("invalid controls and method names fail clearly", {
  d <- il6_fixture()[1:4, ]
  expect_error(fast_mr(d, methods = "not_a_method"), "unknown MR method")
  expect_error(fast_mr(d, nboot = -1), "nboot")
  expect_error(fast_mr(d, threads = 0), "threads")
  expect_error(fast_mr(d, seed = 1.5), "seed")
  expect_error(fast_mr(d[, -1], nboot = 0), "SNP")
  bad <- d
  bad$beta.outcome[[1]] <- "not numeric"
  expect_error(fast_mr(bad, nboot = 0), "beta.outcome")
})

test_that("invalid grid shapes and values fail before native compute", {
  g <- grid_fixture(2, 2)
  expect_error(fast_mr_grid(g$exposure_beta, g$outcome_beta[, -1, drop = FALSE],
                             g$exposure_se, g$outcome_se, nboot = 0), "matching")
  g$outcome_se[[1]] <- 0
  expect_error(fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se,
                            nboot = 0), "positive standard errors")
  expect_error(fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se,
                            matrix("bad", nrow = 2, ncol = 82), nboot = 0), "numeric")
})

test_that("mr_keep false rows produce NA method results rather than errors", {
  d <- il6_fixture()[1:4, ]
  d$mr_keep <- FALSE
  result <- fast_mr(d, methods = "ivw", nboot = 0)
  expect_equal(result$nsnp, 0)
  expect_true(is.na(result$b))
})
