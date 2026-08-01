test_that("grid ordering is exposure-major and outcome-minor", {
  g <- grid_fixture(2, 3)
  rownames(g$exposure_beta) <- c("E1", "E2")
  rownames(g$outcome_beta) <- c("O1", "O2", "O3")
  result <- fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se,
                          methods = "ivw", nboot = 0, threads = 1)
  expect_equal(paste(result$exposure_index, result$outcome_index),
               rep(c("1 1", "1 2", "1 3", "2 1", "2 2", "2 3"), each = 1))
  expect_equal(result$id.exposure, rep(c("E1", "E1", "E1", "E2", "E2", "E2"), each = 1))
  expect_equal(result$id.outcome, rep(c("O1", "O2", "O3", "O1", "O2", "O3"), each = 1))
})

test_that("seeded grid results are deterministic across thread counts", {
  g <- grid_fixture(3, 3)
  a <- fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se,
                    nboot = 7, seed = 20260801, threads = 1)
  b <- fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se,
                    nboot = 7, seed = 20260801, threads = 4)
  d <- fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se,
                    nboot = 7, seed = 20260801, threads = 10)
  expect_equal(a, b, tolerance = 0)
  expect_equal(a, d, tolerance = 0)
  c <- fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se,
                    nboot = 7, seed = 20260802, threads = 1)
  expect_true(any(abs(a$se - c$se) > 0, na.rm = TRUE))
})

test_that("grid returns every method for every pair", {
  g <- grid_fixture(2, 2)
  result <- fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se,
                          methods = c("ivw", "egger", "weighted_median"), nboot = 3, seed = 1)
  expect_equal(nrow(result), 12)
  expect_equal(as.vector(table(result$method_code)), c(4, 4, 4))
})

test_that("simple median grid results are thread deterministic", {
  g <- grid_fixture(2, 2)
  a <- fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se,
                    methods = c("simple_median", "weighted_median"),
                    nboot = 5, seed = 20260804, threads = 1)
  b <- fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se,
                    methods = c("simple_median", "weighted_median"),
                    nboot = 5, seed = 20260804, threads = 5)
  expect_equal(a, b, tolerance = 0)
})

test_that("MR-Egger bootstrap grid results are thread deterministic", {
  g <- grid_fixture(2, 2)
  a <- fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se,
                    methods = "egger_bootstrap", nboot = 5,
                    seed = 20260805, threads = 1)
  b <- fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se,
                    methods = "egger_bootstrap", nboot = 5,
                    seed = 20260805, threads = 5)
  expect_equal(a, b, tolerance = 0)
})

test_that("unweighted regression and sign grid methods cover every pair", {
  g <- grid_fixture(2, 2)
  result <- fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se,
                          methods = c("uwr", "sign"), nboot = 0, threads = 5)
  expect_equal(nrow(result), 8)
  expect_equal(as.vector(table(result$method_code)), c(4, 4))
  expect_true(all(is.finite(result$b)))
})
