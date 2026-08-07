test_that("public masked API matches the internal masked kernel", {
  g <- grid_fixture(2L, 2L)
  exposure_present <- matrix(TRUE, nrow(g$exposure_beta), ncol(g$exposure_beta))
  outcome_present <- matrix(TRUE, nrow(g$outcome_beta), ncol(g$outcome_beta))
  observed <- fast_mr_masked_ivw(
    g$exposure_beta, g$outcome_beta, g$outcome_se,
    exposure_present, outcome_present
  )
  expected <- fastMR:::fastmr_masked_ivw_grid(
    g$exposure_beta, g$outcome_beta, g$outcome_se,
    exposure_present, outcome_present
  )
  for (name in c("nsnp", "beta", "se", "Q", "sigma")) {
    expect_equal(unname(observed[[name]]), unname(expected[[name]]))
  }
  expect_s3_class(observed, "fastmr_masked_ivw_compact")
})

test_that("public sparse API validates CSR and matches the internal kernel", {
  row_ptr <- as.integer(c(0L, 2L, 4L))
  col_index <- as.integer(c(0L, 2L, 1L, 3L))
  exposure_beta <- c(0.1, -0.2, 0.15, 0.08)
  outcome_beta <- matrix(c(0.04, 0.11, -0.03, 0.07,
                           -0.02, 0.08, 0.05, 0.01), 2L, 4L, byrow = TRUE)
  outcome_se <- matrix(0.03, 2L, 4L)
  outcome_present <- matrix(TRUE, 2L, 4L)
  observed <- fast_mr_sparse_ivw(
    row_ptr, col_index, exposure_beta, outcome_beta,
    outcome_se, outcome_present
  )
  expected <- fastMR:::fastmr_sparse_ivw_native(
    row_ptr, col_index, exposure_beta, outcome_beta,
    outcome_se, outcome_present
  )
  for (name in c("nsnp", "beta", "se", "Q", "sigma")) {
  expect_equal(unname(observed[[name]]), unname(expected[[name]]))
  }
  expect_s3_class(observed, "fastmr_sparse_ivw_compact")
  empty_row <- fast_mr_sparse_ivw(
    c(0L, 0L, 2L), c(0L, 2L), c(0.1, -0.2), outcome_beta,
    outcome_se, outcome_present
  )
  expect_equal(unname(empty_row$nsnp[1L, ]), c(0, 0))
  expect_error(fast_mr_sparse_ivw(
    c(0L, 2L), c(0L, 0L), c(0.1, 0.2), outcome_beta,
    outcome_se, outcome_present
  ), "duplicate")
  expect_error(fast_mr_sparse_ivw(
    c(0L, 1L, 4L), c(0L, 4L, 1L, 3L), exposure_beta,
    outcome_beta, outcome_se, outcome_present
  ), "zero-based")
})

test_that("public batch APIs reject invalid masks, SEs, dimensions, and bounds", {
  g <- grid_fixture(2L, 2L)
  exposure_present <- matrix(TRUE, nrow(g$exposure_beta), ncol(g$exposure_beta))
  outcome_present <- matrix(TRUE, nrow(g$outcome_beta), ncol(g$outcome_beta))
  expect_error(fast_mr_masked_ivw(
    g$exposure_beta, g$outcome_beta[, -1L, drop = FALSE], g$outcome_se,
    exposure_present, outcome_present
  ), "dimensions")
  bad_mask <- exposure_present
  bad_mask[1L, 1L] <- NA
  expect_error(fast_mr_masked_ivw(
    g$exposure_beta, g$outcome_beta, g$outcome_se,
    bad_mask, outcome_present
  ), "mask")
  bad_se <- g$outcome_se
  bad_se[1L, 1L] <- 0
  expect_error(fast_mr_masked_ivw(
    g$exposure_beta, g$outcome_beta, bad_se,
    exposure_present, outcome_present
  ), "positive standard errors")
  exposure_present[1L, 1L] <- FALSE
  g$exposure_beta[1L, 1L] <- NA_real_
  expect_s3_class(fast_mr_masked_ivw(
    g$exposure_beta, g$outcome_beta, g$outcome_se,
    exposure_present, outcome_present
  ), "fastmr_masked_ivw_compact")
  expect_error(fast_mr_masked_ivw(
    g$exposure_beta, g$outcome_beta, g$outcome_se,
    exposure_present, outcome_present, max_output_cells = 1
  ), "max_output_cells")
  expect_error(fast_mr_masked_ivw(
    g$exposure_beta, g$outcome_beta, g$outcome_se,
    exposure_present, outcome_present, max_memory_mb = 0.0001
  ), "max_memory_mb")

  outcome_present[1L, 1L] <- NA
  expect_error(fast_mr_sparse_ivw(
    c(0L, 1L), 0L, 0.1, g$outcome_beta, g$outcome_se, outcome_present
  ), "must not contain NA")
  expect_error(fast_mr_sparse_ivw(
    c(0L, 1L), 0L, 0.1, g$outcome_beta, g$outcome_se[, -1L, drop = FALSE],
    matrix(TRUE, nrow(g$outcome_beta), ncol(g$outcome_beta) - 1L)
  ), "dimensions")
})
