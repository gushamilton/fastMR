manual_ivw <- function(x, y, se) {
  w <- 1 / (se * se)
  denominator <- sum(x * x * w)
  numerator <- sum(x * y * w)
  beta <- numerator / denominator
  q <- sum(w * (y - beta * x)^2)
  sigma <- sqrt(q / (length(x) - 1L))
  list(
    nsnp = length(x),
    beta = beta,
    se = sqrt(1 / denominator) * max(1, sigma),
    Q = q,
    sigma = sigma
  )
}

test_that("masked IVW preserves pair-specific instrument membership", {
  exposure_beta <- rbind(
    exposure_a = c(0.10, 0.20, -0.15, 0.08),
    exposure_b = c(0.05, -0.12, 0.09, 0.14)
  )
  outcome_beta <- rbind(
    outcome_a = c(0.04, 0.11, -0.03, 0.07),
    outcome_b = c(-0.02, 0.08, 0.05, 0.01)
  )
  outcome_se <- matrix(c(0.02, 0.03, 0.025, 0.04,
                         0.03, 0.02, 0.04, 0.025), nrow = 2L, byrow = TRUE)
  exposure_present <- matrix(c(TRUE, TRUE, FALSE, TRUE,
                               TRUE, FALSE, TRUE, TRUE), nrow = 2L, byrow = TRUE)
  outcome_present <- matrix(c(TRUE, FALSE, TRUE, TRUE,
                              TRUE, TRUE, FALSE, TRUE), nrow = 2L, byrow = TRUE)

  observed <- fastMR:::fastmr_masked_ivw_grid(
    exposure_beta, outcome_beta, outcome_se,
    exposure_present, outcome_present
  )
  for (i in seq_len(nrow(exposure_beta))) {
    for (j in seq_len(nrow(outcome_beta))) {
      keep <- exposure_present[i, ] & outcome_present[j, ]
      expected <- manual_ivw(exposure_beta[i, keep], outcome_beta[j, keep], outcome_se[j, keep])
      expect_equal(observed$nsnp[i, j], expected$nsnp)
      expect_equal(observed$beta[i, j], expected$beta, tolerance = 1e-12)
      expect_equal(observed$se[i, j], expected$se, tolerance = 1e-12)
      expect_equal(observed$Q[i, j], expected$Q, tolerance = 1e-12)
      expect_equal(observed$sigma[i, j], expected$sigma, tolerance = 1e-12)
    }
  }
  expect_error(
    fastMR:::fastmr_masked_ivw_grid(
      exposure_beta, outcome_beta, outcome_se,
      {x <- exposure_present; x[1, 1] <- NA; x}, outcome_present
    ),
    "must not contain NA"
  )
})

test_that("sparse IVW agrees with manual IVW on a small CSR panel", {
  row_ptr <- as.integer(c(0L, 3L, 5L))
  col_index <- as.integer(c(0L, 2L, 4L, 1L, 3L))
  exposure_beta <- c(0.10, -0.15, 0.08, 0.20, -0.12)
  outcome_beta <- rbind(
    outcome_a = c(0.04, 0.11, -0.03, 0.07, 0.02),
    outcome_b = c(-0.02, 0.08, 0.05, 0.01, -0.04)
  )
  outcome_se <- matrix(c(0.02, 0.03, 0.025, 0.04, 0.03,
                         0.03, 0.02, 0.04, 0.025, 0.05), nrow = 2L, byrow = TRUE)
  outcome_present <- matrix(c(TRUE, TRUE, TRUE, TRUE, FALSE,
                              TRUE, TRUE, TRUE, TRUE, TRUE), nrow = 2L, byrow = TRUE)

  observed <- fastMR:::fastmr_sparse_ivw_native(
    row_ptr, col_index, exposure_beta, outcome_beta, outcome_se,
    outcome_present, threads = 1L
  )
  for (i in seq_len(length(row_ptr) - 1L)) {
    positions <- seq.int(row_ptr[i] + 1L, row_ptr[i + 1L])
    snps <- col_index[positions] + 1L
    for (j in seq_len(nrow(outcome_beta))) {
      keep <- outcome_present[j, snps]
      expected <- manual_ivw(
        exposure_beta[positions][keep], outcome_beta[j, snps][keep], outcome_se[j, snps][keep]
      )
      expect_equal(observed$nsnp[i, j], expected$nsnp)
      expect_equal(observed$beta[i, j], expected$beta, tolerance = 1e-12)
      expect_equal(observed$se[i, j], expected$se, tolerance = 1e-12)
      expect_equal(observed$Q[i, j], expected$Q, tolerance = 1e-12)
      expect_equal(observed$sigma[i, j], expected$sigma, tolerance = 1e-12)
    }
  }
  expect_error(
    fastMR:::fastmr_sparse_ivw_native(
      as.integer(c(0L, 2L)), as.integer(c(0L, 0L)), c(0.1, 0.2),
      outcome_beta, outcome_se, outcome_present, threads = 1L
    ),
    "duplicate SNP indices"
  )
})
