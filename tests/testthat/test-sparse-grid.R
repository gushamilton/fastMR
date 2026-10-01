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

test_that("sparse pair masks preserve NULL and all-TRUE results exactly", {
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

  baseline <- fast_mr_sparse_ivw(
    row_ptr, col_index, exposure_beta, outcome_beta, outcome_se,
    outcome_present
  )
  explicit_null <- fast_mr_sparse_ivw(
    row_ptr, col_index, exposure_beta, outcome_beta, outcome_se,
    outcome_present, pair_snp_keep = NULL
  )
  all_true <- fast_mr_sparse_ivw(
    row_ptr, col_index, exposure_beta, outcome_beta, outcome_se,
    outcome_present,
    pair_snp_keep = matrix(TRUE, nrow(outcome_beta), length(col_index))
  )

  expect_identical(explicit_null, baseline)
  expect_identical(all_true, baseline)
})

test_that("sparse pair masks use outcome rows and concatenated CSR entry columns", {
  # Exposure two is empty. The remaining row boundaries deliberately do not
  # match SNP-column boundaries, so index leakage across CSR rows is visible.
  row_ptr <- as.integer(c(0L, 3L, 3L, 7L))
  col_index <- as.integer(c(5L, 0L, 3L, 1L, 5L, 2L, 4L))
  exposure_beta <- c(0.13, -0.08, 0.21, -0.15, 0.07, 0.18, -0.11)
  outcome_beta <- rbind(
    outcome_a = c(0.04, -0.01, 0.06, 0.09, -0.03, 0.08),
    outcome_b = c(-0.02, 0.05, 0.03, -0.04, 0.07, 0.01),
    outcome_c = c(0.06, 0.02, -0.05, 0.03, 0.04, -0.02)
  )
  outcome_se <- matrix(c(0.02, 0.03, 0.04, 0.025, 0.035, 0.03,
                         0.03, 0.025, 0.02, 0.04, 0.03, 0.035,
                         0.025, 0.04, 0.03, 0.02, 0.035, 0.03),
                       nrow = 3L, byrow = TRUE)
  outcome_present <- matrix(c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                              TRUE, TRUE, FALSE, TRUE, TRUE, TRUE,
                              TRUE, FALSE, TRUE, TRUE, TRUE, TRUE),
                            nrow = 3L, byrow = TRUE)
  pair_snp_keep <- matrix(c(
    TRUE, FALSE, TRUE,  TRUE, TRUE, FALSE, TRUE,
    FALSE, TRUE, TRUE,  TRUE, FALSE, TRUE, TRUE,
    TRUE, TRUE, FALSE,  FALSE, TRUE, TRUE, TRUE
  ), nrow = 3L, byrow = TRUE)

  observed <- fast_mr_sparse_ivw(
    row_ptr, col_index, exposure_beta, outcome_beta, outcome_se,
    outcome_present, threads = 2L, pair_snp_keep = pair_snp_keep
  )
  expect_equal(unname(observed$nsnp[2L, ]), c(0, 0, 0))

  for (i in seq_len(length(row_ptr) - 1L)) {
    positions <- if (row_ptr[i] < row_ptr[i + 1L]) {
      seq.int(row_ptr[i] + 1L, row_ptr[i + 1L])
    } else integer()
    snps <- col_index[positions] + 1L
    for (j in seq_len(nrow(outcome_beta))) {
      keep <- pair_snp_keep[j, positions] & outcome_present[j, snps]
      if (sum(keep) < 2L) {
        expect_equal(unname(observed$nsnp[i, j]), sum(keep))
        expect_true(is.na(observed$beta[i, j]))
      } else {
        expected <- manual_ivw(
          exposure_beta[positions][keep], outcome_beta[j, snps][keep],
          outcome_se[j, snps][keep]
        )
        expect_equal(unname(observed$nsnp[i, j]), expected$nsnp)
        expect_equal(unname(observed$beta[i, j]), expected$beta, tolerance = 1e-12)
        expect_equal(unname(observed$se[i, j]), expected$se, tolerance = 1e-12)
        expect_equal(unname(observed$Q[i, j]), expected$Q, tolerance = 1e-12)
        expect_equal(unname(observed$sigma[i, j]), expected$sigma, tolerance = 1e-12)
      }
    }
  }

  # Processing an outcome block uses only the corresponding mask rows.
  outcome_block <- c(3L, 1L)
  blocked <- fast_mr_sparse_ivw(
    row_ptr, col_index, exposure_beta,
    outcome_beta[outcome_block, , drop = FALSE],
    outcome_se[outcome_block, , drop = FALSE],
    outcome_present[outcome_block, , drop = FALSE],
    pair_snp_keep = pair_snp_keep[outcome_block, , drop = FALSE]
  )
  for (name in c("nsnp", "beta", "se", "Q", "sigma")) {
    expect_identical(unname(blocked[[name]]),
                     unname(observed[[name]][, outcome_block, drop = FALSE]))
  }
})

test_that("sparse pair-specific masks agree with tidy fast_mr IVW", {
  row_ptr <- as.integer(c(0L, 4L, 8L))
  col_index <- as.integer(c(0L, 2L, 3L, 5L, 1L, 2L, 4L, 5L))
  exposure_beta <- c(0.13, -0.08, 0.21, 0.09,
                     -0.15, 0.07, 0.18, -0.11)
  outcome_beta <- rbind(
    outcome_a = c(0.04, -0.01, 0.06, 0.09, -0.03, 0.08),
    outcome_b = c(-0.02, 0.05, 0.03, -0.04, 0.07, 0.01)
  )
  outcome_se <- matrix(c(0.02, 0.03, 0.04, 0.025, 0.035, 0.03,
                         0.03, 0.025, 0.02, 0.04, 0.03, 0.035),
                       nrow = 2L, byrow = TRUE)
  outcome_present <- matrix(TRUE, 2L, 6L)
  pair_snp_keep <- matrix(c(
    TRUE, FALSE, TRUE, TRUE,  TRUE, TRUE, FALSE, TRUE,
    TRUE, TRUE, FALSE, TRUE,  FALSE, TRUE, TRUE, TRUE
  ), nrow = 2L, byrow = TRUE)

  observed <- fast_mr_sparse_ivw(
    row_ptr, col_index, exposure_beta, outcome_beta, outcome_se,
    outcome_present, pair_snp_keep = pair_snp_keep
  )
  for (i in 1:2) {
    positions <- seq.int(row_ptr[i] + 1L, row_ptr[i + 1L])
    snps <- col_index[positions] + 1L
    for (j in 1:2) {
      keep <- pair_snp_keep[j, positions]
      tidy <- fast_mr(data.frame(
        SNP = paste0("rs", snps[keep]),
        beta.exposure = exposure_beta[positions][keep],
        beta.outcome = outcome_beta[j, snps][keep],
        se.exposure = rep(0.01, sum(keep)),
        se.outcome = outcome_se[j, snps][keep]
      ), methods = "ivw", nboot = 0L)
      expect_equal(unname(observed$nsnp[i, j]), tidy$nsnp)
      expect_equal(unname(observed$beta[i, j]), tidy$b, tolerance = 1e-12)
      expect_equal(unname(observed$se[i, j]), tidy$se, tolerance = 1e-12)
      expect_equal(unname(observed$Q[i, j]), tidy$Q, tolerance = 1e-12)
    }
  }
})

test_that("sparse pair masks reject ambiguous or missing entries", {
  row_ptr <- as.integer(c(0L, 2L))
  col_index <- as.integer(c(0L, 2L))
  exposure_beta <- c(0.1, -0.2)
  outcome_beta <- matrix(c(0.04, 0.11, -0.03,
                           -0.02, 0.08, 0.05), 2L, 3L, byrow = TRUE)
  outcome_se <- matrix(0.03, 2L, 3L)
  outcome_present <- matrix(TRUE, 2L, 3L)

  expect_error(fast_mr_sparse_ivw(
    row_ptr, col_index, exposure_beta, outcome_beta, outcome_se,
    outcome_present, pair_snp_keep = matrix(1, 2L, 2L)
  ), "logical matrix")
  expect_error(fast_mr_sparse_ivw(
    row_ptr, col_index, exposure_beta, outcome_beta, outcome_se,
    outcome_present, pair_snp_keep = matrix(TRUE, 2L, 3L)
  ), "one row per outcome")
  expect_error(fast_mr_sparse_ivw(
    row_ptr, col_index, exposure_beta, outcome_beta, outcome_se,
    outcome_present,
    pair_snp_keep = matrix(c(TRUE, NA, TRUE, TRUE), 2L, 2L)
  ), "must not contain NA")
  expect_error(fastMR:::fastmr_sparse_ivw_native(
    row_ptr, col_index, exposure_beta, outcome_beta, outcome_se,
    outcome_present, pair_snp_keep = matrix(TRUE, 2L, 3L)
  ), "one row per outcome")
  expect_error(fastMR:::fastmr_sparse_ivw_native(
    row_ptr, col_index, exposure_beta, outcome_beta, outcome_se,
    outcome_present,
    pair_snp_keep = matrix(c(TRUE, NA, TRUE, TRUE), 2L, 2L)
  ), "must not contain NA")
})
