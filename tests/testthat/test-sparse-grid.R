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

sparse_steiger_fixture <- function(E = 6L, O = 4L, S = 12L, seed = 41L) {
  set.seed(seed)
  col_index <- unlist(lapply(seq_len(E), function(e) sort(sample.int(S, 5L) - 1L)))
  list(
    row_ptr = as.integer(c(0L, cumsum(rep(5L, E)))),
    col_index = as.integer(col_index),
    exposure_beta = rnorm(length(col_index), 0, 0.1),
    outcome_beta = matrix(rnorm(O * S, 0, 0.05), O, S),
    outcome_se = matrix(runif(O * S, 0.02, 0.05), O, S),
    outcome_present = matrix(runif(O * S) > 0.1, O, S),
    rsq_exp = runif(length(col_index), 0, 0.02),
    rsq_out = matrix(runif(O * S, 0, 0.02), O, S)
  )
}

test_that("steiger arguments equal the equivalent pair_snp_keep matrix", {
  x <- sparse_steiger_fixture()
  x$rsq_out[2, 3] <- NA
  x$rsq_exp[4] <- NA
  keep <- outer(seq_len(nrow(x$outcome_beta)), seq_along(x$col_index),
                Vectorize(function(o, k) isTRUE(x$rsq_exp[k] > x$rsq_out[o, x$col_index[k] + 1L])))
  base <- list(x$row_ptr, x$col_index, x$exposure_beta, x$outcome_beta,
               x$outcome_se, x$outcome_present)
  by_steiger <- do.call(fast_mr_sparse_ivw, c(base, list(
    steiger_exposure_rsq = x$rsq_exp, steiger_outcome_rsq = x$rsq_out)))
  by_keep <- do.call(fast_mr_sparse_ivw, c(base, list(pair_snp_keep = keep)))
  expect_identical(by_steiger[names(by_steiger) != "nsnp_prefilter"] |> `class<-`(class(by_steiger)),
                   by_keep)
  expect_true("nsnp_prefilter" %in% names(by_steiger))
  expect_true(all(by_steiger$nsnp <= by_steiger$nsnp_prefilter))
  plain <- do.call(fast_mr_sparse_ivw, base)
  expect_identical(names(plain), c("nsnp", "beta", "se", "Q", "sigma"))
  expect_identical(by_steiger$nsnp_prefilter, plain$nsnp)
  expect_null(do.call(fast_mr_sparse_ivw, c(base, list(pair_snp_keep = keep)))$nsnp_prefilter)

  # pair_snp_drop equals complement keep; combos AND together.
  drops <- which(!keep, arr.ind = TRUE)
  drop <- list(outcome = drops[, 1], entry = drops[, 2])
  set.seed(1)
  extra <- matrix(runif(length(keep)) > 0.2, nrow(keep))
  by_drop <- do.call(fast_mr_sparse_ivw, c(base, list(pair_snp_drop = drop)))
  expect_identical(by_drop, by_steiger)
  combo <- do.call(fast_mr_sparse_ivw, c(base, list(
    steiger_exposure_rsq = x$rsq_exp, steiger_outcome_rsq = x$rsq_out,
    pair_snp_keep = extra)))
  expect_identical(combo[names(combo) != "nsnp_prefilter"] |> `class<-`(class(combo)),
                   do.call(fast_mr_sparse_ivw, c(base, list(pair_snp_keep = keep & extra))))
  rand_drop <- which(!extra, arr.ind = TRUE)
  combo2 <- do.call(fast_mr_sparse_ivw, c(base, list(
    steiger_exposure_rsq = x$rsq_exp, steiger_outcome_rsq = x$rsq_out,
    pair_snp_drop = list(outcome = rand_drop[, 1], entry = rand_drop[, 2]))))
  expect_identical(combo2, combo)

  expect_error(do.call(fast_mr_sparse_ivw, c(base, list(steiger_exposure_rsq = x$rsq_exp))),
               "supplied together")
  expect_error(do.call(fast_mr_sparse_ivw, c(base, list(pair_snp_drop = list(outcome = 5L, entry = 1L)))),
               "1-based")
})

test_that("steiger sparse IVW matches fast_mr on Steiger-filtered data", {
  x <- sparse_steiger_fixture(E = 5L, O = 3L, S = 14L, seed = 7L)
  x$outcome_present[] <- TRUE
  res <- fast_mr_sparse_ivw(x$row_ptr, x$col_index, x$exposure_beta, x$outcome_beta,
    x$outcome_se, x$outcome_present,
    steiger_exposure_rsq = x$rsq_exp, steiger_outcome_rsq = x$rsq_out)
  for (e in seq_len(nrow(res$beta))) for (o in seq_len(ncol(res$beta))) {
    pos <- seq.int(x$row_ptr[e] + 1L, x$row_ptr[e + 1L])
    snp <- x$col_index[pos] + 1L
    keep <- x$rsq_exp[pos] > x$rsq_out[o, snp]
    expect_identical(unname(res$nsnp[e, o]), as.numeric(sum(keep)))
    if (sum(keep) < 2L) next
    d <- data.frame(
      SNP = paste0("s", snp[keep]), id.exposure = "e", id.outcome = "o",
      exposure = "e", outcome = "o",
      beta.exposure = x$exposure_beta[pos][keep], se.exposure = 0.01,
      beta.outcome = x$outcome_beta[o, snp[keep]], se.outcome = x$outcome_se[o, snp[keep]],
      mr_keep = TRUE)
    ref <- fast_mr(d, methods = "ivw")
    expect_equal(res$beta[e, o], ref$b[1], tolerance = 1e-12, ignore_attr = TRUE)
  }
})

test_that("steiger path does not allocate an outcome by entry matrix", {
  # O x N = 2000 x 2e5 logical (1.6 GB) would be prohibitive as a dense mask.
  E <- 2000L; O <- 2000L; S <- 5000L; per <- 100L
  set.seed(3)
  col_index <- unlist(lapply(seq_len(E), function(e) sort(sample.int(S, per) - 1L)))
  N <- length(col_index)
  rsq_exp <- runif(N, 0, 0.01)
  outcome_beta <- matrix(0.01, O, S)
  outcome_se <- matrix(0.02, O, S)
  outcome_present <- matrix(TRUE, O, S)
  rsq_out <- matrix(0.005, O, S)
  exposure_beta <- rnorm(N)
  row_ptr <- as.integer(c(0L, cumsum(rep(per, E))))
  col_index <- as.integer(col_index)
  gc(reset = TRUE)
  before <- sum(gc(reset = TRUE)[, 2L])
  res <- fast_mr_sparse_ivw(row_ptr, col_index, exposure_beta, outcome_beta,
    outcome_se, outcome_present,
    steiger_exposure_rsq = rsq_exp, steiger_outcome_rsq = rsq_out,
    max_output_cells = 1e7)
  g <- gc()
  extra_mb <- sum(g[, ncol(g)]) - before
  expect_equal(dim(res$nsnp), c(E, O))
  expect_true(all(res$nsnp <= res$nsnp_prefilter))
  expect_lt(extra_mb, 400)
})

test_that("fast_mr_steiger_rsq_matrix keeps shape and NA semantics", {
  b <- matrix(c(0.1, NA, 0.2, 0.3), 2L)
  r <- fast_mr_steiger_rsq_matrix(b, se = 0.05, n = 1e4)
  expect_identical(dim(r), dim(b))
  expect_true(is.na(r[2, 1]))
  expect_equal(c(r[1, 1], r[1, 2]),
               fast_mr_steiger_r2(c(0.1, 0.2), 0.05, 1e4)$rsq)
})
