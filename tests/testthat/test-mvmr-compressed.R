mvmr_store_world <- function() {
  identity <- compressor_canonical_fixture()
  keys <- CompreSSoR::compressor_variant_key(
    identity$chromosome, identity$base_pair_location,
    identity$other_allele, identity$effect_allele
  )
  n <- length(keys)
  set.seed(42)
  make <- function(beta, se) {
    d <- identity
    d$beta <- beta
    d$standard_error <- se
    path <- tempfile("fastmr-mvmr-store-")
    CompreSSoR::compress_sumstats(d, path, overwrite = TRUE)
    path
  }
  # Stores may round on encoding: references use the values read back.
  stored <- function(path) {
    got <- fast_read_compressed(path, keys, columns = c("beta", "standard_error"))
    got <- got[match(keys, got$variant_key), ]
    list(beta = got$beta, se = got$standard_error)
  }
  bx <- matrix(rnorm(n * 2, 0, 0.1), n, 2)
  bc <- rnorm(n, 0, 0.1)
  se_x <- matrix(runif(n * 2, 0.005, 0.01), n, 2)
  se_c <- runif(n, 0.005, 0.01)
  by <- sapply(1:3, function(k) 0.2 * k * bx[, 1] - 0.1 * bx[, 2] + 0.4 * bc +
                 rnorm(n, 0, 0.01))
  se_y <- matrix(runif(n * 3, 0.008, 0.02), n, 3)
  exposures <- c(E1 = make(bx[, 1], se_x[, 1]), E2 = make(bx[, 2], se_x[, 2]))
  covariate <- c(PLT = make(bc, se_c))
  outcomes <- c(O1 = make(by[, 1], se_y[, 1]), O2 = make(by[, 2], se_y[, 2]),
                O3 = make(by[, 3], se_y[, 3]))
  ex <- lapply(exposures, stored)
  cv <- stored(covariate)
  out <- lapply(outcomes, stored)
  list(
    keys = keys,
    bx = sapply(ex, `[[`, "beta"), se_x = sapply(ex, `[[`, "se"),
    bc = cv$beta, se_c = cv$se,
    by = sapply(out, `[[`, "beta"), se_y = sapply(out, `[[`, "se"),
    exposures = exposures, covariate = covariate, outcomes = outcomes
  )
}

test_that("fast_mvmr_compressed (mvmr) equals the in-memory kernel on the same rows", {
  skip_if_compressor_unavailable()
  w <- mvmr_store_world()
  inst <- list(E1 = w$keys[1:20], E2 = w$keys[21:40])
  cinst <- w$keys[41:60]
  R <- matrix(c(1, 0.2, 0.2, 1), 2, dimnames = list(c("exposure", "PLT"), c("exposure", "PLT")))
  res <- fast_mvmr_compressed(w$exposures, w$outcomes, w$covariate, inst,
                              covariate_instruments = cinst, method = "mvmr",
                              exposure_cor = R, covariate_estimates = TRUE,
                              threads = 2, io_threads = 2, weak_f = -Inf)
  expect_equal(nrow(res), 6L)
  expect_identical(res$id.exposure, rep(c("E1", "E2"), each = 3))
  expect_identical(res$id.outcome, rep(c("O1", "O2", "O3"), 2))
  for (e in 1:2) {
    rows <- c(if (e == 1) 1:20 else 21:40, 41:60)
    X <- cbind(exposure = w$bx[rows, e], PLT = w$bc[rows])
    ref <- fast_mvmr_ivw(X, w$by[rows, ], w$se_y[rows, ],
                         exposure_se = cbind(w$se_x[rows, e], w$se_c[rows]),
                         exposure_cor = R, weak_f = -Inf)
    got <- res[res$id.exposure == paste0("E", e), ]
    expect_equal(got$b, unname(ref$b["exposure", ]), tolerance = 1e-13)
    expect_equal(got$se, unname(ref$se["exposure", ]), tolerance = 1e-13)
    expect_equal(got$b_PLT, unname(ref$b["PLT", ]), tolerance = 1e-13)
    expect_equal(got$Q_A, unname(ref$Q_A), tolerance = 1e-13)
    expect_equal(got$conditional_F, rep(ref$conditional_F[["exposure"]], 3))
    expect_true(all(got$nsnp == 40))
  }
  meta <- attr(res, "mvmr_input")
  expect_setequal(names(meta$timing), c("io_seconds", "design_seconds",
                                        "estimator_seconds", "total_seconds",
                                        "source_bytes_read"))
  expect_equal(meta$diagnostics$nsnp_design, c(40L, 40L))

  # Matrix output carries the same numbers.
  mat <- fast_mvmr_compressed(w$exposures, w$outcomes, w$covariate, inst,
                              covariate_instruments = cinst, exposure_cor = R,
                              output_format = "matrix", weak_f = -Inf)
  expect_equal(as.vector(t(mat$b)), res$b)
  expect_equal(as.vector(t(mat$se)), res$se)

  # The covariate as an external data frame, with some keys allele-swapped.
  f <- fastMR:::fastmr_mvmr_key_fields(w$keys)
  swap <- seq(1, 80, by = 3)
  df_keys <- w$keys
  df_keys[swap] <- paste(f[swap, 1], f[swap, 2], f[swap, 4], f[swap, 3], sep = ":")
  df_beta <- w$bc
  df_beta[swap] <- -df_beta[swap]
  frame <- data.frame(SNP = df_keys, beta = df_beta, se = w$se_c)
  from_frame <- fast_mvmr_compressed(w$exposures, w$outcomes, list(PLT = frame), inst,
                                     covariate_instruments = cinst, exposure_cor = R,
                                     weak_f = -Inf)
  expect_equal(from_frame$b, res$b, tolerance = 1e-14)
})

test_that("fast_mvmr_compressed (residualised) propagates the covariate uncertainty", {
  skip_if_compressor_unavailable()
  w <- mvmr_store_world()
  inst <- list(E1 = w$keys[1:20], E2 = w$keys[21:40])
  cinst <- w$keys[41:60]
  res <- fast_mvmr_compressed(w$exposures, w$outcomes, w$covariate, inst,
                              covariate_instruments = cinst, method = "residualised",
                              weak_f = -Inf)
  expect_identical(unique(res$method_code), "residualised_ivw")
  gamma <- fast_mvmr_ivw(w$bc[41:60], w$by[41:60, ], w$se_y[41:60, ],
                         se_model = "multiplicative_floored")
  fit <- attr(res, "mvmr_input")$covariate_fit
  expect_equal(unname(fit$b[1, ]), unname(gamma$b[1, ]), tolerance = 1e-13)
  expect_equal(unname(fit$se[1, ]), unname(gamma$se[1, ]), tolerance = 1e-13)
  for (e in 1:2) {
    rows <- if (e == 1) 1:20 else 21:40
    adj <- w$by[rows, ] - outer(w$bc[rows], gamma$b[1, ])
    sd <- sqrt(w$se_y[rows, ]^2 + outer(w$bc[rows]^2, gamma$se[1, ]^2))
    ref <- fast_mvmr_ivw(w$bx[rows, e], adj, sd, se_model = "multiplicative_floored")
    got <- res[res$id.exposure == paste0("E", e), ]
    expect_equal(got$b, unname(ref$b[1, ]), tolerance = 1e-13)
    expect_equal(got$se, unname(ref$se[1, ]), tolerance = 1e-13)
  }
})

test_that("outcome_exclude, shared exposure/outcome stores and LD pairs", {
  skip_if_compressor_unavailable()
  w <- mvmr_store_world()
  inst <- list(E1 = w$keys[1:20])
  cinst <- w$keys[41:60]
  pos <- 100001L + 0:79
  exclude <- data.frame(id.outcome = "O2", chromosome = "1",
                        start = pos[5], end = pos[45])
  res <- fast_mvmr_compressed(w$exposures[1], w$outcomes, w$covariate, inst,
                              covariate_instruments = cinst, outcome_exclude = exclude,
                              weak_f = -Inf, exposure_cor = diag(2))
  expect_equal(res$nsnp, c(40, 40 - 21, 40))
  rows <- c(1:20, 41:60)
  B <- w$by[rows, ]
  B[rows >= 5 & rows <= 45, 2] <- NA
  ref <- fast_mvmr_ivw(cbind(w$bx[rows, 1], w$bc[rows]), B, w$se_y[rows, ])
  expect_equal(res$b, unname(ref$b[1, ]), tolerance = 1e-13)

  # An exposure that is also an outcome is read once, before masking.
  both <- c(w$outcomes, E1 = unname(w$exposures[1]))
  shared <- suppressWarnings(fast_mvmr_compressed(c(E1 = unname(w$exposures[1])), both, w$covariate, inst,
                                 covariate_instruments = cinst,
                                 outcome_exclude = data.frame(id.outcome = "E1",
                                                              chromosome = "1",
                                                              start = pos[1], end = pos[80]),
                                 weak_f = -Inf, strict = FALSE, exposure_cor = diag(2)))
  expect_equal(shared$b[shared$id.outcome == "O1"], res$b[1], tolerance = 1e-13)
  expect_false("E1" %in% shared$id.outcome)  # fully masked: 0 SNPs, omitted

  # LD between instrument 1 and covariate instrument 41: the weaker is dropped.
  ld <- data.frame(a = w$keys[1], b = w$keys[41])
  pruned <- fast_mvmr_compressed(w$exposures[1], w$outcomes, w$covariate, inst,
                                 covariate_instruments = cinst, ld_pairs = ld,
                                 weak_f = -Inf, exposure_cor = diag(2))
  kept <- attr(pruned, "mvmr_input")$instruments$E1
  expect_length(kept, 39L)
  z1 <- max(abs(w$bx[1, 1] / w$se_x[1, 1]), abs(w$bc[1] / w$se_c[1]))
  z41 <- max(abs(w$bx[41, 1] / w$se_x[41, 1]), abs(w$bc[41] / w$se_c[41]))
  expect_identical(setdiff(w$keys[rows], kept), w$keys[if (z1 >= z41) 41 else 1])
})

test_that("strict mode, missing instruments and weak-instrument warnings", {
  skip_if_compressor_unavailable()
  w <- mvmr_store_world()
  inst <- list(E1 = c(w$keys[1:20], "1:999999:A:C"))
  cinst <- w$keys[41:60]
  expect_error(
    fast_mvmr_compressed(w$exposures[1], w$outcomes, w$covariate, inst,
                         covariate_instruments = cinst, weak_f = -Inf,
                         exposure_cor = diag(2)),
    "missing or invalid"
  )
  warned <- capture_warnings(
    res <- fast_mvmr_compressed(w$exposures[1], w$outcomes, w$covariate, inst,
                                covariate_instruments = cinst, strict = FALSE,
                                weak_f = -Inf, exposure_cor = diag(2))
  )
  expect_true(any(grepl("missing or invalid outcome values", warned)))
  expect_true(any(grepl("dropped for missing or invalid exposure values", warned)))
  expect_true(all(res$nsnp == 40))
  expect_error(
    fast_mvmr_compressed(w$exposures[1], w$outcomes, w$covariate, list(E1 = w$keys[1:20]),
                         method = "residualised"),
    "covariate_instruments"
  )
  # Covariate instruments only weakly separate the exposures here, so ask for
  # an F threshold that is certainly not met.
  expect_warning(
    fast_mvmr_compressed(w$exposures[1], w$outcomes, w$covariate, list(E1 = w$keys[1:20]),
                         covariate_instruments = cinst, weak_f = 1e6,
                         exposure_cor = diag(2)),
    "conditional F statistic below"
  )
  expect_warning(
    fast_mvmr_compressed(w$exposures[1], w$outcomes, w$covariate, list(E1 = w$keys[1:20]),
                         covariate_instruments = cinst, weak_f = -Inf),
    "exposure_cor not supplied"
  )
})
