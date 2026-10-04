# One instrument (nsnp = 1).  TwoSampleMR::mr() then reports only
# mr_wald_ratio: b = by / bx, se = se_y / |bx| (first order), and every other
# method is NA.  fastMR's IVW estimators return exactly that Wald ratio on
# every path (fast_mr, the grid, the masked and sparse kernels, compressed
# input); the sparse kernel used to return no estimate.  Q and sigma are
# undefined (0 degrees of freedom).  An exact fit with k >= 2 is a different
# case and keeps the fixed-effect se (test-ivw-selfpair.R).

single_wald <- function(bx, by, sy) {
  b <- by / bx
  se <- sy / abs(bx)
  list(b = b, se = se, pval = stats::pnorm(abs(b) / se, lower.tail = FALSE) * 2)
}

single_data <- function() {
  data.frame(SNP = c("rs1", "rs2", "rs3", "rs4"), id.exposure = c("a", "b", "b", "c"),
             id.outcome = "o", beta.exposure = c(0.2, 0.13, -0.3, -0.07),
             se.exposure = 0.01, beta.outcome = c(0.05, 0.02, 0.01, 0.031),
             se.outcome = c(0.02, 0.015, 0.02, 0.011), stringsAsFactors = FALSE)
}

all_methods <- c("ivw", "ivw_fe", "ivw_mre", "egger", "egger_bootstrap", "uwr", "sign",
                 "simple_median", "weighted_median", "penalised_weighted_median",
                 "simple_mode", "weighted_mode", "wald_ratio")

test_that("fast_mr: IVW at nsnp = 1 is TwoSampleMR's Wald ratio; other methods are NA", {
  d <- single_data()
  res <- fast_mr(d, methods = all_methods, nboot = 20, seed = 1)
  for (id in c("a", "c")) {
    row <- d[d$id.exposure == id, ]
    w <- single_wald(row$beta.exposure, row$beta.outcome, row$se.outcome)
    r <- res[res$id.exposure == id, ]
    ivw <- r[r$method_code %in% c("ivw", "ivw_fe", "ivw_mre", "wald_ratio"), ]
    expect_identical(ivw$b, rep(w$b, 4))
    expect_identical(ivw$se, rep(w$se, 4))
    expect_identical(ivw$pval, rep(w$pval, 4))
    expect_identical(ivw$nsnp, rep(1, 4))
    expect_true(all(is.na(r$Q) & !is.nan(r$Q)))
    expect_true(all(is.na(r$sigma)))
    others <- r[!r$method_code %in% c("ivw", "ivw_fe", "ivw_mre", "wald_ratio"), ]
    for (col in c("b", "se", "pval")) {
      expect_true(all(is.na(others[[col]])), info = col)
      expect_false(any(is.nan(others[[col]])), info = col)
    }
  }
  # k = 2: IVW is the ordinary estimator and the Wald ratio is NA, as TwoSampleMR
  b <- res[res$id.exposure == "b", ]
  expect_true(is.finite(b$b[b$method_code == "ivw"]))
  expect_true(is.na(b$b[b$method_code == "wald_ratio"]))
  # a zero exposure effect has no Wald ratio
  z <- transform(d[1, ], beta.exposure = 0)
  expect_true(is.na(fast_mr(z, methods = "ivw")$b))
})

test_that("nsnp = 1 matches TwoSampleMR::mr_wald_ratio exactly", {
  skip_if_not_installed("TwoSampleMR")
  d <- single_data()[1, ]
  t <- TwoSampleMR::mr_wald_ratio(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome)
  f <- fast_mr(d, methods = c("ivw", "ivw_fe", "ivw_mre"))
  expect_identical(f$b, rep(t$b, 3))
  expect_identical(f$se, rep(t$se, 3))
  expect_identical(f$pval, rep(t$pval, 3))
  for (m in c("mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_ivw")) {
    expect_true(is.na(get(m, asNamespace("TwoSampleMR"))(d$beta.exposure, d$beta.outcome, d$se.exposure,
                                                         d$se.outcome, TwoSampleMR::default_parameters())$b))
  }
})

test_that("fast_mr_grid with one SNP column returns the Wald ratio (tidy, compact, mixed methods)", {
  eb <- matrix(c(0.2, -0.07, 0.4), 3, 1); es <- matrix(0.01, 3, 1)
  ob <- matrix(c(0.05, -0.02), 2, 1); os <- matrix(c(0.02, 0.011), 2, 1)
  expected <- function(e, o) single_wald(eb[e, 1], ob[o, 1], os[o, 1])
  check <- function(res) {
    for (e in 1:3) for (o in 1:2) {
      r <- res[res$id.exposure == as.character(e) & res$id.outcome == as.character(o) &
                 res$method_code %in% c("ivw", "ivw_fe"), ]
      w <- expected(e, o)
      expect_identical(r$b, rep(w$b, nrow(r)))
      expect_identical(r$se, rep(w$se, nrow(r)))
      expect_identical(r$pval, rep(w$pval, nrow(r)))
      expect_true(all(is.na(r$Q)))
    }
  }
  check(fast_mr_grid(eb, ob, es, os, methods = c("ivw", "ivw_fe")))
  check(as.data.frame(fast_mr_grid(eb, ob, es, os, methods = c("ivw", "ivw_fe"), return = "compact")))
  mixed <- fast_mr_grid(eb, ob, es, os, methods = c("ivw", "ivw_fe", "egger", "weighted_median"),
                        nboot = 10, seed = 1)
  check(mixed)
  expect_true(all(is.na(mixed$b[mixed$method_code %in% c("egger", "weighted_median")])))
  # identical to fast_mr on the same pair
  one <- fast_mr(data.frame(SNP = "s", beta.exposure = eb[3, 1], se.exposure = 0.01,
                            beta.outcome = ob[2, 1], se.outcome = os[2, 1]), methods = "ivw")
  g <- fast_mr_grid(eb, ob, es, os, methods = "ivw")
  expect_identical(g$b[g$id.exposure == "3" & g$id.outcome == "2"], one$b)
})

test_that("the sparse and masked IVW kernels return the Wald ratio for one shared instrument", {
  # exposure 1: SNPs 1-3; exposure 2: SNP 2 only; exposure 3: SNPs 1 and 4,
  # of which outcome 2 has only SNP 4.
  ob <- rbind(c(0.05, -0.02, 0.03, 0.01), c(0.04, 0.02, -0.01, 0.06))
  os <- rbind(c(0.02, 0.01, 0.015, 0.012), c(0.03, 0.02, 0.011, 0.014))
  present <- rbind(c(TRUE, TRUE, TRUE, TRUE), c(FALSE, TRUE, TRUE, TRUE))
  row_ptr <- c(0L, 3L, 4L, 6L)
  col_index <- c(0L, 1L, 2L, 1L, 0L, 3L)
  ebeta <- c(0.2, 0.15, -0.1, 0.3, 0.25, -0.4)
  sp <- fast_mr_sparse_ivw(row_ptr, col_index, ebeta, ob, os, present)
  w21 <- single_wald(0.3, ob[1, 2], os[1, 2]); w22 <- single_wald(0.3, ob[2, 2], os[2, 2])
  w32 <- single_wald(-0.4, ob[2, 4], os[2, 4])
  expect_identical(unname(sp$nsnp[2, ]), c(1, 1))
  expect_identical(unname(sp$beta[2, ]), c(w21$b, w22$b))
  expect_identical(unname(sp$se[2, ]), c(w21$se, w22$se))
  expect_identical(unname(sp$nsnp[3, 2]), 1)
  expect_identical(unname(sp$beta[3, 2]), w32$b)
  expect_identical(unname(sp$se[3, 2]), w32$se)
  expect_true(all(is.na(sp$Q[2, ])) && all(is.na(sp$sigma[2, ])))
  expect_false(any(is.nan(sp$beta)))
  expect_true(all(is.finite(sp$beta[1, ])))
  # masked dense form of the same problem
  eb <- matrix(0, 3, 4); ep <- matrix(FALSE, 3, 4)
  for (e in 1:3) {
    idx <- (row_ptr[e] + 1L):row_ptr[e + 1L]
    eb[e, col_index[idx] + 1L] <- ebeta[idx]
    ep[e, col_index[idx] + 1L] <- TRUE
  }
  mk <- fast_mr_masked_ivw(eb, ob, os, ep, present)
  expect_identical(unname(mk$nsnp), unname(sp$nsnp))
  expect_identical(unname(mk$beta[2, ]), unname(sp$beta[2, ]))
  expect_identical(unname(mk$se[2, ]), unname(sp$se[2, ]))
  expect_identical(unname(mk$beta[3, 2]), w32$b)
  expect_identical(unname(mk$se[3, 2]), w32$se)
})

test_that("compressed input: single-instrument exposures get the Wald ratio on every path", {
  skip_if_compressor_unavailable()
  skip_on_os("windows")
  stores <- vapply(c(1, 1.3, 0.7), function(m) {
    p <- tempfile("fm-single-")
    CompreSSoR::compress_sumstats(compressor_canonical_fixture(m), p, overwrite = TRUE)
    p
  }, character(1))
  id <- compressor_canonical_fixture()
  keys <- CompreSSoR::compressor_variant_key(id$chromosome, id$base_pair_location, id$other_allele,
                                             id$effect_allele)
  ex <- stats::setNames(stores[1:2], c("ea", "eb"))
  out <- c(oa = stores[[3]])
  sets <- list(ea = keys[1:5], eb = keys[7])
  auto <- fast_mr_compressed(ex, out, sets, methods = "ivw")
  expect_identical(attr(auto, "compressed_input")$estimator_path, "sparse_ivw")
  pair <- fast_mr_compressed(ex, out, sets, methods = "ivw", estimator = "pairwise")
  expect_identical(attr(pair, "compressed_input")$estimator_path, "pairwise")
  e <- fast_read_compressed(stores[[2]], variants = keys[7], columns = c("beta", "standard_error"))
  o <- fast_read_compressed(stores[[3]], variants = keys[7], columns = c("beta", "standard_error"))
  w <- single_wald(e$beta, o$beta, o$standard_error)
  for (res in list(auto, pair)) {
    r <- res[res$id.exposure == "eb", ]
    expect_identical(r$nsnp, 1)
    expect_identical(r$b, w$b)
    expect_identical(r$se, w$se)
    expect_equal(r$pval, w$pval, tolerance = 1e-14)
    expect_true(is.na(r$Q) && is.na(r$Q_df) && is.na(r$sigma))
  }
  strip <- function(x) { attr(x, "compressed_input") <- NULL; x }
  expect_equal(strip(auto), strip(pair), tolerance = 1e-14)
  # a single shared instrument takes the shared-grid path, with the same answer
  grid <- fast_mr_compressed(ex, out, keys[7], methods = "ivw")
  expect_identical(attr(grid, "compressed_input")$estimator_path, "shared_instrument_grid")
  expect_identical(grid$b[grid$id.exposure == "eb"], w$b)
})
