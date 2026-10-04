# Regression: a self-pair (outcome == exposure) is an exact fit, so the
# residual standard error sigma is exactly 0.  TwoSampleMR's IVW se is
# se_lm / min(1, sigma) (mr_ivw) or se_lm / sigma (mr_ivw_fe) with
# se_lm = base_se * sigma, whose limit is the fixed-effect se base_se; fastMR
# returned se = 0 and p = NA (ivw) or se = NA (ivw_fe).  The sparse IVW kernel
# gave base_se up to 0.1.9 and adopted the 0-se formula in d53afb8.

selfpair_data <- function(n = 12L, seed = 3L) {
  set.seed(seed)
  x <- runif(n, 0.2, 0.45)
  sy <- runif(n, 0.02, 0.04)
  data.frame(SNP = paste0("rs", seq_len(n)), id.exposure = "a", id.outcome = "a",
             exposure = "a", outcome = "a", beta.exposure = x, se.exposure = sy,
             beta.outcome = x, se.outcome = sy, mr_keep = TRUE)
}

test_that("IVW on a self-pair gives the fixed-effect se, as TwoSampleMR", {
  h <- selfpair_data()
  fe <- sqrt(1 / sum(h$beta.exposure^2 / h$se.outcome^2))
  got <- fast_mr(h, methods = c("ivw", "ivw_fe"))
  expect_equal(got$b, c(1, 1))
  expect_equal(got$se, c(fe, fe), tolerance = 1e-12)
  expect_true(all(is.finite(got$pval)))
  expect_equal(got$pval, rep(2 * stats::pnorm(-1 / fe), 2), tolerance = 1e-10)

  # Sparse CSR kernel (fast_mr_compressed's IVW-only path) agrees.
  n <- nrow(h)
  sp <- fast_mr_sparse_ivw(c(0L, n), seq_len(n) - 1L, h$beta.exposure,
                           matrix(h$beta.outcome, 1L), matrix(h$se.outcome, 1L),
                           matrix(TRUE, 1L, n))
  expect_identical(as.numeric(sp$beta), 1)
  expect_equal(as.numeric(sp$se), got$se[1L], tolerance = 1e-14)

  # The shared-grid batch path agrees too.
  grid <- fast_mr(rbind(h, transform(h, id.outcome = "b", outcome = "b",
                                     beta.outcome = beta.outcome * 0.5 + 0.01)),
                  methods = c("ivw", "ivw_fe"))
  self <- grid[grid$id.outcome == "a", ]
  expect_equal(self$se, c(fe, fe), tolerance = 1e-12)

  # Under-dispersed and over-dispersed fits are unchanged.
  h2 <- h
  set.seed(9)
  h2$beta.outcome <- h2$beta.exposure * 0.8 + stats::rnorm(n, 0, 0.002)
  h3 <- h
  h3$beta.outcome <- h3$beta.exposure * 0.8 + stats::rnorm(n, 0, 0.2)
  for (d in list(h2, h3)) {
    r <- fast_mr(d, methods = c("ivw", "ivw_fe"))
    m <- stats::lm(d$beta.outcome ~ -1 + d$beta.exposure, weights = 1 / d$se.outcome^2)
    s <- summary(m)
    expect_equal(r$se[1L], s$coefficients[1L, 2L] / min(1, s$sigma), tolerance = 1e-10)
    expect_equal(r$se[2L], s$coefficients[1L, 2L] / s$sigma, tolerance = 1e-10)
  }

  skip_if_not_installed("TwoSampleMR")
  tsmr <- suppressMessages(TwoSampleMR::mr(h, method_list = c("mr_ivw", "mr_ivw_fe")))
  expect_equal(got$se, tsmr$se, tolerance = 1e-8)
})

test_that("fast_mr_compressed IVW on a self-pair store gives the fixed-effect se", {
  skip_if_compressor_unavailable()
  skip_on_os("windows")
  set.seed(4)
  V <- 300L
  z <- rnorm(V)
  z[sample.int(V, 10L)] <- runif(10L, 7, 15)
  d <- data.frame(chromosome = "1", base_pair_location = seq.int(100001L, by = 7L, length.out = V),
                  reference_allele = "A", alternate_allele = "C", effect_allele = "C",
                  other_allele = "A", beta = z * 0.03, standard_error = 0.03,
                  effect_allele_frequency = 0.3)
  p <- tempfile("fm-self-")
  CompreSSoR::compress_sumstats(d, p, qc = "none", overwrite = TRUE)
  x <- CompreSSoR::read_sumstats(p, columns = c("chromosome", "base_pair_location", "other_allele",
                                                "effect_allele", "beta", "standard_error", "p_value"))
  x <- x[x$p_value <= 5e-8, ]
  keys <- CompreSSoR::compressor_variant_key(x$chromosome, x$base_pair_location,
                                             x$other_allele, x$effect_allele)
  fe <- sqrt(1 / sum(x$beta^2 / x$standard_error^2))
  r <- as.data.frame(fast_mr_compressed(c(a = p), c(a = p), list(a = keys), methods = "ivw"))
  expect_equal(r$b, 1)
  expect_equal(r$se, fe, tolerance = 1e-10)
  expect_true(is.finite(r$pval))
})
