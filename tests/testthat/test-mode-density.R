# The mode density has two evaluation paths: direct convolution of the binned
# weights (small pairs) and the FFT (large pairs, and any direct draw whose top
# two grid densities are too close to be sure of the argmax). Every result must
# be identical whichever path runs.

mode_pairs <- function(sizes, seed = 5L) {
  set.seed(seed)
  g <- rep(seq_along(sizes), sizes)
  n <- length(g)
  bx <- sample(c(-1, 1), n, TRUE) * runif(n, 0.03, 0.12)
  sy <- runif(n, 0.008, 0.02)
  data.frame(
    SNP = paste0("rs", seq_len(n)), id.exposure = paste0("E", g), id.outcome = "O",
    beta.exposure = bx, se.exposure = runif(n, 0.004, 0.01),
    beta.outcome = rnorm(length(sizes), 0, 0.3)[g] * bx + rnorm(n, 0, sy), se.outcome = sy,
    stringsAsFactors = FALSE
  )
}

# Ratios mirrored around 0.37 with equal SEs: the two density peaks tie in
# exact arithmetic, and with tiny SEs every bootstrap draw stays a near-tie.
mode_ties <- function(k, se) {
  a <- seq(0.1, 1, length.out = k %/% 2)
  r <- c(-a, a) + 0.37
  data.frame(SNP = paste0("rs", seq_along(r)), id.exposure = "E", id.outcome = "O",
             beta.exposure = 0.1, se.exposure = se, beta.outcome = 0.1 * r,
             se.outcome = se, stringsAsFactors = FALSE)
}

with_direct_max <- function(value, f) {
  previous <- fastMR:::fastmr_set_mode_direct_max_native(value)
  on.exit(fastMR:::fastmr_set_mode_direct_max_native(previous), add = TRUE)
  f()
}

modes <- c("simple_mode", "weighted_mode")

test_that("direct and FFT mode densities give identical results", {
  d <- mode_pairs(c(3L, 4L, 5L, 10L, 25L, 60L))
  for (m in list(modes, "simple_mode", "weighted_mode")) {
    run <- function() fast_mr(d, methods = m, nboot = 60, seed = 2, phi = 0.7)
    fft <- with_direct_max(0, run)
    expect_identical(with_direct_max(Inf, run), fft)
    expect_identical(run(), fft)
  }
})

test_that("near-ties take the guard and still match the FFT path", {
  fastMR:::fastmr_mode_path_counts_native(TRUE)
  for (se in c(1e-12, 1e-9, 1e-6)) {
    d <- mode_ties(10L, se)
    run <- function() fast_mr(d, methods = modes, nboot = 50, seed = 4)
    expect_identical(with_direct_max(Inf, run), with_direct_max(0, run))
  }
  counts <- fastMR:::fastmr_mode_path_counts_native(TRUE)
  expect_gt(counts[["direct"]], 0)
  expect_gt(counts[["guard"]], 0)
  expect_lte(counts[["guard"]], counts[["direct"]])
})

test_that("the direct-path hook validates and round-trips", {
  previous <- fastMR:::fastmr_set_mode_direct_max_native(7)
  expect_identical(fastMR:::fastmr_set_mode_direct_max_native(previous), 7)
  expect_error(fastMR:::fastmr_set_mode_direct_max_native(-1))
  expect_error(fastMR:::fastmr_set_mode_direct_max_native(NA_real_))
})

test_that("threaded bootstrap normals match for every normal.kind", {
  old <- RNGkind()
  on.exit(RNGkind(old[1], old[2], old[3]), add = TRUE)
  d <- mode_pairs(c(4L, 6L, 9L, 3L, 12L))
  for (kind in c("Inversion", "Box-Muller", "Kinderman-Ramage", "Ahrens-Dieter")) {
    RNGkind(normal.kind = kind)
    for (seed in list(NULL, 11)) {
      set.seed(1)
      ref <- fast_mr(d, nboot = 30, seed = seed, threads = 1)
      ref_state <- .Random.seed
      for (threads in c(2L, 4L)) {
        set.seed(1)
        expect_identical(fast_mr(d, nboot = 30, seed = seed, threads = threads), ref)
        expect_identical(.Random.seed, ref_state)
      }
    }
  }
})
