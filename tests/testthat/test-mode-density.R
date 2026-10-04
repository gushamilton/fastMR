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

with_hull <- function(ratios, recurrence, f) {
  previous <- fastMR:::fastmr_set_mode_hull_native(ratios, recurrence)
  on.exit(fastMR:::fastmr_set_mode_hull_native(previous[["ratios"]], previous[["recurrence"]] == 1),
          add = TRUE)
  f()
}

test_that("the hull path, with and without the recurrence kernel, matches the FFT", {
  sizes <- c(2L, 3L, 4L, 5L, 7L, 10L, 16L, 17L, 33L, 64L, 65L, 100L)
  clean <- mode_pairs(sizes, seed = 8L)
  outliers <- clean
  hit <- seq(1L, nrow(outliers), by = 3L)
  outliers$beta.outcome[hit] <- outliers$beta.outcome[hit] + 0.3
  for (d in list(clean, outliers, mode_ties(12L, 1e-9), mode_ties(7L, 1e-6))) {
    run <- function() fast_mr(d, methods = modes, nboot = 40, seed = 6, phi = 0.8)
    fft <- with_direct_max(0, run)
    fastMR:::fastmr_mode_path_counts_native(TRUE)
    for (recurrence in c(TRUE, FALSE)) {
      expect_identical(with_hull(Inf, recurrence, run), fft)
      expect_identical(with_hull(0, recurrence, run), fft)
    }
    expect_identical(run(), fft)
    expect_gt(fastMR:::fastmr_mode_path_counts_native(TRUE)[["hull"]], 0)
  }
})

test_that("the hull hook validates and round-trips", {
  previous <- fastMR:::fastmr_set_mode_hull_native(9, FALSE)
  restored <- fastMR:::fastmr_set_mode_hull_native(previous[["ratios"]], previous[["recurrence"]] == 1)
  expect_identical(unname(restored), c(9, 0))
  expect_error(fastMR:::fastmr_set_mode_hull_native(-1))
  expect_error(fastMR:::fastmr_set_mode_hull_native(NA_real_))
})
