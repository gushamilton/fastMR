# threads = 8 must give identical output whether the work-based cap keeps a
# small input on one worker or lets a large input use several.
thread_cap_data <- function(pairs, snps, seed = 7L) {
  set.seed(seed)
  g <- rep(seq_len(pairs), each = snps)
  n <- length(g)
  data.frame(
    SNP = paste0("rs", seq_len(n)),
    id.exposure = paste0("E", g), id.outcome = paste0("O", g),
    beta.exposure = rnorm(n, 0.2, 0.05), beta.outcome = rnorm(n, 0.1, 0.05),
    se.exposure = runif(n, 0.01, 0.05), se.outcome = runif(n, 0.01, 0.05),
    mr_keep = TRUE, stringsAsFactors = FALSE
  )
}

test_that("threads = 8 is identical to threads = 1 under the real work cap", {
  previous <- fastMR:::fastmr_set_work_scale_native(1)
  on.exit(fastMR:::fastmr_set_work_scale_native(previous), add = TRUE)
  for (shape in list(c(3L, 6L), c(600L, 20L))) {
    d <- thread_cap_data(shape[1], shape[2])
    expect_identical(fast_mr(d, nboot = 0, threads = 8),
                     fast_mr(d, nboot = 0, threads = 1))
    expect_identical(fast_mr(d, nboot = 20, seed = 3, threads = 8),
                     fast_mr(d, nboot = 20, seed = 3, threads = 1))
    expect_identical(fast_mr_heterogeneity(d, threads = 8),
                     fast_mr_heterogeneity(d, threads = 1))
    expect_identical(fast_mr_leaveoneout(d, method = "egger", threads = 8),
                     fast_mr_leaveoneout(d, method = "egger", threads = 1))
  }
})

test_that("the work scale hook validates its input and round-trips", {
  previous <- fastMR:::fastmr_set_work_scale_native(2)
  expect_identical(fastMR:::fastmr_set_work_scale_native(previous), 2)
  expect_error(fastMR:::fastmr_set_work_scale_native(-1))
  expect_error(fastMR:::fastmr_set_work_scale_native(NA_real_))
})
