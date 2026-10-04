# fastmr_compressed_candidate_data() must not trust a batched p-value flag
# read that returns fewer (or more) rows than a store has flagged: CompreSSoR
# 0.7.0's read_candidates_batch() could silently drop rows in mixed batches,
# and an earlier fast path only fell back when the batch call errored.

guard_stores <- function() {
  mk <- function(seed, hits, offset) {
    set.seed(seed)
    V <- 800L
    z <- pmax(pmin(rnorm(V), 3), -3)
    if (hits) z[sample.int(V, hits)] <- sample(c(-1, 1), hits, TRUE) * runif(hits, 7, 20)
    data.frame(chromosome = "1", base_pair_location = seq.int(100001L + offset, by = 5L, length.out = V),
               reference_allele = "A", alternate_allele = "C", effect_allele = "C",
               other_allele = "A", beta = z * 0.04, standard_error = 0.04,
               effect_allele_frequency = 0.3)
  }
  paths <- c(a = tempfile("fm-guard-a-"), b = tempfile("fm-guard-b-"), c = tempfile("fm-guard-c-"))
  CompreSSoR::compress_sumstats(mk(1, 15, 0), paths[["a"]], overwrite = TRUE, pvalue_flag = TRUE)
  CompreSSoR::compress_sumstats(mk(2, 9, 50000), paths[["b"]], overwrite = TRUE, pvalue_flag = TRUE)
  CompreSSoR::compress_sumstats(mk(3, 0, 90000), paths[["c"]], overwrite = TRUE, pvalue_flag = TRUE)
  paths
}

test_that("a batched flag read that drops rows falls back to the per-store path", {
  skip_if_compressor_unavailable()
  skip_on_os("windows")
  skip_if_not(fastMR:::fastmr_have_flag_candidates_batch(), "batched pvalue_flag reader unavailable")
  paths <- guard_stores()
  labels <- names(paths)
  cand <- function() fastMR:::fastmr_compressed_candidate_data(paths, labels, 5e-8, "pvalue_flag",
                                                              "reconstructed", 1L)
  stores <- lapply(paths, CompreSSoR::open_compressor)
  expect_identical(unname(vapply(stores, fastMR:::fastmr_store_flag_count, numeric(1))), c(15, 9, 0))
  reference <- local({
    testthat::local_mocked_bindings(fastmr_have_flag_candidates_batch = function() FALSE,
                                    .package = "fastMR")
    cand()
  })
  expect_identical(nrow(reference$data), 24L)
  # The real batch is accepted silently and agrees with the per-store path.
  expect_no_warning(fast <- cand())
  expect_identical(fast$data, reference$data)

  real <- CompreSSoR::read_candidates_batch
  # A batch reader that silently drops one of store b's rows.
  testthat::local_mocked_bindings(
    fastmr_read_candidates_batch = function(...) {
      got <- real(...)
      got[[2L]] <- got[[2L]][-1L, , drop = FALSE]
      got
    }, .package = "fastMR")
  expect_warning(dropped <- cand(), "store 'b' returned 8 of 9 flagged rows")
  expect_identical(dropped$data, reference$data)
})

test_that("a batched flag read that errors or loses a store also falls back", {
  skip_if_compressor_unavailable()
  skip_on_os("windows")
  skip_if_not(fastMR:::fastmr_have_flag_candidates_batch(), "batched pvalue_flag reader unavailable")
  paths <- guard_stores()
  labels <- names(paths)
  cand <- function() fastMR:::fastmr_compressed_candidate_data(paths, labels, 5e-8, "pvalue_flag",
                                                              "reconstructed", 1L)
  reference <- local({
    testthat::local_mocked_bindings(fastmr_have_flag_candidates_batch = function() FALSE,
                                    .package = "fastMR")
    cand()
  })
  real <- CompreSSoR::read_candidates_batch
  local({
    testthat::local_mocked_bindings(fastmr_read_candidates_batch = function(...) stop("boom"),
                                    .package = "fastMR")
    expect_warning(x <- cand(), "failed: boom")
    expect_identical(x$data, reference$data)
  })
  local({
    testthat::local_mocked_bindings(fastmr_read_candidates_batch = function(...) real(...)[1:2],
                                    .package = "fastMR")
    expect_warning(x <- cand(), "returned 2 tables for 3 stores")
    expect_identical(x$data, reference$data)
  })
  local({
    # Extra rows are as wrong as missing ones.
    testthat::local_mocked_bindings(fastmr_read_candidates_batch = function(...) {
      got <- real(...)
      got[[1L]] <- rbind(got[[1L]], got[[1L]][1L, , drop = FALSE])
      got
    }, .package = "fastMR")
    expect_warning(x <- cand(), "store 'a' returned 16 of 15 flagged rows")
    expect_identical(x$data, reference$data)
  })
})
