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

# Review finding 6: the count-only guard accepted a batch with the right number
# of rows but the wrong rows or keys, and the per-store fallback reader did not
# check its own row count.

test_that("a batched flag read with the right count but wrong rows or keys is rejected", {
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
  expect_true(".fastmr_abs_z" %in% names(reference$data))
  expect_true(all(is.finite(reference$data$.fastmr_abs_z) & reference$data$.fastmr_abs_z > 5))
  real <- CompreSSoR::read_candidates_batch
  flagged_b <- fastMR:::fastmr_store_flag_rows(CompreSSoR::open_compressor(paths[["b"]]))
  local({
    # same count, but one row id is not a flagged row of store b. fastMR checks
    # row ids itself only on a CompreSSoR build without
    # "candidates_batch_rows_checked" (a build with it compares every store's
    # decoded rows with its own flag selection and stops on a mismatch), so
    # this case is exercised with the capability switched off.
    testthat::local_mocked_bindings(fastmr_have_checked_candidates_batch = function() FALSE,
                                    .package = "fastMR")
    testthat::local_mocked_bindings(fastmr_read_candidates_batch = function(...) {
      got <- real(...)
      got[[2L]]$row[1L] <- setdiff(0:799, flagged_b)[1L]
      got
    }, .package = "fastMR")
    expect_warning(x <- cand(), "store 'b' returned row ids that are not its flagged rows")
    expect_identical(x$data, reference$data)
  })
  local({
    # right rows, but a key decoded against another panel (position disagrees)
    testthat::local_mocked_bindings(fastmr_read_candidates_batch = function(...) {
      got <- real(...)
      got[[2L]]$key[1L] <- got[[1L]]$key[1L]
      got
    }, .package = "fastMR")
    expect_warning(x <- cand(), "store 'b' returned 1 key\\(s\\) whose position does not match")
    expect_identical(x$data, reference$data)
  })
})

test_that("the per-store flagged-row reader stops on a short read", {
  skip_if_compressor_unavailable()
  skip_on_os("windows")
  paths <- guard_stores()
  real <- CompreSSoR::read_sumstats
  testthat::local_mocked_bindings(fastmr_have_flag_candidates_batch = function() FALSE, .package = "fastMR")
  testthat::local_mocked_bindings(
    read_sumstats = function(path, variants = NULL, ...) {
      x <- real(path, variants = variants, ...)
      if (!is.null(variants) && nrow(x) > 1L) x[-1L, , drop = FALSE] else x
    }, .package = "CompreSSoR")
  expect_error(fastMR:::fastmr_compressed_candidate_data(paths, names(paths), 5e-8, "pvalue_flag",
                                                        "reconstructed", 1L),
               "flagged-row read returned 14 of 15 rows")
})

test_that("with a rows-checked CompreSSoR the flag stream is not decoded a second time", {
  skip_if_compressor_unavailable()
  skip_on_os("windows")
  skip_if_not(fastMR:::fastmr_have_flag_candidates_batch(), "batched pvalue_flag reader unavailable")
  skip_if_not(fastMR:::fastmr_have_checked_candidates_batch(), "CompreSSoR without candidates_batch_rows_checked")
  paths <- guard_stores()
  labels <- names(paths)
  cand <- function() fastMR:::fastmr_compressed_candidate_data(paths, labels, 5e-8, "pvalue_flag",
                                                              "reconstructed", 1L)
  reference <- local({
    testthat::local_mocked_bindings(fastmr_have_flag_candidates_batch = function() FALSE,
                                    .package = "fastMR")
    cand()
  })
  unchecked <- local({
    testthat::local_mocked_bindings(fastmr_have_checked_candidates_batch = function() FALSE,
                                    .package = "fastMR")
    cand()
  })
  testthat::local_mocked_bindings(fastmr_store_flag_rows = function(...) stop("flag rows pre-read"),
                                  .package = "fastMR")
  expect_no_warning(fast <- cand())
  expect_identical(fast$data, reference$data)
  expect_identical(fast$data, unchecked$data)
  # the count check (manifest flagged-row count) still guards the batch
  real <- CompreSSoR::read_candidates_batch
  testthat::local_mocked_bindings(
    fastmr_read_candidates_batch = function(...) {
      got <- real(...)
      got[[2L]] <- got[[2L]][-1L, , drop = FALSE]
      got
    }, .package = "fastMR")
  expect_warning(dropped <- cand(), "store 'b' returned 8 of 9 flagged rows")
  expect_identical(dropped$data, reference$data)
})
