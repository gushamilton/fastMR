# fastmr_request_groups() replaces pairwise identical() scans over the
# requested key sets; it must give the same grouping.

naive_request_groups <- function(requests) {
  vapply(seq_along(requests), function(i) {
    for (u in seq_len(i)) if (identical(requests[[u]], requests[[i]])) return(u)
    NA_integer_
  }, integer(1))
}

test_that("fastmr_request_groups matches a pairwise identical() scan", {
  shared <- c("1:100:A:C", "1:200:C:G", "1:300:G:T")
  requests <- list(
    shared, NULL, shared, c("1:100:A:C", "1:250:C:G", "1:300:G:T"),
    c("1:100:A:C", "1:200:C:G", "1:300:G:T"), character(), NULL, character(),
    c(0, 5, 7), c(-0, 5, 7), c(0L, 5L, 7L), c(NA, "1:1:A:C"), c(NA, "1:1:A:C"),
    c(NaN, 1), c(NA_real_, 1), c("NA", "1:1:A:C"), shared[1:2], rev(shared),
    c(a = "1:100:A:C"), "1:100:A:C", list(1, "a"), list(1, "a"), list("a", 1)
  )
  expect_identical(fastMR:::fastmr_request_groups(requests), naive_request_groups(requests))
  expect_identical(fastMR:::fastmr_request_groups(list()), integer())
  set.seed(7)
  pool <- sprintf("1:%d:A:C", 1:500)
  many <- c(lapply(1:400, function(i) sort(sample(pool, sample(1:6, 1)))),
            rep(list(pool), 50L))
  many <- many[sample.int(length(many))]
  expect_identical(fastMR:::fastmr_request_groups(many), naive_request_groups(many))
})

test_that("fastmr_request_index_usable checks every distinct request", {
  skip_if_not(fastMR:::fastmr_compressor_has("request_index"),
              "this CompreSSoR build has no request_index")
  set.seed(3)
  pool <- sprintf("%d:%d:%s:%s", sample(1:22, 300, TRUE), sample(1e3:1e7, 300),
                  "A", sample(c("C", "G", "T"), 300, TRUE))
  distinct <- lapply(1:200, function(i) sample(pool, sample(1:10, 1)))
  union <- unique(unlist(distinct))
  keys <- c(distinct, rep(list(union), 30L))
  paths <- rep("store", length(keys))
  codecs <- rep(list(list(codec = TRUE)), length(keys))
  expect_true(fastMR:::fastmr_request_index_usable(paths, keys, codecs))
  # One non-canonical key in one distinct request, or in a repeat of one.
  bad <- keys
  bad[[150L]] <- c(bad[[150L]], "chr1:5:A:G")
  expect_false(fastMR:::fastmr_request_index_usable(paths, bad, codecs))
  bad <- c(keys, list(c(union[1:3], "1:05:A:G")), list(c(union[1:3], "1:05:A:G")))
  expect_false(fastMR:::fastmr_request_index_usable(rep("store", length(bad)), bad,
                                                    rep(codecs[1L], length(bad))))
  # A non-character request after identical character ones.
  bad <- c(keys, list(1:3))
  expect_false(fastMR:::fastmr_request_index_usable(rep("store", length(bad)), bad,
                                                    rep(codecs[1L], length(bad))))
})

test_that("request index, coded and string extraction agree for many requests", {
  skip_if_compressor_unavailable()
  stores <- vapply(c(1, 1.3, 0.7), function(multiplier) {
    path <- tempfile("fastmr-request-groups-")
    CompreSSoR::compress_sumstats(compressor_canonical_fixture(multiplier), path,
                                  overwrite = TRUE)
    path
  }, character(1))
  identity <- compressor_canonical_fixture()
  keys <- CompreSSoR::compressor_variant_key(
    identity$chromosome, identity$base_pair_location,
    identity$other_allele, identity$effect_allele
  )
  set.seed(5)
  absent <- "2:200000000:A:C"
  exposure_keys <- lapply(1:60, function(i) c(sample(keys, sample(1:8, 1)), absent))
  union_keys <- unique(unlist(exposure_keys))
  requests <- c(exposure_keys, rep(list(union_keys), 9L), exposure_keys[1:4])
  paths <- rep_len(stores, length(requests))
  codecs <- unname(fastMR:::fastmr_compressed_validate_stores(unique(paths), 1L)[paths])
  columns <- c("beta", "standard_error")
  strings <- fastMR:::fastmr_io_map(paths, requests, columns, 2L)
  coded <- fastMR:::fastmr_io_map(paths, requests, columns, 2L, codecs = codecs,
                                  use_request_index = FALSE)
  expect_identical(coded, strings)
  if (fastMR:::fastmr_compressor_has("request_index")) {
    expect_true(fastMR:::fastmr_request_index_usable(paths, requests, codecs))
    indexed <- fastMR:::fastmr_io_map(paths, requests, columns, 2L, codecs = codecs)
    expect_identical(indexed, strings)
  }
})
