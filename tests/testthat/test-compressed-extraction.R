# Numeric-identity extraction and parallel store validation in
# fast_mr_compressed() must reproduce the string-key extraction exactly.

extraction_stores <- function(multipliers = c(1, 1.3, 0.7, -0.4, 0.2)) {
  vapply(multipliers, function(multiplier) {
    path <- tempfile("fastmr-extraction-")
    CompreSSoR::compress_sumstats(compressor_canonical_fixture(multiplier), path,
                                  overwrite = TRUE)
    path
  }, character(1))
}

extraction_keys <- function() {
  identity <- compressor_canonical_fixture()
  CompreSSoR::compressor_variant_key(
    identity$chromosome, identity$base_pair_location,
    identity$other_allele, identity$effect_allele
  )
}

test_that("manifest identity codes match the codes stores return", {
  skip_if_compressor_unavailable()
  path <- extraction_stores(1)
  store <- CompreSSoR::open_compressor(path)
  codec <- fastMR:::fastmr_compressed_identity_codec(store$manifest)
  expect_false(is.null(codec))
  keys <- extraction_keys()
  got <- tryCatch(CompreSSoR::read_sumstats_batch(
    path, keys,
    columns = c("global_position", "substitution", "chromosome",
                "base_pair_location", "effect_allele", "other_allele")
  )[[1L]], error = function(e) NULL)
  skip_if(is.null(got), "this CompreSSoR build does not return identity codes from key reads")
  row_keys <- CompreSSoR::compressor_variant_key(
    got$chromosome, got$base_pair_location, got$other_allele, got$effect_allele
  )
  expect_identical(
    fastMR:::fastmr_compressed_key_codes(row_keys, codec),
    as.numeric(got$global_position) * 16 + as.integer(got$substitution)
  )
  expect_identical(
    fastMR:::fastmr_compressed_decode_keys(
      as.numeric(got$global_position), as.integer(got$substitution), codec
    ),
    row_keys
  )
  # Keys that are not canonical single-nucleotide keys on the table get no code.
  expect_true(all(is.na(fastMR:::fastmr_compressed_key_codes(
    c("chr1:100:A:G", "1:0100:A:G", "1:100:a:g", "1:100:AT:G", "1:100:A:A",
      "1:999999999:A:G", "MT:100:A:G"), codec
  ))))
  broken <- store$manifest
  broken$identity$substitution_encoding <- "something_else"
  expect_null(fastMR:::fastmr_compressed_identity_codec(broken))
})

test_that("coded extraction is identical to string-key extraction", {
  skip_if_compressor_unavailable()
  # Either path (coded, or the string fallback on builds without identity
  # codes) must give the string-key result.
  stores <- extraction_stores()
  keys <- extraction_keys()
  absent <- c("2:200000000:A:C", "1:100:A:G")
  exposure_keys <- list(keys[c(2L, 7L, 14L, 25L, 40L, 61L)], c(keys[c(3L, 8L)], absent[[1L]]))
  union_keys <- unique(c(unlist(exposure_keys), absent))
  paths <- c(stores[1:2], stores[3:5], stores[3L])
  requests <- c(exposure_keys, rep(list(union_keys), 4L))
  codecs <- fastMR:::fastmr_compressed_validate_stores(unique(paths), 1L)
  for (columns in list(c("beta", "standard_error"),
                       c("beta", "standard_error", "effect_allele_frequency", "p_value"))) {
    for (io_threads in c(1L, 2L)) {
      strings <- fastMR:::fastmr_io_map(paths, requests, columns, io_threads)
      coded <- fastMR:::fastmr_io_map(paths, requests, columns, io_threads,
                                      codecs = unname(codecs[paths]))
      expect_identical(coded, strings, info = paste(length(columns), io_threads))
    }
  }
})

test_that("fast_mr_compressed results are unchanged by coded extraction", {
  skip_if_compressor_unavailable()
  stores <- extraction_stores()
  keys <- extraction_keys()
  exposures <- setNames(stores[1:2], c("exposure_a", "exposure_b"))
  outcomes <- setNames(stores[3:5], c("outcome_a", "outcome_b", "outcome_c"))
  instruments <- list(exposure_a = keys[c(2L, 7L, 14L, 25L, 40L, 61L)],
                      exposure_b = keys[c(3L, 8L, 19L)])
  run <- function(...) {
    out <- fast_mr_compressed(exposures, outcomes, instruments, ...)
    meta <- attr(out, "compressed_input")
    meta$timing <- NULL
    attr(out, "compressed_input") <- meta
    out
  }
  # Force the string path by withholding the codecs.
  string_run <- function(...) {
    testthat::with_mocked_bindings(
      run(...),
      fastmr_compressed_identity_codec = function(manifest) NULL,
      .package = "fastMR"
    )
  }
  args <- list(
    list(methods = "ivw", io_threads = 2),
    list(methods = c("wald_ratio", "egger", "weighted_median", "ivw", "weighted_mode"),
         nboot = 50, seed = 3, steiger = TRUE,
         samplesize_exposure = 5e4, samplesize_outcome = 1e5)
  )
  for (a in args) {
    expect_identical(do.call(run, a), do.call(string_run, a))
  }
})

test_that("parallel store validation raises the serial loop's first error", {
  skip_if_compressor_unavailable()
  good <- extraction_stores(c(1, 0.5))
  bad <- tempfile("fastmr-incompatible-")
  CompreSSoR::compress_sumstats(compressor_canonical_fixture(), bad, overwrite = TRUE)
  manifest_path <- file.path(bad, "manifest.json")
  manifest <- CompreSSoR:::read_manifest(manifest_path)
  manifest$identity$effect_allele_is_alt <- FALSE
  CompreSSoR:::write_manifest(manifest, manifest_path)
  CompreSSoR:::seal_pcodec_manifest(manifest_path)
  tampered <- extraction_stores(0.3)
  cat(" ", file = file.path(tampered, "manifest.json"), append = TRUE)
  paths <- normalizePath(c(good[[1L]], bad, tampered, good[[2L]]))
  serial_error <- function(paths) {
    tryCatch({
      invisible(lapply(paths, function(path) {
        fastMR:::fastmr_validate_compressed_store(CompreSSoR::open_compressor(path))
      }))
      NA_character_
    }, error = conditionMessage)
  }
  for (order in list(paths, rev(paths), paths[c(1L, 3L, 2L, 4L)])) {
    expected <- serial_error(order)
    expect_false(is.na(expected))
    for (io_threads in c(1L, 2L)) {
      expect_error(fastMR:::fastmr_compressed_validate_stores(order, io_threads),
                   expected, fixed = TRUE)
    }
  }
  codecs <- fastMR:::fastmr_compressed_validate_stores(normalizePath(good), 2L)
  expect_named(codecs, normalizePath(good))
  expect_false(any(vapply(codecs, is.null, logical(1))))
})

test_that("a build without identity codes falls back to the string path once", {
  skip_if_compressor_unavailable()
  stores <- extraction_stores(c(1, 0.5))
  keys <- extraction_keys()[c(2L, 7L, 14L)]
  codecs <- unname(fastMR:::fastmr_compressed_validate_stores(stores, 1L))
  strings <- fastMR:::fastmr_io_map(stores, list(keys, keys), c("beta", "standard_error"), 1L)
  state <- get(".fastmr_compressed_state", envir = asNamespace("fastMR"))
  previous <- state$coded_reads
  withr::defer(state$coded_reads <- previous)
  state$coded_reads <- NULL
  calls <- 0L
  real <- CompreSSoR::read_sumstats_batch
  refusing <- function(stores, variants = NULL, columns, threads = 1L, region = NULL) {
    calls <<- calls + 1L
    if (any(c("global_position", "substitution") %in% columns)) {
      stop("requested columns are not present: global_position, substitution")
    }
    real(stores, variants, columns = columns, threads = threads)
  }
  testthat::local_mocked_bindings(read_sumstats_batch = refusing, .package = "CompreSSoR")
  got <- fastMR:::fastmr_io_map(stores, list(keys, keys), c("beta", "standard_error"), 1L,
                                codecs = codecs)
  expect_identical(got, strings)
  expect_identical(calls, 2L)
  expect_false(fastMR:::fastmr_coded_reads_supported())
  got <- fastMR:::fastmr_io_map(stores, list(keys, keys), c("beta", "standard_error"), 1L,
                                codecs = codecs)
  expect_identical(got, strings)
  expect_identical(calls, 3L)
})

test_that("CompreSSoR's request index gives the string-key extraction", {
  skip_if_compressor_unavailable()
  skip_if_not(fastMR:::fastmr_compressor_has("request_index"),
              "this CompreSSoR build has no request_index")
  stores <- extraction_stores()
  keys <- extraction_keys()
  absent <- c("2:200000000:A:C", "1:100:A:G")
  exposure_keys <- list(keys[c(2L, 7L, 14L, 25L, 40L, 61L)], c(keys[c(3L, 8L)], absent[[1L]]))
  union_keys <- unique(c(unlist(exposure_keys), absent))
  paths <- c(stores[1:2], stores[3:5], stores[3L])
  requests <- c(exposure_keys, rep(list(union_keys), 4L))
  codecs <- unname(fastMR:::fastmr_compressed_validate_stores(unique(paths), 1L)[paths])
  expect_true(fastMR:::fastmr_request_index_usable(paths, requests, codecs))
  for (columns in list(c("beta", "standard_error"),
                       c("beta", "standard_error", "effect_allele_frequency", "p_value"))) {
    for (io_threads in c(1L, 2L)) {
      strings <- fastMR:::fastmr_io_map(paths, requests, columns, io_threads)
      coded <- fastMR:::fastmr_io_map(paths, requests, columns, io_threads, codecs = codecs,
                                      use_request_index = FALSE)
      indexed <- fastMR:::fastmr_io_map(paths, requests, columns, io_threads, codecs = codecs)
      expect_identical(indexed, strings, info = paste(length(columns), io_threads))
      expect_identical(coded, strings, info = paste(length(columns), io_threads))
    }
  }
  # Requests that are not strictly canonical keep the manifest-decoding path.
  expect_false(fastMR:::fastmr_request_index_usable(paths[1L], list(c(keys[1L], "chr1:5:A:G")),
                                                    codecs[1L]))
  expect_false(fastMR:::fastmr_request_index_usable(paths[1L], list(keys[1:2]), list(NULL)))
})

test_that("fast_mr_compressed results are unchanged by the request index", {
  skip_if_compressor_unavailable()
  skip_if_not(fastMR:::fastmr_compressor_has("request_index"),
              "this CompreSSoR build has no request_index")
  stores <- extraction_stores()
  keys <- extraction_keys()
  exposures <- setNames(stores[1:2], c("exposure_a", "exposure_b"))
  outcomes <- setNames(stores[3:5], c("outcome_a", "outcome_b", "outcome_c"))
  instruments <- list(exposure_a = keys[c(2L, 7L, 14L, 25L, 40L, 61L)],
                      exposure_b = keys[c(3L, 8L, 19L)])
  run <- function(...) {
    out <- fast_mr_compressed(exposures, outcomes, instruments, ...)
    meta <- attr(out, "compressed_input")
    meta$timing <- NULL
    attr(out, "compressed_input") <- meta
    out
  }
  without_index <- function(...) {
    testthat::with_mocked_bindings(
      run(...),
      fastmr_compressor_has = function(capability) FALSE,
      .package = "fastMR"
    )
  }
  args <- list(
    list(methods = "ivw", io_threads = 2),
    list(methods = c("wald_ratio", "egger", "weighted_median", "ivw", "weighted_mode"),
         nboot = 50, seed = 3, steiger = TRUE,
         samplesize_exposure = 5e4, samplesize_outcome = 1e5)
  )
  for (a in args) {
    expect_identical(do.call(run, a), do.call(without_index, a))
  }
})
