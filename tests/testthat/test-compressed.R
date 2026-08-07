compressed_fixture <- function(multiplier = 1) {
  n <- 80L
  data.frame(
    chromosome = rep("1", n),
    base_pair_location = seq.int(100001L, length.out = n),
    reference_allele = "A",
    alternate_allele = rep(c("C", "G", "T"), length.out = n),
    effect_allele = rep(c("C", "G", "T"), length.out = n),
    other_allele = "A",
    beta = multiplier * seq(-0.2, 0.2, length.out = n),
    standard_error = 0.02 + (seq_len(n) %% 7L) / 1000,
    effect_allele_frequency = seq(0.05, 0.95, length.out = n),
    stringsAsFactors = FALSE
  )
}

test_that("FastMR reads canonical keys from CompreSSoR", {
  skip_if_not_installed("CompreSSoR")

  input <- compressed_fixture()
  store <- tempfile("fastmr-compressed-read-")
  CompreSSoR::compress_sumstats(input, store, overwrite = TRUE)
  keys <- CompreSSoR::compressor_variant_key(
    input$chromosome, input$base_pair_location,
    input$other_allele, input$effect_allele
  )[c(3L, 29L, 70L)]
  got <- fast_read_compressed(store, keys, columns = c("beta", "standard_error"))
  expect_identical(names(got), c("beta", "standard_error", "variant_key"))
  expect_setequal(got$variant_key, keys)
  missing <- fast_read_compressed(
    store, "2:200000000:A:C", columns = c("beta", "standard_error")
  )
  expect_equal(nrow(missing), 0L)
  expect_identical(missing$variant_key, character())
})

test_that("compressed 2 x 2 MR equals the explicitly decoded workflow", {
  skip_if_not_installed("CompreSSoR")

  stores <- vapply(c(1, 1.3, 0.7, -0.4), function(multiplier) {
    path <- tempfile("fastmr-compressed-grid-")
    CompreSSoR::compress_sumstats(compressed_fixture(multiplier), path, overwrite = TRUE)
    path
  }, character(1))
  exposures <- setNames(stores[1:2], c("exposure_a", "exposure_b"))
  outcomes <- setNames(stores[3:4], c("outcome_a", "outcome_b"))
  identity <- compressed_fixture()
  all_keys <- CompreSSoR::compressor_variant_key(
    identity$chromosome, identity$base_pair_location,
    identity$other_allele, identity$effect_allele
  )
  instruments <- list(
    exposure_a = all_keys[c(2L, 7L, 14L, 25L, 40L)],
    exposure_b = all_keys[c(3L, 8L, 19L, 31L, 60L)]
  )

  result_path <- if (requireNamespace("arrow", quietly = TRUE)) {
    tempfile(fileext = ".parquet")
  } else {
    NULL
  }
  compressed <- fast_mr_compressed(
    exposures, outcomes, instruments, methods = "ivw", nboot = 0,
    threads = 1, io_threads = 2, output = result_path
  )
  expect_equal(nrow(compressed), 4L)
  expect_true(all(compressed$nsnp == 5L))
  if (!is.null(result_path)) {
    stored <- compressed
    attr(stored, "compressed_input") <- NULL
    expect_equal(
      as.data.frame(arrow::read_parquet(result_path, as_data_frame = TRUE)),
      stored
    )
  }
  expect_equal(nrow(attr(compressed, "compressed_input")$counts), 4L)
  timing <- attr(compressed, "compressed_input")$timing
  expect_named(
    timing, c(
      "io_seconds", "estimator_seconds", "total_seconds",
      "source_bytes_read"
    ),
    ignore.order = FALSE
  )
  expect_true(all(is.finite(unlist(timing[1:3]))))
  expect_true(all(unlist(timing[1:3]) >= 0))
  expect_true(is.na(timing$source_bytes_read) ||
                (is.finite(timing$source_bytes_read) &&
                   timing$source_bytes_read >= 0))

  exposure <- fast_read_compressed(exposures[[1L]], instruments[[1L]])
  outcome <- fast_read_compressed(outcomes[[1L]], instruments[[1L]])
  hit <- match(exposure$variant_key, outcome$variant_key)
  manual <- fast_mr(data.frame(
    SNP = exposure$variant_key,
    beta.exposure = exposure$beta,
    beta.outcome = outcome$beta[hit],
    se.exposure = exposure$standard_error,
    se.outcome = outcome$standard_error[hit]
  ), methods = "ivw", nboot = 0)
  first <- compressed[compressed$id.exposure == "exposure_a" &
                        compressed$id.outcome == "outcome_a", , drop = FALSE]
  expect_equal(first[c("nsnp", "b", "se", "pval")],
               manual[c("nsnp", "b", "se", "pval")], tolerance = 0)

  serial <- fast_mr_compressed(
    exposures, outcomes, instruments, methods = "ivw", nboot = 0,
    threads = 1, io_threads = 1
  )
  attributes(serial)$compressed_input <- NULL
  attributes(compressed)$compressed_input <- NULL
  expect_identical(serial, compressed)
})

test_that("shared compressed instruments use the exact native grid path", {
  keys <- paste0("1:", 1:6, ":A:C")
  panel <- function(multiplier) data.frame(
    beta = multiplier * seq(0.02, 0.12, length.out = length(keys)),
    standard_error = seq(0.01, 0.015, length.out = length(keys)),
    variant_key = keys,
    stringsAsFactors = FALSE
  )
  exposures <- list(a = panel(1), b = panel(1.2))
  outcomes <- list(x = panel(0.7), y = panel(-0.4))
  instruments <- list(a = keys, b = keys)
  controls <- fastMR:::fastmr_validate_controls(0, 42, 1)
  observed <- fastMR:::fastmr_compressed_grid_fast_path(
    exposures, outcomes, instruments, "ivw", controls, 1L,
    c(a = "a.cpr", b = "b.cpr"), c(x = "x.cpr", y = "y.cpr"),
    1L, list()
  )
  expected <- fast_mr_grid(
    do.call(rbind, lapply(exposures, `[[`, "beta")),
    do.call(rbind, lapply(outcomes, `[[`, "beta")),
    do.call(rbind, lapply(exposures, `[[`, "standard_error")),
    do.call(rbind, lapply(outcomes, `[[`, "standard_error")),
    methods = "ivw", nboot = 0, seed = 42, threads = 1
  )
  expected$exposure_index <- NULL
  expected$outcome_index <- NULL
  attributes(observed)$compressed_input <- NULL
  expect_equal(observed, expected, tolerance = 0)

  bootstrap_controls <- fastMR:::fastmr_validate_controls(100, 42, 1)
  expect_null(fastMR:::fastmr_compressed_grid_fast_path(
    exposures, outcomes, instruments, "weighted_median", bootstrap_controls, 1L,
    c(a = "a.cpr", b = "b.cpr"), c(x = "x.cpr", y = "y.cpr"),
    1L, list()
  ))
  expect_s3_class(fastMR:::fastmr_compressed_grid_fast_path(
    exposures, outcomes, instruments, "ivw", bootstrap_controls, 1L,
    c(a = "a.cpr", b = "b.cpr"), c(x = "x.cpr", y = "y.cpr"),
    1L, list()
  ), "data.frame")

  exposures$b$variant_key[[1L]] <- "1:100:A:C"
  expect_null(fastMR:::fastmr_compressed_grid_fast_path(
    exposures, outcomes, instruments, "ivw", controls, 1L,
    c(a = "a.cpr", b = "b.cpr"), c(x = "x.cpr", y = "y.cpr"),
    1L, list()
  ))
})

test_that("compressed MR rejects malformed contracts and missing instruments", {
  skip_if_not_installed("CompreSSoR")
  expect_error(fastMR:::fastmr_normalize_variant_keys("rs123"),
               "chromosome:position:REF:ALT")
  expect_error(fastMR:::fastmr_normalize_variant_keys(c("1:1:A:C", "1:1:A:C")),
               "must not contain duplicates")
  expect_error(fastMR:::fastmr_normalize_instruments(list(a = "1:1:A:C"), c("a", "b")),
               "one list element per exposure")
  expect_error(fastMR:::fastmr_positive_integer_scalar("2", "io_threads"),
               "positive integer")
  expect_error(fast_mr_compressed(character(), character(), "1:1:A:C"),
               "non-empty")

  fake <- tempfile("fastmr-parquet-store-")
  dir.create(fake)
  writeLines('{"format":"CompreSSoR","backend":"parquet"}',
             file.path(fake, "manifest.json"))
  expect_error(fast_read_compressed(fake), "requires a self-contained")
})

test_that("compressed reads reject incompatible effect-orientation contracts", {
  skip_if_not_installed("CompreSSoR")

  path <- tempfile("fastmr-incompatible-contract-")
  CompreSSoR::compress_sumstats(compressed_fixture(), path, overwrite = TRUE)
  manifest_path <- file.path(path, "manifest.json")
  manifest <- CompreSSoR:::read_manifest(manifest_path)
  manifest$identity$effect_allele_is_alt <- FALSE
  CompreSSoR:::write_manifest(manifest, manifest_path)
  CompreSSoR:::seal_pcodec_manifest(manifest_path)
  expect_error(fast_read_compressed(path), "ALT-oriented")
})

test_that("compressed MR reports missing instruments and parallel read failures", {
  skip_if_not_installed("CompreSSoR")

  outcome_input <- compressed_fixture(0.5)
  outcome <- tempfile("fastmr-valid-outcome-")
  broken <- tempfile("fastmr-broken-exposure-")
  for (path in c(outcome, broken)) {
    CompreSSoR::compress_sumstats(outcome_input, path, overwrite = TRUE)
  }
  keys <- CompreSSoR::compressor_variant_key(
    outcome_input$chromosome, outcome_input$base_pair_location,
    outcome_input$other_allele, outcome_input$effect_allele
  )[1:5]
  absent <- "2:200000000:A:C"
  expect_error(
    fast_mr_compressed(c(exp = outcome), c(out = outcome), c(keys, absent)),
    "missing requested exposure"
  )
  expect_warning(
    missing_non_strict <- fast_mr_compressed(
      c(exp = outcome), c(out = outcome), c(keys, absent), strict = FALSE
    ),
    "omitted missing requested"
  )
  expect_equal(missing_non_strict$nsnp, length(keys))

  broken_store <- CompreSSoR::open_compressor(broken)
  unlink(file.path(broken, broken_store$manifest$files$z))
  if (.Platform$OS.type != "windows") {
    expect_error(
      suppressWarnings(fast_mr_compressed(
        c(good = outcome, broken = broken), c(out = outcome), keys,
        io_threads = 2
      )),
      "compressed read failed|failed to read"
    )
  }
})
