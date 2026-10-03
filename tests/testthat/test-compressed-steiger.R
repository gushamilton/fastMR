compressed_steiger_study <- function() {
  stores <- vapply(c(1, 1.3, 0.7, -0.4), function(multiplier) {
    path <- tempfile("fastmr-compressed-steiger-")
    CompreSSoR::compress_sumstats(
      compressor_canonical_fixture(multiplier), path, overwrite = TRUE
    )
    path
  }, character(1))
  identity <- compressor_canonical_fixture()
  keys <- CompreSSoR::compressor_variant_key(
    identity$chromosome, identity$base_pair_location,
    identity$other_allele, identity$effect_allele
  )
  list(
    exposures = setNames(stores[1:2], c("exposure_a", "exposure_b")),
    outcomes = setNames(stores[3:4], c("outcome_a", "outcome_b")),
    keys = keys,
    instruments = list(
      exposure_a = keys[c(2L, 7L, 14L, 25L, 40L)],
      exposure_b = keys[c(3L, 8L, 19L, 31L, 60L)]
    ),
    n_exposure = c(exposure_a = 50000, exposure_b = 80000),
    n_outcome = c(outcome_a = 120000, outcome_b = 30000)
  )
}

# The two-pass workflow: read beta/SE/EAF/p again from every store, join the
# exposure instruments to the outcome rows, then call fast_mr_steiger_filtering.
compressed_steiger_two_pass <- function(study, instruments, units = NULL) {
  columns <- c("beta", "standard_error", "effect_allele_frequency", "p_value")
  rows <- list()
  for (exposure in names(study$exposures)) {
    ex <- fast_read_compressed(study$exposures[[exposure]],
                               instruments[[exposure]], columns)
    ex <- ex[match(instruments[[exposure]], ex$variant_key), , drop = FALSE]
    for (outcome in names(study$outcomes)) {
      ou <- fast_read_compressed(study$outcomes[[outcome]],
                                 instruments[[exposure]], columns)
      hit <- match(ex$variant_key, ou$variant_key)
      keep <- !is.na(hit)
      rows[[length(rows) + 1L]] <- data.frame(
        SNP = ex$variant_key[keep], id.exposure = exposure, id.outcome = outcome,
        exposure = exposure, outcome = outcome,
        beta.exposure = ex$beta[keep], beta.outcome = ou$beta[hit[keep]],
        se.exposure = ex$standard_error[keep],
        se.outcome = ou$standard_error[hit[keep]],
        eaf.exposure = ex$effect_allele_frequency[keep],
        eaf.outcome = ou$effect_allele_frequency[hit[keep]],
        pval.exposure = ex$p_value[keep], pval.outcome = ou$p_value[hit[keep]],
        samplesize.exposure = unname(study$n_exposure[exposure]),
        samplesize.outcome = unname(study$n_outcome[outcome]),
        units.exposure = "", units.outcome = "", mr_keep = TRUE,
        stringsAsFactors = FALSE
      )
    }
  }
  fast_mr_steiger_filtering(do.call(rbind, rows))
}

test_that("one-pass compressed Steiger equals the two-pass workflow", {
  skip_if_compressor_unavailable()
  study <- compressed_steiger_study()
  for (methods in list("ivw", c("ivw", "egger"))) {
    result <- fast_mr_compressed(
      study$exposures, study$outcomes, study$instruments, methods = methods,
      steiger = TRUE, samplesize_exposure = study$n_exposure,
      samplesize_outcome = study$n_outcome
    )
    steiger <- attr(result, "steiger")
    expected <- compressed_steiger_two_pass(study, study$instruments)
    rownames(steiger) <- rownames(expected) <- NULL
    expect_identical(steiger, expected)
    expect_true(all(c("steiger_dir", "steiger_pval", "rsq.exposure",
                      "rsq.outcome") %in% names(steiger)))
    expect_true(all(steiger$steiger_pval_valid))
    timing <- attr(result, "compressed_input")$timing
    expect_true(is.finite(timing$steiger_seconds) && timing$steiger_seconds >= 0)
  }
})

test_that("Steiger is attached on every estimator path without changing MR", {
  skip_if_compressor_unavailable()
  study <- compressed_steiger_study()
  shared <- study$keys[c(2L, 7L, 14L, 25L, 40L)]
  cases <- list(
    list(instruments = shared, methods = "ivw", estimator = "auto",
         path = "shared_instrument_grid"),
    list(instruments = study$instruments, methods = "ivw", estimator = "auto",
         path = "sparse_ivw"),
    list(instruments = study$instruments, methods = "ivw",
         estimator = "pairwise", path = "pairwise")
  )
  for (case in cases) {
    plain <- fast_mr_compressed(
      study$exposures, study$outcomes, case$instruments,
      methods = case$methods, estimator = case$estimator
    )
    with_steiger <- fast_mr_compressed(
      study$exposures, study$outcomes, case$instruments,
      methods = case$methods, estimator = case$estimator, steiger = TRUE,
      samplesize_exposure = study$n_exposure, samplesize_outcome = 75000
    )
    expect_identical(attr(plain, "compressed_input")$estimator_path, case$path)
    expect_null(attr(plain, "steiger"))
    expect_named(
      attr(plain, "compressed_input")$timing,
      c("io_seconds", "estimator_seconds", "total_seconds", "source_bytes_read")
    )
    strip <- function(x) {
      attr(x, "compressed_input") <- NULL
      attr(x, "steiger") <- NULL
      x
    }
    expect_identical(strip(with_steiger), strip(plain))
    expect_identical(
      attr(with_steiger, "compressed_input")$counts,
      attr(plain, "compressed_input")$counts
    )
    steiger <- attr(with_steiger, "steiger")
    expect_equal(nrow(steiger), 4L * 5L)
    expect_true(all(steiger$samplesize.outcome == 75000))
    expect_identical(
      unique(paste(steiger$id.exposure, steiger$id.outcome)),
      paste(plain$id.exposure, plain$id.outcome)[!duplicated(
        paste(plain$id.exposure, plain$id.outcome))]
    )
  }
})

test_that("compressed Steiger uses the rows MR keeps", {
  skip_if_compressor_unavailable()
  study <- compressed_steiger_study()
  instruments <- study$instruments
  instruments$exposure_b <- c(instruments$exposure_b, "2:200000000:A:C")
  expect_warning(
    result <- fast_mr_compressed(
      study$exposures, study$outcomes, instruments, methods = "ivw",
      strict = FALSE, steiger = TRUE, samplesize_exposure = study$n_exposure,
      samplesize_outcome = unname(study$n_outcome)
    ),
    "missing requested"
  )
  steiger <- attr(result, "steiger")
  expect_false("2:200000000:A:C" %in% steiger$SNP)
  expected <- compressed_steiger_two_pass(study, study$instruments)
  rownames(steiger) <- rownames(expected) <- NULL
  expect_identical(steiger, expected)
  counts <- attr(result, "compressed_input")$counts
  expect_identical(
    as.vector(table(factor(paste(steiger$id.exposure, steiger$id.outcome),
                           levels = unique(paste(counts$id.exposure,
                                                 counts$id.outcome))))),
    counts$matched
  )
})

test_that("compressed Steiger supports binary traits", {
  skip_if_compressor_unavailable()
  study <- compressed_steiger_study()
  binary <- data.frame(id = "outcome_b", ncase = 4000, ncontrol = 26000,
                       prevalence = 0.05)
  result <- fast_mr_compressed(
    study$exposures, study$outcomes, study$instruments, steiger = TRUE,
    samplesize_exposure = study$n_exposure,
    samplesize_outcome = c(outcome_a = 120000), steiger_binary = binary
  )
  steiger <- attr(result, "steiger")
  manual <- compressed_steiger_two_pass(study, study$instruments)
  manual <- manual[, !grepl("^(rsq|effective_n|steiger)|rsq_", names(manual))]
  b <- manual$id.outcome == "outcome_b"
  manual$units.outcome[b] <- "log odds"
  manual$samplesize.outcome[b] <- NA_real_
  manual$ncase.exposure <- NA_real_
  manual$ncontrol.exposure <- NA_real_
  manual$prevalence.exposure <- NA_real_
  manual$ncase.outcome <- ifelse(b, 4000, NA_real_)
  manual$ncontrol.outcome <- ifelse(b, 26000, NA_real_)
  manual$prevalence.outcome <- ifelse(b, 0.05, NA_real_)
  expected <- fast_mr_steiger_filtering(manual)
  rownames(steiger) <- rownames(expected) <- NULL
  expect_identical(steiger, expected)
  expect_true(all(steiger$effective_n.outcome[b] == 2 / (1 / 4000 + 1 / 26000)))
})

test_that("compressed Steiger options are validated", {
  skip_if_compressor_unavailable()
  study <- compressed_steiger_study()
  run <- function(...) {
    fast_mr_compressed(study$exposures, study$outcomes, study$instruments, ...)
  }
  expect_error(run(steiger = TRUE, samplesize_outcome = 1e5),
               "samplesize_exposure is required")
  expect_error(run(steiger = TRUE, samplesize_exposure = c(exposure_a = 1e5),
                   samplesize_outcome = 1e5), "missing sample size")
  expect_error(run(steiger = TRUE, samplesize_exposure = c(1e5, 2e5, 3e5),
                   samplesize_outcome = 1e5), "length one")
  expect_error(run(steiger = TRUE, samplesize_exposure = -1,
                   samplesize_outcome = 1e5), "finite and positive")
  expect_error(run(samplesize_exposure = 1e5), "only applies when steiger = TRUE")
  expect_error(run(steiger = NA), "steiger must be TRUE or FALSE")
  expect_error(run(steiger = TRUE, samplesize_exposure = 1e5,
                   samplesize_outcome = 1e5,
                   steiger_binary = data.frame(id = "nope", ncase = 1,
                                               ncontrol = 1, prevalence = 0.1)),
               "not exposure or outcome stores")
  expect_error(run(steiger = TRUE, samplesize_exposure = 1e5,
                   samplesize_outcome = 1e5,
                   steiger_binary = data.frame(id = "outcome_a", ncase = 1,
                                               ncontrol = 1, prevalence = 1)),
               "prevalence")
})
