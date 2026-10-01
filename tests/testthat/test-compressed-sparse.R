sparse_store <- function(drop = integer(), invalid_se = integer(), invalid_beta = integer(),
                         multiplier = 1, seed = 1L) {
  set.seed(seed)
  input <- compressor_canonical_fixture(multiplier)
  input$beta <- input$beta + rnorm(nrow(input), 0, 0.05)
  input$standard_error[invalid_se] <- 0
  input$beta[invalid_beta] <- NaN
  if (length(drop)) input <- input[-drop, , drop = FALSE]
  path <- tempfile("fastmr-sparse-")
  CompreSSoR::compress_sumstats(input, path, overwrite = TRUE)
  path
}

sparse_keys <- function() {
  input <- compressor_canonical_fixture()
  CompreSSoR::compressor_variant_key(
    input$chromosome, input$base_pair_location,
    input$other_allele, input$effect_allele
  )
}

run_both <- function(exposures, outcomes, instruments, ..., what = c("value", "warnings")) {
  collect <- function(estimator) {
    warns <- character()
    value <- tryCatch(
      withCallingHandlers(
        fast_mr_compressed(exposures, outcomes, instruments, estimator = estimator, ...),
        warning = function(w) {
          warns <<- c(warns, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) structure(list(message = conditionMessage(e)), class = "captured_error")
    )
    list(value = value, warnings = warns)
  }
  list(auto = collect("auto"), pairwise = collect("pairwise"))
}

expect_same_result <- function(both, path = "sparse_ivw") {
  expect_identical(both$auto$warnings, both$pairwise$warnings)
  a <- both$auto$value
  p <- both$pairwise$value
  if (inherits(p, "captured_error")) {
    expect_s3_class(a, "captured_error")
    expect_identical(a$message, p$message)
    return(invisible())
  }
  expect_false(inherits(a, "captured_error"))
  info_a <- attr(a, "compressed_input")
  info_p <- attr(p, "compressed_input")
  expect_identical(info_a$estimator_path, path)
  expect_identical(info_p$estimator_path, "pairwise")
  expect_identical(info_a$counts, info_p$counts)
  expect_identical(info_a$instruments, info_p$instruments)
  expect_identical(names(info_a), names(info_p))
  expect_identical(names(info_a$timing), names(info_p$timing))
  attr(a, "compressed_input") <- NULL
  attr(p, "compressed_input") <- NULL
  expect_identical(names(a), names(p))
  expect_identical(vapply(a, typeof, ""), vapply(p, typeof, ""))
  expect_identical(a$nsnp, p$nsnp)
  expect_identical(a[c("id.exposure", "id.outcome", "method", "method_code")],
                   p[c("id.exposure", "id.outcome", "method", "method_code")])
  expect_identical(attributes(a)$row.names, attributes(p)$row.names)
  for (column in c("b", "se", "pval")) {
    expect_equal(a[[column]], p[[column]], tolerance = 1e-12)
    expect_identical(is.na(a[[column]]), is.na(p[[column]]))
  }
  for (column in setdiff(names(a), c("b", "se", "pval"))) {
    if (is.numeric(a[[column]])) {
      expect_identical(is.na(a[[column]]), is.na(p[[column]]))
      expect_equal(a[[column]], p[[column]], tolerance = 1e-8)
    }
  }
  invisible()
}

test_that("sparse IVW dispatch matches the pairwise path on ragged instruments", {
  skip_if_compressor_unavailable()
  keys <- sparse_keys()
  exposures <- c(zz = sparse_store(seed = 1), aa = sparse_store(multiplier = 2, seed = 2),
                 mm = sparse_store(multiplier = -1, seed = 3))
  outcomes <- c(o2 = sparse_store(seed = 4), o1 = sparse_store(multiplier = 0.5, seed = 5))
  set.seed(10)
  instruments <- list(zz = keys[sample(80, 9)], aa = keys[sample(80, 4)],
                      mm = keys[sample(80, 6)])
  both <- run_both(exposures, outcomes, instruments)
  expect_equal(nrow(both$auto$value), 6L)
  expect_same_result(both)
  expect_identical(both$auto$value$id.exposure, rep(c("zz", "aa", "mm"), each = 2L))
  expect_identical(both$auto$value$id.outcome, rep(c("o2", "o1"), 3L))
})

test_that("sparse IVW dispatch reproduces missing, invalid and small-pair behaviour", {
  skip_if_compressor_unavailable()
  keys <- sparse_keys()
  exposures <- c(
    e1 = sparse_store(seed = 1),
    e2 = sparse_store(drop = c(5L, 6L), seed = 2),
    e3 = sparse_store(invalid_se = 11L, invalid_beta = 12L, seed = 3),
    e4 = sparse_store(seed = 4)
  )
  outcomes <- c(
    o1 = sparse_store(seed = 5),
    o2 = sparse_store(drop = c(20L, 21L, 7L), seed = 6),
    o3 = sparse_store(invalid_se = c(30L, 2L), invalid_beta = 31L, seed = 7)
  )
  instruments <- list(
    e1 = keys[c(1:4, 20, 30, 40)],            # outcome o2/o3 losses
    e2 = keys[c(5, 6, 8, 9, 10)],             # exposure missing 5,6
    e3 = keys[c(10:14)],                      # exposure invalid 11,12
    e4 = keys[c(2, 21)]                       # two SNPs; o2/o3 knock out
  )
  for (strict in c(TRUE, FALSE)) {
    for (minimum in c(1L, 2L, 3L)) {
      both <- run_both(exposures, outcomes, instruments, strict = strict,
                       minimum_snps = minimum)
      if (strict) {
        expect_s3_class(both$auto$value, "captured_error")
      }
      expect_same_result(both)
    }
  }
  # Subsets that make strict mode succeed.
  clean_exp <- exposures[c("e1", "e4")]
  clean_out <- outcomes["o1"]
  expect_same_result(run_both(clean_exp, clean_out, instruments[c("e1", "e4")],
                              strict = TRUE, minimum_snps = 2L))
})

test_that("sparse IVW dispatch handles pairs with 0, 1 and 2 SNPs", {
  skip_if_compressor_unavailable()
  keys <- sparse_keys()
  exposures <- c(e1 = sparse_store(seed = 1), e2 = sparse_store(seed = 2),
                 e3 = sparse_store(seed = 3))
  outcomes <- c(o1 = sparse_store(drop = 1:10, seed = 4), o2 = sparse_store(seed = 5))
  instruments <- list(e1 = keys[c(1:5)], e2 = keys[c(3, 40)], e3 = keys[50])
  for (minimum in 1:3) {
    both <- run_both(exposures, outcomes, instruments, strict = FALSE,
                     minimum_snps = minimum)
    expect_same_result(both)
  }
  both <- run_both(exposures, outcomes, instruments, strict = FALSE, minimum_snps = 1L)
  expect_true(any(both$auto$value$nsnp == 1))
  expect_true(any(both$auto$value$nsnp == 2))
  expect_true(any(attr(both$auto$value, "compressed_input")$counts$matched == 0L))
  # No pair retained at all.
  none <- run_both(exposures[2:3], outcomes[1], instruments[2:3], strict = FALSE,
                   minimum_snps = 2L)
  expect_s3_class(none$auto$value, "captured_error")
  expect_same_result(none)
})

test_that("non-IVW methods and estimator = 'pairwise' keep the pairwise path", {
  skip_if_compressor_unavailable()
  keys <- sparse_keys()
  exposures <- c(a = sparse_store(seed = 1), b = sparse_store(seed = 2))
  outcomes <- c(x = sparse_store(seed = 3), y = sparse_store(seed = 4))
  instruments <- list(a = keys[1:6], b = keys[4:12])
  for (methods in list(c("ivw", "ivw_fe"), "ivw_mre", c("ivw", "egger"))) {
    res <- fast_mr_compressed(exposures, outcomes, instruments, methods = methods)
    expect_identical(attr(res, "compressed_input")$estimator_path, "pairwise")
  }
  res <- fast_mr_compressed(exposures, outcomes, instruments, estimator = "pairwise")
  expect_identical(attr(res, "compressed_input")$estimator_path, "pairwise")
  expect_identical(
    attr(fast_mr_compressed(exposures, outcomes, instruments),
         "compressed_input")$estimator_path, "sparse_ivw"
  )
  expect_error(fast_mr_compressed(exposures, outcomes, instruments, estimator = "x"),
               "should be one of")
  shared <- list(a = keys[1:6], b = keys[1:6])
  expect_identical(
    attr(fast_mr_compressed(exposures, outcomes, shared), "compressed_input")$estimator_path,
    "shared_instrument_grid"
  )
})
