issue18_data <- function() {
  set.seed(1)
  n <- 10
  bx <- rnorm(n, 0.2, 0.05)
  data.frame(
    SNP = paste0("s", seq_len(n)), id.exposure = "e", id.outcome = "o",
    exposure = "e", outcome = "o", beta.exposure = bx, se.exposure = 0.01,
    beta.outcome = 0.1 * bx + rnorm(n, 0, 0.01), se.outcome = 0.01,
    mr_keep = TRUE
  )
}

test_that("registry flags exactly the bootstrap-dependent methods", {
  registry <- fastmr_method_registry()
  expect_type(registry$bootstrap, "logical")
  expect_setequal(
    registry$code[registry$bootstrap],
    c("egger_bootstrap", "simple_median", "weighted_median",
      "penalised_weighted_median", "simple_mode", "weighted_mode")
  )
})

test_that("fast_mr warns when bootstrap methods are requested with nboot = 0", {
  d <- issue18_data()
  boot <- fastMR:::fastmr_bootstrap_methods()
  for (method in boot) {
    expect_warning(
      result <- fast_mr(d, methods = method, nboot = 0),
      class = "fastmr_nboot_warning"
    )
    # egger_bootstrap has no point estimate without draws; the others keep b.
    if (method != "egger_bootstrap") expect_true(is.finite(result$b), info = method)
    expect_true(is.na(result$se), info = method)
    expect_true(is.na(result$pval), info = method)
  }
  expect_warning(
    fast_mr(d, methods = c("ivw", "weighted_median", "simple_mode"), nboot = 0),
    "weighted_median, simple_mode"
  )
})

test_that("bootstrap methods give finite inference without a warning when nboot > 0", {
  d <- issue18_data()
  expect_no_warning(
    result <- fast_mr(d, methods = "weighted_median", nboot = 200, seed = 1)
  )
  expect_true(is.finite(result$se))
  expect_true(is.finite(result$pval))
})

test_that("non-bootstrap methods do not warn with nboot = 0", {
  d <- issue18_data()
  expect_no_warning(fast_mr(d, methods = c("ivw", "egger", "uwr", "sign"), nboot = 0))
  g <- grid_fixture()
  expect_no_warning(fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se,
                                 g$outcome_se, methods = "ivw", nboot = 0))
})

test_that("fast_mr_grid warns when bootstrap methods are requested with nboot = 0", {
  g <- grid_fixture()
  expect_warning(
    fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se,
                 methods = c("ivw", "weighted_median"), nboot = 0),
    class = "fastmr_nboot_warning"
  )
})

test_that("fast_mr_compressed warns once for bootstrap methods at its default nboot", {
  skip_if_compressor_unavailable()
  stores <- vapply(c(1, 0.7), function(multiplier) {
    path <- tempfile("fastmr-compressed-nboot-")
    CompreSSoR::compress_sumstats(compressor_canonical_fixture(multiplier), path,
                                  overwrite = TRUE)
    path
  }, character(1))
  identity <- compressor_canonical_fixture()
  keys <- CompreSSoR::compressor_variant_key(
    identity$chromosome, identity$base_pair_location,
    identity$other_allele, identity$effect_allele
  )[c(2L, 7L, 14L, 25L, 40L)]
  exposures <- c(exposure_a = stores[[1L]])
  outcomes <- c(outcome_a = stores[[2L]])

  warnings <- list()
  result <- withCallingHandlers(
    fast_mr_compressed(exposures, outcomes, keys,
                       methods = c("ivw", "weighted_median")),
    warning = function(w) {
      warnings[[length(warnings) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  nboot_warnings <- Filter(function(w) inherits(w, "fastmr_nboot_warning"), warnings)
  expect_length(nboot_warnings, 1L)
  expect_match(conditionMessage(nboot_warnings[[1L]]), "weighted_median")
  expect_true(is.finite(result$b[result$method_code == "weighted_median"]))

  expect_no_warning(fast_mr_compressed(exposures, outcomes, keys,
                                       methods = c("ivw", "weighted_median"),
                                       nboot = 50, seed = 1))
  expect_no_warning(fast_mr_compressed(exposures, outcomes, keys, methods = "ivw"))
})
