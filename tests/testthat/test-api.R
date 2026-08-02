test_that("fast_mr returns tidy TwoSampleMR-style results", {
  d <- il6_fixture()
  result <- fast_mr(d, nboot = 0, seed = 20260801)
  expect_s3_class(result, "data.frame")
  expect_equal(result$method_code, c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode"))
  expect_equal(result$nsnp, rep(82, 5))
  expect_named(result, c(
    "id.exposure", "id.outcome", "method", "method_code", "nsnp", "b", "se", "pval",
    "Q", "Q_df", "Q_pval", "sigma", "intercept", "intercept_se", "intercept_pval",
    "ratio_se_mean", "bootstrap", "phi", "flipped", "se_exposure_mean"
  ))
  expect_equal(result$b[result$method_code == "ivw"], 0.0012105753, tolerance = 1e-7)
  expect_equal(result$b[result$method_code == "egger"], 0.0034024061, tolerance = 1e-7)
})

test_that("dispatcher accepts common TwoSampleMR method spellings", {
  d <- il6_fixture()[1:4, ]
  result <- fast_mr(d, methods = c("mr_ivw", "mr_wald_ratio"), nboot = 0)
  expect_equal(result$method_code, c("ivw", "wald_ratio"))
  expect_equal(result$nsnp[[2]], 4)
  expect_true(is.na(result$b[[2]]))
  expect_true(all(c("code", "method", "description") %in% names(fastmr_method_registry())))
})

test_that("simple median is available through the native kernel", {
  d <- il6_fixture()[1:8, ]
  result <- fast_mr(d, methods = "mr_simple_median", nboot = 0)
  expect_equal(result$method_code, "simple_median")
  expect_equal(result$nsnp, 8)
  expect_true(is.finite(result$b))
  expect_true(is.na(result$ratio_se_mean))
})

test_that("MR-Egger bootstrap matches the native seeded method", {
  d <- il6_fixture()[1:8, ]
  result <- fast_mr(d, methods = "mr_egger_regression_bootstrap",
                    nboot = 20, seed = 20260805)
  expect_equal(result$method_code, "egger_bootstrap")
  expect_equal(result$nsnp, 8)
  expect_true(is.finite(result$b))
  expect_true(is.finite(result$se))
  expect_true(is.finite(result$intercept))
})

test_that("MR-Egger bootstrap handles zero and negative exposure effects", {
  d <- il6_fixture()[1:6, ]
  d$beta.exposure[c(1, 2)] <- c(0, -abs(d$beta.exposure[[2]]))
  result <- fast_mr(d, methods = "egger_bootstrap", nboot = 20, seed = 20260806)
  expect_true(is.finite(result$b))
  expect_true(is.finite(result$se))
  expect_true(is.finite(result$pval))
})

test_that("MR-Egger bootstrap retains three instruments with an exact zero", {
  d <- il6_fixture()[1:3, ]
  d$beta.exposure[[1]] <- 0
  result <- fast_mr(d, methods = "egger_bootstrap", nboot = 20, seed = 20260807)
  expect_equal(result$nsnp, 3)
  expect_true(is.finite(result$b))
  expect_true(is.finite(result$se))
})

test_that("duplicate methods fail and duplicate SNPs are counted once", {
  d <- il6_fixture()[1:4, ]
  expect_error(fast_mr(d, methods = c("ivw", "mr_ivw"), nboot = 0), "unique")
  d$SNP[[2]] <- d$SNP[[1]]
  d$pval.exposure[[2]] <- d$pval.exposure[[1]]
  result <- fast_mr(d, methods = "ivw", nboot = 0)
  expect_equal(result$nsnp, 3)
})

test_that("mode diagnostics report the requested phi", {
  d <- il6_fixture()[1:8, ]
  result <- fast_mr(d, methods = "simple_mode", phi = 0.25, nboot = 0)
  expect_equal(result$phi, 0.25)
})

test_that("tidy IDs remain IDs when trait labels differ", {
  d <- il6_fixture()[1:4, ]
  d$id.exposure <- "ieu-a-123"
  d$id.outcome <- "ieu-b-456"
  d$exposure <- "Human-readable exposure"
  d$outcome <- "Human-readable outcome"
  result <- fast_mr(d, methods = "ivw", nboot = 0)
  expect_equal(result$id.exposure, "ieu-a-123")
  expect_equal(result$id.outcome, "ieu-b-456")
})

test_that("unweighted regression and sign concordance are available", {
  d <- il6_fixture()[1:8, ]
  result <- fast_mr(d, methods = c("mr_uwr", "mr_sign"), nboot = 0)
  expect_equal(result$method_code, c("uwr", "sign"))
  expect_true(all(is.finite(result$b)))
  expect_true(all(is.finite(result$pval)))
})

test_that("grouping preserves tidy pair metadata", {
  d <- il6_fixture()[1:8, ]
  d$id.exposure <- rep(c("A", "B"), each = 4)
  d$id.outcome <- rep(c("Y", "Z"), 4)
  result <- fast_mr(d, methods = "ivw", nboot = 0)
  expect_equal(nrow(result), 4)
  expect_equal(paste(result$id.exposure, result$id.outcome), c("A Y", "A Z", "B Y", "B Z"))
})
