diagnostic_fixture <- function() {
  d <- il6_fixture()
  d$id.exposure <- "E"
  d$id.outcome <- "O"
  d$exposure <- "E"
  d$outcome <- "O"
  d$mr_keep <- TRUE
  d
}

test_that("heterogeneity and Egger pleiotropy use compiled diagnostics", {
  d <- diagnostic_fixture()
  mr <- fast_mr(d, methods = c("ivw", "egger"), nboot = 0)
  heterogeneity <- fast_mr_heterogeneity(d)
  expect_equal(heterogeneity$method, c("Inverse variance weighted", "MR Egger"))
  expect_equal(heterogeneity$Q, mr$Q, tolerance = 0)
  expect_equal(heterogeneity$Q_df, mr$Q_df, tolerance = 0)
  expect_equal(heterogeneity$Q_pval, mr$Q_pval, tolerance = 0)

  pleiotropy <- fast_mr_pleiotropy_test(d)
  egger <- mr[mr$method_code == "egger", , drop = FALSE]
  expect_equal(pleiotropy$egger_intercept, egger$intercept, tolerance = 0)
  expect_equal(pleiotropy$se, egger$intercept_se, tolerance = 0)
  expect_equal(pleiotropy$pval, egger$intercept_pval, tolerance = 0)
})

test_that("single-SNP and leave-one-out results preserve native tidy shapes", {
  d <- diagnostic_fixture()
  single <- fast_mr_singlesnp(d)
  expect_equal(nrow(single), length(unique(d$SNP)) + 2L)
  expect_equal(single$SNP[1], d$SNP[1])
  expect_equal(tail(single$SNP, 2), c("All - Inverse variance weighted", "All - MR Egger"))
  expect_true(all(c("exposure", "outcome", "SNP", "b", "se", "p") %in% names(single)))

  loo <- fast_mr_leaveoneout(d)
  expect_equal(nrow(loo), length(unique(d$SNP)) + 1L)
  expect_equal(tail(loo$SNP, 1), "All")
  expect_equal(loo$b[loo$SNP == "All"], fast_mr(d, methods = "ivw", nboot = 0)$b,
               tolerance = 0)

  duplicated <- rbind(d, d[1, , drop = FALSE])
  expect_equal(nrow(fast_mr_singlesnp(duplicated)), length(unique(d$SNP)) + 2L)
})

test_that("diagnostic method validation is explicit", {
  d <- diagnostic_fixture()
  expect_error(fast_mr_heterogeneity(d, methods = "weighted_median"),
               "heterogeneity is not defined")
  expect_error(fast_mr_singlesnp(d, single_method = "ivw"),
               "single_method")
  expect_error(fast_mr_leaveoneout(d, method = "weighted_median"),
               "heterogeneity-capable")
})
