test_that("Parquet roundtrip works when Arrow is installed", {
  skip_if_not_installed("arrow")
  d <- il6_fixture()[1:6, c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
  path <- tempfile(fileext = ".parquet")
  arrow::write_parquet(d, path)
  expect_equal(fast_read_parquet(path)$SNP, d$SNP)
  result <- fast_mr_parquet(path, methods = "ivw", nboot = 0)
  expect_equal(result$nsnp, 6)
})

test_that("missing Arrow gives a targeted optional-dependency error", {
  skip_if(requireNamespace("arrow", quietly = TRUE), "Arrow is installed in this test environment")
  expect_error(fast_read_parquet("missing.parquet"), "optional 'arrow'")
})
