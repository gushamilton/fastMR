test_that("Parquet roundtrip works when Arrow is installed", {
  skip_if_not_installed("arrow")
  d <- il6_fixture()[1:6, c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
  path <- tempfile(fileext = ".parquet")
  arrow::write_parquet(d, path)
  expect_equal(fast_read_parquet(path)$SNP, d$SNP)
  result <- fast_mr_parquet(path, methods = "ivw", nboot = 0)
  expect_equal(result$nsnp, 6)
})

test_that("fastMR writes safe compressed Parquet results", {
  skip_if_not_installed("arrow")
  d <- il6_fixture()[1:6, c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
  expected <- fast_mr(d, methods = "ivw", nboot = 0)
  path <- tempfile(fileext = ".parquet")

  returned_path <- fast_write_parquet(expected, path)
  expect_identical(returned_path, normalizePath(path, mustWork = TRUE))
  expect_equal(as.data.frame(arrow::read_parquet(path, as_data_frame = TRUE)), expected)
  expect_error(fast_write_parquet(expected, path), "destination already exists")
  expect_silent(fast_write_parquet(expected, path, overwrite = TRUE))

  direct_path <- tempfile(fileext = ".parquet")
  observed <- fast_mr(d, methods = "ivw", nboot = 0, output = direct_path)
  expect_equal(as.data.frame(arrow::read_parquet(direct_path, as_data_frame = TRUE)), observed)

  grid_path <- tempfile(fileext = ".parquet")
  grid <- fast_mr_grid(
    exposure_beta = rbind(exposure_a = c(0.1, 0.2)),
    outcome_beta = rbind(outcome_a = c(0.05, 0.1)),
    exposure_se = rbind(exposure_a = c(0.02, 0.02)),
    outcome_se = rbind(outcome_a = c(0.02, 0.02)),
    methods = "ivw", nboot = 0, output = grid_path
  )
  expect_equal(as.data.frame(arrow::read_parquet(grid_path, as_data_frame = TRUE)), grid)
})

test_that("missing Arrow gives a targeted optional-dependency error", {
  skip_if(requireNamespace("arrow", quietly = TRUE), "Arrow is installed in this test environment")
  expect_error(fast_read_parquet("missing.parquet"), "optional 'arrow'")
  expect_error(fast_write_parquet(data.frame(), "results.parquet"), "optional 'arrow'")
})
