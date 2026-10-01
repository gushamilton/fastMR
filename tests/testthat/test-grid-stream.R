rand_grid <- function(E = 5L, O = 4L, S = 12L, seed = 7L) {
  set.seed(seed)
  list(
    exposure_beta = matrix(rnorm(E * S, .05, .02), E, S, dimnames = list(paste0("e", seq_len(E)), NULL)),
    outcome_beta = matrix(rnorm(O * S, 0, .02), O, S, dimnames = list(paste0("o", seq_len(O)), NULL)),
    exposure_se = matrix(runif(E * S, .005, .02), E, S),
    outcome_se = matrix(runif(O * S, .005, .02), O, S)
  )
}
run_grid <- function(g, ...) {
  fast_mr_grid(g$exposure_beta, g$outcome_beta, g$exposure_se, g$outcome_se, ...)
}
read_streamed <- function(path) as.data.frame(arrow::read_parquet(path, as_data_frame = TRUE))

test_that("default return is unchanged tidy data frame", {
  g <- rand_grid()
  r <- run_grid(g, methods = "ivw", nboot = 0)
  expect_s3_class(r, "data.frame")
  expect_equal(nrow(r), 20L)
  expect_identical(run_grid(g, methods = "ivw", nboot = 0, return = "tidy"), r)
})

test_that("compact object converts to the identical tidy data frame", {
  g <- rand_grid()
  for (m in list("ivw", c("ivw", "egger", "weighted_median", "weighted_mode"))) {
    tidy <- run_grid(g, methods = m, nboot = 20, seed = 11)
    cmp <- run_grid(g, methods = m, nboot = 20, seed = 11, return = "compact")
    expect_s3_class(cmp, "fastmr_compact_grid")
    expect_identical(as.data.frame(cmp), tidy)
    expect_identical(fastmr_grid_chunk(cmp, 1, 20), tidy)
    rows <- ((6 - 1) * length(m) + 1):(9 * length(m))
    expect_identical(fastmr_grid_chunk(cmp, 6, 9), tidy[rows, ], ignore_attr = TRUE)
  }
  expect_output(print(cmp), "fastmr_compact_grid")
  expect_error(fastmr_grid_chunk(cmp, 0, 3), "first and last")
})

test_that("streamed Parquet matches tidy output for IVW and mixed methods", {
  skip_if_not_installed("arrow")
  g <- rand_grid()
  cases <- list(ivw = list(methods = "ivw", nboot = 0),
                mixed = list(methods = c("ivw", "egger", "weighted_median", "weighted_mode"),
                             nboot = 15, seed = 99))
  for (cs in cases) {
    tidy <- do.call(run_grid, c(list(g), cs))
    ref <- tempfile(fileext = ".parquet")
    fast_write_parquet(tidy, ref)
    expected <- read_streamed(ref)
    expect_identical(expected, tidy)  # Arrow round trip itself is lossless
    for (cp in c(1, 7, 1e6)) {
      for (ret in c("compact", "none")) {
        path <- tempfile(fileext = ".parquet")
        res <- do.call(run_grid, c(list(g), cs, list(output = path, return = ret, chunk_pairs = cp)))
        if (ret == "none") expect_identical(res, normalizePath(path)) else expect_s3_class(res, "fastmr_compact_grid")
        expect_identical(read_streamed(path), tidy)
        if (cp == 7 && ret == "compact") {
          expect_equal(arrow::ParquetFileReader$create(path)$num_row_groups,
                       ceiling(20 / 7) * 1L)
        }
      }
    }
  }
  # tidy + output still writes the single-table copy
  path <- tempfile(fileext = ".parquet")
  r <- run_grid(g, methods = "ivw", nboot = 0, output = path)
  expect_identical(read_streamed(path), r)
})

test_that("return = none needs output, and existing outputs are protected", {
  skip_if_not_installed("arrow")
  g <- rand_grid()
  expect_error(run_grid(g, methods = "ivw", nboot = 0, return = "none"), "requires output")
  path <- tempfile(fileext = ".parquet")
  run_grid(g, methods = "ivw", nboot = 0, return = "none", output = path)
  expect_error(run_grid(g, methods = "ivw", nboot = 0, return = "none", output = path), "already exists")
  expect_error(run_grid(g, methods = "ivw", nboot = 0, chunk_pairs = 0), "chunk_pairs")
})

test_that("streamed output requires arrow", {
  skip_if(requireNamespace("arrow", quietly = TRUE), "Arrow is installed in this test environment")
  g <- rand_grid()
  expect_error(run_grid(g, methods = "ivw", nboot = 0, return = "none",
                        output = tempfile(fileext = ".parquet")), "optional 'arrow'")
})
