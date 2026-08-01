#' Read a tidy Parquet file with Arrow
#'
#' Arrow is optional so that the core package remains small. Install it with
#' `install.packages("arrow")` when Parquet input is needed.
#' @param path Path to a Parquet file.
#' @return A data frame returned by `arrow::read_parquet()`.
#' @export
fast_read_parquet <- function(path) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("fast_read_parquet() requires the optional 'arrow' package; install.packages('arrow')", call. = FALSE)
  }
  if (length(path) != 1L || !is.character(path) || !file.exists(path)) {
    stop("path must identify an existing Parquet file", call. = FALSE)
  }
  arrow::read_parquet(path, as_data_frame = TRUE)
}

#' Run fastMR on tidy Parquet input
#'
#' @param path Path to a Parquet file with TwoSampleMR-style columns.
#' @param ... Arguments forwarded to [fast_mr()].
#' @return A tidy fastMR result data frame.
#' @export
fast_mr_parquet <- function(path, ...) {
  fast_mr(fast_read_parquet(path), ...)
}
