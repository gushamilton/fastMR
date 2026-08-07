#' Read a tidy Parquet file with Arrow
#'
#' Arrow is optional so that the core package remains small. Install it with
#' `install.packages("arrow")` when Parquet input or output is needed.
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

#' Save a fastMR result as a compressed Parquet file
#'
#' Results are written with Zstandard compression by default. This is useful
#' for large MR grids, while retaining column types and allowing selective
#' reads with Arrow-compatible tools.
#'
#' @param x A data frame, normally a result returned by a `fast_mr*()`
#'   function.
#' @param path Destination Parquet file. Existing files are rejected unless
#'   `overwrite = TRUE`.
#' @param compression Parquet compression codec passed to
#'   [arrow::write_parquet()]. Defaults to `"zstd"`.
#' @param compression_level Optional codec compression level.
#' @param overwrite Whether an existing destination may be replaced.
#' @return The normalized destination path, invisibly.
#' @export
fast_write_parquet <- function(
    x,
    path,
    compression = "zstd",
    compression_level = NULL,
    overwrite = FALSE) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("fast_write_parquet() requires the optional 'arrow' package; install.packages('arrow')", call. = FALSE)
  }
  if (!is.data.frame(x)) {
    stop("x must be a data.frame", call. = FALSE)
  }
  if (length(path) != 1L || !is.character(path) || is.na(path) || !nzchar(path)) {
    stop("path must be one non-empty file path", call. = FALSE)
  }
  if (length(compression) != 1L || !is.character(compression) ||
      is.na(compression) || !nzchar(compression)) {
    stop("compression must be one non-empty codec name", call. = FALSE)
  }
  if (!is.null(compression_level) &&
      (!is.numeric(compression_level) || length(compression_level) != 1L ||
       is.na(compression_level) || !is.finite(compression_level))) {
    stop("compression_level must be NULL or one finite number", call. = FALSE)
  }
  if (length(overwrite) != 1L || !is.logical(overwrite) || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE", call. = FALSE)
  }
  path <- path.expand(path)
  parent <- dirname(path)
  if (!dir.exists(parent)) {
    stop("the destination directory does not exist: ", parent, call. = FALSE)
  }
  path <- file.path(normalizePath(parent, mustWork = TRUE), basename(path))
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop("destination already exists; set overwrite = TRUE to replace it", call. = FALSE)
  }
  # `compressed_input` is useful on the returned object, but Arrow otherwise
  # serializes it as R-specific schema metadata. Keep result files portable
  # and avoid duplicating potentially large instrument lists in every output.
  output_data <- x
  attr(output_data, "compressed_input") <- NULL
  arrow::write_parquet(
    output_data,
    sink = path,
    compression = compression,
    compression_level = compression_level
  )
  invisible(path)
}

fastmr_write_result <- function(result, output) {
  if (!is.null(output)) fast_write_parquet(result, output)
  result
}

#' Run fastMR on tidy Parquet input
#'
#' @param path Path to a Parquet file with TwoSampleMR-style columns.
#' @param ... Arguments forwarded to [fast_mr()], including its optional
#'   `output` Parquet result path.
#' @return A tidy fastMR result data frame.
#' @export
fast_mr_parquet <- function(path, ...) {
  fast_mr(fast_read_parquet(path), ...)
}
