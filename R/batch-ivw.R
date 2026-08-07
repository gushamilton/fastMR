# Validation and public wrappers for the low-memory all-pairs IVW kernels.

fastmr_validate_batch_limit <- function(value, name, integer = FALSE) {
  if (length(value) != 1L || is.na(value) || !is.finite(value) || value <= 0 ||
      (integer && value != floor(value))) {
    suffix <- if (integer) " positive integer" else " positive finite number"
    stop(name, " must be one", suffix, call. = FALSE)
  }
  as.numeric(value)
}

fastmr_batch_memory_bytes <- function(exposure_count, outcome_count, snp_count,
                                      sparse) {
  pair_count <- as.double(exposure_count) * as.double(outcome_count)
  if (sparse) {
    # Five double result matrices are materialised by the native wrapper.
    40 * pair_count
  } else {
    # Five result matrices, four pair accumulators, and the dense transformed
    # exposure/outcome work arrays used by the BLAS path.
    72 * pair_count +
      24 * as.double(exposure_count) * as.double(snp_count) +
      32 * as.double(outcome_count) * as.double(snp_count)
  }
}

fastmr_check_batch_bounds <- function(exposure_count, outcome_count, snp_count,
                                      max_output_cells, max_memory_mb, sparse) {
  max_output_cells <- fastmr_validate_batch_limit(
    max_output_cells, "max_output_cells", integer = TRUE
  )
  max_memory_mb <- fastmr_validate_batch_limit(max_memory_mb, "max_memory_mb")
  output_cells <- as.double(exposure_count) * as.double(outcome_count)
  if (!is.finite(output_cells) || output_cells > max_output_cells) {
    stop("exposure-by-outcome output has ", format(output_cells, scientific = FALSE),
         " cells; increase max_output_cells or process outcomes in batches",
         call. = FALSE)
  }
  estimated_bytes <- fastmr_batch_memory_bytes(
    exposure_count, outcome_count, snp_count, sparse
  )
  if (!is.finite(estimated_bytes) || estimated_bytes > max_memory_mb * 1024^2) {
    stop("estimated native IVW workspace is ",
         format(estimated_bytes / 1024^2, digits = 6),
         " MiB; increase max_memory_mb or process outcomes in batches",
         call. = FALSE)
  }
  invisible(list(max_output_cells = max_output_cells,
                 max_memory_mb = max_memory_mb))
}

fastmr_validate_batch_masks <- function(exposure_beta, outcome_beta, outcome_se,
                                        exposure_present, outcome_present) {
  exposure_present <- as.matrix(exposure_present)
  outcome_present <- as.matrix(outcome_present)
  if (!is.logical(exposure_present) || !is.logical(outcome_present)) {
    stop("presence masks must be logical matrices", call. = FALSE)
  }
  if (anyNA(exposure_present) || anyNA(outcome_present)) {
    stop("presence masks must not contain NA", call. = FALSE)
  }
  if (!identical(dim(exposure_beta), dim(exposure_present)) ||
      !identical(dim(outcome_beta), dim(outcome_present)) ||
      !identical(dim(outcome_beta), dim(outcome_se))) {
    stop("masked IVW inputs must have matching dimensions", call. = FALSE)
  }
  if (nrow(exposure_beta) == 0L || nrow(outcome_beta) == 0L ||
      ncol(exposure_beta) == 0L) {
    stop("masked IVW inputs must have non-empty dimensions", call. = FALSE)
  }
  exposure_snps <- colnames(exposure_beta)
  outcome_snps <- colnames(outcome_beta)
  if (xor(is.null(exposure_snps), is.null(outcome_snps)) ||
      (!is.null(exposure_snps) && !identical(exposure_snps, outcome_snps))) {
    stop("exposure and outcome matrices must use the same SNP column names and order",
         call. = FALSE)
  }
  if (any(exposure_present & !is.finite(exposure_beta))) {
    stop("present exposure_beta values must be finite", call. = FALSE)
  }
  valid_outcome <- outcome_present &
    (!is.finite(outcome_beta) | !is.finite(outcome_se) | outcome_se <= 0)
  if (any(valid_outcome)) {
    stop("present outcome values must be finite with positive standard errors",
         call. = FALSE)
  }
  invisible(list(exposure_present = exposure_present,
                 outcome_present = outcome_present))
}

fastmr_batch_dimnames <- function(result, exposure_beta, outcome_beta) {
  exposure_names <- rownames(exposure_beta)
  outcome_names <- rownames(outcome_beta)
  if (is.null(exposure_names)) exposure_names <- seq_len(nrow(exposure_beta))
  if (is.null(outcome_names)) outcome_names <- seq_len(nrow(outcome_beta))
  for (name in c("nsnp", "beta", "se", "Q", "sigma")) {
    dimnames(result[[name]]) <- list(exposure_names, outcome_names)
  }
  result
}

fastmr_validate_csr_integer <- function(x, name, allow_empty = TRUE) {
  if (!is.atomic(x) || is.object(x) || is.complex(x) ||
      (length(x) == 0L && !allow_empty)) {
    stop(name, " must be an integer vector", call. = FALSE)
  }
  if (!is.integer(x)) {
    if (!is.numeric(x) || any(!is.finite(x)) || any(x != floor(x)) ||
        any(x < -.Machine$integer.max - 1) || any(x > .Machine$integer.max)) {
      stop(name, " must contain finite 32-bit integer values", call. = FALSE)
    }
    x <- as.integer(x)
  }
  if (anyNA(x)) stop(name, " must not contain NA", call. = FALSE)
  x
}

fastmr_validate_csr <- function(row_ptr, col_index, exposure_beta,
                                outcome_beta, outcome_se, outcome_present) {
  row_ptr <- fastmr_validate_csr_integer(row_ptr, "row_ptr", allow_empty = FALSE)
  col_index <- fastmr_validate_csr_integer(col_index, "col_index")
  if (length(row_ptr) < 2L) stop("row_ptr must have at least two entries", call. = FALSE)
  if (row_ptr[[1L]] != 0L || any(row_ptr < 0L) || any(diff(row_ptr) < 0L)) {
    stop("row_ptr must start at zero and be non-decreasing", call. = FALSE)
  }
  if (row_ptr[[length(row_ptr)]] != length(col_index) ||
      length(exposure_beta) != length(col_index)) {
    stop("row_ptr, col_index, and exposure_beta have incompatible lengths",
         call. = FALSE)
  }
  if (any(!is.finite(exposure_beta))) {
    stop("exposure_beta must contain finite values", call. = FALSE)
  }
  snp_count <- ncol(outcome_beta)
  if (any(col_index < 0L) || any(col_index >= snp_count)) {
    stop("col_index must contain zero-based SNP indices in the outcome range",
         call. = FALSE)
  }
  for (i in seq_len(length(row_ptr) - 1L)) {
    first <- row_ptr[[i]] + 1L
    last <- row_ptr[[i + 1L]]
    indices <- if (first <= last) col_index[seq.int(first, last)] else integer()
    if (anyDuplicated(indices)) {
      stop("CSR rows must not contain duplicate SNP indices", call. = FALSE)
    }
  }
  valid_outcome <- outcome_present &
    (!is.finite(outcome_beta) | !is.finite(outcome_se) | outcome_se <= 0)
  if (any(valid_outcome)) {
    stop("present outcome values must be finite with positive standard errors",
         call. = FALSE)
  }
  list(row_ptr = row_ptr, col_index = col_index,
       exposure_beta = as.numeric(exposure_beta))
}

#' Run low-memory IVW over masked exposure and outcome matrices
#'
#' @param exposure_beta Numeric `E x S` matrix of exposure SNP effects.
#' @param outcome_beta Numeric `O x S` matrix of outcome SNP effects.
#' @param outcome_se Numeric `O x S` matrix of outcome standard errors.
#' @param exposure_present Logical `E x S` matrix. Only `TRUE` entries are
#'   used; beta values at `FALSE` entries may be `NA`.
#' @param outcome_present Logical `O x S` matrix. Only `TRUE` entries are used;
#'   corresponding beta values must be finite and standard errors positive.
#' @param threads Positive integer retained for interface consistency. The
#'   dense path delegates matrix multiplication to the configured BLAS.
#' @param max_output_cells Maximum allowed `E * O` result cells. The default
#'   prevents accidentally materialising an unbounded all-by-all result.
#' @param max_memory_mb Conservative bound for native temporary workspace and
#'   result matrices, excluding the input matrices already held by R.
#' @return A `fastmr_masked_ivw_compact` list containing `nsnp`, `beta`, `se`,
#'   `Q`, and `sigma`, each an `E x O` matrix.
#' @export
fast_mr_masked_ivw <- function(exposure_beta, outcome_beta, outcome_se,
                               exposure_present, outcome_present,
                               threads = 1L, max_output_cells = 1e8,
                               max_memory_mb = 2048) {
  exposure_beta <- fastmr_matrix_numeric(exposure_beta, "exposure_beta")
  outcome_beta <- fastmr_matrix_numeric(outcome_beta, "outcome_beta")
  outcome_se <- fastmr_matrix_numeric(outcome_se, "outcome_se")
  masks <- fastmr_validate_batch_masks(
    exposure_beta, outcome_beta, outcome_se, exposure_present, outcome_present
  )
  controls <- fastmr_validate_controls(0L, NULL, threads)
  fastmr_check_batch_bounds(
    nrow(exposure_beta), nrow(outcome_beta), ncol(exposure_beta),
    max_output_cells, max_memory_mb, sparse = FALSE
  )
  result <- fastmr_native_call(
    fastmr_masked_ivw_native,
    list(exposure_beta = exposure_beta, outcome_beta = outcome_beta,
         outcome_se = outcome_se,
         exposure_present = masks$exposure_present,
         outcome_present = masks$outcome_present,
         threads = controls$threads),
    NULL
  )
  fastmr_batch_dimnames(result, exposure_beta, outcome_beta)
}

#' Run low-memory IVW over a sparse CSR exposure panel
#'
#' @param row_ptr Integer CSR row offsets of length `E + 1`, starting at zero.
#' @param col_index Integer zero-based SNP column indices for each stored
#'   exposure effect. Indices must be unique within each exposure row.
#' @param exposure_beta Numeric vector of stored exposure effects, parallel to
#'   `col_index`.
#' @param outcome_beta Numeric `O x S` matrix of outcome SNP effects.
#' @param outcome_se Numeric `O x S` matrix of outcome standard errors.
#' @param outcome_present Logical `O x S` matrix. Only `TRUE` entries are used;
#'   corresponding beta values must be finite and standard errors positive.
#' @param threads Positive integer native worker count over exposure rows.
#' @param max_output_cells Maximum allowed `E * O` result cells. The default
#'   prevents accidentally materialising an unbounded all-by-all result.
#' @param max_memory_mb Conservative bound for native result workspace,
#'   excluding input matrices already held by R.
#' @return A `fastmr_sparse_ivw_compact` list containing `nsnp`, `beta`, `se`,
#'   `Q`, and `sigma`, each an `E x O` matrix.
#' @export
fast_mr_sparse_ivw <- function(row_ptr, col_index, exposure_beta,
                               outcome_beta, outcome_se, outcome_present,
                               threads = 1L, max_output_cells = 1e8,
                               max_memory_mb = 2048) {
  outcome_beta <- fastmr_matrix_numeric(outcome_beta, "outcome_beta")
  outcome_se <- fastmr_matrix_numeric(outcome_se, "outcome_se")
  outcome_present <- as.matrix(outcome_present)
  if (!is.logical(outcome_present)) {
    stop("outcome_present must be a logical matrix", call. = FALSE)
  }
  if (anyNA(outcome_present)) {
    stop("outcome_present must not contain NA", call. = FALSE)
  }
  if (!identical(dim(outcome_beta), dim(outcome_se)) ||
      !identical(dim(outcome_beta), dim(outcome_present))) {
    stop("sparse IVW outcome matrices must have matching dimensions", call. = FALSE)
  }
  if (nrow(outcome_beta) == 0L || ncol(outcome_beta) == 0L) {
    stop("sparse IVW outcome matrices must have non-empty dimensions", call. = FALSE)
  }
  if (is.factor(exposure_beta)) exposure_beta <- as.numeric(as.character(exposure_beta))
  if (!is.numeric(exposure_beta)) {
    converted <- suppressWarnings(as.numeric(exposure_beta))
    if (anyNA(converted) && any(!is.na(exposure_beta))) {
      stop("exposure_beta must be numeric", call. = FALSE)
    }
    exposure_beta <- converted
  }
  csr <- fastmr_validate_csr(
    row_ptr, col_index, exposure_beta, outcome_beta, outcome_se,
    outcome_present
  )
  controls <- fastmr_validate_controls(0L, NULL, threads)
  fastmr_check_batch_bounds(
    length(csr$row_ptr) - 1L, nrow(outcome_beta), ncol(outcome_beta),
    max_output_cells, max_memory_mb, sparse = TRUE
  )
  result <- fastmr_native_call(
    fastmr_sparse_ivw_native,
    list(row_ptr = csr$row_ptr, col_index = csr$col_index,
         exposure_beta = csr$exposure_beta, outcome_beta = outcome_beta,
         outcome_se = outcome_se, outcome_present = outcome_present,
         threads = controls$threads),
    NULL
  )
  exposure_names <- NULL
  outcome_names <- rownames(outcome_beta)
  if (is.null(outcome_names)) outcome_names <- seq_len(nrow(outcome_beta))
  for (name in c("nsnp", "beta", "se", "Q", "sigma")) {
    dimnames(result[[name]]) <- list(exposure_names, outcome_names)
  }
  result
}
