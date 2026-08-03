fastmr_require_compressor <- function() {
  if (!requireNamespace("CompreSSoR", quietly = TRUE)) {
    stop(
      "compressed GWAS input requires CompreSSoR; install it from github.com/gushamilton/CompreSSoR",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

fastmr_positive_integer_scalar <- function(value, argument) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 1 || value != floor(value)) {
    stop(argument, " must be one positive integer", call. = FALSE)
  }
  as.integer(value)
}

fastmr_normalize_compressed_files <- function(paths, argument) {
  if (!is.character(paths) || !length(paths) || anyNA(paths) || any(!nzchar(paths))) {
    stop(argument, " must be a non-empty character vector of CompreSSoR stores", call. = FALSE)
  }
  missing <- paths[!dir.exists(paths)]
  if (length(missing)) {
    stop(argument, " contains missing store(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  labels <- names(paths)
  if (is.null(labels)) labels <- rep("", length(paths))
  generated <- !nzchar(labels)
  labels[generated] <- basename(sub("[\\/]+$", "", paths[generated]))
  if (any(!nzchar(labels)) || anyDuplicated(labels)) {
    stop(argument, " must have unique non-empty names (or unique directory basenames)", call. = FALSE)
  }
  normalized <- normalizePath(paths, mustWork = TRUE)
  names(normalized) <- labels
  normalized
}

fastmr_normalize_variant_keys <- function(keys) {
  fastmr_require_compressor()
  if (!is.character(keys) || !length(keys) || anyNA(keys) || any(!nzchar(trimws(keys)))) {
    stop("instrument keys must be non-empty chromosome:position:REF:ALT strings", call. = FALSE)
  }
  fields <- strsplit(trimws(keys), ":", fixed = TRUE)
  if (any(lengths(fields) != 4L)) {
    stop("instrument keys must use chromosome:position:REF:ALT", call. = FALSE)
  }
  chromosome <- vapply(fields, `[[`, character(1), 1L)
  position <- suppressWarnings(as.numeric(vapply(fields, `[[`, character(1), 2L)))
  ref <- vapply(fields, `[[`, character(1), 3L)
  alt <- vapply(fields, `[[`, character(1), 4L)
  normalized <- CompreSSoR::compressor_variant_key(chromosome, position, ref, alt)
  if (anyDuplicated(normalized)) {
    stop("instrument keys must not contain duplicates", call. = FALSE)
  }
  normalized
}

fastmr_normalize_instruments <- function(instruments, exposure_labels) {
  if (is.character(instruments)) {
    shared <- fastmr_normalize_variant_keys(instruments)
    out <- rep(list(shared), length(exposure_labels))
    names(out) <- exposure_labels
    return(out)
  }
  if (!is.list(instruments) || length(instruments) != length(exposure_labels)) {
    stop("instruments must be a canonical-key vector or one list element per exposure", call. = FALSE)
  }
  if (!is.null(names(instruments)) && all(nzchar(names(instruments)))) {
    missing <- setdiff(exposure_labels, names(instruments))
    if (length(missing)) {
      stop("named instruments are missing exposure(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
    instruments <- instruments[exposure_labels]
  } else if (!is.null(names(instruments)) && any(nzchar(names(instruments)))) {
    stop("instruments must be either fully named or completely unnamed", call. = FALSE)
  }
  out <- lapply(instruments, fastmr_normalize_variant_keys)
  names(out) <- exposure_labels
  out
}

fastmr_finalize_compressed_read <- function(out, columns) {
  identity <- c("chromosome", "base_pair_location", "effect_allele", "other_allele")
  if (!nrow(out)) {
    out$variant_key <- character()
  } else if (all(identity %in% names(out))) {
    out$variant_key <- CompreSSoR::compressor_variant_key(
      out$chromosome, out$base_pair_location, out$other_allele, out$effect_allele
    )
  }
  if (anyDuplicated(out$variant_key)) {
    stop("compressed store returned duplicate canonical variant keys", call. = FALSE)
  }
  keep <- unique(c(columns, if ("variant_key" %in% names(out)) "variant_key"))
  out[keep]
}

fastmr_io_map <- function(paths, keys, columns, io_threads) {
  exports <- getNamespaceExports("CompreSSoR")
  if ("read_sumstats_batch" %in% exports) {
    batch_reader <- getExportedValue("CompreSSoR", "read_sumstats_batch")
    identity <- c("chromosome", "base_pair_location", "effect_allele", "other_allele")
    requested <- unique(c(columns, identity))
    result <- tryCatch(
      batch_reader(
        unname(paths), unname(keys), columns = requested, threads = io_threads
      ),
      error = function(error) {
        stop("batched compressed read failed: ", conditionMessage(error), call. = FALSE)
      }
    )
    return(lapply(result, fastmr_finalize_compressed_read, columns = columns))
  }
  reader <- function(index) {
    tryCatch(
      fast_read_compressed(paths[[index]], variants = keys[[index]], columns = columns),
      error = function(error) {
        stop("failed to read ", paths[[index]], ": ", conditionMessage(error), call. = FALSE)
      }
    )
  }
  indexes <- seq_along(paths)
  if (.Platform$OS.type != "windows" && io_threads > 1L && length(indexes) > 1L) {
    result <- parallel::mclapply(
      indexes, reader, mc.cores = min(io_threads, length(indexes)),
      mc.preschedule = TRUE
    )
    failed <- vapply(result, inherits, logical(1), "try-error")
    if (any(failed)) {
      first <- which(failed)[1L]
      condition <- attr(result[[first]], "condition")
      detail <- if (inherits(condition, "condition")) {
        conditionMessage(condition)
      } else {
        as.character(result[[first]])
      }
      stop(detail, call. = FALSE)
    }
    return(result)
  }
  lapply(indexes, reader)
}

#' Read selected variants from a Pcodec CompreSSoR GWAS
#'
#' This is the FastMR-facing reader for a self-contained CompreSSoR store. It
#' accepts canonical `chromosome:position:REF:ALT` keys and adds those keys to
#' the returned data. No rsID or shared variant dictionary is required.
#'
#' @param path Pcodec CompreSSoR store directory.
#' @param variants Optional canonical keys or zero-based CompreSSoR row IDs.
#' @param columns Columns requested from CompreSSoR. Identity columns are added
#'   internally when needed to construct `variant_key`.
#' @return A data frame containing the requested columns and, when identity is
#'   available, `variant_key`.
#' @export
fast_read_compressed <- function(
    path,
    variants = NULL,
    columns = c("chromosome", "base_pair_location", "effect_allele",
                "other_allele", "beta", "standard_error")) {
  fastmr_require_compressor()
  if (length(path) != 1L || !is.character(path) || !dir.exists(path)) {
    stop("path must identify one existing CompreSSoR store", call. = FALSE)
  }
  store <- CompreSSoR::open_compressor(path)
  if (!identical(store$manifest$backend, "pcodec")) {
    stop("canonical-key FastMR input currently requires a Pcodec CompreSSoR store",
         call. = FALSE)
  }
  columns <- unique(as.character(columns))
  if (!length(columns) || anyNA(columns) || any(!nzchar(columns))) {
    stop("columns must contain at least one non-empty column name", call. = FALSE)
  }
  identity <- c("chromosome", "base_pair_location", "effect_allele", "other_allele")
  requested <- unique(c(columns, identity))
  out <- CompreSSoR::read_sumstats(store, variants = variants, columns = requested)
  fastmr_finalize_compressed_read(out, columns)
}

#' Run FastMR directly from compressed GWAS files
#'
#' Each exposure is read only for its supplied canonical instruments. Each
#' outcome is read once for the union of all instruments, after which FastMR
#' performs every exposure-outcome analysis. CompreSSoR effects are already
#' aligned to the ALT allele encoded in the canonical key, so matching keys are
#' directly comparable without rsID lookup or another allele-harmonisation
#' pass.
#'
#' @param exposure_files Named character vector of Pcodec CompreSSoR stores.
#' @param outcome_files Named character vector of Pcodec CompreSSoR stores.
#' @param instruments A canonical-key character vector shared by every exposure,
#'   or a named/positional list with one key vector per exposure.
#' @param methods FastMR method codes.
#' @param nboot Number of bootstrap draws.
#' @param seed Optional FastMR seed.
#' @param threads Native FastMR worker count.
#' @param io_threads Number of stores decoded concurrently inside the shared
#'   CompreSSoR process.
#' @param minimum_snps Minimum matched instruments required for every pair.
#' @param strict If `TRUE`, fail on invalid beta/standard-error values and when
#'   a pair has fewer than `minimum_snps`; otherwise omit them with warnings.
#' @param ... Additional options passed to [fast_mr()].
#' @return A tidy FastMR result with extraction metadata in the
#'   `compressed_input` attribute.
#' @export
fast_mr_compressed <- function(
    exposure_files,
    outcome_files,
    instruments,
    methods = "ivw",
    nboot = 0,
    seed = NULL,
    threads = 1,
    io_threads = 1,
    minimum_snps = 1L,
    strict = TRUE,
    ...) {
  fastmr_require_compressor()
  exposure_files <- fastmr_normalize_compressed_files(exposure_files, "exposure_files")
  outcome_files <- fastmr_normalize_compressed_files(outcome_files, "outcome_files")
  controls <- fastmr_validate_controls(nboot, seed, threads)
  io_threads <- fastmr_positive_integer_scalar(io_threads, "io_threads")
  minimum_snps <- fastmr_positive_integer_scalar(minimum_snps, "minimum_snps")
  if (length(strict) != 1L || !is.logical(strict) || is.na(strict)) {
    stop("strict must be TRUE or FALSE", call. = FALSE)
  }
  methods <- fastmr_normalize_methods(methods)
  instrument_sets <- fastmr_normalize_instruments(instruments, names(exposure_files))
  union_keys <- unique(unlist(instrument_sets, use.names = FALSE))
  columns <- c("chromosome", "base_pair_location", "effect_allele", "other_allele",
               "beta", "standard_error")
  outcome_keys <- rep(list(union_keys), length(outcome_files))
  all_data <- fastmr_io_map(
    c(unname(exposure_files), unname(outcome_files)),
    c(unname(instrument_sets), unname(outcome_keys)),
    columns,
    as.integer(io_threads)
  )
  exposure_count <- length(exposure_files)
  exposure_data <- all_data[seq_len(exposure_count)]
  outcome_data <- all_data[exposure_count + seq_along(outcome_files)]
  names(exposure_data) <- names(exposure_files)
  names(outcome_data) <- names(outcome_files)

  rows <- list()
  counts <- list()
  skipped <- character()
  invalid_omitted <- character()
  row_index <- 0L
  count_index <- 0L
  for (exposure_name in names(exposure_data)) {
    exposure <- exposure_data[[exposure_name]]
    wanted <- instrument_sets[[exposure_name]]
    exposure <- exposure[match(wanted, exposure$variant_key, nomatch = 0L), , drop = FALSE]
    exposure_valid <- is.finite(exposure$beta) & is.finite(exposure$standard_error) &
      exposure$standard_error > 0
    for (outcome_name in names(outcome_data)) {
      outcome <- outcome_data[[outcome_name]]
      matched <- match(exposure$variant_key, outcome$variant_key, nomatch = 0L)
      found <- matched > 0L
      outcome_valid <- logical(nrow(exposure))
      if (any(found)) {
        outcome_rows <- matched[found]
        outcome_valid[found] <- is.finite(outcome$beta[outcome_rows]) &
          is.finite(outcome$standard_error[outcome_rows]) &
          outcome$standard_error[outcome_rows] > 0
      }
      invalid_exposure <- sum(!exposure_valid)
      invalid_outcome <- sum(found & !outcome_valid)
      if (isTRUE(strict) && (invalid_exposure || invalid_outcome)) {
        stop(
          "invalid beta/standard_error for ", exposure_name, " -> ", outcome_name,
          " (exposure=", invalid_exposure, ", outcome=", invalid_outcome, ")",
          call. = FALSE
        )
      }
      if (!isTRUE(strict) && (invalid_exposure || invalid_outcome)) {
        invalid_omitted <- c(
          invalid_omitted,
          paste0(exposure_name, " -> ", outcome_name, " (exposure=",
                 invalid_exposure, ", outcome=", invalid_outcome, ")")
        )
      }
      keep <- found & exposure_valid & outcome_valid
      outcome_rows <- matched[keep]
      count_index <- count_index + 1L
      counts[[count_index]] <- data.frame(
        id.exposure = exposure_name, id.outcome = outcome_name,
        requested = length(wanted), exposure_found = nrow(exposure),
        outcome_found = sum(found), invalid_exposure = invalid_exposure,
        invalid_outcome = invalid_outcome, matched = sum(keep),
        stringsAsFactors = FALSE
      )
      if (sum(keep) < minimum_snps) {
        label <- paste0(exposure_name, " -> ", outcome_name, " (", sum(keep), " matched)")
        if (isTRUE(strict)) {
          stop("fewer than minimum_snps for ", label, call. = FALSE)
        }
        skipped <- c(skipped, label)
        next
      }
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        SNP = exposure$variant_key[keep],
        beta.exposure = exposure$beta[keep],
        beta.outcome = outcome$beta[outcome_rows],
        se.exposure = exposure$standard_error[keep],
        se.outcome = outcome$standard_error[outcome_rows],
        id.exposure = exposure_name,
        id.outcome = outcome_name,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) {
    stop("no exposure-outcome pair retained enough matched instruments", call. = FALSE)
  }
  if (length(skipped)) {
    warning("omitted pair(s) below minimum_snps: ", paste(skipped, collapse = "; "), call. = FALSE)
  }
  if (length(invalid_omitted)) {
    warning("omitted invalid instrument value(s): ",
            paste(unique(invalid_omitted), collapse = "; "), call. = FALSE)
  }
  result <- fast_mr(
    do.call(rbind, rows), methods = methods, nboot = controls$nboot,
    seed = controls$seed, threads = controls$threads, ...
  )
  attr(result, "compressed_input") <- list(
    exposure_files = exposure_files,
    outcome_files = outcome_files,
    instruments = instrument_sets,
    counts = do.call(rbind, counts),
    io_threads = as.integer(io_threads)
  )
  result
}
