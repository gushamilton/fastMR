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

fastmr_validate_compressed_store <- function(store) {
  manifest <- store$manifest
  identity <- manifest$identity
  compatible <- identical(manifest$backend, "pcodec") &&
    identical(manifest$genome_build, "GRCh38") &&
    identical(manifest$variant_storage, "self_contained_identity_key") &&
    is.list(identity) &&
    identical(identity$effect_allele_is_alt, TRUE) &&
    identical(identity$other_allele_is_ref, TRUE) &&
    identical(identity$external_reference_required, FALSE)
  if (!isTRUE(compatible)) {
    stop(
      "FastMR requires a self-contained GRCh38 Pcodec store whose effects are ALT-oriented",
      call. = FALSE
    )
  }
  invisible(store)
}

fastmr_finalize_compressed_read <- function(out, columns) {
  if (!is.data.frame(out)) {
    detail <- if (inherits(out, "condition")) conditionMessage(out) else as.character(out)
    stop("compressed reader returned an invalid result: ", detail, call. = FALSE)
  }
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
    failed <- !vapply(result, is.data.frame, logical(1))
    if (any(failed)) {
      first <- which(failed)[1L]
      failed_result <- result[[first]]
      condition <- attr(failed_result, "condition", exact = TRUE)
      detail <- if (inherits(condition, "condition")) {
        conditionMessage(condition)
      } else {
        as.character(failed_result)
      }
      stop(
        "batched compressed read failed for ", paths[[first]], ": ", detail,
        call. = FALSE
      )
    }
    source_bytes_read <- attr(result, "source_bytes_read", exact = TRUE)
    result <- lapply(result, fastmr_finalize_compressed_read, columns = columns)
    attr(result, "source_bytes_read") <- source_bytes_read
    return(result)
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

fastmr_compressed_grid_fast_path <- function(
    exposure_data, outcome_data, instrument_sets, methods, controls,
    minimum_snps, exposure_files, outcome_files, io_threads, dots) {
  bootstrap_methods <- c(
    "egger_bootstrap", "simple_median", "weighted_median",
    "penalised_weighted_median", "simple_mode", "weighted_mode"
  )
  # fast_mr() advances a supplied seed independently for every pair, whereas
  # fast_mr_grid() deliberately shares bootstrap layouts across its grid.  Use
  # the shortcut only when no requested result depends on bootstrap draws so
  # compressed input preserves the established seeded-result contract.
  if (controls$nboot > 0L && any(methods %in% bootstrap_methods)) return(NULL)
  wanted <- instrument_sets[[1L]]
  if (length(wanted) < minimum_snps ||
      !all(vapply(instrument_sets, identical, logical(1), wanted))) {
    return(NULL)
  }
  reorder_complete <- function(data) {
    matched <- match(wanted, data$variant_key)
    if (anyNA(matched)) return(NULL)
    out <- data[matched, c("beta", "standard_error"), drop = FALSE]
    valid <- is.finite(out$beta) & is.finite(out$standard_error) &
      out$standard_error > 0
    if (!all(valid)) return(NULL)
    out
  }
  exposures <- lapply(exposure_data, reorder_complete)
  outcomes <- lapply(outcome_data, reorder_complete)
  if (any(vapply(c(exposures, outcomes), is.null, logical(1)))) return(NULL)
  as_grid <- function(data, column, labels) {
    matrix <- do.call(rbind, lapply(data, `[[`, column))
    rownames(matrix) <- labels
    colnames(matrix) <- wanted
    matrix
  }
  exposure_labels <- names(exposure_data)
  outcome_labels <- names(outcome_data)
  call <- c(list(
    exposure_beta = as_grid(exposures, "beta", exposure_labels),
    outcome_beta = as_grid(outcomes, "beta", outcome_labels),
    exposure_se = as_grid(exposures, "standard_error", exposure_labels),
    outcome_se = as_grid(outcomes, "standard_error", outcome_labels),
    methods = methods, nboot = controls$nboot, seed = controls$seed,
    threads = controls$threads
  ), dots)
  result <- do.call(fast_mr_grid, call)
  result$exposure_index <- NULL
  result$outcome_index <- NULL
  pair <- seq_len(length(exposure_labels) * length(outcome_labels))
  exposure_index <- ((pair - 1L) %/% length(outcome_labels)) + 1L
  outcome_index <- ((pair - 1L) %% length(outcome_labels)) + 1L
  counts <- data.frame(
    id.exposure = exposure_labels[exposure_index],
    id.outcome = outcome_labels[outcome_index],
    requested = length(wanted), exposure_found = length(wanted),
    outcome_found = length(wanted), invalid_exposure = 0L,
    invalid_outcome = 0L, matched = length(wanted),
    stringsAsFactors = FALSE
  )
  attr(result, "compressed_input") <- list(
    exposure_files = exposure_files,
    outcome_files = outcome_files,
    instruments = instrument_sets,
    counts = counts,
    io_threads = as.integer(io_threads),
    estimator_path = "shared_instrument_grid"
  )
  result
}

# Sparse CSR dispatch for IVW-only compressed input.  Reproduces the pairwise
# loop's counts, strict-mode errors, warnings, minimum_snps drops and row order
# without building a data frame per exposure-outcome pair.  Returns NULL when
# the problem is outside the kernel's contract so the caller falls back to the
# pairwise path.
fastmr_compressed_sparse_ivw <- function(
    exposure_data, outcome_data, instrument_sets, union_keys, minimum_snps,
    strict, controls) {
  exposure_names <- names(exposure_data)
  outcome_names <- names(outcome_data)
  exposure_count <- length(exposure_names)
  outcome_count <- length(outcome_names)
  snp_count <- length(union_keys)
  if (!exposure_count || !outcome_count || !snp_count) return(NULL)
  bounds_ok <- tryCatch({
    fastmr_check_batch_bounds(
      exposure_count, outcome_count, snp_count, 1e8, 2048, sparse = TRUE
    )
    TRUE
  }, error = function(e) FALSE)
  if (!bounds_ok) return(NULL)
  valid_values <- function(data) {
    is.finite(data$beta) & is.finite(data$standard_error) &
      data$standard_error > 0
  }
  # Outcome union matrices (outcomes x union instruments).
  outcome_beta <- matrix(0, outcome_count, snp_count)
  outcome_se <- matrix(1, outcome_count, snp_count)
  outcome_found <- matrix(FALSE, outcome_count, snp_count)
  outcome_present <- matrix(FALSE, outcome_count, snp_count)
  for (o in seq_len(outcome_count)) {
    data <- outcome_data[[o]]
    hit <- match(union_keys, data$variant_key, nomatch = 0L)
    found <- hit > 0L
    outcome_found[o, ] <- found
    if (any(found)) {
      rows <- hit[found]
      valid <- valid_values(data)[rows]
      present <- logical(snp_count)
      present[found] <- valid
      outcome_present[o, ] <- present
      cols <- which(found)[valid]
      outcome_beta[o, cols] <- data$beta[rows[valid]]
      outcome_se[o, cols] <- data$standard_error[rows[valid]]
    }
  }
  # CSR over exposure-found and exposure-valid instruments, in instrument
  # order; per-pair outcome counts over all exposure-found instruments.
  requested <- integer(exposure_count)
  exposure_found <- integer(exposure_count)
  invalid_exposure <- integer(exposure_count)
  outcome_found_count <- matrix(0L, exposure_count, outcome_count)
  invalid_outcome <- matrix(0L, exposure_count, outcome_count)
  row_ptr <- integer(exposure_count + 1L)
  col_parts <- vector("list", exposure_count)
  beta_parts <- vector("list", exposure_count)
  for (e in seq_len(exposure_count)) {
    data <- exposure_data[[e]]
    wanted <- instrument_sets[[exposure_names[[e]]]]
    requested[[e]] <- length(wanted)
    hit <- match(wanted, data$variant_key, nomatch = 0L)
    present <- hit > 0L
    exposure_found[[e]] <- sum(present)
    cols <- match(wanted[present], union_keys)
    rows <- hit[present]
    valid <- valid_values(data)[rows]
    invalid_exposure[[e]] <- sum(!valid)
    if (length(cols)) {
      found_block <- outcome_found[, cols, drop = FALSE]
      present_block <- outcome_present[, cols, drop = FALSE]
      outcome_found_count[e, ] <- as.integer(rowSums(found_block))
      invalid_outcome[e, ] <- as.integer(rowSums(found_block & !present_block))
    }
    col_parts[[e]] <- cols[valid] - 1L
    beta_parts[[e]] <- data$beta[rows[valid]]
    row_ptr[[e + 1L]] <- row_ptr[[e]] + sum(valid)
  }
  col_index <- as.integer(unlist(col_parts, use.names = FALSE))
  if (!length(col_index)) return(NULL)
  sparse <- fast_mr_sparse_ivw(
    row_ptr, col_index, as.numeric(unlist(beta_parts, use.names = FALSE)),
    outcome_beta, outcome_se, outcome_present, threads = controls$threads
  )
  matched <- matrix(as.integer(sparse$nsnp), exposure_count, outcome_count)

  # Everything below is exposure-major (exposure outer, outcome inner).
  flat <- function(m) as.vector(t(m))
  exposure_of <- rep(seq_len(exposure_count), each = outcome_count)
  outcome_of <- rep(seq_len(outcome_count), times = exposure_count)
  pair_exposure <- exposure_names[exposure_of]
  pair_outcome <- outcome_names[outcome_of]
  pair_requested <- requested[exposure_of]
  pair_exposure_found <- exposure_found[exposure_of]
  pair_outcome_found <- flat(outcome_found_count)
  pair_invalid_exposure <- invalid_exposure[exposure_of]
  pair_invalid_outcome <- flat(invalid_outcome)
  pair_matched <- flat(matched)

  missing_exposure <- exposure_found < requested
  missing_outcome <- pair_outcome_found < pair_requested
  invalid_pair <- (pair_invalid_exposure > 0L) | (pair_invalid_outcome > 0L)
  too_few <- pair_matched < minimum_snps

  if (isTRUE(strict)) {
    pair_error <- missing_outcome | invalid_pair | too_few
    first_pair <- if (any(pair_error)) which(pair_error)[[1L]] else NA_integer_
    first_exposure <- if (any(missing_exposure)) which(missing_exposure)[[1L]] else NA_integer_
    if (!is.na(first_exposure) &&
        (is.na(first_pair) || first_exposure <= exposure_of[[first_pair]])) {
      stop(
        "missing requested exposure instrument(s) for ",
        exposure_names[[first_exposure]], " (found ",
        exposure_found[[first_exposure]], " of ", requested[[first_exposure]],
        ")", call. = FALSE
      )
    }
    if (!is.na(first_pair)) {
      i <- first_pair
      if (missing_outcome[[i]]) {
        stop(
          "missing requested outcome instrument(s) for ", pair_exposure[[i]],
          " -> ", pair_outcome[[i]], " (found ", pair_outcome_found[[i]],
          " of ", pair_requested[[i]], ")", call. = FALSE
        )
      }
      if (invalid_pair[[i]]) {
        stop(
          "invalid beta/standard_error for ", pair_exposure[[i]], " -> ",
          pair_outcome[[i]], " (exposure=", pair_invalid_exposure[[i]],
          ", outcome=", pair_invalid_outcome[[i]], ")", call. = FALSE
        )
      }
      stop("fewer than minimum_snps for ", pair_exposure[[i]], " -> ",
           pair_outcome[[i]], " (", pair_matched[[i]], " matched)", call. = FALSE)
    }
  }

  retained <- !too_few
  if (!any(retained)) {
    stop("no exposure-outcome pair retained enough matched instruments", call. = FALSE)
  }
  if (!isTRUE(strict)) {
    skipped <- if (any(too_few)) paste0(
      pair_exposure[too_few], " -> ", pair_outcome[too_few],
      " (", pair_matched[too_few], " matched)"
    ) else character()
    invalid_omitted <- if (any(invalid_pair)) paste0(
      pair_exposure[invalid_pair], " -> ", pair_outcome[invalid_pair],
      " (exposure=", pair_invalid_exposure[invalid_pair], ", outcome=",
      pair_invalid_outcome[invalid_pair], ")"
    ) else character()
    # Per exposure: the exposure message, then its outcome messages.
    messages <- matrix(NA_character_, outcome_count + 1L, exposure_count)
    if (any(missing_exposure)) {
      messages[1L, missing_exposure] <- paste0(
        exposure_names[missing_exposure], " exposure (found ",
        exposure_found[missing_exposure], " of ", requested[missing_exposure], ")"
      )
    }
    outcome_messages <- matrix(NA_character_, outcome_count, exposure_count)
    if (any(missing_outcome)) {
      outcome_messages[matrix(missing_outcome, outcome_count, exposure_count)] <- paste0(
        pair_exposure[missing_outcome], " -> ", pair_outcome[missing_outcome],
        " outcome (found ", pair_outcome_found[missing_outcome], " of ",
        pair_requested[missing_outcome], ")"
      )
    }
    messages[-1L, ] <- outcome_messages
    missing_omitted <- as.vector(messages)
    missing_omitted <- missing_omitted[!is.na(missing_omitted)]
    if (length(skipped)) {
      warning("omitted pair(s) below minimum_snps: ",
              paste(skipped, collapse = "; "), call. = FALSE)
    }
    if (length(invalid_omitted)) {
      warning("omitted invalid instrument value(s): ",
              paste(unique(invalid_omitted), collapse = "; "), call. = FALSE)
    }
    if (length(missing_omitted)) {
      warning("omitted missing requested instrument(s): ",
              paste(unique(missing_omitted), collapse = "; "), call. = FALSE)
    }
  }

  counts <- data.frame(
    id.exposure = pair_exposure, id.outcome = pair_outcome,
    requested = pair_requested, exposure_found = pair_exposure_found,
    outcome_found = pair_outcome_found,
    invalid_exposure = pair_invalid_exposure,
    invalid_outcome = pair_invalid_outcome, matched = pair_matched,
    stringsAsFactors = FALSE
  )

  n <- pair_matched[retained]
  b <- flat(sparse$beta)[retained]
  se <- flat(sparse$se)[retained]
  q <- flat(sparse$Q)[retained]
  sigma <- flat(sparse$sigma)[retained]
  has_fit <- is.finite(b)
  usable <- has_fit & is.finite(se) & se != 0
  pval <- rep(NA_real_, length(b))
  pval[usable] <- 2 * stats::pnorm(abs(b[usable] / se[usable]), lower.tail = FALSE)
  q_df <- ifelse(has_fit, n - 1, NA_real_)
  q_pval <- rep(NA_real_, length(b))
  ok_q <- has_fit & is.finite(q)
  q_pval[ok_q] <- stats::pchisq(q[ok_q], q_df[ok_q], lower.tail = FALSE)
  na <- rep(NA_real_, length(b))
  registry <- fastmr_method_registry()
  result <- data.frame(
    id.exposure = pair_exposure[retained],
    id.outcome = pair_outcome[retained],
    method = registry$method[match("ivw", registry$code)],
    method_code = "ivw",
    nsnp = as.numeric(n),
    b = ifelse(has_fit, b, NA_real_),
    se = ifelse(has_fit & is.finite(se), se, NA_real_),
    pval = pval,
    Q = ifelse(has_fit & is.finite(q), q, NA_real_),
    Q_df = q_df,
    Q_pval = q_pval,
    sigma = ifelse(has_fit & is.finite(sigma), sigma, NA_real_),
    intercept = na, intercept_se = na, intercept_pval = na,
    ratio_se_mean = na, bootstrap = na, phi = na, flipped = na,
    se_exposure_mean = na,
    stringsAsFactors = FALSE
  )
  list(result = result, counts = counts)
}

# Per-trait sample sizes for Steiger filtering of compressed input. Stores carry
# no sample size, so the caller supplies one per exposure/outcome label: a
# scalar shared by every trait, a vector named by label (extra names are
# ignored), or an unnamed vector in store order. Labels listed in `binary` take
# their effective N from case/control counts and may be omitted here.
fastmr_compressed_samplesize <- function(value, labels, binary, argument) {
  needed <- setdiff(labels, binary)
  if (is.null(value)) {
    if (length(needed)) {
      stop(argument, " is required when steiger = TRUE (stores carry no sample size)",
           call. = FALSE)
    }
    return(stats::setNames(rep(NA_real_, length(labels)), labels))
  }
  if (!is.numeric(value) || !length(value)) {
    stop(argument, " must be a numeric vector of sample sizes", call. = FALSE)
  }
  if (!is.null(names(value)) && all(nzchar(names(value)))) {
    if (anyDuplicated(names(value))) {
      stop(argument, " must not contain duplicate names", call. = FALSE)
    }
    missing <- setdiff(needed, names(value))
    if (length(missing)) {
      stop(argument, " is missing sample size(s) for: ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    out <- unname(as.numeric(value)[match(labels, names(value))])
  } else if (!is.null(names(value)) && any(nzchar(names(value)))) {
    stop(argument, " must be either fully named or completely unnamed", call. = FALSE)
  } else if (length(value) == 1L) {
    out <- rep(as.numeric(value), length(labels))
  } else if (length(value) == length(labels)) {
    out <- as.numeric(value)
  } else {
    stop(argument, " must have length one, one value per store, or names ",
         "matching the store labels", call. = FALSE)
  }
  names(out) <- labels
  bad <- labels %in% needed & !(is.finite(out) & out > 0)
  if (any(bad)) {
    stop(argument, " must be finite and positive for: ",
         paste(labels[bad], collapse = ", "), call. = FALSE)
  }
  out
}

fastmr_compressed_steiger_binary <- function(binary, labels) {
  if (is.null(binary)) return(NULL)
  required <- c("id", "ncase", "ncontrol", "prevalence")
  if (!is.data.frame(binary) || !all(required %in% names(binary))) {
    stop("steiger_binary must be a data frame with columns ",
         paste(required, collapse = ", "), call. = FALSE)
  }
  id <- as.character(binary$id)
  if (anyNA(id) || any(!nzchar(id)) || anyDuplicated(id)) {
    stop("steiger_binary$id must contain unique non-empty trait labels", call. = FALSE)
  }
  unknown <- setdiff(id, labels)
  if (length(unknown)) {
    stop("steiger_binary$id contains label(s) that are not exposure or outcome stores: ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }
  out <- data.frame(
    id = id,
    ncase = fastmr_numeric(binary$ncase, "steiger_binary$ncase"),
    ncontrol = fastmr_numeric(binary$ncontrol, "steiger_binary$ncontrol"),
    prevalence = fastmr_numeric(binary$prevalence, "steiger_binary$prevalence"),
    stringsAsFactors = FALSE
  )
  bad <- !(is.finite(out$ncase) & out$ncase > 0 &
             is.finite(out$ncontrol) & out$ncontrol > 0 &
             is.finite(out$prevalence) & out$prevalence > 0 & out$prevalence < 1)
  if (any(bad)) {
    stop("steiger_binary needs positive ncase/ncontrol and 0 < prevalence < 1 for: ",
         paste(out$id[bad], collapse = ", "), call. = FALSE)
  }
  out
}

# Validates the Steiger options of fast_mr_compressed(); NULL when disabled.
fastmr_compressed_steiger_options <- function(steiger, samplesize_exposure,
                                              samplesize_outcome, steiger_binary,
                                              exposure_labels, outcome_labels) {
  if (length(steiger) != 1L || !is.logical(steiger) || is.na(steiger)) {
    stop("steiger must be TRUE or FALSE", call. = FALSE)
  }
  if (!steiger) {
    supplied <- c(
      samplesize_exposure = !is.null(samplesize_exposure),
      samplesize_outcome = !is.null(samplesize_outcome),
      steiger_binary = !is.null(steiger_binary)
    )
    if (any(supplied)) {
      stop(paste(names(supplied)[supplied], collapse = ", "),
           " only applies when steiger = TRUE", call. = FALSE)
    }
    return(NULL)
  }
  binary <- fastmr_compressed_steiger_binary(
    steiger_binary, union(exposure_labels, outcome_labels)
  )
  binary_ids <- if (is.null(binary)) character() else binary$id
  list(
    samplesize_exposure = fastmr_compressed_samplesize(
      samplesize_exposure, exposure_labels, binary_ids, "samplesize_exposure"
    ),
    samplesize_outcome = fastmr_compressed_samplesize(
      samplesize_outcome, outcome_labels, binary_ids, "samplesize_outcome"
    ),
    binary = binary
  )
}

# Steiger filtering on the rows already read for MR. A pair's rows are its
# exposure instruments found in both stores with finite beta and positive
# standard error (the rows MR uses), in instrument order; pairs with fewer than
# `minimum_snps` such rows are dropped, as MR drops them. Rows are
# exposure-major, outcomes in store order.
fastmr_compressed_steiger <- function(exposure_data, outcome_data, instrument_sets,
                                      minimum_snps, options) {
  valid_values <- function(data) {
    is.finite(data$beta) & is.finite(data$standard_error) & data$standard_error > 0
  }
  outcome_valid <- lapply(outcome_data, valid_values)
  parts <- list()
  for (exposure_name in names(exposure_data)) {
    exposure <- exposure_data[[exposure_name]]
    rows <- match(instrument_sets[[exposure_name]], exposure$variant_key, nomatch = 0L)
    rows <- rows[rows > 0L]
    rows <- rows[valid_values(exposure)[rows]]
    keys <- exposure$variant_key[rows]
    for (outcome_name in names(outcome_data)) {
      hit <- match(keys, outcome_data[[outcome_name]]$variant_key, nomatch = 0L)
      keep <- hit > 0L
      keep[keep] <- outcome_valid[[outcome_name]][hit[keep]]
      if (sum(keep) < minimum_snps) next
      parts[[length(parts) + 1L]] <- list(
        exposure = exposure_name, outcome = outcome_name,
        exposure_rows = rows[keep], outcome_rows = hit[keep]
      )
    }
  }
  if (!length(parts)) return(data.frame())
  size <- vapply(parts, function(part) length(part$exposure_rows), integer(1))
  id_exposure <- rep(vapply(parts, `[[`, character(1), "exposure"), size)
  id_outcome <- rep(vapply(parts, `[[`, character(1), "outcome"), size)
  gather <- function(side, column) {
    source <- if (side == "exposure") exposure_data else outcome_data
    rows_name <- paste0(side, "_rows")
    unlist(lapply(parts, function(part) {
      source[[part[[side]]]][[column]][part[[rows_name]]]
    }), use.names = FALSE)
  }
  data <- data.frame(
    SNP = gather("exposure", "variant_key"),
    id.exposure = id_exposure,
    id.outcome = id_outcome,
    exposure = id_exposure,
    outcome = id_outcome,
    beta.exposure = gather("exposure", "beta"),
    beta.outcome = gather("outcome", "beta"),
    se.exposure = gather("exposure", "standard_error"),
    se.outcome = gather("outcome", "standard_error"),
    eaf.exposure = gather("exposure", "effect_allele_frequency"),
    eaf.outcome = gather("outcome", "effect_allele_frequency"),
    pval.exposure = gather("exposure", "p_value"),
    pval.outcome = gather("outcome", "p_value"),
    samplesize.exposure = unname(options$samplesize_exposure[id_exposure]),
    samplesize.outcome = unname(options$samplesize_outcome[id_outcome]),
    units.exposure = "",
    units.outcome = "",
    mr_keep = TRUE,
    stringsAsFactors = FALSE
  )
  binary <- options$binary
  if (!is.null(binary)) {
    for (side in c("exposure", "outcome")) {
      index <- match(data[[paste0("id.", side)]], binary$id)
      data[[paste0("units.", side)]][!is.na(index)] <- "log odds"
      data[[paste0("ncase.", side)]] <- binary$ncase[index]
      data[[paste0("ncontrol.", side)]] <- binary$ncontrol[index]
      data[[paste0("prevalence.", side)]] <- binary$prevalence[index]
    }
  }
  fast_mr_steiger_filtering(data)
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
  fastmr_validate_compressed_store(store)
  columns <- unique(as.character(columns))
  if (!length(columns) || anyNA(columns) || any(!nzchar(columns))) {
    stop("columns must contain at least one non-empty column name", call. = FALSE)
  }
  identity <- c("chromosome", "base_pair_location", "effect_allele", "other_allele")
  requested <- unique(c(columns, identity))
  out <- CompreSSoR::read_sumstats(store, variants = variants, columns = requested)
  source_bytes_read <- attr(out, "source_bytes_read", exact = TRUE)
  out <- fastmr_finalize_compressed_read(out, columns)
  attr(out, "source_bytes_read") <- source_bytes_read
  out
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
#' @param strict If `TRUE`, fail when any requested instrument is missing, on
#'   invalid beta/standard-error values, and when a pair has fewer than
#'   `minimum_snps`; otherwise omit unavailable or invalid rows with warnings.
#' @param output Optional path for a Zstandard-compressed Parquet copy of the
#'   result. The path must not already exist; use [fast_write_parquet()] when
#'   an overwrite or another compression codec is required.
#' @param estimator `"auto"` (default) or `"pairwise"`. With `"auto"`, runs
#'   requesting only `methods = "ivw"` (and no extra `...` options) whose
#'   instrument sets are not shared by all exposures use the sparse CSR kernel
#'   [fast_mr_sparse_ivw()] instead of an R loop calling [fast_mr()] per pair.
#'   Counts, errors, warnings, `minimum_snps` handling and row order are
#'   identical, but estimates may differ from `"pairwise"` by rounding
#'   (agreement is within 1e-14 relative on beta, se and Q; the two kernels use
#'   the same two-pass residual formula and summation order, so results are
#'   usually bit-identical). `"pairwise"`
#'   (and every other method set) always uses the per-pair path. The path used
#'   is reported as `estimator_path` in the `compressed_input` attribute
#'   (`"sparse_ivw"`, `"pairwise"` or `"shared_instrument_grid"`).
#' @param steiger If `TRUE`, also run per-SNP Steiger filtering
#'   ([fast_mr_steiger_filtering()]) for every retained exposure-outcome pair.
#'   `effect_allele_frequency` and `p_value` are read in the same pass as
#'   beta/standard error, so no store is read twice. The result is returned in
#'   the `steiger` attribute; the MR result itself is unchanged. Steiger rows
#'   are the rows MR uses: instruments found in both stores with finite beta
#'   and positive standard error, for pairs meeting `minimum_snps`.
#' @param samplesize_exposure,samplesize_outcome Sample sizes for Steiger
#'   (stores carry none): a vector named by store label (as in
#'   `exposure_files`/`outcome_files`), an unnamed vector in store order, or
#'   one value for every store. Required when `steiger = TRUE`, except for
#'   traits listed in `steiger_binary`.
#' @param steiger_binary Optional data frame for binary (log-odds) traits with
#'   columns `id` (an exposure or outcome store label), `ncase`, `ncontrol` and
#'   `prevalence`. Those traits use the log-odds R-squared model with the
#'   harmonic effective sample size, as in [fast_mr_steiger_filtering()]; a
#'   label shared by an exposure and an outcome applies to both.
#' @param ... Additional options passed to [fast_mr()].
#' @return A tidy FastMR result with extraction metadata in the
#'   `compressed_input` attribute. With `steiger = TRUE`, the `steiger`
#'   attribute holds the per-SNP [fast_mr_steiger_filtering()] output
#'   (exposure-major, instruments in supplied order) and
#'   `compressed_input$timing$steiger_seconds` its compute time. The
#'   `steiger` attribute is not written to `output`.
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
    output = NULL,
    estimator = c("auto", "pairwise"),
    steiger = FALSE,
    samplesize_exposure = NULL,
    samplesize_outcome = NULL,
    steiger_binary = NULL,
    ...) {
  estimator <- match.arg(estimator)
  total_started <- unname(proc.time()[["elapsed"]])
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
  steiger_options <- fastmr_compressed_steiger_options(
    steiger, samplesize_exposure, samplesize_outcome, steiger_binary,
    names(exposure_files), names(outcome_files)
  )
  dots <- list(...)
  instrument_sets <- fastmr_normalize_instruments(instruments, names(exposure_files))
  union_keys <- unique(unlist(instrument_sets, use.names = FALSE))
  columns <- c("chromosome", "base_pair_location", "effect_allele", "other_allele",
               "beta", "standard_error")
  if (!is.null(steiger_options)) {
    columns <- c(columns, "effect_allele_frequency", "p_value")
  }
  # Attach Steiger (when requested) after estimation, so it never runs on a
  # study that strict mode or minimum_snps rejected.
  # The Parquet copy is written before the `steiger` attribute is attached.
  finish <- function(result) {
    if (is.null(steiger_options)) return(fastmr_write_result(result, output))
    steiger_started <- unname(proc.time()[["elapsed"]])
    steiger_result <- fastmr_compressed_steiger(
      exposure_data, outcome_data, instrument_sets, minimum_snps, steiger_options
    )
    metadata <- attr(result, "compressed_input")
    metadata$timing$steiger_seconds <- unname(proc.time()[["elapsed"]]) - steiger_started
    attr(result, "compressed_input") <- metadata
    result <- fastmr_write_result(result, output)
    attr(result, "steiger") <- steiger_result
    result
  }
  outcome_keys <- rep(list(union_keys), length(outcome_files))
  io_started <- unname(proc.time()[["elapsed"]])
  invisible(lapply(
    unique(c(unname(exposure_files), unname(outcome_files))),
    function(path) fastmr_validate_compressed_store(
      CompreSSoR::open_compressor(path)
    )
  ))
  all_data <- fastmr_io_map(
    c(unname(exposure_files), unname(outcome_files)),
    c(unname(instrument_sets), unname(outcome_keys)),
    columns,
    as.integer(io_threads)
  )
  source_bytes_read <- attr(all_data, "source_bytes_read", exact = TRUE)
  source_bytes_read <- if (is.null(source_bytes_read)) {
    NA_real_
  } else {
    as.numeric(source_bytes_read)
  }
  io_seconds <- unname(proc.time()[["elapsed"]]) - io_started
  exposure_count <- length(exposure_files)
  exposure_data <- all_data[seq_len(exposure_count)]
  outcome_data <- all_data[exposure_count + seq_along(outcome_files)]
  names(exposure_data) <- names(exposure_files)
  names(outcome_data) <- names(outcome_files)

  estimator_started <- unname(proc.time()[["elapsed"]])
  grid_result <- fastmr_compressed_grid_fast_path(
    exposure_data, outcome_data, instrument_sets, methods, controls,
    minimum_snps, exposure_files, outcome_files, io_threads, dots
  )
  if (!is.null(grid_result)) {
    metadata <- attr(grid_result, "compressed_input")
    metadata$timing <- list(
      io_seconds = io_seconds,
      estimator_seconds = unname(proc.time()[["elapsed"]]) - estimator_started,
      total_seconds = unname(proc.time()[["elapsed"]]) - total_started,
      source_bytes_read = source_bytes_read
    )
    attr(grid_result, "compressed_input") <- metadata
    return(finish(grid_result))
  }

  if (identical(estimator, "auto") && identical(methods, "ivw") && !length(dots)) {
    sparse <- fastmr_compressed_sparse_ivw(
      exposure_data, outcome_data, instrument_sets, union_keys, minimum_snps,
      strict, controls
    )
    if (!is.null(sparse)) {
      result <- sparse$result
      attr(result, "compressed_input") <- list(
        exposure_files = exposure_files,
        outcome_files = outcome_files,
        instruments = instrument_sets,
        counts = sparse$counts,
        io_threads = as.integer(io_threads),
        estimator_path = "sparse_ivw",
        timing = list(
          io_seconds = io_seconds,
          estimator_seconds = unname(proc.time()[["elapsed"]]) - estimator_started,
          total_seconds = unname(proc.time()[["elapsed"]]) - total_started,
          source_bytes_read = source_bytes_read
        )
      )
      return(finish(result))
    }
  }

  rows <- list()
  counts <- list()
  skipped <- character()
  invalid_omitted <- character()
  missing_omitted <- character()
  row_index <- 0L
  count_index <- 0L
  for (exposure_name in names(exposure_data)) {
    exposure <- exposure_data[[exposure_name]]
    wanted <- instrument_sets[[exposure_name]]
    exposure_match <- match(wanted, exposure$variant_key, nomatch = 0L)
    exposure_found <- sum(exposure_match > 0L)
    if (isTRUE(strict) && exposure_found < length(wanted)) {
      stop(
        "missing requested exposure instrument(s) for ", exposure_name,
        " (found ", exposure_found, " of ", length(wanted), ")",
        call. = FALSE
      )
    }
    if (!isTRUE(strict) && exposure_found < length(wanted)) {
      missing_omitted <- c(
        missing_omitted,
        paste0(exposure_name, " exposure (found ", exposure_found, " of ",
               length(wanted), ")")
      )
    }
    exposure <- exposure[exposure_match[exposure_match > 0L], , drop = FALSE]
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
      outcome_found <- sum(found)
      if (isTRUE(strict) && outcome_found < length(wanted)) {
        stop(
          "missing requested outcome instrument(s) for ", exposure_name,
          " -> ", outcome_name, " (found ", outcome_found, " of ",
          length(wanted), ")", call. = FALSE
        )
      }
      if (!isTRUE(strict) && outcome_found < length(wanted)) {
        missing_omitted <- c(
          missing_omitted,
          paste0(exposure_name, " -> ", outcome_name, " outcome (found ",
                 outcome_found, " of ", length(wanted), ")")
        )
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
        requested = length(wanted), exposure_found = exposure_found,
        outcome_found = outcome_found, invalid_exposure = invalid_exposure,
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
  if (length(missing_omitted)) {
    warning(
      "omitted missing requested instrument(s): ",
      paste(unique(missing_omitted), collapse = "; "), call. = FALSE
    )
  }
  result <- do.call(fast_mr, c(list(
    data = do.call(rbind, rows), methods = methods, nboot = controls$nboot,
    seed = controls$seed, threads = controls$threads
  ), dots))
  attr(result, "compressed_input") <- list(
    exposure_files = exposure_files,
    outcome_files = outcome_files,
    instruments = instrument_sets,
    counts = do.call(rbind, counts),
    io_threads = as.integer(io_threads),
    estimator_path = "pairwise",
    timing = list(
      io_seconds = io_seconds,
      estimator_seconds = unname(proc.time()[["elapsed"]]) - estimator_started,
      total_seconds = unname(proc.time()[["elapsed"]]) - total_started,
      source_bytes_read = source_bytes_read
    )
  )
  finish(result)
}
