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
  # An exposure with no instruments is kept as character() here and handled
  # by fastmr_compressed_nonempty_instruments() according to `strict`.
  out <- lapply(instruments, function(keys) {
    if (is.null(keys) || (is.character(keys) && !length(keys))) character()
    else fastmr_normalize_variant_keys(keys)
  })
  names(out) <- exposure_labels
  out
}

# Exposures whose instrument set is non-empty.  An empty set is an error with
# strict = TRUE; with strict = FALSE the exposure is dropped with a warning
# (an error only when no exposure has instruments).
fastmr_compressed_nonempty_instruments <- function(instrument_sets, strict) {
  empty <- lengths(instrument_sets) == 0L
  if (!any(empty)) return(names(instrument_sets))
  labels <- names(instrument_sets)[empty]
  shown <- paste(utils::head(labels, 20L), collapse = ", ")
  if (length(labels) > 20L) shown <- paste0(shown, ", ... (", length(labels), " in total)")
  if (all(empty)) stop("no exposure has any instruments", call. = FALSE)
  if (strict) {
    stop("instrument set is empty for exposure(s): ", shown,
         "; use strict = FALSE to drop them", call. = FALSE)
  }
  warning("dropping exposure(s) with an empty instrument set: ", shown, call. = FALSE)
  names(instrument_sets)[!empty]
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

# Numeric variant identity from a store manifest.  Self-contained Pcodec
# stores document their identity encoding (schema
# compressor_variant_identity_v1): a zero-based global position (chromosome
# offset + position - 1) and a substitution code 4 * REF + ALT over A, C, G, T,
# combined as global_position * 16 + substitution.  Returns NULL when the
# manifest does not declare exactly that encoding.
fastmr_compressed_identity_codec <- function(manifest) {
  identity <- manifest$identity
  table <- identity$chromosome_table
  chromosomes <- as.character(unlist(table$chromosomes))
  offsets <- as.numeric(unlist(identity$chromosome_offsets))
  lengths <- as.numeric(unlist(identity$chromosome_lengths))
  ok <- is.list(identity) &&
    identical(identity$schema, "compressor_variant_identity_v1") &&
    identical(identity$encoding, "global_position_plus_directed_ref_alt_substitution") &&
    identical(identity$position_encoding, "zero_based_global_position") &&
    identical(identity$substitution_encoding, "uint8_4_times_ref_plus_alt") &&
    length(chromosomes) > 0L && length(chromosomes) == length(offsets) &&
    length(lengths) == length(offsets) && !anyNA(offsets) && !anyNA(lengths) &&
    !anyDuplicated(chromosomes) && !is.unsorted(offsets, strictly = TRUE)
  if (!isTRUE(ok)) return(NULL)
  list(chromosomes = chromosomes, offsets = offsets, lengths = lengths)
}

fastmr_compressed_identity_bases <- c("A", "C", "G", "T")

# Identity codes of canonical keys; NA for any key that is not a canonical
# single-nucleotide key within the codec's chromosome table (such a key can never
# be paired with a row by code, and is resolved by the string fallback).
fastmr_compressed_key_codes <- function(keys, codec) {
  pattern <- "^([0-9A-Za-z]+):([1-9][0-9]*):([ACGT]):([ACGT])$"
  canonical <- grepl(pattern, keys, perl = TRUE)
  code <- rep(NA_real_, length(keys))
  if (!any(canonical)) return(code)
  k <- keys[canonical]
  chromosome <- match(sub(pattern, "\\1", k, perl = TRUE), codec$chromosomes)
  position <- as.numeric(sub(pattern, "\\2", k, perl = TRUE))
  ref <- match(sub(pattern, "\\3", k, perl = TRUE), fastmr_compressed_identity_bases) - 1L
  alt <- match(sub(pattern, "\\4", k, perl = TRUE), fastmr_compressed_identity_bases) - 1L
  global <- codec$offsets[chromosome] + position - 1
  in_range <- !is.na(chromosome) & position <= codec$lengths[chromosome] & ref != alt
  code[canonical] <- ifelse(in_range, global * 16 + 4 * ref + alt, NA_real_)
  code
}

# Canonical keys of rows from their global positions and substitution codes
# (the fallback for rows no requested key's code accounts for).
fastmr_compressed_decode_keys <- function(global_position, substitution, codec) {
  chromosome <- findInterval(global_position, codec$offsets)
  position <- global_position - codec$offsets[chromosome] + 1
  CompreSSoR::compressor_variant_key(
    codec$chromosomes[chromosome], position,
    fastmr_compressed_identity_bases[substitution %/% 4L + 1L],
    fastmr_compressed_identity_bases[substitution %% 4L + 1L]
  )
}

# Opens and validates every store (io_threads at a time, in forked workers on
# Unix) and returns each store's identity codec (NULL when its manifest does
# not declare the standard encoding).  The first failing store, in `paths`
# order, raises the error the serial open/validate loop raised.
fastmr_compressed_validate_stores <- function(paths, io_threads = 1L) {
  check <- function(path) {
    tryCatch({
      store <- CompreSSoR::open_compressor(path)
      fastmr_validate_compressed_store(store)
      list(ok = TRUE, codec = fastmr_compressed_identity_codec(store$manifest))
    }, error = function(error) list(ok = FALSE, message = conditionMessage(error)))
  }
  if (.Platform$OS.type != "windows" && io_threads > 1L && length(paths) > 1L) {
    result <- parallel::mclapply(paths, check, mc.cores = min(io_threads, length(paths)),
                                 mc.preschedule = TRUE)
  } else {
    result <- lapply(paths, check)
  }
  for (i in seq_along(result)) {
    r <- result[[i]]
    if (inherits(r, "try-error")) stop(as.character(r), call. = FALSE)
    if (!isTRUE(r$ok)) stop(r$message, call. = FALSE)
  }
  stats::setNames(lapply(result, `[[`, "codec"), paths)
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

# Coded counterpart of fastmr_finalize_compressed_read(): adds `variant_key`
# from identity codes.  Codes of each distinct request (one object shared by
# every outcome store, for example) under each distinct codec are computed
# once.
fastmr_finalize_coded_reads <- function(result, keys, codecs, columns) {
  key_group <- fastmr_request_groups(keys)
  codec_group <- fastmr_request_groups(codecs)
  cache <- new.env(hash = TRUE, parent = emptyenv())
  request_codes <- function(i) {
    slot <- paste(key_group[i], codec_group[i])
    codes <- cache[[slot]]
    if (is.null(codes)) {
      codes <- fastmr_compressed_key_codes(as.character(keys[[i]]), codecs[[i]])
      cache[[slot]] <- codes
    }
    codes
  }
  lapply(seq_along(result), function(i) {
    out <- result[[i]]
    if (!is.data.frame(out)) {
      detail <- if (inherits(out, "condition")) conditionMessage(out) else as.character(out)
      stop("compressed reader returned an invalid result: ", detail, call. = FALSE)
    }
    if (!nrow(out)) {
      out$variant_key <- character()
    } else {
      global <- as.numeric(out$global_position)
      substitution <- as.integer(out$substitution)
      hit <- match(global * 16 + substitution, request_codes(i))
      variant_key <- as.character(keys[[i]])[hit]
      missing <- is.na(hit)
      if (any(missing)) {
        variant_key[missing] <- fastmr_compressed_decode_keys(
          global[missing], substitution[missing], codecs[[i]]
        )
      }
      out$variant_key <- variant_key
    }
    if (anyDuplicated(out$variant_key)) {
      stop("compressed store returned duplicate canonical variant keys", call. = FALSE)
    }
    out[unique(c(columns, "variant_key"))]
  })
}

# Whether the installed CompreSSoR returns numeric identity columns from key
# reads: NULL until known, then TRUE/FALSE for this session and build.
.fastmr_compressed_state <- new.env(parent = emptyenv())
fastmr_coded_reads_supported <- function(value) {
  build <- paste(utils::packageVersion("CompreSSoR"),
                 find.package("CompreSSoR", quiet = TRUE))
  if (!missing(value)) {
    .fastmr_compressed_state$coded_reads <- list(build = build, value = value)
    return(invisible(value))
  }
  known <- .fastmr_compressed_state$coded_reads
  if (is.null(known) || !identical(known$build, build)) NULL else known$value
}

# Whether the installed CompreSSoR reports a capability.
fastmr_compressor_has <- function(capability) {
  exports <- getNamespaceExports("CompreSSoR")
  if (!"compressor_capabilities" %in% exports) return(FALSE)
  capabilities <- tryCatch(getExportedValue("CompreSSoR", "compressor_capabilities")(),
                           error = function(error) character())
  capability %in% capabilities
}

# The request-index path applies when CompreSSoR offers it, every store
# declares the standard identity encoding (`codecs` from
# fastmr_compressed_validate_stores()), and every requested key is a strictly
# canonical single-nucleotide key, so the requested string is exactly the key
# the string path would rebuild for its row.
fastmr_request_index_usable <- function(paths, keys, codecs) {
  if (is.null(codecs) || length(codecs) != length(paths) ||
      any(vapply(codecs, is.null, logical(1))) ||
      !fastmr_compressor_has("request_index")) {
    return(FALSE)
  }
  pattern <- "^([1-9]|1[0-9]|2[0-2]|X|Y):[1-9][0-9]*:[ACGT]:[ACGT]$"
  # Check each distinct request once (the shared outcome request is one).
  group <- fastmr_request_groups(keys)
  for (i in which(group == seq_along(group))) {
    request <- keys[[i]]
    if (!is.character(request)) return(FALSE)
    if (!all(grepl(pattern, request, perl = TRUE))) return(FALSE)
  }
  TRUE
}

# Batched read with CompreSSoR's request index: each row's variant key is the
# requested key it answers, with no identity column read.  NULL when any row
# lacks an index (the caller then uses the manifest-decoding path).
fastmr_io_map_indexed <- function(batch_reader, paths, keys, columns, io_threads) {
  result <- tryCatch(
    batch_reader(unname(paths), unname(keys), columns = unique(columns),
                 threads = io_threads, request_index = TRUE),
    error = function(error) {
      stop("batched compressed read failed: ", conditionMessage(error), call. = FALSE)
    }
  )
  failed <- !vapply(result, is.data.frame, logical(1))
  if (any(failed)) {
    first <- which(failed)[1L]
    condition <- attr(result[[first]], "condition", exact = TRUE)
    detail <- if (inherits(condition, "condition")) conditionMessage(condition) else
      as.character(result[[first]])
    stop("batched compressed read failed for ", paths[[first]], ": ", detail, call. = FALSE)
  }
  source_bytes_read <- attr(result, "source_bytes_read", exact = TRUE)
  out <- vector("list", length(result))
  for (i in seq_along(result)) {
    data <- result[[i]]
    index <- data$request_index
    if (is.null(index) || anyNA(index)) return(NULL)
    data$request_index <- NULL
    data$variant_key <- as.character(keys[[i]])[index]
    if (anyDuplicated(data$variant_key)) {
      stop("compressed store returned duplicate canonical variant keys", call. = FALSE)
    }
    out[[i]] <- data[unique(c(columns, "variant_key"))]
  }
  attr(out, "source_bytes_read") <- source_bytes_read
  out
}

# Reads `keys[[i]]` from `paths[[i]]`.  `codecs` (from
# fastmr_compressed_validate_stores()) enables the numeric identity path: the
# reader returns global positions and substitution codes instead of decoded
# chromosome/allele strings, and each row takes the variant key of the
# requested key with the same identity code, so no key string is rebuilt per
# row.  Rows no requested code accounts for (non-canonical requests) are
# decoded to their canonical key as before.
fastmr_io_map <- function(paths, keys, columns, io_threads, codecs = NULL,
                          use_request_index = TRUE) {
  exports <- getNamespaceExports("CompreSSoR")
  if ("read_sumstats_batch" %in% exports) {
    batch_reader <- getExportedValue("CompreSSoR", "read_sumstats_batch")
    if (isTRUE(use_request_index) && fastmr_request_index_usable(paths, keys, codecs)) {
      indexed <- fastmr_io_map_indexed(batch_reader, paths, keys, columns, io_threads)
      if (!is.null(indexed)) return(indexed)
      return(fastmr_io_map(paths, keys, columns, io_threads, codecs = codecs,
                           use_request_index = FALSE))
    }
    identity <- c("chromosome", "base_pair_location", "effect_allele", "other_allele")
    numeric_identity <- !is.null(codecs) && length(codecs) == length(paths) &&
      !any(vapply(codecs, is.null, logical(1))) &&
      !identical(fastmr_coded_reads_supported(), FALSE)
    requested <- if (numeric_identity) {
      unique(c(columns, "global_position", "substitution"))
    } else {
      unique(c(columns, identity))
    }
    # Older CompreSSoR builds do not return global_position/substitution from
    # key reads; the first such refusal switches this session to the string
    # path and the read is repeated there.
    unsupported <- function(message) {
      numeric_identity && grepl("requested columns are not present", message, fixed = TRUE)
    }
    retry <- FALSE
    result <- tryCatch(
      batch_reader(
        unname(paths), unname(keys), columns = requested, threads = io_threads
      ),
      error = function(error) {
        if (unsupported(conditionMessage(error))) {
          retry <<- TRUE
          return(NULL)
        }
        stop("batched compressed read failed: ", conditionMessage(error), call. = FALSE)
      }
    )
    if (!retry && numeric_identity) {
      details <- vapply(result, function(r) {
        if (is.data.frame(r)) return("")
        condition <- attr(r, "condition", exact = TRUE)
        if (inherits(condition, "condition")) conditionMessage(condition) else
          paste(as.character(r), collapse = " ")
      }, character(1))
      retry <- any(vapply(details, unsupported, logical(1)))
    }
    if (retry) {
      fastmr_coded_reads_supported(FALSE)
      return(fastmr_io_map(paths, keys, columns, io_threads, codecs = NULL,
                           use_request_index = FALSE))
    }
    if (numeric_identity) fastmr_coded_reads_supported(TRUE)
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
    result <- if (numeric_identity) {
      fastmr_finalize_coded_reads(result, keys, codecs, columns)
    } else {
      lapply(result, fastmr_finalize_compressed_read, columns = columns)
    }
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

fastmr_compressed_valid_values <- function(data) {
  is.finite(data$beta) & is.finite(data$standard_error) & data$standard_error > 0
}

# Vectorised exposure-outcome pair index for compressed input.  An "instrument"
# is one requested exposure key found in that exposure's store, in supplied
# order; an "entry" is an instrument kept for one outcome (found in the outcome
# store, finite beta and positive standard error on both sides).  Per-pair
# counts are exposure-major (exposure outer, outcome inner) and entries are
# ordered exposure-major, then by instrument, which is exactly the row order of
# the former per-pair loop.  Work and memory are linear in
# instruments x outcomes (one match() per store, no per-pair R objects).
fastmr_compressed_pair_index <- function(exposure_data, outcome_data, instrument_sets) {
  exposure_names <- names(exposure_data)
  outcome_names <- names(outcome_data)
  exposure_count <- length(exposure_names)
  outcome_count <- length(outcome_names)
  requested <- integer(exposure_count)
  exposure_found <- integer(exposure_count)
  invalid_exposure <- integer(exposure_count)
  rows_list <- vector("list", exposure_count)
  valid_list <- vector("list", exposure_count)
  key_list <- vector("list", exposure_count)
  for (e in seq_len(exposure_count)) {
    data <- exposure_data[[e]]
    wanted <- instrument_sets[[exposure_names[[e]]]]
    hit <- match(wanted, data$variant_key, nomatch = 0L)
    rows <- hit[hit > 0L]
    valid <- fastmr_compressed_valid_values(data)[rows]
    requested[[e]] <- length(wanted)
    exposure_found[[e]] <- length(rows)
    invalid_exposure[[e]] <- sum(!valid)
    rows_list[[e]] <- rows
    valid_list[[e]] <- valid
    key_list[[e]] <- data$variant_key[rows]
  }
  instrument_exposure <- rep.int(seq_len(exposure_count), exposure_found)
  instrument_valid <- as.logical(unlist(valid_list, use.names = FALSE))
  instrument_key <- as.character(unlist(key_list, use.names = FALSE))
  keys <- unique(instrument_key)
  instrument_unique <- match(instrument_key, keys)
  # Contiguous per-exposure group sums via cumulative sums.
  group_end <- cumsum(exposure_found)
  group_start <- group_end - exposure_found
  group_sum <- function(x) {
    cs <- c(0L, cumsum(as.integer(x)))
    cs[group_end + 1L] - cs[group_start + 1L]
  }
  outcome_found <- matrix(0L, exposure_count, outcome_count)
  invalid_outcome <- matrix(0L, exposure_count, outcome_count)
  matched <- matrix(0L, exposure_count, outcome_count)
  entry_parts <- vector("list", outcome_count)
  row_parts <- vector("list", outcome_count)
  for (o in seq_len(outcome_count)) {
    data <- outcome_data[[o]]
    hit <- match(keys, data$variant_key, nomatch = 0L)
    ok <- logical(length(keys))
    found_key <- hit > 0L
    ok[found_key] <- fastmr_compressed_valid_values(data)[hit[found_key]]
    hit_i <- hit[instrument_unique]
    ok_i <- ok[instrument_unique]
    found_i <- hit_i > 0L
    keep <- ok_i & instrument_valid
    outcome_found[, o] <- group_sum(found_i)
    invalid_outcome[, o] <- group_sum(found_i & !ok_i)
    matched[, o] <- group_sum(keep)
    kept <- which(keep)
    entry_parts[[o]] <- kept
    row_parts[[o]] <- hit_i[kept]
  }
  entry_instrument <- as.integer(unlist(entry_parts, use.names = FALSE))
  entry_outcome <- rep.int(seq_len(outcome_count), lengths(entry_parts))
  entry_outcome_row <- as.integer(unlist(row_parts, use.names = FALSE))
  rm(entry_parts, row_parts)
  entry_pair <- (instrument_exposure[entry_instrument] - 1) * outcome_count +
    entry_outcome
  if (length(entry_pair) &&
      exposure_count * outcome_count <= .Machine$integer.max) {
    entry_pair <- as.integer(entry_pair)
  }
  # Stable radix order: exposure-major pairs, instruments in supplied order.
  ord <- order(entry_pair, method = "radix")
  flat <- function(m) as.vector(t(m))
  list(
    exposure_names = exposure_names,
    outcome_names = outcome_names,
    requested = requested,
    exposure_found = exposure_found,
    invalid_exposure = invalid_exposure,
    rows = rows_list,
    instrument_exposure = instrument_exposure,
    instrument_key = instrument_key,
    outcome_found = flat(outcome_found),
    invalid_outcome = flat(invalid_outcome),
    matched = flat(matched),
    entry_instrument = entry_instrument[ord],
    entry_outcome = entry_outcome[ord],
    entry_outcome_row = entry_outcome_row[ord],
    entry_pair = entry_pair[ord]
  )
}

# Every decoded store must carry `column` with one value per row: unlist()
# silently drops a NULL (absent) column, which would misalign the gather.
fastmr_compressed_check_column <- function(data_list, values, column, side) {
  bad <- vapply(values, is.null, logical(1)) |
    lengths(values) != vapply(data_list, NROW, integer(1))
  if (any(bad)) {
    stop("compressed ", side, " store(s) lack a complete '", column, "' column: ",
         paste(names(data_list)[bad], collapse = ", "), call. = FALSE)
  }
}

# One column of the exposure stores, per instrument of a pair index.
fastmr_compressed_instrument_column <- function(exposure_data, index, column) {
  values <- lapply(exposure_data, `[[`, column)
  fastmr_compressed_check_column(exposure_data, values, column, "exposure")
  unlist(lapply(seq_along(values), function(e) {
    values[[e]][index$rows[[e]]]
  }), use.names = FALSE)
}

# One column of the outcome stores at (outcome, row) positions.
fastmr_compressed_outcome_column <- function(outcome_data, column, outcome, row) {
  values <- lapply(outcome_data, `[[`, column)
  fastmr_compressed_check_column(outcome_data, values, column, "outcome")
  offset <- c(0, cumsum(as.numeric(lengths(values))))
  unlist(values, use.names = FALSE)[offset[outcome] + row]
}

# Maximum number of affected pairs/exposures listed in each omission warning
# of fast_mr_compressed(strict = FALSE).  Longer lists are truncated with the
# total count; the per-pair detail is always in the `compressed_input$counts`
# attribute.  Set options(fastMR.warning_pairs = Inf) for the full listing.
fastmr_warning_pair_limit <- function() {
  limit <- getOption("fastMR.warning_pairs", 20)
  if (length(limit) != 1L || !is.numeric(limit) || is.na(limit) || limit < 0) {
    stop("option fastMR.warning_pairs must be a single non-negative number", call. = FALSE)
  }
  limit
}

# Emits one warning listing the first `limit` of `total` messages; `render(k)`
# builds only the messages that are shown.
fastmr_omission_warning <- function(prefix, total, render, limit) {
  if (!total) return(invisible(NULL))
  shown <- if (limit >= total) total else as.integer(limit)
  listed <- if (shown) paste(render(shown), collapse = "; ") else ""
  message <- paste0(prefix, listed)
  if (shown < total) {
    message <- paste0(
      message, if (shown) "; " else "", "... and ", total - shown,
      " more (", total, " in total; per-pair counts are in ",
      "attr(result, \"compressed_input\")$counts)"
    )
  }
  warning(message, call. = FALSE)
}

# Strict-mode errors, the no-retained-pair error and strict = FALSE omission
# warnings shared by the sparse IVW and pairwise paths.  Exposure-level inputs
# have one value per exposure; pair-level inputs are exposure-major.  Errors
# name the first failure in the order the original per-pair loop met it.
# Returns the logical vector of retained pairs.
fastmr_compressed_screen_pairs <- function(
    exposure_names, outcome_names, requested, exposure_found, invalid_exposure,
    outcome_found, invalid_outcome, matched, minimum_snps, strict) {
  exposure_count <- length(exposure_names)
  outcome_count <- length(outcome_names)
  exposure_of <- function(pair) (pair - 1L) %/% outcome_count + 1L
  outcome_of <- function(pair) (pair - 1L) %% outcome_count + 1L
  pair_requested <- rep(requested, each = outcome_count)
  pair_invalid_exposure <- rep(invalid_exposure, each = outcome_count)
  missing_exposure <- exposure_found < requested
  missing_outcome <- outcome_found < pair_requested
  invalid_pair <- (pair_invalid_exposure > 0L) | (invalid_outcome > 0L)
  too_few <- matched < minimum_snps
  pair_label <- function(pair) {
    paste0(exposure_names[exposure_of(pair)], " -> ", outcome_names[outcome_of(pair)])
  }

  if (isTRUE(strict)) {
    pair_error <- missing_outcome | invalid_pair | too_few
    first_pair <- if (any(pair_error)) which(pair_error)[[1L]] else NA_integer_
    first_exposure <- if (any(missing_exposure)) which(missing_exposure)[[1L]] else NA_integer_
    if (!is.na(first_exposure) &&
        (is.na(first_pair) || first_exposure <= exposure_of(first_pair))) {
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
          "missing requested outcome instrument(s) for ", pair_label(i),
          " (found ", outcome_found[[i]], " of ", pair_requested[[i]], ")",
          call. = FALSE
        )
      }
      if (invalid_pair[[i]]) {
        stop(
          "invalid beta/standard_error for ", pair_label(i),
          " (exposure=", pair_invalid_exposure[[i]],
          ", outcome=", invalid_outcome[[i]], ")", call. = FALSE
        )
      }
      stop("fewer than minimum_snps for ", pair_label(i),
           " (", matched[[i]], " matched)", call. = FALSE)
    }
  }

  retained <- !too_few
  if (!any(retained)) {
    stop("no exposure-outcome pair retained enough matched instruments", call. = FALSE)
  }
  if (!isTRUE(strict)) {
    limit <- fastmr_warning_pair_limit()
    skipped <- which(too_few)
    fastmr_omission_warning(
      "omitted pair(s) below minimum_snps: ", length(skipped),
      function(k) {
        i <- skipped[seq_len(k)]
        paste0(pair_label(i), " (", matched[i], " matched)")
      },
      limit
    )
    invalid <- which(invalid_pair)
    fastmr_omission_warning(
      "omitted invalid instrument value(s): ", length(invalid),
      function(k) {
        i <- invalid[seq_len(k)]
        paste0(pair_label(i), " (exposure=", pair_invalid_exposure[i],
               ", outcome=", invalid_outcome[i], ")")
      },
      limit
    )
    # Per exposure: the exposure message, then its outcome messages.  Labels
    # are unique, so every message is distinct.
    missing_e <- which(missing_exposure)
    missing_p <- which(missing_outcome)
    order_key <- c(
      (missing_e - 1) * (outcome_count + 1),
      (exposure_of(missing_p) - 1) * (outcome_count + 1) + outcome_of(missing_p)
    )
    is_pair <- c(rep(FALSE, length(missing_e)), rep(TRUE, length(missing_p)))
    what <- c(missing_e, missing_p)
    fastmr_omission_warning(
      "omitted missing requested instrument(s): ", length(what),
      function(k) {
        pick <- utils::head(order(order_key, method = "radix"), k)
        i <- what[pick]
        pair <- is_pair[pick]
        out <- character(length(pick))
        p <- i[pair]
        out[pair] <- paste0(pair_label(p), " outcome (found ", outcome_found[p],
                            " of ", pair_requested[p], ")")
        e <- i[!pair]
        out[!pair] <- paste0(exposure_names[e], " exposure (found ",
                             exposure_found[e], " of ", requested[e], ")")
        out
      },
      limit
    )
  }
  retained
}

# Pairwise estimator input for compressed stores: the harmonised long table
# passed to one batched fast_mr() call, plus per-pair extraction counts.
fastmr_compressed_pairwise_data <- function(exposure_data, outcome_data,
                                            instrument_sets, minimum_snps, strict) {
  index <- fastmr_compressed_pair_index(exposure_data, outcome_data, instrument_sets)
  outcome_count <- length(index$outcome_names)
  retained <- fastmr_compressed_screen_pairs(
    index$exposure_names, index$outcome_names, index$requested,
    index$exposure_found, index$invalid_exposure, index$outcome_found,
    index$invalid_outcome, index$matched, minimum_snps, strict
  )
  pair_exposure <- rep(index$exposure_names, each = outcome_count)
  pair_outcome <- rep(index$outcome_names, times = length(index$exposure_names))
  counts <- data.frame(
    id.exposure = pair_exposure, id.outcome = pair_outcome,
    requested = rep(index$requested, each = outcome_count),
    exposure_found = rep(index$exposure_found, each = outcome_count),
    outcome_found = index$outcome_found,
    invalid_exposure = rep(index$invalid_exposure, each = outcome_count),
    invalid_outcome = index$invalid_outcome, matched = index$matched,
    stringsAsFactors = FALSE
  )
  rows <- fastmr_compressed_entry_table(
    exposure_data, outcome_data, index, retained[index$entry_pair],
    c("beta", "standard_error")
  )
  data <- data.frame(
    SNP = rows$SNP,
    beta.exposure = rows$values$beta$exposure,
    beta.outcome = rows$values$beta$outcome,
    se.exposure = rows$values$standard_error$exposure,
    se.outcome = rows$values$standard_error$outcome,
    id.exposure = rows$id.exposure,
    id.outcome = rows$id.outcome,
    stringsAsFactors = FALSE
  )
  list(data = data, counts = counts, index = index)
}

# Harmonised rows of the selected entries.  `columns` names the store columns
# gathered from both sides; the result has TwoSampleMR-style names.
fastmr_compressed_entry_table <- function(exposure_data, outcome_data, index,
                                          selected, columns) {
  instrument <- index$entry_instrument[selected]
  outcome <- index$entry_outcome[selected]
  outcome_row <- index$entry_outcome_row[selected]
  exposure <- index$instrument_exposure[instrument]
  side <- function(column) {
    list(
      exposure = fastmr_compressed_instrument_column(
        exposure_data, index, column
      )[instrument],
      outcome = fastmr_compressed_outcome_column(
        outcome_data, column, outcome, outcome_row
      )
    )
  }
  values <- stats::setNames(lapply(columns, side), columns)
  list(
    SNP = index$instrument_key[instrument],
    id.exposure = index$exposure_names[exposure],
    id.outcome = index$outcome_names[outcome],
    values = values
  )
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
  }, error = function(e) FALSE) &&
    fastmr_sparse_ivw_fits(exposure_count, outcome_count, snp_count)
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

  retained <- fastmr_compressed_screen_pairs(
    exposure_names, outcome_names, requested, exposure_found, invalid_exposure,
    pair_outcome_found, pair_invalid_outcome, pair_matched, minimum_snps, strict
  )

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
  # One instrument: the native Wald ratio, with Q undefined (as fast_mr()).
  q_df <- ifelse(has_fit & n >= 2, n - 1, NA_real_)
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
# `index` may be the pair index the pairwise path already built for the same
# stores and instruments; otherwise it is built here.
fastmr_compressed_steiger <- function(exposure_data, outcome_data, instrument_sets,
                                      minimum_snps, options, index = NULL) {
  if (is.null(index)) {
    index <- fastmr_compressed_pair_index(exposure_data, outcome_data, instrument_sets)
  }
  selected <- (index$matched >= minimum_snps)[index$entry_pair]
  if (!any(selected)) return(data.frame())
  rows <- fastmr_compressed_entry_table(
    exposure_data, outcome_data, index, selected,
    c("beta", "standard_error", "effect_allele_frequency", "p_value")
  )
  id_exposure <- rows$id.exposure
  id_outcome <- rows$id.outcome
  values <- rows$values
  data <- data.frame(
    SNP = rows$SNP,
    id.exposure = id_exposure,
    id.outcome = id_outcome,
    exposure = id_exposure,
    outcome = id_outcome,
    beta.exposure = values$beta$exposure,
    beta.outcome = values$beta$outcome,
    se.exposure = values$standard_error$exposure,
    se.outcome = values$standard_error$outcome,
    eaf.exposure = values$effect_allele_frequency$exposure,
    eaf.outcome = values$effect_allele_frequency$outcome,
    pval.exposure = values$p_value$exposure,
    pval.outcome = values$p_value$outcome,
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
      hit <- match(data[[paste0("id.", side)]], binary$id)
      data[[paste0("units.", side)]][!is.na(hit)] <- "log odds"
      data[[paste0("ncase.", side)]] <- binary$ncase[hit]
      data[[paste0("ncontrol.", side)]] <- binary$ncontrol[hit]
      data[[paste0("prevalence.", side)]] <- binary$prevalence[hit]
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
#' When every exposure has the same instrument set, no requested method uses
#' bootstrap draws and every pair has all instruments with valid values, the
#' run uses the shared-grid kernel [fast_mr_grid()] whatever `estimator` says
#' (`estimator = "pairwise"` does not disable it; `estimator_path` in the
#' `compressed_input` attribute is then `"shared_instrument_grid"`).  The
#' sparse IVW path is used only when its estimated memory (native results,
#' the outcome-by-union-instrument matrices it builds and its per-pair count
#' matrices) fits `getOption("fastMR.sparse_ivw_max_memory_mb", 8192)` MiB.
#' An exposure with an empty instrument set is an error with `strict = TRUE`
#' and is dropped with a warning with `strict = FALSE`.
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
#'   Each warning lists at most `getOption("fastMR.warning_pairs", 20)`
#'   affected pairs/exposures followed by the total count; the per-pair detail
#'   is always in the `counts` element of the `compressed_input` attribute.
#'   Set `options(fastMR.warning_pairs = Inf)` for the full listing.
#' @param output Optional path for a Zstandard-compressed Parquet copy of the
#'   result. The path must not already exist; use [fast_write_parquet()] when
#'   an overwrite or another compression codec is required.
#' @param estimator `"auto"` (default) or `"pairwise"`. With `"auto"`, runs
#'   requesting only `methods = "ivw"` (and no extra `...` options) whose
#'   instrument sets are not shared by all exposures use the sparse CSR kernel
#'   [fast_mr_sparse_ivw()] instead of the pairwise path, which assembles one
#'   harmonised long table for every pair and makes one batched [fast_mr()]
#'   call.
#'   Counts, errors, warnings, `minimum_snps` handling and row order are
#'   identical, but estimates may differ from `"pairwise"` by rounding
#'   (agreement is within 1e-14 relative on beta, se and Q; the two kernels use
#'   the same two-pass residual formula and summation order, so results are
#'   usually bit-identical). `"pairwise"`
#'   (and every other method set) always uses the pairwise path. The path used
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
  with_instruments <- fastmr_compressed_nonempty_instruments(instrument_sets, strict)
  exposure_files <- exposure_files[with_instruments]
  instrument_sets <- instrument_sets[with_instruments]
  union_keys <- unique(unlist(instrument_sets, use.names = FALSE))
  # Matching uses variant_key only, so no decoded identity strings are read.
  columns <- c("beta", "standard_error")
  if (!is.null(steiger_options)) {
    columns <- c(columns, "effect_allele_frequency", "p_value")
  }
  # Attach Steiger (when requested) after estimation, so it never runs on a
  # study that strict mode or minimum_snps rejected.
  # The Parquet copy is written before the `steiger` attribute is attached.
  finish <- function(result, index = NULL) {
    if (is.null(steiger_options)) return(fastmr_write_result(result, output))
    steiger_started <- unname(proc.time()[["elapsed"]])
    steiger_result <- fastmr_compressed_steiger(
      exposure_data, outcome_data, instrument_sets, minimum_snps, steiger_options,
      index = index
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
  store_paths <- c(unname(exposure_files), unname(outcome_files))
  codecs <- fastmr_compressed_validate_stores(unique(store_paths), as.integer(io_threads))
  all_data <- fastmr_io_map(
    store_paths,
    c(unname(instrument_sets), unname(outcome_keys)),
    columns,
    as.integer(io_threads),
    codecs = unname(codecs[store_paths])
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

  pairwise <- fastmr_compressed_pairwise_data(
    exposure_data, outcome_data, instrument_sets, minimum_snps, strict
  )
  result <- do.call(fast_mr, c(list(
    data = pairwise$data, methods = methods, nboot = controls$nboot,
    seed = controls$seed, threads = controls$threads
  ), dots))
  attr(result, "compressed_input") <- list(
    exposure_files = exposure_files,
    outcome_files = outcome_files,
    instruments = instrument_sets,
    counts = pairwise$counts,
    io_threads = as.integer(io_threads),
    estimator_path = "pairwise",
    timing = list(
      io_seconds = io_seconds,
      estimator_seconds = unname(proc.time()[["elapsed"]]) - estimator_started,
      total_seconds = unname(proc.time()[["elapsed"]]) - total_started,
      source_bytes_read = source_bytes_read
    )
  )
  finish(result, index = pairwise$index)
}
