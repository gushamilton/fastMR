# Verbatim copies of the per-pair compressed implementations from main at
# 5c87455, before the vectorised pair assembly.  Used as the reference in
# test-compressed-pairwise-assembly.R.  The pairwise loop body is unchanged;
# it returns the rbind()-ed long table it passed to fast_mr() and the counts,
# emitting the same errors and warnings as the old fast_mr_compressed().
old_compressed_pairwise_data <- function(exposure_data, outcome_data,
                                         instrument_sets, minimum_snps, strict) {
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
  list(data = do.call(rbind, rows), counts = do.call(rbind, counts))
}

old_compressed_steiger <- function(exposure_data, outcome_data, instrument_sets,
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
