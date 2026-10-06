# Multivariable and residualised MR directly from CompreSSoR stores.

fastmr_mvmr_key_fields <- function(keys) {
  fields <- strsplit(keys, ":", fixed = TRUE)
  ok <- lengths(fields) == 4L
  out <- matrix(NA_character_, length(keys), 4L)
  if (any(ok)) out[ok, ] <- do.call(rbind, fields[ok])
  out
}

# A covariate data frame -> data.frame(key, beta, se) with canonical keys.
fastmr_mvmr_covariate_frame <- function(data, label) {
  names(data) <- tolower(names(data))
  pick <- function(choices) {
    hit <- intersect(choices, names(data))
    if (!length(hit)) {
      stop("covariate '", label, "' needs one of the columns ",
           paste(choices, collapse = "/"), call. = FALSE)
    }
    data[[hit[[1L]]]]
  }
  key <- trimws(as.character(pick(c("variant_key", "snp", "key", "id"))))
  beta <- fastmr_numeric(pick(c("beta", "b")), paste0("covariate ", label, " beta"))
  se <- fastmr_numeric(pick(c("standard_error", "se")), paste0("covariate ", label, " se"))
  keep <- !is.na(key) & nzchar(key)
  out <- data.frame(key = key[keep], beta = beta[keep], se = se[keep],
                    stringsAsFactors = FALSE)
  if (anyDuplicated(out$key)) {
    stop("covariate '", label, "' has duplicated variant keys", call. = FALSE)
  }
  out
}

# Normalises `covariates` to a named list whose elements are either a store
# path (class "store") or a key/beta/se data frame.
fastmr_mvmr_normalize_covariates <- function(covariates) {
  if (is.data.frame(covariates)) {
    lower <- tolower(names(covariates))
    trait_col <- intersect(c("trait", "id.covariate", "covariate"), lower)
    if (length(trait_col)) {
      traits <- as.character(covariates[[match(trait_col[[1L]], lower)]])
      covariates <- split(covariates, factor(traits, levels = unique(traits)))
    } else {
      covariates <- list(covariate = covariates)
    }
  } else if (is.character(covariates)) {
    covariates <- as.list(covariates)
    if (is.null(names(covariates))) {
      names(covariates) <- basename(sub("[\\/]+$", "", unlist(covariates)))
    }
  }
  if (!is.list(covariates) || !length(covariates)) {
    stop("covariates must be a data frame, a list of data frames/stores or a ",
         "character vector of stores", call. = FALSE)
  }
  labels <- names(covariates)
  if (is.null(labels) || any(!nzchar(labels)) || anyDuplicated(labels)) {
    stop("covariates must have unique non-empty names", call. = FALSE)
  }
  if ("exposure" %in% labels) {
    stop("'exposure' is reserved for the primary exposure; rename the covariate",
         call. = FALSE)
  }
  out <- lapply(labels, function(label) {
    value <- covariates[[label]]
    if (is.character(value) && length(value) == 1L) {
      path <- fastmr_normalize_compressed_files(stats::setNames(value, label),
                                                "covariates")
      return(structure(unname(path), class = "fastmr_store"))
    }
    if (!is.data.frame(value)) {
      stop("covariate '", label, "' must be a data frame or one store path",
           call. = FALSE)
    }
    fastmr_mvmr_covariate_frame(value, label)
  })
  names(out) <- labels
  out
}

# Values of a key/beta/se frame at canonical keys; a key found only with its
# alleles swapped (chr:pos:ALT:REF) has its beta negated.
fastmr_mvmr_frame_lookup <- function(frame, keys) {
  hit <- match(keys, frame$key)
  beta <- frame$beta[hit]
  se <- frame$se[hit]
  missing <- which(is.na(hit))
  if (length(missing)) {
    f <- fastmr_mvmr_key_fields(keys[missing])
    swapped <- paste(f[, 1L], f[, 2L], f[, 4L], f[, 3L], sep = ":")
    swap_hit <- match(swapped, frame$key)
    ok <- !is.na(swap_hit) & !is.na(f[, 1L])
    beta[missing[ok]] <- -frame$beta[swap_hit[ok]]
    se[missing[ok]] <- frame$se[swap_hit[ok]]
  }
  list(beta = beta, se = se)
}

fastmr_mvmr_store_lookup <- function(data, keys) {
  hit <- match(keys, data$variant_key)
  list(beta = data$beta[hit], se = data$standard_error[hit])
}

fastmr_mvmr_valid <- function(beta, se) {
  is.finite(beta) & is.finite(se) & se > 0
}

# Covariate instrument keys per trait.
fastmr_mvmr_normalize_covariate_instruments <- function(value, traits) {
  if (is.null(value)) return(NULL)
  if (is.character(value)) {
    if (length(traits) != 1L) {
      stop("covariate_instruments must be a list named by covariate when there ",
           "is more than one covariate", call. = FALSE)
    }
    value <- stats::setNames(list(value), traits)
  }
  if (!is.list(value) || is.null(names(value)) || !all(traits %in% names(value))) {
    stop("covariate_instruments must be a list with one key vector per covariate (",
         paste(traits, collapse = ", "), ")", call. = FALSE)
  }
  lapply(value[traits], function(keys) {
    if (!length(keys)) character() else fastmr_normalize_variant_keys(keys)
  })
}

# Greedy LD resolution: keep candidates in order of decreasing strength,
# dropping any candidate in LD with an already kept one.  `index` are union
# indices, `strength` their max |z|, `neighbours` a list of union indices.
fastmr_mvmr_greedy <- function(index, strength, neighbours, dropped) {
  order_index <- index[order(-strength, seq_along(strength), method = "radix")]
  kept <- integer()
  for (i in order_index) {
    if (dropped[i]) next
    kept <- c(kept, i)
    nb <- neighbours[[i]]
    if (length(nb)) dropped[nb] <- TRUE
  }
  list(kept = kept)
}

fastmr_mvmr_neighbours <- function(ld_pairs, union_keys) {
  if (is.null(ld_pairs)) return(NULL)
  if (!is.data.frame(ld_pairs) && !is.matrix(ld_pairs)) {
    stop("ld_pairs must be a two-column data frame of variant keys", call. = FALSE)
  }
  ld_pairs <- as.data.frame(ld_pairs, stringsAsFactors = FALSE)
  if (ncol(ld_pairs) < 2L) stop("ld_pairs must have two key columns", call. = FALSE)
  a <- match(as.character(ld_pairs[[1L]]), union_keys)
  b <- match(as.character(ld_pairs[[2L]]), union_keys)
  ok <- !is.na(a) & !is.na(b) & a != b
  a <- a[ok]
  b <- b[ok]
  from <- c(a, b)
  to <- c(b, a)
  out <- vector("list", length(union_keys))
  if (length(from)) {
    groups <- split(to, factor(from, levels = seq_along(union_keys)))
    out <- lapply(groups, unique)
  }
  out
}

# Logical U x K-cell exclusions from per-outcome genomic windows.
fastmr_mvmr_apply_exclusions <- function(panel_beta, union_keys, outcome_exclude,
                                         outcome_labels) {
  if (is.null(outcome_exclude)) return(list(beta = panel_beta, masked = 0))
  if (!is.data.frame(outcome_exclude)) {
    stop("outcome_exclude must be a data frame", call. = FALSE)
  }
  needed <- c("id.outcome", "chromosome", "start", "end")
  if (!all(needed %in% names(outcome_exclude))) {
    stop("outcome_exclude needs columns ", paste(needed, collapse = ", "), call. = FALSE)
  }
  f <- fastmr_mvmr_key_fields(union_keys)
  chromosome <- f[, 1L]
  position <- as.numeric(f[, 2L])
  k <- match(as.character(outcome_exclude$id.outcome), outcome_labels)
  ex_chr <- as.character(outcome_exclude$chromosome)
  ex_chr <- sub("^chr", "", ex_chr, ignore.case = TRUE)
  ex_start <- as.numeric(outcome_exclude$start)
  ex_end <- as.numeric(outcome_exclude$end)
  masked <- 0
  by_chr <- split(seq_along(chromosome), chromosome)
  for (r in which(!is.na(k))) {
    rows <- by_chr[[ex_chr[r]]]
    if (!length(rows)) next
    hit <- rows[position[rows] >= ex_start[r] & position[rows] <= ex_end[r]]
    if (length(hit)) {
      masked <- masked + sum(is.finite(panel_beta[hit, k[r]]))
      panel_beta[hit, k[r]] <- NA_real_
    }
  }
  list(beta = panel_beta, masked = masked)
}

fastmr_mvmr_drop_message <- function(what, total, labels, limit = 10L) {
  if (!total) return(invisible(NULL))
  warning(what, ": ", total, if (length(labels)) paste0(
    " (", paste(utils::head(labels, limit), collapse = "; "),
    if (length(labels) > limit) "; ..." else "", ")"
  ), call. = FALSE)
}

#' Multivariable or residualised MR from compressed GWAS stores
#'
#' Adjusts every exposure -> outcome MR estimate for one or more covariate
#' traits (for example platelet count), reading exposures and outcomes from
#' CompreSSoR stores as [fast_mr_compressed()] does, and fitting every
#' exposure against every outcome with the batched kernel of
#' [fast_mvmr_ivw_batch()].
#'
#' **`method = "mvmr"`** fits classical multivariable IVW (TwoSampleMR
#' `mv_ivw` / `mv_multiple` with shared instruments) with design columns
#' `exposure` and each covariate, on the exposure's instruments plus the
#' covariates' instruments (`covariate_instruments`), optionally LD-resolved
#' with `ld_pairs`: candidates are taken in decreasing order of their
#' strongest |z| across the design columns and any candidate in LD with an
#' already kept one is dropped. The Sanderson-Windmeijer conditional F of
#' each design column is reported, with a warning when it is below `weak_f`.
#' In a one-sample setting, where the exposure, covariate and outcome GWAS
#' share participants, weak conditional instruments bias the estimates toward
#' the observational (confounded) associations and can inflate them; check the
#' conditional F before interpreting the estimates.
#'
#' **`method = "residualised"`** keeps each exposure's own instruments and
#' removes the covariate-mediated part of every outcome effect,
#' `b_y* = b_y - sum_t b_x,t * Gamma_t`, where `Gamma_t` is the covariate's
#' effect on the outcome estimated by IVW (joint multivariable IVW for several
#' covariates) from the covariates' own instruments, on the same outcome
#' stores. Its variance is propagated: `var(b_y*) = se_y^2 + b_x,T' V_Gamma
#' b_x,T`, with `V_Gamma` the covariance of the `Gamma` estimates. The
#' exposure is then fitted by IVW on `b_y*`. This is an approximation, not a
#' formal multivariable MR: it assumes one homogeneous covariate -> outcome
#' effect that applies equally at every exposure instrument (instruments that
#' act on the covariate through another pathway leave residue), and it ignores
#' the covariance between `Gamma` and `b_y` (shared outcome sample) and
#' between the covariate and exposure estimates. It was motivated by an
#' analysis of UKB-PPP trans protein MR in which the classical multivariable
#' fit had a median conditional F of about 3.5 while the residualised
#' adjustment behaved as expected; that analysis is a motivation, not a
#' validation of the estimator.
#'
#' **Covariates** may be CompreSSoR stores or external summary statistics:
#' a data frame with columns `SNP` (or `variant_key`) holding canonical
#' `chromosome:position:REF:ALT` keys (GRCh38, as the stores), `beta` and `se`
#' (or `standard_error`) for the ALT allele, plus an optional `trait` column
#' for several covariates; or a named list of such data frames and/or store
#' paths. A data-frame key found only with its alleles swapped is used with
#' its beta negated.
#'
#' **Missing data.** With `strict = TRUE` any requested instrument missing
#' from, or invalid in, an exposure, covariate or outcome store is an error;
#' with `strict = FALSE` design SNPs without valid exposure/covariate values
#' are dropped and missing outcome values are left out of that outcome's fit,
#' with warnings. `outcome_exclude` masks SNPs for single outcomes (e.g. each
#' outcome's cis window); masked cells are not "missing". Pairs with fewer than
#' `max(minimum_snps, 1)` SNPs are dropped from the long result (an error when
#' `strict = TRUE`); in the matrix result they are `NA`.
#'
#' An exposure store that is also an outcome store is not read twice: its
#' values come from the outcome read, before `outcome_exclude` is applied.
#'
#' @param exposure_files Named character vector of exposure stores.
#' @param outcome_files Named character vector of outcome stores.
#' @param covariates Covariate traits (see Details).
#' @param instruments Canonical-key vector shared by every exposure, or a
#'   list with one key vector per exposure (each exposure's own instruments).
#' @param covariate_instruments Keys of each covariate's own instruments: a
#'   character vector for one covariate or a list named by covariate.
#'   Required for `"residualised"`; for `"mvmr"` they are added to every
#'   exposure's design (omit them when `instruments` already holds the union).
#' @param method `"mvmr"` or `"residualised"`.
#' @param exposure_cor Phenotypic correlation of the design columns (named
#'   `exposure` and the covariate labels), scaled by sample overlap, for the
#'   conditional F and `Q_A`: one matrix for every exposure, or a list named by
#'   exposure. `NULL` assumes zero correlation, with a warning.
#' @param ld_pairs Optional two-column data frame of canonical keys in LD,
#'   used to LD-resolve unions of instruments (see Details).
#' @param outcome_exclude Optional data frame with `id.outcome`,
#'   `chromosome`, `start` and `end`: SNPs in the window are not used for that
#'   outcome (in designs and in the covariate fits).
#' @param se_model Standard-error model of [fast_mvmr_ivw()]. Defaults to
#'   `"multiplicative"` (TwoSampleMR `mv_ivw`) for `"mvmr"` and
#'   `"multiplicative_floored"` (TwoSampleMR / fastMR univariable IVW) for
#'   `"residualised"`.
#' @param covariate_se_model Standard-error model of the `Gamma` fits
#'   (residualised only); default `"multiplicative_floored"`.
#' @param weights,shared_tolerance As in [fast_mvmr_ivw()].
#' @param threads Native worker count.
#' @param io_threads Stores decoded concurrently by CompreSSoR.
#' @param minimum_snps Minimum SNPs per exposure-outcome fit.
#' @param strict See Details.
#' @param output_format `"long"` (a tidy data frame, exposure-major) or
#'   `"matrix"` (a list of exposure x outcome matrices; better for large
#'   grids).
#' @param covariate_estimates For `"mvmr"`, also return the covariates'
#'   coefficients (`b_<covariate>`, `se_<covariate>`, `pval_<covariate>`).
#' @param weak_f Conditional F threshold for the weak-instrument warning.
#' @param output Optional Parquet path for the long result.
#' @return With `output_format = "long"`, a data frame with `id.exposure`,
#'   `id.outcome`, `method`, `method_code`, `nsnp`, `b`, `se`, `pval`, `Q`,
#'   `Q_df`, `Q_pval`, `Q_A`, `Q_A_pval`, `sigma` and `conditional_F` (of the
#'   exposure column). The `mvmr_input` attribute holds the stores, design
#'   SNP sets, per-exposure diagnostics (`diagnostics`: design SNP count and
#'   conditional F of every column), the `Gamma` fits (`covariate_fit`,
#'   residualised), counts and `timing` (`io_seconds`, `design_seconds`,
#'   `estimator_seconds`, `total_seconds`, `source_bytes_read`). With
#'   `"matrix"`, a list with the same quantities as matrices and the same
#'   metadata in `mvmr_input`.
#' @seealso [fast_mvmr_ivw()], [fast_mvmr_ivw_batch()], [fast_mr_compressed()]
#' @export
fast_mvmr_compressed <- function(
    exposure_files,
    outcome_files,
    covariates,
    instruments,
    covariate_instruments = NULL,
    method = c("mvmr", "residualised"),
    exposure_cor = NULL,
    ld_pairs = NULL,
    outcome_exclude = NULL,
    se_model = NULL,
    covariate_se_model = "multiplicative_floored",
    weights = c("exact", "shared"),
    shared_tolerance = 1e-8,
    threads = 1,
    io_threads = 1,
    minimum_snps = 1L,
    strict = TRUE,
    output_format = c("long", "matrix"),
    covariate_estimates = FALSE,
    weak_f = 10,
    output = NULL) {
  total_started <- unname(proc.time()[["elapsed"]])
  method <- match.arg(method)
  weights <- match.arg(weights)
  output_format <- match.arg(output_format)
  if (is.null(se_model)) {
    se_model <- if (method == "mvmr") "multiplicative" else "multiplicative_floored"
  }
  se_model <- match.arg(se_model, fastmr_mvmr_se_models)
  covariate_se_model <- match.arg(covariate_se_model, fastmr_mvmr_se_models)
  fastmr_require_compressor()
  threads <- fastmr_mvmr_validate_common(se_model, weights, shared_tolerance,
                                         threads, weak_f)
  io_threads <- fastmr_positive_integer_scalar(io_threads, "io_threads")
  minimum_snps <- fastmr_positive_integer_scalar(minimum_snps, "minimum_snps")
  strict <- fastmr_mvmr_check_flag(strict, "strict")
  covariate_estimates <- fastmr_mvmr_check_flag(covariate_estimates,
                                                "covariate_estimates")
  if (!is.null(output) && output_format != "long") {
    stop("output requires output_format = \"long\"", call. = FALSE)
  }
  exposure_files <- fastmr_normalize_compressed_files(exposure_files, "exposure_files")
  outcome_files <- fastmr_normalize_compressed_files(outcome_files, "outcome_files")
  covariate_sources <- fastmr_mvmr_normalize_covariates(covariates)
  traits <- names(covariate_sources)
  if (any(traits %in% names(exposure_files))) {
    stop("covariate labels must differ from exposure labels", call. = FALSE)
  }
  instrument_sets <- fastmr_normalize_instruments(instruments, names(exposure_files))
  with_instruments <- fastmr_compressed_nonempty_instruments(instrument_sets, strict)
  exposure_files <- exposure_files[with_instruments]
  instrument_sets <- instrument_sets[with_instruments]
  cov_inst <- fastmr_mvmr_normalize_covariate_instruments(covariate_instruments, traits)
  if (method == "residualised" && is.null(cov_inst)) {
    stop("method = \"residualised\" needs covariate_instruments to estimate the ",
         "covariate -> outcome effects", call. = FALSE)
  }
  exposure_labels <- names(exposure_files)
  outcome_labels <- names(outcome_files)
  E <- length(exposure_labels)
  K <- length(outcome_labels)
  cov_union <- unique(unlist(cov_inst, use.names = FALSE))
  candidates <- if (method == "mvmr") {
    lapply(instrument_sets, function(keys) unique(c(keys, cov_union)))
  } else {
    instrument_sets
  }
  union_keys <- unique(c(unlist(candidates, use.names = FALSE), cov_union))
  U <- length(union_keys)

  # ---- one batched read of every store ------------------------------------
  io_started <- unname(proc.time()[["elapsed"]])
  shared_store <- exposure_files %in% outcome_files
  read_exposures <- exposure_labels[!shared_store]
  store_traits <- traits[vapply(covariate_sources, inherits, logical(1), "fastmr_store")]
  paths <- c(unname(outcome_files), unname(exposure_files[read_exposures]),
             vapply(covariate_sources[store_traits], unclass, character(1)))
  keys <- c(rep(list(union_keys), K), unname(candidates[read_exposures]),
            rep(list(union_keys), length(store_traits)))
  codecs <- fastmr_compressed_validate_stores(unique(paths), io_threads)
  all_data <- fastmr_io_map(paths, keys, c("beta", "standard_error"), io_threads,
                            codecs = unname(codecs[paths]))
  source_bytes_read <- attr(all_data, "source_bytes_read", exact = TRUE)
  source_bytes_read <- if (is.null(source_bytes_read)) NA_real_ else
    as.numeric(source_bytes_read)
  io_seconds <- unname(proc.time()[["elapsed"]]) - io_started

  design_started <- unname(proc.time()[["elapsed"]])
  # ---- outcome panel ------------------------------------------------------
  panel_beta <- matrix(NA_real_, U, K, dimnames = list(union_keys, outcome_labels))
  panel_se <- panel_beta
  outcome_found <- integer(K)
  outcome_valid <- integer(K)
  for (k in seq_len(K)) {
    v <- fastmr_mvmr_store_lookup(all_data[[k]], union_keys)
    ok <- fastmr_mvmr_valid(v$beta, v$se)
    outcome_found[k] <- sum(!is.na(v$beta) | !is.na(v$se))
    outcome_valid[k] <- sum(ok)
    panel_beta[ok, k] <- v$beta[ok]
    panel_se[ok, k] <- v$se[ok]
  }
  exposure_data <- vector("list", E)
  names(exposure_data) <- exposure_labels
  offset <- K
  for (e in seq_len(E)) {
    label <- exposure_labels[e]
    wanted <- candidates[[label]]
    if (shared_store[e]) {
      k <- match(exposure_files[[e]], outcome_files)
      rows <- match(wanted, union_keys)
      exposure_data[[e]] <- list(beta = panel_beta[rows, k], se = panel_se[rows, k])
    } else {
      offset <- offset + 1L
      exposure_data[[e]] <- fastmr_mvmr_store_lookup(all_data[[offset]], wanted)
    }
  }
  trait_beta <- matrix(NA_real_, U, length(traits), dimnames = list(union_keys, traits))
  trait_se <- trait_beta
  for (t in traits) {
    source <- covariate_sources[[t]]
    v <- if (inherits(source, "fastmr_store")) {
      offset <- offset + 1L
      fastmr_mvmr_store_lookup(all_data[[offset]], union_keys)
    } else {
      fastmr_mvmr_frame_lookup(source, union_keys)
    }
    ok <- fastmr_mvmr_valid(v$beta, v$se)
    trait_beta[ok, t] <- v$beta[ok]
    trait_se[ok, t] <- v$se[ok]
  }
  rm(all_data)

  # Outcome completeness (masked cells excluded below are not "missing").
  missing_outcome <- U - outcome_valid
  if (any(missing_outcome > 0L)) {
    bad <- which(missing_outcome > 0L)
    if (strict) {
      stop("missing or invalid requested instrument(s) in outcome store ",
           outcome_labels[bad[1L]], " (", missing_outcome[bad[1L]], " of ", U,
           "); use strict = FALSE to leave them out of that outcome's fits",
           call. = FALSE)
    }
    fastmr_mvmr_drop_message(
      "missing or invalid outcome values left out of the fits (outcome stores affected)",
      length(bad), paste0(outcome_labels[bad], " (", missing_outcome[bad], ")"))
  }
  excluded <- fastmr_mvmr_apply_exclusions(panel_beta, union_keys, outcome_exclude,
                                           outcome_labels)
  panel_beta <- excluded$beta

  # ---- designs --------------------------------------------------------------
  neighbours <- fastmr_mvmr_neighbours(ld_pairs, union_keys)
  no_ld <- logical(U)
  trait_valid <- if (length(traits)) rowSums(!is.finite(trait_beta)) == 0 else rep(TRUE, U)
  designs <- vector("list", E)
  names(designs) <- exposure_labels
  design_keys <- vector("list", E)
  names(design_keys) <- exposure_labels
  lost_exposure <- integer(E)
  lost_covariate <- integer(E)
  for (e in seq_len(E)) {
    wanted <- candidates[[e]]
    rows <- match(wanted, union_keys)
    bx <- exposure_data[[e]]$beta
    sx <- exposure_data[[e]]$se
    ok_x <- fastmr_mvmr_valid(bx, sx)
    ok_t <- trait_valid[rows]
    lost_exposure[e] <- sum(!ok_x)
    lost_covariate[e] <- sum(ok_x & !ok_t)
    keep <- ok_x & ok_t
    rows <- rows[keep]
    bx <- bx[keep]
    sx <- sx[keep]
    if (method == "mvmr") {
      beta <- cbind(exposure = bx, trait_beta[rows, , drop = FALSE])
      se <- cbind(exposure = sx, trait_se[rows, , drop = FALSE])
      if (!is.null(neighbours) && length(rows)) {
        strength <- apply(abs(beta / se), 1L, max)
        g <- fastmr_mvmr_greedy(rows, strength, neighbours, no_ld)
        pick <- match(g$kept, rows)
        rows <- rows[pick]
        beta <- beta[pick, , drop = FALSE]
        se <- se[pick, , drop = FALSE]
      }
    } else {
      beta <- matrix(bx, ncol = 1L, dimnames = list(NULL, "exposure"))
      se <- matrix(sx, ncol = 1L, dimnames = list(NULL, "exposure"))
    }
    designs[[e]] <- list(rows = rows, beta = beta, se = se)
    design_keys[[e]] <- union_keys[rows]
  }
  if (any(lost_exposure > 0L)) {
    bad <- which(lost_exposure > 0L)
    if (strict) {
      stop("missing or invalid requested instrument(s) for exposure ",
           exposure_labels[bad[1L]], " (", lost_exposure[bad[1L]], " of ",
           length(candidates[[bad[1L]]]), ")", call. = FALSE)
    }
    fastmr_mvmr_drop_message("design SNPs dropped for missing or invalid exposure values",
                             sum(lost_exposure),
                             paste0(exposure_labels[bad], " (", lost_exposure[bad], ")"))
  }
  if (any(lost_covariate > 0L)) {
    bad <- which(lost_covariate > 0L)
    if (strict) {
      stop("covariate value(s) missing or invalid at instrument(s) of exposure ",
           exposure_labels[bad[1L]], " (", lost_covariate[bad[1L]], ")", call. = FALSE)
    }
    fastmr_mvmr_drop_message("design SNPs dropped for missing or invalid covariate values",
                             sum(lost_covariate),
                             paste0(exposure_labels[bad], " (", lost_covariate[bad], ")"))
  }
  empty <- vapply(designs, function(d) !length(d$rows), logical(1))
  if (all(empty)) stop("no exposure has a usable design SNP", call. = FALSE)

  # ---- residualisation ------------------------------------------------------
  covariate_fit <- NULL
  fit_beta <- panel_beta
  fit_se <- panel_se
  if (method == "residualised") {
    gamma_rows <- match(cov_union, union_keys)
    gamma_rows <- gamma_rows[trait_valid[gamma_rows]]
    if (length(cov_union) > length(gamma_rows)) {
      lost <- length(cov_union) - length(gamma_rows)
      if (strict) {
        stop(lost, " covariate instrument(s) lack valid values for every covariate",
             call. = FALSE)
      }
      warning(lost, " covariate instrument(s) without valid values for every ",
              "covariate were dropped from the covariate -> outcome fits", call. = FALSE)
    }
    if (!is.null(neighbours) && length(gamma_rows)) {
      strength <- apply(abs(trait_beta[gamma_rows, , drop = FALSE] /
                              trait_se[gamma_rows, , drop = FALSE]), 1L, max)
      gamma_rows <- fastmr_mvmr_greedy(gamma_rows, strength, neighbours, no_ld)$kept
    }
    if (length(gamma_rows) < length(traits)) {
      stop("too few covariate instruments to estimate the covariate -> outcome ",
           "effects", call. = FALSE)
    }
    gamma <- fastmr_mvmr_engine(
      c(0L, length(gamma_rows)), gamma_rows, trait_beta[gamma_rows, , drop = FALSE],
      panel_beta, panel_se, covariate_se_model, threads, return_vcov = TRUE
    )
    n_traits <- length(traits)
    G <- matrix(gamma$beta, K, n_traits)                  # outcome x trait
    V <- array(gamma$vcov, c(K, n_traits, n_traits))      # outcome x trait x trait
    adj <- trait_beta %*% t(G)                  # U x K
    extra <- matrix(0, U, K)
    for (a in seq_len(n_traits)) for (b in seq_len(n_traits)) {
      extra <- extra + outer(trait_beta[, a] * trait_beta[, b], V[, a, b])
    }
    fit_beta <- panel_beta - adj
    fit_se <- sqrt(panel_se^2 + extra)
    covariate_fit <- list(
      b = matrix(t(G), n_traits, K, dimnames = list(traits, outcome_labels)),
      se = matrix(sqrt(apply(V, 1L, diag)), n_traits, K,
                  dimnames = list(traits, outcome_labels)),
      nsnp = stats::setNames(gamma$nsnp, outcome_labels),
      instruments = union_keys[gamma_rows],
      se_model = covariate_se_model
    )
    covariate_fit$pval <- fastmr_mvmr_pval(covariate_fit$b, covariate_fit$se)
    rm(adj, extra)
  }
  design_seconds <- unname(proc.time()[["elapsed"]]) - design_started

  # ---- batched fit ----------------------------------------------------------
  estimator_started <- unname(proc.time()[["elapsed"]])
  fit_designs <- designs[!empty]
  cor_list <- NULL
  terms <- colnames(fit_designs[[1L]]$beta)
  if (is.list(exposure_cor) && !is.data.frame(exposure_cor)) {
    missing_cor <- setdiff(names(fit_designs), names(exposure_cor))
    if (length(missing_cor)) {
      stop("exposure_cor list lacks exposure(s): ",
           paste(utils::head(missing_cor, 10L), collapse = ", "), call. = FALSE)
    }
    for (label in names(fit_designs)) fit_designs[[label]]$cor <- exposure_cor[[label]]
    exposure_cor <- NULL
  }
  batch <- fast_mvmr_ivw_batch(
    fit_designs, fit_beta, fit_se, exposure_cor = exposure_cor, se_model = se_model,
    weights = weights, shared_tolerance = shared_tolerance, threads = threads,
    weak_f = weak_f
  )
  estimator_seconds <- unname(proc.time()[["elapsed"]]) - estimator_started

  # ---- assemble -------------------------------------------------------------
  full <- function(m) {
    out <- matrix(NA_real_, E, K, dimnames = list(exposure_labels, outcome_labels))
    out[names(fit_designs), ] <- m
    out
  }
  full_term <- function(name, j) full(batch[[name]][, , j])
  nsnp <- full(batch$nsnp)
  nsnp[is.na(nsnp)] <- 0
  b <- full_term("b", 1L)
  se <- full_term("se", 1L)
  pval <- full_term("pval", 1L)
  cond <- matrix(NA_real_, E, length(terms), dimnames = list(exposure_labels, terms))
  cond[names(fit_designs), ] <- batch$conditional_F
  too_few <- nsnp < minimum_snps
  if (any(too_few)) {
    if (strict) {
      first <- which(too_few, arr.ind = TRUE)[1L, ]
      stop("fewer than minimum_snps for ", exposure_labels[first[1L]], " -> ",
           outcome_labels[first[2L]], " (", nsnp[first[1L], first[2L]], " SNPs)",
           call. = FALSE)
    }
  }
  diagnostics <- data.frame(
    id.exposure = exposure_labels,
    nsnp_design = vapply(designs, function(d) length(d$rows), integer(1)),
    stringsAsFactors = FALSE
  )
  for (j in seq_along(terms)) {
    diagnostics[[paste0("conditional_F_", terms[j])]] <- cond[, j]
  }
  metadata <- list(
    exposure_files = exposure_files, outcome_files = outcome_files,
    covariates = lapply(covariate_sources, function(s) {
      if (inherits(s, "fastmr_store")) unclass(s) else "data.frame"
    }),
    method = method, se_model = se_model, weights = weights,
    instruments = design_keys, diagnostics = diagnostics,
    covariate_fit = covariate_fit,
    counts = list(union_snps = U, masked_cells = excluded$masked,
                  outcome_missing = stats::setNames(missing_outcome, outcome_labels),
                  exposure_dropped = stats::setNames(lost_exposure, exposure_labels),
                  covariate_dropped = stats::setNames(lost_covariate, exposure_labels)),
    shared_deviation = batch$shared_deviation,
    shared_fits = sum(batch$shared, na.rm = TRUE),
    io_threads = io_threads, threads = threads,
    timing = list(io_seconds = io_seconds, design_seconds = design_seconds,
                  estimator_seconds = estimator_seconds,
                  total_seconds = unname(proc.time()[["elapsed"]]) - total_started,
                  source_bytes_read = source_bytes_read)
  )
  Q <- full(batch$Q)
  Q_df <- full(batch$Q_df)
  QA <- full(batch$Q_A)
  method_label <- if (method == "mvmr") "Multivariable IVW" else "Residualised IVW"
  method_code <- if (method == "mvmr") "mv_ivw" else "residualised_ivw"
  if (output_format == "matrix") {
    blank <- function(m) { m[too_few] <- NA_real_; m }
    out <- list(
      b = blank(b), se = blank(se), pval = blank(pval), nsnp = nsnp,
      Q = blank(Q), Q_df = blank(Q_df), Q_pval = blank(full(batch$Q_pval)),
      Q_A = blank(QA), Q_A_pval = blank(full(batch$Q_A_pval)),
      sigma = blank(full(batch$sigma)), conditional_F = cond,
      method = method_label, method_code = method_code
    )
    if (method == "mvmr" && covariate_estimates && length(traits)) {
      out$covariates <- lapply(stats::setNames(seq_along(traits) + 1L, traits), function(j) {
        list(b = blank(full_term("b", j)), se = blank(full_term("se", j)),
             pval = blank(full_term("pval", j)))
      })
    }
    metadata$timing$total_seconds <- unname(proc.time()[["elapsed"]]) - total_started
    attr(out, "mvmr_input") <- metadata
    return(out)
  }
  if (any(too_few) && !strict) {
    pairs <- which(too_few, arr.ind = TRUE)
    fastmr_omission_warning(
      "omitted pair(s) below minimum_snps: ", nrow(pairs),
      function(n) paste0(exposure_labels[pairs[seq_len(n), 1L]], " -> ",
                         outcome_labels[pairs[seq_len(n), 2L]]),
      fastmr_warning_pair_limit())
  }
  flat <- function(m) as.vector(t(m))
  keep <- !flat(too_few)
  result <- data.frame(
    id.exposure = rep(exposure_labels, each = K)[keep],
    id.outcome = rep(outcome_labels, times = E)[keep],
    method = method_label, method_code = method_code,
    nsnp = flat(nsnp)[keep], b = flat(b)[keep], se = flat(se)[keep],
    pval = flat(pval)[keep], Q = flat(Q)[keep], Q_df = flat(Q_df)[keep],
    Q_pval = flat(full(batch$Q_pval))[keep], Q_A = flat(QA)[keep],
    Q_A_pval = flat(full(batch$Q_A_pval))[keep],
    sigma = flat(full(batch$sigma))[keep],
    conditional_F = rep(cond[, 1L], each = K)[keep],
    stringsAsFactors = FALSE
  )
  if (method == "mvmr" && covariate_estimates) {
    for (j in seq_along(traits)) {
      t <- traits[j]
      result[[paste0("b_", t)]] <- flat(full_term("b", j + 1L))[keep]
      result[[paste0("se_", t)]] <- flat(full_term("se", j + 1L))[keep]
      result[[paste0("pval_", t)]] <- flat(full_term("pval", j + 1L))[keep]
    }
  }
  rownames(result) <- NULL
  metadata$timing$total_seconds <- unname(proc.time()[["elapsed"]]) - total_started
  attr(result, "mvmr_input") <- metadata
  fastmr_write_result(result, output)
}
