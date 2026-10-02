# Group number (first-appearance order of each id.exposure/id.outcome pair)
# of every row, the first row of each group, and the per-group ids and labels.
# Missing ids are "".
fastmr_diagnostic_index <- function(data) {
  n <- nrow(data)
  id.exp <- if ("id.exposure" %in% names(data)) as.character(data$id.exposure) else rep("", n)
  id.out <- if ("id.outcome" %in% names(data)) as.character(data$id.outcome) else rep("", n)
  id.exp[is.na(id.exp)] <- ""
  id.out[is.na(id.out)] <- ""
  code.exp <- match(id.exp, unique(id.exp))
  code.out <- match(id.out, unique(id.out))
  key <- (code.out - 1) * (max(code.exp, 0L) + 1) + code.exp
  group <- match(key, unique(key))
  # Groups are numbered by first appearance, so first rows are in group order.
  starts <- which(!duplicated(group))
  label <- function(column, ids) {
    out <- ids[starts]
    if (column %in% names(data)) {
      value <- as.character(data[[column]][starts])
      keep <- !is.na(value)
      out[keep] <- value[keep]
    }
    out
  }
  list(
    group = group,
    count = length(starts),
    starts = starts,
    id.exposure = id.exp[starts],
    id.outcome = id.out[starts],
    exposure = label("exposure", id.exp),
    outcome = label("outcome", id.out)
  )
}

# Row indices (in first-appearance order) of each id.exposure/id.outcome pair,
# with one split() instead of a which() scan per pair.
fastmr_diagnostic_group_index <- function(data) {
  g <- fastmr_diagnostic_index(data)
  rows <- unname(split(seq_len(nrow(data)), factor(g$group, levels = seq_len(g$count))))
  c(list(rows = rows), g[c("id.exposure", "id.outcome", "exposure", "outcome")])
}

fastmr_diagnostic_groups <- function(data) {
  if (!is.data.frame(data)) stop("data must be a data.frame", call. = FALSE)
  # Validate the same required columns and numeric conversions as fast_mr.
  fastmr_prepare_vectors(data)
  g <- fastmr_diagnostic_group_index(data)
  groups <- vector("list", length(g$rows))
  for (i in seq_along(groups)) {
    groups[[i]] <- list(
      data = data[g$rows[[i]], , drop = FALSE],
      id.exposure = g$id.exposure[[i]],
      id.outcome = g$id.outcome[[i]],
      exposure = g$exposure[[i]],
      outcome = g$outcome[[i]]
    )
  }
  groups
}

# Shared set-up for the batched diagnostics: validates `data` and `threads`
# raising the same error the former per-group fast_mr() loop raised first,
# and returns the group index plus the per-row keep flags, numeric vectors and
# SNP ids (NA as ""). Returns NULL when there are no rows. The per-group loop
# checked each group's rows in turn; with `rows_first` (single-SNP and
# leave-one-out) the first group's rows were checked before `threads`.
fastmr_diagnostic_setup <- function(data, threads, rows_first = FALSE) {
  if (!is.data.frame(data)) stop("data must be a data.frame", call. = FALSE)
  prepared <- fastmr_prepare_vectors(data)
  g <- fastmr_diagnostic_index(data)
  if (!g$count) return(NULL)
  n <- nrow(data)
  keep <- if ("mr_keep" %in% names(data)) {
    !is.na(data$mr_keep) & as.logical(data$mr_keep)
  } else {
    rep(TRUE, n)
  }
  valid <- is.finite(prepared$beta.exposure) & is.finite(prepared$beta.outcome) &
    is.finite(prepared$se.exposure) & is.finite(prepared$se.outcome) &
    prepared$se.exposure > 0 & prepared$se.outcome > 0
  snp <- as.character(data$SNP)
  snp[is.na(snp)] <- ""
  bad_value <- keep & !valid
  bad_snp <- keep & !nzchar(snp)
  first_value <- if (any(bad_value)) min(g$group[bad_value]) else Inf
  first_snp <- if (any(bad_snp)) min(g$group[bad_snp]) else Inf
  fail <- function() {
    if (first_value <= first_snp) {
      stop("kept rows must have finite beta values and positive standard errors", call. = FALSE)
    }
    stop("kept rows must have non-empty SNP identifiers", call. = FALSE)
  }
  if (rows_first && min(first_value, first_snp) == 1) fail()
  fastmr_validate_controls(0, NULL, threads)
  if (is.finite(min(first_value, first_snp))) fail()
  g$keep <- keep
  g$prepared <- prepared
  g$snp <- snp
  g
}

# The first SNP row of each id/SNP pair within its group, over all rows
# (kept or not), as the per-group `keep & !duplicated(snp)` selection, in row
# order.
fastmr_diagnostic_selected <- function(g) {
  code <- match(g$snp, unique(g$snp))
  which(g$keep & !duplicated(g$group + g$count * as.numeric(code)))
}

# fastmr_diagnostic_sample_size() of the frames starting at `rows` (NA rows
# give NA).
fastmr_diagnostic_sample_sizes <- function(data, rows) {
  candidates <- c("samplesize.outcome", "samplesize", "sample_size")
  present <- candidates[candidates %in% names(data)]
  if (!length(present)) return(rep(NA_real_, length(rows)))
  value <- suppressWarnings(as.numeric(data[[present[[1L]]]][rows]))
  value[!is.finite(value)] <- NA_real_
  value
}

# A data frame laid out as do.call(rbind, <one-row frames>) returns it:
# attributes in the order names, row.names, class.
fastmr_rbind_layout <- function(columns, row.names = NULL) {
  if (is.null(row.names)) row.names <- .set_row_names(length(columns[[1L]]))
  structure(columns, row.names = row.names, class = "data.frame")
}

# Method display names as the tidy fast_mr() output uses them.
fastmr_method_names <- function(methods) {
  registry <- fastmr_method_registry()
  registry$method[match(methods, registry$code)]
}

#' Calculate TwoSampleMR-compatible heterogeneity statistics
#'
#' @param data A harmonised TwoSampleMR-style data frame.
#' @param methods Heterogeneity-capable fastMR methods. The default matches
#'   native `TwoSampleMR::mr_heterogeneity()` (`ivw` and `egger`).
#' @param threads Maximum native worker count passed to [fast_mr()].
#' @return A tidy data frame with `Q`, `Q_df`, and `Q_pval` per method and pair.
#' @export
fast_mr_heterogeneity <- function(data, methods = c("ivw", "egger"), threads = 1) {
  methods <- fastmr_normalize_methods(methods)
  supported <- c("ivw", "ivw_fe", "ivw_mre", "egger", "uwr")
  unsupported <- setdiff(methods, supported)
  if (length(unsupported)) {
    stop("heterogeneity is not defined for method(s): ",
         paste(unsupported, collapse = ", "), call. = FALSE)
  }
  g <- fastmr_diagnostic_setup(data, threads)
  if (is.null(g)) return(data.frame())
  # One batched (threaded) call; rows are group-major, method-minor.
  result <- fast_mr(data, methods = methods, nboot = 0, threads = threads)
  pair <- rep(seq_len(g$count), each = length(methods))
  fastmr_rbind_layout(list(
    id.exposure = g$id.exposure[pair],
    id.outcome = g$id.outcome[pair],
    outcome = g$outcome[pair],
    exposure = g$exposure[pair],
    method = rep(fastmr_method_names(methods), g$count),
    Q = result$Q,
    Q_df = result$Q_df,
    Q_pval = result$Q_pval
  ))
}

#' Calculate the MR-Egger intercept pleiotropy test
#'
#' @param data A harmonised TwoSampleMR-style data frame.
#' @param threads Maximum native worker count passed to [fast_mr()].
#' @return A tidy data frame compatible with
#'   `TwoSampleMR::mr_pleiotropy_test()`.
#' @export
fast_mr_pleiotropy_test <- function(data, threads = 1) {
  g <- fastmr_diagnostic_setup(data, threads)
  if (is.null(g)) return(data.frame())
  result <- fast_mr(data, methods = "egger", nboot = 0, threads = threads)
  fastmr_rbind_layout(list(
    id.exposure = g$id.exposure,
    id.outcome = g$id.outcome,
    outcome = g$outcome,
    exposure = g$exposure,
    egger_intercept = result$intercept,
    se = result$intercept_se,
    pval = result$intercept_pval
  ))
}

# Assemble the single-SNP / leave-one-out frame. `group`, `samplesize`, `SNP`,
# `b`, `se` and `p` hold the per-SNP rows followed by the per-group summary
# rows; a stable sort by group puts each group's SNP rows (in row order) before
# its summary rows, the order of the former per-group rbind(). `row_label`, if
# given, maps that ordering to explicit row names.
fastmr_diagnostic_snp_frame <- function(g, group, samplesize, SNP, b, se, p,
                                        row_label = NULL) {
  o <- order(group, method = "radix")
  group <- group[o]
  fastmr_rbind_layout(list(
    exposure = g$exposure[group],
    outcome = g$outcome[group],
    id.exposure = g$id.exposure[group],
    id.outcome = g$id.outcome[group],
    samplesize = samplesize[o],
    SNP = SNP[o],
    b = b[o],
    se = se[o],
    p = p[o]
  ), if (is.null(row_label)) NULL else row_label(o))
}

#' Calculate single-SNP MR estimates and aggregate estimates
#'
#' @param data A harmonised TwoSampleMR-style data frame.
#' @param single_method The single-SNP method; defaults to `wald_ratio`.
#' @param all_method Methods used for the aggregate `All - ...` rows.
#' @param threads Maximum native worker count passed to [fast_mr()].
#' @return A tidy data frame compatible with `TwoSampleMR::mr_singlesnp()`.
#' @export
fast_mr_singlesnp <- function(data, single_method = "wald_ratio",
                              all_method = c("ivw", "egger"), threads = 1) {
  single_method <- fastmr_normalize_methods(single_method)
  all_method <- fastmr_normalize_methods(all_method)
  if (length(single_method) != 1L || single_method != "wald_ratio") {
    stop("single_method must be the wald_ratio method", call. = FALSE)
  }
  g <- fastmr_diagnostic_setup(data, threads, rows_first = TRUE)
  if (is.null(g)) return(data.frame())
  selected <- fastmr_diagnostic_selected(g)
  x <- g$prepared$beta.exposure[selected]
  y <- g$prepared$beta.outcome[selected]
  sy <- g$prepared$se.outcome[selected]
  beta <- y / x
  # Match both TwoSampleMR::mr_wald_ratio() and fastMR's native Wald path:
  # the standard error treats the exposure estimate as fixed.
  se <- sy / abs(x)
  p <- rep(NA_real_, length(beta))
  valid <- is.finite(beta) & is.finite(se) & se > 0
  p[valid] <- 2 * stats::pnorm(abs(beta[valid] / se[valid]), lower.tail = FALSE)
  single_group <- g$group[selected]
  # Each group's SNP rows take the sample size of its first selected row.
  first_selected <- selected[match(single_group, single_group)]
  aggregate <- fast_mr(data, methods = all_method, nboot = 0, threads = threads)
  method_count <- length(all_method)
  all_group <- rep(seq_len(g$count), each = method_count)
  fastmr_diagnostic_snp_frame(
    g,
    group = c(single_group, all_group),
    samplesize = c(fastmr_diagnostic_sample_sizes(data, first_selected),
                   fastmr_diagnostic_sample_sizes(data, g$starts)[all_group]),
    SNP = c(g$snp[selected],
            rep(paste("All -", fastmr_method_names(all_method)), g$count)),
    b = c(beta, aggregate$b),
    se = c(se, aggregate$se),
    p = c(p, aggregate$pval))
}

# IVW-family leave-one-out estimates for the selected rows (sorted by group,
# `offsets` delimiting the groups), as the former per-group closed form.
fastmr_leaveoneout_closed_form <- function(prepared, rows, group, offsets, method) {
  x <- prepared$beta.exposure[rows]
  y <- prepared$beta.outcome[rows]
  sy <- prepared$se.outcome[rows]
  if (method == "uwr") {
    w <- rep(1, length(rows))
  } else {
    w <- 1 / (sy * sy)
  }
  wxx <- w * x * x
  wxy <- w * x * y
  wyy <- w * y * y
  # One in-order long-double sum per group, identical to sum() on the group.
  n <- diff(offsets)[group]
  denominator <- fastmr_group_sum_native(offsets, wxx, FALSE)[group] - wxx
  numerator <- fastmr_group_sum_native(offsets, wxy, FALSE)[group] - wxy
  y_sum <- fastmr_group_sum_native(offsets, wyy, FALSE)[group] - wyy
  beta <- rep(NA_real_, length(rows))
  se <- rep(NA_real_, length(rows))
  p <- rep(NA_real_, length(rows))
  valid <- n > 2L & is.finite(denominator) & denominator > 0
  beta[valid] <- numerator[valid] / denominator[valid]
  rss <- y_sum - numerator * numerator / denominator
  df <- n - 2L
  sigma <- sqrt(pmax(0, rss / df))
  base_se <- sqrt(1 / denominator)
  residual_se <- base_se * sigma
  if (method == "ivw_fe") {
    se[valid] <- base_se[valid]
  } else if (method == "ivw_mre") {
    se[valid] <- residual_se[valid]
  } else {
    correction <- pmin(1, sigma)
    se[valid] <- residual_se[valid] / correction[valid]
  }
  valid_p <- valid & is.finite(beta) & is.finite(se) & se > 0
  p[valid_p] <- 2 * stats::pnorm(abs(beta[valid_p] / se[valid_p]), lower.tail = FALSE)
  list(b = beta, se = se, p = p)
}

# Leave-one-out fits of a fast_mr() method for every selected row (sorted by
# group): one batched native call over drop-one jobs on fast_mr()'s own kept,
# de-duplicated layout, so no expanded copy is built.
fastmr_leaveoneout_refit <- function(data, g, rows, method, threads, first_rest) {
  group <- g$group
  kept <- which(g$keep)
  snp_code <- match(g$snp[kept], unique(g$snp[kept]))
  kept <- kept[!duplicated(group[kept] + g$count * as.numeric(snp_code))]
  kept <- kept[order(group[kept], method = "radix")]
  offsets <- c(0L, cumsum(tabulate(group[kept], nbins = g$count)))
  job_group <- group[rows]
  # A selected row is its SNP's first kept row in the group, so it is in `kept`.
  job_drop <- match(rows, kept) - 1L - offsets[job_group]
  # The former code dropped the SNP with `!(SNP == snp)`, which turns rows with
  # an NA SNP into all-NA rows (ids ""). When such a row came first and the
  # pair's ids were not both "", fast_mr() put it in its own empty group and
  # that empty fit was reported.
  raw_snp <- as.character(data$SNP)
  blank <- !nzchar(g$id.exposure) & !nzchar(g$id.outcome)
  empty <- !is.na(first_rest) & is.na(raw_snp[first_rest]) & !blank[job_group]
  job_group[empty] <- NA_integer_
  job_drop[empty] <- -1L
  native <- fastmr_run_groups_drop_native(
    offsets, g$prepared$beta.exposure[kept], g$prepared$beta.outcome[kept],
    g$prepared$se.exposure[kept], g$prepared$se.outcome[kept],
    job_group, as.integer(job_drop), method, threads = as.integer(threads))
  list(b = native$beta, se = native$se, p = native$pval)
}

# First row of each selected row's group once rows with that row's SNP are
# removed (NA when none remain).
fastmr_leaveoneout_first_rest <- function(data, g, rows) {
  raw_snp <- as.character(data$SNP)
  group <- g$group
  first_snp <- raw_snp[g$starts]
  differs <- which(is.na(raw_snp) | raw_snp != first_snp[group])
  first_other <- differs[match(seq_len(g$count), group[differs])]
  job_group <- group[rows]
  drops_first <- !is.na(first_snp[job_group]) & first_snp[job_group] == raw_snp[rows]
  ifelse(drops_first, first_other[job_group], g$starts[job_group])
}

#' Calculate leave-one-SNP-out MR estimates
#'
#' @param data A harmonised TwoSampleMR-style data frame.
#' @param method A leave-one-out-capable method, normally `ivw` or `egger`.
#' @param threads Maximum native worker count passed to [fast_mr()].
#' @return A tidy data frame compatible with `TwoSampleMR::mr_leaveoneout()`.
#' @export
fast_mr_leaveoneout <- function(data, method = "ivw", threads = 1) {
  method <- fastmr_normalize_methods(method)
  if (length(method) != 1L || !method %in% c("ivw", "ivw_fe", "ivw_mre", "egger", "uwr")) {
    stop("method must be one heterogeneity-capable regression method", call. = FALSE)
  }
  g <- fastmr_diagnostic_setup(data, threads, rows_first = TRUE)
  if (is.null(g)) return(data.frame())
  selected <- fastmr_diagnostic_selected(g)
  rows <- selected[order(g$group[selected], method = "radix")]
  row_group <- g$group[rows]
  closed_form <- method != "egger"
  if (closed_form) {
    offsets <- c(0L, cumsum(tabulate(row_group, nbins = g$count)))
    leave <- fastmr_leaveoneout_closed_form(g$prepared, rows, row_group, offsets, method)
    # Each group's rows take the sample size of its first selected row.
    leave_size <- fastmr_diagnostic_sample_sizes(data, rows[offsets[row_group] + 1L])
  } else {
    first_rest <- fastmr_leaveoneout_first_rest(data, g, rows)
    leave <- fastmr_leaveoneout_refit(data, g, rows, method, threads, first_rest)
    # With no rows left the former fast_mr() call returned no rows at all.
    none <- is.na(first_rest)
    leave$b[none] <- NA_real_
    leave$se[none] <- NA_real_
    leave$p[none] <- NA_real_
    # The former code took the sample size of the first remaining row; an NA
    # SNP row there had become all-NA.
    size_row <- first_rest
    size_row[!is.na(first_rest) & is.na(as.character(data$SNP)[first_rest])] <- NA_integer_
    leave_size <- fastmr_diagnostic_sample_sizes(data, size_row)
  }
  all <- fast_mr(data, methods = method, nboot = 0, threads = threads)
  row_label <- NULL
  if (closed_form) {
    # The former code rbind()-ed one-row slices of each group's frame, whose
    # row names were their positions in it, with fresh one-row "All" frames;
    # reproduce rbind.data.frame()'s labels for that sequence.
    position <- c(seq_along(rows) - offsets[row_group], rep(1L, g$count))
    row_label <- function(o) {
      label <- position[o]
      first <- match(TRUE, label != 1L)
      if (is.na(first)) return(NULL)
      label <- c(seq_len(first - 1L), label[first:length(label)])
      make.unique(as.character(label), sep = "")
    }
  }
  fastmr_diagnostic_snp_frame(
    g,
    group = c(row_group, seq_len(g$count)),
    samplesize = c(leave_size, fastmr_diagnostic_sample_sizes(data, g$starts)),
    SNP = c(g$snp[rows], rep("All", g$count)),
    b = c(leave$b, all$b),
    se = c(leave$se, all$se),
    p = c(leave$p, all$pval),
    row_label = row_label)
}
