#' Run exact summary-statistics Mendelian randomization
#'
#' @param data A data frame with `beta.exposure`, `beta.outcome`,
#'   `se.exposure`, `se.outcome`, and `SNP`, with optional `id.exposure` and
#'   `id.outcome` columns.
#' @param methods Character vector of method codes. See
#'   [fastmr_method_registry()].
#' @param nboot Number of normal bootstrap draws for median and mode methods.
#' @param seed Optional integer seed. Seeded median/mode methods share one
#'   ratio bootstrap layout per pair.
#' @param threads Maximum native worker count. It is most useful for
#'   [fast_mr_grid()]; single-pair calls remain bounded and deterministic.
#' @param output Optional path for a Zstandard-compressed Parquet copy of the
#'   result. The path must not already exist; use [fast_write_parquet()] when
#'   an overwrite or another compression codec is required.
#' @param ... Optional `phi` bandwidth multiplier for mode methods and `penk`
#'   penalty multiplier for penalised weighted median (default 20).
#' @return A tidy data frame using TwoSampleMR-compatible result columns.
#' @export
fast_mr <- function(data,
                    methods = c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode"),
                    nboot = 1000,
                    seed = NULL,
                    threads = 1,
                    output = NULL,
                    ...) {
  if (!is.data.frame(data)) stop("data must be a data.frame", call. = FALSE)
  controls <- fastmr_validate_controls(nboot, seed, threads)
  methods <- fastmr_normalize_methods(methods)
  dots <- list(...)
  unknown_dots <- setdiff(names(dots), c("phi", "penk"))
  if (length(unknown_dots)) stop("unknown option(s): ", paste(unknown_dots, collapse = ", "), call. = FALSE)
  phi <- if (is.null(dots$phi)) 1 else dots$phi
  if (length(phi) != 1L || !is.finite(phi) || phi <= 0) stop("phi must be positive and finite", call. = FALSE)
  penk <- if (is.null(dots$penk)) 20 else dots$penk
  if (length(penk) != 1L || !is.finite(penk) || penk <= 0) stop("penk must be positive and finite", call. = FALSE)
  prepared <- fastmr_prepare_vectors(data)
  n <- nrow(prepared)
  keep <- if ("mr_keep" %in% names(data)) !is.na(data$mr_keep) & as.logical(data$mr_keep) else rep(TRUE, n)
  valid <- is.finite(prepared$beta.exposure) & is.finite(prepared$beta.outcome) &
    is.finite(prepared$se.exposure) & is.finite(prepared$se.outcome) &
    prepared$se.exposure > 0 & prepared$se.outcome > 0
  if (any(keep & !valid)) {
    stop("kept rows must have finite beta values and positive standard errors", call. = FALSE)
  }
  snp <- as.character(data$SNP)
  snp[is.na(snp)] <- ""
  if (any(keep & !nzchar(snp))) stop("kept rows must have non-empty SNP identifiers", call. = FALSE)
  id.exp <- if ("id.exposure" %in% names(data)) as.character(data$id.exposure) else rep("", n)
  id.out <- if ("id.outcome" %in% names(data)) as.character(data$id.outcome) else rep("", n)
  id.exp[is.na(id.exp)] <- ""
  id.out[is.na(id.out)] <- ""
  groups <- unique(data.frame(id.exposure = id.exp, id.outcome = id.out,
                              stringsAsFactors = FALSE))
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    group_index <- which(id.exp == groups$id.exposure[[i]] &
                         id.out == groups$id.outcome[[i]])
    index <- group_index[keep[group_index]]
    # Joins and multi-study exports often repeat the same SNP row. Count each
    # SNP once per MR pair; retain the first row deterministically. Repeated
    # p-values are metadata and do not affect this rule.
    index <- index[!duplicated(snp[index])]
    representative <- group_index[[1L]]
    native <- fastmr_native_call(
      fastmr_run_native,
      list(
        exposure_beta = prepared[["beta.exposure"]][index],
        outcome_beta = prepared[["beta.outcome"]][index],
        exposure_se = prepared[["se.exposure"]][index],
        outcome_se = prepared[["se.outcome"]][index],
        methods = methods, nboot = controls[["nboot"]], seed = NULL,
        threads = controls[["threads"]], phi = phi, penk = penk
      ),
      if (is.null(controls[["seed"]])) NULL else controls[["seed"]] + i - 1L
    )
    label.exp <- if ("exposure" %in% names(data)) as.character(data$exposure[representative]) else id.exp[representative]
    label.out <- if ("outcome" %in% names(data)) as.character(data$outcome[representative]) else id.out[representative]
    rows[[i]] <- fastmr_tidy_native(native, methods, id.exp[representative], id.out[representative],
                                     exposure_label = label.exp, outcome_label = label.out)
  }
  if (!length(rows)) return(fastmr_write_result(fastmr_tidy_native(list(), methods), output))
  fastmr_write_result(do.call(rbind, rows), output)
}

#' Run every exposure/outcome pair in a shared exact grid
#'
#' Matrix rows are exposures/outcomes and columns are shared SNPs. Results are
#' returned in exposure-major, outcome-minor order. R's column-major matrices
#' are copied once at the C++ boundary into contiguous row-major pair layouts.
#' @param exposure_beta Exposure effect matrix, exposures by SNP.
#' @param outcome_beta Outcome effect matrix, outcomes by SNP.
#' @param exposure_se Exposure standard-error matrix.
#' @param outcome_se Outcome standard-error matrix.
#' @param methods Character vector of method codes.
#' @param nboot Number of bootstrap draws.
#' @param seed Optional integer seed.
#' @param threads Maximum native worker count.
#' @param output Optional path for a Zstandard-compressed Parquet copy of the
#'   result. The path must not already exist; use [fast_write_parquet()] when
#'   an overwrite or another compression codec is required.
#' @param ... Optional `phi` bandwidth multiplier for mode methods and `penk`
#'   penalty multiplier for penalised weighted median (default 20).
#' @return A tidy data frame with one row per method and grid pair.
#' @export
fast_mr_grid <- function(exposure_beta, outcome_beta, exposure_se, outcome_se,
                         methods = c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode"),
                         nboot = 1000, seed = NULL, threads = 1, output = NULL, ...) {
  controls <- fastmr_validate_controls(nboot, seed, threads)
  methods <- fastmr_normalize_methods(methods)
  dots <- list(...)
  unknown_dots <- setdiff(names(dots), c("phi", "penk"))
  if (length(unknown_dots)) stop("unknown option(s): ", paste(unknown_dots, collapse = ", "), call. = FALSE)
  phi <- if (is.null(dots$phi)) 1 else dots$phi
  if (length(phi) != 1L || !is.finite(phi) || phi <= 0) stop("phi must be positive and finite", call. = FALSE)
  penk <- if (is.null(dots$penk)) 20 else dots$penk
  if (length(penk) != 1L || !is.finite(penk) || penk <= 0) stop("penk must be positive and finite", call. = FALSE)
  arrays <- Map(fastmr_matrix_numeric,
                list(exposure_beta, outcome_beta, exposure_se, outcome_se),
                c("exposure_beta", "outcome_beta", "exposure_se", "outcome_se"))
  names(arrays) <- c("exposure_beta", "outcome_beta", "exposure_se", "outcome_se")
  if (nrow(arrays$exposure_beta) == 0L || nrow(arrays$outcome_beta) == 0L ||
      ncol(arrays$exposure_beta) == 0L ||
      any(!is.finite(arrays$exposure_beta)) || any(!is.finite(arrays$outcome_beta)) ||
      any(!is.finite(arrays$exposure_se)) || any(!is.finite(arrays$outcome_se)) ||
      any(arrays$exposure_se <= 0) || any(arrays$outcome_se <= 0)) {
    stop("grid inputs must be non-empty with finite beta values and positive standard errors", call. = FALSE)
  }
  exp.snps <- colnames(arrays$exposure_beta)
  out.snps <- colnames(arrays$outcome_beta)
  if (xor(is.null(exp.snps), is.null(out.snps)) ||
      (!is.null(exp.snps) && !identical(exp.snps, out.snps))) {
    stop("exposure and outcome matrices must use the same SNP column names and order", call. = FALSE)
  }
  native <- fastmr_native_call(
    fastmr_grid_native,
    list(
      exposure_beta = arrays[["exposure_beta"]],
      outcome_beta = arrays[["outcome_beta"]],
      exposure_se = arrays[["exposure_se"]],
      outcome_se = arrays[["outcome_se"]],
      methods = methods, nboot = controls[["nboot"]], seed = NULL,
      threads = controls[["threads"]], phi = phi, penk = penk
    ),
    controls[["seed"]]
  )
  exp.labels <- rownames(arrays$exposure_beta)
  out.labels <- rownames(arrays$outcome_beta)
  if (is.null(exp.labels)) exp.labels <- as.character(seq_len(nrow(arrays$exposure_beta)))
  if (is.null(out.labels)) out.labels <- as.character(seq_len(nrow(arrays$outcome_beta)))
  fastmr_write_result(fastmr_tidy_grid_native(native, methods, exp.labels, out.labels), output)
}
