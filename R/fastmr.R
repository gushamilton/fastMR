#' Run exact summary-statistics Mendelian randomization
#'
#' @param data A data frame with `beta.exposure`, `beta.outcome`,
#'   `se.exposure`, `se.outcome`, and optionally `SNP`, `id.exposure`, and
#'   `id.outcome` columns.
#' @param methods Character vector of method codes. See
#'   [fastmr_method_registry()].
#' @param nboot Number of normal bootstrap draws for median and mode methods.
#' @param seed Optional integer seed. Seeded median/mode methods share one
#'   ratio bootstrap layout per pair.
#' @param threads Maximum native worker count. It is most useful for
#'   [fast_mr_grid()]; single-pair calls remain bounded and deterministic.
#' @param ... Optional `phi` bandwidth multiplier for mode methods.
#' @return A tidy data frame using TwoSampleMR-compatible result columns.
#' @export
fast_mr <- function(data,
                    methods = c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode"),
                    nboot = 1000,
                    seed = NULL,
                    threads = 1,
                    ...) {
  if (!is.data.frame(data)) stop("data must be a data.frame", call. = FALSE)
  controls <- fastmr_validate_controls(nboot, seed, threads)
  methods <- fastmr_normalize_methods(methods)
  dots <- list(...)
  unknown_dots <- setdiff(names(dots), "phi")
  if (length(unknown_dots)) stop("unknown option(s): ", paste(unknown_dots, collapse = ", "), call. = FALSE)
  phi <- if (is.null(dots$phi)) 1 else dots$phi
  if (length(phi) != 1L || !is.finite(phi) || phi <= 0) stop("phi must be positive and finite", call. = FALSE)
  prepared <- fastmr_prepare_vectors(data)
  n <- nrow(prepared)
  id.exp <- if ("id.exposure" %in% names(data)) as.character(data$id.exposure) else rep("", n)
  id.out <- if ("id.outcome" %in% names(data)) as.character(data$id.outcome) else rep("", n)
  id.exp[is.na(id.exp)] <- ""
  id.out[is.na(id.out)] <- ""
  keys <- paste(id.exp, id.out, sep = "\r")
  groups <- unique(keys)
  rows <- vector("list", length(groups))
  keep <- if ("mr_keep" %in% names(data)) !is.na(data$mr_keep) & as.logical(data$mr_keep) else rep(TRUE, n)
  for (i in seq_along(groups)) {
    index <- which(keys == groups[[i]] & keep)
    representative <- which(keys == groups[[i]])[[1L]]
    native <- fastmr_native_call(
      fastmr_run_native,
      list(
        exposure_beta = prepared[["beta.exposure"]][index],
        outcome_beta = prepared[["beta.outcome"]][index],
        exposure_se = prepared[["se.exposure"]][index],
        outcome_se = prepared[["se.outcome"]][index],
        methods = methods, nboot = controls[["nboot"]], seed = NULL,
        threads = controls[["threads"]], phi = phi
      ),
      controls[["seed"]]
    )
    label.exp <- if ("exposure" %in% names(data)) as.character(data$exposure[representative]) else id.exp[representative]
    label.out <- if ("outcome" %in% names(data)) as.character(data$outcome[representative]) else id.out[representative]
    rows[[i]] <- fastmr_tidy_native(native, methods, id.exp[representative], id.out[representative],
                                     exposure_label = label.exp, outcome_label = label.out)
  }
  if (!length(rows)) return(fastmr_tidy_native(list(), methods))
  do.call(rbind, rows)
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
#' @param ... Optional `phi` bandwidth multiplier for mode methods.
#' @return A tidy data frame with one row per method and grid pair.
#' @export
fast_mr_grid <- function(exposure_beta, outcome_beta, exposure_se, outcome_se,
                         methods = c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode"),
                         nboot = 1000, seed = NULL, threads = 1, ...) {
  controls <- fastmr_validate_controls(nboot, seed, threads)
  methods <- fastmr_normalize_methods(methods)
  dots <- list(...)
  unknown_dots <- setdiff(names(dots), "phi")
  if (length(unknown_dots)) stop("unknown option(s): ", paste(unknown_dots, collapse = ", "), call. = FALSE)
  phi <- if (is.null(dots$phi)) 1 else dots$phi
  if (length(phi) != 1L || !is.finite(phi) || phi <= 0) stop("phi must be positive and finite", call. = FALSE)
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
  native <- fastmr_native_call(
    fastmr_grid_native,
    list(
      exposure_beta = arrays[["exposure_beta"]],
      outcome_beta = arrays[["outcome_beta"]],
      exposure_se = arrays[["exposure_se"]],
      outcome_se = arrays[["outcome_se"]],
      methods = methods, nboot = controls[["nboot"]], seed = NULL,
      threads = controls[["threads"]], phi = phi
    ),
    controls[["seed"]]
  )
  exp.labels <- rownames(arrays$exposure_beta)
  out.labels <- rownames(arrays$outcome_beta)
  if (is.null(exp.labels)) exp.labels <- as.character(seq_len(nrow(arrays$exposure_beta)))
  if (is.null(out.labels)) out.labels <- as.character(seq_len(nrow(arrays$outcome_beta)))
  fastmr_tidy_grid_native(native, methods, exp.labels, out.labels)
}
