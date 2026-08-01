#' List the fastMR method registry
#'
#' @return A data frame mapping short method codes to tidy result names and
#'   descriptions.
#' @export
fastmr_method_registry <- function() {
  data.frame(
    code = c("ivw", "ivw_fe", "ivw_mre", "egger", "weighted_median",
             "simple_mode", "weighted_mode", "wald_ratio"),
    method = c(
      "Inverse variance weighted",
      "Inverse variance weighted (fixed effects)",
      "Inverse variance weighted (multiplicative random effects)",
      "MR Egger",
      "Weighted median",
      "Simple mode",
      "Weighted mode",
      "Wald ratio"
    ),
    description = c(
      "Multiplicative random-effects IVW with under-dispersion correction",
      "Fixed-effects IVW standard error",
      "Multiplicative random-effects IVW without under-dispersion correction",
      "Weighted Egger regression with an intercept",
      "Weighted median of delta-method Wald ratios",
      "Unweighted kernel mode of Wald ratios",
      "Ratio-SE-weighted kernel mode of Wald ratios",
      "Single-SNP Wald ratio"
    ),
    stringsAsFactors = FALSE
  )
}

fastmr_normalize_methods <- function(methods) {
  if (length(methods) == 0L) stop("methods must contain at least one method", call. = FALSE)
  if (!is.character(methods)) stop("methods must be character names", call. = FALSE)
  aliases <- c(
    mr_ivw = "ivw",
    mr_ivw_fe = "ivw_fe",
    mr_ivw_mre = "ivw_mre",
    mr_egger_regression = "egger",
    mr_weighted_median = "weighted_median",
    mr_simple_mode = "simple_mode",
    mr_weighted_mode = "weighted_mode",
    mr_wald_ratio = "wald_ratio",
    `Inverse variance weighted` = "ivw",
    `MR Egger` = "egger",
    `Weighted median` = "weighted_median",
    `Simple mode` = "simple_mode",
    `Weighted mode` = "weighted_mode",
    `Wald ratio` = "wald_ratio"
  )
  normalized <- unname(ifelse(methods %in% names(aliases), aliases[methods], methods))
  allowed <- fastmr_method_registry()$code
  unsupported <- setdiff(normalized, allowed)
  if (length(unsupported)) {
    stop("unknown MR method(s): ", paste(unsupported, collapse = ", "), call. = FALSE)
  }
  normalized
}

fastmr_validate_controls <- function(nboot, seed, threads) {
  if (length(nboot) != 1L || is.na(nboot) || !is.finite(nboot) || nboot < 0 || nboot != floor(nboot)) {
    stop("nboot must be one non-negative integer", call. = FALSE)
  }
  if (length(threads) != 1L || is.na(threads) || !is.finite(threads) || threads < 1 || threads != floor(threads)) {
    stop("threads must be one positive integer", call. = FALSE)
  }
  if (!is.null(seed)) {
    if (length(seed) != 1L || is.na(seed) || !is.finite(seed) || seed != floor(seed)) {
      stop("seed must be NULL or one finite integer", call. = FALSE)
    }
  }
  invisible(list(nboot = as.integer(nboot), seed = if (is.null(seed)) NULL else as.numeric(seed),
                 threads = as.integer(threads)))
}

fastmr_numeric <- function(x, name) {
  if (is.factor(x)) x <- as.character(x)
  converted <- suppressWarnings(as.numeric(x))
  invalid <- is.na(converted) & !is.na(x)
  if (any(invalid)) stop(name, " must be numeric", call. = FALSE)
  converted
}

fastmr_prepare_vectors <- function(data) {
  required <- c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  data.frame(
    beta.exposure = fastmr_numeric(data$beta.exposure, "beta.exposure"),
    beta.outcome = fastmr_numeric(data$beta.outcome, "beta.outcome"),
    se.exposure = fastmr_numeric(data$se.exposure, "se.exposure"),
    se.outcome = fastmr_numeric(data$se.outcome, "se.outcome"),
    stringsAsFactors = FALSE
  )
}

fastmr_matrix_numeric <- function(x, name) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) stop(name, " must be a matrix", call. = FALSE)
  if (!is.numeric(x)) {
    original <- x
    converted <- suppressWarnings(as.numeric(original))
    if (any(is.na(converted) & !is.na(original))) stop(name, " must be numeric", call. = FALSE)
    dim(converted) <- dim(original)
    dimnames(converted) <- dimnames(original)
    x <- converted
  }
  storage.mode(x) <- "double"
  x
}

fastmr_scalar <- function(x, name, default = NA_real_) {
  if (is.null(x[[name]]) || length(x[[name]]) == 0L) return(default)
  x[[name]][[1L]]
}

fastmr_tidy_native <- function(native_results, methods, id.exposure = "", id.outcome = "",
                               exposure_index = NULL, outcome_index = NULL,
                               exposure_label = id.exposure, outcome_label = id.outcome) {
  registry <- fastmr_method_registry()
  code <- vapply(native_results, fastmr_scalar, character(1), name = "method", default = "")
  display <- registry$method[match(code, registry$code)]
  n <- vapply(native_results, fastmr_scalar, numeric(1), name = "n", default = NA_real_)
  out <- data.frame(
    id.exposure = rep(exposure_label, length(native_results)),
    id.outcome = rep(outcome_label, length(native_results)),
    method = display,
    method_code = code,
    nsnp = n,
    b = vapply(native_results, fastmr_scalar, numeric(1), name = "beta"),
    se = vapply(native_results, fastmr_scalar, numeric(1), name = "se"),
    pval = vapply(native_results, fastmr_scalar, numeric(1), name = "pval"),
    Q = vapply(native_results, fastmr_scalar, numeric(1), name = "Q"),
    Q_df = vapply(native_results, fastmr_scalar, numeric(1), name = "Q_df"),
    Q_pval = vapply(native_results, fastmr_scalar, numeric(1), name = "Q_pval"),
    sigma = vapply(native_results, fastmr_scalar, numeric(1), name = "sigma"),
    intercept = vapply(native_results, fastmr_scalar, numeric(1), name = "intercept"),
    intercept_se = vapply(native_results, fastmr_scalar, numeric(1), name = "intercept_se"),
    intercept_pval = vapply(native_results, fastmr_scalar, numeric(1), name = "intercept_pval"),
    ratio_se_mean = vapply(native_results, fastmr_scalar, numeric(1), name = "ratio_se_mean"),
    bootstrap = vapply(native_results, fastmr_scalar, numeric(1), name = "bootstrap"),
    phi = vapply(native_results, fastmr_scalar, numeric(1), name = "phi"),
    flipped = vapply(native_results, fastmr_scalar, numeric(1), name = "flipped"),
    se_exposure_mean = vapply(native_results, fastmr_scalar, numeric(1), name = "se_exposure_mean"),
    stringsAsFactors = FALSE
  )
  if (!is.null(exposure_index)) out$exposure_index <- exposure_index
  if (!is.null(outcome_index)) out$outcome_index <- outcome_index
  out
}


fastmr_native_call <- function(native, args, seed) {
  if (is.null(seed)) return(do.call(native, args))
  state_env <- .GlobalEnv
  had_state <- exists(".Random.seed", envir = state_env, inherits = FALSE)
  old_state <- if (had_state) get(".Random.seed", envir = state_env, inherits = FALSE) else NULL
  on.exit({
    if (had_state) {
      assign(".Random.seed", old_state, envir = state_env)
    } else if (exists(".Random.seed", envir = state_env, inherits = FALSE)) {
      rm(".Random.seed", envir = state_env)
    }
  }, add = TRUE)
  set.seed(seed)
  args[["seed"]] <- NULL
  do.call(native, args)
}
