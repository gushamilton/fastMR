#' List the fastMR method registry
#'
#' @return A data frame mapping short method codes to tidy result names and
#'   descriptions.
#' @export
fastmr_method_registry <- function() {
  data.frame(
    code = c("ivw", "ivw_fe", "ivw_mre", "egger", "egger_bootstrap", "uwr",
             "sign", "simple_median", "weighted_median", "penalised_weighted_median", "simple_mode",
             "weighted_mode", "wald_ratio"),
    method = c(
      "Inverse variance weighted",
      "Inverse variance weighted (fixed effects)",
      "Inverse variance weighted (multiplicative random effects)",
      "MR Egger",
      "MR Egger (bootstrap)",
      "Unweighted regression",
      "Sign concordance test",
      "Simple median",
      "Weighted median",
      "Penalised weighted median",
      "Simple mode",
      "Weighted mode",
      "Wald ratio"
    ),
    description = c(
      "Multiplicative random-effects IVW with under-dispersion correction",
      "Fixed-effects IVW standard error",
      "Multiplicative random-effects IVW without under-dispersion correction",
      "Weighted Egger regression with an intercept",
      "Parametric bootstrap MR-Egger regression",
      "Unweighted no-intercept regression",
      "Exact binomial sign concordance test",
      "Unweighted median of Wald ratios",
      "Weighted median of delta-method Wald ratios",
      "Penalised weighted median with chi-square down-weighting",
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
    mr_egger_regression_bootstrap = "egger_bootstrap",
    mr_uwr = "uwr",
    mr_sign = "sign",
    mr_simple_median = "simple_median",
    mr_weighted_median = "weighted_median",
    mr_penalised_weighted_median = "penalised_weighted_median",
    mr_simple_mode = "simple_mode",
    mr_weighted_mode = "weighted_mode",
    mr_wald_ratio = "wald_ratio",
    `Inverse variance weighted` = "ivw",
    `MR Egger` = "egger",
    `MR Egger (bootstrap)` = "egger_bootstrap",
    `Unweighted regression` = "uwr",
    `Sign concordance test` = "sign",
    `Simple median` = "simple_median",
    `Weighted median` = "weighted_median",
    `Penalised weighted median` = "penalised_weighted_median",
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
  if (anyDuplicated(normalized)) {
    stop("methods must be unique after alias normalization", call. = FALSE)
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
    if (abs(seed) > .Machine$integer.max) {
      stop("seed must lie in [-", .Machine$integer.max, ", ", .Machine$integer.max, "]", call. = FALSE)
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
    id.exposure = rep(id.exposure, length(native_results)),
    id.outcome = rep(id.outcome, length(native_results)),
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

# Tidy rows for pairs `first:last` (global, exposure-major pair numbers) of a
# compact native result whose first column is global pair `pair_offset + 1`.
# With the defaults this converts the whole result.
fastmr_tidy_compact <- function(native_results, methods, exposure_labels, outcome_labels,
                                first = 1L, last = NULL, pair_offset = 0) {
  method_codes <- as.character(native_results$methods)
  method_count <- length(method_codes)
  outcome_count <- length(outcome_labels)
  pair_count <- length(exposure_labels) * outcome_count
  if (is.null(last)) last <- pair_count
  whole <- identical(first, 1L) && last == pair_count && pair_offset == 0
  npairs <- last - first + 1
  total <- method_count * npairs
  pairs <- seq.int(as.integer(first), as.integer(last))
  pair_index <- rep(pairs, each = method_count)
  exposure_index <- ((pair_index - 1L) %/% outcome_count) + 1L
  outcome_index <- ((pair_index - 1L) %% outcome_count) + 1L
  registry <- fastmr_method_registry()
  code <- rep(method_codes, times = npairs)
  cols <- pairs - as.integer(pair_offset)
  values <- function(name, default = NA_real_) {
    if (is.null(native_results[[name]])) return(rep(default, total))
    m <- native_results[[name]]
    if (whole) as.vector(m) else as.vector(m[, cols, drop = FALSE])
  }
  nsnp <- if (!is.null(native_results$nsnp)) values("nsnp") else
    rep(as.numeric(native_results$n)[1L], total)
  out <- data.frame(
    id.exposure = exposure_labels[exposure_index],
    id.outcome = outcome_labels[outcome_index],
    method = registry$method[match(code, registry$code)],
    method_code = code,
    nsnp = nsnp,
    b = values("beta"),
    se = values("se"),
    pval = values("pval"),
    Q = values("Q"),
    Q_df = values("Q_df"),
    Q_pval = values("Q_pval"),
    sigma = values("sigma"),
    intercept = values("intercept"),
    intercept_se = values("intercept_se"),
    intercept_pval = values("intercept_pval"),
    ratio_se_mean = values("ratio_se_mean"),
    bootstrap = values("bootstrap"),
    phi = values("phi"),
    flipped = values("flipped"),
    se_exposure_mean = values("se_exposure_mean"),
    exposure_index = exposure_index,
    outcome_index = outcome_index,
    stringsAsFactors = FALSE
  )
  out
}

fastmr_tidy_grid_native <- function(native_results, methods, exposure_labels, outcome_labels) {
  if (!length(native_results)) return(fastmr_tidy_native(list(), methods))
  if (inherits(native_results, "fastmr_ivw_compact") ||
      inherits(native_results, "fastmr_grid_compact")) {
    return(fastmr_tidy_compact(native_results, methods, exposure_labels, outcome_labels))
  }
  flat <- unlist(native_results, recursive = FALSE, use.names = FALSE)
  method_count <- length(methods)
  pair_index <- rep(seq_along(native_results), each = method_count)
  exposure_count <- length(exposure_labels)
  outcome_count <- length(outcome_labels)
  exposure_index <- ((pair_index - 1L) %/% outcome_count) + 1L
  outcome_index <- ((pair_index - 1L) %% outcome_count) + 1L
  registry <- fastmr_method_registry()
  code <- vapply(flat, fastmr_scalar, character(1), name = "method", default = "")
  display <- registry$method[match(code, registry$code)]
  out <- data.frame(
    id.exposure = exposure_labels[exposure_index],
    id.outcome = outcome_labels[outcome_index],
    method = display,
    method_code = code,
    nsnp = vapply(flat, fastmr_scalar, numeric(1), name = "n", default = NA_real_),
    b = vapply(flat, fastmr_scalar, numeric(1), name = "beta"),
    se = vapply(flat, fastmr_scalar, numeric(1), name = "se"),
    pval = vapply(flat, fastmr_scalar, numeric(1), name = "pval"),
    Q = vapply(flat, fastmr_scalar, numeric(1), name = "Q"),
    Q_df = vapply(flat, fastmr_scalar, numeric(1), name = "Q_df"),
    Q_pval = vapply(flat, fastmr_scalar, numeric(1), name = "Q_pval"),
    sigma = vapply(flat, fastmr_scalar, numeric(1), name = "sigma"),
    intercept = vapply(flat, fastmr_scalar, numeric(1), name = "intercept"),
    intercept_se = vapply(flat, fastmr_scalar, numeric(1), name = "intercept_se"),
    intercept_pval = vapply(flat, fastmr_scalar, numeric(1), name = "intercept_pval"),
    ratio_se_mean = vapply(flat, fastmr_scalar, numeric(1), name = "ratio_se_mean"),
    bootstrap = vapply(flat, fastmr_scalar, numeric(1), name = "bootstrap"),
    phi = vapply(flat, fastmr_scalar, numeric(1), name = "phi"),
    flipped = vapply(flat, fastmr_scalar, numeric(1), name = "flipped"),
    se_exposure_mean = vapply(flat, fastmr_scalar, numeric(1), name = "se_exposure_mean"),
    exposure_index = exposure_index,
    outcome_index = outcome_index,
    stringsAsFactors = FALSE
  )
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

# Integer group codes for (id.exp, id.out) pairs, numbered by first appearance;
# attr "n" holds the number of groups.
fastmr_group_ids <- function(id.exp, id.out) {
  a <- match(id.exp, unique(id.exp))
  b <- match(id.out, unique(id.out))
  key <- a + as.numeric(max(a, 0L)) * (b - 1)
  gid <- match(key, unique(key))
  attr(gid, "n") <- length(unique(key))
  gid
}

# Methods whose native code draws random numbers when nboot > 0.
fastmr_methods_use_rng <- function(methods, nboot) {
  nboot > 0 && any(methods %in% c("simple_median", "weighted_median",
                                  "penalised_weighted_median", "egger_bootstrap",
                                  "simple_mode", "weighted_mode"))
}

# Vectorised equivalent of rbind-ing fastmr_tidy_native() per group, from the
# flat group-major / method-minor output of fastmr_run_groups_native().
fastmr_tidy_groups_native <- function(native, method_count, id.exposure, id.outcome) {
  registry <- fastmr_method_registry()
  code <- native$method
  out <- data.frame(
    id.exposure = rep(id.exposure, each = method_count),
    id.outcome = rep(id.outcome, each = method_count),
    method = registry$method[match(code, registry$code)],
    method_code = code,
    nsnp = native$n,
    b = native$beta,
    se = native$se,
    pval = native$pval,
    Q = native$Q,
    Q_df = native$Q_df,
    Q_pval = native$Q_pval,
    sigma = native$sigma,
    intercept = native$intercept,
    intercept_se = native$intercept_se,
    intercept_pval = native$intercept_pval,
    ratio_se_mean = native$ratio_se_mean,
    bootstrap = native$bootstrap,
    phi = native$phi,
    flipped = native$flipped,
    se_exposure_mean = native$se_exposure_mean,
    stringsAsFactors = FALSE
  )
  # rbind() of per-group frames yields attributes in the order names,
  # row.names, class; reproduce it so serialize() output is byte-identical.
  attrs <- attributes(out)
  attributes(out) <- NULL
  names(out) <- attrs$names
  attr(out, "row.names") <- .set_row_names(length(attrs$row.names))  # compact c(NA, -n), as rbind
  class(out) <- attrs$class
  out
}

# ---- compact (lightweight) grid results and chunked tidy conversion --------

fastmr_new_compact_grid <- function(native, methods, exposure_labels, outcome_labels) {
  structure(
    list(native = native, methods = methods,
         exposure_labels = exposure_labels, outcome_labels = outcome_labels),
    class = "fastmr_compact_grid"
  )
}

fastmr_compact_pairs <- function(x) length(x$exposure_labels) * length(x$outcome_labels)

# Iterator over tidy chunks of `chunk_pairs` pairs of a compact grid.
fastmr_compact_chunks <- function(x, chunk_pairs) {
  total <- fastmr_compact_pairs(x)
  start <- 1
  function() {
    if (start > total) return(NULL)
    end <- min(total, start + chunk_pairs - 1)
    out <- fastmr_grid_chunk(x, start, end)
    start <<- end + 1
    out
  }
}

# Iterator for IVW-only grids: native kernel per block of whole exposures
# (about `chunk_pairs` pairs), tidied in sub-chunks of `chunk_pairs`.
fastmr_ivw_blocks <- function(arrays, run_native, exposure_labels, outcome_labels,
                              methods, chunk_pairs) {
  O <- length(outcome_labels)
  E <- length(exposure_labels)
  block <- max(1, ceiling(chunk_pairs / O))
  next_exp <- 1
  cur <- NULL; cur_offset <- 0; cur_next <- 1; cur_last <- 0
  function() {
    if (is.null(cur) || cur_next > cur_last) {
      if (next_exp > E) return(NULL)
      e1 <- min(E, next_exp + block - 1)
      idx <- next_exp:e1
      cur <<- run_native(arrays$exposure_beta[idx, , drop = FALSE],
                         arrays$exposure_se[idx, , drop = FALSE])
      cur_offset <<- (next_exp - 1) * O
      cur_next <<- cur_offset + 1
      cur_last <<- e1 * O
      next_exp <<- e1 + 1
    }
    end <- min(cur_last, cur_next + chunk_pairs - 1)
    out <- fastmr_tidy_compact(cur, methods, exposure_labels, outcome_labels,
                               first = cur_next, last = end, pair_offset = cur_offset)
    cur_next <<- end + 1
    out
  }
}

#' Tidy rows for a range of pairs of a compact grid result
#'
#' Accessor for [fast_mr_grid()] results created with `return = "compact"`.
#' Converts only pairs `first:last` (exposure-major pair numbers) to the tidy
#' layout, identical to the corresponding rows of `as.data.frame(x)`.
#' @param x A `fastmr_compact_grid` object.
#' @param first,last First and last pair number (1-based, inclusive).
#' @return A tidy data frame with `length(methods) * (last - first + 1)` rows.
#' @export
fastmr_grid_chunk <- function(x, first = 1, last = fastmr_compact_pairs(x)) {
  if (!inherits(x, "fastmr_compact_grid")) stop("x must be a fastmr_compact_grid", call. = FALSE)
  total <- fastmr_compact_pairs(x)
  if (!is.numeric(first) || !is.numeric(last) || length(first) != 1L || length(last) != 1L ||
      is.na(first) || is.na(last) || first < 1 || last > total || first > last) {
    stop("first and last must satisfy 1 <= first <= last <= number of pairs", call. = FALSE)
  }
  fastmr_tidy_compact(x$native, x$methods, x$exposure_labels, x$outcome_labels,
                      first = if (first == 1) 1L else first, last = last)
}

#' @export
as.data.frame.fastmr_compact_grid <- function(x, ...) {
  fastmr_tidy_grid_native(x$native, x$methods, x$exposure_labels, x$outcome_labels)
}

#' @export
print.fastmr_compact_grid <- function(x, ...) {
  cat(sprintf("<fastmr_compact_grid> %d exposures x %d outcomes = %.0f pairs, methods: %s\n",
              length(x$exposure_labels), length(x$outcome_labels), fastmr_compact_pairs(x),
              paste(as.character(x$native$methods), collapse = ", ")))
  cat("Use as.data.frame() for the tidy table or fastmr_grid_chunk() for a pair range.\n")
  invisible(x)
}
