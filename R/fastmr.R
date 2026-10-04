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
#' @param threads Maximum native worker count. Exposure/outcome pairs run in
#'   parallel, including bootstrap methods: their normal draws are made on the
#'   main thread in the serial order, so results (seeded or unseeded) and the
#'   RNG state afterwards are identical for every thread count.
#' @param output Optional path for a Zstandard-compressed Parquet copy of the
#'   result. The path must not already exist; use [fast_write_parquet()] when
#'   an overwrite or another compression codec is required.
#' @param ... Optional `phi` bandwidth multiplier for mode methods and `penk`
#'   penalty multiplier for penalised weighted median (default 20).
#' @details With a single instrument (`nsnp = 1`) the IVW estimators (`"ivw"`,
#'   `"ivw_fe"`, `"ivw_mre"`) return the Wald ratio, as TwoSampleMR's `mr()`
#'   reports it through `mr_wald_ratio()`: `b = beta.outcome / beta.exposure`
#'   and `se = se.outcome / |beta.exposure|` (first order); `Q` and `sigma` are
#'   `NA`. Every other method needs at least two (`"uwr"`), three or six
#'   (`"sign"`) instruments and returns `NA`, as TwoSampleMR does. The same
#'   holds for [fast_mr_grid()], [fast_mr_sparse_ivw()],
#'   [fast_mr_masked_ivw()] and [fast_mr_compressed()]. An exact fit with two
#'   or more instruments (residual standard error 0) is a different case and
#'   keeps the fixed-effect standard error. Missing results are `NA`.
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
  # Group rows once (first-appearance order of each exposure/outcome pair).
  # Integer codes avoid any separator collision between ids.
  gid <- fastmr_group_ids(id.exp, id.out)
  group_count <- attr(gid, "n")
  if (!group_count) return(fastmr_write_result(fastmr_tidy_native(list(), methods), output))
  # Joins and multi-study exports often repeat the same SNP row. Count each
  # SNP once per MR pair; retain the first kept row deterministically. Repeated
  # p-values are metadata and do not affect this rule.
  kept <- which(keep)
  snp_code <- match(snp[kept], unique(snp[kept]))
  kept <- kept[!duplicated(gid[kept] + group_count * as.numeric(snp_code))]
  kept <- kept[order(gid[kept], method = "radix")]
  counts <- tabulate(gid[kept], nbins = group_count)
  offsets <- c(0L, cumsum(counts))
  # Groups are numbered by first appearance, so first rows are in group order.
  first <- which(!duplicated(gid))
  args <- list(
    offsets = offsets,
    exposure_beta = prepared[["beta.exposure"]][kept],
    outcome_beta = prepared[["beta.outcome"]][kept],
    exposure_se = prepared[["se.exposure"]][kept],
    outcome_se = prepared[["se.outcome"]][kept],
    methods = methods, nboot = controls[["nboot"]],
    threads = controls[["threads"]], phi = phi, penk = penk
  )
  native <- if (fastmr_methods_use_rng(methods, controls[["nboot"]])) {
    fastmr_check_seed_streams(controls[["seed"]], group_count)
    fastmr_run_bootstrap_groups(args, controls[["seed"]])
  } else {
    fastmr_native_call(fastmr_run_groups_native, args, NULL)
  }
  fastmr_write_result(
    fastmr_tidy_groups_native(native, length(methods), id.exp[first], id.out[first]),
    output)
}

# Bootstrap groups, threaded and bit-identical to one fastmr_run_native() call
# per group. All RNG draws happen natively on the main thread in the order those
# calls consumed the RNG: one continuous stream from the caller's .Random.seed
# when `seed` is NULL, otherwise a fresh set.seed(seed + i - 1) stream per group
# i (the caller's RNG state is then restored). With one worker, draws stream
# straight into each group's bootstrap layout as before; with several, groups
# run in batches whose normals (at most getOption("fastMR.bootstrap_batch_draws"),
# default 2^23 = 64 MB, or a single group streamed serially) are drawn into a
# reused native buffer. Batching does not change any result.
# Group i reseeds with set.seed(seed + i - 1), so the last group's seed must
# still be a valid R integer. Check before any group runs rather than failing
# part-way through the native loop.
fastmr_check_seed_streams <- function(seed, n_groups) {
  if (is.null(seed) || n_groups < 1) return(invisible(NULL))
  last <- seed + n_groups - 1
  if (last > .Machine$integer.max) {
    stop("seed + (number of exposure/outcome pairs) - 1 = ", format(last, scientific = FALSE),
         " exceeds the largest R integer (", .Machine$integer.max, "); each pair i is seeded with ",
         "seed + i - 1, so use seed <= ", format(.Machine$integer.max - n_groups + 1, scientific = FALSE),
         call. = FALSE)
  }
  invisible(NULL)
}

fastmr_run_bootstrap_groups <- function(args, seed) {
  args$batch_draws <- as.numeric(getOption("fastMR.bootstrap_batch_draws", 2^23))
  if (!is.null(seed)) {
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
    args$reseed <- function(i) set.seed(seed + i - 1L)
  }
  # The native call's RNGScope creates .Random.seed when absent, as each
  # per-group call did.
  do.call(fastmr_run_groups_boot_native, args)
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
#' @param return What to return. `"tidy"` (default) is the full tidy data frame.
#'   `"compact"` returns a light `fastmr_compact_grid` object (see
#'   [fastmr_grid_chunk()]) that converts to the identical tidy data frame via
#'   `as.data.frame()`. `"none"` requires `output` and returns the path
#'   invisibly. Streamed IVW-only output (`return = "none"`) is computed in
#'   exposure blocks and equals the tidy result to about 1e-15 relative (a
#'   blocked BLAS `dgemm`, e.g. OpenBLAS, can differ in the last bits); it is
#'   identical when the grid is not blocked.
#' @param chunk_pairs Number of grid pairs converted to tidy form per Parquet
#'   row group when a non-`"tidy"` `return` is combined with `output`
#'   (default 1e6). Only the tidy conversion is chunked, so results are
#'   identical for every chunk size. For IVW-only grids the native kernel is
#'   also run in exposure blocks of about `chunk_pairs` pairs, which is
#'   deterministic and pair-independent (blocked results agree with the
#'   unblocked ones to about 1e-15 relative, see `return`); other methods
#'   (including seeded bootstraps) always use one native call.
#' @param ... Optional `phi` bandwidth multiplier for mode methods and `penk`
#'   penalty multiplier for penalised weighted median (default 20).
#' @return With `return = "tidy"`, a tidy data frame with one row per method
#'   and grid pair; see `return` for the other modes. With `return != "tidy"`
#'   and `output`, the Parquet file is written in row groups of `chunk_pairs`
#'   pairs (Zstandard compressed, schema identical to the tidy output).
#' @export
fast_mr_grid <- function(exposure_beta, outcome_beta, exposure_se, outcome_se,
                         methods = c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode"),
                         nboot = 1000, seed = NULL, threads = 1, output = NULL,
                         return = c("tidy", "compact", "none"), chunk_pairs = 1e6, ...) {
  return <- match.arg(return)
  if (!is.numeric(chunk_pairs) || length(chunk_pairs) != 1L || !is.finite(chunk_pairs) || chunk_pairs < 1) {
    stop("chunk_pairs must be one number >= 1", call. = FALSE)
  }
  chunk_pairs <- floor(chunk_pairs)
  if (return == "none" && is.null(output)) stop("return = \"none\" requires output", call. = FALSE)
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
  exp.labels <- rownames(arrays$exposure_beta)
  out.labels <- rownames(arrays$outcome_beta)
  if (is.null(exp.labels)) exp.labels <- as.character(seq_len(nrow(arrays$exposure_beta)))
  if (is.null(out.labels)) out.labels <- as.character(seq_len(nrow(arrays$outcome_beta)))
  run_native <- function(eb, es) {
    fastmr_native_call(
      fastmr_grid_native,
      list(
        exposure_beta = eb,
        outcome_beta = arrays[["outcome_beta"]],
        exposure_se = es,
        outcome_se = arrays[["outcome_se"]],
        methods = methods, nboot = controls[["nboot"]], seed = NULL,
        threads = controls[["threads"]], phi = phi, penk = penk
      ),
      controls[["seed"]]
    )
  }
  if (return == "none" && identical(methods, "ivw")) {
    # RNG-free: run the kernel in exposure blocks and stream each to Parquet.
    path <- fastmr_stream_parquet(
      fastmr_ivw_blocks(arrays, run_native, exp.labels, out.labels, methods, chunk_pairs),
      output)
    return(invisible(path))
  }
  native <- run_native(arrays[["exposure_beta"]], arrays[["exposure_se"]])
  if (return == "tidy") {
    return(fastmr_write_result(fastmr_tidy_grid_native(native, methods, exp.labels, out.labels), output))
  }
  res <- fastmr_new_compact_grid(native, methods, exp.labels, out.labels)
  if (!is.null(output)) {
    path <- fast_write_parquet(res, output, chunk_pairs = chunk_pairs)
    if (return == "none") return(invisible(path))
  }
  res
}
