# Reference copy of the pre-optimisation fast_mr() (commit 58ac74f) used to
# check that the split/native-batched implementation is identical().
fast_mr_reference <- local({
f <- function(data,
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
  environment(f) <- asNamespace("fastMR")
  f
})
