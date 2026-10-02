# Verbatim copies of the pre-vectorisation implementations (commit 58ac74f).
old_diagnostic_groups <- function(data) {
  if (!is.data.frame(data)) stop("data must be a data.frame", call. = FALSE)
  # Validate the same required columns and numeric conversions as fast_mr.
  fastmr_prepare_vectors(data)
  n <- nrow(data)
  id.exp <- if ("id.exposure" %in% names(data)) as.character(data$id.exposure) else rep("", n)
  id.out <- if ("id.outcome" %in% names(data)) as.character(data$id.outcome) else rep("", n)
  id.exp[is.na(id.exp)] <- ""
  id.out[is.na(id.out)] <- ""
  pairs <- unique(data.frame(id.exposure = id.exp, id.outcome = id.out,
                             stringsAsFactors = FALSE))
  groups <- vector("list", nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    index <- which(id.exp == pairs$id.exposure[[i]] & id.out == pairs$id.outcome[[i]])
    label.exp <- if ("exposure" %in% names(data)) as.character(data$exposure[index[[1L]]]) else pairs$id.exposure[[i]]
    label.out <- if ("outcome" %in% names(data)) as.character(data$outcome[index[[1L]]]) else pairs$id.outcome[[i]]
    if (is.na(label.exp)) label.exp <- pairs$id.exposure[[i]]
    if (is.na(label.out)) label.out <- pairs$id.outcome[[i]]
    groups[[i]] <- list(
      data = data[index, , drop = FALSE],
      id.exposure = pairs$id.exposure[[i]],
      id.outcome = pairs$id.outcome[[i]],
      exposure = label.exp,
      outcome = label.out
    )
  }
  groups
}


old_steiger_filtering <- function(data) {
  if (!is.data.frame(data)) stop("data must be a data.frame", call. = FALSE)
  groups <- old_diagnostic_groups(data)
  rows <- lapply(groups, function(group) {
    x <- group$data
    if (!"units.exposure" %in% names(x)) x$units.exposure <- NA_character_
    if (!"units.outcome" %in% names(x)) x$units.outcome <- NA_character_
    if (!fastmr_steiger_unique(x$exposure) || !fastmr_steiger_unique(x$outcome) ||
        !fastmr_steiger_unique(x$units.exposure) ||
        !fastmr_steiger_unique(x$units.outcome)) {
      stop("each exposure/outcome pair must have unique labels and units",
           call. = FALSE)
    }
    x <- fastmr_steiger_add_rsq_one(x, "exposure")
    x <- fastmr_steiger_add_rsq_one(x, "outcome")
    if (!"effective_n.exposure" %in% names(x)) {
      x$effective_n.exposure <- NA_real_
    }
    if (!"effective_n.outcome" %in% names(x)) {
      x$effective_n.outcome <- NA_real_
    }
    x$steiger_dir <- x$rsq.exposure > x$rsq.outcome
    x$steiger_pval <- fastmr_steiger_rtest_p_vector(
      sqrt(x$rsq.exposure), sqrt(x$rsq.outcome),
      x$effective_n.exposure, x$effective_n.outcome)
    p_reason <- rep("ok", nrow(x))
    p_reason <- fastmr_steiger_set_reason(
      p_reason, !x$rsq_valid.exposure,
      paste0("exposure_", x$rsq_reason.exposure))
    p_reason <- fastmr_steiger_set_reason(
      p_reason, !x$rsq_valid.outcome,
      paste0("outcome_", x$rsq_reason.outcome))
    p_reason <- fastmr_steiger_set_reason(
      p_reason, is.na(x$effective_n.exposure), "missing_exposure_n")
    p_reason <- fastmr_steiger_set_reason(
      p_reason, !is.finite(x$effective_n.exposure) |
        x$effective_n.exposure <= 3, "invalid_exposure_n")
    p_reason <- fastmr_steiger_set_reason(
      p_reason, is.na(x$effective_n.outcome), "missing_outcome_n")
    p_reason <- fastmr_steiger_set_reason(
      p_reason, !is.finite(x$effective_n.outcome) |
        x$effective_n.outcome <= 3, "invalid_outcome_n")
    p_reason <- fastmr_steiger_set_reason(
      p_reason, !is.finite(x$steiger_pval), "invalid_steiger_pval")
    x$steiger_pval_valid <- p_reason == "ok"
    x$steiger_pval_reason <- p_reason
    x
  })
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}
