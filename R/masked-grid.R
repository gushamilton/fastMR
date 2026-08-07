# Fast IVW grid with variable instrument membership.
#
# The dense beta matrices are zero-filled outside each exposure's instrument
# set, while the logical masks preserve the pair-specific SNP intersection for
# nsnp, Q, and multiplicative random-effects SE.
fastmr_masked_ivw_grid <- function(exposure_beta, outcome_beta, outcome_se,
                                   exposure_present, outcome_present,
                                   threads = 1L) {
  exposure_beta <- fastmr_matrix_numeric(exposure_beta, "exposure_beta")
  outcome_beta <- fastmr_matrix_numeric(outcome_beta, "outcome_beta")
  outcome_se <- fastmr_matrix_numeric(outcome_se, "outcome_se")
  exposure_present <- as.matrix(exposure_present)
  outcome_present <- as.matrix(outcome_present)
  if (!is.logical(exposure_present) || !is.logical(outcome_present)) {
    stop("presence masks must be logical matrices", call. = FALSE)
  }
  if (anyNA(exposure_present) || anyNA(outcome_present)) {
    stop("presence masks must not contain NA", call. = FALSE)
  }
  if (!identical(dim(exposure_beta), dim(exposure_present)) ||
      !identical(dim(outcome_beta), dim(outcome_present)) ||
      !identical(dim(outcome_beta), dim(outcome_se))) {
    stop("masked IVW inputs must have matching dimensions", call. = FALSE)
  }
  fastmr_native_call(
    fastmr_masked_ivw_native,
    list(exposure_beta = exposure_beta, outcome_beta = outcome_beta,
         outcome_se = outcome_se, exposure_present = exposure_present,
         outcome_present = outcome_present, threads = as.integer(threads)),
    NULL
  )
}
