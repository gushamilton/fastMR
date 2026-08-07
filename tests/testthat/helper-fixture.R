library(fastMR)

il6_fixture <- function() {
  path <- system.file("extdata", "il6_crp_primary_100.tsv", package = "fastMR")
  if (!nzchar(path)) path <- file.path("inst", "extdata", "il6_crp_primary_100.tsv")
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

grid_fixture <- function(n_exposure = 3L, n_outcome = 4L) {
  d <- il6_fixture()
  list(
    exposure_beta = matrix(rep(d$beta.exposure, n_exposure), nrow = n_exposure, byrow = TRUE),
    outcome_beta = matrix(rep(d$beta.outcome, n_outcome), nrow = n_outcome, byrow = TRUE),
    exposure_se = matrix(rep(d$se.exposure, n_exposure), nrow = n_exposure, byrow = TRUE),
    outcome_se = matrix(rep(d$se.outcome, n_outcome), nrow = n_outcome, byrow = TRUE)
  )
}

diagnostic_fixture <- function() {
  d <- il6_fixture()
  d$id.exposure <- "E"
  d$id.outcome <- "O"
  d$exposure <- "E"
  d$outcome <- "O"
  d$mr_keep <- TRUE
  d
}
skip_if_compressor_unconfigured <- function() {
  skip_if_not_installed("CompreSSoR")
  reference <- Sys.getenv("COMPRESSOR_CANONICAL_REFERENCE", unset = "")
  skip_if(
    !nzchar(reference) || !file.exists(reference),
    "CompreSSoR native tests require COMPRESSOR_CANONICAL_REFERENCE"
  )
}
