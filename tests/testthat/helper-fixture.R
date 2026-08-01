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
