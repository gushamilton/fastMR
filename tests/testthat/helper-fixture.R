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
skip_if_compressor_unavailable <- function() {
  skip_if_not_installed("CompreSSoR")
  skip_if(
    utils::packageVersion("CompreSSoR") < "0.5.0",
    "CompreSSoR native tests require version 0.5.0 or newer"
  )
}

# This is deliberately a prepared, explicit-identity GRCh38 fixture.  In
# CompreSSoR 0.5 the canonical identity is stored in the native store itself;
# no dbSNP/EBI dictionary or COMPRESSOR_CANONICAL_REFERENCE path is needed.
compressor_canonical_fixture <- function(multiplier = 1) {
  n <- 80L
  data.frame(
    chromosome = rep("1", n),
    base_pair_location = seq.int(100001L, length.out = n),
    reference_allele = "A",
    alternate_allele = rep(c("C", "G", "T"), length.out = n),
    effect_allele = rep(c("C", "G", "T"), length.out = n),
    other_allele = "A",
    beta = multiplier * seq(-0.2, 0.2, length.out = n),
    standard_error = 0.02 + (seq_len(n) %% 7L) / 1000,
    effect_allele_frequency = seq(0.05, 0.95, length.out = n),
    stringsAsFactors = FALSE
  )
}
