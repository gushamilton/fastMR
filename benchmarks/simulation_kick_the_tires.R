args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]]) else getwd()
.libPaths(c(file.path(root, ".local", "Rlib"),
            "/Users/fergushamilton/projects/twosamplemr-fast/.local/Rlib",
            .libPaths()))
suppressPackageStartupMessages({
  library(fastMR)
  library(TwoSampleMR)
})

max_delta <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (any(ok)) max(abs(a[ok] - b[ok])) else 0
}

complement <- c(A = "T", T = "A", C = "G", G = "C")

simulate_summary_stats <- function(n = 60L, seed = 20260812L) {
  set.seed(seed)
  snp <- paste0("rs_sim_", seq_len(n))
  allele_pairs <- rep(c("AG", "CT", "AC", "GT", "AT", "CG"), length.out = n)
  effect <- substr(allele_pairs, 1L, 1L)
  other <- substr(allele_pairs, 2L, 2L)
  effect[5L] <- "A"; other[5L] <- "T"
  effect[6L] <- "C"; other[6L] <- "G"
  effect[7L] <- "A"; other[7L] <- "C"
  beta_exposure <- rnorm(n, 0, 0.08)
  se_exposure <- runif(n, 0.01, 0.04)
  beta_outcome <- 0.35 * beta_exposure + rnorm(n, 0, 0.03)
  se_outcome <- runif(n, 0.01, 0.04)
  eaf_exposure <- runif(n, 0.1, 0.9)
  eaf_outcome <- eaf_exposure + rnorm(n, 0, 0.03)
  eaf_outcome <- pmin(0.95, pmax(0.05, eaf_outcome))

  # Deliberate orientation errors in the outcome data.
  out_effect <- effect
  out_other <- other
  out_beta <- beta_outcome
  out_eaf <- eaf_outcome
  out_effect[2L] <- other[2L]; out_other[2L] <- effect[2L]
  out_beta[2L] <- -beta_outcome[2L]; out_eaf[2L] <- 1 - eaf_outcome[2L]
  out_effect[3L] <- unname(complement[effect[3L]])
  out_other[3L] <- unname(complement[other[3L]])
  out_effect[4L] <- unname(complement[other[4L]])
  out_other[4L] <- unname(complement[effect[4L]])
  out_beta[4L] <- -beta_outcome[4L]; out_eaf[4L] <- 1 - eaf_outcome[4L]
  eaf_exposure[5L] <- 0.20; out_eaf[5L] <- 0.80
  eaf_exposure[6L] <- 0.45; out_eaf[6L] <- 0.55
  out_effect[7L] <- "A"; out_other[7L] <- "G"

  exposure <- data.frame(
    SNP = snp, beta.exposure = beta_exposure, se.exposure = se_exposure,
    effect_allele.exposure = effect, other_allele.exposure = other,
    eaf.exposure = eaf_exposure, id.exposure = "sim_exposure",
    exposure = "simulated exposure", stringsAsFactors = FALSE
  )
  outcome <- data.frame(
    SNP = snp, beta.outcome = out_beta, se.outcome = se_outcome,
    effect_allele.outcome = out_effect, other_allele.outcome = out_other,
    eaf.outcome = out_eaf, id.outcome = "sim_outcome",
    outcome = "simulated outcome", stringsAsFactors = FALSE
  )
  list(exposure = exposure, outcome = outcome,
       beta_exposure = beta_exposure, beta_outcome = beta_outcome,
       se_exposure = se_exposure, se_outcome = se_outcome, snp = snp)
}

sim <- simulate_summary_stats()
fast_h <- suppressWarnings(fast_harmonise_data(sim$exposure, sim$outcome, action = 2L))
native_h <- suppressMessages(suppressWarnings(
  TwoSampleMR::harmonise_data(sim$exposure, sim$outcome, action = 2L)))
common <- intersect(fast_h$SNP, native_h$SNP)
fast_h <- fast_h[match(common, fast_h$SNP), , drop = FALSE]
native_h <- native_h[match(common, native_h$SNP), , drop = FALSE]
harmony <- data.frame(
  common_snps = length(common), fast_rows = nrow(fast_h), native_rows = nrow(native_h),
  keep_mismatches = sum(fast_h$mr_keep != native_h$mr_keep, na.rm = TRUE),
  remove_mismatches = sum(fast_h$remove != native_h$remove, na.rm = TRUE),
  palindrome_mismatches = sum(fast_h$palindromic != native_h$palindromic, na.rm = TRUE),
  ambiguous_mismatches = sum(fast_h$ambiguous != native_h$ambiguous, na.rm = TRUE),
  max_beta_delta = max_delta(fast_h$beta.outcome, native_h$beta.outcome),
  stringsAsFactors = FALSE
)
stopifnot(harmony$fast_rows == harmony$native_rows,
          harmony$keep_mismatches == 0L,
          harmony$remove_mismatches == 0L,
          harmony$palindrome_mismatches == 0L,
          harmony$ambiguous_mismatches == 0L,
          harmony$max_beta_delta <= 1e-12)

harmony_actions <- do.call(rbind, lapply(1:3, function(action) {
  fh <- suppressWarnings(fast_harmonise_data(sim$exposure, sim$outcome, action = action))
  nh <- suppressMessages(suppressWarnings(
    TwoSampleMR::harmonise_data(sim$exposure, sim$outcome, action = action)))
  ids <- intersect(fh$SNP, nh$SNP)
  fh <- fh[match(ids, fh$SNP), , drop = FALSE]
  nh <- nh[match(ids, nh$SNP), , drop = FALSE]
  data.frame(
    action = action, fast_rows = nrow(fh), native_rows = nrow(nh),
    keep_mismatches = sum(fh$mr_keep != nh$mr_keep, na.rm = TRUE),
    remove_mismatches = sum(fh$remove != nh$remove, na.rm = TRUE),
    ambiguous_mismatches = sum(fh$ambiguous != nh$ambiguous, na.rm = TRUE),
    max_beta_delta = max_delta(fh$beta.outcome, nh$beta.outcome),
    stringsAsFactors = FALSE
  )
}))
stopifnot(all(harmony_actions$fast_rows == harmony_actions$native_rows),
          all(harmony_actions$keep_mismatches == 0L),
          all(harmony_actions$remove_mismatches == 0L),
          all(harmony_actions$ambiguous_mismatches == 0L),
          max(harmony_actions$max_beta_delta) <= 1e-12)

# Unequal exposure/outcome SNP sets: only the intersection should survive.
set.seed(20260813L)
exposure_index <- sort(sample(seq_len(nrow(sim$exposure)), 47L))
outcome_index <- sort(sample(seq_len(nrow(sim$outcome)), 34L))
unequal_h <- suppressWarnings(fast_harmonise_data(
  sim$exposure[exposure_index, , drop = FALSE],
  sim$outcome[outcome_index, , drop = FALSE], action = 2L))
expected_overlap <- length(intersect(sim$exposure$SNP[exposure_index],
                                     sim$outcome$SNP[outcome_index]))
unequal <- data.frame(
  exposure_snps = length(exposure_index), outcome_snps = length(outcome_index),
  expected_overlap = expected_overlap, harmonised_rows = nrow(unequal_h),
  stringsAsFactors = FALSE
)
stopifnot(unequal$harmonised_rows == unequal$expected_overlap)

# Compare several MR estimators on the uneven, harmonised simulation.
fast_unequal <- suppressWarnings(fast_mr(
  unequal_h, methods = c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode"),
  nboot = 0L, seed = 20260814L))
native_unequal <- suppressMessages(suppressWarnings(TwoSampleMR::mr(
  unequal_h,
  method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median",
                  "mr_simple_mode", "mr_weighted_mode"),
  parameters = modifyList(TwoSampleMR::default_parameters(), list(nboot = 0L)))))
native_map <- c("Inverse variance weighted" = "ivw", "MR Egger" = "egger",
                "Weighted median" = "weighted_median", "Simple mode" = "simple_mode",
                "Weighted mode" = "weighted_mode")
native_unequal$method_code <- unname(native_map[native_unequal$method])
method_order <- c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode")
method_parity <- do.call(rbind, lapply(method_order, function(method) {
  a <- fast_unequal[fast_unequal$method_code == method, , drop = FALSE]
  b <- native_unequal[native_unequal$method_code == method, , drop = FALSE]
  data.frame(method = method, fast_rows = nrow(a), native_rows = nrow(b),
             beta_delta = max_delta(a$b, b$b), se_delta = max_delta(a$se, b$se),
             pval_delta = max_delta(a$pval, b$pval), stringsAsFactors = FALSE)
}))
stopifnot(all(method_parity$fast_rows == 1L), all(method_parity$native_rows == 1L),
          max(method_parity$beta_delta) <= 1e-10)

# A rectangular grid on the same simulated SNP panel, with native checks on
# representative exposure/outcome pairs.
grid_snps <- intersect(sim$exposure$SNP, sim$outcome$SNP)
grid_n_exp <- 7L; grid_n_out <- 11L
exp_scale <- 1 + (seq_len(grid_n_exp) - 4L) * 0.01
out_scale <- 1 + (seq_len(grid_n_out) - 6L) * 0.008
grid <- list(
  exposure_beta = matrix(rep(sim$beta_exposure[match(grid_snps, sim$snp)], grid_n_exp),
                         nrow = grid_n_exp, byrow = TRUE) * exp_scale,
  outcome_beta = matrix(rep(sim$beta_outcome[match(grid_snps, sim$snp)], grid_n_out),
                        nrow = grid_n_out, byrow = TRUE) * out_scale,
  exposure_se = matrix(rep(sim$se_exposure[match(grid_snps, sim$snp)], grid_n_exp),
                       nrow = grid_n_exp, byrow = TRUE),
  outcome_se = matrix(rep(sim$se_outcome[match(grid_snps, sim$snp)], grid_n_out),
                      nrow = grid_n_out, byrow = TRUE)
)
grid_methods <- c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode")
fast_grid <- fast_mr_grid(grid$exposure_beta, grid$outcome_beta,
                          grid$exposure_se, grid$outcome_se,
                          methods = grid_methods, nboot = 0L, threads = 5L)
representatives <- data.frame(exposure = c(1L, 4L, 7L), outcome = c(1L, 8L, 11L))
grid_rows <- list()
for (i in seq_len(nrow(representatives))) {
  e <- representatives$exposure[i]; o <- representatives$outcome[i]
  dat <- data.frame(
    SNP = grid_snps, beta.exposure = grid$exposure_beta[e, ],
    beta.outcome = grid$outcome_beta[o, ], se.exposure = grid$exposure_se[e, ],
    se.outcome = grid$outcome_se[o, ], id.exposure = paste0("E", e),
    id.outcome = paste0("O", o), exposure = paste0("E", e), outcome = paste0("O", o),
    mr_keep = TRUE, stringsAsFactors = FALSE
  )
  native <- suppressMessages(suppressWarnings(TwoSampleMR::mr(
    dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median",
                         "mr_simple_mode", "mr_weighted_mode"),
    parameters = modifyList(TwoSampleMR::default_parameters(), list(nboot = 0L)))))
  native$method_code <- unname(native_map[native$method])
  for (method in grid_methods) {
    a <- fast_grid[fast_grid$exposure_index == e & fast_grid$outcome_index == o &
                     fast_grid$method_code == method, , drop = FALSE]
    b <- native[native$method_code == method, , drop = FALSE]
    grid_rows[[length(grid_rows) + 1L]] <- data.frame(
      exposure = e, outcome = o, method = method,
      beta_delta = max_delta(a$b, b$b), se_delta = max_delta(a$se, b$se),
      pval_delta = max_delta(a$pval, b$pval), stringsAsFactors = FALSE)
  }
}
grid_parity <- do.call(rbind, grid_rows)
stopifnot(max(grid_parity$beta_delta) <= 1e-10)

# Empty intersection and a deliberately short instrument set should return
# safely rather than error.
no_overlap <- fast_harmonise_data(
  sim$exposure[1:4, , drop = FALSE],
  transform(sim$outcome[5:8, , drop = FALSE], SNP = paste0("no_", SNP)), action = 2L)
short <- suppressWarnings(fast_mr(no_overlap, methods = "ivw", nboot = 0L))
edge <- data.frame(no_overlap_rows = nrow(no_overlap), short_result_rows = nrow(short),
                   short_beta_all_na = all(is.na(short$b)), stringsAsFactors = FALSE)
stopifnot(edge$no_overlap_rows == 0L, edge$short_result_rows == 0L,
          edge$short_beta_all_na)

dir.create(file.path(root, "outputs"), showWarnings = FALSE, recursive = TRUE)
write.csv(harmony, file.path(root, "outputs", "simulation_harmonisation_summary.csv"), row.names = FALSE)
write.csv(harmony_actions, file.path(root, "outputs", "simulation_harmonisation_actions.csv"), row.names = FALSE)
write.csv(unequal, file.path(root, "outputs", "simulation_unequal_snps.csv"), row.names = FALSE)
write.csv(method_parity, file.path(root, "outputs", "simulation_unequal_mr_parity.csv"), row.names = FALSE)
write.csv(grid_parity, file.path(root, "outputs", "simulation_grid_parity.csv"), row.names = FALSE)
write.csv(edge, file.path(root, "outputs", "simulation_edge_cases.csv"), row.names = FALSE)
writeLines(c(
  "# Simulation and harmonisation tyre-kick", "",
  sprintf("Simulation: %d SNPs, causal effect 0.35; deliberately swapped, complemented, reverse-complemented, palindromic, and incompatible alleles.", nrow(sim$exposure)),
  sprintf("Harmonisation: %d common SNPs; fast/native keep mismatches %d; max beta delta %.3e.", harmony$common_snps, harmony$keep_mismatches, harmony$max_beta_delta),
  sprintf("Harmonisation actions 1/2/3: maximum keep mismatch %d; maximum beta delta %.3e.", max(harmony_actions$keep_mismatches), max(harmony_actions$max_beta_delta)),
  sprintf("Unequal SNP sets: %d exposure, %d outcome, %d expected/retained overlap.", unequal$exposure_snps, unequal$outcome_snps, unequal$harmonised_rows),
  sprintf("Unequal-set MR: maximum beta delta %.3e across five methods.", max(method_parity$beta_delta)),
  sprintf("Grid: %dx%d pairs on %d SNPs; representative native parity maximum beta delta %.3e.", grid_n_exp, grid_n_out, length(grid_snps), max(grid_parity$beta_delta)),
  sprintf("Empty-intersection edge case: %d rows; safe short result with all-NA beta: %s.", edge$no_overlap_rows, edge$short_beta_all_na)
), file.path(root, "outputs", "simulation_kick_the_tires.md"))
print(harmony, row.names = FALSE)
print(harmony_actions, row.names = FALSE)
print(unequal, row.names = FALSE)
print(method_parity, row.names = FALSE)
print(data.frame(grid_pairs = grid_n_exp * grid_n_out, grid_snps = length(grid_snps),
                 representative_rows = nrow(grid_parity),
                 max_beta_delta = max(grid_parity$beta_delta)), row.names = FALSE)
print(edge, row.names = FALSE)
