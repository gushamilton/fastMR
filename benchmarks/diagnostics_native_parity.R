args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]]) else getwd()
.libPaths(c(file.path(root, ".local", "Rlib"),
            "/Users/fergushamilton/projects/twosamplemr-fast/.local/Rlib",
            .libPaths()))

suppressPackageStartupMessages({
  library(fastMR)
  library(TwoSampleMR)
})

d <- read.delim(file.path(root, "inst", "extdata", "il6_crp_primary_100.tsv"),
                check.names = FALSE, stringsAsFactors = FALSE)
d$id.exposure <- "E"
d$id.outcome <- "O"
d$exposure <- "E"
d$outcome <- "O"
d$mr_keep <- TRUE
steiger_d <- d
steiger_d$r.exposure <- seq_len(nrow(steiger_d)) / 1000
steiger_d$r.outcome <- seq_len(nrow(steiger_d)) / 1200
steiger_d$samplesize.exposure <- 10000
steiger_d$samplesize.outcome <- 12000
filtering_d <- d
filtering_d$units.exposure <- "SD"
filtering_d$units.outcome <- "SD"
filtering_d$samplesize.exposure <- 10000
filtering_d$samplesize.outcome <- 12000

fast_repeats <- 20L
native_repeats <- 5L
threads <- 5L

timed <- function(fun, repeats) {
  invisible(fun())
  elapsed <- numeric(repeats)
  result <- NULL
  for (i in seq_len(repeats)) {
    started <- proc.time()[["elapsed"]]
    result <- fun()
    elapsed[i] <- proc.time()[["elapsed"]] - started
  }
  list(result = result, seconds = median(elapsed))
}

compare <- function(component, fast_fun, native_fun, key = NULL,
                    primary = "b", secondary = "se", p_value = "pval") {
  fast <- timed(fast_fun, fast_repeats)
  native <- timed(native_fun, native_repeats)
  a <- fast$result
  b <- native$result
  if (!is.null(key)) {
    a <- a[match(key, a$SNP), , drop = FALSE]
    b <- b[match(key, b$SNP), , drop = FALSE]
  }
  finite_max <- function(x) if (any(is.finite(x))) max(abs(x), na.rm = TRUE) else NA_real_
  data.frame(
    component = component,
    fastMR_seconds = fast$seconds,
    TwoSampleMR_seconds = native$seconds,
    speedup = native$seconds / fast$seconds,
    fast_rows = nrow(a),
    native_rows = nrow(b),
    row_count_match = nrow(a) == nrow(b),
    key_match = if (is.null(key)) NA else identical(as.character(a$SNP), as.character(b$SNP)),
    max_abs_primary_delta = finite_max(a[[primary]] - b[[primary]]),
    max_abs_secondary_delta = finite_max(a[[secondary]] - b[[secondary]]),
    max_abs_pval_delta = finite_max(a[[p_value]] - b[[p_value]]),
    stringsAsFactors = FALSE
  )
}

native_heterogeneity <- function() suppressMessages(suppressWarnings(
  TwoSampleMR::mr_heterogeneity(d, method_list = c("mr_ivw", "mr_egger_regression"))))
native_pleiotropy <- function() suppressMessages(suppressWarnings(
  TwoSampleMR::mr_pleiotropy_test(d)))
native_singlesnp <- function() suppressMessages(suppressWarnings(
  TwoSampleMR::mr_singlesnp(d, single_method = "mr_wald_ratio",
                            all_method = c("mr_ivw", "mr_egger_regression"))))
native_leaveoneout <- function() suppressMessages(suppressWarnings(
  TwoSampleMR::mr_leaveoneout(d, method = TwoSampleMR::mr_ivw)))
native_directionality <- function() suppressMessages(suppressWarnings(
  TwoSampleMR::directionality_test(steiger_d)))
native_steiger_filtering <- function() suppressMessages(suppressWarnings(
  TwoSampleMR::steiger_filtering(filtering_d)))

fast_heterogeneity <- function() fast_mr_heterogeneity(d, threads = threads)
fast_pleiotropy <- function() fast_mr_pleiotropy_test(d, threads = threads)
fast_singlesnp <- function() fast_mr_singlesnp(d, threads = threads)
fast_leaveoneout <- function() fast_mr_leaveoneout(d, method = "ivw", threads = threads)
fast_directionality <- function() fast_mr_directionality_test(steiger_d)
fast_steiger_filtering <- function() fast_mr_steiger_filtering(filtering_d)

singlesnp_native_once <- native_singlesnp()
singlesnp_key <- singlesnp_native_once$SNP
leaveoneout_native_once <- native_leaveoneout()
leaveoneout_key <- leaveoneout_native_once$SNP
steiger_filtering_key <- native_steiger_filtering()$SNP

rows <- list(
  compare("heterogeneity", fast_heterogeneity, native_heterogeneity,
          primary = "Q", secondary = "Q_df", p_value = "Q_pval"),
  compare("egger_pleiotropy", fast_pleiotropy, native_pleiotropy,
          primary = "egger_intercept"),
  compare("single_snp", fast_singlesnp, native_singlesnp, singlesnp_key,
          primary = "b", secondary = "se", p_value = "p"),
  compare("leave_one_out", fast_leaveoneout, native_leaveoneout, leaveoneout_key,
          primary = "b", secondary = "se", p_value = "p"),
  compare("directionality", fast_directionality, native_directionality,
          primary = "snp_r2.exposure", secondary = "snp_r2.outcome",
          p_value = "steiger_pval"),
  compare("steiger_filtering", fast_steiger_filtering, native_steiger_filtering,
          steiger_filtering_key, primary = "rsq.exposure",
          secondary = "rsq.outcome", p_value = "steiger_pval")
)
result <- do.call(rbind, rows)
dir.create(file.path(root, "outputs"), showWarnings = FALSE, recursive = TRUE)
write.csv(result, file.path(root, "outputs", "diagnostics_native_parity.csv"), row.names = FALSE)

lines <- c(
  "# fastMR diagnostics versus native TwoSampleMR", "",
  sprintf("TwoSampleMR 0.7.9; IL6 fixture: %d SNPs; fastMR threads=%d; timings are medians of %d fast and %d native calls.",
          nrow(d), threads, fast_repeats, native_repeats), "",
  "| component | fastMR s | TwoSampleMR s | speedup | rows | row/key match | max primary delta | max secondary delta | max p delta |",
  "|---|---:|---:|---:|---:|---|---:|---:|---:|"
)
for (i in seq_len(nrow(result))) {
  x <- result[i, ]
  match_text <- if (is.na(x$key_match)) as.character(x$row_count_match) else
    paste(x$row_count_match, x$key_match, sep = "/")
  lines <- c(lines, sprintf(
    "| %s | %.6f | %.6f | %.2fx | %d/%d | %s | %.3e | %.3e | %.3e |",
    x$component, x$fastMR_seconds, x$TwoSampleMR_seconds, x$speedup,
    x$fast_rows, x$native_rows, match_text, x$max_abs_primary_delta,
    x$max_abs_secondary_delta, x$max_abs_pval_delta))
}
writeLines(lines, file.path(root, "outputs", "diagnostics_native_parity.md"))
print(result, row.names = FALSE, digits = 6)
