args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]]) else getwd()
.libPaths(c(file.path(root, ".local", "Rlib"), .libPaths()))
suppressPackageStartupMessages(library(fastMR))

d <- read.delim(file.path(root, "inst/extdata", "il6_crp_primary_100.tsv"),
                check.names = FALSE, stringsAsFactors = FALSE)
sizes <- c(50L, 100L, 250L, 500L, 1000L)
rows <- vector("list", length(sizes))
for (i in seq_along(sizes)) {
  n <- sizes[[i]]
  x <- matrix(rep(d$beta.exposure, n), nrow = n, byrow = TRUE)
  y <- matrix(rep(d$beta.outcome, n), nrow = n, byrow = TRUE)
  sx <- matrix(rep(d$se.exposure, n), nrow = n, byrow = TRUE)
  sy <- matrix(rep(d$se.outcome, n), nrow = n, byrow = TRUE)
  invisible(fastMR:::fastmr_grid_native(x, y, sx, sy, "ivw", 0L, NULL, 5L))
  repeats <- if (n <= 100L) 100L else if (n <= 250L) 10L else if (n <= 500L) 5L else 3L
  started <- proc.time()[["elapsed"]]
  for (repeat_index in seq_len(repeats)) {
    result <- fastMR:::fastmr_grid_native(x, y, sx, sy, "ivw", 0L, NULL, 5L)
  }
  elapsed <- (proc.time()[["elapsed"]] - started) / repeats
  rows[[i]] <- data.frame(
    exposures = n, outcomes = n, pairs = n * n,
    repeats = repeats, seconds = elapsed, pairs_per_second = n * n / elapsed,
    result_layout = paste(dim(result$beta), collapse = "x"),
    stringsAsFactors = FALSE
  )
  rm(x, y, sx, sy, result)
  gc(FALSE)
}
result <- do.call(rbind, rows)
dir.create(file.path(root, "outputs"), showWarnings = FALSE, recursive = TRUE)
write.csv(result, file.path(root, "outputs", "ivw_large_grid_scaling.csv"), row.names = FALSE)
writeLines(c(
  "# Raw compact IVW grid scaling", "",
  "IL6/CRP 82-SNP fixture; nboot=0; five requested native workers.",
  "This measures the compact native grid before constructing the tidy R data frame.",
  "", "| exposures | outcomes | pairs | repeats | seconds | pairs/s | layout |",
  "|---:|---:|---:|---:|---:|---:|---|",
  vapply(seq_len(nrow(result)), function(i) sprintf(
    "| %d | %d | %d | %d | %.6f | %.0f | %s |", result$exposures[i],
    result$outcomes[i], result$pairs[i], result$repeats[i], result$seconds[i],
    result$pairs_per_second[i], result$result_layout[i]), character(1))),
  file.path(root, "outputs", "ivw_large_grid_scaling.md"))
print(result, row.names = FALSE)
