#!/usr/bin/env Rscript

# End-to-end 5 x 5 IVW benchmark using exact GRCh38 REF/ALT keys. The same
# real FinnGen GWAS may stand in for each study so that storage access—not
# biological heterogeneity—is isolated. Every path can be overridden by an
# environment variable; large inputs and stores remain outside the repository.

suppressPackageStartupMessages({
  library(CompreSSoR)
  library(data.table)
  library(fastMR)
  library(jsonlite)
})

compressed_path <- Sys.getenv(
  "FASTMR_COMPRESSED_STORE",
  "/Volumes/crucial_x9/CompreSSoR-benchmarks/finngen-full-pcodec-canonical.cpr"
)
tsv_path <- Sys.getenv(
  "FASTMR_TSV_GZ",
  "/Volumes/crucial_x9/CompreSSoR-bench-reset/finngen-snp/core-canonical.tsv.gz"
)
vcf_path <- Sys.getenv(
  "FASTMR_VCF_GZ",
  "/Volumes/crucial_x9/CompreSSoR-bench-reset/finngen-snp/vcf/core.vcf.gz"
)
tabix <- Sys.getenv("FASTMR_TABIX", "/opt/homebrew/bin/tabix")
output_dir <- Sys.getenv("FASTMR_BENCH_OUTPUT", file.path(getwd(), "outputs"))
runs <- as.integer(Sys.getenv("FASTMR_BENCH_RUNS", "5"))
stopifnot(runs >= 5L, dir.exists(compressed_path), file.exists(tsv_path),
          file.exists(vcf_path), file.exists(paste0(vcf_path, ".tbi")),
          file.exists(tabix))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
compressed_files_on_disk <- list.files(
  compressed_path, recursive = TRUE, full.names = TRUE,
  all.files = TRUE, no.. = TRUE
)
compressed_bytes <- sum(file.info(compressed_files_on_disk)$size)

identity_columns <- c(
  "chromosome", "base_pair_location", "effect_allele", "other_allele"
)
n <- open_compressor(compressed_path)$manifest$n_rows
row_ids <- unique(as.integer(floor(seq(0, n - 1, length.out = 25L))))
identity <- read_sumstats(
  compressed_path, variants = row_ids, columns = identity_columns
)
keys <- compressor_variant_key(
  identity$chromosome, identity$base_pair_location,
  identity$other_allele, identity$effect_allele
)
key_fields <- tstrsplit(keys, ":", fixed = TRUE)
wanted <- data.table(
  chrom = key_fields[[1L]],
  pos = as.integer(key_fields[[2L]]),
  want_ref = key_fields[[3L]],
  want_alt = key_fields[[4L]],
  variant_key = keys
)

order_hits <- function(hits) {
  hits <- hits[a2 == want_ref & a1 == want_alt]
  hits[, variant_key := compressor_variant_key(chrom, pos, a2, a1)]
  hits <- hits[match(keys, variant_key)]
  stopifnot(nrow(hits) == length(keys), identical(unname(hits$variant_key), unname(keys)),
            all(is.finite(hits$beta)), all(is.finite(hits$se)), all(hits$se > 0))
  hits[, .(variant_key, beta = as.numeric(beta), se = as.numeric(se))]
}

read_tsv <- function() {
  x <- fread(
    tsv_path, select = c("chrom", "pos", "a1", "a2", "beta", "se"),
    showProgress = FALSE
  )
  x[, chrom := as.character(chrom)]
  order_hits(x[wanted, on = .(chrom, pos), nomatch = 0L, allow.cartesian = TRUE])
}

info_value <- function(info, field) {
  prefix <- paste0(field, "=")
  vapply(strsplit(info, ";", fixed = TRUE), function(parts) {
    value <- parts[startsWith(parts, prefix)]
    if (!length(value)) return(NA_real_)
    as.numeric(sub(prefix, "", value[[1L]], fixed = TRUE))
  }, numeric(1))
}

regions_path <- tempfile("fastmr-tabix-regions-", fileext = ".tsv")
on.exit(unlink(regions_path, force = TRUE), add = TRUE)
fwrite(unique(wanted[, .(chrom, pos, end = pos)]), regions_path,
       sep = "\t", col.names = FALSE)

read_tabix <- function() {
  lines <- system2(tabix, c("-R", shQuote(regions_path), shQuote(vcf_path)),
                   stdout = TRUE, stderr = TRUE)
  status <- attr(lines, "status")
  if (!is.null(status) && status != 0L) stop(paste(lines, collapse = "\n"))
  x <- fread(text = lines, sep = "\t", header = FALSE, showProgress = FALSE)
  stopifnot(ncol(x) >= 8L)
  setnames(x, 1:8, c("chrom", "pos", "id", "a2", "a1", "qual", "filter", "info"))
  x[, `:=`(
    chrom = as.character(chrom),
    beta = info_value(info, "BETA"),
    se = info_value(info, "SE")
  )]
  order_hits(x[wanted, on = .(chrom, pos), nomatch = 0L, allow.cartesian = TRUE])
}

run_ivw_grid <- function(exposures, outcomes) {
  rows <- vector("list", length(exposures) * length(outcomes))
  at <- 0L
  for (i in seq_along(exposures)) {
    for (j in seq_along(outcomes)) {
      at <- at + 1L
      rows[[at]] <- data.frame(
        SNP = keys,
        beta.exposure = exposures[[i]]$beta,
        beta.outcome = outcomes[[j]]$beta,
        se.exposure = exposures[[i]]$se,
        se.outcome = outcomes[[j]]$se,
        id.exposure = paste0("ex", i),
        id.outcome = paste0("out", j),
        stringsAsFactors = FALSE
      )
    }
  }
  fast_mr(rbindlist(rows), methods = "ivw", threads = 1L)
}

run_repeated_reader <- function(reader) {
  studies <- lapply(seq_len(10L), function(index) reader())
  run_ivw_grid(studies[1:5], studies[6:10])
}

compressed_files <- setNames(rep(compressed_path, 5L), paste0("ex", 1:5))
outcome_files <- setNames(rep(compressed_path, 5L), paste0("out", 1:5))
workloads <- list(
  `CompreSSoR Pcodec batch` = function() fast_mr_compressed(
    compressed_files, outcome_files, keys, methods = "ivw",
    threads = 1L, io_threads = 1L
  ),
  `VCF.gz + Tabix` = function() run_repeated_reader(read_tabix),
  `TSV.gz full scan` = function() run_repeated_reader(read_tsv)
)

# Warm only the two indexed paths. The five recorded gzip runs include their
# first full scan so cache treatment is visible rather than hidden.
invisible(workloads[[1L]]())
invisible(workloads[[2L]]())
records <- list()
for (label in names(workloads)) {
  for (run in seq_len(runs)) {
    gc()
    result <- NULL
    elapsed <- system.time(result <- workloads[[label]]())[["elapsed"]]
    stopifnot(nrow(result) == 25L, all(result$nsnp == 25L),
              max(abs(result$b - 1)) < 1e-12)
    records[[length(records) + 1L]] <- data.table(
      format = label,
      run = run,
      seconds = as.numeric(elapsed),
      pairs = nrow(result),
      instruments_per_pair = unique(result$nsnp),
      estimate = result$b[[1L]]
    )
    message(sprintf("%s run %d/%d: %.3f s", label, run, runs, elapsed))
  }
}
records <- rbindlist(records)
summary <- records[, .(
  runs = .N,
  median_seconds = median(seconds),
  min_seconds = min(seconds),
  max_seconds = max(seconds)
), by = format]
summary[, speedup_vs_tsv_gz :=
          summary[format == "TSV.gz full scan", median_seconds] / median_seconds]

csv_path <- file.path(output_dir, "compressed_io_benchmark.csv")
json_path <- file.path(output_dir, "compressed_io_benchmark.json")
fwrite(records, csv_path)
write_json(list(
  schema_version = "1.0.0",
  measured_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  machine = Sys.info()[c("nodename", "sysname", "release", "machine")],
  rows_in_gwas = n,
  instrument_keys = keys,
  study_reads = 10L,
  mr_pairs = 25L,
  inputs = list(
    compressed = list(file = basename(compressed_path), bytes = unname(compressed_bytes)),
    tsv_gz = list(file = basename(tsv_path), bytes = unname(file.info(tsv_path)$size)),
    vcf_gz = list(
      file = basename(vcf_path),
      bytes = unname(file.info(vcf_path)$size),
      index_bytes = unname(file.info(paste0(vcf_path, ".tbi"))$size)
    )
  ),
  summary = summary,
  runs = records
), json_path, auto_unbox = TRUE, pretty = TRUE)
print(summary)
