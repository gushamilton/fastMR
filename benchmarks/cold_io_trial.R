#!/usr/bin/env Rscript

# Execute one fresh-process trial for the cache-controlled compressed-I/O
# benchmark. Fixture staging and trial randomisation are handled by
# run_cold_io_suite.py; this worker times only access, matching, and MR.

suppressPackageStartupMessages({
  library(CompreSSoR)
  library(data.table)
  library(fastMR)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: cold_io_trial.R SPEC.json RESULT.json", call. = FALSE)
}
spec <- fromJSON(args[[1L]], simplifyVector = FALSE)
result_path <- args[[2L]]

elapsed <- function() unname(proc.time()[["elapsed"]])
as_character_vector <- function(value) unname(unlist(value, use.names = FALSE))
paths <- as_character_vector(spec$paths)
keys <- as_character_vector(spec$keys)
format <- as.character(spec$format)
workload <- as.character(spec$workload)
tabix <- as.character(spec$tabix)
bcftools <- as.character(spec$bcftools)
exposures <- as.integer(spec$exposures)
outcomes <- as.integer(spec$outcomes)
io_threads <- as.integer(spec$io_threads)
stopifnot(length(io_threads) == 1L, is.finite(io_threads), io_threads >= 1L)
stopifnot(length(paths) == exposures + outcomes || workload == "full_load")

options(
  CompreSSoR.persistent_worker = TRUE,
  CompreSSoR.native_bridge = TRUE,
  CompreSSoR.coalesce_batch_reads = FALSE,
  CompreSSoR.tempdir = as.character(spec$tempdir)
)
Sys.setenv(
  COMPRESSOR_PAGE_CACHE_PAGES = "0",
  COMPRESSOR_RELOAD_EXCEPTIONS = "1"
)

key_parts <- tstrsplit(keys, ":", fixed = TRUE)
wanted <- data.table(
  chromosome = key_parts[[1L]],
  base_pair_location = as.integer(key_parts[[2L]]),
  reference_allele = key_parts[[3L]],
  alternate_allele = key_parts[[4L]],
  variant_key = keys
)

normalize_hits <- function(data) {
  setDT(data)
  data[, chromosome := sub("^chr", "", as.character(chromosome),
                           ignore.case = TRUE)]
  hits <- data[wanted, on = .(chromosome, base_pair_location),
               nomatch = 0L, allow.cartesian = TRUE]
  hits <- hits[
    reference_allele == i.reference_allele &
      alternate_allele == i.alternate_allele
  ]
  hits[, variant_key := compressor_variant_key(
    chromosome, base_pair_location, reference_allele, alternate_allele
  )]
  hits <- hits[match(keys, variant_key)]
  stopifnot(
    nrow(hits) == length(keys),
    identical(unname(hits$variant_key), unname(keys)),
    all(is.finite(hits$beta)), all(is.finite(hits$standard_error)),
    all(hits$standard_error > 0)
  )
  hits[, .(
    variant_key,
    beta = as.numeric(beta),
    standard_error = as.numeric(standard_error)
  )]
}

read_tsv_sparse <- function(path) {
  data <- fread(
    path,
    select = c("chrom", "pos", "a1", "a2", "beta", "se"),
    showProgress = FALSE
  )
  setnames(
    data,
    c("chrom", "pos", "a1", "a2", "beta", "se"),
    c("chromosome", "base_pair_location", "alternate_allele",
      "reference_allele", "beta", "standard_error")
  )
  normalize_hits(data)
}

info_value <- function(info, field) {
  prefix <- paste0(field, "=")
  vapply(strsplit(info, ";", fixed = TRUE), function(parts) {
    value <- parts[startsWith(parts, prefix)]
    if (!length(value)) return(NA_real_)
    as.numeric(sub(prefix, "", value[[1L]], fixed = TRUE))
  }, numeric(1L))
}

regions_path <- file.path(as.character(spec$tempdir), "regions.tsv")
fwrite(
  unique(wanted[, .(
    chromosome, base_pair_location, end = base_pair_location
  )]),
  regions_path, sep = "\t", col.names = FALSE
)

read_vcf_sparse <- function(path) {
  lines <- system2(
    tabix, c("-R", shQuote(regions_path), shQuote(path)),
    stdout = TRUE, stderr = TRUE
  )
  status <- attr(lines, "status")
  if (!is.null(status) && status != 0L) stop(paste(lines, collapse = "\n"))
  data <- fread(text = lines, sep = "\t", header = FALSE,
                showProgress = FALSE)
  stopifnot(ncol(data) >= 8L)
  setnames(
    data, 1:8,
    c("chromosome", "base_pair_location", "id", "reference_allele",
      "alternate_allele", "qual", "filter", "info")
  )
  data[, `:=`(
    beta = info_value(info, "BETA"),
    standard_error = info_value(info, "SE")
  )]
  normalize_hits(data)
}

read_pcodec_sparse <- function(path) {
  data <- fast_read_compressed(
    path, variants = keys, columns = c("beta", "standard_error")
  )
  data <- as.data.table(data)
  data <- data[match(keys, variant_key)]
  stopifnot(
    nrow(data) == length(keys),
    identical(unname(data$variant_key), unname(keys))
  )
  data[, .(
    variant_key,
    beta = as.numeric(beta),
    standard_error = as.numeric(standard_error)
  )]
}

run_grid <- function(studies) {
  exposure_data <- studies[seq_len(exposures)]
  outcome_data <- studies[exposures + seq_len(outcomes)]
  make_matrix <- function(data, column, labels) {
    output <- do.call(rbind, lapply(data, `[[`, column))
    rownames(output) <- labels
    colnames(output) <- keys
    output
  }
  fast_mr_grid(
    exposure_beta = make_matrix(
      exposure_data, "beta", paste0("ex", seq_len(exposures))
    ),
    outcome_beta = make_matrix(
      outcome_data, "beta", paste0("out", seq_len(outcomes))
    ),
    exposure_se = make_matrix(
      exposure_data, "standard_error", paste0("ex", seq_len(exposures))
    ),
    outcome_se = make_matrix(
      outcome_data, "standard_error", paste0("out", seq_len(outcomes))
    ),
    methods = "ivw", nboot = 0L, threads = 1L
  )
}

serialize_mr <- function(result) {
  result <- as.data.table(result)
  setorderv(result, c("id.exposure", "id.outcome", "method"))
  result[, .(
    id_exposure = id.exposure,
    id_outcome = id.outcome,
    method,
    nsnp = as.integer(nsnp),
    b = as.numeric(b),
    se = as.numeric(se),
    pval = as.numeric(pval)
  )]
}

chromosome_number <- function(value) {
  match(sub("^chr", "", as.character(value), ignore.case = TRUE),
        c(as.character(1:22), "X", "Y"))
}
base_number <- function(value) match(as.character(value), c("A", "C", "G", "T"))

touch_full <- function(data) {
  stopifnot(
    all(c(
      "chromosome", "base_pair_location", "reference_allele",
      "alternate_allele", "beta", "standard_error",
      "effect_allele_frequency"
    ) %in% names(data))
  )
  z <- data$beta / data$standard_error
  list(
    rows = nrow(data),
    identity_sum = sum(
      chromosome_number(data$chromosome) * 31 +
        as.numeric(data$base_pair_location) * 17 +
        base_number(data$reference_allele) * 5 +
        base_number(data$alternate_allele),
      na.rm = TRUE
    ),
    beta_sum = sum(data$beta, na.rm = TRUE),
    beta_sum_squares = sum(data$beta * data$beta, na.rm = TRUE),
    se_sum = sum(data$standard_error, na.rm = TRUE),
    eaf_sum = sum(data$effect_allele_frequency, na.rm = TRUE),
    z_sum = sum(z, na.rm = TRUE),
    first_key = compressor_variant_key(
      data$chromosome[[1L]], data$base_pair_location[[1L]],
      data$reference_allele[[1L]], data$alternate_allele[[1L]]
    ),
    last_key = compressor_variant_key(
      data$chromosome[[nrow(data)]], data$base_pair_location[[nrow(data)]],
      data$reference_allele[[nrow(data)]],
      data$alternate_allele[[nrow(data)]]
    )
  )
}

read_full_pcodec <- function(path) {
  data <- read_sumstats(path, columns = c(
    "chromosome", "base_pair_location", "effect_allele", "other_allele",
    "beta", "standard_error", "effect_allele_frequency"
  ))
  setnames(
    data, c("effect_allele", "other_allele"),
    c("alternate_allele", "reference_allele")
  )
  data
}

read_full_tsv <- function(path) {
  data <- fread(
    path,
    select = c("chrom", "pos", "a1", "a2", "beta", "se", "eaf"),
    showProgress = FALSE
  )
  setnames(
    data,
    c("chrom", "pos", "a1", "a2", "beta", "se", "eaf"),
    c("chromosome", "base_pair_location", "alternate_allele",
      "reference_allele", "beta", "standard_error",
      "effect_allele_frequency")
  )
  data
}

read_full_vcf <- function(path) {
  query_format <- paste0(
    "%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/BETA\\t%INFO/SE",
    "\\t%INFO/EAF\\n"
  )
  command <- paste(
    shQuote(bcftools), "query -f", shQuote(query_format), shQuote(path)
  )
  data <- fread(
    cmd = command,
    sep = "\t", header = FALSE,
    col.names = c(
      "chromosome", "base_pair_location", "reference_allele",
      "alternate_allele", "beta", "standard_error",
      "effect_allele_frequency"
    ),
    showProgress = FALSE
  )
  data
}

gc()
total_started <- elapsed()
if (workload == "full_load") {
  read_started <- elapsed()
  data <- switch(
    format,
    pcodec_full = read_full_pcodec(paths[[1L]]),
    tsv_full = read_full_tsv(paths[[1L]]),
    vcf_full = read_full_vcf(paths[[1L]]),
    stop("unsupported full-load format: ", format, call. = FALSE)
  )
  read_seconds <- elapsed() - read_started
  touch_started <- elapsed()
  checksum <- touch_full(data)
  touch_seconds <- elapsed() - touch_started
  result <- list(
    format = format,
    workload = workload,
    total_seconds = elapsed() - total_started,
    io_seconds = read_seconds,
    estimator_seconds = 0,
    touch_seconds = touch_seconds,
    object_bytes = as.numeric(object.size(data)),
    checksum = checksum
  )
} else if (format == "pcodec_direct") {
  exposure_paths <- setNames(
    paths[seq_len(exposures)], paste0("ex", seq_len(exposures))
  )
  outcome_paths <- setNames(
    paths[exposures + seq_len(outcomes)], paste0("out", seq_len(outcomes))
  )
  mr <- fast_mr_compressed(
    exposure_paths, outcome_paths, keys,
    methods = "ivw", nboot = 0L, threads = 1L, io_threads = io_threads
  )
  internal <- attr(mr, "compressed_input")$timing
  result <- list(
    format = format,
    workload = workload,
    total_seconds = elapsed() - total_started,
    io_seconds = internal$io_seconds,
    estimator_seconds = internal$estimator_seconds,
    io_threads = io_threads,
    touch_seconds = 0,
    pairs = nrow(mr),
    mr = serialize_mr(mr)
  )
} else {
  reader <- switch(
    format,
    pcodec_explicit = read_pcodec_sparse,
    tsv_gz = read_tsv_sparse,
    vcf_tabix = read_vcf_sparse,
    stop("unsupported MR format: ", format, call. = FALSE)
  )
  io_started <- elapsed()
  studies <- lapply(paths, reader)
  io_seconds <- elapsed() - io_started
  estimator_started <- elapsed()
  mr <- run_grid(studies)
  estimator_seconds <- elapsed() - estimator_started
  result <- list(
    format = format,
    workload = workload,
    total_seconds = elapsed() - total_started,
    io_seconds = io_seconds,
    estimator_seconds = estimator_seconds,
    touch_seconds = 0,
    pairs = nrow(mr),
    mr = serialize_mr(mr)
  )
}

result$measured_utc <- base::format(Sys.time(), tz = "UTC", usetz = TRUE)
result$package_versions <- list(
  CompreSSoR = as.character(packageVersion("CompreSSoR")),
  fastMR = as.character(packageVersion("fastMR")),
  data.table = as.character(packageVersion("data.table"))
)
write_json(result, result_path, auto_unbox = TRUE, pretty = TRUE, digits = NA)
