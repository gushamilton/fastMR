#!/usr/bin/env Rscript

# Bounded ordinary two-sample MR benchmark for two moderate TSV/CSV gzip files.
#
# The script deliberately separates:
#   1. gzip parsing and column normalisation;
#   2. an explicit SNP-only inner join;
#   3. fastMR's allele harmonisation (which performs its own join);
#   4. one IVW calculation on the harmonised pair.
#
# It does not read or write any production fastMR result. The output directory
# is supplied explicitly by the Slurm wrapper and defaults to the disjoint
# results/two_sample_mr_gzip_benchmark directory.

suppressPackageStartupMessages({
  library(data.table)
  library(fastMR)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = "") {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}
arg_integer <- function(name, default) {
  value <- suppressWarnings(as.integer(arg_value(name, as.character(default))))
  if (is.na(value) || value < 1L) stop("--", name, " must be a positive integer", call. = FALSE)
  value
}

output_dir <- normalizePath(
  arg_value("output-dir", Sys.getenv("OUTPUT_DIR", "results/two_sample_mr_gzip_benchmark")),
  mustWork = FALSE
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

exposure_path <- arg_value("exposure", Sys.getenv("EXPOSURE_GZ", ""))
outcome_path <- arg_value("outcome", Sys.getenv("OUTCOME_GZ", ""))
data_root <- arg_value("data-root", Sys.getenv("DATA_ROOT", "/user/work/fh6520"))
min_bytes <- as.numeric(arg_value("min-bytes", Sys.getenv("MIN_BYTES", "52428800")))
max_bytes <- as.numeric(arg_value("max-bytes", Sys.getenv("MAX_BYTES", "786432000")))
if (!is.finite(min_bytes) || !is.finite(max_bytes) || min_bytes <= 0 || max_bytes <= min_bytes) {
  stop("min/max byte bounds are invalid", call. = FALSE)
}

stop_if_bad_file <- function(path, label) {
  if (!nzchar(path) || !file.exists(path)) stop(label, " file does not exist: ", path, call. = FALSE)
  if (!grepl("\\.(tsv|txt|csv)\\.gz$", tolower(path))) {
    stop(label, " must be a TSV/TXT/CSV gzip file: ", path, call. = FALSE)
  }
  bytes <- file.info(path)$size
  if (!is.finite(bytes)) stop("could not stat ", label, " file: ", path, call. = FALSE)
  invisible(bytes)
}

# Header aliases cover the common GWAS exports used by TwoSampleMR and the
# canonical FinnGen-style TSV.gz exports used in the fastMR benchmarks.
normalise_name <- function(x) tolower(gsub("[^a-z0-9]", "", x))
pick_column <- function(names_in, candidates, label, optional = FALSE) {
  normal <- normalise_name(names_in)
  wanted <- normalise_name(candidates)
  hit <- match(wanted, normal)
  hit <- hit[!is.na(hit)]
  if (!length(hit)) {
    if (optional) return(NA_character_)
    stop("could not find ", label, " column; available columns: ",
         paste(names_in, collapse = ", "), call. = FALSE)
  }
  names_in[[hit[[1L]]]]
}

column_map <- function(names_in, role) {
  suffix <- if (role == "exposure") ".exposure" else ".outcome"
  list(
    SNP = pick_column(names_in,
                      c("SNP", "rsid", "rs_id", "variant_id", "variant", "markername", "marker"),
                      "SNP"),
    beta = pick_column(names_in,
                       c(paste0("beta", suffix), "beta", "b", "effect", "effect_size",
                         "beta_hat", "es", "BETA"), "beta"),
    se = pick_column(names_in,
                     c(paste0("se", suffix), "se", "stderr", "standard_error",
                       "standarderror", "sebeta", "SE"), "standard error"),
    effect_allele = pick_column(names_in,
                                c(paste0("effect_allele", suffix), "effect_allele",
                                  "ea", "a1", "alt", "ALT"), "effect allele"),
    other_allele = pick_column(names_in,
                               c(paste0("other_allele", suffix), "other_allele",
                                 "oa", "a2", "ref", "REF"), "other allele"),
    eaf = pick_column(names_in,
                      c(paste0("eaf", suffix), "eaf", "eaf01", "effect_allele_frequency",
                        "effectallelefrequency", "af", "maf"), "effect allele frequency", optional = TRUE)
  )
}

discover_files <- function(root, min_bytes, max_bytes) {
  if (!dir.exists(root)) stop("DATA_ROOT does not exist: ", root, call. = FALSE)
  paths <- list.files(root, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
  paths <- paths[grepl("\\.(tsv|txt|csv)\\.gz$", tolower(paths))]
  paths <- paths[!grepl("(^|/)(tar|fat|huge)([^/]*)$|\\.tar\\.gz$", tolower(paths))]
  paths <- paths[!grepl("(^|/)results(/|$)|(^|/)logs(/|$)", paths)]
  if (!length(paths)) stop("no candidate TSV/TXT/CSV.gz files found below ", root, call. = FALSE)
  sizes <- file.info(paths)$size
  paths <- paths[is.finite(sizes) & sizes >= min_bytes & sizes <= max_bytes]
  if (length(paths) < 2L) {
    stop("fewer than two moderate gzip files below byte bounds under ", root,
         " (", format(min_bytes, scientific = FALSE), "-", format(max_bytes, scientific = FALSE), ")",
         call. = FALSE)
  }
  # Check headers before selecting. This is intentionally nrows=0: it avoids
  # decompressing the bodies during discovery and rejects unrelated gzip data.
  eligible <- vapply(paths, function(path) {
    header <- tryCatch(fread(path, nrows = 0L, showProgress = FALSE), error = function(e) NULL)
    if (is.null(header)) return(FALSE)
    ok <- tryCatch({
      column_map(names(header), "exposure")
      TRUE
    }, error = function(e) FALSE)
    isTRUE(ok)
  }, logical(1))
  paths <- paths[eligible]
  if (length(paths) < 2L) stop("fewer than two eligible summary-statistics gzip files after header checks", call. = FALSE)
  sizes <- file.info(paths)$size
  # Choose two different files near the middle of the bounded size range so
  # the comparison is representative rather than dominated by a tiny fixture.
  order_index <- order(abs(log(pmax(sizes, 1) / sqrt(min_bytes * max_bytes))), paths)
  paths[order_index[seq_len(2L)]]
}

if (!nzchar(exposure_path) || !nzchar(outcome_path)) {
  selected <- discover_files(data_root, min_bytes, max_bytes)
  if (!nzchar(exposure_path)) exposure_path <- selected[[1L]]
  if (!nzchar(outcome_path)) outcome_path <- selected[[2L]]
}
if (normalizePath(exposure_path, mustWork = FALSE) == normalizePath(outcome_path, mustWork = FALSE)) {
  stop("exposure and outcome must be two different gzip files", call. = FALSE)
}
exposure_bytes <- stop_if_bad_file(exposure_path, "exposure")
outcome_bytes <- stop_if_bad_file(outcome_path, "outcome")

read_gwas <- function(path, role) {
  header <- fread(path, nrows = 0L, showProgress = FALSE)
  map <- column_map(names(header), role)
  select <- unname(unlist(map, use.names = FALSE))
  select <- select[!is.na(select)]
  dat <- fread(path, select = select, showProgress = FALSE, nThread = 1L)
  setnames(dat, select, c(
    "SNP", paste0("beta.", role), paste0("se.", role),
    paste0("effect_allele.", role), paste0("other_allele.", role),
    if (!is.na(map$eaf)) paste0("eaf.", role)
  ))
  dat[, SNP := as.character(SNP)]
  dat <- dat[!is.na(SNP) & nzchar(SNP)]
  dat[, (paste0("beta.", role)) := suppressWarnings(as.numeric(get(paste0("beta.", role))))]
  dat[, (paste0("se.", role)) := suppressWarnings(as.numeric(get(paste0("se.", role))))]
  dat[, (paste0("effect_allele.", role)) := toupper(as.character(get(paste0("effect_allele.", role))))]
  dat[, (paste0("other_allele.", role)) := toupper(as.character(get(paste0("other_allele.", role))))]
  eaf_col <- paste0("eaf.", role)
  if (eaf_col %in% names(dat)) dat[, (eaf_col) := suppressWarnings(as.numeric(get(eaf_col)))]
  dat <- dat[!duplicated(SNP)]
  dat
}

overall_start <- proc.time()[["elapsed"]]
read_start <- proc.time()[["elapsed"]]
exposure <- read_gwas(exposure_path, "exposure")
outcome <- read_gwas(outcome_path, "outcome")
read_seconds <- proc.time()[["elapsed"]] - read_start

exposure$id.exposure <- basename(exposure_path)
outcome$id.outcome <- basename(outcome_path)

join_start <- proc.time()[["elapsed"]]
joined <- merge(exposure, outcome, by = "SNP", sort = FALSE, all = FALSE)
join_seconds <- proc.time()[["elapsed"]] - join_start

harm_start <- proc.time()[["elapsed"]]
harmonised <- suppressWarnings(suppressMessages(
  fast_harmonise_data(exposure, outcome, action = 2L)
))
harmonisation_seconds <- proc.time()[["elapsed"]] - harm_start
if (!nrow(harmonised)) stop("harmonisation retained zero SNPs", call. = FALSE)

ivw_start <- proc.time()[["elapsed"]]
ivw <- fast_mr(harmonised, methods = "ivw", nboot = 0L, threads = 1L)
ivw_seconds <- proc.time()[["elapsed"]] - ivw_start
total_seconds <- proc.time()[["elapsed"]] - overall_start

log_attr <- attr(harmonised, "log", exact = TRUE)
harm_summary <- data.table(
  candidate_exposure_rows = nrow(exposure),
  candidate_outcome_rows = nrow(outcome),
  joined_rows = nrow(joined),
  harmonised_rows = nrow(harmonised),
  mr_keep_rows = if ("mr_keep" %in% names(harmonised)) sum(harmonised$mr_keep %in% TRUE, na.rm = TRUE) else NA_integer_,
  removed_rows = if ("mr_keep" %in% names(harmonised)) sum(harmonised$mr_keep %in% FALSE, na.rm = TRUE) else NA_integer_
)
if (is.data.frame(log_attr) && nrow(log_attr)) {
  for (nm in intersect(names(log_attr), c("candidate_variants", "variants_absent_from_reference", "total_variants", "total_variants_for_mr"))) {
    harm_summary[[paste0("fastmr_", nm)]] <- log_attr[[nm]][[1L]]
  }
}

metrics <- data.table(
  measured_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  exposure_file = normalizePath(exposure_path, mustWork = TRUE),
  outcome_file = normalizePath(outcome_path, mustWork = TRUE),
  exposure_bytes = as.numeric(exposure_bytes),
  outcome_bytes = as.numeric(outcome_bytes),
  exposure_rows = nrow(exposure),
  outcome_rows = nrow(outcome),
  joined_rows = nrow(joined),
  harmonised_rows = nrow(harmonised),
  mr_keep_rows = harm_summary$mr_keep_rows,
  read_seconds = read_seconds,
  join_seconds = join_seconds,
  harmonisation_seconds = harmonisation_seconds,
  ivw_seconds = ivw_seconds,
  total_seconds = total_seconds,
  ivw_pairs = 1L,
  mr_rows_per_second = 1 / total_seconds,
  package_fastMR = as.character(packageVersion("fastMR")),
  package_data_table = as.character(packageVersion("data.table")),
  hostname = Sys.info()[["nodename"]]
)

fwrite(metrics, file.path(output_dir, "metrics.tsv"), sep = "\t")
fwrite(harm_summary, file.path(output_dir, "harmonisation_summary.tsv"), sep = "\t")
fwrite(ivw, file.path(output_dir, "mr_result.tsv"), sep = "\t", na = "NA")
fwrite(data.table(
  role = c("exposure", "outcome"),
  path = c(normalizePath(exposure_path, mustWork = TRUE), normalizePath(outcome_path, mustWork = TRUE)),
  bytes = c(as.numeric(exposure_bytes), as.numeric(outcome_bytes)),
  sha256 = c(unname(tools::sha256sum(exposure_path)), unname(tools::sha256sum(outcome_path)))
), file.path(output_dir, "selected_inputs.tsv"), sep = "\t")

message(sprintf(
  paste0("completed: read=%.3fs join=%.3fs harmonisation=%.3fs IVW=%.3fs total=%.3fs ",
         "joined=%d harmonised=%d mr_keep=%d"),
  read_seconds, join_seconds, harmonisation_seconds, ivw_seconds, total_seconds,
  nrow(joined), nrow(harmonised), harm_summary$mr_keep_rows
))
