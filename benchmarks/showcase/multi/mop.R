#!/usr/bin/env Rscript
# One multi-trait cell in a fresh process (under /usr/bin/time -v from run_mop()).
# op=instruments|extract fmt=tsv_gz|vcf|cpr method=... threads=T ntraits=N size=K out=file
source(file.path(Sys.getenv("MULTI_SCRIPTS", "/user/work/fh6520/showcase/storage/multi/scripts"), "mcommon.R"))
a <- parse_args(commandArgs(TRUE))
op <- a$op; fmt <- a$fmt; method <- a$method; T <- as.integer(a$threads); N <- as.integer(a$ntraits); K <- as.integer(a$size)
now <- function() proc.time()[[3]]
t_load <- now()
if (fmt == "cpr") suppressPackageStartupMessages(library(CompreSSoR))
suppressPackageStartupMessages(library(parallel))
load_s <- now() - t_load
setDTthreads(T)
P <- trait_paths(fmt, N)
stopifnot(all(file.exists(P)))
THR <- 5e-8
# Per-file parallelism for the tools without a native multi-file reader: 8 cores = up to 8 files at once.
mc <- function(f) {
  if (T == 1L || N == 1L) return(lapply(seq_along(P), f))
  mclapply(seq_along(P), f, mc.cores = min(T, N), mc.preschedule = FALSE)
}
inner_T <- if (N == 1L) T else max(1L, T %/% min(T, N))
keys <- NULL
if (op == "extract") {
  keys <- fread(file.path(MROOT, "lists", sprintf("keys_%d.tsv", K)), colClasses = list(character = c("chrom", "ref", "alt")))
  rf <- tempfile(fileext = ".tsv"); u <- unique(keys[, .(chrom, pos)]); setorder(u, chrom, pos)
  fwrite(u, rf, sep = "\t", col.names = FALSE)
}
join_keys <- function(x) x[keys, on = c("chrom", "pos", "ref", "alt"), nomatch = NULL]
t0 <- now()
out <- if (fmt == "tsv_gz") {
  mc(function(i) { x <- read_sim_tsv(P[i], inner_T); if (op == "instruments") x[p < THR] else join_keys(x) })
} else if (fmt == "vcf") {
  flag <- if (method == "bcftools_T") "-T" else "-R"
  mc(function(i) {
    cmd <- if (op == "instruments") vcf_cmd(P[i], sprintf("-i 'FORMAT/LP>%.10f'", -log10(THR)))
           else vcf_cmd(P[i], sprintf("%s %s", flag, shQuote(rf)))
    x <- read_vcf_fread(cmd, 1L); if (op == "extract" && nrow(x)) join_keys(x) else x })
} else if (fmt == "cpr") {
  cols <- CPR_COLS
  if (op == "instruments") {
    if (method == "batch") {
      r <- read_candidates_batch(P, THR, columns = cols, threads = T, strategy = "pvalue_flag")
      lapply(r, cpr_to_logical)
    } else lapply(P, function(s) cpr_to_logical(read_candidates(s, THR, columns = cols, threads = T, strategy = "pvalue_flag")))
  } else {
    v <- mk_key(keys$chrom, keys$pos, keys$ref, keys$alt)
    if (method == "batch") lapply(read_sumstats_batch(P, variants = v, columns = cols, threads = T), cpr_to_logical)
    else lapply(P, function(s) cpr_to_logical(read_sumstats(s, variants = v, columns = cols, threads = T)))
  }
} else stop("fmt")
secs <- now() - t0
bad <- vapply(out, function(x) inherits(x, "try-error") || inherits(x, "error") || !is.data.frame(x), logical(1))
if (any(bad)) stop("worker failed: ", paste(head(as.character(out[bad]), 2), collapse = " | "))
n_rows <- sum(vapply(out, nrow, 0)); chk <- sum(vapply(out, function(x) sum(as.numeric(x$pos)), 0))
# Per-trait result tables for cross-format agreement (only the largest N, to keep disk small).
if (!is.null(a$keep) && nzchar(a$keep)) saveRDS(setNames(out, TRAIT_IDS[seq_len(N)]), a$keep)
fwrite(data.table(load_s = load_s, seconds = secs, n_rows = n_rows, chk_s = chk,
                  inner_threads = if (fmt == "cpr") T else inner_T,
                  note = sprintf("sum_beta=%s", signif(sum(vapply(out, function(x) sum(x$beta), 0)), 10))), a$out)
