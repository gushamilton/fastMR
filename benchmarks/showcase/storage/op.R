#!/usr/bin/env Rscript
# One benchmark operation in a fresh process. Invoked by run_op() under /usr/bin/time -v.
# Args: op=write|fullread|region|lookup fmt=... method=... threads=T size=N out=file [dest=dir] [perq=file]
args <- commandArgs(TRUE)
SCRIPT_DIR <- Sys.getenv("STORAGE_SCRIPTS", file.path(Sys.getenv("STORAGE_ROOT", "/user/work/fh6520/showcase/storage"), "scripts"))
source(file.path(SCRIPT_DIR, "common.R"))
a <- parse_args(args)
op <- a$op; fmt <- a$fmt; method <- a$method %||% ""; T <- as.integer(a$threads %||% 1L)
size <- as.integer(a$size %||% 0L); out <- a$out
now <- function() proc.time()[[3]]
res <- list(load_s = NA_real_, seconds = NA_real_, n_rows = NA_real_, chk_s = NA_real_, n_q = NA_real_,
            perq_median_s = NA_real_, perq_max_s = NA_real_, bytes = NA_real_,
            index_bytes = NA_real_, note = "")
# Preload only the packages this op needs, OUTSIDE the timed region; load time is its own column.
needs_arrow <- fmt %in% c("parquet", "qparquet") || (fmt == "tsv_gz" && method == "arrow_csv")
t_load <- now()
if (needs_arrow) suppressPackageStartupMessages({library(arrow); library(dplyr)})
if (fmt == "cpr") suppressPackageStartupMessages(library(CompreSSoR))
if (method == "gwasvcf") suppressPackageStartupMessages(library(gwasvcf))
res$load_s <- now() - t_load
nthreads_set(T)
emit <- function() { fwrite(as.data.table(res), out) }
timed_queries <- function(n, f, deadline = Inf, open_s = 0) {
  tt <- numeric(0); rows <- 0; ssum <- 0; t_start <- now()
  for (i in seq_len(n)) {
    t0 <- now(); x <- f(i); tt <- c(tt, now() - t0)
    if (!is.null(x) && nrow(x)) { rows <- rows + nrow(x); ssum <- ssum + sum(as.numeric(x$pos)) }
    if (now() - t_start > deadline) { res$note <<- paste0("partial: ", i, "/", n, " queries before in-process deadline"); break }
  }
  res$seconds <<- open_s + sum(tt); res$n_rows <<- rows; res$chk_s <<- ssum
  res$n_q <<- length(tt); res$perq_median_s <<- median(tt); res$perq_max_s <<- max(tt)
  if (!is.null(a$perq)) { dir.create(dirname(a$perq), FALSE, TRUE); fwrite(data.table(i = seq_along(tt), seconds = tt), a$perq) }
}

# ---------------------------------------------------------------- write
if (op == "write") {
  dest <- a$dest %||% STORES
  dir.create(dest, FALSE, TRUE)
  p <- file.path(dest, basename(store_path(fmt)))
  unlink(c(p, paste0(p, ".tbi")), recursive = TRUE)
  if (fmt == "vcf") {   # bcftools path from the TSV.gz store (no in-memory data)
    tsv <- store_path("tsv_gz")
    hdr <- tempfile(fileext = ".hdr")
    writeLines(c("##fileformat=VCFv4.2", "##FILTER=<ID=PASS,Description=\"All filters passed\">",
      sprintf("##contig=<ID=%d>", 1:22),
      "##FORMAT=<ID=ES,Number=A,Type=Float,Description=\"Effect size estimate relative to the alternative allele\">",
      "##FORMAT=<ID=SE,Number=A,Type=Float,Description=\"Standard error of effect size estimate\">",
      "##FORMAT=<ID=LP,Number=A,Type=Float,Description=\"-log10 p-value for effect estimate\">",
      "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Alternate allele frequency in the association study\">",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFINNGEN"), hdr)
    awkf <- tempfile(fileext = ".awk")
    writeLines("NR>1{lp=($8>0)?-log($8)/log(10):999; printf(\"%s\\t%s\\t.\\t%s\\t%s\\t.\\tPASS\\t.\\tES:SE:LP:AF\\t%s:%s:%.8g:%s\\n\",$1,$2,$3,$4,$5,$6,lp,$7)}", awkf)
    awk <- sprintf("awk -F'\\t' -f %s", awkf)
    n_lim <- if (N_ROWS > 0L) sprintf(" | head -n %d", N_ROWS + 1L) else ""
    cmd <- sprintf("bash -c \"set -o pipefail; (cat %s; pigz -dc -p 2 %s %s | %s) | bcftools view -Oz --threads %d -o %s && bcftools index -t -f %s\"",
                   hdr, shQuote(tsv), n_lim, awk, T, shQuote(p), shQuote(p))
    t0 <- now(); st <- system(cmd); res$seconds <- now() - t0
    if (st != 0L) stop("vcf write failed")
    res$note <- "from TSV.gz via awk|bcftools view -Oz; LP printed %.8g"
  } else {
    d <- load_source(T)
    if (fmt == "tsv_gz") {
      fa <- if ("compressLevel" %in% names(formals(fwrite))) list(compressLevel = 6L) else list()
      t0 <- now(); do.call(fwrite, c(list(d, p, sep = "\t", compress = "gzip", nThread = T), fa)); res$seconds <- now() - t0
      res$note <- "fwrite(compress='gzip', level 6)"
    } else if (fmt == "tsv_bgzip_tabix") {
      tmp <- paste0(p, ".plain.tmp")
      t0 <- now()
      fwrite(d, tmp, sep = "\t", nThread = T)
      system(sprintf("bgzip -@ %d -c %s > %s", T, shQuote(tmp), shQuote(p)))
      system(sprintf("tabix -f -s1 -b2 -e2 -S1 %s", shQuote(p)))
      res$seconds <- now() - t0; unlink(tmp)
      res$note <- "fwrite plain + bgzip -@T + tabix (timed together)"
    } else if (fmt == "parquet" || fmt == "qparquet") {
      suppressPackageStartupMessages(library(arrow))
      dict <- c(chrom = TRUE, pos = FALSE, ref = TRUE, alt = TRUE, beta = FALSE, se = FALSE, eaf = FALSE, p = FALSE)
      if (fmt == "parquet") {
        t0 <- now()
        write_parquet(d, p, compression = "zstd", compression_level = 9L, chunk_size = 1e6L,
                      use_dictionary = dict, write_statistics = TRUE)
        res$seconds <- now() - t0
        res$note <- "zstd-9, 1M-row groups, dict only on chrom/ref/alt, statistics on; no page index (R arrow 22 does not expose it)"
      } else {
        dict <- c(chrom = TRUE, pos = FALSE, ref = TRUE, alt = TRUE, beta = FALSE, se = FALSE, eaf = FALSE, nlp = FALSE)
        t0 <- now()
        tb <- arrow_table(chrom = d$chrom, pos = Array$create(d$pos, type = int32()),
                          ref = d$ref, alt = d$alt,
                          beta = Array$create(d$beta, type = float32()),
                          se = Array$create(d$se, type = float32()),
                          eaf = Array$create(d$eaf, type = float32()),
                          nlp = Array$create(-log10(d$p), type = float32()))
        write_parquet(tb, p, compression = "zstd", compression_level = 9L, chunk_size = 1e6L,
                      use_dictionary = dict, write_statistics = TRUE)
        res$seconds <- now() - t0
        res$note <- "float32 beta/se/eaf, p as float32 -log10, zstd-9, 1M-row groups"
      }
    } else if (fmt == "cpr") {
      suppressPackageStartupMessages(library(CompreSSoR))
      cd <- data.frame(chromosome = d$chrom, base_pair_location = d$pos, reference_allele = d$ref,
                       alternate_allele = d$alt, effect_allele = d$alt, other_allele = d$ref,
                       beta = d$beta, standard_error = d$se, effect_allele_frequency = d$eaf,
                       p_value = d$p)
      rm(d); invisible(gc())
      t0 <- now()
      st <- compress_sumstats(cd, p, input_build = "GRCh38", store_build = "GRCh38", threads = T,
                              pvalue_order = TRUE, overwrite = TRUE)
      res$seconds <- now() - t0
      res$note <- paste0("compress_sumstats defaults, pvalue_order=TRUE, threads=", T)
    } else stop("unknown fmt")
  }
  res$bytes <- path_bytes(p); res$index_bytes <- path_bytes(paste0(p, ".tbi"))
  if (!is.na(res$index_bytes)) res$bytes <- res$bytes
  emit(); quit(status = 0)
}

# ---------------------------------------------------------------- full read
if (op == "fullread") {
  if (fmt == "cpr") suppressPackageStartupMessages(library(CompreSSoR))
  t0 <- now(); x <- read_full(fmt, method, T); res$seconds <- now() - t0
  res$n_rows <- nrow(x); res$chk_s <- sum(as.numeric(x$pos))
  res$note <- paste0("sum_beta=", signif(sum(x$beta), 8)); emit(); quit(status = 0)
}

# ---------------------------------------------------------------- region
if (op == "region") {
  reg <- fread(file.path(LISTS, "regions.tsv"), colClasses = list(character = "chrom"))
  nthreads_set(T)
  rstr <- sprintf("%s:%d-%d", reg$chrom, reg$start, reg$end)
  P <- store_path(fmt)
  if (fmt == "tsv_gz") {
    if (method == "scan_per_query") {
      timed_queries(nrow(reg), function(i) { x <- scan_tsv(P, T); r <- x[chrom == reg$chrom[i] & pos >= reg$start[i] & pos <= reg$end[i]]; rm(x); r },
                    deadline = 540)
    } else if (method == "scan_single") {
      t0 <- now(); x <- scan_tsv(P, T)
      rows <- 0; ssum <- 0
      for (i in seq_len(nrow(reg))) { r <- x[chrom == reg$chrom[i] & pos >= reg$start[i] & pos <= reg$end[i]]; rows <- rows + nrow(r); ssum <- ssum + sum(as.numeric(r$pos)) }
      res$seconds <- now() - t0; res$n_rows <- rows; res$chk_s <- ssum; res$n_q <- nrow(reg)
      res$perq_median_s <- res$seconds / nrow(reg); res$note <- "one full scan serving all 100 windows; per-query = total/100"
    } else stop("method")
  } else if (fmt == "tsv_bgzip_tabix") {
    timed_queries(nrow(reg), function(i) {
      x <- fread(cmd = sprintf("tabix %s %s", shQuote(P), rstr[i]), header = FALSE, showProgress = FALSE, nThread = 1,
                 colClasses = list(character = c(1L, 3L, 4L)))
      if (nrow(x)) setnames(x, LOGICAL); x })
  } else if (fmt == "vcf" && method == "bcftools") {
    timed_queries(nrow(reg), function(i) {
      read_vcf_fread(vcf_cmd(P, sprintf("-r %s", rstr[i])), 1L) })
  } else if (fmt == "vcf" && method == "gwasvcf") {
    suppressPackageStartupMessages(library(gwasvcf))
    timed_queries(nrow(reg), function(i) {
      v <- query_gwas(P, chrompos = rstr[i]); x <- as.data.table(vcf_to_tibble(v))
      if (nrow(x)) x[, pos := start]; x })
    res$note <- paste(res$note, "gwasvcf::query_gwas + vcf_to_tibble (readVcf)")
  } else if (fmt %in% c("parquet", "qparquet")) {
    suppressPackageStartupMessages({library(arrow); library(dplyr)})
    t0 <- now(); ds <- open_dataset(P); open_s <- now() - t0
    timed_queries(nrow(reg), function(i) {
      ch <- reg$chrom[i]; s <- reg$start[i]; e <- reg$end[i]
      as.data.table(ds |> filter(chrom == ch, pos >= s, pos <= e) |> collect()) }, open_s = open_s)
  } else if (fmt == "cpr") {
    suppressPackageStartupMessages(library(CompreSSoR))
    timed_queries(nrow(reg), function(i) {
      x <- read_sumstats(P, region = sprintf("chr%s", rstr[i]), columns = CPR_COLS, threads = T)
      cpr_to_logical(x) })
  } else stop("unknown combination")
  emit(); quit(status = 0)
}

# ---------------------------------------------------------------- lookup
if (op == "lookup") {
  nthreads_set(T)
  keys <- fread(file.path(LISTS, sprintf("keys_%d.tsv", size)), colClasses = list(character = c("chrom", "ref", "alt")))
  P <- store_path(fmt)
  finish <- function(x) { res$n_rows <<- nrow(x); res$chk_s <<- sum(as.numeric(x$pos)) }
  join_keys <- function(x) {  # keep only rows whose full key is requested
    if (!nrow(x)) return(x)
    setnames(x, LOGICAL[1:4], c("chrom", "pos", "ref", "alt"), skip_absent = TRUE)
    x[keys[, .(chrom, pos, ref, alt)], on = c("chrom", "pos", "ref", "alt"), nomatch = NULL]
  }
  sorted_pos <- function() { u <- unique(keys[, .(chrom, pos)]); setorder(u, chrom, pos); u }
  t0 <- now()
  if (fmt == "tsv_gz" && method == "scan_join") {
    x <- scan_tsv(P, T); r <- join_keys(x); finish(r)
  } else if ((fmt == "tsv_bgzip_tabix" || (fmt == "vcf" && method %in% c("bcftools_R", "bcftools_T"))) ) {
    u <- sorted_pos(); rf <- tempfile(fileext = ".tsv"); fwrite(u, rf, sep = "\t", col.names = FALSE)
    flag <- if (grepl("_T$", method)) "-T" else "-R"
    if (fmt == "tsv_bgzip_tabix") {
      x <- fread(cmd = sprintf("tabix %s %s %s", flag, shQuote(rf), shQuote(P)), header = FALSE, showProgress = FALSE,
                 colClasses = list(character = c(1L, 3L, 4L)))
      if (nrow(x)) setnames(x, LOGICAL)
    } else {
      x <- read_vcf_fread(vcf_cmd(P, sprintf("%s %s", flag, shQuote(rf))), 1L)
    }
    r <- join_keys(x); finish(r)
  } else if (fmt == "vcf" && method == "gwasvcf") {
    suppressPackageStartupMessages(library(gwasvcf))
    u <- sorted_pos()
    v <- query_gwas(P, chrompos = sprintf("%s:%d-%d", u$chrom, u$pos, u$pos))
    x <- as.data.table(vcf_to_tibble(v))
    if (nrow(x)) { x[, `:=`(chrom = as.character(seqnames), pos = start)]; setnames(x, c("REF", "ALT"), c("ref", "alt"), skip_absent = TRUE)
      x[, alt := as.character(alt)]; x[, ref := as.character(ref)] }
    r <- join_keys(x); finish(r)
  } else if (fmt %in% c("parquet", "qparquet")) {
    suppressPackageStartupMessages({library(arrow); library(dplyr)})
    ds <- open_dataset(P)
    kt <- arrow_table(chrom = keys$chrom, pos = Array$create(keys$pos, type = int32()), ref = keys$ref, alt = keys$alt)
    r <- as.data.table(ds |> semi_join(kt, by = c("chrom", "pos", "ref", "alt")) |> collect()); finish(r)
  } else if (fmt == "cpr") {
    suppressPackageStartupMessages(library(CompreSSoR))
    r <- cpr_to_logical(read_sumstats(P, variants = mk_key(keys$chrom, keys$pos, keys$ref, keys$alt),
                                      columns = CPR_COLS, threads = T)); finish(r)
  } else stop("unknown combination")
  res$seconds <- now() - t0; res$n_q <- size
  emit(); quit(status = 0)
}
stop("unknown op")
