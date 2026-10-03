#!/usr/bin/env Rscript
# One replicate: every operation in a fresh Rscript under /usr/bin/time -v, formats in random order.
source(file.path(Sys.getenv("STORAGE_SCRIPTS", file.path(Sys.getenv("STORAGE_ROOT", "/user/work/fh6520/showcase/storage"), "scripts")), "common.R"))
REP <- Sys.getenv("REP", Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))
TAG <- Sys.getenv("RESULT_TAG", "rep")                       # "dry" for the dry run
CAP <- as.integer(Sys.getenv("OP_CAP", "600"))
SIZES <- as.integer(strsplit(Sys.getenv("LOOKUP_SIZES", "25,1000,100000"), ",")[[1]])
CPUS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "NA"))
if (is.na(CPUS) || CPUS != 8L) stop("SLURM_CPUS_PER_TASK must be 8 (got ", CPUS, ")")
RES <- file.path(ROOT, "results"); dir.create(RES, FALSE, TRUE)
resfile <- file.path(RES, sprintf("%s%s_ops.csv", TAG, REP))
unlink(resfile)
host <- system("hostname", intern = TRUE)
sh <- function(cmd) tryCatch(paste(system(cmd, intern = TRUE, ignore.stderr = TRUE)[1]), error = function(e) NA)
bi <- tryCatch(paste(capture.output(print(CompreSSoR::compressor_build_info())), collapse = " "), error = function(e) NA)
envrow <- data.table(rep = REP, tag = TAG, host = host, cpu = sh("lscpu | grep 'Model name' | sed 's/.*: *//'"),
  nproc = as.integer(sh("nproc")), slurm_cpus_per_task = CPUS, slurm_job_id = Sys.getenv("SLURM_JOB_ID"),
  slurm_array_task = Sys.getenv("SLURM_ARRAY_TASK_ID"), partition = Sys.getenv("SLURM_JOB_PARTITION"),
  R = R.version.string, arrow = as.character(packageVersion("arrow")), data.table = as.character(packageVersion("data.table")),
  CompreSSoR = as.character(packageVersion("CompreSSoR")), CompreSSoR_lib = LIB_CS,
  gwasvcf = tryCatch(as.character(packageVersion("gwasvcf")), error = function(e) NA),
  bcftools = sh("bcftools --version | head -1"), tabix = sh("tabix --version | head -1"),
  pigz = sh("pigz --version 2>&1 | head -1"), bgzip = sh("bgzip --version | head -1"),
  build_info = bi, date = as.character(Sys.time()))
fwrite(envrow, file.path(RES, sprintf("%s%s_env.csv", TAG, REP)))

ops_for <- function(fmt) {
  S <- list(); add <- function(...) S[[length(S) + 1L]] <<- list(...)
  if (fmt != "vcf") add(op = "write", fmt = fmt, method = "", threads = 8L, size = 0L)
  fr <- switch(fmt, tsv_gz = c("fread_pigz", "arrow_csv"), tsv_bgzip_tabix = "fread_bgzip",
               vcf = "bcftools_fread", parquet = "read_parquet", qparquet = "read_parquet", cpr = "read_sumstats")
  for (m in fr) for (T in c(1L, 8L)) add(op = "fullread", fmt = fmt, method = m, threads = T, size = 0L)
  reg <- switch(fmt, tsv_gz = list(c("scan_per_query", 1L), c("scan_per_query", 8L), c("scan_single", 1L), c("scan_single", 8L)),
                tsv_bgzip_tabix = list(c("tabix", 1L)), vcf = list(c("bcftools", 1L), c("gwasvcf", 1L)),
                parquet = , qparquet = , cpr = list(c("dataset_filter_or_region", 1L), c("dataset_filter_or_region", 8L)))
  for (z in reg) add(op = "region", fmt = fmt, method = if (fmt == "cpr") "read_sumstats_region" else if (z[1] == "dataset_filter_or_region") "open_dataset_filter" else z[1],
                     threads = as.integer(z[2]), size = 0L)
  lk <- switch(fmt, tsv_gz = list(c("scan_join", 1L), c("scan_join", 8L)),
               tsv_bgzip_tabix = list(c("tabix_R", 1L), c("tabix_T", 1L)),
               vcf = list(c("bcftools_R", 1L), c("bcftools_T", 1L), c("gwasvcf", 1L)),
               parquet = , qparquet = list(c("semi_join", 1L), c("semi_join", 8L)),
               cpr = list(c("read_sumstats_variants", 1L), c("read_sumstats_variants", 8L)))
  for (z in lk) for (s in SIZES) add(op = "lookup", fmt = fmt, method = z[1], threads = as.integer(z[2]), size = s)
  S
}
# op.R decides by (fmt, method); map labels to what op.R expects
fix_method <- function(sp) {
  if (sp$op == "region" && sp$fmt %in% c("parquet", "qparquet")) sp$method <- "dataset_filter"
  if (sp$op == "region" && sp$fmt == "cpr") sp$method <- "read_sumstats"
  if (sp$op == "lookup" && sp$fmt == "cpr") sp$method <- "read_sumstats"
  sp
}
set.seed(as.integer(REP) * 7919L + 13L)
order_fmt <- sample(FORMATS)
cat("format order:", order_fmt, "\n"); fwrite(data.table(rep = REP, order = paste(order_fmt, collapse = ",")), file.path(RES, sprintf("%s%s_order.csv", TAG, REP)))
capped <- character(); ALL <- list()
dest <- file.path(ROOT, "work", sprintf("%s%s", TAG, REP)); dir.create(dest, FALSE, TRUE)
for (fmt in order_fmt) {
  for (sp in ops_for(fmt)) {
    label <- sp$method; sp <- fix_method(sp)
    key <- paste(sp$op, sp$fmt, sp$method, sp$threads)
    if (sp$op == "lookup" && key %in% capped) {
      r <- data.table(op = sp$op, format = sp$fmt, method = sp$method, threads = sp$threads, size = sp$size,
                      status = "> cap (skipped: smaller size capped)")
    } else {
      sp$perq <- file.path(RES, "perq", sprintf("%s%s_%s_%s_%s_T%d_%d.csv", TAG, REP, sp$op, sp$fmt, sp$method, sp$threads, sp$size))
      if (sp$op == "write") sp$dest <- dest
      r <- run_op(sp, cap = CAP)
      if (sp$op == "write") unlink(file.path(dest, "*"), recursive = TRUE)
    }
    r[, `:=`(rep = REP, host = host)]
    if (grepl("^> cap", r$status)) capped <- c(capped, key)
    cat(sprintf("[%s] %-9s %-16s %-24s T%d size=%s -> %s  %.2fs rss=%.0fMB\n", format(Sys.time(), "%H:%M:%S"), sp$op, sp$fmt, sp$method, sp$threads,
                sp$size, r$status, r$seconds %||% NA_real_, r$maxrss_mb %||% NA_real_))
    ALL[[length(ALL) + 1L]] <- r
    fwrite(rbindlist(ALL, fill = TRUE), resfile)
  }
}
unlink(dest, recursive = TRUE)
cat("replicate done\n")
