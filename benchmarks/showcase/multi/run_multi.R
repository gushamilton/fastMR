#!/usr/bin/env Rscript
# One replicate of the multi-trait suite. Cells in fresh processes, formats in random order, 25 min cap.
# A (op, format, method, threads[, size]) series stops once a smaller N hits the cap or is predicted to.
source(file.path(Sys.getenv("MULTI_SCRIPTS", "/user/work/fh6520/showcase/storage/multi/scripts"), "mcommon.R"))
REP <- Sys.getenv("REP", Sys.getenv("SLURM_ARRAY_TASK_ID", "0")); TAG <- Sys.getenv("RESULT_TAG", "rep")
CAP <- as.integer(Sys.getenv("CELL_CAP", "1500"))
NS <- as.integer(strsplit(Sys.getenv("NTRAITS", "1,5,10,20"), ",")[[1]])
KS <- as.integer(strsplit(Sys.getenv("KEY_SIZES", "1000,100000"), ",")[[1]])
FMTS <- strsplit(Sys.getenv("MULTI_FORMATS", "tsv_gz,vcf,cpr"), ",")[[1]]   # e.g. "vcf" to re-run one format
CPUS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "NA"))
if (is.na(CPUS) || CPUS != 8L) stop("SLURM_CPUS_PER_TASK must be 8 (got ", CPUS, ")")
RES <- file.path(MROOT, "results"); dir.create(RES, FALSE, TRUE)
resfile <- file.path(RES, sprintf("%s%s_multi.csv", TAG, REP)); unlink(resfile)
host <- system("hostname", intern = TRUE)
sh <- function(cmd) tryCatch(paste(system(cmd, intern = TRUE, ignore.stderr = TRUE)[1]), error = function(e) NA)
fwrite(data.table(rep = REP, tag = TAG, host = host, cpu = sh("lscpu | grep 'Model name' | sed 's/.*: *//'"),
  slurm_cpus_per_task = CPUS, slurm_job_id = Sys.getenv("SLURM_JOB_ID"), R = R.version.string,
  data.table = as.character(packageVersion("data.table")), CompreSSoR = as.character(packageVersion("CompreSSoR")),
  CompreSSoR_lib = LIB_CS, cpr_dir = CPR_DIR, bcftools = sh("bcftools --version | head -1"),
  pigz = sh("pigz --version 2>&1 | head -1"), date = as.character(Sys.time())), file.path(RES, sprintf("%s%s_multi_env.csv", TAG, REP)))
# Footprint: bytes on disk for N traits per format (deterministic; rep 1 only).
if (REP == "1") {
  fp <- rbindlist(lapply(FMTS, function(f) rbindlist(lapply(NS, function(n) {
    p <- trait_paths(f, n); idx <- if (f == "vcf") paste0(p, ".tbi") else character()
    data.table(format = f, ntraits = n, bytes = sum(vapply(c(p, idx), path_bytes, 0))) }))))
  fp[, bytes_per_trait := bytes / ntraits]; fwrite(fp, file.path(RES, sprintf("%s_multi_footprint.csv", TAG)))
}
series <- function(fmt) {
  S <- list(); add <- function(...) S[[length(S) + 1L]] <<- list(...)
  im <- switch(fmt, tsv_gz = "fread_pigz", vcf = "bcftools_query_i", cpr = c("batch", "loop"))
  for (m in im) for (T in c(1L, 8L)) add(op = "instruments", fmt = fmt, method = m, threads = T, size = 0L)
  em <- switch(fmt, tsv_gz = "fread_join", vcf = c("bcftools_R", "bcftools_T"), cpr = c("batch", "loop"))
  for (k in KS) for (m in em) for (T in c(1L, 8L)) {
    if (m == "bcftools_T" && k < 10000L) next          # -T streams the whole file; -R wins for small K
    add(op = "extract", fmt = fmt, method = m, threads = T, size = k)
  }
  S
}
set.seed(as.integer(REP) * 7919L + 29L)
order_fmt <- sample(FMTS); cat("format order:", order_fmt, "\n")
ALL <- list()
for (fmt in order_fmt) for (sp in series(fmt)) {
  prev <- NULL
  for (n in NS) {
    sp$ntraits <- n
    sp$keep <- if (n == max(NS) && REP == "1") file.path(MROOT, "work", sprintf("%s_%s_%s_%s_T%d_K%d.rds", TAG, sp$op, sp$fmt, sp$method, sp$threads, sp$size)) else ""
    if (nzchar(sp$keep)) dir.create(dirname(sp$keep), FALSE, TRUE)
    skip <- NULL
    if (!is.null(prev)) {
      if (grepl("^> cap", prev$status)) skip <- "> cap (skipped: smaller N capped)"
      else if (prev$status == "ok" && prev$wall_s * n / prev$ntraits > CAP) skip <- "> cap (predicted)"
    }
    r <- if (!is.null(skip)) data.table(op = sp$op, format = sp$fmt, method = sp$method, threads = sp$threads,
                                        ntraits = n, size = sp$size, status = skip)
         else run_mop(sp, cap = CAP)
    r[, `:=`(rep = REP, host = host)]
    cat(sprintf("[%s] %-11s %-6s %-16s T%d N=%-2d K=%-6d -> %s %.2fs rss=%.0fMB\n", format(Sys.time(), "%H:%M:%S"), sp$op, sp$fmt,
                sp$method, sp$threads, n, sp$size, r$status, r$wall_s %||% NA_real_, r$maxrss_mb %||% NA_real_))
    ALL[[length(ALL) + 1L]] <- r; fwrite(rbindlist(ALL, fill = TRUE), resfile)
    prev <- r
  }
}
cat("replicate done\n")
