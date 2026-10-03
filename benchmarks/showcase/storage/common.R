# Shared helpers for the storage-format benchmark (showcase section 1).
ROOT <- Sys.getenv("STORAGE_ROOT", "/user/work/fh6520/showcase/storage")
SRC <- Sys.getenv("SRC_TSV", "/user/work/fh6520/CompreSSoR-bp-thread-test/external/benchmark-10m/finngen_10m_snps.tsv.gz")
LIB_CS <- Sys.getenv("LIB_CS", "/user/work/fh6520/regress-audit/lib-cs-final")
LIB_ST <- Sys.getenv("LIB_STORAGE", "/user/work/fh6520/showcase/lib-storage")
N_ROWS <- as.integer(Sys.getenv("N_ROWS", "0"))      # 0 = all rows (testing only)
.libPaths(c(LIB_ST, LIB_CS, .libPaths()))
STORES <- file.path(ROOT, "stores")
LISTS <- file.path(ROOT, "stores", "lists")
SCRIPTS <- Sys.getenv("STORAGE_SCRIPTS", file.path(ROOT, "scripts"))
LOGICAL <- c("chrom", "pos", "ref", "alt", "beta", "se", "eaf", "p")
FORMATS <- c("tsv_gz", "tsv_bgzip_tabix", "vcf", "parquet", "qparquet", "cpr")
suppressPackageStartupMessages(library(data.table))

store_path <- function(fmt) file.path(STORES, switch(fmt,
  tsv_gz = "finngen.tsv.gz", tsv_bgzip_tabix = "finngen.bgz.tsv.gz",
  vcf = "finngen.vcf.gz", parquet = "finngen.parquet",
  qparquet = "finngen.q.parquet", cpr = "finngen.cpr"))
index_path <- function(fmt) switch(fmt,
  tsv_bgzip_tabix = paste0(store_path(fmt), ".tbi"),
  vcf = paste0(store_path(fmt), ".tbi"), NA_character_)
path_bytes <- function(p) {
  if (is.na(p) || !file.exists(p)) return(NA_real_)
  if (dir.exists(p)) sum(file.info(list.files(p, full.names = TRUE, recursive = TRUE))$size)
  else file.info(p)$size
}
nthreads_set <- function(T) {
  data.table::setDTthreads(T)
  if (isNamespaceLoaded("arrow")) {   # never load arrow here: callers preload it outside the timer
    # Arrow hangs with a 1-thread pool for datasets/CSV ("use num_threads >= 2"), so for T = 1
    # the pools keep 2 threads but all scans/reads run with use_threads = FALSE (serial).
    arrow::set_cpu_count(max(T, 2L)); arrow::set_io_thread_count(max(T, 2L))
    options(arrow.use_threads = T > 1L)
  }
}
load_source <- function(T = 8L) {
  data.table::setDTthreads(T)
  d <- fread(cmd = paste("pigz -dc", shQuote(SRC)), nThread = T, showProgress = FALSE,
             colClasses = list(character = c("chrom", "ref", "alt")),
             nrows = if (N_ROWS > 0L) N_ROWS else Inf)
  setnames(d, LOGICAL)
  # canonical genomic order (numeric chromosome, pos, ref, alt) for every format; the source file is
  # in lexicographic chromosome order, which is not the order a .cpr store can hold.
  d[, .ci := as.integer(chrom)]; setorderv(d, c(".ci", "pos", "ref", "alt")); d[, .ci := NULL]
  d
}
mk_key <- function(chrom, pos, ref, alt) paste(chrom, pos, ref, alt, sep = ":")
chk_of <- function(x) {            # cheap content checksum on identity columns
  if (is.null(x) || !nrow(x)) return(c(n = 0, s = 0))
  c(n = nrow(x), s = sum(as.numeric(x$pos)))
}
CPR_COLS <- c("chromosome", "base_pair_location", "other_allele", "effect_allele",
              "beta", "standard_error", "effect_allele_frequency", "p_value")
cpr_to_logical <- function(x) {
  setDT(x)
  setnames(x, CPR_COLS, c("chrom", "pos", "ref", "alt", "beta", "se", "eaf", "p"))
  x[, c(LOGICAL), with = FALSE]
}
TSV_CLS <- list(character = c("chrom", "ref", "alt"))
VCF_FMT <- "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%ES\\t%SE\\t%AF\\t%LP]\\n"
vcf_cmd <- function(path, extra = "") sprintf("bcftools query %s -f '%s' %s", extra, VCF_FMT, shQuote(path))
vcf_to_logical <- function(x) {
  setnames(x, c("chrom", "pos", "ref", "alt", "beta", "se", "eaf", "lp"))
  x[, p := 10^(-lp)][, lp := NULL]
  setcolorder(x, LOGICAL)
  x
}
read_vcf_fread <- function(cmd, T) {
  x <- fread(cmd = cmd, header = FALSE, nThread = T, showProgress = FALSE,
             colClasses = list(character = c(1L, 3L, 4L)))
  vcf_to_logical(x)
}
arrow_tsv_schema <- function() {
  arrow::schema(chrom = arrow::string(), pos = arrow::int32(), ref = arrow::string(),
                alt = arrow::string(), beta = arrow::float64(), se = arrow::float64(),
                eaf = arrow::float64(), p = arrow::float64())
}
scan_tsv <- function(path, T, how = "pigz") {
  cmd <- if (how == "bgzip") sprintf("bgzip -@ %d -dc %s", T, shQuote(path))
         else sprintf("pigz -dc -p %d %s", T, shQuote(path))
  x <- fread(cmd = cmd, nThread = T, showProgress = FALSE, colClasses = TSV_CLS)
  setnames(x, LOGICAL); x
}

# Full read into a data.table with the logical columns.
read_full <- function(fmt, method, T) {
  nthreads_set(T)
  if (fmt == "tsv_gz" && method == "fread_pigz") return(scan_tsv(store_path(fmt), T, "pigz"))
  if (fmt == "tsv_gz" && method == "arrow_csv") {
    x <- arrow::read_csv_arrow(store_path(fmt), schema = arrow_tsv_schema(), skip = 1,
                               parse_options = arrow::csv_parse_options(delimiter = "\t"),
                               as_data_frame = TRUE)
    return(setDT(x))
  }
  if (fmt == "tsv_bgzip_tabix" && method == "fread_bgzip") return(scan_tsv(store_path(fmt), T, "bgzip"))
  if (fmt == "vcf" && method == "bcftools_fread") return(read_vcf_fread(vcf_cmd(store_path(fmt)), T))
  if (fmt %in% c("parquet", "qparquet") && method == "read_parquet") {
    x <- setDT(arrow::read_parquet(store_path(fmt), as_data_frame = TRUE))
    if (fmt == "qparquet") {
      x[, p := 10^(-as.numeric(nlp))][, nlp := NULL]
      x[, `:=`(beta = as.numeric(beta), se = as.numeric(se), eaf = as.numeric(eaf))]
      setcolorder(x, LOGICAL)
    }
    return(x)
  }
  if (fmt == "cpr" && method == "read_sumstats") {
    x <- CompreSSoR::read_sumstats(store_path(fmt), columns = CPR_COLS, threads = T)
    return(cpr_to_logical(x))
  }
  stop("unknown read ", fmt, "/", method)
}

# ---- orchestration: run one op in a fresh Rscript under /usr/bin/time -v -------------
parse_time_v <- function(f) {
  l <- if (file.exists(f)) readLines(f) else character()
  g <- function(pat) { m <- grep(pat, l, value = TRUE); if (length(m)) sub("^.*: ", "", m[1]) else NA }
  el <- g("Elapsed \\(wall")
  wall <- NA_real_
  if (!is.na(el)) { p <- as.numeric(strsplit(el, ":")[[1]]); wall <- sum(p * c(3600, 60, 1)[(4 - length(p)):3]) }
  list(wall = wall, user = as.numeric(g("User time")), sys = as.numeric(g("System time")),
       maxrss_kb = as.numeric(g("Maximum resident set size")),
       exit = as.numeric(g("Exit status")))
}
run_op <- function(spec, cap = 600L, rscript = "Rscript", logdir = file.path(ROOT, "logs"),
                   env_extra = "") {
  dir.create(logdir, showWarnings = FALSE, recursive = TRUE)
  base <- file.path(logdir, paste0("op_", Sys.getpid(), "_", as.integer(runif(1) * 1e9)))
  outf <- paste0(base, ".csv"); tf <- paste0(base, ".time"); lf <- paste0(base, ".log")
  spec$out <- outf
  args <- paste0(names(spec), "=", vapply(spec, as.character, ""), collapse = " ")
  cmd <- sprintf("env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 %s /usr/bin/time -v -o %s timeout -k 10 %d %s --no-save %s/op.R %s > %s 2>&1",
                 env_extra, tf, cap, rscript, SCRIPTS, args, lf)
  st <- system(cmd)
  tv <- parse_time_v(tf)
  row <- if (file.exists(outf)) fread(outf, colClasses = list(character = "note")) else data.table()
  capped <- (st == 124L || (!is.na(tv$exit) && tv$exit == 124))
  status <- if (capped) "> cap" else if (st != 0L || !nrow(row)) "error" else "ok"
  res <- data.table(op = spec$op, format = spec$fmt, method = spec$method %||% "",
                    threads = as.integer(spec$threads %||% 1L), size = as.integer(spec$size %||% NA),
                    status = status, wall_s = tv$wall, user_s = tv$user, sys_s = tv$sys,
                    maxrss_mb = tv$maxrss_kb / 1024)
  if (nrow(row)) res <- cbind(res, row[, setdiff(names(row), names(res)), with = FALSE])
  if (status == "error") { res$note <- paste(tail(readLines(lf), 3), collapse = " | ") }
  if (status == "ok") unlink(c(tf, lf, outf))
  res
}
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
parse_args <- function(a) { kv <- strsplit(a, "=", fixed = TRUE); setNames(lapply(kv, function(z) paste(z[-1], collapse = "=")), vapply(kv, `[`, "", 1)) }
write_result <- function(out, ...) {
  r <- as.data.table(list(...)); fwrite(r, out)
}
