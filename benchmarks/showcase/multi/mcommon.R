# Multi-trait storage suite: N simulated GWAS (same 8.97M-variant set) per format.
source(file.path(Sys.getenv("STORAGE_SCRIPTS", "/user/work/fh6520/showcase/storage/scripts"), "common.R"))
MROOT <- Sys.getenv("MULTI_ROOT", "/user/work/fh6520/showcase/storage/multi")
TR <- Sys.getenv("TRAITS", "/user/work/fh6520/showcase/prep/traits")
CPR_DIR <- Sys.getenv("CPR_DIR", file.path(TR, "cpr_new"))
TRAIT_IDS <- c(sprintf("exp%02d", 1:10), sprintf("out%02d", 1:10))
trait_paths <- function(fmt, n) {
  id <- TRAIT_IDS[seq_len(n)]
  switch(fmt, tsv_gz = file.path(TR, "tsv", paste0(id, ".tsv.gz")),
         vcf = file.path(TR, "vcf", paste0(id, ".vcf.gz")),
         cpr = file.path(CPR_DIR, paste0(id, ".cpr")))
}
SIM_TSV_COLS <- c("chromosome", "base_pair_location", "effect_allele", "other_allele", "beta",
                  "standard_error", "effect_allele_frequency", "p_value")
read_sim_tsv <- function(path, T) {   # -> logical columns chrom pos ref alt beta se eaf p
  x <- fread(cmd = sprintf("pigz -dc -p %d %s", max(1L, T), shQuote(path)), nThread = T, showProgress = FALSE,
             select = SIM_TSV_COLS, colClasses = list(character = c("chromosome", "effect_allele", "other_allele")))
  setnames(x, SIM_TSV_COLS, c("chrom", "pos", "alt", "ref", "beta", "se", "eaf", "p"))
  setcolorder(x, LOGICAL); x
}
run_mop <- function(spec, cap = 1500L, logdir = file.path(MROOT, "logs")) {
  dir.create(logdir, showWarnings = FALSE, recursive = TRUE)
  base <- file.path(logdir, paste0("mop_", Sys.getpid(), "_", as.integer(runif(1) * 1e9)))
  outf <- paste0(base, ".csv"); tf <- paste0(base, ".time"); lf <- paste0(base, ".log")
  spec$out <- outf
  args <- paste0(names(spec), "=", vapply(spec, as.character, ""), collapse = " ")
  cmd <- sprintf("env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 /usr/bin/time -v -o %s timeout -k 10 %d Rscript --no-save %s/mop.R %s > %s 2>&1",
                 tf, cap, Sys.getenv("MULTI_SCRIPTS", file.path(MROOT, "scripts")), args, lf)
  st <- system(cmd); tv <- parse_time_v(tf)
  row <- if (file.exists(outf)) fread(outf, colClasses = list(character = "note")) else data.table()
  capped <- (st == 124L || (!is.na(tv$exit) && tv$exit == 124))
  status <- if (capped) "> cap" else if (st != 0L || !nrow(row)) "error" else "ok"
  res <- data.table(op = spec$op, format = spec$fmt, method = spec$method, threads = as.integer(spec$threads),
                    ntraits = as.integer(spec$ntraits), size = as.integer(spec$size), status = status,
                    wall_s = tv$wall, user_s = tv$user, sys_s = tv$sys, maxrss_mb = tv$maxrss_kb / 1024)
  if (nrow(row)) res <- cbind(res, row[, setdiff(names(row), names(res)), with = FALSE])
  if (status == "error") res$note <- paste(tail(readLines(lf), 3), collapse = " | ")
  if (status == "ok") unlink(c(tf, lf, outf))
  res
}
