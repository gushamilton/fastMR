#!/usr/bin/env Rscript
# MR-compute benchmark worker: ONE scenario x ONE arm x ONE size in a fresh R process.
# Writes one CSV row (--cellcsv) and a standardised result RDS (--rds dir) used by compare.R.
# Libraries come from LIB_FASTMR / LIB_TSMR env vars. Data are SYNTHETIC and on the FORWARD strand
# (harmonisation uses action = 1) except the explicitly labelled harmonise "action2_strandflips" variant.
args <- commandArgs(trailingOnly = TRUE)
o <- list(scen = "single", size = "10", arm = "tsmr1", set = "default", rep = "1", rds = tempdir(),
          cellcsv = NA, ref = NA, pref = NA, p1bin = NA, p2bin = NA)
for (a in args) { kv <- regmatches(a, regexec("^--([^=]+)=(.*)$", a))[[1]]; if (length(kv) == 3) o[[kv[2]]] <- kv[3] }
size <- as.numeric(o$size)

# ---- arm table: arm -> package, threads, bootstrap seed
ARMS <- list(tsmr1 = list(pkg = "tsmr", thr = 1L, seed = 1L), tsmr8 = list(pkg = "tsmr", thr = 8L, seed = 1L),
             tsmr1_seedB = list(pkg = "tsmr", thr = 1L, seed = 2L),
             p2c1 = list(pkg = "plink2", thr = 1L, seed = 1L),
             fast1 = list(pkg = "fastmr", thr = 1L, seed = 1L), fast8 = list(pkg = "fastmr", thr = 8L, seed = 1L))
arm <- ARMS[[o$arm]]; if (is.null(arm)) stop("unknown arm")
threads <- arm$thr; is_tsmr <- arm$pkg == "tsmr"; is_fast <- arm$pkg == "fastmr"

# ---- thread control (explicit everywhere) and allocation assertion
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK", "")
if (slurm_cpus != "8") stop("SLURM_CPUS_PER_TASK must be 8, got '", slurm_cpus, "'")
nproc_avail <- as.integer(system("env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc", intern = TRUE))
if (nproc_avail < threads) stop("nproc (", nproc_avail, ") < arm threads")
# ---- libraries
lib_fm <- Sys.getenv("LIB_FASTMR", ""); lib_ts <- Sys.getenv("LIB_TSMR", "")
if (is_fast && nzchar(lib_fm)) .libPaths(c(lib_fm, .libPaths()))
if (!is_fast && nzchar(lib_ts)) .libPaths(c(lib_ts, .libPaths()))
suppressPackageStartupMessages({ library(data.table); if (is_tsmr) library(TwoSampleMR) else if (is_fast) library(fastMR) })
data.table::setDTthreads(if (is_fast) threads else 1L)   # TSMR arms: 1 (x8 parallelism is mclapply forks)
pkg_label <- c(tsmr = "TwoSampleMR", plink2 = "plink2", fastmr = "fastMR")[[arm$pkg]]
pkg_version <- switch(arm$pkg, tsmr = as.character(packageVersion("TwoSampleMR")),
                      fastmr = paste0(packageVersion("fastMR"), "+", Sys.getenv("FASTMR_SHA", "?")), plink2 = "plink2")
cpu <- local({ l <- grep("^model name", readLines("/proc/cpuinfo", warn = FALSE), value = TRUE); trimws(sub("^[^:]*:", "", l[1])) })
host <- Sys.info()[["nodename"]]
resf <- file.path(o$rds, paste0(paste(o$scen, o$size, o$set, o$arm, sep = "_"), ".rds"))

emit <- function(wall, cpu_s, nreps, notes = "") {
  row <- data.frame(scenario = o$scen, size = size, method_set = o$set, arm = o$arm, package = pkg_label, threads = threads,
    replicate = as.integer(o$rep), hostname = host, cpu_model = cpu, nproc = nproc_avail, slurm_cpus = slurm_cpus,
    status = "ok", wall_s = wall, cpu_s = cpu_s, timing_reps = nreps, notes = notes, pkg_version = pkg_version,
    boot_seed = arm$seed, stringsAsFactors = FALSE)
  write.csv(row, o$cellcsv, row.names = FALSE)
  cat(sprintf("[%s %s %s %s] wall=%.4g cpu=%.4g reps=%d %s\n", o$scen, o$size, o$set, o$arm, wall, cpu_s, nreps, notes))
}
cpu_now <- function() { p <- proc.time(); unname(p[1] + p[2] + p[4] + p[5]) }   # user+sys, self + reaped children (mclapply, plink)
# Fast calls repeated to >= 2 s (>= 3 reps, <= 300); calls whose first run exceeds 30 s are timed once. Median wall.
time_calls <- function(f, min_total = 2, min_reps = 3, max_reps = 300, once_above = 30) {
  ts <- numeric(); val <- NULL; c0 <- cpu_now()
  repeat {
    t0 <- proc.time()[["elapsed"]]; val <- f(); ts <- c(ts, proc.time()[["elapsed"]] - t0)
    if (ts[1] > once_above || length(ts) >= max_reps || (sum(ts) >= min_total && length(ts) >= min_reps)) break
  }
  list(wall = median(ts), cpu = (cpu_now() - c0) / length(ts), n = length(ts), val = val)
}
q <- function(f) suppressWarnings(suppressMessages(f()))
dothread <- function(d, key, nch = 8L) {   # interleave pairs over nch chunks
  u <- unique(key); split(d, factor((match(key, u) - 1L) %% nch, levels = 0:(nch - 1L)), drop = TRUE)
}
# TSMR x8: mclapply over chunks of pairs (split + fork + rbind all inside the timed region).
par_apply <- function(d, key, fun, seed) {
  if (threads == 1L) { set.seed(seed); return(fun(d)) }
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  ch <- dothread(d, key, threads)
  res <- parallel::mclapply(ch, fun, mc.cores = threads, mc.preschedule = TRUE)
  bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), NA)
  if (any(bad)) stop("mclapply failure: ", paste(res[bad][[1]], collapse = " "))
  do.call(rbind, res)
}

# ------------------------------------------------------------------ generators (deterministic, shared across arms)
gen_pairs <- function(P, seed = 11L, mean_k = 20L, fixed_k = NULL) {
  set.seed(seed)
  k <- if (is.null(fixed_k)) pmax(5L, rpois(P, mean_k)) else rep(as.integer(fixed_k), P)
  n <- sum(k); g <- rep(seq_len(P), k)
  theta <- rnorm(P, 0, 0.3)
  bx <- sample(c(-1, 1), n, TRUE) * runif(n, 0.03, 0.12)
  sx <- runif(n, 0.004, 0.01); sy <- runif(n, 0.008, 0.02)
  by <- theta[g] * bx + rnorm(n, 0, 0.004) + rnorm(n, 0, sy)   # mild balanced pleiotropy
  data.frame(SNP = paste0("rs", sample.int(1e9L, n)),
    beta.exposure = bx, se.exposure = sx, beta.outcome = by, se.outcome = sy,
    id.exposure = paste0("exp", g), id.outcome = paste0("out", g), exposure = paste0("exp", g), outcome = paste0("out", g),
    mr_keep = TRUE, pval.exposure = 2 * pnorm(-abs(bx / sx)), pval.outcome = 2 * pnorm(-abs(by / sy)),
    samplesize.exposure = 50000, samplesize.outcome = 300000, units.exposure = "SD", units.outcome = "SD",
    eaf.exposure = runif(n, 0.1, 0.9), eaf.outcome = runif(n, 0.1, 0.9), stringsAsFactors = FALSE)
}
pkey <- function(d) paste(d$id.exposure, d$id.outcome, sep = "|")
std_mr <- function(r) { r <- as.data.frame(r); data.frame(key = paste(r$id.exposure, r$id.outcome, sep = "|"), method = r$method,
  b = r$b, se = r$se, pval = r$pval, stringsAsFactors = FALSE) }
save_std <- function(s) saveRDS(s, resf)
default_fast <- c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode")

# ------------------------------------------------------------------ MR (single pair / many pairs)
run_mr <- function(d) {
  force(d)   # must be generated BEFORE any RNGkind change in par_apply
  isivw <- o$set == "ivw"
  if (is_tsmr) {
    fun <- if (isivw) function(x) q(function() mr(x, method_list = "mr_ivw")) else function(x) q(function() mr(x))
    f <- function() par_apply(d, pkey(d), fun, arm$seed)
  } else {
    f <- if (isivw) function() q(function() fast_mr(d, methods = "ivw", nboot = 0, threads = threads))
         else function() q(function() fast_mr(d, methods = default_fast, nboot = 1000, seed = arm$seed, threads = threads))
  }
  tm <- time_calls(f); save_std(std_mr(tm$val))
  emit(tm$wall, tm$cpu, tm$n, sprintf("pairs=%d;rows=%d", length(unique(pkey(d))), nrow(d)))
}
kfix <- function() { m <- regmatches(o$set, regexec("^default_k([0-9]+)$", o$set))[[1]]; if (length(m)) as.integer(m[2]) else NULL }
sc_single <- function() run_mr(gen_pairs(1L, fixed_k = size))
sc_many <- function() run_mr(gen_pairs(size, fixed_k = kfix()))

# ------------------------------------------------------------------ harmonisation
# forward_only=TRUE (action = 1 runs): same / swapped alleles only, no strand flips. Palindromic SNPs (A/T, C/G) are kept
# (action = 1 assumes forward strand, so palindromes are aligned by allele identity, as in an on-strand pipeline).
# forward_only=FALSE (labelled variant, action = 2): adds strand flips (+ swap) with allele-frequency palindrome inference.
gen_harm <- function(n, seed = 5L, forward_only = TRUE, id_e = "E1", id_o = "O1") {
  set.seed(seed)
  pairs <- list(c("A","G"), c("A","C"), c("C","T"), c("G","T"), c("A","T"), c("C","G"))
  pi <- sample(6, n, TRUE, prob = c(.22,.22,.22,.22,.06,.06))
  ea <- vapply(pairs, `[`, "", 1)[pi]; oa <- vapply(pairs, `[`, "", 2)[pi]
  comp <- c(A="T", T="A", C="G", G="C")
  eaf <- runif(n, 0.05, 0.95); pal <- pi >= 5; eaf[pal] <- ifelse(runif(sum(pal)) < .7, runif(sum(pal), .05, .3), runif(sum(pal), .4, .6))
  snp <- paste0("rs", seq_len(n)); b <- rnorm(n, 0, .05); se <- runif(n, .005, .02)
  list(snp = snp, ea = ea, oa = oa, eaf = eaf, b = b, se = se, comp = comp, n = n, forward_only = forward_only, id_e = id_e, id_o = id_o)
}
mk_exp <- function(u, idx, id, seed) {
  set.seed(seed); b <- rnorm(length(idx), 0, .05); se <- runif(length(idx), .005, .02)
  data.frame(SNP = u$snp[idx], beta.exposure = b, se.exposure = se, effect_allele.exposure = u$ea[idx], other_allele.exposure = u$oa[idx],
    eaf.exposure = u$eaf[idx], pval.exposure = 2*pnorm(-abs(b/se)), samplesize.exposure = 50000, id.exposure = id, exposure = id,
    mr_keep.exposure = TRUE, stringsAsFactors = FALSE)
}
mk_out <- function(u, keep, bx_by_snp, seed) {   # outcome rows for SNP indices `keep`; beta from bx_by_snp
  set.seed(seed); m <- length(keep)
  mode <- if (u$forward_only) sample(2, m, TRUE, prob = c(.65, .35)) else sample(4, m, TRUE, prob = c(.55, .3, .1, .05))
  oea <- u$ea[keep]; ooa <- u$oa[keep]; sgn <- rep(1, m); oeaf <- u$eaf[keep] + rnorm(m, 0, .02)
  sw <- mode %in% c(2, 4); tmp <- oea[sw]; oea[sw] <- ooa[sw]; ooa[sw] <- tmp; sgn[sw] <- -1; oeaf[sw] <- 1 - oeaf[sw]
  st <- mode %in% c(3, 4); oea[st] <- u$comp[oea[st]]; ooa[st] <- u$comp[ooa[st]]
  ob <- sgn * bx_by_snp * 0.5 + rnorm(m, 0, .01); ose <- runif(m, .005, .02)
  data.frame(SNP = u$snp[keep], beta.outcome = ob, se.outcome = ose, effect_allele.outcome = oea, other_allele.outcome = ooa,
    eaf.outcome = pmin(pmax(oeaf, .001), .999), pval.outcome = 2*pnorm(-abs(ob/ose)), samplesize.outcome = 300000,
    id.outcome = u$id_o, outcome = u$id_o, mr_keep.outcome = TRUE, stringsAsFactors = FALSE)
}
harm_action <- function() if (o$set == "action2_strandflips") 2L else 1L
harm_fun <- function(action) if (is_tsmr) function(e, ot) q(function() harmonise_data(e, ot, action = action)) else
  function(e, ot) q(function() fast_harmonise_data(e, ot, action = action))
std_harm <- function(r) { r <- as.data.frame(r); data.frame(key = paste(r$id.exposure, r$id.outcome, r$SNP, sep = "|"), method = "all",
  beta.outcome = r$beta.outcome, beta.exposure = r$beta.exposure, eaf.outcome = r$eaf.outcome, flag_mr_keep = as.numeric(r$mr_keep), stringsAsFactors = FALSE) }
sc_harmonise <- function() {   # single pair, `size` rows
  fo <- o$set != "action2_strandflips"
  u <- gen_harm(size, forward_only = fo); e <- mk_exp(u, seq_len(size), "E1", 6L)
  keep <- { set.seed(7L); sample(size, round(.9 * size)) }
  ot <- mk_out(u, keep, e$beta.exposure[keep], 8L); set.seed(9L); ot <- ot[sample(nrow(ot)), ]
  h <- harm_fun(harm_action()); tm <- time_calls(function() h(e, ot))
  save_std(std_harm(tm$val)); r <- tm$val
  emit(tm$wall, tm$cpu, tm$n, sprintf("action=%d;rows_in=%d;rows_out=%d;kept=%d", harm_action(), size, nrow(r), sum(r$mr_keep)))
}
sc_harmonise_many <- function() {   # `size` pairs (exposures) x 30 SNPs against ONE outcome, forward strand, action = 1
  P <- size; Uu <- 50000L; u <- gen_harm(Uu, seed = 5L)
  set.seed(12L); idx <- lapply(seq_len(P), function(i) sort(sample(Uu, 30L)))
  e <- do.call(rbind, lapply(seq_len(P), function(i) mk_exp(u, idx[[i]], paste0("E", i), 100L + i)))
  keep <- { set.seed(13L); sample(Uu, round(.9 * Uu)) }
  ot <- mk_out(u, keep, rep(0.03, length(keep)), 14L)
  h <- harm_fun(1L)
  f <- if (is_tsmr) function() par_apply(e, e$id.exposure, function(x) h(x, ot), 1L) else function() h(e, ot)
  tm <- time_calls(f); save_std(std_harm(tm$val)); r <- tm$val
  emit(tm$wall, tm$cpu, tm$n, sprintf("action=1;pairs=%d;snps_per_pair=30;rows_in=%d;rows_out=%d;kept=%d", P, nrow(e), nrow(r), sum(r$mr_keep)))
}

# ------------------------------------------------------------------ Steiger / directionality / heterogeneity / pleiotropy
sc_steiger <- function() {   # `size` = rows (pairs of K=20 SNPs)
  d <- gen_pairs(size / 20, seed = 21L, fixed_k = 20L)
  fun <- if (is_tsmr) function(x) q(function() steiger_filtering(x)) else function(x) q(function() fast_mr_steiger_filtering(x))
  tm <- time_calls(function() if (is_tsmr) par_apply(d, pkey(d), fun, 1L) else fun(d)); r <- as.data.frame(tm$val)
  save_std(data.frame(key = paste(r$id.exposure, r$id.outcome, r$SNP, sep = "|"), method = "all", rsq.exposure = r$rsq.exposure,
    rsq.outcome = r$rsq.outcome, steiger_pval = r$steiger_pval, flag_steiger_dir = as.numeric(r$steiger_dir), stringsAsFactors = FALSE))
  emit(tm$wall, tm$cpu, tm$n, sprintf("rows=%d;pairs=%d", nrow(d), size / 20))
}
sc_diag <- function(kind) {   # 1000 pairs, K ~ Poisson(20)
  d <- gen_pairs(size, seed = 31L)
  fun <- if (kind == "het") { if (is_tsmr) function(x) q(function() mr_heterogeneity(x)) else function(x) q(function() fast_mr_heterogeneity(d = x, methods = c("ivw", "egger"), threads = threads)) }
         else { if (is_tsmr) function(x) q(function() mr_pleiotropy_test(x)) else function(x) q(function() fast_mr_pleiotropy_test(x, threads = threads)) }
  if (!is_tsmr && kind == "het") fun <- function(x) q(function() fast_mr_heterogeneity(x, methods = c("ivw", "egger"), threads = threads))
  tm <- time_calls(function() if (is_tsmr) par_apply(d, pkey(d), fun, 1L) else fun(d)); r <- as.data.frame(tm$val)
  if (kind == "het") save_std(data.frame(key = pkey(r), method = r$method, Q = r$Q, Q_df = r$Q_df, Q_pval = r$Q_pval, stringsAsFactors = FALSE))
  else save_std(data.frame(key = pkey(r), method = "egger_intercept", egger_intercept = r$egger_intercept, se = r$se, pval = r$pval, stringsAsFactors = FALSE))
  emit(tm$wall, tm$cpu, tm$n, sprintf("pairs=%d;rows=%d", size, nrow(d)))
}

# ------------------------------------------------------------------ clumping (SNP ids drawn from the LD reference's own variants)
gen_clump <- function(E, seed = 41L) {
  set.seed(seed)
  bim <- data.table::fread(paste0(o$ref, ".bim"), header = FALSE, select = 1:4, col.names = c("chr", "SNP", "cm", "bp"))
  bad <- bim$SNP == "." | duplicated(bim$SNP)
  nloci <- 300L; pool <- which(!bad); centers <- sort(sample(seq_len(length(pool) - 200L), nloci))
  loci <- lapply(centers, function(c0) pool[c0 + 0:99])
  rows <- lapply(seq_len(E), function(e) {
    L <- sample(nloci, 15); idx <- unique(unlist(lapply(L, function(l) sample(loci[[l]], 20))))
    p <- 10^(-(7.5 + rexp(length(idx), 1/3))); p[sample(length(idx), 15)] <- 10^(-runif(15, 12, 40))   # all candidates are genome-wide-significant (instrument-like): plink1.9 default --clump-p2 0.01 is then moot
    data.frame(SNP = bim$SNP[idx], pval.exposure = p, id.exposure = paste0("E", e), exposure = paste0("E", e),
               chr_name = as.character(bim$chr[idx]), chrom_start = bim$bp[idx], stringsAsFactors = FALSE) })
  do.call(rbind, rows)
}
plink2_clump_one <- function(x) {   # per-exposure plink2 --clump
  fn <- tempfile(); write.table(data.frame(SNP = x$SNP, P = x$pval.exposure), fn, row.names = FALSE, quote = FALSE)
  system2(o$p2bin, c("--pfile", o$pref, "--clump", fn, "--clump-id-field", "SNP", "--clump-p-field", "P", "--clump-p1", "1", "--clump-p2", "1",
    "--clump-r2", "0.001", "--clump-kb", "10000", "--threads", "1", "--out", fn), stdout = FALSE, stderr = FALSE)
  cf <- paste0(fn, ".clumps"); ids <- if (file.exists(cf)) data.table::fread(cf, select = "ID")$ID else character()
  unlink(paste0(fn, "*")); x[x$SNP %in% ids, ]
}
sc_clump <- function() {
  d <- gen_clump(size)
  ex <- d$id.exposure
  if (is_tsmr) {
    fun <- function(x) do.call(rbind, lapply(split(x, x$id.exposure), function(y)
      q(function() ieugwasr::ld_clump(data.frame(rsid = y$SNP, pval = y$pval.exposure, id = y$id.exposure), clump_kb = 10000,
        clump_r2 = 0.001, clump_p = 1, bfile = o$ref, plink_bin = o$p1bin)) |> (\(r) y[y$SNP %in% r$rsid, ])()))
    f <- function() par_apply(d, ex, fun, 1L)
  } else if (arm$pkg == "plink2") {
    f <- function() do.call(rbind, lapply(split(d, d$id.exposure), plink2_clump_one))
  } else {
    f <- function() q(function() { a <- list(clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, plink2_bin = o$p2bin, threads = threads)
      if (!is.na(o$pref) && file.exists(paste0(o$pref, ".pgen"))) a$pfile <- o$pref else a$bfile <- o$ref
      do.call(fast_clump_data_graph, c(list(d), a))$data })
  }
  tm <- time_calls(f); r <- as.data.frame(tm$val)
  save_std(data.frame(key = paste(r$id.exposure, r$SNP, sep = "|"), method = "instruments", stringsAsFactors = FALSE))
  emit(tm$wall, tm$cpu, tm$n, sprintf("exposures=%d;rows_in=%d;retained=%d;ref=%s", size, nrow(d), nrow(r), basename(o$ref)))
}

switch(o$scen, single = sc_single(), many = sc_many(), harmonise = sc_harmonise(), harmonise_many = sc_harmonise_many(),
       steiger = sc_steiger(), heterogeneity = sc_diag("het"), pleiotropy = sc_diag("pleio"), clump = sc_clump(),
       stop("unknown scenario"))
