#!/usr/bin/env Rscript
# Self-contained fastMR end-to-end performance benchmark.
#
# Usage:
#   Rscript fastmr_perf_bench.R --scenario=<name> [--scale=small|full]
#     [--label=baseline] [--sha=<git sha>] [--repo=<path to fastMR checkout>]
#     [--lib=<R library containing the fastMR build to test>]
#     [--output=<csv path, appended>] [--replicate=1] [--threads=1]
#     [--profile=<Rprof output path>] [--clump_exposures=N]
#   (clump_mock runs the 'global' frontier only for N <= 100: at N = 300 it
#    exceeded 23 min / 6.3 GB RSS on an M5 Pro before being stopped.)
#
# Scenarios: grid_ivw, grid_ivw_parquet, grid_mixed, fast_mr_groups,
#   fast_mr_groups_mixed, sparse_ivw, sparse_ivw_mask, steiger_filtering,
#   clump_mock, compressed_pairwise, all
#
# One CSV row is emitted per scenario phase with host/CPU/R/git metadata,
# wall time, R-heap peak (gc), Linux VmHWM (process peak RSS), and two
# correctness checksums (exact md5 of the result, and md5 after signif(.,10)
# so cross-BLAS/cross-node comparisons are robust).  Inputs are generated with
# a fixed seed so baseline and optimised SHAs see identical data.

args <- commandArgs(trailingOnly = TRUE)
opt <- list(scenario = "all", scale = "small", label = "unlabelled", sha = NA_character_,
            repo = NA_character_, lib = NA_character_, output = NA_character_,
            replicate = "1", threads = "1", profile = NA_character_)
for (a in args) {
  kv <- regmatches(a, regexec("^--([^=]+)=(.*)$", a))[[1L]]
  if (length(kv) == 3L) opt[[kv[2L]]] <- kv[3L]
}
if (!is.na(opt$lib)) .libPaths(c(normalizePath(opt$lib), .libPaths()))
suppressPackageStartupMessages(library(fastMR))
threads <- as.integer(opt$threads)
full <- identical(opt$scale, "full")

sha <- opt$sha
if (is.na(sha) && !is.na(opt$repo)) {
  sha <- tryCatch(system2("git", c("-C", shQuote(opt$repo), "rev-parse", "--short=12", "HEAD"),
                          stdout = TRUE, stderr = FALSE)[1L], error = function(e) NA_character_)
}
cpu_model <- local({
  if (file.exists("/proc/cpuinfo")) {
    l <- grep("^model name", readLines("/proc/cpuinfo", warn = FALSE), value = TRUE)
    if (length(l)) return(trimws(sub("^[^:]*:", "", l[1L])))
  }
  out <- tryCatch(system2("sysctl", c("-n", "machdep.cpu.brand_string"), stdout = TRUE, stderr = FALSE),
                  error = function(e) character())
  if (length(out)) out[1L] else NA_character_
})
nproc <- parallel::detectCores()
vmhwm_mb <- function() {
  if (!file.exists("/proc/self/status")) return(NA_real_)
  l <- grep("^VmHWM:", readLines("/proc/self/status", warn = FALSE), value = TRUE)
  if (!length(l)) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", l)) / 1024
}
checksums <- function(x) {
  strip <- function(v) {
    if (is.data.frame(v)) { attr(v, "compressed_input") <- NULL; rownames(v) <- NULL }
    v
  }
  x <- strip(x)
  round_obj <- function(v) {
    if (is.list(v)) { v[] <- lapply(v, round_obj); return(v) }
    if (is.double(v)) return(signif(v, 10))
    v
  }
  f1 <- tempfile(); f2 <- tempfile()
  saveRDS(x, f1, compress = FALSE, version = 3)
  saveRDS(round_obj(x), f2, compress = FALSE, version = 3)
  out <- unname(tools::md5sum(c(f1, f2)))
  unlink(c(f1, f2))
  out
}
rows <- list()
emit <- function(scenario, phase, params, seconds, r_peak_mb, result) {
  cs <- if (is.null(result)) c(NA, NA) else checksums(result)
  row <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), label = opt$label, git_sha = sha,
    hostname = Sys.info()[["nodename"]], cpu_model = cpu_model, nproc = nproc,
    r_version = R.version.string, fastmr_version = as.character(utils::packageVersion("fastMR")),
    scenario = scenario, phase = phase, scale = opt$scale, replicate = opt$replicate,
    threads = threads, params = params, wall_s = seconds, r_heap_peak_mb = r_peak_mb,
    vmhwm_mb = vmhwm_mb(), checksum_exact = cs[1L], checksum_signif10 = cs[2L],
    stringsAsFactors = FALSE)
  rows[[length(rows) + 1L]] <<- row
  cat(sprintf("%-22s %-16s %10.3f s  heap %8.1f MB  %s\n", scenario, phase, seconds, r_peak_mb, cs[2L]))
}
timed <- function(expr) {
  invisible(gc(reset = TRUE, full = TRUE))
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  el <- proc.time()[["elapsed"]] - t0
  g <- gc(full = FALSE)
  list(value = value, seconds = el, peak_mb = sum(g[, ncol(g)]))
}

# ---------------------------------------------------------------- generators
gen_grid <- function(E, O, S, seed = 1L) {
  set.seed(seed)
  bx <- matrix(rnorm(E * S, 0.05, 0.02), E, S)
  by <- matrix(rnorm(O * S, 0, 0.02), O, S)
  sx <- matrix(runif(E * S, 0.005, 0.02), E, S)
  sy <- matrix(runif(O * S, 0.005, 0.02), O, S)
  rownames(bx) <- rownames(sx) <- paste0("exp", seq_len(E))
  rownames(by) <- rownames(sy) <- paste0("out", seq_len(O))
  list(bx = bx, by = by, sx = sx, sy = sy)
}
gen_long <- function(groups, snps_per_group, seed = 2L, ragged = TRUE, steiger = FALSE) {
  set.seed(seed)
  k <- if (ragged) pmax(3L, rpois(groups, snps_per_group)) else rep(snps_per_group, groups)
  n <- sum(k)
  g <- rep(seq_len(groups), k)
  n_exp <- ceiling(sqrt(groups))
  d <- data.frame(
    SNP = paste0("rs", sample.int(50L * n, n)),
    id.exposure = paste0("exp", (g - 1L) %% n_exp + 1L),
    id.outcome = paste0("out", (g - 1L) %/% n_exp + 1L),
    beta.exposure = rnorm(n, 0.05, 0.02), beta.outcome = rnorm(n, 0, 0.02),
    se.exposure = runif(n, 0.005, 0.02), se.outcome = runif(n, 0.005, 0.02),
    stringsAsFactors = FALSE)
  d$exposure <- d$id.exposure; d$outcome <- d$id.outcome
  if (steiger) {
    d$samplesize.exposure <- 30000; d$samplesize.outcome <- 300000
    d$pval.exposure <- 2 * pnorm(-abs(d$beta.exposure / d$se.exposure))
    d$pval.outcome <- 2 * pnorm(-abs(d$beta.outcome / d$se.outcome))
  }
  d
}
gen_csr <- function(E, O, S, per_exp, seed = 3L) {
  set.seed(seed)
  k <- pmax(1L, rpois(E, per_exp))
  col <- unlist(lapply(k, function(m) sort(sample.int(S, m)) - 1L))
  list(row_ptr = c(0L, cumsum(k)), col_index = as.integer(col),
       exposure_beta = rnorm(length(col), 0.05, 0.02),
       outcome_beta = matrix(rnorm(O * S, 0, 0.02), O, S,
                             dimnames = list(paste0("out", seq_len(O)), NULL)),
       outcome_se = matrix(runif(O * S, 0.005, 0.02), O, S),
       outcome_present = matrix(runif(O * S) < 0.97, O, S))
}

# ---------------------------------------------------------------- scenarios
sc <- list()
sc$grid_ivw <- function() {
  E <- if (full) 1000L else 300L; S <- 100L
  g <- gen_grid(E, E, S)
  r <- timed(fast_mr_grid(g$bx, g$by, g$sx, g$sy, methods = "ivw", nboot = 0, threads = threads))
  emit("grid_ivw", "total", sprintf("E=%d;O=%d;S=%d;methods=ivw", E, E, S), r$seconds, r$peak_mb, r$value)
  if (!exists("fastmr_grid_native", asNamespace("fastMR"))) return(invisible())
  n <- timed(get("fastmr_grid_native", asNamespace("fastMR"))(g$bx, g$by, g$sx, g$sy, "ivw", 0L, NULL, threads, 1, 20))
  emit("grid_ivw", "native_only", sprintf("E=%d;O=%d;S=%d", E, E, S), n$seconds, n$peak_mb, NULL)
}
sc$grid_ivw_parquet <- function() {
  if (!requireNamespace("arrow", quietly = TRUE)) return(invisible())
  E <- if (full) 1000L else 300L; S <- 100L
  g <- gen_grid(E, E, S)
  path <- tempfile(fileext = ".parquet")
  r <- timed(fast_mr_grid(g$bx, g$by, g$sx, g$sy, methods = "ivw", nboot = 0, threads = threads, output = path))
  emit("grid_ivw_parquet", "total", sprintf("E=%d;O=%d;S=%d;bytes=%d", E, E, S, file.size(path)),
       r$seconds, r$peak_mb, r$value)
  unlink(path)
}
sc$grid_mixed <- function() {
  E <- if (full) 200L else 60L; S <- 60L
  g <- gen_grid(E, E, S)
  m <- c("ivw", "egger", "weighted_median", "weighted_mode")
  r <- timed(fast_mr_grid(g$bx, g$by, g$sx, g$sy, methods = m, nboot = 100, seed = 1, threads = threads))
  emit("grid_mixed", "total", sprintf("E=%d;O=%d;S=%d;nboot=100;4methods", E, E, S), r$seconds, r$peak_mb, r$value)
}
sc$fast_mr_groups <- function() {
  G <- if (full) 20000L else 4000L
  d <- gen_long(G, 20L)
  r <- timed(fast_mr(d, methods = "ivw", nboot = 0, threads = threads))
  emit("fast_mr_groups", "total", sprintf("groups=%d;rows=%d;methods=ivw", G, nrow(d)), r$seconds, r$peak_mb, r$value)
}
sc$fast_mr_groups_mixed <- function() {
  G <- if (full) 2000L else 400L
  d <- gen_long(G, 20L)
  m <- c("ivw", "egger", "weighted_median", "weighted_mode")
  r <- timed(fast_mr(d, methods = m, nboot = 100, seed = 1, threads = threads))
  emit("fast_mr_groups_mixed", "total", sprintf("groups=%d;rows=%d;nboot=100;4methods", G, nrow(d)),
       r$seconds, r$peak_mb, r$value)
}
sc$sparse_ivw <- function() {
  E <- if (full) 20000L else 4000L; O <- if (full) 100L else 50L; S <- 20000L
  x <- gen_csr(E, O, S, 15)
  r <- timed(fast_mr_sparse_ivw(x$row_ptr, x$col_index, x$exposure_beta, x$outcome_beta,
                                x$outcome_se, x$outcome_present, threads = threads,
                                max_output_cells = 1e9, max_memory_mb = 1e5))
  emit("sparse_ivw", "total", sprintf("E=%d;O=%d;S=%d;nnz=%d", E, O, S, length(x$col_index)),
       r$seconds, r$peak_mb, r$value)
}
sc$sparse_ivw_mask <- function() {
  E <- if (full) 20000L else 4000L; O <- if (full) 100L else 50L; S <- 20000L
  x <- gen_csr(E, O, S, 15)
  set.seed(9)
  keep <- matrix(runif(O * length(x$col_index)) < 0.9, O, length(x$col_index))
  r <- timed(fast_mr_sparse_ivw(x$row_ptr, x$col_index, x$exposure_beta, x$outcome_beta,
                                x$outcome_se, x$outcome_present, threads = threads,
                                max_output_cells = 1e9, max_memory_mb = 1e5, pair_snp_keep = keep))
  emit("sparse_ivw_mask", "total", sprintf("E=%d;O=%d;nnz=%d;mask_MB=%.0f", E, O, length(x$col_index),
       object.size(keep) / 2^20), r$seconds, r$peak_mb, r$value)
}
sc$steiger_filtering <- function() {
  G <- if (full) 20000L else 4000L
  d <- gen_long(G, 10L, steiger = TRUE)
  r <- timed(fast_mr_steiger_filtering(d))
  emit("steiger_filtering", "per_pair_R", sprintf("pairs=%d;rows=%d", G, nrow(d)), r$seconds, r$peak_mb,
       r$value[order(r$value$id.exposure, r$value$id.outcome, r$value$SNP),
               c("SNP", "id.exposure", "id.outcome", "steiger_dir", "steiger_pval")])
  v <- timed({
    rx <- fast_mr_steiger_r2(d$beta.exposure, d$se.exposure, d$samplesize.exposure)
    ry <- fast_mr_steiger_r2(d$beta.outcome, d$se.outcome, d$samplesize.outcome)
    list(dir = rx$rsq > ry$rsq)
  })
  emit("steiger_filtering", "vector_r2_floor", sprintf("rows=%d", nrow(d)), v$seconds, v$peak_mb, NULL)
}
# PLINK is replaced by a deterministic synthetic LD oracle so the R frontier
# bookkeeping can be measured without a reference panel.  The oracle returns
# r2 = exp(-|dbp| / 25kb) for same-chromosome pairs within the window.
sc$clump_mock <- function() {
  # Cis-like candidate sets (cf. issue #8: ~274 candidates/exposure within
  # gene +/-300 kb, ~12 retained).  Exposures draw gene centres from a shared
  # pool so lead/target overlap across exposures is realistic.
  E <- if (!is.null(opt$clump_exposures)) as.integer(opt$clump_exposures) else if (full) 1000L else 100L
  set.seed(5)
  genes <- data.frame(chr = sample(as.character(1:22), 1500L, TRUE), centre = sample.int(2e8, 1500L) + 1e6)
  rows <- lapply(seq_len(E), function(e) {
    g <- genes[sample.int(nrow(genes), 1L), ]
    pos <- seq.int(g$centre - 300000L, g$centre + 300000L, by = 1000L)
    pos <- pos[runif(length(pos)) < 0.45]
    data.frame(SNP = sprintf("%s:%d:A:C", g$chr, pos), id.exposure = paste0("e", e),
               pval.exposure = runif(length(pos))^3 * 0.01, chr_name = g$chr, chrom_start = pos,
               stringsAsFactors = FALSE)
  })
  dat <- do.call(rbind, rows)
  u <- unique(dat$SNP)
  parts <- do.call(rbind, strsplit(u, ":", fixed = TRUE))
  snp <- u; chr <- parts[, 1L]; bp <- as.numeric(parts[, 2L])
  calls <- 0L; oracle_s <- 0
  oracle <- function(leads, targets, reference_args, plink2_bin, clump_kb, clump_r2, threads, workdir, round) {
    calls <<- calls + 1L
    t0 <- proc.time()[["elapsed"]]
    on.exit(oracle_s <<- oracle_s + proc.time()[["elapsed"]] - t0)
    li <- match(unique(leads), snp); ti <- match(unique(targets), snp)
    edges <- ld_edges(li, ti, clump_kb, clump_r2)
    data.frame(lead = snp[edges$l], target = snp[edges$t], stringsAsFactors = FALSE)
  }
  # All pairs (l in li, t in ti) on one chromosome with r2 >= thr (distance
  # bounded by the r2 decay), found via sorted positions, not a cross-join.
  ld_edges <- function(li, ti, clump_kb, thr) {
    dmax <- min(clump_kb * 1000, -20000 * log(thr))
    key <- paste(chr[ti]); o <- order(key, bp[ti]); ti <- ti[o]
    out_l <- list(); out_t <- list()
    for (cc in unique(chr[li])) {
      tt <- ti[chr[ti] == cc]; if (!length(tt)) next
      ll <- li[chr[li] == cc]
      lo <- findInterval(bp[ll] - dmax, bp[tt], left.open = TRUE) + 1L
      hi <- findInterval(bp[ll] + dmax, bp[tt])
      n <- pmax(0L, hi - lo + 1L)
      if (!sum(n)) next
      L <- rep(ll, n); T <- tt[unlist(mapply(function(a, b) if (b >= a) a:b else integer(), lo, hi, SIMPLIFY = FALSE))]
      k <- exp(-abs(bp[L] - bp[T]) / 20000) >= thr
      out_l[[cc]] <- L[k]; out_t[[cc]] <- T[k]
    }
    list(l = unlist(out_l, use.names = FALSE), t = unlist(out_t, use.names = FALSE))
  }
  # Prototype of the proposed exact strategy: one LD-graph query over the
  # candidate union, then per-exposure greedy clumping against the in-memory
  # graph.  Must reproduce the PLINK-frontier instruments exactly.
  graph_clump <- function(dat, clump_kb, clump_r2, clump_p1) {
    ui <- match(dat$SNP, snp)
    all <- unique(ui)
    e <- ld_edges(all, all, clump_kb, clump_r2)
    adj <- split(c(e$t, e$l), factor(c(e$l, e$t), levels = seq_along(snp)))
    p <- dat$pval.exposure
    res <- lapply(split(seq_len(nrow(dat)), dat$id.exposure), function(ii) {
      ii <- ii[is.finite(p[ii]) & p[ii] <= clump_p1]
      ii <- ii[order(p[ii], dat$SNP[ii], method = "radix")]
      v <- ui[ii]; dead <- logical(length(v)); keep <- integer()
      pos <- match(seq_along(snp), v)  # local index lookup
      for (k in seq_along(v)) {
        if (dead[k]) next
        keep <- c(keep, k)
        nb <- pos[adj[[v[k]]]]; nb <- nb[!is.na(nb)]
        dead[nb] <- TRUE
      }
      dat$SNP[ii[keep]]
    })
    res
  }
  # One all-pairs call per chromosome (the 'graph' partition's oracle).
  graph_oracle <- function(snps, reference_args, plink2_bin, clump_kb, clump_r2, threads, workdir, tag, ...) {
    calls <<- calls + 1L
    t0 <- proc.time()[["elapsed"]]
    on.exit(oracle_s <<- oracle_s + proc.time()[["elapsed"]] - t0)
    i <- match(unique(snps), snp)
    edges <- ld_edges(i, i, clump_kb, clump_r2)
    data.frame(lead = snp[edges$l], target = snp[edges$t], stringsAsFactors = FALSE)
  }
  ns <- asNamespace("fastMR")
  if (!exists("fastmr_clump_run_frontier", ns, inherits = FALSE)) stop("clump_mock needs fastmr_clump_run_frontier")
  swap <- function(name, fun) {
    orig <- get(name, ns); unlockBinding(name, ns); assign(name, fun, ns)
    function() { assign(name, orig, ns); lockBinding(name, ns) }
  }
  restore <- list(swap("fastmr_clump_run_frontier", oracle))
  if (exists("fastmr_clump_run_graph", ns, inherits = FALSE)) restore <- c(restore, swap("fastmr_clump_run_graph", graph_oracle))
  on.exit(for (f in restore) f())
  strategies <- c(if (E <= 100L) "global", "lead_row", if (exists("fast_clump_data_graph")) "graph")
  for (strategy in strategies) {
    fun <- switch(strategy, global = fast_clump_data_batched, lead_row = fast_clump_data_lead_rows,
                  graph = fast_clump_data_graph)
    calls <- 0L; oracle_s <- 0
    r <- timed(fun(dat, clump_kb = 10000, clump_r2 = 0.01, clump_p1 = 0.01, plink2_bin = "/bin/true", bfile = "mock"))
    emit("clump_mock", strategy, sprintf("E=%d;rows=%d;retained=%d;rounds=%d;oracle_calls=%d;oracle_s=%.2f;logical_pairs=%.0f",
         E, nrow(dat), nrow(r$value$data), r$value$diagnostics$rounds, calls, oracle_s,
         r$value$diagnostics$logical_pairs),
         r$seconds, r$peak_mb, lapply(r$value$instruments, sort))
  }
  g <- timed(graph_clump(dat, 10000, 0.01, 0.01))
  emit("clump_mock", "graph_prototype_R", sprintf("E=%d;rows=%d;retained=%d", E, nrow(dat), sum(lengths(g$value))),
       g$seconds, g$peak_mb, lapply(g$value, sort))
}
sc$compressed_pairwise <- function() {
  if (!requireNamespace("CompreSSoR", quietly = TRUE)) return(invisible())
  E <- if (full) 100L else 30L; O <- if (full) 100L else 30L; V <- 4000L
  dir <- tempfile("cmp_"); dir.create(dir)
  set.seed(7)
  make <- function(name, mult) {
    df <- data.frame(chromosome = "1", base_pair_location = seq.int(100001L, length.out = V),
                     reference_allele = "A", alternate_allele = "C", effect_allele = "C", other_allele = "A",
                     beta = rnorm(V, 0, 0.02 * mult), standard_error = runif(V, 0.005, 0.02),
                     effect_allele_frequency = runif(V, 0.05, 0.95), stringsAsFactors = FALSE)
    p <- file.path(dir, name)
    invisible(CompreSSoR::compress_sumstats(df, p, overwrite = TRUE))
    p
  }
  ef <- vapply(seq_len(E), function(i) make(paste0("e", i), 2), ""); names(ef) <- paste0("e", seq_len(E))
  of <- vapply(seq_len(O), function(i) make(paste0("o", i), 1), ""); names(of) <- paste0("o", seq_len(O))
  keys <- CompreSSoR::compressor_variant_key("1", seq.int(100001L, length.out = V), "A", "C")
  inst <- lapply(seq_len(E), function(i) sample(keys, 25L)); names(inst) <- names(ef)
  r <- timed(fast_mr_compressed(ef, of, inst, methods = "ivw", nboot = 0, threads = threads))
  tm <- attr(r$value, "compressed_input")$timing
  emit("compressed_pairwise", "total", sprintf("E=%d;O=%d;io_s=%.3f;estimator_s=%.3f;path=%s", E, O,
       tm$io_seconds, tm$estimator_seconds, attr(r$value, "compressed_input")$estimator_path),
       r$seconds, r$peak_mb, r$value)
  unlink(dir, recursive = TRUE)
}

todo <- if (identical(opt$scenario, "all")) names(sc) else strsplit(opt$scenario, ",")[[1L]]
if (!is.na(opt$profile)) Rprof(opt$profile, interval = 0.005, memory.profiling = TRUE)
for (s in todo) sc[[s]]()
if (!is.na(opt$profile)) Rprof(NULL)
out <- do.call(rbind, rows)
if (!is.na(opt$output)) {
  dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(out, opt$output, sep = ",", row.names = FALSE, qmethod = "double",
                     col.names = !file.exists(opt$output), append = file.exists(opt$output))
}
