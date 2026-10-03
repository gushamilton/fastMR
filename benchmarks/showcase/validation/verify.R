a <- commandArgs(TRUE); lib <- a[1]; bref <- a[2]; pref <- a[3]; r2 <- as.numeric(a[4]); kb <- as.numeric(a[5]); E <- as.integer(a[6]); seed <- as.integer(a[7]); mode <- if (length(a) >= 8) a[8] else "pfile"
.libPaths(c(lib, "/user/work/fh6520/regress-audit/lib-fm-final", .libPaths()))
suppressPackageStartupMessages({library(data.table); library(fastMR)})
p2b <- Sys.which("plink2"); p1b <- Sys.which("plink")
gen_clump <- function(E, seed) {
  set.seed(seed)
  bim <- fread(paste0(bref, ".bim"), header = FALSE, select = 1:4, col.names = c("chr", "SNP", "cm", "bp"))
  bad <- bim$SNP == "." | duplicated(bim$SNP)
  nloci <- 300L; pool <- which(!bad); centers <- sort(sample(seq_len(length(pool) - 200L), nloci))
  loci <- lapply(centers, function(c0) pool[c0 + 0:99])
  rows <- lapply(seq_len(E), function(e) {
    L <- sample(nloci, 15); idx <- unique(unlist(lapply(L, function(l) sample(loci[[l]], 20))))
    p <- 10^(-(7.5 + rexp(length(idx), 1/3))); p[sample(length(idx), 15)] <- 10^(-runif(15, 12, 40))
    data.frame(SNP = bim$SNP[idx], pval.exposure = p, id.exposure = paste0("E", e),
               chr_name = as.character(bim$chr[idx]), chrom_start = bim$bp[idx], stringsAsFactors = FALSE) })
  do.call(rbind, rows)
}
d <- gen_clump(E, seed)
ra <- if (mode == "pfile") list(pfile = pref) else list(bfile = bref)
fc <- function(f) do.call(f, c(list(d, clump_kb = kb, clump_r2 = r2, plink2_bin = p2b, threads = 4), ra))$data
wd <- tempfile("v_", tmpdir = getwd()); dir.create(wd)
key <- function(x) sort(paste(x$id.exposure, x$SNP))
per <- function(fun) do.call(rbind, lapply(split(d, d$id.exposure), function(x) { fn <- file.path(wd, x$id.exposure[1])
  write.table(data.frame(SNP = x$SNP, P = x$pval.exposure), fn, row.names = FALSE, quote = FALSE); x[x$SNP %in% fun(fn), ] }))
res <- list(
  plink19 = per(function(fn) { system2(p1b, c("--bfile", bref, "--clump", fn, "--clump-p1 1 --clump-p2 1 --clump-r2", r2, "--clump-kb", kb, "--out", fn), stdout = FALSE, stderr = FALSE)
    cf <- paste0(fn, ".clumped"); if (file.exists(cf)) fread(cf)$SNP else character() }),
  plink2 = per(function(fn) { system2(p2b, c("--pfile", pref, "--clump", fn, "--clump-p1 1 --clump-p2 1 --clump-r2", r2, "--clump-kb", kb, "--out", fn), stdout = FALSE, stderr = FALSE)
    cf <- paste0(fn, ".clumps"); if (file.exists(cf)) fread(cf)$ID else character() }),
  graph = fc(fast_clump_data_graph),
  lead_row = fc(fast_clump_data_lead_rows),
  global = fc(fast_clump_data_batched))
k <- lapply(res, key)
cmp <- function(z, r) if (identical(z, r)) "=" else sprintf("d%d", length(union(setdiff(z, r), setdiff(r, z))))
cat(sprintf("ref=%s fm=%s r2=%g kb=%g E=%d seed=%d rows=%d | %s\n", basename(bref), mode, r2, kb, E, seed, nrow(d),
  paste(sprintf("%s=%d[vs19:%s vs2:%s]", names(k), lengths(k), vapply(k, cmp, "", k$plink19), vapply(k, cmp, "", k$plink2)), collapse = " ")))
unlink(wd, recursive = TRUE)
