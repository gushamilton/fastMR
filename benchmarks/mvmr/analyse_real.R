# Mac-side: compare the BP batched-MVMR runs with F-platelet's M1/R1/R2 matrices and
# summarise the residualised "how close to null" progression (F's 09 metrics).
# Usage: Rscript analyse_real.R <results_dir> <inputs.rds> <out_prefix>
suppressMessages(library(data.table))
args <- commandArgs(TRUE); RD <- args[[1]]; inp <- readRDS(args[[2]]); OUT <- args[[3]]
Fw <- path.expand("~/agent-runs/cistrans/F-platelet/work")
A <- path.expand("~/agent-runs/cistrans")
gene <- inp$gene; prot <- inp$excl$id.outcome
mask_same <- function(m) { ex <- rownames(m); same <- outer(gene[match(ex, prot)], gene, "=="); m[same] <- NA; m }
ours <- function(run) { r <- readRDS(file.path(RD, paste0(run, ".rds"))); r$b <- mask_same(r$b); r$se <- mask_same(r$se); r }
fref <- function(v) readRDS(file.path(Fw, paste0("mvmr_trans_maf_", v, ".rds")))
cmp <- function(o, f, label) {
  ex <- intersect(rownames(o$b), rownames(f$b)); ob <- o$b[ex, prot]; fb <- f$b[ex, prot]; os <- o$se[ex, prot]; fs <- f$se[ex, prot]
  both <- is.finite(ob) & is.finite(fb); only_o <- is.finite(ob) & !is.finite(fb); only_f <- !is.finite(ob) & is.finite(fb)
  rb <- abs(ob - fb)[both] / pmax(abs(fb[both]), 1e-12); rs <- abs(os - fs)[both] / fs[both]
  data.table(comparison = label, exposures_ours = nrow(o$b), exposures_F = nrow(f$b), common_exposures = length(ex),
             pairs_both = sum(both), only_ours = sum(only_o), only_F = sum(only_f),
             max_abs_db = max(abs(ob - fb)[both]), median_rel_db = median(rb), q99_rel_db = quantile(rb, .99),
             max_rel_db = max(rb), frac_rel_db_lt_1e10 = mean(rb < 1e-10),
             max_rel_dse = max(rs), frac_rel_dse_lt_1e10 = mean(rs < 1e-10), cor_b = cor(ob[both], fb[both]))
}
res <- list()
M1 <- ours("M1"); R1 <- ours("R1")
res[[1]] <- cmp(M1, fref("M1"), "M1 (F's union sets, floored se) vs F M1")
res[[2]] <- cmp(R1, fref("R1"), "R1 residualised PLT vs F R1")
if (file.exists(file.path(RD, "R2.rds"))) res[[3]] <- cmp(ours("R2"), fref("R2"), "R2 residualised PLT+MPV vs F R2")
if (file.exists(file.path(RD, "M1_greedy.rds"))) {
  G <- ours("M1_greedy"); res[[4]] <- cmp(G, fref("M1"), "M1 package LD-greedy (max |z|) vs F M1")
  same_set <- mapply(function(a, b) setequal(a, b), G$instruments[names(inp$m1)], inp$m1)
  cat("package greedy union sets identical to F's for", sum(same_set), "of", length(same_set), "exposures\n")
}
if (file.exists(file.path(RD, "M1_tsmr.rds"))) {
  T <- ours("M1_tsmr"); ratio <- (T$se / M1$se)[is.finite(T$se) & is.finite(M1$se)]
  cat("M1 TwoSampleMR (unfloored) se / floored se: median", median(ratio), "frac < 1:", mean(ratio < 1 - 1e-12), "\n")
}
if (file.exists(file.path(RD, "M1_shared.rds"))) {
  S <- ours("M1_shared"); res[[5]] <- cmp(S, M1, "M1 shared-weight path (tol 0.25) vs exact M1")
  cat("shared fits:", S$shared_fits, "of", sum(is.finite(M1$b)), "\n")
}
cmpT <- rbindlist(res, fill = TRUE); print(cmpT); fwrite(cmpT, paste0(OUT, "_reproduction.tsv"), sep = "\t")
# conditional F: package (MVMR::strength_mvmr delta, cor 0, L-(p-1)) vs F (optimised delta)
fM1 <- fref("M1"); cf <- merge(data.table(e = rownames(M1$conditional_F), ours = M1$conditional_F[, "exposure"]),
                             fM1$info[, .(e, F = condF)], by = "e")
cat("conditional F (exposure): median ours", median(cf$ours), "F", median(cf$F), "| ours>=F:", mean(cf$ours >= cf$F - 1e-9),
    "| n<10 ours", sum(cf$ours < 10), "F", sum(cf$F < 10), "| cor", cor(cf$ours, cf$F), "\n")
# ---- progression metrics (F's 09_mvmr_summary.R definitions) ----
L <- fread(file.path(A, "F-platelet/out/factor_loadings_oriented.tsv")); stopifnot(identical(L$id, prot))
v1 <- L$C_v1 / sqrt(sum(L$C_v1^2))
hb <- fread(file.path(A, "C-network/out/hotspot_bins_gwas_labelled.tsv"))
sf <- fread(file.path(A, "D-distributions/out/signflip_by_nloci.tsv"))
binof <- function(s) paste0(sub(":.*", "", s), ":", floor(as.numeric(sub("^[^:]+:([0-9]+):.*", "\\1", s)) / 1e6))
plt_exp <- names(inp$own)[vapply(inp$own, function(s) any(grepl("latelet", hb$gwas_label[match(binof(s), hb$bin)])), logical(1))]
nloc <- lengths(inp$own)
p_null <- function(n) { br <- c(1, 2, 3, 5, 10, 20, 50, 1e4); i <- findInterval(n, br, left.open = TRUE)
  out <- sf$p_z4_null[pmax(i, 1)]; out[n <= 1] <- 2 * pnorm(-4); out }
versions <- list(M0 = fref("M0"), R1_F = fref("R1"), R2_F = fref("R2"), R1 = R1)
for (v in c("R2", "R3")) if (file.exists(file.path(RD, paste0(v, ".rds")))) versions[[v]] <- ours(v)
versions$M0$b <- mask_same(versions$M0$b); versions$M0$se <- mask_same(versions$M0$se)
ex <- Reduce(intersect, lapply(versions, function(r) rownames(r$b)))
fin <- Reduce(`&`, lapply(versions, function(r) is.finite(r$b[ex, prot]) & is.finite(r$se[ex, prot]) & r$se[ex, prot] > 0))
pe <- ex %in% plt_exp
summ <- rbindlist(lapply(names(versions), function(v) {
  r <- versions[[v]]; Z <- r$b[ex, prot] / r$se[ex, prot]; Z[!fin] <- NA
  Zc <- pmin(pmax(Z, -15), 15); Zc[is.na(Zc)] <- 0; proj <- as.vector(Zc %*% v1)
  Zp <- Z[pe, , drop = FALSE]; nl <- matrix(nloc[ex[pe]], nrow(Zp), ncol(Zp)); ok <- is.finite(Zp)
  obs <- mean(abs(Zp[ok]) > 4); nul <- mean(p_null(nl[ok]))
  tm <- if (!is.null(r$timing)) r$timing else NULL
  data.table(version = v, n_exp = length(ex), n_pairs = sum(fin), share_Z2_on_v1 = sum(proj^2) / sum(Zc^2),
             pct_pos_absZ4 = 100 * mean(Z[abs(Z) > 4] > 0, na.rm = TRUE), mean_Z2 = mean(Z^2, na.rm = TRUE),
             plt_exposures = sum(pe), plt_pairs = sum(ok), plt_frac_absZ4 = obs, plt_null_frac_absZ4 = nul,
             plt_excess_over_null = obs - nul, plt_obs_over_null = obs / nul,
             io_s = if (is.null(tm)) NA else tm$io_seconds, design_s = if (is.null(tm)) NA else tm$design_seconds,
             estimator_s = if (is.null(tm)) NA else tm$estimator_seconds, total_s = if (is.null(tm)) NA else tm$total_seconds)
}))
print(summ[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)])
fwrite(summ, paste0(OUT, "_progression.tsv"), sep = "\t")
tim <- rbindlist(lapply(list.files(RD, pattern = "^[MR].*\\.rds$", full.names = TRUE), function(f) {
  r <- readRDS(f); data.table(run = sub("\\.rds$", "", basename(f)), threads = r$threads, cpu = r$cpu,
    io_s = r$timing$io_seconds, design_s = r$timing$design_seconds, estimator_s = r$timing$estimator_seconds,
    total_s = r$timing$total_seconds, union_snps = r$counts$union_snps, GB_read = r$timing$source_bytes_read / 1e9)
}))
print(tim); fwrite(tim, paste0(OUT, "_timing.tsv"), sep = "\t")
