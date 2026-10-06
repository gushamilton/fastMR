# GP6-type component (E-mediation section 6 w_k approach) as an extra covariate trait:
# w = rank-1 direction of the 19:55 bin's SNP z-fingerprints (E's WB[["19:55"]], over the
# 2,940 outcomes); each panel SNP's "effect" on the component is s_i = sum_k z_ik w_k with
# outcomes within +/-1 Mb of the SNP masked (as E did). GP6 instruments: the 19:55 bin's
# trans-instrument SNPs, thinned greedily by |s| to >= 500 kb apart.
SET <- "trans_maf"
commandArgs <- function(...) SET
source(path.expand("~/agent-runs/fastmr-mvmr/bench/f06_defs.R"))  # F-platelet/code/06_mvmr.R up to "versions <-"
fm <- readRDS(file.path(A, "E-mediation/work/factor_model.rds"))
w <- fm$WB[["19:55"]]; stopifnot(length(w) == length(prot))
Z <- PB / PS
for (ch in unique(pos$chr)) {               # mask outcomes whose gene is within 1 Mb of the SNP
  r <- which(pos$chr == ch); cand <- which(meta$chr_name == as.character(ch))
  for (k in cand) { hit <- r[pos$pos[r] >= meta$gene_start[k] - 1e6 & pos$pos[r] <= meta$gene_end[k] + 1e6]; Z[hit, k] <- NA }
}
Z[!is.finite(Z)] <- 0
s <- as.vector(Z %*% w)
frame <- data.frame(SNP = rownames(PB), beta = s, se = 1)
inst_all <- unique(unlist(instA))
bin <- pos$chr == 19 & floor(pos$pos / 1e6) == 55 & pos$SNP %in% inst_all
cand <- pos$SNP[bin]; o <- order(-abs(s[match(cand, rownames(PB))])); cand <- cand[o]
kept <- character(); kp <- numeric()
for (c in cand) { p <- pos$pos[match(c, pos$SNP)]; if (all(abs(kp - p) >= 5e5)) { kept <- c(kept, c); kp <- c(kp, p) } }
cat("19:55 bin trans-instrument SNPs:", length(cand), "kept:", length(kept), "\n"); print(data.frame(SNP = kept, s = s[match(kept, rownames(PB))]))
inp <- readRDS(path.expand("~/agent-runs/fastmr-mvmr/data/real_inputs.rds"))
inp$gp6 <- list(frame = frame, instruments = kept, w = w)
saveRDS(inp, path.expand("~/agent-runs/fastmr-mvmr/data/real_inputs.rds"))
cat("saved\n")
