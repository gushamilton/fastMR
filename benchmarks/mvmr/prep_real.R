# Mac-side preparation of the real-data MVMR benchmark inputs (F-platelet).
# Sources F's 06_mvmr.R definitions (panel, LD neighbours, greedy union) unchanged.
SET <- "trans_maf"
commandArgs <- function(...) SET
source(path.expand("~/agent-runs/fastmr-mvmr/bench/f06_defs.R"))  # F-platelet/code/06_mvmr.R up to "versions <-"
OUT <- path.expand("~/agent-runs/fastmr-mvmr/data")
exps <- names(instA)[names(instA) %in% prot]
cat("exposures", length(exps), "\n")
# F's own-instrument sets (as fit_one: in panel, finite exposure values)
own <- lapply(exps, function(e) { ec <- match(e, prot); sA <- intersect(instA[[e]], rownames(PB)); sA[is.finite(PB[sA, ec]) & is.finite(PS[sA, ec])] })
names(own) <- exps
# F's M1 union sets (A + PLT, greedy by smallest p, then finite rows)
m1 <- lapply(exps, function(e) {
  ec <- match(e, prot); sA <- own[[e]]
  snps <- union_snps(sA, ec, "PLT")
  X <- cbind(PB[snps, ec], tb("PLT", snps)); Sx <- cbind(PS[snps, ec], tse("PLT", snps))
  snps[rowSums(!is.finite(X)) == 0 & rowSums(!is.finite(Sx)) == 0]
})
names(m1) <- exps
# Trait instrument sets used by F's trait_fit (Gamma)
gsets <- lapply(c(PLT = "PLT", MPV = "MPV"), function(t) {
  s <- union_snps(character(0), NA, t); s[is.finite(tb(t, s))] })
gPM <- { s <- union_snps(character(0), NA, c("PLT", "MPV")); s[is.finite(tb("PLT", s)) & is.finite(tb("MPV", s))] }
cat("trait sets", lengths(gsets), length(gPM), "\n")
covs <- lapply(c(PLT = "PLT", MPV = "MPV"), function(t) {
  d <- data.frame(SNP = lk$SNP, beta = lk[[paste0(t, "_b")]], se = lk[[paste0(t, "_se")]])
  d[is.finite(d$beta) & is.finite(d$se), ] })
excl <- data.frame(id.outcome = prot, chromosome = meta$chr_name,
                   start = meta$gene_start - 1e6, end = meta$gene_end + 1e6)
saveRDS(list(own = own, m1 = m1, gsets = gsets, gPM = gPM, covs = covs, excl = excl,
             gene = gene_of, ld = as.data.frame(ld[a < b]), trait_inst = tinst[trait %in% c("PLT", "MPV")]),
        file.path(OUT, "real_inputs.rds"))
cat("saved\n")
