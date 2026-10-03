#!/usr/bin/env Rscript
# Deterministic, seeded simulated-GWAS generator (plan revision 15).
#
# Generator in one paragraph.  Null background: one `plink2 --glm` pass over the EUR
# MAF>=0.01 panel (~633 samples) with 20 independent N(0,1) phenotype columns gives
# LD-correlated null test statistics; each T is mapped to an exact N(0,1) Z per variant
# (qnorm(pt(T))).  Signal: each exposure gets K~U{20,100} causal SNPs (MAF>=0.05) with
# standardised effects b_j (lambda_j = sqrt(N)*b_j ~ +-U(4.5,12)); the marginal Z within
# +-1 Mb is Z_null + sqrt(N) * sum_j r_ij b_j with r from the same genotypes
# (plink2 --r-unphased).  Each outcome shares every exposure's causal SNPs with a true
# theta (30% exactly 0, else +-U(0.2,0.8)), outcome effect = theta*b_ej plus mild balanced
# pleiotropy (20% of the shared SNPs get an extra N(0,(0.2|b|)^2) direct effect).
# SE = 1/sqrt(2N f(1-f)) * exp(0.02*N(0,1)); beta = Z*SE; p = 2*Phi(-|Z|); forward
# strand, effect allele = ALT, other = REF.  Values are rounded to 4 significant digits
# (as in the FinnGen R1 file) before p/LP are written.
#
# Env: LIB_CS LIB_FASTMR LIB_TSMR SEED NTHREADS REF TRAITS  (see common.R)
args <- commandArgs(TRUE)
src_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(src_dir, "common.R"))
setup_libs(c("cs"))
suppressMessages(library(data.table))
SEED <- as.integer(Sys.getenv("SEED", "20261002")); NT <- as.integer(Sys.getenv("NTHREADS", "16"))
setDTthreads(NT)
NEXP <- 10L; NOUT <- 10L; NT_TOT <- NEXP + NOUT
TMP <- file.path(TRAITS, "tmp"); dir.create(TMP, recursive = TRUE, showWarnings = FALSE)
log <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")
set.seed(SEED)

# ---- design (all random draws up front: deterministic) ---------------------------------
pv <- fread(cmd = paste0("grep -v '^##' ", REF, ".pvar"), select = 1:5, nThread = NT)
setnames(pv, c("chr", "pos", "id", "ref", "alt")); nv <- nrow(pv)
psam <- fread(paste0(REF, ".psam")); ns <- nrow(psam)
log("variants", nv, "samples", ns)
Ns <- c(round(runif(NEXP, 50e3, 400e3), -3), round(runif(NOUT, 100e3, 500e3), -3))
K <- sample(20:100, NEXP, replace = TRUE)
trait_id <- c(sprintf("exp%02d", 1:NEXP), sprintf("out%02d", 1:NOUT))
theta <- matrix(0, NEXP, NOUT, dimnames = list(trait_id[1:NEXP], trait_id[NEXP + 1:NOUT]))
nz <- matrix(runif(NEXP * NOUT) > 0.3, NEXP, NOUT)
theta[nz] <- sample(c(-1, 1), sum(nz), TRUE) * runif(sum(nz), 0.2, 0.8)
pheno <- matrix(rnorm(ns * NT_TOT), ns, NT_TOT)

# ---- stage 1: plink2 passes (cached) ---------------------------------------------------
glm_done <- file.path(TMP, "glm.done")
if (!file.exists(glm_done)) {
  ph <- data.table(IID = psam[[1]]); ph <- cbind(ph, as.data.table(pheno))
  setnames(ph, c("#IID", trait_id)); fwrite(ph, file.path(TMP, "pheno.txt"), sep = "\t")
  log("plink2 --glm")
  rc <- system2("plink2", c("--pfile", REF, "--pheno", file.path(TMP, "pheno.txt"), "--glm", "allow-no-covars",
                "cols=chrom,pos,tz", "--threads", NT, "--memory", 60000, "--out", file.path(TMP, "glm")))
  stopifnot(rc == 0)
  rc <- system2("plink2", c("--pfile", REF, "--freq", "--threads", NT, "--out", file.path(TMP, "freq")))
  stopifnot(rc == 0); file.create(glm_done)
}
fq <- fread(file.path(TMP, "freq.afreq")); stopifnot(nrow(fq) == nv, all(fq$ID == pv$id))
f <- fq$ALT_FREQS

# ---- causal SNPs, effects --------------------------------------------------------------
eligible <- which(pmin(f, 1 - f) >= 0.05 & pv$chr %in% 1:22)
causal <- lapply(1:NEXP, function(e) sort(sample(eligible, K[e])))
used <- unique(unlist(causal)); stopifnot(!anyDuplicated(unlist(causal)))
bstd <- lapply(1:NEXP, function(e) sample(c(-1, 1), K[e], TRUE) * runif(K[e], 4.5, 12) / sqrt(Ns[e]))
ld_done <- file.path(TMP, "ld.done")
if (!file.exists(ld_done)) {
  writeLines(pv$id[used], file.path(TMP, "causal.txt"))
  log("plink2 --r-unphased around", length(used), "causal SNPs")
  rc <- system2("plink2", c("--pfile", REF, "--ld-snp-list", file.path(TMP, "causal.txt"), "--r-unphased",
        "cols=id,maj,nonmaj", "--ld-window", 1e9, "--ld-window-kb", 1000, "--ld-window-r2", 0,
        "--threads", NT, "--out", file.path(TMP, "ld")))
  stopifnot(rc == 0); file.create(ld_done)
}
ld <- fread(file.path(TMP, "ld.vcor"))
cat("LD columns:", names(ld), "\n")
# plink2 reports r for the NON-MAJOR alleles; flip sign per variant so r refers to ALT dosages
nmA <- grep("^NONMAJ_A", names(ld), value = TRUE); nmB <- grep("^NONMAJ_B", names(ld), value = TRUE)
stopifnot(length(nmA) == 1, length(nmB) == 1)
setnames(ld, "#ID_A", "ID_A")
ld <- ld[, .(a = match(ID_A, pv$id), b = match(ID_B, pv$id), r = UNPHASED_R, nmA = get(nmA), nmB = get(nmB))]
ld[, r := r * ifelse(nmA == pv$alt[a], 1, -1) * ifelse(nmB == pv$alt[b], 1, -1)][, c("nmA", "nmB") := NULL]
isc <- logical(nv); isc[used] <- TRUE
rev <- ld[isc[b]]; setnames(rev, c("a", "b", "r"), c("b", "a", "r"))
ld <- unique(rbind(ld[isc[a]], rev[isc[a]], data.table(a = used, b = used, r = 1)), by = c("a", "b"))
ld <- ld[abs(pv$pos[a] - pv$pos[b]) <= 1e6 & pv$chr[a] == pv$chr[b]]
log("LD pairs", nrow(ld))

# standardised effect vector per trait over all variants (sparse over `used`)
beta_vec <- vector("list", NT_TOT)
for (e in 1:NEXP) { v <- numeric(nv); v[causal[[e]]] <- bstd[[e]]; beta_vec[[e]] <- v }
truth <- list()
for (o in 1:NOUT) {
  v <- numeric(nv)
  for (e in 1:NEXP) {
    be <- bstd[[e]]; eff <- theta[e, o] * be
    pl <- runif(K[e]) < 0.2
    eff[pl] <- eff[pl] + rnorm(sum(pl), 0, 0.2 * abs(be[pl]))
    v[causal[[e]]] <- eff
    truth[[length(truth) + 1]] <- data.table(exposure = trait_id[e], outcome = trait_id[NEXP + o],
                                             theta = theta[e, o], n_causal = K[e], n_pleio = sum(pl))
  }
  beta_vec[[NEXP + o]] <- v
}

# design tables (written by every run; identical content)
fwrite(rbindlist(truth), file.path(TRAITS, "truth_theta.csv"))
fwrite(data.table(trait = trait_id, role = rep(c("exposure", "outcome"), each = 10), N = Ns,
                  K = c(K, rep(NA, NOUT))), file.path(TRAITS, "traits.csv"))
saveRDS(list(causal_ids = lapply(causal, function(i) pv$id[i]), seed = SEED), file.path(TRAITS, "causal.rds"))
if (nzchar(Sys.getenv("STAGE1_ONLY"))) { log("stage 1 done"); quit(save = "no") }
# traits to write: TRAIT_SET="1,2,3" or TASK_CHUNK=<n traits per array task> with SLURM_ARRAY_TASK_ID
todo <- seq_len(NT_TOT)
if (nzchar(Sys.getenv("TRAIT_SET"))) todo <- as.integer(strsplit(Sys.getenv("TRAIT_SET"), ",")[[1]])
if (nzchar(Sys.getenv("TASK_CHUNK"))) { ch <- as.integer(Sys.getenv("TASK_CHUNK")); k <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  todo <- intersect(todo, ((k - 1) * ch + 1):(k * ch)) }
dir.create(file.path(TRAITS, "conv"), showWarnings = FALSE)
# ---- per-trait output ------------------------------------------------------------------
dir.create(file.path(TRAITS, "tsv"), showWarnings = FALSE); dir.create(file.path(TRAITS, "vcf"), showWarnings = FALSE)
dir.create(file.path(TRAITS, "cpr"), showWarnings = FALSE)
conv <- list()
for (t in todo) {
  tid <- trait_id[t]; N <- Ns[t]
  if (file.exists(file.path(TRAITS, "conv", paste0(tid, ".csv")))) { log(tid, "already done"); next }
  tt <- fread(file.path(TMP, sprintf("glm.%s.glm.linear", tid)), select = c("ID", "T_STAT"), nThread = NT)
  tstat <- tt$T_STAT[match(pv$id, tt$ID)]; rm(tt)
  z <- qnorm(pt(tstat, df = ns - 2, log.p = TRUE), log.p = TRUE)  # exact N(0,1) per variant
  z[!is.finite(z)] <- 0
  b <- beta_vec[[t]]; nzc <- which(b != 0)
  sig <- numeric(nv)
  # signal in Z units: sqrt(N) * sum_j r_ij b_j
  l <- ld[a %in% nzc]
  add <- rowsum(l$r * b[l$a], l$b, reorder = FALSE)
  sig[as.integer(rownames(add))] <- add[, 1]
  z <- z + sqrt(N) * sig
  se <- 1 / sqrt(2 * N * f * (1 - f)) * exp(0.02 * rnorm(nv))
  bb <- signif(z * se, 4); se <- signif(se, 4)
  zr <- bb / se
  p <- signif(2 * pnorm(-abs(z)), 4)
  d <- data.table(chromosome = pv$chr, base_pair_location = pv$pos, effect_allele = pv$alt, other_allele = pv$ref,
                  beta = bb, standard_error = se, effect_allele_frequency = signif(f, 4), p_value = p,
                  rsid = pv$id, n = N)
  tsv <- file.path(TRAITS, "tsv", paste0(tid, ".tsv.gz"))
  tw <- system.time({ fwrite(d, file.path(TMP, paste0(tid, ".tsv")), sep = "\t", nThread = NT)
    system2("gzip", c("-6", "-f", file.path(TMP, paste0(tid, ".tsv")))); file.rename(file.path(TMP, paste0(tid, ".tsv.gz")), tsv) })[["elapsed"]]
  rm(d, z, zr, se, bb, p, sig)
  cc <- convert_tsv_to_cpr(tsv, file.path(TRAITS, "cpr", paste0(tid, ".cpr")), threads = NT, tmpdir = TMP)
  vt <- convert_tsv_to_vcf(tsv, file.path(TRAITS, "vcf", paste0(tid, ".vcf.gz")), tid, threads = NT, tmpdir = TMP)
  cv <- data.table(trait = tid, N = N, tsv_write_s = tw, cpr_prepare_s = cc[["prepare"]], cpr_encode_s = cc[["encode"]],
                          cpr_total_s = cc[["total"]], vcf_total_s = vt)
  fwrite(cv, file.path(TRAITS, "conv", paste0(tid, ".csv")))
  log(tid, "N", N, "done; cpr s", round(cc[["total"]], 1), "vcf s", round(vt, 1))
}
log("done")
