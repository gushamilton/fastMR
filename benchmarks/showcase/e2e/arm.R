#!/usr/bin/env Rscript
# One end-to-end arm in a fresh R process.
# Usage: Rscript arm.R ARM E O THREADS OUTDIR     ARM in A, A8, B, C
#   A  : TSV.gz + TwoSampleMR (MR serial);  A8: same, MR via parallel::mclapply over pairs (8 cores)
#   B  : GWAS-VCF + bcftools + TwoSampleMR (MR serial)
#   C  : .cpr + fastMR (THREADS = fastMR threads / io_threads / plink2 threads)
# Env: LIB_CS LIB_FASTMR LIB_TSMR REF TRAITS SEED  (see common.R); PLINK1 PLINK2 optional paths
a <- commandArgs(TRUE); ARM <- a[1]; E <- as.integer(a[2]); O <- as.integer(a[3]); TH <- as.integer(a[4]); OUT <- a[5]
src_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(src_dir, "common.R"))
setup_libs(if (ARM == "C") c("cs", "fastmr", "tsmr") else c("tsmr"))
suppressMessages(library(data.table)); setDTthreads(TH)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SEED <- as.integer(Sys.getenv("SEED", "1")); P1 <- 0.001; KB <- 10000; PTHR <- 5e-8
PLINK1 <- Sys.getenv("PLINK1", Sys.which("plink")); PLINK2 <- Sys.getenv("PLINK2", Sys.which("plink2"))
tr <- fread(file.path(TRAITS, "traits.csv")); expo <- sprintf("exp%02d", seq_len(E)); outc <- sprintf("out%02d", seq_len(O))
Nof <- setNames(tr$N, tr$trait)
cpuline <- system("grep -m1 'model name' /proc/cpuinfo", intern = TRUE)
meta <- data.table(arm = ARM, E = E, O = O, threads = TH, host = Sys.info()[["nodename"]], cpu = sub(".*: ", "", cpuline),
                   slurm_cpus = Sys.getenv("SLURM_CPUS_PER_TASK"), r = R.version.string,
                   fastmr = tryCatch(as.character(packageVersion("fastMR")), error = function(e) NA),
                   compressor = tryCatch(as.character(packageVersion("CompreSSoR")), error = function(e) NA),
                   tsmr = tryCatch(as.character(packageVersion("TwoSampleMR")), error = function(e) NA))
fwrite(meta, file.path(OUT, "meta.csv"))
if (nzchar(Sys.getenv("SLURM_CPUS_PER_TASK")) && as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) < TH) stop("threads exceed allocation")
stages <- file.path(OUT, "stages.csv"); if (file.exists(stages)) file.remove(stages)
cpu0 <- proc.time()
stage <- function(name, expr) {
  t0 <- proc.time(); on.exit(NULL); v <- force(expr); t1 <- proc.time() - t0
  fwrite(data.table(stage = name, wall_s = t1[["elapsed"]], cpu_s = t1[["user.self"]] + t1[["sys.self"]] + t1[["user.child"]] + t1[["sys.child"]]),
         stages, append = file.exists(stages)); invisible(v)
}
rd <- function(f, cols = NULL) {   # best TSV.gz reader on this box: pigz -dc | fread
  fread(cmd = paste("pigz -dc -p", TH, shQuote(f)), nThread = TH, select = cols)
}
tsmr_fmt <- function(d, id, type) {
  d <- as.data.frame(d); d$Phenotype <- id; d$id <- id
  TwoSampleMR::format_data(d, type = type, phenotype_col = "Phenotype", id_col = "id", snp_col = "rsid", beta_col = "beta",
    se_col = "standard_error", eaf_col = "effect_allele_frequency", effect_allele_col = "effect_allele",
    other_allele_col = "other_allele", pval_col = "p_value", samplesize_col = "n", chr_col = "chromosome", pos_col = "base_pair_location")
}
vcf_query <- function(vcf, extra) {   # GWAS-VCF -> data.table in the same logical columns
  fmt <- "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t[%ES\\t%SE\\t%LP\\t%AF\\t%SS]\\n"  # literal backslash-t for bcftools
  d <- fread(cmd = paste("bcftools query", extra, "-f", shQuote(fmt), shQuote(vcf)),
             header = FALSE, col.names = c("chromosome", "base_pair_location", "rsid", "other_allele", "effect_allele",
                                          "beta", "standard_error", "lp", "effect_allele_frequency", "n"))
  d[, p_value := 10^-lp]; d[, lp := NULL]; d
}
res <- list(); pair_key <- function(d) paste(d$id.exposure, d$id.outcome, sep = "|")

if (ARM %in% c("A", "A8", "B")) {
  suppressMessages({library(TwoSampleMR); library(ieugwasr)})
  # ---- select: instruments p < 5e-8 per exposure ----
  expdat <- stage("select", {
    L <- lapply(expo, function(id) {
      d <- if (ARM == "B") vcf_query(file.path(TRAITS, "vcf", paste0(id, ".vcf.gz")), "-i 'FORMAT/LP>7.30103'")
           else { x <- rd(file.path(TRAITS, "tsv", paste0(id, ".tsv.gz"))); x[p_value < PTHR] }
      tsmr_fmt(d, id, "exposure")
    })
    do.call(rbind, L)
  })
  # NB: TwoSampleMR::format_data() lower-cases SNP ids (chr:pos:a:g); the plink reference and files use upper-case alleles,
  # so ids are upper-cased again for clumping/extraction and for the agreement sets.
  res$instruments_pre <- split(toupper(expdat$SNP), expdat$id.exposure)
  # ---- clump: local plink 1.9 via ieugwasr::ld_clump ----
  expc <- stage("clump", {
    cl <- ieugwasr::ld_clump(dplyr::tibble(rsid = toupper(expdat$SNP), pval = expdat$pval.exposure, id = expdat$id.exposure),
                             clump_kb = KB, clump_r2 = P1, clump_p = 1, plink_bin = PLINK1, bfile = REF)
    expdat[paste(toupper(expdat$SNP), expdat$id.exposure) %in% paste(cl$rsid, cl$id), ]
  })
  res$instruments <- split(toupper(expc$SNP), expc$id.exposure)
  # ---- extract: outcome rows at the instrument union ----
  outdat <- stage("extract", {
    snps <- toupper(unique(expc$SNP))
    L <- lapply(outc, function(id) {
      d <- if (ARM == "B") {
        pos <- unique(expc[, c("chr.exposure", "pos.exposure")]); pos <- pos[order(as.integer(pos[[1]]), pos[[2]]), ]
        pf <- tempfile(fileext = ".tsv"); fwrite(pos, pf, sep = "\t", col.names = FALSE)
        vcf_query(file.path(TRAITS, "vcf", paste0(id, ".vcf.gz")), paste("-R", pf))
      } else { x <- rd(file.path(TRAITS, "tsv", paste0(id, ".tsv.gz"))); x }
      d <- d[rsid %in% snps]
      tsmr_fmt(d, id, "outcome")
    })
    do.call(rbind, L)
  })
  # ---- harmonise (forward strand: action = 1) ----
  dat <- stage("harmonise", TwoSampleMR::harmonise_data(expc, outdat, action = 1))
  res$harm <- split(toupper(dat$SNP[dat$mr_keep]), pair_key(dat)[dat$mr_keep])
  saveRDS(dat, file.path(OUT, "harmonised.rds"))
  # ---- MR, default method set, default nboot (1000) ----
  mrres <- stage("mr", {
    set.seed(SEED)
    if (ARM == "A8") {
      pk <- pair_key(dat); up <- unique(pk)
      parts <- parallel::mclapply(seq_along(up), function(i) { set.seed(SEED + i); TwoSampleMR::mr(dat[pk == up[i], ]) },
                                  mc.cores = TH, mc.preschedule = FALSE)
      bad <- vapply(parts, function(x) inherits(x, "try-error") || is.null(x), logical(1)); if (any(bad)) warning(sum(bad), " pairs failed")
      rbindlist(parts[!bad], fill = TRUE)
    } else as.data.table(TwoSampleMR::mr(dat))
  })
  stg <- stage("steiger", as.data.table(TwoSampleMR::steiger_filtering(dat)))
  stage("write", { fwrite(mrres, file.path(OUT, "mr.tsv"), sep = "\t"); fwrite(stg, file.path(OUT, "steiger.tsv"), sep = "\t") })
  res$mr <- as.data.frame(mrres)
}

if (ARM == "C") {
  suppressMessages({library(CompreSSoR); library(fastMR); loadNamespace("arrow")})   # loaded before any timed stage, like the other arms
  ef <- setNames(file.path(TRAITS, "cpr", paste0(expo, ".cpr")), expo)
  of <- setNames(file.path(TRAITS, "cpr", paste0(outc, ".cpr")), outc)
  # select + clump are one public call (candidate read from the p-value flag/order domains, then graph clump)
  cl <- stage("clump", fast_clump_compressed(ef, pvalue_threshold = PTHR, candidate_source = "pvalue_flag",
          pvalue_order = "require_exact", partition = "graph", pfile = REF, plink2_bin = PLINK2,
          clump_kb = KB, clump_r2 = P1, threads = TH, io_threads = TH))
  inst <- cl$instruments[expo]
  res$instruments <- lapply(inst, as.character)
  methods <- c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode")
  mrres <- stage("mr_total_incl_extract", fast_mr_compressed(ef, of, inst, methods = methods, nboot = 1000, seed = SEED,
                                                threads = TH, io_threads = TH))
  tim <- attr(mrres, "compressed_input")$timing
  fwrite(data.table(stage = c("extract", "harmonise", "mr"), wall_s = c(tim$io_seconds, 0, tim$estimator_seconds), cpu_s = NA_real_),
         stages, append = TRUE)
  res$counts <- attr(mrres, "compressed_input")$counts
  # Steiger: the store carries no sample size, so N comes from the simulation design table.
  stg <- stage("steiger", {
    # Stores are read in parallel (one per worker) and joined with data.table: the
    # reads are I/O the TSMR arms already did in their extract stage.
    rd1 <- function(path, keys) as.data.table(fast_read_compressed(path, variants = keys,
              columns = c("beta", "standard_error", "effect_allele_frequency", "p_value")))
    allk <- unique(unlist(inst, use.names = FALSE))
    par <- function(ids, f) if (TH > 1L) parallel::mclapply(ids, f, mc.cores = min(TH, length(ids))) else lapply(ids, f)
    ex <- rbindlist(par(expo, function(id) rd1(ef[[id]], inst[[id]])[, id.exposure := id]))
    ou <- rbindlist(par(outc, function(id) rd1(of[[id]], allk)[, id.outcome := id]))
    setnames(ex, c("beta", "standard_error", "effect_allele_frequency", "p_value"),
             c("beta.exposure", "se.exposure", "eaf.exposure", "pval.exposure"))
    setnames(ou, c("beta", "standard_error", "effect_allele_frequency", "p_value"),
             c("beta.outcome", "se.outcome", "eaf.outcome", "pval.outcome"))
    d <- ex[ou, on = "variant_key", nomatch = NULL, allow.cartesian = TRUE]
    d[, `:=`(SNP = variant_key, exposure = id.exposure, outcome = id.outcome, units.exposure = "", units.outcome = "",
             samplesize.exposure = Nof[id.exposure], samplesize.outcome = Nof[id.outcome], mr_keep = TRUE)]
    setorder(d, id.exposure, id.outcome)
    d <- as.data.frame(d); res$harm <- split(d$SNP, paste(d$id.exposure, d$id.outcome, sep = "|"))
    fast_mr_steiger_filtering(d)
  })
  stage("write", { fast_write_parquet(as.data.frame(mrres), file.path(OUT, "mr.parquet"), overwrite = TRUE)
                   fast_write_parquet(as.data.frame(stg), file.path(OUT, "steiger.parquet"), overwrite = TRUE) })
  fwrite(as.data.table(mrres), file.path(OUT, "mr_for_agreement.tsv"), sep = "\t")   # untimed copy for the agreement script
  res$mr <- as.data.frame(mrres)
}
tot <- proc.time() - cpu0
res$total_cpu_s <- tot[["user.self"]] + tot[["sys.self"]] + tot[["user.child"]] + tot[["sys.child"]]
saveRDS(res, file.path(OUT, "result.rds"))
cat("ARM", ARM, "OK; stages:\n"); print(fread(stages))
