# Shared helpers for the showcase prep + e2e suites.
# Library locations come from environment variables (final libs are swapped in later):
#   LIB_CS (CompreSSoR), LIB_FASTMR (fastMR), LIB_TSMR (TwoSampleMR + ieugwasr)
SHOW <- Sys.getenv("SHOWCASE", "/user/work/fh6520/showcase")
setup_libs <- function(need = c("cs", "fastmr", "tsmr")) {
  base <- "/software/local/languages/miniforge3/envs/r-4.5.1/lib/R/library"
  lib <- c(if ("tsmr" %in% need) Sys.getenv("LIB_TSMR", file.path(SHOW, "lib-tsmr")),
           if ("cs" %in% need) Sys.getenv("LIB_CS", "/user/work/fh6520/regress-audit/lib-cs-final"),
           if ("fastmr" %in% need) Sys.getenv("LIB_FASTMR", "/user/work/fh6520/regress-audit/lib-fm-final"),
           "/user/work/fh6520/r_packages", base)
  # NB: r_packages (old TSMR 0.5.7) is placed AFTER lib-tsmr, so lib-tsmr wins.
  .libPaths(unique(lib[nzchar(lib)]))
  invisible(.libPaths())
}
REF <- Sys.getenv("REF", file.path(SHOW, "prep/ref/EUR_maf01"))
TRAITS <- Sys.getenv("TRAITS", file.path(SHOW, "prep/traits"))

# Prepared table for CompreSSoR: explicit REF/ALT columns, ALT-oriented effects.
# Returns seconds spent (prepare + encode) as list(prepare=, encode=, total=).
convert_tsv_to_cpr <- function(tsv_gz, cpr, threads = 4L, tmpdir = tempdir(), ...) {
  stopifnot(requireNamespace("data.table"), requireNamespace("CompreSSoR"))
  t0 <- proc.time()[["elapsed"]]
  d <- data.table::fread(tsv_gz, nThread = threads)
  d[, `:=`(reference_allele = other_allele, alternate_allele = effect_allele)]
  prep <- file.path(tmpdir, paste0(basename(cpr), ".prepared.tsv.gz"))
  data.table::fwrite(d, prep, sep = "\t", compress = "gzip", nThread = threads)
  rm(d); t1 <- proc.time()[["elapsed"]]
  CompreSSoR::compress_sumstats(prep, cpr, input_build = "GRCh38", store_build = "GRCh38",
                                threads = threads, overwrite = TRUE, ...)
  t2 <- proc.time()[["elapsed"]]
  unlink(prep)
  c(prepare = t1 - t0, encode = t2 - t1, total = t2 - t0)
}

# GWAS-VCF (IEU-style): FORMAT = ES:SE:LP:AF:SS, effect allele = ALT, ID = chr:pos:ref:alt.
# Simplification vs OpenGWAS: no ##META/contig lengths, ID is not an rsID, ES/SE are the
# same rounded values as in the TSV.  Needs bcftools on PATH.
convert_tsv_to_vcf <- function(tsv_gz, vcf_gz, sample_id, threads = 4L, tmpdir = tempdir()) {
  t0 <- proc.time()[["elapsed"]]
  d <- data.table::fread(tsv_gz, nThread = threads)
  lp <- -(pnorm(-abs(d$beta / d$standard_error), log.p = TRUE) + log(2)) / log(10)
  p_ok <- d$p_value > 0
  lp[p_ok] <- -log10(d$p_value[p_ok])
  out <- data.table::data.table(
    `#CHROM` = d$chromosome, POS = d$base_pair_location, ID = d$rsid,
    REF = d$other_allele, ALT = d$effect_allele, QUAL = ".", FILTER = "PASS", INFO = ".",
    FORMAT = "ES:SE:LP:AF:SS",
    S = paste(d$beta, d$standard_error, signif(lp, 4), d$effect_allele_frequency, d$n, sep = ":"))
  data.table::setnames(out, "S", sample_id)
  hdr <- c("##fileformat=VCFv4.2",
    '##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect size estimate relative to the alternative allele">',
    '##FORMAT=<ID=SE,Number=A,Type=Float,Description="Standard error of effect size estimate">',
    '##FORMAT=<ID=LP,Number=A,Type=Float,Description="-log10 p-value for effect estimate">',
    '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency in the association study">',
    '##FORMAT=<ID=SS,Number=A,Type=Float,Description="Sample size used to estimate genetic effect">',
    paste0("##contig=<ID=", 1:22, ">"), paste0("##source=showcase-sim"))
  plain <- file.path(tmpdir, paste0(sample_id, ".vcf"))
  writeLines(hdr, plain)
  data.table::fwrite(out, plain, sep = "\t", append = TRUE, col.names = TRUE, nThread = threads, quote = FALSE)
  rm(out, d)
  rc <- system2("bcftools", c("view", "-Oz", "--threads", threads, "-o", vcf_gz, plain))
  if (rc != 0) stop("bcftools view failed")
  rc <- system2("bcftools", c("index", "-t", "-f", vcf_gz))
  unlink(plain)
  if (rc != 0) stop("bcftools index failed")
  proc.time()[["elapsed"]] - t0
}
