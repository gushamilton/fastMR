#!/usr/bin/env Rscript
# Fidelity of every store against the source (read via the standard read path, 8 threads).
source(file.path(Sys.getenv("STORAGE_SCRIPTS", file.path(Sys.getenv("STORAGE_ROOT", "/user/work/fh6520/showcase/storage"), "scripts")), "common.R"))
suppressPackageStartupMessages(library(CompreSSoR))
RES <- file.path(ROOT, "results"); dir.create(RES, FALSE, TRUE)
nthreads_set(8L)
src <- load_source(8L)
n <- nrow(src)
src_nlp <- -log10(src$p); src_z <- src$beta / src$se
err1 <- function(fmt, field, ref, new) {
  ad <- abs(new - ref); ok <- is.finite(ad)
  rel <- ad / abs(ref); rel[!is.finite(rel) | ref == 0] <- NA
  if (all(is.na(rel))) rel[1] <- NA_real_
  data.table(format = fmt, field = field, n = sum(ok), max_abs = max(ad[ok]),
             p999_abs = quantile(ad[ok], 0.999, names = FALSE),
             max_rel = if (all(is.na(rel))) NA_real_ else max(rel, na.rm = TRUE),
             p999_rel = if (all(is.na(rel))) NA_real_ else quantile(rel, 0.999, na.rm = TRUE, names = FALSE),
             n_na_new = sum(!is.finite(new)))
}
errs <- list(); flips <- list(); exact <- list()
flip_tab <- function(fmt, p_new) {
  rbindlist(lapply(c(5e-8, 1e-5, 0.01), function(th) {
    s <- src$p <= th; r <- p_new <= th; r[is.na(r)] <- FALSE
    data.table(format = fmt, threshold = th, n_na_recon_p = sum(is.na(p_new)), n_source = sum(s), n_recon = sum(r),
               lost = sum(s & !r), gained = sum(!s & r), flips = sum(s != r),
               flip_pct_of_source_hits = 100 * sum(s != r) / max(1, sum(s)))
  }))
}

# ---- CompreSSoR
st <- store_path("cpr")
x <- read_sumstats(st, columns = c(CPR_COLS, "z"), threads = 8L); setDT(x)
stopifnot(nrow(x) == n)
kx <- mk_key(x$chromosome, x$base_pair_location, x$other_allele, x$effect_allele)
ks <- mk_key(src$chrom, src$pos, src$ref, src$alt)
store_order_same <- identical(kx, ks)
mi <- match(ks, kx); stopifnot(!anyNA(mi), !anyDuplicated(mi))
x <- x[mi]; rm(kx, ks, mi)       # align store rows to source rows by full key
errs <- c(errs, list(err1("cpr", "beta", src$beta, x$beta), err1("cpr", "se", src$se, x$standard_error),
  err1("cpr", "z", src_z, x$z), err1("cpr", "eaf", src$eaf, x$effect_allele_frequency),
  err1("cpr", "nlp", src_nlp, -log10(x$p_value)), err1("cpr", "p", src$p, x$p_value),
  err1("cpr", "beta_in_se_units", rep(0, n), (x$beta - src$beta) / src$se)))
flips <- c(flips, list(flip_tab("cpr", x$p_value)))
bad <- which(!is.finite(x$p_value) | !is.finite(x$z) | !is.finite(x$beta))
fwrite(cbind(src[bad], data.table(cpr_z = x$z[bad], cpr_beta = x$beta[bad], cpr_se = x$standard_error[bad], cpr_p = x$p_value[bad])),
       file.path(RES, "fidelity_cpr_nonfinite_rows.csv"))
cat("cpr non-finite reconstructed rows:", length(bad), "\n")
# exact p-order domain
ord <- read_pvalue_order(st)
# the order domain refers to store rows; recompute expected order in store-row coordinates
x_store_p <- NULL
xs <- read_sumstats(st, columns = c("chromosome", "base_pair_location", "other_allele", "effect_allele", "p_value"), threads = 8L)
ksto <- mk_key(xs$chromosome, xs$base_pair_location, xs$other_allele, xs$effect_allele)
sp <- src$p[match(ksto, mk_key(src$chrom, src$pos, src$ref, src$alt))]
hit <- which(sp <= 0.01); exp_ord <- hit[order(sp[hit], hit)] - 1L
order_row <- data.table(store_row_order_equals_source_order = store_order_same, n_source_p_le_0.01 = length(hit), n_order = length(ord),
                        identical_to_source_order = identical(as.integer(ord), as.integer(exp_ord)),
                        n_mismatch = if (length(ord) == length(exp_ord)) sum(as.integer(ord) != exp_ord) else NA_integer_)
rh <- which(xs$p_value <= 0.01); rec_ord <- rh[order(xs$p_value[rh], rh)] - 1L
order_row[, `:=`(reconstructed_order_n = length(rec_ord),
  reconstructed_order_identical = identical(as.integer(rec_ord), as.integer(exp_ord)),
  reconstructed_order_mismatch_pct = if (length(rec_ord) == length(exp_ord)) 100 * mean(as.integer(rec_ord) != exp_ord) else NA_real_)]
rm(xs, ksto, sp)
fwrite(order_row, file.path(RES, "fidelity_pvalue_order.csv"))
rm(x); invisible(gc())

# ---- quantised parquet
q <- read_full("qparquet", "read_parquet", 8L)
stopifnot(identical(q$pos, src$pos), identical(q$chrom, src$chrom), identical(q$ref, src$ref), identical(q$alt, src$alt))
errs <- c(errs, list(err1("qparquet", "beta", src$beta, q$beta), err1("qparquet", "se", src$se, q$se),
  err1("qparquet", "z", src_z, q$beta / q$se), err1("qparquet", "eaf", src$eaf, q$eaf),
  err1("qparquet", "nlp", src_nlp, -log10(q$p)), err1("qparquet", "p", src$p, q$p)))
flips <- c(flips, list(flip_tab("qparquet", q$p)))
rm(q); invisible(gc())

# ---- formats that should be exact
cmp_exact <- function(fmt, method) {
  y <- read_full(fmt, method, 8L)
  setcolorder(y, LOGICAL)
  rows_same <- nrow(y) == n
  r <- data.table(format = fmt, method = method, n_rows = nrow(y), same_rows = rows_same)
  if (rows_same) {
    for (cn in c("chrom", "pos", "ref", "alt", "beta", "se", "eaf"))
      r[[paste0("identical_", cn)]] <- identical(y[[cn]], src[[cn]])
    r$max_abs_p <- max(abs(y$p - src$p))
    r$max_rel_p <- max(abs(y$p - src$p) / src$p)
    r$max_abs_nlp <- max(abs(-log10(y$p) - src_nlp))
    r$identical_p <- identical(y$p, src$p)
  }
  r
}
exact <- list(cmp_exact("tsv_gz", "fread_pigz"), cmp_exact("tsv_gz", "arrow_csv"),
              cmp_exact("tsv_bgzip_tabix", "fread_bgzip"), cmp_exact("vcf", "bcftools_fread"),
              cmp_exact("parquet", "read_parquet"))
fwrite(rbindlist(errs), file.path(RES, "fidelity_errors.csv"))
fwrite(rbindlist(flips), file.path(RES, "fidelity_flips.csv"))
fwrite(rbindlist(exact, fill = TRUE), file.path(RES, "fidelity_exactness.csv"))
print(rbindlist(errs)); print(rbindlist(flips)); print(order_row); print(rbindlist(exact, fill = TRUE))

# size of a .cpr written without the exact p-order domain (explains difference vs the README headline)
suppressPackageStartupMessages(library(CompreSSoR))
cd <- data.frame(chromosome = src$chrom, base_pair_location = src$pos, reference_allele = src$ref, alternate_allele = src$alt,
                 effect_allele = src$alt, other_allele = src$ref, beta = src$beta, standard_error = src$se,
                 effect_allele_frequency = src$eaf, p_value = src$p)
tmp <- file.path(ROOT, "work", "cpr_noorder.cpr")
compress_sumstats(cd, tmp, input_build = "GRCh38", store_build = "GRCh38", threads = 8L, pvalue_order = FALSE, overwrite = TRUE)
fwrite(data.table(variant = c("default (pvalue_order=TRUE)", "pvalue_order=FALSE"), bytes = c(path_bytes(store_path("cpr")), path_bytes(tmp))),
       file.path(RES, "cpr_size_order_domain.csv"))
unlink(tmp, recursive = TRUE)
