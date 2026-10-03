source(file.path(Sys.getenv("STORAGE_SCRIPTS", "/user/work/fh6520/showcase/storage/scripts"), "common.R"))
suppressPackageStartupMessages({library(CompreSSoR); library(arrow)})
P <- store_path("cpr"); out <- list()
sets <- list(z_only = "z", beta_se_eaf = c("beta", "standard_error", "effect_allele_frequency"),
  numeric_all = c("base_pair_location", "beta", "standard_error", "effect_allele_frequency"),
  numeric_all_plus_p = c("base_pair_location", "beta", "standard_error", "effect_allele_frequency", "p_value"),
  identity_only = c("chromosome", "base_pair_location", "effect_allele", "other_allele"),
  all_logical = CPR_COLS)
for (T in c(1L, 8L)) for (nm in names(sets)) for (rep in 1:3) {
  t <- system.time(x <- read_sumstats(P, columns = sets[[nm]], threads = T))[["elapsed"]]
  out[[length(out) + 1L]] <- data.table(T = T, columns = nm, rep = rep, seconds = t, nrow = nrow(x)); rm(x); invisible(gc())
}
nthreads_set(8L)
for (rep in 1:3) { t <- system.time(x <- arrow::read_parquet(store_path("parquet"), col_select = c("pos", "beta", "se", "eaf", "p")))[["elapsed"]]
  out[[length(out) + 1L]] <- data.table(T = 8L, columns = "parquet_numeric_only", rep = rep, seconds = t, nrow = nrow(x)); rm(x) }
for (rep in 1:3) { t <- system.time(x <- arrow::read_parquet(store_path("parquet")))[["elapsed"]]
  out[[length(out) + 1L]] <- data.table(T = 8L, columns = "parquet_all", rep = rep, seconds = t, nrow = nrow(x)); rm(x) }
r <- rbindlist(out); fwrite(r, file.path(ROOT, "results", "diag_cpr_read_columns.csv"))
print(r[, .(median_s = median(seconds), min_s = min(seconds)), by = .(T, columns)])
