#!/usr/bin/env Rscript
# One-off prep job: create the stores used for reads (timed first writes recorded), query lists, fidelity.
source(file.path(Sys.getenv("STORAGE_SCRIPTS", file.path(Sys.getenv("STORAGE_ROOT", "/user/work/fh6520/showcase/storage"), "scripts")), "common.R"))
RES <- file.path(ROOT, "results"); dir.create(RES, FALSE, TRUE)
dir.create(STORES, FALSE, TRUE)
host <- system("hostname", intern = TRUE)
PF <- Sys.getenv("PREP_FORMATS", "")
if (!nzchar(PF)) stopifnot(system(sprintf("Rscript %s/make_lists.R", SCRIPTS)) == 0L)
out <- list()
for (fmt in if (nzchar(PF)) strsplit(PF, ",")[[1]] else c("tsv_gz", "tsv_bgzip_tabix", "parquet", "qparquet", "cpr", "vcf")) {
  cat("writing", fmt, "\n")
  r <- run_op(list(op = "write", fmt = fmt, method = "", threads = 8L, size = 0L), cap = 7200L)
  r[, `:=`(host = host, stage = "prep")]
  out[[fmt]] <- r; print(r)
  fwrite(rbindlist(out, fill = TRUE), file.path(RES, if (nzchar(PF)) paste0("prep_writes_", gsub(",", "-", PF), ".csv") else "prep_writes.csv"))
}
sz <- rbindlist(lapply(FORMATS, function(f) data.table(format = f, data_bytes = path_bytes(store_path(f)),
      index_bytes = path_bytes(index_path(f)))))
sz[, total_bytes := data_bytes + fifelse(is.na(index_bytes), 0, index_bytes)]
fwrite(sz, file.path(RES, "sizes.csv")); print(sz)
cat("fidelity\n")
st <- system(sprintf("/usr/bin/time -v -o %s/logs/fidelity.time Rscript %s/fidelity.R > %s/logs/fidelity.log 2>&1", ROOT, SCRIPTS, ROOT))
cat("fidelity exit", st, "\n")
