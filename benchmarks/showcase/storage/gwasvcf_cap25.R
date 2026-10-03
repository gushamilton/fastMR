# Re-run the one storage cell censored at the old 10-min cap with the 25-min cap: 3 replicates.
source("/user/work/fh6520/showcase/storage/scripts/common.R")
out <- rbindlist(lapply(1:3, function(r) {
  x <- run_op(list(op = "lookup", fmt = "vcf", method = "gwasvcf", threads = 1L, size = 100000L), cap = 1500L)
  x[, `:=`(rep = r, host = Sys.info()[["nodename"]])]; print(x); x }), fill = TRUE)
fwrite(out, file.path(ROOT, "results", "rep_gwasvcf_lookup100k_cap1500.csv"))
