# After CompreSSoR #49 (byte-identical stores): re-time only the .cpr write paths. One replicate per array task.
source("/user/work/fh6520/showcase/storage/scripts/common.R")
REP <- Sys.getenv("SLURM_ARRAY_TASK_ID"); host <- Sys.info()[["nodename"]]
dest <- file.path(ROOT, "work", paste0("cprw", REP)); dir.create(dest, FALSE, TRUE)
w <- run_op(list(op = "write", fmt = "cpr", method = "", threads = 8L, size = 0L, dest = dest), cap = 1500L)
w[, `:=`(rep = REP, host = host, CompreSSoR_lib = LIB_CS)]; unlink(file.path(dest, "*"), recursive = TRUE)
fwrite(w, file.path(ROOT, "results", sprintf("rep%s_cprwrite_a27b32d.csv", REP)))
# e2e amortised conversion line: TSV.gz -> .cpr for one simulated trait (prepare + encode), same code as prep/common.R
source("/user/work/fh6520/showcase/prep/common.R")
t <- convert_tsv_to_cpr("/user/work/fh6520/showcase/prep/traits/tsv/exp01.tsv.gz", file.path(dest, "exp01.cpr"), threads = 8L, tmpdir = dest)
fwrite(data.table(rep = REP, host = host, trait = "exp01", prepare_s = t[["prepare"]], encode_s = t[["encode"]], total_s = t[["total"]]),
       sprintf("/user/work/fh6520/showcase/e2e/results/conversion_cpr_a27b32d_rep%s.csv", REP))
unlink(dest, recursive = TRUE)
