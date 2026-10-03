# Rebuild prep/traits cpr stores with the final CompreSSoR (writer z-bug fix). Writes to cpr_new/.
.libPaths(c(Sys.getenv("LIB_CS"), .libPaths()))
suppressMessages({ library(data.table); library(CompreSSoR) })
source("/user/work/fh6520/showcase/prep/common.R")
TR <- "/user/work/fh6520/showcase/prep/traits"
ids <- sort(sub("\\.tsv\\.gz$", "", list.files(file.path(TR, "tsv"), "\\.tsv\\.gz$")))
k <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")); ch <- 5L
ids <- ids[intersect(seq_along(ids), ((k - 1) * ch + 1):(k * ch))]
dir.create(file.path(TR, "cpr_new"), showWarnings = FALSE); dir.create(file.path(TR, "conv_new"), showWarnings = FALSE)
NT <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
for (tid in ids) {
  out <- file.path(TR, "cpr_new", paste0(tid, ".cpr"))
  cc <- convert_tsv_to_cpr(file.path(TR, "tsv", paste0(tid, ".tsv.gz")), out, threads = NT, tmpdir = Sys.getenv("TMPDIR"))
  fwrite(data.table(trait = tid, cpr_prepare_s = cc[["prepare"]], cpr_encode_s = cc[["encode"]], cpr_total_s = cc[["total"]],
                    compressor = as.character(packageVersion("CompreSSoR")), node = Sys.info()[["nodename"]]),
         file.path(TR, "conv_new", paste0(tid, ".csv")))
  cat(tid, round(cc[["total"]], 1), "s\n")
}
