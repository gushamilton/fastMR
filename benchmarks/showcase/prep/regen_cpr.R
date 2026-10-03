#!/usr/bin/env Rscript
# Re-encode only the .cpr stores from the existing TSV.gz (e.g. after swapping LIB_CS).  Task chunking as in gen_traits.R.
src_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(src_dir, "common.R")); setup_libs("cs"); suppressMessages(library(data.table))
NT <- as.integer(Sys.getenv("NTHREADS", "16")); setDTthreads(NT)
tr <- fread(file.path(TRAITS, "traits.csv"))$trait
ch <- as.integer(Sys.getenv("TASK_CHUNK", "3")); k <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
todo <- tr[intersect(seq_along(tr), ((k - 1) * ch + 1):(k * ch))]
dir.create(file.path(TRAITS, "conv_cpr"), showWarnings = FALSE)
cat("CompreSSoR", as.character(packageVersion("CompreSSoR")), "from", find.package("CompreSSoR"), "\n")
for (id in todo) {
  out <- file.path(TRAITS, "cpr", paste0(id, ".cpr")); unlink(out, recursive = TRUE)
  cc <- convert_tsv_to_cpr(file.path(TRAITS, "tsv", paste0(id, ".tsv.gz")), out, threads = NT, tmpdir = tempdir())
  fwrite(data.table(trait = id, cpr_prepare_s = cc[["prepare"]], cpr_encode_s = cc[["encode"]], cpr_total_s = cc[["total"]],
                    compressor = as.character(packageVersion("CompreSSoR"))), file.path(TRAITS, "conv_cpr", paste0(id, ".csv")))
  cat(id, "done", round(cc[["total"]], 1), "s\n")
}
