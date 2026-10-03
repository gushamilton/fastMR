#!/usr/bin/env Rscript
# Rebuild only the GWAS-VCFs from the existing TSV.gz (LP now written at 7 s.f.).  Task chunking as in gen_traits.R.
src_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(src_dir, "common.R")); suppressMessages(library(data.table))
NT <- as.integer(Sys.getenv("NTHREADS", "16")); setDTthreads(NT)
tr <- fread(file.path(TRAITS, "traits.csv"))$trait
ch <- as.integer(Sys.getenv("TASK_CHUNK", "3")); k <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
todo <- tr[intersect(seq_along(tr), ((k - 1) * ch + 1):(k * ch))]
dir.create(file.path(TRAITS, "conv_vcf"), showWarnings = FALSE)
for (id in todo) {
  out <- file.path(TRAITS, "vcf", paste0(id, ".vcf.gz"))
  vt <- convert_tsv_to_vcf(file.path(TRAITS, "tsv", paste0(id, ".tsv.gz")), out, id, threads = NT, tmpdir = tempdir())
  fwrite(data.table(trait = id, vcf_total_s = vt, lp_signif = 7L), file.path(TRAITS, "conv_vcf", paste0(id, ".csv")))
  cat(id, "done", round(vt, 1), "s
")
}
