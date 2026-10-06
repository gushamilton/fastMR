# Re-encode the simulated trait stores with the final CompreSSoR (0.7.0, z10/eaf8/se8+xse), straight from
# the GWAS-SSF TSV via allele_columns (no prepared file). Times each conversion for the e2e amortised line.
# STORE_TAG (default v7) names the output sets traits_<tag>, traits25_<tag> and their conv_<tag> timing dirs.
.libPaths(c(Sys.getenv("LIB_CS"), .libPaths()))
suppressMessages({ library(data.table); library(CompreSSoR) })
S <- "/user/work/fh6520/showcase/prep"
TAG <- Sys.getenv("STORE_TAG", "v7")
sets <- setNames(list(file.path(S, "traits"), file.path(S, "traits25")), paste0(c("traits_", "traits25_"), TAG))
jobs <- rbindlist(lapply(names(sets), function(d) {
  src <- sets[[d]]; dst <- file.path(S, d)
  data.table(set = d, src = src, dst = dst, id = fread(file.path(src, "traits.csv"))$trait) }))
k <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")); nt <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
mine <- jobs[seq_len(.N) %% nt == (k - 1)]
NT <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
for (i in seq_len(nrow(mine))) {
  j <- mine[i]; dir.create(file.path(j$dst, "cpr"), FALSE, TRUE); dir.create(file.path(j$dst, paste0("conv_", TAG)), FALSE, TRUE)
  out <- file.path(j$dst, "cpr", paste0(j$id, ".cpr")); unlink(out, recursive = TRUE)
  t <- system.time(compress_sumstats(file.path(j$src, "tsv", paste0(j$id, ".tsv.gz")), out, input_build = "GRCh38",
        store_build = "GRCh38", allele_columns = c(ref = "other_allele", alt = "effect_allele"), threads = NT))[["elapsed"]]
  fwrite(data.table(trait = j$id, set = j$set, cpr_total_s = t, bytes = sum(file.info(list.files(out, full.names = TRUE))$size),
                    compressor = as.character(packageVersion("CompreSSoR")), node = Sys.info()[["nodename"]]),
         file.path(j$dst, paste0("conv_", TAG), paste0(j$id, ".csv")))
  cat(j$set, j$id, round(t, 1), "s\n")
}
