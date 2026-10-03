#!/usr/bin/env Rscript
# Fixed, seeded query lists shared by every format.
source(file.path(Sys.getenv("STORAGE_SCRIPTS", file.path(Sys.getenv("STORAGE_ROOT", "/user/work/fh6520/showcase/storage"), "scripts")), "common.R"))
dir.create(LISTS, FALSE, TRUE)
d <- load_source(8L)
n <- nrow(d); set.seed(20260801L)
ci <- sort(sample.int(n, 100L))
reg <- d[ci, .(chrom, centre = pos)]
reg[, `:=`(start = pmax(1L, centre - 500000L), end = centre + 500000L)]
# clip to GRCh38 chromosome length: a CompreSSoR region running past the chromosome end leaks into the next
# chromosome's rows (reported separately), and tabix/parquet/VCF would not, so windows stay inside the chromosome
grch38_len <- c(248956422L, 242193529L, 198295559L, 190214555L, 181538259L, 170805979L, 159345973L, 145138636L,
  138394717L, 133797422L, 135086622L, 133275309L, 114364328L, 107043718L, 101991189L, 90338345L, 83257441L,
  80373285L, 58617616L, 64444167L, 46709983L, 50818468L)
reg[, end := pmin(end, grch38_len[as.integer(chrom)])]
fwrite(reg, file.path(LISTS, "regions.tsv"), sep = "\t")
perm <- sample.int(n)                      # nested: keys_25 subset of keys_1000 subset of keys_100000
for (k in c(25L, 1000L, 100000L)) {
  idx <- perm[seq_len(min(k, n))]
  fwrite(d[idx, .(chrom, pos, ref, alt)], file.path(LISTS, sprintf("keys_%d.tsv", k)), sep = "\t")
}
cat("rows", n, "\n")
