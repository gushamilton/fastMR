# Seeded random variant keys from the shared 1000G EUR panel (= the simulated traits' variant set).
source(file.path(Sys.getenv("MULTI_SCRIPTS", "/user/work/fh6520/showcase/storage/multi/scripts"), "mcommon.R"))
d <- file.path(MROOT, "lists"); dir.create(d, FALSE, TRUE)
sizes <- as.integer(strsplit(Sys.getenv("KEY_SIZES", "1000,100000"), ",")[[1]])
if (all(file.exists(file.path(d, sprintf("keys_%d.tsv", sizes))))) quit(save = "no")
pv <- fread(cmd = "grep -v '^##' /user/work/fh6520/showcase/prep/ref/EUR_maf01.pvar", select = 1:5,
            colClasses = list(character = c(1L, 3L, 4L, 5L)))
setnames(pv, c("chrom", "pos", "id", "ref", "alt"))
for (k in sizes) {
  set.seed(20261002L + k)
  x <- pv[sort(sample.int(nrow(pv), k)), .(chrom, pos, ref, alt)]
  tmp <- tempfile(tmpdir = d); fwrite(x, tmp, sep = "\t"); file.rename(tmp, file.path(d, sprintf("keys_%d.tsv", k)))
}
