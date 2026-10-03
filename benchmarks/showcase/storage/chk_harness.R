source(file.path(Sys.getenv("STORAGE_SCRIPTS"), "common.R"))
for (s in list(list(op="fullread",fmt="cpr",method="read_sumstats",threads=8L),
               list(op="fullread",fmt="parquet",method="read_parquet",threads=8L),
               list(op="fullread",fmt="tsv_gz",method="fread_pigz",threads=1L),
               list(op="region",fmt="cpr",method="read_sumstats",threads=1L),
               list(op="lookup",fmt="parquet",method="semi_join",threads=8L,size=1000L))) {
  print(run_op(s)[, .(op, format, method, threads, status, load_s, seconds, wall_s, n_rows)]) }
