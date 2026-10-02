.libPaths(c("/Users/fh6520/projects/fastMR/.local/Rlib", .libPaths())); suppressMessages(library(fastMR))
E <- as.integer(commandArgs(TRUE)[1]); S <- 100L; set.seed(1)
bx <- matrix(rnorm(E*S,.05,.02),E,S); by <- matrix(rnorm(E*S,0,.02),E,S)
sx <- matrix(runif(E*S,.005,.02),E,S); sy <- matrix(runif(E*S,.005,.02),E,S)
rownames(bx) <- paste0("e",1:E); rownames(by) <- paste0("o",1:E)
gc(reset=TRUE); t <- system.time(n <- fastMR:::fastmr_grid_native(bx,by,sx,sy,"ivw",0L,NULL,1L,1,20))[["elapsed"]]
g1 <- gc(); cat(sprintf("E=O=%d native %.2fs heap_peak %.0fMB native_obj %.0fMB\n", E, t, sum(g1[,ncol(g1)]), object.size(n)/2^20))
gc(reset=TRUE); t <- system.time(r <- fast_mr_grid(bx,by,sx,sy,methods="ivw",nboot=0))[["elapsed"]]
g2 <- gc(); cat(sprintf("E=O=%d fast_mr_grid tidy %.2fs heap_peak %.0fMB tidy_obj %.0fMB\n", E, t, sum(g2[,ncol(g2)]), object.size(r)/2^20))
if (requireNamespace("arrow", quietly=TRUE)) { p <- tempfile(fileext=".parquet"); t <- system.time(fast_write_parquet(r, p))[["elapsed"]]; cat(sprintf("parquet write %.2fs %0.fMB\n", t, file.size(p)/2^20)) }
