args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1L]] else getwd()
.libPaths(c(file.path(root, ".local", "Rlib"), "/Users/fergushamilton/projects/twosamplemr-fast/.local/Rlib", .libPaths()))
suppressPackageStartupMessages({ library(fastMR); library(TwoSampleMR) })
d <- read.delim(file.path(root, "inst/extdata/il6_crp_primary_100.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
n <- 50L
g <- list(exposure_beta=matrix(rep(d$beta.exposure,n),nrow=n,byrow=TRUE), outcome_beta=matrix(rep(d$beta.outcome,n),nrow=n,byrow=TRUE), exposure_se=matrix(rep(d$se.exposure,n),nrow=n,byrow=TRUE), outcome_se=matrix(rep(d$se.outcome,n),nrow=n,byrow=TRUE))
methods <- c("ivw","egger","weighted_median","simple_mode","weighted_mode")
timed <- function(fun, reps=3L) unname(vapply(seq_len(reps), function(i) {
  gc(FALSE)
  system.time(fun())[["elapsed"]]
}, numeric(1)))
rows <- do.call(rbind, lapply(c(0L,100L), function(nboot) do.call(rbind, lapply(c(1L,2L,5L,10L), function(threads) { f <- function() fast_mr_grid(g$exposure_beta,g$outcome_beta,g$exposure_se,g$outcome_se,methods=methods,nboot=nboot,seed=20260803,threads=threads); x <- timed(f, if(nboot) 3L else 5L); data.frame(component="MR grid",nboot=nboot,threads=threads,median_seconds=median(x),min_seconds=min(x),pairs=n*n) }))))
set.seed(4102); nc <- 1500L
cd <- data.frame(SNP=paste0("rs",seq_len(nc)),id.exposure="E",pval.exposure=sort(runif(nc)),chr_name=1L,chrom_start=seq(1L,by=1000L,length.out=nc))
ld <- diag(nc); for(i in seq(1L,nc-1L,by=2L)) ld[i,i+1L] <- ld[i+1L,i] <- .8
rownames(ld) <- colnames(ld) <- cd$SNP
f <- function() fast_clump_data(cd,clump_kb=10,clump_r2=.5,ld_matrix=ld)
x <- timed(f,3L)
rows <- rbind(rows,data.frame(component="local LD-matrix clump",nboot=NA,threads=NA,median_seconds=median(x),min_seconds=min(x),pairs=nc))
dir.create(file.path(root,"outputs"),showWarnings=FALSE)
write.csv(rows,file.path(root,"outputs/iteration_scaling.csv"),row.names=FALSE)
md <- c("# fastMR implementation and scaling iteration","","Prior native TwoSampleMR 50x50 reference: 338.919 seconds at nboot=100.","","| component | nboot | threads | median seconds | minimum seconds | rows/pairs |","|---|---:|---:|---:|---:|---:|")
for(i in seq_len(nrow(rows))) md <- c(md,sprintf("| %s | %s | %s | %.4f | %.4f | %d |",rows$component[i],rows$nboot[i],rows$threads[i],rows$median_seconds[i],rows$min_seconds[i],rows$pairs[i]))
writeLines(md,file.path(root,"outputs/iteration_scaling.md")); print(rows,row.names=FALSE)
