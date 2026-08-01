args <- commandArgs(trailingOnly = TRUE)
root <- args[[1L]]
.libPaths(c(file.path(root, ".local", "Rlib"), "/Users/fergushamilton/projects/twosamplemr-fast/.local/Rlib", .libPaths()))
suppressPackageStartupMessages({ library(fastMR); library(TwoSampleMR) })
d <- read.delim(file.path(root, "inst/extdata/il6_crp_primary_100.tsv"), check.names=FALSE, stringsAsFactors=FALSE)
scenarios <- data.frame(name=c("balanced_50x50","one_exposure_250_outcomes","one_outcome_250_exposures","ten_exposures_100_outcomes","hundred_exposures_ten_outcomes"), ne=c(50L,1L,250L,10L,100L), no=c(50L,250L,1L,100L,10L))
make_grid <- function(ne,no) list(exposure_beta=matrix(rep(d$beta.exposure,ne),nrow=ne,byrow=TRUE),outcome_beta=matrix(rep(d$beta.outcome,no),nrow=no,byrow=TRUE),exposure_se=matrix(rep(d$se.exposure,ne),nrow=ne,byrow=TRUE),outcome_se=matrix(rep(d$se.outcome,no),nrow=no,byrow=TRUE))
run_native <- function(g, ne, no, method) {
  out <- vector("list", ne*no); k <- 0L
  for (i in seq_len(ne)) for (j in seq_len(no)) {
    k <- k + 1L
    z <- data.frame(SNP=d$SNP,beta.exposure=g$exposure_beta[i,],beta.outcome=g$outcome_beta[j,],se.exposure=g$exposure_se[i,],se.outcome=g$outcome_se[j,],id.exposure="E",id.outcome="O",exposure="E",outcome="O",mr_keep=TRUE)
    out[[k]] <- suppressMessages(suppressWarnings(TwoSampleMR::mr(z,method_list=method)))
  }
  do.call(rbind,out)
}
methods <- data.frame(code=c("uwr","sign"),native=c("mr_uwr","mr_sign"),stringsAsFactors=FALSE)
rows <- list(); k <- 0L
for (m in seq_len(nrow(methods))) for (s in seq_len(nrow(scenarios))) {
  sc <- scenarios[s,]; g <- make_grid(sc$ne,sc$no); pairs <- sc$ne*sc$no; code <- methods$code[m]; native_method <- methods$native[m]
  invisible(fast_mr_grid(g$exposure_beta,g$outcome_beta,g$exposure_se,g$outcome_se,methods=code,nboot=0,threads=5))
  a <- system.time(fast <- fast_mr_grid(g$exposure_beta,g$outcome_beta,g$exposure_se,g$outcome_se,methods=code,nboot=0,threads=5))[["elapsed"]]
  b <- system.time(native <- run_native(g,sc$ne,sc$no,native_method))[["elapsed"]]
  beta <- fast$b-native$b; p <- fast$pval-native$pval; k <- k+1L
  rows[[k]] <- data.frame(method=code,scenario=sc$name,pairs=pairs,fastMR_seconds=a,TwoSampleMR_seconds=b,speedup=b/a,max_abs_beta_difference=max(abs(beta),na.rm=TRUE),max_abs_pval_difference=max(abs(p),na.rm=TRUE))
}
out <- do.call(rbind,rows); write.csv(out,file.path(root,"outputs/native_tsmr_uwr_sign_benchmark.csv"),row.names=FALSE)
md <- c("# Native TwoSampleMR deterministic method benchmarks","","IL6 fixture: 82 SNPs; nboot=0; fastMR threads=5.","","| method | scenario | pairs | fastMR s | TwoSampleMR s | speedup | max beta delta | max p-value delta |","|---|---|---:|---:|---:|---:|---:|---:|")
for(i in seq_len(nrow(out))) md <- c(md,sprintf("| %s | %s | %d | %.3f | %.3f | %.2fx | %.3e | %.3e |",out$method[i],out$scenario[i],out$pairs[i],out$fastMR_seconds[i],out$TwoSampleMR_seconds[i],out$speedup[i],out$max_abs_beta_difference[i],out$max_abs_pval_difference[i]))
writeLines(md,file.path(root,"outputs/native_tsmr_uwr_sign_benchmark.md"))
print(out,row.names=FALSE)
