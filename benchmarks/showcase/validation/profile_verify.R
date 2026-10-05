#!/usr/bin/env Rscript
# Cost of CompreSSoR's shared-identity verification (d158622: hash one store's position + substitution
# streams, byte-compare every other store's with it) for the 25 exposure and 25 outcome stores of the
# 25x25 study, in a fresh process. Rscript profile_verify.R LIB TRAITS_DIR
a <- commandArgs(TRUE); .libPaths(c(a[1], .libPaths())); suppressMessages(library(CompreSSoR))
f <- list.files(file.path(a[2], "cpr"), pattern = "\\.cpr$", full.names = TRUE)
t_open <- system.time(st <- lapply(f, CompreSSoR::open_compressor))[["elapsed"]]
bytes <- sum(vapply(st, function(s) sum(file.size(file.path(s$path, CompreSSoR:::pcodec_identity_file_names(s)))), 0))
t1 <- system.time(CompreSSoR:::pcodec_verify_identity_files(st, threads = 8L))[["elapsed"]]
t2 <- system.time(CompreSSoR:::pcodec_verify_identity_files(st, threads = 8L))[["elapsed"]]
cat(sprintf("stores=%d identity_bytes=%.0f MB open=%.2fs verify_first=%.2fs verify_cached=%.3fs host=%s\n",
            length(st), bytes / 1e6, t_open, t1, t2, Sys.info()[["nodename"]]))
