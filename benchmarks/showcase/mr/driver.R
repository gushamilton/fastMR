#!/usr/bin/env Rscript
# MR-compute benchmark driver: runs every (scenario, set, size, arm) cell in a fresh R process under /usr/bin/time -v with a
# per-cell cap, arms in RANDOMISED order within each (scenario, set, size) on the same node, and predicted-skip logic.
# usage: Rscript driver.R --rep=1 --out=/path/results --work=/path/work [--dry=1] [--cap=1500]
a <- commandArgs(TRUE); o <- list(rep = "1", out = NA, work = NA, dry = "0", cap = "1500", seed = NA)
for (x in a) { kv <- regmatches(x, regexec("^--([^=]+)=(.*)$", x))[[1]]; if (length(kv) == 3) o[[kv[2]]] <- kv[3] }
REP <- as.integer(o$rep); CAP <- as.numeric(o$cap); DRY <- o$dry == "1"
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
dir.create(o$out, showWarnings = FALSE, recursive = TRUE); dir.create(o$work, showWarnings = FALSE, recursive = TRUE)
stopifnot(Sys.getenv("SLURM_CPUS_PER_TASK") == "8")
set.seed(if (is.na(o$seed)) 1000L + REP else as.integer(o$seed))   # randomisation of arm order (reproducible per replicate)
resf <- file.path(o$out, sprintf("rep%d%s.csv", REP, if (DRY) "_dry" else ""))
agf <- file.path(o$out, sprintf("agreement_rep%d%s.csv", REP, if (DRY) "_dry" else ""))
bootf <- file.path(o$out, sprintf("bootse_rep%d%s.csv", REP, if (DRY) "_dry" else ""))
for (f in c(resf, agf, bootf)) if (file.exists(f)) file.remove(f)
REF <- Sys.getenv("CLUMP_REF"); PREF <- Sys.getenv("CLUMP_PREF")   # bed prefix / pgen prefix
P1 <- Sys.getenv("PLINK1"); P2 <- Sys.getenv("PLINK2")
CPU <- local({ l <- grep("^model name", readLines("/proc/cpuinfo", warn = FALSE), value = TRUE); trimws(sub("^[^:]*:", "", l[1])) })
HOST <- Sys.info()[["nodename"]]; NPROC <- as.integer(system("env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc", intern = TRUE))

# ------------------------------------------------------------------ the plan
plan <- data.frame(scen = character(), set = character(), arm = character(), size = numeric(), exp = numeric(),
                   force_skip = character(), stringsAsFactors = FALSE)
add <- function(scen, set, arm, sizes, exp = 1, skip = numeric(), why = "") {
  plan <<- rbind(plan, data.frame(scen = scen, set = set, arm = arm, size = sizes, exp = exp,
    force_skip = ifelse(sizes %in% skip, why, ""), stringsAsFactors = FALSE)) }
if (!DRY) {
  for (s in c("default", "ivw")) { add("single", s, "tsmr1", c(10, 100, 1000)); add("single", s, "fast1", c(10, 100, 1000)) }
  add("single", "default", "fast8", c(10, 100, 1000))          # fastMR 8 threads only in bootstrap cells
  I <- c(100, 1e3, 1e4, 3e4, 1e5)
  add("many", "ivw", "tsmr1", I)   # TSMR 0.7.11 IVW is ~4x faster than 0.5.7: 100k predicted ~300 s, kept
  add("many", "ivw", "tsmr8", I); add("many", "ivw", "fast1", I)
  for (K in c(10, 50)) { s <- paste0("default_k", K)
    add("many", s, "tsmr1", c(10, 30, 100)); add("many", s, "tsmr8", c(100, 300, 1000))
    add("many", s, "fast1", c(10, 30, 100, 1e3, 1e4)); add("many", s, "fast8", c(10, 30, 100, 1e3, 1e4)) }
  add("many", "default_k10", "tsmr1_seedB", c(10, 30))         # TSMR seed A vs seed B null control for bootstrap SEs
  H <- c(1e3, 1e4, 1e5)
  add("harmonise", "action1", "tsmr1", H, exp = 2); add("harmonise", "action1", "fast1", H)
  add("harmonise", "action2_strandflips", "tsmr1", 1e4, exp = 2); add("harmonise", "action2_strandflips", "fast1", 1e4)
  add("harmonise_many", "action1", "tsmr1", 1000); add("harmonise_many", "action1", "tsmr8", 1000); add("harmonise_many", "action1", "fast1", 1000)
  for (arm in c("tsmr1", "tsmr8", "fast1")) add("steiger", "all", arm, c(1e3, 2e4))
  for (sc in c("heterogeneity", "pleiotropy")) for (arm in c("tsmr1", "tsmr8", "fast1")) add(sc, "all", arm, 1000)
  for (arm in c("tsmr1", "tsmr8", "p2c1")) add("clump", "all", arm, c(10, 100))
  for (arm in c("fast1", "fast8")) add("clump", "all", arm, c(10, 100, 1000))
} else {   # dry run: smallest sizes everywhere, plus fast-only probes at larger sizes for the budget table
  for (s in c("default", "ivw")) { add("single", s, "tsmr1", 10); add("single", s, "fast1", 10) }
  add("single", "default", "fast8", 10)
  add("many", "ivw", "tsmr1", 100); add("many", "ivw", "tsmr8", 100); add("many", "ivw", "fast1", c(100, 1e4)); add("many", "ivw", "tsmr1", 1000)
  add("many", "default_k10", "tsmr1", 10); add("many", "default_k10", "tsmr8", 10); add("many", "default_k10", "tsmr1_seedB", 10)
  add("many", "default_k10", "fast1", c(10, 100, 1000)); add("many", "default_k10", "fast8", c(10, 100, 1000))
  add("many", "default_k50", "tsmr1", 10); add("many", "default_k50", "fast1", 10)
  add("harmonise", "action1", "tsmr1", 1e3); add("harmonise", "action1", "fast1", 1e3)
  add("harmonise", "action2_strandflips", "tsmr1", 1e3); add("harmonise", "action2_strandflips", "fast1", 1e3)
  for (arm in c("tsmr1", "tsmr8", "fast1")) add("harmonise_many", "action1", arm, 100)
  for (arm in c("tsmr1", "tsmr8", "fast1")) add("steiger", "all", arm, 1e3)
  for (sc in c("heterogeneity", "pleiotropy")) for (arm in c("tsmr1", "tsmr8", "fast1", "fast8")) add(sc, "all", arm, if (arm == "fast8") 1000 else 1000)
  for (arm in c("tsmr1", "tsmr8", "p2c1", "fast1", "fast8")) add("clump", "all", arm, 10)
}
ARM_PKG <- c(tsmr1 = "TwoSampleMR", tsmr8 = "TwoSampleMR", tsmr1_seedB = "TwoSampleMR", p2c1 = "plink2", fast1 = "fastMR", fast8 = "fastMR")
ARM_THR <- c(tsmr1 = 1, tsmr8 = 8, tsmr1_seedB = 1, p2c1 = 1, fast1 = 1, fast8 = 8)
COLS <- c("scenario", "size", "method_set", "arm", "package", "threads", "replicate", "hostname", "cpu_model", "nproc", "slurm_cpus",
          "status", "wall_s", "cpu_s", "timing_reps", "proc_wall_s", "proc_cpu_s", "peak_rss_mb", "notes", "pkg_version", "boot_seed")
append_row <- function(r) {
  full <- setNames(rep(list(NA), length(COLS)), COLS); for (n in intersect(names(r), COLS)) full[[n]] <- r[[n]]
  write.table(as.data.frame(full, stringsAsFactors = FALSE), resf, sep = ",", row.names = FALSE, col.names = !file.exists(resf),
              append = file.exists(resf), qmethod = "double") }
ts_tag <- function() format(Sys.time(), "%H:%M:%S")
state <- new.env()   # per (scen,set,arm): last_size, last_wall, dead
rdsdir <- file.path(o$work, "rds"); cdir <- file.path(o$work, "cells"); dir.create(rdsdir, FALSE, TRUE); dir.create(cdir, FALSE, TRUE)
t_task0 <- proc.time()[["elapsed"]]

run_cell <- function(scen, set, arm, size, exp, force_skip) {
  base <- list(scenario = scen, size = size, method_set = set, arm = arm, package = ARM_PKG[[arm]], threads = ARM_THR[[arm]], replicate = REP,
               hostname = HOST, cpu_model = CPU, nproc = NPROC, slurm_cpus = Sys.getenv("SLURM_CPUS_PER_TASK"))
  key <- paste(scen, set, arm, sep = "|"); st <- state[[key]]
  skip <- function(status, note) { r <- base; r$status <- status; r$notes <- note; append_row(r); cat(sprintf("%s SKIP %s\n", ts_tag(), paste(key, size, status))) }
  if (nzchar(force_skip)) return(skip(sub(" .*", "", force_skip), force_skip))
  if (!is.null(st) && st$dead) return(skip("> cap (predicted)", "a smaller size was already censored"))
  if (!is.null(st) && !is.na(st$wall)) {
    pred <- st$wall * (size / st$size)^exp
    if (pred > CAP) return(skip("> cap (predicted)", sprintf("predicted %.0f s = t(%g)=%.3g s x (%g/%g)^%g > cap %g s", pred, st$size, st$wall, size, st$size, exp, CAP))) }
  cellcsv <- file.path(cdir, "cell.csv"); tfile <- file.path(cdir, "time.txt"); logf <- file.path(cdir, "run.log")
  file.remove(Filter(file.exists, c(cellcsv, tfile)))
  cmd <- sprintf("/usr/bin/time -v -o %s timeout -k 10 %d Rscript %s --scen=%s --size=%s --set=%s --arm=%s --rep=%d --rds=%s --cellcsv=%s --ref=%s --pref=%s --p1bin=%s --p2bin=%s > %s 2>&1",
                 tfile, CAP, file.path(HERE, "worker.R"), scen, format(size, scientific = FALSE), set, arm, REP, rdsdir, cellcsv, REF, PREF, P1, P2, logf)
  rc <- system(cmd)
  tl <- if (file.exists(tfile)) readLines(tfile) else character()
  num <- function(p) { l <- grep(p, tl, value = TRUE); if (length(l)) as.numeric(sub("^.*: ", "", l[1])) else NA_real_ }
  pw <- { l <- grep("Elapsed \\(wall", tl, value = TRUE); if (length(l)) { s <- sub("^.*: ", "", l[1]); p <- as.numeric(strsplit(s, ":")[[1]]); sum(p * c(60^(rev(seq_along(p)) - 1))) } else NA_real_ }
  extra <- list(proc_wall_s = pw, proc_cpu_s = sum(num("User time"), num("System time"), na.rm = TRUE), peak_rss_mb = num("Maximum resident") / 1024)
  if (rc %in% c(124, 137) || (!file.exists(cellcsv) && !is.na(pw) && pw >= CAP - 1)) {
    r <- c(base, extra); r$status <- "> cap"; r$notes <- sprintf("process killed at %d s cap (includes data generation)", CAP); append_row(r)
    assign(key, list(size = size, wall = NA, dead = TRUE), envir = state); cat(sprintf("%s CENSORED %s %s\n", ts_tag(), key, size)); return(invisible()) }
  if (rc != 0 || !file.exists(cellcsv)) {
    r <- c(base, extra); r$status <- "error"; r$notes <- gsub("[\",\n]", " ", paste(tail(readLines(logf), 3), collapse = " | ")); append_row(r)
    assign(key, list(size = size, wall = NA, dead = TRUE), envir = state); cat(sprintf("%s ERROR %s %s: %s\n", ts_tag(), key, size, r$notes)); return(invisible()) }
  r <- as.list(read.csv(cellcsv, stringsAsFactors = FALSE)); r[names(extra)] <- extra; append_row(r)
  assign(key, list(size = size, wall = r$wall_s, dead = FALSE), envir = state)
  cat(sprintf("%s ok %s %s wall=%.4g cpu=%.4g rss=%.0fMB (task elapsed %.0f s)\n", ts_tag(), key, size, r$wall_s, r$cpu_s, extra$peak_rss_mb, proc.time()[["elapsed"]] - t_task0))
}

groups <- unique(plan[, c("scen", "set")])
for (gi in seq_len(nrow(groups))) {
  g <- plan[plan$scen == groups$scen[gi] & plan$set == groups$set[gi], ]
  for (sz in sort(unique(g$size))) {
    cells <- g[g$size == sz, ]; cells <- cells[sample(nrow(cells)), ]   # randomised arm order, same node
    for (ci in seq_len(nrow(cells))) with(cells[ci, ], run_cell(scen, set, arm, size, exp, force_skip))
    system(sprintf("Rscript %s --rds=%s --scen=%s --size=%s --set=%s --rep=%d --out=%s --boot=%s", file.path(HERE, "compare.R"), rdsdir,
                   groups$scen[gi], format(sz, scientific = FALSE), groups$set[gi], REP, agf, bootf))
    unlink(list.files(rdsdir, full.names = TRUE))
  }
}
cat(sprintf("TASK_DONE rep=%d dry=%s elapsed=%.0f s\n", REP, DRY, proc.time()[["elapsed"]] - t_task0))
