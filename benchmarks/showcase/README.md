# Showcase benchmarks: fastMR + CompreSSoR against the standard MR toolchain

Scripts behind the published comparison of fastMR + CompreSSoR with TwoSampleMR,
PLINK, TSV.gz, bgzip + tabix, GWAS-VCF and Parquet. Everything ran on BluePebble
(University of Bristol HPC, Slurm account `sscm013902`) under
`/user/work/fh6520/showcase/`. The scripts use those absolute paths; edit the
roots at the top of each `common.R` to run elsewhere. `PLAN.md` is the design,
including the revisions made after review (section 6).

| Folder | On BluePebble | What it does |
|---|---|---|
| `prep/` | `showcase/prep/` | Builds the EUR MAF ≥ 0.01 LD reference (bed + pgen) and 20 simulated GWAS (TSV.gz, GWAS-VCF, `.cpr`); `regen_*` rebuild one format |
| `storage/` | `showcase/storage/scripts/` | Single real GWAS (FinnGen, 10M variants) in each format: write, full read, regional query, lookups, fidelity |
| `multi/` | `showcase/storage/multi/scripts/` | 1–20 simulated GWAS per format: instrument selection and variant extraction across files |
| `mr/` | `showcase/mr/` | MR compute: fastMR vs TwoSampleMR on harmonised data, 1 and 8 cores, plus clumping |
| `e2e/` | `showcase/e2e/` | End-to-end 1×1 and 10×10 studies, files on disk to results on disk, per-stage timing and agreement |
| `build/` | `showcase/final/` | Offline builds of the package versions under test |
| `validation/` | `showcase/e2e/`, `showcase/clumpdebug/` | One-off checks: clumping against PLINK 1.9/2, the single-call clumping change, profiling |
| `results_v8/` | `showcase/summary_v8/` | Old-vs-new summary CSVs from `e2e/compare_v8.R` and the scaling model |

Each suite is one Slurm array with one replicate per task (3 tasks, different nodes),
cells in fresh R processes under `/usr/bin/time -v`, arms in randomised order, a
25-minute cap per cell and no extrapolation. Aggregate with
`mr/aggregate.R`, `storage/aggregate.R`, `multi/aggregate_multi.R` and
`e2e/aggregate_e2e.R` (after `e2e/agree.R` per replicate and size).

Versions used for the published numbers (v8 rerun, 2026-10-05): fastMR 0.2.0 at `543cae2` and CompreSSoR 0.7.2
at `d158622` (library `final/lib7`, built by `build/final_build7.sbatch`; stores re-encoded by `prep/reencode_v8.sbatch`,
store format unchanged), TwoSampleMR 0.7.11 (`f492045`), ieugwasr 1.2.0, PLINK 1.9 b7.7, PLINK2 2.00a6.8,
bcftools/htslib 1.19, R 4.5.1, all on Intel Xeon Gold 6226R (Slurm constraint `Cascadelake`), 8 Slurm CPUs
(4 physical cores with hyperthreading). The fastMR arms come from `e2e/e2e_v8.sbatch` (1x1 to 50x50); the
TwoSampleMR and GWAS-VCF arms do not depend on fastMR or CompreSSoR and are reused from the v7 runs on the same
node type (`e2e/assemble_v8.sh` links them next to the new fastMR arms). `e2e/agree_v8.sbatch` then runs
`agree.R` (instrument, harmonised-set and estimate agreement), `wald_check.R` (single-instrument Wald ratio
pairs built from the strongest instrument per exposure, since the showcase designs have none),
`scaling_model.R` (per-stage a + b*n + c*n^2 model fitted on 1x1 and 10x10, checked against the measured
25x25 and 50x50) and `compare_v8.R` (old-vs-new tables). The storage and multi-trait cpr rows were re-timed
with `storage/prep_v8.sbatch` followed by `storage/rep.sbatch` (`RESULT_TAG=v8 STORAGE_FORMATS=cpr`) and
`multi/multi.sbatch` (`RESULT_TAG=repv8 MULTI_FORMATS=cpr CPR_DIR=prep/traits_v8/cpr`).

The previous build (v7): fastMR 0.2.0 at `fb62057`, CompreSSoR 0.7.0 at `f91bb8d` (store format 0.4.6; stores
re-encoded by `prep/reencode_v7.*`), arms from `e2e/e2e_v7.sbatch`. `storage/run_rep.R` and `multi/run_multi.R`
take `STORAGE_FORMATS` / `MULTI_FORMATS` to re-time a single format.

## v8 / v9 rerun against v7

- v7: fastMR `fb62057` + CompreSSoR `f91bb8d` (lib6), the previous published numbers.
- v8: fastMR `543cae2` + CompreSSoR `d158622` (lib7), merged main.
- v9: fastMR `671c85f` (PR #37, branch `fix/clump-flag-preread`) + CompreSSoR `d158622` (lib8). It is v8 plus the
  fix for the clump-stage regression v8 introduced.

All cells ran on Intel Xeon Gold 6226R with 8 Slurm CPUs (4 physical cores with hyperthreading). Values are medians
of 3 replicates in fresh processes; min-max are in `results_v8/` and `results_v9/` (`e2e/compare_v8.R`,
`NEW_TAG=v8|v9`). The TwoSampleMR and GWAS-VCF arms are the v7 measurements, since they do not use this code. v9's
MR outputs are byte-identical to v8's in all 18 C1/C8 cells.

End to end, wall seconds (fastMR C8 = 8 CPUs):

| Study | v7 C8 | v8 C8 | v9 C8 | TSMR, MR serial (1 core) | TSMR, MR over 8 CPUs | GWAS-VCF + TSMR |
|---|---|---|---|---|---|---|
| 1x1 | 7.38 | 8.80 | 7.62 | 48.97 | - | 23.0 |
| 10x10 | 10.87 | 16.44 | 7.74 | 902.5 | 518.6 | 649.6 |
| 25x25 | 23.17 | 28.53 | 20.87 | 4,346 | 1,613 | 3,562.5 |
| 50x50 | 109.70 | 108.42 | 87.33 | not run | not run | not run |

fastMR C1 (1 thread), v7 / v9 in seconds: 1x1 4.84 / 3.53, 10x10 21.04 / 16.94, 25x25 62.63 / 64.23,
50x50 300.85 / 267.79.

v9 C8 speedups:
- 25x25: 208x vs TSMR on 1 core, 77x vs TSMR on 8 CPUs, 171x vs GWAS-VCF (v7: 188x, 70x, 154x).
- 10x10: 117x, 67x and 84x (v7: 83x, 48x, 60x).

Per stage at 25x25 C8 (v7 / v9): clump 7.46 / 7.21 s, extract 4.70 / 4.41 s, MR 6.44 / 3.90 s, Steiger
0.16 / 0.24 s. At 50x50: clump 28.8 / 27.3 s, MR 55.1 / 35.7 s.

Same-node interleaved A/B (`e2e/e2e_ab.sbatch`, `BUILDS="v7 v8 v9"`), C8, 6 cells per build and size:

| Study | v7 | v8 | v9 |
|---|---|---|---|
| 10x10 | 13.21 s (clump 5.20) | 16.33 s (clump 8.35) | 11.18 s (clump 5.17) |
| 25x25 | 23.56 s (clump 7.30) | 29.49 s (clump 16.82) | 19.78 s (clump 7.92) |

The v8 clump regression came from fastMR's batched-candidate guard (#33). It decoded every store's flagged row
ids, one store at a time, before the batch read (`validation/prof_clump_ab.R`: about 7.3 s of a 9.8 s call at 25
exposures). CompreSSoR's identity verification costs only 0.2-0.7 s for 50 stores (`validation/profile_verify.R`).
v9 skips that pre-read when CompreSSoR reports `candidates_batch_rows_checked`.

**Agreement (v7 = v8 = v9)**
- Instruments are identical to TSMR: 1,259/1,259 at 25x25 and 532/532 at 10x10, all with Jaccard 1.
- Matched estimate rows: 3,125 at 25x25 and 500 at 10x10.
- IVW quantisation effect: median 0.0014 SE (max 0.0127) at 25x25.
- The showcase designs have no single-instrument pairs, so `e2e/wald_check.R` builds them from the strongest
  TSMR instrument per exposure. All 625/625 and 100/100 Wald ratio rows match TSMR's `mr_wald_ratio`. The match
  is exact for `fast_mr()` on the same harmonised rows. For `fast_mr_compressed()` on the stores it is within
  0.009 SE (median 2.4e-6 SE).

**Storage** (FinnGen, 10M variants, cpr; v7 -> v8 seconds; CompreSSoR is the same in v8 and v9)

| Operation | v7 | v8 |
|---|---|---|
| Write | 15.07 | 13.71 |
| Full read (8 threads) | 1.74 | 1.62 |
| Region | 1.21 | 0.88 |
| Lookup, 25 / 1,000 / 100,000 keys | 0.060 / 0.105 / 0.335 | 0.052 / 0.097 / 0.289 |

Every stream file is byte-identical to v7, both here and in all 70 trait stores. Only the manifest and
`native.index.json` differ, so the store goes from 42,654,358 to 42,653,378 bytes (43 MB).

**Multi-trait** (20 simulated GWAS, 8 threads; cpr v7 -> v8 wall seconds; TSV.gz from the full run)

| Operation | v7 | v8 | TSV.gz |
|---|---|---|---|
| Instrument selection, batch | 6.84 | 1.51 | 26.35 |
| Instrument selection, loop | 13.27 | 2.64 | 26.35 |
| Extract 100,000 keys, batch | 2.85 | 2.82 | 28.68 |
| Extract 1,000 keys, batch | 1.34 | 1.51 | 27.86 |

Row counts match TSV.gz in every cell, and the extracted values equal v7's. The batch extract is slower than v7
in 8 of its 16 cells, by 3-34%. `validation/prof_batch_extract.R` puts this on CompreSSoR #55 (`58d6eb0`), which
made same-panel key sharing in `read_sumstats_batch()` opt-in (`options(CompreSSoR.batch_share_panels = TRUE)`).
By default each of the 20 stores now resolves the same keys itself. Identity verification is not the cause: it
was never entered. With sharing switched on, d158622 beats v7 in every profiled cell
(`results_v9/batch_extract_profile.csv`):
- 1,000 keys, 1 thread: about 3.3 s, against 3.3-7.7 s on v7 and 5.3-6.4 s by default.
- 8 threads: 0.9 s against 1.3 s.

**Scaling model** (a + b*n + c*n^2 per stage, fit on 1x1 and 10x10; a model, not a measurement)
- On v7 it predicts the measured 25x25 to within -5.6% to +0.4% (C8 -3.5%).
- On v8 C8 it over-predicts 25x25 by +15.7%.
- On v9 it under-predicts 25x25 by 26% (C8 15.4 vs 20.9 s) and 24% (C1). With the clump fix, 10x10 (7.74 s) is
  barely slower than the noisy 1x1 cell (2.99-7.95 s), so the two-point fit is unstable.
- 50x50 is under-predicted by 38-62% on every build. The 50x50 set reuses the 25x25 traits as both exposures
  and outcomes (`e2e/make_traits50b.sh`), so it is not a scaled copy of the 10x10 design.

**Caveats**
- Most v8/v9 cells shared bp1-compute053, a busy node; the storage replicates ran on bp1-compute184. v7 used
  several nodes of the same CPU model.
- The `Cascadelake` feature also matches Gold 6226 (2.7 GHz) nodes (bp1-compute066, bp1-compute124). The
  e2e, A/B and multi scripts therefore refuse any CPU other than Gold 6226R (`REQUIRE_CPU`).
