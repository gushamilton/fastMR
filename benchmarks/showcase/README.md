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

## v8 rerun (fastMR `543cae2`, CompreSSoR `d158622`) against v7

All cells on Intel Xeon Gold 6226R with 8 Slurm CPUs; medians of 3 replicates in fresh processes (min-max in
`results_v8/`). TwoSampleMR and GWAS-VCF arms are the v7 measurements (they do not use this code). Tables are
written by `e2e/compare_v8.R`; CSVs are in `results_v8/`.

End to end, wall seconds (fastMR arm C8 = 8 CPUs):

| Study | v7 C8 | v8 C8 | change | TSMR, MR serial (1 core) | TSMR, MR over 8 CPUs | GWAS-VCF + TSMR |
|---|---|---|---|---|---|---|
| 10x10 | 10.87 | 16.44 | +51% | 902.5 | 518.6 | 649.6 |
| 25x25 | 23.17 | 28.53 | +23% | 4,346 | 1,613 | 3,562.5 |
| 50x50 | 109.70 | 108.42 | -1% | not run | not run | not run |

Speedups of v8 C8 (8 CPUs): 25x25 152x vs TSMR 1 core and 57x vs TSMR 8 CPUs (v7: 188x, 70x); 10x10 55x and
32x (v7: 83x, 48x). Per stage at 25x25 C8 (v7 -> v8): clump 7.46 -> 15.93 s, extract 4.70 -> 4.26 s,
MR 6.44 -> 3.90 s, Steiger 0.16 -> 0.13 s. In a same-node interleaved A/B (`e2e/e2e_ab.sbatch`, 6 cells per
build and size) 25x25 C8 took 21.85 s on v7 and 25.95 s on v8 (clump 6.97 -> 14.72 s, MR 5.84 -> 3.90 s).
`validation/prof_clump_ab.R` puts the clump regression in fastMR's batched-candidate guard: it reads every
exposure store's flagged row ids with `CompreSSoR::read_pvalue_flag(as = "row_ids")`, one store at a time,
before the batch read (about 7.3 s of 9.8 s profiled at 25 exposures). CompreSSoR's identity verification
costs 0.2-0.7 s for the 50 stores (`validation/profile_verify.R`).

Agreement (v7 = v8): instruments identical to TSMR (1,259/1,259 at 25x25, 532/532 at 10x10, all Jaccard 1);
3,125 and 500 matched estimate rows; IVW quantisation effect median 0.0014 SE (max 0.0127) at 25x25. The
showcase designs have no single-instrument pairs, so `e2e/wald_check.R` builds them (strongest TSMR instrument
per exposure): 625/625 and 100/100 Wald ratio rows match TSMR's `mr_wald_ratio`, exactly for `fast_mr()` on the
same harmonised rows and within 0.009 SE (median 2.4e-6 SE) for `fast_mr_compressed()` on the stores.

Storage (FinnGen, 10M variants, cpr, v7 -> v8 seconds): write 15.07 -> 13.71; full read (8 threads)
1.74 -> 1.62; region 1.21 -> 0.88; lookups 25 / 1,000 / 100,000 keys 0.060 / 0.105 / 0.335 -> 0.052 / 0.097 /
0.289. Every stream file is byte-identical to v7 (here and in all 70 trait stores); only the manifest and
`native.index.json` differ, so the store is 42,654,358 -> 42,653,378 bytes (43 MB).

Multi-trait (20 simulated GWAS, 8 threads, cpr v7 -> v8 wall seconds; TSV.gz from the full run): instrument
selection, batch 6.84 -> 1.51 and loop 13.27 -> 2.64 (TSV.gz 26.35); extract 100,000 keys, batch 2.85 -> 2.82
(TSV.gz 28.68); extract 1,000 keys, batch 1.34 -> 1.51 (TSV.gz 27.86). Row counts match TSV.gz in every cell.
Instrument selection and the loop extract are faster in every cell; the batch extract is slower than v7 in 8
of its 16 cells, by 3-34% (worst: 1,000 keys from 20 stores on 1 thread, 3.64 -> 4.88 s).

Scaling model (a + b*n + c*n^2 per stage, fit on 1x1 and 10x10; a model, not a measurement), prediction vs
measured 25x25: v8 C8 33.0 vs 28.5 s (+15.7%), C1 63.9 vs 68.0 s (-6.1%); the same script on v7 gives C8 22.4 vs
23.2 s (-3.5%), C1 59.1 vs 62.6 s (-5.6%), TSMR 4,246 vs 4,346 s (-2.3%), GWAS-VCF 3,575 vs 3,562.5 s (+0.4%).
At 50x50 it under-predicts both builds by 38-53%: the 50x50 set reuses the 25x25 traits as both exposures
and outcomes (`e2e/make_traits50b.sh`), so it is not a scaled copy of the 10x10 design.

Caveats: the storage replicates all landed on one node (bp1-compute184), the v8 fastMR end-to-end and
multi-trait replicates on bp1-compute053; v7 used several nodes of the same model. The `Cascadelake` feature
also matches Gold 6226 (2.7 GHz) nodes, so `e2e/e2e_v8.sbatch`, `e2e/e2e_ab.sbatch` and `multi/multi.sbatch`
(`REQUIRE_CPU`) now refuse any CPU other than Gold 6226R.
