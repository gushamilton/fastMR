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

Each suite is one Slurm array with one replicate per task (3 tasks, different nodes),
cells in fresh R processes under `/usr/bin/time -v`, arms in randomised order, a
25-minute cap per cell and no extrapolation. Aggregate with
`mr/aggregate.R`, `storage/aggregate.R`, `multi/aggregate_multi.R` and
`e2e/aggregate_e2e.R` (after `e2e/agree.R` per replicate and size).

Versions used for the published numbers: fastMR 0.2.0 at `5a189f7`, CompreSSoR 0.6.0
(`d84e1aa`; `.cpr` write re-timed on `a27b32d`), TwoSampleMR 0.7.11 (`f492045`),
ieugwasr 1.2.0, PLINK 1.9 b7.7, PLINK2 2.00a6.8, bcftools/htslib 1.19, R 4.5.1.
