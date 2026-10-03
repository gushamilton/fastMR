# Showcase benchmark: CompreSSoR + fastMR versus the standard MR toolchain

Audience: a newcomer seeing these packages for the first time. They do not care
about speedups over our own earlier commits. They care how long a real MR study
takes, and how much disk it needs, with the tools they already use, compared
with CompreSSoR + fastMR, and whether the answers agree.

All compute runs on BluePebble (Slurm account sscm013902, under
/user/work/fh6520/showcase/). Final package SHAs come from `main` of both repos
after the perf PRs merge. Every timed cell is a fresh R process under
`/usr/bin/time -v`, with 3 replicates (best/median of three) on distinct nodes. Record the
host, CPU model, cores allocated, package versions/SHAs, wall time and peak
RSS. Report medians with min–max ranges.

## 0. Data

### Real inputs (storage component)
- **R1:** FinnGen R2 ANTIDEPRESSANTS, normalised to 10,000,000 biallelic SNPs.
  Already on BP: `CompreSSoR-bp-thread-test/external/benchmark-10m/finngen_10m_snps.tsv.gz`
  plus its metadata JSON.
- **R2:** UKB-PPP C3 protein, about 23.7M variants. Already on BP under
  `PPPMR/results/C3_.../all_variants.cpr`, plus its prepared source if still
  present; otherwise skip R2.

### Simulated studies (end-to-end and scaling)
- **Variant panel:** 1000 Genomes EUR GRCh38 biallelic SNPs with MAF ≥ 0.01,
  from `1kg_ref/grch38` (pgen) or `grch38_superpops/EUR.*`. That is about
  8–9M SNPs, a realistic imputed-GWAS row count. Use the same panel and the
  same 1000G EUR genotypes as the clumping LD reference, so clumping is
  realistic.
- **Trait generator:** this must be deterministic and seeded, and written once
  as the shared source of truth.
  - Each exposure gets K ~ U{20, 300} causal SNPs with N = 50k–400k.
  - Marginal Z within ±1 Mb of each causal SNP is induced through real LD (r
    from the 1000G EUR genotypes). Elsewhere, Z is iid N(0,1).
  - Outcomes share causal variants with the exposures (true θ drawn per pair,
    including θ = 0 nulls) plus balanced pleiotropy.
  - Columns are GWAS-SSF style: chromosome, base_pair_location,
    effect_allele, other_allele, beta, standard_error, effect_allele_frequency,
    p_value, rsid, n.
  - Document the generator in one paragraph (LD-induced signal regions, iid
    null background).
- **Trait counts:** up to 1000 exposures (proteome-like) and 100 outcomes.
  Materialise only what each format needs. TSV.gz for 1000 traits would be
  about 130 GB, so the standard pipeline only gets the sizes it can finish
  (see §2).

## 1. Storage component (CompreSSoR versus formats people actually use)

### Formats
Each format stores the same rows and row order:

| # | Format | How it is written |
|---|---|---|
| 1 | TSV.gz | GWAS Catalog harmonised / GWAS-SSF, gzip -6 |
| 2 | bgzip TSV + tabix index | |
| 3 | GWAS-VCF + .tbi | OpenGWAS/IEU format (`bcftools/1.19`, `htslib/1.19.1` modules) |
| 4 | Parquet | arrow, zstd, row-group sorted by position, no dictionary on numerics |
| 5 | CompreSSoR .cpr | defaults, including the exact p-order domain |

DuckDB is not installed. Add it only if it's a cheap install into a private
lib on /user/work, and label it as optional.

### Metrics per format, on R1 (and R2 if available)
- On-disk bytes and bytes per variant; write/convert time.
- Full read into an R data.frame with the standard logical columns.
- **Regional query:** a cis window of ±500 kb around 100 fixed loci. Report the
  total and the per-query time.
- **Variant lookup:** 25, 1,000 and 100,000 random keys (chr:pos:ref:alt),
  including file open.
- Peak RSS for each operation; cold first query versus warm repeat, reported
  separately.
- **Fidelity:** CompreSSoR is lossy, so report max/99.9th-percentile
  abs/relative error of beta, SE, Z, EAF and −log10 p versus the source. Also
  report how many p ≤ 5e-8 hits and p ≤ 1e-5 / 0.01 rows change membership.
  The other formats should be exact; verify that.
- **Projected footprint:** bytes for 1,000 / 3,000 traits (UKB-PPP scale) for
  each format, as a computed column.

### Plots
- Size versus full-read time, as a Pareto plot.
- Grouped bars for region and lookup latency.
- Projected footprint.

## 2. End-to-end MR study (small, illustrative)

Scope decision (user): the weight of evidence comes from the simulated component
benchmarks (sections 1 and 3), each with a time cutoff and censored points marked.
The end-to-end run is only 1×1 and 10×10, on synthetic data, to show the stage
breakdown and the agreement between arms.

### Pipeline
For E exposures × O outcomes, from files on disk to results on disk:

1. Instrument selection at p < 5e-8 from each exposure GWAS.
2. LD clumping on the 1000G EUR reference (r² < 0.001, 10 Mb), using local
   plink.
3. Extract those SNPs from every outcome.
4. Harmonise.
5. MR with TwoSampleMR's default method set (IVW, MR-Egger, weighted median,
   simple mode, weighted mode); nboot as TwoSampleMR's default.
6. Steiger filtering/directionality.
7. Write a results table.

### Standard arm
This should be the best realistic usage, not a straw man.
- `data.table::fread` the TSV.gz (nThread = cores), filter p, then
  `TwoSampleMR::format_data`.
- `ieugwasr::ld_clump(plink_bin = genetics.binaRies::get_plink_binary(), bfile = EUR)`,
  or `TwoSampleMR::clump_data` pointed at local plink.
- `TwoSampleMR::read_outcome_data(snps = ...)` per outcome.
- `harmonise_data`, `mr`, `steiger_filtering`, `fwrite`.
- **Second standard variant (OpenGWAS-style):** GWAS-VCF + tabix, with outcome
  SNPs extracted via `bcftools view -R` or `Rsamtools::scanTabix`. This is the
  strongest existing indexed alternative.

### New arm
- `.cpr` stores, `fastMR::fast_clump_compressed(partition = "graph")` and
  `fast_mr_compressed()`, or the documented recommended path.
- Use the same thresholds and the same reference.
- Use the same default method set and nboot.
- Run fastMR harmonisation/Steiger as the API exposes them.
- Write Parquet output.

### Sizes
1×1 and 10×10 only.
- Cap each component run at 20–25 min (as in the existing harness): censored
points are shown as "> cap", with no extrapolation. The standard arm gets 8 cores: it may use data.table/plink threads.
- The new arm is reported at 1 and 8 threads.

### Report
- Wall time per stage (select/clump/extract/harmonise/MR/Steiger/write), so the
  bottleneck in each arm is visible.
- Total wall time and peak RSS, plus total disk footprint of the inputs.
- **Agreement:** instrument sets identical (Jaccard); harmonised SNP sets;
  max |Δb| for IVW and Egger; median-method |Δb|; bootstrap-SE
  difference distributions (RNG differs, which is expected). Any mismatch
  must be explained.

## 3. MR compute component (fastMR versus TwoSampleMR, harmonised data in memory)

Re-use `/user/work/fh6520/fastmr-vs-tsmr/worker.R` and `run_task2.sh`. Their
scenarios are already sound: single pair 10/100/1000 SNPs; many pairs
100…100k for IVW and default methods; harmonise 1e4–1e6; Steiger 1k/20k;
directionality; heterogeneity; pleiotropy; clumping (TwoSampleMR/plink1.9
versus fastMR graph/plink2).

Point them at the final `main` build and drop the `fastMR_baseline` column.
Add fastMR at 8 threads.

### Plots
Time versus number of pairs (log–log), one line per tool and method set, with
censored points marked.

## 4. Deliverables
- Raw CSVs, per-replicate, plus aggregated summaries.
- Scripts in `fastMR/benchmarks/showcase/` (generator, writers, arms,
  slurm arrays, aggregate, plots), committed on a branch and PR'd.
- Figures (SVG/PNG) for slides.
- Update the published results page: TSMR section, storage, end-to-end.
- Clean up regenerable data on BP when done. Keep results, logs and scripts.

## 5. Job layout and time budget (user: few jobs, 3 replicates, think about timings)

- Few jobs: one Slurm array per suite (storage, MR compute, end-to-end), 3 tasks each
  (one per replicate, on different nodes; `--exclude` earlier nodes if needed), so about
  9 tasks in total. There are no per-cell jobs. All three suites can run concurrently.
- Each task runs its cells sequentially in fresh R processes with a per-cell cap.
  A size sweep stops for that tool once a smaller size hits the cap, so the slow
  tool never burns hours on larger sizes.
- Per-cell cap: 15 min for MR-compute cells, 10 min for storage operations. The cap
  is set so that each task's worst case stays under about 4 h (the short partition is
  fine, and nothing runs overnight).
- Before submitting, write a budget table from the earlier TwoSampleMR timings on
  BP (TwoSampleMR default methods ≈ 5.7 s/pair; IVW ≈ 12 ms/pair; harmonise 1e5
  rows = 183 s) and choose sizes so the expected censored count is small and
  deliberate. Example: TwoSampleMR default methods 10/100 pairs run, 1000 censored.
- Timing repetition within a cell: fast calls are repeated to ≥ 2 s total and the
  median is used; slow calls (> 30 s) are timed once per replicate.
- Do a dry run of one replicate at the smallest sizes first (≤ 15 min) to check
  the plumbing before the arrays go in.

## 6. Revisions after the Opus review (these supersede the earlier text)

1. Use one cap everywhere: 25 min per cell, with no extrapolation. A cell
   predicted to exceed the cap from the previous size is marked "> cap
   (predicted)" and skipped.
2. **TwoSampleMR at equal cores:** add a TSMR × 8 arm (`parallel::mclapply`
   over pairs) wherever fastMR runs with 8 threads. Headline both 1 versus 1
   core and 8 versus 8 cores.
3. Install the current TwoSampleMR from MRCIEU GitHub (0.6.x) into a private
   lib on /user/work, and pin ieugwasr. Record the versions and SHAs.
4. Report agreement in two parts:
   - (a) fastMR versus TSMR on the same exact TSV input, which should be
     ≈ 1e-12;
   - (b) fastMR on `.cpr` versus fastMR on TSV, in SE units, as the
     quantisation effect.
   Never mix the two.
5. Make instrument sets identical by construction: `candidate_source =
   "pvalue_flag"` plus `pvalue_order = "require_exact"`. Report
   reconstructed-p membership flips at 5e-8, 1e-5 and 0.01 only as a storage
   fidelity metric.
6. Harmonisation: synthetic data is on the forward strand, and TSMR uses
   `action = 1`. State this explicitly.
7. Build one LD reference: EUR, MAF ≥ 0.01, biallelic SNPs, about 9M, from
   `1kg_ref/grch38_superpops/EUR`. Write it once as both bed (plink 1.9) and
   pgen (plink2), from the same samples and variants. Use it for every clump
   cell; `EUR_sub` (903k SNPs) is no longer used.
8. Drop R2 (C3): it is 14.4M rows, written with the legacy python backend,
   and its source is not on BP. Use the real R1 FinnGen 10M for storage, and
   base footprint projections on R1's bytes per variant.
9. **Thread control:** set `setDTthreads`, `arrow::set_cpu_count`,
   OMP/OPENBLAS=1, plink `--threads`, and fastMR `threads × io_threads ≤
   allocation`. Assert that `SLURM_CPUS_PER_TASK` matches the threads label.
   Record user+sys CPU.
10. **Storage:**
    - Add a quantised-Parquet control (float32 beta/se, zstd).
    - Write Parquet sorted, with about 1M-row row groups, and query it with
      `open_dataset` filter pushdown.
    - Read TSV.gz with the best of `fread(cmd = "pigz -dc")` and
      `arrow::read_csv_arrow`.
    - Report CompreSSoR and Parquet at 1 thread and 8 threads.
    - VCF via bcftools; for 100k lookups use the best of `-R` and `-T`.
      gwasvcf is not installed, so skip it.
    - Label all timings "warm"; there is no cold/warm split.
11. **Bootstrap SEs:** use a TSMR seed A versus seed B control as the null
    distribution for fastMR-versus-TSMR SE differences.
12. Vary K in many-pairs default runs: K = 10 and 50.
13. **Cuts:** single-pair harmonise at 1e6 (replaced by many-pairs harmonise,
    1k pairs × 30 SNPs), DuckDB, and the cold/warm split.
14. **Clump:** TSMR (plink 1.9) per exposure, an optional plink2 `--clump` per
    exposure, and fastMR graph. Use 10 and 100 exposures, plus 1000 for fastMR
    only. Verify instrument identity.
15. **Simulated traits:**
    - Null Z comes from one `plink2 --glm` pass over the EUR panel with random
      phenotype columns, which gives LD-correlated null Z.
    - Add signal √N·Rβ within ±1 Mb windows.
    - SE = 1/√(2N·f(1−f)) × (1 + small noise); β = Z·SE.
    - Report the simulated compression ratio against R1 as a sanity check.
16. **End-to-end:** include the one-off conversion cost (prepare + encode) as an
    amortised line. The arms are TSV.gz + TSMR (×1 and ×8), VCF + tabix +
    TSMR, and `.cpr` + fastMR.
17. **Replicates:** 3 per suite (user decision), with arms paired within the
    same job on the same node, in randomised order. Prefer Gold 6226 nodes.
    Expected wall time per task is ≤ 2 h. All three arrays run concurrently,
    for about 9 tasks in total.
