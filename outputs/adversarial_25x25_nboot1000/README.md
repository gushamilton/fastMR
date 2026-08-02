# Adversarial 25 x 25 fastMR validation

Simulation: 25 exposures x 25 outcomes, 82 SNPs, five main methods, nboot=1000.
Rows: 3125 fastMR and 3125 native TwoSampleMR; expected 3125.
Timings: fastMR 10 threads 3.630s; fastMR 1 thread 20.177s; native TwoSampleMR 824.585s.
Speedup: 227.16x versus native TwoSampleMR.
Seeded thread reproducibility maximum delta: 0.000e+00.
Maximum native point-estimate beta delta: 4.526e-12.
Edge checks: 7/7 passed.

Bootstrap SE and p-value differences are not treated as point-estimate failures: fastMR and TwoSampleMR use independent bootstrap implementations/streams. The method_summary.csv file reports their per-method distributions.

The input includes deliberately varied signs, scales, heavy-tailed effects, and heteroskedastic standard errors. Separate edge checks cover 1/2/3/7-SNP inputs, rejected non-finite and mismatched matrices, and flipped alleles through fast_harmonise_data(action = 2).
