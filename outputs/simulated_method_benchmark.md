# Independent simulated method benchmark

Causal simulation seed=20260824; nboot=100; fastMR threads=5; fastMR repeats=50; native repeats=10.
The 400-SNP simulation uses a causal effect, heavy-tailed pleiotropic noise,
heterogeneous standard errors, and mixed-sign instrument effects. Wald ratio
is tested separately on a one-SNP simulation, as required by the method.

| scenario | method | SNPs | fastMR s | TwoSampleMR s | speedup | abs beta delta | abs SE delta | abs p-value delta |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| causal_400_snp | ivw | 400 | 0.001000 | 0.006500 | 6.50x | 1.665e-16 | 6.939e-18 | 9.124e-33 |
| causal_400_snp | ivw_fe | 400 | 0.001000 | 0.006000 | 6.00x | 1.665e-16 | 3.469e-18 | 5.223e-39 |
| causal_400_snp | ivw_mre | 400 | 0.001000 | 0.006000 | 6.00x | 1.665e-16 | 6.939e-18 | 9.124e-33 |
| causal_400_snp | egger | 400 | 0.001000 | 0.006000 | 6.00x | 2.776e-17 | 2.082e-17 | 3.151e-19 |
| causal_400_snp | egger_bootstrap | 400 | 0.002000 | 0.010000 | 5.00x | 0.000e+00 | 2.082e-17 | 0.000e+00 |
| causal_400_snp | uwr | 400 | 0.001000 | 0.006000 | 6.00x | 8.327e-17 | 4.441e-16 | 0.000e+00 |
| causal_400_snp | sign | 400 | 0.001000 | 0.006000 | 6.00x | 0.000e+00 | NA | 0.000e+00 |
| causal_400_snp | simple_median | 400 | 0.003000 | 0.010000 | 3.33x | 3.275e-15 | 4.163e-17 | 1.544e-22 |
| causal_400_snp | weighted_median | 400 | 0.003000 | 0.010000 | 3.33x | 0.000e+00 | 6.939e-18 | 1.165e-21 |
| causal_400_snp | penalised_weighted_median | 400 | 0.004000 | 0.014000 | 3.50x | 5.551e-17 | 0.000e+00 | 7.031e-23 |
| causal_400_snp | simple_mode | 400 | 0.004000 | 0.088500 | 22.12x | 0.000e+00 | 5.551e-16 | 1.176e-17 |
| causal_400_snp | weighted_mode | 400 | 0.004000 | 0.090500 | 22.62x | 0.000e+00 | 1.041e-16 | 1.059e-20 |
| causal_1_snp | wald_ratio | 1 | 0.001000 | 0.006000 | 6.00x | 0.000e+00 | 0.000e+00 | 0.000e+00 |
