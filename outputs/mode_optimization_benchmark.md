# Exact mode-kernel optimisation benchmark

IL6 fixture: 82 SNPs; 50x50 grid; nboot=100; threads=10; five warm repeats.
The pre-optimisation baseline is the measured implementation before reusable FFT workspaces and plans.
The estimator and 512-point density grid are unchanged.

| method | median fastMR s | min fastMR s | pre-optimisation s | speedup |
|---|---:|---:|---:|---:|
| simple_mode | 1.760 | 1.746 | 1.886 | 1.07x |
| weighted_mode | 1.770 | 1.740 | 1.910 | 1.08x |
