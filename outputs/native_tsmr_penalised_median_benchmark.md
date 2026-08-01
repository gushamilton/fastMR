# Native TwoSampleMR penalised weighted median benchmark

IL6 fixture: 82 SNPs; nboot=100; fastMR threads=5; seed=20260904.
The native comparator is the installed TwoSampleMR R workflow.

| scenario | pairs | fastMR s | TwoSampleMR s | speedup | max beta delta | median SE delta |
|---|---:|---:|---:|---:|---:|---:|
| balanced_50x50 | 2500 | 0.143 | 20.544 | 143.66x | 4.337e-19 | 1.603e-04 |
| one_exposure_250_outcomes | 250 | 0.017 | 2.083 | 122.53x | 4.337e-19 | 1.662e-04 |
| one_outcome_250_exposures | 250 | 0.016 | 2.068 | 129.25x | 4.337e-19 | 1.750e-04 |
| ten_exposures_100_outcomes | 1000 | 0.054 | 8.261 | 152.98x | 6.505e-19 | 1.448e-04 |
| hundred_exposures_ten_outcomes | 1000 | 0.052 | 8.294 | 159.50x | 6.505e-19 | 1.642e-04 |
