# Native TwoSampleMR deterministic method benchmarks

IL6 fixture: 82 SNPs; nboot=0; fastMR threads=5.

| method | scenario | pairs | fastMR s | TwoSampleMR s | speedup | max beta delta | max p-value delta |
|---|---|---:|---:|---:|---:|---:|---:|
| uwr | balanced_50x50 | 2500 | 1.307 | 13.910 | 10.64x | 4.337e-19 | 0.000e+00 |
| uwr | one_exposure_250_outcomes | 250 | 0.132 | 1.397 | 10.58x | 4.337e-19 | 0.000e+00 |
| uwr | one_outcome_250_exposures | 250 | 0.133 | 1.390 | 10.45x | 4.337e-19 | 0.000e+00 |
| uwr | ten_exposures_100_outcomes | 1000 | 0.532 | 5.669 | 10.66x | 4.337e-19 | 0.000e+00 |
| uwr | hundred_exposures_ten_outcomes | 1000 | 0.530 | 5.627 | 10.62x | 4.337e-19 | 0.000e+00 |
| sign | balanced_50x50 | 2500 | 1.240 | 13.708 | 11.05x | 0.000e+00 | 0.000e+00 |
| sign | one_exposure_250_outcomes | 250 | 0.129 | 1.457 | 11.29x | 0.000e+00 | 0.000e+00 |
| sign | one_outcome_250_exposures | 250 | 0.138 | 1.495 | 10.83x | 0.000e+00 | 0.000e+00 |
| sign | ten_exposures_100_outcomes | 1000 | 0.509 | 5.522 | 10.85x | 0.000e+00 | 0.000e+00 |
| sign | hundred_exposures_ten_outcomes | 1000 | 0.508 | 5.513 | 10.85x | 0.000e+00 | 0.000e+00 |
