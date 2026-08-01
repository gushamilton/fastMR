# Native TwoSampleMR exact mode benchmark

IL6 fixture: 82 SNPs; nboot=100; fastMR threads=5; seed=20260818.

| scenario | method | pairs | fastMR s | TwoSampleMR s | speedup | max beta delta | max SE delta |
|---|---|---:|---:|---:|---:|---:|---:|
| balanced_50x50 | simple_mode | 2500 | 1.352 | 177.175 | 131.05x | 2.455e-16 | 7.050e-03 |
| balanced_50x50 | weighted_mode | 2500 | 1.476 | 178.290 | 120.79x | 2.520e-16 | 1.016e-03 |
| one_exposure_250_outcomes | simple_mode | 250 | 0.154 | 17.745 | 115.23x | 1.448e-16 | 5.910e-03 |
| one_exposure_250_outcomes | weighted_mode | 250 | 0.160 | 17.986 | 112.41x | 2.290e-16 | 9.642e-04 |
| one_outcome_250_exposures | simple_mode | 250 | 0.143 | 17.571 | 122.87x | 2.463e-16 | 5.372e-03 |
| one_outcome_250_exposures | weighted_mode | 250 | 0.148 | 17.510 | 118.31x | 2.515e-16 | 7.125e-04 |
| ten_exposures_100_outcomes | simple_mode | 1000 | 0.541 | 69.751 | 128.93x | 2.663e-16 | 6.395e-03 |
| ten_exposures_100_outcomes | weighted_mode | 1000 | 0.540 | 71.459 | 132.33x | 2.177e-16 | 1.035e-03 |
| hundred_exposures_ten_outcomes | simple_mode | 1000 | 0.596 | 70.148 | 117.70x | 2.472e-16 | 6.086e-03 |
| hundred_exposures_ten_outcomes | weighted_mode | 1000 | 0.535 | 69.696 | 130.27x | 2.520e-16 | 7.662e-04 |
