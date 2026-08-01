# Native TwoSampleMR method benchmarks

IL6 fixture: 82 SNPs; nboot=100; fastMR threads=5; seed=20260802.
The native comparator is the installed TwoSampleMR R workflow, timed one method at a time.
Bootstrap SE differences are Monte Carlo differences when comparing a shared fastMR grid layout with native per-pair streams.

| scenario | method | pairs | fastMR s | TwoSampleMR s | speedup | max beta delta | median SE delta |
|---|---|---:|---:|---:|---:|---:|---:|
| balanced_50x50 | ivw | 2500 | 1.298 | 13.661 | 10.52x | 3.469e-18 | 4.337e-19 |
| balanced_50x50 | egger | 2500 | 1.436 | 14.118 | 9.83x | 8.674e-18 | 4.337e-19 |
| balanced_50x50 | weighted_median | 2500 | 1.369 | 16.840 | 12.30x | 4.337e-19 | 1.108e-04 |
| balanced_50x50 | simple_mode | 2500 | 3.598 | 169.190 | 47.02x | 2.455e-16 | 2.147e-03 |
| balanced_50x50 | weighted_mode | 2500 | 3.640 | 170.680 | 46.89x | 2.520e-16 | 2.310e-04 |
| one_exposure_250_outcomes | ivw | 250 | 0.136 | 1.437 | 10.57x | 2.819e-18 | 8.674e-19 |
| one_exposure_250_outcomes | egger | 250 | 0.153 | 1.460 | 9.54x | 7.373e-18 | 1.735e-18 |
| one_exposure_250_outcomes | weighted_median | 250 | 0.143 | 1.724 | 12.06x | 4.337e-19 | 1.041e-04 |
| one_exposure_250_outcomes | simple_mode | 250 | 0.379 | 17.101 | 45.12x | 1.448e-16 | 1.573e-03 |
| one_exposure_250_outcomes | weighted_mode | 250 | 0.365 | 17.014 | 46.61x | 2.290e-16 | 2.224e-04 |
| one_outcome_250_exposures | ivw | 250 | 0.136 | 1.471 | 10.82x | 3.253e-18 | 4.337e-19 |
| one_outcome_250_exposures | egger | 250 | 0.168 | 1.457 | 8.67x | 7.373e-18 | 4.337e-19 |
| one_outcome_250_exposures | weighted_median | 250 | 0.142 | 1.748 | 12.31x | 4.337e-19 | 1.098e-04 |
| one_outcome_250_exposures | simple_mode | 250 | 0.363 | 17.080 | 47.05x | 2.463e-16 | 1.866e-03 |
| one_outcome_250_exposures | weighted_mode | 250 | 0.358 | 17.058 | 47.65x | 2.515e-16 | 1.744e-04 |
| ten_exposures_100_outcomes | ivw | 1000 | 0.548 | 5.737 | 10.47x | 2.602e-18 | 4.337e-19 |
| ten_exposures_100_outcomes | egger | 1000 | 0.597 | 5.881 | 9.85x | 6.939e-18 | 8.674e-19 |
| ten_exposures_100_outcomes | weighted_median | 1000 | 0.554 | 6.909 | 12.47x | 4.337e-19 | 1.091e-04 |
| ten_exposures_100_outcomes | simple_mode | 1000 | 1.463 | 71.668 | 48.99x | 2.663e-16 | 2.135e-03 |
| ten_exposures_100_outcomes | weighted_mode | 1000 | 1.575 | 70.589 | 44.82x | 2.177e-16 | 2.177e-04 |
| hundred_exposures_ten_outcomes | ivw | 1000 | 0.585 | 5.899 | 10.08x | 2.819e-18 | 4.337e-19 |
| hundred_exposures_ten_outcomes | egger | 1000 | 0.615 | 5.995 | 9.75x | 7.806e-18 | 4.337e-19 |
| hundred_exposures_ten_outcomes | weighted_median | 1000 | 0.569 | 7.518 | 13.21x | 4.337e-19 | 1.043e-04 |
| hundred_exposures_ten_outcomes | simple_mode | 1000 | 1.453 | 69.027 | 47.51x | 2.472e-16 | 1.852e-03 |
| hundred_exposures_ten_outcomes | weighted_mode | 1000 | 1.434 | 69.238 | 48.28x | 2.520e-16 | 1.867e-04 |
