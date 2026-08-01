# Native TwoSampleMR IVW benchmark

IL6 fixture: 82 SNPs; nboot=100; fastMR threads=5; seed=20260802; fastMR timings averaged over 20-100 repeats.

| scenario | pairs | fastMR s | TwoSampleMR s | speedup | max beta delta | max SE delta |
|---|---:|---:|---:|---:|---:|---:|
| balanced_50x50 | 2500 | 0.001 | 13.948 | 14682.11x | 3.253e-18 | 3.036e-18 |
| one_exposure_250_outcomes | 250 | 0.001 | 1.471 | 1987.84x | 1.952e-18 | 2.168e-18 |
| one_outcome_250_exposures | 250 | 0.001 | 1.419 | 1867.11x | 1.735e-18 | 1.301e-18 |
| ten_exposures_100_outcomes | 1000 | 0.001 | 5.616 | 7488.00x | 2.819e-18 | 3.036e-18 |
| hundred_exposures_ten_outcomes | 1000 | 0.001 | 5.692 | 7589.33x | 2.602e-18 | 2.168e-18 |
