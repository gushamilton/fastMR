# fastMR optimization history

Fixture: 82 rows; grid: 2500 pairs; nboot: 100; seed: 20260801.

| loop | wall seconds | pairs/s | memory bytes | correctness |
|---|---:|---:|---:|:---:|
| single_pair_warm | 0.008000 | 125.000 | NA | TRUE |
| single_pair_ivw_egger | 0.001000 | 1000.000 | NA | TRUE |
| grid_point_threads_1 | 1.242000 | 2012.882 | NA | TRUE |
| grid_point_threads_4 | 1.122000 | 2228.164 | NA | TRUE |
| grid_point_threads_10 | 1.117000 | 2238.138 | NA | TRUE |
| grid_boot_threads_1 | 16.601000 | 150.593 | NA | TRUE |
| grid_boot_threads_4 | 5.279000 | 473.575 | NA | TRUE |
| grid_boot_threads_10 | 3.971000 | 629.564 | NA | TRUE |
| parquet_read | 0.001000 | 82000.000 | NA | TRUE |
| grid_correctness_threads | 4.136000 | 604.449 | NA | TRUE |
| prototype_native_reference | 1.663021 | 1503.288 | 139493376 | TRUE |
| prototype_python_reference | 95.032412 | 26.307 | NA | TRUE |
| TwoSampleMR_0.7.9_warm | 0.141000 | 7.092 | NA | TRUE |

Maximum seeded thread correctness difference: 0e+00.

Each row's hypothesis and exact command are retained in optimization_history.csv/json.
