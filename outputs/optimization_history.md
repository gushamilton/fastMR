# fastMR optimization history

Fixture: 82 rows; grid: 2500 pairs; nboot: 100; seed: 20260801.

| loop | wall seconds | pairs/s | memory bytes | correctness |
|---|---:|---:|---:|:---:|
| single_pair_warm | 0.007000 | 142.857 | NA | TRUE |
| single_pair_ivw_egger | 0.000000 | NA | NA | TRUE |
| grid_point_threads_1 | 0.547000 | 4570.384 | NA | TRUE |
| grid_point_threads_4 | 0.441000 | 5668.934 | NA | TRUE |
| grid_point_threads_10 | 0.423000 | 5910.165 | NA | TRUE |
| grid_boot_threads_1 | 15.294000 | 163.463 | NA | TRUE |
| grid_boot_threads_4 | 4.345000 | 575.374 | NA | TRUE |
| grid_boot_threads_10 | 2.973000 | 840.901 | NA | TRUE |
| parquet_read | 0.001000 | 82000.000 | NA | TRUE |
| grid_correctness_threads | 2.969000 | 842.034 | NA | TRUE |
| prototype_native_reference | 1.663021 | 1503.288 | 139493376 | TRUE |
| prototype_python_reference | 95.032412 | 26.307 | NA | TRUE |
| TwoSampleMR_0.7.9_warm | 0.138000 | 7.246 | NA | TRUE |

Maximum seeded thread correctness difference: 0e+00.

Each row's hypothesis and exact command are retained in optimization_history.csv/json.
