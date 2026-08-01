# fastMR optimization history

Fixture: 82 rows; grid: 2500 pairs; nboot: 100; seed: 20260801.

| loop | wall seconds | pairs/s | memory bytes | correctness |
|---|---:|---:|---:|:---:|
| single_pair_warm | 0.008000 | 125.000 | NA | TRUE |
| single_pair_ivw_egger | 0.001000 | 1000.000 | NA | TRUE |
| grid_point_threads_1 | 2.367000 | 1056.189 | NA | TRUE |
| grid_point_threads_4 | 2.349000 | 1064.283 | NA | TRUE |
| grid_point_threads_10 | 2.287000 | 1093.135 | NA | TRUE |
| grid_boot_threads_1 | 16.725000 | 149.477 | NA | TRUE |
| grid_boot_threads_4 | 6.543000 | 382.088 | NA | TRUE |
| grid_boot_threads_10 | 5.089000 | 491.256 | NA | TRUE |
| parquet_read | 0.001000 | 82000.000 | NA | TRUE |
| grid_correctness_threads | 5.103000 | 489.908 | NA | TRUE |
| prototype_native_reference | 1.663021 | 1503.288 | 139493376 | TRUE |
| prototype_python_reference | 95.032412 | 26.307 | NA | TRUE |
| TwoSampleMR_0.7.9_warm | 0.139000 | 7.194 | NA | TRUE |

Maximum seeded thread correctness difference: 0e+00.

Each row's hypothesis and exact command are retained in optimization_history.csv/json.
