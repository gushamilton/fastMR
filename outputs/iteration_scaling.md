# fastMR implementation and scaling iteration

Prior native TwoSampleMR 50x50 reference: 338.919 seconds at nboot=100.

| component | nboot | threads | median seconds | minimum seconds | rows/pairs |
|---|---:|---:|---:|---:|---:|
| MR grid | 0 | 1 | 2.4280 | 2.4120 | 2500 |
| MR grid | 0 | 2 | 2.4020 | 2.3820 | 2500 |
| MR grid | 0 | 5 | 2.3690 | 2.3320 | 2500 |
| MR grid | 0 | 10 | 2.3330 | 2.3200 | 2500 |
| MR grid | 100 | 1 | 18.8970 | 18.4480 | 2500 |
| MR grid | 100 | 2 | 10.8400 | 10.7540 | 2500 |
| MR grid | 100 | 5 | 6.7610 | 6.4920 | 2500 |
| MR grid | 100 | 10 | 5.2970 | 5.2850 | 2500 |
| local LD-matrix clump | NA | NA | 0.0380 | 0.0350 | 1500 |
