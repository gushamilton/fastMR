# fastMR implementation and scaling iteration

Prior native TwoSampleMR 50x50 reference: 338.919 seconds at nboot=100.

| component | nboot | threads | median seconds | minimum seconds | rows/pairs |
|---|---:|---:|---:|---:|---:|
| MR grid | 0 | 1 | 0.4780 | 0.4780 | 2500 |
| MR grid | 0 | 2 | 0.4350 | 0.4340 | 2500 |
| MR grid | 0 | 5 | 0.4510 | 0.4250 | 2500 |
| MR grid | 0 | 10 | 0.4420 | 0.4390 | 2500 |
| MR grid | 100 | 1 | 8.4260 | 8.4250 | 2500 |
| MR grid | 100 | 2 | 4.5320 | 4.4800 | 2500 |
| MR grid | 100 | 5 | 2.5160 | 2.5140 | 2500 |
| MR grid | 100 | 10 | 2.0170 | 1.9660 | 2500 |
| local LD-matrix clump | NA | NA | 0.0360 | 0.0360 | 1500 |
