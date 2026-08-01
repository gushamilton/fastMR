# fastMR implementation and scaling iteration

Prior native TwoSampleMR 50x50 reference: 338.919 seconds at nboot=100.

| component | nboot | threads | median seconds | minimum seconds | rows/pairs |
|---|---:|---:|---:|---:|---:|
| MR grid | 0 | 1 | 1.4700 | 1.4120 | 2500 |
| MR grid | 0 | 2 | 1.2760 | 1.2540 | 2500 |
| MR grid | 0 | 5 | 1.2080 | 1.1960 | 2500 |
| MR grid | 0 | 10 | 1.1970 | 1.1820 | 2500 |
| MR grid | 100 | 1 | 17.1060 | 17.0260 | 2500 |
| MR grid | 100 | 2 | 9.3080 | 9.2820 | 2500 |
| MR grid | 100 | 5 | 4.9080 | 4.8770 | 2500 |
| MR grid | 100 | 10 | 3.8870 | 3.8840 | 2500 |
| local LD-matrix clump | NA | NA | 0.0350 | 0.0350 | 1500 |
