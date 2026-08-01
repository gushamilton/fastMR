# Five-round fastMR optimization cycle

Workload: IL6/CRP fixture, 50x50 exposure/outcome grid, five methods, nboot=100, ten requested threads.
The rejected round is retained to show that the slower scheduler was not kept.

| stage | change | seconds | delta vs baseline | accepted |
|---|---|---:|---:|:---:|
| baseline | Pre-cycle O2 build, nested mixed-grid results, original mode scan | 1.544 | 0.00% | TRUE |
| round_1 | Compiler -O3 (rejected: R CMD check portability warning) | 1.523 | -1.36% | FALSE |
| round_2 | Static fallback thread chunks (rejected) | 1.728 | 11.92% | FALSE |
| round_3 | Flat mixed-grid Result storage | 1.538 | -0.39% | TRUE |
| round_4 | Linearized exact-mode interpolation positions | 1.534 | -0.65% | TRUE |
| round_5 | Fused simple/weighted exact-mode max scans | 1.498 | -2.98% | TRUE |
| final_validation | Final retained implementation | 1.503 | -2.66% | TRUE |

Final validation median: 1.503 seconds (-2.66% versus baseline).
