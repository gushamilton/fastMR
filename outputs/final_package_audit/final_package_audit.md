# Final fastMR package audit

Local BMI -> CRP: 104 BMI rows, 104 matched local CRP rows, 5 repeats.
Median BMI -> CRP pipeline: fastMR 0.005000s; TwoSampleMR 0.127000s; speedup 25.40x.
BMI -> CRP harmonisation parity was exact: TRUE.

Different simulation 50x50 grid: 400 SNPs, five methods, nboot=100, fastMR threads=10.
Grid: fastMR 2.030000s; TwoSampleMR 391.764000s; speedup 192.99x.
Grid maximum absolute deltas: beta 6.924e-11, SE 3.623e+04, p-value 1.282e-01.
Grid beta parity is exact to numerical precision for IVW/Egger (max beta delta 1.554e-15); all methods have max beta delta 6.924e-11.
Bootstrap SE/p-value variation by method (median SE delta; 95th percentile SE delta; median p-value delta): egger: median SE 4.163e-17, p95 SE 1.665e-16, median p 4.718e-16; ivw: median SE 2.082e-17, p95 SE 6.245e-17, median p 2.220e-16; simple_mode: median SE 4.339e+00, p95 SE 6.592e+02, median p 5.417e-04; weighted_median: median SE 3.098e-03, p95 SE 9.258e-03, median p 1.452e-02; weighted_mode: median SE 4.486e+00, p95 SE 7.355e+02, median p 5.591e-04.
The large raw SE maximum is a stress-test property of the heavy-tailed simulation: near-zero exposure draws make ratio bootstraps unstable. It is not a point-estimate disagreement; bootstrap streams are independently sampled between fastMR and TwoSampleMR.

Different simulation 1x1: fastMR 0.008000s; TwoSampleMR 0.165000s; speedup 20.63x.
1x1 maximum absolute deltas: beta 2.776e-16, SE 8.466e+00, p-value 4.618e-03.

Interpretation: the main weakness for 1x1 workloads is R/data-frame and harmonisation overhead; the large-grid path is dominated by bootstrap-heavy modes/Egger rather than IVW. Local PLINK timing remains unavailable because no PLINK executable is installed on the Mac mini; the wrapper was tested separately with a contract-validating stub.
