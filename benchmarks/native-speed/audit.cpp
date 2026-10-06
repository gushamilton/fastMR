// Audit of the hull path and the recurrence-kernel error bound.
#include "fastmr_copy.h"
// Max relative error of the recurrence kernel vs mode_kernel() over all
// distances used, in units of u = 2^-53, for each bandwidth ratio a.
// [[Rcpp::export]]
Rcpp::NumericVector recurrence_error_ulps(Rcpp::NumericVector deltas, Rcpp::NumericVector bandwidths) {
  Rcpp::NumericVector out(deltas.size());
  const double u = DBL_EPSILON / 2.0;
  for (R_xlen_t i = 0; i < deltas.size(); ++i) {
    const double delta = deltas[i], bandwidth = bandwidths[i];
    const double a = delta / bandwidth;
    const double radius_cells = kModeKernelCutoff * bandwidth / delta;
    const int radius = radius_cells < 510.0 ? static_cast<int>(radius_cells) + 1 : 511;
    const double aa = a * a, r = std::exp(-aa);
    double worst = 0.0;
    for (int d0 = 0; d0 <= radius; d0 += kModeRecurrenceBlock) {
      double value = mode_kernel(d0, delta, bandwidth);
      double q = std::exp(-0.5 * aa * (2.0 * static_cast<double>(d0) + 1.0));
      const int end = std::min(radius, d0 + kModeRecurrenceBlock - 1);
      for (int d = d0; d <= end; ++d) {
        const double exact = mode_kernel(d, delta, bandwidth);
        worst = std::max(worst, std::abs(value / exact - 1.0) / u);
        value *= q; q *= r;
      }
    }
    out[i] = worst;
  }
  return out;
}
// For each draw (column): run the hull path; when it certifies, check its
// index against the FFT path and record max |d_hull - d_fft| / max(y) over the
// scanned points. Returns counts and the worst relative discrepancy.
// [[Rcpp::export]]
Rcpp::NumericVector hull_audit(Rcpp::NumericMatrix ratios, Rcpp::NumericMatrix weights, double phi) {
  const int k = ratios.nrow(), draws = ratios.ncol(), n = MODE_GRID_SIZE;
  long certified = 0, mismatches = 0, fallbacks = 0; double worst = 0.0;
  for (int d = 0; d < draws; ++d) {
    const double* v = &ratios(0, d); const double* w = &weights(0, d);
    ModeGrid g = mode_grid(v, k, phi, mode_workspace().scratch);
    int idx = -1;
    const double* wp[1] = {w};
    mode_workspace().hull_density[0].assign(n, 0.0);
    const bool ok = mode_index_hull(v, wp, 1, k, g, &idx);
    if (!ok) { ++fallbacks; continue; }
    ++certified;
    std::vector<double> yh(mode_workspace().hull_density[0]);
    const int f = mode_index_fft(v, w, k, g);
    if (f != idx) ++mismatches;
    const std::vector<std::complex<double>>& yb = mode_workspace().binned;
    // compare over cells touched by the hull (non-zero hull density)
    double maxy = 0.0, diff = 0.0;
    for (int j = 1; j < n - 1; ++j) if (yh[j] != 0.0) { maxy = std::max(maxy, yh[j]); diff = std::max(diff, std::abs(yh[j] - yb[j].real())); }
    if (maxy > 0) worst = std::max(worst, diff / maxy);
  }
  return Rcpp::NumericVector::create(Rcpp::_["certified"] = certified, Rcpp::_["fallbacks"] = fallbacks,
    Rcpp::_["mismatches"] = mismatches, Rcpp::_["worst_rel_discrepancy"] = worst);
}
