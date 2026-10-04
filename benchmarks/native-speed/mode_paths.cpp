// Microbenchmark: includes the package source verbatim and times each mode path per draw.
#include "fastmr_copy.h"
#include <time.h>
static double tcpu(){timespec t; clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t); return t.tv_sec+1e-9*t.tv_nsec;}
// [[Rcpp::export]]
Rcpp::NumericVector mode_paths(Rcpp::NumericMatrix ratios, Rcpp::NumericVector w, Rcpp::NumericVector w2, int reps) {
  const int k = ratios.nrow(), draws = ratios.ncol();
  const double* wp[2] = {w.begin(), w2.begin()};
  double best[5] = {1e300,1e300,1e300,1e300,1e300}; int mism = 0, fb = 0;
  for (int r = 0; r < reps; ++r) {
    for (int path = 0; path < 5; ++path) {
      mode_hull_recurrence.store(path == 2 || path == 4);
      double t0 = tcpu();
      for (int d = 0; d < draws; ++d) {
        const double* v = &ratios(0, d);
        ModeGrid g = mode_grid(v, k, 1.0, mode_workspace().scratch);
        int idx[2] = {0, 0};
        if (path == 0) { idx[0] = mode_index_fft(v, wp[0], k, g); }
        else if (path == 1 || path == 2) {
          if (!mode_index_hull(v, wp, 1, k, g, idx)) { if (r == 0 && path == 2) ++fb; if (!mode_index_direct(v, wp, 1, k, g, idx)) idx[0] = mode_index_fft(v, wp[0], k, g); }
          if (r == 0 && path == 2 && idx[0] != mode_index_fft(v, wp[0], k, g)) ++mism;
        } else if (path == 3) {
          if (!mode_index_direct(v, wp, 1, k, g, idx)) idx[0] = mode_index_fft(v, wp[0], k, g);
        } else {
          if (!mode_index_hull(v, wp, 2, k, g, idx)) { if (!mode_index_direct(v, wp, 2, k, g, idx)) idx[0] = mode_index_fft(v, wp[0], k, g); }
        }
      }
      best[path] = std::min(best[path], tcpu() - t0);
    }
  }
  // pair direct for reference
  double bp = 1e300;
  for (int r = 0; r < reps; ++r) { double t0 = tcpu();
    for (int d = 0; d < draws; ++d) { const double* v = &ratios(0, d); ModeGrid g = mode_grid(v, k, 1.0, mode_workspace().scratch); int idx[2];
      if (!mode_index_direct(v, wp, 2, k, g, idx)) idx[0] = mode_index_fft(v, wp[0], k, g); }
    bp = std::min(bp, tcpu() - t0); }
  const double per = 1e9 / draws;
  return Rcpp::NumericVector::create(Rcpp::_["fft"] = best[0]*per, Rcpp::_["direct"] = best[3]*per, Rcpp::_["hull"] = best[1]*per,
    Rcpp::_["hull_rec"] = best[2]*per, Rcpp::_["pair_direct"] = bp*per, Rcpp::_["pair_hull_rec"] = best[4]*per,
    Rcpp::_["fallbacks"] = fb, Rcpp::_["mismatch"] = mism);
}
