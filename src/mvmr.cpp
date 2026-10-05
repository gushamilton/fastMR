// Batched multivariable IVW (MVMR) kernel.
//
// A batch holds E designs.  Design e is an n_e x p matrix of SNP effects on p
// exposures (column 1 usually the primary exposure, the rest covariate traits)
// whose rows address rows of a shared U x K outcome panel (SNPs x outcomes).
// Every (design, outcome) pair is an independent weighted least-squares fit
// through the origin with outcome-specific weights 1 / se_y^2, so the p x p
// normal equations differ by outcome.  p is small, so each pair accumulates
// X'WX and X'Wy in one pass over the design rows, solves by Cholesky and
// takes a second pass for the residual sum of squares (the two-pass form; the
// normal-equation identity y'Wy - b'X'Wy cancels badly for good fits).
//
// Optional shared-weight path: when the outcome standard errors factor as
// se[i, k] = row_se[i] * outcome_scale[k], every outcome's normal matrix is
// X' diag(1 / row_se^2) X / outcome_scale[k]^2, so one inverse per design
// serves every outcome whose design rows are all present.  The caller decides
// per outcome whether the factorisation holds (shared_outcome[k]).
//
// Results are written by index, so they do not depend on the thread count or
// the schedule.  Workers never touch the R API.

#include <Rcpp.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <thread>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

const double MV_NA = NA_REAL;

// Relative Cholesky pivot below which the design is treated as collinear
// (|R_jj| / ||W^(1/2) x_j|| < 1e-6 in QR terms).
constexpr double kPivotTolerance = 1e-12;

template <typename Job>
void mv_parallel(std::size_t jobs, int threads, Job job) {
  int thread_count = std::max(1, threads);
  if (static_cast<std::size_t>(thread_count) > jobs) {
    thread_count = static_cast<int>(std::max<std::size_t>(1, jobs));
  }
  if (thread_count == 1) {
    for (std::size_t index = 0; index < jobs; ++index) job(index);
    return;
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(thread_count)
  for (long long index = 0; index < static_cast<long long>(jobs); ++index) {
    job(static_cast<std::size_t>(index));
  }
#else
  std::atomic<std::size_t> next(0);
  std::vector<std::thread> pool;
  pool.reserve(static_cast<std::size_t>(thread_count));
  for (int worker = 0; worker < thread_count; ++worker) {
    pool.emplace_back([&]() {
      while (true) {
        const std::size_t index = next.fetch_add(1, std::memory_order_relaxed);
        if (index >= jobs) break;
        job(index);
      }
    });
  }
  for (std::thread& worker : pool) worker.join();
#endif
}

// In-place Cholesky of the lower triangle of the p x p row-major matrix `a`,
// then its inverse into `inv` (full, symmetric).  False when a pivot is not
// positive or falls below kPivotTolerance of its original diagonal.
// `scratch` must hold p * p + p doubles; no allocation per call.
bool spd_inverse(std::vector<double>& a, std::vector<double>& inv, int p,
                 std::vector<double>& scratch) {
  double* diag = scratch.data();
  double* linv = scratch.data() + p;
  for (int j = 0; j < p; ++j) diag[j] = a[j * p + j];
  for (int j = 0; j < p; ++j) {
    if (!(diag[j] > 0.0) || !std::isfinite(diag[j])) return false;
    double d = a[j * p + j];
    for (int k = 0; k < j; ++k) d -= a[j * p + k] * a[j * p + k];
    if (!(d > kPivotTolerance * diag[j])) return false;
    const double l = std::sqrt(d);
    a[j * p + j] = l;
    for (int i = j + 1; i < p; ++i) {
      double s = a[i * p + j];
      for (int k = 0; k < j; ++k) s -= a[i * p + k] * a[j * p + k];
      a[i * p + j] = s / l;
    }
  }
  // L^{-1} (lower) into `linv`, then inv = L^{-T} L^{-1}.
  std::fill(linv, linv + static_cast<std::size_t>(p) * p, 0.0);
  for (int j = 0; j < p; ++j) {
    linv[j * p + j] = 1.0 / a[j * p + j];
    for (int i = j + 1; i < p; ++i) {
      double s = 0.0;
      for (int k = j; k < i; ++k) s -= a[i * p + k] * linv[k * p + j];
      linv[i * p + j] = s / a[i * p + i];
    }
  }
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j <= i; ++j) {
      double s = 0.0;
      for (int k = i; k < p; ++k) s += linv[k * p + i] * linv[k * p + j];
      inv[i * p + j] = s;
      inv[j * p + i] = s;
    }
  }
  return true;
}

} // namespace

// row_ptr: E + 1 offsets into `rows`/`design`; rows: zero-based panel rows;
// design: N x p (stacked designs); outcome_beta/outcome_se: U x K.
// se_model: 0 = multiplicative random effects (se = se_fixed * sigma, as
// TwoSampleMR mv_ivw / mv_multiple), 1 = multiplicative floored at one
// (se_fixed * max(1, sigma), as TwoSampleMR mr_ivw), 2 = fixed effects.
// design_se (N x p) and design_cor (p x p x E) enable Q_A.  shared_row_se (U),
// shared_outcome_scale (K) and shared_outcome (K, logical) enable the
// shared-weight path.  return_vcov adds the p x p coefficient covariance.
// [[Rcpp::export]]
Rcpp::List fastmr_mvmr_batch_native(
    Rcpp::IntegerVector row_ptr, Rcpp::IntegerVector rows,
    Rcpp::NumericMatrix design, Rcpp::NumericMatrix outcome_beta,
    Rcpp::NumericMatrix outcome_se, int se_model, int threads,
    Rcpp::Nullable<Rcpp::NumericMatrix> design_se = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> design_cor = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> shared_row_se = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> shared_outcome_scale = R_NilValue,
    Rcpp::Nullable<Rcpp::LogicalVector> shared_outcome = R_NilValue,
    bool return_vcov = false) {
  const int design_count = row_ptr.size() - 1;
  const int p = design.ncol();
  const R_xlen_t total_rows = design.nrow();
  const int panel_rows = outcome_beta.nrow();
  const int outcome_count = outcome_beta.ncol();
  if (design_count < 1 || p < 1 || outcome_count < 1) {
    Rcpp::stop("MVMR inputs must have at least one design, exposure and outcome");
  }
  if (outcome_se.nrow() != panel_rows || outcome_se.ncol() != outcome_count) {
    Rcpp::stop("outcome_se must have the dimensions of outcome_beta");
  }
  if (rows.size() != total_rows || row_ptr[0] != 0 || row_ptr[design_count] != total_rows) {
    Rcpp::stop("invalid MVMR design offsets");
  }
  if (se_model < 0 || se_model > 2) Rcpp::stop("invalid se_model");
  for (int e = 0; e < design_count; ++e) {
    if (row_ptr[e + 1] < row_ptr[e]) Rcpp::stop("design offsets must be non-decreasing");
  }
  for (R_xlen_t i = 0; i < total_rows; ++i) {
    if (rows[i] == NA_INTEGER || rows[i] < 0 || rows[i] >= panel_rows) {
      Rcpp::stop("design rows must address the outcome panel");
    }
  }
  // Row-major copy of the design so each SNP's p effects are contiguous.
  std::vector<double> x(static_cast<std::size_t>(total_rows) * p);
  for (R_xlen_t i = 0; i < total_rows; ++i) {
    for (int j = 0; j < p; ++j) {
      const double value = design(i, j);
      if (!std::isfinite(value)) Rcpp::stop("design effects must be finite");
      x[static_cast<std::size_t>(i) * p + j] = value;
    }
  }
  const bool has_qa = design_se.isNotNull();
  std::vector<double> sx;
  std::vector<double> cor;
  if (has_qa) {
    Rcpp::NumericMatrix se_matrix = Rcpp::as<Rcpp::NumericMatrix>(design_se);
    if (se_matrix.nrow() != total_rows || se_matrix.ncol() != p) {
      Rcpp::stop("design_se must have the dimensions of design");
    }
    sx.resize(static_cast<std::size_t>(total_rows) * p);
    for (R_xlen_t i = 0; i < total_rows; ++i) {
      for (int j = 0; j < p; ++j) {
        const double value = se_matrix(i, j);
        if (!std::isfinite(value) || value <= 0.0) {
          Rcpp::stop("design standard errors must be finite and positive");
        }
        sx[static_cast<std::size_t>(i) * p + j] = value;
      }
    }
    if (design_cor.isNull()) Rcpp::stop("design_cor is required with design_se");
    Rcpp::NumericVector c = Rcpp::as<Rcpp::NumericVector>(design_cor);
    if (c.size() != static_cast<R_xlen_t>(p) * p * design_count) {
      Rcpp::stop("design_cor must hold one p x p matrix per design");
    }
    cor.assign(c.begin(), c.end());
  }
  const bool has_shared = shared_row_se.isNotNull();
  std::vector<double> row_se;
  std::vector<double> scale;
  std::vector<int> shared_ok;
  std::vector<double> shared_weight;
  if (has_shared) {
    Rcpp::NumericVector r = Rcpp::as<Rcpp::NumericVector>(shared_row_se);
    Rcpp::NumericVector s = Rcpp::as<Rcpp::NumericVector>(shared_outcome_scale);
    Rcpp::LogicalVector o = Rcpp::as<Rcpp::LogicalVector>(shared_outcome);
    if (r.size() != panel_rows || s.size() != outcome_count || o.size() != outcome_count) {
      Rcpp::stop("shared-weight inputs must have one value per panel row / outcome");
    }
    row_se.assign(r.begin(), r.end());
    scale.assign(s.begin(), s.end());
    shared_ok.assign(o.begin(), o.end());
    shared_weight.resize(row_se.size());
    for (std::size_t i = 0; i < row_se.size(); ++i) {
      shared_weight[i] = row_se[i] > 0.0 ? 1.0 / (row_se[i] * row_se[i]) : MV_NA;
    }
    for (int k = 0; k < outcome_count; ++k) {
      if (shared_ok[k] == NA_LOGICAL) shared_ok[k] = 0;
      if (shared_ok[k] && !(std::isfinite(scale[k]) && scale[k] > 0.0)) shared_ok[k] = 0;
    }
  }

  const std::size_t cells = static_cast<std::size_t>(design_count) * outcome_count;
  Rcpp::NumericVector out_beta(cells * p, MV_NA);
  Rcpp::NumericVector out_se(cells * p, MV_NA);
  Rcpp::NumericVector out_nsnp(cells, 0.0);
  Rcpp::NumericVector out_q(cells, MV_NA);
  Rcpp::NumericVector out_sigma(cells, MV_NA);
  Rcpp::NumericVector out_qa(cells, MV_NA);
  Rcpp::IntegerVector out_path(cells, 0);
  Rcpp::NumericVector out_vcov(return_vcov ? cells * p * p : 0, MV_NA);
  double* beta_ptr = out_beta.begin();
  double* se_ptr = out_se.begin();
  double* nsnp_ptr = out_nsnp.begin();
  double* q_ptr = out_q.begin();
  double* sigma_ptr = out_sigma.begin();
  double* qa_ptr = out_qa.begin();
  int* path_ptr = out_path.begin();
  double* vcov_ptr = return_vcov ? out_vcov.begin() : nullptr;
  const double* yb = outcome_beta.begin();
  const double* ys = outcome_se.begin();
  const int* row_data = rows.begin();
  const int* ptr_data = row_ptr.begin();

  // Shared-path inverse per design (weights 1 / row_se^2); empty when the
  // design is singular or has too few rows.
  std::vector<std::vector<double>> shared_inverse(has_shared ? design_count : 0);
  if (has_shared) {
    mv_parallel(static_cast<std::size_t>(design_count), threads, [&](std::size_t e) {
      const int first = ptr_data[e];
      const int last = ptr_data[e + 1];
      if (last - first < p) return;
      std::vector<double> a(static_cast<std::size_t>(p) * p, 0.0);
      for (int idx = first; idx < last; ++idx) {
        const double r = row_se[row_data[idx]];
        if (!(std::isfinite(r) && r > 0.0)) return;
        const double w = 1.0 / (r * r);
        const double* xi = &x[static_cast<std::size_t>(idx) * p];
        for (int u = 0; u < p; ++u) {
          const double wx = w * xi[u];
          for (int v = 0; v <= u; ++v) a[u * p + v] += wx * xi[v];
        }
      }
      std::vector<double> inv(static_cast<std::size_t>(p) * p);
      std::vector<double> scratch(static_cast<std::size_t>(p) * p + p);
      if (spd_inverse(a, inv, p, scratch)) shared_inverse[e].swap(inv);
    });
  }

  // Jobs: (block of outcomes, design), designs the fast index, so concurrent
  // workers read the same outcome columns (rows shared between designs, such
  // as covariate instruments, stay in cache).
  const int block = 32;
  const int blocks = (outcome_count + block - 1) / block;
  const std::size_t jobs = static_cast<std::size_t>(design_count) * blocks;
  mv_parallel(jobs, threads, [&](std::size_t job) {
    const int e = static_cast<int>(job % design_count);
    const int k_first = static_cast<int>(job / design_count) * block;
    const int k_last = std::min(outcome_count, k_first + block);
    const int first = ptr_data[e];
    const int last = ptr_data[e + 1];
    const int n_rows = last - first;
    std::vector<double> a(static_cast<std::size_t>(p) * p);
    std::vector<double> inv(static_cast<std::size_t>(p) * p);
    std::vector<double> c(static_cast<std::size_t>(p));
    std::vector<double> b(static_cast<std::size_t>(p));
    std::vector<double> scratch(static_cast<std::size_t>(p) * p + p);
    const double* cor_e = has_qa ? &cor[static_cast<std::size_t>(e) * p * p] : nullptr;
    for (int k = k_first; k < k_last; ++k) {
      const std::size_t cell = static_cast<std::size_t>(e) +
                               static_cast<std::size_t>(design_count) * k;
      const double* yk = yb + static_cast<std::size_t>(panel_rows) * k;
      const double* sk = ys + static_cast<std::size_t>(panel_rows) * k;
      bool use_shared = has_shared && shared_ok[k] && !shared_inverse[e].empty();
      double n = 0.0;
      double scale2 = 1.0;  // outcome_scale^2 on the shared path
      std::fill(c.begin(), c.end(), 0.0);
      bool solved = false;
      if (use_shared) {
        // One pass; any missing outcome value sends the pair to the exact path.
        for (int idx = first; idx < last; ++idx) {
          const int i = row_data[idx];
          const double y = yk[i];
          if (!std::isfinite(y)) { use_shared = false; break; }
          const double wy = shared_weight[i] * y;
          const double* xi = &x[static_cast<std::size_t>(idx) * p];
          for (int u = 0; u < p; ++u) c[u] += xi[u] * wy;
        }
        if (!use_shared) std::fill(c.begin(), c.end(), 0.0);
      }
      if (use_shared) {
        n = static_cast<double>(n_rows);
        scale2 = scale[k] * scale[k];
        const std::vector<double>& si = shared_inverse[e];
        for (int u = 0; u < p; ++u) {
          double s = 0.0;
          for (int v = 0; v < p; ++v) s += si[u * p + v] * c[v];
          b[u] = s;
        }
        // Covariance of the fixed-effect fit with weights 1 / (row_se * scale)^2.
        for (int u = 0; u < p * p; ++u) inv[u] = si[u] * scale2;
        solved = true;
        path_ptr[cell] = 1;
      } else {
        std::fill(a.begin(), a.end(), 0.0);
        for (int idx = first; idx < last; ++idx) {
          const int i = row_data[idx];
          const double y = yk[i];
          const double s = sk[i];
          if (!std::isfinite(y) || !std::isfinite(s) || !(s > 0.0)) continue;
          const double w = 1.0 / (s * s);
          const double* xi = &x[static_cast<std::size_t>(idx) * p];
          for (int u = 0; u < p; ++u) {
            const double wx = w * xi[u];
            c[u] += wx * y;
            for (int v = 0; v <= u; ++v) a[u * p + v] += wx * xi[v];
          }
          n += 1.0;
        }
        if (n >= p && spd_inverse(a, inv, p, scratch)) {
          for (int u = 0; u < p; ++u) {
            double s = 0.0;
            for (int v = 0; v < p; ++v) s += inv[u * p + v] * c[v];
            b[u] = s;
          }
          solved = true;
        }
      }
      nsnp_ptr[cell] = n;
      if (!solved) continue;
      bool finite = true;
      for (int u = 0; u < p; ++u) finite = finite && std::isfinite(b[u]);
      if (!finite) continue;
      // Second pass: weighted residual sum of squares (Cochran's Q) and Q_A.
      double q = 0.0;
      double qa = 0.0;
      for (int idx = first; idx < last; ++idx) {
        const int i = row_data[idx];
        const double y = yk[i];
        double var_y;
        if (use_shared) {
          const double r = row_se[i];
          var_y = r * r * scale2;
        } else {
          const double s = sk[i];
          if (!std::isfinite(y) || !std::isfinite(s) || !(s > 0.0)) continue;
          var_y = s * s;
        }
        const double* xi = &x[static_cast<std::size_t>(idx) * p];
        double fit = 0.0;
        for (int u = 0; u < p; ++u) fit += xi[u] * b[u];
        const double residual = y - fit;
        q += residual * residual / var_y;
        if (has_qa) {
          const double* si = &sx[static_cast<std::size_t>(idx) * p];
          double extra = 0.0;
          for (int u = 0; u < p; ++u) {
            for (int v = 0; v < p; ++v) {
              extra += b[u] * b[v] * cor_e[u + p * v] * si[u] * si[v];
            }
          }
          qa += residual * residual / (var_y + extra);
        }
      }
      const double df = n - p;
      double multiplier = 1.0;
      bool se_defined = true;
      if (df > 0.0) {
        const double sigma = std::sqrt(q / df);
        sigma_ptr[cell] = sigma;
        q_ptr[cell] = q;
        if (has_qa) qa_ptr[cell] = qa;
        if (se_model == 0) multiplier = sigma;
        else if (se_model == 1) multiplier = std::max(1.0, sigma);
      } else if (se_model == 0) {
        // n == p: exactly determined; the residual variance is undefined.
        se_defined = false;
      }
      for (int u = 0; u < p; ++u) {
        const std::size_t out = cell + cells * static_cast<std::size_t>(u);
        beta_ptr[out] = b[u];
        if (se_defined) se_ptr[out] = std::sqrt(std::max(inv[u * p + u], 0.0)) * multiplier;
      }
      if (return_vcov && se_defined) {
        const double m2 = multiplier * multiplier;
        for (int u = 0; u < p; ++u) {
          for (int v = 0; v < p; ++v) {
            vcov_ptr[cell + cells * (static_cast<std::size_t>(u) + static_cast<std::size_t>(p) * v)] =
              inv[u * p + v] * m2;
          }
        }
      }
    }
  });

  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["beta"] = out_beta, Rcpp::_["se"] = out_se, Rcpp::_["nsnp"] = out_nsnp,
    Rcpp::_["Q"] = out_q, Rcpp::_["sigma"] = out_sigma, Rcpp::_["Q_A"] = out_qa,
    Rcpp::_["shared"] = out_path);
  if (return_vcov) out["vcov"] = out_vcov;
  return out;
}
