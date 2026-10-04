#include <Rcpp.h>
#include <R_ext/BLAS.h>

#include <algorithm>
#include <atomic>
#include <cfloat>
#include <cstddef>
#include <cmath>
#include <complex>
#include <exception>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <thread>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// The implementation below is a direct R/C++ port of the validated exact
// shared-grid kernel in twosamplemr-fast/native_backend.cpp. R matrices are
// column-major at the API boundary; they are copied once into row-major
// vectors so each pair can stream contiguous SNP columns.

namespace {

constexpr int MODE_GRID_SIZE = 512;
// R's NA_real_ (not a bare NaN), so missing results reach R as NA as in
// TwoSampleMR.  Only ever tested with std::isfinite()/std::isnan() here.
const double NA_VALUE = NA_REAL;
std::atomic<bool> defer_r_math(false);
// Per-thread deferral, for workers that run while the main thread keeps
// computing R-backed p-values for an earlier batch (see
// fastmr_run_groups_boot_native()).
thread_local bool defer_r_math_here = false;

bool r_math_deferred() {
  return defer_r_math_here || defer_r_math.load(std::memory_order_relaxed);
}

double finite_or_na(double x) {
  return std::isfinite(x) ? x : NA_VALUE;
}

double z_pvalue(double statistic) {
  if (r_math_deferred()) return NA_VALUE;
  if (!std::isfinite(statistic)) return NA_VALUE;
  return 2.0 * R::pnorm5(std::abs(statistic), 0.0, 1.0, false, false);
}

double t_pvalue(double statistic, int df) {
  if (r_math_deferred()) return NA_VALUE;
  if (!std::isfinite(statistic) || df <= 0) return NA_VALUE;
  return 2.0 * R::pt(std::abs(statistic), static_cast<double>(df), false, false);
}

double chi_square_pvalue(double q, int df) {
  if (r_math_deferred()) return NA_VALUE;
  if (!std::isfinite(q) || df <= 0) return NA_VALUE;
  return R::pchisq(q, static_cast<double>(df), false, false);
}

double chi_square1_survival(double q) {
  if (!std::isfinite(q) || q < 0.0) return NA_VALUE;
  return std::erfc(std::sqrt(0.5 * q));
}

double safe_statistic(double numerator, double denominator) {
  if (!std::isfinite(numerator) || !std::isfinite(denominator) || denominator == 0.0) {
    return NA_VALUE;
  }
  return numerator / denominator;
}

double sample_std(const std::vector<double>& values) {
  if (values.size() < 2) return NA_VALUE;
  double mean = std::accumulate(values.begin(), values.end(), 0.0) /
                static_cast<double>(values.size());
  double ss = 0.0;
  for (double value : values) {
    const double d = value - mean;
    ss += d * d;
  }
  return std::sqrt(ss / static_cast<double>(values.size() - 1));
}

double median_inplace(std::vector<double>& values) {
  if (values.empty()) return NA_VALUE;
  const std::size_t middle = values.size() / 2;
  if (values.size() % 2 != 0) {
    std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(middle), values.end());
    return values[middle];
  }
  std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(middle - 1), values.end());
  const double lower = values[middle - 1];
  std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(middle), values.end());
  return 0.5 * (lower + values[middle]);
}

double mad(const std::vector<double>& values) {
  if (values.empty()) return NA_VALUE;
  std::vector<double> sorted(values);
  const double center = median_inplace(sorted);
  std::vector<double> deviations;
  deviations.reserve(values.size());
  for (double value : values) deviations.push_back(std::abs(value - center));
  const double result = median_inplace(deviations);
  return 1.4826 * result;
}

double weighted_median_ordered(const double* values, const double* weights,
                               std::size_t count,
                               const std::vector<std::size_t>& order) {
  if (count == 0 || values == nullptr || weights == nullptr || order.size() != count) {
    return NA_VALUE;
  }
  double total = 0.0;
  for (std::size_t index : order) total += weights[index];
  if (!std::isfinite(total) || total <= 0.0) return NA_VALUE;

  double cumulative = 0.0;
  std::size_t last_below = order.size();
  for (std::size_t i = 0; i < order.size(); ++i) {
    const double midpoint = (cumulative + 0.5 * weights[order[i]]) / total;
    cumulative += weights[order[i]];
    if (midpoint < 0.5) last_below = i;
  }
  if (last_below == order.size()) return values[order.front()];
  if (last_below + 1 >= order.size()) return values[order[last_below]];

  cumulative = 0.0;
  for (std::size_t i = 0; i <= last_below; ++i) cumulative += weights[order[i]];
  const double left = (cumulative - 0.5 * weights[order[last_below]]) / total;
  const double right = (cumulative + 0.5 * weights[order[last_below + 1]]) / total;
  const double gap = right - left;
  if (gap <= 0.0) return values[order[last_below]];
  return values[order[last_below]] +
         (values[order[last_below + 1]] - values[order[last_below]]) *
         (0.5 - left) / gap;
}

double weighted_median_point_ptr(const double* values, const double* weights,
                                 std::size_t count, std::vector<std::size_t>& order) {
  if (count == 0 || values == nullptr || weights == nullptr) return NA_VALUE;
  order.resize(count);
  std::iota(order.begin(), order.end(), static_cast<std::size_t>(0));
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    return values[a] < values[b];
  });
  return weighted_median_ordered(values, weights, count, order);
}

double sample_std_ptr(const double* values, std::size_t count) {
  if (count < 2) return NA_VALUE;
  double mean = 0.0;
  for (std::size_t i = 0; i < count; ++i) mean += values[i];
  mean /= static_cast<double>(count);
  double ss = 0.0;
  for (std::size_t i = 0; i < count; ++i) {
    const double d = values[i] - mean;
    ss += d * d;
  }
  return std::sqrt(ss / static_cast<double>(count - 1));
}

double mad_ptr(const double* values, std::size_t count, std::vector<double>& scratch) {
  if (count == 0) return NA_VALUE;
  scratch.assign(values, values + count);
  const double center = median_inplace(scratch);
  for (std::size_t i = 0; i < count; ++i) scratch[i] = std::abs(values[i] - center);
  return 1.4826 * median_inplace(scratch);
}

struct FFTPlan {
  std::size_t n;
  std::vector<std::size_t> bit_reverse;
  std::vector<std::vector<std::complex<double>>> forward_factors;
  std::vector<std::vector<std::complex<double>>> inverse_factors;
  // Largest distance of a stored twiddle factor from cos/sin evaluated
  // directly (the recurrence below accumulates rounding). Only used to bound
  // the transform's rounding error; it does not change the factors.
  double twiddle_error = 0.0;

  explicit FFTPlan(std::size_t size) : n(size), bit_reverse(size) {
    for (std::size_t i = 1, j = 0; i < n; ++i) {
      std::size_t bit = n >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
      j ^= bit;
      bit_reverse[i] = j;
    }
    for (std::size_t length = 2; length <= n; length <<= 1) {
      const double angle = 2.0 * 3.14159265358979323846 / static_cast<double>(length);
      const std::size_t half = length >> 1;
      std::vector<std::complex<double>> forward(half);
      std::vector<std::complex<double>> inverse(half);
      const std::complex<double> forward_step(std::cos(-angle), std::sin(-angle));
      const std::complex<double> inverse_step(std::cos(angle), std::sin(angle));
      std::complex<double> forward_factor(1.0, 0.0);
      std::complex<double> inverse_factor(1.0, 0.0);
      for (std::size_t i = 0; i < half; ++i) {
        forward[i] = forward_factor;
        inverse[i] = inverse_factor;
        const double theta = angle * static_cast<double>(i);
        const std::complex<double> direct(std::cos(theta), std::sin(theta));
        twiddle_error = std::max(twiddle_error, std::abs(inverse_factor - direct));
        twiddle_error = std::max(twiddle_error, std::abs(forward_factor - std::conj(direct)));
        forward_factor *= forward_step;
        inverse_factor *= inverse_step;
      }
      forward_factors.emplace_back(std::move(forward));
      inverse_factors.emplace_back(std::move(inverse));
    }
  }
};

const FFTPlan& fft_plan(std::size_t n) {
  static const FFTPlan plan_1024(1024);
  if (n == plan_1024.n) return plan_1024;
  static thread_local std::unique_ptr<FFTPlan> fallback;
  if (!fallback || fallback->n != n) fallback = std::make_unique<FFTPlan>(n);
  return *fallback;
}

// Bound on ||fft(x) - FFT(x)||_2 / ||FFT(x)||_2 for fft_inplace() with this
// plan (Higham 2002, Accuracy and Stability of Numerical Algorithms, Thm 24.2:
// log2(n) eta / (1 - log2(n) eta), eta = mu + gamma_4 (sqrt(2) + mu), with mu
// the twiddle-factor error). mu is the measured error plus 16 unit roundoffs
// for the error of the cos/sin reference itself.
double fft_relative_error_bound(const FFTPlan& plan) {
  const double u = DBL_EPSILON / 2.0;
  const double mu = plan.twiddle_error + 16.0 * u;
  const double gamma4 = 4.0 * u / (1.0 - 4.0 * u);
  const double eta = mu + gamma4 * (std::sqrt(2.0) + mu);
  const double stages = std::log2(static_cast<double>(plan.n));
  return stages * eta / (1.0 - stages * eta);
}

// Complex product in the same operation order as the compiler's complex
// multiply (libgcc / compiler-rt __muldc3) for finite operands, without its
// NaN-recovery branch, which blocks vectorisation. Mode inputs are finite, so
// the result is bit-identical. Separate statements keep compilers that only
// contract within one expression (clang's default -ffp-contract=on) from
// fusing a multiply-add here.
inline std::complex<double> cmul(const std::complex<double>& a, const std::complex<double>& b) {
  const double ac = a.real() * b.real();
  const double bd = a.imag() * b.imag();
  const double ad = a.real() * b.imag();
  const double bc = a.imag() * b.real();
  return std::complex<double>(ac - bd, ad + bc);
}

void fft_inplace(std::vector<std::complex<double>>& values, bool inverse) {
  const std::size_t n = values.size();
  const FFTPlan& plan = fft_plan(n);
  for (std::size_t i = 1; i < n; ++i) {
    const std::size_t j = plan.bit_reverse[i];
    if (i < j) std::swap(values[i], values[j]);
  }
  for (std::size_t stage = 0, length = 2; length <= n; ++stage, length <<= 1) {
    const std::vector<std::complex<double>>& factors = inverse
      ? plan.inverse_factors[stage] : plan.forward_factors[stage];
    for (std::size_t start = 0; start < n; start += length) {
      const std::size_t half = length >> 1;
      for (std::size_t i = 0; i < half; ++i) {
        const std::complex<double> even = values[start + i];
        const std::complex<double> odd = cmul(factors[i], values[start + i + half]);
        values[start + i] = even + odd;
        values[start + i + half] = even - odd;
      }
    }
  }
  if (inverse) {
    const double scale = 1.0 / static_cast<double>(n);
    for (std::complex<double>& value : values) value *= scale;
  }
}

struct ModeDensityWorkspace {
  std::vector<double> scratch;
  std::vector<std::complex<double>> binned;
  std::vector<std::complex<double>> simple;
  std::vector<std::complex<double>> weighted;
  std::vector<std::complex<double>> kernel;
  // Direct path: binned weights and densities per weight vector, and the
  // truncated symmetric kernel.
  std::vector<double> direct_bins[2];
  std::vector<double> direct_density[2];
  std::vector<double> direct_kernel;
  // Hull path: bins are kept all-zero between calls (only the occupied cells
  // are cleared after use); densities are only written and read on the hull.
  std::vector<double> hull_bins[2];
  std::vector<double> hull_density[2];
  std::vector<double> hull_kernel;
};

ModeDensityWorkspace& mode_workspace() {
  static thread_local ModeDensityWorkspace workspace;
  return workspace;
}

// Bandwidth and grids shared by every mode-density path: density() on n = 512
// points with a bw.nrd0-style bandwidth (times phi) and cut = 3, linearly
// binned on [lo, up], then approx() back onto the output grid [from, to].
struct ModeGrid {
  double bandwidth;
  double from;
  double output_step;
  double lo;
  double delta;
  double position_start;
  double position_step;
};

ModeGrid mode_grid(const double* values, std::size_t count, double phi,
                   std::vector<double>& scratch) {
  scratch.clear();
  scratch.reserve(count);
  const double raw_bandwidth = 0.9 * std::min(sample_std_ptr(values, count),
                                                mad_ptr(values, count, scratch)) /
                               std::pow(static_cast<double>(count), 0.2);
  double bandwidth = std::isfinite(raw_bandwidth) ? std::max(1e-8, raw_bandwidth) : 1e-8;
  bandwidth *= phi;
  double minimum = values[0], maximum = values[0];
  for (std::size_t i = 1; i < count; ++i) {
    minimum = std::min(minimum, values[i]);
    maximum = std::max(maximum, values[i]);
  }
  const double from = minimum - 3.0 * bandwidth;
  const double to = maximum + 3.0 * bandwidth;
  const int n = MODE_GRID_SIZE;
  const double lo = from - 4.0 * bandwidth;
  const double up = to + 4.0 * bandwidth;
  ModeGrid grid;
  grid.bandwidth = bandwidth;
  grid.from = from;
  grid.lo = lo;
  grid.delta = (up - lo) / static_cast<double>(n - 1);
  grid.output_step = (to - from) / static_cast<double>(n - 1);
  grid.position_start = (from - lo) / grid.delta;
  grid.position_step = grid.output_step / grid.delta;
  return grid;
}

// Linear binning of `weights` onto the n-point grid (density()'s BinDist).
// Calls add(cell, amount) in the same order for every path, so all paths see
// bit-identical bins.
template <typename Add>
void mode_bin(const double* values, const double* weights, std::size_t count,
              const ModeGrid& grid, Add add) {
  const int n = MODE_GRID_SIZE;
  for (std::size_t i = 0; i < count; ++i) {
    const double xpos = (values[i] - grid.lo) / grid.delta;
    if (!std::isfinite(xpos) || xpos > static_cast<double>(std::numeric_limits<int>::max()) ||
        xpos < static_cast<double>(std::numeric_limits<int>::min())) continue;
    const int index = static_cast<int>(std::floor(xpos));
    const double fraction = xpos - static_cast<double>(index);
    if (0 <= index && index <= n - 2) {
      add(index, (1.0 - fraction) * weights[i]);
      add(index + 1, fraction * weights[i]);
    } else if (index == -1) {
      add(0, fraction * weights[i]);
    } else if (index == n - 1) {
      add(index, (1.0 - fraction) * weights[i]);
    }
  }
}

// Gaussian kernel at distance d * delta, as density() evaluates it.
inline double mode_kernel(int d, double delta, double bandwidth) {
  const double distance = static_cast<double>(d) * delta;
  const double z = distance / bandwidth;
  return std::exp(-0.5 * z * z) / (bandwidth * std::sqrt(2.0 * 3.14159265358979323846));
}

// Argmax of the interpolated density on the output grid (first maximum wins).
// `density(cell)` returns the convolved density at a binning cell. Also
// reports the largest density at any other output point (`second`).
template <typename Density>
int mode_argmax(const ModeGrid& grid, Density density_at, double& best, double& second) {
  const int n = MODE_GRID_SIZE;
  int best_index = 0;
  best = -std::numeric_limits<double>::infinity();
  second = best;
  for (int i = 0; i < n; ++i) {
    const double position = grid.position_start + grid.position_step * static_cast<double>(i);
    const int left = static_cast<int>(std::floor(position));
    double density = 0.0;
    if (left < 0) density = std::max(0.0, density_at(0));
    else if (left >= n - 1) density = std::max(0.0, density_at(n - 1));
    else {
      const double fraction = position - static_cast<double>(left);
      density = (1.0 - fraction) * density_at(left) + fraction * density_at(left + 1);
      density = std::max(0.0, density);
    }
    if (density > best) {
      second = best;
      best = density;
      best_index = i;
    } else if (density > second) {
      second = density;
    }
  }
  return best_index;
}

// FFT path: zero-padded circular convolution (length 2n) of the bins with the
// kernel, as density() computes it. Returns the argmax output-grid index.
int mode_index_fft(const double* values, const double* weights, std::size_t count,
                   const ModeGrid& grid) {
  const int n = MODE_GRID_SIZE;
  const int length = 2 * n;
  ModeDensityWorkspace& workspace = mode_workspace();
  std::vector<std::complex<double>>& binned = workspace.binned;
  binned.assign(length, std::complex<double>(0.0, 0.0));
  mode_bin(values, weights, count, grid, [&](int cell, double amount) { binned[cell] += amount; });
  std::vector<std::complex<double>>& kernel = workspace.kernel;
  kernel.resize(length);
  // The kernel is even: kernel[length - i] == kernel[i] bit for bit (the
  // distance only changes sign), so evaluate exp() once per distance.
  for (int i = 0; i <= n; ++i) kernel[i] = mode_kernel(i, grid.delta, grid.bandwidth);
  for (int i = n + 1; i < length; ++i) kernel[i] = kernel[length - i];
  fft_inplace(binned, false);
  fft_inplace(kernel, false);
  for (int i = 0; i < length; ++i) binned[i] = cmul(binned[i], std::conj(kernel[i]));
  fft_inplace(binned, true);
  double best = 0.0, second = 0.0;
  return mode_argmax(grid, [&](int cell) { return binned[cell].real(); }, best, second);
}

// FFT path for two weight vectors on the same ratios (one kernel transform).
std::pair<int, int> mode_index_fft_pair(const double* values, const double* simple_weights,
                                        const double* weighted_weights, std::size_t count,
                                        const ModeGrid& grid) {
  const int n = MODE_GRID_SIZE;
  const int length = 2 * n;
  ModeDensityWorkspace& workspace = mode_workspace();
  std::vector<std::complex<double>>& simple = workspace.simple;
  std::vector<std::complex<double>>& weighted = workspace.weighted;
  simple.assign(length, std::complex<double>(0.0, 0.0));
  weighted.assign(length, std::complex<double>(0.0, 0.0));
  mode_bin(values, simple_weights, count, grid, [&](int cell, double amount) { simple[cell] += amount; });
  mode_bin(values, weighted_weights, count, grid, [&](int cell, double amount) { weighted[cell] += amount; });
  std::vector<std::complex<double>>& kernel = workspace.kernel;
  kernel.resize(length);
  for (int i = 0; i <= n; ++i) kernel[i] = mode_kernel(i, grid.delta, grid.bandwidth);
  for (int i = n + 1; i < length; ++i) kernel[i] = kernel[length - i];
  fft_inplace(simple, false);
  fft_inplace(weighted, false);
  fft_inplace(kernel, false);
  for (int i = 0; i < length; ++i) {
    simple[i] = cmul(simple[i], std::conj(kernel[i]));
    weighted[i] = cmul(weighted[i], std::conj(kernel[i]));
  }
  fft_inplace(simple, true);
  fft_inplace(weighted, true);
  double best = 0.0, second = 0.0;
  const int simple_index = mode_argmax(grid, [&](int cell) { return simple[cell].real(); }, best, second);
  const int weighted_index = mode_argmax(grid, [&](int cell) { return weighted[cell].real(); }, best, second);
  return std::make_pair(simple_index, weighted_index);
}

// ---- Direct path ---------------------------------------------------------
// The bins have at most 2 * count non-zero cells, so for small count the
// linear convolution y[j] = sum_m b[m] k(|j - m| delta) is far cheaper than
// the FFT path's five 1024-point transforms (there is no wrap-around: the FFT
// is zero-padded to 2n). Same bandwidth, grids, bins, kernel values and
// argmax rule; the densities differ from the FFT's only by rounding, so the
// chosen grid point can only differ on a near-tie. A guard detects every such
// draw and recomputes it on the FFT path, so the selected mode is always the
// FFT path's (and main's) bit for bit.

// Ratios up to which the direct path is tried. Its convolution work is bounded
// by the grid (at most n non-zero bins, kernel radius at most ~366 cells), so
// it was faster than the FFT at every size measured: 11-13x at 3 ratios, 2-3x
// at 1000, 1.1x at 20000 (both paths are then dominated by the shared
// bandwidth/binning work). Above the largest size measured, keep the FFT.
constexpr double kModeDirectMaxRatios = 20000.0;
// The kernel is truncated beyond this many bandwidths: exp(-50) ~ 2e-22 of
// its peak, far below the guard's tolerance (bounded explicitly below).
constexpr double kModeKernelCutoff = 10.0;
// Floor of the guard's relative tolerance (see mode_index_direct()).
constexpr double kModeGuardRelative = 1e-9;

std::atomic<double> mode_direct_max(kModeDirectMaxRatios);
std::atomic<unsigned long long> mode_direct_total(0), mode_guard_total(0);
std::atomic<unsigned long long> mode_hull_total(0), mode_hull_fallback_total(0);
// Per-thread draw counts, flushed into the totals once per mode fit.
thread_local unsigned long long mode_direct_draws = 0, mode_guard_draws = 0;
thread_local unsigned long long mode_hull_draws = 0, mode_hull_fallback_draws = 0;

void flush_mode_counters() {
  if (mode_direct_draws) mode_direct_total.fetch_add(mode_direct_draws, std::memory_order_relaxed);
  if (mode_guard_draws) mode_guard_total.fetch_add(mode_guard_draws, std::memory_order_relaxed);
  if (mode_hull_draws) mode_hull_total.fetch_add(mode_hull_draws, std::memory_order_relaxed);
  if (mode_hull_fallback_draws) {
    mode_hull_fallback_total.fetch_add(mode_hull_fallback_draws, std::memory_order_relaxed);
  }
  mode_direct_draws = 0;
  mode_guard_draws = 0;
  mode_hull_draws = 0;
  mode_hull_fallback_draws = 0;
}

#ifdef _OPENMP
#define FASTMR_SIMD _Pragma("omp simd")
#else
#define FASTMR_SIMD
#endif

// Direct densities for `vectors` (1 or 2) weight vectors; stores each argmax
// output-grid index in index[]. Returns false, leaving the draw to the FFT
// path, unless every argmax is certain to equal the FFT path's.
//
// The guard. For each density, |d_fft(i) - d_direct(i)| <= D at every output
// point i, so if best - second > 2 D the FFT has the same unique argmax. With
// b the bins, W = sum(b) and k the FFT path's length-2n kernel vector, the
// error analysis of an FFT convolution (forward transforms of b and k, one
// complex product, inverse transform; each transform within relative 2-norm
// error tau, Higham Thm 24.2, from the plan's measured twiddle error) gives
//   ||y_fft - y||_inf <= E_fft = tau (2 ||b||_2 ||k||_1 + W ||k||_2)
//                                + 3 u ||b||_2 ||k||_1   (+ O(tau^2)),
// while the direct sums of at most `nonzero` non-negative terms give
//   ||y_direct - y||_inf <= E_direct = (nonzero + 2) u max(y) + truncation,
// truncation <= W k(cutoff) (zero when nothing is truncated). The linear
// interpolation adds at most 8 u max(y). D = 2 E_fft + E_direct + 8 u max(y)
// (E_fft doubled to absorb the second-order terms), and the guard falls back
// when best - second <= max(2 D, 1e-9 best). The 1e-9 floor is a fixed margin
// on top of the bound: measured over 135k random draws (3-1000 ratios), D was
// 1e-11 to 8e-11 of best and the actual discrepancy at most 3.8e-14 of best
// (at most 1.3% of D), so the floor is what binds, ~2.6e4 times the largest
// discrepancy; gaps that small occurred in ~1e-5 of those draws, so the
// fallback costs nothing measurable.
bool mode_index_direct(const double* values, const double* const* weights, int vectors,
                       std::size_t count, const ModeGrid& grid, int* index) {
  const int n = MODE_GRID_SIZE;
  const double delta = grid.delta;
  const double bandwidth = grid.bandwidth;
  const double radius_cells = kModeKernelCutoff * bandwidth / delta;
  if (!std::isfinite(radius_cells) || !(delta > 0.0) || !(bandwidth > 0.0)) return false;
  ModeDensityWorkspace& workspace = mode_workspace();
  int touched_lo = n, touched_hi = -1;
  for (int v = 0; v < vectors; ++v) {
    std::vector<double>& bins = workspace.direct_bins[v];
    bins.assign(n, 0.0);
    mode_bin(values, weights[v], count, grid, [&](int cell, double amount) {
      bins[cell] += amount;
      touched_lo = std::min(touched_lo, cell);
      touched_hi = std::max(touched_hi, cell);
    });
  }
  // Truncated symmetric kernel: kernel[radius + d] = kernel[radius - d] = k(d).
  const int radius = radius_cells < static_cast<double>(n - 2)
    ? static_cast<int>(radius_cells) + 1 : n - 1;
  std::vector<double>& kernel = workspace.direct_kernel;
  kernel.resize(2 * static_cast<std::size_t>(radius) + 1);
  double kernel_l1 = 0.0, kernel_l2 = 0.0;  // norms of the FFT path's kernel vector
  for (int d = 0; d <= radius; ++d) {
    const double value = mode_kernel(d, delta, bandwidth);
    kernel[radius + d] = value;
    kernel[radius - d] = value;
    const double copies = d == 0 ? 1.0 : 2.0;  // entries d and 2n - d
    kernel_l1 += copies * value;
    kernel_l2 += copies * value * value;
  }
  // Entries the direct path does not use: distances radius+1..n (n once).
  // Beyond the cutoff each is at most k(cutoff); with no truncation only
  // distance n remains, which no output cell needs (|j - m| <= n - 1).
  const double tail_value = radius < n - 1
    ? mode_kernel(0, delta, bandwidth) * std::exp(-0.5 * kModeKernelCutoff * kModeKernelCutoff)
    : mode_kernel(n, delta, bandwidth);
  const double tail_count = static_cast<double>(2 * n - (2 * radius + 1));
  kernel_l1 += tail_count * tail_value;
  kernel_l2 += tail_count * tail_value * tail_value;
  // Cells the output grid reads (positions are increasing in i).
  const double first_position = grid.position_start;
  const double last_position = grid.position_start + grid.position_step * static_cast<double>(n - 1);
  const int cell_lo = static_cast<int>(std::max(0.0, std::min(static_cast<double>(n - 1), std::floor(first_position))));
  const int cell_hi = static_cast<int>(std::max(0.0, std::min(static_cast<double>(n - 1), std::floor(last_position) + 1.0)));
  double weight_sum[2] = {0.0, 0.0}, weight_ss[2] = {0.0, 0.0};
  int nonzero = 0;
  for (int v = 0; v < vectors; ++v) workspace.direct_density[v].assign(n, 0.0);
  double* const y0 = workspace.direct_density[0].data();
  double* const y1 = workspace.direct_density[vectors > 1 ? 1 : 0].data();
  const double* const b0 = workspace.direct_bins[0].data();
  const double* const b1 = workspace.direct_bins[vectors > 1 ? 1 : 0].data();
  for (int m = touched_lo; m <= touched_hi; ++m) {
    const double c0 = b0[m];
    const double c1 = vectors > 1 ? b1[m] : 0.0;
    if (c0 == 0.0 && c1 == 0.0) continue;
    ++nonzero;
    weight_sum[0] += c0;
    weight_ss[0] += c0 * c0;
    weight_sum[1] += c1;
    weight_ss[1] += c1 * c1;
    const int first = std::max(cell_lo, m - radius);
    const int last = std::min(cell_hi, m + radius);
    if (first > last) continue;
    const double* const k = kernel.data() + (first - m + radius);
    const int width = last - first + 1;
    double* const t0 = y0 + first;
    if (vectors == 1) {
      FASTMR_SIMD
      for (int t = 0; t < width; ++t) t0[t] += c0 * k[t];
    } else {
      double* const t1 = y1 + first;
      FASTMR_SIMD
      for (int t = 0; t < width; ++t) {
        t0[t] += c0 * k[t];
        t1[t] += c1 * k[t];
      }
    }
  }
  static const double tau = fft_relative_error_bound(fft_plan(2 * static_cast<std::size_t>(n)));
  const double u = DBL_EPSILON / 2.0;
  for (int v = 0; v < vectors; ++v) {
    const double* const y = workspace.direct_density[v].data();
    double best = 0.0, second = 0.0;
    const int best_index = mode_argmax(grid, [&](int cell) { return y[cell]; }, best, second);
    double max_y = 0.0;
    for (int j = cell_lo; j <= cell_hi; ++j) max_y = std::max(max_y, y[j]);
    const double b_l2 = std::sqrt(weight_ss[v]);
    const double s1 = b_l2 * kernel_l1;
    const double s2 = weight_sum[v] * std::sqrt(kernel_l2);
    const double fft_error = tau * (2.0 * s1 + s2) + 3.0 * u * s1;
    const double truncation = radius < n - 1 ? weight_sum[v] * tail_value : 0.0;
    const double direct_error = static_cast<double>(nonzero + 2) * u * max_y + truncation;
    const double discrepancy = 2.0 * fft_error + direct_error + 8.0 * u * max_y;
    const double tolerance = std::max(2.0 * discrepancy, kModeGuardRelative * best);
    if (!(best - second > tolerance)) return false;
    index[v] = best_index;
  }
  return true;
}

// ---- Hull path -------------------------------------------------------------
// For small pairs the direct path's cost is the 512-cell grid, not the data:
// it evaluates ~10 bandwidths of kernel (one exp() per distance, up to ~366),
// convolves every occupied bin over that whole radius and scans all 512 output
// points. The hull path restricts all three to the cells spanning the occupied
// bins, and (optionally) evaluates the kernel by a recurrence with an exact
// exp() every 16th distance. It returns false whenever it cannot certify the
// FFT path's argmax; the caller then runs mode_index_direct() (and, failing
// that, the FFT), so every selected mode is still the FFT path's bit for bit.
//
// Why the hull suffices. Let [A, B] be the occupied bins (every bin outside is
// exactly 0) and y*(j) = sum_m b[m] g(|j - m|) the exact density with the true
// Gaussian g. For j < j' <= A every term has |j' - m| < |j - m|, so y* is
// non-decreasing on cells [0, A]; likewise non-increasing on [B, n - 1]. The
// output point at position p interpolates cells floor(p) and floor(p) + 1, so
// its exact density is non-decreasing in p on [0, A] and non-increasing on
// [B, n - 1]. We scan exactly the output points whose left cell lies in
// [A - 1, B] (positions [A - 1, B + 1)), which read cells [A - 1, B + 1] only.
// The first scanned point (position < A) bounds every earlier point and the
// last (position >= B) every later one. If the best scanned point is strictly
// inside the scanned range, both boundary points are counted in `second`, so
// the guard below also separates best from every unscanned point. The FFT and
// this path use the floating-point kernel K_f, not g; for distances within
// 11 bandwidths |log(K_f / g)| <= 309 u (see the recurrence bound below) and
// beyond 10 bandwidths K_f <= 2 tail with tail = g(0) exp(-50). Hence for an
// unscanned point o and its boundary point q,
//   y_K(o) <= (1 + 310 u) / (1 - 310 u) y_K(q) + 4 W tail,
// so the guard adds hull_term = 1024 u max(y) + 4 W tail (generous; it also
// absorbs the FFT's 8 u interpolation rounding at o) to D before doubling.
//
// The FFT kernel-vector norms the guard needs are replaced by upper bounds
// (g is unimodal, so sum_d g(d delta) <= 1/delta + g(0) and
// sum_d g(d delta)^2 <= 1/(2 h sqrt(pi) delta) + g(0)^2), inflated by 1e-12
// relative for K_f vs g and by 2n entries of 2 tail for the far tail.
//
// The recurrence kernel. With a = delta / h, g(d + 1) = g(d) q(d) where
// q(d) = exp(-a^2 (2d + 1) / 2) and q(d + 1) = q(d) r, r = exp(-a^2). In each
// block of 16 distances d0..d0+15, K_rec(d0) = K_f(d0) exactly (mode_kernel())
// and K_rec(d0 + t) = fl(K_rec(d0 + t - 1) q^(t - 1)), with q^0 = fl(exp(-aa
// (2 d0 + 1) / 2)), q^j = fl(q^(j - 1) r^), aa = fl(a^ a^), a^ = fl(delta / h),
// r^ = fl(exp(-aa)). It is used only when a <= 1 and 1e-250 <= g(0) <= 1e300,
// so every distance used has z = d a <= 10 + a <= 11 and every intermediate
// is a normal double (no underflow, so fl(x op y) = (x op y)(1 + e), |e| <= u,
// u = 2^-53). Assume exp() is within 2 ulps (relative error <= 4 u; glibc,
// Apple libm and the UCRT are all within 1 ulp). Then, to first order (the
// constants below already include a 1% margin for higher-order terms):
//  * K_f(d) = fl(fl(exp(-0.5 fl(z^ z^))) / c), z^ = fl(fl(d delta) / h): the
//    exponent has relative error <= 5 u, i.e. absolute <= 5 u z^2 / 2 <= 303 u
//    at z <= 11, plus 4 u (exp) and u (division): |log(K_f / g)| <= 309 u.
//  * |log(q^0 / q(d0))| <= 4 u * 11.5 + 4 u = 50.3 u (exponent a^2 (2 d0 + 1)/2
//    <= a z + a^2 / 2 <= 11.5, four roundings), |log(r^ / r)| <= 3 u a^2 + 4 u
//    <= 7.1 u, so |log(q^j / q(d0 + j))| <= 50.3 u + j (7.1 u + 1.01 u).
//  * Summing j = 0..t-1 for t <= 15 gives <= 754.5 u + 851.6 u = 1606.1 u;
//    the t products add <= 15.2 u and K_f(d0) itself 309 u, so
//    |log(K_rec / g)| <= 1930.3 u and |log(K_rec / K_f)| <= 2239.3 u.
// So |K_rec(d) - K_f(d)| <= 2240 u K_f(d) for every distance, and the direct
// sums over the hull differ from the same sums with K_f by at most
// 2240 u y_K <= 2240 u max(y) (+ second order). The guard uses
// kernel_rel = 4096 u (about 9.1e-13), added to the direct-path error term.
// Every quantity in D remains a proven upper bound, as before, and the 1e-9
// relative floor stays on top.

// Ratios up to which the hull path is tried first (finding: it wins at small
// k, where the 512-cell grid dominates; above this the direct path's
// per-bin work dominates and the hull saves little). 0 disables it.
constexpr double kModeHullMaxRatios = 64.0;
// Kernel blocks: one exact exp() per this many distances.
constexpr int kModeRecurrenceBlock = 16;
std::atomic<double> mode_hull_max(kModeHullMaxRatios);
std::atomic<bool> mode_hull_recurrence(true);

bool mode_index_hull(const double* values, const double* const* weights, int vectors,
                     std::size_t count, const ModeGrid& grid, int* index) {
  const int n = MODE_GRID_SIZE;
  const double delta = grid.delta;
  const double bandwidth = grid.bandwidth;
  const double radius_cells = kModeKernelCutoff * bandwidth / delta;
  if (!std::isfinite(radius_cells) || !(delta > 0.0) || !(bandwidth > 0.0)) return false;
  ModeDensityWorkspace& workspace = mode_workspace();
  int touched_lo = n, touched_hi = -1;
  for (int v = 0; v < vectors; ++v) {
    std::vector<double>& bins = workspace.hull_bins[v];
    if (bins.size() != static_cast<std::size_t>(n)) bins.assign(n, 0.0);
    double* const b = bins.data();
    mode_bin(values, weights[v], count, grid, [&](int cell, double amount) {
      b[cell] += amount;
      touched_lo = std::min(touched_lo, cell);
      touched_hi = std::max(touched_hi, cell);
    });
  }
  // Restore the all-zero invariant on every exit.
  struct ClearBins {
    ModeDensityWorkspace& ws; int vectors; const int& lo; const int& hi;
    ~ClearBins() {
      for (int v = 0; v < vectors; ++v) {
        double* const b = ws.hull_bins[v].data();
        for (int m = lo; m <= hi; ++m) b[m] = 0.0;
      }
    }
  } clear_bins{workspace, vectors, touched_lo, touched_hi};
  if (touched_hi < 0) return false;
  const int A = touched_lo, B = touched_hi;
  const int lo = A - 1, hi = B + 1;
  // Cells the full output grid reads; the hull must sit strictly inside, so
  // every scanned point interpolates two in-range cells and has unscanned
  // neighbours handled by the monotonicity argument.
  const double first_position = grid.position_start;
  const double last_position = grid.position_start + grid.position_step * static_cast<double>(n - 1);
  const int cell_lo = static_cast<int>(std::max(0.0, std::min(static_cast<double>(n - 1), std::floor(first_position))));
  const int cell_hi = static_cast<int>(std::max(0.0, std::min(static_cast<double>(n - 1), std::floor(last_position) + 1.0)));
  if (lo < cell_lo + 1 || hi > cell_hi - 1) return false;
  const int span = hi - lo;  // >= every |j - m| needed (j in [lo, hi], m in [A, B])
  const int full_radius = radius_cells < static_cast<double>(n - 2)
    ? static_cast<int>(radius_cells) + 1 : n - 1;
  const int radius = std::min(span, full_radius);
  const bool truncated = radius < span;
  const double k0 = mode_kernel(0, delta, bandwidth);
  if (!std::isfinite(k0) || !(k0 > 0.0)) return false;
  const double tail_value = k0 * std::exp(-0.5 * kModeKernelCutoff * kModeKernelCutoff);
  std::vector<double>& kernel = workspace.hull_kernel;
  kernel.resize(2 * static_cast<std::size_t>(radius) + 1);
  double* const kc = kernel.data() + radius;  // kc[d] = kc[-d] = K(d)
  double kernel_rel = 0.0;
  const double a = delta / bandwidth;
  const bool recurrence = mode_hull_recurrence.load(std::memory_order_relaxed) &&
    a <= 1.0 && k0 >= 1e-250 && k0 <= 1e300 && radius >= kModeRecurrenceBlock;
  if (recurrence) {
    const double aa = a * a;
    const double r = std::exp(-aa);
    for (int d0 = 0; d0 <= radius; d0 += kModeRecurrenceBlock) {
      double value = mode_kernel(d0, delta, bandwidth);
      double q = std::exp(-0.5 * aa * (2.0 * static_cast<double>(d0) + 1.0));
      const int end = std::min(radius, d0 + kModeRecurrenceBlock - 1);
      for (int d = d0; d <= end; ++d) {
        kc[d] = value;
        kc[-d] = value;
        value *= q;
        q *= r;
      }
    }
    kernel_rel = 4096.0 * (DBL_EPSILON / 2.0);
  } else {
    for (int d = 0; d <= radius; ++d) {
      const double value = mode_kernel(d, delta, bandwidth);
      kc[d] = value;
      kc[-d] = value;
    }
  }
  // Upper bounds of the FFT path's length-2n kernel-vector norms (see above).
  const double far = 4.0 * static_cast<double>(n) * tail_value;
  const double kernel_l1 = (1.0 / delta + k0) * (1.0 + 1e-12) + far;
  const double kernel_l2 = std::sqrt(
    (1.0 / (2.0 * bandwidth * std::sqrt(3.14159265358979323846) * delta) + k0 * k0) * (1.0 + 1e-12) +
    2.0 * far * tail_value);
  double weight_sum[2] = {0.0, 0.0}, weight_ss[2] = {0.0, 0.0};
  int nonzero = 0;
  for (int v = 0; v < vectors; ++v) {
    std::vector<double>& density = workspace.hull_density[v];
    if (density.size() != static_cast<std::size_t>(n)) density.assign(n, 0.0);
    std::fill(density.begin() + lo, density.begin() + hi + 1, 0.0);
  }
  double* const y0 = workspace.hull_density[0].data();
  double* const y1 = workspace.hull_density[vectors > 1 ? 1 : 0].data();
  const double* const b0 = workspace.hull_bins[0].data();
  const double* const b1 = workspace.hull_bins[vectors > 1 ? 1 : 0].data();
  for (int m = A; m <= B; ++m) {
    const double c0 = b0[m];
    const double c1 = vectors > 1 ? b1[m] : 0.0;
    if (c0 == 0.0 && c1 == 0.0) continue;
    ++nonzero;
    weight_sum[0] += c0;
    weight_ss[0] += c0 * c0;
    weight_sum[1] += c1;
    weight_ss[1] += c1 * c1;
    const int first = std::max(lo, m - radius);
    const int last = std::min(hi, m + radius);
    const double* const k = kc + (first - m);
    const int width = last - first + 1;
    double* const t0 = y0 + first;
    if (vectors == 1) {
      FASTMR_SIMD
      for (int t = 0; t < width; ++t) t0[t] += c0 * k[t];
    } else {
      double* const t1 = y1 + first;
      FASTMR_SIMD
      for (int t = 0; t < width; ++t) {
        t0[t] += c0 * k[t];
        t1[t] += c1 * k[t];
      }
    }
  }
  // Output points with left cell in [lo, B], computed exactly as mode_argmax().
  const double ps = grid.position_start, step = grid.position_step;
  auto left_of = [&](int i) {
    return static_cast<int>(std::floor(ps + step * static_cast<double>(i)));
  };
  int i_first = static_cast<int>(std::max(0.0, std::min(static_cast<double>(n - 1),
    std::ceil((static_cast<double>(lo) - ps) / step) - 2.0)));
  while (i_first > 0 && left_of(i_first) >= lo) --i_first;
  while (i_first < n && left_of(i_first) < lo) ++i_first;
  if (i_first >= n || left_of(i_first) > B) return false;
  int i_last = i_first;
  while (i_last + 1 < n && left_of(i_last + 1) <= B) ++i_last;
  if (i_last - i_first < 2) return false;
  static const double tau = fft_relative_error_bound(fft_plan(2 * static_cast<std::size_t>(n)));
  const double u = DBL_EPSILON / 2.0;
  for (int v = 0; v < vectors; ++v) {
    const double* const y = workspace.hull_density[v].data();
    int best_index = i_first;
    double best = -std::numeric_limits<double>::infinity();
    double second = best;
    for (int i = i_first; i <= i_last; ++i) {
      const double position = ps + step * static_cast<double>(i);
      const int left = static_cast<int>(std::floor(position));
      const double fraction = position - static_cast<double>(left);
      double density = (1.0 - fraction) * y[left] + fraction * y[left + 1];
      density = std::max(0.0, density);
      if (density > best) {
        second = best;
        best = density;
        best_index = i;
      } else if (density > second) {
        second = density;
      }
    }
    if (best_index == i_first || best_index == i_last) return false;
    double max_y = 0.0;
    for (int j = lo; j <= hi; ++j) max_y = std::max(max_y, y[j]);
    const double s1 = std::sqrt(weight_ss[v]) * kernel_l1;
    const double s2 = weight_sum[v] * kernel_l2;
    const double fft_error = tau * (2.0 * s1 + s2) + 3.0 * u * s1;
    const double truncation = truncated ? 2.0 * weight_sum[v] * tail_value : 0.0;
    const double direct_error = static_cast<double>(nonzero + 2) * u * max_y + truncation +
                                kernel_rel * max_y;
    const double hull_term = 1024.0 * u * max_y + 4.0 * weight_sum[v] * tail_value;
    const double discrepancy = 2.0 * fft_error + direct_error + 8.0 * u * max_y + hull_term;
    const double tolerance = std::max(2.0 * discrepancy, kModeGuardRelative * best);
    if (!(best - second > tolerance)) return false;
    index[v] = best_index;
  }
  return true;
}

bool mode_try_direct(std::size_t count) {
  return static_cast<double>(count) <= mode_direct_max.load(std::memory_order_relaxed);
}

bool mode_try_hull(std::size_t count) {
  return static_cast<double>(count) <= mode_hull_max.load(std::memory_order_relaxed);
}

double mode_point_r_density(const double* values, const double* weights,
                            std::size_t count, double phi) {
  if (count == 0) return NA_VALUE;
  for (std::size_t i = 0; i < count; ++i) {
    if (!std::isfinite(values[i]) || !std::isfinite(weights[i]) || weights[i] < 0.0) return NA_VALUE;
  }
  const ModeGrid grid = mode_grid(values, count, phi, mode_workspace().scratch);
  int index = 0;
  bool done = false;
  if (mode_try_direct(count)) {
    ++mode_direct_draws;
    if (mode_try_hull(count)) {
      ++mode_hull_draws;
      done = mode_index_hull(values, &weights, 1, count, grid, &index);
      if (!done) ++mode_hull_fallback_draws;
    }
    if (!done) done = mode_index_direct(values, &weights, 1, count, grid, &index);
    if (!done) ++mode_guard_draws;
  }
  if (!done) index = mode_index_fft(values, weights, count, grid);
  return grid.from + grid.output_step * static_cast<double>(index);
}

std::pair<double, double> mode_point_r_density_pair(const double* values,
                                                      const double* simple_weights,
                                                      const double* weighted_weights,
                                                      std::size_t count,
                                                      double phi) {
  if (count == 0) return std::make_pair(NA_VALUE, NA_VALUE);
  for (std::size_t i = 0; i < count; ++i) {
    if (!std::isfinite(values[i]) || !std::isfinite(simple_weights[i]) ||
        !std::isfinite(weighted_weights[i]) || simple_weights[i] < 0.0 ||
        weighted_weights[i] < 0.0) return std::make_pair(NA_VALUE, NA_VALUE);
  }
  const ModeGrid grid = mode_grid(values, count, phi, mode_workspace().scratch);
  int index[2] = {0, 0};
  bool done = false;
  if (mode_try_direct(count)) {
    const double* const weights[2] = {simple_weights, weighted_weights};
    ++mode_direct_draws;
    if (mode_try_hull(count)) {
      ++mode_hull_draws;
      done = mode_index_hull(values, weights, 2, count, grid, index);
      if (!done) ++mode_hull_fallback_draws;
    }
    if (!done) done = mode_index_direct(values, weights, 2, count, grid, index);
    if (!done) ++mode_guard_draws;
  }
  if (!done) {
    const std::pair<int, int> fft = mode_index_fft_pair(values, simple_weights, weighted_weights, count, grid);
    index[0] = fft.first;
    index[1] = fft.second;
  }
  return std::make_pair(grid.from + grid.output_step * static_cast<double>(index[0]),
                        grid.from + grid.output_step * static_cast<double>(index[1]));
}

struct Result {
  std::string method;
  int n = 0;
  double beta = NA_VALUE;
  double se = NA_VALUE;
  double pval = NA_VALUE;
  bool ratio_se_mean = false;
  double ratio_se_mean_value = NA_VALUE;
  bool bootstrap = false;
  int bootstrap_value = 0;
  bool phi = false;
  double phi_value = 1.0;
  bool q = false;
  double q_value = NA_VALUE;
  int q_df = 0;
  double q_pval = NA_VALUE;
  bool sigma = false;
  double sigma_value = NA_VALUE;
  bool intercept = false;
  double intercept_value = NA_VALUE;
  double intercept_se = NA_VALUE;
  double intercept_pval = NA_VALUE;
  int flipped = 0;
  double se_exposure_mean = NA_VALUE;
};

Result empty_result(const std::string& method, int n) {
  Result result;
  result.method = method;
  result.n = n;
  return result;
}

void populate_result_pvalues(Result& result) {
  if (result.method == "egger_bootstrap") {
    return;
  }
  if (result.method == "sign") {
    if (std::isfinite(result.beta) && result.n > 0) {
      const int concordant = static_cast<int>(std::llround(
        0.5 * (result.beta + 1.0) * static_cast<double>(result.n)));
      const int lower = std::min(concordant, result.n - concordant);
      double pval = 2.0 * R::pbinom(static_cast<double>(lower),
                                    static_cast<double>(result.n), 0.5, true, false);
      result.pval = std::min(1.0, pval);
    }
  } else if (result.method == "egger") {
    result.pval = t_pvalue(safe_statistic(result.beta, result.se), result.n - 2);
  } else if (result.method == "simple_mode" || result.method == "weighted_mode") {
    result.pval = t_pvalue(safe_statistic(result.beta, result.se), result.n - 1);
  } else {
    result.pval = z_pvalue(safe_statistic(result.beta, result.se));
  }
  if (result.q) result.q_pval = chi_square_pvalue(result.q_value, result.q_df);
  if (result.intercept) {
    result.intercept_pval = t_pvalue(
      safe_statistic(result.intercept_value, result.intercept_se), result.n - 2);
  }
}

struct Prepared {
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> sx;
  std::vector<double> sy;
  std::vector<double> ratio;
  std::vector<double> ratio_se;
  std::vector<double> bootstrap;
  std::vector<double> penalised_bootstrap;
  std::vector<double> mode_bootstrap;
  std::vector<double> egger_x_bootstrap;
  std::vector<double> egger_y_bootstrap;
};

void prepare_ratios(Prepared& p) {
  p.ratio.clear();
  p.ratio_se.clear();
  p.ratio.reserve(p.x.size());
  p.ratio_se.reserve(p.x.size());
  for (std::size_t i = 0; i < p.x.size(); ++i) {
    if (p.x[i] == 0.0) continue;
    const double x = p.x[i];
    const double y = p.y[i];
    const double sx = p.sx[i];
    const double sy = p.sy[i];
    p.ratio.push_back(y / x);
    const double outcome_part = sy / x;
    const double exposure_part = y * sx / (x * x);
    p.ratio_se.push_back(std::sqrt(outcome_part * outcome_part +
                                   exposure_part * exposure_part));
  }
}

// IVW standard errors as TwoSampleMR: multiplicative random effects
// se / min(1, sigma) (mr_ivw) and fixed effects se / sigma (mr_ivw_fe), where
// se = base_se * sigma is the residual-scaled lm() standard error.  Both
// reduce to the fixed-effect base_se when the residual variance is under-
// dispersed; for an exact fit (sigma == 0, e.g. a self-pair with outcome ==
// exposure) the ratio is 0/0, and its limit, base_se, is returned rather than 0.
inline double ivw_mre_se(double base_se, double residual_se, double sigma) {
  if (!std::isfinite(sigma)) return residual_se;
  return sigma > 0.0 ? residual_se / std::min(1.0, sigma) : base_se;
}

inline double ivw_fe_se(double base_se, double residual_se, double sigma) {
  if (!std::isfinite(sigma)) return NA_VALUE;
  return sigma > 0.0 ? residual_se / sigma : base_se;
}

// One instrument (k = 1).  TwoSampleMR::mr() then reports only
// mr_wald_ratio: b = b_out / b_exp, se = se_out / |b_exp| (first order), with
// every other method NA.  The IVW estimators (ivw, ivw_fe, ivw_mre) return
// exactly that Wald ratio; Q and sigma are undefined (0 degrees of freedom)
// and stay NA.  This is distinct from an exact fit with k >= 2 (sigma == 0),
// which keeps the fixed-effect se (ivw_mre_se()/ivw_fe_se()).  NA when the
// exposure effect is 0 or a value is not finite.
inline bool single_snp_wald(double x, double y, double sy, double& beta, double& se) {
  beta = NA_VALUE;
  se = NA_VALUE;
  if (!std::isfinite(x) || x == 0.0 || !std::isfinite(y) || !std::isfinite(sy) || !(sy > 0.0)) {
    return false;
  }
  beta = y / x;
  se = sy / std::abs(x);
  if (!std::isfinite(beta) || !std::isfinite(se)) {
    beta = NA_VALUE;
    se = NA_VALUE;
    return false;
  }
  return true;
}

Result single_snp_ivw_result(const std::string& method, double x, double y, double sy) {
  Result result = empty_result(method, 1);
  double beta = NA_VALUE;
  double se = NA_VALUE;
  if (!single_snp_wald(x, y, sy, beta, se)) return result;
  result.beta = beta;
  result.se = se;
  result.pval = z_pvalue(safe_statistic(beta, se));
  return result;
}

Result compute_ivw(const Prepared& p, const std::string& method) {
  const int n = static_cast<int>(p.x.size());
  if (n == 1) return single_snp_ivw_result(method, p.x[0], p.y[0], p.sy[0]);
  if (n < 2) return empty_result(method, n);
  double denominator = 0.0;
  double numerator = 0.0;
  for (int i = 0; i < n; ++i) {
    const double weight = 1.0 / (p.sy[i] * p.sy[i]);
    denominator += weight * p.x[i] * p.x[i];
    numerator += weight * p.x[i] * p.y[i];
  }
  if (!(denominator > 0.0) || !std::isfinite(denominator)) return empty_result(method, n);
  const double beta = numerator / denominator;
  double rss = 0.0;
  for (int i = 0; i < n; ++i) {
    const double weight = 1.0 / (p.sy[i] * p.sy[i]);
    const double residual = p.y[i] - beta * p.x[i];
    rss += weight * residual * residual;
  }
  const int df = n - 1;
  const double sigma = std::sqrt(rss / static_cast<double>(df));
  const double base_se = std::sqrt(1.0 / denominator);
  const double residual_se = base_se * sigma;
  double se = residual_se;
  if (method == "ivw") {
    se = ivw_mre_se(base_se, residual_se, sigma);
  } else if (method == "ivw_fe") {
    se = ivw_fe_se(base_se, residual_se, sigma);
  }
  Result result = empty_result(method, n);
  result.beta = beta;
  result.se = se;
  result.pval = z_pvalue(safe_statistic(beta, se));
  result.q = true;
  result.q_value = rss;
  result.q_df = df;
  result.q_pval = chi_square_pvalue(rss, df);
  result.sigma = true;
  result.sigma_value = sigma;
  return result;
}

Result compute_uwr(const Prepared& p) {
  const int n = static_cast<int>(p.x.size());
  if (n < 2) return empty_result("uwr", n);
  double denominator = 0.0;
  double numerator = 0.0;
  for (int i = 0; i < n; ++i) {
    denominator += p.x[i] * p.x[i];
    numerator += p.x[i] * p.y[i];
  }
  if (!(denominator > 0.0) || !std::isfinite(denominator)) return empty_result("uwr", n);
  const double beta = numerator / denominator;
  double rss = 0.0;
  for (int i = 0; i < n; ++i) {
    const double residual = p.y[i] - beta * p.x[i];
    rss += residual * residual;
  }
  const int df = n - 1;
  const double sigma = std::sqrt(rss / static_cast<double>(df));
  const double residual_se = std::sqrt(1.0 / denominator) * sigma;
  const double se = sigma > 0.0 && std::isfinite(sigma)
    ? residual_se / std::min(1.0, sigma) : residual_se;
  Result result = empty_result("uwr", n);
  result.beta = beta;
  result.se = se;
  result.pval = z_pvalue(safe_statistic(beta, se));
  result.q = true;
  result.q_value = rss;
  result.q_df = df;
  result.q_pval = chi_square_pvalue(rss, df);
  result.sigma = true;
  result.sigma_value = sigma;
  return result;
}

Result compute_sign(const Prepared& p) {
  int n = 0;
  int concordant = 0;
  for (std::size_t i = 0; i < p.x.size(); ++i) {
    if (p.x[i] == 0.0 || p.y[i] == 0.0) continue;
    ++n;
    if ((p.x[i] > 0.0) == (p.y[i] > 0.0)) ++concordant;
  }
  if (n < 6) return empty_result("sign", n);
  const double beta = (2.0 * static_cast<double>(concordant) / static_cast<double>(n)) - 1.0;
  const int lower = std::min(concordant, n - concordant);
  double pval = r_math_deferred()
    ? NA_VALUE
    : 2.0 * R::pbinom(static_cast<double>(lower), static_cast<double>(n), 0.5, true, false);
  if (std::isfinite(pval)) pval = std::min(1.0, pval);
  Result result = empty_result("sign", n);
  result.beta = beta;
  result.pval = pval;
  return result;
}

Result compute_egger(const Prepared& p) {
  const int n = static_cast<int>(p.x.size());
  if (n < 3) return empty_result("egger", n);
  double sw = 0.0, swx = 0.0, swxx = 0.0, swy = 0.0, swxy = 0.0;
  double sx_sum = 0.0;
  int flipped = 0;
  std::vector<double> x(n), y(n), weights(n);
  for (int i = 0; i < n; ++i) {
    const double sign = p.x[i] == 0.0 || p.x[i] > 0.0 ? 1.0 : -1.0;
    if (sign < 0.0) ++flipped;
    x[i] = std::abs(p.x[i]);
    y[i] = p.y[i] * sign;
    weights[i] = 1.0 / (p.sy[i] * p.sy[i]);
    sw += weights[i];
    swx += weights[i] * x[i];
    swxx += weights[i] * x[i] * x[i];
    swy += weights[i] * y[i];
    swxy += weights[i] * x[i] * y[i];
    sx_sum += p.sx[i];
  }
  const double determinant = sw * swxx - swx * swx;
  if (determinant == 0.0 || !std::isfinite(determinant)) return empty_result("egger", n);
  const double intercept = (swxx * swy - swx * swxy) / determinant;
  const double beta = (sw * swxy - swx * swy) / determinant;
  const double cov00 = swxx / determinant;
  const double cov11 = sw / determinant;
  double rss = 0.0;
  for (int i = 0; i < n; ++i) {
    const double residual = y[i] - intercept - beta * x[i];
    rss += weights[i] * residual * residual;
  }
  const int df = n - 2;
  const double sigma = std::sqrt(rss / static_cast<double>(df));
  const double correction = std::isfinite(sigma) ? std::min(1.0, sigma) : NA_VALUE;
  const double intercept_se = correction > 0.0 ? std::sqrt(cov00) * sigma / correction : NA_VALUE;
  const double beta_se = correction > 0.0 ? std::sqrt(cov11) * sigma / correction : NA_VALUE;
  Result result = empty_result("egger", n);
  result.beta = beta;
  result.se = beta_se;
  result.pval = t_pvalue(safe_statistic(beta, beta_se), df);
  result.intercept = true;
  result.intercept_value = intercept;
  result.intercept_se = intercept_se;
  result.intercept_pval = t_pvalue(safe_statistic(intercept, intercept_se), df);
  result.flipped = flipped;
  result.se_exposure_mean = sx_sum / static_cast<double>(n);
  result.q = true;
  result.q_value = rss;
  result.q_df = df;
  result.q_pval = chi_square_pvalue(rss, df);
  result.sigma = true;
  result.sigma_value = sigma;
  return result;
}

Result compute_egger_bootstrap(const Prepared& p, int nboot) {
  const int n = static_cast<int>(p.x.size());
  if (n < 3 || nboot <= 0 || p.egger_x_bootstrap.empty()) {
    return empty_result("egger_bootstrap", n);
  }
  std::vector<double> betas(nboot), intercepts(nboot);
  std::vector<double> weights(n);
  for (int i = 0; i < n; ++i) weights[i] = 1.0 / (p.sy[i] * p.sy[i]);
  for (int draw = 0; draw < nboot; ++draw) {
    const double* xs = p.egger_x_bootstrap.data() + static_cast<std::size_t>(draw) * n;
    const double* ys = p.egger_y_bootstrap.data() + static_cast<std::size_t>(draw) * n;
    double mean_x = 0.0;
    double mean_y = 0.0;
    double mean_xw = 0.0;
    double mean_yw = 0.0;
    for (int i = 0; i < n; ++i) {
      const double sign = xs[i] >= 0.0 ? 1.0 : -1.0;
      const double x = std::abs(xs[i]);
      const double y = ys[i] * sign;
      mean_x += x;
      mean_y += y;
      mean_xw += x * weights[i];
      mean_yw += y * weights[i];
    }
    mean_x /= static_cast<double>(n);
    mean_y /= static_cast<double>(n);
    mean_xw /= static_cast<double>(n);
    mean_yw /= static_cast<double>(n);
    double covariance = 0.0;
    double variance = 0.0;
    for (int i = 0; i < n; ++i) {
      const double sign = xs[i] >= 0.0 ? 1.0 : -1.0;
      const double xw = std::abs(xs[i]) * weights[i];
      const double yw = ys[i] * sign * weights[i];
      covariance += (xw - mean_xw) * (yw - mean_yw);
      variance += (xw - mean_xw) * (xw - mean_xw);
    }
    const double beta = variance > 0.0 ? covariance / variance : NA_VALUE;
    betas[draw] = beta;
    intercepts[draw] = std::isfinite(beta) ? mean_y - mean_x * beta : NA_VALUE;
  }
  double beta_mean = 0.0;
  double intercept_mean = 0.0;
  int finite_count = 0;
  for (int draw = 0; draw < nboot; ++draw) {
    if (std::isfinite(betas[draw])) {
      beta_mean += betas[draw];
      ++finite_count;
    }
    if (std::isfinite(intercepts[draw])) intercept_mean += intercepts[draw];
  }
  if (finite_count == 0) return empty_result("egger_bootstrap", n);
  beta_mean /= static_cast<double>(finite_count);
  int intercept_count = 0;
  for (double value : intercepts) if (std::isfinite(value)) ++intercept_count;
  if (intercept_count > 0) intercept_mean /= static_cast<double>(intercept_count);
  std::vector<double> finite_betas;
  std::vector<double> finite_intercepts;
  finite_betas.reserve(betas.size());
  finite_intercepts.reserve(intercepts.size());
  for (double value : betas) if (std::isfinite(value)) finite_betas.push_back(value);
  for (double value : intercepts) if (std::isfinite(value)) finite_intercepts.push_back(value);
  const double beta_se = sample_std(finite_betas);
  const double intercept_se = sample_std(finite_intercepts);
  const int sign = beta_mean > 0.0 ? 1 : (beta_mean < 0.0 ? -1 : 0);
  int opposite = 0;
  for (double value : betas) if (std::isfinite(value) && sign * value < 0.0) ++opposite;
  Result result = empty_result("egger_bootstrap", n);
  result.beta = beta_mean;
  result.se = beta_se;
  result.pval = static_cast<double>(opposite) / static_cast<double>(nboot);
  result.intercept = true;
  result.intercept_value = intercept_mean;
  result.intercept_se = intercept_se;
  const int intercept_sign = intercept_mean > 0.0 ? 1 : (intercept_mean < 0.0 ? -1 : 0);
  int intercept_opposite = 0;
  for (double value : intercepts) if (std::isfinite(value) && intercept_sign * value < 0.0) ++intercept_opposite;
  result.intercept_pval = static_cast<double>(intercept_opposite) / static_cast<double>(nboot);
  result.bootstrap = true;
  result.bootstrap_value = nboot;
  return result;
}

Result compute_median(const Prepared& p, const std::string& method,
                      int nboot, bool weighted) {
  const int n = static_cast<int>(p.ratio.size());
  if (n < 3) return empty_result(method, n);
  std::vector<double> weights(n);
  if (weighted) {
    for (int i = 0; i < n; ++i) weights[i] = 1.0 / (p.ratio_se[i] * p.ratio_se[i]);
  } else {
    std::fill(weights.begin(), weights.end(), 1.0);
  }
  std::vector<std::size_t> order;
  const double beta = weighted_median_point_ptr(
    p.ratio.data(), weights.data(), p.ratio.size(), order);
  double se = NA_VALUE;
  if (nboot > 0 && !p.bootstrap.empty()) {
    std::vector<double> estimates(nboot);
    for (int draw = 0; draw < nboot; ++draw) {
      const double* row = p.bootstrap.data() + static_cast<std::size_t>(draw) * n;
      estimates[draw] = weighted_median_point_ptr(row, weights.data(), n, order);
    }
    se = sample_std(estimates);
  }
  Result result = empty_result(method, n);
  result.beta = beta;
  result.se = se;
  result.pval = z_pvalue(safe_statistic(beta, se));
  if (weighted) {
    result.ratio_se_mean = true;
    result.ratio_se_mean_value = std::accumulate(p.ratio_se.begin(), p.ratio_se.end(), 0.0) /
                                 static_cast<double>(p.ratio_se.size());
  }
  result.bootstrap = true;
  result.bootstrap_value = nboot;
  return result;
}

void compute_both_medians(const Prepared& p, int nboot,
                          Result& simple, Result& weighted) {
  const int n = static_cast<int>(p.ratio.size());
  if (n < 3) {
    simple = empty_result("simple_median", n);
    weighted = empty_result("weighted_median", n);
    return;
  }
  std::vector<double> simple_weights(n, 1.0);
  std::vector<double> weighted_weights(n);
  for (int i = 0; i < n; ++i) {
    weighted_weights[i] = 1.0 / (p.ratio_se[i] * p.ratio_se[i]);
  }
  std::vector<std::size_t> order;
  const double simple_beta = weighted_median_point_ptr(
    p.ratio.data(), simple_weights.data(), p.ratio.size(), order);
  const double weighted_beta = weighted_median_ordered(
    p.ratio.data(), weighted_weights.data(), p.ratio.size(), order);
  double simple_se = NA_VALUE;
  double weighted_se = NA_VALUE;
  if (nboot > 0 && !p.bootstrap.empty()) {
    std::vector<double> simple_estimates(nboot), weighted_estimates(nboot);
    for (int draw = 0; draw < nboot; ++draw) {
      const double* row = p.bootstrap.data() + static_cast<std::size_t>(draw) * n;
      // Both estimators sort the same bootstrap ratios; only their weights
      // differ, so reuse the permutation for the second weighted pass.
      const double simple_value = weighted_median_point_ptr(
        row, simple_weights.data(), n, order);
      simple_estimates[draw] = simple_value;
      weighted_estimates[draw] = weighted_median_ordered(
        row, weighted_weights.data(), n, order);
    }
    simple_se = sample_std(simple_estimates);
    weighted_se = sample_std(weighted_estimates);
  }
  simple = empty_result("simple_median", n);
  simple.beta = simple_beta;
  simple.se = simple_se;
  simple.pval = z_pvalue(safe_statistic(simple_beta, simple_se));
  simple.bootstrap = true;
  simple.bootstrap_value = nboot;
  weighted = empty_result("weighted_median", n);
  weighted.beta = weighted_beta;
  weighted.se = weighted_se;
  weighted.pval = z_pvalue(safe_statistic(weighted_beta, weighted_se));
  weighted.ratio_se_mean = true;
  weighted.ratio_se_mean_value = std::accumulate(
    p.ratio_se.begin(), p.ratio_se.end(), 0.0) / static_cast<double>(p.ratio_se.size());
  weighted.bootstrap = true;
  weighted.bootstrap_value = nboot;
}

Result compute_penalised_median(const Prepared& p, int nboot, double penk) {
  const int n = static_cast<int>(p.ratio.size());
  if (n < 3) return empty_result("penalised_weighted_median", n);
  std::vector<double> weights(n), penalised_weights(n);
  std::vector<std::size_t> order;
  for (int i = 0; i < n; ++i) weights[i] = 1.0 / (p.ratio_se[i] * p.ratio_se[i]);
  const double weighted_beta = weighted_median_point_ptr(
    p.ratio.data(), weights.data(), p.ratio.size(), order);
  for (int i = 0; i < n; ++i) {
    const double statistic = weights[i] * (p.ratio[i] - weighted_beta) *
                             (p.ratio[i] - weighted_beta);
    const double penalty = chi_square1_survival(statistic);
    penalised_weights[i] = weights[i] * std::min(1.0, penalty * penk);
  }
  const double beta = weighted_median_point_ptr(
    p.ratio.data(), penalised_weights.data(), p.ratio.size(), order);
  double se = NA_VALUE;
  const std::vector<double>& bootstrap = p.penalised_bootstrap.empty()
    ? p.bootstrap : p.penalised_bootstrap;
  if (nboot > 0 && !bootstrap.empty()) {
    std::vector<double> estimates(nboot);
    for (int draw = 0; draw < nboot; ++draw) {
      const double* row = bootstrap.data() + static_cast<std::size_t>(draw) * n;
      estimates[draw] = weighted_median_point_ptr(
        row, penalised_weights.data(), n, order);
    }
    se = sample_std(estimates);
  }
  Result result = empty_result("penalised_weighted_median", n);
  result.beta = beta;
  result.se = se;
  result.pval = z_pvalue(safe_statistic(beta, se));
  result.bootstrap = true;
  result.bootstrap_value = nboot;
  return result;
}

Result compute_mode(const Prepared& p, const std::string& method, int nboot, double phi) {
  const int n = static_cast<int>(p.ratio.size());
  if (n < 3) return empty_result(method, n);
  std::vector<double> point_se(n), weights(n);
  if (method == "simple_mode") {
    std::fill(point_se.begin(), point_se.end(), 1.0);
    std::fill(weights.begin(), weights.end(), 1.0);
  } else {
    point_se = p.ratio_se;
    for (int i = 0; i < n; ++i) weights[i] = 1.0 / (p.ratio_se[i] * p.ratio_se[i]);
  }
  const double beta = mode_point_r_density(p.ratio.data(), weights.data(), p.ratio.size(), phi);
  double se = NA_VALUE;
  if (nboot > 0 && !p.mode_bootstrap.empty()) {
    std::vector<double> estimates(nboot);
    for (int draw = 0; draw < nboot; ++draw) {
      estimates[draw] = mode_point_r_density(p.mode_bootstrap.data() + static_cast<std::size_t>(draw) * n, weights.data(), p.ratio.size(), phi);
    }
    se = mad(estimates);
  }
  flush_mode_counters();
  Result result = empty_result(method, n);
  result.beta = beta;
  result.se = se;
  result.pval = t_pvalue(safe_statistic(beta, se), n - 1);
  result.bootstrap = true;
  result.bootstrap_value = nboot;
  result.phi = true;
  result.phi_value = phi;
  return result;
}

void compute_both_modes(const Prepared& p, int nboot, double phi,
                        Result& simple, Result& weighted) {
  const int n = static_cast<int>(p.ratio.size());
  if (n < 3) {
    simple = empty_result("simple_mode", n);
    weighted = empty_result("weighted_mode", n);
    return;
  }
  std::vector<double> simple_weights(n, 1.0);
  std::vector<double> weighted_weights(n);
  for (int i = 0; i < n; ++i) weighted_weights[i] = 1.0 / (p.ratio_se[i] * p.ratio_se[i]);
  const std::pair<double, double> point_modes = mode_point_r_density_pair(
    p.ratio.data(), simple_weights.data(), weighted_weights.data(), p.ratio.size(), phi);
  const double simple_beta = point_modes.first;
  const double weighted_beta = point_modes.second;
  double simple_se = NA_VALUE;
  double weighted_se = NA_VALUE;
  if (nboot > 0 && !p.mode_bootstrap.empty()) {
    std::vector<double> simple_estimates(nboot), weighted_estimates(nboot);
    for (int draw = 0; draw < nboot; ++draw) {
      const double* row = p.mode_bootstrap.data() + static_cast<std::size_t>(draw) * n;
      const std::pair<double, double> modes = mode_point_r_density_pair(
        row, simple_weights.data(), weighted_weights.data(), p.ratio.size(), phi);
      simple_estimates[draw] = modes.first;
      weighted_estimates[draw] = modes.second;
    }
    simple_se = mad(simple_estimates);
    weighted_se = mad(weighted_estimates);
  }
  flush_mode_counters();
  simple = empty_result("simple_mode", n);
  simple.beta = simple_beta;
  simple.se = simple_se;
  simple.pval = t_pvalue(safe_statistic(simple_beta, simple_se), n - 1);
  simple.bootstrap = true;
  simple.bootstrap_value = nboot;
  simple.phi = true;
  simple.phi_value = phi;
  weighted = empty_result("weighted_mode", n);
  weighted.beta = weighted_beta;
  weighted.se = weighted_se;
  weighted.pval = t_pvalue(safe_statistic(weighted_beta, weighted_se), n - 1);
  weighted.bootstrap = true;
  weighted.bootstrap_value = nboot;
  weighted.phi = true;
  weighted.phi_value = phi;
}

Result compute_wald(const Prepared& p) {
  const int n = static_cast<int>(p.x.size());
  if (n != 1 || p.x.front() == 0.0) return empty_result("wald_ratio", n);
  const double beta = p.y.front() / p.x.front();
  const double se = p.sy.front() / std::abs(p.x.front());
  Result result = empty_result("wald_ratio", 1);
  result.beta = beta;
  result.se = se;
  result.pval = z_pvalue(safe_statistic(beta, se));
  return result;
}

// Standard-normal draw source. The default reads R's RNG (main thread only);
// PredrawnNormals replays draws made beforehand on the main thread in exactly the
// order the default would consume them, so worker threads never touch R.
struct RNormals {
  double operator()() const { return R::rnorm(0.0, 1.0); }
};

struct PredrawnNormals {
  const double* next;
  double operator()() { return *next++; }
};

// Number of standard-normal draws make_bootstrap() + make_mode_bootstrap()
// consume for a pair with `snps` usable rows and `ratios` non-zero exposures.
double bootstrap_draw_count(std::size_t snps, std::size_t ratios, int nboot,
                            bool needs_median, bool needs_egger,
                            bool needs_penalised, bool needs_mode) {
  if (nboot <= 0) return 0.0;
  double count = 0.0;
  if (snps >= 3 && (needs_median || needs_egger || needs_penalised)) {
    const double stream = 2.0 * static_cast<double>(nboot) * static_cast<double>(snps);
    count += stream;
    if (needs_penalised) count += stream;
  }
  if (needs_mode && ratios >= 3) {
    count += static_cast<double>(nboot) * static_cast<double>(ratios);
  }
  return count;
}

// Serial cost of a mode-bootstrap normal relative to a median/Egger stream
// normal (one density per draw vs one sort/regression): measured 12-45x at
// 3-30 ratios with the direct density.
constexpr double kModeDrawCost = 16.0;

// bootstrap_draw_count() with mode normals weighted by kModeDrawCost: the
// work estimate the thread cap uses for bootstrap batches.
double bootstrap_draw_work(std::size_t snps, std::size_t ratios, int nboot,
                           bool needs_median, bool needs_egger,
                           bool needs_penalised, bool needs_mode) {
  const double stream = bootstrap_draw_count(snps, ratios, nboot, needs_median,
                                             needs_egger, needs_penalised, false);
  const double mode = bootstrap_draw_count(snps, ratios, nboot, false, false,
                                           false, needs_mode);
  return stream + kModeDrawCost * mode;
}

template <typename Draw>
void make_bootstrap(Prepared& p, int nboot, Draw& next_normal,
                    bool needs_median, bool needs_egger,
                    bool needs_penalised) {
  if (nboot <= 0 || p.x.size() < 3 ||
      (!needs_median && !needs_egger && !needs_penalised)) return;
  const std::size_t n = p.ratio.size();
  auto draw_stream = [&](bool penalised) {
    if (penalised) {
      p.penalised_bootstrap.resize(static_cast<std::size_t>(nboot) * n);
    } else {
      if (needs_median) p.bootstrap.resize(static_cast<std::size_t>(nboot) * n);
      if (needs_egger) {
        p.egger_x_bootstrap.resize(static_cast<std::size_t>(nboot) * p.x.size());
        p.egger_y_bootstrap.resize(static_cast<std::size_t>(nboot) * p.x.size());
      }
    }
    std::vector<double> exp_z(static_cast<std::size_t>(nboot) * p.x.size());
    std::vector<double> out_z(static_cast<std::size_t>(nboot) * p.x.size());
    // R fills matrix(rnorm(nboot*n, mean=rep(x, each=nboot)), nrow=nboot)
    // column by column, so consume each RNG stream in SNP-major order.
    for (std::size_t snp = 0; snp < p.x.size(); ++snp) {
      for (int draw = 0; draw < nboot; ++draw) {
        exp_z[static_cast<std::size_t>(draw) * p.x.size() + snp] = next_normal();
      }
    }
    for (std::size_t snp = 0; snp < p.x.size(); ++snp) {
      for (int draw = 0; draw < nboot; ++draw) {
        out_z[static_cast<std::size_t>(draw) * p.x.size() + snp] = next_normal();
      }
    }
    for (int draw = 0; draw < nboot; ++draw) {
      std::size_t ratio_index = 0;
      for (std::size_t snp = 0; snp < p.x.size(); ++snp) {
        const double exp_draw = p.x[snp] + p.sx[snp] * exp_z[static_cast<std::size_t>(draw) * p.x.size() + snp];
        const double out_draw = p.y[snp] + p.sy[snp] * out_z[static_cast<std::size_t>(draw) * p.x.size() + snp];
        if (!penalised && needs_egger) {
          p.egger_x_bootstrap[static_cast<std::size_t>(draw) * p.x.size() + snp] = exp_draw;
          p.egger_y_bootstrap[static_cast<std::size_t>(draw) * p.x.size() + snp] = out_draw;
        }
        if (p.x[snp] == 0.0) continue;
        if (penalised) {
          p.penalised_bootstrap[static_cast<std::size_t>(draw) * n + ratio_index] =
            exp_draw == 0.0 ? NA_VALUE : out_draw / exp_draw;
        } else if (needs_median) {
          p.bootstrap[static_cast<std::size_t>(draw) * n + ratio_index] =
            exp_draw == 0.0 ? NA_VALUE : out_draw / exp_draw;
        }
        ++ratio_index;
      }
    }
  };
  // Native penalised weighted median first calls weighted median (and consumes
  // its full bootstrap stream), then bootstraps the penalised weights again.
  draw_stream(false);
  if (needs_penalised) draw_stream(true);
}

template <typename Draw>
void make_mode_bootstrap(Prepared& p, int nboot, Draw& next_normal) {
  if (nboot <= 0 || p.ratio.size() < 3) return;
  const std::size_t n = p.ratio.size();
  p.mode_bootstrap.resize(static_cast<std::size_t>(nboot) * n);
  // TwoSampleMR::mr_mode draws delta-method Wald ratios directly.
  for (std::size_t snp = 0; snp < n; ++snp) {
    for (int draw = 0; draw < nboot; ++draw) {
      p.mode_bootstrap[static_cast<std::size_t>(draw) * n + snp] =
        p.ratio[snp] + p.ratio_se[snp] * next_normal();
    }
  }
}

struct BootstrapNeeds {
  bool median = false;
  bool penalised = false;
  bool mode = false;
  bool egger = false;
};

BootstrapNeeds bootstrap_needs(const std::vector<std::string>& methods) {
  BootstrapNeeds needs;
  for (const std::string& method : methods) {
    needs.median = needs.median || method == "simple_median" ||
                   method == "weighted_median";
    needs.penalised = needs.penalised || method == "penalised_weighted_median";
    needs.mode = needs.mode || method == "simple_mode" || method == "weighted_mode";
    needs.egger = needs.egger || method == "egger_bootstrap";
  }
  return needs;
}

// `draws` == nullptr draws from R's RNG (main thread only); otherwise the
// bootstrap consumes pre-drawn standard normals from `draws` in the same order.
std::vector<Result> compute_pair(Prepared p,
                                 const std::vector<std::string>& methods,
                                 int nboot, SEXP seed,
                                 bool prepare_bootstrap = true,
                                 double phi = 1.0,
                                 double penk = 20.0,
                                 const double* draws = nullptr) {
  (void) seed;
  prepare_ratios(p);
  if (prepare_bootstrap) {
    const BootstrapNeeds needs = bootstrap_needs(methods);
    auto fill = [&](auto& source) {
      if (needs.median || needs.egger || needs.penalised) {
        make_bootstrap(p, nboot, source, needs.median, needs.egger, needs.penalised);
      }
      if (needs.mode) make_mode_bootstrap(p, nboot, source);
    };
    if (draws == nullptr) {
      RNormals source;
      fill(source);
    } else {
      PredrawnNormals source{draws};
      fill(source);
    }
  }
  std::vector<Result> result;
  result.resize(methods.size());
  bool has_simple = false;
  bool has_weighted = false;
  std::size_t simple_index = 0;
  std::size_t weighted_index = 0;
  bool has_simple_median = false;
  bool has_weighted_median = false;
  std::size_t simple_median_index = 0;
  std::size_t weighted_median_index = 0;
  for (std::size_t i = 0; i < methods.size(); ++i) {
    if (methods[i] == "simple_mode") { has_simple = true; simple_index = i; }
    if (methods[i] == "weighted_mode") { has_weighted = true; weighted_index = i; }
    if (methods[i] == "simple_median") {
      has_simple_median = true;
      simple_median_index = i;
    }
    if (methods[i] == "weighted_median") {
      has_weighted_median = true;
      weighted_median_index = i;
    }
  }
  for (std::size_t i = 0; i < methods.size(); ++i) {
    const std::string& method = methods[i];
    if ((method == "simple_mode" || method == "weighted_mode") && has_simple && has_weighted) continue;
    if ((method == "simple_median" || method == "weighted_median") &&
        has_simple_median && has_weighted_median) continue;
    if (method == "ivw" || method == "ivw_fe" || method == "ivw_mre") {
      result[i] = compute_ivw(p, method);
    } else if (method == "uwr") {
      result[i] = compute_uwr(p);
    } else if (method == "sign") {
      result[i] = compute_sign(p);
    } else if (method == "egger") {
      result[i] = compute_egger(p);
    } else if (method == "egger_bootstrap") {
      result[i] = compute_egger_bootstrap(p, nboot);
    } else if (method == "simple_median" || method == "weighted_median") {
      result[i] = compute_median(p, method, nboot, method == "weighted_median");
    } else if (method == "penalised_weighted_median") {
      result[i] = compute_penalised_median(p, nboot, penk);
    } else if (method == "simple_mode" || method == "weighted_mode") {
      result[i] = compute_mode(p, method, nboot, phi);
    } else if (method == "wald_ratio") {
      result[i] = compute_wald(p);
    }
  }
  if (has_simple && has_weighted) {
    compute_both_modes(p, nboot, phi, result[simple_index], result[weighted_index]);
  }
  if (has_simple_median && has_weighted_median) {
    compute_both_medians(p, nboot, result[simple_median_index], result[weighted_median_index]);
  }
  return result;
}

Rcpp::List result_to_list(const Result& result) {
  Rcpp::List out;
  out["method"] = result.method;
  out["n"] = result.n;
  out["beta"] = finite_or_na(result.beta);
  out["se"] = finite_or_na(result.se);
  out["pval"] = finite_or_na(result.pval);
  if (result.ratio_se_mean) out["ratio_se_mean"] = finite_or_na(result.ratio_se_mean_value);
  if (result.bootstrap) out["bootstrap"] = result.bootstrap_value;
  if (result.phi) out["phi"] = finite_or_na(result.phi_value);
  if (result.q) {
    out["Q"] = finite_or_na(result.q_value);
    out["Q_df"] = result.q_df;
    out["Q_pval"] = finite_or_na(result.q_pval);
  }
  if (result.sigma) out["sigma"] = finite_or_na(result.sigma_value);
  if (result.intercept) {
    out["intercept"] = finite_or_na(result.intercept_value);
    out["intercept_se"] = finite_or_na(result.intercept_se);
    out["intercept_pval"] = finite_or_na(result.intercept_pval);
    out["flipped"] = result.flipped;
    out["se_exposure_mean"] = finite_or_na(result.se_exposure_mean);
  }
  return out;
}

std::vector<std::string> parse_methods(Rcpp::CharacterVector methods) {
  const std::vector<std::string> allowed = {
    "ivw", "ivw_fe", "ivw_mre", "egger", "egger_bootstrap", "uwr", "sign",
    "simple_median", "weighted_median", "penalised_weighted_median",
    "simple_mode", "weighted_mode", "wald_ratio"
  };
  std::vector<std::string> parsed;
  parsed.reserve(methods.size());
  for (R_xlen_t i = 0; i < methods.size(); ++i) {
    if (methods[i] == NA_STRING) Rcpp::stop("methods cannot contain NA");
    const std::string value = Rcpp::as<std::string>(methods[i]);
    if (std::find(allowed.begin(), allowed.end(), value) == allowed.end()) {
      Rcpp::stop("unknown MR method: " + value);
    }
    parsed.push_back(value);
  }
  if (parsed.empty()) Rcpp::stop("methods must contain at least one method");
  return parsed;
}

void validate_controls(int nboot, int threads, double phi) {
  if (nboot < 0) Rcpp::stop("nboot must be a non-negative integer");
  if (threads < 1) Rcpp::stop("threads must be at least 1");
  if (!std::isfinite(phi) || phi <= 0.0) Rcpp::stop("phi must be positive and finite");
}

Prepared one_pair_from_vectors(Rcpp::NumericVector x, Rcpp::NumericVector y,
                               Rcpp::NumericVector sx, Rcpp::NumericVector sy) {
  if (x.size() != y.size() || x.size() != sx.size() || x.size() != sy.size()) {
    Rcpp::stop("MR vectors must have equal lengths");
  }
  Prepared p;
  p.x.reserve(x.size()); p.y.reserve(x.size()); p.sx.reserve(x.size()); p.sy.reserve(x.size());
  for (R_xlen_t i = 0; i < x.size(); ++i) {
    if (Rcpp::NumericVector::is_na(x[i]) || Rcpp::NumericVector::is_na(y[i]) ||
        Rcpp::NumericVector::is_na(sx[i]) || Rcpp::NumericVector::is_na(sy[i])) continue;
    if (!std::isfinite(x[i]) || !std::isfinite(y[i]) || !std::isfinite(sx[i]) ||
        !std::isfinite(sy[i]) || sx[i] <= 0.0 || sy[i] <= 0.0) continue;
    p.x.push_back(x[i]); p.y.push_back(y[i]); p.sx.push_back(sx[i]); p.sy.push_back(sy[i]);
  }
  return p;
}

void validate_grid_shapes(Rcpp::NumericMatrix exp_beta, Rcpp::NumericMatrix out_beta,
                          Rcpp::NumericMatrix exp_se, Rcpp::NumericMatrix out_se) {
  if (exp_beta.nrow() == 0 || out_beta.nrow() == 0 || exp_beta.ncol() == 0) {
    Rcpp::stop("grid matrices must have matching non-empty SNP dimensions");
  }
  if (exp_beta.nrow() != exp_se.nrow() || exp_beta.ncol() != exp_se.ncol() ||
      out_beta.nrow() != out_se.nrow() || out_beta.ncol() != out_se.ncol() ||
      exp_beta.ncol() != out_beta.ncol()) {
    Rcpp::stop("grid matrices must have matching beta/se and SNP dimensions");
  }
}

struct GridData {
  int exposure_count = 0;
  int outcome_count = 0;
  int snp_count = 0;
  std::vector<double> exp_beta, out_beta, exp_se, out_se;
  std::vector<double> exp_draws;
  std::vector<double> exp_inverse, out_draws;
  std::vector<double> exp_inverse_penalised, out_draws_penalised;
  std::vector<double> mode_z;
};

double stable_grid_rss(const GridData& grid, int exposure, int outcome,
                       double beta) {
  const std::size_t exposure_offset =
    static_cast<std::size_t>(exposure) * grid.snp_count;
  const std::size_t outcome_offset =
    static_cast<std::size_t>(outcome) * grid.snp_count;
  long double sum = 0.0L;
  long double correction = 0.0L;
  const long double coefficient = static_cast<long double>(beta);
  for (int snp = 0; snp < grid.snp_count; ++snp) {
    const std::size_t exposure_index = exposure_offset + snp;
    const std::size_t outcome_index = outcome_offset + snp;
    const long double residual =
      static_cast<long double>(grid.out_beta[outcome_index]) -
      coefficient * static_cast<long double>(grid.exp_beta[exposure_index]);
    const long double standard_error =
      static_cast<long double>(grid.out_se[outcome_index]);
    const long double term =
      (residual / standard_error) * (residual / standard_error);
    const long double adjusted = term - correction;
    const long double updated = sum + adjusted;
    correction = (updated - sum) - adjusted;
    sum = updated;
  }
  return static_cast<double>(sum);
}

bool only_ivw_methods(const std::vector<std::string>& methods) {
  if (methods.empty()) return false;
  for (const std::string& method : methods) {
    if (method != "ivw" && method != "ivw_fe" && method != "ivw_mre") return false;
  }
  return true;
}

std::vector<Result> compute_ivw_grid_blas(
    const GridData& grid, const std::vector<std::string>& methods) {
  const int exposure_count = grid.exposure_count;
  const int outcome_count = grid.outcome_count;
  const int snp_count = grid.snp_count;
  const std::size_t pair_count = static_cast<std::size_t>(exposure_count) *
                                 static_cast<std::size_t>(outcome_count);
  const std::size_t method_count = methods.size();
  std::vector<Result> results(pair_count * method_count);
  if (snp_count == 1) {
    // One instrument: the Wald ratio (see single_snp_wald()).
    for (int exposure = 0; exposure < exposure_count; ++exposure) {
      for (int outcome = 0; outcome < outcome_count; ++outcome) {
        const std::size_t pair = static_cast<std::size_t>(exposure) * outcome_count + outcome;
        for (std::size_t i = 0; i < methods.size(); ++i) {
          results[pair * method_count + i] = single_snp_ivw_result(
            methods[i], grid.exp_beta[static_cast<std::size_t>(exposure)],
            grid.out_beta[static_cast<std::size_t>(outcome)],
            grid.out_se[static_cast<std::size_t>(outcome)]);
        }
      }
    }
    return results;
  }
  if (snp_count < 2) {
    for (std::size_t pair = 0; pair < pair_count; ++pair) {
      for (std::size_t i = 0; i < methods.size(); ++i) {
        results[pair * method_count + i] = empty_result(methods[i], snp_count);
      }
    }
    return results;
  }

  // Store the exposure matrices as nSNP x nExposure and the outcome matrices
  // as nSNP x nOutcome in memory. The row-major pair layout already has this
  // byte order, so BLAS can consume the contiguous vectors without another
  // transpose or R-to-C++ conversion.
  const std::size_t exp_size = static_cast<std::size_t>(snp_count) * exposure_count;
  const std::size_t out_size = static_cast<std::size_t>(snp_count) * outcome_count;
  std::vector<double> exp_squared(exp_size);
  std::vector<double> out_weight(out_size);
  std::vector<double> out_weighted(out_size);
  std::vector<double> out_y2(static_cast<std::size_t>(outcome_count), 0.0);
  for (int exposure = 0; exposure < exposure_count; ++exposure) {
    const std::size_t offset = static_cast<std::size_t>(exposure) * snp_count;
    for (int snp = 0; snp < snp_count; ++snp) {
      const double x = grid.exp_beta[offset + static_cast<std::size_t>(snp)];
      exp_squared[offset + static_cast<std::size_t>(snp)] = x * x;
    }
  }
  for (int outcome = 0; outcome < outcome_count; ++outcome) {
    const std::size_t offset = static_cast<std::size_t>(outcome) * snp_count;
    double y2 = 0.0;
    for (int snp = 0; snp < snp_count; ++snp) {
      const std::size_t index = offset + static_cast<std::size_t>(snp);
      const double weight = 1.0 / (grid.out_se[index] * grid.out_se[index]);
      const double y = grid.out_beta[index];
      out_weight[index] = weight;
      out_weighted[index] = weight * y;
      y2 += weight * y * y;
    }
    out_y2[static_cast<std::size_t>(outcome)] = y2;
  }

  // C is exposure-major in its mathematical shape but column-major in memory:
  // C[exposure + exposure_count * outcome]. This is the same order as the
  // public grid result list (exposure-major, outcome-minor) after indexing.
  std::vector<double> numerator(pair_count, 0.0);
  std::vector<double> denominator(pair_count, 0.0);
  const char transposed = 'T';
  const char normal = 'N';
  const int m = exposure_count;
  const int n = outcome_count;
  const int k = snp_count;
  const int lda = snp_count;
  const int ldb = snp_count;
  const int ldc = exposure_count;
  const double alpha = 1.0;
  const double beta_zero = 0.0;
  F77_CALL(dgemm)(&transposed, &normal, &m, &n, &k, &alpha,
                  grid.exp_beta.data(), &lda, out_weighted.data(), &ldb,
                  &beta_zero, numerator.data(), &ldc FCONE FCONE);
  F77_CALL(dgemm)(&transposed, &normal, &m, &n, &k, &alpha,
                  exp_squared.data(), &lda, out_weight.data(), &ldb,
                  &beta_zero, denominator.data(), &ldc FCONE FCONE);

  for (int exposure = 0; exposure < exposure_count; ++exposure) {
    for (int outcome = 0; outcome < outcome_count; ++outcome) {
      const std::size_t matrix_index = static_cast<std::size_t>(exposure) +
                                       static_cast<std::size_t>(exposure_count) * outcome;
      const double den = denominator[matrix_index];
      const double num = numerator[matrix_index];
      const double beta_value = den > 0.0 && std::isfinite(den) ? num / den : NA_VALUE;
      double rss = NA_VALUE;
      double sigma = NA_VALUE;
      double base_se = NA_VALUE;
      if (std::isfinite(beta_value)) {
        const double y2 = out_y2[static_cast<std::size_t>(outcome)];
        rss = y2 - 2.0 * beta_value * num + beta_value * beta_value * den;
        // The normal-equation subtraction is fast but loses precision for
        // nearly exact fits. Recompute only cancellation-prone pairs from
        // their residuals, using compensated long-double accumulation.
        const double cancellation_scale = std::max({
          1.0, std::abs(y2), std::abs(2.0 * beta_value * num),
          std::abs(beta_value * beta_value * den)
        });
        if (rss < 0.0 || std::abs(rss) <= 1e-10 * cancellation_scale) {
          rss = stable_grid_rss(grid, exposure, outcome, beta_value);
        }
        if (rss >= 0.0 && std::isfinite(rss)) {
          sigma = std::sqrt(rss / static_cast<double>(snp_count - 1));
          base_se = std::sqrt(1.0 / den);
        }
      }
      const double residual_se = std::isfinite(sigma) && std::isfinite(base_se)
        ? base_se * sigma : NA_VALUE;
      const std::size_t result_offset =
        (static_cast<std::size_t>(exposure) * outcome_count + outcome) * method_count;
      for (std::size_t i = 0; i < methods.size(); ++i) {
        const std::string& method = methods[i];
        Result result = empty_result(method, snp_count);
        result.beta = beta_value;
        if (method == "ivw") {
          result.se = ivw_mre_se(base_se, residual_se, sigma);
        } else if (method == "ivw_fe") {
          result.se = ivw_fe_se(base_se, residual_se, sigma);
        } else {
          result.se = residual_se;
        }
        result.pval = z_pvalue(safe_statistic(result.beta, result.se));
        result.q = true;
        result.q_value = rss;
        result.q_df = snp_count - 1;
        result.q_pval = chi_square_pvalue(rss, snp_count - 1);
        result.sigma = true;
        result.sigma_value = sigma;
        results[result_offset + i] = result;
      }
    }
  }
  return results;
}

Rcpp::List compute_ivw_grid_compact(const GridData& grid,
                                    const std::vector<std::string>& methods) {
  const int exposure_count = grid.exposure_count;
  const int outcome_count = grid.outcome_count;
  const int snp_count = grid.snp_count;
  const std::size_t pair_count = static_cast<std::size_t>(exposure_count) *
                                 static_cast<std::size_t>(outcome_count);
  const int method_count = static_cast<int>(methods.size());
  Rcpp::NumericMatrix beta(method_count, pair_count);
  Rcpp::NumericMatrix se(method_count, pair_count);
  Rcpp::NumericMatrix pval(method_count, pair_count);
  Rcpp::NumericMatrix q(method_count, pair_count);
  Rcpp::NumericMatrix q_df(method_count, pair_count);
  Rcpp::NumericMatrix q_pval(method_count, pair_count);
  Rcpp::NumericMatrix sigma(method_count, pair_count);
  std::fill(beta.begin(), beta.end(), NA_VALUE);
  std::fill(se.begin(), se.end(), NA_VALUE);
  std::fill(pval.begin(), pval.end(), NA_VALUE);
  std::fill(q.begin(), q.end(), NA_VALUE);
  std::fill(q_df.begin(), q_df.end(), NA_VALUE);
  std::fill(q_pval.begin(), q_pval.end(), NA_VALUE);
  std::fill(sigma.begin(), sigma.end(), NA_VALUE);

  if (snp_count == 1) {
    // One instrument: the Wald ratio (see single_snp_wald()); Q and sigma NA.
    for (int exposure = 0; exposure < exposure_count; ++exposure) {
      for (int outcome = 0; outcome < outcome_count; ++outcome) {
        const std::size_t pair = static_cast<std::size_t>(exposure) * outcome_count + outcome;
        double beta_value = NA_VALUE;
        double se_value = NA_VALUE;
        single_snp_wald(grid.exp_beta[static_cast<std::size_t>(exposure)],
                        grid.out_beta[static_cast<std::size_t>(outcome)],
                        grid.out_se[static_cast<std::size_t>(outcome)], beta_value, se_value);
        const double p_value = z_pvalue(safe_statistic(beta_value, se_value));
        for (int method_index = 0; method_index < method_count; ++method_index) {
          beta(method_index, pair) = beta_value;
          se(method_index, pair) = se_value;
          pval(method_index, pair) = p_value;
        }
      }
    }
  }

  if (snp_count >= 2) {
    const std::size_t exp_size = static_cast<std::size_t>(snp_count) * exposure_count;
    const std::size_t out_size = static_cast<std::size_t>(snp_count) * outcome_count;
    std::vector<double> exp_squared(exp_size);
    std::vector<double> out_weight(out_size);
    std::vector<double> out_weighted(out_size);
    std::vector<double> out_y2(static_cast<std::size_t>(outcome_count), 0.0);
    for (int exposure = 0; exposure < exposure_count; ++exposure) {
      const std::size_t offset = static_cast<std::size_t>(exposure) * snp_count;
      for (int snp = 0; snp < snp_count; ++snp) {
        const double x = grid.exp_beta[offset + static_cast<std::size_t>(snp)];
        exp_squared[offset + static_cast<std::size_t>(snp)] = x * x;
      }
    }
    for (int outcome = 0; outcome < outcome_count; ++outcome) {
      const std::size_t offset = static_cast<std::size_t>(outcome) * snp_count;
      double y2 = 0.0;
      for (int snp = 0; snp < snp_count; ++snp) {
        const std::size_t index = offset + static_cast<std::size_t>(snp);
        const double weight = 1.0 / (grid.out_se[index] * grid.out_se[index]);
        const double y = grid.out_beta[index];
        out_weight[index] = weight;
        out_weighted[index] = weight * y;
        y2 += weight * y * y;
      }
      out_y2[static_cast<std::size_t>(outcome)] = y2;
    }

    std::vector<double> numerator(pair_count, 0.0);
    std::vector<double> denominator(pair_count, 0.0);
    const char transposed = 'T';
    const char normal = 'N';
    const int m = exposure_count;
    const int n = outcome_count;
    const int k = snp_count;
    const int lda = snp_count;
    const int ldb = snp_count;
    const int ldc = exposure_count;
    const double alpha = 1.0;
    const double beta_zero = 0.0;
    F77_CALL(dgemm)(&transposed, &normal, &m, &n, &k, &alpha,
                    grid.exp_beta.data(), &lda, out_weighted.data(), &ldb,
                    &beta_zero, numerator.data(), &ldc FCONE FCONE);
    F77_CALL(dgemm)(&transposed, &normal, &m, &n, &k, &alpha,
                    exp_squared.data(), &lda, out_weight.data(), &ldb,
                    &beta_zero, denominator.data(), &ldc FCONE FCONE);

    for (int exposure = 0; exposure < exposure_count; ++exposure) {
      for (int outcome = 0; outcome < outcome_count; ++outcome) {
        const std::size_t pair = static_cast<std::size_t>(exposure) * outcome_count + outcome;
        const std::size_t matrix_index = static_cast<std::size_t>(exposure) +
                                         static_cast<std::size_t>(exposure_count) * outcome;
        const double den = denominator[matrix_index];
        const double num = numerator[matrix_index];
        const double beta_value = den > 0.0 && std::isfinite(den) ? num / den : NA_VALUE;
        double rss = NA_VALUE;
        double sigma_value = NA_VALUE;
        double base_se = NA_VALUE;
        if (std::isfinite(beta_value)) {
          const double y2 = out_y2[static_cast<std::size_t>(outcome)];
          rss = y2 - 2.0 * beta_value * num + beta_value * beta_value * den;
          const double cancellation_scale = std::max({
            1.0, std::abs(y2), std::abs(2.0 * beta_value * num),
            std::abs(beta_value * beta_value * den)
          });
          if (rss < 0.0 || std::abs(rss) <= 1e-10 * cancellation_scale) {
            rss = stable_grid_rss(grid, exposure, outcome, beta_value);
          }
          if (rss >= 0.0 && std::isfinite(rss)) {
            sigma_value = std::sqrt(rss / static_cast<double>(snp_count - 1));
            base_se = std::sqrt(1.0 / den);
          }
        }
        const double residual_se = std::isfinite(sigma_value) && std::isfinite(base_se)
          ? base_se * sigma_value : NA_VALUE;
        for (int method_index = 0; method_index < method_count; ++method_index) {
          const std::string& method = methods[static_cast<std::size_t>(method_index)];
          double method_se = residual_se;
          if (method == "ivw") {
            method_se = ivw_mre_se(base_se, residual_se, sigma_value);
          } else if (method == "ivw_fe") {
            method_se = ivw_fe_se(base_se, residual_se, sigma_value);
          }
          beta(method_index, pair) = beta_value;
          se(method_index, pair) = method_se;
          pval(method_index, pair) = z_pvalue(safe_statistic(beta_value, method_se));
          q(method_index, pair) = rss;
          q_df(method_index, pair) = snp_count - 1;
          q_pval(method_index, pair) = chi_square_pvalue(rss, snp_count - 1);
          sigma(method_index, pair) = sigma_value;
        }
      }
    }
  }

  Rcpp::CharacterVector method_codes(method_count);
  for (int i = 0; i < method_count; ++i) method_codes[i] = methods[static_cast<std::size_t>(i)];
  Rcpp::List output = Rcpp::List::create(
    Rcpp::_["methods"] = method_codes,
    Rcpp::_["n"] = snp_count,
    Rcpp::_["beta"] = beta,
    Rcpp::_["se"] = se,
    Rcpp::_["pval"] = pval,
    Rcpp::_["Q"] = q,
    Rcpp::_["Q_df"] = q_df,
    Rcpp::_["Q_pval"] = q_pval,
    Rcpp::_["sigma"] = sigma
  );
  output.attr("class") = "fastmr_ivw_compact";
  return output;
}

Rcpp::List compute_masked_ivw_grid(
    Rcpp::NumericMatrix exp_beta, Rcpp::NumericMatrix out_beta,
    Rcpp::NumericMatrix out_se, Rcpp::LogicalMatrix exp_present,
    Rcpp::LogicalMatrix out_present) {
  const int exposure_count = exp_beta.nrow();
  const int outcome_count = out_beta.nrow();
  const int snp_count = exp_beta.ncol();
  if (exposure_count == 0 || outcome_count == 0 || snp_count == 0 ||
      exp_present.nrow() != exposure_count || exp_present.ncol() != snp_count ||
      out_present.nrow() != outcome_count || out_present.ncol() != snp_count ||
      out_se.nrow() != outcome_count || out_se.ncol() != snp_count) {
    Rcpp::stop("masked IVW matrices must have matching non-empty dimensions");
  }

  const std::size_t exp_size = static_cast<std::size_t>(exposure_count) * snp_count;
  const std::size_t out_size = static_cast<std::size_t>(outcome_count) * snp_count;
  std::vector<double> exp_values(exp_size, 0.0);
  std::vector<double> exp_squared(exp_size, 0.0);
  std::vector<double> exp_presence(exp_size, 0.0);
  std::vector<double> out_weight(out_size, 0.0);
  std::vector<double> out_weighted(out_size, 0.0);
  std::vector<double> out_y2_weighted(out_size, 0.0);
  std::vector<double> out_valid(out_size, 0.0);

  for (int exposure = 0; exposure < exposure_count; ++exposure) {
    for (int snp = 0; snp < snp_count; ++snp) {
      const std::size_t index = static_cast<std::size_t>(exposure) * snp_count + snp;
      if (!exp_present(exposure, snp)) continue;
      const double value = exp_beta(exposure, snp);
      if (!std::isfinite(value)) Rcpp::stop("present exposure values must be finite");
      exp_values[index] = value;
      exp_squared[index] = value * value;
      exp_presence[index] = 1.0;
    }
  }
  for (int outcome = 0; outcome < outcome_count; ++outcome) {
    for (int snp = 0; snp < snp_count; ++snp) {
      const std::size_t index = static_cast<std::size_t>(outcome) * snp_count + snp;
      if (!out_present(outcome, snp)) continue;
      const double value = out_beta(outcome, snp);
      const double se = out_se(outcome, snp);
      if (!std::isfinite(value) || !std::isfinite(se) || se <= 0.0) {
        Rcpp::stop("present outcome values must be finite with positive standard errors");
      }
      const double weight = 1.0 / (se * se);
      out_weight[index] = weight;
      out_weighted[index] = weight * value;
      out_y2_weighted[index] = weight * value * value;
      out_valid[index] = 1.0;
    }
  }

  const std::size_t pair_count = static_cast<std::size_t>(exposure_count) * outcome_count;
  std::vector<double> numerator(pair_count, 0.0);
  std::vector<double> denominator(pair_count, 0.0);
  std::vector<double> yy(pair_count, 0.0);
  std::vector<double> nsnp(pair_count, 0.0);
  const char transposed = 'T';
  const char normal = 'N';
  const int m = exposure_count;
  const int n = outcome_count;
  const int k = snp_count;
  const int lda = snp_count;
  const int ldb = snp_count;
  const int ldc = exposure_count;
  const double alpha = 1.0;
  const double beta_zero = 0.0;
  F77_CALL(dgemm)(&transposed, &normal, &m, &n, &k, &alpha,
                  exp_values.data(), &lda, out_weighted.data(), &ldb,
                  &beta_zero, numerator.data(), &ldc FCONE FCONE);
  F77_CALL(dgemm)(&transposed, &normal, &m, &n, &k, &alpha,
                  exp_squared.data(), &lda, out_weight.data(), &ldb,
                  &beta_zero, denominator.data(), &ldc FCONE FCONE);
  F77_CALL(dgemm)(&transposed, &normal, &m, &n, &k, &alpha,
                  exp_presence.data(), &lda, out_y2_weighted.data(), &ldb,
                  &beta_zero, yy.data(), &ldc FCONE FCONE);
  F77_CALL(dgemm)(&transposed, &normal, &m, &n, &k, &alpha,
                  exp_presence.data(), &lda, out_valid.data(), &ldb,
                  &beta_zero, nsnp.data(), &ldc FCONE FCONE);

  Rcpp::NumericMatrix result_beta(exposure_count, outcome_count);
  Rcpp::NumericMatrix result_se(exposure_count, outcome_count);
  Rcpp::NumericMatrix result_q(exposure_count, outcome_count);
  Rcpp::NumericMatrix result_sigma(exposure_count, outcome_count);
  Rcpp::NumericMatrix result_nsnp(exposure_count, outcome_count);
  std::fill(result_beta.begin(), result_beta.end(), NA_VALUE);
  std::fill(result_se.begin(), result_se.end(), NA_VALUE);
  std::fill(result_q.begin(), result_q.end(), NA_VALUE);
  std::fill(result_sigma.begin(), result_sigma.end(), NA_VALUE);
  std::fill(result_nsnp.begin(), result_nsnp.end(), 0.0);
  // Present SNPs per exposure, built only if some pair has one instrument.
  std::vector<std::vector<int>> exposure_snps;
  for (int exposure = 0; exposure < exposure_count; ++exposure) {
    for (int outcome = 0; outcome < outcome_count; ++outcome) {
      const std::size_t index = static_cast<std::size_t>(exposure) +
                                static_cast<std::size_t>(exposure_count) * outcome;
      const double count = nsnp[index];
      result_nsnp(exposure, outcome) = count;
      if (count == 1.0) {
        // One instrument: the Wald ratio (see single_snp_wald()); Q, sigma NA.
        if (exposure_snps.empty()) {
          exposure_snps.resize(static_cast<std::size_t>(exposure_count));
          for (int e = 0; e < exposure_count; ++e) {
            for (int snp = 0; snp < snp_count; ++snp) {
              if (exp_presence[static_cast<std::size_t>(e) * snp_count + snp] != 0.0) {
                exposure_snps[static_cast<std::size_t>(e)].push_back(snp);
              }
            }
          }
        }
        for (int snp : exposure_snps[static_cast<std::size_t>(exposure)]) {
          const std::size_t out_index = static_cast<std::size_t>(outcome) * snp_count + snp;
          if (out_valid[out_index] == 0.0) continue;
          double beta_value = NA_VALUE;
          double se_value = NA_VALUE;
          if (single_snp_wald(exp_values[static_cast<std::size_t>(exposure) * snp_count + snp],
                              out_beta(outcome, snp), out_se(outcome, snp), beta_value, se_value)) {
            result_beta(exposure, outcome) = beta_value;
            result_se(exposure, outcome) = se_value;
          }
          break;
        }
        continue;
      }
      if (count < 2.0 || !std::isfinite(denominator[index]) || denominator[index] <= 0.0) continue;
      const double beta_value = numerator[index] / denominator[index];
      if (!std::isfinite(beta_value)) continue;
      double q_value = yy[index] - 2.0 * beta_value * numerator[index] +
                       beta_value * beta_value * denominator[index];
      if (q_value < 0.0) q_value = 0.0;
      const double sigma = std::sqrt(q_value / (count - 1.0));
      const double base_se = std::sqrt(1.0 / denominator[index]);
      if (!std::isfinite(sigma) || !std::isfinite(base_se)) continue;
      result_beta(exposure, outcome) = beta_value;
      result_sigma(exposure, outcome) = sigma;
      result_se(exposure, outcome) = base_se * std::max(1.0, sigma);
      result_q(exposure, outcome) = q_value;
    }
  }
  Rcpp::List output = Rcpp::List::create(
    Rcpp::_["nsnp"] = result_nsnp,
    Rcpp::_["beta"] = result_beta,
    Rcpp::_["se"] = result_se,
    Rcpp::_["Q"] = result_q,
    Rcpp::_["sigma"] = result_sigma
  );
  output.attr("class") = "fastmr_masked_ivw_compact";
  return output;
}

// Fused sparse exposure-by-dense outcome IVW.  The exposure instruments are
// supplied as a CSR matrix: row_ptr has E + 1 entries and col_index contains
// zero-based SNP column indices.  Outcome matrices are B x U in ordinary R
// layout (outcome rows, SNP columns).  An optional B x N mask addresses the
// concatenated N CSR entries: outcome o can independently retain entry k,
// while each exposure reads only its row_ptr-delimited entries.  This avoids
// materialising the very sparse E x U exposure panel and performs the four
// IVW accumulations in one pass over the non-zero exposure instruments.
Rcpp::List compute_sparse_ivw_grid(
    Rcpp::IntegerVector row_ptr, Rcpp::IntegerVector col_index,
    Rcpp::NumericVector exposure_beta, Rcpp::NumericMatrix outcome_beta,
    Rcpp::NumericMatrix outcome_se, Rcpp::LogicalMatrix outcome_present,
    int threads,
    Rcpp::Nullable<Rcpp::LogicalMatrix> pair_snp_keep = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> steiger_exposure_rsq = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> steiger_outcome_rsq = R_NilValue,
    Rcpp::Nullable<Rcpp::IntegerVector> drop_outcome = R_NilValue,
    Rcpp::Nullable<Rcpp::IntegerVector> drop_entry = R_NilValue) {
  const int exposure_count = row_ptr.size() - 1;
  const int outcome_count = outcome_beta.nrow();
  const int snp_count = outcome_beta.ncol();
  if (exposure_count <= 0 || outcome_count <= 0 || snp_count <= 0) {
    Rcpp::stop("sparse IVW inputs must have positive dimensions");
  }
  if (outcome_se.nrow() != outcome_count || outcome_se.ncol() != snp_count ||
      outcome_present.nrow() != outcome_count || outcome_present.ncol() != snp_count) {
    Rcpp::stop("sparse IVW outcome matrices must have matching dimensions");
  }
  if (row_ptr[0] != 0 || row_ptr[exposure_count] != col_index.size() ||
      exposure_beta.size() != col_index.size()) {
    Rcpp::stop("invalid sparse IVW CSR offsets or values");
  }
  const bool has_pair_snp_keep = pair_snp_keep.isNotNull();
  Rcpp::LogicalMatrix pair_snp_keep_matrix(0, 0);
  if (has_pair_snp_keep) {
    pair_snp_keep_matrix = Rcpp::as<Rcpp::LogicalMatrix>(pair_snp_keep);
    if (pair_snp_keep_matrix.nrow() != outcome_count ||
        pair_snp_keep_matrix.ncol() != col_index.size()) {
      Rcpp::stop("pair_snp_keep must have one row per outcome and one column per CSR entry");
    }
    for (R_xlen_t index = 0; index < pair_snp_keep_matrix.size(); ++index) {
      if (pair_snp_keep_matrix[index] == NA_LOGICAL) {
        Rcpp::stop("pair_snp_keep must not contain NA");
      }
    }
  }
  // Steiger masks are evaluated on the fly (no O x N matrix): keep iff
  // rsq_exposure[entry] > rsq_outcome[outcome, snp]; NA/NaN drops.
  const bool has_steiger = steiger_exposure_rsq.isNotNull() || steiger_outcome_rsq.isNotNull();
  Rcpp::NumericVector steiger_exp(0);
  Rcpp::NumericMatrix steiger_out(0, 0);
  if (has_steiger) {
    if (steiger_exposure_rsq.isNull() || steiger_outcome_rsq.isNull()) {
      Rcpp::stop("steiger_exposure_rsq and steiger_outcome_rsq must be supplied together");
    }
    steiger_exp = Rcpp::as<Rcpp::NumericVector>(steiger_exposure_rsq);
    steiger_out = Rcpp::as<Rcpp::NumericMatrix>(steiger_outcome_rsq);
    if (steiger_exp.size() != col_index.size()) {
      Rcpp::stop("steiger_exposure_rsq must have one value per CSR entry");
    }
    if (steiger_out.nrow() != outcome_count || steiger_out.ncol() != snp_count) {
      Rcpp::stop("steiger_outcome_rsq must have the same dimensions as outcome_beta");
    }
  }
  // Sparse drop list, sorted by (outcome, entry) with per-outcome offsets.
  const bool has_drop = drop_outcome.isNotNull() || drop_entry.isNotNull();
  std::vector<int> drop_offset;
  std::vector<int> drop_sorted;
  if (has_drop) {
    if (drop_outcome.isNull() || drop_entry.isNull()) {
      Rcpp::stop("pair_snp_drop must contain both outcome and entry");
    }
    Rcpp::IntegerVector d_out = Rcpp::as<Rcpp::IntegerVector>(drop_outcome);
    Rcpp::IntegerVector d_ent = Rcpp::as<Rcpp::IntegerVector>(drop_entry);
    if (d_out.size() != d_ent.size()) {
      Rcpp::stop("pair_snp_drop outcome and entry must have equal length");
    }
    drop_offset.assign(static_cast<std::size_t>(outcome_count) + 1, 0);
    for (R_xlen_t i = 0; i < d_out.size(); ++i) {
      if (d_out[i] == NA_INTEGER || d_ent[i] == NA_INTEGER ||
          d_out[i] < 1 || d_out[i] > outcome_count ||
          d_ent[i] < 1 || d_ent[i] > col_index.size()) {
        Rcpp::stop("pair_snp_drop indices must be 1-based and within the outcome and CSR entry ranges");
      }
      ++drop_offset[static_cast<std::size_t>(d_out[i])];
    }
    for (int o = 0; o < outcome_count; ++o) drop_offset[o + 1] += drop_offset[o];
    drop_sorted.resize(static_cast<std::size_t>(d_out.size()));
    std::vector<int> cursor(drop_offset.begin(), drop_offset.end() - 1);
    for (R_xlen_t i = 0; i < d_out.size(); ++i) {
      drop_sorted[cursor[d_out[i] - 1]++] = d_ent[i] - 1;
    }
    for (int o = 0; o < outcome_count; ++o) {
      std::sort(drop_sorted.begin() + drop_offset[o], drop_sorted.begin() + drop_offset[o + 1]);
    }
  }
  const bool report_prefilter = has_steiger || has_drop;
  for (int exposure = 0; exposure < exposure_count; ++exposure) {
    if (row_ptr[exposure] < 0 || row_ptr[exposure + 1] < row_ptr[exposure]) {
      Rcpp::stop("CSR row_ptr must be non-decreasing and non-negative");
    }
    std::unordered_set<int> seen;
    seen.reserve(static_cast<std::size_t>(row_ptr[exposure + 1] - row_ptr[exposure]));
    for (int index = row_ptr[exposure]; index < row_ptr[exposure + 1]; ++index) {
      if (!seen.insert(col_index[index]).second) {
        Rcpp::stop("sparse IVW CSR rows must not contain duplicate SNP indices");
      }
    }
  }
  for (int index = 0; index < col_index.size(); ++index) {
    if (col_index[index] < 0 || col_index[index] >= snp_count ||
        !std::isfinite(exposure_beta[index])) {
      Rcpp::stop("sparse IVW exposure entries must have valid SNP indices and finite betas");
    }
  }
  for (int outcome = 0; outcome < outcome_count; ++outcome) {
    for (int snp = 0; snp < snp_count; ++snp) {
      if (!outcome_present(outcome, snp)) continue;
      const double beta = outcome_beta(outcome, snp);
      const double se = outcome_se(outcome, snp);
      if (!std::isfinite(beta) || !std::isfinite(se) || se <= 0.0) {
        Rcpp::stop("present outcome values must be finite with positive standard errors");
      }
    }
  }

  Rcpp::NumericMatrix result_beta(exposure_count, outcome_count);
  Rcpp::NumericMatrix result_se(exposure_count, outcome_count);
  Rcpp::NumericMatrix result_q(exposure_count, outcome_count);
  Rcpp::NumericMatrix result_sigma(exposure_count, outcome_count);
  Rcpp::NumericMatrix result_nsnp(exposure_count, outcome_count);
  std::fill(result_beta.begin(), result_beta.end(), NA_VALUE);
  std::fill(result_se.begin(), result_se.end(), NA_VALUE);
  std::fill(result_q.begin(), result_q.end(), NA_VALUE);
  std::fill(result_sigma.begin(), result_sigma.end(), NA_VALUE);
  std::fill(result_nsnp.begin(), result_nsnp.end(), 0.0);
  Rcpp::NumericMatrix result_prefilter(report_prefilter ? exposure_count : 0,
                                       report_prefilter ? outcome_count : 0);

#ifdef _OPENMP
  const int thread_count = std::max(1, std::min(threads, exposure_count));
#pragma omp parallel for schedule(static) num_threads(thread_count)
#endif
  for (int exposure = 0; exposure < exposure_count; ++exposure) {
    const int first = row_ptr[exposure];
    const int last = row_ptr[exposure + 1];
    for (int outcome = 0; outcome < outcome_count; ++outcome) {
      double numerator = 0.0;
      double denominator = 0.0;
      double count = 0.0;
      double prefilter = 0.0;
      double single_x = NA_VALUE;
      double single_y = NA_VALUE;
      double single_se = NA_VALUE;
      const int* drop_it = nullptr;
      const int* drop_end = nullptr;
      if (has_drop) {
        const int* base = drop_sorted.data();
        drop_it = std::lower_bound(base + drop_offset[outcome],
                                   base + drop_offset[outcome + 1], first);
        drop_end = base + drop_offset[outcome + 1];
      }
      const int* drop_begin = drop_it;
      // Same filters for both passes; `cursor` walks the sorted drop list.
      auto kept = [&](int index, int snp, const int*& cursor) -> bool {
        if (has_pair_snp_keep && !pair_snp_keep_matrix(outcome, index)) return false;
        if (has_steiger &&
            !(steiger_exp[index] > steiger_out(outcome, snp))) return false;
        if (has_drop) {
          while (cursor != drop_end && *cursor < index) ++cursor;
          if (cursor != drop_end && *cursor == index) return false;
        }
        return true;
      };
      for (int index = first; index < last; ++index) {
        const int snp = col_index[index];
        if (!outcome_present(outcome, snp)) continue;
        prefilter += 1.0;
        if (!kept(index, snp, drop_it)) continue;
        const double x = exposure_beta[index];
        const double y = outcome_beta(outcome, snp);
        const double se = outcome_se(outcome, snp);
        const double weight = 1.0 / (se * se);
        // Same association order as compute_ivw() so results are bit-identical.
        denominator += weight * x * x;
        numerator += weight * x * y;
        count += 1.0;
        single_x = x;
        single_y = y;
        single_se = se;
      }
      const std::size_t result_index = static_cast<std::size_t>(exposure) +
                                       static_cast<std::size_t>(exposure_count) * outcome;
      result_nsnp[ result_index ] = count;
      if (report_prefilter) result_prefilter[ result_index ] = prefilter;
      if (count == 1.0) {
        // One instrument: the Wald ratio (see single_snp_wald()); Q, sigma NA.
        double beta_value = NA_VALUE;
        double se_value = NA_VALUE;
        if (single_snp_wald(single_x, single_y, single_se, beta_value, se_value)) {
          result_beta[ result_index ] = beta_value;
          result_se[ result_index ] = se_value;
        }
        continue;
      }
      if (count < 2.0 || !std::isfinite(denominator) || denominator <= 0.0) continue;
      const double beta_value = numerator / denominator;
      if (!std::isfinite(beta_value)) continue;
      // Two-pass residual sum (as compute_ivw); the normal-equation form
      // yy - 2*b*num + b^2*den cancels catastrophically for good fits.
      double q_value = 0.0;
      const int* cursor = drop_begin;
      for (int index = first; index < last; ++index) {
        const int snp = col_index[index];
        if (!outcome_present(outcome, snp)) continue;
        if (!kept(index, snp, cursor)) continue;
        const double se = outcome_se(outcome, snp);
        const double weight = 1.0 / (se * se);
        const double residual = outcome_beta(outcome, snp) - beta_value * exposure_beta[index];
        q_value += weight * residual * residual;
      }
      const double sigma_value = std::sqrt(q_value / (count - 1.0));
      const double base_se = std::sqrt(1.0 / denominator);
      if (!std::isfinite(sigma_value) || !std::isfinite(base_se)) continue;
      const double residual_se = base_se * sigma_value;
      result_beta[ result_index ] = beta_value;
      result_sigma[ result_index ] = sigma_value;
      result_se[ result_index ] = ivw_mre_se(base_se, residual_se, sigma_value);
      result_q[ result_index ] = q_value;
    }
  }
  Rcpp::List output = Rcpp::List::create(
    Rcpp::_["nsnp"] = result_nsnp,
    Rcpp::_["beta"] = result_beta,
    Rcpp::_["se"] = result_se,
    Rcpp::_["Q"] = result_q,
    Rcpp::_["sigma"] = result_sigma
  );
  if (report_prefilter) output["nsnp_prefilter"] = result_prefilter;
  output.attr("class") = "fastmr_sparse_ivw_compact";
  return output;
}

GridData copy_grid(Rcpp::NumericMatrix exp_beta, Rcpp::NumericMatrix out_beta,
                  Rcpp::NumericMatrix exp_se, Rcpp::NumericMatrix out_se) {
  validate_grid_shapes(exp_beta, out_beta, exp_se, out_se);
  GridData grid;
  grid.exposure_count = exp_beta.nrow();
  grid.outcome_count = out_beta.nrow();
  grid.snp_count = exp_beta.ncol();
  const std::size_t exp_size = static_cast<std::size_t>(grid.exposure_count) * grid.snp_count;
  const std::size_t out_size = static_cast<std::size_t>(grid.outcome_count) * grid.snp_count;
  grid.exp_beta.resize(exp_size); grid.exp_se.resize(exp_size);
  grid.out_beta.resize(out_size); grid.out_se.resize(out_size);
  for (int i = 0; i < grid.exposure_count; ++i) for (int j = 0; j < grid.snp_count; ++j) {
    grid.exp_beta[static_cast<std::size_t>(i) * grid.snp_count + j] = exp_beta(i, j);
    grid.exp_se[static_cast<std::size_t>(i) * grid.snp_count + j] = exp_se(i, j);
  }
  for (int i = 0; i < grid.outcome_count; ++i) for (int j = 0; j < grid.snp_count; ++j) {
    grid.out_beta[static_cast<std::size_t>(i) * grid.snp_count + j] = out_beta(i, j);
    grid.out_se[static_cast<std::size_t>(i) * grid.snp_count + j] = out_se(i, j);
  }
  return grid;
}

void fill_grid_bootstrap_layout(GridData& grid, int nboot, SEXP seed,
                                 bool needs_median, bool needs_mode,
                                 bool needs_egger, bool needs_penalised) {
  if (nboot <= 0 || (!needs_median && !needs_mode && !needs_egger && !needs_penalised)) return;
  (void) seed;
  const std::size_t n = static_cast<std::size_t>(grid.snp_count);
  const std::size_t block = static_cast<std::size_t>(nboot) * n;

  if (needs_median || needs_egger || needs_penalised) {
    if (needs_median) grid.exp_inverse.resize(static_cast<std::size_t>(grid.exposure_count) * block);
    if (needs_egger) grid.exp_draws.resize(static_cast<std::size_t>(grid.exposure_count) * block);
    grid.out_draws.resize(static_cast<std::size_t>(grid.outcome_count) * block);
    std::vector<double> ze(static_cast<std::size_t>(nboot) * n);
    std::vector<double> zo(static_cast<std::size_t>(nboot) * n);
    for (std::size_t snp = 0; snp < n; ++snp) {
      for (int draw = 0; draw < nboot; ++draw) {
        ze[static_cast<std::size_t>(draw) * n + snp] = R::rnorm(0.0, 1.0);
      }
    }
    for (std::size_t snp = 0; snp < n; ++snp) {
      for (int draw = 0; draw < nboot; ++draw) {
        zo[static_cast<std::size_t>(draw) * n + snp] = R::rnorm(0.0, 1.0);
      }
    }
    for (int draw = 0; draw < nboot; ++draw) {
      for (int exposure = 0; exposure < grid.exposure_count; ++exposure) {
        const std::size_t source = static_cast<std::size_t>(exposure) * n;
        const std::size_t target = static_cast<std::size_t>(exposure) * block +
                                   static_cast<std::size_t>(draw) * n;
        for (std::size_t snp = 0; snp < n; ++snp) {
          const double value = grid.exp_beta[source + snp] +
            grid.exp_se[source + snp] * ze[static_cast<std::size_t>(draw) * n + snp];
          if (needs_egger) grid.exp_draws[target + snp] = value;
          if (needs_median) grid.exp_inverse[target + snp] = value == 0.0 ? NA_VALUE : 1.0 / value;
        }
      }
      for (int outcome = 0; outcome < grid.outcome_count; ++outcome) {
        const std::size_t source = static_cast<std::size_t>(outcome) * n;
        const std::size_t target = static_cast<std::size_t>(outcome) * block +
                                   static_cast<std::size_t>(draw) * n;
        for (std::size_t snp = 0; snp < n; ++snp) {
          grid.out_draws[target + snp] = grid.out_beta[source + snp] +
            grid.out_se[source + snp] * zo[static_cast<std::size_t>(draw) * n + snp];
        }
      }
    }
  }

  if (needs_penalised) {
    grid.exp_inverse_penalised.resize(static_cast<std::size_t>(grid.exposure_count) * block);
    grid.out_draws_penalised.resize(static_cast<std::size_t>(grid.outcome_count) * block);
    std::vector<double> ze(static_cast<std::size_t>(nboot) * n);
    std::vector<double> zo(static_cast<std::size_t>(nboot) * n);
    for (std::size_t snp = 0; snp < n; ++snp) {
      for (int draw = 0; draw < nboot; ++draw) {
        ze[static_cast<std::size_t>(draw) * n + snp] = R::rnorm(0.0, 1.0);
      }
    }
    for (std::size_t snp = 0; snp < n; ++snp) {
      for (int draw = 0; draw < nboot; ++draw) {
        zo[static_cast<std::size_t>(draw) * n + snp] = R::rnorm(0.0, 1.0);
      }
    }
    for (int draw = 0; draw < nboot; ++draw) {
      for (int exposure = 0; exposure < grid.exposure_count; ++exposure) {
        const std::size_t source = static_cast<std::size_t>(exposure) * n;
        const std::size_t target = static_cast<std::size_t>(exposure) * block +
                                   static_cast<std::size_t>(draw) * n;
        for (std::size_t snp = 0; snp < n; ++snp) {
          const double value = grid.exp_beta[source + snp] +
            grid.exp_se[source + snp] * ze[static_cast<std::size_t>(draw) * n + snp];
          grid.exp_inverse_penalised[target + snp] = value == 0.0 ? NA_VALUE : 1.0 / value;
        }
      }
      for (int outcome = 0; outcome < grid.outcome_count; ++outcome) {
        const std::size_t source = static_cast<std::size_t>(outcome) * n;
        const std::size_t target = static_cast<std::size_t>(outcome) * block +
                                   static_cast<std::size_t>(draw) * n;
        for (std::size_t snp = 0; snp < n; ++snp) {
          grid.out_draws_penalised[target + snp] = grid.out_beta[source + snp] +
            grid.out_se[source + snp] * zo[static_cast<std::size_t>(draw) * n + snp];
        }
      }
    }
  }

  if (needs_mode) {
    grid.mode_z.resize(block);
    for (std::size_t snp = 0; snp < n; ++snp) {
      for (int draw = 0; draw < nboot; ++draw) {
        grid.mode_z[static_cast<std::size_t>(draw) * n + snp] = R::rnorm(0.0, 1.0);
      }
    }
  }
}

Prepared pair_from_grid(const GridData& grid, int exposure, int outcome,
                        int nboot) {
  const std::size_t n = static_cast<std::size_t>(grid.snp_count);
  const std::size_t exp_offset = static_cast<std::size_t>(exposure) * n;
  const std::size_t out_offset = static_cast<std::size_t>(outcome) * n;
  Prepared p;
  p.x.resize(n); p.y.resize(n); p.sx.resize(n); p.sy.resize(n);
  for (std::size_t snp = 0; snp < n; ++snp) {
    p.x[snp] = grid.exp_beta[exp_offset + snp];
    p.sx[snp] = grid.exp_se[exp_offset + snp];
    p.y[snp] = grid.out_beta[out_offset + snp];
    p.sy[snp] = grid.out_se[out_offset + snp];
  }
  prepare_ratios(p);
  if (nboot > 0 && p.ratio.size() >= 3 && (!grid.exp_inverse.empty() || !grid.exp_draws.empty())) {
    const std::size_t block = static_cast<std::size_t>(nboot) * n;
    std::size_t ratio_count = p.ratio.size();
    if (!grid.exp_inverse.empty()) p.bootstrap.resize(static_cast<std::size_t>(nboot) * ratio_count);
    if (!grid.exp_draws.empty()) {
      p.egger_x_bootstrap.resize(static_cast<std::size_t>(nboot) * n);
      p.egger_y_bootstrap.resize(static_cast<std::size_t>(nboot) * n);
    }
    for (int draw = 0; draw < nboot; ++draw) {
      const std::size_t raw = static_cast<std::size_t>(draw) * n;
      const std::size_t exp_raw = static_cast<std::size_t>(exposure) * block + raw;
      const std::size_t out_raw = static_cast<std::size_t>(outcome) * block + raw;
      std::size_t ratio_index = 0;
      for (std::size_t snp = 0; snp < n; ++snp) {
        if (!grid.exp_draws.empty()) {
          p.egger_x_bootstrap[static_cast<std::size_t>(draw) * n + snp] = grid.exp_draws[exp_raw + snp];
          p.egger_y_bootstrap[static_cast<std::size_t>(draw) * n + snp] = grid.out_draws[out_raw + snp];
        }
        if (p.x[snp] == 0.0) continue;
        if (!grid.exp_inverse.empty()) {
          p.bootstrap[static_cast<std::size_t>(draw) * ratio_count + ratio_index] =
            grid.out_draws[out_raw + snp] * grid.exp_inverse[exp_raw + snp];
        }
        ++ratio_index;
      }
    }
  }
  if (nboot > 0 && p.ratio.size() >= 3 && !grid.exp_inverse_penalised.empty()) {
    const std::size_t block = static_cast<std::size_t>(nboot) * n;
    const std::size_t ratio_count = p.ratio.size();
    p.penalised_bootstrap.resize(static_cast<std::size_t>(nboot) * ratio_count);
    for (int draw = 0; draw < nboot; ++draw) {
      const std::size_t raw = static_cast<std::size_t>(draw) * n;
      const std::size_t exp_raw = static_cast<std::size_t>(exposure) * block + raw;
      const std::size_t out_raw = static_cast<std::size_t>(outcome) * block + raw;
      std::size_t ratio_index = 0;
      for (std::size_t snp = 0; snp < n; ++snp) {
        if (p.x[snp] == 0.0) continue;
        p.penalised_bootstrap[static_cast<std::size_t>(draw) * ratio_count + ratio_index] =
          grid.out_draws_penalised[out_raw + snp] * grid.exp_inverse_penalised[exp_raw + snp];
        ++ratio_index;
      }
    }
  }
  if (nboot > 0 && p.ratio.size() >= 3 && !grid.mode_z.empty()) {
    const std::size_t ratio_count = p.ratio.size();
    p.mode_bootstrap.resize(static_cast<std::size_t>(nboot) * ratio_count);
    for (int draw = 0; draw < nboot; ++draw) {
      const std::size_t raw = static_cast<std::size_t>(draw) * n;
      std::size_t ratio_index = 0;
      for (std::size_t snp = 0; snp < n; ++snp) {
        if (p.x[snp] == 0.0) continue;
        p.mode_bootstrap[static_cast<std::size_t>(draw) * ratio_count + ratio_index] =
          p.ratio[ratio_index] + p.ratio_se[ratio_index] * grid.mode_z[raw + snp];
        ++ratio_index;
      }
    }
  }
  return p;
}

Rcpp::List results_to_list(const std::vector<Result>& results) {
  Rcpp::List output(results.size());
  for (std::size_t i = 0; i < results.size(); ++i) output[i] = result_to_list(results[i]);
  return output;
}

Rcpp::List results_to_compact_grid(const std::vector<Result>& results,
                                   const std::vector<std::string>& methods) {
  const int method_count = static_cast<int>(methods.size());
  const std::size_t pair_count = method_count == 0 ? 0 : results.size() / static_cast<std::size_t>(method_count);
  Rcpp::NumericMatrix nsnp(method_count, pair_count);
  Rcpp::NumericMatrix beta(method_count, pair_count);
  Rcpp::NumericMatrix se(method_count, pair_count);
  Rcpp::NumericMatrix pval(method_count, pair_count);
  Rcpp::NumericMatrix ratio_se_mean(method_count, pair_count);
  Rcpp::NumericMatrix bootstrap(method_count, pair_count);
  Rcpp::NumericMatrix phi(method_count, pair_count);
  Rcpp::NumericMatrix q(method_count, pair_count);
  Rcpp::NumericMatrix q_df(method_count, pair_count);
  Rcpp::NumericMatrix q_pval(method_count, pair_count);
  Rcpp::NumericMatrix sigma(method_count, pair_count);
  Rcpp::NumericMatrix intercept(method_count, pair_count);
  Rcpp::NumericMatrix intercept_se(method_count, pair_count);
  Rcpp::NumericMatrix intercept_pval(method_count, pair_count);
  Rcpp::NumericMatrix flipped(method_count, pair_count);
  Rcpp::NumericMatrix se_exposure_mean(method_count, pair_count);
  const std::vector<Rcpp::NumericMatrix*> fields = {
    &nsnp, &beta, &se, &pval, &ratio_se_mean, &bootstrap, &phi, &q, &q_df,
    &q_pval, &sigma, &intercept, &intercept_se, &intercept_pval, &flipped,
    &se_exposure_mean
  };
  for (Rcpp::NumericMatrix* field : fields) std::fill(field->begin(), field->end(), NA_VALUE);
  for (std::size_t pair = 0; pair < pair_count; ++pair) {
    for (int method_index = 0; method_index < method_count; ++method_index) {
      const Result& result = results[pair * static_cast<std::size_t>(method_count) +
                                      static_cast<std::size_t>(method_index)];
      nsnp(method_index, pair) = result.n;
      beta(method_index, pair) = finite_or_na(result.beta);
      se(method_index, pair) = finite_or_na(result.se);
      pval(method_index, pair) = finite_or_na(result.pval);
      if (result.ratio_se_mean) ratio_se_mean(method_index, pair) = finite_or_na(result.ratio_se_mean_value);
      if (result.bootstrap) bootstrap(method_index, pair) = result.bootstrap_value;
      if (result.phi) phi(method_index, pair) = finite_or_na(result.phi_value);
      if (result.q) {
        q(method_index, pair) = finite_or_na(result.q_value);
        q_df(method_index, pair) = result.q_df;
        q_pval(method_index, pair) = finite_or_na(result.q_pval);
      }
      if (result.sigma) sigma(method_index, pair) = finite_or_na(result.sigma_value);
      if (result.intercept) {
        intercept(method_index, pair) = finite_or_na(result.intercept_value);
        intercept_se(method_index, pair) = finite_or_na(result.intercept_se);
        intercept_pval(method_index, pair) = finite_or_na(result.intercept_pval);
        flipped(method_index, pair) = result.flipped;
        se_exposure_mean(method_index, pair) = finite_or_na(result.se_exposure_mean);
      }
    }
  }
  Rcpp::CharacterVector method_codes(method_count);
  for (int i = 0; i < method_count; ++i) method_codes[i] = methods[static_cast<std::size_t>(i)];
  Rcpp::List output = Rcpp::List::create(
    Rcpp::_["methods"] = method_codes,
    Rcpp::_["nsnp"] = nsnp,
    Rcpp::_["beta"] = beta,
    Rcpp::_["se"] = se,
    Rcpp::_["pval"] = pval,
    Rcpp::_["ratio_se_mean"] = ratio_se_mean,
    Rcpp::_["bootstrap"] = bootstrap,
    Rcpp::_["phi"] = phi,
    Rcpp::_["Q"] = q,
    Rcpp::_["Q_df"] = q_df,
    Rcpp::_["Q_pval"] = q_pval,
    Rcpp::_["sigma"] = sigma,
    Rcpp::_["intercept"] = intercept,
    Rcpp::_["intercept_se"] = intercept_se,
    Rcpp::_["intercept_pval"] = intercept_pval,
    Rcpp::_["flipped"] = flipped,
    Rcpp::_["se_exposure_mean"] = se_exposure_mean
  );
  output.attr("class") = "fastmr_grid_compact";
  return output;
}

int bounded_threads(int requested, std::size_t jobs) {
  const int bounded = std::max(1, std::min<int>(requested, static_cast<int>(std::max<std::size_t>(1, jobs))));
  return bounded;
}

// Minimum work per worker before an extra thread pays for its spawn/join and
// serial overheads. Work is counted in SNP rows (RNG-free group fits, grid
// pairs x SNPs) or, for bootstrap batches, in median/Egger stream normals
// (~20 ns each) plus kModeDrawCost per mode normal (bootstrap_draw_work()).
// The closed-form methods (ivw, egger, ...) cost ~10-100 ns per row and are
// dominated by serial p-value work, so they need ~100k rows per worker.
// Calibrated (scripts in the mode-density PR) as the smallest work at which
// 2/4/8 threads stopped being slower than 1, on a Mac mini (std::thread
// fallback, spawns threads on every call) and on 8 Slurm CPUs (4 cores x 2
// hyperthreads, OpenMP): stream-only bootstraps (weighted median, Egger
// bootstrap) broke even at 16-32k normals per worker on the Mac and 32-64k on
// the cluster; mode bootstraps at < 1000 units; RNG-free mode fits at
// 200-500 rows.
constexpr double kMinRowsPerWorkerHeavy = 512.0;
constexpr double kMinRowsPerWorkerCheap = 100000.0;
constexpr double kMinDrawsPerWorker = 65536.0;

// Test hook: scales the minimum work per worker (0 forces the parallel paths
// on tiny inputs so the thread-equivalence tests still exercise them).
std::atomic<double> min_work_scale(1.0);

// Only use extra threads when each would get at least `min_work_per_worker`
// units of work. Never changes results: they are identical for any count.
int worthwhile_threads(int requested, std::size_t jobs, double work,
                       double min_work_per_worker) {
  const int bounded = bounded_threads(requested, jobs);
  const double min_work = min_work_per_worker * min_work_scale.load(std::memory_order_relaxed);
  if (!(min_work > 0.0)) return bounded;
  const double by_work = std::floor(work / min_work);
  if (!(by_work >= 1.0)) return 1;
  return std::min(bounded, static_cast<int>(std::min(by_work, 1.0e6)));
}

double min_rows_per_worker(const BootstrapNeeds& needs) {
  return needs.median || needs.mode || needs.penalised || needs.egger
    ? kMinRowsPerWorkerHeavy : kMinRowsPerWorkerCheap;
}

} // namespace

// Internal test hook: set the multiplier on the minimum work per worker and
// return the previous value (1 by default; 0 uses all requested threads).
// [[Rcpp::export]]
double fastmr_set_work_scale_native(double scale) {
  if (!std::isfinite(scale) || scale < 0.0) Rcpp::stop("scale must be non-negative and finite");
  return min_work_scale.exchange(scale);
}

// Internal test hooks for the mode-density paths. Set the largest ratio count
// that tries the direct path (0 = always FFT, Inf = always direct) and return
// the previous value; results are identical for every setting.
// [[Rcpp::export]]
double fastmr_set_mode_direct_max_native(double ratios) {
  if (std::isnan(ratios) || ratios < 0.0) Rcpp::stop("ratios must be non-negative");
  return mode_direct_max.exchange(ratios);
}

// Mode-density draws that tried the direct path, and those of them the guard
// sent to the FFT path, since the last reset.
// [[Rcpp::export]]
Rcpp::NumericVector fastmr_mode_path_counts_native(bool reset = false) {
  Rcpp::NumericVector out = Rcpp::NumericVector::create(
    Rcpp::_["direct"] = static_cast<double>(mode_direct_total.load()),
    Rcpp::_["guard"] = static_cast<double>(mode_guard_total.load()),
    Rcpp::_["hull"] = static_cast<double>(mode_hull_total.load()),
    Rcpp::_["hull_fallback"] = static_cast<double>(mode_hull_fallback_total.load()));
  if (reset) {
    mode_direct_total.store(0);
    mode_guard_total.store(0);
    mode_hull_total.store(0);
    mode_hull_fallback_total.store(0);
  }
  return out;
}

// Internal test hook for the hull path: set the largest ratio count that
// tries it first (0 disables it) and whether its kernel uses the recurrence;
// returns the previous settings. Results are identical for every setting.
// [[Rcpp::export]]
Rcpp::NumericVector fastmr_set_mode_hull_native(double ratios, bool recurrence = true) {
  if (std::isnan(ratios) || ratios < 0.0) Rcpp::stop("ratios must be non-negative");
  const double previous_max = mode_hull_max.exchange(ratios);
  const bool previous_recurrence = mode_hull_recurrence.exchange(recurrence);
  return Rcpp::NumericVector::create(Rcpp::_["ratios"] = previous_max,
                                     Rcpp::_["recurrence"] = previous_recurrence ? 1.0 : 0.0);
}

// [[Rcpp::export]]
Rcpp::List fastmr_run_native(Rcpp::NumericVector exposure_beta,
                             Rcpp::NumericVector outcome_beta,
                             Rcpp::NumericVector exposure_se,
                             Rcpp::NumericVector outcome_se,
                             Rcpp::CharacterVector methods,
                             int nboot = 1000,
                             SEXP seed = R_NilValue,
                             int threads = 1,
                             double phi = 1.0,
                             double penk = 20.0) {
  Rcpp::RNGScope scope;
  validate_controls(nboot, threads, phi);
  if (!std::isfinite(penk) || penk <= 0.0) Rcpp::stop("penk must be positive and finite");
  const std::vector<std::string> parsed_methods = parse_methods(methods);
  Prepared prepared = one_pair_from_vectors(exposure_beta, outcome_beta, exposure_se, outcome_se);
  return results_to_list(compute_pair(std::move(prepared), parsed_methods, nboot, seed, true, phi, penk));
}

// [[Rcpp::export]]
Rcpp::List fastmr_grid_native(Rcpp::NumericMatrix exposure_beta,
                              Rcpp::NumericMatrix outcome_beta,
                              Rcpp::NumericMatrix exposure_se,
                              Rcpp::NumericMatrix outcome_se,
                              Rcpp::CharacterVector methods,
                              int nboot = 1000,
                              SEXP seed = R_NilValue,
                              int threads = 1,
                              double phi = 1.0,
                              double penk = 20.0) {
  Rcpp::RNGScope scope;
  validate_controls(nboot, threads, phi);
  if (!std::isfinite(penk) || penk <= 0.0) Rcpp::stop("penk must be positive and finite");
  const std::vector<std::string> parsed_methods = parse_methods(methods);
  GridData grid = copy_grid(exposure_beta, outcome_beta, exposure_se, outcome_se);
  const std::size_t pair_count = static_cast<std::size_t>(grid.exposure_count) *
                                 static_cast<std::size_t>(grid.outcome_count);
  if (only_ivw_methods(parsed_methods)) {
    return compute_ivw_grid_compact(grid, parsed_methods);
  }
  std::vector<std::string> ivw_methods;
  std::vector<std::string> other_methods;
  ivw_methods.reserve(parsed_methods.size());
  other_methods.reserve(parsed_methods.size());
  for (std::size_t i = 0; i < parsed_methods.size(); ++i) {
    const std::string& method = parsed_methods[i];
    if (method == "ivw" || method == "ivw_fe" || method == "ivw_mre") {
      ivw_methods.push_back(method);
    } else {
      other_methods.push_back(method);
    }
  }
  // Compute every IVW flavour through the same BLAS batch used by an IVW-only
  // request. The remaining methods still run pairwise, but never redo the
  // dominant cross-product for the IVW rows.
  std::vector<Result> ivw_results;
  if (!ivw_methods.empty()) ivw_results = compute_ivw_grid_blas(grid, ivw_methods);
  bool needs_median = false;
  bool needs_penalised = false;
  bool needs_mode = false;
  bool needs_egger = false;
  for (const std::string& method : parsed_methods) {
    needs_median = needs_median || method == "simple_median" ||
                   method == "weighted_median";
    needs_penalised = needs_penalised || method == "penalised_weighted_median";
    needs_mode = needs_mode || method == "simple_mode" || method == "weighted_mode";
    needs_egger = needs_egger || method == "egger_bootstrap";
  }
  fill_grid_bootstrap_layout(grid, nboot, seed, needs_median, needs_mode, needs_egger, needs_penalised);
  std::vector<Result> results(pair_count * parsed_methods.size());
  // Resampling pairs cost what the same bootstrap costs in a batch fit
  // (bootstrap_draw_work() per pair); otherwise count rows.
  const bool resampling =
    (needs_median || needs_penalised || needs_mode || needs_egger) && nboot > 0;
  const std::size_t grid_snps = static_cast<std::size_t>(grid.snp_count);
  const double pair_rows = static_cast<double>(pair_count) * static_cast<double>(grid_snps);
  const int thread_count = resampling
    ? worthwhile_threads(threads, pair_count,
        static_cast<double>(pair_count) *
          bootstrap_draw_work(grid_snps, grid_snps, nboot, needs_median, needs_egger,
                              needs_penalised, needs_mode) + pair_rows,
        kMinDrawsPerWorker)
    : worthwhile_threads(threads, pair_count, pair_rows,
        min_rows_per_worker(BootstrapNeeds{needs_median, needs_penalised, needs_mode, needs_egger}));
  defer_r_math.store(true, std::memory_order_relaxed);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(thread_count)
  // Signed 64-bit index: pair_count = exposures x outcomes can exceed INT_MAX.
  for (std::ptrdiff_t index = 0; index < static_cast<std::ptrdiff_t>(pair_count); ++index) {
    const std::size_t pair = static_cast<std::size_t>(index);
    const int exposure = static_cast<int>(pair / static_cast<std::size_t>(grid.outcome_count));
    const int outcome = static_cast<int>(pair % static_cast<std::size_t>(grid.outcome_count));
    std::vector<Result> other_results = compute_pair(
      pair_from_grid(grid, exposure, outcome, nboot), other_methods, nboot,
      R_NilValue, false, phi, penk);
    std::vector<Result> full_results(parsed_methods.size());
    std::size_t other_index = 0;
    std::size_t ivw_index = 0;
    for (std::size_t method_index = 0; method_index < parsed_methods.size(); ++method_index) {
      if (parsed_methods[method_index] == "ivw" ||
          parsed_methods[method_index] == "ivw_fe" ||
          parsed_methods[method_index] == "ivw_mre") {
        full_results[method_index] = ivw_results[pair * ivw_methods.size() + ivw_index++];
      } else {
        full_results[method_index] = std::move(other_results[other_index++]);
      }
    }
    std::move(full_results.begin(), full_results.end(),
              results.begin() + static_cast<std::ptrdiff_t>(pair * parsed_methods.size()));
  }
#else
  if (thread_count == 1) {
    for (std::size_t index = 0; index < pair_count; ++index) {
      const int exposure = static_cast<int>(index / static_cast<std::size_t>(grid.outcome_count));
      const int outcome = static_cast<int>(index % static_cast<std::size_t>(grid.outcome_count));
      std::vector<Result> other_results = compute_pair(
        pair_from_grid(grid, exposure, outcome, nboot), other_methods, nboot,
        R_NilValue, false, phi, penk);
      std::vector<Result> full_results(parsed_methods.size());
      std::size_t other_index = 0;
      std::size_t ivw_index = 0;
      for (std::size_t method_index = 0; method_index < parsed_methods.size(); ++method_index) {
        if (parsed_methods[method_index] == "ivw" ||
            parsed_methods[method_index] == "ivw_fe" ||
            parsed_methods[method_index] == "ivw_mre") {
          full_results[method_index] = ivw_results[index * ivw_methods.size() + ivw_index++];
        } else {
          full_results[method_index] = std::move(other_results[other_index++]);
        }
      }
      std::move(full_results.begin(), full_results.end(),
                results.begin() + static_cast<std::ptrdiff_t>(index * parsed_methods.size()));
    }
  } else {
    std::atomic<std::size_t> next(0);
    std::vector<std::thread> pool;
    pool.reserve(static_cast<std::size_t>(thread_count));
    for (int worker = 0; worker < thread_count; ++worker) {
      pool.emplace_back([&]() {
        while (true) {
          const std::size_t index = next.fetch_add(1, std::memory_order_relaxed);
          if (index >= pair_count) break;
          const int exposure = static_cast<int>(index / static_cast<std::size_t>(grid.outcome_count));
          const int outcome = static_cast<int>(index % static_cast<std::size_t>(grid.outcome_count));
          std::vector<Result> other_results = compute_pair(
            pair_from_grid(grid, exposure, outcome, nboot), other_methods, nboot,
            R_NilValue, false, phi, penk);
          std::vector<Result> full_results(parsed_methods.size());
          std::size_t other_index = 0;
          std::size_t ivw_index = 0;
          for (std::size_t method_index = 0; method_index < parsed_methods.size(); ++method_index) {
            if (parsed_methods[method_index] == "ivw" ||
                parsed_methods[method_index] == "ivw_fe" ||
                parsed_methods[method_index] == "ivw_mre") {
              full_results[method_index] = ivw_results[index * ivw_methods.size() + ivw_index++];
            } else {
              full_results[method_index] = std::move(other_results[other_index++]);
            }
          }
          std::move(full_results.begin(), full_results.end(),
                    results.begin() + static_cast<std::ptrdiff_t>(index * parsed_methods.size()));
        }
      });
    }
    for (std::thread& worker : pool) worker.join();
  }
#endif

  defer_r_math.store(false, std::memory_order_relaxed);
  for (Result& result : results) populate_result_pvalues(result);
  return results_to_compact_grid(results, parsed_methods);
}

// [[Rcpp::export]]
Rcpp::List fastmr_masked_ivw_native(
    Rcpp::NumericMatrix exposure_beta, Rcpp::NumericMatrix outcome_beta,
    Rcpp::NumericMatrix outcome_se, Rcpp::LogicalMatrix exposure_present,
    Rcpp::LogicalMatrix outcome_present, int threads = 1) {
  validate_controls(0, threads, 1.0);
  return compute_masked_ivw_grid(
    exposure_beta, outcome_beta, outcome_se, exposure_present, outcome_present);
}

// [[Rcpp::export]]
Rcpp::List fastmr_sparse_ivw_native(
    Rcpp::IntegerVector row_ptr, Rcpp::IntegerVector col_index,
    Rcpp::NumericVector exposure_beta, Rcpp::NumericMatrix outcome_beta,
    Rcpp::NumericMatrix outcome_se, Rcpp::LogicalMatrix outcome_present,
    int threads = 1,
    Rcpp::Nullable<Rcpp::LogicalMatrix> pair_snp_keep = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> steiger_exposure_rsq = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> steiger_outcome_rsq = R_NilValue,
    Rcpp::Nullable<Rcpp::IntegerVector> drop_outcome = R_NilValue,
    Rcpp::Nullable<Rcpp::IntegerVector> drop_entry = R_NilValue) {
  validate_controls(0, threads, 1.0);
  return compute_sparse_ivw_grid(row_ptr, col_index, exposure_beta,
                                 outcome_beta, outcome_se, outcome_present,
                                 threads, pair_snp_keep, steiger_exposure_rsq,
                                 steiger_outcome_rsq, drop_outcome, drop_entry);
}

namespace {

struct GroupInputs {
  std::vector<R_xlen_t> offsets;
  const double* x;
  const double* y;
  const double* sx;
  const double* sy;
  R_xlen_t groups;
};

GroupInputs check_group_inputs(Rcpp::IntegerVector offsets,
                               Rcpp::NumericVector exposure_beta,
                               Rcpp::NumericVector outcome_beta,
                               Rcpp::NumericVector exposure_se,
                               Rcpp::NumericVector outcome_se) {
  if (offsets.size() < 1 || offsets[0] != 0) Rcpp::stop("offsets must start at 0");
  const R_xlen_t groups = offsets.size() - 1;
  const R_xlen_t total_rows = exposure_beta.size();
  if (outcome_beta.size() != total_rows || exposure_se.size() != total_rows ||
      outcome_se.size() != total_rows || offsets[groups] != total_rows) {
    Rcpp::stop("MR vectors must have equal lengths matching offsets");
  }
  for (R_xlen_t g = 0; g < groups; ++g) {
    if (offsets[g + 1] < offsets[g]) Rcpp::stop("offsets must be non-decreasing");
  }
  GroupInputs in;
  in.offsets.assign(offsets.begin(), offsets.end());
  in.x = REAL(exposure_beta); in.y = REAL(outcome_beta);
  in.sx = REAL(exposure_se); in.sy = REAL(outcome_se);
  in.groups = groups;
  return in;
}

// Same row filter as one_pair_from_vectors() (NA is never finite), without
// touching R objects so worker threads may call it.
Prepared group_prepared(const GroupInputs& in, R_xlen_t g) {
  Prepared p;
  const R_xlen_t begin = in.offsets[static_cast<std::size_t>(g)];
  const R_xlen_t end = in.offsets[static_cast<std::size_t>(g) + 1];
  p.x.reserve(end - begin); p.y.reserve(end - begin);
  p.sx.reserve(end - begin); p.sy.reserve(end - begin);
  for (R_xlen_t i = begin; i < end; ++i) {
    const double x = in.x[i], y = in.y[i], sx = in.sx[i], sy = in.sy[i];
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(sx) ||
        !std::isfinite(sy) || sx <= 0.0 || sy <= 0.0) continue;
    p.x.push_back(x); p.y.push_back(y); p.sx.push_back(sx); p.sy.push_back(sy);
  }
  return p;
}

// Flat group-major, method-minor list matching rbind(fastmr_tidy_native()).
Rcpp::List group_results_to_flat(const std::vector<Result>& results) {
  const std::size_t total = results.size();
  Rcpp::CharacterVector method_out(total);
  Rcpp::NumericVector n_out(total, NA_REAL), beta(total, NA_REAL), se(total, NA_REAL),
    pval(total, NA_REAL), q(total, NA_REAL), q_df(total, NA_REAL), q_pval(total, NA_REAL),
    sigma(total, NA_REAL), intercept(total, NA_REAL), intercept_se(total, NA_REAL),
    intercept_pval(total, NA_REAL), ratio_se_mean(total, NA_REAL), boot(total, NA_REAL),
    phi_out(total, NA_REAL), flipped(total, NA_REAL), se_exposure_mean(total, NA_REAL);
  for (std::size_t k = 0; k < total; ++k) {
    const Result& r = results[k];
    method_out[k] = r.method;
    n_out[k] = r.n;
    beta[k] = finite_or_na(r.beta);
    se[k] = finite_or_na(r.se);
    pval[k] = finite_or_na(r.pval);
    if (r.ratio_se_mean) ratio_se_mean[k] = finite_or_na(r.ratio_se_mean_value);
    if (r.bootstrap) boot[k] = r.bootstrap_value;
    if (r.phi) phi_out[k] = finite_or_na(r.phi_value);
    if (r.q) {
      q[k] = finite_or_na(r.q_value);
      q_df[k] = r.q_df;
      q_pval[k] = finite_or_na(r.q_pval);
    }
    if (r.sigma) sigma[k] = finite_or_na(r.sigma_value);
    if (r.intercept) {
      intercept[k] = finite_or_na(r.intercept_value);
      intercept_se[k] = finite_or_na(r.intercept_se);
      intercept_pval[k] = finite_or_na(r.intercept_pval);
      flipped[k] = r.flipped;
      se_exposure_mean[k] = finite_or_na(r.se_exposure_mean);
    }
  }
  return Rcpp::List::create(
    Rcpp::_["method"] = method_out, Rcpp::_["n"] = n_out, Rcpp::_["beta"] = beta,
    Rcpp::_["se"] = se, Rcpp::_["pval"] = pval, Rcpp::_["Q"] = q, Rcpp::_["Q_df"] = q_df,
    Rcpp::_["Q_pval"] = q_pval, Rcpp::_["sigma"] = sigma, Rcpp::_["intercept"] = intercept,
    Rcpp::_["intercept_se"] = intercept_se, Rcpp::_["intercept_pval"] = intercept_pval,
    Rcpp::_["ratio_se_mean"] = ratio_se_mean, Rcpp::_["bootstrap"] = boot,
    Rcpp::_["phi"] = phi_out, Rcpp::_["flipped"] = flipped,
    Rcpp::_["se_exposure_mean"] = se_exposure_mean);
}

// Defers R-backed p-values for the guard's lifetime (restored on unwind).
struct DeferRMath {
  DeferRMath() { defer_r_math.store(true, std::memory_order_relaxed); }
  ~DeferRMath() { defer_r_math.store(false, std::memory_order_relaxed); }
};

// Run job(index) for index in [0, jobs) on `threads` workers. job must not
// touch the R API. Dynamic scheduling; results are written by index, so the
// output does not depend on the schedule.
template <typename Job>
void run_parallel(std::size_t jobs, int threads, Job job) {
  const int thread_count = bounded_threads(threads, jobs);
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

} // namespace

namespace {

// Run compute_pair() on `jobs` RNG-free inputs, job-major / method-minor.
// prepared(job) must build the job's Prepared without touching R. With one
// worker every p-value is computed inline, exactly as one fastmr_run_native()
// call per job; with several, R-backed p-values are deferred during the
// parallel section and filled serially afterwards (as in
// fastmr_run_groups_boot_native()), which gives bit-identical results for
// every thread count. Requests that would draw random numbers always run on
// one worker so R's RNG is only touched from the main thread.
template <typename MakePrepared>
std::vector<Result> run_rng_free_jobs(std::size_t jobs, int threads,
                                      double work_rows,
                                      const std::vector<std::string>& methods,
                                      int nboot, double phi, double penk,
                                      MakePrepared prepared) {
  const std::size_t method_count = methods.size();
  std::vector<Result> all(jobs * method_count);
  auto run = [&](std::size_t job) {
    std::vector<Result> results =
      compute_pair(prepared(job), methods, nboot, R_NilValue, true, phi, penk);
    std::move(results.begin(), results.end(),
              all.begin() + static_cast<std::ptrdiff_t>(job * method_count));
  };
  const BootstrapNeeds needs = bootstrap_needs(methods);
  const bool draws = nboot > 0 &&
    (needs.median || needs.egger || needs.penalised || needs.mode);
  const int workers = worthwhile_threads(threads, jobs, work_rows,
                                         min_rows_per_worker(needs));
  if (draws || workers == 1) {
    for (std::size_t job = 0; job < jobs; ++job) run(job);
    return all;
  }
  {
    DeferRMath deferred;
    run_parallel(jobs, workers, run);
  }
  for (Result& result : all) populate_result_pvalues(result);
  return all;
}

} // namespace

// Batched equivalent of calling fastmr_run_native() once per group. Groups are
// CSR-style slices [offsets[g], offsets[g+1]) of the concatenated vectors. Each
// group goes through exactly the same one_pair_from_vectors() filtering and
// compute_pair() code path, so results are bit-identical to the per-group
// calls. Only used when no requested method draws random numbers (bootstrap
// methods use fastmr_run_groups_boot_native()); groups run on up to `threads`
// workers with identical results for every thread count.
// Returns flat vectors in group-major, method-minor order.
// [[Rcpp::export]]
Rcpp::List fastmr_run_groups_native(Rcpp::IntegerVector offsets,
                                    Rcpp::NumericVector exposure_beta,
                                    Rcpp::NumericVector outcome_beta,
                                    Rcpp::NumericVector exposure_se,
                                    Rcpp::NumericVector outcome_se,
                                    Rcpp::CharacterVector methods,
                                    int nboot = 0,
                                    int threads = 1,
                                    double phi = 1.0,
                                    double penk = 20.0) {
  validate_controls(nboot, threads, phi);
  if (!std::isfinite(penk) || penk <= 0.0) Rcpp::stop("penk must be positive and finite");
  const std::vector<std::string> parsed_methods = parse_methods(methods);
  const GroupInputs in = check_group_inputs(offsets, exposure_beta, outcome_beta,
                                            exposure_se, outcome_se);
  const std::vector<Result> all = run_rng_free_jobs(
    static_cast<std::size_t>(in.groups), threads,
    static_cast<double>(in.offsets[in.groups] - in.offsets[0]),
    parsed_methods, nboot, phi, penk,
    [&](std::size_t g) { return group_prepared(in, static_cast<R_xlen_t>(g)); });
  return group_results_to_flat(all);
}

// Leave-one-out fits without materialising the drop-one layouts. Job j fits
// group job_group[j] (1-based) with its job_drop[j]-th row (0-based, -1 for
// none) omitted; an NA group is an empty fit. Each job is bit-identical to
// fastmr_run_groups_native() on the group with that row removed. Only
// RNG-free requests are valid (nboot is 0). Returns flat vectors in
// job-major, method-minor order.
// [[Rcpp::export]]
Rcpp::List fastmr_run_groups_drop_native(Rcpp::IntegerVector offsets,
                                         Rcpp::NumericVector exposure_beta,
                                         Rcpp::NumericVector outcome_beta,
                                         Rcpp::NumericVector exposure_se,
                                         Rcpp::NumericVector outcome_se,
                                         Rcpp::IntegerVector job_group,
                                         Rcpp::IntegerVector job_drop,
                                         Rcpp::CharacterVector methods,
                                         int threads = 1,
                                         double phi = 1.0,
                                         double penk = 20.0) {
  validate_controls(0, threads, phi);
  if (!std::isfinite(penk) || penk <= 0.0) Rcpp::stop("penk must be positive and finite");
  const std::vector<std::string> parsed_methods = parse_methods(methods);
  const GroupInputs in = check_group_inputs(offsets, exposure_beta, outcome_beta,
                                            exposure_se, outcome_se);
  const R_xlen_t jobs = job_group.size();
  if (job_drop.size() != jobs) Rcpp::stop("job_group and job_drop must have equal lengths");
  for (R_xlen_t j = 0; j < jobs; ++j) {
    const int g = job_group[j];
    if (g == NA_INTEGER) continue;
    if (g < 1 || g > in.groups) Rcpp::stop("job_group out of range");
    const R_xlen_t size = in.offsets[g] - in.offsets[g - 1];
    if (job_drop[j] == NA_INTEGER || job_drop[j] < -1 || job_drop[j] >= size) {
      Rcpp::stop("job_drop out of range");
    }
  }
  double drop_rows = 0.0;  // rows fitted across all jobs (work estimate)
  for (R_xlen_t j = 0; j < jobs; ++j) {
    const int g = job_group[j];
    if (g != NA_INTEGER) drop_rows += static_cast<double>(in.offsets[g] - in.offsets[g - 1]);
  }
  const std::vector<int> groups(job_group.begin(), job_group.end());
  const std::vector<int> drops(job_drop.begin(), job_drop.end());
  const std::vector<Result> all = run_rng_free_jobs(
    static_cast<std::size_t>(jobs), threads,
    drop_rows, parsed_methods, 0, phi, penk,
    [&](std::size_t j) {
      Prepared p;
      if (groups[j] == NA_INTEGER) return p;
      const std::size_t g = static_cast<std::size_t>(groups[j] - 1);
      const R_xlen_t begin = in.offsets[g], end = in.offsets[g + 1];
      const R_xlen_t skip = begin + drops[j];
      p.x.reserve(end - begin); p.y.reserve(end - begin);
      p.sx.reserve(end - begin); p.sy.reserve(end - begin);
      for (R_xlen_t i = begin; i < end; ++i) {
        if (drops[j] >= 0 && i == skip) continue;
        const double x = in.x[i], y = in.y[i], sx = in.sx[i], sy = in.sy[i];
        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(sx) ||
            !std::isfinite(sy) || sx <= 0.0 || sy <= 0.0) continue;  // group_prepared() filter
        p.x.push_back(x); p.y.push_back(y); p.sx.push_back(sx); p.sy.push_back(sy);
      }
      return p;
    });
  return group_results_to_flat(all);
}

// sum(x[[g]], na.rm = narm) for each CSR group, with R's own accumulation
// (long double, in order, clamped to +-Inf), so values are identical to the
// per-group sum() calls.
// [[Rcpp::export]]
Rcpp::NumericVector fastmr_group_sum_native(Rcpp::IntegerVector offsets,
                                            Rcpp::NumericVector x, bool narm) {
  const R_xlen_t groups = offsets.size() - 1;
  if (groups < 0 || offsets[0] != 0 || offsets[groups] != x.size()) {
    Rcpp::stop("offsets must start at 0 and end at length(x)");
  }
  Rcpp::NumericVector out(groups);
  for (R_xlen_t g = 0; g < groups; ++g) {
    if (offsets[g + 1] < offsets[g]) Rcpp::stop("offsets must be non-decreasing");
    long double s = 0.0;
    for (R_xlen_t i = offsets[g]; i < offsets[g + 1]; ++i) {
      if (!narm || !ISNAN(x[i])) s += x[i];
    }
    if (s > DBL_MAX) out[g] = R_PosInf;
    else if (s < -DBL_MAX) out[g] = R_NegInf;
    else out[g] = static_cast<double>(s);
  }
  return out;
}

// mean(x[[g]], na.rm = TRUE) for each CSR group, replicating R's real_mean()
// (long double sum, overflow-safe fallback, one refinement pass).
// [[Rcpp::export]]
Rcpp::NumericVector fastmr_group_mean_native(Rcpp::IntegerVector offsets,
                                             Rcpp::NumericVector x) {
  const R_xlen_t groups = offsets.size() - 1;
  if (groups < 0 || offsets[0] != 0 || offsets[groups] != x.size()) {
    Rcpp::stop("offsets must start at 0 and end at length(x)");
  }
  Rcpp::NumericVector out(groups);
  std::vector<double> values;
  for (R_xlen_t g = 0; g < groups; ++g) {
    if (offsets[g + 1] < offsets[g]) Rcpp::stop("offsets must be non-decreasing");
    values.clear();
    for (R_xlen_t i = offsets[g]; i < offsets[g + 1]; ++i) {
      if (!ISNAN(x[i])) values.push_back(x[i]);
    }
    const R_xlen_t n = static_cast<R_xlen_t>(values.size());
    long double s = 0.0;
    for (double v : values) s += v;
    if (R_FINITE(static_cast<double>(s))) {
      s /= n;
    } else {
      long double t = 0.0;
      for (double v : values) t += v / n;
      s = t;
    }
    if (R_FINITE(static_cast<double>(s))) {
      long double t = 0.0;
      for (double v : values) t += (v - s);
      s += t / n;
    }
    out[g] = static_cast<double>(s);
  }
  return out;
}

namespace {

// True when R's normal generator is "Inversion" (the default), whose norm_rand()
// is (R 4.5 src/nmath/snorm.c, unchanged since R 1.7.0):
//   u = unif_rand(); u = (int)(BIG * u) + unif_rand();
//   return qnorm5(u / BIG, 0.0, 1.0, 1, 0);      with BIG = 134217728 (2^27)
// and rnorm(0, 1) returns 0.0 + 1.0 * norm_rand().
bool normal_kind_is_inversion() {
  Rcpp::Function rng_kind("RNGkind", R_BaseNamespace);
  const Rcpp::CharacterVector kinds = rng_kind();
  return kinds.size() >= 2 && Rcpp::as<std::string>(kinds[1]) == "Inversion";
}

// Main-thread part of R::rnorm(0, 1) under "Inversion": draws the two uniforms
// per normal in stream order (consuming the RNG exactly as R::rnorm() would)
// and stores u = (int)(BIG * u1) + u2. inversion_normals_finish() completes
// them; together they produce R::rnorm(0, 1)'s values bit for bit.
void inversion_normals_rng_part(double* out, std::size_t count) {
  for (std::size_t j = 0; j < count; ++j) {
    const double u = unif_rand();
    out[j] = static_cast<int>(134217728.0 * u) + unif_rand();
  }
}

// out[j] = 0 + 1 * qnorm(out[j] / BIG, 0, 1), as rnorm(0, 1) and norm_rand()
// compute it, in parallel. R::qnorm() is Rf_qnorm5() from nmath, the function
// norm_rand() calls: it is pure arithmetic on its arguments (no globals, no
// allocation; its only warning path, ML_WARN_return_NAN with ME_DOMAIN, never
// calls warning(), and p is in (0, 1) here anyway), so it is safe off the main
// thread.
void inversion_normals_finish(double* data, std::size_t count, int workers) {
  const std::size_t chunk = 65536;
  const std::size_t chunks = (count + chunk - 1) / chunk;
  run_parallel(chunks, workers, [&](std::size_t c) {
    const std::size_t end = std::min(count, (c + 1) * chunk);
    for (std::size_t j = c * chunk; j < end; ++j) {
      data[j] = 0.0 + 1.0 * R::qnorm(data[j] / 134217728.0, 0.0, 1.0, 1, 0);
    }
  });
}

} // namespace

// Equivalent of calling fastmr_run_native() once per group with bootstrap
// methods, threaded across groups with bit-identical results. All of R's RNG
// is consumed on the main thread in exactly the serial order: one continuous
// stream from the caller's state, or, when `reseed` is an R function, a fresh
// stream started by reseed(g) (1-based group number) before each group that
// draws. Groups are processed in batches of at most `batch_draws` standard
// normals:
//  * a batch that would run on one worker (threads = 1, or a single group,
//    e.g. one group larger than the budget) streams R::rnorm() straight into
//    the bootstrap layouts with inline p-values, exactly as the per-group
//    calls did, so it needs no draw buffer at all;
//  * otherwise the batch's normals are drawn into one of two reused native
//    buffers, workers rebuild the identical layouts from it (PredrawnNormals)
//    with R-backed p-values deferred, and those p-values are filled serially.
//    Drawing batch b + 1 (main thread) overlaps computing batch b (workers).
// Each buffer never exceeds `batch_draws` doubles. Batching never changes any
// result. Returns flat group-major, method-minor vectors.
// [[Rcpp::export]]
Rcpp::List fastmr_run_groups_boot_native(Rcpp::IntegerVector offsets,
                                         Rcpp::NumericVector exposure_beta,
                                         Rcpp::NumericVector outcome_beta,
                                         Rcpp::NumericVector exposure_se,
                                         Rcpp::NumericVector outcome_se,
                                         Rcpp::CharacterVector methods,
                                         int nboot = 1000,
                                         int threads = 1,
                                         double phi = 1.0,
                                         double penk = 20.0,
                                         SEXP reseed = R_NilValue,
                                         double batch_draws = 8388608.0) {
  validate_controls(nboot, threads, phi);
  if (!std::isfinite(penk) || penk <= 0.0) Rcpp::stop("penk must be positive and finite");
  if (!Rf_isNull(reseed) && !Rf_isFunction(reseed)) Rcpp::stop("reseed must be NULL or a function");
  const std::vector<std::string> parsed_methods = parse_methods(methods);
  const GroupInputs in = check_group_inputs(offsets, exposure_beta, outcome_beta,
                                            exposure_se, outcome_se);
  const std::size_t groups = static_cast<std::size_t>(in.groups);
  const std::size_t method_count = parsed_methods.size();
  const BootstrapNeeds needs = bootstrap_needs(parsed_methods);
  // Per-group draw counts (the number of normals compute_pair() consumes).
  std::vector<double> counts(groups), work(groups);
  for (std::size_t g = 0; g < groups; ++g) {
    const R_xlen_t begin = in.offsets[g], end = in.offsets[g + 1];
    std::size_t snps = 0, ratios = 0;
    for (R_xlen_t i = begin; i < end; ++i) {
      const double x = in.x[i], y = in.y[i], sx = in.sx[i], sy = in.sy[i];
      if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(sx) ||
          !std::isfinite(sy) || sx <= 0.0 || sy <= 0.0) continue;  // group_prepared() filter
      ++snps;
      if (x != 0.0) ++ratios;
    }
    counts[g] = bootstrap_draw_count(snps, ratios, nboot, needs.median,
                                     needs.egger, needs.penalised, needs.mode);
    work[g] = bootstrap_draw_work(snps, ratios, nboot, needs.median,
                                  needs.egger, needs.penalised, needs.mode);
  }
  auto start_stream = [&](std::size_t g) {
    if (Rf_isNull(reseed) || counts[g] <= 0.0) return;
    Rcpp::Function fn(reseed);
    fn(static_cast<int>(g + 1));
  };
  std::vector<Result> all(groups * method_count);
  auto store = [&](std::size_t g, std::vector<Result>& results) {
    std::move(results.begin(), results.end(),
              all.begin() + static_cast<std::ptrdiff_t>(g * method_count));
  };
  // Under normal.kind = "Inversion" (R's default) a normal is qnorm() of two
  // uniforms, and only the uniforms touch the RNG state; see
  // inversion_normals_rng_part().
  int inversion = -1;  // unknown until the first multi-worker batch
  // Multi-worker batches are double-buffered: while workers compute batch b
  // (on a coordinator thread, from buffer b % 2), the main thread draws batch
  // b + 1's normals into the other buffer, then fills batch b - 1's p-values.
  // R's RNG is still consumed only on the main thread, in group order, so the
  // draws, every result and the final RNG state are exactly the serial ones.
  // Workers never touch R: they read pre-drawn normals, run R::qnorm() (pure
  // nmath, see inversion_normals_finish()) and defer R-backed p-values through
  // the thread-local flag, which leaves the main thread free to compute p-values
  // for the previous batch meanwhile.
  struct DrawBuffer {
    std::vector<double> draws;
    std::vector<std::size_t> start;
  };
  DrawBuffer buffers[2];
  struct InFlight {
    std::thread thread;
    std::exception_ptr error;
    std::size_t first = 0, last = 0;
    ~InFlight() { if (thread.joinable()) thread.join(); }  // on unwind only
  } in_flight;
  bool has_in_flight = false;
  // Groups [pending_first, pending_last) are computed but still need their
  // R-backed p-values, which only the main thread fills.
  std::size_t pending_first = 0, pending_last = 0;
  auto join_in_flight = [&]() {
    if (!has_in_flight) return;
    in_flight.thread.join();
    has_in_flight = false;
    if (in_flight.error) std::rethrow_exception(in_flight.error);
    pending_first = in_flight.first;
    pending_last = in_flight.last;
  };
  auto fill_pending = [&]() {
    for (std::size_t k = pending_first * method_count; k < pending_last * method_count; ++k) {
      populate_result_pvalues(all[k]);
    }
    pending_first = pending_last = 0;
  };
  // Batch budget: never above batch_draws, and small enough that a large run
  // has several batches to overlap (batching never changes any result).
  double total_draws = 0.0;
  for (std::size_t g = 0; g < groups; ++g) total_draws += counts[g];
  const double mean_draws = groups > 0 ? total_draws / static_cast<double>(groups) : 0.0;
  const double budget = threads > 1
    ? std::min(batch_draws, std::max({total_draws / 8.0, 1048576.0,
                                      32.0 * static_cast<double>(threads) * mean_draws}))
    : batch_draws;
  std::size_t parity = 0;
  std::size_t first = 0;
  while (first < groups) {
    // Greedy batch: at least one group, then add groups while within budget.
    double batch_total = counts[first], batch_work = work[first];
    std::size_t last = first + 1;
    while (last < groups && batch_total + counts[last] <= budget) {
      batch_total += counts[last];
      batch_work += work[last];
      ++last;
    }
    const std::size_t batch_groups = last - first;
    const double batch_rows = static_cast<double>(in.offsets[last] - in.offsets[first]);
    const int workers = batch_total > 0.0
      ? worthwhile_threads(threads, batch_groups, batch_work + batch_rows, kMinDrawsPerWorker)
      : worthwhile_threads(threads, batch_groups, batch_rows, min_rows_per_worker(needs));
    if (workers == 1) {
      join_in_flight();
      fill_pending();
      for (std::size_t g = first; g < last; ++g) {
        start_stream(g);
        std::vector<Result> results = compute_pair(
          group_prepared(in, static_cast<R_xlen_t>(g)), parsed_methods, nboot,
          R_NilValue, true, phi, penk);
        store(g, results);
      }
    } else {
      DrawBuffer& buffer = buffers[parity];
      parity ^= 1;
      buffer.draws.resize(static_cast<std::size_t>(batch_total));
      buffer.start.resize(batch_groups);
      std::size_t pos = 0;
      for (std::size_t g = first; g < last; ++g) {
        buffer.start[g - first] = pos;
        start_stream(g);
        const std::size_t k = static_cast<std::size_t>(counts[g]);
        double* out = buffer.draws.data() + pos;
        if (inversion < 0) inversion = normal_kind_is_inversion() ? 1 : 0;
        if (inversion) {
          inversion_normals_rng_part(out, k);
        } else {
          for (std::size_t j = 0; j < k; ++j) out[j] = R::rnorm(0.0, 1.0);
        }
        pos += k;
      }
      // At most one batch computes at a time: the previous one must finish
      // before this one starts (its buffer is reused by the next batch).
      join_in_flight();
      in_flight.first = first;
      in_flight.last = last;
      in_flight.error = nullptr;
      const bool finish_normals = inversion == 1;
      DrawBuffer* const data = &buffer;
      in_flight.thread = std::thread([&, data, first, pos, workers, batch_groups, finish_normals]() {
        try {
          if (finish_normals) inversion_normals_finish(data->draws.data(), pos, workers);
          const double* draw_data = data->draws.data();
          run_parallel(batch_groups, workers, [&](std::size_t index) {
            struct DeferHere {
              bool previous = defer_r_math_here;
              DeferHere() { defer_r_math_here = true; }
              ~DeferHere() { defer_r_math_here = previous; }
            } deferred;
            const std::size_t g = first + index;
            std::vector<Result> results = compute_pair(
              group_prepared(in, static_cast<R_xlen_t>(g)), parsed_methods, nboot,
              R_NilValue, true, phi, penk, draw_data + data->start[index]);
            store(g, results);
          });
        } catch (...) {
          in_flight.error = std::current_exception();
        }
      });
      has_in_flight = true;
      fill_pending();  // the previous batch's p-values, while this one runs
    }
    first = last;
  }
  join_in_flight();
  fill_pending();
  return group_results_to_flat(all);
}
