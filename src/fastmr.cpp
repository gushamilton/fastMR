#include <Rcpp.h>
#include <R_ext/BLAS.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <thread>
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
const double NA_VALUE = std::numeric_limits<double>::quiet_NaN();

double finite_or_na(double x) {
  return std::isfinite(x) ? x : NA_VALUE;
}

double z_pvalue(double statistic) {
  if (!std::isfinite(statistic)) return NA_VALUE;
  return 2.0 * R::pnorm5(std::abs(statistic), 0.0, 1.0, false, false);
}

double t_pvalue(double statistic, int df) {
  if (!std::isfinite(statistic) || df <= 0) return NA_VALUE;
  return 2.0 * R::pt(std::abs(statistic), static_cast<double>(df), false, false);
}

double chi_square_pvalue(double q, int df) {
  if (!std::isfinite(q) || df <= 0) return NA_VALUE;
  return R::pchisq(q, static_cast<double>(df), false, false);
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

double weighted_median_point(const std::vector<double>& values,
                             const std::vector<double>& weights) {
  if (values.empty() || values.size() != weights.size()) return NA_VALUE;
  std::vector<std::size_t> order(values.size());
  std::iota(order.begin(), order.end(), static_cast<std::size_t>(0));
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    return values[a] < values[b];
  });
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
  std::vector<std::complex<double>> forward_steps;
  std::vector<std::complex<double>> inverse_steps;

  explicit FFTPlan(std::size_t size) : n(size), bit_reverse(size) {
    for (std::size_t i = 1, j = 0; i < n; ++i) {
      std::size_t bit = n >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
      j ^= bit;
      bit_reverse[i] = j;
    }
    for (std::size_t length = 2; length <= n; length <<= 1) {
      const double angle = 2.0 * 3.14159265358979323846 / static_cast<double>(length);
      forward_steps.emplace_back(std::cos(-angle), std::sin(-angle));
      inverse_steps.emplace_back(std::cos(angle), std::sin(angle));
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

void fft_inplace(std::vector<std::complex<double>>& values, bool inverse) {
  const std::size_t n = values.size();
  const FFTPlan& plan = fft_plan(n);
  for (std::size_t i = 1; i < n; ++i) {
    const std::size_t j = plan.bit_reverse[i];
    if (i < j) std::swap(values[i], values[j]);
  }
  for (std::size_t stage = 0, length = 2; length <= n; ++stage, length <<= 1) {
    const std::complex<double> step = inverse ? plan.inverse_steps[stage]
                                              : plan.forward_steps[stage];
    for (std::size_t start = 0; start < n; start += length) {
      std::complex<double> factor(1.0, 0.0);
      const std::size_t half = length >> 1;
      for (std::size_t i = 0; i < half; ++i) {
        const std::complex<double> even = values[start + i];
        const std::complex<double> odd = factor * values[start + i + half];
        values[start + i] = even + odd;
        values[start + i + half] = even - odd;
        factor *= step;
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
};

ModeDensityWorkspace& mode_workspace() {
  static thread_local ModeDensityWorkspace workspace;
  return workspace;
}

double mode_point_r_density(const double* values, const double* weights,
                            std::size_t count, double phi) {
  if (count == 0) return NA_VALUE;
  for (std::size_t i = 0; i < count; ++i) {
    if (!std::isfinite(values[i]) || !std::isfinite(weights[i]) || weights[i] < 0.0) return NA_VALUE;
  }
  ModeDensityWorkspace& workspace = mode_workspace();
  std::vector<double>& scratch = workspace.scratch;
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
  const int length = 2 * n;
  const double lo = from - 4.0 * bandwidth;
  const double up = to + 4.0 * bandwidth;
  const double delta = (up - lo) / static_cast<double>(n - 1);
  std::vector<std::complex<double>>& binned = workspace.binned;
  binned.assign(length, std::complex<double>(0.0, 0.0));
  for (std::size_t i = 0; i < count; ++i) {
    const double xpos = (values[i] - lo) / delta;
    if (!std::isfinite(xpos) || xpos > static_cast<double>(std::numeric_limits<int>::max()) ||
        xpos < static_cast<double>(std::numeric_limits<int>::min())) continue;
    const int index = static_cast<int>(std::floor(xpos));
    const double fraction = xpos - static_cast<double>(index);
    if (0 <= index && index <= n - 2) {
      binned[index] += (1.0 - fraction) * weights[i];
      binned[index + 1] += fraction * weights[i];
    } else if (index == -1) {
      binned[0] += fraction * weights[i];
    } else if (index == n - 1) {
      binned[index] += (1.0 - fraction) * weights[i];
    }
  }
  std::vector<std::complex<double>>& kernel = workspace.kernel;
  kernel.resize(length);
  for (int i = 0; i < length; ++i) {
    const double distance = (i <= n) ? static_cast<double>(i) * delta
                                     : -static_cast<double>(length - i) * delta;
    const double z = distance / bandwidth;
    kernel[i] = std::exp(-0.5 * z * z) / (bandwidth * std::sqrt(2.0 * 3.14159265358979323846));
  }
  fft_inplace(binned, false);
  fft_inplace(kernel, false);
  for (int i = 0; i < length; ++i) binned[i] *= std::conj(kernel[i]);
  fft_inplace(binned, true);
  int best_index = 0;
  double best_density = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < n; ++i) {
    const double x = from + (to - from) * static_cast<double>(i) / static_cast<double>(n - 1);
    const double position = (x - lo) / delta;
    int left = static_cast<int>(std::floor(position));
    double density = 0.0;
    if (left < 0) density = std::max(0.0, binned[0].real());
    else if (left >= n - 1) density = std::max(0.0, binned[n - 1].real());
    else {
      const double fraction = position - static_cast<double>(left);
      density = (1.0 - fraction) * binned[left].real() + fraction * binned[left + 1].real();
      density = std::max(0.0, density);
    }
    if (density > best_density) {
      best_density = density;
      best_index = i;
    }
  }
  return from + (to - from) * static_cast<double>(best_index) / static_cast<double>(n - 1);
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
  ModeDensityWorkspace& workspace = mode_workspace();
  std::vector<double>& scratch = workspace.scratch;
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
  const int length = 2 * n;
  const double lo = from - 4.0 * bandwidth;
  const double up = to + 4.0 * bandwidth;
  const double delta = (up - lo) / static_cast<double>(n - 1);
  std::vector<std::complex<double>>& simple = workspace.simple;
  std::vector<std::complex<double>>& weighted = workspace.weighted;
  simple.assign(length, std::complex<double>(0.0, 0.0));
  weighted.assign(length, std::complex<double>(0.0, 0.0));
  auto bin = [&](std::vector<std::complex<double>>& target, const double* source) {
    for (std::size_t i = 0; i < count; ++i) {
      const double xpos = (values[i] - lo) / delta;
      if (!std::isfinite(xpos) || xpos > static_cast<double>(std::numeric_limits<int>::max()) ||
          xpos < static_cast<double>(std::numeric_limits<int>::min())) continue;
      const int index = static_cast<int>(std::floor(xpos));
      const double fraction = xpos - static_cast<double>(index);
      if (0 <= index && index <= n - 2) {
        target[index] += (1.0 - fraction) * source[i];
        target[index + 1] += fraction * source[i];
      } else if (index == -1) {
        target[0] += fraction * source[i];
      } else if (index == n - 1) {
        target[index] += (1.0 - fraction) * source[i];
      }
    }
  };
  bin(simple, simple_weights);
  bin(weighted, weighted_weights);
  std::vector<std::complex<double>>& kernel = workspace.kernel;
  kernel.resize(length);
  for (int i = 0; i < length; ++i) {
    const double distance = (i <= n) ? static_cast<double>(i) * delta
                                     : -static_cast<double>(length - i) * delta;
    const double z = distance / bandwidth;
    kernel[i] = std::exp(-0.5 * z * z) / (bandwidth * std::sqrt(2.0 * 3.14159265358979323846));
  }
  fft_inplace(simple, false);
  fft_inplace(weighted, false);
  fft_inplace(kernel, false);
  for (int i = 0; i < length; ++i) {
    simple[i] *= std::conj(kernel[i]);
    weighted[i] *= std::conj(kernel[i]);
  }
  fft_inplace(simple, true);
  fft_inplace(weighted, true);
  auto find_max = [&](const std::vector<std::complex<double>>& density_values) {
    int best_index = 0;
    double best_density = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < n; ++i) {
      const double x = from + (to - from) * static_cast<double>(i) / static_cast<double>(n - 1);
      const double position = (x - lo) / delta;
      const int left = static_cast<int>(std::floor(position));
      double density = 0.0;
      if (left < 0) density = std::max(0.0, density_values[0].real());
      else if (left >= n - 1) density = std::max(0.0, density_values[n - 1].real());
      else {
        const double fraction = position - static_cast<double>(left);
        density = (1.0 - fraction) * density_values[left].real() +
                  fraction * density_values[left + 1].real();
        density = std::max(0.0, density);
      }
      if (density > best_density) {
        best_density = density;
        best_index = i;
      }
    }
    return from + (to - from) * static_cast<double>(best_index) / static_cast<double>(n - 1);
  };
  return std::make_pair(find_max(simple), find_max(weighted));
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

Result compute_ivw(const Prepared& p, const std::string& method) {
  const int n = static_cast<int>(p.x.size());
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
    se = sigma > 0.0 && std::isfinite(sigma)
      ? residual_se / std::min(1.0, sigma)
      : residual_se;
  } else if (method == "ivw_fe") {
    se = sigma > 0.0 && std::isfinite(sigma) ? residual_se / sigma : NA_VALUE;
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
  double pval = 2.0 * R::pbinom(static_cast<double>(lower), static_cast<double>(n), 0.5, true, false);
  pval = std::min(1.0, pval);
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
  const double beta = weighted_median_point(p.ratio, weights);
  double se = NA_VALUE;
  if (nboot > 0 && !p.bootstrap.empty()) {
    std::vector<double> estimates(nboot);
    for (int draw = 0; draw < nboot; ++draw) {
      std::vector<double> row(p.bootstrap.begin() + static_cast<std::size_t>(draw) * n,
                              p.bootstrap.begin() + static_cast<std::size_t>(draw + 1) * n);
      estimates[draw] = weighted_median_point(row, weights);
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

Result compute_penalised_median(const Prepared& p, int nboot, double penk) {
  const int n = static_cast<int>(p.ratio.size());
  if (n < 3) return empty_result("penalised_weighted_median", n);
  std::vector<double> weights(n), penalised_weights(n);
  for (int i = 0; i < n; ++i) weights[i] = 1.0 / (p.ratio_se[i] * p.ratio_se[i]);
  const double weighted_beta = weighted_median_point(p.ratio, weights);
  for (int i = 0; i < n; ++i) {
    const double statistic = weights[i] * (p.ratio[i] - weighted_beta) *
                             (p.ratio[i] - weighted_beta);
    const double penalty = R::pchisq(statistic, 1.0, false, false);
    penalised_weights[i] = weights[i] * std::min(1.0, penalty * penk);
  }
  const double beta = weighted_median_point(p.ratio, penalised_weights);
  double se = NA_VALUE;
  const std::vector<double>& bootstrap = p.penalised_bootstrap.empty()
    ? p.bootstrap : p.penalised_bootstrap;
  if (nboot > 0 && !bootstrap.empty()) {
    std::vector<double> estimates(nboot);
    for (int draw = 0; draw < nboot; ++draw) {
      const double* row = bootstrap.data() + static_cast<std::size_t>(draw) * n;
      estimates[draw] = weighted_median_point(
        std::vector<double>(row, row + n), penalised_weights);
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
  Result result = empty_result(method, n);
  result.beta = beta;
  result.se = se;
  result.pval = t_pvalue(safe_statistic(beta, se), n - 1);
  result.bootstrap = true;
  result.bootstrap_value = nboot;
  result.phi = true;
  result.phi_value = 1.0;
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

void make_bootstrap(Prepared& p, int nboot, SEXP seed,
                    bool needs_median, bool needs_egger,
                    bool needs_penalised) {
  if (nboot <= 0 || p.ratio.size() < 3 ||
      (!needs_median && !needs_egger && !needs_penalised)) return;
  (void) seed;
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
        exp_z[static_cast<std::size_t>(draw) * p.x.size() + snp] = R::rnorm(0.0, 1.0);
      }
    }
    for (std::size_t snp = 0; snp < p.x.size(); ++snp) {
      for (int draw = 0; draw < nboot; ++draw) {
        out_z[static_cast<std::size_t>(draw) * p.x.size() + snp] = R::rnorm(0.0, 1.0);
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

void make_mode_bootstrap(Prepared& p, int nboot) {
  if (nboot <= 0 || p.ratio.size() < 3) return;
  const std::size_t n = p.ratio.size();
  p.mode_bootstrap.resize(static_cast<std::size_t>(nboot) * n);
  // TwoSampleMR::mr_mode draws delta-method Wald ratios directly.
  for (std::size_t snp = 0; snp < n; ++snp) {
    for (int draw = 0; draw < nboot; ++draw) {
      p.mode_bootstrap[static_cast<std::size_t>(draw) * n + snp] =
        p.ratio[snp] + p.ratio_se[snp] * R::rnorm(0.0, 1.0);
    }
  }
}

std::vector<Result> compute_pair(Prepared p,
                                 const std::vector<std::string>& methods,
                                 int nboot, SEXP seed,
                                 bool prepare_bootstrap = true,
                                 double phi = 1.0,
                                 double penk = 20.0) {
  prepare_ratios(p);
  if (prepare_bootstrap) {
    bool needs_median = false;
    bool needs_penalised = false;
    bool needs_mode = false;
    bool needs_egger = false;
    for (const std::string& method : methods) {
      needs_median = needs_median || method == "simple_median" ||
                     method == "weighted_median";
      needs_penalised = needs_penalised || method == "penalised_weighted_median";
      needs_mode = needs_mode || method == "simple_mode" || method == "weighted_mode";
      needs_egger = needs_egger || method == "egger_bootstrap";
    }
    if (needs_median || needs_egger || needs_penalised) {
      make_bootstrap(p, nboot, seed, needs_median, needs_egger, needs_penalised);
    }
    if (needs_mode) make_mode_bootstrap(p, nboot);
  }
  std::vector<Result> result;
  result.resize(methods.size());
  bool has_simple = false;
  bool has_weighted = false;
  std::size_t simple_index = 0;
  std::size_t weighted_index = 0;
  for (std::size_t i = 0; i < methods.size(); ++i) {
    if (methods[i] == "simple_mode") { has_simple = true; simple_index = i; }
    if (methods[i] == "weighted_mode") { has_weighted = true; weighted_index = i; }
  }
  for (std::size_t i = 0; i < methods.size(); ++i) {
    const std::string& method = methods[i];
    if ((method == "simple_mode" || method == "weighted_mode") && has_simple && has_weighted) continue;
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

bool only_ivw_methods(const std::vector<std::string>& methods) {
  if (methods.empty()) return false;
  for (const std::string& method : methods) {
    if (method != "ivw" && method != "ivw_fe" && method != "ivw_mre") return false;
  }
  return true;
}

std::vector<std::vector<Result>> compute_ivw_grid_blas(
    const GridData& grid, const std::vector<std::string>& methods) {
  const int exposure_count = grid.exposure_count;
  const int outcome_count = grid.outcome_count;
  const int snp_count = grid.snp_count;
  const std::size_t pair_count = static_cast<std::size_t>(exposure_count) *
                                 static_cast<std::size_t>(outcome_count);
  std::vector<std::vector<Result>> results(pair_count);
  if (snp_count < 2) {
    for (std::size_t pair = 0; pair < pair_count; ++pair) {
      results[pair].resize(methods.size());
      for (std::size_t i = 0; i < methods.size(); ++i) {
        results[pair][i] = empty_result(methods[i], snp_count);
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
        // BLAS reduction order can leave a tiny negative residue at an exact
        // fit. Preserve the scalar path's non-negative RSS contract.
        if (rss < 0.0 && rss > -1e-12 * std::max(1.0, y2)) rss = 0.0;
        if (rss >= 0.0 && std::isfinite(rss)) {
          sigma = std::sqrt(rss / static_cast<double>(snp_count - 1));
          base_se = std::sqrt(1.0 / den);
        }
      }
      const double residual_se = std::isfinite(sigma) && std::isfinite(base_se)
        ? base_se * sigma : NA_VALUE;
      results[static_cast<std::size_t>(exposure) * outcome_count + outcome].resize(methods.size());
      for (std::size_t i = 0; i < methods.size(); ++i) {
        const std::string& method = methods[i];
        Result result = empty_result(method, snp_count);
        result.beta = beta_value;
        if (method == "ivw") {
          result.se = sigma > 0.0 && std::isfinite(sigma)
            ? residual_se / std::min(1.0, sigma) : residual_se;
        } else if (method == "ivw_fe") {
          result.se = sigma > 0.0 && std::isfinite(sigma) ? residual_se / sigma : NA_VALUE;
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
        results[static_cast<std::size_t>(exposure) * outcome_count + outcome][i] = result;
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
          if (rss < 0.0 && rss > -1e-12 * std::max(1.0, y2)) rss = 0.0;
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
            method_se = sigma_value > 0.0 && std::isfinite(sigma_value)
              ? residual_se / std::min(1.0, sigma_value) : residual_se;
          } else if (method == "ivw_fe") {
            method_se = sigma_value > 0.0 && std::isfinite(sigma_value)
              ? residual_se / sigma_value : NA_VALUE;
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

int bounded_threads(int requested, std::size_t jobs) {
  const int bounded = std::max(1, std::min<int>(requested, static_cast<int>(std::max<std::size_t>(1, jobs))));
  return bounded;
}

} // namespace

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
  std::vector<std::vector<Result>> results(pair_count);
  const int thread_count = bounded_threads(threads, pair_count);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(thread_count)
  for (int index = 0; index < static_cast<int>(pair_count); ++index) {
    const int exposure = index / grid.outcome_count;
    const int outcome = index % grid.outcome_count;
    results[static_cast<std::size_t>(index)] = compute_pair(
      pair_from_grid(grid, exposure, outcome, nboot), parsed_methods, nboot,
      R_NilValue, false, phi, penk);
  }
#else
  if (thread_count == 1) {
    for (std::size_t index = 0; index < pair_count; ++index) {
      const int exposure = static_cast<int>(index / static_cast<std::size_t>(grid.outcome_count));
      const int outcome = static_cast<int>(index % static_cast<std::size_t>(grid.outcome_count));
      results[index] = compute_pair(pair_from_grid(grid, exposure, outcome, nboot), parsed_methods, nboot, R_NilValue, false, phi, penk);
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
          results[index] = compute_pair(pair_from_grid(grid, exposure, outcome, nboot), parsed_methods, nboot, R_NilValue, false, phi, penk);
        }
      });
    }
    for (std::thread& worker : pool) worker.join();
  }
#endif

  Rcpp::List output(pair_count);
  for (std::size_t i = 0; i < pair_count; ++i) output[i] = results_to_list(results[i]);
  return output;
}
