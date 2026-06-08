#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace py = pybind11;

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

struct ProposalResult {
  std::vector<double> tess;
  int n_centres = 0;
  std::vector<int> dim;
  std::string modification = "Change";
};

struct ReducedMetric {
  std::vector<int> metric;
  std::vector<int> member_counts;
};

bool in_vector(int value, const std::vector<int>& values) {
  return std::find(values.begin(), values.end(), value) != values.end();
}

double period_shift(double value, double limit) {
  while (value >= limit) {
    value -= 2.0 * limit;
  }
  while (value < -limit) {
    value += 2.0 * limit;
  }
  return value;
}

bool is_last_member_column(int index, const std::vector<int>& members) {
  return index == static_cast<int>(members.size()) - 1 || members[index + 1] != members[index];
}

double uniform01(std::mt19937_64& rng) {
  static constexpr double min_open = std::numeric_limits<double>::min();
  std::uniform_real_distribution<double> dist(min_open, 1.0);
  return dist(rng);
}

double log_binomial_pmf(int k, int n, double p) {
  if (k < 0 || k > n) {
    return -std::numeric_limits<double>::infinity();
  }
  return std::lgamma(n + 1.0) - std::lgamma(k + 1.0) - std::lgamma(n - k + 1.0) +
         k * std::log(p) + (n - k) * std::log1p(-p);
}

double log_poisson_pmf(int k, double lambda) {
  if (k < 0) {
    return -std::numeric_limits<double>::infinity();
  }
  return k * std::log(lambda) - lambda - std::lgamma(k + 1.0);
}

double euclidean_distance(const std::vector<double>& first,
                          const std::vector<double>& second,
                          int offset,
                          int size) {
  double total = 0.0;
  for (int idx = 0; idx < size; ++idx) {
    const double diff = first[offset + idx] - second[offset + idx];
    total += diff * diff;
  }
  return total;
}

double spherical_distance(const std::vector<double>& first,
                          const std::vector<double>& second,
                          int offset,
                          int size) {
  if (size == 1) {
    const double a1 = std::abs(first[offset] - second[offset]);
    const double a2 = 2.0 * kPi - a1;
    return std::min(a1, a2) * std::min(a1, a2);
  }

  double angle_diff = std::cos(first[offset + size - 1] - second[offset + size - 1]);
  for (int idx = size - 2; idx >= 0; --idx) {
    double internal = std::sin(first[offset + idx]) * std::sin(second[offset + idx]) +
                      std::cos(first[offset + idx]) * std::cos(second[offset + idx]) * angle_diff;
    internal = std::clamp(internal, -1.0, 1.0);
    angle_diff = (idx == 0) ? std::acos(internal) : internal;
  }
  return angle_diff * angle_diff;
}

double calc_distance(const std::vector<double>& first,
                     const std::vector<double>& second,
                     const std::vector<int>& member_counts,
                     const std::vector<int>& metric) {
  int offset = 0;
  double total = 0.0;
  for (int group = 0; group < static_cast<int>(member_counts.size()); ++group) {
    const int size = member_counts[group];
    if (metric[group] == 0) {
      total += euclidean_distance(first, second, offset, size);
    } else if (metric[group] == 1) {
      total += spherical_distance(first, second, offset, size);
    }
    offset += size;
  }
  return total;
}

ReducedMetric make_reduced_metric(const std::vector<int>& metric, const std::vector<int>& members) {
  if (metric.size() != members.size()) {
    throw std::invalid_argument("metric and members must have the same length.");
  }

  ReducedMetric reduced;
  int idx = 0;
  while (idx < static_cast<int>(members.size())) {
    const int current_member = members[idx];
    int count = 0;
    while (idx + count < static_cast<int>(members.size()) && members[idx + count] == current_member) {
      ++count;
    }
    reduced.metric.push_back(metric[idx]);
    reduced.member_counts.push_back(count);
    idx += count;
  }
  return reduced;
}

std::vector<int> knn1_internal(const double* x,
                               int n,
                               int p,
                               const std::vector<double>& centres,
                               int n_centres,
                               int d,
                               const std::vector<int>& dim,
                               const std::vector<int>& metric_red,
                               const std::vector<int>& member_red) {
  std::vector<int> result(n, 0);
  if (n_centres == 1) {
    return result;
  }

  std::vector<double> query_point(p);
  std::vector<double> tess_point(p);
  for (int obs = 0; obs < n; ++obs) {
    for (int col = 0; col < p; ++col) {
      query_point[col] = x[obs * p + col];
      tess_point[col] = query_point[col];
    }

    double best_distance = std::numeric_limits<double>::infinity();
    int best_centre = 0;
    for (int centre = 0; centre < n_centres; ++centre) {
      for (int local_dim = 0; local_dim < d; ++local_dim) {
        tess_point[dim[local_dim]] = centres[centre * d + local_dim];
      }
      const double distance = calc_distance(query_point, tess_point, member_red, metric_red);
      if (distance < best_distance) {
        best_distance = distance;
        best_centre = centre;
      }
    }
    result[obs] = best_centre;
  }

  return result;
}

void aggregate_residuals(const std::vector<double>& residuals,
                         const std::vector<int>& idx_old,
                         int old_centres,
                         const std::vector<int>& idx_new,
                         int new_centres,
                         std::vector<double>& r_old,
                         std::vector<int>& n_old,
                         std::vector<double>& r_new,
                         std::vector<int>& n_new) {
  r_old.assign(old_centres, 0.0);
  n_old.assign(old_centres, 0);
  r_new.assign(new_centres, 0.0);
  n_new.assign(new_centres, 0);

  for (int obs = 0; obs < static_cast<int>(residuals.size()); ++obs) {
    r_old[idx_old[obs]] += residuals[obs];
    n_old[idx_old[obs]] += 1;
    r_new[idx_new[obs]] += residuals[obs];
    n_new[idx_new[obs]] += 1;
  }
}

double log_acceptance_probability(const std::vector<double>& r_old,
                                  const std::vector<int>& n_old,
                                  const std::vector<double>& r_new,
                                  const std::vector<int>& n_new,
                                  int d_new,
                                  int n_centres_new,
                                  double sigma_squared,
                                  double sigma_squared_mu,
                                  double omega,
                                  double lambda_rate,
                                  int p,
                                  const std::string& modification) {
  double sum_log_old = 0.0;
  double sum_r_old = 0.0;
  for (int idx = 0; idx < static_cast<int>(n_old.size()); ++idx) {
    const double den = n_old[idx] * sigma_squared_mu + sigma_squared;
    sum_log_old += std::log(den);
    sum_r_old += (r_old[idx] * r_old[idx]) / den;
  }

  double sum_log_new = 0.0;
  double sum_r_new = 0.0;
  for (int idx = 0; idx < static_cast<int>(n_new.size()); ++idx) {
    const double den = n_new[idx] * sigma_squared_mu + sigma_squared;
    sum_log_new += std::log(den);
    sum_r_new += (r_new[idx] * r_new[idx]) / den;
  }

  const double log_lik = 0.5 * (sum_log_old - sum_log_new) +
                         (sigma_squared_mu / (2.0 * sigma_squared)) * (sum_r_new - sum_r_old);
  const double prob = std::min(1.0 - 1e-10, std::max(0.0, omega / static_cast<double>(p)));
  double out = log_lik;

  if (modification == "AD") {
    out += log_binomial_pmf(d_new - 1, p - 1, prob) -
           log_binomial_pmf(d_new - 2, p - 1, prob) -
           std::log(static_cast<double>(d_new));
    if (d_new == 1) {
      out += std::log(0.5);
    } else if (d_new == p - 1) {
      out += std::log(2.0);
    }
  } else if (modification == "RD") {
    out += log_binomial_pmf(d_new - 1, p, prob) -
           log_binomial_pmf(d_new, p, prob) +
           std::log(static_cast<double>(d_new + 1));
    if (d_new == p) {
      out += std::log(0.5);
    } else if (d_new == 2) {
      out += std::log(2.0);
    }
  } else if (modification == "AC") {
    out += log_poisson_pmf(n_centres_new - 1, lambda_rate) -
           log_poisson_pmf(n_centres_new - 2, lambda_rate) -
           std::log(static_cast<double>(n_centres_new)) +
           0.5 * std::log(sigma_squared);
    if (n_centres_new == 1) {
      out += std::log(0.5);
    }
  } else if (modification == "RC") {
    out += log_poisson_pmf(n_centres_new - 1, lambda_rate) -
           log_poisson_pmf(n_centres_new, lambda_rate) +
           std::log(static_cast<double>(n_centres_new + 1)) -
           0.5 * std::log(sigma_squared);
    if (n_centres_new == 2) {
      out += std::log(2.0);
    }
  }

  return out;
}

std::vector<double> sample_mu(const std::vector<double>& r_cell,
                              const std::vector<int>& n_cell,
                              double sigma_squared_mu,
                              double sigma_squared,
                              std::mt19937_64& rng) {
  std::vector<double> out(r_cell.size());
  std::normal_distribution<double> normal(0.0, 1.0);
  for (int idx = 0; idx < static_cast<int>(r_cell.size()); ++idx) {
    const double den = sigma_squared_mu * n_cell[idx] + sigma_squared;
    const double mean = (sigma_squared_mu * r_cell[idx]) / den;
    const double sd = std::sqrt((sigma_squared * sigma_squared_mu) / den);
    out[idx] = mean + normal(rng) * sd;
  }
  return out;
}

ProposalResult propose(const std::vector<double>& tess,
                       int n_centres,
                       int d,
                       const std::vector<int>& dim,
                       int p,
                       const std::vector<double>& proposal_sd,
                       const std::vector<double>& proposal_mu,
                       const std::vector<int>& metric,
                       const std::vector<int>& members,
                       std::mt19937_64& rng) {
  ProposalResult result;
  result.tess = tess;
  result.n_centres = n_centres;
  result.dim = dim;
  result.modification = "Change";

  const double choice = uniform01(rng);
  std::normal_distribution<double> normal(0.0, 1.0);

  auto sample_coordinate = [&](int global_dim) {
    double value = proposal_mu[global_dim] + normal(rng) * proposal_sd[global_dim];
    if (metric[global_dim] == 1 && is_last_member_column(global_dim, members)) {
      value = period_shift(value, kPi);
    }
    return value;
  };

  if ((choice < 0.2 && d != p) || (d == 1 && d != p && choice < 0.4)) {
    result.modification = "AD";
    int new_dim = 0;
    std::uniform_int_distribution<int> dim_dist(0, p - 1);
    do {
      new_dim = dim_dist(rng);
    } while (in_vector(new_dim, result.dim));

    result.dim.push_back(new_dim);
    result.tess.assign(n_centres * (d + 1), 0.0);
    for (int row = 0; row < n_centres; ++row) {
      for (int col = 0; col < d; ++col) {
        result.tess[row * (d + 1) + col] = tess[row * d + col];
      }
      result.tess[row * (d + 1) + d] = sample_coordinate(new_dim);
    }
  } else if (choice < 0.4 && d > 1) {
    result.modification = "RD";
    std::uniform_int_distribution<int> remove_dist(0, d - 1);
    const int remove_idx = remove_dist(rng);
    result.dim.erase(result.dim.begin() + remove_idx);
    result.tess.assign(n_centres * (d - 1), 0.0);
    for (int row = 0; row < n_centres; ++row) {
      int out_col = 0;
      for (int col = 0; col < d; ++col) {
        if (col == remove_idx) {
          continue;
        }
        result.tess[row * (d - 1) + out_col] = tess[row * d + col];
        ++out_col;
      }
    }
  } else if (choice < 0.6 || (choice < 0.8 && n_centres == 1)) {
    result.modification = "AC";
    result.n_centres = n_centres + 1;
    result.tess.assign(result.n_centres * d, 0.0);
    for (int row = 0; row < n_centres; ++row) {
      for (int col = 0; col < d; ++col) {
        result.tess[row * d + col] = tess[row * d + col];
      }
    }
    for (int col = 0; col < d; ++col) {
      result.tess[n_centres * d + col] = sample_coordinate(dim[col]);
    }
  } else if (choice < 0.8 && n_centres > 1) {
    result.modification = "RC";
    std::uniform_int_distribution<int> remove_dist(0, n_centres - 1);
    const int remove_row = remove_dist(rng);
    result.n_centres = n_centres - 1;
    result.tess.reserve(result.n_centres * d);
    for (int row = 0; row < n_centres; ++row) {
      if (row == remove_row) {
        continue;
      }
      for (int col = 0; col < d; ++col) {
        result.tess.push_back(tess[row * d + col]);
      }
    }
  } else if (choice < 0.9 || d == p) {
    std::uniform_int_distribution<int> centre_dist(0, n_centres - 1);
    const int centre = centre_dist(rng);
    for (int col = 0; col < d; ++col) {
      result.tess[centre * d + col] = sample_coordinate(dim[col]);
    }
  } else {
    result.modification = "Swap";
    std::uniform_int_distribution<int> local_dim_dist(0, d - 1);
    std::uniform_int_distribution<int> global_dim_dist(0, p - 1);
    const int local_dim = local_dim_dist(rng);
    int new_dim = 0;
    do {
      new_dim = global_dim_dist(rng);
    } while (in_vector(new_dim, result.dim));
    result.dim[local_dim] = new_dim;
    for (int row = 0; row < n_centres; ++row) {
      result.tess[row * d + local_dim] = sample_coordinate(new_dim);
    }
  }

  return result;
}

std::vector<double> copy_double_array(py::handle handle, int expected_ndim, std::vector<ssize_t>* shape = nullptr) {
  py::array_t<double, py::array::c_style | py::array::forcecast> arr = py::cast<py::array_t<double>>(handle);
  py::buffer_info info = arr.request();
  if (info.ndim != expected_ndim) {
    throw std::invalid_argument("Unexpected array dimensionality.");
  }
  if (shape != nullptr) {
    shape->assign(info.shape.begin(), info.shape.end());
  }
  const auto* ptr = static_cast<const double*>(info.ptr);
  return std::vector<double>(ptr, ptr + info.size);
}

std::vector<int> copy_int_array(py::handle handle) {
  py::array_t<int, py::array::c_style | py::array::forcecast> arr = py::cast<py::array_t<int>>(handle);
  py::buffer_info info = arr.request();
  if (info.ndim != 1) {
    throw std::invalid_argument("Expected a one-dimensional integer array.");
  }
  const auto* ptr = static_cast<const int*>(info.ptr);
  return std::vector<int>(ptr, ptr + info.size);
}

py::array_t<double> make_matrix(const std::vector<double>& values, int rows, int cols) {
  py::array_t<double> arr({rows, cols});
  auto out = arr.mutable_unchecked<2>();
  for (int row = 0; row < rows; ++row) {
    for (int col = 0; col < cols; ++col) {
      out(row, col) = values[row * cols + col];
    }
  }
  return arr;
}

py::array_t<double> make_vector(const std::vector<double>& values) {
  py::array_t<double> arr(static_cast<py::ssize_t>(values.size()));
  auto out = arr.mutable_unchecked<1>();
  for (int idx = 0; idx < static_cast<int>(values.size()); ++idx) {
    out(idx) = values[idx];
  }
  return arr;
}

py::array_t<int> make_int_vector(const std::vector<int>& values) {
  py::array_t<int> arr(static_cast<py::ssize_t>(values.size()));
  auto out = arr.mutable_unchecked<1>();
  for (int idx = 0; idx < static_cast<int>(values.size()); ++idx) {
    out(idx) = values[idx];
  }
  return arr;
}

}  // namespace

py::dict run_mcmc(py::array_t<double, py::array::c_style | py::array::forcecast> x_scaled_arr,
                  py::array_t<double, py::array::c_style | py::array::forcecast> y_scaled_arr,
                  py::array_t<int, py::array::c_style | py::array::forcecast> metric_arr,
                  py::array_t<int, py::array::c_style | py::array::forcecast> member_arr,
                  int m,
                  int total_iter,
                  int burn_in,
                  int thinning,
                  double nu,
                  double lambda,
                  double sigma_squared_mu,
                  double omega,
                  double lambda_rate,
                  py::array_t<double, py::array::c_style | py::array::forcecast> proposal_sd_arr,
                  py::array_t<double, py::array::c_style | py::array::forcecast> proposal_mu_arr,
                  py::list init_tess,
                  py::list init_dim,
                  py::list init_pred,
                  py::array_t<int, py::array::c_style | py::array::forcecast> binary_cols_arr,
                  double cat_scaling,
                  std::uint64_t seed,
                  bool verbose) {
  py::buffer_info x_info = x_scaled_arr.request();
  py::buffer_info y_info = y_scaled_arr.request();
  if (x_info.ndim != 2 || y_info.ndim != 1) {
    throw std::invalid_argument("x_scaled must be 2D and y_scaled must be 1D.");
  }
  const int n = static_cast<int>(x_info.shape[0]);
  const int p = static_cast<int>(x_info.shape[1]);
  if (static_cast<int>(y_info.shape[0]) != n) {
    throw std::invalid_argument("x_scaled and y_scaled have incompatible row counts.");
  }

  const auto* x_scaled = static_cast<const double*>(x_info.ptr);
  const auto* y_scaled = static_cast<const double*>(y_info.ptr);
  std::vector<int> metric = copy_int_array(metric_arr);
  std::vector<int> members = copy_int_array(member_arr);
  std::vector<double> proposal_sd = copy_double_array(proposal_sd_arr, 1);
  std::vector<double> proposal_mu = copy_double_array(proposal_mu_arr, 1);
  std::vector<int> binary_cols = copy_int_array(binary_cols_arr);
  const ReducedMetric reduced = make_reduced_metric(metric, members);

  if (static_cast<int>(metric.size()) != p || static_cast<int>(members.size()) != p ||
      static_cast<int>(proposal_sd.size()) != p || static_cast<int>(proposal_mu.size()) != p) {
    throw std::invalid_argument("Per-feature arrays must match x_scaled column count.");
  }
  if (static_cast<int>(init_tess.size()) != m || static_cast<int>(init_dim.size()) != m ||
      static_cast<int>(init_pred.size()) != m) {
    throw std::invalid_argument("Initial state lists must have length n_tessellations.");
  }

  std::vector<std::vector<double>> tess(m);
  std::vector<int> tess_n_centres(m);
  std::vector<int> tess_dim_count(m);
  std::vector<std::vector<int>> dim(m);
  std::vector<std::vector<double>> pred(m);

  for (int idx = 0; idx < m; ++idx) {
    std::vector<ssize_t> tess_shape;
    tess[idx] = copy_double_array(init_tess[idx], 2, &tess_shape);
    tess_n_centres[idx] = static_cast<int>(tess_shape[0]);
    tess_dim_count[idx] = static_cast<int>(tess_shape[1]);
    dim[idx] = copy_int_array(init_dim[idx]);
    pred[idx] = copy_double_array(init_pred[idx], 1);
    if (static_cast<int>(dim[idx].size()) != tess_dim_count[idx] ||
        static_cast<int>(pred[idx].size()) != tess_n_centres[idx]) {
      throw std::invalid_argument("Initial tessellation, dimension, and prediction shapes are incompatible.");
    }
  }

  std::mt19937_64 rng(seed);
  std::vector<std::vector<int>> current_indices(m);
  std::vector<double> sum_all_tess(n, 0.0);
  for (int idx = 0; idx < m; ++idx) {
    current_indices[idx] = knn1_internal(
        x_scaled,
        n,
        p,
        tess[idx],
        tess_n_centres[idx],
        tess_dim_count[idx],
        dim[idx],
        reduced.metric,
        reduced.member_counts);
    for (int obs = 0; obs < n; ++obs) {
      sum_all_tess[obs] += pred[idx][current_indices[idx][obs]];
    }
  }

  const int num_samples = total_iter > burn_in ? (total_iter - burn_in) / thinning : 0;
  py::list posterior_tess;
  py::list posterior_dim;
  py::list posterior_pred;
  std::vector<double> posterior_sigma(num_samples);
  std::vector<double> prediction_matrix(n * num_samples, 0.0);

  double sigma_squared = 1.0;
  std::vector<double> last_tess_pred(n, 0.0);
  int storage_idx = 0;
  const int progress_step = std::max(1, total_iter / 10);

  for (int iter = 1; iter <= total_iter; ++iter) {
    if (verbose && (iter % progress_step == 0 || iter == total_iter)) {
      py::print("  Iteration", iter, "/", total_iter);
    }

    double sum_sq = 0.0;
    for (int obs = 0; obs < n; ++obs) {
      const double residual = y_scaled[obs] - sum_all_tess[obs];
      sum_sq += residual * residual;
    }
    const double shape = (nu + n) / 2.0;
    const double rate = (nu * lambda + sum_sq) / 2.0;
    std::gamma_distribution<double> gamma(shape, 1.0 / rate);
    sigma_squared = 1.0 / gamma(rng);

    for (int j = 0; j < m; ++j) {
      if (j == 0) {
        for (int obs = 0; obs < n; ++obs) {
          sum_all_tess[obs] -= pred[j][current_indices[j][obs]];
        }
      } else {
        for (int obs = 0; obs < n; ++obs) {
          sum_all_tess[obs] += last_tess_pred[obs] - pred[j][current_indices[j][obs]];
        }
      }

      std::vector<double> residuals(n);
      for (int obs = 0; obs < n; ++obs) {
        residuals[obs] = y_scaled[obs] - sum_all_tess[obs];
      }

      ProposalResult proposal = propose(
          tess[j],
          tess_n_centres[j],
          tess_dim_count[j],
          dim[j],
          p,
          proposal_sd,
          proposal_mu,
          metric,
          members,
          rng);

      if (!binary_cols.empty()) {
        for (int local_dim = 0; local_dim < static_cast<int>(proposal.dim.size()); ++local_dim) {
          if (!in_vector(proposal.dim[local_dim], binary_cols)) {
            continue;
          }
          const int d_new = static_cast<int>(proposal.dim.size());
          for (int row = 0; row < proposal.n_centres; ++row) {
            double& value = proposal.tess[row * d_new + local_dim];
            value = std::clamp(value, 0.0, cat_scaling);
          }
        }
      }

      std::vector<int> proposed_indices = knn1_internal(
          x_scaled,
          n,
          p,
          proposal.tess,
          proposal.n_centres,
          static_cast<int>(proposal.dim.size()),
          proposal.dim,
          reduced.metric,
          reduced.member_counts);

      std::vector<double> r_old;
      std::vector<double> r_new;
      std::vector<int> n_old;
      std::vector<int> n_new;
      aggregate_residuals(
          residuals,
          current_indices[j],
          tess_n_centres[j],
          proposed_indices,
          proposal.n_centres,
          r_old,
          n_old,
          r_new,
          n_new);

      const bool has_empty = std::any_of(n_new.begin(), n_new.end(), [](int value) { return value == 0; });
      bool accepted = false;
      if (!has_empty) {
        const double log_alpha = log_acceptance_probability(
            r_old,
            n_old,
            r_new,
            n_new,
            static_cast<int>(proposal.dim.size()),
            proposal.n_centres,
            sigma_squared,
            sigma_squared_mu,
            omega,
            lambda_rate,
            p,
            proposal.modification);
        accepted = std::log(uniform01(rng)) < log_alpha;
      }

      if (accepted) {
        tess[j] = std::move(proposal.tess);
        tess_n_centres[j] = proposal.n_centres;
        tess_dim_count[j] = static_cast<int>(proposal.dim.size());
        dim[j] = std::move(proposal.dim);
        current_indices[j] = std::move(proposed_indices);
        pred[j] = sample_mu(r_new, n_new, sigma_squared_mu, sigma_squared, rng);
        for (int obs = 0; obs < n; ++obs) {
          last_tess_pred[obs] = pred[j][current_indices[j][obs]];
        }
      } else {
        pred[j] = sample_mu(r_old, n_old, sigma_squared_mu, sigma_squared, rng);
        for (int obs = 0; obs < n; ++obs) {
          last_tess_pred[obs] = pred[j][current_indices[j][obs]];
        }
      }

      if (j == m - 1) {
        for (int obs = 0; obs < n; ++obs) {
          sum_all_tess[obs] += last_tess_pred[obs];
        }
      }
    }

    if (iter > burn_in && (iter - burn_in) % thinning == 0) {
      for (int obs = 0; obs < n; ++obs) {
        prediction_matrix[obs * num_samples + storage_idx] = sum_all_tess[obs];
      }
      posterior_sigma[storage_idx] = sigma_squared;

      py::list sample_tess;
      py::list sample_dim;
      py::list sample_pred;
      for (int j = 0; j < m; ++j) {
        sample_tess.append(make_matrix(tess[j], tess_n_centres[j], tess_dim_count[j]));
        sample_dim.append(make_int_vector(dim[j]));
        sample_pred.append(make_vector(pred[j]));
      }
      posterior_tess.append(sample_tess);
      posterior_dim.append(sample_dim);
      posterior_pred.append(sample_pred);
      ++storage_idx;
    }
  }

  py::array_t<double> pred_matrix_arr({n, num_samples});
  auto pred_matrix = pred_matrix_arr.mutable_unchecked<2>();
  for (int obs = 0; obs < n; ++obs) {
    for (int sample = 0; sample < num_samples; ++sample) {
      pred_matrix(obs, sample) = prediction_matrix[obs * num_samples + sample];
    }
  }

  py::dict result;
  result["posterior_tess"] = posterior_tess;
  result["posterior_dim"] = posterior_dim;
  result["posterior_pred"] = posterior_pred;
  result["posterior_sigma"] = make_vector(posterior_sigma);
  result["prediction_matrix"] = pred_matrix_arr;
  return result;
}

py::array_t<int> cell_indices(py::array_t<double, py::array::c_style | py::array::forcecast> query_arr,
                              py::array_t<double, py::array::c_style | py::array::forcecast> centres_arr,
                              py::array_t<int, py::array::c_style | py::array::forcecast> dim_arr,
                              py::array_t<int, py::array::c_style | py::array::forcecast> metric_red_arr,
                              py::array_t<int, py::array::c_style | py::array::forcecast> member_red_arr) {
  py::buffer_info query_info = query_arr.request();
  py::buffer_info centres_info = centres_arr.request();
  if (query_info.ndim != 2 || centres_info.ndim != 2) {
    throw std::invalid_argument("query and centres must be two-dimensional arrays.");
  }

  const int n = static_cast<int>(query_info.shape[0]);
  const int p = static_cast<int>(query_info.shape[1]);
  const int n_centres = static_cast<int>(centres_info.shape[0]);
  const int d = static_cast<int>(centres_info.shape[1]);
  const auto* query = static_cast<const double*>(query_info.ptr);
  const auto* centres_ptr = static_cast<const double*>(centres_info.ptr);
  std::vector<double> centres(centres_ptr, centres_ptr + centres_info.size);
  std::vector<int> dim = copy_int_array(dim_arr);
  std::vector<int> metric_red = copy_int_array(metric_red_arr);
  std::vector<int> member_red = copy_int_array(member_red_arr);

  if (static_cast<int>(dim.size()) != d) {
    throw std::invalid_argument("dim length must match centres column count.");
  }
  const int member_total = std::accumulate(member_red.begin(), member_red.end(), 0);
  if (member_total != p) {
    throw std::invalid_argument("member_red must sum to query column count.");
  }

  const std::vector<int> indices = knn1_internal(query, n, p, centres, n_centres, d, dim, metric_red, member_red);
  return make_int_vector(indices);
}

PYBIND11_MODULE(_core, m) {
  m.doc() = "Standalone C++ backend for the AddiVortes Python package.";
  m.def("run_mcmc", &run_mcmc, py::arg("x_scaled"), py::arg("y_scaled"), py::arg("metric"),
        py::arg("members"), py::arg("n_tessellations"), py::arg("total_iter"), py::arg("burn_in"),
        py::arg("thinning"), py::arg("nu"), py::arg("lambda_value"), py::arg("sigma_squared_mu"),
        py::arg("omega"), py::arg("lambda_rate"), py::arg("proposal_sd"), py::arg("proposal_mu"),
        py::arg("init_tess"), py::arg("init_dim"), py::arg("init_pred"), py::arg("binary_cols"),
        py::arg("cat_scaling"), py::arg("seed"), py::arg("verbose"));
  m.def("cell_indices", &cell_indices, py::arg("query"), py::arg("centres"), py::arg("dim"),
        py::arg("metric_red"), py::arg("member_red"));
}
