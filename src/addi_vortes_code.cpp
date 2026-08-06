// [[Rcpp::plugins(cpp20)]

// Hot-path optimisations ported from lessons in the Rust addivortes engine:
// 1. Incremental cell reassignment with cached winning distance keys
// 2. Active-dimension-only Euclidean distance
// 3. Row-major (observation-major) packing of X for the NN loop
// 4. Specialised all-Euclidean assign path
// 5. Single-pass residual aggregation reused for MH and mu redraw
// 6. Preallocated scratch buffers outside the j/iter loops
// 8. Deferred posterior packaging (compact C++ store, R lists at end)
// 10/11. Same NN kernel + flattened posterior traversal for predict
// 13. Tightened helpers (masks, no per-proposal which_elem allocations)

#include <vector>
#include <string>
#include <algorithm>
#include <span>
#include <cmath>
#include <cstring>
#include <limits>

#define R_NO_REMAP

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Random.h>
#include <R_ext/Utils.h>

// ---------------------------------------------------------------------------
// Small helpers
// ---------------------------------------------------------------------------

static inline bool in_vector(int value, const std::vector<int>& vec) {
  return std::find(vec.begin(), vec.end(), value) != vec.end();
}

static inline int n_elem(int value, const std::vector<int>& vec) {
  int total = 0;
  for (int i = 0; i < static_cast<int>(vec.size()); ++i) {
    if (vec[i] == value) total++;
  }
  return total;
}

static inline double period_shift(double val, double lim) {
  while (val >= lim) val -= 2 * lim;
  while (val < -lim) val += 2 * lim;
  return val;
}

static inline void pack_row_major(const double* col_major, int n, int p,
                                  std::vector<double>& row_major) {
  row_major.resize(static_cast<size_t>(n) * static_cast<size_t>(p));
  for (int r = 0; r < n; ++r) {
    for (int c = 0; c < p; ++c) {
      row_major[static_cast<size_t>(r) * p + c] = col_major[r + c * n];
    }
  }
}

static inline bool all_euclidean_metric(const std::vector<int>& metric) {
  for (int m : metric) {
    if (m != 0) return false;
  }
  return true;
}

// ---------------------------------------------------------------------------
// Distance kernels (general / spherical / categorical)
// ---------------------------------------------------------------------------

static double euclidean_distance(std::span<const double> p1,
                                 std::span<const double> p2) {
  if (p1.size() != p2.size()) {
    Rf_error("Points have incompatible dimensions.");
  }
  double dist = 0.0;
  for (size_t i = 0; i < p1.size(); ++i) {
    const double diff = p1[i] - p2[i];
    dist += diff * diff;
  }
  return dist;
}

static double spherical_distance(std::span<const double> p1,
                                 std::span<const double> p2) {
  if (p1.size() != p2.size()) {
    Rf_error("Points have incompatible dimensions.");
  }
  if (p1.size() == 1) {
    double a1 = std::fabs(p1[0] - p2[0]);
    double a2 = 2 * M_PI - a1;
    return (a1 < a2) ? (a1 * a1) : (a2 * a2);
  }
  double angle_diff = cos(p1[p1.size() - 1] - p2[p2.size() - 1]);
  for (int i = static_cast<int>(p1.size()) - 2; i >= 0; --i) {
    double internal = sin(p1[i]) * sin(p2[i]) +
                      cos(p1[i]) * cos(p2[i]) * angle_diff;
    if (internal > 1) internal = 1;
    if (internal < -1) internal = -1.0;
    if (i == 0) angle_diff = acos(internal);
    else angle_diff = internal;
  }
  return angle_diff * angle_diff;
}

static double categorical_distance(std::span<const double> p1,
                                   std::span<const double> p2,
                                   const std::vector<int>& ncat) {
  if (p1.size() != p2.size()) {
    Rf_error("Points have incompatible dimensions.");
  }
  if (p1.size() != ncat.size()) {
    Rf_error("Point dimension does not match stated number of categorical variables.");
  }
  double dist = 0;
  for (size_t i = 0; i < p1.size(); ++i) {
    if (floor(p1[i]) != p1[i]) Rf_error("Not all coordinates in p1 are integer.");
    if (floor(p2[i]) != p2[i]) Rf_error("Not all coordinates in p2 are integer.");
    if (p1[i] != p2[i]) dist += 2.0 / (ncat[i] * ncat[i]);
  }
  return dist;
}

// Distance over full-width synthesised points (inactive dims equal => 0 contrib).
static double calc_distance(const double* vec1, const double* vec2, int p,
                            const std::vector<int>& nvals,
                            const std::vector<int>& type,
                            const std::vector<int>& cats) {
  int idx = 0;
  double tot = 0;
  int cat_idx = 0;
  for (size_t i = 0; i < nvals.size(); ++i) {
    int these_vals = nvals[i];
    std::span<const double> subvec1(vec1 + idx, these_vals);
    std::span<const double> subvec2(vec2 + idx, these_vals);
    if (type[i] == 0) {
      tot += euclidean_distance(subvec1, subvec2);
    } else if (type[i] == 1) {
      tot += spherical_distance(subvec1, subvec2);
    } else if (type[i] == 2) {
      std::vector<int> these_cats(cats.begin() + cat_idx,
                                  cats.begin() + cat_idx + these_vals);
      tot += categorical_distance(subvec1, subvec2, these_cats);
      cat_idx += these_vals;
    }
    idx += these_vals;
  }
  (void)p;
  return tot;
}

// ---------------------------------------------------------------------------
// Assignment cache + incremental reassignment
// ---------------------------------------------------------------------------

enum class AssignmentDelta {
  CentreAdded,
  CentreRemoved,
  CentreMoved,
  FullRecompute
};

struct AssignmentCache {
  std::vector<int> assignment;   // 0-based centre index per observation
  std::vector<double> best_keys; // winning comparison key per observation
};

static inline double euclidean_key_active(const double* active,
                                         const double* centres_cm,
                                         int nC, int c, int d) {
  double key = 0.0;
  for (int di = 0; di < d; ++di) {
    const double diff = active[di] - centres_cm[c + di * nC];
    key += diff * diff;
  }
  return key;
}

// Full n x nC Euclidean assign over active dimensions only.
static void assign_full_euclidean(const double* x_row, int n, int p,
                                  const double* centres, int nC, int d,
                                  const int* dim0,
                                  AssignmentCache& out,
                                  std::vector<double>& active_scratch) {
  out.assignment.resize(n);
  out.best_keys.resize(n);
  active_scratch.resize(d);
  for (int obs = 0; obs < n; ++obs) {
    const double* row = x_row + static_cast<size_t>(obs) * p;
    for (int di = 0; di < d; ++di) active_scratch[di] = row[dim0[di]];
    double best = std::numeric_limits<double>::infinity();
    int best_c = 0;
    for (int c = 0; c < nC; ++c) {
      const double key = euclidean_key_active(active_scratch.data(),
                                              centres, nC, c, d);
      if (key < best) {
        best = key;
        best_c = c;
      }
    }
    out.assignment[obs] = best_c;
    out.best_keys[obs] = best;
  }
}

static void reassign_added_euclidean(const double* x_row, int n, int p,
                                     const double* centres, int nC, int d,
                                     const int* dim0,
                                     const AssignmentCache& prev,
                                     AssignmentCache& out,
                                     std::vector<double>& active_scratch) {
  out.assignment = prev.assignment;
  out.best_keys = prev.best_keys;
  active_scratch.resize(d);
  const int added = nC - 1;
  for (int obs = 0; obs < n; ++obs) {
    const double* row = x_row + static_cast<size_t>(obs) * p;
    for (int di = 0; di < d; ++di) active_scratch[di] = row[dim0[di]];
    const double key = euclidean_key_active(active_scratch.data(),
                                            centres, nC, added, d);
    if (key < out.best_keys[obs]) {
      out.best_keys[obs] = key;
      out.assignment[obs] = added;
    }
  }
}

static void reassign_moved_euclidean(const double* x_row, int n, int p,
                                     const double* centres, int nC, int d,
                                     const int* dim0,
                                     int moved,
                                     const AssignmentCache& prev,
                                     AssignmentCache& out,
                                     std::vector<double>& active_scratch) {
  out.assignment = prev.assignment;
  out.best_keys = prev.best_keys;
  active_scratch.resize(d);
  for (int obs = 0; obs < n; ++obs) {
    const double* row = x_row + static_cast<size_t>(obs) * p;
    for (int di = 0; di < d; ++di) active_scratch[di] = row[dim0[di]];
    const int incumbent = prev.assignment[obs];
    if (incumbent == moved) {
      double best = std::numeric_limits<double>::infinity();
      int best_c = 0;
      for (int c = 0; c < nC; ++c) {
        const double key = euclidean_key_active(active_scratch.data(),
                                                centres, nC, c, d);
        if (key < best) {
          best = key;
          best_c = c;
        }
      }
      out.assignment[obs] = best_c;
      out.best_keys[obs] = best;
    } else {
      const double key = euclidean_key_active(active_scratch.data(),
                                              centres, nC, moved, d);
      if (key < out.best_keys[obs] ||
          (key == out.best_keys[obs] && moved < incumbent)) {
        out.best_keys[obs] = key;
        out.assignment[obs] = moved;
      }
    }
  }
}

static void reassign_removed_euclidean(const double* x_row, int n, int p,
                                       const double* centres, int nC, int d,
                                       const int* dim0,
                                       int removed,
                                       const AssignmentCache& prev,
                                       AssignmentCache& out,
                                       std::vector<double>& active_scratch) {
  out.assignment.resize(n);
  out.best_keys.resize(n);
  active_scratch.resize(d);
  for (int obs = 0; obs < n; ++obs) {
    const int incumbent = prev.assignment[obs];
    if (incumbent == removed) {
      const double* row = x_row + static_cast<size_t>(obs) * p;
      for (int di = 0; di < d; ++di) active_scratch[di] = row[dim0[di]];
      double best = std::numeric_limits<double>::infinity();
      int best_c = 0;
      for (int c = 0; c < nC; ++c) {
        const double key = euclidean_key_active(active_scratch.data(),
                                                centres, nC, c, d);
        if (key < best) {
          best = key;
          best_c = c;
        }
      }
      out.assignment[obs] = best_c;
      out.best_keys[obs] = best;
    } else {
      out.assignment[obs] = (incumbent > removed) ? incumbent - 1 : incumbent;
      out.best_keys[obs] = prev.best_keys[obs];
    }
  }
}

// General (non-all-Euclidean) full assign via synthesised centre rows.
static void assign_full_general(const double* x_row, int n, int p,
                                const double* centres, int nC, int d,
                                const int* dim0,
                                const std::vector<int>& metric,
                                const std::vector<int>& members,
                                const std::vector<int>& cats,
                                AssignmentCache& out,
                                std::vector<double>& synth) {
  out.assignment.resize(n);
  out.best_keys.resize(n);
  synth.resize(p);
  for (int obs = 0; obs < n; ++obs) {
    const double* row = x_row + static_cast<size_t>(obs) * p;
    std::memcpy(synth.data(), row, static_cast<size_t>(p) * sizeof(double));
    double best = std::numeric_limits<double>::infinity();
    int best_c = 0;
    for (int c = 0; c < nC; ++c) {
      for (int di = 0; di < d; ++di) synth[dim0[di]] = centres[c + di * nC];
      const double key = calc_distance(row, synth.data(), p, members, metric, cats);
      if (key < best) {
        best = key;
        best_c = c;
      }
    }
    out.assignment[obs] = best_c;
    out.best_keys[obs] = best;
  }
}

static void reassign_added_general(const double* x_row, int n, int p,
                                   const double* centres, int nC, int d,
                                   const int* dim0,
                                   const std::vector<int>& metric,
                                   const std::vector<int>& members,
                                   const std::vector<int>& cats,
                                   const AssignmentCache& prev,
                                   AssignmentCache& out,
                                   std::vector<double>& synth) {
  out.assignment = prev.assignment;
  out.best_keys = prev.best_keys;
  synth.resize(p);
  const int added = nC - 1;
  for (int obs = 0; obs < n; ++obs) {
    const double* row = x_row + static_cast<size_t>(obs) * p;
    std::memcpy(synth.data(), row, static_cast<size_t>(p) * sizeof(double));
    for (int di = 0; di < d; ++di) synth[dim0[di]] = centres[added + di * nC];
    const double key = calc_distance(row, synth.data(), p, members, metric, cats);
    if (key < out.best_keys[obs]) {
      out.best_keys[obs] = key;
      out.assignment[obs] = added;
    }
  }
}

static void reassign_moved_general(const double* x_row, int n, int p,
                                   const double* centres, int nC, int d,
                                   const int* dim0,
                                   int moved,
                                   const std::vector<int>& metric,
                                   const std::vector<int>& members,
                                   const std::vector<int>& cats,
                                   const AssignmentCache& prev,
                                   AssignmentCache& out,
                                   std::vector<double>& synth) {
  out.assignment = prev.assignment;
  out.best_keys = prev.best_keys;
  synth.resize(p);
  for (int obs = 0; obs < n; ++obs) {
    const double* row = x_row + static_cast<size_t>(obs) * p;
    std::memcpy(synth.data(), row, static_cast<size_t>(p) * sizeof(double));
    const int incumbent = prev.assignment[obs];
    if (incumbent == moved) {
      double best = std::numeric_limits<double>::infinity();
      int best_c = 0;
      for (int c = 0; c < nC; ++c) {
        for (int di = 0; di < d; ++di) synth[dim0[di]] = centres[c + di * nC];
        const double key = calc_distance(row, synth.data(), p, members, metric, cats);
        if (key < best) {
          best = key;
          best_c = c;
        }
      }
      out.assignment[obs] = best_c;
      out.best_keys[obs] = best;
    } else {
      for (int di = 0; di < d; ++di) synth[dim0[di]] = centres[moved + di * nC];
      const double key = calc_distance(row, synth.data(), p, members, metric, cats);
      if (key < out.best_keys[obs] ||
          (key == out.best_keys[obs] && moved < incumbent)) {
        out.best_keys[obs] = key;
        out.assignment[obs] = moved;
      }
    }
  }
}

static void reassign_removed_general(const double* x_row, int n, int p,
                                     const double* centres, int nC, int d,
                                     const int* dim0,
                                     int removed,
                                     const std::vector<int>& metric,
                                     const std::vector<int>& members,
                                     const std::vector<int>& cats,
                                     const AssignmentCache& prev,
                                     AssignmentCache& out,
                                     std::vector<double>& synth) {
  out.assignment.resize(n);
  out.best_keys.resize(n);
  synth.resize(p);
  for (int obs = 0; obs < n; ++obs) {
    const int incumbent = prev.assignment[obs];
    if (incumbent == removed) {
      const double* row = x_row + static_cast<size_t>(obs) * p;
      std::memcpy(synth.data(), row, static_cast<size_t>(p) * sizeof(double));
      double best = std::numeric_limits<double>::infinity();
      int best_c = 0;
      for (int c = 0; c < nC; ++c) {
        for (int di = 0; di < d; ++di) synth[dim0[di]] = centres[c + di * nC];
        const double key = calc_distance(row, synth.data(), p, members, metric, cats);
        if (key < best) {
          best = key;
          best_c = c;
        }
      }
      out.assignment[obs] = best_c;
      out.best_keys[obs] = best;
    } else {
      out.assignment[obs] = (incumbent > removed) ? incumbent - 1 : incumbent;
      out.best_keys[obs] = prev.best_keys[obs];
    }
  }
}

struct AssignScratch {
  std::vector<double> active;
  std::vector<double> synth;
  std::vector<int> dim0;
};

static void reassign(const double* x_row, int n, int p,
                     const double* centres, int nC, int d,
                     const std::vector<int>& dim1, // 1-based
                     AssignmentDelta delta, int touched,
                     bool euclidean,
                     const std::vector<int>& metric,
                     const std::vector<int>& members,
                     const std::vector<int>& cats,
                     const AssignmentCache& prev,
                     AssignmentCache& out,
                     AssignScratch& scratch) {
  scratch.dim0.resize(d);
  for (int i = 0; i < d; ++i) scratch.dim0[i] = dim1[i] - 1;

  const bool cold = prev.assignment.size() != static_cast<size_t>(n) ||
                    prev.best_keys.size() != static_cast<size_t>(n);
  if (cold || delta == AssignmentDelta::FullRecompute) {
    if (euclidean) {
      assign_full_euclidean(x_row, n, p, centres, nC, d, scratch.dim0.data(),
                            out, scratch.active);
    } else {
      assign_full_general(x_row, n, p, centres, nC, d, scratch.dim0.data(),
                          metric, members, cats, out, scratch.synth);
    }
    return;
  }

  if (euclidean) {
    switch (delta) {
      case AssignmentDelta::CentreAdded:
        reassign_added_euclidean(x_row, n, p, centres, nC, d, scratch.dim0.data(),
                                 prev, out, scratch.active);
        break;
      case AssignmentDelta::CentreMoved:
        reassign_moved_euclidean(x_row, n, p, centres, nC, d, scratch.dim0.data(),
                                 touched, prev, out, scratch.active);
        break;
      case AssignmentDelta::CentreRemoved:
        reassign_removed_euclidean(x_row, n, p, centres, nC, d, scratch.dim0.data(),
                                   touched, prev, out, scratch.active);
        break;
      case AssignmentDelta::FullRecompute:
        assign_full_euclidean(x_row, n, p, centres, nC, d, scratch.dim0.data(),
                              out, scratch.active);
        break;
    }
  } else {
    switch (delta) {
      case AssignmentDelta::CentreAdded:
        reassign_added_general(x_row, n, p, centres, nC, d, scratch.dim0.data(),
                               metric, members, cats, prev, out, scratch.synth);
        break;
      case AssignmentDelta::CentreMoved:
        reassign_moved_general(x_row, n, p, centres, nC, d, scratch.dim0.data(),
                               touched, metric, members, cats, prev, out,
                               scratch.synth);
        break;
      case AssignmentDelta::CentreRemoved:
        reassign_removed_general(x_row, n, p, centres, nC, d, scratch.dim0.data(),
                                 touched, metric, members, cats, prev, out,
                                 scratch.synth);
        break;
      case AssignmentDelta::FullRecompute:
        assign_full_general(x_row, n, p, centres, nC, d, scratch.dim0.data(),
                            metric, members, cats, out, scratch.synth);
        break;
    }
  }
}

// Aggregate residuals into per-cell sums/counts (one pass). Optionally only
// the "new" side, or both old and new.
static void aggregate_residuals_both(
    const std::vector<double>& R_j,
    const std::vector<int>& idx_old, int nC_old,
    const std::vector<int>& idx_new, int nC_new,
    std::vector<double>& R_old, std::vector<int>& n_old,
    std::vector<double>& R_new, std::vector<int>& n_new) {
  R_old.assign(nC_old, 0.0);
  n_old.assign(nC_old, 0);
  R_new.assign(nC_new, 0.0);
  n_new.assign(nC_new, 0);
  const int n = static_cast<int>(R_j.size());
  for (int obs = 0; obs < n; ++obs) {
    R_old[idx_old[obs]] += R_j[obs];
    n_old[idx_old[obs]]++;
    R_new[idx_new[obs]] += R_j[obs];
    n_new[idx_new[obs]]++;
  }
}

// ---------------------------------------------------------------------------
// Progress helpers
// ---------------------------------------------------------------------------

static void progress_bar(const char* label, int current, int total,
                         int width = 40) {
  if (total <= 0) return;
  if (current < 0) current = 0;
  if (current > total) current = total;
  const double frac = static_cast<double>(current) / static_cast<double>(total);
  int filled = static_cast<int>(frac * width + 1e-12);
  if (filled > width) filled = width;
  Rprintf("\r%s [", label);
  for (int i = 0; i < width; ++i) Rprintf("%c", i < filled ? '=' : ' ');
  Rprintf("] %3.0f%%  %d/%d", 100.0 * frac, current, total);
  if (current >= total) Rprintf("\n");
  R_FlushConsole();
}

static void maybe_progress(const char* label, int current, int total,
                           int width, int& last_filled, bool showProgress) {
  if (showProgress) {
    int filled = 0;
    if (total > 0) {
      filled = static_cast<int>(
        static_cast<double>(current) / static_cast<double>(total) * width + 1e-12);
      if (filled > width) filled = width;
    }
    if (current == 1 || current == total || filled != last_filled) {
      progress_bar(label, current, total, width);
      last_filled = filled;
      R_CheckUserInterrupt();
    }
  } else if (current == 1 || current % 64 == 0 || current == total) {
    R_CheckUserInterrupt();
  }
}

// ---------------------------------------------------------------------------
// MH / mu helpers
// ---------------------------------------------------------------------------

static double tessellation_log_likelihood_component(
    const std::vector<double>& R, const std::vector<int>& n_counts,
    double sigmaSquared, double sigmaSquaredMu) {
  double sum_log = 0.0, sum_R = 0.0;
  for (int k = 0; k < static_cast<int>(n_counts.size()); ++k) {
    double den = n_counts[k] * sigmaSquaredMu + sigmaSquared;
    sum_log += log(den);
    sum_R += (R[k] * R[k]) / den;
  }
  return -0.5 * sum_log + (sigmaSquaredMu / (2.0 * sigmaSquared)) * sum_R;
}

struct AcceptanceComponents {
  double logAlpha;
};

static AcceptanceComponents log_acceptance_components(
    const std::vector<double>& R_old, const std::vector<int>& n_old,
    const std::vector<double>& R_new, const std::vector<int>& n_new,
    int d_new, int nC_new,
    double sigmaSquared, double sigmaSquaredMu,
    double omega, double lambdaRate, int p,
    const std::string& mod) {
  double old_log_lik = tessellation_log_likelihood_component(
    R_old, n_old, sigmaSquared, sigmaSquaredMu);
  double new_log_lik = tessellation_log_likelihood_component(
    R_new, n_new, sigmaSquared, sigmaSquaredMu);
  double acc = new_log_lik - old_log_lik;

  if (mod == "AD") {
    double log_ts_tr = log(static_cast<double>(p - d_new + 1))
                     - log(static_cast<double>(d_new - 1))
                     + log(omega)
                     - log(p - omega);
    acc += log_ts_tr;
    if (d_new == 2) acc += -log(2);
    if (d_new == p) acc += log(2);
  } else if (mod == "RD") {
    double log_ts_tr = log(static_cast<double>(d_new))
                     - log(static_cast<double>(p - d_new))
                     + log(p - omega)
                     - log(omega);
    acc += log_ts_tr;
    if (d_new == (p - 1)) acc += -log(2);
    if (d_new == 1) acc += log(2);
  } else if (mod == "AC") {
    double log_ts_tr = log(lambdaRate) - log(static_cast<double>(nC_new - 1));
    acc += log_ts_tr + 0.5 * log(sigmaSquared);
    if (nC_new == 2) acc += -log(2);
  } else if (mod == "RC") {
    double log_ts_tr = log(static_cast<double>(nC_new)) - log(lambdaRate);
    acc += log_ts_tr - 0.5 * log(sigmaSquared);
    if (nC_new == 1) acc += log(2);
  }
  return {acc};
}

static void sample_mu_into(const std::vector<double>& R_ij,
                           const std::vector<int>& n_ij,
                           double sigmaSquaredMu, double sigmaSquared,
                           std::vector<double>& result) {
  const int N = static_cast<int>(R_ij.size());
  result.resize(N);
  for (int k = 0; k < N; ++k) {
    double den = sigmaSquaredMu * n_ij[k] + sigmaSquared;
    double mean = (sigmaSquaredMu * R_ij[k]) / den;
    double sd = sqrt((sigmaSquared * sigmaSquaredMu) / den);
    result[k] = mean + norm_rand() * sd;
  }
}

// ---------------------------------------------------------------------------
// Tessellation proposals (with AssignmentDelta)
// ---------------------------------------------------------------------------

struct ProposalResult {
  std::vector<double> tess; // column-major nC x d
  int nC;
  std::vector<int> dim; // 1-based
  std::string mod;
  AssignmentDelta delta;
  int touched; // moved/removed centre index (0-based); unused otherwise
};

// Precomputed categorical column -> ncats index (-1 if not categorical).
static ProposalResult propose_internal(
    const std::vector<double>& tess_j, int nC, int d_j,
    const std::vector<int>& dim_j,
    int p,
    const double* sd, const double* mus,
    const std::vector<int>& metric,
    const std::vector<int>& members,
    const std::vector<int>& cats,
    const std::vector<int>& cat_index_of_col) {
  ProposalResult r;
  r.tess = tess_j;
  r.nC = nC;
  r.dim = dim_j;
  r.mod = "Change";
  r.delta = AssignmentDelta::CentreMoved;
  r.touched = 0;

  const double prand = unif_rand();
  double new_val;

  auto draw_coord = [&](int global0) {
    double v = mus[global0] + norm_rand() * sd[global0];
    if (metric[global0] == 1) {
      if (global0 == static_cast<int>(members.size()) - 1 ||
          members[global0 + 1] != members[global0]) {
        v = period_shift(v, M_PI);
      }
    } else if (metric[global0] == 2) {
      const int ci = cat_index_of_col[global0];
      v = 1.0 + floor(unif_rand() * cats[ci]);
    }
    return v;
  };

  if ((prand < 0.2 && d_j != p) || (d_j == 1 && d_j != p && prand < 0.4)) {
    r.mod = "AD";
    r.delta = AssignmentDelta::FullRecompute;
    int new_dim;
    do { new_dim = static_cast<int>(unif_rand() * p) + 1; }
    while (in_vector(new_dim, r.dim));
    r.dim.push_back(new_dim);

    std::vector<double> new_tess(static_cast<size_t>(nC) * (d_j + 1));
    for (int row = 0; row < nC; ++row) {
      for (int col = 0; col < d_j; ++col)
        new_tess[row + col * nC] = tess_j[row + col * nC];
      new_tess[row + d_j * nC] = draw_coord(new_dim - 1);
    }
    r.tess = std::move(new_tess);
    r.nC = nC;

  } else if (prand < 0.4 && d_j > 1) {
    r.mod = "RD";
    r.delta = AssignmentDelta::FullRecompute;
    int rm_idx = static_cast<int>(unif_rand() * d_j);
    r.dim.erase(r.dim.begin() + rm_idx);
    std::vector<double> new_tess(static_cast<size_t>(nC) * (d_j - 1));
    int cur_col = 0;
    for (int col = 0; col < d_j; ++col) {
      if (col != rm_idx) {
        for (int row = 0; row < nC; ++row)
          new_tess[row + cur_col * nC] = tess_j[row + col * nC];
        cur_col++;
      }
    }
    r.tess = std::move(new_tess);
    r.nC = nC;

  } else if (prand < 0.6 || (prand < 0.8 && nC == 1)) {
    r.mod = "AC";
    r.delta = AssignmentDelta::CentreAdded;
    r.tess = tess_j;
    for (int i = 0; i < d_j; ++i) {
      new_val = draw_coord(dim_j[i] - 1);
      r.tess.insert(r.tess.begin() + (i * (nC + 1)) + nC, new_val);
    }
    r.nC = nC + 1;

  } else if (prand < 0.8 && nC > 1) {
    r.mod = "RC";
    r.delta = AssignmentDelta::CentreRemoved;
    int rm_row = static_cast<int>(unif_rand() * nC);
    r.touched = rm_row;
    std::vector<double> new_tess;
    new_tess.reserve(static_cast<size_t>(nC - 1) * d_j);
    for (int col = 0; col < d_j; ++col)
      for (int row = 0; row < nC; ++row)
        if (row != rm_row) new_tess.push_back(tess_j[row + col * nC]);
    r.tess = std::move(new_tess);
    r.nC = nC - 1;

  } else if (prand < 0.9 || d_j == p) {
    r.mod = "Change";
    r.delta = AssignmentDelta::CentreMoved;
    int ci = static_cast<int>(unif_rand() * nC);
    r.touched = ci;
    for (int col = 0; col < d_j; ++col) {
      r.tess[ci + col * nC] = draw_coord(dim_j[col] - 1);
    }

  } else {
    r.mod = "Swap";
    r.delta = AssignmentDelta::FullRecompute;
    int swap_idx = static_cast<int>(unif_rand() * d_j);
    int new_dim;
    do { new_dim = static_cast<int>(unif_rand() * p) + 1; }
    while (in_vector(new_dim, r.dim));
    r.dim[swap_idx] = new_dim;
    for (int row = 0; row < nC; ++row) {
      r.tess[row + swap_idx * nC] = draw_coord(new_dim - 1);
    }
  }

  return r;
}

// ---------------------------------------------------------------------------
// Compact posterior store (deferred R packaging)
// ---------------------------------------------------------------------------

struct StoredTess {
  std::vector<double> centres; // column-major nC x d
  int nC = 0;
  int d = 0;
  std::vector<int> dim;        // 1-based
  std::vector<double> mu;
};

struct StoredDraw {
  std::vector<StoredTess> tessellations;
  double sigma = 0.0;
};

static SEXP pack_posterior_lists(const std::vector<StoredDraw>& draws, int m) {
  const int numSamples = static_cast<int>(draws.size());
  SEXP outTess = PROTECT(Rf_allocVector(VECSXP, numSamples));
  SEXP outDim = PROTECT(Rf_allocVector(VECSXP, numSamples));
  SEXP outPred = PROTECT(Rf_allocVector(VECSXP, numSamples));

  for (int s = 0; s < numSamples; ++s) {
    SEXP sampleTess = PROTECT(Rf_allocVector(VECSXP, m));
    SEXP sampleDim = PROTECT(Rf_allocVector(VECSXP, m));
    SEXP samplePred = PROTECT(Rf_allocVector(VECSXP, m));
    for (int j = 0; j < m; ++j) {
      const StoredTess& t = draws[s].tessellations[j];
      SEXP rTess = PROTECT(Rf_allocMatrix(REALSXP, t.nC, t.d));
      if (t.nC * t.d > 0) {
        std::memcpy(REAL(rTess), t.centres.data(),
                    static_cast<size_t>(t.nC) * t.d * sizeof(double));
      }
      SET_VECTOR_ELT(sampleTess, j, rTess);
      UNPROTECT(1);

      SEXP rDim = PROTECT(Rf_allocVector(INTSXP, t.d));
      if (t.d > 0) {
        std::memcpy(INTEGER(rDim), t.dim.data(),
                    static_cast<size_t>(t.d) * sizeof(int));
      }
      SET_VECTOR_ELT(sampleDim, j, rDim);
      UNPROTECT(1);

      SEXP rPred = PROTECT(Rf_allocVector(REALSXP, t.nC));
      if (t.nC > 0) {
        std::memcpy(REAL(rPred), t.mu.data(),
                    static_cast<size_t>(t.nC) * sizeof(double));
      }
      SET_VECTOR_ELT(samplePred, j, rPred);
      UNPROTECT(1);
    }
    SET_VECTOR_ELT(outTess, s, sampleTess);
    SET_VECTOR_ELT(outDim, s, sampleDim);
    SET_VECTOR_ELT(outPred, s, samplePred);
    UNPROTECT(3);
  }

  SEXP result = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(result, 0, outTess);
  SET_VECTOR_ELT(result, 1, outDim);
  SET_VECTOR_ELT(result, 2, outPred);
  UNPROTECT(4); // outTess, outDim, outPred, result (caller re-protects)
  return result;
}

// Flatten nested posterior lists once for predict.
struct FlatPosterior {
  int numSamples = 0;
  int m = 0;
  std::vector<double> centres; // concatenated column-major blocks
  std::vector<double> mus;
  std::vector<int> dims;       // 1-based, concatenated
  std::vector<int> centre_off; // index into centres for (s,j)
  std::vector<int> mu_off;
  std::vector<int> dim_off;
  std::vector<int> nC;
  std::vector<int> d;
};

static FlatPosterior flatten_posterior(SEXP posteriorTess_sexp,
                                       SEXP posteriorDim_sexp,
                                       SEXP posteriorPred_sexp) {
  FlatPosterior flat;
  flat.numSamples = Rf_length(posteriorTess_sexp);
  if (flat.numSamples == 0) return flat;
  flat.m = Rf_length(VECTOR_ELT(posteriorTess_sexp, 0));
  const int total = flat.numSamples * flat.m;
  flat.centre_off.resize(total);
  flat.mu_off.resize(total);
  flat.dim_off.resize(total);
  flat.nC.resize(total);
  flat.d.resize(total);

  int c_off = 0, m_off = 0, d_off = 0;
  for (int s = 0; s < flat.numSamples; ++s) {
    SEXP sampleTess = VECTOR_ELT(posteriorTess_sexp, s);
    SEXP sampleDim = VECTOR_ELT(posteriorDim_sexp, s);
    SEXP samplePred = VECTOR_ELT(posteriorPred_sexp, s);
    if (Rf_length(sampleTess) != flat.m ||
        Rf_length(sampleDim) != flat.m ||
        Rf_length(samplePred) != flat.m) {
      Rf_error("Posterior sample %d has inconsistent tessellation counts", s + 1);
    }
    for (int j = 0; j < flat.m; ++j) {
      const int idx = s * flat.m + j;
      SEXP tess_j = VECTOR_ELT(sampleTess, j);
      SEXP dim_j = VECTOR_ELT(sampleDim, j);
      SEXP pred_j = VECTOR_ELT(samplePred, j);
      const int nC = Rf_nrows(tess_j);
      const int d = Rf_ncols(tess_j);
      if (Rf_length(dim_j) != d) {
        Rf_error("Tessellation %d in sample %d has mismatched dim length",
                 j + 1, s + 1);
      }
      if (Rf_length(pred_j) != nC) {
        Rf_error("Tessellation %d in sample %d has mismatched pred length",
                 j + 1, s + 1);
      }
      flat.nC[idx] = nC;
      flat.d[idx] = d;
      flat.centre_off[idx] = c_off;
      flat.mu_off[idx] = m_off;
      flat.dim_off[idx] = d_off;

      flat.centres.resize(c_off + nC * d);
      if (nC * d > 0) {
        std::memcpy(flat.centres.data() + c_off, REAL(tess_j),
                    static_cast<size_t>(nC) * d * sizeof(double));
      }
      flat.mus.resize(m_off + nC);
      if (nC > 0) {
        std::memcpy(flat.mus.data() + m_off, REAL(pred_j),
                    static_cast<size_t>(nC) * sizeof(double));
      }
      flat.dims.resize(d_off + d);
      if (d > 0) {
        std::memcpy(flat.dims.data() + d_off, INTEGER(dim_j),
                    static_cast<size_t>(d) * sizeof(int));
      }
      c_off += nC * d;
      m_off += nC;
      d_off += d;
    }
  }
  return flat;
}

extern "C" {

  // ---------------------------------------------------------------------------
  // knnx_index_cpp — uses the same active-dim / row-major assign path
  // ---------------------------------------------------------------------------
  SEXP knnx_index_cpp(SEXP tess_sexp, SEXP query_sexp, SEXP dim_sexp,
                      SEXP dist_sexp, SEXP member_sexp) {
    double* p_tess = REAL(tess_sexp);
    double* p_query = REAL(query_sexp);
    int* dim_p = INTEGER(dim_sexp);
    int* member_ptr = INTEGER(member_sexp);
    int* metric_ptr = INTEGER(dist_sexp);

    std::vector<int> dim_p_temp(dim_p, dim_p + Rf_length(dim_sexp));
    const std::vector<int> metric(metric_ptr, metric_ptr + Rf_length(dist_sexp));
    const std::vector<int> members(member_ptr, member_ptr + Rf_length(member_sexp));

    int tess_rows = Rf_nrows(tess_sexp);
    int tess_cols = Rf_ncols(tess_sexp);
    int query_rows = Rf_nrows(query_sexp);
    int query_cols = Rf_ncols(query_sexp);

    int mem_sum = 0;
    for (int i = 0; i < static_cast<int>(members.size()); ++i) mem_sum += members[i];
    if (mem_sum != query_cols) {
      Rf_error("Length of metric must match number of columns in query/data matrices");
    }
    if (tess_cols != query_cols) {
      Rf_error("Dimensions of tess and query matrices must match");
    }
    if (tess_rows <= 0) {
      Rf_error("Reference set must contain at least one point");
    }

    std::vector<char> active_dim_mask(query_cols, 0);
    std::vector<int> active_dim_idx;
    active_dim_idx.reserve(dim_p_temp.size());
    for (int i = 0; i < static_cast<int>(dim_p_temp.size()); ++i) {
      const int d0 = dim_p_temp[i] - 1;
      if (d0 < 0 || d0 >= query_cols) {
        Rf_error("Values in dim must be valid 1-based column indices");
      }
      if (!active_dim_mask[d0]) {
        active_dim_mask[d0] = 1;
        active_dim_idx.push_back(d0);
      }
    }

    std::vector<int> metric_aug;
    for (int i = 0; i < static_cast<int>(metric.size()); ++i) {
      for (int j = 0; j < members[i]; ++j) metric_aug.push_back(metric[i]);
    }
    std::vector<int> ncats;
    if (in_vector(2, metric)) {
      for (int i = 0; i < static_cast<int>(metric_aug.size()); ++i) {
        if (metric_aug[i] == 2) {
          int max_val = 0;
          for (int j = 0; j < query_rows; ++j) {
            if (p_query[i * query_rows + j] > max_val)
              max_val = static_cast<int>(p_query[i * query_rows + j]);
          }
          ncats.push_back(max_val);
        }
      }
    }

    // Build reduced centres (active dims only, column-major) and dim1.
    const int d = static_cast<int>(active_dim_idx.size());
    std::vector<double> centres(static_cast<size_t>(tess_rows) * d);
    std::vector<int> dim1(d);
    for (int di = 0; di < d; ++di) {
      dim1[di] = active_dim_idx[di] + 1;
      for (int t = 0; t < tess_rows; ++t) {
        centres[t + di * tess_rows] = p_tess[t + active_dim_idx[di] * tess_rows];
      }
    }

    std::vector<double> x_row;
    pack_row_major(p_query, query_rows, query_cols, x_row);

    AssignmentCache cache;
    AssignScratch scratch;
    const bool eucl = all_euclidean_metric(metric);
    reassign(x_row.data(), query_rows, query_cols,
             centres.data(), tess_rows, d, dim1,
             AssignmentDelta::FullRecompute, 0, eucl,
             metric, members, ncats,
             AssignmentCache{}, cache, scratch);

    SEXP result = PROTECT(Rf_allocMatrix(INTSXP, query_rows, 1));
    int* p_result = INTEGER(result);
    for (int q = 0; q < query_rows; ++q) p_result[q] = cache.assignment[q] + 1;
    UNPROTECT(1);
    return result;
  }

  // ---------------------------------------------------------------------------
  // calculate_residuals_cpp (kept for registration compatibility)
  // ---------------------------------------------------------------------------
  SEXP calculate_residuals_cpp(SEXP R_j_sexp, SEXP indexes_sexp,
                               SEXP indexesStar_sexp, SEXP num_levels_old_sexp,
                               SEXP num_centres_new_sexp) {
    double* p_R_j = REAL(R_j_sexp);
    int* p_indexes = INTEGER(indexes_sexp);
    int* p_indexesStar = INTEGER(indexesStar_sexp);
    int num_levels_old = INTEGER(num_levels_old_sexp)[0];
    int num_centres_new = INTEGER(num_centres_new_sexp)[0];
    int n_obs = Rf_length(R_j_sexp);

    std::vector<double> R_ijOld(num_levels_old, 0.0);
    std::vector<int> n_ijOld(num_levels_old, 0);
    std::vector<double> R_ijNew(num_centres_new, 0.0);
    std::vector<int> n_ijNew(num_centres_new, 0);

    for (int i = 0; i < n_obs; ++i) {
      int idx_old = p_indexes[i] - 1;
      if (idx_old >= 0 && idx_old < num_levels_old) {
        R_ijOld[idx_old] += p_R_j[i];
        n_ijOld[idx_old]++;
      }
      int idx_new = p_indexesStar[i] - 1;
      if (idx_new >= 0 && idx_new < num_centres_new) {
        R_ijNew[idx_new] += p_R_j[i];
        n_ijNew[idx_new]++;
      }
    }

    SEXP res_R_ijOld, res_n_ijOld, res_R_ijNew, res_n_ijNew, result_list, list_names;
    PROTECT(res_R_ijOld = Rf_allocVector(REALSXP, num_levels_old));
    std::memcpy(REAL(res_R_ijOld), R_ijOld.data(), num_levels_old * sizeof(double));
    PROTECT(res_n_ijOld = Rf_allocVector(INTSXP, num_levels_old));
    std::memcpy(INTEGER(res_n_ijOld), n_ijOld.data(), num_levels_old * sizeof(int));
    PROTECT(res_R_ijNew = Rf_allocVector(REALSXP, num_centres_new));
    std::memcpy(REAL(res_R_ijNew), R_ijNew.data(), num_centres_new * sizeof(double));
    PROTECT(res_n_ijNew = Rf_allocVector(INTSXP, num_centres_new));
    std::memcpy(INTEGER(res_n_ijNew), n_ijNew.data(), num_centres_new * sizeof(int));
    PROTECT(result_list = Rf_allocVector(VECSXP, 4));
    SET_VECTOR_ELT(result_list, 0, res_R_ijOld);
    SET_VECTOR_ELT(result_list, 1, res_n_ijOld);
    SET_VECTOR_ELT(result_list, 2, res_R_ijNew);
    SET_VECTOR_ELT(result_list, 3, res_n_ijNew);
    PROTECT(list_names = Rf_allocVector(STRSXP, 4));
    SET_STRING_ELT(list_names, 0, Rf_mkChar("R_ijOld"));
    SET_STRING_ELT(list_names, 1, Rf_mkChar("n_ijOld"));
    SET_STRING_ELT(list_names, 2, Rf_mkChar("R_ijNew"));
    SET_STRING_ELT(list_names, 3, Rf_mkChar("n_ijNew"));
    Rf_setAttrib(result_list, R_NamesSymbol, list_names);
    UNPROTECT(6);
    return result_list;
  }

  // ---------------------------------------------------------------------------
  // addi_vortes_mcmc_cpp
  // ---------------------------------------------------------------------------
  SEXP addi_vortes_mcmc_cpp(
      SEXP xScaled_sexp,
      SEXP yScaled_sexp,
      SEXP metric_sexp,
      SEXP member_sexp,
      SEXP m_sexp,
      SEXP totalMCMCIter_sexp,
      SEXP mcmcBurnIn_sexp,
      SEXP thinning_sexp,
      SEXP nu_sexp,
      SEXP lambda_sexp,
      SEXP sigmaSquaredMu_sexp,
      SEXP omega_sexp,
      SEXP lambdaRate_sexp,
      SEXP sd_sexp,
      SEXP mus_sexp,
      SEXP init_tess_sexp,
      SEXP init_dim_sexp,
      SEXP init_pred_sexp,
      SEXP binaryCols_sexp,
      SEXP catScaling_sexp,
      SEXP showProgress_sexp) {

    const double* xScaled = REAL(xScaled_sexp);
    const double* yScaled = REAL(yScaled_sexp);
    int n = Rf_nrows(xScaled_sexp);
    int p = Rf_ncols(xScaled_sexp);
    int m = INTEGER(m_sexp)[0];
    int totalIter = INTEGER(totalMCMCIter_sexp)[0];
    int burnIn = INTEGER(mcmcBurnIn_sexp)[0];
    int thinning = INTEGER(thinning_sexp)[0];
    double nu = REAL(nu_sexp)[0];
    double lambda = REAL(lambda_sexp)[0];
    double sigSqMu = REAL(sigmaSquaredMu_sexp)[0];
    double omega = REAL(omega_sexp)[0];
    double lambdaRate = REAL(lambdaRate_sexp)[0];
    const double* sd = REAL(sd_sexp);
    const double* mus = REAL(mus_sexp);
    double catScaling = REAL(catScaling_sexp)[0];
    bool showProgress = LOGICAL(showProgress_sexp)[0];

    std::vector<int> metric(INTEGER(metric_sexp), INTEGER(metric_sexp) + p);
    std::vector<int> members(INTEGER(member_sexp), INTEGER(member_sexp) + p);

    // Binary column mask (option 13)
    std::vector<char> is_binary(p, 0);
    if (!Rf_isNull(binaryCols_sexp)) {
      int* bc = INTEGER(binaryCols_sexp);
      int nb = Rf_length(binaryCols_sexp);
      for (int i = 0; i < nb; ++i) {
        int g0 = bc[i] - 1;
        if (g0 >= 0 && g0 < p) is_binary[g0] = 1;
      }
    }

    // Reduced metric / membership for distance grouping
    std::vector<int> metric_red, member_red;
    {
      int i = 0;
      while (i < static_cast<int>(metric.size())) {
        int this_elem = members[i];
        int this_metric = metric[i];
        int how_many = n_elem(this_elem, members);
        member_red.push_back(how_many);
        metric_red.push_back(this_metric);
        i += how_many;
      }
    }

    std::vector<int> ncats;
    std::vector<int> cat_index_of_col(p, -1);
    if (in_vector(2, metric_red)) {
      int cat_i = 0;
      for (int i = 0; i < p; ++i) {
        if (metric[i] == 2) {
          int max_val = 0;
          for (int j = 0; j < n; ++j) {
            if (xScaled[i * n + j] > max_val)
              max_val = static_cast<int>(xScaled[i * n + j]);
          }
          ncats.push_back(max_val);
          cat_index_of_col[i] = cat_i++;
        }
      }
    }

    const bool euclidean = all_euclidean_metric(metric_red);

    // Row-major X (option 3)
    std::vector<double> x_row;
    pack_row_major(xScaled, n, p, x_row);

    // Initial state
    std::vector<std::vector<double>> tess(m);
    std::vector<int> tess_nC(m), tess_d(m);
    std::vector<std::vector<int>> dim_j(m);
    std::vector<std::vector<double>> pred(m);
    for (int j = 0; j < m; ++j) {
      SEXP t_j = VECTOR_ELT(init_tess_sexp, j);
      int rows = Rf_nrows(t_j), cols = Rf_ncols(t_j);
      tess_nC[j] = rows;
      tess_d[j] = cols;
      double* pt = REAL(t_j);
      tess[j].assign(pt, pt + rows * cols);
      SEXP d_j = VECTOR_ELT(init_dim_sexp, j);
      dim_j[j].assign(INTEGER(d_j), INTEGER(d_j) + Rf_length(d_j));
      SEXP p_j = VECTOR_ELT(init_pred_sexp, j);
      pred[j].assign(REAL(p_j), REAL(p_j) + Rf_length(p_j));
    }

    // Assignment caches (option 1)
    std::vector<AssignmentCache> caches(m);
    AssignScratch assign_scratch;
    AssignmentCache prop_cache;
    for (int j = 0; j < m; ++j) {
      reassign(x_row.data(), n, p,
               tess[j].data(), tess_nC[j], tess_d[j], dim_j[j],
               AssignmentDelta::FullRecompute, 0, euclidean,
               metric_red, member_red, ncats,
               AssignmentCache{}, caches[j], assign_scratch);
    }

    std::vector<double> sumAllTess(n, 0.0);
    for (int j = 0; j < m; ++j) {
      for (int obs = 0; obs < n; ++obs)
        sumAllTess[obs] += pred[j][caches[j].assignment[obs]];
    }

    int numSamples = 0;
    if (totalIter > burnIn) numSamples = (totalIter - burnIn) / thinning;
    if (numSamples < 0) numSamples = 0;

    // Deferred compact posterior store (option 8)
    std::vector<StoredDraw> stored;
    stored.reserve(numSamples);
    std::vector<double> predictionMatrix(static_cast<size_t>(n) * numSamples, 0.0);

    SEXP outTraceIteration = PROTECT(Rf_allocVector(INTSXP, totalIter));
    SEXP outTraceBurnIn = PROTECT(Rf_allocVector(LGLSXP, totalIter));
    SEXP outTraceAvgCenters = PROTECT(Rf_allocVector(REALSXP, totalIter));
    SEXP outTraceSdCenters = PROTECT(Rf_allocVector(REALSXP, totalIter));
    SEXP outTraceAvgDims = PROTECT(Rf_allocVector(REALSXP, totalIter));
    SEXP outTraceLogLik = PROTECT(Rf_allocVector(REALSXP, totalIter));

    // Scratch buffers (option 6)
    std::vector<double> R_j(n);
    std::vector<double> lastTessPred(n, 0.0);
    std::vector<double> R_old, R_new;
    std::vector<int> n_old, n_new;

    GetRNGstate();
    double sigmaSquared = 1.0;
    int storageIdx = 0;
    const int progressWidth = 40;
    int lastFilled = -1;

    for (int iter = 1; iter <= totalIter; ++iter) {
      maybe_progress("MCMC", iter, totalIter, progressWidth, lastFilled, showProgress);

      double sum_sq = 0.0;
      for (int obs = 0; obs < n; ++obs) {
        double r = yScaled[obs] - sumAllTess[obs];
        sum_sq += r * r;
      }
      double shape = (nu + n) / 2.0;
      double rate = (nu * lambda + sum_sq) / 2.0;
      sigmaSquared = 1.0 / rgamma(shape, 1.0 / rate);

      for (int j = 0; j < m; ++j) {
        if (j == 0) {
          for (int obs = 0; obs < n; ++obs)
            sumAllTess[obs] -= pred[j][caches[j].assignment[obs]];
        } else {
          for (int obs = 0; obs < n; ++obs)
            sumAllTess[obs] += lastTessPred[obs] - pred[j][caches[j].assignment[obs]];
        }

        for (int obs = 0; obs < n; ++obs)
          R_j[obs] = yScaled[obs] - sumAllTess[obs];

        ProposalResult prop = propose_internal(
          tess[j], tess_nC[j], tess_d[j], dim_j[j],
          p, sd, mus, metric, members, ncats, cat_index_of_col);

        // Clamp binary columns via mask (option 13)
        if (!Rf_isNull(binaryCols_sexp)) {
          int d_star = static_cast<int>(prop.dim.size());
          for (int di = 0; di < d_star; ++di) {
            int g0 = prop.dim[di] - 1;
            if (g0 >= 0 && g0 < p && is_binary[g0]) {
              for (int row = 0; row < prop.nC; ++row) {
                double v = prop.tess[row + di * prop.nC];
                if (v < 0.0) v = 0.0;
                if (v > catScaling) v = catScaling;
                prop.tess[row + di * prop.nC] = v;
              }
            }
          }
        }

        reassign(x_row.data(), n, p,
                 prop.tess.data(), prop.nC, static_cast<int>(prop.dim.size()),
                 prop.dim, prop.delta, prop.touched, euclidean,
                 metric_red, member_red, ncats,
                 caches[j], prop_cache, assign_scratch);

        // Option 5: one pass builds both aggregates; reused for MH and mu draw
        aggregate_residuals_both(R_j,
          caches[j].assignment, tess_nC[j],
          prop_cache.assignment, prop.nC,
          R_old, n_old, R_new, n_new);

        bool hasEmpty = false;
        for (int k = 0; k < prop.nC; ++k) {
          if (n_new[k] == 0) { hasEmpty = true; break; }
        }

        bool accepted = false;
        if (!hasEmpty) {
          AcceptanceComponents acc = log_acceptance_components(
            R_old, n_old, R_new, n_new,
            static_cast<int>(prop.dim.size()), prop.nC,
            sigmaSquared, sigSqMu, omega, lambdaRate, p, prop.mod);
          accepted = (log(unif_rand()) < acc.logAlpha);
        }

        if (accepted) {
          tess[j] = std::move(prop.tess);
          tess_nC[j] = prop.nC;
          tess_d[j] = static_cast<int>(prop.dim.size());
          dim_j[j] = std::move(prop.dim);
          caches[j] = std::move(prop_cache);
          sample_mu_into(R_new, n_new, sigSqMu, sigmaSquared, pred[j]);
          for (int obs = 0; obs < n; ++obs)
            lastTessPred[obs] = pred[j][caches[j].assignment[obs]];
        } else {
          sample_mu_into(R_old, n_old, sigSqMu, sigmaSquared, pred[j]);
          for (int obs = 0; obs < n; ++obs)
            lastTessPred[obs] = pred[j][caches[j].assignment[obs]];
        }

        if (j == m - 1) {
          for (int obs = 0; obs < n; ++obs)
            sumAllTess[obs] += lastTessPred[obs];
        }
      }

      double meanCenters = 0.0;
      double meanDims = 0.0;
      double retainedLogLikSum = 0.0;
      for (int j = 0; j < m; ++j) {
        meanCenters += tess_nC[j];
        meanDims += tess_d[j];

        R_old.assign(tess_nC[j], 0.0);
        n_old.assign(tess_nC[j], 0);
        for (int obs = 0; obs < n; ++obs) {
          int cell = caches[j].assignment[obs];
          double tessContribution = pred[j][cell];
          double r = yScaled[obs] - (sumAllTess[obs] - tessContribution);
          R_old[cell] += r;
          n_old[cell]++;
        }
        retainedLogLikSum += tessellation_log_likelihood_component(
          R_old, n_old, sigmaSquared, sigSqMu);
      }
      meanCenters /= m;
      meanDims /= m;

      double sdCenters = 0.0;
      if (m > 1) {
        for (int j = 0; j < m; ++j) {
          double diff = tess_nC[j] - meanCenters;
          sdCenters += diff * diff;
        }
        sdCenters = sqrt(sdCenters / (m - 1));
      }

      int traceIdx = iter - 1;
      INTEGER(outTraceIteration)[traceIdx] = iter;
      LOGICAL(outTraceBurnIn)[traceIdx] = iter <= burnIn;
      REAL(outTraceAvgCenters)[traceIdx] = meanCenters;
      REAL(outTraceSdCenters)[traceIdx] = sdCenters;
      REAL(outTraceAvgDims)[traceIdx] = meanDims;
      REAL(outTraceLogLik)[traceIdx] = retainedLogLikSum / m;

      if (iter > burnIn && (iter - burnIn) % thinning == 0) {
        for (int obs = 0; obs < n; ++obs)
          predictionMatrix[obs + storageIdx * n] = sumAllTess[obs];

        StoredDraw draw;
        draw.sigma = sigmaSquared;
        draw.tessellations.resize(m);
        for (int j = 0; j < m; ++j) {
          StoredTess& st = draw.tessellations[j];
          st.nC = tess_nC[j];
          st.d = tess_d[j];
          st.centres = tess[j];
          st.dim = dim_j[j];
          st.mu = pred[j];
        }
        stored.push_back(std::move(draw));
        storageIdx++;
      }
    }

    PutRNGstate();

    // Pack deferred posterior into R lists (option 8)
    SEXP packed = PROTECT(pack_posterior_lists(stored, m));
    SEXP outTess = VECTOR_ELT(packed, 0);
    SEXP outDim = VECTOR_ELT(packed, 1);
    SEXP outPred = VECTOR_ELT(packed, 2);

    SEXP outSigma = PROTECT(Rf_allocVector(REALSXP, numSamples));
    for (int s = 0; s < numSamples; ++s) REAL(outSigma)[s] = stored[s].sigma;

    SEXP outPredMatrix = PROTECT(Rf_allocMatrix(REALSXP, n, numSamples));
    if (numSamples > 0) {
      std::memcpy(REAL(outPredMatrix), predictionMatrix.data(),
                  static_cast<size_t>(n) * numSamples * sizeof(double));
    }

    SEXP outTraceStats = PROTECT(Rf_allocVector(VECSXP, 6));
    SEXP traceNames = PROTECT(Rf_allocVector(STRSXP, 6));
    SET_VECTOR_ELT(outTraceStats, 0, outTraceIteration);
    SET_VECTOR_ELT(outTraceStats, 1, outTraceBurnIn);
    SET_VECTOR_ELT(outTraceStats, 2, outTraceAvgCenters);
    SET_VECTOR_ELT(outTraceStats, 3, outTraceSdCenters);
    SET_VECTOR_ELT(outTraceStats, 4, outTraceAvgDims);
    SET_VECTOR_ELT(outTraceStats, 5, outTraceLogLik);
    SET_STRING_ELT(traceNames, 0, Rf_mkChar("iteration"));
    SET_STRING_ELT(traceNames, 1, Rf_mkChar("isBurnIn"));
    SET_STRING_ELT(traceNames, 2, Rf_mkChar("averageCentresPerTessellation"));
    SET_STRING_ELT(traceNames, 3, Rf_mkChar("sdCentresPerTessellation"));
    SET_STRING_ELT(traceNames, 4, Rf_mkChar("averageDimensionsPerTessellation"));
    SET_STRING_ELT(traceNames, 5, Rf_mkChar("logLikelihood"));
    Rf_setAttrib(outTraceStats, R_NamesSymbol, traceNames);
    SEXP dataFrameClass = PROTECT(Rf_mkString("data.frame"));
    Rf_setAttrib(outTraceStats, R_ClassSymbol, dataFrameClass);
    SEXP rowNames = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(rowNames)[0] = NA_INTEGER;
    INTEGER(rowNames)[1] = -totalIter;
    Rf_setAttrib(outTraceStats, R_RowNamesSymbol, rowNames);

    SEXP result = PROTECT(Rf_allocVector(VECSXP, 6));
    SEXP listNames = PROTECT(Rf_allocVector(STRSXP, 6));
    SET_VECTOR_ELT(result, 0, outTess);
    SET_VECTOR_ELT(result, 1, outDim);
    SET_VECTOR_ELT(result, 2, outPred);
    SET_VECTOR_ELT(result, 3, outSigma);
    SET_VECTOR_ELT(result, 4, outPredMatrix);
    SET_VECTOR_ELT(result, 5, outTraceStats);
    SET_STRING_ELT(listNames, 0, Rf_mkChar("posteriorTess"));
    SET_STRING_ELT(listNames, 1, Rf_mkChar("posteriorDim"));
    SET_STRING_ELT(listNames, 2, Rf_mkChar("posteriorPred"));
    SET_STRING_ELT(listNames, 3, Rf_mkChar("posteriorSigma"));
    SET_STRING_ELT(listNames, 4, Rf_mkChar("predictionMatrix"));
    SET_STRING_ELT(listNames, 5, Rf_mkChar("traceStats"));
    Rf_setAttrib(result, R_NamesSymbol, listNames);

    UNPROTECT(15); // 6 traces + packed + sigma + predMat + traceStats + names
                   // + class + rownames + result + listNames
    return result;
  }

  // ---------------------------------------------------------------------------
  // addi_vortes_predict_cpp — flattened posterior + shared NN kernel
  // ---------------------------------------------------------------------------
  SEXP addi_vortes_predict_cpp(
      SEXP xNew_sexp,
      SEXP posteriorTess_sexp,
      SEXP posteriorDim_sexp,
      SEXP posteriorPred_sexp,
      SEXP metric_sexp,
      SEXP member_sexp,
      SEXP showProgress_sexp) {

    const double* xNew = REAL(xNew_sexp);
    int n = Rf_nrows(xNew_sexp);
    int p = Rf_ncols(xNew_sexp);
    bool showProgress = LOGICAL(showProgress_sexp)[0];

    FlatPosterior flat = flatten_posterior(
      posteriorTess_sexp, posteriorDim_sexp, posteriorPred_sexp);
    const int numSamples = flat.numSamples;
    if (numSamples == 0) {
      SEXP empty = PROTECT(Rf_allocMatrix(REALSXP, n, 0));
      UNPROTECT(1);
      return empty;
    }
    const int m = flat.m;

    int nMetric = Rf_length(metric_sexp);
    if (nMetric != Rf_length(member_sexp)) {
      Rf_error("metric and member vectors must have the same length");
    }
    std::vector<int> metric(INTEGER(metric_sexp), INTEGER(metric_sexp) + nMetric);
    std::vector<int> members(INTEGER(member_sexp), INTEGER(member_sexp) + nMetric);

    int mem_sum = 0;
    for (int i = 0; i < nMetric; ++i) mem_sum += members[i];
    if (mem_sum != p) {
      Rf_error("Sum of member counts must match number of columns in newdata");
    }

    std::vector<int> metric_aug;
    metric_aug.reserve(p);
    for (int i = 0; i < nMetric; ++i) {
      for (int j = 0; j < members[i]; ++j) metric_aug.push_back(metric[i]);
    }
    std::vector<int> ncats;
    if (in_vector(2, metric)) {
      for (int i = 0; i < p; ++i) {
        if (metric_aug[i] == 2) {
          int max_val = 0;
          for (int j = 0; j < n; ++j) {
            if (xNew[i * n + j] > max_val)
              max_val = static_cast<int>(xNew[i * n + j]);
          }
          ncats.push_back(max_val);
        }
      }
    }

    const bool euclidean = all_euclidean_metric(metric);
    std::vector<double> x_row;
    pack_row_major(xNew, n, p, x_row);

    SEXP outPredMatrix = PROTECT(Rf_allocMatrix(REALSXP, n, numSamples));
    double* p_out = REAL(outPredMatrix);

    AssignScratch assign_scratch;
    AssignmentCache cache;
    std::vector<double> drawPred(n, 0.0);
    std::vector<int> dim1;

    const int progressWidth = 40;
    int lastFilled = -1;

    for (int s = 0; s < numSamples; ++s) {
      maybe_progress("Predict", s + 1, numSamples, progressWidth,
                     lastFilled, showProgress);
      std::fill(drawPred.begin(), drawPred.end(), 0.0);

      for (int j = 0; j < m; ++j) {
        const int idx = s * m + j;
        const int nC = flat.nC[idx];
        const int d = flat.d[idx];
        const double* centres = flat.centres.data() + flat.centre_off[idx];
        const double* mu = flat.mus.data() + flat.mu_off[idx];
        dim1.assign(flat.dims.begin() + flat.dim_off[idx],
                    flat.dims.begin() + flat.dim_off[idx] + d);

        reassign(x_row.data(), n, p, centres, nC, d, dim1,
                 AssignmentDelta::FullRecompute, 0, euclidean,
                 metric, members, ncats,
                 AssignmentCache{}, cache, assign_scratch);

        for (int obs = 0; obs < n; ++obs)
          drawPred[obs] += mu[cache.assignment[obs]];
      }

      for (int obs = 0; obs < n; ++obs)
        p_out[obs + s * n] = drawPred[obs];
    }

    UNPROTECT(1);
    return outPredMatrix;
  }

} // extern "C"
