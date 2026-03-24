// 1. C++ Standard Library headers first
#include <vector>    // For std::vector
#include <string>    // For std::string
#include <numeric>   // For std::accumulate
#include <algorithm> // For std::find
#include <cmath>     // For sqrt
#include <cstring>   // For memcpy

// 2. Add this to prevent R from creating problematic macros
#define R_NO_REMAP

// 3. R headers last
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Random.h> // For unif_rand() and norm_rand()

// Helper function to check if a value is in a vector
bool in_vector(int value, const std::vector<int>& vec) {
  return std::find(vec.begin(), vec.end(), value) != vec.end();
}

int n_elem(int value, const std::vector<int>& vec) {
  int total = 0;
  for (int i = 0; i < vec.size(); ++i) {
    if (vec[i] == value) {
      total++;
    }
  }
  return total;
}

std::vector<int> which_elem(int value, const std::vector<int>& vec) {
  std::vector<int> indexes(n_elem(value, vec));
  int indexes_index = 0;
  for (int i = 0; i < vec.size(); ++i) {
    if (vec[i] == value) {
      indexes[indexes_index] = i;
      indexes_index++;
    }
  }
  return indexes;
}

double period_shift(double val, double lim) {
  while (val >= lim) {
    val -= 2*lim;
  }
  while (val < -lim) {
    val += 2*lim;
  }
  return val;
}

// Computes Euclidean distance between two points
double euclidean_distance(const std::vector<double>& p1, const std::vector<double>& p2) {
  if (p1.size() != p2.size()) {
    Rf_error("Points have incompatible dimensions.");
  }
  double dist = 0.0;
  for (int i = 0; i < static_cast<int>(p1.size()); ++i) {
    const double diff = p1[i] - p2[i];
    dist += diff * diff;
  }
  return dist;
}

// Computes great circle distance on an n-sphere.
// It assumes that the last dimension is the 'azimuthal' distance; i.e. the one with range [-pi, pi].
// Optimisation not great: needs some work with the trigonometric functions.
// The issue is that one needs all the covariates to calculate: you can't just cull the covariates not in the tessellation.
double spherical_distance(const std::vector<double>& p1, const std::vector<double>& p2) {
  if (p1.size() != p2.size()) {
    Rf_error("Points have incompatible dimensions.");
  }
  if (p1.size() == 1) {
    double a1 = abs(p1[0]-p2[0]);
    double a2 = 2*M_PI-a1;
    if (a1 < a2) return(a1*a1);
    return(a2*a2);
  }
  double angle_diff = cos(p1[p1.size()-1]-p2[p2.size()-1]);
  for (int i = p1.size()-2; i >= 0; --i) {
    double internal = sin(p1[i]) * sin(p2[i]) + cos(p1[i]) * cos(p2[i]) * angle_diff;
    if (internal > 1) {
      internal = 1;
    }
    if (internal < -1) {
      internal = -1.0;
    }
    if (i == 0) {
      angle_diff = acos(internal);
    }
    else {
      angle_diff = internal;
    }
  }
  return(angle_diff * angle_diff);
}

extern "C" {
  
  // ---------------------------------------------------------------------------
  // 0. knnx_index_cpp
  // ---------------------------------------------------------------------------
  // This function implements k-nearest neighbors index search to replace FNN::knnx.index
  // The R wrapper function is `knnx.index`.
  SEXP knnx_index_cpp(SEXP tess_sexp, SEXP query_sexp, SEXP k_sexp, SEXP dim_sexp, SEXP dist_sexp) {
    // --- Unpack arguments ---
    double* p_tess = REAL(tess_sexp);
    double* p_query = REAL(query_sexp);
    int k = INTEGER(k_sexp)[0];
    int* dim_p = INTEGER(dim_sexp);
    int* metric_ptr = INTEGER(dist_sexp);

    std::vector<int> dim_p_temp(dim_p, dim_p + Rf_length(dim_sexp));
    std::vector<int> metric(metric_ptr, metric_ptr + Rf_length(dist_sexp));
    
    int tess_rows = Rf_nrows(tess_sexp);
    int tess_cols = Rf_ncols(tess_sexp);
    int query_rows = Rf_nrows(query_sexp);
    int query_cols = Rf_ncols(query_sexp);

    if (static_cast<int>(metric.size()) != query_cols) {
      Rf_error("Length of metric must match number of columns in query/data matrices");
    }
    
    // Check dimensions match
    if (tess_cols != query_cols) {
      Rf_error("Dimensions of tess and query matrices must match");
    }
    
    // Check k is valid
    if (k <= 0 || k > tess_rows) {
      Rf_error("k must be positive and not greater than number of reference points");
    }
    
    // --- Create result matrix ---
    SEXP result;
    PROTECT(result = Rf_allocMatrix(INTSXP, query_rows, k));
    int* p_result = INTEGER(result);

    // Define variables to be used in the loop
    int how_many_E = n_elem(0, metric);
    const bool has_euclidean = how_many_E > 0;
    std::vector<double> t_pt_E(how_many_E);
    std::vector<double> q_pt_E(how_many_E);
    int how_many_S = n_elem(1, metric);
    const bool has_spherical = how_many_S > 0;
    std::vector<double> t_pt_S(how_many_S);
    std::vector<double> q_pt_S(how_many_S);
    std::vector<int> coind(metric.size());
    int count_S = 0;
    int count_E = 0;
    for (int i = 0; i < metric.size(); i++) {
      if (metric[i] == 0) {
        coind[i] = count_E;
        count_E++;
      }
      if (metric[i] == 1) {
        coind[i] = count_S;
        count_S++;
      }
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

    double dval;
    
    // --- Main Logic: For each query point, find k nearest neighbors ---
    for (int q = 0; q < query_rows; ++q) {
      for (int d = 0; d < query_cols; ++d) {
        if (metric[d] == 0) {
          q_pt_E[coind[d]] = p_query[q + d * query_rows];
        }
        if (metric[d] == 1) {
          q_pt_S[coind[d]] = p_query[q + d * query_rows];
        }
      }

      // Calculate distances to all tessellation points
      std::vector<std::pair<double, int>> distances(tess_rows);

      for (int t = 0; t < tess_rows; ++t) {
        dval = 0.0;

        if (has_euclidean) {
          t_pt_E = q_pt_E;
        }
        if (has_spherical) {
          t_pt_S = q_pt_S;
        }

        for (int i = 0; i < static_cast<int>(active_dim_idx.size()); ++i) {
          const int d = active_dim_idx[i];
          if (metric[d] == 0) {
            t_pt_E[coind[d]] = p_tess[t + d * tess_rows];
          }
          if (metric[d] == 1) {
            t_pt_S[coind[d]] = p_tess[t + d * tess_rows];
          }
        }

        if (has_spherical) {
          dval += spherical_distance(q_pt_S, t_pt_S);
        }
        if (has_euclidean) {
          dval += euclidean_distance(q_pt_E, t_pt_E);
        }
        distances[t] = std::make_pair(dval, t + 1); // +1 for R 1-based indexing
      }
      
      // Sort by squared distance to get k nearest neighbors
      std::partial_sort(distances.begin(), distances.begin() + k, distances.end());
      
      // Store the indices of k nearest neighbors
      for (int i = 0; i < k; ++i) {
        p_result[q + i * query_rows] = distances[i].second;
      }
    }
    
    UNPROTECT(1);
    return result;
  }
  
  
  // ---------------------------------------------------------------------------
  // 1. calculate_residuals_cpp
  // ---------------------------------------------------------------------------
  // This function calculates residuals for new centres based on the 
  // provided indexes.
  // The R wrapper function is `calculateResiduals`.
  SEXP calculate_residuals_cpp(SEXP R_j_sexp, SEXP indexes_sexp,
                               SEXP indexesStar_sexp, SEXP num_levels_old_sexp,
                               SEXP num_centres_new_sexp) {
    
    // --- Unpack arguments ---
    double* p_R_j = REAL(R_j_sexp);
    int* p_indexes = INTEGER(indexes_sexp);
    int* p_indexesStar = INTEGER(indexesStar_sexp);
    int num_levels_old = INTEGER(num_levels_old_sexp)[0];
    int num_centres_new = INTEGER(num_centres_new_sexp)[0];
    int n_obs = Rf_length(R_j_sexp);
    
    // --- Create C++ vectors for processing ---
    std::vector<double> R_ijOld(num_levels_old, 0.0);
    std::vector<int> n_ijOld(num_levels_old, 0);
    std::vector<double> R_ijNew(num_centres_new, 0.0);
    std::vector<int> n_ijNew(num_centres_new, 0);
    
    // --- Main Logic ---
    for(int i = 0; i < n_obs; ++i) {
      int idx_old = p_indexes[i] - 1; // R is 1-based
      if (idx_old >= 0 && idx_old < num_levels_old) {
        R_ijOld[idx_old] += p_R_j[i];
        n_ijOld[idx_old]++;
      }
      
      int idx_new = p_indexesStar[i] - 1; // R is 1-based
      if (idx_new >= 0 && idx_new < num_centres_new) {
        R_ijNew[idx_new] += p_R_j[i];
        n_ijNew[idx_new]++;
      }
    }
    
    // --- Pack results into a named list for R ---
    SEXP res_R_ijOld, res_n_ijOld, res_R_ijNew, res_n_ijNew, result_list, list_names;
    
    PROTECT(res_R_ijOld = Rf_allocVector(REALSXP, num_levels_old));
    memcpy(REAL(res_R_ijOld), R_ijOld.data(), num_levels_old * sizeof(double));
    
    PROTECT(res_n_ijOld = Rf_allocVector(INTSXP, num_levels_old));
    memcpy(INTEGER(res_n_ijOld), n_ijOld.data(), num_levels_old * sizeof(int));
    
    PROTECT(res_R_ijNew = Rf_allocVector(REALSXP, num_centres_new));
    memcpy(REAL(res_R_ijNew), R_ijNew.data(), num_centres_new * sizeof(double));
    
    PROTECT(res_n_ijNew = Rf_allocVector(INTSXP, num_centres_new));
    memcpy(INTEGER(res_n_ijNew), n_ijNew.data(), num_centres_new * sizeof(int));
    
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
    
    UNPROTECT(6); // 4 result vectors + list + names
    return result_list;
  }
  
  
  // ---------------------------------------------------------------------------
  // 2. propose_tessellation_cpp
  // ---------------------------------------------------------------------------
  // This function proposes a new tessellation based on the current one,
  // modifying it according to a set of rules and a random number generator.
  // The R wrapper function is `proposeTessellation`.
  SEXP propose_tessellation_cpp(SEXP tess_j_sexp, SEXP dim_j_sexp, SEXP sd_sexp, SEXP mu_sexp, SEXP num_cov_sexp, SEXP metric_sexp) {
    
    // --- Unpack arguments ---
    double* p_tess_j = REAL(tess_j_sexp);
    int* p_dim_j = INTEGER(dim_j_sexp);
    double* sd = REAL(sd_sexp);
    double* mu = REAL(mu_sexp);
    int numCovariates = INTEGER(num_cov_sexp)[0];
    int* metric_ptr = INTEGER(metric_sexp);

    std::vector<int> metric(metric_ptr, metric_ptr + Rf_length(metric_sexp));
    
    int tess_j_rows = Rf_nrows(tess_j_sexp);
    int d_j_length = Rf_length(dim_j_sexp);
    
    // Get R's random number generator state
    GetRNGstate();
    
    // --- Create C++ copies for manipulation ---
    std::vector<int> dim_j_star(p_dim_j, p_dim_j + d_j_length);
    std::vector<double> tess_j_star(p_tess_j, p_tess_j + (tess_j_rows * d_j_length));
    std::string modification = "Change";

    double new_val = 0.0;
    
    double p = unif_rand();

    std::vector<int> sphere_index;
    if (in_vector(1, metric)) {
      sphere_index = which_elem(1, metric);
    }
    
    // --- Main Logic ---
    // Add Dimension (AD): ensure we don't try to add a dimension when all covariates are already selected
    if ((p < 0.2 && d_j_length != numCovariates) || (d_j_length == 1 && d_j_length != numCovariates && p < 0.4)) {
      modification = "AD";
      int new_dim;
      do {
        new_dim = floor(unif_rand() * numCovariates) + 1;
      } while (in_vector(new_dim, dim_j_star));
      
      dim_j_star.push_back(new_dim);
      
      std::vector<double> new_tess(tess_j_rows * (d_j_length + 1));
      for (int r = 0; r < tess_j_rows; ++r) {
        for (int c = 0; c < d_j_length; ++c) {
          new_tess[r + c * tess_j_rows] = tess_j_star[r + c * tess_j_rows];
        }
        new_val = mu[new_dim-1] + norm_rand() * sd[new_dim-1];
        if (metric[new_dim-1] == 1) {
          if (new_dim - 1 == sphere_index[sphere_index.size()-1]) {
            new_val = period_shift(new_val, M_PI);
          }
          // else {
          //   new_val = period_shift(new_val, M_PI_2);
          // }
        }
        new_tess[r + d_j_length * tess_j_rows] = new_val;
      }
      tess_j_star = new_tess;
      
    } else if (p < 0.4 && d_j_length > 1) {
      modification = "RD";
      int removed_dim_idx = floor(unif_rand() * d_j_length);
      dim_j_star.erase(dim_j_star.begin() + removed_dim_idx);
      
      std::vector<double> new_tess(tess_j_rows * (d_j_length - 1));
      int current_col = 0;
      for (int c = 0; c < d_j_length; ++c) {
        if (c != removed_dim_idx) {
          for (int r = 0; r < tess_j_rows; ++r) {
            new_tess[r + current_col * tess_j_rows] = tess_j_star[r + c * tess_j_rows];
          }
          current_col++;
        }
      }
      tess_j_star = new_tess;
      
    } else if (p < 0.6 || (p < 0.8 && tess_j_rows == 1)) {
      modification = "AC";
      for (int i = 0; i < d_j_length; ++i) {
        new_val = mu[i] + norm_rand() * sd[i];
        if (metric[i] == 1) {
          if (i == sphere_index[sphere_index.size()-1]) {
            new_val = period_shift(new_val, M_PI);
          }
          // else {
          //   new_val = period_shift(new_val, M_PI_2);
          // }
        }
        tess_j_star.insert(tess_j_star.begin() + (i * (tess_j_rows + 1)) + tess_j_rows, new_val);
      }
      
    } else if (p < 0.8 && tess_j_rows > 1) {
      modification = "RC";
      int removed_row_idx = floor(unif_rand() * tess_j_rows);
      std::vector<double> new_tess;
      new_tess.reserve((tess_j_rows - 1) * d_j_length);
      for(int c = 0; c < d_j_length; ++c){
        for(int r = 0; r < tess_j_rows; ++r){
          if(r != removed_row_idx) {
            new_tess.push_back(tess_j_star[r + c * tess_j_rows]);
          }
        }
      }
      tess_j_star = new_tess;
      
    } else if (p < 0.9 || d_j_length == numCovariates) {
      int centre_to_change_idx = floor(unif_rand() * tess_j_rows);
      for (int c = 0; c < d_j_length; ++c) {
        new_val = mu[c] + norm_rand() * sd[c];
        if (metric[c] == 1) {
          if (c == sphere_index[sphere_index.size()-1]) {
            new_val = period_shift(new_val, M_PI);
          }
          // else {
          //   new_val = period_shift(new_val, M_PI_2);
          // }
        }
        tess_j_star[centre_to_change_idx + c * tess_j_rows] = new_val;
      }
      
    } else {
      modification = "Swap";
      int dim_to_change_idx = floor(unif_rand() * d_j_length);
      int new_dim;
      do {
        new_dim = floor(unif_rand() * numCovariates) + 1;
      } while (in_vector(new_dim, dim_j_star));
      
      dim_j_star[dim_to_change_idx] = new_dim;
      for (int r = 0; r < tess_j_rows; ++r) {
        new_val = mu[dim_to_change_idx] + norm_rand() * sd[dim_to_change_idx];
        if (metric[dim_to_change_idx] == 1) {
          if (dim_to_change_idx == sphere_index[sphere_index.size()-1]) {
            new_val = period_shift(new_val, M_PI);
          }
          // else {
          //   new_val = period_shift(new_val, M_PI_2);
          // }
        }
        tess_j_star[r + dim_to_change_idx * tess_j_rows] = new_val;
      }
    }
    
    // Update R's random number generator state
    PutRNGstate();
    
    // --- Pack results into a named list for R ---
    SEXP res_tess_j_star, res_dim_j_star, res_mod, result_list, list_names;
    
    int new_rows = tess_j_star.size() / dim_j_star.size();
    PROTECT(res_tess_j_star = Rf_allocMatrix(REALSXP, new_rows, dim_j_star.size()));
    memcpy(REAL(res_tess_j_star), tess_j_star.data(), tess_j_star.size() * sizeof(double));
    
    PROTECT(res_dim_j_star = Rf_allocVector(INTSXP, dim_j_star.size()));
    memcpy(INTEGER(res_dim_j_star), dim_j_star.data(), dim_j_star.size() * sizeof(int));
    
    PROTECT(res_mod = Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(res_mod, 0, Rf_mkChar(modification.c_str()));
    
    PROTECT(result_list = Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(result_list, 0, res_tess_j_star);
    SET_VECTOR_ELT(result_list, 1, res_dim_j_star);
    SET_VECTOR_ELT(result_list, 2, res_mod);
    
    PROTECT(list_names = Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(list_names, 0, Rf_mkChar("tess_j_star"));
    SET_STRING_ELT(list_names, 1, Rf_mkChar("dim_j_star"));
    SET_STRING_ELT(list_names, 2, Rf_mkChar("Modification"));
    Rf_setAttrib(result_list, R_NamesSymbol, list_names);
    
    UNPROTECT(5); // 3 result vectors + list + names
    return result_list;
  }
  
} // extern "C"

// =============================================================================
// Internal helpers for the unified MCMC function
// =============================================================================

// Find the nearest centre (k=1 NN) for every observation.
// obs_data : n x p  column-major double array
// centres  : nC x d column-major double array (active dims only)
// dim1     : d active dimension indices, 1-BASED
// Returns  : n vector of 0-based centre indices
static std::vector<int> knn1_internal(
    const double* obs_data, int n, int p,
    const double* centres, int nC, int d,
    const std::vector<int>& dim1,
    const std::vector<int>& metric,
    const std::vector<int>& coind,
    int nE, int nS) {

  std::vector<int> result(n, 0);
  if (nC == 1) return result;

  const bool hasE = (nE > 0), hasS = (nS > 0);
  std::vector<double> q_E(nE), q_S(nS), t_E(nE), t_S(nS);

  for (int obs = 0; obs < n; obs++) {
    for (int g = 0; g < p; g++) {
      double v = obs_data[obs + g * n];
      if (metric[g] == 0) q_E[coind[g]] = v;
      else                 q_S[coind[g]] = v;
    }
    double best = 1e300;
    int best_c = 0;
    for (int c = 0; c < nC; c++) {
      if (hasE) t_E = q_E;
      if (hasS) t_S = q_S;
      for (int di = 0; di < d; di++) {
        int g = dim1[di] - 1;  // 0-based global index
        if (metric[g] == 0) t_E[coind[g]] = centres[c + di * nC];
        else                 t_S[coind[g]] = centres[c + di * nC];
      }
      double dist = 0.0;
      if (hasS) dist += spherical_distance(q_S, t_S);
      if (hasE) dist += euclidean_distance(q_E, t_E);
      if (dist < best) { best = dist; best_c = c; }
    }
    result[obs] = best_c;
  }
  return result;
}

// Aggregate partial residuals R_j into per-cell sums and counts.
static void aggregate_residuals(
    const std::vector<double>& R_j,
    const std::vector<int>& idx_old, int nC_old,
    const std::vector<int>& idx_new, int nC_new,
    std::vector<double>& R_old, std::vector<int>& n_old,
    std::vector<double>& R_new, std::vector<int>& n_new) {

  R_old.assign(nC_old, 0.0); n_old.assign(nC_old, 0);
  R_new.assign(nC_new, 0.0); n_new.assign(nC_new, 0);
  for (int obs = 0; obs < (int)R_j.size(); obs++) {
    R_old[idx_old[obs]] += R_j[obs]; n_old[idx_old[obs]]++;
    R_new[idx_new[obs]] += R_j[obs]; n_new[idx_new[obs]]++;
  }
}

// Compute the log acceptance probability.
static double log_acceptance_prob(
    const std::vector<double>& R_old, const std::vector<int>& n_old,
    const std::vector<double>& R_new, const std::vector<int>& n_new,
    int d_new, int nC_new,
    double sigmaSquared, double sigmaSquaredMu,
    double omega, double lambdaRate, int p,
    const std::string& mod) {

  double sum_log_old = 0.0, sum_R_old = 0.0;
  for (int k = 0; k < (int)n_old.size(); k++) {
    double den = n_old[k] * sigmaSquaredMu + sigmaSquared;
    sum_log_old += log(den);
    sum_R_old   += (R_old[k] * R_old[k]) / den;
  }
  double sum_log_new = 0.0, sum_R_new = 0.0;
  for (int k = 0; k < (int)n_new.size(); k++) {
    double den = n_new[k] * sigmaSquaredMu + sigmaSquared;
    sum_log_new += log(den);
    sum_R_new   += (R_new[k] * R_new[k]) / den;
  }

  double log_lik = 0.5 * (sum_log_old - sum_log_new) +
    (sigmaSquaredMu / (2.0 * sigmaSquared)) * (sum_R_new - sum_R_old);

  const double prob_eps = 1e-10;
  double prob = std::min(1.0 - prob_eps, std::max(0.0, omega / p));

  double acc = log_lik;

  if (mod == "AD") {
    double log_ts_tr = Rf_dbinom(d_new - 1, p - 1, prob, 1)
                     - Rf_dbinom(d_new - 2, p - 1, prob, 1)
                     - log((double)d_new);
    acc += log_ts_tr;
    if (d_new == 1)     acc += log(0.5);
    else if (d_new == p - 1) acc += log(2.0);

  } else if (mod == "RD") {
    double log_ts_tr = Rf_dbinom(d_new - 1, p, prob, 1)
                     - Rf_dbinom(d_new,     p, prob, 1)
                     + log((double)(d_new + 1));
    acc += log_ts_tr;
    if (d_new == p) acc += log(0.5);
    else if (d_new == 2) acc += log(2.0);

  } else if (mod == "AC") {
    double log_ts_tr = Rf_dpois(nC_new - 1, lambdaRate, 1)
                     - Rf_dpois(nC_new - 2, lambdaRate, 1)
                     - log((double)nC_new);
    acc += log_ts_tr + 0.5 * log(sigmaSquared);
    if (nC_new == 1) acc += log(0.5);

  } else if (mod == "RC") {
    double log_ts_tr = Rf_dpois(nC_new - 1, lambdaRate, 1)
                     - Rf_dpois(nC_new,     lambdaRate, 1)
                     + log((double)(nC_new + 1));
    acc += log_ts_tr - 0.5 * log(sigmaSquared);
    if (nC_new == 2) acc += log(2.0);
  }
  // "Change" and "Swap": log(TessStructure * TransitionRatio) = 0

  return acc;
}

// Sample mu values for all centres of tessellation j.
static std::vector<double> sample_mu_internal(
    const std::vector<double>& R_ij, const std::vector<int>& n_ij,
    double sigmaSquaredMu, double sigmaSquared) {

  int N = (int)R_ij.size();
  std::vector<double> result(N);
  for (int k = 0; k < N; k++) {
    double den  = sigmaSquaredMu * n_ij[k] + sigmaSquared;
    double mean = (sigmaSquaredMu * R_ij[k]) / den;
    double sd   = sqrt((sigmaSquared * sigmaSquaredMu) / den);
    result[k]   = mean + norm_rand() * sd;
  }
  return result;
}

// Propose a new tessellation. Mirrors propose_tessellation_cpp exactly,
// operating on C++ vectors to avoid R↔C++ crossing.
struct ProposalResult {
  std::vector<double> tess;
  int nC;
  std::vector<int> dim;  // 1-based covariate indices
  std::string mod;
};

static ProposalResult propose_internal(
    const std::vector<double>& tess_j, int nC, int d_j,
    const std::vector<int>& dim_j,   // 1-based
    int p,
    const double* sd, const double* mus,
    const std::vector<int>& metric,
    const std::vector<int>& sphere_index) {  // 0-based spherical dim indices

  ProposalResult r;
  r.tess = tess_j; r.nC = nC; r.dim = dim_j; r.mod = "Change";

  double prand = unif_rand();
  double new_val;

  if ((prand < 0.2 && d_j != p) || (d_j == 1 && d_j != p && prand < 0.4)) {
    // Add Dimension
    r.mod = "AD";
    int new_dim;
    do { new_dim = (int)(unif_rand() * p) + 1; }
    while (in_vector(new_dim, r.dim));
    r.dim.push_back(new_dim);

    std::vector<double> new_tess(nC * (d_j + 1));
    for (int row = 0; row < nC; row++) {
      for (int col = 0; col < d_j; col++)
        new_tess[row + col * nC] = tess_j[row + col * nC];
      new_val = mus[new_dim - 1] + norm_rand() * sd[new_dim - 1];
      if (metric[new_dim - 1] == 1)
        if ((new_dim - 1) == sphere_index.back())
          new_val = period_shift(new_val, M_PI);
      new_tess[row + d_j * nC] = new_val;
    }
    r.tess = new_tess; r.nC = nC;

  } else if (prand < 0.4 && d_j > 1) {
    // Remove Dimension
    r.mod = "RD";
    int rm_idx = (int)(unif_rand() * d_j);
    r.dim.erase(r.dim.begin() + rm_idx);

    std::vector<double> new_tess(nC * (d_j - 1));
    int cur_col = 0;
    for (int col = 0; col < d_j; col++) {
      if (col != rm_idx) {
        for (int row = 0; row < nC; row++)
          new_tess[row + cur_col * nC] = tess_j[row + col * nC];
        cur_col++;
      }
    }
    r.tess = new_tess; r.nC = nC;

  } else if (prand < 0.6 || (prand < 0.8 && nC == 1)) {
    // Add Centre
    r.mod = "AC";
    r.tess = tess_j;
    for (int i = 0; i < d_j; i++) {
      // NOTE: mus/sd are indexed by local position i, and metric is checked
      // using i rather than the global covariate index dim_j[i]-1.
      // This mirrors the original propose_tessellation_cpp behaviour exactly.
      new_val = mus[i] + norm_rand() * sd[i];
      if (metric[i] == 1)
        if (i == (int)sphere_index.back())
          new_val = period_shift(new_val, M_PI);
      r.tess.insert(r.tess.begin() + (i * (nC + 1)) + nC, new_val);
    }
    r.nC = nC + 1;

  } else if (prand < 0.8 && nC > 1) {
    // Remove Centre
    r.mod = "RC";
    int rm_row = (int)(unif_rand() * nC);
    std::vector<double> new_tess;
    new_tess.reserve((nC - 1) * d_j);
    for (int col = 0; col < d_j; col++)
      for (int row = 0; row < nC; row++)
        if (row != rm_row) new_tess.push_back(tess_j[row + col * nC]);
    r.tess = new_tess; r.nC = nC - 1;

  } else if (prand < 0.9 || d_j == p) {
    // Change Centre — same local-index convention as propose_tessellation_cpp
    int ci = (int)(unif_rand() * nC);
    for (int col = 0; col < d_j; col++) {
      new_val = mus[col] + norm_rand() * sd[col];
      if (metric[col] == 1)
        if (col == (int)sphere_index.back())
          new_val = period_shift(new_val, M_PI);
      r.tess[ci + col * nC] = new_val;
    }

  } else {
    // Swap Dimension — same local-index convention as propose_tessellation_cpp
    r.mod = "Swap";
    int swap_idx = (int)(unif_rand() * d_j);
    int new_dim;
    do { new_dim = (int)(unif_rand() * p) + 1; }
    while (in_vector(new_dim, r.dim));
    r.dim[swap_idx] = new_dim;
    for (int row = 0; row < nC; row++) {
      new_val = mus[swap_idx] + norm_rand() * sd[swap_idx];
      if (metric[swap_idx] == 1)
        if (swap_idx == (int)sphere_index.back())
          new_val = period_shift(new_val, M_PI);
      r.tess[row + swap_idx * nC] = new_val;
    }
  }

  return r;
}

extern "C" {

  // ---------------------------------------------------------------------------
  // 4. addi_vortes_mcmc_cpp
  // ---------------------------------------------------------------------------
  // Runs the complete MCMC loop for model fitting.  All per-iteration R↔C++
  // crossings are eliminated: a single call returns all posterior samples.
  //
  // Arguments
  //   xScaled_sexp        n x p double matrix (column-major)
  //   yScaled_sexp        n double vector
  //   metric_sexp         p integer vector  (0=Euclidean, 1=Spherical)
  //   m_sexp              integer — number of tessellations
  //   totalMCMCIter_sexp  integer
  //   mcmcBurnIn_sexp     integer
  //   thinning_sexp       integer
  //   nu_sexp             double — degrees of freedom for sigma prior
  //   lambda_sexp         double — scale for sigma prior
  //   sigmaSquaredMu_sexp double — prior variance on mu values
  //   omega_sexp          double — dimension inclusion prior (Omega)
  //   lambdaRate_sexp     double — Poisson rate for number of centres
  //   sd_sexp             p double vector — proposal s.d. per covariate
  //   mus_sexp            n double vector — proposal mean per covariate
  //   init_tess_sexp      R list of m matrices (nC_j x d_j, 1-based dims)
  //   init_dim_sexp       R list of m integer vectors (1-based)
  //   init_pred_sexp      R list of m double vectors (mu values)
  //   binaryCols_sexp     integer vector of binary column indices (1-based),
  //                       or R_NilValue when there are no categorical covariates
  //   catScaling_sexp     double — upper bound for binary column clamping
  //
  // Returns a named R list:
  //   posteriorTess    — list[numSamples] of list[m] of matrices
  //   posteriorDim     — list[numSamples] of list[m] of integer vectors
  //   posteriorPred    — list[numSamples] of list[m] of double vectors
  //   posteriorSigma   — double vector[numSamples]
  //   predictionMatrix — n x numSamples double matrix
  SEXP addi_vortes_mcmc_cpp(
      SEXP xScaled_sexp,
      SEXP yScaled_sexp,
      SEXP metric_sexp,
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
      SEXP catScaling_sexp) {

    // -------------------------------------------------------------------------
    // 1. Unpack scalar / vector inputs
    // -------------------------------------------------------------------------
    const double* xScaled  = REAL(xScaled_sexp);
    const double* yScaled  = REAL(yScaled_sexp);
    int n  = Rf_nrows(xScaled_sexp);
    int p  = Rf_ncols(xScaled_sexp);
    int m  = INTEGER(m_sexp)[0];
    int totalIter  = INTEGER(totalMCMCIter_sexp)[0];
    int burnIn     = INTEGER(mcmcBurnIn_sexp)[0];
    int thinning   = INTEGER(thinning_sexp)[0];
    double nu           = REAL(nu_sexp)[0];
    double lambda       = REAL(lambda_sexp)[0];
    double sigSqMu      = REAL(sigmaSquaredMu_sexp)[0];
    double omega        = REAL(omega_sexp)[0];
    double lambdaRate   = REAL(lambdaRate_sexp)[0];
    const double* sd    = REAL(sd_sexp);
    const double* mus   = REAL(mus_sexp);
    double catScaling   = REAL(catScaling_sexp)[0];

    int* metric_ptr = INTEGER(metric_sexp);
    std::vector<int> metric(metric_ptr, metric_ptr + p);

    // Binary column indices (0-based internally)
    std::vector<int> binaryCols;
    if (!Rf_isNull(binaryCols_sexp)) {
      int* bc = INTEGER(binaryCols_sexp);
      int nb  = Rf_length(binaryCols_sexp);
      for (int i = 0; i < nb; i++) binaryCols.push_back(bc[i] - 1);
    }

    // Precompute coind: global dim index -> position in E or S sub-vector
    std::vector<int> coind(p, 0);
    int cntE = 0, cntS = 0;
    for (int g = 0; g < p; g++) {
      if (metric[g] == 0) coind[g] = cntE++;
      else                 coind[g] = cntS++;
    }
    int nE = cntE, nS = cntS;

    // Precompute 0-based spherical dimension indices (mirrors which_elem in C++)
    std::vector<int> sphere_index;
    for (int g = 0; g < p; g++)
      if (metric[g] == 1) sphere_index.push_back(g);

    // -------------------------------------------------------------------------
    // 2. Unpack initial tessellation state from R lists
    // -------------------------------------------------------------------------
    std::vector<std::vector<double>> tess(m);
    std::vector<int>                 tess_nC(m), tess_d(m);
    std::vector<std::vector<int>>    dim_j(m);
    std::vector<std::vector<double>> pred(m);

    for (int j = 0; j < m; j++) {
      SEXP t_j = VECTOR_ELT(init_tess_sexp, j);
      int rows = Rf_nrows(t_j), cols = Rf_ncols(t_j);
      tess_nC[j] = rows; tess_d[j] = cols;
      double* pt = REAL(t_j);
      tess[j].assign(pt, pt + rows * cols);

      SEXP d_j = VECTOR_ELT(init_dim_sexp, j);
      int* pd = INTEGER(d_j);
      int  nd = Rf_length(d_j);
      dim_j[j].assign(pd, pd + nd);

      SEXP p_j = VECTOR_ELT(init_pred_sexp, j);
      double* pp = REAL(p_j);
      int    np  = Rf_length(p_j);
      pred[j].assign(pp, pp + np);
    }

    // -------------------------------------------------------------------------
    // 3. Initial sumOfAllTess and cell indices
    // -------------------------------------------------------------------------
    std::vector<double> sumAllTess(n, 0.0);
    // sumOfAllTess starts as the sum of pred[[j]][indices[[j]]] for all j
    std::vector<std::vector<int>> curIdx(m);
    for (int j = 0; j < m; j++) {
      curIdx[j] = knn1_internal(
        xScaled, n, p,
        tess[j].data(), tess_nC[j], tess_d[j], dim_j[j],
        metric, coind, nE, nS);
      for (int obs = 0; obs < n; obs++)
        sumAllTess[obs] += pred[j][curIdx[j][obs]];
    }

    // -------------------------------------------------------------------------
    // 4. Allocate output storage
    // -------------------------------------------------------------------------
    int numSamples = 0;
    if (totalIter > burnIn)
      numSamples = (totalIter - burnIn) / thinning;
    if (numSamples < 0) numSamples = 0;

    SEXP outTess  = PROTECT(Rf_allocVector(VECSXP,  numSamples));
    SEXP outDim   = PROTECT(Rf_allocVector(VECSXP,  numSamples));
    SEXP outPred  = PROTECT(Rf_allocVector(VECSXP,  numSamples));
    SEXP outSigma = PROTECT(Rf_allocVector(REALSXP, numSamples));
    SEXP outPredMatrix = PROTECT(Rf_allocMatrix(REALSXP, n, numSamples));
    double* p_outPredMatrix = REAL(outPredMatrix);

    // -------------------------------------------------------------------------
    // 5. MCMC loop
    // -------------------------------------------------------------------------
    GetRNGstate();

    double sigmaSquared = 1.0;
    std::vector<double> lastTessPred(n, 0.0);
    int storageIdx = 0;

    for (int iter = 1; iter <= totalIter; iter++) {

      // Sample sigma squared from inverse-gamma
      double sum_sq = 0.0;
      for (int obs = 0; obs < n; obs++) {
        double r = yScaled[obs] - sumAllTess[obs];
        sum_sq += r * r;
      }
      double shape = (nu + n) / 2.0;
      double rate  = (nu * lambda + sum_sq) / 2.0;
      sigmaSquared = 1.0 / Rf_rgamma(shape, 1.0 / rate);

      for (int j = 0; j < m; j++) {

        // Update sumAllTess to exclude tessellation j
        if (j == 0) {
          for (int obs = 0; obs < n; obs++)
            sumAllTess[obs] -= pred[j][curIdx[j][obs]];
        } else {
          for (int obs = 0; obs < n; obs++)
            sumAllTess[obs] += lastTessPred[obs] - pred[j][curIdx[j][obs]];
        }

        // Partial residuals: R_j = y - sumAllTess (excluding j)
        std::vector<double> R_j(n);
        for (int obs = 0; obs < n; obs++)
          R_j[obs] = yScaled[obs] - sumAllTess[obs];

        // Propose new tessellation
        ProposalResult prop = propose_internal(
          tess[j], tess_nC[j], tess_d[j], dim_j[j],
          p, sd, mus, metric, sphere_index);

        // Clamp binary columns to [0, catScaling] in the proposal
        if (!binaryCols.empty()) {
          int d_star = (int)prop.dim.size();
          for (int di = 0; di < d_star; di++) {
            int g0 = prop.dim[di] - 1;  // 0-based global index
            if (in_vector(g0, binaryCols)) {
              for (int row = 0; row < prop.nC; row++) {
                double v = prop.tess[row + di * prop.nC];
                if (v < 0.0)        v = 0.0;
                if (v > catScaling) v = catScaling;
                prop.tess[row + di * prop.nC] = v;
              }
            }
          }
        }

        // Cell assignments for proposed tessellation
        std::vector<int> idxStar = knn1_internal(
          xScaled, n, p,
          prop.tess.data(), prop.nC, (int)prop.dim.size(), prop.dim,
          metric, coind, nE, nS);

        // Aggregate residuals for old and new tessellations
        std::vector<double> R_old, R_new;
        std::vector<int>    n_old, n_new;
        aggregate_residuals(R_j,
          curIdx[j], tess_nC[j],
          idxStar,   prop.nC,
          R_old, n_old, R_new, n_new);

        // Accept only if no empty cells in proposal
        bool hasEmpty = false;
        for (int k = 0; k < prop.nC; k++)
          if (n_new[k] == 0) { hasEmpty = true; break; }

        bool accepted = false;
        if (!hasEmpty) {
          double logAlpha = log_acceptance_prob(
            R_old, n_old, R_new, n_new,
            (int)prop.dim.size(), prop.nC,
            sigmaSquared, sigSqMu,
            omega, lambdaRate, p,
            prop.mod);
          accepted = (log(unif_rand()) < logAlpha);
        }

        if (accepted) {
          tess[j]    = prop.tess;
          tess_nC[j] = prop.nC;
          tess_d[j]  = (int)prop.dim.size();
          dim_j[j]   = prop.dim;
          curIdx[j]  = idxStar;
          pred[j]    = sample_mu_internal(R_new, n_new, sigSqMu, sigmaSquared);
          for (int obs = 0; obs < n; obs++)
            lastTessPred[obs] = pred[j][idxStar[obs]];
        } else {
          pred[j] = sample_mu_internal(R_old, n_old, sigSqMu, sigmaSquared);
          for (int obs = 0; obs < n; obs++)
            lastTessPred[obs] = pred[j][curIdx[j][obs]];
        }

        // After last tessellation, restore sumAllTess to full sum
        if (j == m - 1) {
          for (int obs = 0; obs < n; obs++)
            sumAllTess[obs] += lastTessPred[obs];
        }
      } // end j loop

      // Store posterior sample (post burn-in, respecting thinning)
      if (iter > burnIn && (iter - burnIn) % thinning == 0) {
        // predictionMatrix column
        for (int obs = 0; obs < n; obs++)
          p_outPredMatrix[obs + storageIdx * n] = sumAllTess[obs];

        REAL(outSigma)[storageIdx] = sigmaSquared;

        // Build sample-level list[m] for tess, dim, pred
        SEXP sampleTess = PROTECT(Rf_allocVector(VECSXP, m));
        SEXP sampleDim  = PROTECT(Rf_allocVector(VECSXP, m));
        SEXP samplePred = PROTECT(Rf_allocVector(VECSXP, m));

        for (int j = 0; j < m; j++) {
          int nCj = tess_nC[j], dj = tess_d[j];

          SEXP rTess = PROTECT(Rf_allocMatrix(REALSXP, nCj, dj));
          memcpy(REAL(rTess), tess[j].data(), nCj * dj * sizeof(double));
          SET_VECTOR_ELT(sampleTess, j, rTess);
          UNPROTECT(1);

          SEXP rDim = PROTECT(Rf_allocVector(INTSXP, dj));
          memcpy(INTEGER(rDim), dim_j[j].data(), dj * sizeof(int));
          SET_VECTOR_ELT(sampleDim, j, rDim);
          UNPROTECT(1);

          SEXP rPred = PROTECT(Rf_allocVector(REALSXP, nCj));
          memcpy(REAL(rPred), pred[j].data(), nCj * sizeof(double));
          SET_VECTOR_ELT(samplePred, j, rPred);
          UNPROTECT(1);
        }

        SET_VECTOR_ELT(outTess, storageIdx, sampleTess);
        SET_VECTOR_ELT(outDim,  storageIdx, sampleDim);
        SET_VECTOR_ELT(outPred, storageIdx, samplePred);
        UNPROTECT(3); // sampleTess, sampleDim, samplePred — protected by outer lists

        storageIdx++;
      }
    } // end iter loop

    PutRNGstate();

    // -------------------------------------------------------------------------
    // 6. Build and return named result list
    // -------------------------------------------------------------------------
    SEXP result    = PROTECT(Rf_allocVector(VECSXP, 5));
    SEXP listNames = PROTECT(Rf_allocVector(STRSXP, 5));

    SET_VECTOR_ELT(result, 0, outTess);
    SET_VECTOR_ELT(result, 1, outDim);
    SET_VECTOR_ELT(result, 2, outPred);
    SET_VECTOR_ELT(result, 3, outSigma);
    SET_VECTOR_ELT(result, 4, outPredMatrix);

    SET_STRING_ELT(listNames, 0, Rf_mkChar("posteriorTess"));
    SET_STRING_ELT(listNames, 1, Rf_mkChar("posteriorDim"));
    SET_STRING_ELT(listNames, 2, Rf_mkChar("posteriorPred"));
    SET_STRING_ELT(listNames, 3, Rf_mkChar("posteriorSigma"));
    SET_STRING_ELT(listNames, 4, Rf_mkChar("predictionMatrix"));
    Rf_setAttrib(result, R_NamesSymbol, listNames);

    UNPROTECT(7); // result, listNames, outTess, outDim, outPred, outSigma, outPredMatrix
    return result;
  }

} // extern "C" (second block)