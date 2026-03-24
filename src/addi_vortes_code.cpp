// 1. C++ Standard Library headers first
#include <vector>    // For std::vector
#include <string>    // For std::string
#include <numeric>   // For std::accumulate
#include <algorithm> // For std::find
#include <cmath>     // For sqrt
#include <cstring>   // For memcpy
#include <set>       // For std::set

// 2. Add this to prevent R from creating problematic macros
#define R_NO_REMAP

// 3. R headers last
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Random.h> // For unif_rand() and norm_rand()
#include <Rmath.h>        // For rgamma

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

// ============================================================
// Internal helpers used only by addi_vortes_mcmc_cpp
// ============================================================

// Modification type constants
static const int MOD_AD     = 0;
static const int MOD_RD     = 1;
static const int MOD_AC     = 2;
static const int MOD_RC     = 3;
static const int MOD_CHANGE = 4;
static const int MOD_SWAP   = 5;

// knn1_internal: assign each of the n observations (rows of xScaled, n x p column-major)
// to its nearest centre in tess_j (nCentres x nDims column-major, with dim_j giving the
// 0-based global column indices).  Writes 0-based centre indices into result (length n).
static void knn1_internal(
  const double* xScaled, int n, int p,
  const double* tess_j, int nCentres, int nDims,
  const std::vector<int>& dim_j,
  const std::vector<int>& metric,
  int nE, int nS,
  const std::vector<int>& coind,
  std::vector<int>& result
) {
  if (nCentres == 1) {
    std::fill(result.begin(), result.end(), 0);
    return;
  }

  const bool hasE = nE > 0;
  const bool hasS = nS > 0;
  std::vector<double> q_pt_E(nE), q_pt_S(nS);
  std::vector<double> t_pt_E(nE), t_pt_S(nS);

  for (int obs = 0; obs < n; obs++) {
    // Extract observation coordinates into E and S sub-vectors
    for (int d = 0; d < p; d++) {
      if (metric[d] == 0) {
        q_pt_E[coind[d]] = xScaled[obs + d * n];
      } else {
        q_pt_S[coind[d]] = xScaled[obs + d * n];
      }
    }

    double best_dist = 1e300;
    int best_centre = 0;

    for (int c = 0; c < nCentres; c++) {
      // Start from the query point, then overwrite the active dimensions
      if (hasE) t_pt_E = q_pt_E;
      if (hasS) t_pt_S = q_pt_S;

      for (int d_local = 0; d_local < nDims; d_local++) {
        int d_global = dim_j[d_local];
        double tval = tess_j[c + d_local * nCentres];
        if (metric[d_global] == 0) {
          t_pt_E[coind[d_global]] = tval;
        } else {
          t_pt_S[coind[d_global]] = tval;
        }
      }

      double dval = 0.0;
      if (hasE) dval += euclidean_distance(q_pt_E, t_pt_E);
      if (hasS) dval += spherical_distance(q_pt_S, t_pt_S);

      if (dval < best_dist) {
        best_dist = dval;
        best_centre = c;
      }
    }
    result[obs] = best_centre;
  }
}

// log(dbinom(x, n, p)) – implemented without Rmath to avoid potential naming conflicts
static double log_dbinom(int x, int n, double p) {
  if (x < 0 || x > n) return -1e300;
  if (p <= 0.0) return (x == 0) ? 0.0 : -1e300;
  if (p >= 1.0) return (x == n) ? 0.0 : -1e300;
  return lgamma(n + 1.0) - lgamma(x + 1.0) - lgamma(n - x + 1.0)
         + x * log(p) + (n - x) * log(1.0 - p);
}

// log(dpois(x, lambda))
static double log_dpois(int x, double lambda) {
  if (x < 0 || lambda <= 0.0) return -1e300;
  return x * log(lambda) - lambda - lgamma(x + 1.0);
}

// Log acceptance probability for the Metropolis-Hastings step.
// d and cStar refer to the NEW tessellation (after the proposed modification).
// modification is one of MOD_*.
static double acceptance_prob_internal(
  const std::vector<double>& R_ijOld, const std::vector<int>& n_ijOld,
  const std::vector<double>& R_ijNew, const std::vector<int>& n_ijNew,
  int cStar, int d,
  double SigmaSquared, int modification,
  double SigmaSquaredMu, double Omega, double LambdaRate, int NumCovariates
) {
  const double prob_eps = 1e-10;
  double prob = Omega / NumCovariates;
  if (prob < 0.0) prob = 0.0;
  if (prob > 1.0 - prob_eps) prob = 1.0 - prob_eps;

  // Log-likelihood ratio
  double logLR = 0.0;
  for (int k = 0; k < (int)n_ijOld.size(); k++) {
    logLR += 0.5 * log((double)n_ijOld[k] * SigmaSquaredMu + SigmaSquared);
  }
  for (int k = 0; k < (int)n_ijNew.size(); k++) {
    logLR -= 0.5 * log((double)n_ijNew[k] * SigmaSquaredMu + SigmaSquared);
  }
  double sum_new = 0.0, sum_old = 0.0;
  for (int k = 0; k < (int)n_ijNew.size(); k++) {
    double denom = (double)n_ijNew[k] * SigmaSquaredMu + SigmaSquared;
    sum_new += R_ijNew[k] * R_ijNew[k] / denom;
  }
  for (int k = 0; k < (int)n_ijOld.size(); k++) {
    double denom = (double)n_ijOld[k] * SigmaSquaredMu + SigmaSquared;
    sum_old += R_ijOld[k] * R_ijOld[k] / denom;
  }
  logLR += (SigmaSquaredMu / (2.0 * SigmaSquared)) * (sum_new - sum_old);

  double acc = logLR;

  if (modification == MOD_AD) {
    double tessStr = log_dbinom(d - 1, NumCovariates - 1, prob)
                   - log_dbinom(d - 2, NumCovariates - 1, prob)
                   - log((double)(NumCovariates - d + 1));
    double transR  = log((double)(NumCovariates - d + 1)) - log((double)d);
    acc += tessStr + transR;
    if (d == 1)                    acc += log(0.5);
    else if (d == NumCovariates - 1) acc += log(2.0);
  } else if (modification == MOD_RD) {
    double tessStr = log_dbinom(d - 1, NumCovariates, prob)
                   + log((double)(NumCovariates - d))
                   - log_dbinom(d, NumCovariates, prob);
    double transR  = log((double)(d + 1)) - log((double)(NumCovariates - d));
    acc += tessStr + transR;
    if (d == NumCovariates) acc += log(0.5);
    else if (d == 2)        acc += log(2.0);
  } else if (modification == MOD_AC) {
    double tessStr = log_dpois(cStar - 1, LambdaRate) - log_dpois(cStar - 2, LambdaRate);
    double transR  = -log((double)cStar);
    acc += tessStr + transR + 0.5 * log(SigmaSquared);
    if (cStar == 1) acc += log(0.5);
  } else if (modification == MOD_RC) {
    double tessStr = log_dpois(cStar - 1, LambdaRate) - log_dpois(cStar, LambdaRate);
    double transR  = log((double)(cStar + 1));
    acc += tessStr + transR - 0.5 * log(SigmaSquared);
    if (cStar == 2) acc += log(2.0);
  }
  // MOD_CHANGE and MOD_SWAP: TessStructure = 1, TransitionRatio = 1 → add 0

  return acc;
}

// Sample mu values for each cell.  Writes norm_rand()-based samples into pred_j.
static void sample_mu_internal(
  const std::vector<double>& R_ij,
  const std::vector<int>& n_ij,
  double SigmaSquaredMu, double SigmaSquared,
  std::vector<double>& pred_j
) {
  int N = (int)R_ij.size();
  pred_j.resize(N);
  for (int k = 0; k < N; k++) {
    double denom    = SigmaSquaredMu * (double)n_ij[k] + SigmaSquared;
    double mean_k   = SigmaSquaredMu * R_ij[k] / denom;
    double sd_k     = sqrt(SigmaSquared * SigmaSquaredMu / denom);
    pred_j[k]       = mean_k + norm_rand() * sd_k;
  }
}

// Propose a new tessellation.  Mirrors propose_tessellation_cpp exactly (same RNG calls)
// but operates on C++ vectors.  modification is set to one of MOD_*.
static void propose_tess_internal(
  const std::vector<double>& tess_j, int nCentres, int nDims,
  const std::vector<int>& dim_j,
  const double* sd_ptr, const double* mu_ptr,
  int numCovariates,
  const std::vector<int>& metric,
  const std::vector<int>& sphere_index,
  std::vector<double>& new_tess,
  std::vector<int>& new_dim,
  int& modification
) {
  new_dim  = dim_j;
  new_tess = tess_j;
  modification = MOD_CHANGE;

  double p = unif_rand();

  // Add Dimension (AD)
  if ((p < 0.2 && nDims != numCovariates) ||
      (nDims == 1 && nDims != numCovariates && p < 0.4)) {
    modification = MOD_AD;
    int new_d_0; // 0-based global dim index
    do {
      new_d_0 = (int)(unif_rand() * numCovariates); // 0-based: [0, numCovariates)
    } while (in_vector(new_d_0, new_dim));
    new_dim.push_back(new_d_0); // store 0-based

    int new_nd = nDims + 1;
    std::vector<double> tmp(nCentres * new_nd);
    for (int r = 0; r < nCentres; r++) {
      for (int c = 0; c < nDims; c++) {
        tmp[r + c * nCentres] = new_tess[r + c * nCentres];
      }
      double nv = mu_ptr[new_d_0] + norm_rand() * sd_ptr[new_d_0];
      if (metric[new_d_0] == 1) {
        if (!sphere_index.empty() && new_d_0 == sphere_index.back()) {
          nv = period_shift(nv, M_PI);
        }
      }
      tmp[r + nDims * nCentres] = nv;
    }
    new_tess = tmp;

  // Remove Dimension (RD)
  } else if (p < 0.4 && nDims > 1) {
    modification = MOD_RD;
    int rem = (int)(unif_rand() * nDims);
    new_dim.erase(new_dim.begin() + rem);
    std::vector<double> tmp;
    tmp.reserve(nCentres * (nDims - 1));
    int cur_col = 0;
    for (int c = 0; c < nDims; c++) {
      if (c != rem) {
        for (int r = 0; r < nCentres; r++) {
          tmp.push_back(new_tess[r + c * nCentres]);
        }
        cur_col++;
      }
    }
    (void)cur_col;
    new_tess = tmp;

  // Add Centre (AC)
  } else if (p < 0.6 || (p < 0.8 && nCentres == 1)) {
    modification = MOD_AC;
    // Insert a new row at the end of each column, using the global dimension's proposal params
    for (int i = 0; i < nDims; i++) {
      int gd = new_dim[i]; // 0-based global dimension index
      double nv = mu_ptr[gd] + norm_rand() * sd_ptr[gd];
      if (metric[gd] == 1) {
        if (!sphere_index.empty() && gd == (int)sphere_index.back()) {
          nv = period_shift(nv, M_PI);
        }
      }
      new_tess.insert(new_tess.begin() + (i * (nCentres + 1)) + nCentres, nv);
    }

  // Remove Centre (RC)
  } else if (p < 0.8 && nCentres > 1) {
    modification = MOD_RC;
    int rem_row = (int)(unif_rand() * nCentres);
    std::vector<double> tmp;
    tmp.reserve((nCentres - 1) * nDims);
    for (int c = 0; c < nDims; c++) {
      for (int r = 0; r < nCentres; r++) {
        if (r != rem_row) tmp.push_back(new_tess[r + c * nCentres]);
      }
    }
    new_tess = tmp;

  // Change one centre (default)
  } else if (p < 0.9 || nDims == numCovariates) {
    modification = MOD_CHANGE;
    int ctr = (int)(unif_rand() * nCentres);
    for (int c = 0; c < nDims; c++) {
      int gd = new_dim[c]; // 0-based global dimension index
      double nv = mu_ptr[gd] + norm_rand() * sd_ptr[gd];
      if (metric[gd] == 1) {
        if (!sphere_index.empty() && gd == (int)sphere_index.back()) {
          nv = period_shift(nv, M_PI);
        }
      }
      new_tess[ctr + c * nCentres] = nv;
    }

  // Swap one dimension
  } else {
    modification = MOD_SWAP;
    int dim_chg = (int)(unif_rand() * nDims);
    int new_d_0; // 0-based global dim index for the replacement
    do {
      new_d_0 = (int)(unif_rand() * numCovariates); // 0-based
    } while (in_vector(new_d_0, new_dim));
    new_dim[dim_chg] = new_d_0; // store 0-based
    for (int r = 0; r < nCentres; r++) {
      double nv = mu_ptr[new_d_0] + norm_rand() * sd_ptr[new_d_0];
      if (metric[new_d_0] == 1) {
        if (!sphere_index.empty() && new_d_0 == (int)sphere_index.back()) {
          nv = period_shift(nv, M_PI);
        }
      }
      new_tess[r + dim_chg * nCentres] = nv;
    }
  }
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
  
  // ---------------------------------------------------------------------------
  // 3. addi_vortes_mcmc_cpp
  // ---------------------------------------------------------------------------
  // Runs the full MCMC back-fitting loop for the AddiVortes model in a single
  // C++ call, eliminating repeated R <-> C++ round-trips.
  //
  // Arguments (18):
  //   xScaled_sexp       – n x p covariate matrix (column-major double)
  //   yScaled_sexp       – n-vector of scaled responses
  //   m_sexp             – number of tessellations
  //   totalMCMCIter_sexp – total MCMC iterations
  //   mcmcBurnIn_sexp    – burn-in count
  //   thinning_sexp      – thinning rate for posterior storage
  //   nu_sexp            – prior degrees of freedom
  //   lambda_sexp        – prior scale for sigma^2
  //   SigmaSquaredMu_sexp– prior variance for cell means
  //   Omega_sexp         – dimension inclusion prior rate
  //   LambdaRate_sexp    – Poisson rate for number of centres
  //   metric_sexp        – p-integer vector (0=Euclidean, 1=Spherical)
  //   sd_prop_sexp       – p-vector of proposal SDs
  //   mus_sexp           – p-vector of proposal means
  //   init_tess_sexp     – list of m initial tessellation matrices
  //   init_dim_sexp      – list of m initial dimension integer vectors (1-based)
  //   binaryCols_sexp    – integer vector of 1-based binary column indices (or 0-length)
  //   catScaling_sexp    – scalar: upper clamp for binary tessellation centres
  //
  // Returns a named list:
  //   posteriorTess      – list[numPosteriorSamples] of list[m] of matrices
  //   posteriorDim       – list[numPosteriorSamples] of list[m] of int vectors
  //   posteriorPred      – list[numPosteriorSamples] of list[m] of real vectors
  //   posteriorSigma     – real vector[numPosteriorSamples]
  //   predictionMatrix   – n x numPosteriorSamples matrix
  SEXP addi_vortes_mcmc_cpp(
    SEXP xScaled_sexp, SEXP yScaled_sexp,
    SEXP m_sexp, SEXP totalMCMCIter_sexp,
    SEXP mcmcBurnIn_sexp, SEXP thinning_sexp,
    SEXP nu_sexp, SEXP lambda_sexp,
    SEXP SigmaSquaredMu_sexp,
    SEXP Omega_sexp, SEXP LambdaRate_sexp,
    SEXP metric_sexp, SEXP sd_prop_sexp, SEXP mus_sexp,
    SEXP init_tess_sexp, SEXP init_dim_sexp,
    SEXP binaryCols_sexp, SEXP catScaling_sexp
  ) {
    // ---- Unpack scalars ----
    int n           = Rf_nrows(xScaled_sexp);
    int p           = Rf_ncols(xScaled_sexp);
    int m           = INTEGER(m_sexp)[0];
    int totalIter   = INTEGER(totalMCMCIter_sexp)[0];
    int burnIn      = INTEGER(mcmcBurnIn_sexp)[0];
    int thinning    = INTEGER(thinning_sexp)[0];
    double nu       = REAL(nu_sexp)[0];
    double lambda   = REAL(lambda_sexp)[0];
    double sigMu    = REAL(SigmaSquaredMu_sexp)[0];
    double Omega    = REAL(Omega_sexp)[0];
    double lambdaR  = REAL(LambdaRate_sexp)[0];
    double catScale = REAL(catScaling_sexp)[0];

    // ---- Unpack arrays ----
    const double* xScaled   = REAL(xScaled_sexp);
    const double* yScaled   = REAL(yScaled_sexp);
    const double* sd_ptr    = REAL(sd_prop_sexp);
    const double* mu_ptr    = REAL(mus_sexp);

    std::vector<int> metric(INTEGER(metric_sexp),
                            INTEGER(metric_sexp) + Rf_length(metric_sexp));

    // Binary columns (0-based)
    int nBin = Rf_length(binaryCols_sexp);
    std::set<int> binaryColsSet;
    if (nBin > 0) {
      int* bcp = INTEGER(binaryCols_sexp);
      for (int i = 0; i < nBin; i++) binaryColsSet.insert(bcp[i] - 1);
    }

    // ---- Posterior sample count ----
    int numPS = 0;
    if (totalIter > burnIn) {
      numPS = (totalIter - burnIn) / thinning;
    }

    // ---- Precompute metric helpers (mirrors knnx_index_cpp) ----
    int nE = 0, nS = 0;
    std::vector<int> coind(p);
    for (int d = 0; d < p; d++) {
      if (metric[d] == 0) { coind[d] = nE++; }
      else                  { coind[d] = nS++; }
    }

    // Sphere index (0-based global dim indices where metric==1)
    std::vector<int> sphere_index;
    for (int d = 0; d < p; d++) {
      if (metric[d] == 1) sphere_index.push_back(d);
    }

    // ---- Unpack initial tessellations ----
    // tess_data[j]: flat column-major matrix (nCentres_j x nDims_j)
    // dim_data[j]:  0-based active global column indices
    std::vector<std::vector<double>> tess_data(m);
    std::vector<std::vector<int>>    dim_data(m);
    std::vector<int>                 nCentres_v(m), nDims_v(m);

    for (int j = 0; j < m; j++) {
      SEXP tj   = VECTOR_ELT(init_tess_sexp, j);
      SEXP dj   = VECTOR_ELT(init_dim_sexp,  j);
      int  nr   = Rf_nrows(tj);
      int  nd   = Rf_length(dj);
      nCentres_v[j] = nr;
      nDims_v[j]    = nd;
      tess_data[j].assign(REAL(tj), REAL(tj) + nr * nd);
      // Convert from 1-based R to 0-based C++
      int* dp = INTEGER(dj);
      dim_data[j].resize(nd);
      for (int k = 0; k < nd; k++) dim_data[j][k] = dp[k] - 1;
    }

    // ---- Compute mean_y (used for initial pred and sumOfAllTess) ----
    double sum_y = 0.0;
    for (int i = 0; i < n; i++) sum_y += yScaled[i];
    double mean_y = (n > 0) ? sum_y / n : 0.0;

    // ---- Initialise state ----
    // pred_data[j]: predicted value per centre (initially mean_y/m for each)
    // idx_data[j]:  0-based centre assignment per obs (initially 0, since 1 centre)
    // sumOfAllTess: n-vector, initially mean_y everywhere
    std::vector<std::vector<double>> pred_data(m);
    std::vector<std::vector<int>>    idx_data(m);
    for (int j = 0; j < m; j++) {
      pred_data[j].assign(nCentres_v[j], mean_y / m);
      idx_data[j].assign(n, 0);
    }
    std::vector<double> sumOfAllTess(n, mean_y);
    std::vector<double> last_tess_pred(n, 0.0);

    // ---- Posterior storage (C++ side, converted to SEXP after the loop) ----
    // Outer index: sample index (0 .. numPS-1)
    // Inner index: tessellation j (0 .. m-1)
    std::vector<std::vector<std::vector<double>>> post_tess(numPS, std::vector<std::vector<double>>(m));
    std::vector<std::vector<std::vector<int>>>    post_dim( numPS, std::vector<std::vector<int>>(m));
    std::vector<std::vector<std::vector<double>>> post_pred(numPS, std::vector<std::vector<double>>(m));
    std::vector<std::vector<int>>                 post_nrows(numPS, std::vector<int>(m));
    std::vector<double>                           post_sigma(numPS, 0.0);
    std::vector<double>                           predMatrix((size_t)n * numPS, 0.0);

    int storeIdx = 0; // next slot in posterior arrays

    // ---- Temporaries reused across iterations ----
    std::vector<int>    new_dim_prop;
    std::vector<double> new_tess_prop;
    std::vector<int>    new_idx(n);
    std::vector<double> R_ijOld, R_ijNew;
    std::vector<int>    n_ijOld, n_ijNew;
    std::vector<double> new_pred_j;
    std::vector<double> R_j(n);

    // ---- MCMC loop ----
    GetRNGstate();

    for (int i = 0; i < totalIter; i++) {

      // -- Sample sigma^2 from InvGamma --
      double ss_resid = 0.0;
      for (int obs = 0; obs < n; obs++) {
        double diff = yScaled[obs] - sumOfAllTess[obs];
        ss_resid += diff * diff;
      }
      double shape    = (nu + n) / 2.0;
      double rate_val = (nu * lambda + ss_resid) / 2.0;
      double sigSq    = 1.0 / rgamma(shape, 1.0 / rate_val);

      // -- Back-fitting over tessellations --
      for (int j = 0; j < m; j++) {

        int nC = nCentres_v[j];
        int nD = nDims_v[j];

        // Propose new tessellation
        int mod_type;
        propose_tess_internal(
          tess_data[j], nC, nD, dim_data[j],
          sd_ptr, mu_ptr, p, metric, sphere_index,
          new_tess_prop, new_dim_prop, mod_type
        );

        // Clamp binary-column coordinates to [0, catScaling]
        if (!binaryColsSet.empty()) {
          int new_nC  = (int)new_dim_prop.size() > 0 ?
                        (int)new_tess_prop.size() / (int)new_dim_prop.size() : 0;
          for (int d_local = 0; d_local < (int)new_dim_prop.size(); d_local++) {
            int d_global = new_dim_prop[d_local];
            if (binaryColsSet.count(d_global)) {
              for (int r = 0; r < new_nC; r++) {
                double& v = new_tess_prop[r + d_local * new_nC];
                if (v < 0.0)        v = 0.0;
                if (v > catScale)   v = catScale;
              }
            }
          }
        }

        int new_nC = (int)new_dim_prop.size() > 0 ?
                     (int)new_tess_prop.size() / (int)new_dim_prop.size() : 0;
        int new_nD = (int)new_dim_prop.size();

        // Compute new cell assignments
        new_idx.resize(n);
        knn1_internal(xScaled, n, p,
                      new_tess_prop.data(), new_nC, new_nD, new_dim_prop,
                      metric, nE, nS, coind, new_idx);

        // Calculate partial residuals R_j, then aggregate per cell
        // (mirrors calculateResiduals in R)
        if (j == 0) {
          for (int obs = 0; obs < n; obs++) {
            sumOfAllTess[obs] -= pred_data[j][idx_data[j][obs]];
          }
        } else {
          for (int obs = 0; obs < n; obs++) {
            sumOfAllTess[obs] += last_tess_pred[obs] - pred_data[j][idx_data[j][obs]];
          }
        }
        for (int obs = 0; obs < n; obs++) R_j[obs] = yScaled[obs] - sumOfAllTess[obs];

        // Aggregate residuals for old tessellation
        int old_nLevels = nC;
        R_ijOld.assign(old_nLevels, 0.0);
        n_ijOld.assign(old_nLevels, 0);
        for (int obs = 0; obs < n; obs++) {
          int k = idx_data[j][obs];
          R_ijOld[k] += R_j[obs];
          n_ijOld[k]++;
        }

        // Aggregate residuals for proposed tessellation
        R_ijNew.assign(new_nC, 0.0);
        n_ijNew.assign(new_nC, 0);
        for (int obs = 0; obs < n; obs++) {
          int k = new_idx[obs];
          R_ijNew[k] += R_j[obs];
          n_ijNew[k]++;
        }

        // Check for empty cells in proposal
        bool any_empty = false;
        for (int k = 0; k < new_nC; k++) {
          if (n_ijNew[k] == 0) { any_empty = true; break; }
        }

        bool accepted = false;
        if (!any_empty) {
          double logAcc = acceptance_prob_internal(
            R_ijOld, n_ijOld, R_ijNew, n_ijNew,
            new_nC, new_nD, sigSq, mod_type,
            sigMu, Omega, lambdaR, p
          );
          if (log(unif_rand()) < logAcc) {
            accepted = true;
          }
        }

        if (accepted) {
          tess_data[j]  = new_tess_prop;
          dim_data[j]   = new_dim_prop;
          nCentres_v[j] = new_nC;
          nDims_v[j]    = new_nD;
          idx_data[j]   = new_idx;
          sample_mu_internal(R_ijNew, n_ijNew, sigMu, sigSq, pred_data[j]);
          for (int obs = 0; obs < n; obs++) {
            last_tess_pred[obs] = pred_data[j][new_idx[obs]];
          }
        } else {
          sample_mu_internal(R_ijOld, n_ijOld, sigMu, sigSq, pred_data[j]);
          for (int obs = 0; obs < n; obs++) {
            last_tess_pred[obs] = pred_data[j][idx_data[j][obs]];
          }
        }

        // After the last tessellation, add its prediction to sumOfAllTess
        if (j == m - 1) {
          for (int obs = 0; obs < n; obs++) {
            sumOfAllTess[obs] += last_tess_pred[obs];
          }
        }
      } // end for j

      // Store posterior sample (after burn-in, respecting thinning)
      if (i >= burnIn && ((i - burnIn) % thinning == 0)) {
        int s = storeIdx++;
        for (int j = 0; j < m; j++) {
          post_tess[s][j]  = tess_data[j];
          post_dim[s][j]   = dim_data[j];  // 0-based; converted to 1-based when building SEXP
          post_pred[s][j]  = pred_data[j];
          post_nrows[s][j] = nCentres_v[j];
        }
        post_sigma[s] = sigSq;
        for (int obs = 0; obs < n; obs++) {
          predMatrix[(size_t)obs + (size_t)s * n] = sumOfAllTess[obs];
        }
      }
    } // end for i (MCMC)

    PutRNGstate();

    // ============================================================
    // Build the R return value
    // ============================================================
    // We protect the 5 top-level output objects + the result list + names
    SEXP res_postTess, res_postDim, res_postPred, res_postSigma, res_predMatrix;
    SEXP result_list, list_names;

    PROTECT(res_postTess  = Rf_allocVector(VECSXP, numPS));
    PROTECT(res_postDim   = Rf_allocVector(VECSXP, numPS));
    PROTECT(res_postPred  = Rf_allocVector(VECSXP, numPS));
    PROTECT(res_postSigma = Rf_allocVector(REALSXP, numPS));
    PROTECT(res_predMatrix = Rf_allocMatrix(REALSXP, n, numPS));

    // Fill posteriorSigma and predictionMatrix
    if (numPS > 0) {
      memcpy(REAL(res_postSigma), post_sigma.data(), numPS * sizeof(double));
      memcpy(REAL(res_predMatrix), predMatrix.data(), (size_t)n * numPS * sizeof(double));
    }

    // Fill posteriorTess / posteriorDim / posteriorPred
    for (int s = 0; s < numPS; s++) {
      SEXP inner_tess, inner_dim, inner_pred;

      // Protect each inner list so GC can't collect them while we fill them
      PROTECT(inner_tess = Rf_allocVector(VECSXP, m));
      PROTECT(inner_dim  = Rf_allocVector(VECSXP, m));
      PROTECT(inner_pred = Rf_allocVector(VECSXP, m));

      for (int j = 0; j < m; j++) {
        int nR = post_nrows[s][j];
        int nD = (int)post_dim[s][j].size();

        // tessellation matrix (nR x nD)
        SEXP mat_j, dim_j_sexp_r, pred_j_sexp;
        PROTECT(mat_j = Rf_allocMatrix(REALSXP, nR, nD));
        memcpy(REAL(mat_j), post_tess[s][j].data(), (size_t)nR * nD * sizeof(double));
        SET_VECTOR_ELT(inner_tess, j, mat_j);
        UNPROTECT(1); // mat_j now owned by inner_tess

        // dim vector (1-based for R)
        PROTECT(dim_j_sexp_r = Rf_allocVector(INTSXP, nD));
        for (int d = 0; d < nD; d++) {
          INTEGER(dim_j_sexp_r)[d] = post_dim[s][j][d] + 1; // 0-based → 1-based
        }
        SET_VECTOR_ELT(inner_dim, j, dim_j_sexp_r);
        UNPROTECT(1); // dim_j_sexp_r now owned by inner_dim

        // pred vector
        PROTECT(pred_j_sexp = Rf_allocVector(REALSXP, nR));
        memcpy(REAL(pred_j_sexp), post_pred[s][j].data(), nR * sizeof(double));
        SET_VECTOR_ELT(inner_pred, j, pred_j_sexp);
        UNPROTECT(1); // pred_j_sexp now owned by inner_pred
      }

      SET_VECTOR_ELT(res_postTess, s, inner_tess);
      SET_VECTOR_ELT(res_postDim,  s, inner_dim);
      SET_VECTOR_ELT(res_postPred, s, inner_pred);
      UNPROTECT(3); // inner_tess, inner_dim, inner_pred now owned by res_post*
    }

    // Assemble result list
    PROTECT(result_list = Rf_allocVector(VECSXP, 5));
    SET_VECTOR_ELT(result_list, 0, res_postTess);
    SET_VECTOR_ELT(result_list, 1, res_postDim);
    SET_VECTOR_ELT(result_list, 2, res_postPred);
    SET_VECTOR_ELT(result_list, 3, res_postSigma);
    SET_VECTOR_ELT(result_list, 4, res_predMatrix);

    PROTECT(list_names = Rf_allocVector(STRSXP, 5));
    SET_STRING_ELT(list_names, 0, Rf_mkChar("posteriorTess"));
    SET_STRING_ELT(list_names, 1, Rf_mkChar("posteriorDim"));
    SET_STRING_ELT(list_names, 2, Rf_mkChar("posteriorPred"));
    SET_STRING_ELT(list_names, 3, Rf_mkChar("posteriorSigma"));
    SET_STRING_ELT(list_names, 4, Rf_mkChar("predictionMatrix"));
    Rf_setAttrib(result_list, R_NamesSymbol, list_names);

    UNPROTECT(7); // res_postTess, res_postDim, res_postPred, res_postSigma,
                  // res_predMatrix, result_list, list_names
    return result_list;
  }

} // extern "C"