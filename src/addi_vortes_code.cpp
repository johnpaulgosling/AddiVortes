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
double euclidean_distance(std::vector<double>& p1, std::vector<double>& p2) {
  if (p1.size() != p2.size()) {
    Rf_error("Points have incompatible dimensions.");
  }
  double dist = 0.0;
  for (int i = 0; i < p1.size(); ++i) {
    dist += pow(p1[i]-p2[i], 2);
  }
  return dist;
}

// Computes great circle distance on an n-sphere.
// It assumes that the last dimension is the 'azimuthal' distance; i.e. the one with range [-pi, pi].
// Optimisation not great: needs some work with the trigonometric functions.
// The issue is that one needs all the covariates to calculate: you can't just cull the covariates not in the tessellation.
double spherical_distance(std::vector<double>& p1, std::vector<double>& p2) {
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
    std::vector<double> t_pt_E(how_many_E);
    std::vector<double> q_pt_E(how_many_E);
    int how_many_S = n_elem(1, metric);
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
    double dval;
    //std::vector<double> q_pt(query_cols);
    //std::vector<double> t_pt(tess_cols);
    
    // --- Main Logic: For each query point, find k nearest neighbors ---
    for (int q = 0; q < query_rows; ++q) {
      
      // Calculate distances to all tessellation points
      std::vector<std::pair<double, int>> distances(tess_rows);
      
      /// Getting ref errors here when building the q_pt_S and t_pt_S
      for (int t = 0; t < tess_rows; ++t) {
        dval = 0.0;
        for (int d = 0; d < query_cols; ++d) {
          if (metric[d] == 0) {
            q_pt_E[coind[d]] = p_query[q + d * query_rows];
          }
          if (metric[d] == 1) {
            q_pt_S[coind[d]] = p_query[q + d * query_rows];
          }
        }
        for (int d = 0; d < tess_cols; ++d) {
          if (in_vector(d+1, dim_p_temp)) {
            if (metric[d] == 0) {
              t_pt_E[coind[d]] = p_tess[t + d * tess_rows];
            }
            if (metric[d] == 1) {
              t_pt_S[coind[d]] = p_tess[t + d * tess_rows];
            }
          }
          else {
            if (metric[d] == 0) {
              t_pt_E[coind[d]] = p_query[q + d * query_rows];
            }
            if (metric[d] == 1) {
              t_pt_S[coind[d]] = p_query[q + d * query_rows];
            }
          }
        }
        if (in_vector(1, metric)) {
          dval += spherical_distance(q_pt_S, t_pt_S);
        }
        if (in_vector(0, metric)) {
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