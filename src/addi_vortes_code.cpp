#include <R.h>
#include <Rinternals.h>
#include <R_ext/Random.h> // For unif_rand() and norm_rand()

#include <vector>    // For std::vector
#include <string>    // For std::string
#include <numeric>   // For std::accumulate
#include <algorithm> // For std::find
#include <cmath>     // For sqrt

// Helper function to check if a value is in a vector
bool in_vector(int value, const std::vector<int>& vec) {
  return std::find(vec.begin(), vec.end(), value) != vec.end();
}


extern "C" {
  
  // ---------------------------------------------------------------------------
  // 0. knnx_index_cpp
  // ---------------------------------------------------------------------------
  // This function implements k-nearest neighbors index search to replace FNN::knnx.index
  // The R wrapper function is `knnx.index`.
  SEXP knnx_index_cpp(SEXP tess_sexp, SEXP query_sexp, SEXP k_sexp) {
    
    // --- Unpack arguments ---
    double* p_tess = REAL(tess_sexp);
    double* p_query = REAL(query_sexp);
    int k = INTEGER(k_sexp)[0];
    
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
    
    // --- Main Logic: For each query point, find k nearest neighbors ---
    for (int q = 0; q < query_rows; ++q) {
      
      // Calculate distances to all tessellation points
      std::vector<std::pair<double, int>> distances(tess_rows);
      
      for (int t = 0; t < tess_rows; ++t) {
        double dist_sq = 0.0;
        for (int d = 0; d < tess_cols; ++d) {
          double diff = p_query[q + d * query_rows] - p_tess[t + d * tess_rows];
          dist_sq += diff * diff;
        }
        distances[t] = std::make_pair(dist_sq, t + 1); // +1 for R 1-based indexing
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
  SEXP propose_tessellation_cpp(SEXP tess_j_sexp, SEXP dim_j_sexp, SEXP var_sexp, SEXP num_cov_sexp) {
    
    // --- Unpack arguments ---
    double* p_tess_j = REAL(tess_j_sexp);
    int* p_dim_j = INTEGER(dim_j_sexp);
    double var = REAL(var_sexp)[0];
    int numCovariates = INTEGER(num_cov_sexp)[0];
    
    int tess_j_rows = Rf_nrows(tess_j_sexp);
    int d_j_length = Rf_length(dim_j_sexp);
    
    // Get R's random number generator state
    GetRNGstate();
    
    // --- Create C++ copies for manipulation ---
    std::vector<int> dim_j_star(p_dim_j, p_dim_j + d_j_length);
    std::vector<double> tess_j_star(p_tess_j, p_tess_j + (tess_j_rows * d_j_length));
    std::string modification = "Change";
    
    double p = unif_rand();
    
    // --- Main Logic ---
    if ((p < 0.2 && d_j_length != numCovariates) || (d_j_length == 1 && p < 0.4)) {
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
        new_tess[r + d_j_length * tess_j_rows] = norm_rand() * sqrt(var);
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
        tess_j_star.insert(tess_j_star.begin() + (i * (tess_j_rows + 1)) + tess_j_rows, norm_rand() * sqrt(var));
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
        tess_j_star[centre_to_change_idx + c * tess_j_rows] = norm_rand() * sqrt(var);
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
        tess_j_star[r + dim_to_change_idx * tess_j_rows] = norm_rand() * sqrt(var);
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