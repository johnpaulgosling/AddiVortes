#include <R.h>
#include <Rinternals.h>
#include <R_ext/Random.h> // For unif_rand() and norm_rand()
#include <Rmath.h>        // For Rf_dbinom, Rf_dpois

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

// ---------------------------------------------------------------------------
// 3. sample_mu_values_cpp
// ---------------------------------------------------------------------------
// C++ equivalent of the R function `sampleMuValues`
std::vector<double> sample_mu_values_cpp(int n_centres,
                                         const std::vector<double>& R_ij,
                                         const std::vector<int>& n_ij,
                                         double sigma_sq_mu,
                                         double sigma_sq) {
    std::vector<double> mu_values(n_centres);
    for (int i = 0; i < n_centres; ++i) {
        double denominator = sigma_sq_mu * n_ij[i] + sigma_sq;
        double mean = (sigma_sq_mu * R_ij[i]) / denominator;
        double sd = sqrt((sigma_sq * sigma_sq_mu) / denominator);
        mu_values[i] = norm_rand() * sd + mean;
    }
    return mu_values;
}


// ---------------------------------------------------------------------------
// 4. acceptance_prob_cpp
// ---------------------------------------------------------------------------
// C++ equivalent of the R function `acceptanceProbability`
double acceptance_prob_cpp(const std::vector<double>& R_ijOld, const std::vector<int>& n_ijOld,
                             const std::vector<double>& R_ijNew, const std::vector<int>& n_ijNew,
                             int cStar, int d_star_len,
                             const std::string& modification,
                             double sigma_sq, double sigma_sq_mu,
                             double omega, double lambda_rate, int num_cov) {

    double log_like_ratio = 0.0;
    double term_old = 0.0;
    double term_new = 0.0;
    double sum_old = 0.0;
    double sum_new = 0.0;

    for (size_t i = 0; i < n_ijOld.size(); ++i) {
        term_old += log(n_ijOld[i] * sigma_sq_mu + sigma_sq);
        sum_old += (R_ijOld[i] * R_ijOld[i]) / (n_ijOld[i] * sigma_sq_mu + sigma_sq);
    }
    for (size_t i = 0; i < n_ijNew.size(); ++i) {
        term_new += log(n_ijNew[i] * sigma_sq_mu + sigma_sq);
        sum_new += (R_ijNew[i] * R_ijNew[i]) / (n_ijNew[i] * sigma_sq_mu + sigma_sq);
    }
    log_like_ratio = 0.5 * (term_old - term_new) + (sigma_sq_mu / (2.0 * sigma_sq)) * (sum_new - sum_old);

    double tess_structure = 1.0;
    double trans_ratio = 1.0;
    double acceptance_prob = log_like_ratio;

    if (modification == "AD") {
        tess_structure = (Rf_dbinom(d_star_len - 2, num_cov - 1, omega / num_cov, 0)) /
                         (Rf_dbinom(d_star_len - 3, num_cov - 1, omega / num_cov, 0) * (num_cov - d_star_len + 2));
        trans_ratio = (double)(num_cov - d_star_len + 2) / (d_star_len-1);
        acceptance_prob += log(tess_structure) + log(trans_ratio);
        if (d_star_len -1 == 1) acceptance_prob += log(0.5);
        else if (d_star_len-1 == num_cov - 1) acceptance_prob += log(2.0);

    } else if (modification == "RD") {
        tess_structure = (Rf_dbinom(d_star_len - 1, num_cov, omega / num_cov, 0) * (num_cov - d_star_len)) /
                         (Rf_dbinom(d_star_len, num_cov, omega / num_cov, 0));
        trans_ratio = (double)(d_star_len + 1) / (num_cov - d_star_len);
        acceptance_prob += log(tess_structure) + log(trans_ratio);
        if (d_star_len == num_cov) acceptance_prob += log(0.5);
        else if (d_star_len == 2) acceptance_prob += log(2.0);

    } else if (modification == "AC") {
        tess_structure = Rf_dpois(cStar - 1, lambda_rate, 0) / Rf_dpois(cStar - 2, lambda_rate, 0);
        trans_ratio = 1.0 / cStar;
        acceptance_prob += log(tess_structure) + log(trans_ratio) + 0.5 * log(sigma_sq);
        if (cStar == 1) acceptance_prob += log(0.5);

    } else if (modification == "RC") {
        tess_structure = Rf_dpois(cStar - 1, lambda_rate, 0) / Rf_dpois(cStar, lambda_rate, 0);
        trans_ratio = (double)cStar + 1.0;
        acceptance_prob += log(tess_structure) + log(trans_ratio) - 0.5 * log(sigma_sq);
        if (cStar == 2) acceptance_prob += log(2.0);
    }
    // For "Change" and "Swap", ratios are 1, so we just use the log_like_ratio.

    return acceptance_prob;
}


// ---------------------------------------------------------------------------
// 5. mcmc_step_j_cpp (THE MAIN NEW FUNCTION)
// ---------------------------------------------------------------------------
extern "C" {
SEXP mcmc_step_j_cpp(SEXP y_scaled_sexp, SEXP x_scaled_sexp, SEXP j_sexp,
                       SEXP sum_of_all_tess_sexp,
                       SEXP pred_list_sexp, SEXP tess_list_sexp, SEXP dim_list_sexp,
                       SEXP last_tess_pred_sexp, SEXP indexes_j_sexp,
                       SEXP sigma_sq_sexp, SEXP sigma_sq_mu_sexp,
                       SEXP sd_sexp, SEXP num_cov_sexp,
                       SEXP omega_sexp, SEXP lambda_rate_sexp) {

    // --- 1. Unpack Arguments ---
    double* p_y = REAL(y_scaled_sexp);
    double* p_x = REAL(x_scaled_sexp);
    int j = INTEGER(j_sexp)[0] - 1; // Convert to 0-based index
    int n_obs = Rf_length(y_scaled_sexp);
    int num_cov = INTEGER(num_cov_sexp)[0];

    // Get current state for tessellation j
    SEXP tess_j = VECTOR_ELT(tess_list_sexp, j);
    SEXP dim_j = VECTOR_ELT(dim_list_sexp, j);
    SEXP pred_j = VECTOR_ELT(pred_list_sexp, j);

    // Other scalar parameters
    double sigma_sq = REAL(sigma_sq_sexp)[0];
    double sigma_sq_mu = REAL(sigma_sq_mu_sexp)[0];
    double sd_val = REAL(sd_sexp)[0];
    double omega = REAL(omega_sexp)[0];
    double lambda_rate = REAL(lambda_rate_sexp)[0];

    // --- 2. Calculate Partial Residuals ---
    std::vector<double> sum_of_all_tess(REAL(sum_of_all_tess_sexp), REAL(sum_of_all_tess_sexp) + n_obs);
    std::vector<double> current_tess_pred(n_obs);
    int* p_indexes_j = INTEGER(indexes_j_sexp);
    double* p_pred_j = REAL(pred_j);

    for(int i = 0; i < n_obs; ++i) {
        current_tess_pred[i] = p_pred_j[p_indexes_j[i] - 1];
    }
    
    if (j == 0) { // First tessellation (j was decremented)
        for(int i=0; i<n_obs; ++i) sum_of_all_tess[i] -= current_tess_pred[i];
    } else {
        double* p_last_tess_pred = REAL(last_tess_pred_sexp);
        for(int i=0; i<n_obs; ++i) sum_of_all_tess[i] += p_last_tess_pred[i] - current_tess_pred[i];
    }

    std::vector<double> R_j(n_obs);
    for(int i=0; i<n_obs; ++i) R_j[i] = p_y[i] - sum_of_all_tess[i];

    // --- 3. Propose New Tessellation ---
    SEXP proposal_list;
    PROTECT(proposal_list = propose_tessellation_cpp(tess_j, dim_j, sd_sexp, num_cov_sexp));
    SEXP tess_j_star_sexp = VECTOR_ELT(proposal_list, 0);
    SEXP dim_j_star_sexp = VECTOR_ELT(proposal_list, 1);
    std::string modification = CHAR(STRING_ELT(VECTOR_ELT(proposal_list, 2), 0));

    // --- 4. Find New Indices (k-NN) ---
    int num_centres_new = Rf_nrows(tess_j_star_sexp);
    SEXP indexes_star_sexp;
    if (num_centres_new == 1) {
        PROTECT(indexes_star_sexp = Rf_allocVector(INTSXP, n_obs));
        int* p_indexes_star = INTEGER(indexes_star_sexp);
        for(int i=0; i<n_obs; ++i) p_indexes_star[i] = 1;
    } else {
        int d_star_len = Rf_length(dim_j_star_sexp);
        SEXP query_matrix;
        PROTECT(query_matrix = Rf_allocMatrix(REALSXP, n_obs, d_star_len));
        double* p_query = REAL(query_matrix);
        int* p_dim_star = INTEGER(dim_j_star_sexp);
        for(int c=0; c<d_star_len; ++c) {
            for(int r=0; r<n_obs; ++r) {
                p_query[r + c*n_obs] = p_x[r + (p_dim_star[c]-1)*n_obs];
            }
        }
        PROTECT(indexes_star_sexp = knnx_index_cpp(tess_j_star_sexp, query_matrix, Rf_ScalarInteger(1)));
        UNPROTECT(1); // query_matrix
    }

    // --- 5. Calculate Residual Summaries ---
    int num_centres_old = Rf_nrows(tess_j);
    std::vector<double> R_ijOld(num_centres_old, 0.0), R_ijNew(num_centres_new, 0.0);
    std::vector<int> n_ijOld(num_centres_old, 0), n_ijNew(num_centres_new, 0);

    for(int i = 0; i < n_obs; ++i) {
        R_ijOld[INTEGER(indexes_j_sexp)[i] - 1] += R_j[i];
        n_ijOld[INTEGER(indexes_j_sexp)[i] - 1]++;
        R_ijNew[INTEGER(indexes_star_sexp)[i] - 1] += R_j[i];
        n_ijNew[INTEGER(indexes_star_sexp)[i] - 1]++;
    }

    // --- 6. Acceptance Step ---
    bool accepted = false;
    bool any_empty_cells = false;
    for(int count : n_ijNew) {
        if (count == 0) {
            any_empty_cells = true;
            break;
        }
    }

    if (!any_empty_cells) {
        double log_acceptance_prob = acceptance_prob_cpp(R_ijOld, n_ijOld, R_ijNew, n_ijNew,
                                                         num_centres_new, Rf_length(dim_j_star_sexp),
                                                         modification, sigma_sq, sigma_sq_mu,
                                                         omega, lambda_rate, num_cov);
        if (log(unif_rand()) < log_acceptance_prob) {
            accepted = true;
        }
    }

    // --- 7. Sample New Mu Values & Update State ---
    SEXP final_tess, final_dim, final_pred, final_indexes;
    if(accepted) {
        std::vector<double> pred_new_vals = sample_mu_values_cpp(num_centres_new, R_ijNew, n_ijNew, sigma_sq_mu, sigma_sq);
        PROTECT(final_pred = Rf_allocVector(REALSXP, num_centres_new));
        memcpy(REAL(final_pred), pred_new_vals.data(), num_centres_new * sizeof(double));
        final_tess = tess_j_star_sexp;
        final_dim = dim_j_star_sexp;
        final_indexes = indexes_star_sexp;
    } else {
        std::vector<double> pred_old_vals = sample_mu_values_cpp(num_centres_old, R_ijOld, n_ijOld, sigma_sq_mu, sigma_sq);
        PROTECT(final_pred = Rf_allocVector(REALSXP, num_centres_old));
        memcpy(REAL(final_pred), pred_old_vals.data(), num_centres_old * sizeof(double));
        final_tess = tess_j;
        final_dim = dim_j;
        final_indexes = indexes_j_sexp;
    }

    // --- 8. Calculate lastTessPred & update sumOfAllTess for the final time ---
    SEXP last_tess_pred_new;
    PROTECT(last_tess_pred_new = Rf_allocVector(REALSXP, n_obs));
    double* p_last_tess_pred_new = REAL(last_tess_pred_new);
    double* p_final_pred = REAL(final_pred);
    int* p_final_indexes = INTEGER(final_indexes);

    for(int i=0; i<n_obs; ++i) {
        p_last_tess_pred_new[i] = p_final_pred[p_final_indexes[i]-1];
    }
    
    // This is only needed for the *final* tessellation of the loop (j=m-1)
    for(int i=0; i<n_obs; ++i) {
        sum_of_all_tess[i] += p_last_tess_pred_new[i];
    }

    // --- 9. Package and Return Results ---
    SEXP result_list, list_names;
    PROTECT(result_list = Rf_allocVector(VECSXP, 6));

    // Updated sum of all tessellations
    SEXP sum_tess_out = Rf_allocVector(REALSXP, n_obs);
    memcpy(REAL(sum_tess_out), sum_of_all_tess.data(), n_obs * sizeof(double));
    SET_VECTOR_ELT(result_list, 0, sum_tess_out);

    // Final state of tess, dim, pred, indexes for j
    SET_VECTOR_ELT(result_list, 1, final_tess);
    SET_VECTOR_ELT(result_list, 2, final_dim);
    SET_VECTOR_ELT(result_list, 3, final_pred);
    SET_VECTOR_ELT(result_list, 4, final_indexes);
    SET_VECTOR_ELT(result_list, 5, last_tess_pred_new);
    
    PROTECT(list_names = Rf_allocVector(STRSXP, 6));
    SET_STRING_ELT(list_names, 0, Rf_mkChar("sumOfAllTess"));
    SET_STRING_ELT(list_names, 1, Rf_mkChar("tess_j"));
    SET_STRING_ELT(list_names, 2, Rf_mkChar("dim_j"));
    SET_STRING_ELT(list_names, 3, Rf_mkChar("pred_j"));
    SET_STRING_ELT(list_names, 4, Rf_mkChar("indexes_j"));
    SET_STRING_ELT(list_names, 5, Rf_mkChar("lastTessPred"));
    Rf_setAttrib(result_list, R_NamesSymbol, list_names);

    UNPROTECT(6); // proposal_list, indexes_star_sexp (if created), final_pred, last_tess_pred_new, result_list, list_names
    return result_list;
}
} // extern "C"
