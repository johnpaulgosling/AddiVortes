#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <utility>

#define R_NO_REMAP

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Random.h>
#include <R_ext/Print.h>
#include <R_ext/Utils.h>

namespace {

bool in_vector(int value, const std::vector<int>& vec) {
  return std::find(vec.begin(), vec.end(), value) != vec.end();
}
  
  double calc_acceptance_cpp(const std::vector<double>& rIjOld, const std::vector<double>& denomOld,
                             const std::vector<double>& rIjNew, const std::vector<double>& denomNew,
                             int old_centres, int new_centres, int d_old, int d_new,
                             double sigma2, const std::string& mod, double sigma2mu, double omega, double lambda, int p,
                             double sum_notin_old, double sum_notin_new) {
    
    double prob_eps = 1e-10;
    double prob = fmin(1.0 - prob_eps, fmax(0.0, omega / p));
    
    double sum_log_old = 0.0, sum_log_new = 0.0;
    double sum_R2_old = 0.0, sum_R2_new = 0.0;
    
    const double* __restrict__ p_rIjOld = rIjOld.data();
    const double* __restrict__ p_dOld = denomOld.data();
    
    for(size_t i = 0; i < denomOld.size(); ++i) {
      sum_log_old += log(p_dOld[i]);
      sum_R2_old += (p_rIjOld[i] * p_rIjOld[i]) / p_dOld[i];
    }
    
    const double* __restrict__ p_rIjNew = rIjNew.data();
    const double* __restrict__ p_dNew = denomNew.data();
    
    for(size_t i = 0; i < denomNew.size(); ++i) {
      sum_log_new += log(p_dNew[i]);
      sum_R2_new += (p_rIjNew[i] * p_rIjNew[i]) / p_dNew[i];
    }
    
    double LogLikelihoodRatio = 0.5 * (sum_log_old - sum_log_new) +
      (sigma2mu / (2.0 * sigma2)) * (sum_R2_new - sum_R2_old);
    
    double AcceptanceProb = LogLikelihoodRatio;
    
    if (mod == "AD") {
      double tess_num = dbinom(d_new - 1, p - 1, prob, 0);
      double tess_den = dbinom(d_old - 1, p - 1, prob, 0);
      double trans = sum_notin_old / d_new;
      AcceptanceProb += log(tess_num / tess_den) + log(trans);
      if (d_new == 2) AcceptanceProb += log(0.5);
      else if (d_new == p) AcceptanceProb += log(2.0);
    } else if (mod == "RD") {
      double tess_num = dbinom(d_new - 1, p - 1, prob, 0);
      double tess_den = dbinom(d_old - 1, p - 1, prob, 0);
      double trans = (double)(d_old) / sum_notin_new;
      AcceptanceProb += log(tess_num / tess_den) + log(trans);
      if (d_old == p) AcceptanceProb += log(0.5);
      else if (d_old == 2) AcceptanceProb += log(2.0);
    } else if (mod == "AC") {
      double tess_num = dpois(new_centres - 1, lambda, 0);
      double tess_den = dpois(new_centres - 2, lambda, 0);
      double trans = 1.0 / new_centres;
      AcceptanceProb += log(tess_num / tess_den) + log(trans) + 0.5 * log(sigma2);
      if (new_centres == 1) AcceptanceProb += log(0.5);
    } else if (mod == "RC") {
      double tess_num = dpois(new_centres - 1, lambda, 0);
      double tess_den = dpois(new_centres, lambda, 0);
      double trans = new_centres + 1.0;
      AcceptanceProb += log(tess_num / tess_den) + log(trans) - 0.5 * log(sigma2);
      if (new_centres == 2) AcceptanceProb += log(2.0);
    } else if (mod == "Swap") {
      double trans = sum_notin_old / sum_notin_new;
      AcceptanceProb += log(trans);
    }
    
    return AcceptanceProb;
  }
  
  std::vector<double> sample_mu_cpp(int new_centres, const std::vector<double>& rIjNew, 
                                    const std::vector<double>& denomNew, double sigma2mu, double sigma2) {
    
    std::vector<double> new_mu(new_centres, 0.0);
    const double* __restrict__ p_rIjNew = rIjNew.data();
    const double* __restrict__ p_dNew = denomNew.data();
    double* __restrict__ p_new_mu = new_mu.data();
    
    for(int i = 0; i < new_centres; ++i) {
      double mean = (sigma2mu * p_rIjNew[i]) / p_dNew[i];
      double var = (sigma2 * sigma2mu) / p_dNew[i];
      double sd = sqrt(var);
      p_new_mu[i] = norm_rand() * sd + mean;
    }
    return new_mu;
  }
  
  double sample_sigma_squared_cpp(int n, double nu, double lambda, const double* __restrict__ p_y, const double* __restrict__ p_sum) {
    double sse = 0.0;
    for (int i = 0; i < n; ++i) {
      double diff = p_y[i] - p_sum[i];
      sse += diff * diff;
    }
    double shape = (nu + n) / 2.0;
    double rate = (nu * lambda + sse) / 2.0;
    return 1.0 / rgamma(shape, 1.0 / rate);
  }
  
}

extern "C" {
  
  SEXP knnx_index_cpp(SEXP data_sexp, SEXP tess_sexp, SEXP dim_sexp, SEXP tessstar_sexp, SEXP k_sexp, SEXP dimstar_sexp, SEXP dist_sexp, SEXP modification_sexp, SEXP row_column_modified_sexp, SEXP old_idx_sexp) {
    
    const double* __restrict__ p_data = REAL(data_sexp);
    const double* __restrict__ p_tess = REAL(tess_sexp);
    const double* __restrict__ p_tessstar = REAL(tessstar_sexp);
    const double* __restrict__ p_dist = REAL(dist_sexp);
    const int* __restrict__ p_old_idx = INTEGER(old_idx_sexp);
    
    const int* __restrict__ p_dim = INTEGER(dim_sexp);
    const int* __restrict__ p_dimstar = INTEGER(dimstar_sexp);
    int k = INTEGER(k_sexp)[0];
    
    const char* mod_str = CHAR(STRING_ELT(modification_sexp, 0));
    std::string modification(mod_str);
    int mod_idx = INTEGER(row_column_modified_sexp)[0] - 1; 
    
    int n_obs = Rf_nrows(data_sexp);
    int n_centres_old = Rf_nrows(tess_sexp);
    int n_centres_new = Rf_nrows(tessstar_sexp);
    int n_dims_new = Rf_length(dimstar_sexp);
    
    if (k <= 0 || k > n_centres_new) {
      Rf_error("k must be positive and not greater than number of reference points");
    }
    
    SEXP res_dist, res_nn, result_list, list_names;
    PROTECT(res_dist = Rf_allocMatrix(REALSXP, n_obs, n_centres_new));
    double* __restrict__ p_new_dist = REAL(res_dist);
    
    if (modification == "RC") {
      if (mod_idx > 0) {
        size_t bytes_before = mod_idx * n_obs * sizeof(double);
        std::memcpy(p_new_dist, p_dist, bytes_before);
      }
      if (mod_idx < n_centres_old - 1) {
        size_t cols_after = (n_centres_old - 1) - mod_idx;
        size_t bytes_after = cols_after * n_obs * sizeof(double);
        std::memcpy(p_new_dist + (mod_idx * n_obs), 
                    p_dist + ((mod_idx + 1) * n_obs), 
                    bytes_after);
      }
    } else if (modification == "AC") {
      size_t total_old_bytes = n_centres_old * n_obs * sizeof(double);
      std::memcpy(p_new_dist, p_dist, total_old_bytes);
      int c = n_centres_new - 1; 
      for (int i = 0; i < n_obs; ++i) {
        p_new_dist[i + c * n_obs] = 0.0;
      }
      for (int d = 0; d < n_dims_new; ++d) {
        int cov_idx = p_dimstar[d] - 1;
        int cov_offset = cov_idx * n_obs;
        double tess_val = p_tessstar[c + d * n_centres_new];
        for (int i = 0; i < n_obs; ++i) {
          double diff = p_data[i + cov_offset] - tess_val;
          p_new_dist[i + c * n_obs] += diff * diff;
        }
      }
    } else if (modification == "Change") {
      size_t total_bytes = n_centres_new * n_obs * sizeof(double);
      std::memcpy(p_new_dist, p_dist, total_bytes);
      int c = mod_idx;
      for (int i = 0; i < n_obs; ++i) {
        p_new_dist[i + c * n_obs] = 0.0;
      }
      for (int d = 0; d < n_dims_new; ++d) {
        int cov_idx = p_dimstar[d] - 1;
        int cov_offset = cov_idx * n_obs;
        double tess_val = p_tessstar[c + d * n_centres_new];
        for (int i = 0; i < n_obs; ++i) {
          double diff = p_data[i + cov_offset] - tess_val;
          p_new_dist[i + c * n_obs] += diff * diff;
        }
      }
    } else if (modification == "AD") {
      int new_dim_cov_idx = p_dimstar[mod_idx] - 1;
      int cov_offset = new_dim_cov_idx * n_obs;
      for (int c = 0; c < n_centres_new; ++c) {
        double tess_val = p_tessstar[c + mod_idx * n_centres_new];
        int c_offset = c * n_obs;
        for (int i = 0; i < n_obs; ++i) {
          double diff = p_data[i + cov_offset] - tess_val;
          p_new_dist[i + c_offset] = p_dist[i + c_offset] + (diff * diff);
        }
      }
    } else if (modification == "RD") {
      int old_dim_cov_idx = p_dim[mod_idx] - 1;
      int cov_offset = old_dim_cov_idx * n_obs;
      for (int c = 0; c < n_centres_new; ++c) {
        double tess_val = p_tess[c + mod_idx * n_centres_old];
        int c_offset = c * n_obs;
        for (int i = 0; i < n_obs; ++i) {
          double diff = p_data[i + cov_offset] - tess_val;
          double new_val = p_dist[i + c_offset] - (diff * diff);
          p_new_dist[i + c_offset] = (new_val < 0.0) ? 0.0 : new_val;
        }
      }
    } else if (modification == "Swap") {
      int old_dim_cov_idx = p_dim[mod_idx] - 1;
      int new_dim_cov_idx = p_dimstar[mod_idx] - 1;
      int old_cov_offset = old_dim_cov_idx * n_obs;
      int new_cov_offset = new_dim_cov_idx * n_obs;
      for (int c = 0; c < n_centres_new; ++c) {
        double old_tess_val = p_tess[c + mod_idx * n_centres_old];
        double new_tess_val = p_tessstar[c + mod_idx * n_centres_new];
        int c_offset = c * n_obs;
        for (int i = 0; i < n_obs; ++i) {
          double old_diff = p_data[i + old_cov_offset] - old_tess_val;
          double new_diff = p_data[i + new_cov_offset] - new_tess_val;
          double new_val = p_dist[i + c_offset] - (old_diff * old_diff) + (new_diff * new_diff);
          p_new_dist[i + c_offset] = (new_val < 0.0) ? 0.0 : new_val;
        }
      }
    } else {
      for (int c = 0; c < n_centres_new; ++c) {
        for (int i = 0; i < n_obs; ++i) {
          double dval = 0.0;
          for (int d = 0; d < n_dims_new; ++d) {
            int cov_idx = p_dimstar[d] - 1;
            double diff = p_data[i + cov_idx * n_obs] - p_tessstar[c + d * n_centres_new];
            dval += diff * diff;
          }
          p_new_dist[i + c * n_obs] = dval;
        }
      }
    }
    
    PROTECT(res_nn = Rf_allocMatrix(INTSXP, n_obs, k));
    int* __restrict__ p_result = INTEGER(res_nn);
    
    if (k == 1) {
      if (modification == "Change") {
        int mod_offset = mod_idx * n_obs;
        for (int q = 0; q < n_obs; ++q) {
          int old_best = p_old_idx[q] - 1;
          if (old_best == mod_idx) {
            double min_dist = p_new_dist[q]; 
            int best_idx = 0;
            for (int t = 1; t < n_centres_new; ++t) {
              double dist = p_new_dist[q + t * n_obs];
              if (dist < min_dist) {
                min_dist = dist;
                best_idx = t;
              }
            }
            p_result[q] = best_idx + 1;
          } else {
            double old_min_dist = p_dist[q + old_best * n_obs];
            double new_dist = p_new_dist[q + mod_offset];
            p_result[q] = (new_dist < old_min_dist) ? (mod_idx + 1) : (old_best + 1);
          }
        }
      } else if (modification == "AC") {
        int mod_offset = mod_idx * n_obs;
        for (int q = 0; q < n_obs; ++q) {
          int old_best = p_old_idx[q] - 1;
          double old_min_dist = p_dist[q + old_best * n_obs];
          double new_dist = p_new_dist[q + mod_offset];
          p_result[q] = (new_dist < old_min_dist) ? (mod_idx + 1) : (old_best + 1);
        }
      } else if (modification == "RC") {
        for (int q = 0; q < n_obs; ++q) {
          int old_best = p_old_idx[q] - 1;
          if (old_best == mod_idx) {
            double min_dist = p_new_dist[q]; 
            int best_idx = 0;
            for (int t = 1; t < n_centres_new; ++t) {
              double dist = p_new_dist[q + t * n_obs];
              if (dist < min_dist) {
                min_dist = dist;
                best_idx = t;
              }
            }
            p_result[q] = best_idx + 1;
          } else {
            p_result[q] = (old_best > mod_idx) ? old_best : (old_best + 1);
          }
        }
      } else {
        for (int q = 0; q < n_obs; ++q) {
          double min_dist = p_new_dist[q]; 
          int best_idx = 0;
          for (int t = 1; t < n_centres_new; ++t) {
            double current_dist = p_new_dist[q + t * n_obs];
            if (current_dist < min_dist) {
              min_dist = current_dist;
              best_idx = t;
            }
          }
          p_result[q] = best_idx + 1;
        }
      }
    } else {
      std::vector<std::pair<double, int>> distances(n_centres_new);
      for (int q = 0; q < n_obs; ++q) {
        for (int t = 0; t < n_centres_new; ++t) {
          distances[t].first = p_new_dist[q + t * n_obs];
          distances[t].second = t + 1; 
        }
        std::partial_sort(distances.begin(), distances.begin() + k, distances.end());
        for (int i = 0; i < k; ++i) {
          p_result[q + i * n_obs] = distances[i].second;
        }
      }
    }
    
    PROTECT(result_list = Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(result_list, 0, res_nn);
    SET_VECTOR_ELT(result_list, 1, res_dist);
    
    PROTECT(list_names = Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(list_names, 0, Rf_mkChar("nn_index"));
    SET_STRING_ELT(list_names, 1, Rf_mkChar("new_distance_matrix"));
    Rf_setAttrib(result_list, R_NamesSymbol, list_names);
    
    UNPROTECT(4);
    return result_list;
  }
  
  SEXP propose_tessellation_cpp(SEXP tess_j_sexp, SEXP dim_j_sexp, SEXP sd_sexp, SEXP mu_sexp, SEXP num_cov_sexp, SEXP s_weights_sexp) {
    
    const double* __restrict__ p_tess_j = REAL(tess_j_sexp);
    const int* __restrict__ p_dim_j = INTEGER(dim_j_sexp);
    const double* __restrict__ sd = REAL(sd_sexp);
    const double* __restrict__ mu = REAL(mu_sexp);
    int numCovariates = INTEGER(num_cov_sexp)[0];
    const double* __restrict__ p_s_weights = REAL(s_weights_sexp);
    
    int tess_j_rows = Rf_nrows(tess_j_sexp);
    int d_j_length = Rf_length(dim_j_sexp);
    
    std::vector<int> dim_j_star(p_dim_j, p_dim_j + d_j_length);
    std::vector<double> tess_j_star(p_tess_j, p_tess_j + (tess_j_rows * d_j_length));
    std::string modification = "Change";
    
    double new_val = 0.0;
    int modified_idx = 0;
    
    double p = unif_rand();
    
    if ((p < 0.2 && d_j_length != numCovariates) || (d_j_length == 1 && d_j_length != numCovariates && p < 0.4)) {
      modification = "AD";
      
      std::vector<int> eligible;
      double sum_w = 0.0;
      for (int i = 1; i <= numCovariates; ++i) {
        if (!in_vector(i, dim_j_star)) {
          eligible.push_back(i);
          sum_w += p_s_weights[i - 1];
        }
      }
      int new_dim = eligible.back();
      if (eligible.size() > 1) {
        double u = unif_rand() * sum_w;
        double cum = 0.0;
        for (int e : eligible) {
          cum += p_s_weights[e - 1];
          if (u <= cum) {
            new_dim = e;
            break;
          }
        }
      }
      
      dim_j_star.push_back(new_dim);
      modified_idx = d_j_length + 1; 
      
      std::vector<double> new_tess(tess_j_rows * (d_j_length + 1));
      for (int r = 0; r < tess_j_rows; ++r) {
        for (int c = 0; c < d_j_length; ++c) {
          new_tess[r + c * tess_j_rows] = tess_j_star[r + c * tess_j_rows];
        }
        new_val = mu[new_dim-1] + norm_rand() * sd[0];
        new_tess[r + d_j_length * tess_j_rows] = new_val;
      }
      tess_j_star = new_tess;
      
    } else if (p < 0.4 && d_j_length > 1) {
      modification = "RD";
      int removed_dim_idx = floor(unif_rand() * d_j_length);
      dim_j_star.erase(dim_j_star.begin() + removed_dim_idx);
      modified_idx = removed_dim_idx + 1;
      
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
      modified_idx = tess_j_rows + 1;
      for (int i = 0; i < d_j_length; ++i) {
        new_val = mu[i] + norm_rand() * sd[0];
        tess_j_star.insert(tess_j_star.begin() + (i * (tess_j_rows + 1)) + tess_j_rows, new_val);
      }
      
    } else if (p < 0.8 && tess_j_rows > 1) {
      modification = "RC";
      int removed_row_idx = floor(unif_rand() * tess_j_rows);
      modified_idx = removed_row_idx + 1;
      
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
      modified_idx = centre_to_change_idx + 1;
      for (int c = 0; c < d_j_length; ++c) {
        new_val = mu[c] + norm_rand() * sd[0];
        tess_j_star[centre_to_change_idx + c * tess_j_rows] = new_val;
      }
      
    } else {
      modification = "Swap";
      int dim_to_change_idx = floor(unif_rand() * d_j_length);
      modified_idx = dim_to_change_idx + 1;
      
      std::vector<int> eligible;
      double sum_w = 0.0;
      for (int i = 1; i <= numCovariates; ++i) {
        if (!in_vector(i, dim_j_star)) {
          eligible.push_back(i);
          sum_w += p_s_weights[i - 1];
        }
      }
      int new_dim = eligible.back();
      if (eligible.size() > 1) {
        double u = unif_rand() * sum_w;
        double cum = 0.0;
        for (int e : eligible) {
          cum += p_s_weights[e - 1];
          if (u <= cum) {
            new_dim = e;
            break;
          }
        }
      }
      
      dim_j_star[dim_to_change_idx] = new_dim;
      for (int r = 0; r < tess_j_rows; ++r) {
        new_val = mu[dim_to_change_idx] + norm_rand() * sd[0];
        tess_j_star[r + dim_to_change_idx * tess_j_rows] = new_val;
      }
    }
    
    SEXP res_tess_j_star, res_dim_j_star, res_mod, res_row_column_modified, result_list, list_names;
    
    int new_rows = tess_j_star.size() / dim_j_star.size();
    PROTECT(res_tess_j_star = Rf_allocMatrix(REALSXP, new_rows, dim_j_star.size()));
    memcpy(REAL(res_tess_j_star), tess_j_star.data(), tess_j_star.size() * sizeof(double));
    
    PROTECT(res_dim_j_star = Rf_allocVector(INTSXP, dim_j_star.size()));
    memcpy(INTEGER(res_dim_j_star), dim_j_star.data(), dim_j_star.size() * sizeof(int));
    
    PROTECT(res_mod = Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(res_mod, 0, Rf_mkChar(modification.c_str()));
    
    PROTECT(res_row_column_modified = Rf_allocVector(INTSXP, 1));
    INTEGER(res_row_column_modified)[0] = modified_idx;
    
    PROTECT(result_list = Rf_allocVector(VECSXP, 4));
    SET_VECTOR_ELT(result_list, 0, res_tess_j_star);
    SET_VECTOR_ELT(result_list, 1, res_dim_j_star);
    SET_VECTOR_ELT(result_list, 2, res_mod);
    SET_VECTOR_ELT(result_list, 3, res_row_column_modified);
    
    PROTECT(list_names = Rf_allocVector(STRSXP, 4));
    SET_STRING_ELT(list_names, 0, Rf_mkChar("tess_j_star"));
    SET_STRING_ELT(list_names, 1, Rf_mkChar("dim_j_star"));
    SET_STRING_ELT(list_names, 2, Rf_mkChar("Modification"));
    SET_STRING_ELT(list_names, 3, Rf_mkChar("row_column_modified"));
    Rf_setAttrib(result_list, R_NamesSymbol, list_names);
    
    UNPROTECT(6);
    return result_list;
  }
  
  SEXP super_mcmc_loop_cpp(SEXP x_sexp, SEXP y_sexp, SEXP sum_sexp, SEXP tess_list, SEXP dim_list, 
                           SEXP indices_list, SEXP sqdist_list, SEXP pred_list, SEXP m_sexp, SEXP p_sexp, 
                           SEXP sd_sexp, SEXP mu_sexp, SEXP sigma2mu_sexp, SEXP omega_sexp, SEXP poisson_lambda_sexp,
                           SEXP total_iter_sexp, SEXP burn_in_sexp, SEXP thinning_sexp,
                           SEXP nu_sexp, SEXP invgamma_lambda_sexp, SEXP show_progress_sexp,
                           SEXP alpha_sexp, SEXP a_alpha_sexp, SEXP b_alpha_sexp, SEXP rho_alpha_sexp, SEXP dirichlet_warmup_sexp,
                           SEXP adapt_boost_sexp, SEXP adapt_penalty_sexp, SEXP momentum_decay_sexp) {
    
    int m = INTEGER(m_sexp)[0];
    int p = INTEGER(p_sexp)[0];
    int n_obs = Rf_length(y_sexp);
    double sigma2mu = REAL(sigma2mu_sexp)[0];
    double omega = REAL(omega_sexp)[0];
    double lambda_poisson = REAL(poisson_lambda_sexp)[0];
    
    int total_iter = INTEGER(total_iter_sexp)[0];
    int burn_in = INTEGER(burn_in_sexp)[0];
    int thinning = INTEGER(thinning_sexp)[0];
    double nu = REAL(nu_sexp)[0];
    double lambda_invgamma = REAL(invgamma_lambda_sexp)[0];
    int show_progress = INTEGER(show_progress_sexp)[0];
    
    double alpha = REAL(alpha_sexp)[0];
    double a_alpha = REAL(a_alpha_sexp)[0];
    double b_alpha = REAL(b_alpha_sexp)[0];
    double rho_alpha = REAL(rho_alpha_sexp)[0];
    int dirichlet_warmup = INTEGER(dirichlet_warmup_sexp)[0];
    
    double adapt_boost = REAL(adapt_boost_sexp)[0];
    double adapt_penalty = REAL(adapt_penalty_sexp)[0];
    double momentum_decay = REAL(momentum_decay_sexp)[0];
    
    // [The rest of the initialisation code remains exactly the same]
    
    int num_stored = std::max(0, (total_iter - burn_in) / thinning);
    
    SEXP return_tess, return_dim, return_indices, return_sqdist, return_pred;
    PROTECT(return_tess = Rf_duplicate(tess_list));
    PROTECT(return_dim = Rf_duplicate(dim_list));
    PROTECT(return_indices = Rf_duplicate(indices_list));
    PROTECT(return_sqdist = Rf_duplicate(sqdist_list));
    PROTECT(return_pred = Rf_duplicate(pred_list));
    
    SEXP return_sum;
    PROTECT(return_sum = Rf_duplicate(sum_sexp));
    double* __restrict__ p_sum = REAL(return_sum);
    const double* __restrict__ p_y = REAL(y_sexp);
    
    SEXP out_tess = PROTECT(Rf_allocVector(VECSXP, num_stored));
    SEXP out_dim = PROTECT(Rf_allocVector(VECSXP, num_stored));
    SEXP out_pred = PROTECT(Rf_allocVector(VECSXP, num_stored));
    SEXP out_sigma = PROTECT(Rf_allocVector(REALSXP, num_stored));
    SEXP out_pred_mat = PROTECT(Rf_allocMatrix(REALSXP, n_obs, num_stored));
    double* __restrict__ p_out_mat = REAL(out_pred_mat);
    
    SEXP k_sexp = PROTECT(Rf_allocVector(INTSXP, 1));
    INTEGER(k_sexp)[0] = 1;
    
    SEXP s_weights_sexp = PROTECT(Rf_allocVector(REALSXP, p));
    double* __restrict__ p_s_weights = REAL(s_weights_sexp);
    for(int i = 0; i < p; ++i) p_s_weights[i] = 1.0 / p;
    
    SEXP out_dirichlet = PROTECT(Rf_allocMatrix(REALSXP, p, num_stored));
    SEXP out_var_sel = PROTECT(Rf_allocMatrix(INTSXP, p, num_stored));
    SEXP out_alpha = PROTECT(Rf_allocVector(REALSXP, num_stored));
    
    GetRNGstate();
    
    std::vector<double> adaptive_momentum(p, 0.0);
    std::vector<double> R_j(n_obs);
    std::vector<double> rIjOld;
    std::vector<int> nIjOld;
    std::vector<double> denomOld;
    std::vector<double> rIjNew;
    std::vector<int> nIjNew;
    std::vector<double> denomNew;
    
    int store_idx = 0;
    
    for (int iter = 0; iter < total_iter; ++iter) {
      
      if (show_progress && (iter + 1) % std::max(1, total_iter / 100) == 0) {
        Rprintf("\rProcessing MCMC Iteration %d of %d...", iter + 1, total_iter);
        R_FlushConsole();
        R_CheckUserInterrupt();
      }
      
      double current_sigma2 = sample_sigma_squared_cpp(n_obs, nu, lambda_invgamma, p_y, p_sum);
      
      for (int j = 0; j < m; ++j) {
        
        SEXP old_pred_sexp = VECTOR_ELT(return_pred, j);
        const double* __restrict__ p_old_pred = REAL(old_pred_sexp);
        SEXP old_idx_sexp = VECTOR_ELT(return_indices, j);
        const int* __restrict__ p_old_idx = INTEGER(old_idx_sexp);
        double* __restrict__ p_R_j = R_j.data();
        
        for(int i = 0; i < n_obs; ++i) {
          p_sum[i] -= p_old_pred[p_old_idx[i] - 1];
          p_R_j[i] = p_y[i] - p_sum[i];
        }
        
        SEXP newTessOutput = propose_tessellation_cpp(VECTOR_ELT(return_tess, j), VECTOR_ELT(return_dim, j), sd_sexp, mu_sexp, p_sexp, s_weights_sexp);
        PROTECT(newTessOutput);
        SEXP tess_j_star_sexp = VECTOR_ELT(newTessOutput, 0);
        SEXP dim_j_star_sexp = VECTOR_ELT(newTessOutput, 1);
        SEXP mod_sexp = VECTOR_ELT(newTessOutput, 2);
        SEXP rcm_sexp = VECTOR_ELT(newTessOutput, 3);
        
        SEXP dist_calc = knnx_index_cpp(x_sexp, VECTOR_ELT(return_tess, j), VECTOR_ELT(return_dim, j), tess_j_star_sexp, k_sexp, dim_j_star_sexp, VECTOR_ELT(return_sqdist, j), mod_sexp, rcm_sexp, old_idx_sexp);
        PROTECT(dist_calc);
        SEXP indexesStar_sexp = VECTOR_ELT(dist_calc, 0);
        SEXP sqdist_j_star_sexp = VECTOR_ELT(dist_calc, 1);
        
        int num_levels_old = Rf_nrows(VECTOR_ELT(return_tess, j));
        int num_centres_new = Rf_nrows(tess_j_star_sexp);
        
        rIjOld.assign(num_levels_old, 0.0);
        nIjOld.assign(num_levels_old, 0);
        rIjNew.assign(num_centres_new, 0.0);
        nIjNew.assign(num_centres_new, 0);
        
        const int* __restrict__ p_new_idx = INTEGER(indexesStar_sexp);
        for(int i = 0; i < n_obs; ++i) {
          int idx_old = p_old_idx[i] - 1;
          rIjOld[idx_old] += p_R_j[i];
          nIjOld[idx_old]++;
          
          int idx_new = p_new_idx[i] - 1;
          rIjNew[idx_new] += p_R_j[i];
          nIjNew[idx_new]++;
        }
        
        denomOld.assign(num_levels_old, 0.0);
        denomNew.assign(num_centres_new, 0.0);
        
        double* __restrict__ p_denomOld = denomOld.data();
        const int* __restrict__ p_nIjOld_calc = nIjOld.data();
        for(int c = 0; c < num_levels_old; ++c) {
          p_denomOld[c] = sigma2mu * p_nIjOld_calc[c] + current_sigma2;
        }
        
        double* __restrict__ p_denomNew = denomNew.data();
        const int* __restrict__ p_nIjNew_calc = nIjNew.data();
        for(int c = 0; c < num_centres_new; ++c) {
          p_denomNew[c] = sigma2mu * p_nIjNew_calc[c] + current_sigma2;
        }
        
        int old_dims = Rf_length(VECTOR_ELT(return_dim, j));
        int new_dims = Rf_length(dim_j_star_sexp);
        std::string mod = CHAR(STRING_ELT(mod_sexp, 0));
        
        double sum_notin_old = 1.0;
        const int* p_old_dim = INTEGER(VECTOR_ELT(return_dim, j));
        for (int k = 0; k < old_dims; ++k) {
          sum_notin_old -= p_s_weights[p_old_dim[k] - 1];
        }
        sum_notin_old = fmax(sum_notin_old, 1e-10);
        
        double sum_notin_new = 1.0;
        const int* p_new_dim = INTEGER(dim_j_star_sexp);
        for (int k = 0; k < new_dims; ++k) {
          sum_notin_new -= p_s_weights[p_new_dim[k] - 1];
        }
        sum_notin_new = fmax(sum_notin_new, 1e-10);
        
        double log_prob = calc_acceptance_cpp(rIjOld, denomOld, rIjNew, denomNew, num_levels_old, num_centres_new, old_dims, new_dims, current_sigma2, mod, sigma2mu, omega, lambda_poisson, p, sum_notin_old, sum_notin_new);
        
        int mod_idx_cpp = INTEGER(rcm_sexp)[0] - 1;
        double acceptance_prob = exp(fmin(0.0, log_prob));
        
        if (mod == "AD") {
          int added_cov = p_new_dim[mod_idx_cpp] - 1;
          adaptive_momentum[added_cov] += adapt_boost * acceptance_prob;
        } else if (mod == "RD") {
          int removed_cov = p_old_dim[mod_idx_cpp] - 1;
          adaptive_momentum[removed_cov] -= adapt_penalty * acceptance_prob;
        } else if (mod == "Swap") {
          int added_cov = p_new_dim[mod_idx_cpp] - 1;
          int removed_cov = p_old_dim[mod_idx_cpp] - 1;
          adaptive_momentum[added_cov] += adapt_boost * acceptance_prob;
          adaptive_momentum[removed_cov] -= adapt_penalty * acceptance_prob;
        }
        
        if(log(unif_rand()) < log_prob) {
          SET_VECTOR_ELT(return_tess, j, tess_j_star_sexp);
          SET_VECTOR_ELT(return_dim, j, dim_j_star_sexp);
          SET_VECTOR_ELT(return_indices, j, indexesStar_sexp);
          SET_VECTOR_ELT(return_sqdist, j, sqdist_j_star_sexp);
          
          std::vector<double> new_mu = sample_mu_cpp(num_centres_new, rIjNew, denomNew, sigma2mu, current_sigma2);
          SEXP new_pred_sexp = PROTECT(Rf_allocMatrix(REALSXP, num_centres_new, 1));
          memcpy(REAL(new_pred_sexp), new_mu.data(), num_centres_new * sizeof(double));
          SET_VECTOR_ELT(return_pred, j, new_pred_sexp);
          UNPROTECT(1);
        } else {
          std::vector<double> old_mu = sample_mu_cpp(num_levels_old, rIjOld, denomOld, sigma2mu, current_sigma2);
          SEXP new_pred_sexp = PROTECT(Rf_allocMatrix(REALSXP, num_levels_old, 1));
          memcpy(REAL(new_pred_sexp), old_mu.data(), num_levels_old * sizeof(double));
          SET_VECTOR_ELT(return_pred, j, new_pred_sexp);
          UNPROTECT(1);
        }
        
        SEXP active_pred_sexp = VECTOR_ELT(return_pred, j);
        const double* __restrict__ p_active_pred = REAL(active_pred_sexp);
        SEXP active_idx_sexp = VECTOR_ELT(return_indices, j);
        const int* __restrict__ p_active_idx = INTEGER(active_idx_sexp);
        
        for(int i = 0; i < n_obs; ++i) {
          p_sum[i] += p_active_pred[p_active_idx[i] - 1];
        }
        
        UNPROTECT(2);
      }
      
      if (iter > dirichlet_warmup) {
        std::vector<double> m_counts(p, 0.0);
        for (int j = 0; j < m; ++j) {
          SEXP dim_j_sexp = VECTOR_ELT(return_dim, j);
          int d_len = Rf_length(dim_j_sexp);
          const int* p_dim_j = INTEGER(dim_j_sexp);
          
          for(int k = 0; k < d_len; ++k) {
            m_counts[p_dim_j[k] - 1] += 1.0;
          }
          
          if (d_len > 1) {
            std::vector<int> ordered_dims(p_dim_j, p_dim_j + d_len);
            for (int k = d_len - 1; k > 0; --k) {
              int r = floor(unif_rand() * (k + 1));
              std::swap(ordered_dims[k], ordered_dims[r]);
            }
            
            for (int k = 1; k < d_len; ++k) {
              std::vector<int> ineligible;
              double sum_inel = 0.0;
              for (int e = 0; e < k; ++e) {
                ineligible.push_back(ordered_dims[e]);
                sum_inel += p_s_weights[ordered_dims[e] - 1];
              }
              
              double prob_success = fmax(1e-10, fmin(1.0 - 1e-10, 1.0 - sum_inel));
              double failures = rgeom(prob_success);
              if (failures > 0) {
                if (ineligible.size() == 1) {
                  m_counts[ineligible[0] - 1] += failures;
                } else {
                  // Mathematically distribute failures instantly to avoid looping billions of times
                  double sum_w = fmax(sum_inel, 1e-16); 
                  for (int e : ineligible) {
                    m_counts[e - 1] += failures * (p_s_weights[e - 1] / sum_w);
                  }
                }
              }
            }
          }
        }
        
        double sum_s = 0.0;
        for (int i = 0; i < p; ++i) {
          double raw_shape = (alpha / p) + m_counts[i] + adaptive_momentum[i];
          double final_shape = fmax(raw_shape, 1e-10);
          
          p_s_weights[i] = fmax(rgamma(final_shape, 1.0), 1e-16);
          sum_s += p_s_weights[i];
        }
        
        for (int i = 0; i < p; ++i) {
          p_s_weights[i] /= sum_s;
          adaptive_momentum[i] *= momentum_decay; 
        }
        
        double alpha_star = exp(norm_rand() * 0.5 + log(alpha));
        
        auto ll_alpha = [&](double a_val) {
          double sum_log_s = 0.0;
          for (int i = 0; i < p; ++i) {
            sum_log_s += log(fmax(p_s_weights[i], 1e-16));
          }
          return lgammafn(a_val) - p * lgammafn(a_val / p) + 
            (a_val / p) * sum_log_s + 
            (a_alpha - 1.0) * log(a_val) - 
            (a_alpha + b_alpha) * log(a_val + rho_alpha);
        };
        
        double log_accept_ratio = ll_alpha(alpha_star) - ll_alpha(alpha) + log(alpha_star) - log(alpha);
        if (log(unif_rand()) < log_accept_ratio) {
          alpha = alpha_star;
        }
      } else {
        for(int i = 0; i < p; ++i) {
          p_s_weights[i] = 1.0 / p;
        }
      }
      
      if (iter >= burn_in && (iter + 1 - burn_in) % thinning == 0) {
        SET_VECTOR_ELT(out_tess, store_idx, Rf_duplicate(return_tess));
        SET_VECTOR_ELT(out_dim, store_idx, Rf_duplicate(return_dim));
        SET_VECTOR_ELT(out_pred, store_idx, Rf_duplicate(return_pred));
        REAL(out_sigma)[store_idx] = current_sigma2;
        
        std::memcpy(p_out_mat + store_idx * n_obs, p_sum, n_obs * sizeof(double));
        
        for(int i = 0; i < p; ++i) {
          REAL(out_dirichlet)[i + store_idx * p] = p_s_weights[i];
        }
        REAL(out_alpha)[store_idx] = alpha;
        
        std::vector<int> current_covariate_counts(p, 0);
        for (int j = 0; j < m; ++j) {
          SEXP dim_j_sexp = VECTOR_ELT(return_dim, j);
          int d_len = Rf_length(dim_j_sexp);
          const int* p_dim_j = INTEGER(dim_j_sexp);
          for(int k = 0; k < d_len; ++k) {
            current_covariate_counts[p_dim_j[k] - 1] += 1;
          }
        }
        for(int i = 0; i < p; ++i) {
          INTEGER(out_var_sel)[i + store_idx * p] = current_covariate_counts[i];
        }
        
        store_idx++;
      }
    }
    
    if (show_progress) {
      Rprintf("\n");
      R_FlushConsole();
    }
    
    PutRNGstate();
    
    SEXP return_list, list_names;
    PROTECT(return_list = Rf_allocVector(VECSXP, 8));
    SET_VECTOR_ELT(return_list, 0, out_tess);
    SET_VECTOR_ELT(return_list, 1, out_dim);
    SET_VECTOR_ELT(return_list, 2, out_pred);
    SET_VECTOR_ELT(return_list, 3, out_sigma);
    SET_VECTOR_ELT(return_list, 4, out_pred_mat);
    SET_VECTOR_ELT(return_list, 5, out_dirichlet);
    SET_VECTOR_ELT(return_list, 6, out_var_sel);
    SET_VECTOR_ELT(return_list, 7, out_alpha);
    
    PROTECT(list_names = Rf_allocVector(STRSXP, 8));
    SET_STRING_ELT(list_names, 0, Rf_mkChar("posteriorTess"));
    SET_STRING_ELT(list_names, 1, Rf_mkChar("posteriorDim"));
    SET_STRING_ELT(list_names, 2, Rf_mkChar("posteriorPred"));
    SET_STRING_ELT(list_names, 3, Rf_mkChar("posteriorSigma"));
    SET_STRING_ELT(list_names, 4, Rf_mkChar("predictionMatrix"));
    SET_STRING_ELT(list_names, 5, Rf_mkChar("posteriorDirichletWeights"));
    SET_STRING_ELT(list_names, 6, Rf_mkChar("posteriorVariableSelection"));
    SET_STRING_ELT(list_names, 7, Rf_mkChar("posteriorAlpha"));
    Rf_setAttrib(return_list, R_NamesSymbol, list_names);
    
    UNPROTECT(18);
    return return_list;
  }
  
  SEXP knnx_index_predict_cpp(SEXP tess_sexp, SEXP query_sexp, SEXP k_sexp, SEXP dim_sexp) {
    
    const double* __restrict__ p_tess = REAL(tess_sexp);
    const double* __restrict__ p_query = REAL(query_sexp);
    int k = INTEGER(k_sexp)[0];
    const int* __restrict__ dim_p = INTEGER(dim_sexp);
    
    int tess_rows = Rf_nrows(tess_sexp);
    int tess_cols = Rf_ncols(tess_sexp);
    int query_rows = Rf_nrows(query_sexp);
    int query_cols = Rf_ncols(query_sexp);
    
    if (tess_cols != query_cols) {
      Rf_error("Dimensions of tess and query matrices must match");
    }
    if (k <= 0 || k > tess_rows) {
      Rf_error("k must be positive and not greater than number of reference points");
    }
    
    int n_dims = Rf_length(dim_sexp);
    std::vector<int> active_dim_idx;
    active_dim_idx.reserve(n_dims);
    std::vector<char> active_dim_mask(query_cols, 0);
    
    for (int i = 0; i < n_dims; ++i) {
      int d0 = dim_p[i] - 1;
      if (d0 < 0 || d0 >= query_cols) {
        Rf_error("Values in dim must be valid 1-based column indices");
      }
      if (!active_dim_mask[d0]) {
        active_dim_mask[d0] = 1;
        active_dim_idx.push_back(d0);
      }
    }
    
    SEXP result;
    PROTECT(result = Rf_allocMatrix(INTSXP, query_rows, k));
    int* __restrict__ p_result = INTEGER(result);
    
    if (k == 1) {
      for (int q = 0; q < query_rows; ++q) {
        double min_dist = R_PosInf;
        int best_idx = 0;
        
        for (int t = 0; t < tess_rows; ++t) {
          double dval = 0.0;
          for (size_t i = 0; i < active_dim_idx.size(); ++i) {
            int d = active_dim_idx[i];
            double diff = p_query[q + d * query_rows] - p_tess[t + d * tess_rows];
            dval += diff * diff;
            if (dval >= min_dist) {
              break;
            }
          }
          
          if (dval < min_dist) {
            min_dist = dval;
            best_idx = t;
          }
        }
        p_result[q] = best_idx + 1;
      }
    } else {
      std::vector<std::pair<double, int>> distances(tess_rows);
      for (int q = 0; q < query_rows; ++q) {
        for (int t = 0; t < tess_rows; ++t) {
          double dval = 0.0;
          for (size_t i = 0; i < active_dim_idx.size(); ++i) {
            int d = active_dim_idx[i];
            double diff = p_query[q + d * query_rows] - p_tess[t + d * tess_rows];
            dval += diff * diff;
          }
          distances[t].first = dval;
          distances[t].second = t + 1; 
        }
        
        std::partial_sort(distances.begin(), distances.begin() + k, distances.end());
        for (int i = 0; i < k; ++i) {
          p_result[q + i * query_rows] = distances[i].second;
        }
      }
    }
    
    UNPROTECT(1);
    return result;
  }
}