#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#include <RcppEigen.h>
#pragma GCC diagnostic pop

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
  
  struct Hypers {
    double sigma2; double sigma2mu; double omega; double lambda_poisson;
    double p_shape; double p_rate; double p_sd; double nu; double lambda_invgamma;
    int num_covariates;
    void UpdateSigma(const Eigen::VectorXd& res) {
      double sse = res.squaredNorm(); double n = res.size();
      double shape = (nu + n) / 2.0; double rate = (nu * lambda_invgamma + sse) / 2.0;
      sigma2 = 1.0 / R::rgamma(shape, 1.0 / rate);
    }
  };
  
  struct PosteriorStats {
    Eigen::VectorXd mu_post; Eigen::MatrixXd L_chol;
    double quad_form; double log_det_Q; bool valid;
  };
  
  struct Tessellation {
    std::vector<double> tess; std::vector<int> dim;
    int n_centres; int n_dims; double p_val; Eigen::VectorXd mu;
    Eigen::MatrixXd cached_sqdist; Eigen::MatrixXd cached_weights; PosteriorStats current_stats;
  };
  
  double calc_acceptance_hard(const std::vector<double>& rIjOld, const std::vector<double>& denomOld,
                              const std::vector<double>& rIjNew, const std::vector<double>& denomNew,
                              int old_centres, int new_centres, int d_old, int d_new,
                              double sigma2, const std::string& mod, double sigma2mu, double omega, double lambda, int p,
                              double sum_notin_old, double sum_notin_new, int var_sel_mode) {
    double prob = fmin(1.0 - 1e-10, fmax(1e-10, omega / p)); 
    double sum_log_old = 0.0, sum_log_new = 0.0, sum_R2_old = 0.0, sum_R2_new = 0.0;
    for(size_t i = 0; i < denomOld.size(); ++i) {
      sum_log_old += log(fmax(denomOld[i], 1e-16)); 
      sum_R2_old += (rIjOld[i] * rIjOld[i]) / fmax(denomOld[i], 1e-16); 
    }
    for(size_t i = 0; i < denomNew.size(); ++i) {
      sum_log_new += log(fmax(denomNew[i], 1e-16)); 
      sum_R2_new += (rIjNew[i] * rIjNew[i]) / fmax(denomNew[i], 1e-16); 
    }
    double safe_sigma2 = fmax(sigma2, 1e-10);
    double ratio = 0.5 * (sum_log_old - sum_log_new) + (sigma2mu / (2.0 * safe_sigma2)) * (sum_R2_new - sum_R2_old);
    
    if (mod == "AD") {
      double t_num = R::dbinom(d_new - 1, p - 1, prob, 0); double t_den = R::dbinom(d_old - 1, p - 1, prob, 0);
      ratio += log(fmax(t_num / fmax(t_den, 1e-16), 1e-16)) + log((var_sel_mode == 0) ? 1.0 : (sum_notin_old / d_new));
      if (d_new == 2) ratio += log(0.5); else if (d_new == p) ratio += log(2.0);
    } else if (mod == "RD") {
      double t_num = R::dbinom(d_new - 1, p - 1, prob, 0); double t_den = R::dbinom(d_old - 1, p - 1, prob, 0);
      ratio += log(fmax(t_num / fmax(t_den, 1e-16), 1e-16)) + log((var_sel_mode == 0) ? 1.0 : ((double)(d_old) / sum_notin_new));
      if (d_old == p) ratio += log(0.5); else if (d_old == 2) ratio += log(2.0);
    } else if (mod == "AC") {
      double t_num = R::dpois(new_centres - 1, lambda, 0); double t_den = R::dpois(new_centres - 2, lambda, 0);
      ratio += log(fmax(t_num / fmax(t_den, 1e-16), 1e-16)) + log(1.0 / new_centres) + 0.5 * log(safe_sigma2);
      if (new_centres == 2) ratio += log(0.5); 
    } else if (mod == "RC") {
      double t_num = R::dpois(new_centres - 1, lambda, 0); double t_den = R::dpois(new_centres, lambda, 0);
      ratio += log(fmax(t_num / fmax(t_den, 1e-16), 1e-16)) + log(new_centres + 1.0) - 0.5 * log(safe_sigma2);
      if (new_centres == 1) ratio += log(2.0); 
    } else if (mod == "Swap") {
      ratio += log((var_sel_mode == 0) ? 1.0 : (sum_notin_old / sum_notin_new));
    }
    return std::isnan(ratio) ? -1e10 : ratio;
  }
  
  double calc_acceptance_soft(const PosteriorStats& StatsOld, const PosteriorStats& StatsNew,
                              int old_c, int new_c, int old_d, int new_d,
                              const std::string& mod, const Hypers& hypers, int p,
                              double sum_notin_old, double sum_notin_new, int var_sel_mode) {
    double prob = fmin(1.0 - 1e-10, fmax(0.0, hypers.omega / p));
    double ratio = -0.5 * (new_c - old_c) * log(hypers.sigma2mu) + 0.5 * (StatsOld.log_det_Q - StatsNew.log_det_Q) + 0.5 * (StatsNew.quad_form - StatsOld.quad_form);    
    
    if (mod == "AD") {
      double t_num = R::dbinom(new_d - 1, p - 1, prob, 0); double t_den = R::dbinom(old_d - 1, p - 1, prob, 0);
      ratio += log(fmax(t_num / fmax(t_den, 1e-16), 1e-16)) + log((var_sel_mode == 0) ? 1.0 : (sum_notin_old / new_d));
      if (new_d == 2) ratio += log(0.5); else if (new_d == p) ratio += log(2.0);
    } else if (mod == "RD") {
      double t_num = R::dbinom(new_d - 1, p - 1, prob, 0); double t_den = R::dbinom(old_d - 1, p - 1, prob, 0);
      ratio += log(fmax(t_num / fmax(t_den, 1e-16), 1e-16)) + log((var_sel_mode == 0) ? 1.0 : ((double)(old_d) / sum_notin_new));
      if (old_d == p) ratio += log(0.5); else if (old_d == 2) ratio += log(2.0);
    } else if (mod == "AC") {
      double t_num = R::dpois(new_c - 1, hypers.lambda_poisson, 0); double t_den = R::dpois(new_c - 2, hypers.lambda_poisson, 0);
      ratio += log(fmax(t_num / fmax(t_den, 1e-16), 1e-16)) + log(1.0 / new_c);
      if (new_c == 2) ratio += log(0.5); 
    } else if (mod == "RC") {
      double t_num = R::dpois(new_c - 1, hypers.lambda_poisson, 0); double t_den = R::dpois(new_c, hypers.lambda_poisson, 0);
      ratio += log(fmax(t_num / fmax(t_den, 1e-16), 1e-16)) + log(new_c + 1.0);
      if (new_c == 1) ratio += log(2.0); 
    } else if (mod == "Swap") {
      ratio += log((var_sel_mode == 0) ? 1.0 : (sum_notin_old / sum_notin_new));
    }
    return std::isnan(ratio) ? -1e10 : ratio;
  }
  
  std::vector<double> sample_mu_hard(int new_centres, const std::vector<double>& rIjNew, const std::vector<int>& nIjNew, double sigma2mu, double sigma2) {
    std::vector<double> new_mu(new_centres, 0.0);
    double safe_sigma2 = fmax(sigma2, 1e-10); 
    for(int i = 0; i < new_centres; ++i) {
      if (nIjNew[i] == 0.0) { new_mu[i] = norm_rand() * sqrt(sigma2mu); continue; }
      double denom = fmax(nIjNew[i] * sigma2mu + safe_sigma2, 1e-10); 
      new_mu[i] = norm_rand() * sqrt(fmax((safe_sigma2 * sigma2mu) / denom, 0.0)) + ((sigma2mu * rIjNew[i]) / denom);
    }
    return new_mu;
  }
  
  double sample_sigma_squared_cpp(int n, double nu, double lambda, const double* p_y, const double* p_sum) {
    double sse = 0.0;
    for (int i = 0; i < n; ++i) { double diff = p_y[i] - p_sum[i]; sse += diff * diff; }
    return fmax(1.0 / R::rgamma((nu + n) / 2.0, 1.0 / ((nu * lambda + sse) / 2.0)), 1e-10); 
  }
  
  void calculate_distances_hard(const double* p_data, const double* p_tess, const double* p_tessstar, const double* p_dist,
                                const int* p_old_idx, const int* p_dim, const int* p_dimstar, int k, const std::string& modification,
                                int mod_idx, int n_obs, int n_centres_old, int n_centres_new, int n_dims_new, double* p_new_dist, int* p_result) {
    if (modification == "RC") {
      if (mod_idx > 0) std::memcpy(p_new_dist, p_dist, mod_idx * n_obs * sizeof(double));
      if (mod_idx < n_centres_old - 1) std::memcpy(p_new_dist + (mod_idx * n_obs), p_dist + ((mod_idx + 1) * n_obs), ((n_centres_old - 1) - mod_idx) * n_obs * sizeof(double));
    } else if (modification == "AC") {
      std::memcpy(p_new_dist, p_dist, n_centres_old * n_obs * sizeof(double));
      int c = n_centres_new - 1; int c_offset = c * n_obs;
      for (int i = 0; i < n_obs; ++i) p_new_dist[i + c_offset] = 0.0;
      for (int d = 0; d < n_dims_new; ++d) {
        int cov_offset = (p_dimstar[d] - 1) * n_obs; double tess_val = p_tessstar[c + d * n_centres_new];
        for (int i = 0; i < n_obs; ++i) { double diff = p_data[i + cov_offset] - tess_val; p_new_dist[i + c_offset] += diff * diff; }
      }
    } else if (modification == "Change") {
      std::memcpy(p_new_dist, p_dist, n_centres_new * n_obs * sizeof(double));
      int c = mod_idx; int c_offset = c * n_obs;
      for (int i = 0; i < n_obs; ++i) p_new_dist[i + c_offset] = 0.0;
      for (int d = 0; d < n_dims_new; ++d) {
        int cov_offset = (p_dimstar[d] - 1) * n_obs; double tess_val = p_tessstar[c + d * n_centres_new];
        for (int i = 0; i < n_obs; ++i) { double diff = p_data[i + cov_offset] - tess_val; p_new_dist[i + c_offset] += diff * diff; }
      }
    } else if (modification == "AD") {
      int cov_offset = (p_dimstar[mod_idx] - 1) * n_obs;
      for (int c = 0; c < n_centres_new; ++c) {
        double tess_val = p_tessstar[c + mod_idx * n_centres_new]; int c_offset = c * n_obs;
        for (int i = 0; i < n_obs; ++i) { double diff = p_data[i + cov_offset] - tess_val; p_new_dist[i + c_offset] = p_dist[i + c_offset] + (diff * diff); }
      }
    } else if (modification == "RD") {
      int cov_offset = (p_dim[mod_idx] - 1) * n_obs;
      for (int c = 0; c < n_centres_new; ++c) {
        double tess_val = p_tess[c + mod_idx * n_centres_old]; int c_offset = c * n_obs;
        for (int i = 0; i < n_obs; ++i) { double diff = p_data[i + cov_offset] - tess_val; double val = p_dist[i + c_offset] - (diff * diff); p_new_dist[i + c_offset] = (val < 0.0) ? 0.0 : val; }
      }
    } else if (modification == "Swap") {
      int old_cov_offset = (p_dim[mod_idx] - 1) * n_obs; int new_cov_offset = (p_dimstar[mod_idx] - 1) * n_obs;
      for (int c = 0; c < n_centres_new; ++c) {
        double old_tess_val = p_tess[c + mod_idx * n_centres_old]; double new_tess_val = p_tessstar[c + mod_idx * n_centres_new]; int c_offset = c * n_obs;
        for (int i = 0; i < n_obs; ++i) { double o_diff = p_data[i + old_cov_offset] - old_tess_val; double n_diff = p_data[i + new_cov_offset] - new_tess_val; double val = p_dist[i + c_offset] - (o_diff * o_diff) + (n_diff * n_diff); p_new_dist[i + c_offset] = (val < 0.0) ? 0.0 : val; }
      }
    } else {
      for (int c = 0; c < n_centres_new; ++c) {
        int c_offset = c * n_obs; for (int i = 0; i < n_obs; ++i) p_new_dist[i + c_offset] = 0.0;
        for (int d = 0; d < n_dims_new; ++d) {
          int cov_offset = (p_dimstar[d] - 1) * n_obs; double tess_val = p_tessstar[c + d * n_centres_new];
          for (int i = 0; i < n_obs; ++i) { double diff = p_data[i + cov_offset] - tess_val; p_new_dist[i + c_offset] += diff * diff; }
        }
      }
    }
    
    if (modification == "Change") {
      int mod_offset = mod_idx * n_obs;
      for (int q = 0; q < n_obs; ++q) {
        int old_best = p_old_idx[q] - 1;
        if (old_best != mod_idx) { double old_min_dist = p_dist[q + old_best * n_obs]; double new_dist = p_new_dist[q + mod_offset]; p_result[q] = (new_dist < old_min_dist) ? (mod_idx + 1) : (old_best + 1); } 
        else { p_result[q] = 1; }
      }
      std::vector<double> current_min_dist(n_obs, R_PosInf);
      for(int q = 0; q < n_obs; ++q) if(p_result[q] == 1) current_min_dist[q] = p_new_dist[q];
      for (int t = 1; t < n_centres_new; ++t) {
        int t_offset = t * n_obs;
        for (int q = 0; q < n_obs; ++q) { if (p_old_idx[q] - 1 == mod_idx) { double dist = p_new_dist[q + t_offset]; if (dist < current_min_dist[q]) { current_min_dist[q] = dist; p_result[q] = t + 1; } } }
      }
    } else if (modification == "AC") {
      int mod_offset = mod_idx * n_obs;
      for (int q = 0; q < n_obs; ++q) { int old_best = p_old_idx[q] - 1; double old_min_dist = p_dist[q + old_best * n_obs]; double new_dist = p_new_dist[q + mod_offset]; p_result[q] = (new_dist < old_min_dist) ? (mod_idx + 1) : (old_best + 1); }
    } else if (modification == "RC") {
      for (int q = 0; q < n_obs; ++q) { int old_best = p_old_idx[q] - 1; if (old_best == mod_idx) { p_result[q] = 1; } else { p_result[q] = (old_best > mod_idx) ? old_best : (old_best + 1); } }
      std::vector<double> current_min_dist(n_obs, R_PosInf);
      for(int q = 0; q < n_obs; ++q) if(p_result[q] == 1) current_min_dist[q] = p_new_dist[q];
      for (int t = 1; t < n_centres_new; ++t) {
        int t_offset = t * n_obs;
        for (int q = 0; q < n_obs; ++q) { if (p_old_idx[q] - 1 == mod_idx) { double dist = p_new_dist[q + t_offset]; if (dist < current_min_dist[q]) { current_min_dist[q] = dist; p_result[q] = t + 1; } } }
      }
    } else {
      std::vector<double> current_min_dist(p_new_dist, p_new_dist + n_obs);
      for (int q = 0; q < n_obs; ++q) p_result[q] = 1;
      for (int t = 1; t < n_centres_new; ++t) {
        int t_offset = t * n_obs;
        for (int q = 0; q < n_obs; ++q) { double dist = p_new_dist[q + t_offset]; if (dist < current_min_dist[q]) { current_min_dist[q] = dist; p_result[q] = t + 1; } }
      }
    }
  }
  
  void build_sqdist(const double* p_data, int n_obs, const std::vector<double>& old_tess, int old_c, const std::vector<double>& new_tess, int new_c,
                    const std::vector<int>& old_dim, int old_d, const std::vector<int>& new_dim, int new_d, const Eigen::MatrixXd& old_sqdist, const std::string& mod, int mod_idx, Eigen::MatrixXd& out_sqdist) {
    out_sqdist.resize(n_obs, new_c); double* p_new = out_sqdist.data(); const double* p_old = old_sqdist.data();
    if (mod == "RC") {
      if (mod_idx > 0) std::memcpy(p_new, p_old, mod_idx * n_obs * sizeof(double));
      if (mod_idx < old_c - 1) std::memcpy(p_new + (mod_idx * n_obs), p_old + ((mod_idx + 1) * n_obs), ((old_c - 1) - mod_idx) * n_obs * sizeof(double));
    } else if (mod == "AC") {
      std::memcpy(p_new, p_old, old_c * n_obs * sizeof(double)); int c = new_c - 1; for (int i = 0; i < n_obs; ++i) p_new[i + c * n_obs] = 0.0;
      for (int d = 0; d < new_d; ++d) {
        int offset = (new_dim[d] - 1) * n_obs; double t_val = new_tess[c + d * new_c];
        for (int i = 0; i < n_obs; ++i) { double diff = p_data[i + offset] - t_val; p_new[i + c * n_obs] += diff * diff; }
      }
    } else if (mod == "Change") {
      std::memcpy(p_new, p_old, new_c * n_obs * sizeof(double)); int c = mod_idx; for (int i = 0; i < n_obs; ++i) p_new[i + c * n_obs] = 0.0;
      for (int d = 0; d < new_d; ++d) {
        int offset = (new_dim[d] - 1) * n_obs; double t_val = new_tess[c + d * new_c];
        for (int i = 0; i < n_obs; ++i) { double diff = p_data[i + offset] - t_val; p_new[i + c * n_obs] += diff * diff; }
      }
    } else if (mod == "AD") {
      int offset = (new_dim[mod_idx] - 1) * n_obs;
      for (int c = 0; c < new_c; ++c) {
        double t_val = new_tess[c + mod_idx * new_c]; int c_offset = c * n_obs;
        for (int i = 0; i < n_obs; ++i) { double diff = p_data[i + offset] - t_val; p_new[i + c_offset] = p_old[i + c_offset] + (diff * diff); }
      }
    } else if (mod == "RD") {
      int offset = (old_dim[mod_idx] - 1) * n_obs;
      for (int c = 0; c < new_c; ++c) {
        double t_val = old_tess[c + mod_idx * old_c]; int c_offset = c * n_obs;
        for (int i = 0; i < n_obs; ++i) { double diff = p_data[i + offset] - t_val; double val = p_old[i + c_offset] - (diff * diff); p_new[i + c_offset] = (val < 0.0) ? 0.0 : val; }
      }
    } else if (mod == "Swap") {
      int old_off = (old_dim[mod_idx] - 1) * n_obs; int new_off = (new_dim[mod_idx] - 1) * n_obs;
      for (int c = 0; c < new_c; ++c) {
        double old_t = old_tess[c + mod_idx * old_c]; double new_t = new_tess[c + mod_idx * new_c]; int c_off = c * n_obs;
        for (int i = 0; i < n_obs; ++i) { double o_diff = p_data[i + old_off] - old_t; double n_diff = p_data[i + new_off] - new_t; double val = p_old[i + c_off] - (o_diff * o_diff) + (n_diff * n_diff); p_new[i + c_off] = (val < 0.0) ? 0.0 : val; }
      }
    } else {
      for (int c = 0; c < new_c; ++c) {
        for (int i = 0; i < n_obs; ++i) {
          double dval = 0.0;
          for (int d = 0; d < new_d; ++d) { double diff = p_data[i + (new_dim[d] - 1) * n_obs] - new_tess[c + d * new_c]; dval += diff * diff; }
          p_new[i + c * n_obs] = dval;
        }
      }
    }
  }
  
  void build_weights(const Eigen::MatrixXd& sqdist, int n_dims, double p_val, Eigen::MatrixXd& out_weights, double epsilon = 1e-10) {
    int n_obs = sqdist.rows(); int n_centres = sqdist.cols();
    out_weights.resize(n_obs, n_centres);
    double half_power = p_val / 2.0; double eps_sq = epsilon * epsilon;
    bool fast_path = std::abs(half_power - 1.0) < 1e-7;
    double dim_scale = 1.0 / static_cast<double>(n_dims);
    
    for (int i = 0; i < n_obs; ++i) {
      bool is_zero = false; int z_idx = -1;
      for (int c = 0; c < n_centres; ++c) { if (sqdist(i, c) * dim_scale < eps_sq) { is_zero = true; z_idx = c; break; } }
      if (is_zero) {
        for (int c = 0; c < n_centres; ++c) out_weights(i, c) = (c == z_idx) ? 1.0 : 0.0;
      } else {
        double rsum = 0.0;
        for (int c = 0; c < n_centres; ++c) {
          double sqd = sqdist(i, c) * dim_scale;
          double w = fast_path ? 1.0 / (sqd + epsilon) : 1.0 / (std::pow(sqd, half_power) + epsilon);
          out_weights(i, c) = w; rsum += w;
        }
        double inv_sum = 1.0 / rsum;
        for (int c = 0; c < n_centres; ++c) out_weights(i, c) *= inv_sum;
      }
    }
  }
  
  PosteriorStats build_stats(const Eigen::MatrixXd& w, const Eigen::VectorXd& res, const Hypers& hypers) {
    int k = w.cols(); PosteriorStats stats; Eigen::MatrixXd Q(k, k);
    Q.noalias() = w.transpose() * w; Q /= hypers.sigma2; Q.diagonal().array() += (1.0 / hypers.sigma2mu);
    Eigen::LLT<Eigen::MatrixXd> llt(Q);
    if (llt.info() != Eigen::NumericalIssue) {
      stats.L_chol = llt.matrixL(); Eigen::VectorXd z(k); z.noalias() = (1.0 / hypers.sigma2) * w.transpose() * res;
      stats.mu_post = llt.solve(z); stats.quad_form = stats.mu_post.dot(z); stats.log_det_Q = 2.0 * stats.L_chol.diagonal().array().log().sum(); stats.valid = true;
    } else { stats.valid = false; }
    return stats;
  }
  
  void update_mu_soft(Tessellation& t) {
    Eigen::VectorXd z(t.n_centres); for(int i = 0; i < t.n_centres; ++i) z(i) = norm_rand();
    t.mu = t.current_stats.mu_post + t.current_stats.L_chol.triangularView<Eigen::Lower>().transpose().solve(z);
  }
  
  SEXP propose_tessellation_cpp(SEXP tess_j_sexp, SEXP dim_j_sexp, SEXP sd_sexp, SEXP mu_sexp, SEXP num_cov_sexp, SEXP s_weights_sexp) {
    const double* p_tess_j = REAL(tess_j_sexp); const int* p_dim_j = INTEGER(dim_j_sexp);
    const double* sd = REAL(sd_sexp); const double* mu = REAL(mu_sexp); int numCovariates = INTEGER(num_cov_sexp)[0]; const double* p_s_weights = REAL(s_weights_sexp);
    int mu_len = Rf_length(mu_sexp); int sd_len = Rf_length(sd_sexp);
    auto get_mu = [&](int idx) { return mu_len > 1 ? mu[idx] : mu[0]; }; auto get_sd = [&](int idx) { return sd_len > 1 ? sd[idx] : sd[0]; };
    
    int tess_j_rows = Rf_nrows(tess_j_sexp); int d_j_length = Rf_length(dim_j_sexp);
    std::vector<int> dim_j_star(p_dim_j, p_dim_j + d_j_length); std::vector<double> tess_j_star(p_tess_j, p_tess_j + (tess_j_rows * d_j_length));
    std::string modification = "Change"; int modified_idx = 0; double p = unif_rand();
    
    if ((p < 0.2 && d_j_length != numCovariates) || (d_j_length == 1 && d_j_length != numCovariates && p < 0.4)) {
      modification = "AD"; std::vector<int> eligible; double sum_w = 0.0;
      for (int i = 1; i <= numCovariates; ++i) { if (!in_vector(i, dim_j_star)) { eligible.push_back(i); sum_w += p_s_weights[i - 1]; } }
      int new_dim = eligible.back();
      if (eligible.size() > 1) { double u = unif_rand() * sum_w; double cum = 0.0; for (int e : eligible) { cum += p_s_weights[e - 1]; if (u <= cum) { new_dim = e; break; } } }
      dim_j_star.push_back(new_dim); modified_idx = d_j_length + 1; 
      std::vector<double> new_tess(tess_j_rows * (d_j_length + 1));
      for (int r = 0; r < tess_j_rows; ++r) {
        for (int c = 0; c < d_j_length; ++c) new_tess[r + c * tess_j_rows] = tess_j_star[r + c * tess_j_rows];
        new_tess[r + d_j_length * tess_j_rows] = get_mu(new_dim - 1) + norm_rand() * get_sd(new_dim - 1); 
      }
      tess_j_star = new_tess;
    } else if (p < 0.4 && d_j_length > 1) {
      modification = "RD"; int removed_dim_idx = floor(unif_rand() * d_j_length);
      dim_j_star.erase(dim_j_star.begin() + removed_dim_idx); modified_idx = removed_dim_idx + 1;
      std::vector<double> new_tess(tess_j_rows * (d_j_length - 1)); int current_col = 0;
      for (int c = 0; c < d_j_length; ++c) { if (c != removed_dim_idx) { for (int r = 0; r < tess_j_rows; ++r) new_tess[r + current_col * tess_j_rows] = tess_j_star[r + c * tess_j_rows]; current_col++; } }
      tess_j_star = new_tess;
    } else if (p < 0.6 || (p < 0.8 && tess_j_rows == 1)) {
      modification = "AC"; modified_idx = tess_j_rows + 1;
      for (int i = 0; i < d_j_length; ++i) { int cov_idx = dim_j_star[i] - 1; tess_j_star.insert(tess_j_star.begin() + (i * (tess_j_rows + 1)) + tess_j_rows, get_mu(cov_idx) + norm_rand() * get_sd(cov_idx)); }
    } else if (p < 0.8 && tess_j_rows > 1) {
      modification = "RC"; modified_idx = floor(unif_rand() * tess_j_rows);
      std::vector<double> new_tess; new_tess.reserve((tess_j_rows - 1) * d_j_length);
      for(int c = 0; c < d_j_length; ++c) { for(int r = 0; r < tess_j_rows; ++r) { if(r != modified_idx) new_tess.push_back(tess_j_star[r + c * tess_j_rows]); } }
      tess_j_star = new_tess; modified_idx += 1;
    } else if (p < 0.9 || d_j_length == numCovariates) {
      modification = "Change"; modified_idx = floor(unif_rand() * tess_j_rows);
      for (int c = 0; c < d_j_length; ++c) { int cov_idx = dim_j_star[c] - 1; tess_j_star[modified_idx + c * tess_j_rows] = get_mu(cov_idx) + norm_rand() * get_sd(cov_idx); }
      modified_idx += 1;
    } else {
      modification = "Swap"; modified_idx = floor(unif_rand() * d_j_length);
      std::vector<int> eligible; double sum_w = 0.0;
      for (int i = 1; i <= numCovariates; ++i) { if (!in_vector(i, dim_j_star)) { eligible.push_back(i); sum_w += p_s_weights[i - 1]; } }
      int new_dim = eligible.back();
      if (eligible.size() > 1) { double u = unif_rand() * sum_w; double cum = 0.0; for (int e : eligible) { cum += p_s_weights[e - 1]; if (u <= cum) { new_dim = e; break; } } }
      dim_j_star[modified_idx] = new_dim;
      for (int r = 0; r < tess_j_rows; ++r) { tess_j_star[r + modified_idx * tess_j_rows] = get_mu(new_dim - 1) + norm_rand() * get_sd(new_dim - 1); }
      modified_idx += 1;
    }
    
    SEXP res_tess_j_star, res_dim_j_star, res_mod, res_row_column_modified, result_list, list_names;
    int new_rows = tess_j_star.size() / dim_j_star.size();
    PROTECT(res_tess_j_star = Rf_allocMatrix(REALSXP, new_rows, dim_j_star.size())); memcpy(REAL(res_tess_j_star), tess_j_star.data(), tess_j_star.size() * sizeof(double));
    PROTECT(res_dim_j_star = Rf_allocVector(INTSXP, dim_j_star.size())); memcpy(INTEGER(res_dim_j_star), dim_j_star.data(), dim_j_star.size() * sizeof(int));
    PROTECT(res_mod = Rf_allocVector(STRSXP, 1)); SET_STRING_ELT(res_mod, 0, Rf_mkChar(modification.c_str()));
    PROTECT(res_row_column_modified = Rf_allocVector(INTSXP, 1)); INTEGER(res_row_column_modified)[0] = modified_idx; 
    PROTECT(result_list = Rf_allocVector(VECSXP, 4)); SET_VECTOR_ELT(result_list, 0, res_tess_j_star); SET_VECTOR_ELT(result_list, 1, res_dim_j_star); SET_VECTOR_ELT(result_list, 2, res_mod); SET_VECTOR_ELT(result_list, 3, res_row_column_modified);
    PROTECT(list_names = Rf_allocVector(STRSXP, 4)); SET_STRING_ELT(list_names, 0, Rf_mkChar("tess_j_star")); SET_STRING_ELT(list_names, 1, Rf_mkChar("dim_j_star")); SET_STRING_ELT(list_names, 2, Rf_mkChar("Modification")); SET_STRING_ELT(list_names, 3, Rf_mkChar("row_column_modified"));
    Rf_setAttrib(result_list, R_NamesSymbol, list_names); UNPROTECT(6); return result_list;
  }
}

extern "C" {
  
  SEXP super_mcmc_loop_cpp(SEXP x_sexp, SEXP y_sexp, SEXP sum_sexp, SEXP tess_list, SEXP dim_list, 
                           SEXP indices_list, SEXP sqdist_list, SEXP pred_list, SEXP m_sexp, SEXP p_sexp, 
                           SEXP sd_sexp, SEXP mu_sexp, SEXP sigma2mu_sexp, SEXP omega_sexp, SEXP poisson_lambda_sexp,
                           SEXP total_iter_sexp, SEXP burn_in_sexp, SEXP thinning_sexp,
                           SEXP nu_sexp, SEXP invgamma_lambda_sexp, SEXP show_progress_sexp,
                           SEXP alpha_sexp, SEXP a_alpha_sexp, SEXP b_alpha_sexp, SEXP rho_alpha_sexp, SEXP dirichlet_warmup_sexp,
                           SEXP adapt_boost_sexp, SEXP adapt_penalty_sexp, SEXP momentum_decay_sexp, SEXP kappa_sexp, 
                           SEXP var_sel_mode_sexp, SEXP split_mode_sexp, SEXP power_sexp, SEXP p_shape_sexp, SEXP p_rate_sexp, SEXP p_sd_sexp,
                           SEXP is_class_sexp) { // Added parameter for classification
    
    Hypers hypers;
    int var_sel_mode = INTEGER(var_sel_mode_sexp)[0];
    int split_mode = INTEGER(split_mode_sexp)[0]; 
    int is_classification = INTEGER(is_class_sexp)[0]; // Classification flag
    
    hypers.num_covariates = INTEGER(p_sexp)[0];
    hypers.sigma2mu = REAL(sigma2mu_sexp)[0]; hypers.omega = REAL(omega_sexp)[0]; hypers.lambda_poisson = REAL(poisson_lambda_sexp)[0];
    hypers.nu = REAL(nu_sexp)[0]; hypers.lambda_invgamma = REAL(invgamma_lambda_sexp)[0];
    if (split_mode == 0) { hypers.p_shape = REAL(p_shape_sexp)[0]; hypers.p_rate = REAL(p_rate_sexp)[0]; hypers.p_sd = REAL(p_sd_sexp)[0]; }
    
    int m = INTEGER(m_sexp)[0]; int p = INTEGER(p_sexp)[0]; int n_obs = Rf_length(y_sexp);
    int total_iter = INTEGER(total_iter_sexp)[0]; int burn_in = INTEGER(burn_in_sexp)[0];
    int thinning = INTEGER(thinning_sexp)[0]; int show_progress = INTEGER(show_progress_sexp)[0];
    
    double alpha = REAL(alpha_sexp)[0]; double a_alpha = REAL(a_alpha_sexp)[0]; double b_alpha = REAL(b_alpha_sexp)[0]; double rho_alpha = REAL(rho_alpha_sexp)[0];
    int dirichlet_warmup = INTEGER(dirichlet_warmup_sexp)[0];
    double adapt_boost = REAL(adapt_boost_sexp)[0]; double adapt_penalty = REAL(adapt_penalty_sexp)[0];
    double momentum_decay = REAL(momentum_decay_sexp)[0]; double kappa = REAL(kappa_sexp)[0];
    int num_stored = std::max(0, (total_iter - burn_in) / thinning);
    
    SEXP return_tess, return_dim, return_sqdist, return_pred, return_indices;
    PROTECT(return_tess = Rf_duplicate(tess_list)); PROTECT(return_dim = Rf_duplicate(dim_list)); PROTECT(return_sqdist = Rf_duplicate(sqdist_list)); PROTECT(return_pred = Rf_duplicate(pred_list));
    if (split_mode == 1) { PROTECT(return_indices = Rf_duplicate(indices_list)); } else { PROTECT(return_indices = Rf_allocVector(VECSXP, 1)); }
    
    SEXP return_sum; PROTECT(return_sum = Rf_duplicate(sum_sexp));
    double* p_sum = REAL(return_sum); const double* p_x = REAL(x_sexp); 
    
    // Allocate latent continuous response vector for Probit augmentation
    std::vector<double> Y_latent(n_obs);
    std::memcpy(Y_latent.data(), REAL(y_sexp), n_obs * sizeof(double));
    const double* original_y = REAL(y_sexp);
    
    SEXP out_tess = PROTECT(Rf_allocVector(VECSXP, num_stored)); SEXP out_dim = PROTECT(Rf_allocVector(VECSXP, num_stored)); SEXP out_pred = PROTECT(Rf_allocVector(VECSXP, num_stored)); SEXP out_sigma = PROTECT(Rf_allocVector(REALSXP, num_stored)); SEXP out_pred_mat = PROTECT(Rf_allocMatrix(REALSXP, n_obs, num_stored));
    double* p_out_mat = REAL(out_pred_mat);
    
    SEXP s_weights_sexp = PROTECT(Rf_allocVector(REALSXP, p)); double* p_s_weights = REAL(s_weights_sexp);
    for(int i = 0; i < p; ++i) p_s_weights[i] = 1.0 / p;
    SEXP out_dirichlet = PROTECT(Rf_allocMatrix(REALSXP, p, num_stored)); SEXP out_var_sel = PROTECT(Rf_allocMatrix(INTSXP, p, num_stored)); SEXP out_alpha = PROTECT(Rf_allocVector(REALSXP, num_stored));
    
    SEXP out_power, out_mu;
    if (split_mode == 0) { out_power = PROTECT(Rf_allocMatrix(REALSXP, m, num_stored)); out_mu = PROTECT(Rf_allocVector(VECSXP, num_stored)); } 
    else { out_power = PROTECT(Rf_allocVector(REALSXP, 1)); out_mu = PROTECT(Rf_allocVector(VECSXP, 1)); }
    
    GetRNGstate();
    std::vector<double> adaptive_momentum(p, 0.0); std::vector<double> R_j(n_obs);
    
    std::vector<double> rIjOld, denomOld, rIjNew, denomNew; std::vector<int> nIjOld, nIjNew;
    std::vector<double> workspace_dist; std::vector<int> workspace_idx(n_obs);
    if (split_mode == 1) workspace_dist.reserve(n_obs * 50); 
    
    Eigen::MatrixXd workspace_sqdist; Eigen::MatrixXd workspace_weights; std::vector<Tessellation> ensemble;
    Eigen::Map<Eigen::VectorXd> Y(Y_latent.data(), n_obs); Eigen::VectorXd Y_hat(n_obs); std::vector<Eigen::VectorXd> Forest_Preds;
    
    if (split_mode == 0) {
      workspace_sqdist.resize(n_obs, 50); workspace_weights.resize(n_obs, 50); ensemble.resize(m); 
      double* p_vec = REAL(power_sexp);
      for(int j = 0; j < m; ++j) {
        SEXP rtess = VECTOR_ELT(tess_list, j); SEXP rdim = VECTOR_ELT(dim_list, j);
        ensemble[j].n_centres = Rf_nrows(rtess); ensemble[j].n_dims = Rf_length(rdim);
        ensemble[j].tess.assign(REAL(rtess), REAL(rtess) + (ensemble[j].n_centres * ensemble[j].n_dims)); 
        ensemble[j].dim.assign(INTEGER(rdim), INTEGER(rdim) + ensemble[j].n_dims); 
        ensemble[j].p_val = p_vec[j];
        
        Eigen::MatrixXd empty_sq; 
        build_sqdist(p_x, n_obs, ensemble[j].tess, ensemble[j].n_centres, ensemble[j].tess, ensemble[j].n_centres, 
                     ensemble[j].dim, ensemble[j].n_dims, ensemble[j].dim, ensemble[j].n_dims, empty_sq, "Init", 0, ensemble[j].cached_sqdist);
        build_weights(ensemble[j].cached_sqdist, ensemble[j].n_dims, ensemble[j].p_val, ensemble[j].cached_weights); 
        
        ensemble[j].mu = Eigen::VectorXd::Zero(ensemble[j].n_centres);
        SEXP initial_mu_sexp = VECTOR_ELT(pred_list, j);
        for(int c = 0; c < ensemble[j].n_centres; ++c) {
          ensemble[j].mu(c) = REAL(initial_mu_sexp)[c];
        }
      }
      
      std::memcpy(Y_hat.data(), REAL(sum_sexp), n_obs * sizeof(double)); 
      Forest_Preds.resize(m);
      for(int j = 0; j < m; ++j) {
        Forest_Preds[j] = ensemble[j].cached_weights * ensemble[j].mu;
      }
    }
    
    int store_idx = 0; 
    
    for (int iter = 0; iter < total_iter; ++iter) {
      if (show_progress && (iter + 1) % std::max(1, total_iter / 100) == 0) { Rprintf("\rProcessing MCMC Iteration %d of %d...", iter + 1, total_iter); R_FlushConsole(); R_CheckUserInterrupt(); }
      
      // Albert-Chib Truncated Normal Data Augmentation for Probit Model
      if (is_classification == 1) {
        for(int i = 0; i < n_obs; ++i) {
          double current_mean = (split_mode == 1) ? p_sum[i] : Y_hat(i);
          double u = unif_rand();
          if (original_y[i] == 1.0) {
            double prob_lower = R::pnorm5(0.0, current_mean, 1.0, 1, 0);
            Y_latent[i] = R::qnorm5(prob_lower + u * (1.0 - prob_lower), current_mean, 1.0, 1, 0);
          } else {
            double prob_upper = R::pnorm5(0.0, current_mean, 1.0, 1, 0);
            Y_latent[i] = R::qnorm5(u * prob_upper, current_mean, 1.0, 1, 0);
          }
          if (std::isnan(Y_latent[i]) || std::isinf(Y_latent[i])) {
            Y_latent[i] = (original_y[i] == 1.0) ? std::max(0.001, current_mean) : std::min(-0.001, current_mean);
          }
        }
        hypers.sigma2 = 1.0;
      } else {
        if (split_mode == 1) { hypers.sigma2 = sample_sigma_squared_cpp(n_obs, hypers.nu, hypers.lambda_invgamma, Y_latent.data(), p_sum); } 
        else { hypers.UpdateSigma(Y - Y_hat); }
      }
      
      for (int j = 0; j < m; ++j) {
        SEXP newTessOutput = propose_tessellation_cpp(VECTOR_ELT(return_tess, j), VECTOR_ELT(return_dim, j), sd_sexp, mu_sexp, p_sexp, s_weights_sexp); PROTECT(newTessOutput);
        SEXP tess_j_star_sexp = VECTOR_ELT(newTessOutput, 0); SEXP dim_j_star_sexp = VECTOR_ELT(newTessOutput, 1); SEXP mod_sexp = VECTOR_ELT(newTessOutput, 2); SEXP rcm_sexp = VECTOR_ELT(newTessOutput, 3);
        int num_levels_old = Rf_nrows(VECTOR_ELT(return_tess, j)); int num_centres_new = Rf_nrows(tess_j_star_sexp); int new_dims = Rf_length(dim_j_star_sexp);
        std::string mod = CHAR(STRING_ELT(mod_sexp, 0)); int mod_idx_cpp = INTEGER(rcm_sexp)[0] - 1; int old_dims = Rf_length(VECTOR_ELT(return_dim, j));
        
        double sum_notin_old = 1.0; const int* p_old_dim = INTEGER(VECTOR_ELT(return_dim, j));
        for (int k = 0; k < old_dims; ++k) sum_notin_old -= p_s_weights[p_old_dim[k] - 1]; sum_notin_old = fmax(sum_notin_old, 1e-10);
        double sum_notin_new = 1.0; const int* p_new_dim = INTEGER(dim_j_star_sexp);
        for (int k = 0; k < new_dims; ++k) sum_notin_new -= p_s_weights[p_new_dim[k] - 1]; sum_notin_new = fmax(sum_notin_new, 1e-10);
        
        double log_prob = 0.0;
        
        if (split_mode == 1) {
          SEXP old_pred_sexp = VECTOR_ELT(return_pred, j); const double* p_old_pred = REAL(old_pred_sexp); SEXP old_idx_sexp = VECTOR_ELT(return_indices, j); const int* p_old_idx = INTEGER(old_idx_sexp); double* p_R_j = R_j.data();
          for(int i = 0; i < n_obs; ++i) { p_sum[i] -= p_old_pred[p_old_idx[i] - 1]; p_R_j[i] = Y_latent[i] - p_sum[i]; }
          if (workspace_dist.capacity() < static_cast<size_t>(n_obs * num_centres_new)) { workspace_dist.reserve(n_obs * num_centres_new * 2); } workspace_dist.resize(n_obs * num_centres_new);
          
          calculate_distances_hard(p_x, REAL(VECTOR_ELT(return_tess, j)), REAL(tess_j_star_sexp), REAL(VECTOR_ELT(return_sqdist, j)), p_old_idx, INTEGER(VECTOR_ELT(return_dim, j)), INTEGER(dim_j_star_sexp), 1, mod, mod_idx_cpp, n_obs, num_levels_old, num_centres_new, new_dims, workspace_dist.data(), workspace_idx.data());
          rIjOld.assign(num_levels_old, 0.0); nIjOld.assign(num_levels_old, 0); rIjNew.assign(num_centres_new, 0.0); nIjNew.assign(num_centres_new, 0); const int* p_new_idx = workspace_idx.data();
          for(int i = 0; i < n_obs; ++i) { int idx_old = p_old_idx[i] - 1; rIjOld[idx_old] += p_R_j[i]; nIjOld[idx_old]++; int idx_new = p_new_idx[i] - 1; rIjNew[idx_new] += p_R_j[i]; nIjNew[idx_new]++; }
          denomOld.assign(num_levels_old, 0.0); denomNew.assign(num_centres_new, 0.0);
          for(int c = 0; c < num_levels_old; ++c) denomOld[c] = fmax(hypers.sigma2mu * nIjOld[c] + hypers.sigma2, 1e-10);
          for(int c = 0; c < num_centres_new; ++c) denomNew[c] = fmax(hypers.sigma2mu * nIjNew[c] + hypers.sigma2, 1e-10);
          
          log_prob = calc_acceptance_hard(rIjOld, denomOld, rIjNew, denomNew, num_levels_old, num_centres_new, old_dims, new_dims, hypers.sigma2, mod, hypers.sigma2mu, hypers.omega, hypers.lambda_poisson, p, sum_notin_old, sum_notin_new, var_sel_mode);
        } else {
          Eigen::VectorXd res = Y - (Y_hat - Forest_Preds[j]); ensemble[j].current_stats = build_stats(ensemble[j].cached_weights, res, hypers);
          std::vector<double> prop_tess(REAL(tess_j_star_sexp), REAL(tess_j_star_sexp) + (num_centres_new * new_dims)); std::vector<int> prop_dim(INTEGER(dim_j_star_sexp), INTEGER(dim_j_star_sexp) + new_dims);
          build_sqdist(p_x, n_obs, ensemble[j].tess, num_levels_old, prop_tess, num_centres_new, ensemble[j].dim, old_dims, prop_dim, new_dims, ensemble[j].cached_sqdist, mod, mod_idx_cpp, workspace_sqdist);
          build_weights(workspace_sqdist, new_dims, ensemble[j].p_val, workspace_weights); PosteriorStats prop_stats = build_stats(workspace_weights, res, hypers);
          if (prop_stats.valid && ensemble[j].current_stats.valid) { log_prob = calc_acceptance_soft(ensemble[j].current_stats, prop_stats, num_levels_old, num_centres_new, old_dims, new_dims, mod, hypers, p, sum_notin_old, sum_notin_new, var_sel_mode); } else { log_prob = -1e10; }
          if (log_prob > -1e10) { ensemble[j].current_stats = std::move(prop_stats); }
        }
        
        double acceptance_prob = exp(fmin(0.0, log_prob));
        if (var_sel_mode == 2) {
          double adaptation_scale = (iter > dirichlet_warmup) ? pow(iter - dirichlet_warmup + 1.0, -kappa) : 1.0;
          double current_boost = adapt_boost * adaptation_scale; double current_penalty = adapt_penalty * adaptation_scale;
          const int* p_old_dim = INTEGER(VECTOR_ELT(return_dim, j)); const int* p_new_dim = INTEGER(dim_j_star_sexp);
          if (mod == "AD") { adaptive_momentum[p_new_dim[mod_idx_cpp] - 1] += current_boost * acceptance_prob; } 
          else if (mod == "RD") { adaptive_momentum[p_old_dim[mod_idx_cpp] - 1] -= current_penalty * acceptance_prob; } 
          else if (mod == "Swap") { adaptive_momentum[p_new_dim[mod_idx_cpp] - 1] += current_boost * acceptance_prob; adaptive_momentum[p_old_dim[mod_idx_cpp] - 1] -= current_penalty * acceptance_prob; }
        }
        
        if(log(unif_rand()) < log_prob) {
          SET_VECTOR_ELT(return_tess, j, tess_j_star_sexp); SET_VECTOR_ELT(return_dim, j, dim_j_star_sexp);
          if (split_mode == 1) {
            SEXP new_indices_sexp = PROTECT(Rf_allocMatrix(INTSXP, n_obs, 1)); std::memcpy(INTEGER(new_indices_sexp), workspace_idx.data(), n_obs * sizeof(int)); SET_VECTOR_ELT(return_indices, j, new_indices_sexp);
            SEXP new_sqdist_sexp = PROTECT(Rf_allocMatrix(REALSXP, n_obs, num_centres_new)); std::memcpy(REAL(new_sqdist_sexp), workspace_dist.data(), n_obs * num_centres_new * sizeof(double)); SET_VECTOR_ELT(return_sqdist, j, new_sqdist_sexp);
            std::vector<double> new_mu = sample_mu_hard(num_centres_new, rIjNew, nIjNew, hypers.sigma2mu, hypers.sigma2);
            SEXP new_pred_sexp = PROTECT(Rf_allocMatrix(REALSXP, num_centres_new, 1)); memcpy(REAL(new_pred_sexp), new_mu.data(), num_centres_new * sizeof(double)); SET_VECTOR_ELT(return_pred, j, new_pred_sexp);
            UNPROTECT(3);
          } else {
            ensemble[j].tess.assign(REAL(tess_j_star_sexp), REAL(tess_j_star_sexp) + (num_centres_new * new_dims));
            ensemble[j].dim.assign(INTEGER(dim_j_star_sexp), INTEGER(dim_j_star_sexp) + new_dims);
            ensemble[j].n_centres = num_centres_new; ensemble[j].n_dims = new_dims;
            ensemble[j].cached_sqdist = std::move(workspace_sqdist); ensemble[j].cached_weights = std::move(workspace_weights);
            update_mu_soft(ensemble[j]);
          }
        } else {
          if (split_mode == 1) {
            std::vector<double> old_mu = sample_mu_hard(num_levels_old, rIjOld, nIjOld, hypers.sigma2mu, hypers.sigma2);
            SEXP new_pred_sexp = PROTECT(Rf_allocMatrix(REALSXP, num_levels_old, 1)); memcpy(REAL(new_pred_sexp), old_mu.data(), num_levels_old * sizeof(double)); SET_VECTOR_ELT(return_pred, j, new_pred_sexp);
            UNPROTECT(1);
          } else {
            Eigen::VectorXd res = Y - (Y_hat - Forest_Preds[j]);
            ensemble[j].current_stats = build_stats(ensemble[j].cached_weights, res, hypers);
            update_mu_soft(ensemble[j]);
          }
        }
        
        if (split_mode == 1) {
          SEXP active_pred_sexp = VECTOR_ELT(return_pred, j); const double* p_active_pred = REAL(active_pred_sexp); SEXP active_idx_sexp = VECTOR_ELT(return_indices, j); const int* p_active_idx = INTEGER(active_idx_sexp);
          for(int i = 0; i < n_obs; ++i) p_sum[i] += p_active_pred[p_active_idx[i] - 1];
        } else {
          Eigen::VectorXd res = Y - (Y_hat - Forest_Preds[j]);
          double log_p_star = log(ensemble[j].p_val) + norm_rand() * hypers.p_sd; double p_star = exp(log_p_star);
          build_weights(ensemble[j].cached_sqdist, ensemble[j].n_dims, p_star, workspace_weights); PosteriorStats prop_stats = build_stats(workspace_weights, res, hypers);
          if (prop_stats.valid && ensemble[j].current_stats.valid) {
            double ll_old = 0.5 * ensemble[j].current_stats.quad_form - 0.5 * ensemble[j].current_stats.log_det_Q; double ll_new = 0.5 * prop_stats.quad_form - 0.5 * prop_stats.log_det_Q;
            double prior_old = R::dexp(ensemble[j].p_val, 1.0 / hypers.p_rate, 1); double prior_new = R::dexp(p_star, 1.0 / hypers.p_rate, 1);
            if (log(unif_rand()) < (ll_new - ll_old + prior_new - prior_old + log_p_star - log(ensemble[j].p_val))) {
              ensemble[j].p_val = p_star; ensemble[j].cached_weights = std::move(workspace_weights); ensemble[j].current_stats = std::move(prop_stats);
            }
          }
          Forest_Preds[j] = ensemble[j].cached_weights * ensemble[j].mu; Y_hat = Y - res + Forest_Preds[j];
        }
        UNPROTECT(1);
      }
      
      // Removed dynamic updating of boost and penalty terms during burn-in
      /*
       if (var_sel_mode == 2 && iter < burn_in && (iter + 1) % tune_window == 0 && jump_count > 0) {
       double avg_acc = running_acc / jump_count; double error = target_acceptance - avg_acc;
       adapt_boost = fmax(0.001, adapt_boost + learning_rate * error); adapt_penalty = fmax(0.001, adapt_penalty - learning_rate * error);
       running_acc = 0.0; jump_count = 0;
       }
       */
      
      if (var_sel_mode > 0 && iter > dirichlet_warmup) {
        std::vector<double> m_counts(p, 0.0);
        for (int j = 0; j < m; ++j) {
          SEXP dim_j_sexp = VECTOR_ELT(return_dim, j); int d_len = Rf_length(dim_j_sexp); const int* p_dim_j = INTEGER(dim_j_sexp);
          for(int k = 0; k < d_len; ++k) m_counts[p_dim_j[k] - 1] += 1.0;
          if (d_len > 1) {
            std::vector<int> ordered_dims(p_dim_j, p_dim_j + d_len);
            for (int k = d_len - 1; k > 0; --k) { int r = floor(unif_rand() * (k + 1)); std::swap(ordered_dims[k], ordered_dims[r]); }
            for (int k = 1; k < d_len; ++k) {
              std::vector<int> ineligible; double sum_inel = 0.0;
              for (int e = 0; e < k; ++e) { ineligible.push_back(ordered_dims[e]); sum_inel += p_s_weights[ordered_dims[e] - 1]; }
              double prob_success = fmax(1e-5, fmin(1.0 - 1e-10, 1.0 - sum_inel)); double raw_failures = R::rgeom(prob_success); int failures = static_cast<int>(fmin(raw_failures, 10000.0));
              if (failures > 0) {
                if (ineligible.size() == 1) { m_counts[ineligible[0] - 1] += failures; } 
                else {
                  int k_classes = ineligible.size(); std::vector<double> probs(k_classes); double sum_w = 0.0;
                  for (int i = 0; i < k_classes; ++i) { probs[i] = p_s_weights[ineligible[i] - 1]; sum_w += probs[i]; } sum_w = fmax(sum_w, 1e-16); double running_sum = 0.0;
                  for (int i = 0; i < k_classes - 1; ++i) { probs[i] /= sum_w; running_sum += probs[i]; } probs[k_classes - 1] = fmax(0.0, 1.0 - running_sum); 
                  std::vector<int> alloc_counts(k_classes, 0); Rf_rmultinom(failures, probs.data(), k_classes, alloc_counts.data());
                  for (int i = 0; i < k_classes; ++i) m_counts[ineligible[i] - 1] += alloc_counts[i];
                }
              }
            }
          }
        }
        double sum_s = 0.0;
        for (int i = 0; i < p; ++i) {
          double raw_shape = (alpha / p) + m_counts[i]; if (var_sel_mode == 2) raw_shape += adaptive_momentum[i];
          p_s_weights[i] = fmax(R::rgamma(fmax(raw_shape, 1e-10), 1.0), 1e-16); sum_s += p_s_weights[i];
        }
        double sum_log_s = 0.0;
        for (int i = 0; i < p; ++i) {
          p_s_weights[i] /= sum_s; sum_log_s += log(fmax(p_s_weights[i], 1e-16));
          if (var_sel_mode == 2) adaptive_momentum[i] *= momentum_decay; 
        }
        double alpha_star = exp(norm_rand() * 0.5 + log(alpha));
        auto eval_log_post_alpha = [&](double a_val) { return R::lgammafn(a_val) - p * R::lgammafn(a_val / p) + (a_val / p) * sum_log_s + (a_alpha - 1.0) * log(a_val) - (a_alpha + b_alpha) * log(a_val + rho_alpha); };
        if (log(unif_rand()) < (eval_log_post_alpha(alpha_star) - eval_log_post_alpha(alpha) + log(alpha_star) - log(alpha))) { alpha = alpha_star; }
      } else {
        for(int i = 0; i < p; ++i) { p_s_weights[i] = 1.0 / p; if (var_sel_mode == 2) adaptive_momentum[i] *= momentum_decay; }
      }
      
      if (iter >= burn_in && (iter + 1 - burn_in) % thinning == 0) {
        if (split_mode == 0) {
          SEXP tl = PROTECT(Rf_allocVector(VECSXP, m)); SEXP dl = PROTECT(Rf_allocVector(VECSXP, m)); SEXP pl = PROTECT(Rf_allocVector(VECSXP, m)); SEXP ml = PROTECT(Rf_allocVector(VECSXP, m));
          for(int j = 0; j < m; ++j) {
            SEXP out_t = PROTECT(Rf_allocMatrix(REALSXP, ensemble[j].n_centres, ensemble[j].n_dims)); std::memcpy(REAL(out_t), ensemble[j].tess.data(), ensemble[j].tess.size() * sizeof(double));
            SEXP out_d = PROTECT(Rf_allocVector(INTSXP, ensemble[j].n_dims)); std::memcpy(INTEGER(out_d), ensemble[j].dim.data(), ensemble[j].dim.size() * sizeof(int));
            SEXP out_p = PROTECT(Rf_allocMatrix(REALSXP, n_obs, 1)); std::memcpy(REAL(out_p), Forest_Preds[j].data(), n_obs * sizeof(double));
            SEXP out_m = PROTECT(Rf_allocMatrix(REALSXP, ensemble[j].n_centres, 1)); std::memcpy(REAL(out_m), ensemble[j].mu.data(), ensemble[j].n_centres * sizeof(double));
            SET_VECTOR_ELT(tl, j, out_t); SET_VECTOR_ELT(dl, j, out_d); SET_VECTOR_ELT(pl, j, out_p); SET_VECTOR_ELT(ml, j, out_m); UNPROTECT(4); REAL(out_power)[j + store_idx * m] = ensemble[j].p_val;
          }
          SET_VECTOR_ELT(out_tess, store_idx, tl); SET_VECTOR_ELT(out_dim, store_idx, dl); SET_VECTOR_ELT(out_pred, store_idx, pl); SET_VECTOR_ELT(out_mu, store_idx, ml); UNPROTECT(4);
        } else {
          SET_VECTOR_ELT(out_tess, store_idx, Rf_duplicate(return_tess)); SET_VECTOR_ELT(out_dim, store_idx, Rf_duplicate(return_dim)); SET_VECTOR_ELT(out_pred, store_idx, Rf_duplicate(return_pred));
        }
        REAL(out_sigma)[store_idx] = hypers.sigma2;
        if (split_mode == 1) std::memcpy(p_out_mat + store_idx * n_obs, p_sum, n_obs * sizeof(double)); else std::memcpy(p_out_mat + store_idx * n_obs, Y_hat.data(), n_obs * sizeof(double));
        for(int i = 0; i < p; ++i) REAL(out_dirichlet)[i + store_idx * p] = p_s_weights[i];
        REAL(out_alpha)[store_idx] = alpha;
        std::vector<int> current_covariate_counts(p, 0);
        for (int j = 0; j < m; ++j) { SEXP dim_j_sexp = VECTOR_ELT(return_dim, j); int d_len = Rf_length(dim_j_sexp); const int* p_dim_j = INTEGER(dim_j_sexp); for(int k = 0; k < d_len; ++k) current_covariate_counts[p_dim_j[k] - 1] += 1; }
        for(int i = 0; i < p; ++i) INTEGER(out_var_sel)[i + store_idx * p] = current_covariate_counts[i];
        store_idx++;
      }
    }
    
    if (show_progress) { Rprintf("\n"); R_FlushConsole(); } PutRNGstate();
    
    SEXP return_list, list_names;
    PROTECT(return_list = Rf_allocVector(VECSXP, 10));
    SET_VECTOR_ELT(return_list, 0, out_tess); SET_VECTOR_ELT(return_list, 1, out_dim); SET_VECTOR_ELT(return_list, 2, out_pred); SET_VECTOR_ELT(return_list, 3, out_sigma); SET_VECTOR_ELT(return_list, 4, out_pred_mat); SET_VECTOR_ELT(return_list, 5, out_dirichlet); SET_VECTOR_ELT(return_list, 6, out_var_sel); SET_VECTOR_ELT(return_list, 7, out_alpha);
    SET_VECTOR_ELT(return_list, 8, out_power); SET_VECTOR_ELT(return_list, 9, out_mu);
    
    PROTECT(list_names = Rf_allocVector(STRSXP, 10));
    SET_STRING_ELT(list_names, 0, Rf_mkChar("posteriorTess")); SET_STRING_ELT(list_names, 1, Rf_mkChar("posteriorDim")); SET_STRING_ELT(list_names, 2, Rf_mkChar("posteriorPred")); SET_STRING_ELT(list_names, 3, Rf_mkChar("posteriorSigma")); SET_STRING_ELT(list_names, 4, Rf_mkChar("predictionMatrix")); SET_STRING_ELT(list_names, 5, Rf_mkChar("posteriorDirichletWeights")); SET_STRING_ELT(list_names, 6, Rf_mkChar("posteriorVariableSelection")); SET_STRING_ELT(list_names, 7, Rf_mkChar("posteriorAlpha"));
    SET_STRING_ELT(list_names, 8, Rf_mkChar("posteriorPower")); SET_STRING_ELT(list_names, 9, Rf_mkChar("posteriorMu"));
    
    Rf_setAttrib(return_list, R_NamesSymbol, list_names);
    
    UNPROTECT(19); return return_list;
  }
  // ---------------------------------------------------------
  // PREDICTION: SOFT IDW
  // ---------------------------------------------------------
  SEXP soft_predict_cpp(SEXP tess_sexp, SEXP query_sexp, SEXP dim_sexp, SEXP mu_sexp, SEXP power_sexp) {
    const double* p_tess = REAL(tess_sexp); const double* p_query = REAL(query_sexp);
    const int* dim_p = INTEGER(dim_sexp); const double* p_mu = REAL(mu_sexp);
    double p_val = REAL(power_sexp)[0];
    
    int tess_rows = Rf_nrows(tess_sexp); int query_rows = Rf_nrows(query_sexp); int n_dims = Rf_length(dim_sexp);
    
    std::vector<double> dists(query_rows * tess_rows, 0.0);
    for (int d = 0; d < n_dims; ++d) {
      int cov_idx = dim_p[d] - 1; int query_offset = cov_idx * query_rows; int tess_offset = d * tess_rows;
      for (int c = 0; c < tess_rows; ++c) {
        double tess_val = p_tess[c + tess_offset]; int dist_offset = c * query_rows;
        for (int q = 0; q < query_rows; ++q) { double diff = p_query[q + query_offset] - tess_val; dists[q + dist_offset] += diff * diff; }
      }
    }
    
    SEXP result = PROTECT(Rf_allocMatrix(REALSXP, query_rows, 1)); double* p_result = REAL(result);
    double half_power = p_val / 2.0; double epsilon = 1e-10; double eps_sq = epsilon * epsilon;
    bool fast_path = std::abs(half_power - 1.0) < 1e-7;
    double dim_scale = 1.0 / static_cast<double>(n_dims);
    
    for (int q = 0; q < query_rows; ++q) {
      bool is_zero = false; int z_idx = -1;
      for (int c = 0; c < tess_rows; ++c) { if (dists[q + c * query_rows] * dim_scale < eps_sq) { is_zero = true; z_idx = c; break; } }
      if (is_zero) {
        p_result[q] = p_mu[z_idx];
      } else {
        double rsum = 0.0; double pred_val = 0.0;
        for (int c = 0; c < tess_rows; ++c) {
          double sqd = dists[q + c * query_rows] * dim_scale;
          double w = fast_path ? 1.0 / (sqd + epsilon) : 1.0 / (std::pow(sqd, half_power) + epsilon);
          pred_val += w * p_mu[c]; rsum += w;
        }
        p_result[q] = pred_val / rsum;
      }
    }
    UNPROTECT(1); return result;
  }
  
  // ---------------------------------------------------------
  // PREDICTION: HARD VORONOI
  // ---------------------------------------------------------
  SEXP knnx_index_predict_cpp(SEXP tess_sexp, SEXP query_sexp, SEXP k_sexp, SEXP dim_sexp) {
    const double* p_tess = REAL(tess_sexp); const double* p_query = REAL(query_sexp);
    int k = INTEGER(k_sexp)[0]; const int* dim_p = INTEGER(dim_sexp);
    int tess_rows = Rf_nrows(tess_sexp); int query_rows = Rf_nrows(query_sexp); int n_dims = Rf_length(dim_sexp);
    
    if (k <= 0 || k > tess_rows) { Rf_error("k must be positive and not greater than number of reference points"); }
    std::vector<double> dists(query_rows * tess_rows, 0.0);
    for (int d = 0; d < n_dims; ++d) {
      int cov_idx = dim_p[d] - 1; int query_offset = cov_idx * query_rows; int tess_offset = d * tess_rows;
      for (int c = 0; c < tess_rows; ++c) {
        double tess_val = p_tess[c + tess_offset]; int dist_offset = c * query_rows;
        for (int q = 0; q < query_rows; ++q) { double diff = p_query[q + query_offset] - tess_val; dists[q + dist_offset] += diff * diff; }
      }
    }
    
    SEXP result; PROTECT(result = Rf_allocMatrix(INTSXP, query_rows, k)); int* p_result = INTEGER(result);
    if (k == 1) {
      std::vector<double> current_min_dist(dists.begin(), dists.begin() + query_rows);
      for (int q = 0; q < query_rows; ++q) p_result[q] = 1;
      for (int c = 1; c < tess_rows; ++c) {
        int c_offset = c * query_rows;
        for (int q = 0; q < query_rows; ++q) { double dist = dists[q + c_offset]; if (dist < current_min_dist[q]) { current_min_dist[q] = dist; p_result[q] = c + 1; } }
      }
    } else {
      std::vector<std::pair<double, int>> distances(tess_rows);
      for (int q = 0; q < query_rows; ++q) {
        for (int c = 0; c < tess_rows; ++c) { distances[c].first = dists[q + c * query_rows]; distances[c].second = c + 1; }
        std::partial_sort(distances.begin(), distances.begin() + k, distances.end());
        for (int i = 0; i < k; ++i) { p_result[q + i * query_rows] = distances[i].second; }
      }
    }
    UNPROTECT(1); return result;
  }
}