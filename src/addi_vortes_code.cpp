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
    double sigma2;
    double sigma2mu;
    double omega;
    double lambda_poisson;
    double p_shape;
    double p_rate;
    double p_sd;
    double nu;
    double lambda_invgamma;
    int num_covariates;
    
    void UpdateSigma(const Eigen::VectorXd& res) {
      double sse = res.squaredNorm();
      double n = res.size();
      double shape = (nu + n) / 2.0;
      double rate = (nu * lambda_invgamma + sse) / 2.0;
      sigma2 = 1.0 / R::rgamma(shape, 1.0 / rate);
    }
  };
  
  struct PosteriorStats {
    Eigen::VectorXd mu_post;
    Eigen::MatrixXd L_chol;
    double quad_form;
    double log_det_Q;
    bool valid;
  };
  
  // Pure data container, no internal allocations
  struct Tessellation {
    std::vector<double> tess;
    std::vector<int> dim;
    int n_centres;
    int n_dims;
    double p_val;
    Eigen::VectorXd mu;
    
    Eigen::MatrixXd cached_sqdist;
    Eigen::MatrixXd cached_weights;
    PosteriorStats current_stats;
  };
  
  double calc_acceptance_cpp(const PosteriorStats& StatsOld, const PosteriorStats& StatsNew,
                             int old_c, int new_c, int old_d, int new_d,
                             const std::string& mod, const Hypers& hypers, int p) {
    
    double prob_eps = 1e-10;
    double prob = fmin(1.0 - prob_eps, fmax(0.0, hypers.omega / p));
    double LogLikelihoodRatio = 0.5 * (StatsOld.log_det_Q - StatsNew.log_det_Q) + 0.5 * (StatsNew.quad_form - StatsOld.quad_form);    
    double ratio = LogLikelihoodRatio;
    
    if (mod == "AD") {
      ratio += log(R::dbinom(new_d - 1, p - 1, prob, 0) / R::dbinom(new_d - 2, p - 1, prob, 0)) + log((double)(p - new_d + 1) / new_d);
      if (old_d == 1) ratio += log(0.5); else if (old_d == p - 1) ratio += log(2.0);
    } else if (mod == "RD") {
      ratio += log(R::dbinom(new_d - 1, p, prob, 0) / R::dbinom(new_d, p, prob, 0)) + log((double)(new_d + 1) / (p - new_d));
      if (old_d == p) ratio += log(0.5); else if (old_d == 2) ratio += log(2.0);
    } else if (mod == "AC") {
      ratio += log(R::dpois(new_c - 1, hypers.lambda_poisson, 0) / R::dpois(new_c - 2, hypers.lambda_poisson, 0)) + log(1.0 / new_c);
      if (old_c == 1) ratio += log(0.5);
    } else if (mod == "RC") {
      ratio += log(R::dpois(new_c - 1, hypers.lambda_poisson, 0) / R::dpois(new_c, hypers.lambda_poisson, 0)) + log(new_c + 1.0);
      if (old_c == 2) ratio += log(2.0);
    }
    return ratio;
  }
  
  double sample_sigma_squared_cpp(int n, double nu, double lambda, const double* __restrict__ p_y, const double* __restrict__ p_sum) {
    double sse = 0.0;
    for (int i = 0; i < n; ++i) {
      double diff = p_y[i] - p_sum[i];
      sse += diff * diff;
    }
    return 1.0 / R::rgamma((nu + n) / 2.0, 1.0 / ((nu * lambda + sse) / 2.0));
  }
  
  // Memory builder functions (NRVO ensures zero-copy return)
  Eigen::MatrixXd build_sqdist(const double* p_data, int n_obs, 
                               const std::vector<double>& old_tess, int old_c, 
                               const std::vector<double>& new_tess, int new_c,
                               const std::vector<int>& old_dim, int old_d,
                               const std::vector<int>& new_dim, int new_d,
                               const Eigen::MatrixXd& old_sqdist,
                               const std::string& mod, int mod_idx) {
    
    Eigen::MatrixXd new_sq(n_obs, new_c);
    double* p_new = new_sq.data();
    const double* p_old = old_sqdist.data();
    
    if (mod == "RC") {
      if (mod_idx > 0) std::memcpy(p_new, p_old, mod_idx * n_obs * sizeof(double));
      if (mod_idx < old_c - 1) {
        std::memcpy(p_new + (mod_idx * n_obs), p_old + ((mod_idx + 1) * n_obs), ((old_c - 1) - mod_idx) * n_obs * sizeof(double));
      }
    } else if (mod == "AC") {
      std::memcpy(p_new, p_old, old_c * n_obs * sizeof(double));
      int c = new_c - 1; 
      for (int i = 0; i < n_obs; ++i) p_new[i + c * n_obs] = 0.0;
      for (int d = 0; d < new_d; ++d) {
        int offset = (new_dim[d] - 1) * n_obs;
        double t_val = new_tess[c + d * new_c];
        for (int i = 0; i < n_obs; ++i) {
          double diff = p_data[i + offset] - t_val;
          p_new[i + c * n_obs] += diff * diff;
        }
      }
    } else if (mod == "Change") {
      std::memcpy(p_new, p_old, new_c * n_obs * sizeof(double));
      int c = mod_idx;
      for (int i = 0; i < n_obs; ++i) p_new[i + c * n_obs] = 0.0;
      for (int d = 0; d < new_d; ++d) {
        int offset = (new_dim[d] - 1) * n_obs;
        double t_val = new_tess[c + d * new_c];
        for (int i = 0; i < n_obs; ++i) {
          double diff = p_data[i + offset] - t_val;
          p_new[i + c * n_obs] += diff * diff;
        }
      }
    } else if (mod == "AD") {
      int offset = (new_dim[mod_idx] - 1) * n_obs;
      for (int c = 0; c < new_c; ++c) {
        double t_val = new_tess[c + mod_idx * new_c];
        int c_offset = c * n_obs;
        for (int i = 0; i < n_obs; ++i) {
          double diff = p_data[i + offset] - t_val;
          p_new[i + c_offset] = p_old[i + c_offset] + (diff * diff);
        }
      }
    } else if (mod == "RD") {
      int offset = (old_dim[mod_idx] - 1) * n_obs;
      for (int c = 0; c < new_c; ++c) {
        double t_val = old_tess[c + mod_idx * old_c];
        int c_offset = c * n_obs;
        for (int i = 0; i < n_obs; ++i) {
          double diff = p_data[i + offset] - t_val;
          double val = p_old[i + c_offset] - (diff * diff);
          p_new[i + c_offset] = (val < 0.0) ? 0.0 : val;
        }
      }
    } else if (mod == "Swap") {
      int old_off = (old_dim[mod_idx] - 1) * n_obs;
      int new_off = (new_dim[mod_idx] - 1) * n_obs;
      for (int c = 0; c < new_c; ++c) {
        double old_t = old_tess[c + mod_idx * old_c];
        double new_t = new_tess[c + mod_idx * new_c];
        int c_off = c * n_obs;
        for (int i = 0; i < n_obs; ++i) {
          double o_diff = p_data[i + old_off] - old_t;
          double n_diff = p_data[i + new_off] - new_t;
          double val = p_old[i + c_off] - (o_diff * o_diff) + (n_diff * n_diff);
          p_new[i + c_off] = (val < 0.0) ? 0.0 : val;
        }
      }
    } else { // Init
      for (int c = 0; c < new_c; ++c) {
        for (int i = 0; i < n_obs; ++i) {
          double dval = 0.0;
          for (int d = 0; d < new_d; ++d) {
            double diff = p_data[i + (new_dim[d] - 1) * n_obs] - new_tess[c + d * new_c];
            dval += diff * diff;
          }
          p_new[i + c * n_obs] = dval;
        }
      }
    }
    return new_sq;
  }
  
  Eigen::MatrixXd build_weights(const Eigen::MatrixXd& sqdist, double p_val, double epsilon = 1e-10) {
    int n_obs = sqdist.rows();
    int n_centres = sqdist.cols();
    Eigen::MatrixXd weights(n_obs, n_centres);
    double half_power = p_val / 2.0;
    double eps_sq = epsilon * epsilon;
    bool fast_path = std::abs(half_power - 1.0) < 1e-7;
    
    for (int i = 0; i < n_obs; ++i) {
      bool is_zero = false;
      int z_idx = -1;
      for (int c = 0; c < n_centres; ++c) {
        if (sqdist(i, c) < eps_sq) { is_zero = true; z_idx = c; break; }
      }
      if (is_zero) {
        for (int c = 0; c < n_centres; ++c) weights(i, c) = (c == z_idx) ? 1.0 : 0.0;
      } else {
        double rsum = 0.0;
        for (int c = 0; c < n_centres; ++c) {
          double w = fast_path ? 1.0 / (sqdist(i, c) + epsilon) : 1.0 / (std::pow(sqdist(i, c), half_power) + epsilon);
          weights(i, c) = w;
          rsum += w;
        }
        double inv_sum = 1.0 / rsum;
        for (int c = 0; c < n_centres; ++c) weights(i, c) *= inv_sum;
      }
    }
    return weights;
  }
  
  PosteriorStats build_stats(const Eigen::MatrixXd& w, const Eigen::VectorXd& res, const Hypers& hypers) {
    int k = w.cols();
    PosteriorStats stats;
    Eigen::MatrixXd Q(k, k);
    Q.noalias() = w.transpose() * w; 
    Q /= hypers.sigma2;
    Q.diagonal().array() += (1.0 / hypers.sigma2mu);
    
    Eigen::LLT<Eigen::MatrixXd> llt(Q);
    if (llt.info() != Eigen::NumericalIssue) {
      stats.L_chol = llt.matrixL();
      Eigen::VectorXd z(k);
      z.noalias() = (1.0 / hypers.sigma2) * w.transpose() * res;
      stats.mu_post = llt.solve(z);
      stats.quad_form = stats.mu_post.dot(z);
      stats.log_det_Q = 2.0 * stats.L_chol.diagonal().array().log().sum();
      stats.valid = true;
    } else {
      stats.valid = false;
    }
    return stats;
  }
  
  void update_mu(Tessellation& t) {
    Eigen::VectorXd z(t.n_centres);
    for(int i = 0; i < t.n_centres; ++i) z(i) = norm_rand();
    t.mu = t.current_stats.mu_post + t.current_stats.L_chol.triangularView<Eigen::Lower>().transpose().solve(z);
  }
  
  void update_p(Tessellation& t, const Eigen::VectorXd& res, const Hypers& hypers) {
    double log_p_star = log(t.p_val) + norm_rand() * hypers.p_sd;
    double p_star = exp(log_p_star);
    
    Eigen::MatrixXd prop_weights = build_weights(t.cached_sqdist, p_star);
    PosteriorStats prop_stats = build_stats(prop_weights, res, hypers);
    
    if (prop_stats.valid && t.current_stats.valid) {
      double ll_old = 0.5 * t.current_stats.quad_form - 0.5 * t.current_stats.log_det_Q;
      double ll_new = 0.5 * prop_stats.quad_form - 0.5 * prop_stats.log_det_Q;
      double prior_old = R::dexp(t.p_val, 1.0 / hypers.p_rate, 1);
      double prior_new = R::dexp(p_star, 1.0 / hypers.p_rate, 1);
      
      if (log(unif_rand()) < (ll_new - ll_old + prior_new - prior_old + log_p_star - log(t.p_val))) {
        t.p_val = p_star;
        t.cached_weights = std::move(prop_weights);
        t.current_stats = std::move(prop_stats);
      }
    }
  }
  
  void perturb_tessellation(Tessellation& t, const double* p_data, int n_obs, const Eigen::VectorXd& res, const Hypers& hypers, double sd_val, const double* mu_val) {
    std::string mod = "Change";
    int mod_idx = 0;
    double p = unif_rand();
    
    std::vector<double> prop_tess = t.tess;
    std::vector<int> prop_dim = t.dim;
    int prop_c = t.n_centres;
    int prop_d = t.n_dims;
    
    if ((p < 0.2 && prop_d != hypers.num_covariates) || (prop_d == 1 && p < 0.4)) {
      mod = "AD";
      int new_d; do { new_d = floor(unif_rand() * hypers.num_covariates) + 1; } while (in_vector(new_d, prop_dim));
      prop_dim.push_back(new_d);
      mod_idx = prop_d; 
      std::vector<double> n_tess(prop_c * (prop_d + 1));
      for (int r = 0; r < prop_c; ++r) {
        for (int c = 0; c < prop_d; ++c) n_tess[r + c * prop_c] = prop_tess[r + c * prop_c];
        n_tess[r + prop_d * prop_c] = mu_val[new_d-1] + norm_rand() * sd_val;
      }
      prop_tess = std::move(n_tess); prop_d++;
    } else if (p < 0.4 && prop_d > 1) {
      mod = "RD"; mod_idx = floor(unif_rand() * prop_d);
      prop_dim.erase(prop_dim.begin() + mod_idx);
      std::vector<double> n_tess(prop_c * (prop_d - 1));
      int cur_c = 0;
      for (int c = 0; c < prop_d; ++c) {
        if (c != mod_idx) {
          for (int r = 0; r < prop_c; ++r) n_tess[r + cur_c * prop_c] = prop_tess[r + c * prop_c];
          cur_c++;
        }
      }
      prop_tess = std::move(n_tess); prop_d--;
    } else if (p < 0.6 || (p < 0.8 && prop_c == 1)) {
      mod = "AC"; mod_idx = prop_c;
      for (int i = 0; i < prop_d; ++i) prop_tess.insert(prop_tess.begin() + (i * (prop_c + 1)) + prop_c, mu_val[i] + norm_rand() * sd_val);
      prop_c++;
    } else if (p < 0.8 && prop_c > 1) {
      mod = "RC"; mod_idx = floor(unif_rand() * prop_c);
      std::vector<double> n_tess; n_tess.reserve((prop_c - 1) * prop_d);
      for(int c = 0; c < prop_d; ++c) {
        for(int r = 0; r < prop_c; ++r) {
          if(r != mod_idx) n_tess.push_back(prop_tess[r + c * prop_c]);
        }
      }
      prop_tess = std::move(n_tess); prop_c--;
    } else if (p < 0.9 || prop_d == hypers.num_covariates) {
      mod = "Change"; mod_idx = floor(unif_rand() * prop_c);
      for (int c = 0; c < prop_d; ++c) prop_tess[mod_idx + c * prop_c] = mu_val[c] + norm_rand() * sd_val;
    } else {
      mod = "Swap"; mod_idx = floor(unif_rand() * prop_d);
      int new_d; do { new_d = floor(unif_rand() * hypers.num_covariates) + 1; } while (in_vector(new_d, prop_dim));
      prop_dim[mod_idx] = new_d;
      for (int r = 0; r < prop_c; ++r) prop_tess[r + mod_idx * prop_c] = mu_val[mod_idx] + norm_rand() * sd_val;
    }
    
    // Generate the Ghost Matrices
    Eigen::MatrixXd prop_sqdist = build_sqdist(p_data, n_obs, t.tess, t.n_centres, prop_tess, prop_c, t.dim, t.n_dims, prop_dim, prop_d, t.cached_sqdist, mod, mod_idx);
    Eigen::MatrixXd prop_weights = build_weights(prop_sqdist, t.p_val);
    PosteriorStats prop_stats = build_stats(prop_weights, res, hypers);
    
    if (prop_stats.valid && t.current_stats.valid) {
      double ll_old = 0.5 * t.current_stats.quad_form - 0.5 * t.current_stats.log_det_Q;
      double ll_new = 0.5 * prop_stats.quad_form - 0.5 * prop_stats.log_det_Q;
      double ratio = calc_acceptance_cpp(t.current_stats, prop_stats, t.n_centres, prop_c, t.n_dims, prop_d, mod, hypers, hypers.num_covariates);
      
      if (log(unif_rand()) < (ll_new - ll_old + ratio)) {
        // INSTANT ZERO-COPY POINTER SWAP
        t.tess = std::move(prop_tess);
        t.dim = std::move(prop_dim);
        t.n_centres = prop_c;
        t.n_dims = prop_d;
        t.cached_sqdist = std::move(prop_sqdist);
        t.cached_weights = std::move(prop_weights);
        t.current_stats = std::move(prop_stats);
      }
    }
  }
  
} // End Namespace

extern "C" {
  
  SEXP super_mcmc_loop_cpp(SEXP x_sexp, SEXP y_sexp, SEXP sum_sexp, SEXP tess_list, SEXP dim_list, 
                           SEXP pred_list, SEXP m_sexp, SEXP p_sexp, 
                           SEXP sd_sexp, SEXP mu_sexp, SEXP sigma2mu_sexp, SEXP omega_sexp, SEXP poisson_lambda_sexp,
                           SEXP total_iter_sexp, SEXP burn_in_sexp, SEXP thinning_sexp,
                           SEXP nu_sexp, SEXP invgamma_lambda_sexp, SEXP show_progress_sexp, 
                           SEXP power_sexp, SEXP p_shape_sexp, SEXP p_rate_sexp, SEXP p_sd_sexp) {
    
    Hypers hypers;
    hypers.num_covariates = INTEGER(p_sexp)[0];
    hypers.sigma2mu = REAL(sigma2mu_sexp)[0];
    hypers.omega = REAL(omega_sexp)[0];
    hypers.lambda_poisson = REAL(poisson_lambda_sexp)[0];
    hypers.p_shape = REAL(p_shape_sexp)[0];
    hypers.p_rate = REAL(p_rate_sexp)[0];
    hypers.p_sd = REAL(p_sd_sexp)[0];
    hypers.nu = REAL(nu_sexp)[0];
    hypers.lambda_invgamma = REAL(invgamma_lambda_sexp)[0];
    
    int m = INTEGER(m_sexp)[0];
    int n_obs = Rf_length(y_sexp);
    const double* p_x = REAL(x_sexp);
    double sd_val = REAL(sd_sexp)[0];
    const double* mu_val = REAL(mu_sexp);
    double* p_vec = REAL(power_sexp);
    
    int total_iter = INTEGER(total_iter_sexp)[0];
    int burn_in = INTEGER(burn_in_sexp)[0];
    int thinning = INTEGER(thinning_sexp)[0];
    int show_progress = INTEGER(show_progress_sexp)[0];
    int num_stored = std::max(0, (total_iter - burn_in) / thinning);
    
    // The Ensemble is now a contiguous block in memory!
    std::vector<Tessellation> ensemble(m);
    for(int j = 0; j < m; ++j) {
      SEXP rtess = VECTOR_ELT(tess_list, j);
      SEXP rdim = VECTOR_ELT(dim_list, j);
      ensemble[j].n_centres = Rf_nrows(rtess);
      ensemble[j].n_dims = Rf_length(rdim);
      ensemble[j].tess.assign(REAL(rtess), REAL(rtess) + (ensemble[j].n_centres * ensemble[j].n_dims));
      ensemble[j].dim.assign(INTEGER(rdim), INTEGER(rdim) + ensemble[j].n_dims);
      ensemble[j].p_val = p_vec[j];
      
      Eigen::MatrixXd empty_sq;
      ensemble[j].cached_sqdist = build_sqdist(p_x, n_obs, ensemble[j].tess, ensemble[j].n_centres, ensemble[j].tess, ensemble[j].n_centres, 
                                               ensemble[j].dim, ensemble[j].n_dims, ensemble[j].dim, ensemble[j].n_dims, empty_sq, "Init", 0);
      ensemble[j].cached_weights = build_weights(ensemble[j].cached_sqdist, ensemble[j].p_val);
      ensemble[j].mu = Eigen::VectorXd::Zero(ensemble[j].n_centres); 
    }
    
    Eigen::Map<const Eigen::VectorXd> Y(REAL(y_sexp), n_obs);
    Eigen::VectorXd Y_hat(n_obs);
    std::memcpy(Y_hat.data(), REAL(sum_sexp), n_obs * sizeof(double));
    std::vector<Eigen::VectorXd> Forest_Preds(m);
    for(int j = 0; j < m; ++j) Forest_Preds[j] = Eigen::Map<const Eigen::VectorXd>(REAL(VECTOR_ELT(pred_list, j)), n_obs);
    
    SEXP out_tess = PROTECT(Rf_allocVector(VECSXP, num_stored));
    SEXP out_dim = PROTECT(Rf_allocVector(VECSXP, num_stored));
    SEXP out_pred = PROTECT(Rf_allocVector(VECSXP, num_stored));
    SEXP out_sigma = PROTECT(Rf_allocVector(REALSXP, num_stored));
    SEXP out_power = PROTECT(Rf_allocMatrix(REALSXP, m, num_stored));
    SEXP out_mu = PROTECT(Rf_allocVector(VECSXP, num_stored));
    SEXP out_pred_mat = PROTECT(Rf_allocMatrix(REALSXP, n_obs, num_stored));
    
    GetRNGstate();
    int store_idx = 0;
    
    for (int iter = 0; iter < total_iter; ++iter) {
      if (show_progress && (iter + 1) % std::max(1, total_iter / 100) == 0) {
        Rprintf("\rProcessing MCMC Iteration %d of %d...", iter + 1, total_iter);
        R_FlushConsole(); R_CheckUserInterrupt();
      }
      
      hypers.UpdateSigma(Y - Y_hat);
      
      for(int j = 0; j < m; ++j) {
        Eigen::VectorXd res = Y - (Y_hat - Forest_Preds[j]);
        ensemble[j].current_stats = build_stats(ensemble[j].cached_weights, res, hypers);
        
        perturb_tessellation(ensemble[j], p_x, n_obs, res, hypers, sd_val, mu_val);
        update_p(ensemble[j], res, hypers);
        update_mu(ensemble[j]);
        
        Forest_Preds[j] = ensemble[j].cached_weights * ensemble[j].mu;
        Y_hat = Y - res + Forest_Preds[j];
      }
      
      if (iter >= burn_in && (iter + 1 - burn_in) % thinning == 0) {
        SEXP tl = PROTECT(Rf_allocVector(VECSXP, m));
        SEXP dl = PROTECT(Rf_allocVector(VECSXP, m));
        SEXP pl = PROTECT(Rf_allocVector(VECSXP, m));
        SEXP ml = PROTECT(Rf_allocVector(VECSXP, m));
        
        for(int j = 0; j < m; ++j) {
          SEXP out_t = PROTECT(Rf_allocMatrix(REALSXP, ensemble[j].n_centres, ensemble[j].n_dims));
          std::memcpy(REAL(out_t), ensemble[j].tess.data(), ensemble[j].tess.size() * sizeof(double));
          SEXP out_d = PROTECT(Rf_allocVector(INTSXP, ensemble[j].n_dims));
          std::memcpy(INTEGER(out_d), ensemble[j].dim.data(), ensemble[j].dim.size() * sizeof(int));
          
          SEXP out_p = PROTECT(Rf_allocMatrix(REALSXP, n_obs, 1));
          std::memcpy(REAL(out_p), Forest_Preds[j].data(), n_obs * sizeof(double));
          SEXP out_m = PROTECT(Rf_allocMatrix(REALSXP, ensemble[j].n_centres, 1));
          std::memcpy(REAL(out_m), ensemble[j].mu.data(), ensemble[j].n_centres * sizeof(double));
          
          SET_VECTOR_ELT(tl, j, out_t); SET_VECTOR_ELT(dl, j, out_d);
          SET_VECTOR_ELT(pl, j, out_p); SET_VECTOR_ELT(ml, j, out_m);
          UNPROTECT(4);
          REAL(out_power)[j + store_idx * m] = ensemble[j].p_val;
        }
        
        SET_VECTOR_ELT(out_tess, store_idx, tl); SET_VECTOR_ELT(out_dim, store_idx, dl);
        SET_VECTOR_ELT(out_pred, store_idx, pl); SET_VECTOR_ELT(out_mu, store_idx, ml);
        UNPROTECT(4);
        
        REAL(out_sigma)[store_idx] = hypers.sigma2;
        std::memcpy(REAL(out_pred_mat) + store_idx * n_obs, Y_hat.data(), n_obs * sizeof(double));
        store_idx++;
      }
    }
    
    if (show_progress) { Rprintf("\n"); R_FlushConsole(); }
    PutRNGstate();
    
    SEXP return_list, list_names;
    PROTECT(return_list = Rf_allocVector(VECSXP, 7));
    SET_VECTOR_ELT(return_list, 0, out_tess); SET_VECTOR_ELT(return_list, 1, out_dim);
    SET_VECTOR_ELT(return_list, 2, out_pred); SET_VECTOR_ELT(return_list, 3, out_sigma);
    SET_VECTOR_ELT(return_list, 4, out_power); SET_VECTOR_ELT(return_list, 5, out_mu);
    SET_VECTOR_ELT(return_list, 6, out_pred_mat);
    
    PROTECT(list_names = Rf_allocVector(STRSXP, 7));
    SET_STRING_ELT(list_names, 0, Rf_mkChar("posteriorTess")); SET_STRING_ELT(list_names, 1, Rf_mkChar("posteriorDim"));
    SET_STRING_ELT(list_names, 2, Rf_mkChar("posteriorPred")); SET_STRING_ELT(list_names, 3, Rf_mkChar("posteriorSigma"));
    SET_STRING_ELT(list_names, 4, Rf_mkChar("posteriorPower")); SET_STRING_ELT(list_names, 5, Rf_mkChar("posteriorMu"));
    SET_STRING_ELT(list_names, 6, Rf_mkChar("predictionMatrix"));
    Rf_setAttrib(return_list, R_NamesSymbol, list_names);
    
    UNPROTECT(9);
    return return_list;
  }
  
  SEXP soft_predict_cpp(SEXP tess_sexp, SEXP query_sexp, SEXP dim_sexp, SEXP mu_sexp, SEXP power_sexp) {
    Tessellation t;
    t.n_centres = Rf_nrows(tess_sexp); t.n_dims = Rf_length(dim_sexp);
    t.tess.assign(REAL(tess_sexp), REAL(tess_sexp) + (t.n_centres * t.n_dims));
    t.dim.assign(INTEGER(dim_sexp), INTEGER(dim_sexp) + t.n_dims);
    
    int query_rows = Rf_nrows(query_sexp);
    Eigen::MatrixXd empty_sq;
    Eigen::MatrixXd SqDist = build_sqdist(REAL(query_sexp), query_rows, t.tess, t.n_centres, t.tess, t.n_centres, t.dim, t.n_dims, t.dim, t.n_dims, empty_sq, "Init", 0);
    Eigen::MatrixXd Weights = build_weights(SqDist, REAL(power_sexp)[0]);
    
    Eigen::Map<const Eigen::VectorXd> mu_vec(REAL(mu_sexp), t.n_centres);
    Eigen::VectorXd preds = Weights * mu_vec;
    
    SEXP result = PROTECT(Rf_allocMatrix(REALSXP, query_rows, 1));
    std::memcpy(REAL(result), preds.data(), query_rows * sizeof(double));
    UNPROTECT(1); return result;
  }
}