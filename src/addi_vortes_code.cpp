#include <Rcpp.h>
// [[Rcpp::depends(Rcpp)]]
using namespace Rcpp;

// [[Rcpp::export]]
List calculate_residuals_cpp(NumericVector R_j, IntegerVector indexes,
                             IntegerVector indexesStar, int num_levels_old,
                             int num_centres_new) {
  // For the 'Old' tessellation
  NumericVector R_ijOld(num_levels_old, 0.0);
  IntegerVector n_ijOld(num_levels_old, 0);
  
  for(int i = 0; i < R_j.length(); ++i) {
    int idx = indexes[i] - 1; // R is 1-based, C++ is 0-based
    if (idx >= 0 && idx < num_levels_old) {
      R_ijOld[idx] += R_j[i];
      n_ijOld[idx]++;
    }
  }
  
  // For the 'New' (star) tessellation
  NumericVector R_ijNew(num_centres_new, 0.0);
  IntegerVector n_ijNew(num_centres_new, 0);
  
  for(int i = 0; i < R_j.length(); ++i) {
    int idx = indexesStar[i] - 1; // R is 1-based, C++ is 0-based
    if (idx >= 0 && idx < num_centres_new) {
      R_ijNew[idx] += R_j[i];
      n_ijNew[idx]++;
    }
  }
  
  return List::create(
    _["R_ijOld"] = R_ijOld,
    _["n_ijOld"] = n_ijOld,
    _["R_ijNew"] = R_ijNew,
    _["n_ijNew"] = n_ijNew
  );
}

// [[Rcpp::export]]
List propose_tessellation_cpp(NumericMatrix tess_j, IntegerVector dim_j, double var, int numCovariates) {
  double p = R::runif(0, 1);
  int d_j_length = dim_j.size();
  int tess_j_rows = tess_j.nrow();
  
  IntegerVector dim_j_star = clone(dim_j);
  NumericMatrix tess_j_star = clone(tess_j);
  std::string modification = "Change";
  
  if ((p < 0.2 && d_j_length != numCovariates) || (d_j_length == 1 && p < 0.4)) {
    modification = "AD";
    int new_dim;
    do {
      new_dim = floor(R::runif(0, 1) * numCovariates) + 1;
    } while (is_true(any(dim_j == new_dim)));
    
    dim_j_star.push_back(new_dim);
    NumericVector new_col = Rcpp::rnorm(tess_j_rows, 0, var);
    tess_j_star = Rcpp::cbind(tess_j_star, new_col);
    
  } else if (p < 0.4 && d_j_length > 1) {
    modification = "RD";
    int removed_dim_idx = floor(R::runif(0, 1) * d_j_length);
    dim_j_star.erase(removed_dim_idx);
    
    NumericMatrix temp_tess(tess_j_rows, d_j_length - 1);
    int current_col = 0;
    for (int i = 0; i < d_j_length; ++i) {
      if (i != removed_dim_idx) {
        temp_tess(_, current_col++) = tess_j_star(_, i);
      }
    }
    tess_j_star = temp_tess;
    
  } else if (p < 0.6 || (p < 0.8 && tess_j_rows == 1)) {
    modification = "AC";
    NumericVector new_row = Rcpp::rnorm(d_j_length, 0, var);
    NumericMatrix temp_tess(tess_j_rows + 1, d_j_length);
    for(int i = 0; i < tess_j_rows; ++i) temp_tess(i, _) = tess_j_star(i, _);
    temp_tess(tess_j_rows, _) = new_row;
    tess_j_star = temp_tess;
    
  } else if (p < 0.8 && tess_j_rows > 1) {
    modification = "RC";
    int removed_row_idx = floor(R::runif(0, 1) * tess_j_rows);
    NumericMatrix temp_tess(tess_j_rows - 1, d_j_length);
    int current_row = 0;
    for (int i = 0; i < tess_j_rows; ++i) {
      if (i != removed_row_idx) {
        temp_tess(current_row++, _) = tess_j_star(i, _);
      }
    }
    tess_j_star = temp_tess;
    
  } else if (p < 0.9 || d_j_length == numCovariates) {
    int centre_to_change_idx = floor(R::runif(0, 1) * tess_j_rows);
    tess_j_star(centre_to_change_idx, _) = Rcpp::rnorm(d_j_length, 0, var);
    
  } else {
    modification = "Swap";
    int dim_to_change_idx = floor(R::runif(0, 1) * d_j_length);
    int new_dim;
    do {
      new_dim = floor(R::runif(0, 1) * numCovariates) + 1;
    } while (is_true(any(dim_j == new_dim)));
    
    dim_j_star[dim_to_change_idx] = new_dim;
    tess_j_star(_, dim_to_change_idx) = Rcpp::rnorm(tess_j_rows, 0, var);
  }
  
  return List::create(
    _["tess_j_star"] = tess_j_star,
    _["dim_j_star"] = dim_j_star,
    _["Modification"] = modification
  );
}