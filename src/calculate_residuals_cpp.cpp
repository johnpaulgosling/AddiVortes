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