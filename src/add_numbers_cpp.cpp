#include <Rcpp.h>

using namespace Rcpp;

//' Add Two Numbers in C++
 //'
 //' This function adds two numbers together.
 //'
 //' @param x The first number.
 //' @param y The second number.
 //'
 //' @return The sum of `x` and `y`.
 //'
 //' @examples
 //' add_numbers_cpp(2, 3)
 // [[Rcpp::export]]
 float add_numbers_cpp(float x, float y) {
   return x + y;
 }
