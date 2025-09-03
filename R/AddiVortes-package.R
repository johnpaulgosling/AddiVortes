#' AddiVortes: Bayesian Additive Voronoi Tessellations for Machine Learning
#' 
#' @description 
#' AddiVortes implements Bayesian Additive Voronoi Tessellation models for 
#' machine learning regression and non-parametric statistical modeling. This 
#' package provides a flexible alternative to BART (Bayesian Additive Regression 
#' Trees), using Voronoi tessellations instead of trees for spatial partitioning.
#' The method is particularly effective for spatial data analysis, complex 
#' function approximation, and Bayesian regression modeling.
#' 
#' @details
#' Key features include:
#' \itemize{
#'   \item Machine learning regression with Bayesian inference
#'   \item Alternative to BART using Voronoi tessellations
#'   \item Spatial data analysis and modeling
#'   \item Non-parametric regression capabilities
#'   \item Complex function approximation
#'   \item Uncertainty quantification through posterior inference
#' }
#' 
#' @keywords package machine-learning bayesian regression BART spatial tessellation
#' @references 
#' Stone, A. and Gosling, J.P. (2025). AddiVortes: (Bayesian) additive Voronoi 
#' tessellations. Journal of Computational and Graphical Statistics.
#' 
#' @seealso 
#' \url{https://johnpaulgosling.github.io/AddiVortes/}
"_PACKAGE"

## usethis namespace: start
#' @useDynLib AddiVortes, .registration = TRUE
## usethis namespace: end
NULL