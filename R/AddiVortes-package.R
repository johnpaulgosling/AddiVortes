#' AddiVortes: Bayesian Additive Voronoi Tessellations for Machine Learning
#'
#' @description
#' AddiVortes implements Bayesian Additive Voronoi Tessellation models for
#' machine learning regression, classification and non-parametric statistical
#' modelling. This package provides a flexible alternative to BART (Bayesian
#' Additive Regression Trees), using Voronoi tessellations instead of trees for
#' spatial partitioning. The method is particularly effective for spatial data
#' analysis, complex function approximation, and Bayesian regression and
#' classification.
#'
#' @details
#' Key features include:
#' \itemize{
#'   \item Machine learning regression and classification with Bayesian inference
#'   \item Binary and multinomial probit models with latent-variable data
#'     augmentation
#'   \item Alternative to BART using Voronoi tessellations
#'   \item Spatial data analysis and modelling
#'   \item Non-parametric regression capabilities
#'   \item Complex function approximation
#'   \item Uncertainty quantification through posterior inference
#' }
#'
#' @keywords package machine-learning bayesian regression classification BART spatial tessellation
#' @references
#' Stone, A. and Gosling, J.P. (2025). AddiVortes: (Bayesian) additive Voronoi
#' tessellations. Journal of Computational and Graphical Statistics.
#'
#' Stone, A.J., Ogundimu, E. and Gosling, J.P. (2026). Binary AddiVortes:
#' (Bayesian) Additive Voronoi Tessellations for Binary Classification with an
#' application to Predicting Home Mortgage Application Outcomes.
#'
#' Albert, J.H. and Chib, S. (1993). Bayesian analysis of binary and
#' polychotomous response data. Journal of the American Statistical Association.
#'
#' Kindo, B.P., Wang, H. and Peña, E.A. (2016). Multinomial probit Bayesian
#' additive regression trees. Stat.
#'
#' @seealso
#' \url{https://johnpaulgosling.github.io/AddiVortes/}
"_PACKAGE"

## usethis namespace: start
#' @useDynLib AddiVortes, .registration = TRUE
## usethis namespace: end
NULL
