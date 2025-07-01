#' @title Fitting function
#'
#' @description A function that calculates the squared difference between sigma squared hat
#' and the inverse gamma function.
#'
#' @param lambda The lambda value.
#' @param q The q value.
#' @param nu The nu value.
#' @param sigmaSquared_hat The sigma squared hat value.
#'
#' @return The squared difference between sigma squared hat and the inverse gamma function.
Fitting_Function <- function(lambda,
                             q,
                             nu,
                             sigmaSquared_hat) {
  return((sigmaSquared_hat - qinvgamma_internal(q,
                                                shape = nu / 2,
                                                rate = nu * lambda / 2))^2)
}
