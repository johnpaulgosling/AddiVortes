#' @title Fitting function
#'
#' @description A function that calculates the squared difference between sigma squared hat
#' and the inverse gamma function.
#'
#' @param lambda The lambda value.
#' @param q The q value.
#' @param nu The nu value.
#' @param sigmaSquaredHat The sigma squared hat value.
#'
#' @return The squared difference between sigma squared hat and the inverse gamma function.
fittingFunction <- function(lambda,
                            q,
                            nu,
                            sigmaSquaredHat) {
  return((sigmaSquaredHat - qinvgamma_internal(q,
                                               shape = nu / 2,
                                               rate = nu * lambda / 2))^2)
}
