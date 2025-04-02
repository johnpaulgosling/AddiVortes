#' @title Sample Sigma Squared
#'
#' @description This function samples sigma squared from the inverse gamma distribution.
#'
#' @param yScaled A vector of scaled data.
#' @param nu The prior degrees of freedom.
#' @param lambda The prior scale parameter.
#' @param SumOfAllTess The sum of all the tessellations.
#'
#' @return SigmaSquared The sampled sigma squared.
#'
#' @export
Sample_Sigma_Squared <- function(yScaled, nu, lambda, SumOfAllTess) {
  n <- length(yScaled)
  SigmaSquared <- rinvgamma(1, shape = (nu + n) / 2, rate = (nu * lambda + sum((yScaled - SumOfAllTess)^2)) / 2)

  return(SigmaSquared)
}
