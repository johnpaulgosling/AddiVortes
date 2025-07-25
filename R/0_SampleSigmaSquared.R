#' @title Sample Sigma Squared
#'
#' @description This function samples sigma squared from the inverse gamma distribution.
#'
#' @param yScaled A vector of scaled data.
#' @param nu The prior degrees of freedom.
#' @param lambda The prior scale parameter.
#' @param SumOfAllTess The sum of all the tessellations.
#'
#' @return sigmaSquared The sampled sigma squared.
#'
#' @keywords internal
#' @noRd
#'
sampleSigmaSquared <- function(yScaled, nu,
                               lambda, SumOfAllTess) {
  n <- length(yScaled)
  sigmaSquared <- rinvgamma_internal(1,
    shape = (nu + n) / 2,
    rate = (nu * lambda +
      sum((yScaled - SumOfAllTess)^2)) / 2
  )

  return(sigmaSquared)
}
