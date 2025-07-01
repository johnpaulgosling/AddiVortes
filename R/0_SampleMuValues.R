#' @title Sample ,u Values
#'
#' @description This function samples the new mean output values for
#' the proposed tessellation.
#'
#' @param j The index of the tessellation.
#' @param Tess The tessellation.
#' @param R_ijNew The new R_ij values.
#' @param n_ijNew The new n_ij values.
#' @param sigmaSquaredMu The sigma squared mu value.
#' @param sigmaSquared The Sigma squared value.
#'
#' @return The new output values for the new tessellation.
#'
#' @keywords internal
#' @noRd
#'
sampleMuValues <- function(j, Tess,
                           R_ijNew, n_ijNew,
                           sigmaSquaredMu, sigmaSquared) {
  # 1. Get the number of samples needed
  N <- length(Tess[[j]][, 1]) # Or nrow(Tess[[j]])

  # 2. Calculate the vector of means for rnorm
  denominatorVec <- sigmaSquaredMu * n_ijNew + sigmaSquared
  meanVec <- (sigmaSquaredMu * R_ijNew) / denominatorVec

  # 3. Calculate the vector of standard deviations for rnorm
  #    Calculate variance first for clarity, then sqrt
  #    Ensure SigmaSquared, sigmaSquaredMu are non-negative as expected for variances
  varianceVec <- (sigmaSquared * sigmaSquaredMu) / denominatorVec
  # Handle potential negative variance if inputs aren't guaranteed non-negative,
  # though typically variance terms are >= 0. sqrt() will produce NaN for negative input.
  sdVec <- sqrt(varianceVec)

  # 4. Call rnorm once with vector arguments
  return(rnorm(n = N,
               mean = meanVec,
               sd = sdVec))
}
