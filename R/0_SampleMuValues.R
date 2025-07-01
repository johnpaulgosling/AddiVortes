#' @title Sample_mu_values
#'
#' @description This function samples the new mean output values for
#' the proposed tessellation.
#'
#' @param j The index of the tessellation.
#' @param Tess The tessellation.
#' @param R_ijNew The new R_ij values.
#' @param n_ijNew The new n_ij values.
#' @param sigmaSquaredMu The sigma squared mu value.
#' @param SigmaSquared The Sigma squared value.
#'
#' @return The new output values for the new tessellation.
Sample_mu_values <- function(j, Tess,
                             R_ijNew, n_ijNew,
                             sigmaSquaredMu, SigmaSquared) {
  # 1. Get the number of samples needed
  N <- length(Tess[[j]][, 1]) # Or nrow(Tess[[j]])

  # 2. Calculate the vector of means for rnorm
  denominator_vec <- sigmaSquaredMu * n_ijNew + SigmaSquared
  mean_vec <- (sigmaSquaredMu * R_ijNew) / denominator_vec

  # 3. Calculate the vector of standard deviations for rnorm
  #    Calculate variance first for clarity, then sqrt
  #    Ensure SigmaSquared, sigmaSquaredMu are non-negative as expected for variances
  variance_vec <- (SigmaSquared * sigmaSquaredMu) / denominator_vec
  # Handle potential negative variance if inputs aren't guaranteed non-negative,
  # though typically variance terms are >= 0. sqrt() will produce NaN for negative input.
  sd_vec <- sqrt(variance_vec)

  # 4. Call rnorm once with vector arguments
  return(rnorm(n = N,
               mean = mean_vec,
               sd = sd_vec))
}
