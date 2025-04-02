#' @title Sample_mu_values
#'
#' @description This function samples the new output values for the new tessellation.
#'
#' @param j The index of the tessellation.
#' @param Tess The tessellation.
#' @param R_ijNew The new R_ij values.
#' @param n_ijNew The new n_ij values.
#' @param sigmaSquaredMu The sigma squared mu value.
#' @param SigmaSquared The Sigma squared value.
#'
#' @return The new output values for the new tessellation.
#'
#' @export
Sample_mu_values <- function(j, Tess, R_ijNew, n_ijNew, sigmaSquaredMu, SigmaSquared) {
  PredSet <- rep(0, length(Tess[[j]][, 1]))
  for (i in 1:length(Tess[[j]][, 1])) {
    PredSet[i] <- rnorm(
      1, sigmaSquaredMu * R_ijNew[i] /
        (sigmaSquaredMu * n_ijNew[i] + SigmaSquared),
      ((SigmaSquared * sigmaSquaredMu) /
        (n_ijNew[i] * sigmaSquaredMu + SigmaSquared))^0.5
    )
  }
  return(PredSet)
}
