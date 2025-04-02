#' @title Acceptance Probability
#'
#' @description Calculates the acceptance rate of the proposed tessellation.
#'
#' @param x The covariate matrix.
#' @param Tess The tessellation.
#' @param Dim The dimensions of the tessellation.
#' @param j The index of the tessellation.
#' @param R_ijOld The residuals of the old tessellation.
#' @param n_ijOld The number of observations in each cell of the old tessellation.
#' @param R_ijNew The residuals of the new tessellation.
#' @param n_ijNew The number of observations in each cell of the new tessellation.
#' @param SigmaSquared The variance of the residuals.
#' @param Modification The type of modification.
#' @param SigmaSquaredMu The variance of the random effects.
#' @param Omega The prior probability of adding a dimension.
#' @param lambda_rate The rate of the Poisson distribution for the number of centers.
#'
#' @return The acceptance probability.
#'
#' @export
Acceptance_Probability <- function(x, Tess, Dim, j, R_ijOld, n_ijOld, R_ijNew, n_ijNew, SigmaSquared,
                                   Modification, SigmaSquaredMu, Omega, lambda_rate) {
  d <- length(Dim[[j]])
  NumCovariates <- length(x[1, ])
  cStar <- length(Tess[[j]][, 1])

  # The Log Likelihood Ratio in the acceptance ratio
  LOGlikelihoodRatio <- 0.5 * (log(prod(n_ijOld * SigmaSquaredMu + SigmaSquared)) - log(prod(n_ijNew * SigmaSquaredMu + SigmaSquared))) + ((SigmaSquaredMu / (2 * SigmaSquared)) * (sum((R_ijNew^2) / (n_ijNew * SigmaSquaredMu + SigmaSquared)) - sum((R_ijOld^2) / (n_ijOld * SigmaSquaredMu + SigmaSquared))))

  # Calculating the acceptance probability for
  # "AD"=Adding a dimension, "RD"=Removing a dimension,
  # "AC"=Adding a center, "RC"=Removing a center,
  # "Change"=Changing the coordinates of a center and Swopping a dimension.
  if (Modification == "AD") {
    TessStructure <- (dbinom(d - 1, NumCovariates - 1, Omega / NumCovariates)) / (dbinom(d - 2, NumCovariates - 1, Omega / NumCovariates) * (NumCovariates - d + 1))
    TransitionRatio <- (NumCovariates - d + 1) / d
    AcceptanceProb <- LOGlikelihoodRatio + log(TessStructure) + log(TransitionRatio)

    # Adjustments.
    if (length(Dim[[j]]) == 1) {
      AcceptanceProb <- AcceptanceProb + log(1 / 2)
    } else if (length(Dim[[j]]) == NumCovariates - 1) {
      AcceptanceProb <- AcceptanceProb + log(2)
    }
  } else if (Modification == "RD") {
    TessStructure <- (dbinom(d - 1, NumCovariates, Omega / NumCovariates) * (NumCovariates - d)) / (dbinom(d, NumCovariates, Omega / NumCovariates))
    TransitionRatio <- (d + 1) / (NumCovariates - d)
    AcceptanceProb <- LOGlikelihoodRatio + log(TessStructure) + log(TransitionRatio)

    # Adjustments.
    if (length(Dim[[j]]) == NumCovariates) {
      AcceptanceProb <- AcceptanceProb + log(1 / 2)
    } else if (length(Dim[[j]]) == 2) {
      AcceptanceProb <- AcceptanceProb + log(2)
    }
  } else if (Modification == "AC") {
    TessStructure <- dpois(cStar - 1, lambda_rate) / dpois(cStar - 2, lambda_rate)
    TransitionRatio <- 1 / cStar
    AcceptanceProb <- LOGlikelihoodRatio + log(TessStructure) + log(TransitionRatio) + 0.5 * log(SigmaSquared)

    # Adjustments.
    if (cStar == 1) {
      AcceptanceProb <- AcceptanceProb + log(1 / 2)
    }
  } else if (Modification == "RC") {
    TessStructure <- dpois(cStar - 1, lambda_rate) / dpois(cStar, lambda_rate)
    TransitionRatio <- cStar + 1
    AcceptanceProb <- LOGlikelihoodRatio + log(TessStructure) + log(TransitionRatio) - 0.5 * log(SigmaSquared)

    # Adjustments.
    if (cStar == 2) {
      AcceptanceProb <- AcceptanceProb + log(2)
    }
  } else if (Modification == "Change") {
    TessStructure <- 1
    TransitionRatio <- 1
    AcceptanceProb <- LOGlikelihoodRatio + log(TessStructure) + log(TransitionRatio)
  } else {
    TessStructure <- 1
    TransitionRatio <- 1
    AcceptanceProb <- LOGlikelihoodRatio + log(TessStructure) + log(TransitionRatio)
  }
  return(AcceptanceProb)
}
