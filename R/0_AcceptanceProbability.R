#' @title Acceptance Probability
#'
#' @description Calculates the acceptance rate of the proposed tessellation.
#'
#' @param x The covariate matrix.
#' @param Tess The tessellation.
#' @param Dim The dimensions of the tessellation.
#' @param j The index of the tessellation.
#' @param R_ijOld The residuals of the old tessellation.
#' @param n_ijOld The number of observations in each cell of the old
#' tessellation.
#' @param R_ijNew The residuals of the new tessellation.
#' @param n_ijNew The number of observations in each cell of the new
#' tessellation.
#' @param SigmaSquared The variance of the residuals.
#' @param Modification The type of modification.
#' @param SigmaSquaredMu The variance of the random effects.
#' @param Omega The prior probability of adding a dimension.
#' @param LambdaRate The rate of the Poisson distribution for the number of
#' centres.
#' @param NumCovariates The total number of covariates.
#'
#' @return The acceptance probability.
#'
#' @keywords internal
#' @noRd
#'
acceptanceProbability <- function(R_ijOld, n_ijOld,
                                  R_ijNew, n_ijNew,
                                  tess_j_star, dim_j_star, # pass single objects
                                  SigmaSquared,
                                  Modification, SigmaSquaredMu,
                                  Omega, LambdaRate, NumCovariates) {
  
  d <- length(dim_j_star)
  cStar <- nrow(tess_j_star)
  
  LogLikelihoodRatio <- 0.5 * (log(prod(n_ijOld * SigmaSquaredMu + SigmaSquared)) -
                                 log(prod(n_ijNew * SigmaSquaredMu + SigmaSquared))) +
    ((SigmaSquaredMu / (2 * SigmaSquared)) *
       (sum((R_ijNew^2) / (n_ijNew * SigmaSquaredMu + SigmaSquared)) -
          sum((R_ijOld^2) / (n_ijOld * SigmaSquaredMu + SigmaSquared))))
  
  # This part of the logic remains unchanged from your original
  if (Modification == "AD") {
    TessStructure <- (dbinom(d - 1, NumCovariates - 1, Omega / NumCovariates)) /
      (dbinom(d - 2, NumCovariates - 1, Omega / NumCovariates) *
         (NumCovariates - d + 1))
    TransitionRatio <- (NumCovariates - d + 1) / d
    AcceptanceProb <- LogLikelihoodRatio + log(TessStructure) +
      log(TransitionRatio)
    if (length(dim_j_star) == 1) {
      AcceptanceProb <- AcceptanceProb + log(1 / 2)
    } else if (length(dim_j_star) == NumCovariates - 1) {
      AcceptanceProb <- AcceptanceProb + log(2)
    }
  } else if (Modification == "RD") {
    TessStructure <- (dbinom(d - 1, NumCovariates, Omega / NumCovariates) *
                        (NumCovariates - d)) /
      (dbinom(d, NumCovariates, Omega / NumCovariates))
    TransitionRatio <- (d + 1) / (NumCovariates - d)
    AcceptanceProb <- LogLikelihoodRatio + log(TessStructure) +
      log(TransitionRatio)
    if (length(dim_j_star) == NumCovariates) {
      AcceptanceProb <- AcceptanceProb + log(1 / 2)
    } else if (length(dim_j_star) == 2) {
      AcceptanceProb <- AcceptanceProb + log(2)
    }
  } else if (Modification == "AC") {
    TessStructure <- dpois(cStar - 1, LambdaRate) / dpois(cStar - 2, LambdaRate)
    TransitionRatio <- 1 / cStar
    AcceptanceProb <- LogLikelihoodRatio + log(TessStructure) +
      log(TransitionRatio) + 0.5 * log(SigmaSquared)
    if (cStar == 1) {
      AcceptanceProb <- AcceptanceProb + log(1 / 2)
    }
  } else if (Modification == "RC") {
    TessStructure <- dpois(cStar - 1, LambdaRate) / dpois(cStar, LambdaRate)
    TransitionRatio <- cStar + 1
    AcceptanceProb <- LogLikelihoodRatio + log(TessStructure) +
      log(TransitionRatio) - 0.5 * log(SigmaSquared)
    if (cStar == 2) {
      AcceptanceProb <- AcceptanceProb + log(2)
    }
  } else if (Modification == "Change") {
    TessStructure <- 1
    TransitionRatio <- 1
    AcceptanceProb <- LogLikelihoodRatio + log(TessStructure) +
      log(TransitionRatio)
  } else { # Includes "Swap"
    TessStructure <- 1
    TransitionRatio <- 1
    AcceptanceProb <- LogLikelihoodRatio + log(TessStructure) +
      log(TransitionRatio)
  }
  return(AcceptanceProb)
}