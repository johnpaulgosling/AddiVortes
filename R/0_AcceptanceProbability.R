#' @title Calculate Acceptance Probability
#' @description This function calculates the log of the acceptance probability 
#'   for a proposed modification to a tessellation in a Bayesian model. It is 
#'   used within a Metropolis-Hastings algorithm to decide whether to accept a 
#'   new tessellation.
#'
#' @param R_ijOld A numeric vector of residuals for each cell in the current
#'   tessellation.
#' @param n_ijOld A numeric vector containing the number of observations in each
#'   cell of the current tessellation.
#' @param R_ijNew A numeric vector of residuals for each cell in the proposed
#'   new tessellation.
#' @param n_ijNew A numeric vector containing the number of observations in each
#'   cell of the proposed new tessellation.
#' @param tess_j_star The proposed new tessellation structure, which is a matrix
#'   of centres.
#' @param dim_j_star A vector of integers representing the dimensions included in
#'   the proposed new tessellation.
#' @param SigmaSquared The variance of the model's residuals.
#' @param Modification A character string indicating the type of modification
#'   proposed. Accepted values are "AD" (Add Dimension), "RD" (Remove
#'   Dimension), "AC" (Add Centre), "RC" (Remove Centre), "Change" and "Swap".
#' @param SigmaSquaredMu The variance of the random effects component of the
#'   model.
#' @param Omega The prior probability of including a covariate in a
#'   tessellation's dimension.
#' @param LambdaRate The rate parameter for the Poisson prior on the number of
#'   centres in the tessellation.
#' @param NumCovariates The total number of covariates available in the dataset.
#'
#' @return The natural logarithm of the acceptance probability for the proposed
#'   modification.
#'
#' @keywords internal
#' @noRd
acceptanceProbability <- function(R_ijOld, n_ijOld,
                                  R_ijNew, n_ijNew,
                                  tess_j_star, dim_j_star,
                                  SigmaSquared,
                                  Modification, SigmaSquaredMu,
                                  Omega, LambdaRate,
                                  NumCovariates) {
  
  d <- length(dim_j_star)
  cStar <- nrow(tess_j_star)
  
  # Clamp probability to [0, 1] range to avoid NaN in dbinom.
  # When Omega/NumCovariates > 1 (e.g., Omega=3, NumCovariates=2), 
  # clamping to 1.0 represents 100% probability of including dimensions.
  # This allows the model to work with small numbers of covariates.
  prob <- min(1, max(0, Omega / NumCovariates))
  
  LogLikelihoodRatio <- 0.5 * (log(prod(n_ijOld * SigmaSquaredMu +
                                          SigmaSquared)) -
                                 log(prod(n_ijNew * SigmaSquaredMu +
                                            SigmaSquared))) +
    ((SigmaSquaredMu / (2 * SigmaSquared)) *
       (sum((R_ijNew^2) / (n_ijNew * SigmaSquaredMu + SigmaSquared)) -
          sum((R_ijOld^2) / (n_ijOld * SigmaSquaredMu + SigmaSquared))))
  
  if (Modification == "AD") {
    TessStructure <- (dbinom(d - 1, NumCovariates - 1, prob)) /
      (dbinom(d - 2, NumCovariates - 1, prob) *
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
    TessStructure <- (dbinom(d - 1, NumCovariates, prob) *
                        (NumCovariates - d)) /
      (dbinom(d, NumCovariates, prob))
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