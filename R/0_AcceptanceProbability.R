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
#' @param dim_j_star A vector of integers representing the dimensions included
#'   in the proposed new tessellation.
#' @param SigmaSquared The variance of the model's residuals.
#' @param Modification A character string indicating the type of modification
#'   proposed. Accepted values are "AD" (Add Dimension), "RD" (Remove
#'   Dimension), "AC" (Add Centre), "RC" (Remove Centre), "Change" and "Swap".
#' @param SigmaSquaredMu The variance of the random effects component of the
#'   model.
#' @param Omega Omega/(number of covariates) is the prior probability of adding
#' a dimension.
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

  # Clamp probability to [0, 1) range to avoid NaN in dbinom.
  # When Omega/NumCovariates >= 1, clamping to just below 1 represents a very
  # high (but finite) probability of including dimensions. This keeps the
  # acceptance probability numerically well-defined: dbinom(k, n, 1) = 0 for
  # k < n, which causes 0/0 = NaN in the AD/RD ratio formulas.
  prob_eps <- 1e-10 # keep prob strictly below 1 to avoid 0/0 in dbinom
  prob <- min(1 - prob_eps, max(0, Omega / NumCovariates))

  addDimensionProbability <- function(dCurrent) {
    if (dCurrent >= NumCovariates) {
      return(0)
    }
    if (dCurrent == 1) {
      return(0.4)
    }
    0.2
  }

  removeDimensionProbability <- function(dCurrent) {
    if (dCurrent <= 1) {
      return(0)
    }
    if (dCurrent == NumCovariates) {
      return(0.4)
    }
    0.2
  }

  addCentreProbability <- function(cCurrent) {
    if (cCurrent == 1) {
      return(0.4)
    }
    0.2
  }

  removeCentreProbability <- function(cCurrent) {
    if (cCurrent > 1) {
      return(0.2)
    }
    0
  }

  LogLikelihoodRatio <- 0.5 * (log(prod(n_ijOld * SigmaSquaredMu +
    SigmaSquared)) -
    log(prod(n_ijNew * SigmaSquaredMu +
      SigmaSquared))) +
    ((SigmaSquaredMu / (2 * SigmaSquared)) *
      (sum((R_ijNew^2) / (n_ijNew * SigmaSquaredMu + SigmaSquared)) -
        sum((R_ijOld^2) / (n_ijOld * SigmaSquaredMu + SigmaSquared))))

  if (Modification == "AD") {
    dOld <- d - 1
    AcceptanceProb <- LogLikelihoodRatio +
      2 * log(NumCovariates - dOld) -
      log(dOld) -
      log(d) +
      log(prob) -
      log1p(-prob) +
      log(removeDimensionProbability(d)) -
      log(addDimensionProbability(dOld))
  } else if (Modification == "RD") {
    dOld <- d + 1
    AcceptanceProb <- LogLikelihoodRatio +
      log(dOld) +
      log(d) -
      2 * log(NumCovariates - d) +
      log1p(-prob) -
      log(prob) +
      log(addDimensionProbability(d)) -
      log(removeDimensionProbability(dOld))
  } else if (Modification == "AC") {
    cOld <- length(n_ijOld)
    AcceptanceProb <- LogLikelihoodRatio +
      log(LambdaRate) -
      2 * log(cStar) +
      log(removeCentreProbability(cStar)) -
      log(addCentreProbability(cOld)) +
      0.5 * log(SigmaSquared)
  } else if (Modification == "RC") {
    cOld <- length(n_ijOld)
    AcceptanceProb <- LogLikelihoodRatio +
      2 * log(cOld) -
      log(LambdaRate) +
      log(addCentreProbability(cStar)) -
      log(removeCentreProbability(cOld)) -
      0.5 * log(SigmaSquared)
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
