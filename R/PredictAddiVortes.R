#' @title predictAddiVortes
#'
#' @description
#' Utilises a fitted addiVortes model to make predictions on new test data and
#' calculates the Root Mean Squared Error (RMSE). This function aims to
#' reproduce the Out-of-Sample RMSE calculation performed internally by the
#' addiVortes fitting function.
#'
#' @param addiVortesModelFit The 'addiVortesModel' list object returned by the
#' addiVortes function. This contains posterior samples
#' and scaling parameters.
#' @param xNew A matrix of the covariates for the new test set. The number of
#' columns must match the original training data.
#' @param yNew A vector of the true output values for the new test set. Must have
#' the same number of observations as rows in `xNew`.
#' @param quantiles A numeric vector of quantiles to compute for the predictions.
#'
#' @return A list that contains:
#' - `outOfSampleRmse`: The Root Mean Squared Error of the predictions on the
#' new test data.
#' - `meanYhatNewUnscaled`: The mean predictions for the new test data,
#' unscaled to the original scale of the output variable.
#' - `quantileYhatNewUnscaled`: A matrix of quantiles of the predictions for the
#' new test data, unscaled to the original output variable scale.
#'
#' @details
#' This function requires the helper functions `applyScalingInternal` and
#' `testSetPrediction` (used by the original `addiVortes` function) to be
#' available in the R environment.
#'
#' @export
predictAddiVortes <- function(addiVortesModelFit,
                              xNew,
                              yNew,
                              quantiles = c(0.025, 0.975)) {

  # Extract components from the fitted model object
  posteriorTessSamples <- addiVortesModelFit[[1]]
  posteriorDimSamples  <- addiVortesModelFit[[2]]
  posteriorPredSamples <- addiVortesModelFit[[3]]
  # currentStorageIdx indicates (number of stored samples + 1)
  numStoredSamples     <- addiVortesModelFit[[4]] - 1

  xModelCenters <- addiVortesModelFit[[5]]
  xModelRanges  <- addiVortesModelFit[[6]]
  yModelCenter  <- addiVortesModelFit[[7]]
  yModelRange   <- addiVortesModelFit[[8]]

  # Basic input validation
  if (!is.matrix(xNew)) {
    stop("'xNew' must be a matrix.")
  }
  if (ncol(xNew) != length(xModelCenters)) {
    stop("Number of columns in 'xNew' does not match the original training data.")
  }
  if (length(yNew) != nrow(xNew)) {
    stop("Length of 'yNew' must match the number of rows in 'xNew'.")
  }

  # Handle cases with no stored samples
  if (numStoredSamples <= 0) {
    warning("The addiVortes model contains no stored posterior samples. Cannot make predictions.")
    return(NA_real_)
  }

  # Scale the new test covariates using stored scaling parameters
  xNewScaled <- applyScaling_internal(mat = xNew,
                                     centers = xModelCenters,
                                     ranges = xModelRanges)

  # Determine the number of tessellations (m) from the first stored sample
  # Assumes m is constant across all samples
  mTessellations <- length(posteriorTessSamples[[1]])

  # Initialize a matrix to store predictions for each posterior sample
  # Dimensions: [number of new test observations, number of stored posterior samples]
  newTestDataPredictionsMatrix <- array(dim = c(nrow(xNewScaled),
                                                numStoredSamples))

  # Loop through each stored posterior sample to make predictions
  for (sIdx in 1:numStoredSamples) {
    currentTess <- posteriorTessSamples[[sIdx]]
    currentDim  <- posteriorDimSamples[[sIdx]]
    currentPred <- posteriorPredSamples[[sIdx]]

    # Get predictions for the current posterior sample
    # Assumes testSetPrediction(xTestScaled, m, tess, dim, pred) is available
    predictionsForSampleS <- testSetPrediction(xNewScaled,
                                               mTessellations,
                                               currentTess,
                                               currentDim,
                                               currentPred)
    newTestDataPredictionsMatrix[, sIdx] <- predictionsForSampleS
  }

  # Calculate the mean of predictions across all posterior samples (still scaled)
  meanYhatNewScaled <- rowSums(newTestDataPredictionsMatrix) / numStoredSamples

  # Calculate the specified quantiles of the predictions
  quantileYhatNewScaled <- apply(newTestDataPredictionsMatrix, 1, quantile,
                                 probs = quantiles)

  # Unscale the mean predictions
  meanYhatNewUnscaled <- meanYhatNewScaled * yModelRange + yModelCenter

  # Unscale the quantiles of predictions
  quantileYhatNewUnscaled <- quantileYhatNewScaled * yModelRange + yModelCenter

  # Calculate the Out-of-Sample RMSE for the new test data
  outOfSampleRmse <- sqrt(mean((yNew - meanYhatNewUnscaled)^2))

  return(list(outOfSampleRmse,
              meanYhatNewUnscaled,
              quantileYhatNewUnscaled))
}
