#' @title Create an AddiVortesFit Object
#'
#' @description A constructor for the AddiVortesFit class.
#'
#' @param posteriorTess A list of the posterior samples of the tessellations.
#' @param posteriorDim A list of the posterior samples of the dimensions.
#' @param posteriorPred A list of the posterior samples of the predictions.
#' @param xCentres The centres of the covariates.
#' @param xRanges The ranges of the covariates.
#' @param yCentre The centre of the output values.
#' @param yRange The range of the output values.
#' @param inSampleRmse The in-sample RMSE.
#'
#' @return An object of class AddiVortesFit.
#' @export
new_AddiVortesFit <- function(posteriorTess, posteriorDim, posteriorPred,
                              xCentres, xRanges, yCentre, yRange,
                              inSampleRmse) {
  structure(
    list(
      posteriorTess = posteriorTess,
      posteriorDim = posteriorDim,
      posteriorPred = posteriorPred,
      xCentres = xCentres,
      xRanges = xRanges,
      yCentre = yCentre,
      yRange = yRange,
      inSampleRmse = inSampleRmse
    ),
    class = "AddiVortesFit"
  )
}


#' @title Predict Method for AddiVortesFit
#'
#' @description
#' Predicts outcomes for new data using a fitted `AddiVortesFit` model object.
#' It can return mean predictions, quantiles, and optionally calculate the
#' Root Mean Squared Error (RMSE) if true outcomes are provided.
#'
#' @param object An object of class `AddiVortesFit`, typically the result of a
#'   call to `AddiVortes()`.
#' @param newdata A matrix of covariates for the new test set. The number of
#'   columns must match the original training data.
#' @param type The type of prediction required. The default, `"response"`, gives the
#'   mean prediction. The alternative, `"quantile"`, returns the quantiles
#'   specified by the `quantiles` argument.
#' @param quantiles A numeric vector of probabilities with values in [0, 1] to
#'   compute for the predictions when `type = "quantile"`.
#' @param ... Further arguments passed to or from other methods (currently unused).
#'
#' @return
#' If `type = "response"`, a numeric vector of mean predictions.
#' If `type = "quantile"`, a matrix where each row corresponds to an observation
#' in `newdata` and each column to a quantile.
#' If `yNew` is provided, the returned object will have an `rmse` attribute
#' containing the calculated Root Mean Squared Error.
#'
#' @details
#' This function relies on the internal helper functions `applyScaling_internal`
#' and `testSetPrediction` being available in the environment, which are used
#' by the main `AddiVortes` function.
#'
#' @export
#' @method predict AddiVortesFit
predict.AddiVortesFit <- function(object, newdata,
                                  type = c("response", "quantile"),
                                  quantiles = c(0.025, 0.975), ...) {
  # Match the type argument
  type <- match.arg(type)

  # --- Input Validation ---
  if (!inherits(object, "AddiVortesFit")) {
    stop("`object` must be of class 'AddiVortesFit'.")
  }
  if (!is.matrix(newdata)) {
    stop("`newdata` must be a matrix.")
  }
  if (ncol(newdata) != length(object$xCentres)) {
    stop("Number of columns in `newdata` does not match the original training data.")
  }

  # Extract model components from the fitted object
  posteriorTessSamples <- object$posteriorTess
  posteriorDimSamples <- object$posteriorDim
  posteriorPredSamples <- object$posteriorPred
  numStoredSamples <- length(posteriorTessSamples)

  # Handle cases with no stored samples
  if (numStoredSamples == 0) {
    warning("The AddiVortes model contains no posterior samples. Cannot make predictions.")
    return(NA_real_)
  }

  # Scale the new test covariates using stored scaling parameters
  xNewScaled <- applyScaling_internal(
    mat = newdata,
    centres = object$xCentres,
    ranges = object$xRanges
  )

  # Determine the number of tessellations (m) from the first stored sample
  mTessellations <- length(posteriorTessSamples[[1]])

  # Initialize a matrix to store predictions for each posterior sample
  newTestDataPredictionsMatrix <- array(dim = c(nrow(xNewScaled), numStoredSamples))

  # --- Prediction Loop ---
  for (sIdx in 1:numStoredSamples) {
    # Get predictions for the current posterior sample
    predictionsForSampleS <- testSetPrediction(
      xNewScaled,
      mTessellations,
      posteriorTessSamples[[sIdx]],
      posteriorDimSamples[[sIdx]],
      posteriorPredSamples[[sIdx]]
    )
    newTestDataPredictionsMatrix[, sIdx] <- predictionsForSampleS
  }

  # --- Process and Unscale Predictions ---
  if (type == "response") {
    # Calculate the mean of predictions across all posterior samples (still scaled)
    meanYhatNewScaled <- rowMeans(newTestDataPredictionsMatrix)
    # Unscale the mean predictions
    predictions <- meanYhatNewScaled * object$yRange + object$yCentre
  } else if (type == "quantile") {
    # Calculate the specified quantiles of the predictions
    quantileYhatNewScaled <- apply(newTestDataPredictionsMatrix, 1, quantile,
      probs = quantiles, na.rm = TRUE
    )
    # Unscale the quantiles of predictions (transpose for correct dimensions)
    predictions <- t(quantileYhatNewScaled * object$yRange + object$yCentre)
  }

  return(predictions)
}
