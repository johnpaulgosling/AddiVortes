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
#' It can return mean predictions, quantiles and optionally calculate the
#' Root Mean Squared Error (RMSE) if true outcomes are provided.
#'
#' @param object An object of class `AddiVortesFit`, typically the result of a
#'   call to `AddiVortes()`.
#' @param newdata A matrix of covariates for the new test set. The number of
#'   columns must match the original training data.
#' @param type The type of prediction required. The default `"response"` gives the
#'   mean prediction. The alternative `"quantile"` returns the quantiles
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
#' This function relies on the internal helper function `applyScaling_internal`
#' being available in the environment, which is used by the main
#' `AddiVortes` function.
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
    # --- Inlined testSetPrediction logic ---
    current_tess <- posteriorTessSamples[[sIdx]]
    current_dim  <- posteriorDimSamples[[sIdx]]
    current_pred <- posteriorPredSamples[[sIdx]]
    
    # Get predictions for each of the m tessellations in the current posterior sample
    predictionList <- lapply(1:mTessellations, function(j) {
      NewTessIndexes <- cellIndices(xNewScaled, current_tess[[j]], current_dim[[j]])
      current_pred[[j]][NewTessIndexes]
    })
    
    # Sum the predictions from all m tessellations for this posterior sample
    predictionsForSampleS <- rowSums(do.call(cbind, predictionList))
    # --- End of inlined logic ---
    
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

#' @title Plot Method for AddiVortesFit
#'
#' @description
#' Generates diagnostic plots for a fitted `AddiVortesFit` object.
#' This function creates a scatter plot of the true versus predicted values
#' for the training set to help visualise model fit.
#'
#' @param aModel An object of class `AddiVortesFit`, typically the result of a
#'   call to `AddiVortes()`.
#' @param x_train A matrix of the original training covariates.
#' @param y_train A numeric vector of the original training true outcomes.
#' @param ... Additional graphical parameters to be passed to the `plot` function
#'   (e.g., `pch`, `cex`).
#'
#' @return
#' This function is called for its side effect of creating a plot and returns
#' `NULL` invisibly.
#'
#' @details
#' The function internally calls `predict.AddiVortesFit` on the provided
#' training data to get the mean predictions. It then plots these predictions
#' against the true outcome values (`y_train`). A dashed red line with an
#' intercept of 0 and a slope of 1 is added to the plot to represent a perfect
#' prediction, making it easy to assess the model's accuracy.
#'
#' @importFrom graphics plot abline title
#' @export
#' @method plot AddiVortesFit
#'
#' @examples
#' \dontrun{
#' # Assuming 'fit' is a trained object of class AddiVortesFit,
#' # and 'x_train_data', 'y_train_data' are your training datasets.
#'
#' plot(fit, x_train = x_train_data, y_train = y_train_data)
#' }
plot.AddiVortesFit <- function(aModel, x_train, y_train, ...) {
  # --- Input Validation ---
  if (!inherits(aModel, "AddiVortesFit")) {
    stop("`x` must be an object of class 'AddiVortesFit'.")
  }
  if (missing(x_train) || missing(y_train)) {
    stop("`x_train` and `y_train` must be provided to plot true vs. predicted values.")
  }
  if (!is.matrix(x_train)) {
    stop("`x_train` must be a matrix.")
  }
  if (!is.numeric(y_train)) {
    stop("`y_train` must be a numeric vector.")
  }
  if (nrow(x_train) != length(y_train)) {
    stop("The number of rows in `x_train` must match the length of `y_train`.")
  }
  
  # Generate mean predictions for the training set
  preds <- predict(aModel, 
                   newdata = x_train,
                   type = "response")
  
  # --- Create the Plot ---
  # Plot true values vs. predicted values
  plot(y_train,
       preds,
       xlab = "True Values",
       ylab = "Predicted Values",
       main = "AddiVortes Predictions vs True Values",
       xlim = range(c(y_train, preds)),
       ylim = range(c(y_train, preds)),
       pch = 19, col = "darkblue"
  )
  
  # Add the line of equality (y = x) for reference
  abline(a = 0, b = 1, col = "darkred", lwd = 2)
  
  # Get quantile predictions to create error bars/intervals
  preds_quantile <- predict(aModel,
                            x_train,
                            "quantile")
  
  # Add error segments for each prediction
  for (i in 1:nrow(preds_quantile)) {
    segments(y_train, preds_quantile[i, 1],
             y_train, preds_quantile[i, 2],
             col = "darkblue", lwd = 1
    )
  }
  
  # Add legend
  legend("bottomright",
         legend=c("Prediction",
                  "Equality Line"), 
         col=c("darkblue",
               "darkred"),
         lty=1, pch=c(19, NA), lwd=2)
  
  # The function is used for plotting alone, so return NULL invisibly
  invisible(NULL)
}