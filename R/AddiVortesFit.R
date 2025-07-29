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

#' @title Print Method for AddiVortesFit
#'
#' @description
#' Prints a summary of a fitted `AddiVortesFit` object, providing information
#' about the model structure, dimensions, and fit quality similar to the
#' output of a linear model summary.
#'
#' @param x An object of class `AddiVortesFit`, typically the result of a
#'   call to `AddiVortes()`.
#' @param ... Further arguments passed to or from other methods (currently 
#' unused).
#'
#' @return
#' The function is called for its side effect of printing model information
#' and returns the input object `x` invisibly.
#'
#' @details
#' The print method displays:
#' - The model formula representation
#' - Number of covariates and posterior samples
#' - Number of tessellations used
#' - In-sample RMSE
#' - Covariate scaling information
#'
#' @export
#' @method print AddiVortesFit
print.AddiVortesFit <- function(x, ...) {
  # --- Input Validation ---
  if (!inherits(x, "AddiVortesFit")) {
    stop("`x` must be an object of class 'AddiVortesFit'.")
  }
  
  cat("AddiVortes Model\n")
  cat("================\n\n")
  
  # Model equation representation
  num_covariates <- length(x$xCentres)
  covariate_names <- if (is.null(names(x$xCentres))) {
    paste0("X", 1:num_covariates)
  } else {
    names(x$xCentres)
  }
  
  cat("Model Formula:\n")
  if (num_covariates == 1) {
    cat("Y ~ f(", covariate_names[1], ")\n\n")
  } else {
    formula_str <- paste0("Y ~ f(", paste(covariate_names, collapse = ", "), ")")
    cat(formula_str, "\n")
    cat("where f(.) is represented by additive Voronoi tessellations\n\n")
  }
  
  # Model dimensions
  num_samples <- length(x$posteriorTess)
  num_tessellations <- if (num_samples > 0) {
    length(x$posteriorTess[[1]])
  } else {
    0
  }
  
  cat("Model Information:\n")
  cat("Number of covariates:     ", num_covariates, "\n")
  cat("Number of tessellations:  ", num_tessellations, "\n")
  cat("Posterior samples:        ", num_samples, "\n")
  cat("In-sample RMSE:           ", round(x$inSampleRmse, 4), "\n\n")
  
  # Scaling information
  cat("Covariate Scaling:\n")
  scaling_df <- data.frame(
    Covariate = covariate_names,
    Centre = round(x$xCentres, 4),
    Range = round(x$xRanges, 4)
  )
  print(scaling_df, row.names = FALSE)
  
  cat("\nOutput Scaling:\n")
  cat("Centre: ", round(x$yCentre, 4), "\n")
  cat("Range:  ", round(x$yRange, 4), "\n\n")
  
  # Additional model information
  if (num_samples > 0) {
    # Get some statistics about the tessellations
    tess_sizes <- sapply(1:min(5, num_samples), function(i) {
      sapply(x$posteriorTess[[i]], function(tess) nrow(tess))
    })
    
    if (is.matrix(tess_sizes)) {
      avg_tess_size <- round(mean(tess_sizes), 1)
      range_tess_size <- range(tess_sizes)
    } else {
      avg_tess_size <- round(mean(tess_sizes), 1)
      range_tess_size <- range(tess_sizes)
    }
    
    cat("Tessellation Statistics (from first ", min(5, num_samples), " samples):\n")
    cat("Average cells per tessellation: ", avg_tess_size, "\n")
    cat("Range of cells per tessellation: [", range_tess_size[1], ", ", 
        range_tess_size[2], "]\n")
  } else {
    cat("No posterior samples available.\n")
  }
  
  # Return the object invisibly
  invisible(x)
}

#' @title Summary Method for AddiVortesFit
#'
#' @description
#' Provides a detailed summary of a fitted `AddiVortesFit` object, including
#' more comprehensive information than the print method.
#'
#' @param object An object of class `AddiVortesFit`, typically the result of a
#'   call to `AddiVortes()`.
#' @param ... Further arguments passed to or from other methods (currently 
#' unused).
#'
#' @return
#' The function is called for its side effect of printing detailed model 
#' information and returns the input object `object` invisibly.
#'
#' @export
#' @method summary AddiVortesFit
summary.AddiVortesFit <- function(object, ...) {
  # Call the print method first
  print(object)
  
  # Add additional summary information
  if (length(object$posteriorTess) > 0) {
    cat("\nDetailed Posterior Information:\n")
    cat("===============================\n")
    
    # Analyze tessellation complexity across samples
    all_tess_sizes <- sapply(object$posteriorTess, function(sample) {
      sapply(sample, function(tess) nrow(tess))
    })
    
    if (is.matrix(all_tess_sizes)) {
      cat("Tessellation complexity across all samples:\n")
      for (j in 1:nrow(all_tess_sizes)) {
        cat("  Tessellation ", j, ": mean = ", round(mean(all_tess_sizes[j,]), 1),
            ", sd = ", round(sd(all_tess_sizes[j,]), 2), "\n")
      }
    }
    
    # Dimension information if available
    if (length(object$posteriorDim) > 0) {
      dim_info <- sapply(object$posteriorDim, function(sample) {
        sapply(sample, length)
      })
      
      if (is.matrix(dim_info)) {
        cat("\nActive dimensions per tessellation:\n")
        for (j in 1:nrow(dim_info)) {
          cat("  Tessellation ", j, ": mean = ", round(mean(dim_info[j,]), 1),
              " dimensions\n")
        }
      }
    }
  }
  
  invisible(object)
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
#' @param type The type of prediction required. The default `"response"` gives 
#'   the mean prediction. The alternative `"quantile"` returns the quantiles
#'   specified by the `quantiles` argument.
#' @param quantiles A numeric vector of probabilities with values in [0, 1] to
#'   compute for the predictions when `type = "quantile"`.
#' @param ... Further arguments passed to or from other methods (currently 
#' unused).
#'
#' @return
#' If `type = "response"`, a numeric vector of mean predictions.
#' If `type = "quantile"`, a matrix where each row corresponds to an observation
#' in `newdata` and each column to a quantile.
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
  newTestDataPredictionsMatrix <- array(dim = c(nrow(xNewScaled),
                                                numStoredSamples))
  
  # --- Prediction Loop ---
  for (sIdx in 1:numStoredSamples) {
    # --- Inlined testSetPrediction logic ---
    current_tess <- posteriorTessSamples[[sIdx]]
    current_dim  <- posteriorDimSamples[[sIdx]]
    current_pred <- posteriorPredSamples[[sIdx]]
    
    # Get preds for each of the m tessellations in current posterior sample
    predictionList <- lapply(1:mTessellations, function(j) {
      NewTessIndexes <- cellIndices(xNewScaled, current_tess[[j]],
                                    current_dim[[j]])
      current_pred[[j]][NewTessIndexes]
    })
    
    # Sum the predictions from all m tessellations for this posterior sample
    predictionsForSampleS <- rowSums(do.call(cbind, predictionList))
    # --- End of inlined logic ---
    
    newTestDataPredictionsMatrix[, sIdx] <- predictionsForSampleS
  }
  
  # --- Process and Unscale Predictions ---
  if (type == "response") {
    # Calculate the mean of predictions across all posterior samples 
    # (still scaled)
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
#' @param ... Additional graphical parameters to be passed to the `plot` 
#' function (e.g., `pch`, `cex`).
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
    stop(paste0("`x_train` and `y_train` must be provided to plot true",
                " vs. predicted values."))
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