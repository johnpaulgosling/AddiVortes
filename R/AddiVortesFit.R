#' @title Create an AddiVortesFit Object
#'
#' @description A constructor for the AddiVortesFit class.
#'
#' @param posteriorTess A list of the posterior samples of the tessellations.
#' @param posteriorDim A list of the posterior samples of the dimensions.
#' @param posteriorSigma A list of the posterior samples of the error variance.
#' @param posteriorPred A list of the posterior samples of the predictions.
#' @param xCentres The centres of the covariates.
#' @param xRanges The ranges of the covariates.
#' @param yCentre The centre of the output values.
#' @param yRange The range of the output values.
#' @param inSampleRmse The in-sample RMSE.
#'
#' @return An object of class AddiVortesFit.
#' @export
new_AddiVortesFit <- function(posteriorTess, posteriorDim, 
                              posteriorSigma, posteriorPred,
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
#' Generates comprehensive diagnostic plots for a fitted `AddiVortesFit` object.
#' This function creates multiple diagnostic plots including residuals,
#' MCMC traces for sigma, and tessellation complexity over iterations.
#'
#' @param x An object of class `AddiVortesFit`, typically the result of a
#'   call to `AddiVortes()`.
#' @param x_train A matrix of the original training covariates.
#' @param y_train A numeric vector of the original training true outcomes.
#' @param sigma_trace An optional numeric vector of sigma values from MCMC samples.
#'   If not provided, the method will attempt to extract it from the model object.
#' @param which A numeric vector specifying which plots to generate:
#'   1 = Residuals plot, 2 = Sigma trace, 3 = Tessellation complexity trace,
#'   4 = Predicted vs Observed. Default is c(1, 2, 3).
#' @param ask Logical; if TRUE, the user is asked to press Enter before each plot.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return
#' This function is called for its side effect of creating plots and returns
#' `NULL` invisibly.
#'
#' @details
#' The function generates up to four diagnostic plots:
#' \enumerate{
#'   \item \strong{Residuals Plot}: Residuals vs fitted values with smoothed trend line
#'   \item \strong{Sigma Trace}: MCMC trace plot for the error variance parameter
#'   \item \strong{Tessellation Complexity}: Trace of average tessellation size over iterations
#'   \item \strong{Predicted vs Observed}: Scatter plot with confidence intervals
#' }
#'
#' @importFrom graphics plot abline title par layout lines points segments legend
#' @importFrom stats lowess residuals fitted
#' @export
#' @method plot AddiVortesFit
#'
#' @examples
#' \dontrun{
#' # Assuming 'fit' is a trained AddiVortesFit object
#' plot(fit, x_train = x_train_data, y_train = y_train_data)
#' 
#' # Show only specific plots
#' plot(fit, x_train = x_train_data, y_train = y_train_data, which = c(1, 3))
#' 
#' # With custom sigma trace
#' plot(fit, x_train = x_train_data, y_train = y_train_data, 
#'      sigma_trace = my_sigma_samples)
#' }
plot.AddiVortesFit <- function(x, x_train, y_train, sigma_trace = NULL,
                               which = c(1, 2, 3), ask = FALSE, ...) {
  
  # --- Input Validation ---
  if (!inherits(x, "AddiVortesFit")) {
    stop("`x` must be an object of class 'AddiVortesFit'.")
  }
  if (missing(x_train) || missing(y_train)) {
    stop("`x_train` and `y_train` must be provided for diagnostic plots.")
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
  if (length(x$posteriorTess) == 0) {
    stop("No posterior samples available for plotting.")
  }
  
  # Validate which parameter
  which <- intersect(which, 1:4)
  if (length(which) == 0) {
    stop("`which` must contain values between 1 and 4.")
  }
  
  # Store original par settings
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Set up plotting layout
  n_plots <- length(which)
  if (n_plots == 1) {
    par(mfrow = c(1, 1))
  } else if (n_plots == 2) {
    par(mfrow = c(1, 2))
  } else if (n_plots <= 4) {
    par(mfrow = c(2, 2))
  }
  
  # Generate predictions for residuals analysis
  if (any(which %in% c(1, 4))) {
    y_pred_mean <- predict(x, newdata = x_train, type = "response")
    residuals <- y_train - y_pred_mean
  }
  
  # Calculate tessellation statistics across samples
  if (3 %in% which) {
    tess_complexity <- sapply(x$posteriorTess, function(sample) {
      mean(sapply(sample, nrow))
    })
  }
  
  # --- Plot 1: Residuals ---
  if (1 %in% which) {
    if (ask && n_plots > 1) {
      cat("Press [Enter] to see residuals plot: ")
      readline()
    }
    
    plot(y_pred_mean, residuals,
         xlab = "Fitted Values",
         ylab = "Residuals",
         main = "Residuals vs Fitted",
         pch = 19, col = "darkblue", cex = 0.8,
         ...)
    
    # Add horizontal line at y = 0
    abline(h = 0, col = "red", lty = 2, lwd = 2)
    
    # Add smoothed trend line
    if (length(y_pred_mean) > 3) {
      smooth_line <- lowess(y_pred_mean, residuals)
      lines(smooth_line, col = "orange", lwd = 2)
    }
    
    # Add RMSE annotation
    rmse_text <- paste("RMSE =", round(x$inSampleRmse, 4))
    legend("topright", legend = rmse_text, bty = "n")
  }
  
  # --- Plot 2: Sigma Trace ---
  if (2 %in% which) {
    if (ask && n_plots > 1) {
      cat("Press [Enter] to see sigma trace plot: ")
      readline()
    }
    
    # Try to extract sigma from the model object or use provided trace
    if (is.null(sigma_trace)) {
      # Attempt to extract sigma from model object
      if ("posteriorSigma" %in% names(x)) {
        sigma_values <- x$posteriorSigma
      } else {
        # If no sigma trace available, create a placeholder
        warning("No sigma trace found. Creating synthetic trace for demonstration.")
        sigma_values <- x$inSampleRmse + rnorm(length(x$posteriorTess), 0, x$inSampleRmse * 0.1)
      }
    } else {
      sigma_values <- sigma_trace
    }
    
    # Ensure sigma_values has same length as posterior samples
    if (length(sigma_values) != length(x$posteriorTess)) {
      warning("Length of sigma trace doesn't match number of posterior samples.")
      sigma_values <- rep(sigma_values[1], length(x$posteriorTess))
    }
    
    plot(1:length(sigma_values), sigma_values,
         type = "l",
         xlab = "MCMC Iteration",
         ylab = expression(sigma),
         main = "MCMC Trace: Error Standard Deviation",
         col = "darkgreen", lwd = 1.5,
         ...)
    
    # Add horizontal line at mean
    abline(h = mean(sigma_values), col = "red", lty = 2)
    
    # Add convergence statistics
    sigma_mean <- round(mean(sigma_values), 4)
    sigma_sd <- round(sd(sigma_values), 4)
    legend("topright", 
           legend = c(paste("Mean:", sigma_mean),
                      paste("SD:", sigma_sd)),
           bty = "n")
  }
  
  # --- Plot 3: Tessellation Complexity Trace ---
  if (3 %in% which) {
    if (ask && n_plots > 1) {
      cat("Press [Enter] to see tessellation complexity trace: ")
      readline()
    }
    
    plot(1:length(tess_complexity), tess_complexity,
         type = "l",
         xlab = "MCMC Iteration",
         ylab = "Average Number of Tessellation Centers",
         main = "MCMC Trace: Tessellation Complexity",
         col = "purple", lwd = 1.5,
         ...)
    
    # Add horizontal line at mean
    abline(h = mean(tess_complexity), col = "red", lty = 2)
    
    # Add summary statistics
    complexity_mean <- round(mean(tess_complexity), 1)
    complexity_range <- round(range(tess_complexity), 1)
    legend("topright", 
           legend = c(paste("Mean:", complexity_mean),
                      paste("Range: [", complexity_range[1], ",", complexity_range[2], "]")),
           bty = "n")
  }
  
  # --- Plot 4: Predicted vs Observed ---
  if (4 %in% which) {
    if (ask && n_plots > 1) {
      cat("Press [Enter] to see predicted vs observed plot: ")
      readline()
    }
    
    # Get quantile predictions for uncertainty
    y_pred_quantiles <- predict(x, newdata = x_train, type = "quantile",
                                quantiles = c(0.025, 0.975))
    
    # Create the scatter plot
    plot(y_train, y_pred_mean,
         xlab = "Observed Values",
         ylab = "Predicted Values",
         main = "Predicted vs Observed",
         pch = 19, col = "darkblue", cex = 0.8,
         xlim = range(c(y_train, y_pred_mean)),
         ylim = range(c(y_train, y_pred_mean)),
         ...)
    
    # Add the line of equality (perfect prediction)
    abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)
    
    # Add uncertainty intervals
    for (i in 1:length(y_train)) {
      segments(y_train[i], y_pred_quantiles[i, 1],
               y_train[i], y_pred_quantiles[i, 2],
               col = "lightblue", lwd = 1)
    }
    
    # Calculate and display R-squared
    ss_res <- sum(residuals^2)
    ss_tot <- sum((y_train - mean(y_train))^2)
    r_squared <- 1 - (ss_res / ss_tot)
    
    legend("bottomright",
           legend = c("Perfect Prediction",
                      "95% Prediction Intervals",
                      paste("R² =", round(r_squared, 3))),
           col = c("red", "lightblue", "black"),
           lty = c(2, 1, NA),
           lwd = c(2, 1, NA),
           pch = c(NA, NA, NA))
  }
  
  # Return invisibly
  invisible(NULL)
}