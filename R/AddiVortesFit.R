#' @title Create an AddiVortesFit Object
#'
#' @description A constructor for the AddiVortesFit class.
#'
#' @export
new_AddiVortesFit <- function(posteriorTess, posteriorDim, 
                              posteriorSigma, posteriorPred,
                              posteriorDirichletWeights = NULL,
                              posteriorVariableSelection = NULL,
                              posteriorAlpha = NULL,
                              predictionMatrix = NULL,
                              posteriorPower = NULL,
                              splitMode = 1,
                              isClassification = FALSE,
                              xCentres, xRanges, yCentre, yRange,
                              inSampleRmse,
                              posteriorDirichletWeightsMean = NULL,
                              posteriorDirichletWeightsLower = NULL,
                              posteriorDirichletWeightsUpper = NULL,
                              ensembleInclusionProbabilities = NULL) {
  structure(
    list(
      posteriorTess = posteriorTess,
      posteriorDim = posteriorDim,
      posteriorSigma = posteriorSigma,
      posteriorPred = posteriorPred,
      posteriorDirichletWeights = posteriorDirichletWeights,
      posteriorVariableSelection = posteriorVariableSelection,
      posteriorAlpha = posteriorAlpha,
      predictionMatrix = predictionMatrix,
      posteriorPower = posteriorPower,
      splitMode = splitMode,
      isClassification = isClassification,
      xCentres = xCentres,
      xRanges = xRanges,
      yCentre = yCentre,
      yRange = yRange,
      inSampleRmse = inSampleRmse,
      posteriorDirichletWeightsMean = posteriorDirichletWeightsMean,
      posteriorDirichletWeightsLower = posteriorDirichletWeightsLower,
      posteriorDirichletWeightsUpper = posteriorDirichletWeightsUpper,
      ensembleInclusionProbabilities = ensembleInclusionProbabilities
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
#' @param x An object of class `AddiVortesFit`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @export
#' @method print AddiVortesFit
print.AddiVortesFit <- function(x, ...) {
  if (!inherits(x, "AddiVortesFit")) {
    stop("`x` must be an object of class 'AddiVortesFit'.")
  }
  
  cat("Adaptive AddiVortes Model (with Dirichlet Sparsity)\n")
  cat("===================================================\n\n")
  
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
  
  num_samples <- length(x$posteriorTess)
  num_tessellations <- if (num_samples > 0) length(x$posteriorTess[[1]]) else 0
  
  cat("Model Information:\n")
  cat("Number of covariates:     ", num_covariates, "\n")
  cat("Number of tessellations:  ", num_tessellations, "\n")
  cat("Posterior samples:        ", num_samples, "\n")
  
  if (isTRUE(x$isClassification)) {
    cat("In-sample Error Rate:     ", round(x$inSampleRmse, 4), "\n\n")
  } else {
    cat("In-sample RMSE:           ", round(x$inSampleRmse, 4), "\n\n")
  }
  
  invisible(x)
}

#' @title Summary Method for AddiVortesFit
#'
#' @description
#' Provides a detailed summary of a fitted `AddiVortesFit` object, including
#' more comprehensive information than the print method, such as variable importance.
#'
#' @param object An object of class `AddiVortesFit`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @importFrom stats sd
#' @export
#' @method summary AddiVortesFit
summary.AddiVortesFit <- function(object, ...) {
  print(object)
  
  if (length(object$posteriorTess) > 0) {
    cat("\nDetailed Posterior Information:\n")
    cat("===============================\n")
    
    all_tess_sizes <- sapply(object$posteriorTess, function(sample) {
      sapply(sample, function(tess) nrow(tess))
    })
    
    cat(sprintf("Average cells per tessellation: %.2f (SD: %.2f)\n", 
                mean(all_tess_sizes), sd(as.vector(all_tess_sizes))))
    
    if (!is.null(object$posteriorDirichletWeights)) {
      cat("\nVariable Selection (Mean Dirichlet Inclusion Weights):\n")
      cat("----------------------------------------------------\n")
      
      mean_weights <- rowMeans(object$posteriorDirichletWeights)
      cov_names <- if(is.null(names(object$xCentres))) paste0("X", 1:length(mean_weights)) else names(object$xCentres)
      
      weight_df <- data.frame(
        Covariate = cov_names,
        MeanWeight = round(mean_weights, 4),
        RelativeImportance = round(mean_weights / max(mean_weights), 3)
      )
      
      weight_df <- weight_df[order(weight_df$MeanWeight, decreasing = TRUE), ]
      
      if(nrow(weight_df) > 15) {
        print(weight_df[1:15, ], row.names = FALSE)
        cat("... and", nrow(weight_df) - 15, "more covariates omitted.\n")
      } else {
        print(weight_df, row.names = FALSE)
      }
    }
  }
  invisible(object)
}

#' @title Predict Method for AddiVortesFit
#'
#' @description
#' Predicts outcomes for new data using a fitted AddiVortesFit model object.
#'
#' @param object An object of class AddiVortesFit.
#' @param newdata A matrix of covariates for the new test set. 
#' @param type The type of prediction required ("response", "quantile", "prob", "class", or "matrix"). 
#' @param quantiles A numeric vector of probabilities to compute.
#' @param interval The type of interval calculation. 
#' @param showProgress Logical indicator for a progress bar.
#' @param parallel Logical indicator for parallel processing.
#' @param cores The number of CPU cores to use.
#' @param ... Further arguments passed to or from other methods.
#'
#' @importFrom parallel makeCluster stopCluster clusterEvalQ detectCores
#' @importFrom pbapply pblapply
#' @importFrom stats rnorm quantile pnorm
#' @export
#' @method predict AddiVortesFit
predict.AddiVortesFit <- function(object, newdata,
                                  type = c("response", "quantile", "prob", "class", "matrix"),
                                  quantiles = c(0.025, 0.975),
                                  interval = c("credible", "prediction"),
                                  showProgress = TRUE,
                                  parallel = TRUE,
                                  cores = NULL,
                                  ...) {
  type <- match.arg(type)
  interval <- match.arg(interval)
  
  if (!inherits(object, "AddiVortesFit")) stop("object must be of class AddiVortesFit.")
  if (!is.matrix(newdata)) stop("newdata must be a matrix.")
  if (ncol(newdata) != length(object$xCentres)) stop("Number of columns does not match.")
  
  posteriorTessSamples  <- object$posteriorTess
  posteriorDimSamples   <- object$posteriorDim
  posteriorPredSamples  <- object$posteriorPred
  posteriorSigmaSamples <- object$posteriorSigma
  
  splitMode <- object$splitMode
  if (is.null(splitMode)) splitMode <- 1 
  posteriorPowerSamples <- object$posteriorPower 
  isClassification <- isTRUE(object$isClassification)
  
  numStoredSamples <- length(posteriorTessSamples)
  
  if (interval == "prediction" && type == "quantile" && !isClassification) {
    if (is.null(posteriorSigmaSamples) || length(posteriorSigmaSamples) == 0) {
      stop("Prediction intervals require posterior sigma samples.")
    }
  }
  
  if (numStoredSamples == 0) {
    warning("The model contains no posterior samples.")
    return(NA_real_)
  }
  
  xNewScaled <- applyScaling_internal(
    mat = newdata,
    centres = object$xCentres,
    ranges = object$xRanges
  )
  xNewScaled <- as.matrix(xNewScaled)
  storage.mode(xNewScaled) <- "double"
  
  k_int <- as.integer(1)
  nObs <- nrow(xNewScaled)
  mTessellations <- length(posteriorTessSamples[[1]])
  
  max_cores_safe <- 1
  if (is.null(cores)) {
    detected_cores <- parallel::detectCores()
    if (is.na(detected_cores)) detected_cores <- 1
    cores <- max(1, min(detected_cores - 1, max_cores_safe))
  } else {
    cores <- max(1, min(cores, max_cores_safe))
  }
  useParallel <- parallel && (cores > 1)
  cl <- NULL
  
  if (useParallel && .Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, { requireNamespace("AddiVortes", quietly = TRUE) })
  }
  
  if (showProgress) {
    cat("Generating predictions for ", nObs,
        " observations using ", numStoredSamples,
        " posterior samples (Mode: ", ifelse(splitMode == 1, "Hard", "Soft"), ")...\n", sep = "")
  }
  
  prediction_list <- pbapply::pblapply(
    X = 1:numStoredSamples,
    FUN = function(sIdx) {
      current_tess <- posteriorTessSamples[[sIdx]]
      current_dim  <- posteriorDimSamples[[sIdx]]
      current_pred <- posteriorPredSamples[[sIdx]]
      
      pred_list <- lapply(seq_len(mTessellations), function(j) {
        n_tess <- nrow(current_tess[[j]])
        
        if (splitMode == 1) {
          if (n_tess == 1L) {
            NewTessIndexes <- rep.int(1L, nObs)
          } else {
            NewTessIndexes <- .Call("knnx_index_predict_cpp", 
                                    current_tess[[j]], 
                                    xNewScaled, 
                                    k_int, 
                                    as.integer(current_dim[[j]]))
            NewTessIndexes <- as.vector(NewTessIndexes)
          }
          return(current_pred[[j]][NewTessIndexes])
          
        } else {
          if (n_tess == 1L) {
            return(rep(current_pred[[j]][1], nObs))
          } else {
            current_p <- posteriorPowerSamples[j, sIdx] 
            pred_vals <- .Call("soft_predict_cpp", 
                               current_tess[[j]], 
                               xNewScaled, 
                               as.integer(current_dim[[j]]),
                               as.numeric(current_pred[[j]]), 
                               as.numeric(current_p))
            return(as.vector(pred_vals))
          }
        }
      })
      
      # Sum across all trees (Latent Z)
      model_predictions <- rowSums(do.call(cbind, pred_list))
      
      if (isClassification) {
        # Transform the latent sum into the probability of being active (Class 1)
        model_predictions <- stats::pnorm(model_predictions)
      } else {
        model_predictions <- model_predictions * object$yRange + object$yCentre
        if (interval == "prediction" && type == "quantile") {
          current_sigma <- posteriorSigmaSamples[sIdx] 
          model_predictions <- model_predictions + stats::rnorm(nObs, mean = 0, sd = sqrt(current_sigma))
        }
      }
      
      model_predictions
    },
    cl = if (useParallel) (if (.Platform$OS.type == "windows") cl else cores) else NULL
  )
  
  if (showProgress) cat("\nPrediction generation completed.\n\n")
  
  # Because of the transformation above, if isClassification == TRUE, 
  # this matrix is now purely probabilities.
  newTestDataPredictionsMatrix <- do.call(cbind, prediction_list)
  
  if (type == "response" || type == "prob") {
    if (!isClassification && type == "prob") {
      stop("Probability predictions are only valid for classification models.")
    }
    # rowMeans of probabilities gives the posterior mean probability
    predictions <- rowMeans(newTestDataPredictionsMatrix)
    
  } else if (type == "class") {
    if (!isClassification) stop("Class predictions are only valid for classification models.")
    prob_predictions <- rowMeans(newTestDataPredictionsMatrix)
    predictions <- ifelse(prob_predictions > 0.5, 1, 0)
    
  } else if (type == "quantile") {
    predictions <- t(apply(newTestDataPredictionsMatrix, 1, stats::quantile,
                           probs = quantiles, na.rm = TRUE))
    
  } else if (type == "matrix") {
    # Returns the full [N x Samples] matrix of probabilities
    predictions <- newTestDataPredictionsMatrix
  }
  
  return(predictions)
}

#' @title Plot Method for AddiVortesFit
#'
#' @description
#' Generates comprehensive diagnostic plots for a fitted `AddiVortesFit` object.
#'
#' @param x An object of class `AddiVortesFit`.
#' @param x_train A matrix of the original training covariates.
#' @param y_train A numeric vector of the original training true outcomes.
#' @param which A numeric vector specifying which plots to generate:
#'   1 = Residuals plot, 2 = Sigma trace, 3 = Tessellation complexity trace,
#'   4 = Predicted vs Observed, 5 = Variable Importance, 6 = Alpha Trace. Default is c(1, 2, 3).
#' @param ask Logical; if TRUE, the user is asked to press Enter before each plot.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @importFrom graphics plot abline title par layout lines points segments legend barplot text boxplot
#' @importFrom stats lowess residuals fitted predict sd
#' @export
#' @method plot AddiVortesFit
plot.AddiVortesFit <- function(x, x_train = NULL, y_train = NULL, 
                               which = c(1, 2, 3), ask = FALSE, ...) {
  
  if (!inherits(x, "AddiVortesFit")) stop("`x` must be an object of class 'AddiVortesFit'.")
  which <- intersect(which, 1:6)
  if (length(which) == 0) stop("`which` must contain values between 1 and 6.")
  
  if (any(which %in% c(1, 4)) && (is.null(x_train) || is.null(y_train))) {
    stop("`x_train` and `y_train` must be provided to plot residuals or predicted vs observed (Plots 1 and 4).")
  }
  
  isClassification <- isTRUE(x$isClassification)
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  n_plots <- length(which)
  if (n_plots == 1) { par(mfrow = c(1, 1)) } 
  else if (n_plots == 2) { par(mfrow = c(1, 2)) } 
  else if (n_plots <= 4) { par(mfrow = c(2, 2)) }
  else { par(mfrow = c(2, 3)) }
  
  if (any(which %in% c(1, 4))) {
    y_pred_mean <- predict(x, newdata = x_train, type = "response", showProgress = FALSE)
    if (!isClassification) {
      residuals <- y_train - y_pred_mean
    }
  }
  
  if (1 %in% which) {
    if (ask && n_plots > 1) { cat("Press [Enter] to see residuals plot: "); readline() }
    if (isClassification) {
      boxplot(y_pred_mean ~ y_train, col = c("lightcoral", "lightblue"),
              xlab = "Observed Class", ylab = "Predicted Probability",
              main = "Predicted Probabilities by Class")
    } else {
      plot(y_pred_mean, residuals, xlab = "Fitted Values", ylab = "Residuals",
           main = "Residuals vs Fitted", pch = 19, col = "darkblue", cex = 0.8, ...)
      abline(h = 0, col = "red", lty = 2, lwd = 2)
      if (length(y_pred_mean) > 3) lines(lowess(y_pred_mean, residuals), col = "orange", lwd = 2)
      legend("topright", legend = paste("RMSE =", round(x$inSampleRmse, 4)), bty = "n")
    }
  }
  
  if (2 %in% which) {
    if (ask && n_plots > 1) { cat("Press [Enter] to see sigma trace plot: "); readline() }
    if (isClassification) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "Trace: Error Standard Deviation")
      text(1, 1, "Sigma is fixed to 1\n(Probit Classification Mode)")
    } else {
      sigma_values <- if (!is.null(x$posteriorSigma)) sqrt(x$posteriorSigma) else rep(NA, length(x$posteriorTess))
      plot(1:length(sigma_values), sigma_values, type = "l", xlab = "MCMC Iteration",
           ylab = expression(sigma), main = "Trace: Error Standard Deviation",
           col = "darkgreen", lwd = 1.5, ...)
      abline(h = mean(sigma_values, na.rm = TRUE), col = "red", lty = 2)
    }
  }
  
  if (3 %in% which) {
    if (ask && n_plots > 1) { cat("Press [Enter] to see tessellation complexity trace: "); readline() }
    tess_complexity <- sapply(x$posteriorTess, function(sample) mean(sapply(sample, nrow)))
    plot(1:length(tess_complexity), tess_complexity, type = "l", xlab = "MCMC Iteration",
         ylab = "Average Centers", main = "Trace: Tessellation Complexity",
         col = "purple", lwd = 1.5, ...)
    abline(h = mean(tess_complexity), col = "red", lty = 2)
  }
  
  if (4 %in% which) {
    if (ask && n_plots > 1) { cat("Press [Enter] to see predicted vs observed plot: "); readline() }
    if (isClassification) {
      plot(jitter(y_train, factor = 0.5), y_pred_mean, xlab = "Observed Class (Jittered)", 
           ylab = "Predicted Probability", main = "Predicted vs Observed", 
           pch = 19, col = rgb(0, 0, 0.5, 0.3), ...)
      abline(h = 0.5, col = "red", lty = 2)
    } else {
      y_pred_quantiles <- predict(x, newdata = x_train, type = "quantile", quantiles = c(0.025, 0.975), showProgress = FALSE)
      plot(y_train, y_pred_mean, xlab = "Observed Values", ylab = "Predicted Values",
           main = "Predicted vs Observed", pch = 19, col = "darkblue", cex = 0.8,
           xlim = range(c(y_train, y_pred_mean)), ylim = range(c(y_train, y_pred_mean)), ...)
      abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)
      for (i in 1:length(y_train)) {
        segments(y_train[i], y_pred_quantiles[i, 1], y_train[i], y_pred_quantiles[i, 2], col = "lightblue", lwd = 1)
      }
    }
  }
  
  if (5 %in% which) {
    if (ask && n_plots > 1) { cat("Press [Enter] to see variable importance plot: "); readline() }
    
    if (!is.null(x$posteriorDirichletWeightsMean) && 
        !is.null(x$posteriorDirichletWeightsLower) && 
        !is.null(x$posteriorDirichletWeightsUpper)) {
      
      old_mar <- par("mar")
      
      mean_weights <- x$posteriorDirichletWeightsMean
      lower_weights <- x$posteriorDirichletWeightsLower
      upper_weights <- x$posteriorDirichletWeightsUpper
      
      cov_names <- if(!is.null(names(x$xCentres))) names(x$xCentres) else paste0("X", 1:length(mean_weights))
      
      ord <- order(mean_weights, decreasing = FALSE)
      
      max_vars_to_plot <- 20
      if (length(mean_weights) > max_vars_to_plot) {
        ord <- tail(ord, max_vars_to_plot)
        main_title <- paste("Top", max_vars_to_plot, "Variable Importance (95% CI)")
      } else {
        main_title <- "Variable Importance (95% CI)"
      }
      
      max_name_len <- max(nchar(cov_names[ord]))
      par(mar = c(5, max(4.1, max_name_len * 0.6), 4, 2) + 0.1)
      
      bp <- barplot(mean_weights[ord], horiz = TRUE, names.arg = cov_names[ord], las = 1,
                    col = "steelblue", border = NA, xlab = "Posterior Inclusion Weight",
                    main = main_title, cex.names = 0.8, 
                    xlim = c(0, max(upper_weights[ord]) * 1.05))
      
      segments(x0 = lower_weights[ord], y0 = bp, x1 = upper_weights[ord], y1 = bp, 
               col = "black", lwd = 1.5)
      
      epsilon <- 0.15
      segments(x0 = lower_weights[ord], y0 = bp - epsilon, x1 = lower_weights[ord], y1 = bp + epsilon, 
               col = "black", lwd = 1.5)
      segments(x0 = upper_weights[ord], y0 = bp - epsilon, x1 = upper_weights[ord], y1 = bp + epsilon, 
               col = "black", lwd = 1.5)
      
      abline(v = 1/length(mean_weights), col = "red", lty = 2, lwd = 1.5) 
      
      par(mar = old_mar)
      
    } else {
      plot(1, type="n", axes=FALSE, xlab="", ylab="", main="Variable Importance Missing")
      text(1, 1, "Dirichlet Weight Summaries Not Found")
    }
  }
  
  if (6 %in% which) {
    if (ask && n_plots > 1) { cat("Press [Enter] to see alpha trace plot: "); readline() }
    if (!is.null(x$posteriorAlpha)) {
      plot(1:length(x$posteriorAlpha), x$posteriorAlpha, type = "l", xlab = "MCMC Iteration",
           ylab = expression(alpha), main = "Trace: Dirichlet Concentration",
           col = "darkorange", lwd = 1.5, ...)
      abline(h = mean(x$posteriorAlpha, na.rm=TRUE), col = "red", lty = 2)
    } else {
      plot(1, type="n", axes=FALSE, xlab="", ylab="", main="Alpha Trace Missing")
      text(1, 1, "No Alpha Trace Found")
    }
  }
  
  invisible(NULL)
}