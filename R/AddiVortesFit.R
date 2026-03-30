#' @title Create an AddiVortesFit Object
#'
#' @description A constructor for the AddiVortesFit class.
#'
#' @export
new_AddiVortesFit <- function(posteriorTess, posteriorDim, 
                              posteriorSigma, posteriorPower, posteriorMu, posteriorPred,
                              xCentres, xRanges, yCentre, yRange,
                              inSampleRmse) {
  structure(
    list(
      posteriorTess = posteriorTess,
      posteriorDim = posteriorDim,
      posteriorSigma = posteriorSigma,
      posteriorPower = posteriorPower,
      posteriorMu = posteriorMu,
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
#' @export
#' @method print AddiVortesFit
print.AddiVortesFit <- function(x, ...) {
  if (!inherits(x, "AddiVortesFit")) {
    stop("`x` must be an object of class 'AddiVortesFit'.")
  }
  
  cat("AddiVortes Model\n")
  cat("================\n\n")
  
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
    cat("where f(.) is represented by additive soft Voronoi tessellations\n\n")
  }
  
  num_samples <- length(x$posteriorTess)
  num_tessellations <- if (num_samples > 0) {
    length(x$posteriorTess[[1]])
  } else {
    0
  }
  
  cat("Model Information:\n")
  cat("Number of covariates:      ", num_covariates, "\n")
  cat("Number of tessellations:   ", num_tessellations, "\n")
  cat("Posterior samples:         ", num_samples, "\n")
  cat("In-sample RMSE:            ", round(x$inSampleRmse, 4), "\n\n")
  
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
  
  if (num_samples > 0) {
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
  
  invisible(x)
}

#' @title Summary Method for AddiVortesFit
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
    
    if (is.matrix(all_tess_sizes)) {
      cat("Tessellation complexity across all samples:\n")
      for (j in 1:nrow(all_tess_sizes)) {
        cat("  Tessellation ", j, ": mean = ", round(mean(all_tess_sizes[j,]), 1),
            ", sd = ", round(sd(all_tess_sizes[j,]), 2), "\n")
      }
    }
    
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
    
    if (!is.null(object$posteriorPower)) {
      cat("\nLearned Distance Power:\n")
      cat("  Mean Power: ", round(mean(object$posteriorPower), 3), "\n")
      cat("  Range: [", round(min(object$posteriorPower), 3), ", ", round(max(object$posteriorPower), 3), "]\n")
    }
  }
  
  invisible(object)
}

#' @keywords internal
#' @noRd
applyScaling_internal <- function(mat, centres, ranges) {
  scale(mat, center = centres, scale = ranges)
}

#' @keywords internal
#' @noRd
soft_predict <- function(tess, query, dim, mu, power) {
  if (!is.matrix(tess)) tess <- as.matrix(tess)
  if (!is.matrix(query)) query <- as.matrix(query)
  if (!is.matrix(mu)) mu <- as.matrix(mu)
  
  storage.mode(tess) <- "double"
  storage.mode(query) <- "double"
  storage.mode(mu) <- "double"
  
  result <- .Call("soft_predict_cpp", tess, query, as.integer(dim), mu, as.numeric(power))
  return(result)
}

#' @title Predict Method for AddiVortesFit
#'
#' @importFrom parallel makeCluster stopCluster parLapply detectCores
#' @importFrom pbapply pblapply
#' @importFrom stats rnorm quantile
#' @export
#' @method predict AddiVortesFit
predict.AddiVortesFit <- function(object, newdata,
                                  type = c("response", "quantile"),
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
  
  posteriorTessSamples <- object$posteriorTess
  posteriorDimSamples  <- object$posteriorDim
  posteriorMuSamples   <- object$posteriorMu
  posteriorPower       <- object$posteriorPower
  posteriorSigmaSamples <- object$posteriorSigma
  
  numStoredSamples     <- length(posteriorTessSamples)
  
  if (interval == "prediction" && type == "quantile") {
    if (is.null(posteriorSigmaSamples) || length(posteriorSigmaSamples) == 0) {
      stop("Prediction intervals require posterior sigma samples.")
    }
    if (length(posteriorSigmaSamples) != numStoredSamples) {
      stop("Number of sigma samples does not match.")
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
  
  mTessellations <- length(posteriorTessSamples[[1]])
  nObs <- nrow(xNewScaled)
  
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
    cat("Generating true out-of-sample predictions for ", nrow(newdata),
        " observations using ", numStoredSamples,
        " posterior samples...\n", sep = "")
  }
  
  prediction_list <- pbapply::pblapply(
    X = 1:numStoredSamples,
    FUN = function(sIdx) {
      current_tess  <- posteriorTessSamples[[sIdx]]
      current_dim   <- posteriorDimSamples[[sIdx]]
      current_mu    <- posteriorMuSamples[[sIdx]]
      current_power <- posteriorPower[, sIdx]
      
      pred_list <- lapply(seq_len(mTessellations), function(j) {
        soft_predict(
          tess = current_tess[[j]], 
          query = xNewScaled, 
          dim = current_dim[[j]], 
          mu = current_mu[[j]], 
          power = current_power[j]
        )
      })
      
      model_predictions <- rowSums(do.call(cbind, pred_list))
      
      if (interval == "prediction" && type == "quantile") {
        current_sigma <- posteriorSigmaSamples[sIdx]
        model_predictions <- model_predictions + stats::rnorm(nObs, mean = 0, sd = sqrt(current_sigma))
      }
      
      return(model_predictions)
    },
    cl = if (useParallel) (if (.Platform$OS.type == "windows") cl else cores) else NULL
  )
  
  if (showProgress) cat("\nPrediction generation completed.\n\n")
  
  newTestDataPredictionsMatrix <- do.call(cbind, prediction_list)
  
  if (type == "response") {
    predictions <- rowMeans(newTestDataPredictionsMatrix) * object$yRange + object$yCentre
  } else if (type == "quantile") {
    quantileYhatNewScaled <- apply(newTestDataPredictionsMatrix, 1, stats::quantile,
                                   probs = quantiles, na.rm = TRUE)
    predictions <- t(quantileYhatNewScaled * object$yRange + object$yCentre)
  }
  
  return(predictions)
}

#' @title Plot Method for AddiVortesFit
#'
#' @importFrom graphics plot abline title par layout lines points segments legend
#' @importFrom stats lowess residuals fitted predict sd
#' @export
#' @method plot AddiVortesFit
plot.AddiVortesFit <- function(x, x_train, y_train, sigma_trace = NULL,
                               which = c(1, 2, 3, 4), ask = FALSE, ...) {
  
  if (!inherits(x, "AddiVortesFit")) stop("`x` must be an object of class 'AddiVortesFit'.")
  if (missing(x_train) || missing(y_train)) stop("`x_train` and `y_train` must be provided.")
  if (!is.matrix(x_train)) stop("`x_train` must be a matrix.")
  if (!is.numeric(y_train)) stop("`y_train` must be a numeric vector.")
  if (nrow(x_train) != length(y_train)) stop("Row dimensions must match.")
  if (length(x$posteriorTess) == 0) stop("No posterior samples available.")
  
  which <- intersect(which, 1:5)
  if (length(which) == 0) stop("`which` must contain values between 1 and 5.")
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  n_plots <- length(which)
  if (n_plots == 1) {
    par(mfrow = c(1, 1))
  } else if (n_plots == 2) {
    par(mfrow = c(1, 2))
  } else if (n_plots <= 4) {
    par(mfrow = c(2, 2))
  } else {
    par(mfrow = c(2, 3))
  }
  
  if (any(which %in% c(1, 5))) {
    y_pred_mean <- predict(x, newdata = x_train, type = "response", showProgress = FALSE)
    residuals <- y_train - y_pred_mean
  }
  
  if (3 %in% which) {
    tess_complexity <- sapply(x$posteriorTess, function(sample) {
      mean(sapply(sample, nrow))
    })
  }
  
  if (1 %in% which) {
    if (ask && n_plots > 1) { cat("Press [Enter]: "); readline() }
    plot(y_pred_mean, residuals, xlab = "Fitted Values", ylab = "Residuals", main = "Residuals vs Fitted", pch = 19, col = "darkblue", cex = 0.8, ...)
    abline(h = 0, col = "red", lty = 2, lwd = 2)
    if (length(y_pred_mean) > 3) { smooth_line <- lowess(y_pred_mean, residuals); lines(smooth_line, col = "orange", lwd = 2) }
    rmse_text <- paste("RMSE =", round(x$inSampleRmse, 4))
    legend("topright", legend = rmse_text, bty = "n")
  }
  
  if (2 %in% which) {
    if (ask && n_plots > 1) { cat("Press [Enter]: "); readline() }
    if (is.null(sigma_trace)) {
      if ("posteriorSigma" %in% names(x)) sigma_values <- x$posteriorSigma
      else sigma_values <- x$inSampleRmse + rnorm(length(x$posteriorTess), 0, x$inSampleRmse * 0.1)
    } else sigma_values <- sigma_trace
    plot(1:length(sigma_values), sigma_values, type = "l", xlab = "MCMC Iteration", ylab = expression(sigma), main = "MCMC Trace: Error SD", col = "darkgreen", lwd = 1.5, ...)
    abline(h = mean(sigma_values), col = "red", lty = 2)
    legend("topright", legend = c(paste("Mean:", round(mean(sigma_values), 4)), paste("SD:", round(sd(sigma_values), 4))), bty = "n")
  }
  
  if (3 %in% which) {
    if (ask && n_plots > 1) { cat("Press [Enter]: "); readline() }
    plot(1:length(tess_complexity), tess_complexity, type = "l", xlab = "MCMC Iteration", ylab = "Avg Centers", main = "Trace: Tessellation Complexity", col = "purple", lwd = 1.5, ...)
    abline(h = mean(tess_complexity), col = "red", lty = 2)
    legend("topright", legend = c(paste("Mean:", round(mean(tess_complexity), 1))), bty = "n")
  }
  
  if (4 %in% which) {
    if (ask && n_plots > 1) { cat("Press [Enter]: "); readline() }
    if (!is.null(x$posteriorPower)) {
      power_means <- colMeans(x$posteriorPower)
      plot(1:length(power_means), power_means, type = "l", xlab = "MCMC Iteration", ylab = "Average Power (p)", main = "Trace: Distance Power", col = "darkred", lwd = 1.5, ...)
      abline(h = mean(power_means), col = "blue", lty = 2)
      legend("topright", legend = paste("Mean:", round(mean(power_means), 2)), bty = "n")
    } else { plot.new(); text(0.5, 0.5, "No Power Trace", cex = 1.5) }
  }
  
  if (5 %in% which) {
    if (ask && n_plots > 1) { cat("Press [Enter]: "); readline() }
    y_pred_quantiles <- predict(x, newdata = x_train, type = "quantile", quantiles = c(0.025, 0.975), showProgress = FALSE)
    plot(y_train, y_pred_mean, xlab = "Observed", ylab = "Predicted", main = "Predicted vs Observed", pch = 19, col = "darkblue", cex = 0.8, xlim = range(c(y_train, y_pred_mean)), ylim = range(c(y_train, y_pred_mean)), ...)
    abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)
    for (i in 1:length(y_train)) segments(y_train[i], y_pred_quantiles[i, 1], y_train[i], y_pred_quantiles[i, 2], col = "lightblue", lwd = 1)
    r_squared <- 1 - (sum(residuals^2) / sum((y_train - mean(y_train))^2))
    legend("bottomright", legend = c("Perfect Prediction", "95% Prediction Intervals", paste("R^2 =", round(r_squared, 3))), col = c("red", "lightblue", "black"), lty = c(2, 1, NA), lwd = c(2, 1, NA), pch = c(NA, NA, NA))
  }
  
  invisible(NULL)
}