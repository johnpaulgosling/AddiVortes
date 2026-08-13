#' @title Create an AddiVortes Object
#'
#' @description A constructor for the AddiVortes class.
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
#' @param metric The metric used for scaling covariates (default "E" for Euclidean).
#' @param members The membership vector for the covariates
#' @param metric_aug The augmented metric after categorical variables are converted
#' to one-hot
#' @param member_aug The membership vector corresponding to metric_aug
#' @param catEncoding Optional list of categorical encoding metadata returned by
#'   \code{encodeCategories_internal}, or \code{NULL} if no categorical covariates
#'   were present.
#' @param traceStats Optional data frame of per-iteration MCMC trace statistics.
#' @param task The modelling task: `"regression"`, `"binary"` or `"multinomial"`.
#' @param classLevels Character vector of class labels for classification fits,
#'   or `NULL` for regression.
#' @param nLatents Number of latent probit dimensions. 1 for regression and
#'   binary classification; \eqn{K-1}{K-1} for \eqn{K}{K}-class multinomial
#'   models.
#' @param mPerLatent Number of tessellations per latent ensemble.
#' @param inSampleAccuracy In-sample classification accuracy, or `NA` for
#'   regression.
#' @param inSampleBrier In-sample Brier score for binary classification, or `NA`
#'   otherwise.
#'
#' @return An object of class AddiVortes.
#' @export
new_AddiVortes <- function(posteriorTess, posteriorDim,
                           posteriorSigma, posteriorPred,
                           xCentres, xRanges, yCentre, yRange,
                           inSampleRmse, metric = "E",
                           members = rep(1, length(xCentres)),
                           metric_aug = "E",
                           member_aug = rep(1, length(xCentres)),
                           catEncoding = NULL,
                           traceStats = NULL,
                           task = "regression",
                           classLevels = NULL,
                           nLatents = 1L,
                           mPerLatent = NA_integer_,
                           inSampleAccuracy = NA_real_,
                           inSampleBrier = NA_real_) {
  member_length <- sapply(unique(member_aug), function(x) sum(member_aug == x))
  metric_type <- sapply(unique(member_aug), function(x) metric_aug[which(member_aug == x)[1]])
  if (is.na(mPerLatent)) {
    n_tess <- if (length(posteriorTess) > 0) length(posteriorTess[[1]]) else 0L
    mPerLatent <- as.integer(n_tess / max(1L, nLatents))
  }
  structure(
    list(
      posteriorTess = posteriorTess,
      posteriorDim = posteriorDim,
      posteriorSigma = posteriorSigma,
      posteriorPred = posteriorPred,
      xCentres = xCentres,
      xRanges = xRanges,
      yCentre = yCentre,
      yRange = yRange,
      inSampleRmse = inSampleRmse,
      metric = metric,
      members = members,
      metric_red = metric_type,
      member_red = member_length,
      catEncoding = catEncoding,
      traceStats = traceStats,
      task = task,
      classLevels = classLevels,
      nLatents = as.integer(nLatents),
      mPerLatent = as.integer(mPerLatent),
      inSampleAccuracy = inSampleAccuracy,
      inSampleBrier = inSampleBrier
    ),
    class = "AddiVortes"
  )
}

#' @title Print Method for AddiVortes
#'
#' @description
#' Prints a summary of a fitted `AddiVortes` object, providing information
#' about the model structure, dimensions, and fit quality similar to the
#' output of a linear model summary.
#'
#' @param x An object of class `AddiVortes`, typically the result of a
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
#' @method print AddiVortes
print.AddiVortes <- function(x, ...) {
  # --- Input Validation ---
  if (!inherits(x, "AddiVortes")) {
    stop("`x` must be an object of class 'AddiVortes'.")
  }

  is_class <- isClassification_internal(x)
  if (isTRUE(x$task == "binary")) {
    cat("AddiVortes Binary Classification Model\n")
  } else if (isTRUE(x$task == "multinomial")) {
    cat("AddiVortes Multinomial Classification Model\n")
  } else {
    cat("AddiVortes Model\n")
  }
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
  if (is_class) {
    cat("Task:                     ", x$task, "\n")
    cat("Classes:                  ", paste(x$classLevels, collapse = ", "), "\n")
    if (!is.null(x$nLatents) && x$nLatents > 1L) {
      cat("Latent dimensions:        ", x$nLatents, "\n")
    }
    if (!is.na(x$inSampleAccuracy)) {
      cat("In-sample accuracy:       ", round(x$inSampleAccuracy, 4), "\n")
    }
    if (!is.na(x$inSampleBrier)) {
      cat("In-sample Brier score:    ", round(x$inSampleBrier, 4), "\n")
    }
    cat("\n")
  } else {
    cat("In-sample RMSE:           ", round(x$inSampleRmse, 4), "\n\n")
  }

  # Scaling information
  cat("Covariate Scaling:\n")
  scaling_df <- data.frame(
    Covariate = covariate_names,
    Centre = round(x$xCentres, 4),
    Range = round(x$xRanges, 4)
  )
  print(scaling_df, row.names = FALSE)

  if (is_class) {
    cat("\nLatent scale: probit (sigma = 1; response is not scaled)\n\n")
  } else {
    cat("\nOutput Scaling:\n")
    cat("Centre: ", round(x$yCentre, 4), "\n")
    cat("Range:  ", round(x$yRange, 4), "\n\n")
  }

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
    cat(
      "Range of cells per tessellation: [", range_tess_size[1], ", ",
      range_tess_size[2], "]\n"
    )
  } else {
    cat("No posterior samples available.\n")
  }

  # Return the object invisibly
  invisible(x)
}

#' @title Summary Method for AddiVortes
#'
#' @description
#' Provides a detailed summary of a fitted `AddiVortes` object, including
#' more comprehensive information than the print method.
#'
#' @param object An object of class `AddiVortes`, typically the result of a
#'   call to `AddiVortes()`.
#' @param ... Further arguments passed to or from other methods (currently
#' unused).
#'
#' @return
#' The function is called for its side effect of printing detailed model
#' information and returns the input object `object` invisibly.
#'
#' @importFrom stats sd
#' @export
#' @method summary AddiVortes
summary.AddiVortes <- function(object, ...) {
  # --- Input Validation ---
  if (!inherits(object, "AddiVortes")) {
    stop("`object` must be an object of class 'AddiVortes'.")
  }

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
      for (j in seq_len(nrow(all_tess_sizes))) {
        cat(
          "  Tessellation ", j, ": mean = ", round(mean(all_tess_sizes[j, ]), 1),
          ", sd = ", round(sd(all_tess_sizes[j, ]), 2), "\n"
        )
      }
    }

    # Dimension information if available
    if (length(object$posteriorDim) > 0) {
      dim_info <- sapply(object$posteriorDim, function(sample) {
        sapply(sample, length)
      })

      if (is.matrix(dim_info)) {
        cat("\nActive dimensions per tessellation:\n")
        for (j in seq_len(nrow(dim_info))) {
          cat(
            "  Tessellation ", j, ": mean = ", round(mean(dim_info[j, ]), 1),
            " dimensions\n"
          )
        }
      }
    }
  }

  invisible(object)
}

#' @title Predict Method for AddiVortes
#'
#' @description
#' Predicts outcomes for new data using a fitted `AddiVortes` model object.
#' Regression fits return means or quantiles of the response. Classification
#' fits return class probabilities, class labels, latent-scale values, or
#' quantiles of the class probabilities.
#'
#' @param object An object of class `AddiVortes`, typically the result of a
#'   call to `AddiVortes()`.
#' @param newdata A matrix of covariates for the new test set. The number of
#'   columns must match the original training data.
#' @param type The type of prediction required. The default `"response"` gives
#'   the mean prediction (class probabilities for classification). `"quantile"`
#'   returns the quantiles specified by `quantiles`. `"class"` returns predicted
#'   class labels (classification only). `"link"` returns the latent sum of
#'   tessellations \eqn{G(x)}{G(x)}.
#' @param quantiles A numeric vector of probabilities to
#'   compute for the predictions when `type = "quantile"`.
#' @param interval The type of interval calculation. The default `"credible"`
#'   accounts only for uncertainty in the mean (similar to `lm`'s confidence interval).
#'   The alternative `"prediction"` also includes the model's error variance,
#'   producing wider intervals (similar to `lm`'s prediction interval). Not used
#'   for classification models.
#' @param showProgress Logical; if TRUE, a progress bar is shown during prediction.
#' @param ... Further arguments passed to or from other methods (currently
#' unused).
#'
#' @return
#' If `type = "response"`, a numeric vector of mean predictions for regression
#' or binary classification, or an \eqn{n \times K}{n x K} probability matrix
#' for multinomial classification. If `type = "quantile"`, a matrix of quantiles
#' (binary/regression) or a named list of such matrices (multinomial). If
#' `type = "class"`, a factor of predicted labels. If `type = "link"`, the
#' latent function \eqn{G(x)}{G(x)}.
#'
#' @details
#' This function relies on the internal helper function `applyScaling_internal`
#' being available in the environment, which is used by the main
#' `AddiVortes` function.
#'
#' Predictions traverse all retained draws and tessellations in a single C++
#' call, avoiding repeated R/C++ boundary crossings per tessellation.
#'
#' When `interval = "prediction"` and `type = "quantile"`, the function samples
#' additional Gaussian noise with variance equal to the sampled sigma squared
#' from the posterior. This accounts for the inherent variability in individual
#' predictions, not just uncertainty in the mean function. The noise is added
#' in the scaled space before unscaling predictions. Classification uses a
#' probit link with residual variance fixed at 1, so prediction intervals are
#' not defined; use `type = "quantile"` for credible intervals on probabilities.
#'
#' For binary classification, `"response"` is the posterior mean of
#' \eqn{\Phi(G^{(s)}(x))}{Phi(G^(s)(x))}. For multinomial classification, class
#' probabilities are estimated from independent \eqn{N(G, I)}{N(G, I)} latents,
#' with the first class as the reference.
#'
#' @examples
#' \donttest{
#' # Fit a model
#' set.seed(123)
#' X <- matrix(rnorm(100), 20, 5)
#' Y <- rnorm(20)
#' fit <- AddiVortes(Y, X, m = 5, totalMCMCIter = 50, mcmcBurnIn = 10)
#'
#' # New data for prediction
#' X_new <- matrix(rnorm(25), 5, 5)
#'
#' # Mean predictions
#' pred_mean <- predict(fit, X_new, type = "response")
#'
#' # Credible intervals (uncertainty in mean only)
#' pred_conf <- predict(fit, X_new,
#'   type = "quantile",
#'   interval = "credible",
#'   quantiles = c(0.025, 0.975)
#' )
#'
#' # Prediction intervals (includes error variance)
#' pred_pred <- predict(fit, X_new,
#'   type = "quantile",
#'   interval = "prediction",
#'   quantiles = c(0.025, 0.975)
#' )
#'
#' # Prediction intervals are wider than credible intervals
#' mean(pred_pred[, 2] - pred_pred[, 1]) > mean(pred_conf[, 2] - pred_conf[, 1])
#' }
#'
#' @importFrom stats rnorm quantile pnorm
#' @export
#' @method predict AddiVortes
predict.AddiVortes <- function(object, newdata,
                               type = c("response", "quantile", "class", "link"),
                               quantiles = c(0.025, 0.975),
                               interval = c("credible", "prediction"),
                               showProgress = interactive(),
                               ...) {
  type <- match.arg(type)
  interval <- match.arg(interval)

  # --- Input validation ---
  if (!inherits(object, "AddiVortes")) {
    stop("`object` must be of class 'AddiVortes'.")
  }
  if (!is.matrix(newdata) && !is.data.frame(newdata)) {
    stop("`newdata` must be a matrix or data frame.")
  }

  is_class <- isClassification_internal(object)
  if (type == "class" && !is_class) {
    stop("`type = \"class\"` is only valid for classification models.",
         call. = FALSE)
  }
  if (is_class && interval == "prediction") {
    stop(
      "Prediction intervals are not used for classification models; ",
      "use type = \"quantile\" for credible intervals on class probabilities.",
      call. = FALSE
    )
  }

  # Apply categorical encoding if the model was trained with categorical covariates
  if (!is.null(object$catEncoding)) {
    if (ncol(newdata) != object$catEncoding$origNCols) {
      stop("Number of columns in `newdata` does not match the original training data.")
    }
    san_data <- covariateStructure_internal(newdata, object$metric, object$members)
    newdata <- san_data$data
    encResult <- encodeCategories_internal(newdata, encoding = object$catEncoding)
    newdata <- encResult$encoded
  } else {
    if (ncol(newdata) != length(object$xCentres)) {
      stop("Number of columns in `newdata` does not match the original training data.")
    }
    san_data <- covariateStructure_internal(newdata, object$metric, object$members, one.hot = FALSE)
    newdata <- as.matrix(san_data$data)
  }

  posteriorTessSamples <- object$posteriorTess
  posteriorDimSamples <- object$posteriorDim
  posteriorPredSamples <- object$posteriorPred
  posteriorSigmaSamples <- object$posteriorSigma
  numStoredSamples <- length(posteriorTessSamples)

  # Validate sigma samples for prediction intervals
  if (interval == "prediction" && type == "quantile") {
    if (is.null(posteriorSigmaSamples) || length(posteriorSigmaSamples) == 0) {
      stop("Prediction intervals require posterior sigma samples, which are not available in this model object.")
    }
    if (length(posteriorSigmaSamples) != numStoredSamples) {
      stop("Number of sigma samples does not match number of posterior samples.")
    }
  }

  if (numStoredSamples == 0) {
    warning("The AddiVortes model contains no posterior samples. Cannot make predictions.")
    return(NA_real_)
  }

  # Scale new data if required
  xNewScaled <- applyScaling_internal(
    mat = newdata,
    centres = object$xCentres,
    ranges = object$xRanges
  )
  # Leave spherical columns on their original radian scale, exactly as the fit
  # function does. The mask must be aligned to the (potentially reordered and
  # augmented) column order of `newdata`, not the original user column order in
  # `object$metric`: `covariateStructure_internal` may reorder covariates (e.g.
  # move Euclidean covariates ahead of spherical ones). The per-column augmented
  # metric is recovered from the reduced representation used for distances, which
  # is already aligned to the reordered columns.
  metricAug <- rep(object$metric_red, object$member_red)
  xNewScaled[, metricAug != 0] <- newdata[, metricAug != 0]
  # Binary columns from categorical encoding are kept at their encoded values
  # (0 or catScaling) rather than being further scaled
  if (!is.null(object$catEncoding)) {
    binaryCols <- object$catEncoding$encodedBinaryCols
    xNewScaled[, binaryCols] <- newdata[, binaryCols]
  }

  nObs <- nrow(xNewScaled)

  if (showProgress) {
    cat("Generating predictions for ", nrow(newdata),
      " observations using ", numStoredSamples,
      " posterior samples...\n",
      sep = ""
    )
  }

  if (is_class) {
    G_list <- latentLinkMatrices_internal(object, xNewScaled, showProgress)
    if (showProgress) cat("Done.\n\n")
    return(summariseClassificationPredictions_internal(
      object, G_list, type, quantiles
    ))
  }

  # Single compiled pass over all draws and tessellations
  newTestDataPredictionsMatrix <- .Call(
    "addi_vortes_predict_cpp",
    xNewScaled,
    posteriorTessSamples,
    posteriorDimSamples,
    posteriorPredSamples,
    as.integer(object$metric_red),
    as.integer(object$member_red),
    as.logical(showProgress)
  )

  # Add Gaussian noise for prediction intervals when computing quantiles
  if (interval == "prediction" && type == "quantile") {
    for (sIdx in seq_len(numStoredSamples)) {
      current_sigma <- posteriorSigmaSamples[sIdx]
      # Note: sigma is stored as sigma^2 (variance), so we need sqrt for sd
      newTestDataPredictionsMatrix[, sIdx] <-
        newTestDataPredictionsMatrix[, sIdx] +
        rnorm(nObs, mean = 0, sd = sqrt(current_sigma))
    }
  }

  if (showProgress) cat("Done.\n\n")

  # --- Unscale and summarise predictions ---
  if (type == "link" || type == "response") {
    predictions <- rowMeans(newTestDataPredictionsMatrix) * object$yRange + object$yCentre
  } else if (type == "quantile") {
    quantileYhatNewScaled <- apply(newTestDataPredictionsMatrix, 1, quantile,
      probs = quantiles, na.rm = TRUE
    )
    predictions <- t(quantileYhatNewScaled * object$yRange + object$yCentre)
  }

  return(predictions)
}

#' @title Summarise classification predictions from latent G draws
#'
#' @keywords internal
#' @noRd
summariseClassificationPredictions_internal <- function(object, G_list, type,
                                                        quantiles) {
  class_levels <- object$classLevels
  if (object$task == "binary") {
    g_mat <- G_list[[1]]
    p_mat <- stats::pnorm(g_mat)
    mean_p <- rowMeans(p_mat)
    if (type == "link") {
      return(rowMeans(g_mat))
    }
    if (type == "response") {
      return(mean_p)
    }
    if (type == "class") {
      idx <- ifelse(mean_p > 0.5, 2L, 1L)
      return(factor(class_levels[idx], levels = class_levels))
    }
    quantile_p <- apply(p_mat, 1, stats::quantile, probs = quantiles, na.rm = TRUE)
    return(t(quantile_p))
  }

  n_latents <- length(G_list)
  if (type == "link") {
    link_mat <- vapply(G_list, rowMeans, numeric(nrow(G_list[[1]])))
    colnames(link_mat) <- paste0("latent", seq_len(n_latents))
    return(link_mat)
  }

  probs <- multinomialProbabilities_internal(
    G_list,
    nMC = 32L,
    perDraw = type == "quantile"
  )
  colnames(probs$mean) <- class_levels
  if (type == "response") {
    return(probs$mean)
  }
  if (type == "class") {
    idx <- max.col(probs$mean, ties.method = "first")
    return(factor(class_levels[idx], levels = class_levels))
  }

  quantile_list <- vector("list", length(class_levels))
  names(quantile_list) <- class_levels
  n_draw <- dim(probs$draw)[3]
  n_obs <- nrow(probs$mean)
  for (k in seq_along(class_levels)) {
    draw_k <- matrix(probs$draw[, k, ], nrow = n_obs, ncol = n_draw)
    quantile_list[[k]] <- t(apply(draw_k, 1, stats::quantile,
                                  probs = quantiles, na.rm = TRUE))
  }
  quantile_list
}

extractErrorStandardDeviationTrace_internal <- function(x, sigma_trace = NULL,
                                                        expected_length = length(x$posteriorTess)) {
  if (is.null(sigma_trace)) {
    if ("posteriorSigma" %in% names(x) && length(x$posteriorSigma) > 0) {
      sigma_values <- sqrt(as.numeric(x$posteriorSigma))
    } else {
      stop(
        "No sigma trace found. Provide `sigma_trace` or use an AddiVortes ",
        "object with posterior sigma samples."
      )
    }
  } else {
    if (!is.numeric(sigma_trace)) {
      stop("`sigma_trace` must be a numeric vector.")
    }
    sigma_values <- as.numeric(sigma_trace)
  }

  if (length(sigma_values) != expected_length) {
    warning("Length of sigma trace doesn't match number of posterior samples.")
    sigma_values <- rep(sigma_values[1], expected_length)
  }

  sigma_values
}

traceplotData_internal <- function(x) {
  expected_columns <- c(
    "iteration",
    "isBurnIn",
    "averageCentresPerTessellation",
    "sdCentresPerTessellation",
    "averageDimensionsPerTessellation",
    "logLikelihood"
  )

  if (!is.null(x$traceStats)) {
    missing_columns <- setdiff(expected_columns, names(x$traceStats))
    if (length(missing_columns) > 0) {
      stop("Trace statistics are missing required columns: ",
           paste(missing_columns, collapse = ", "))
    }
    return(as.data.frame(x$traceStats[expected_columns]))
  }

  if (length(x$posteriorTess) == 0) {
    stop("No posterior samples available for plotting.")
  }
  if (length(x$posteriorDim) == 0) {
    stop("No posterior dimension samples available for plotting.")
  }
  if (length(x$posteriorDim) != length(x$posteriorTess)) {
    stop("Number of posterior dimension samples does not match number of posterior tessellation samples.")
  }

  num_samples <- length(x$posteriorTess)

  centre_counts <- lapply(x$posteriorTess, function(sample) {
    if (length(sample) == 0) {
      stop("Posterior tessellation samples must contain at least one tessellation.")
    }
    vapply(sample, nrow, numeric(1))
  })

  dimension_counts <- lapply(x$posteriorDim, function(sample) {
    if (length(sample) == 0) {
      stop("Posterior dimension samples must contain at least one tessellation.")
    }
    vapply(sample, length, numeric(1))
  })

  data.frame(
    iteration = seq_len(num_samples),
    isBurnIn = rep(FALSE, num_samples),
    averageCentresPerTessellation = vapply(centre_counts, mean, numeric(1)),
    sdCentresPerTessellation = vapply(centre_counts, function(counts) {
      if (length(counts) <= 1) {
        0
      } else {
        sd(counts)
      }
    }, numeric(1)),
    averageDimensionsPerTessellation = vapply(dimension_counts, mean, numeric(1)),
    logLikelihood = rep(NA_real_, num_samples)
  )
}

plotBurnInTrace_internal <- function(trace_data, y, ylab, main, col, legend_digits = 2,
                                     show_summary = TRUE, ...) {
  burn_in <- trace_data$isBurnIn %in% TRUE
  finite_y <- y[is.finite(y)]
  if (length(finite_y) == 0) {
    plot(trace_data$iteration, rep(0, length(trace_data$iteration)),
      type = "n",
      xlab = "MCMC Iteration",
      ylab = ylab,
      main = main,
      ylim = c(0, 1),
      ...
    )
    text(mean(range(trace_data$iteration)), 0.5, "Trace unavailable", cex = 0.9)
    return(invisible(NULL))
  }

  plot(trace_data$iteration, y,
    type = "l",
    xlab = "MCMC Iteration",
    ylab = ylab,
    main = main,
    col = col, lwd = 1.5,
    ...
  )

  if (any(burn_in)) {
    lines(trace_data$iteration[burn_in], y[burn_in], col = "black", lwd = 1.5)
  }

  summary_y <- y[!burn_in & is.finite(y)]
  if (length(summary_y) == 0) {
    summary_y <- finite_y
  }

  if (show_summary && length(summary_y) > 0) {
    abline(h = mean(summary_y), col = "red", lty = 2)
  }

  legend_entries <- character(0)
  legend_cols <- character(0)
  legend_lty <- numeric(0)
  legend_lwd <- numeric(0)

  if (any(burn_in)) {
    legend_entries <- c(legend_entries, "Burn-in")
    legend_cols <- c(legend_cols, "black")
    legend_lty <- c(legend_lty, 1)
    legend_lwd <- c(legend_lwd, 1.5)
  }

  if (show_summary && length(summary_y) > 0) {
    legend_entries <- c(
      legend_entries,
      paste("Mean:", round(mean(summary_y), legend_digits))
    )
    legend_cols <- c(legend_cols, "red")
    legend_lty <- c(legend_lty, 2)
    legend_lwd <- c(legend_lwd, 1)
  }

  if (length(legend_entries) > 0) {
    legend("bottomright",
      legend = legend_entries,
      col = legend_cols,
      lty = legend_lty,
      lwd = legend_lwd,
      bty = "n", cex = 0.9
    )
  }
}

#' @title Trace Plot Method for AddiVortes
#'
#' @description
#' Displays four MCMC trace plots for a fitted `AddiVortes` object:
#' the average number of centres per tessellation, the standard deviation
#' of the number of centres per tessellation, the average number of dimensions
#' used per tessellation, and the retained-state log-likelihood component whose
#' differences form the likelihood part of the acceptance ratio.
#'
#' @param x An object of class `AddiVortes`, typically the result of a
#'   call to `AddiVortes()`.
#' @param ask Logical; if TRUE, the user is asked to press Enter before each plot.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return
#' This function is called for its side effect of creating plots and returns
#' `NULL` invisibly.
#'
#' @details
#' The four trace plots are:
#' \enumerate{
#'   \item \strong{Average Centres}: Average number of centres per tessellation.
#'   \item \strong{Centre Count Standard Deviation}: Standard deviation of the
#'     number of centres per tessellation.
#'   \item \strong{Average Dimensions}: Average number of active dimensions used
#'     per tessellation.
#'   \item \strong{Log Likelihood}: Average retained-state log-likelihood
#'     component at the end of each MCMC iteration.
#' }
#'
#' @importFrom graphics plot abline legend par text
#' @importFrom stats sd
#' @name traceplots.AddiVortes
#' @aliases traceplots
#' @usage \method{traceplots}{AddiVortes}(x, ask = FALSE, ...)
#' @export
#' @method traceplots AddiVortes
#'
#' @examples
#' \dontrun{
#' # Assuming 'fit' is a trained AddiVortes object
#' traceplots(fit)
#' }
traceplots.AddiVortes <- function(x, ask = FALSE, ...) {
  if (!inherits(x, "AddiVortes")) {
    stop("`x` must be an object of class 'AddiVortes'.")
  }

  trace_data <- traceplotData_internal(x)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(2, 2))

  if (ask) {
    cat("Press [Enter] to see average centres trace plot: ")
    readline()
  }
  plotBurnInTrace_internal(trace_data, trace_data$averageCentresPerTessellation,
    ylab = "Average Number of Tessellation Centers",
    main = "MCMC Trace: Average Centres",
    col = "purple",
    legend_digits = 1,
    ...
  )

  if (ask) {
    cat("Press [Enter] to see centre-count standard deviation trace plot: ")
    readline()
  }
  plotBurnInTrace_internal(trace_data, trace_data$sdCentresPerTessellation,
    ylab = "SD of Tessellation Centers",
    main = "MCMC Trace: Centre Count SD",
    col = "darkorange",
    ...
  )

  if (ask) {
    cat("Press [Enter] to see average dimensions trace plot: ")
    readline()
  }
  plotBurnInTrace_internal(trace_data, trace_data$averageDimensionsPerTessellation,
    ylab = "Average Number of Dimensions",
    main = "MCMC Trace: Average Dimensions",
    col = "darkblue",
    legend_digits = 1,
    ...
  )

  if (ask) {
    cat("Press [Enter] to see log-likelihood trace plot: ")
    readline()
  }
  plotBurnInTrace_internal(trace_data, trace_data$logLikelihood,
    ylab = "Log Likelihood",
    main = "MCMC Trace: Log Likelihood",
    col = "darkgreen",
    legend_digits = 4,
    ...
  )

  invisible(NULL)
}

# S3 generic for dispatch; documented only via the AddiVortes method.
#' @export
#' @rdname traceplots.AddiVortes
#' @usage NULL
traceplots <- function(x, ...) {
  UseMethod("traceplots")
}

#' @title Plot Method for AddiVortes
#'
#' @description
#' Generates comprehensive diagnostic plots for a fitted `AddiVortes` object.
#' This function creates multiple diagnostic plots including residuals,
#' MCMC traces for sigma, and tessellation complexity over iterations.
#'
#' @param x An object of class `AddiVortes`, typically the result of a
#'   call to `AddiVortes()`.
#' @param x_train A matrix or data frame of the original training covariates.
#' @param y_train The original training response. Numeric for regression;
#'   factor, character, logical or 0/1 numeric for classification.
#' @param sigma_trace An optional numeric vector of sigma values from MCMC samples.
#'   If not provided, the method will attempt to extract the posterior error
#'   standard deviation from the model object.
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
#'   \item \strong{Sigma Trace}: MCMC trace plot for the error standard deviation
#'   \item \strong{Tessellation Complexity}: Trace of average tessellation size over iterations
#'   \item \strong{Predicted vs Observed}: Scatter plot with credible intervals
#' }
#'
#' @importFrom graphics plot abline title par layout lines points segments legend
#' @importFrom stats lowess residuals fitted predict sd
#' @export
#' @method plot AddiVortes
#'
#' @examples
#' \dontrun{
#' # Assuming 'fit' is a trained AddiVortes object
#' plot(fit, x_train = x_train_data, y_train = y_train_data)
#'
#' # Show only specific plots
#' plot(fit, x_train = x_train_data, y_train = y_train_data, which = c(1, 3))
#'
#' # With custom sigma trace
#' plot(fit,
#'   x_train = x_train_data, y_train = y_train_data,
#'   sigma_trace = my_sigma_samples
#' )
#' }
plot.AddiVortes <- function(x, x_train, y_train, sigma_trace = NULL,
                            which = c(1, 2, 3), ask = FALSE, ...) {
  # --- Input Validation ---
  if (!inherits(x, "AddiVortes")) {
    stop("`x` must be an object of class 'AddiVortes'.")
  }
  if (missing(x_train) || missing(y_train)) {
    stop("`x_train` and `y_train` must be provided for diagnostic plots.")
  }
  if (!is.matrix(x_train) && !is.data.frame(x_train)) {
    stop("`x_train` must be a matrix or data frame.")
  }
  is_class <- isClassification_internal(x)
  if (!is_class && !is.numeric(y_train)) {
    stop("`y_train` must be a numeric vector.")
  }
  if (NROW(x_train) != length(y_train)) {
    stop("The number of rows in `x_train` must match the length of `y_train`.")
  }
  if (length(x$posteriorTess) == 0) {
    stop("No posterior samples available for plotting.")
  }

  # Validate which parameter
  which <- intersect(which, 1:4)
  if (is_class) {
    which <- setdiff(which, 2L)
  }
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
  y_pred_prob <- NULL
  if (any(which %in% c(1, 4))) {
    y_pred_mean <- predict(x, newdata = x_train, type = "response",
                           showProgress = FALSE)
    if (isTRUE(x$task == "binary")) {
      y01 <- if (is.numeric(y_train)) {
        as.numeric(y_train == 1 | y_train == x$classLevels[2])
      } else {
        as.numeric(as.character(y_train) == x$classLevels[2] |
                     y_train == x$classLevels[2])
      }
      if (is.logical(y_train)) {
        y01 <- as.numeric(y_train)
      }
      y_pred_prob <- y_pred_mean
      residuals <- y01 - y_pred_mean
    } else if (isTRUE(x$task == "multinomial")) {
      y_pred_prob <- y_pred_mean
      true_idx <- match(as.character(y_train), x$classLevels)
      if (anyNA(true_idx)) {
        true_idx <- as.integer(factor(y_train, levels = x$classLevels))
      }
      residuals <- y_pred_mean[cbind(seq_len(nrow(y_pred_mean)), true_idx)]
    } else {
      residuals <- y_train - y_pred_mean
    }
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

    if (isTRUE(x$task == "multinomial")) {
      plot(seq_along(residuals), residuals,
        xlab = "Observation",
        ylab = "Predicted Probability of True Class",
        main = "Probability of Observed Class",
        pch = 19, col = "darkblue", cex = 0.8,
        ylim = c(0, 1),
        ...
      )
      abline(h = 1 / length(x$classLevels), col = "red", lty = 2, lwd = 2)
    } else {
      plot(y_pred_mean, residuals,
        xlab = if (is_class) "Fitted Probability" else "Fitted Values",
        ylab = if (is_class) "Observed minus Probability" else "Residuals",
        main = "Residuals vs Fitted",
        pch = 19, col = "darkblue", cex = 0.8,
        ...
      )

      # Add horizontal line at y = 0
      abline(h = 0, col = "red", lty = 2, lwd = 2)

      # Add smoothed trend line
      if (length(y_pred_mean) > 3) {
        smooth_line <- lowess(y_pred_mean, residuals)
        lines(smooth_line, col = "orange", lwd = 2)
      }
    }

    # Add RMSE or accuracy annotation
    if (is_class) {
      acc_text <- paste("Accuracy =", round(x$inSampleAccuracy, 4))
      legend("topright", legend = acc_text, bty = "n", cex = 0.9)
    } else {
      rmse_text <- paste("RMSE =", round(x$inSampleRmse, 4))
      legend("topright", legend = rmse_text, bty = "n", cex = 0.9)
    }
  }

  # --- Plot 2: Sigma Trace ---
  if (2 %in% which) {
    if (ask && n_plots > 1) {
      cat("Press [Enter] to see sigma trace plot: ")
      readline()
    }

    sigma_values <- extractErrorStandardDeviationTrace_internal(x, sigma_trace)

    plot(seq_along(sigma_values), sigma_values,
      type = "l",
      xlab = "MCMC Iteration",
      ylab = expression(sigma),
      main = "MCMC Trace: Error Standard Deviation",
      col = "darkgreen", lwd = 1.5,
      ...
    )

    # Add horizontal line at mean
    abline(h = mean(sigma_values), col = "red", lty = 2)

    # Add convergence statistics
    sigma_mean <- round(mean(sigma_values), 4)
    sigma_sd <- round(sd(sigma_values), 4)
    legend("topright",
      legend = c(
        paste("Mean:", sigma_mean),
        paste("SD:", sigma_sd)
      ),
      bty = "n", cex = 0.9
    )
  }

  # --- Plot 3: Tessellation Complexity Trace ---
  if (3 %in% which) {
    if (ask && n_plots > 1) {
      cat("Press [Enter] to see tessellation complexity trace: ")
      readline()
    }

    plot(seq_along(tess_complexity), tess_complexity,
      type = "l",
      xlab = "MCMC Iteration",
      ylab = "Average Number of Tessellation Centers",
      main = "MCMC Trace: Tessellation Complexity",
      col = "purple", lwd = 1.5,
      ...
    )

    # Add horizontal line at mean
    abline(h = mean(tess_complexity), col = "red", lty = 2)

    # Add summary statistics
    complexity_mean <- round(mean(tess_complexity), 1)
    complexity_range <- round(range(tess_complexity), 1)
    legend("topright",
      legend = c(
        paste("Mean:", complexity_mean),
        paste(
          "Range: [", complexity_range[1], ",",
          complexity_range[2], "]"
        )
      ),
      bty = "n", cex = 0.9
    )
  }

  # --- Plot 4: Predicted vs Observed ---
  if (4 %in% which) {
    if (ask && n_plots > 1) {
      cat("Press [Enter] to see predicted vs observed plot: ")
      readline()
    }

    if (isTRUE(x$task == "binary")) {
      y01 <- if (is.logical(y_train)) {
        as.numeric(y_train)
      } else if (is.numeric(y_train) && all(y_train %in% c(0, 1))) {
        as.numeric(y_train)
      } else {
        as.numeric(as.character(y_train) == x$classLevels[2])
      }
      plot(y01, y_pred_prob,
        xlab = "Observed Class (0/1)",
        ylab = "Predicted Probability",
        main = "Predicted Probability vs Observed",
        pch = 19, col = "darkblue", cex = 0.8,
        ylim = c(0, 1),
        ...
      )
      abline(h = 0.5, col = "red", lty = 2, lwd = 2)
    } else if (isTRUE(x$task == "multinomial")) {
      true_idx <- match(as.character(y_train), x$classLevels)
      p_true <- y_pred_prob[cbind(seq_len(nrow(y_pred_prob)), true_idx)]
      plot(seq_along(p_true), p_true,
        xlab = "Observation",
        ylab = "Predicted Probability of True Class",
        main = "Probability of Observed Class",
        pch = 19, col = "darkblue", cex = 0.8,
        ylim = c(0, 1),
        ...
      )
      abline(h = 1 / length(x$classLevels), col = "red", lty = 2, lwd = 2)
    } else {
      # Get quantile predictions for uncertainty
      y_pred_quantiles <- predict(x,
        newdata = x_train, type = "quantile",
        quantiles = c(0.025, 0.975),
        showProgress = FALSE
      )

      # Create the scatter plot
      plot(y_train, y_pred_mean,
        xlab = "Observed Values",
        ylab = "Predicted Values",
        main = "Predicted vs Observed",
        pch = 19, col = "darkblue", cex = 0.8,
        xlim = range(c(y_train, y_pred_mean)),
        ylim = range(c(y_train, y_pred_mean)),
        ...
      )

      # Add the line of equality (perfect prediction)
      abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)

      # Add uncertainty intervals
      for (i in seq_along(y_train)) {
        segments(y_train[i], y_pred_quantiles[i, 1],
          y_train[i], y_pred_quantiles[i, 2],
          col = "lightblue", lwd = 1
        )
      }

      # Calculate and display R-squared
      ss_res <- sum(residuals^2)
      ss_tot <- sum((y_train - mean(y_train))^2)
      r_squared <- 1 - (ss_res / ss_tot)

      legend("topleft",
        legend = c(
          "Perfect Prediction",
          "95% Prediction Intervals"
        ),
        col = c("red", "lightblue"),
        lty = c(2, 1),
        lwd = c(2, 1),
        pch = c(NA, NA), cex = 0.9, bty = "n"
      )
      legend("bottomright",
        legend = c(paste("R^2 =", round(r_squared, 3))),
        col = c("black"),
        lty = c(NA),
        lwd = c(NA),
        pch = c(NA), cex = 0.9, bty = "n"
      )
    }
  }

  # Return invisibly
  invisible(NULL)
}
