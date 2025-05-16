#' @title PredictAddiVortes
#'
#' @description
#' Utilises a fitted AddiVortes model to make predictions on new test data and
#' calculates the Root Mean Squared Error (RMSE). This function aims to
#' reproduce the Out-of-Sample RMSE calculation performed internally by the
#' AddiVortes fitting function.
#'
#' @param addivortes_model_fit The 'AddiVortes_model' list object returned by the
#'                             AddiVortes function. This contains posterior samples
#'                             and scaling parameters.
#' @param x_new A matrix of the covariates for the new test set. The number of
#'              columns must match the original training data.
#' @param y_new A vector of the true output values for the new test set. Must have
#'              the same number of observations as rows in `x_new`.
#'
#' @return The RMSE value for the new test samples. Returns NA if the model
#'         contains no posterior samples.
#'
#' @details
#' This function requires the helper functions `apply_scaling_internal` and
#' `TestSet_Prediction` (used by the original `AddiVortes` function) to be
#' available in the R environment.
#'
#' @export
PredictAddiVortes <- function(addivortes_model_fit,
                              x_new,
                              y_new) {

  # Extract components from the fitted model object
  posterior_Tess_samples <- addivortes_model_fit[[1]]
  posterior_Dim_samples  <- addivortes_model_fit[[2]]
  posterior_Pred_samples <- addivortes_model_fit[[3]]
  # current_storage_idx indicates (number of stored samples + 1)
  num_stored_samples     <- addivortes_model_fit[[4]] - 1

  x_model_centers <- addivortes_model_fit[[5]]
  x_model_ranges  <- addivortes_model_fit[[6]]
  y_model_center  <- addivortes_model_fit[[7]]
  y_model_range   <- addivortes_model_fit[[8]]

  # Basic input validation
  if (!is.matrix(x_new)) {
    stop("'x_new' must be a matrix.")
  }
  if (ncol(x_new) != length(x_model_centers)) {
    stop("Number of columns in 'x_new' does not match the original training data.")
  }
  if (length(y_new) != nrow(x_new)) {
    stop("Length of 'y_new' must match the number of rows in 'x_new'.")
  }

  # Handle cases with no stored samples
  if (num_stored_samples <= 0) {
    warning("The AddiVortes model contains no stored posterior samples. Cannot make predictions.")
    return(NA_real_)
  }

  # Scale the new test covariates using stored scaling parameters
  # Assumes apply_scaling_internal(mat, centers, ranges) is available
  x_new_scaled <- apply_scaling_internal(mat = x_new,
                                         centers = x_model_centers,
                                         ranges = x_model_ranges)

  # Determine the number of tessellations (m) from the first stored sample
  # Assumes m is constant across all samples
  m_tessellations <- length(posterior_Tess_samples[[1]])

  # Initialize a matrix to store predictions for each posterior sample
  # Dimensions: [number of new test observations, number of stored posterior samples]
  NewTestDataPredictionsMatrix <- array(dim = c(nrow(x_new_scaled),
                                                num_stored_samples))

  # Loop through each stored posterior sample to make predictions
  for (s_idx in 1:num_stored_samples) {
    current_Tess <- posterior_Tess_samples[[s_idx]]
    current_Dim  <- posterior_Dim_samples[[s_idx]]
    current_Pred <- posterior_Pred_samples[[s_idx]]

    # Get predictions for the current posterior sample
    # Assumes TestSet_Prediction(xTest_scaled, m, Tess, Dim, Pred) is available
    predictions_for_sample_s <- TestSet_Prediction(x_new_scaled,
                                                   m_tessellations,
                                                   current_Tess,
                                                   current_Dim,
                                                   current_Pred)
    NewTestDataPredictionsMatrix[, s_idx] <- predictions_for_sample_s
  }

  # Calculate the mean of predictions across all posterior samples (still scaled)
  mean_yhat_new_scaled <- rowSums(NewTestDataPredictionsMatrix) / num_stored_samples

  # Unscale the mean predictions
  mean_yhat_new_unscaled <- mean_yhat_new_scaled * y_model_range + y_model_center

  # Calculate the Out-of-Sample RMSE for the new test data
  out_of_sample_RMSE <- sqrt(mean((y_new - mean_yhat_new_unscaled)^2))

  return(list(out_of_sample_RMSE,
              mean_yhat_new_unscaled))
}
