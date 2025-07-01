#' @title AddiVortes
#'
#' @description
#' The AddiVortes function is a Bayesian nonparametric regression model that uses a
#' tessellation to model the relationship between the covariates and the output values.
#' The model uses a backfitting algorithm to sample from the posterior distribution of
#' the output values for each tessellation. The function returns the RMSE value for
#' the test samples.
#'
#' @param y A vector of the output values.
#' @param x A matrix of the covariates.
#' @param m The number of tessellations.
#' @param max_iter The number of iterations.
#' @param burn_in The number of burn in iterations.
#' @param nu The degrees of freedom.
#' @param q The quantile.
#' @param k The number of centers.
#' @param sd The standard deviation.
#' @param Omega The prior probability of adding a dimension.
#' @param lambda_rate The rate of the Poisson distribution for the number of centers.
#' @param yTest A vector of the output values for the test set.
#' @param xTest A matrix of the covariates for the test set.
#' @param IntialSigma The method used to calculate the initial variance.
#' @param thinning The thinning rate.
#'
#' @return A list containing the following elements:
#' - `AddiVortes_model`: A list containing the posterior samples of the tessellations,
#'  dimensions, and predictions.
#' - `In_sample_RMSE`: The RMSE value for the training samples.
#'
#' @export
AddiVortes <- function(y, x, m = 200, max_iter = 1200,
                       burn_in = 200, nu = 6, q = 0.85,
                       k = 3, sd = 0.8, Omega = 3, lambda_rate = 25,
                       yTest, xTest, IntialSigma = "Linear",
                       thinning = 1) {
  #### Scaling x and y ---------------------------------------------------------
  y_scaling_result <- scale_data_internal(y)
  yScaled <- y_scaling_result$scaled_data # Vector of values
  y_center <- y_scaling_result$centers
  y_range <- y_scaling_result$ranges

  x_scaling_result <- scale_data_internal(x)
  xScaled <- x_scaling_result$scaled_data # Matrix of values
  x_centers <- x_scaling_result$centers # Vector of values
  x_ranges <- x_scaling_result$ranges   # Vector of values

  #### Initialise predictions --------------------------------------------------
  # Initialise:
  # Prediction Set (A list of vectors with the output values for each tessellation),
  # Dimension set (A list of vectors with the covariates included in the tessellations);
  # and Tessellation Set (A list of matrices that give the
  #                       coordinates of the centres in the tessellations)
  Pred <- rep(list(matrix(mean(yScaled) / m)), m)
  Dim <- sapply(1:m, function(ignored_index) {
    list(sample(1:length(x[1, ]), 1))
  })
  Tess <- sapply(1:m, function(ignored_index) {
    list(matrix(rnorm(1, 0, sd)))
  })

  #### Set-up MCMC -------------------------------------------------------------
  # Prepare some variables used in the backfitting algorithm
  SumOfAllTess <- rep(mean(yScaled),
                      length(yScaled))
  SigmaSquaredMu <- (0.5 / (k * sqrt(m)))^2
  LastTessPred <- matrix

  # Matrices that will hold the samples from the posterior distribution
  # for the training samples and test samples.
  posterior_samples <- floor((max_iter - burn_in) / thinning) + 1
  PredictionMatrix <- array(dim = c(length(y),
                                    posterior_samples))
  TestMatrix <- array(dim = c(length(yTest),
                              posterior_samples))

  # Finding lambda
  if (IntialSigma == "Naive") {
    # Usually used if p is greater then n. Uses Standard deviation of y to predict sigma.
    SigmaSquaredHat <- var(yScaled)
  } else {
    # Default method using residual standard deviation from a least-squared linear
    # regression of y, to predict sigma.
    MultiLinear <- lm(yScaled ~ xScaled)
    SigmaSquaredHat <- sum(MultiLinear$residuals^2) /
      (length(yScaled) - length(xScaled[1, ]) - 1)
  }
  lambda <- optim(
    par = 1,
    Fitting_Function,
    method = "Brent",
    lower = 0.001,
    upper = 100,
    q = q, nu = nu,
    sigmaSquared_hat = SigmaSquaredHat
  )$par

  # Determine number of samples to store
  num_posterior_samples_to_store <- 0
  if (max_iter > burn_in) {
    num_posterior_samples_to_store <- floor((max_iter - burn_in) / thinning)
  }
  if (num_posterior_samples_to_store < 0) num_posterior_samples_to_store <- 0

  # Lists to store the states of Tess, Dim, Pred for the model object output
  output_posterior_Tess <- vector("list", num_posterior_samples_to_store)
  output_posterior_Dim <- vector("list", num_posterior_samples_to_store)
  output_posterior_Pred <- vector("list", num_posterior_samples_to_store)

  current_storage_idx <- 1 # Index for the new output lists

  # Setting up progress bar
  pbar <- utils::txtProgressBar(min = 0, max = max_iter,
                                style = 3, width = 50, char = "=")

  #### MCMC Loop ---------------------------------------------------------------
  for (i in 1:max_iter) {
    # Sample Sigma squared using all tessellations to predict the outcome variables
    SigmaSquared <- Sample_Sigma_Squared(yScaled,
                                         nu,
                                         lambda,
                                         SumOfAllTess)

    for (j in 1:m) {
      # Propose new Tessellation
      NewTessOutput <- Propose_Tessellation(xScaled,
                                            j,
                                            Tess,
                                            Dim,
                                            sd)
      TessStar <- NewTessOutput[[1]]
      DimStar <- NewTessOutput[[2]]
      Modification <- NewTessOutput[[3]]

      # Calculate the n-vector of partial residuals derived from a fitting process
      # that excludes the jth tessellation and the number of observations in each cell.
      ResidualsOutput <- Calculate_Residuals(yScaled,
                                             xScaled,
                                             j,
                                             SumOfAllTess,
                                             Tess,
                                             Dim,
                                             Pred,
                                             TessStar,
                                             DimStar,
                                             LastTessPred)
      # Old and New refer to the original and proposed tessellations
      R_ijOld <- ResidualsOutput[[1]]
      n_ijOld <- ResidualsOutput[[2]]
      R_ijNew <- ResidualsOutput[[3]]
      n_ijNew <- ResidualsOutput[[4]]
      # Keeps track of the prediction for all tessellations to help
      # sample sigma squared.
      SumOfAllTess <- ResidualsOutput[[5]]
      # Gives the row of each observation for the cell it falls in for the
      # proposed tessellation.
      IndexesStar <- ResidualsOutput[[6]]
      # Gives the row of each observation for the cell it falls in for the
      # original tessellation.
      Indexes <- ResidualsOutput[[7]]

      if (!any(n_ijNew == 0)) {
        # Automatically reject proposed tessellation if there exists a cell
        # with no observations in.
        LogAcceptanceProb <- Acceptance_Probability(xScaled,
                                                    TessStar,
                                                    DimStar,
                                                    j,
                                                    R_ijOld, n_ijOld,
                                                    R_ijNew, n_ijNew,
                                                    SigmaSquared,
                                                    Modification,
                                                    SigmaSquaredMu,
                                                    Omega,
                                                    lambda_rate)

        if (log(runif(n = 1)) < LogAcceptanceProb) {
          # Accepts the proposed tessellation is accepted then calculates the new
          # output values for the new tessellation.
          Tess <- TessStar
          Dim <- DimStar
          Pred[[j]] <- Sample_mu_values(j, TessStar,
                                        R_ijNew, n_ijNew,
                                        SigmaSquaredMu,
                                        SigmaSquared)
          LastTessPred <- Pred[[j]][IndexesStar]
        } else {
          # Rejects the proposed tessellation then calculates new output values
          # for the original tessellation.
          Pred[[j]] <- Sample_mu_values(j, Tess, R_ijOld, n_ijOld,
                                        SigmaSquaredMu, SigmaSquared)
          LastTessPred <- Pred[[j]][Indexes]
        }
      } else {
        # Rejects the proposed tessellation then calculates new output values for
        # the original tessellation.
        Pred[[j]] <- Sample_mu_values(j, Tess, R_ijOld, n_ijOld,
                                      SigmaSquaredMu, SigmaSquared)
        LastTessPred <- Pred[[j]][Indexes]
      }
      if (j == m) {
        # If j equals m then adds the last tessellation output values
        # to give a prediction.
        SumOfAllTess <- SumOfAllTess + LastTessPred
      }
    }

    # Update the progress bar
    utils::setTxtProgressBar(pbar, i)

    if (i >= burn_in & (i - burn_in) %% thinning == 0) {
      # vectors that hold the predictions for each iteration after burn in.
      PredictionMatrix[, 1 + (i - burn_in) / thinning] <- SumOfAllTess
    }

    # Store the posterior samples
    if (num_posterior_samples_to_store > 0 && i >= burn_in & (i - burn_in) %% thinning == 0) {
      # Store the current state of Tess, Dim, Pred
      output_posterior_Tess[[current_storage_idx]] <- Tess
      output_posterior_Dim[[current_storage_idx]] <- Dim
      output_posterior_Pred[[current_storage_idx]] <- Pred
      current_storage_idx <- current_storage_idx + 1
      #  }
    }
  } # End of MCMC Loop

  # Close the progress bar
  close(pbar)

  # Finding the mean of the prediction over the iterations and then unscaling
  # the predictions.
  mean_yhat <- (rowSums(PredictionMatrix) / (posterior_samples)) * y_range +
    y_center

  return( # Returns the RMSE value for the test samples.
    list(AddiVortes_model = list(output_posterior_Tess,
                                 output_posterior_Dim,
                                 output_posterior_Pred,
                                 current_storage_idx,
                                 x_centers,
                                 x_ranges,
                                 y_center,
                                 y_range),
         In_sample_RMSE = sqrt(mean((y - mean_yhat)^2))
    )
  )
}
