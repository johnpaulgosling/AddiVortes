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
#' @param totalMCMCIter The number of iterations.
#' @param mcmcBurnIn The number of burn in iterations.
#' @param nu The degrees of freedom.
#' @param q The quantile.
#' @param k The number of centers.
#' @param sd The standard deviation.
#' @param omega The prior probability of adding a dimension.
#' @param lambdaRate The rate of the Poisson distribution for the number of centers.
#' @param yTest A vector of the output values for the test set.
#' @param xTest A matrix of the covariates for the test set.
#' @param IntialSigma The method used to calculate the initial variance.
#' @param thinning The thinning rate.
#'
#' @return A list containing the following elements:
#' - `AddiVortesModel`: A list containing the posterior samples of the tessellations,
#'   dimensions and predictions.
#' - `inSampleRmse`: The RMSE value for the training samples.
#'
#' @export
AddiVortes <- function(y, x, m = 200, totalMCMCIter = 1200,
                       mcmcBurnIn = 200, nu = 6, q = 0.85,
                       k = 3, sd = 0.8, omega = 3, lambdaRate = 25,
                       yTest, xTest, IntialSigma = "Linear",
                       thinning = 1) {
  #### Scaling x and y ---------------------------------------------------------
  yScalingResult <- scaleData_internal(y)
  yScaled <- yScalingResult$scaled_data # Vector of values
  yCenter <- yScalingResult$centers
  yRange <- yScalingResult$ranges

  xScalingResult <- scaleData_internal(x)
  xScaled <- xScalingResult$scaled_data # Matrix of values
  xCenters <- xScalingResult$centers # Vector of values
  xRanges <- xScalingResult$ranges   # Vector of values

  #### Initialise predictions --------------------------------------------------
  # Initialise:
  # Prediction Set (A list of vectors with the output values for each tessellation),
  # Dimension set (A list of vectors with the covariates included in the tessellations);
  # and Tessellation Set (A list of matrices that give the
  #                      coordinates of the centres in the tessellations)
  pred <- rep(list(matrix(mean(yScaled) / m)), m)
  dim <- sapply(1:m, function(ignoredIndex) {
    list(sample(1:length(x[1, ]), 1))
  })
  tess <- sapply(1:m, function(ignoredIndex) {
    list(matrix(rnorm(1, 0, sd)))
  })

  #### Set-up MCMC -------------------------------------------------------------
  # Prepare some variables used in the backfitting algorithm.
  # We start off with the mean of the scaled y values as the prediction for all
  # tessellations.
  sumOfAllTess <- rep(mean(yScaled),
                      length(yScaled))
  # The variance that captures variability around the mean of the scaled y values.
  sigmaSquaredMu <- (0.5 / (k * sqrt(m)))^2
  lastTessPred <- matrix

  # Matrices that will hold the samples from the posterior distribution
  # for the training samples and test samples.
  posteriorSamples <- floor((totalMCMCIter - mcmcBurnIn) / thinning)
  predictionMatrix <- array(dim = c(length(y),
                                    posteriorSamples))
  testMatrix <- array(dim = c(length(yTest),
                              posteriorSamples))

  # Finding lambda
  if (IntialSigma == "Naive") {
    # Usually used if p is greater then n. Uses Standard deviation of y to predict sigma.
    sigmaSquaredHat <- var(yScaled)
  } else {
    # Default method using residual standard deviation from a least-squared linear
    # regression of y, to predict sigma.
    multiLinear <- lm(yScaled ~ xScaled)
    sigmaSquaredHat <- sum(multiLinear$residuals^2) /
      (length(yScaled) - length(xScaled[1, ]) - 1)
  }
  lambda <- optim(
    par = 1,
    fittingFunction,
    method = "Brent",
    lower = 0.001,
    upper = 100,
    q = q, nu = nu,
    sigmaSquaredHat = sigmaSquaredHat
  )$par

  # Determine number of samples to store
  numPosteriorSamplesToStore <- 0
  if (totalMCMCIter > mcmcBurnIn) {
    numPosteriorSamplesToStore <- floor((totalMCMCIter - mcmcBurnIn) / thinning)
  }
  if (numPosteriorSamplesToStore < 0) numPosteriorSamplesToStore <- 0

  # Lists to store the states of tess, dim, pred for the model object output
  outputPosteriorTess <- vector("list",
                                numPosteriorSamplesToStore)
  outputPosteriorDim <- vector("list",
                               numPosteriorSamplesToStore)
  outputPosteriorPred <- vector("list",
                                numPosteriorSamplesToStore)

  currentStorageIdx <- 1 # Index for the new output lists

  # Setting up progress bar
  pbar <- utils::txtProgressBar(min = 0, max = totalMCMCIter,
                                style = 3, width = 50, char = "=")

  #### MCMC Loop ---------------------------------------------------------------
  for (i in 1:totalMCMCIter) {
    # Sample sigma squared using all tessellations to predict the outcome variables
    sigmaSquared <- sampleSigmaSquared(yScaled,
                                       nu,
                                       lambda,
                                       sumOfAllTess)

    for (j in 1:m) {
      # Propose new Tessellation
      newTessOutput <- proposeTessellation(xScaled,
                                           j,
                                           tess,
                                           dim,
                                           sd)
      tessStar <- newTessOutput[[1]]
      dimStar <- newTessOutput[[2]]
      modification <- newTessOutput[[3]]

      # Calculate the n-vector of partial residuals derived from a fitting process
      # that excludes the jth tessellation and the number of observations in each cell.
      residualsOutput <- calculateResiduals(yScaled,
                                            xScaled,
                                            j,
                                            sumOfAllTess,
                                            tess,
                                            dim,
                                            pred,
                                            tessStar,
                                            dimStar,
                                            lastTessPred)
      # Old and New refer to the original and proposed tessellations
      rIjOld <- residualsOutput[[1]]
      nIjOld <- residualsOutput[[2]]
      rIjNew <- residualsOutput[[3]]
      nIjNew <- residualsOutput[[4]]
      # Keeps track of the prediction for all tessellations to help
      # sample sigma squared.
      sumOfAllTess <- residualsOutput[[5]]
      # Gives the row of each observation for the cell it falls in for the
      # proposed tessellation.
      indexesStar <- residualsOutput[[6]]
      # Gives the row of each observation for the cell it falls in for the
      # original tessellation.
      indexes <- residualsOutput[[7]]

      if (!any(nIjNew == 0)) {
        # Automatically reject proposed tessellation if there exists a cell
        # with no observations in.
        logAcceptanceProb <- acceptanceProbability(xScaled,
                                                   tessStar,
                                                   dimStar,
                                                   j,
                                                   rIjOld, nIjOld,
                                                   rIjNew, nIjNew,
                                                   sigmaSquared,
                                                   modification,
                                                   sigmaSquaredMu,
                                                   omega,
                                                   lambdaRate)

        if (log(runif(n = 1)) < logAcceptanceProb) {
          # Accepts the proposed tessellation is accepted then calculates the new
          # output values for the new tessellation.
          tess <- tessStar
          dim <- dimStar
          pred[[j]] <- sampleMuValues(j, tessStar,
                                      rIjNew, nIjNew,
                                      sigmaSquaredMu,
                                      sigmaSquared)
          lastTessPred <- pred[[j]][indexesStar]
        } else {
          # Rejects the proposed tessellation then calculates new output values
          # for the original tessellation.
          pred[[j]] <- sampleMuValues(j, tess, rIjOld, nIjOld,
                                      sigmaSquaredMu, sigmaSquared)
          lastTessPred <- pred[[j]][indexes]
        }
      } else {
        # Rejects the proposed tessellation then calculates new output values for
        # the original tessellation.
        pred[[j]] <- sampleMuValues(j, tess, rIjOld, nIjOld,
                                    sigmaSquaredMu, sigmaSquared)
        lastTessPred <- pred[[j]][indexes]
      }
      if (j == m) {
        # If j equals m then adds the last tessellation output values
        # to give a prediction.
        sumOfAllTess <- sumOfAllTess + lastTessPred
      }
    }

    # Update the progress bar
    utils::setTxtProgressBar(pbar, i)

    if (i > mcmcBurnIn & (i - mcmcBurnIn) %% thinning == 0) {
      # vectors that hold the predictions for each iteration after burn in.
      predictionMatrix[, (i - mcmcBurnIn) / thinning] <- sumOfAllTess
    }

    # Store the posterior samples
    if (numPosteriorSamplesToStore > 0 &&
        i > mcmcBurnIn &
        (i - mcmcBurnIn) %% thinning == 0) {
      # Store the current state of tess, dim, pred
      outputPosteriorTess[[currentStorageIdx]] <- tess
      outputPosteriorDim[[currentStorageIdx]] <- dim
      outputPosteriorPred[[currentStorageIdx]] <- pred
      currentStorageIdx <- currentStorageIdx + 1
    }
  } # End of MCMC Loop

  # Close the progress bar
  close(pbar)

  # Finding the mean of the prediction over the iterations and then unscaling
  # the predictions.
  meanYhat <- (rowSums(predictionMatrix) / (posteriorSamples)) * yRange +
    yCenter

  return( # Returns the RMSE value for the test samples.
    list(AddiVortesModel = list(outputPosteriorTess,
                                outputPosteriorDim,
                                outputPosteriorPred,
                                currentStorageIdx,
                                xCenters,
                                xRanges,
                                yCenter,
                                yRange),
         inSampleRmse = sqrt(mean((y - meanYhat)^2))
    )
  )
}
