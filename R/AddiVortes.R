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
#' @param k The number of centres.
#' @param sd The standard deviation.
#' @param omega omega/(number of covariates) is the prior probability of adding a dimension.
#' @param lambdaRate The rate of the Poisson distribution for the number of centres.
#' @param IntialSigma The method used to calculate the initial variance.
#' @param thinning The thinning rate.
#' @param showProgress Logical; if TRUE (default), progress bars and messages are shown during fitting.
#'
#' @return An AddiVortesFit object containing the posterior samples of the
#' tessellations, dimensions and predictions.
#'
#' @examples
#' \donttest{
#' # Simple example with simulated data
#' set.seed(123)
#' x <- matrix(rnorm(50), 10, 5)
#' y <- rnorm(10)
#' # Fit model with reduced iterations for quick example
#' fit <- AddiVortes(y, x, m = 5, totalMCMCIter = 50, mcmcBurnIn = 10)
#' }
#'
#' @export
AddiVortes <- function(y, x, m = 200, totalMCMCIter = 1200,
                       mcmcBurnIn = 200, nu = 6, q = 0.85,
                       k = 3, sd = 0.8, 
                       omega = min(3, ncol(x)), 
                       lambdaRate = 25,
                       IntialSigma = "Linear",
                       thinning = 1, showProgress = TRUE) {
  #### Scaling x and y ---------------------------------------------------------
  yScalingResult <- scaleData_internal(y)
  yScaled <- yScalingResult$scaledData # Vector of values
  yCentre <- yScalingResult$centres
  yRange <- yScalingResult$ranges
  
  xScalingResult <- scaleData_internal(x)
  xScaled <- xScalingResult$scaledData # Matrix of values
  xCentres <- xScalingResult$centres # Vector of values
  xRanges <- xScalingResult$ranges # Vector of values
  
  #### Check dimensions --------------------------------------------------------
  n <- length(y)
  p <- ncol(xScaled)
  if (p > n) {
    warning(
      "Number of covariates (p = ", p, ") exceeds number of observations (n = ", n, "). ",
      "Model results may not be stable. Consider reducing the number of covariates or ",
      "increasing the sample size.",
      call. = FALSE
    )
  }
  
  if (omega > p) {
    message(
      "Note: omega (", omega, ") exceeds number of covariates (", p, "). ",
      "The dimension inclusion probability will be clamped to 100%."
    )
  }
  
  #### Initialise predictions --------------------------------------------------
  # Initialise:
  # Prediction Set (A list of vectors with the output values for each tessellation),
  # Dimension set (A list of vectors with the covariates included in the tessellations);
  # and Tessellation Set (A list of matrices that give the
  #                       coordinates of the centres in the tessellations)
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
  sumOfAllTess <- rep(
    mean(yScaled),
    length(yScaled)
  )
  # The variance that captures variability around the mean of the scaled y values.
  sigmaSquaredMu <- (0.5 / (k * sqrt(m)))^2
  lastTessPred <- matrix
  
  # Matrices that will hold the samples from the posterior distribution
  # for the training samples and test samples.
  posteriorSamples <- floor((totalMCMCIter - mcmcBurnIn) / thinning)
  predictionMatrix <- array(dim = c(
    length(y),
    posteriorSamples
  ))
  
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
  outputPosteriorTess <- vector(
    "list",
    numPosteriorSamplesToStore
  )
  outputPosteriorDim <- vector(
    "list",
    numPosteriorSamplesToStore
  )
  outputPosteriorPred <- vector(
    "list",
    numPosteriorSamplesToStore
  )
  outputPosteriorSigma <- numeric(numPosteriorSamplesToStore)
  sigmaSquared <- NULL
  
  currentStorageIdx <- 1 # Index for the new output lists
  
  # Some precalculations
  numCovariates <- ncol(xScaled)
  covariateIndices <- seq_len(numCovariates)
  current_indices <- vector("list", m)
  for(k in 1:m) {
    current_indices[[k]] <- cellIndices(xScaled, tess[[k]], dim[[k]])
  }
  
  # Initial message and progress bar setup
  if (showProgress) {
    cat("Fitting AddiVortes model to input data...\n")
    cat("Input dimensions: ", nrow(xScaled),
        " observations, ", ncol(xScaled),
        " covariates\n", sep = "")
    cat("Model configuration: ", m,
        " tessellations, ", totalMCMCIter,
        " total iterations (", mcmcBurnIn,
        " burn-in)\n\n", sep = "")
  }
  
  #### MCMC Loop ---------------------------------------------------------------
  
  # Initialize progress tracking variables
  pbar_burnin <- NULL
  pbar_sampling <- NULL
  
  # Start burn-in phase
  if (showProgress && mcmcBurnIn > 0) {
    cat("Phase 1: Burn-in sampling (", mcmcBurnIn, " iterations)\n", sep = "")
    pbar_burnin <- txtProgressBar(min = 0, max = mcmcBurnIn,
                                 style = 3, width = 50, char = "=")
  }
  
  for (i in 1:totalMCMCIter) {
    # Progress bar management
    if (showProgress) {
      if (i <= mcmcBurnIn && !is.null(pbar_burnin)) {
        setTxtProgressBar(pbar_burnin, i)
      } else if (i == mcmcBurnIn + 1 && mcmcBurnIn > 0) {
        # Close burn-in progress bar and start sampling phase
        if (!is.null(pbar_burnin)) {
          close(pbar_burnin)
          cat("\n")
        }
        if (totalMCMCIter > mcmcBurnIn) {
          cat("Phase 2: Posterior sampling (",
              totalMCMCIter - mcmcBurnIn,
              " iterations)\n",
              sep = "")
          pbar_sampling <- txtProgressBar(min = 0, max = totalMCMCIter - mcmcBurnIn,
                                         style = 3, width = 50, char = "=")
        }
      } else if (i > mcmcBurnIn && !is.null(pbar_sampling)) {
        setTxtProgressBar(pbar_sampling, i - mcmcBurnIn)
      } else if (mcmcBurnIn == 0 && i == 1 && totalMCMCIter > 0) {
        # No burn-in phase, start directly with sampling
        cat("Posterior sampling (", totalMCMCIter, " iterations)\n")
        pbar_sampling <- txtProgressBar(min = 0, max = totalMCMCIter,
                                       style = 3, width = 50, char = "=")
        setTxtProgressBar(pbar_sampling, i)
      }
    }
    # Sample sigma squared using all tessellations to predict the outcome variables
    sigmaSquared[i] <- sampleSigmaSquared(
      yScaled,
      nu,
      lambda,
      sumOfAllTess
    )
    
    for (j in 1:m) {
      # Propose new Tessellation for component j
      newTessOutput <- proposeTessellation(
        tess[[j]],
        dim[[j]],
        sd,
        covariateIndices,
        numCovariates
      )
      tess_j_star <- newTessOutput[[1]]
      dim_j_star <- newTessOutput[[2]]
      modification <- newTessOutput[[3]]
      
      # Retrieve old indices from cache
      indexes <- current_indices[[j]]
      # Calculate new indices for the proposal
      indexesStar <- cellIndices(xScaled, tess_j_star, dim_j_star)
      
      residualsOutput <- calculateResiduals(
        y = yScaled,
        j = j,
        SumOfAllTess = sumOfAllTess,
        Pred = pred,
        lastTessPred = lastTessPred,
        indexes = indexes,
        indexesStar = indexesStar,
        num_centres_new = nrow(tess_j_star) 
      )
      
      rIjOld <- residualsOutput[[1]]
      nIjOld <- residualsOutput[[2]]
      rIjNew <- residualsOutput[[3]]
      nIjNew <- residualsOutput[[4]]
      sumOfAllTess <- residualsOutput[[5]]
      
      if (!any(nIjNew == 0)) {
        # Call the acceptanceProbability function 
        logAcceptanceProb <- acceptanceProbability(
          rIjOld, nIjOld,
          rIjNew, nIjNew,
          tess_j_star, dim_j_star,
          sigmaSquared[i],
          modification,
          sigmaSquaredMu,
          omega,
          lambdaRate,
          numCovariates
        )
        
        if (log(runif(n = 1)) < logAcceptanceProb) {
          # Accept proposal: update lists IN-PLACE
          tess[[j]] <- tess_j_star
          dim[[j]] <- dim_j_star
          current_indices[[j]] <- indexesStar
          
          pred[[j]] <- sampleMuValues(
            j, tess,
            rIjNew, nIjNew,
            sigmaSquaredMu,
            sigmaSquared[i]
          )
          lastTessPred <- pred[[j]][indexesStar]
          
        } else {
          # Reject proposal
          pred[[j]] <- sampleMuValues(
            j, tess, rIjOld, nIjOld,
            sigmaSquaredMu, sigmaSquared[i]
          )
          lastTessPred <- pred[[j]][indexes]
        }
      } else {
        # Reject proposal (empty cell)
        pred[[j]] <- sampleMuValues(
          j, tess, rIjOld, nIjOld,
          sigmaSquaredMu, sigmaSquared[i]
        )
        lastTessPred <- pred[[j]][indexes]
      }
      
      if (j == m) {
        sumOfAllTess <- sumOfAllTess + lastTessPred
      }
    }
    
    if (i > mcmcBurnIn & (i - mcmcBurnIn) %% thinning == 0) {
      # vectors that hold the predictions for each iteration after burn in.
      predictionMatrix[, (i - mcmcBurnIn) / thinning] <- sumOfAllTess
    }
    
    # Store the posterior samples
    if (numPosteriorSamplesToStore > 0 &&
        i > mcmcBurnIn &
        (i - mcmcBurnIn) %% thinning == 0) {
      # Store the current state of tess, dim, pred, sigma
      outputPosteriorTess[[currentStorageIdx]] <- tess
      outputPosteriorDim[[currentStorageIdx]] <- dim
      outputPosteriorPred[[currentStorageIdx]] <- pred
      outputPosteriorSigma[currentStorageIdx] <- sigmaSquared[i]
      currentStorageIdx <- currentStorageIdx + 1
    }
  } # End of MCMC Loop
  
  # Close any remaining progress bar
  if (showProgress) {
    if (!is.null(pbar_sampling)) {
      close(pbar_sampling)
      cat("\n")
    } else if (!is.null(pbar_burnin)) {
      close(pbar_burnin)
      cat("\n")
    }
    cat("MCMC sampling completed.\n\n")
  }
  
  # Finding the mean of the prediction over the iterations and then unscaling
  # the predictions.
  meanYhat <- (rowSums(predictionMatrix) / (posteriorSamples)) * yRange +
    yCentre
  
  # Create and return the AddiVortesFit object
  new_AddiVortesFit(
    posteriorTess = outputPosteriorTess,
    posteriorDim = outputPosteriorDim,
    posteriorSigma = outputPosteriorSigma,
    posteriorPred = outputPosteriorPred,
    xCentres = xCentres,
    xRanges = xRanges,
    yCentre = yCentre,
    yRange = yRange,
    inSampleRmse = sqrt(mean((y - meanYhat)^2))
  )
}
