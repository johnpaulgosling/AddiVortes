#' @title AddiVortes
#'
#' @description
#' The AddiVortes function is a Bayesian nonparametric regression model that uses a
#' tessellation to model the relationship between the covariates and the output values.
#' The model uses a backfitting algorithm to sample from the posterior distribution of
#' the output values for each tessellation. The function returns the RMSE value for
#' the test samples.
#' 
#' For spherical data, it is assumed that the final spherical dimension is the
#' polar angle: i.e. that with range [0, 2*pi].
#'
#' @param y A vector of the output values.
#' @param x A matrix or data frame of the covariates. Character and factor columns
#'   are treated as categorical variables and automatically converted to d-1 binary
#'   indicator variables via one-hot encoding (with the first level as reference).
#' @param m The number of tessellations.
#' @param totalMCMCIter The number of iterations.
#' @param mcmcBurnIn The number of burn in iterations.
#' @param nu The degrees of freedom.
#' @param q The quantile.
#' @param k The number of centres.
#' @param sd The standard deviation used in centre proposals.
#' @param Omega Omega/(number of covariates) is the prior probability of adding a dimension.
#' @param LambdaRate The rate of the Poisson distribution for the number of centres.
#' @param InitialSigma The method used to calculate the initial variance.
#' @param thinning The thinning rate.
#' @param metric Either "E" (Euclidean, default) or "S" (Spherical).
#' @param catScaling Numeric scalar controlling the scale of binary indicator
#'   variables created from categorical covariates. The indicator takes values
#'   0 or \code{catScaling} (default 1). Larger values give categorical differences
#'   more weight in the distance calculations.
#' @param showProgress Logical; if TRUE, progress bars and messages are shown during fitting.
#'
#' @return An AddiVortes object containing the posterior samples of the
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
#' 
#' # Larger example with categorical covariates (d=2 and d=3) and a test set
#' set.seed(456)
#' n_train <- 200
#' n_test  <- 50
#' x_train <- data.frame(
#'   x1   = rnorm(n_train),
#'   x2   = runif(n_train),
#'   grp2 = sample(c("A", "B"), n_train, replace = TRUE),
#'   grp3 = sample(c("low", "mid", "high"), n_train, replace = TRUE)
#' )
#' y_train <- x_train$x1 + ifelse(x_train$grp2 == "B", 1, 0) + rnorm(n_train, sd = 0.5)
#' 
#' fit2 <- AddiVortes(y_train, x_train, m = 10, totalMCMCIter = 200, mcmcBurnIn = 50,
#'                    catScaling = 1, showProgress = FALSE)
#' 
#' x_test <- data.frame(
#'   x1   = rnorm(n_test),
#'   x2   = runif(n_test),
#'   grp2 = sample(c("A", "B"), n_test, replace = TRUE),
#'   grp3 = sample(c("low", "mid", "high"), n_test, replace = TRUE)
#' )
#' y_test <- x_test$x1 + ifelse(x_test$grp2 == "B", 1, 0) + rnorm(n_test, sd = 0.5)
#' 
#' preds <- predict(fit2, x_test, showProgress = FALSE)
#' test_rmse <- sqrt(mean((y_test - preds)^2))
#' }
#'
#' @importFrom stats var lm optim quantile runif dbinom dpois
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
AddiVortes <- function(y, x, m = 200,
                       totalMCMCIter = 1200,
                       mcmcBurnIn = 200,
                       nu = 6, q = 0.85,
                       k = 3, sd = 0.8, 
                       Omega = min(3, ncol(x)), 
                       LambdaRate = 25,
                       InitialSigma = "Linear",
                       thinning = 1,
                       metric = "E",
                       catScaling = 1,
                       showProgress = interactive()) {
  #### Encode categorical covariates -------------------------------------------
  if (!is.numeric(catScaling) || length(catScaling) != 1 || catScaling <= 0)
    stop("'catScaling' must be a single positive number.")
  encResult <- encodeCategories_internal(x, catScaling = catScaling)
  catEncoding <- encResult$encoding
  x <- encResult$encoded

  #### Dealing with choice of metric -------------------------------------------
  if (length(metric) == 1) {
    if (metric == "E" || metric == "Euc" || metric == "Euclidean")
      metric = rep(0, ncol(x))
    else if (metric == "S" || metric == "Sphere" || metric == "Spherical")
      metric = rep(1, ncol(x))
  }
  metric[metric == "E" | metric == "Euc" | metric == "Euclidean"] <- 0
  metric[metric == "S" | metric == "Sphere" | metric == "Spherical"] <- 1
  metric <- as.integer(metric)
  if (1 %in% metric) {
    sphere_ranges <- list()
    for (i in seq_len(sum(metric == 1)-1))
      sphere_ranges[[length(sphere_ranges)+1]] <- c(-pi/2, pi/2)
    sphere_ranges[[length(sphere_ranges)+1]] <- c(-pi, pi)
  }
  else
    sphere_ranges <- NULL
  
  #### Scaling x and y ---------------------------------------------------------
  yScalingResult <- scaleData_internal(y)
  yScaled <- yScalingResult$scaledData # Vector of values
  yCentre <- yScalingResult$centres
  yRange <- yScalingResult$ranges
  
  xScalingResult <- scaleData_internal(x)
  xScaled <- xScalingResult$scaledData # Matrix of values
  xCentres <- xScalingResult$centres # Vector of values
  xRanges <- xScalingResult$ranges # Vector of values
  
  ##### Dealing with unscaled data ---------------------------------------------
  xScaled[,metric != 0] <- x[,metric != 0]
  # Binary columns from categorical encoding keep their {0, catScaling} values
  # rather than being further scaled, so they directly control distance weight
  if (!is.null(catEncoding)) {
    binaryCols <- catEncoding$encodedBinaryCols
    xScaled[, binaryCols] <- x[, binaryCols]
  }
  mus <- rep(0, nrow(x))
  mus[metric != 0] <- xCentres[metric != 0]
  
  #### Handling NULL sigma choice and ensuring it's vectorised
  sd <- sapply(xRanges, function(r) uniroot(function(x) qnorm(0.75,0,x)-r/2, c(0,r))$root)
  sd[metric == 0] <- 0.8
  
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
  
  if (Omega > p) {
    message(
      "Note: Omega (", Omega, ") exceeds number of covariates (", p, "). ",
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
  ## Make sure that tessellation proposals are within the region, if periodic
  if (!is.null(sphere_ranges)) {
    for (i in seq_along(tess)) {
      if (metric[dim[[i]]] == 1) {
        sph_ind <- sum(metric[1:dim[[i]]] == 1)
        while (tess[[i]][1,1] > sphere_ranges[[sph_ind]][2])
          tess[[i]][1,1] <- tess[[i]][1,1] - sphere_ranges[[sph_ind]][2]
        while (tess[[i]][1,1] < sphere_ranges[[sph_ind]][1])
          tess[[i]][1,1] <- tess[[i]][1,1] + sphere_ranges[[sph_ind]][2]
      }
    }
  }
  ## Constrain initial centres for binary (one-hot) dimensions to [0, catScaling]
  ## Initialization always produces single-dimension tessellations, so dim[[i]] is
  ## a scalar. Guard against any future multi-column initial states by iterating
  ## over all active dims.
  if (!is.null(catEncoding) && length(catEncoding$encodedBinaryCols) > 0) {
    binaryColsInit <- catEncoding$encodedBinaryCols
    cs <- catEncoding$catScaling
    for (i in seq_along(tess)) {
      active_dims <- as.integer(unlist(dim[[i]]))
      local_bin_pos <- which(active_dims %in% binaryColsInit)
      if (length(local_bin_pos) > 0) {
        for (lp in local_bin_pos) {
          tess[[i]][, lp] <- runif(nrow(tess[[i]]), 0, cs)
        }
      }
    }
  }
  
  #### Set-up MCMC -------------------------------------------------------------
  # Prepare some variables used in the backfitting algorithm.
  # We start off with the mean of the scaled y values as the prediction for all
  # tessellations.
  sumOfAllTess <- rep(
    mean(yScaled),
    length(yScaled)
  )
  # The variance that captures variability around the mean of the scaled y values.
  SigmaSquaredMu <- (0.5 / (k * sqrt(m)))^2
  lastTessPred <- matrix
  
  # Matrices that will hold the samples from the posterior distribution
  # for the training samples and test samples.
  posteriorSamples <- floor((totalMCMCIter - mcmcBurnIn) / thinning)
  predictionMatrix <- array(dim = c(
    length(y),
    posteriorSamples
  ))
  
  # Finding lambda
  if (InitialSigma == "Naive") {
    # Usually used if p is greater then n. Uses Standard deviation of y to predict sigma.
    SigmaSquaredHat <- var(yScaled)
  } else {
    # Default method using residual standard deviation from a least-squared linear
    # regression of y, to predict sigma.
    multiLinear <- lm(yScaled ~ xScaled)
    SigmaSquaredHat <- sum(multiLinear$residuals^2) /
      (length(yScaled) - length(xScaled[1, ]) - 1)
  }
  lambda <- optim(
    par = 1,
    fittingFunction,
    method = "Brent",
    lower = 0.001,
    upper = 100,
    q = q, nu = nu,
    SigmaSquaredHat = SigmaSquaredHat
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
  SigmaSquared <- NULL
  
  currentStorageIdx <- 1 # Index for the new output lists
  
  # Some precalculations
  NumCovariates <- ncol(xScaled)
  covariateIndices <- seq_len(NumCovariates)
  currentIndices <- vector("list", m)
  for(k in 1:m) {
    currentIndices[[k]] <- cellIndices(xScaled, tess[[k]], dim[[k]], metric)
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
    SigmaSquared[i] <- sampleSigmaSquared(
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
        mus,
        covariateIndices,
        NumCovariates,
        metric
      )
      tess_j_star <- newTessOutput[[1]]
      dim_j_star <- newTessOutput[[2]]
      modification <- newTessOutput[[3]]
      
      ## Clamp proposed centres for binary (one-hot) dimensions to [0, catScaling]
      if (!is.null(catEncoding) && length(catEncoding$encodedBinaryCols) > 0) {
        local_bin_pos <- which(dim_j_star %in% catEncoding$encodedBinaryCols)
        if (length(local_bin_pos) > 0) {
          cs <- catEncoding$catScaling
          for (lp in local_bin_pos) {
            tess_j_star[, lp] <- pmin(pmax(tess_j_star[, lp], 0), cs)
          }
        }
      }
      
      # Retrieve old indices from cache
      indexes <- currentIndices[[j]]
      # Calculate new indices for the proposal
      indexesStar <- cellIndices(xScaled, tess_j_star, dim_j_star, metric)
      
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
          SigmaSquared[i],
          modification,
          SigmaSquaredMu,
          Omega,
          LambdaRate,
          NumCovariates
        )
        
        if (log(runif(n = 1)) < logAcceptanceProb) {
          # Accept proposal: update lists IN-PLACE
          tess[[j]] <- tess_j_star
          dim[[j]] <- dim_j_star
          currentIndices[[j]] <- indexesStar
          
          pred[[j]] <- sampleMuValues(
            j, tess,
            rIjNew, nIjNew,
            SigmaSquaredMu,
            SigmaSquared[i]
          )
          lastTessPred <- pred[[j]][indexesStar]
          
        } else {
          # Reject proposal
          pred[[j]] <- sampleMuValues(
            j, tess, rIjOld, nIjOld,
            SigmaSquaredMu, SigmaSquared[i]
          )
          lastTessPred <- pred[[j]][indexes]
        }
      } else {
        # Reject proposal (empty cell)
        pred[[j]] <- sampleMuValues(
          j, tess, rIjOld, nIjOld,
          SigmaSquaredMu, SigmaSquared[i]
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
      outputPosteriorSigma[currentStorageIdx] <- SigmaSquared[i]
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
  
  # Create and return the AddiVortes object
  new_AddiVortes(
    posteriorTess = outputPosteriorTess,
    posteriorDim = outputPosteriorDim,
    posteriorSigma = outputPosteriorSigma,
    posteriorPred = outputPosteriorPred,
    xCentres = xCentres,
    xRanges = xRanges,
    yCentre = yCentre,
    yRange = yRange,
    inSampleRmse = sqrt(mean((y - meanYhat)^2)),
    metric = metric,
    catEncoding = catEncoding
  )
}
