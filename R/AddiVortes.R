#' @title AddiVortes_Local
#'
#' @description
#' The AddiVortes function is a Bayesian nonparametric regression model that uses a
#' tessellation to model the relationship between the covariates and the output values.
#' The model uses a backfitting algorithm to sample from the posterior distribution of
#' the output values for each tessellation. The function returns the RMSE value for
#' the test samples.
#'
#' @importFrom stats var lm optim quantile runif dbinom dpois
#' @export

qinvgamma_internal <- function(p, shape, rate) {
  1 / stats::qgamma(1 - p, shape = shape, rate = rate)
}

fittingFunction <- function(lambda, q, nu, SigmaSquaredHat) {
  (qinvgamma_internal(q, shape = nu / 2, rate = nu * lambda / 2) - SigmaSquaredHat)^2
}

scaleData_internal <- function(data) {
  if (is.matrix(data) || is.data.frame(data)) {
    centres <- apply(data, 2, mean)
    ranges <- apply(data, 2, function(x) diff(range(x)))
    ranges[ranges == 0] <- 1
    scaledData <- scale(data, center = centres, scale = ranges)
  } else {
    centres <- mean(data)
    ranges <- diff(range(data))
    if (ranges == 0) ranges <- 1
    scaledData <- scale(data, center = centres, scale = ranges)
  }
  list(scaledData = as.matrix(scaledData), centres = centres, ranges = ranges)
}

AddiVortes <- function (y, x, m = 200, totalMCMCIter = 2500, mcmcBurnIn = 1000, 
                        nu = 6, q = 0.85, k = 3, sd = 0.8, Omega = 1, 
                        LambdaRate = 25, IntialSigma = "Linear", thinning = 1, showProgress = TRUE,
                        alpha = 1, a_alpha = 0.5, b_alpha = 1, rho_alpha = ncol(x), dirichletWarmup = NULL,
                        adaptBoost = 1, adaptPenalty = 1, momentumDecay = 0.90, kappa = 0.60,
                        varSelMode = 2, numChains = 4,
                        splitMode = 1, power = 2.0, p_shape = 2.0, p_rate = 2.0, p_sd = 1) 
{
  if (is.null(dirichletWarmup)) {
    dirichletWarmup <- floor(mcmcBurnIn / 2)
  }
  
  unique_y <- unique(y)
  is_classification <- length(unique_y) == 2
  
  if (is_classification) {
    yScaled <- as.numeric(as.factor(y)) - 1
    yCentre <- 0
    yRange <- 1
    if (showProgress) cat("Binary response detected. Engaging Probit classification mode.\n")
  } else {
    yScalingResult <- scaleData_internal(y)
    yScaled <- yScalingResult$scaledData
    yCentre <- yScalingResult$centres
    yRange <- yScalingResult$ranges
  }
  
  xScalingResult <- scaleData_internal(x)
  xScaled <- xScalingResult$scaledData
  xCentres <- xScalingResult$centres
  xRanges <- xScalingResult$ranges
  
  n <- length(y)
  p <- ncol(xScaled)
  
  storage.mode(xScaled) <- "double"
  storage.mode(yScaled) <- "double"
  p_int <- as.integer(p)
  sd_dbl <- as.numeric(sd)
  mus_dbl <- rep(0.0, p)
  
  splitMode_int <- as.integer(splitMode)
  p_vec <- rep(as.numeric(power), m)
  is_class_int <- as.integer(is_classification)
  
  chain_results <- lapply(1:numChains, function(chain) {
    if (showProgress && numChains > 1) {
      cat(sprintf("\nInitialising and running MCMC chain %d of %d...\n", chain, numChains))
    }
    
    pred <- rep(list(matrix(mean(yScaled)/m)), m)
    dim <- lapply(1:m, function(ignoredIndex) as.integer(sample(1:p, 1)))
    tess <- lapply(1:m, function(ignoredIndex) matrix(rnorm(1, 0, sd)))
    
    sqdist <- lapply(1:m, function(i) {
      matrix((xScaled[, dim[[i]]] - as.numeric(tess[[i]]))^2, ncol = 1)
    })
    
    if (is_classification) {
      # Map raw proportion to the latent probit scale, clamping to avoid infinity
      p_hat <- mean(yScaled)
      latent_offset <- qnorm(max(0.01, min(0.99, p_hat))) 
      
      pred <- rep(list(matrix(latent_offset / m)), m)
      sumOfAllTess <- rep(latent_offset, length(yScaled))
      
      # Scale the prior variance to span [-3, 3] in the latent space
      SigmaSquaredMu <- (3.0 / (k * sqrt(m)))^2 
      
      SigmaSquaredHat <- 1.0
      lambda_invgamma <- 1.0 
    } else {
      # Standard continuous initialisation and prior
      pred <- rep(list(matrix(mean(yScaled) / m)), m)
      sumOfAllTess <- rep(mean(yScaled), length(yScaled))
      
      SigmaSquaredMu <- (0.5 / (k * sqrt(m)))^2
      
      if (IntialSigma == "Naive") {
        SigmaSquaredHat <- var(yScaled)
      } else if (IntialSigma == "LASSO") {
        cv_fit <- glmnet::cv.glmnet(xScaled, yScaled, alpha = 1)
        best_idx <- which(cv_fit$lambda == cv_fit$lambda.min)
        SigmaSquaredHat <- cv_fit$cvm[best_idx]
      } else {
        multiLinear <- lm(yScaled ~ xScaled)
        SigmaSquaredHat <- sum(multiLinear$residuals^2) / (length(yScaled) - p - 1)
      }
      lambda_invgamma <- optim(par = 1, fittingFunction, method = "Brent", lower = 0.001, upper = 100, q = q, nu = nu, SigmaSquaredHat = SigmaSquaredHat)$par
    }
    
    if (!is_classification) {
      print(paste("Initial Sigma Squared Estimate for Chain", chain, ":", round(SigmaSquaredHat*yRange^2, 4)))
    }
    
    current_indices <- lapply(1:m, function(idx) rep.int(1L, nrow(xScaled)))
    
    if (showProgress && numChains == 1) {
      cat("Fitting AddiVortes model to input data...\n")
    }
    
    super_call_result <- .Call("super_mcmc_loop_cpp",
                               xScaled, yScaled, sumOfAllTess, tess, dim,
                               current_indices, sqdist, pred, as.integer(m),
                               p_int, sd_dbl, mus_dbl, as.numeric(SigmaSquaredMu),
                               as.numeric(Omega), as.numeric(LambdaRate),
                               as.integer(totalMCMCIter), as.integer(mcmcBurnIn),
                               as.integer(thinning), as.numeric(nu),
                               as.numeric(lambda_invgamma),
                               as.integer(if(showProgress) 1L else 0L),
                               as.numeric(alpha), as.numeric(a_alpha),
                               as.numeric(b_alpha), as.numeric(rho_alpha),
                               as.integer(dirichletWarmup),
                               as.numeric(adaptBoost),
                               as.numeric(adaptPenalty),
                               as.numeric(momentumDecay),
                               as.numeric(kappa),
                               as.integer(varSelMode),
                               splitMode_int,        
                               as.numeric(p_vec),    
                               as.numeric(p_shape),  
                               as.numeric(p_rate),   
                               as.numeric(p_sd),
                               is_class_int)         
    
    return(super_call_result)
  })
  
  posteriorTessCombined <- do.call(c, lapply(chain_results, function(res) res$posteriorTess))
  posteriorDimCombined <- do.call(c, lapply(chain_results, function(res) res$posteriorDim))
  posteriorSigmaCombined <- do.call(c, lapply(chain_results, function(res) res$posteriorSigma))
  posteriorPredCombined <- do.call(c, lapply(chain_results, function(res) res$posteriorPred))
  
  posteriorPowerCombined <- do.call(cbind, lapply(chain_results, function(res) res$posteriorPower))
  posteriorMuCombined <- do.call(c, lapply(chain_results, function(res) res$posteriorMu))
  
  predictionMatrixCombined <- do.call(cbind, lapply(chain_results, function(res) res$predictionMatrix))
  dirichletWeightsCombined <- do.call(cbind, lapply(chain_results, function(res) res$posteriorDirichletWeights))
  variableSelectionCombined <- do.call(cbind, lapply(chain_results, function(res) res$posteriorVariableSelection))
  posteriorAlphaCombined <- do.call(c, lapply(chain_results, function(res) res$posteriorAlpha))
  
  posteriorSamplesPerChain <- floor((totalMCMCIter - mcmcBurnIn)/thinning)
  totalPosteriorSamples <- posteriorSamplesPerChain * numChains
  
  if (is_classification) {
    meanYhat <- pnorm(rowSums(predictionMatrixCombined)/totalPosteriorSamples)
    in_sample_error <- mean((meanYhat > 0.5) != y)
  } else {
    meanYhat <- (rowSums(predictionMatrixCombined)/totalPosteriorSamples) * yRange + yCentre
    in_sample_error <- sqrt(mean((y - meanYhat)^2))
  }
  
  inclusion_indicator_matrix <- variableSelectionCombined > 0
  ensemble_inclusion_probabilities <- rowMeans(inclusion_indicator_matrix)
  
  dirichlet_weights_mean <- rowMeans(dirichletWeightsCombined)
  dirichlet_weights_lower <- apply(dirichletWeightsCombined, 1, quantile, probs = 0.025)
  dirichlet_weights_upper <- apply(dirichletWeightsCombined, 1, quantile, probs = 0.975)
  
  unscaledPosteriorSigmaSquared <- posteriorSigmaCombined * (yRange^2)
  posteriorPowerCombined <- do.call(cbind, lapply(chain_results, function(res) res$posteriorPower))
  
  final_result <- new_AddiVortesFit(
    posteriorTess = posteriorTessCombined, 
    posteriorDim = posteriorDimCombined, 
    posteriorSigma = unscaledPosteriorSigmaSquared, 
    posteriorPred = if(splitMode == 1) posteriorPredCombined else posteriorMuCombined, 
    posteriorDirichletWeights = dirichletWeightsCombined,
    posteriorVariableSelection = variableSelectionCombined,
    posteriorAlpha = posteriorAlphaCombined,
    predictionMatrix = predictionMatrixCombined,
    posteriorPower = posteriorPowerCombined,
    splitMode = splitMode,
    xCentres = xCentres, xRanges = xRanges, yCentre = yCentre, yRange = yRange, 
    inSampleRmse = in_sample_error,
    posteriorDirichletWeightsMean = dirichlet_weights_mean,
    posteriorDirichletWeightsLower = dirichlet_weights_lower,
    posteriorDirichletWeightsUpper = dirichlet_weights_upper,
    ensembleInclusionProbabilities = ensemble_inclusion_probabilities,
    isClassification = is_classification
  )
  class(final_result) <- "AddiVortesFit"
  
  return(final_result)
}