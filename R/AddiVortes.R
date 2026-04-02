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

AddiVortes <- function (y, x, m = 50, totalMCMCIter = 2500, mcmcBurnIn = 1000, 
                        nu = 6, q = 0.85, k = 3, sd = 0.8, Omega = 1, 
                        LambdaRate = 25, IntialSigma = "Linear", thinning = 1, showProgress = TRUE,
                        alpha = 1, a_alpha = 0.5, b_alpha = 1, rho_alpha = ncol(x), dirichletWarmup = NULL,
                        adaptBoost = 1, adaptPenalty = 1, momentumDecay = 0.90) 
{
  if (is.null(dirichletWarmup)) {
    dirichletWarmup <- floor(mcmcBurnIn / 2)
  }
  
  yScalingResult <- scaleData_internal(y)
  yScaled <- yScalingResult$scaledData
  yCentre <- yScalingResult$centres
  yRange <- yScalingResult$ranges
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
  
  pred <- rep(list(matrix(mean(yScaled)/m)), m)
  dim <- lapply(1:m, function(ignoredIndex) as.integer(sample(1:p, 1)))
  tess <- lapply(1:m, function(ignoredIndex) matrix(rnorm(1, 0, sd)))
  
  sqdist <- lapply(1:m, function(i) {
    matrix((xScaled[, dim[[i]]] - as.numeric(tess[[i]]))^2, ncol = 1)
  })
  
  sumOfAllTess <- rep(mean(yScaled), length(yScaled))
  storage.mode(sumOfAllTess) <- "double"
  SigmaSquaredMu <- (0.5/(k * sqrt(m)))^2
  
  if (IntialSigma == "Naive") {
    SigmaSquaredHat <- var(yScaled)
  } else {
    multiLinear <- lm(yScaled ~ xScaled)
    SigmaSquaredHat <- sum(multiLinear$residuals^2)/(length(yScaled) - p - 1)
  }
  
  lambda_invgamma <- optim(par = 1, fittingFunction, method = "Brent", 
                           lower = 0.001, upper = 100, q = q, nu = nu, SigmaSquaredHat = SigmaSquaredHat)$par
  
  current_indices <- lapply(1:m, function(idx) rep.int(1L, nrow(xScaled)))
  
  if (showProgress) {
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
                             as.numeric(momentumDecay))
  
  posteriorSamples <- floor((totalMCMCIter - mcmcBurnIn)/thinning)
  meanYhat <- (rowSums(super_call_result$predictionMatrix)/(posteriorSamples)) * yRange + yCentre
  
  inclusion_indicator_matrix <- super_call_result$posteriorVariableSelection > 0
  ensemble_inclusion_probabilities <- rowMeans(inclusion_indicator_matrix)
  
  dirichlet_weights_mean <- rowMeans(super_call_result$posteriorDirichletWeights)
  dirichlet_weights_lower <- apply(super_call_result$posteriorDirichletWeights, 1, quantile, probs = 0.025)
  dirichlet_weights_upper <- apply(super_call_result$posteriorDirichletWeights, 1, quantile, probs = 0.975)
  
  final_result <- list(
    posteriorTess = super_call_result$posteriorTess, 
    posteriorDim = super_call_result$posteriorDim, 
    posteriorSigma = super_call_result$posteriorSigma, 
    posteriorPred = super_call_result$posteriorPred, 
    xCentres = xCentres, xRanges = xRanges, yCentre = yCentre, yRange = yRange, 
    inSampleRmse = sqrt(mean((y - meanYhat)^2)),
    posteriorDirichletWeightsMean = dirichlet_weights_mean,
    posteriorDirichletWeightsLower = dirichlet_weights_lower,
    posteriorDirichletWeightsUpper = dirichlet_weights_upper,
    ensembleInclusionProbabilities = ensemble_inclusion_probabilities,
    posteriorAlphaTrace = super_call_result$posteriorAlpha
  )
  class(final_result) <- "AddiVortesFit"
  
  return(final_result)
}