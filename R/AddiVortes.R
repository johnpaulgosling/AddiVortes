#' @title AddiVortes_Local
#'
#' @description
#' The AddiVortes function is a Bayesian nonparametric regression model that uses a
#' tessellation to model the relationship between the covariates and the output values.
#' The model uses a backfitting algorithm to sample from the posterior distribution of
#' the output values for each tessellation. The function returns the RMSE value for
#' the test samples.
#'
#' @param y A vector of the output values.
#' @param x A matrix or data frame of the covariates.
#' @param m The number of tessellations.
#' @param totalMCMCIter The number of iterations.
#' @param mcmcBurnIn The number of burn in iterations.
#' @param nu The degrees of freedom.
#' @param q The quantile.
#' @param k The number of centres.
#' @param sd The standard deviation used in centre proposals.
#' @param Omega Omega/(number of covariates) is the prior probability of adding a dimension.
#' @param LambdaRate The rate of the Poisson distribution for the number of centres.
#' @param IntialSigma The method used to calculate the initial variance.
#' @param thinning The thinning rate.
#' @param showProgress Logical; if TRUE, progress bars and messages are shown during fitting.
#'
#' @return An AddiVortes object containing the posterior samples of the
#' tessellations, dimensions and predictions.
#'
#' @importFrom stats var lm optim quantile runif dbinom dpois
#' @export
#' 

#' @keywords internal
#' @noRd
qinvgamma_internal <- function(p, shape, rate) {
  1 / stats::qgamma(1 - p, shape = shape, rate = rate)
}

#' @keywords internal
#' @noRd
fittingFunction <- function(lambda, q, nu, SigmaSquaredHat) {
  (qinvgamma_internal(q, shape = nu / 2, rate = nu * lambda / 2) - SigmaSquaredHat)^2
}

#' @keywords internal
#' @noRd
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

AddiVortes <- function (y, x, m = 200, totalMCMCIter = 1200, mcmcBurnIn = 200, 
                        nu = 6, q = 0.85, k = 3, sd = 0.8, Omega = min(3, ncol(x)), 
                        LambdaRate = 25, IntialSigma = "Linear", thinning = 1, showProgress = TRUE) 
{
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
    cat(sprintf("Input dimensions: %d observations, %d covariates\n", nrow(xScaled), ncol(xScaled)))
    cat(sprintf("Model configuration: %d tessellations, %d total iterations (%d burn-in)\n\n", m, totalMCMCIter, mcmcBurnIn))
  }
  
  super_call_result <- .Call("super_mcmc_loop_cpp",
                             xScaled,
                             yScaled,
                             sumOfAllTess,
                             tess,
                             dim,
                             current_indices,
                             sqdist,
                             pred,
                             as.integer(m),
                             p_int,
                             sd_dbl,
                             mus_dbl,
                             as.numeric(SigmaSquaredMu),
                             as.numeric(Omega),
                             as.numeric(LambdaRate),
                             as.integer(totalMCMCIter),
                             as.integer(mcmcBurnIn),
                             as.integer(thinning),
                             as.numeric(nu),
                             as.numeric(lambda_invgamma),
                             as.integer(if(showProgress) 1L else 0L))
  
  if (showProgress) {
    cat("MCMC sampling completed.\n\n")
  }
  
  posteriorSamples <- floor((totalMCMCIter - mcmcBurnIn)/thinning)
  meanYhat <- (rowSums(super_call_result$predictionMatrix)/(posteriorSamples)) * yRange + yCentre
  
  final_result <- list(
    posteriorTess = super_call_result$posteriorTess, 
    posteriorDim = super_call_result$posteriorDim, 
    posteriorSigma = super_call_result$posteriorSigma, 
    posteriorPred = super_call_result$posteriorPred, 
    xCentres = xCentres, 
    xRanges = xRanges, 
    yCentre = yCentre, 
    yRange = yRange, 
    inSampleRmse = sqrt(mean((y - meanYhat)^2))
  )
  class(final_result) <- "AddiVortesFit"
  
  return(final_result)
}