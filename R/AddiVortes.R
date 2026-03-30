#' @title AddiVortes
#'
#' @description
#' The AddiVortes function is a Bayesian nonparametric regression model that uses a
#' soft tessellation to model the relationship between the covariates and the output values.
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
    centres <- apply(data, 2, function(x) (max(x) + min(x)) / 2)
    ranges <- apply(data, 2, function(x) diff(range(x)))
    ranges[ranges == 0] <- 1
    scaledData <- scale(data, center = centres, scale = ranges)
  } else {
    centres <- (max(data) + min(data)) / 2
    ranges <- diff(range(data))
    if (ranges == 0) ranges <- 1
    scaledData <- scale(data, center = centres, scale = ranges)
  }
  list(scaledData = as.matrix(scaledData), centres = centres, ranges = ranges)
}

AddiVortes <- function (y, x, m = 200, totalMCMCIter = 1200, mcmcBurnIn = 200, 
                        nu = 6, q = 0.85, k = 3, sd = 0.8, Omega = min(3, ncol(x)), 
                        LambdaRate = 25, IntialSigma = "Linear", thinning = 1, showProgress = TRUE, 
                        distancePower = 2.0, p_shape = 2.0, p_rate = 1.0, p_sd = 0.1) 
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
  
  pred <- lapply(1:m, function(ignoredIndex) matrix(mean(yScaled)/m, nrow = n, ncol = 1))
  dim <- lapply(1:m, function(ignoredIndex) as.integer(sample(1:p, 1)))
  tess <- lapply(1:m, function(ignoredIndex) matrix(rnorm(1, 0, sd)))
  
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
  
  p_vec <- rep(as.numeric(distancePower), m)
  
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
                             as.integer(if(showProgress) 1L else 0L),
                             p_vec,
                             as.numeric(p_shape),
                             as.numeric(p_rate),
                             as.numeric(p_sd))
  
  if (showProgress) {
    cat("MCMC sampling completed.\n\n")
  }
  
  posteriorSamples <- floor((totalMCMCIter - mcmcBurnIn)/thinning)
  meanYhat <- (rowSums(super_call_result$predictionMatrix)/(posteriorSamples)) * yRange + yCentre
  
  final_result <- list(
    posteriorTess = super_call_result$posteriorTess, 
    posteriorDim = super_call_result$posteriorDim, 
    posteriorSigma = super_call_result$posteriorSigma, 
    posteriorPower = super_call_result$posteriorPower,
    posteriorMu = super_call_result$posteriorMu,
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