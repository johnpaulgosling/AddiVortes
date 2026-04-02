#' @title AddiVortes
#'
#' @description
#' The AddiVortes function is a Bayesian nonparametric regression model that
#' uses a tessellation to model the relationship between the covariates and
#' the output values. The model uses a backfitting algorithm to sample from
#' the posterior distribution of the output values for each tessellation.
#' The function returns the RMSE value for the test samples.
#'
#' For spherical data, coordinates must be supplied in radians. The last
#' spherical column is treated as the azimuthal (longitude) dimension with
#' range \eqn{[-\pi, \pi]}; all preceding spherical columns are treated as
#' polar (latitude) dimensions with range \eqn{[-\pi/2, \pi/2]}.
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
#' @param Omega Omega/(number of covariates) is the prior probability of
#'   adding a dimension.
#' @param LambdaRate The rate of the Poisson distribution for the number of centres.
#' @param InitialSigma The method used to calculate the initial variance.
#' @param thinning The thinning rate.
#' @param metric The distance metric to use for each covariate column. Either a
#'   single string (\code{"E"} / \code{"Euclidean"} for all columns Euclidean,
#'   \code{"S"} / \code{"Spherical"} for all columns spherical) or a vector of
#'   the same length as the number of covariate columns containing \code{"E"},
#'   \code{"S"}, or the integers \code{0} (Euclidean) / \code{1} (Spherical),
#'   e.g. \code{metric = c("S", "E", "S")} for a dataset where columns 1 and 3
#'   are spherical and column 2 is Euclidean. Spherical coordinates must be in
#'   radians: polar (latitude-type) dimensions in \eqn{[-\pi/2, \pi/2]} and the
#'   last spherical (azimuthal / longitude) dimension in \eqn{[-\pi, \pi]}.
#' @param catScaling Numeric scalar controlling the scale of binary indicator
#'   variables created from categorical covariates. Each binary indicator takes
#'   values 0 (reference level) or \code{catScaling} (non-reference level).
#'   The default value of 1 matches the range of continuous covariates, which are
#'   normalised to \code{[-0.5, 0.5]} (range = 1) during fitting, so categorical
#'   differences receive comparable weight to continuous differences in the distance
#'   calculations. Increase above 1 to give categorical differences more weight;
#'   decrease below 1 to give them less weight. Binary indicator columns are named
#'   \code{<colname>_<level>} (e.g. a column \code{grp} with levels \code{"A"},
#'   \code{"B"}, \code{"C"} produces columns \code{grp_B} and \code{grp_C}, with
#'   \code{"A"} as the reference level).
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
#' n_test <- 50
#' x_train <- data.frame(
#'   x1   = rnorm(n_train),
#'   x2   = runif(n_train),
#'   grp2 = sample(c("A", "B"), n_train, replace = TRUE),
#'   grp3 = sample(c("low", "mid", "high"), n_train, replace = TRUE)
#' )
#' y_train <- x_train$x1 + ifelse(x_train$grp2 == "B", 1, 0) + rnorm(n_train, sd = 0.5)
#'
#' fit2 <- AddiVortes(y_train, x_train,
#'   m = 10, totalMCMCIter = 200, mcmcBurnIn = 50,
#'   catScaling = 1, showProgress = FALSE
#' )
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
#' @importFrom stats var lm optim quantile runif rnorm dbinom dpois
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
  # Force evaluation of Omega using the *original* x before categorical encoding
  # replaces x with the encoded matrix. Without this, R's lazy evaluation would
  # use ncol() of the encoded matrix, potentially making Omega = NumCovariates
  # and causing prob = 1 in acceptanceProbability (which produces 0/0 = NaN).
  force(Omega)
  #### Encode categorical covariates -------------------------------------------
  if (!is.numeric(catScaling) || length(catScaling) != 1 || catScaling <= 0) {
    stop("'catScaling' must be a single positive number.")
  }
  encResult <- encodeCategories_internal(x, catScaling = catScaling)
  catEncoding <- encResult$encoding
  covariateSummary <- formatCovariateSummary_internal(x, metric, catEncoding)
  x <- encResult$encoded

  #### Dealing with choice of metric -------------------------------------------
  if (length(metric) == 1) {
    if (metric == "E" || metric == "Euc" || metric == "Euclidean") {
      metric <- rep(0, ncol(x))
    } else if (metric == "S" || metric == "Sphere" || metric == "Spherical") {
      metric <- rep(1, ncol(x))
    }
  }
  metric[metric == "E" | metric == "Euc" | metric == "Euclidean"] <- 0
  metric[metric == "S" | metric == "Sphere" | metric == "Spherical"] <- 1
  metric <- as.integer(metric)
  # When the user supplies a per-column metric for the original (pre-encoding)
  # columns and categorical expansion has added extra binary columns, extend the
  # metric automatically: binary indicator columns are always Euclidean (0).
  if (!is.null(catEncoding) && length(metric) == catEncoding$origNCols) {
    enc_metric <- integer(0)
    for (j in seq_len(catEncoding$origNCols)) {
      if (j %in% catEncoding$catColIndices) {
        n_bin <- length(catEncoding$colEncodings[[j]]$levels) - 1L
        enc_metric <- c(enc_metric, rep(0L, n_bin))
      } else {
        enc_metric <- c(enc_metric, metric[j])
      }
    }
    metric <- enc_metric
  }
  if (length(metric) != ncol(x)) {
    stop(
      "'metric' must have one element per covariate column (", ncol(x),
      " required, ", length(metric), " supplied)."
    )
  }
  if (!all(metric %in% 0:1)) {
    stop("'metric' must contain only 0 (Euclidean) or 1 (Spherical) values.")
  }
  if (1 %in% metric) {
    sphere_ranges <- list()
    for (i in seq_len(sum(metric == 1) - 1)) {
      sphere_ranges[[length(sphere_ranges) + 1]] <- c(-pi / 2, pi / 2)
    }
    sphere_ranges[[length(sphere_ranges) + 1]] <- c(-pi, pi)
  } else {
    sphere_ranges <- NULL
  }

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
  xScaled[, metric != 0] <- x[, metric != 0]
  # Binary columns from categorical encoding keep their {0, catScaling} values
  # rather than being further scaled, so they directly control distance weight
  if (!is.null(catEncoding)) {
    binaryCols <- catEncoding$encodedBinaryCols
    xScaled[, binaryCols] <- x[, binaryCols]
  }
  mus <- rep(0, ncol(x))
  mus[metric != 0] <- xCentres[metric != 0]

  #### Handling NULL sigma choice and ensuring it's vectorised
  sd <- sapply(xRanges, function(r) uniroot(function(x) qnorm(0.75, 0, x) - r / 2, c(0, r))$root)
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
    list(sample(seq_len(ncol(x)), 1))
  })
  tess <- sapply(1:m, function(ignoredIndex) {
    list(matrix(rnorm(1, 0, sd)))
  })
  ## Make sure that tessellation proposals are within the region, if periodic
  if (!is.null(sphere_ranges)) {
    for (i in seq_along(tess)) {
      if (metric[dim[[i]]] == 1) {
        sph_ind <- sum(metric[1:dim[[i]]] == 1)
        lim <- sphere_ranges[[sph_ind]][2]  # half-period: pi/2 or pi
        val <- tess[[i]][1, 1]
        while (val >= lim) val <- val - 2 * lim
        while (val < -lim) val <- val + 2 * lim
        tess[[i]][1, 1] <- val
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
  # The variance that captures variability around the mean of the scaled y values.
  SigmaSquaredMu <- (0.5 / (k * sqrt(m)))^2

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

  # Normalise tess, dim and pred to plain R objects before passing to C++.
  # sapply may wrap results in lists or simplify to arrays; ensure each
  # element is a plain double matrix / integer vector / double vector.
  init_tess <- lapply(seq_len(m), function(j) {
    t_j <- tess[[j]]
    if (is.list(t_j)) t_j <- t_j[[1]]
    m_j <- as.matrix(t_j)
    storage.mode(m_j) <- "double"
    m_j
  })
  init_dim  <- lapply(dim,  function(d) as.integer(unlist(d)))
  init_pred <- lapply(pred, function(p_j) as.double(p_j))

  # Binary column indices for categorical clamping (NULL when not applicable)
  binaryCols_r <- if (!is.null(catEncoding) && length(catEncoding$encodedBinaryCols) > 0) {
    as.integer(catEncoding$encodedBinaryCols)
  } else {
    NULL
  }
  catScaling_r <- if (!is.null(catEncoding)) catEncoding$catScaling else 0.0

  # Progress message
  if (showProgress) {
    cat("Fitting AddiVortes model to input data...\n")
    if (length(covariateSummary) > 0) {
      cat(paste(covariateSummary, collapse = "\n"), "\n\n", sep = "")
    }
    cat("Input dimensions: ", nrow(xScaled),
      " observations, ", ncol(xScaled),
      " covariates\n",
      sep = ""
    )
    cat("Model configuration: ", m,
      " tessellations, ", totalMCMCIter,
      " total iterations (", mcmcBurnIn,
      " burn-in)\n\n",
      sep = ""
    )
    cat("Running MCMC...\n")
  }

  #### MCMC (single C++ call) --------------------------------------------------
  mcmcResult <- .Call(
    "addi_vortes_mcmc_cpp",
    matrix(as.double(xScaled), nrow = nrow(xScaled), ncol = ncol(xScaled)),
    as.double(yScaled),
    as.integer(metric),
    as.integer(m),
    as.integer(totalMCMCIter),
    as.integer(mcmcBurnIn),
    as.integer(thinning),
    as.double(nu),
    as.double(lambda),
    as.double(SigmaSquaredMu),
    as.double(Omega),
    as.double(LambdaRate),
    as.double(sd),
    as.double(mus),
    init_tess,
    init_dim,
    init_pred,
    binaryCols_r,
    as.double(catScaling_r),
    as.logical(showProgress)
  )

  if (showProgress) cat("MCMC sampling completed.\n\n")

  outputPosteriorTess  <- mcmcResult$posteriorTess
  outputPosteriorDim   <- mcmcResult$posteriorDim
  outputPosteriorPred  <- mcmcResult$posteriorPred
  outputPosteriorSigma <- mcmcResult$posteriorSigma
  predictionMatrix     <- mcmcResult$predictionMatrix

  posteriorSamples <- ncol(predictionMatrix)

  # Finding the mean of the prediction over the iterations and then unscaling
  # the predictions.
  meanYhat <- if (posteriorSamples > 0) {
    (rowSums(predictionMatrix) / posteriorSamples) * yRange + yCentre
  } else {
    rep(yCentre, length(y))
  }

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

countCovariateTypes_internal <- function(x, metric) {
  if (!is.data.frame(x)) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  }

  is_categorical <- vapply(x, function(col) is.character(col) || is.factor(col), logical(1))
  metric_vec <- rep_len(metric, ncol(x))
  metric_vec <- as.character(metric_vec)
  metric_vec[metric_vec %in% c("E", "Euc", "Euclidean")] <- "0"
  metric_vec[metric_vec %in% c("S", "Sphere", "Spherical")] <- "1"
  metric_vec <- suppressWarnings(as.integer(metric_vec))

  list(
    continuous = sum(metric_vec[!is_categorical] == 0L, na.rm = TRUE),
    spherical = sum(metric_vec[!is_categorical] == 1L, na.rm = TRUE),
    categorical = sum(is_categorical)
  )
}

formatCovariateSummary_internal <- function(x, metric, catEncoding = NULL) {
  counts <- countCovariateTypes_internal(x, metric)
  parts <- character(0)

  if (counts$continuous > 0) {
    parts <- c(parts, sprintf("%d continuous", counts$continuous))
  }
  if (counts$spherical > 0) {
    parts <- c(parts, sprintf("%d spherical", counts$spherical))
  }
  if (counts$categorical > 0) {
    parts <- c(parts, sprintf("%d categorical", counts$categorical))
  }

  if (length(parts) == 0) {
    return(character(0))
  }

  lines <- c(sprintf("Covariate summary: %s.", paste(parts, collapse = ", ")))

  if (counts$categorical > 0) {
    if (!is.null(catEncoding) && !is.null(catEncoding$encodedBinaryCols)) {
      binary_cols <- length(catEncoding$encodedBinaryCols)
      lines <- c(
        lines,
        sprintf(
          paste0(
            "Categorical covariates are expanded to %d one-hot encoded ",
            "binary column%s, with the first level of each categorical ",
            "variable used as the reference category."
          ),
          binary_cols,
          if (binary_cols == 1L) "" else "s"
        )
      )
    } else {
      lines <- c(
        lines,
        paste0(
          "Categorical covariates are expanded to one-hot encoded binary ",
          "columns, with the first level of each categorical variable used ",
          "as the reference category."
        )
      )
    }
  }

  lines
}
