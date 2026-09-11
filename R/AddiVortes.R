#' @title Fit an AddiVortes regression or classification model
#'
#' @description
#' The AddiVortes model is a Bayesian nonparametric model that uses additive
#' Voronoi tessellations to relate covariates to a response. For a numeric
#' response the model is Gaussian regression. For a classification response
#' it uses a probit link with Albert-Chib latent variables (binary) or
#' independent multinomial probit latents (three or more classes). The task is
#' chosen automatically from `y`.
#'
#' The function can handle multiple types of covariates, including continuous,
#' spherical and categorical. Categorical covariates are automatically detected.
#' By default (`cat.onehot = TRUE`) they are one-hot encoded, with the first
#' level of each categorical variable used as the reference category; the
#' `catScaling` parameter then controls the weight of categorical differences in
#' distance calculations. Setting `cat.onehot = FALSE` instead keeps each
#' categorical covariate as a single integer-coded column and uses Eskin
#' distance (Eskin et al., 2002). For spherical covariates, the function assumes
#' that the final spherical dimension corresponds to the polar angle, which has
#' a range of 0 to 2*pi. The `metric` parameter can be used to specify the type
#' of each covariate (Euclidean, Spherical, or Categorical), and the `members`
#' parameter can indicate membership of covariates into different subspaces when
#' using multiple spheres in covariate space.
#'
#' @param y A vector of response values. Numeric `y` is treated as regression,
#'   except when it has exactly two unique values in `{0, 1}`, which is
#'   binary classification. Factor, character and logical vectors are treated as
#'   classification: two levels give a binary probit model and three or more
#'   levels give a multinomial probit model. The first factor level (or 0 for
#'   numeric 0/1 responses) is the reference class. Missing values are not
#'   allowed.
#' @param x A matrix or data frame of the covariates. Character and factor columns
#'   are treated as categorical variables and automatically converted to d-1 binary
#'   indicator variables via one-hot encoding (with the first level as reference).
#' @param m The number of tessellations. For multinomial classification this is
#'   the number of tessellations **per latent dimension** (there are
#'   \eqn{K-1}{K-1} latents for \eqn{K}{K} classes).
#' @param totalMCMCIter The number of iterations.
#' @param mcmcBurnIn The number of burn in iterations.
#' @param nu The degrees of freedom for the inverse-gamma prior on the residual
#'   variance. Ignored for classification, where the latent residual variance is
#'   fixed at 1.
#' @param q The quantile used to set the inverse-gamma prior on the residual
#'   variance. Ignored for classification.
#' @param k Prior scale for tessellation output values. For regression,
#'   \eqn{\sigma_\mu = 0.5/(k\sqrt{m})}{sigma_mu = 0.5/(k sqrt(m))} on the scaled
#'   response. For classification,
#'   \eqn{\sigma_\mu = 3/(k\sqrt{m})}{sigma_mu = 3/(k sqrt(m))} on the latent
#'   probit scale.
#' @param sd The standard deviation used in centre proposals.
#' @param Omega Omega/(number of covariates) is the prior probability of
#'   adding a dimension.
#' @param LambdaRate The rate of the Poisson distribution for the number of centres.
#' @param InitialSigma The method used to calculate the initial residual
#'   variance for regression (`"Linear"` or `"Naive"`). Ignored for
#'   classification.
#' @param thinning The thinning rate.
#' @param metric Either "E" (Euclidean, default), "S" (Spherical), or "C" (Categorical).
#' @param members If needed, indicates membership of covariates into different
#'   subspaces (needed if using multiple spheres in covariate space). Default `NULL`.
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
#' @param cat.onehot Should categorical covariates be one-hot encoded? Default
#'   `TRUE`. When `TRUE`, each categorical covariate with *d* levels is expanded
#'   to *d* − 1 binary indicators and distances are Euclidean (weighted by
#'   \code{catScaling}). When `FALSE`, categories are kept as a single
#'   integer-coded column and mismatches use Eskin distance (Eskin et al., 2002),
#'   with squared cost \eqn{2 / d^2} when levels differ and 0 when they match;
#'   \code{catScaling} is then ignored. See the categorical covariates vignette
#'   for a comparison and guidance on which to use.
#' @param showProgress Logical; if TRUE, a progress bar is shown during fitting.
#'
#' @return An AddiVortes object containing the posterior samples of the
#' tessellations, dimensions and predictions, plus per-iteration trace
#' statistics used by `traceplots()`. Classification fits also store `task`,
#' `classLevels`, `nLatents` and in-sample accuracy.
#'
#' @references
#' Stone, A. and Gosling, J.P. (2025). AddiVortes: (Bayesian) additive Voronoi
#' tessellations. *Journal of Computational and Graphical Statistics*.
#'
#' Stone, A.J., Ogundimu, E. and Gosling, J.P. (2026). Binary AddiVortes:
#' (Bayesian) Additive Voronoi Tessellations for Binary Classification with an
#' application to Predicting Home Mortgage Application Outcomes.
#'
#' Albert, J.H. and Chib, S. (1993). Bayesian analysis of binary and
#' polychotomous response data. *Journal of the American Statistical
#' Association*, 88(422), 669–679.
#'
#' Kindo, B.P., Wang, H. and Peña, E.A. (2016). Multinomial probit Bayesian
#' additive regression trees. *Stat*, 5(1), 171–181.
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
#'
#' # Binary classification is selected automatically from a 0/1 or factor y
#' set.seed(789)
#' x_clf <- matrix(runif(80), 40, 2)
#' y_clf <- factor(ifelse(x_clf[, 1] + x_clf[, 2] > 1, "yes", "no"))
#' fit_clf <- AddiVortes(y_clf, x_clf, m = 8, totalMCMCIter = 80,
#'                       mcmcBurnIn = 20, showProgress = FALSE)
#' p_clf <- predict(fit_clf, x_clf, type = "response", showProgress = FALSE)
#' cls_clf <- predict(fit_clf, x_clf, type = "class", showProgress = FALSE)
#' }
#'
#' @importFrom stats var lm optim quantile runif rnorm dbinom dpois qnorm uniroot pnorm
#' @export
AddiVortes <- function(y, x, m = 200,
                       totalMCMCIter = 1200,
                       mcmcBurnIn = 200,
                       nu = 6, q = 0.85,
                       k = 3, sd = 0.8,
                       Omega = min(3, ncol(x)),
                       LambdaRate = 5,
                       InitialSigma = "Linear",
                       thinning = 1,
                       metric = "E",
                       members = NULL,
                       catScaling = 1,
                       cat.onehot = TRUE,
                       showProgress = interactive()) {
  # Force evaluation of Omega using the *original* x before categorical encoding
  # replaces x with the encoded matrix. Without this, R's lazy evaluation would
  # use ncol() of the encoded matrix, potentially making Omega = NumCovariates
  # and causing prob = 1 in the dimension acceptance ratio (which produces 0/0 = NaN).
  force(Omega)
  ### Pre-processing data

  if (NROW(x) != length(y)) {
    stop("'y' must have length equal to the number of rows of 'x'.", call. = FALSE)
  }
  xOriginal <- x
  responseInfo <- infer_response_type_internal(y)
  task <- responseInfo$task
  classLevels <- responseInfo$classLevels
  yClass <- responseInfo$yClass
  nLatents <- responseInfo$nLatents
  mPerLatent <- m
  mTotal <- m * nLatents
  isClassification <- task != "regression"

  #### Encode categorical covariates -------------------------------------------
  if (!is.numeric(catScaling) || length(catScaling) != 1 || catScaling <= 0) {
    stop("'catScaling' must be a single positive number.")
  }
  #### Dealing with choice of metric -------------------------------------------
  if (length(metric) == 1) {
    if (metric == "E" || metric == "Euc" || metric == "Euclidean") {
      metric <- rep("E", ncol(x))
    } else if (metric == "S" || metric == "Sphere" || metric == "Spherical") {
      metric <- rep("S", ncol(x))
    } else if (metric == "C" || metric == "Cat" || metric == "Categorical") {
      metric <- rep("C", ncol(x))
    }
  }

  old_metric <- metric
  old_metric[old_metric == "E" | old_metric == "Euc" |
               old_metric == "Euclidean"] <- 0
  old_metric[old_metric == "S" | old_metric == "Sphere" |
               old_metric == "Spherical"] <- 1
  old_metric[old_metric == "C" | old_metric == "Cat" |
               old_metric == "Categorical"] <- 2
  old_metric <- as.integer(old_metric)
  old_members <- if(is.null(members)) NULL else as.integer(members)

  san_data <- covariateStructure_internal(x,
                                          metric, members, cat.onehot)
  encResult <- encodeCategories_internal(san_data$data, catScaling = catScaling)
  catEncoding <- encResult$encoding
  covariateSummary <- formatCovariateSummary_internal(x,
                                                       metric,
                                                       catEncoding,
                                                       cat.onehot)
  if (cat.onehot)
    x <- encResult$encoded
  else
    x <- as.matrix(san_data$data)
  
  members <- as.integer(san_data$membership)
  
  metric <- san_data$structure
  metric[metric == "E"] <- 0
  metric[metric == "S"] <- 1
  metric[metric == "C"] <- 2
  metric <- as.integer(metric)
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
  if (isClassification) {
    yScaled <- as.double(rep(0, length(y)))
    yCentre <- 0
    yRange <- 1
  } else {
    yScalingResult <- scaleData_internal(y)
    yScaled <- yScalingResult$scaledData # Vector of values
    yCentre <- yScalingResult$centres
    yRange <- yScalingResult$ranges
  }
  
  xScalingResult <- scaleData_internal(x)
  xScaled <- xScalingResult$scaledData # Matrix of values
  xCentres <- xScalingResult$centres # Vector of values
  xRanges <- xScalingResult$ranges # Vector of values
  
  ##### Dealing with unscaled data ---------------------------------------------
  xScaled[, metric != 0] <- x[, metric != 0]
  # Binary columns from categorical encoding keep their {0, catScaling} values
  # rather than being further scaled, so they directly control distance weight
  if (cat.onehot && !is.null(catEncoding)) {
    binaryCols <- catEncoding$encodedBinaryCols
    xScaled[, binaryCols] <- x[, binaryCols]
  }
  mus <- rep(0, ncol(x))
  mus[metric != 0] <- xCentres[metric != 0]
  
  #### Handling NULL sigma choice and ensuring it's vectorised
  sd <- sapply(xRanges,
               function(r) uniroot(function(x) 
                 qnorm(0.75, 0, x) - r / 2, c(0, r))$root)
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
  pred <- rep(list(matrix(if (isClassification) 0 else mean(yScaled) / m)), mTotal)
  dim <- sapply(seq_len(mTotal), function(ignoredIndex) {
    list(sample(seq_len(ncol(x)), 1))
  })
  tess <- sapply(seq_len(mTotal), function(ignoredIndex) {
    list(matrix(rnorm(1, 0, sd)))
  })
  ## Make sure that tessellation proposals are within the region, if periodic
  if (!is.null(sphere_ranges)) {
    for (i in seq_along(tess)) {
      if (metric[dim[[i]]] == 1) {
        sph_ind <- sum(metric[1:dim[[i]]] == 1)
        while (tess[[i]][1, 1] > sphere_ranges[[sph_ind]][2]) {
          tess[[i]][1, 1] <- tess[[i]][1, 1] - sphere_ranges[[sph_ind]][2]
        }
        while (tess[[i]][1, 1] < sphere_ranges[[sph_ind]][1]) {
          tess[[i]][1, 1] <- tess[[i]][1, 1] + sphere_ranges[[sph_ind]][2]
        }
      }
    }
  }
  ## Constrain initial centres for binary (one-hot) dimensions to [0, catScaling]
  ## Initialization always produces single-dimension tessellations, so dim[[i]] is
  ## a scalar. Guard against any future multi-column initial states by iterating
  ## over all active dims.
  if (cat.onehot && !is.null(catEncoding) && length(catEncoding$encodedBinaryCols) > 0) {
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
  } else {
    for (i in seq_along(tess)) {
      if (metric[dim[[i]]] == 2) {
        tess[[i]][1,1] <- sample(unique(xScaled[,dim[[i]]]), 1)
      }
    }
  }
  
  #### Set-up MCMC -------------------------------------------------------------
  # The variance that captures variability around the mean of the scaled y
  # values (regression) or the latent probit scale (classification).
  if (isClassification) {
    SigmaSquaredMu <- (3 / (k * sqrt(m)))^2
    lambda <- 1
  } else {
    SigmaSquaredMu <- (0.5 / (k * sqrt(m)))^2

    # Finding lambda
    if (InitialSigma == "Naive") {
      # Usually used if p is greater then n. Uses std dev of y to predict sigma.
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
  }
  
  # Normalise tess, dim and pred to plain R objects before passing to C++.
  # sapply may wrap results in lists or simplify to arrays; ensure each
  # element is a plain double matrix / integer vector / double vector.
  init_tess <- lapply(seq_len(mTotal), function(j) {
    t_j <- tess[[j]]
    if (is.list(t_j)) t_j <- t_j[[1]]
    m_j <- as.matrix(t_j)
    storage.mode(m_j) <- "double"
    m_j
  })
  init_dim  <- lapply(dim,  function(d) as.integer(unlist(d)))
  init_pred <- lapply(pred, function(p_j) as.double(p_j))
  
  # Binary column indices for categorical clamping (NULL when not applicable)
  binaryCols_r <- if (cat.onehot && !is.null(catEncoding) && 
                      length(catEncoding$encodedBinaryCols) > 0) {
    as.integer(catEncoding$encodedBinaryCols)
  } else {
    NULL
  }
  catScaling_r <- if (cat.onehot && !is.null(catEncoding)) catEncoding$catScaling else 0.0
  
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
    cat("Model configuration: ", mTotal,
        " tessellation", if (mTotal == 1L) "" else "s",
        if (isClassification && nLatents > 1L) {
          paste0(" (", m, " per latent, ", nLatents, " latents)")
        } else {
          ""
        },
        ", ", totalMCMCIter,
        " total iterations (", mcmcBurnIn,
        " burn-in)",
        if (isClassification) {
          paste0("\nTask: ", task, " classification with classes ",
                 paste(classLevels, collapse = ", "))
        } else {
          ""
        },
        "\n\n",
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
    as.integer(members),
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
    as.logical(showProgress),
    as.logical(isClassification),
    if (isClassification) as.integer(yClass) else NULL,
    as.integer(nLatents)
  )
  
  if (showProgress) cat("Done.\n\n")
  
  outputPosteriorTess  <- mcmcResult$posteriorTess
  outputPosteriorDim   <- mcmcResult$posteriorDim
  outputPosteriorPred  <- mcmcResult$posteriorPred
  outputPosteriorSigma <- mcmcResult$posteriorSigma
  predictionMatrix     <- mcmcResult$predictionMatrix
  traceStats           <- mcmcResult$traceStats
  
  posteriorSamples <- ncol(predictionMatrix)

  inSampleRmse <- NA_real_
  inSampleAccuracy <- NA_real_
  inSampleBrier <- NA_real_
  if (task == "regression") {
    meanYhat <- if (posteriorSamples > 0) {
      (rowSums(predictionMatrix) / posteriorSamples) * yRange + yCentre
    } else {
      rep(yCentre, length(y))
    }
    inSampleRmse <- sqrt(mean((y - meanYhat)^2))
  }

  # Create and return the AddiVortes object
  if (cat.onehot)
    this_enc <- catEncoding
  else
    this_enc <- NULL
  fit <- new_AddiVortes(
    posteriorTess = outputPosteriorTess,
    posteriorDim = outputPosteriorDim,
    posteriorSigma = outputPosteriorSigma,
    posteriorPred = outputPosteriorPred,
    xCentres = xCentres,
    xRanges = xRanges,
    yCentre = yCentre,
    yRange = yRange,
    inSampleRmse = inSampleRmse,
    metric = old_metric,
    members = old_members,
    metric_aug = metric,
    member_aug = members,
    catEncoding = this_enc,
    traceStats = traceStats,
    task = task,
    classLevels = classLevels,
    nLatents = nLatents,
    mPerLatent = mPerLatent,
    inSampleAccuracy = inSampleAccuracy,
    inSampleBrier = inSampleBrier
  )

  if (isClassification && posteriorSamples > 0) {
    in_sample <- predict(fit, xOriginal, type = "response", showProgress = FALSE)
    if (task == "binary") {
      y01 <- as.numeric(yClass)
      pred01 <- as.numeric(in_sample > 0.5)
      fit$inSampleAccuracy <- mean(pred01 == y01)
      fit$inSampleBrier <- mean((in_sample - y01)^2)
    } else {
      pred_idx <- max.col(in_sample, ties.method = "first")
      fit$inSampleAccuracy <- mean(pred_idx == (yClass + 1L))
    }
  }

  fit
}

countCovariateTypes_internal <- function(x, metric) {
  if (!is.data.frame(x)) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  }
  
  is_categorical <- vapply(x, function(col) is.character(col) ||
                             is.factor(col), logical(1))
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

formatCovariateSummary_internal <- function(x, metric, catEncoding = NULL,
                                            coh = TRUE) {
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
    if (coh) {
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
    } else {
      lines <- c(
        lines,
        paste0(
          "Categorical covariates use Eskin distance (cat.onehot = FALSE): ",
          "each stays as one integer-coded column, with mismatch cost 2/d^2 ",
          "for a variable with d levels."
        )
      )
    }
  }
  
  lines
}
