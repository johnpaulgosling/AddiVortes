# --- Test Suite for AddiVortes Class Methods Syntax Issues ---
# Tests for checking input validation, error handling, and syntax issues
# for all AddiVortes class methods

# --- Helper: Create a simple test AddiVortes object ---
create_test_object <- function() {
  set.seed(42)
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  
  AddiVortes(Y, X, m = 3, totalMCMCIter = 20, mcmcBurnIn = 5, showProgress = FALSE)
}

# --- Tests for new_AddiVortes() constructor ---

test_that("new_AddiVortes requires all arguments", {
  # Test that constructor fails without required arguments
  expect_error(
    new_AddiVortes(
      posteriorTess = list(),
      posteriorDim = list(),
      posteriorSigma = numeric()
      # Missing arguments
    ),
    "argument .* is missing"
  )
})

test_that("new_AddiVortes creates correct class", {
  # Test that constructor creates object of correct class
  obj <- new_AddiVortes(
    posteriorTess = list(),
    posteriorDim = list(),
    posteriorSigma = numeric(),
    posteriorPred = list(),
    xCentres = c(0, 0),
    xRanges = c(1, 1),
    yCentre = 0,
    yRange = 1,
    inSampleRmse = 0.5,
    mcmcLog = data.frame()
  )
  
  expect_s3_class(obj, "AddiVortes")
  expect_true(inherits(obj, "AddiVortes"))
})

# --- Tests for print.AddiVortes() ---

test_that("print.AddiVortes requires AddiVortes object", {
  # Test that print fails with non-AddiVortes object
  not_addivortes <- list(a = 1, b = 2)
  
  expect_error(
    print.AddiVortes(not_addivortes),
    "must be an object of class 'AddiVortes'"
  )
})

test_that("print.AddiVortes works with valid object", {
  # Test that print works with a valid AddiVortes object
  obj <- create_test_object()
  
  expect_output(print(obj), "AddiVortes Model")
  expect_output(print(obj), "Model Formula:")
  expect_output(print(obj), "In-sample RMSE:")
})

test_that("print.AddiVortes handles empty posterior samples", {
  # Test print with object containing no posterior samples
  obj <- new_AddiVortes(
    posteriorTess = list(),
    posteriorDim = list(),
    posteriorSigma = numeric(),
    posteriorPred = list(),
    xCentres = c(0),
    xRanges = c(1),
    yCentre = 0,
    yRange = 1,
    inSampleRmse = 0.5,
    mcmcLog = data.frame()
  )
  
  expect_output(print(obj), "No posterior samples available")
})

# --- Tests for summary.AddiVortes() ---

test_that("summary.AddiVortes requires AddiVortes object", {
  # Test that summary fails with non-AddiVortes object
  not_addivortes <- list(a = 1, b = 2)
  
  expect_error(
    summary.AddiVortes(not_addivortes),
    "`object` must be an object of class 'AddiVortes'"
  )
})

test_that("summary.AddiVortes works with valid object", {
  # Test that summary works with a valid AddiVortes object
  obj <- create_test_object()
  
  expect_output(summary(obj), "AddiVortes Model")
  expect_output(summary(obj), "Model Information:")
})

# --- Tests for predict.AddiVortes() ---

test_that("predict.AddiVortes requires AddiVortes object", {
  # Test that predict fails with non-AddiVortes object
  not_addivortes <- list(a = 1, b = 2)
  X_new <- matrix(rnorm(25), 5, 5)
  
  expect_error(
    predict.AddiVortes(not_addivortes, X_new),
    "must be of class 'AddiVortes'"
  )
})

test_that("predict.AddiVortes requires matrix newdata", {
  # Test that predict fails if newdata is not a matrix
  obj <- create_test_object()
  
  expect_error(
    predict(obj, newdata = c(1, 2, 3, 4, 5)),
    "must be a matrix"
  )
  
  expect_error(
    predict(obj, newdata = data.frame(a = 1:5, b = 1:5)),
    "must be a matrix"
  )
})

test_that("predict.AddiVortes checks newdata dimensions", {
  # Test that predict fails if newdata has wrong number of columns
  obj <- create_test_object()
  X_wrong <- matrix(rnorm(30), 10, 3)  # Wrong number of columns
  
  expect_error(
    predict(obj, X_wrong),
    "Number of columns.*does not match"
  )
})

test_that("predict.AddiVortes works with valid inputs", {
  # Test that predict works with valid inputs
  obj <- create_test_object()
  X_new <- matrix(rnorm(25), 5, 5)
  
  preds <- predict(obj, X_new, showProgress = FALSE)
  
  expect_type(preds, "double")
  expect_length(preds, 5)
  expect_false(any(is.na(preds)))
})

test_that("predict.AddiVortes type argument works", {
  # Test different type arguments
  obj <- create_test_object()
  X_new <- matrix(rnorm(25), 5, 5)
  
  # Test response type
  preds_response <- predict(obj, X_new, type = "response", showProgress = FALSE)
  expect_type(preds_response, "double")
  expect_length(preds_response, 5)
  
  # Test quantile type
  preds_quantile <- predict(obj, X_new, type = "quantile", 
                           quantiles = c(0.025, 0.975), showProgress = FALSE)
  expect_true(is.matrix(preds_quantile))
  expect_equal(nrow(preds_quantile), 5)
  expect_equal(ncol(preds_quantile), 2)
})

test_that("predict.AddiVortes interval argument works", {
  # Test different interval arguments
  obj <- create_test_object()
  X_new <- matrix(rnorm(25), 5, 5)
  
  # Test credible intervals
  preds_credible <- predict(obj, X_new, type = "quantile",
                           interval = "credible", showProgress = FALSE)
  expect_true(is.matrix(preds_credible))
  
  # Test prediction intervals
  preds_prediction <- predict(obj, X_new, type = "quantile",
                             interval = "prediction", showProgress = FALSE)
  expect_true(is.matrix(preds_prediction))
})

test_that("predict.AddiVortes handles empty posterior samples", {
  # Test predict with object containing no posterior samples
  obj <- new_AddiVortes(
    posteriorTess = list(),
    posteriorDim = list(),
    posteriorSigma = numeric(),
    posteriorPred = list(),
    xCentres = c(0, 0, 0, 0, 0),
    xRanges = c(1, 1, 1, 1, 1),
    yCentre = 0,
    yRange = 1,
    inSampleRmse = 0.5,
    mcmcLog = data.frame()
  )
  X_new <- matrix(rnorm(25), 5, 5)
  
  expect_warning(
    predict(obj, X_new, showProgress = FALSE),
    "contains no posterior samples"
  )
})

# --- Tests for plot.AddiVortes() ---

test_that("plot.AddiVortes requires AddiVortes object", {
  # Test that plot fails with non-AddiVortes object
  not_addivortes <- list(a = 1, b = 2)
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  
  expect_error(
    plot.AddiVortes(not_addivortes, X, Y),
    "must be an object of class 'AddiVortes'"
  )
})

test_that("plot.AddiVortes requires x_train and y_train", {
  # Test that plot fails without required training data
  obj <- create_test_object()
  
  expect_error(
    plot(obj),
    "x_train.*and.*y_train.*must be provided"
  )
})

test_that("plot.AddiVortes checks x_train is a matrix", {
  # Test that plot fails if x_train is not a matrix
  obj <- create_test_object()
  Y <- rnorm(10)
  
  expect_error(
    plot(obj, x_train = c(1, 2, 3), y_train = Y),
    "must be a matrix"
  )
})

test_that("plot.AddiVortes checks y_train is numeric", {
  # Test that plot fails if y_train is not numeric
  obj <- create_test_object()
  X <- matrix(rnorm(50), 10, 5)
  
  expect_error(
    plot(obj, x_train = X, y_train = c("a", "b", "c")),
    "must be a numeric vector"
  )
})

test_that("plot.AddiVortes checks dimensions match", {
  # Test that plot fails if dimensions don't match
  obj <- create_test_object()
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(5)  # Wrong length
  
  expect_error(
    plot(obj, x_train = X, y_train = Y),
    "number of rows.*must match"
  )
})

test_that("plot.AddiVortes handles empty posterior samples", {
  # Test plot with object containing no posterior samples
  obj <- new_AddiVortes(
    posteriorTess = list(),
    posteriorDim = list(),
    posteriorSigma = numeric(),
    posteriorPred = list(),
    xCentres = c(0, 0, 0, 0, 0),
    xRanges = c(1, 1, 1, 1, 1),
    yCentre = 0,
    yRange = 1,
    inSampleRmse = 0.5,
    mcmcLog = data.frame()
  )
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  
  expect_error(
    plot(obj, x_train = X, y_train = Y),
    "No posterior samples available"
  )
})

test_that("plot.AddiVortes validates which parameter", {
  # Test that plot validates the which parameter
  obj <- create_test_object()
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  
  # Invalid which values should be filtered out
  expect_error(
    plot(obj, x_train = X, y_train = Y, which = c(5, 6, 7)),
    "which.*must contain values between 1 and 4"
  )
})

# --- Tests for mcmcLog field ---

test_that("AddiVortes object contains mcmcLog field", {
  obj <- create_test_object()
  
  expect_true("mcmcLog" %in% names(obj))
  expect_s3_class(obj$mcmcLog, "data.frame")
})

test_that("mcmcLog has correct columns", {
  obj <- create_test_object()
  
  expected_cols <- c("iteration", "tessellation", "move_type", "log_alpha", "accepted")
  expect_true(all(expected_cols %in% names(obj$mcmcLog)))
})

test_that("mcmcLog has correct number of rows", {
  set.seed(42)
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  m <- 3
  totalMCMCIter <- 20
  
  obj <- AddiVortes(Y, X, m = m, totalMCMCIter = totalMCMCIter,
                    mcmcBurnIn = 5, showProgress = FALSE)
  
  expect_equal(nrow(obj$mcmcLog), totalMCMCIter * m)
})

test_that("mcmcLog move_type values are valid", {
  obj <- create_test_object()
  
  valid_types <- c("AD", "RD", "AC", "RC", "Change", "Swap")
  expect_true(all(obj$mcmcLog$move_type %in% valid_types))
})

test_that("mcmcLog accepted column is logical", {
  obj <- create_test_object()
  
  expect_type(obj$mcmcLog$accepted, "logical")
})

test_that("mcmcLog log_alpha is numeric", {
  obj <- create_test_object()
  
  expect_type(obj$mcmcLog$log_alpha, "double")
})

test_that("mcmcLog iteration and tessellation indices are correct", {
  set.seed(42)
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  m <- 3
  totalMCMCIter <- 10
  
  obj <- AddiVortes(Y, X, m = m, totalMCMCIter = totalMCMCIter,
                    mcmcBurnIn = 5, showProgress = FALSE)
  
  expect_true(all(obj$mcmcLog$iteration >= 1 & obj$mcmcLog$iteration <= totalMCMCIter))
  expect_true(all(obj$mcmcLog$tessellation >= 1 & obj$mcmcLog$tessellation <= m))
})

test_that("mcmcLog auto-rejected proposals have log_alpha = -Inf and accepted = FALSE", {
  obj <- create_test_object()
  
  # Any auto-rejected row (if any exist) must have log_alpha == -Inf and accepted == FALSE
  auto_rejected <- obj$mcmcLog[is.infinite(obj$mcmcLog$log_alpha) &
                                  obj$mcmcLog$log_alpha < 0, ]
  if (nrow(auto_rejected) > 0) {
    expect_true(all(!auto_rejected$accepted))
  }
})

test_that("mcmcLog variable column exists and is integer", {
  obj <- create_test_object()
  
  expect_true("variable" %in% names(obj$mcmcLog))
  expect_type(obj$mcmcLog$variable, "integer")
})

test_that("mcmcLog variable is NA for non-AD/RD moves", {
  obj <- create_test_object()
  log <- obj$mcmcLog
  
  non_dim_rows <- log[!(log$move_type %in% c("AD", "RD")), ]
  if (nrow(non_dim_rows) > 0) {
    expect_true(all(is.na(non_dim_rows$variable)))
  }
})

test_that("mcmcLog variable is a valid covariate index for AD/RD moves", {
  set.seed(99)
  X <- matrix(rnorm(100), 20, 5)
  Y <- rnorm(20)
  
  obj <- AddiVortes(Y, X, m = 5, totalMCMCIter = 50, mcmcBurnIn = 10,
                    showProgress = FALSE)
  log <- obj$mcmcLog
  
  ad_rd_rows <- log[log$move_type %in% c("AD", "RD"), ]
  if (nrow(ad_rd_rows) > 0) {
    expect_true(all(!is.na(ad_rd_rows$variable)))
    expect_true(all(ad_rd_rows$variable >= 1 & ad_rd_rows$variable <= ncol(X)))
  }
})
