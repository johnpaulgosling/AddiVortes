# --- Test Suite for prediction intervals ---

test_that("Prediction intervals are wider than confidence intervals", {
  set.seed(12345)
  X <- matrix(rnorm(500), 100, 5)
  Y <- rnorm(100, 0, 2)
  X_test <- matrix(rnorm(50), 10, 5)
  
  results <- AddiVortes(Y, X, m = 5,
                       totalMCMCIter = 100, mcmcBurnIn = 20,
                       showProgress = FALSE)
  
  # Get confidence intervals (default)
  pred_conf <- predict(results, X_test, type = "quantile",
                      interval = "confidence",
                      quantiles = c(0.025, 0.975),
                      showProgress = FALSE)
  
  # Get prediction intervals
  pred_pred <- predict(results, X_test, type = "quantile",
                      interval = "prediction",
                      quantiles = c(0.025, 0.975),
                      showProgress = FALSE)
  
  # Calculate interval widths
  conf_width <- pred_conf[, 2] - pred_conf[, 1]
  pred_width <- pred_pred[, 2] - pred_pred[, 1]
  
  # Prediction intervals should generally be wider than confidence intervals
  # (allowing for some statistical variation due to sampling)
  expect_true(mean(pred_width) > mean(conf_width))
})

test_that("Confidence interval is default", {
  set.seed(54321)
  X <- matrix(rnorm(100), 20, 5)
  Y <- rnorm(20)
  X_test <- matrix(rnorm(25), 5, 5)
  
  results <- AddiVortes(Y, X, m = 3,
                       totalMCMCIter = 50, mcmcBurnIn = 10,
                       showProgress = FALSE)
  
  # Get quantiles without specifying interval (should default to confidence)
  pred_default <- predict(results, X_test, type = "quantile",
                         quantiles = c(0.1, 0.9),
                         showProgress = FALSE)
  
  # Get quantiles with explicit confidence interval
  pred_conf <- predict(results, X_test, type = "quantile",
                      interval = "confidence",
                      quantiles = c(0.1, 0.9),
                      showProgress = FALSE)
  
  # They should be identical
  expect_equal(pred_default, pred_conf)
})

test_that("posteriorSigma is stored in model object", {
  set.seed(99999)
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  
  results <- AddiVortes(Y, X, m = 3,
                       totalMCMCIter = 40, mcmcBurnIn = 10,
                       showProgress = FALSE)
  
  # Check that posteriorSigma exists
  expect_true("posteriorSigma" %in% names(results))
  
  # Check that it has the correct length (thinned posterior samples)
  expected_length <- floor((40 - 10) / 1)  # (totalMCMCIter - mcmcBurnIn) / thinning
  expect_equal(length(results$posteriorSigma), expected_length)
  
  # Check that all sigma values are positive
  expect_true(all(results$posteriorSigma > 0))
})

test_that("Prediction intervals require posteriorSigma", {
  set.seed(777)
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  X_test <- matrix(rnorm(25), 5, 5)
  
  results <- AddiVortes(Y, X, m = 2,
                       totalMCMCIter = 30, mcmcBurnIn = 5,
                       showProgress = FALSE)
  
  # Manually remove posteriorSigma to simulate an old model object
  results$posteriorSigma <- NULL
  
  # Attempting to get prediction intervals should fail
  expect_error(
    predict(results, X_test, type = "quantile",
           interval = "prediction",
           quantiles = c(0.025, 0.975),
           showProgress = FALSE),
    "Prediction intervals require posterior sigma samples"
  )
})

test_that("Response type ignores interval parameter", {
  set.seed(333)
  X <- matrix(rnorm(100), 20, 5)
  Y <- rnorm(20)
  X_test <- matrix(rnorm(25), 5, 5)
  
  results <- AddiVortes(Y, X, m = 3,
                       totalMCMCIter = 50, mcmcBurnIn = 10,
                       showProgress = FALSE)
  
  # Get mean predictions with different interval settings
  pred_conf <- predict(results, X_test, type = "response",
                      interval = "confidence",
                      showProgress = FALSE)
  
  pred_pred <- predict(results, X_test, type = "response",
                      interval = "prediction",
                      showProgress = FALSE)
  
  # They should be identical (interval only affects quantiles)
  expect_equal(pred_conf, pred_pred)
})
