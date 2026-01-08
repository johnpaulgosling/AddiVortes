# --- Test Suite for small number of covariates ---
# This test file checks that AddiVortes works correctly with 1, 2, or 3 covariates
# and produces finite, valid results without errors or infinite loops.

test_that("AddiVortes works with 1 covariate", {
  set.seed(42)
  # Create data with only 1 covariate
  X <- matrix(rnorm(20), 20, 1)
  Y <- rnorm(20)
  
  # Fit model - should not crash or produce errors
  expect_no_error(
    results <- AddiVortes(Y, X, 
                         m = 5,
                         Omega = 1,
                         totalMCMCIter = 50, 
                         mcmcBurnIn = 10,
                         showProgress = FALSE)
  )
  
  # Check that results are valid
  expect_true(is.finite(results$inSampleRmse))
  expect_true(results$inSampleRmse >= 0)
  
  # Check that posterior samples exist
  expect_true(length(results$posteriorSigma) > 0)
  expect_true(all(is.finite(results$posteriorSigma)))
  expect_true(all(results$posteriorSigma > 0))
})

test_that("AddiVortes works with 2 covariates", {
  set.seed(123)
  # Create data with only 2 covariates
  X <- matrix(rnorm(40), 20, 2)
  Y <- rnorm(20)
  
  # Fit model - should not crash or produce errors
  expect_no_error(
    results <- AddiVortes(Y, X, 
                         m = 5,
                         Omega = 1,
                         totalMCMCIter = 50, 
                         mcmcBurnIn = 10,
                         showProgress = FALSE)
  )
  
  # Check that results are valid
  expect_true(is.finite(results$inSampleRmse))
  expect_true(results$inSampleRmse >= 0)
  
  # Check that posterior samples exist
  expect_true(length(results$posteriorSigma) > 0)
  expect_true(all(is.finite(results$posteriorSigma)))
  expect_true(all(results$posteriorSigma > 0))
})

test_that("AddiVortes works with 3 covariates", {
  set.seed(456)
  # Create data with only 3 covariates
  X <- matrix(rnorm(60), 20, 3)
  Y <- rnorm(20)
  
  # Fit model - should not crash or produce errors
  expect_no_error(
    results <- AddiVortes(Y, X, 
                         m = 5,
                         Omega = 1,
                         totalMCMCIter = 50, 
                         mcmcBurnIn = 10,
                         showProgress = FALSE)
  )
  
  # Check that results are valid
  expect_true(is.finite(results$inSampleRmse))
  expect_true(results$inSampleRmse >= 0)
  
  # Check that posterior samples exist
  expect_true(length(results$posteriorSigma) > 0)
  expect_true(all(is.finite(results$posteriorSigma)))
  expect_true(all(results$posteriorSigma > 0))
})

test_that("Predictions work with 1 covariate", {
  set.seed(789)
  X <- matrix(rnorm(20), 20, 1)
  Y <- rnorm(20)
  X_test <- matrix(rnorm(10), 10, 1)
  
  results <- AddiVortes(Y, X, 
                       m = 5,
                       Omega = 1,
                       totalMCMCIter = 50, 
                       mcmcBurnIn = 10,
                       showProgress = FALSE)
  
  # Test predictions
  expect_no_error(
    preds <- predict(results, X_test, showProgress = FALSE)
  )
  
  expect_equal(length(preds), nrow(X_test))
  expect_true(all(is.finite(preds)))
})

test_that("Predictions work with 2 covariates", {
  set.seed(101)
  X <- matrix(rnorm(40), 20, 2)
  Y <- rnorm(20)
  X_test <- matrix(rnorm(20), 10, 2)
  
  results <- AddiVortes(Y, X, 
                       m = 5, 
                       Omega = 1,
                       totalMCMCIter = 50, 
                       mcmcBurnIn = 10,
                       showProgress = FALSE)
  
  # Test predictions
  expect_no_error(
    preds <- predict(results, X_test, showProgress = FALSE)
  )
  
  expect_equal(length(preds), nrow(X_test))
  expect_true(all(is.finite(preds)))
})

test_that("Predictions work with 3 covariates", {
  set.seed(202)
  X <- matrix(rnorm(60), 20, 3)
  Y <- rnorm(20)
  X_test <- matrix(rnorm(30), 10, 3)
  
  results <- AddiVortes(Y, X, 
                       m = 5, 
                       Omega = 1,
                       totalMCMCIter = 50, 
                       mcmcBurnIn = 10,
                       showProgress = FALSE)
  
  # Test predictions
  expect_no_error(
    preds <- predict(results, X_test, showProgress = FALSE)
  )
  
  expect_equal(length(preds), nrow(X_test))
  expect_true(all(is.finite(preds)))
})

test_that("No NaN or Inf values in results with small covariates", {
  set.seed(303)
  # Test with 1 covariate - most restrictive case
  X1 <- matrix(rnorm(30), 30, 1)
  Y1 <- rnorm(30)
  
  # Run with more iterations to exercise more of the MCMC
  expect_no_error(
    results1 <- AddiVortes(Y1, X1, 
                          m = 10, 
                          Omega = 1,
                          totalMCMCIter = 100, 
                          mcmcBurnIn = 20,
                          showProgress = FALSE)
  )
  
  # Verify all sigma values are finite and positive
  expect_true(all(is.finite(results1$posteriorSigma)))
  expect_true(all(results1$posteriorSigma > 0))
  
  # Test with 2 covariates
  X2 <- matrix(rnorm(60), 30, 2)
  Y2 <- rnorm(30)
  
  expect_no_error(
    results2 <- AddiVortes(Y2, X2, 
                          m = 10, 
                          Omega = 1,
                          totalMCMCIter = 100, 
                          mcmcBurnIn = 20,
                          showProgress = FALSE)
  )
  
  expect_true(all(is.finite(results2$posteriorSigma)))
  expect_true(all(results2$posteriorSigma > 0))
  
  # Test with 3 covariates
  X3 <- matrix(rnorm(90), 30, 3)
  Y3 <- rnorm(30)
  
  expect_no_error(
    results3 <- AddiVortes(Y3, X3, 
                          m = 10, 
                          Omega = 1,
                          totalMCMCIter = 100, 
                          mcmcBurnIn = 20,
                          showProgress = FALSE)
  )
  
  expect_true(all(is.finite(results3$posteriorSigma)))
  expect_true(all(results3$posteriorSigma > 0))
})
