# --- Test Suite for overall functionality ---

test_that("Simple AddiVortes fit 1", {
  set.seed(111333)
  X <- matrix(rnorm(100), 10, 10)
  Y <- rnorm(10)
  X_test <- matrix(rnorm(100), 10, 10)
  Y_test <- rnorm(10)

  results <- AddiVortes(Y, X, 10,
    90, 10,
    6, 0.85, 3, 0.8, 3, 25,
    InitialSigma = "Linear"
  )

  expect_equal(round(results$inSampleRmse, 3), 0.754)
})

test_that("Simple AddiVortes fit 2", {
  set.seed(1789)
  X <- matrix(runif(500), 100, 5)
  Y <- rnorm(100, -5, 3)
  X_test <- matrix(runif(100), 20, 5)
  Y_test <- rnorm(20, -5, 3)

  results <- AddiVortes(Y, X, 5,
    150, 50,
    6, 0.85, 3, 0.8, 3, 25,
    InitialSigma = "Linear"
  )


  expect_equal(round(results$inSampleRmse, 3), 2.671)
})

test_that("Simple AddiVortes fit 3", {
  set.seed(1234)
  X <- matrix(rnorm(10000), 1000, 10)
  Y <- runif(1000, -1, 3)
  X_test <- matrix(rnorm(1000), 100, 10)
  Y_test <- runif(100, -1, 3)

  results <- AddiVortes(Y, X, 10,
    200, 100,
    6, 0.85, 3, 0.8, 3, 25,
    InitialSigma = "Linear"
  )

  expect_equal(round(results$inSampleRmse, 3), 1.148)
})

test_that("Warning when p > n", {
  set.seed(42)
  # Create data where p > n (5 observations, 10 covariates)
  X <- matrix(rnorm(50), 5, 10)
  Y <- rnorm(5)
  
  # Expect a warning to be issued
  expect_warning(
    AddiVortes(Y, X, m = 2, totalMCMCIter = 10, mcmcBurnIn = 5, showProgress = FALSE),
    "Number of covariates.*exceeds number of observations"
  )
})

test_that("No warning when p <= n", {
  set.seed(42)
  # Create data where p <= n (10 observations, 5 covariates)
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  
  # Expect no warning to be issued
  expect_no_warning(
    AddiVortes(Y, X, m = 2, totalMCMCIter = 10, mcmcBurnIn = 5, showProgress = FALSE)
  )
})

test_that("Thinning parameter works correctly", {
  set.seed(567813)
  X <- matrix(rnorm(500), 50, 10)
  Y <- runif(50, -1, 3)
  X_test <- matrix(rnorm(100), 10, 10)
  
  # Test with thinning = 3
  results <- AddiVortes(Y, X, 10,
    120, 30,
    6, 0.85, 3, 0.8, 3, 25,
    InitialSigma = "Linear",
    thinning = 3,
    showProgress = FALSE
  )
  
  # Check that posterior samples are correctly thinned
  # Expected number of samples: floor((totalMCMCIter - mcmcBurnIn) / thinning)
  expected_samples <- floor((120 - 30) / 3)
  expect_equal(length(results$posteriorSigma), expected_samples)
  
  # Verify predictions work with thinning
  preds <- predict(results, X_test, showProgress = FALSE)
  expect_equal(length(preds), nrow(X_test))
})
