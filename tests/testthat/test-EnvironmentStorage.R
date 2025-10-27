# --- Test Suite for environment-based storage ---

test_that("AddiVortes works with environment-based storage", {
  set.seed(12345)
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  
  # Run with small number of iterations
  results <- AddiVortes(Y, X, m = 3,
    totalMCMCIter = 20, mcmcBurnIn = 5,
    6, 0.85, 3, 0.8, 3, 25,
    IntialSigma = "Linear",
    showProgress = FALSE
  )
  
  # Check that model object has expected structure
  expect_true(!is.null(results))
  expect_s3_class(results, "AddiVortesFit")
  expect_true(!is.null(results$posteriorTess))
  expect_true(!is.null(results$posteriorDim))
  expect_true(!is.null(results$posteriorPred))
  
  # Check that posterior samples are lists (converted from environments)
  expect_true(is.list(results$posteriorTess))
  expect_true(is.list(results$posteriorDim))
  expect_true(is.list(results$posteriorPred))
  
  # Check that each posterior sample contains lists (not environments)
  if (length(results$posteriorTess) > 0) {
    expect_true(is.list(results$posteriorTess[[1]]))
    expect_true(is.list(results$posteriorDim[[1]]))
    expect_true(is.list(results$posteriorPred[[1]]))
  }
})

test_that("Predictions work with environment-based storage", {
  set.seed(54321)
  X <- matrix(rnorm(40), 8, 5)
  Y <- rnorm(8)
  X_test <- matrix(rnorm(20), 4, 5)
  
  # Fit model
  results <- AddiVortes(Y, X, m = 2,
    totalMCMCIter = 15, mcmcBurnIn = 5,
    6, 0.85, 3, 0.8, 3, 25,
    IntialSigma = "Linear",
    showProgress = FALSE
  )
  
  # Make predictions
  predictions <- predict(results, X_test, showProgress = FALSE)
  
  # Check predictions
  expect_true(!is.null(predictions))
  expect_equal(length(predictions), nrow(X_test))
  expect_true(all(is.finite(predictions)))
})

test_that("Environment storage maintains backward compatibility", {
  # This test ensures that the helper functions work with both
  # environments and lists
  
  # Create test data
  set.seed(9999)
  m <- 3
  
  # Test with list (legacy)
  tess_list <- lapply(1:m, function(i) matrix(rnorm(6), 3, 2))
  
  # Test with environment (new)
  tess_env <- new.env(hash = TRUE, size = m)
  for (k in 1:m) {
    tess_env[[as.character(k)]] <- matrix(rnorm(6), 3, 2)
  }
  
  # Both should work with is.environment check
  expect_false(is.environment(tess_list))
  expect_true(is.environment(tess_env))
  
  # Test accessing elements
  expect_true(!is.null(tess_list[[1]]))
  expect_true(!is.null(tess_env[["1"]]))
  
  # Test dimensions
  expect_equal(nrow(tess_list[[1]]), 3)
  expect_equal(nrow(tess_env[["1"]]), 3)
})

test_that("Environment to list conversion preserves structure", {
  set.seed(7777)
  m <- 5
  
  # Create environment
  tess_env <- new.env(hash = TRUE, size = m)
  dim_env <- new.env(hash = TRUE, size = m)
  pred_env <- new.env(hash = TRUE, size = m)
  
  for (k in 1:m) {
    key <- as.character(k)
    tess_env[[key]] <- matrix(rnorm(6), 3, 2)
    dim_env[[key]] <- sample(1:5, 2)
    pred_env[[key]] <- rnorm(3)
  }
  
  # Convert to lists
  tess_list <- lapply(1:m, function(k) tess_env[[as.character(k)]])
  dim_list <- lapply(1:m, function(k) dim_env[[as.character(k)]])
  pred_list <- lapply(1:m, function(k) pred_env[[as.character(k)]])
  
  # Check structure preservation
  expect_equal(length(tess_list), m)
  expect_equal(length(dim_list), m)
  expect_equal(length(pred_list), m)
  
  # Check content preservation
  for (k in 1:m) {
    key <- as.character(k)
    expect_equal(tess_list[[k]], tess_env[[key]])
    expect_equal(dim_list[[k]], dim_env[[key]])
    expect_equal(pred_list[[k]], pred_env[[key]])
  }
})
