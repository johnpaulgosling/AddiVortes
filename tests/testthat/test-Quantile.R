# --- Test Suite for quantile functionality ---

test_that("AddiVortes fit 1 with quantile check", {
  set.seed(111333)
  X <- matrix(rnorm(100),10,10)
  Y <- rnorm(10)
  X_test <- matrix(rnorm(100),10,10)
  Y_test <- rnorm(10)
  
  results <- AddiVortes(Y,X,10,
                        90,10,
                        6,0.85,3,0.8,3,25,
                        Y_test,X_test,
                        IntialSigma = "Linear")
  predictions <- predict(results,
                         X_test,
                         "quantile",
                         c(0.5,0.75))
  
  # Check ordering of all predictions
  expect_true(all(predictions[,1] <= predictions[,2]))
  expect_equal(round(as.numeric(predictions[1,2]),3),
               0.380)
  expect_equal(round(as.numeric(predictions[10,1]),3),
               -0.098)
})

test_that("AddiVortes fit 2 with reverse quantile check", {
  set.seed(1789)
  X <- matrix(runif(500),100, 5)
  Y <- rnorm(100, -5, 3)
  X_test <- matrix(runif(100),20, 5)
  Y_test <- rnorm(20, -5, 3)
  
  results <- AddiVortes(Y,X,5,
                        150,50,
                        6,0.85,3,0.8,3,25,
                        Y_test,X_test,
                        IntialSigma = "Linear")
  predictions <- predict(results,
                         X_test,
                         "quantile",
                         c(0.95,0.75))
  
  # Check ordering of all predictions
  expect_true(all(predictions[,1] >= predictions[,2]))
})

test_that("AddiVortes fit 3 with quantile count", {
  set.seed(1234)
  X <- matrix(rnorm(10000),1000, 10)
  Y <- runif(1000, -1, 3)
  X_test <- matrix(rnorm(1000),100, 10)
  Y_test <- runif(100, -1, 3)
  
  results <- AddiVortes(Y,X,10,
                        100,30,
                        6,0.85,3,0.8,3,25,
                        Y_test,X_test,
                        IntialSigma = "Linear")
  predictions <- predict(results,
                         X_test,
                         "quantile",
                         c(0.01,0.2,0.3,0.4,0.95,0.96,0.97))
  
  expect_equal(ncol(predictions), 7)
})
