# --- Test Suite for overall functionality ---

test_that("Simple AddiVortes fit 1", {
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

  expect_equal(round(results[[1]],3), 0.724)
  expect_equal(round(results[[2]],3), 1.052)
})

test_that("Simple AddiVortes fit 2", {
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


  expect_equal(round(results[[1]],3), 2.641)
  expect_equal(round(results[[2]],3), 3.008)
})

test_that("Simple AddiVortes fit 3", {
  set.seed(1234)
  X <- matrix(rnorm(10000),1000, 10)
  Y <- runif(1000, -1, 3)
  X_test <- matrix(rnorm(1000),100, 10)
  Y_test <- runif(100, -1, 3)

  results <- AddiVortes(Y,X,10,
             200,100,
             6,0.85,3,0.8,3,25,
             Y_test,X_test,
             IntialSigma = "Linear")

  expect_equal(round(results[[1]],3), 1.149)
  expect_equal(round(results[[2]],3), 1.094)
})
