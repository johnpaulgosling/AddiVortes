# Fast smoke tests that run on CRAN — deliberately minimal to stay under time limits.
# All other model-fitting tests use skip_on_cran().

test_that("AddiVortes returns AddiVortes object", {
  withr::local_seed(1)
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  fit <- AddiVortes(Y, X,
    m = 3, totalMCMCIter = 20, mcmcBurnIn = 5,
    showProgress = FALSE
  )
  expect_s3_class(fit, "AddiVortes")
})

test_that("predict returns correct dimensions", {
  withr::local_seed(2)
  X <- matrix(rnorm(50), 10, 5)
  Y <- rnorm(10)
  fit <- AddiVortes(Y, X,
    m = 3, totalMCMCIter = 20, mcmcBurnIn = 5,
    showProgress = FALSE
  )
  X_new <- matrix(rnorm(25), 5, 5)
  preds <- predict(fit, X_new, showProgress = FALSE)
  expect_length(preds, 5)
})

test_that("binary classification smoke fit returns probabilities", {
  withr::local_seed(3)
  X <- matrix(rnorm(50), 10, 5)
  Y <- as.integer(X[, 1] > 0)
  fit <- AddiVortes(Y, X,
    m = 3, totalMCMCIter = 20, mcmcBurnIn = 5,
    showProgress = FALSE
  )
  expect_s3_class(fit, "AddiVortes")
  expect_equal(fit$task, "binary")
  preds <- predict(fit, X, type = "response", showProgress = FALSE)
  expect_length(preds, 10)
  expect_true(all(preds >= 0 & preds <= 1))
})
