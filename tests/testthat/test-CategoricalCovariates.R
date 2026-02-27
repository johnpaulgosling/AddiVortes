# --- Tests for categorical covariate support ---

test_that("AddiVortes fits with one binary categorical covariate (d = 2)", {
  set.seed(42)
  n <- 30
  # Two numeric covariates plus one binary categorical (d = 2: "A", "B")
  x <- data.frame(
    x1 = rnorm(n),
    x2 = runif(n),
    cat2 = sample(c("A", "B"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x, m = 5, totalMCMCIter = 30, mcmcBurnIn = 10,
                    showProgress = FALSE)

  # Model should fit without error and produce finite RMSE
  expect_true(is.finite(fit$inSampleRmse))

  # catEncoding should be stored and reference the correct column
  expect_false(is.null(fit$catEncoding))
  expect_equal(fit$catEncoding$catColIndices, 3L)

  # A categorical variable with d = 2 produces d - 1 = 1 binary column
  # Total encoded columns: 2 numeric + 1 binary = 3
  expect_equal(length(fit$xCentres), 3L)
})

test_that("AddiVortes fits with one 3-level categorical covariate (d = 3)", {
  set.seed(7)
  n <- 40
  # Two numeric covariates plus one categorical with 3 levels (d = 3)
  x <- data.frame(
    x1 = rnorm(n),
    x2 = runif(n),
    cat3 = sample(c("low", "mid", "high"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x, m = 5, totalMCMCIter = 30, mcmcBurnIn = 10,
                    showProgress = FALSE)

  expect_true(is.finite(fit$inSampleRmse))
  expect_false(is.null(fit$catEncoding))

  # A categorical variable with d = 3 produces d - 1 = 2 binary columns
  # Total encoded columns: 2 numeric + 2 binary = 4
  expect_equal(length(fit$xCentres), 4L)
})

test_that("AddiVortes fits with both a d=2 and a d=3 categorical covariate", {
  set.seed(123)
  n <- 50
  x <- data.frame(
    x1    = rnorm(n),
    cat2  = sample(c("yes", "no"), n, replace = TRUE),
    x2    = runif(n),
    cat3  = sample(c("red", "green", "blue"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x, m = 5, totalMCMCIter = 30, mcmcBurnIn = 10,
                    showProgress = FALSE)

  expect_true(is.finite(fit$inSampleRmse))
  expect_false(is.null(fit$catEncoding))

  # cat2 -> 1 binary column; cat3 -> 2 binary columns; 2 numeric -> 2 columns
  # Total encoded columns: 2 + 1 + 2 = 5
  expect_equal(length(fit$xCentres), 5L)

  # catColIndices: columns 2 and 4 in the original data frame
  expect_equal(sort(fit$catEncoding$catColIndices), c(2L, 4L))
})

test_that("predict works after fitting with categorical covariates", {
  set.seed(99)
  n <- 40
  x <- data.frame(
    x1   = rnorm(n),
    cat2 = sample(c("A", "B"), n, replace = TRUE),
    cat3 = sample(c("p", "q", "r"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x, m = 5, totalMCMCIter = 30, mcmcBurnIn = 10,
                    showProgress = FALSE)

  # Predict on new data with the same column structure
  x_new <- data.frame(
    x1   = rnorm(10),
    cat2 = sample(c("A", "B"), 10, replace = TRUE),
    cat3 = sample(c("p", "q", "r"), 10, replace = TRUE),
    stringsAsFactors = FALSE
  )

  preds <- predict(fit, x_new, showProgress = FALSE)
  expect_equal(length(preds), 10L)
  expect_true(all(is.finite(preds)))
})

test_that("catScaling controls the binary indicator value", {
  set.seed(55)
  n <- 30
  x <- data.frame(
    x1   = rnorm(n),
    cat2 = sample(c("A", "B"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit_default <- AddiVortes(y, x, m = 5, totalMCMCIter = 20, mcmcBurnIn = 5,
                             catScaling = 1, showProgress = FALSE)
  fit_half    <- AddiVortes(y, x, m = 5, totalMCMCIter = 20, mcmcBurnIn = 5,
                             catScaling = 0.5, showProgress = FALSE)

  # Both models should store their respective catScaling
  expect_equal(fit_default$catEncoding$catScaling, 1)
  expect_equal(fit_half$catEncoding$catScaling, 0.5)
})

test_that("factor columns are handled like character columns", {
  set.seed(77)
  n <- 30
  x <- data.frame(
    x1   = rnorm(n),
    cat2 = factor(sample(c("A", "B"), n, replace = TRUE)),
    stringsAsFactors = FALSE
  )
  y <- rnorm(n)

  fit <- AddiVortes(y, x, m = 5, totalMCMCIter = 20, mcmcBurnIn = 5,
                    showProgress = FALSE)

  expect_true(is.finite(fit$inSampleRmse))
  # factor with d = 2 levels produces 1 binary column; total = 1 numeric + 1 binary = 2
  expect_equal(length(fit$xCentres), 2L)
})
