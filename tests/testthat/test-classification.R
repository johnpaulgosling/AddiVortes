# Tests for classification: task detection, binary and multinomial probit fits.

test_that("infer_response_type_internal detects regression and classification", {
  expect_equal(infer_response_type_internal(rnorm(10))$task, "regression")
  expect_equal(infer_response_type_internal(c(0, 1, 0, 1))$task, "binary")
  expect_equal(infer_response_type_internal(c(0L, 1L, 0L))$classLevels, c("0", "1"))
  expect_equal(infer_response_type_internal(c(0, 1, 2))$task, "regression")
  expect_equal(infer_response_type_internal(factor(c("a", "b", "a")))$task, "binary")
  expect_equal(infer_response_type_internal(factor(c("a", "b", "c")))$task, "multinomial")
  expect_equal(infer_response_type_internal(c("yes", "no", "yes"))$task, "binary")
  expect_equal(infer_response_type_internal(c(TRUE, FALSE, TRUE))$task, "binary")
  expect_equal(infer_response_type_internal(c(TRUE, FALSE))$nLatents, 1L)
  expect_equal(infer_response_type_internal(factor(c("a", "b", "c")))$nLatents, 2L)
})

test_that("infer_response_type_internal rejects invalid responses", {
  expect_error(infer_response_type_internal(c(1, NA, 0)), "missing")
  expect_error(infer_response_type_internal(factor("only")), "at least two")
  expect_error(infer_response_type_internal(character(0)), "length")
})

test_that("AddiVortes stores the detected task on the fitted object", {
  skip_on_cran()
  withr::local_seed(11)
  x <- matrix(rnorm(80), 20, 4)
  fit_reg <- AddiVortes(rnorm(20), x, m = 3, totalMCMCIter = 20,
                        mcmcBurnIn = 5, showProgress = FALSE)
  expect_equal(fit_reg$task, "regression")

  y_bin <- as.integer(x[, 1] > 0)
  fit_bin <- AddiVortes(y_bin, x, m = 3, totalMCMCIter = 20,
                        mcmcBurnIn = 5, showProgress = FALSE)
  expect_equal(fit_bin$task, "binary")
  expect_equal(fit_bin$classLevels, c("0", "1"))
  expect_equal(fit_bin$nLatents, 1L)
})

test_that("binary classification probabilities are in [0, 1] and beat chance", {
  skip_on_cran()
  withr::local_seed(42)
  n <- 80
  x <- matrix(runif(n * 2), n, 2)
  y <- as.integer(x[, 1] + x[, 2] > 1)
  fit <- AddiVortes(y, x, m = 8, totalMCMCIter = 80, mcmcBurnIn = 20,
                    showProgress = FALSE)
  expect_equal(fit$task, "binary")
  p <- predict(fit, x, type = "response", showProgress = FALSE)
  expect_length(p, n)
  expect_true(all(p >= 0 & p <= 1))
  cls <- predict(fit, x, type = "class", showProgress = FALSE)
  expect_s3_class(cls, "factor")
  expect_equal(levels(cls), c("0", "1"))
  expect_gt(mean(as.character(cls) == as.character(y)), 0.7)
  expect_gt(fit$inSampleAccuracy, 0.7)
  expect_true(is.finite(fit$inSampleBrier))
  q <- predict(fit, x, type = "quantile", quantiles = c(0.1, 0.9),
               showProgress = FALSE)
  expect_equal(dim(q), c(n, 2))
  expect_true(all(q[, 1] <= q[, 2]))
  link <- predict(fit, x, type = "link", showProgress = FALSE)
  expect_length(link, n)
})

test_that("binary classification works with a two-level factor", {
  skip_on_cran()
  withr::local_seed(7)
  n <- 60
  x <- matrix(rnorm(n * 3), n, 3)
  y <- factor(ifelse(x[, 1] > 0, "yes", "no"), levels = c("no", "yes"))
  fit <- AddiVortes(y, x, m = 6, totalMCMCIter = 60, mcmcBurnIn = 15,
                    showProgress = FALSE)
  expect_equal(fit$task, "binary")
  expect_equal(fit$classLevels, c("no", "yes"))
  cls <- predict(fit, x, type = "class", showProgress = FALSE)
  expect_equal(levels(cls), c("no", "yes"))
  expect_gt(mean(cls == y), 0.65)
})

test_that("multinomial classification returns probabilities that sum to 1", {
  skip_on_cran()
  withr::local_seed(99)
  n <- 75
  x <- matrix(rnorm(n * 2), n, 2)
  y <- factor(ifelse(x[, 1] < -0.4, "A", ifelse(x[, 1] > 0.4, "C", "B")))
  fit <- AddiVortes(y, x, m = 5, totalMCMCIter = 50, mcmcBurnIn = 15,
                    showProgress = FALSE)
  expect_equal(fit$task, "multinomial")
  expect_equal(fit$nLatents, 2L)
  expect_equal(fit$mPerLatent, 5L)
  expect_equal(length(fit$posteriorTess[[1]]), 10L)
  p <- predict(fit, x, type = "response", showProgress = FALSE)
  expect_equal(ncol(p), 3)
  expect_equal(colnames(p), levels(y))
  expect_equal(rowSums(p), rep(1, n), tolerance = 1e-8)
  expect_true(all(p >= 0 & p <= 1))
  cls <- predict(fit, x, type = "class", showProgress = FALSE)
  expect_equal(levels(cls), levels(y))
  expect_gt(mean(cls == y), 1 / 3)
  expect_gt(fit$inSampleAccuracy, 1 / 3)
  link <- predict(fit, x, type = "link", showProgress = FALSE)
  expect_equal(dim(link), c(n, 2))
})

test_that("type = 'class' errors on a regression fit", {
  skip_on_cran()
  withr::local_seed(3)
  x <- matrix(rnorm(40), 10, 4)
  fit <- AddiVortes(rnorm(10), x, m = 3, totalMCMCIter = 15,
                    mcmcBurnIn = 5, showProgress = FALSE)
  expect_error(
    predict(fit, x, type = "class", showProgress = FALSE),
    "only valid for classification"
  )
})

test_that("prediction intervals are rejected for classification", {
  skip_on_cran()
  withr::local_seed(4)
  x <- matrix(rnorm(40), 10, 4)
  y <- as.integer(x[, 1] > 0)
  fit <- AddiVortes(y, x, m = 3, totalMCMCIter = 15, mcmcBurnIn = 5,
                    showProgress = FALSE)
  expect_error(
    predict(fit, x, type = "quantile", interval = "prediction",
            showProgress = FALSE),
    "Prediction intervals are not used"
  )
})

test_that("print mentions classification for a binary fit", {
  skip_on_cran()
  withr::local_seed(5)
  x <- matrix(rnorm(40), 10, 4)
  y <- as.integer(x[, 1] > 0)
  fit <- AddiVortes(y, x, m = 3, totalMCMCIter = 15, mcmcBurnIn = 5,
                    showProgress = FALSE)
  expect_output(print(fit), "Binary Classification")
  expect_output(print(fit), "In-sample accuracy")
})
