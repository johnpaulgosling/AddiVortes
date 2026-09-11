# Fixed-seed Friedman benchmark used to guard MCMC/predict hot-path changes.
# The classic Friedman function (Friedman, 1991) with p = 10 (five signal
# covariates and five noise covariates).

friedman_response <- function(x) {
  10 * sin(pi * x[, 1] * x[, 2]) +
    20 * (x[, 3] - 0.5)^2 +
    10 * x[, 4] +
    5 * x[, 5]
}

friedman_data <- function(n, p = 10, sigma = 1, seed) {
  withr::local_seed(seed)
  x <- matrix(runif(n * p), n, p)
  y <- friedman_response(x) + rnorm(n, sd = sigma)
  list(x = x, y = y)
}

# Compact configuration: long enough to exercise proposals and predict, short
# enough for routine testthat runs.
FRIEDMAN_BENCH <- list(
  seed_train = 20260806L,
  seed_test = 20260807L,
  n_train = 500L,
  n_test = 200L,
  p = 10L,
  sigma = 1,
  m = 100L,
  totalMCMCIter = 400L,
  mcmcBurnIn = 100L,
  thinning = 1L,
  # Locked against the pre-optimisation chain; regenerate deliberately if the
  # sampler's RNG contract changes.
  expected_in_sample_rmse = 0.454,
  expected_test_rmse = 1.431,
  rmse_digits = 3L
)

test_that("Friedman fixed-seed fit and predict stay bit-stable", {
  skip_on_cran()

  train <- friedman_data(
    FRIEDMAN_BENCH$n_train,
    FRIEDMAN_BENCH$p,
    FRIEDMAN_BENCH$sigma,
    FRIEDMAN_BENCH$seed_train
  )
  test <- friedman_data(
    FRIEDMAN_BENCH$n_test,
    FRIEDMAN_BENCH$p,
    FRIEDMAN_BENCH$sigma,
    FRIEDMAN_BENCH$seed_test
  )

  withr::local_seed(FRIEDMAN_BENCH$seed_train)
  fit_secs <- system.time({
    fit <- AddiVortes(
      train$y,
      train$x,
      m = FRIEDMAN_BENCH$m,
      totalMCMCIter = FRIEDMAN_BENCH$totalMCMCIter,
      mcmcBurnIn = FRIEDMAN_BENCH$mcmcBurnIn,
      thinning = FRIEDMAN_BENCH$thinning,
      showProgress = FALSE
    )
  })[["elapsed"]]

  expect_s3_class(fit, "AddiVortes")
  expect_equal(
    round(fit$inSampleRmse, FRIEDMAN_BENCH$rmse_digits),
    FRIEDMAN_BENCH$expected_in_sample_rmse
  )

  pred_secs <- system.time({
    preds <- predict(fit, test$x, showProgress = FALSE)
  })[["elapsed"]]
  test_rmse <- sqrt(mean((test$y - preds)^2))
  expect_equal(
    round(test_rmse, FRIEDMAN_BENCH$rmse_digits),
    FRIEDMAN_BENCH$expected_test_rmse
  )

  # Timing is informational: printed when tests are run interactively / via
  # testthat reporters that show messages.
  message(sprintf(
    "Friedman benchmark (n=%d, m=%d, iter=%d): fit=%.3fs predict=%.3fs in=%.3f test=%.3f",
    FRIEDMAN_BENCH$n_train,
    FRIEDMAN_BENCH$m,
    FRIEDMAN_BENCH$totalMCMCIter,
    fit_secs,
    pred_secs,
    fit$inSampleRmse,
    test_rmse
  ))
})

test_that("Friedman ensemble predict matches cellIndices reference", {
  skip_on_cran()

  train <- friedman_data(80L, 10L, 1, 4242L)
  withr::local_seed(4242L)
  fit <- AddiVortes(
    train$y, train$x,
    m = 10L, totalMCMCIter = 60L, mcmcBurnIn = 20L,
    showProgress = FALSE
  )
  test <- friedman_data(40L, 10L, 1, 4243L)

  preds_cpp <- predict(fit, test$x, showProgress = FALSE)

  xNewScaled <- applyScaling_internal(test$x, fit$xCentres, fit$xRanges)
  metricAug <- rep(fit$metric_red, fit$member_red)
  xNewScaled[, metricAug != 0] <- test$x[, metricAug != 0]
  nObs <- nrow(xNewScaled)
  m <- length(fit$posteriorTess[[1]])
  pred_mat <- matrix(0, nObs, length(fit$posteriorTess))
  for (s in seq_along(fit$posteriorTess)) {
    draw_pred <- numeric(nObs)
    for (j in seq_len(m)) {
      idx <- AddiVortes:::cellIndices(
        xNewScaled,
        fit$posteriorTess[[s]][[j]],
        fit$posteriorDim[[s]][[j]],
        fit$metric_red,
        fit$member_red
      )
      draw_pred <- draw_pred + fit$posteriorPred[[s]][[j]][idx]
    }
    pred_mat[, s] <- draw_pred
  }
  preds_ref <- rowMeans(pred_mat) * fit$yRange + fit$yCentre
  expect_equal(preds_cpp, preds_ref)
})
