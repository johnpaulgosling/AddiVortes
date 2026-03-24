#### ---- Unit Tests for qinvgamma_internal ---

test_that("qinvgamma_internal calculates quantiles correctly", {
  # 1. Test against the known relationship with qgamma
  shape <- 5
  rate <- 2
  p <- 0.5 # The median

  # The quantile for invgamma(p) should be 1 / qgamma(1 - p, ...) for the same parameters
  # Or more directly, 1 / qgamma(p, ..., lower.tail = FALSE)
  expected_quantile <- 1 / stats::qgamma(p, shape = shape, rate = rate, lower.tail = FALSE)
  actual_quantile <- qinvgamma_internal(p, shape = shape, rate = rate)

  expect_equal(actual_quantile, expected_quantile)

  # 2. Test vectorised 'p'
  p_vec <- c(0.25, 0.5, 0.75)
  expected_vec <- 1 / stats::qgamma(p_vec, shape = shape, rate = rate, lower.tail = FALSE)
  expect_equal(qinvgamma_internal(p_vec, shape = shape, rate = rate), expected_vec)
})

test_that("qinvgamma_internal handles edge cases for p correctly", {
  # p=0 should result in a quantile of 0
  expect_equal(qinvgamma_internal(0, shape = 2, rate = 1), 0)

  # p=1 should result in a quantile of Inf
  expect_equal(qinvgamma_internal(1, shape = 2, rate = 1), Inf)

  # Test with log.p = TRUE
  expect_equal(qinvgamma_internal(log(0), shape = 2, rate = 1, log.p = TRUE), 0)
  expect_equal(qinvgamma_internal(log(1), shape = 2, rate = 1, log.p = TRUE), Inf)
})

test_that("qinvgamma_internal handles lower.tail and log.p arguments", {
  p <- 0.8
  shape <- 3
  rate <- 2

  # P(X > x) = 1 - p  is the same as P(X <= x) = p with lower.tail = FALSE
  q_lower <- qinvgamma_internal(p,
    shape, rate,
    lower.tail = TRUE
  )
  q_upper <- qinvgamma_internal(1 - p,
    shape, rate,
    lower.tail = FALSE
  )
  expect_equal(q_lower, q_upper)

  # Test with log.p
  log_p <- log(p)
  q_log <- qinvgamma_internal(log_p, shape, rate, log.p = TRUE)
  expect_equal(q_lower, q_log)
})

test_that("qinvgamma_internal handles rate and scale parameters correctly", {
  p <- 0.5
  shape <- 3

  q_from_rate <- qinvgamma_internal(p,
    shape = shape,
    rate = 4
  )
  # scale = 1/4 should give identical results to rate = 4
  q_from_scale <- qinvgamma_internal(p,
    shape = shape,
    scale = 0.25
  )

  expect_equal(
    q_from_rate,
    q_from_scale
  )
})

test_that("qinvgamma_internal validates inputs and throws errors", {
  # Test invalid 'shape'
  expect_error(qinvgamma_internal(0.5, shape = -1), "'shape' must be a single positive finite number.")

  # Test invalid 'rate' and 'scale'
  expect_error(qinvgamma_internal(0.5,
    shape = 1,
    rate = -1
  ), "'rate' must be a single positive finite number.")
  expect_error(qinvgamma_internal(0.5,
    shape = 1,
    scale = -1
  ), "'scale' must be a single positive finite number.")

  # Test invalid 'p'
  # Check if it returns NaN as stats::qgamma would
  suppressWarnings(expect_true(is.nan(qinvgamma_internal(1.1,
    shape = 1
  ))))
})
