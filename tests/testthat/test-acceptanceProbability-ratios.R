log_likelihood_ratio_for_test <- function(R_old, n_old, R_new, n_new,
                                          sigma_sq, sigma_mu) {
  0.5 * (
    log(prod(n_old * sigma_mu + sigma_sq)) -
      log(prod(n_new * sigma_mu + sigma_sq))
  ) +
    (sigma_mu / (2 * sigma_sq)) *
      (
        sum((R_new^2) / (n_new * sigma_mu + sigma_sq)) -
          sum((R_old^2) / (n_old * sigma_mu + sigma_sq))
      )
}

add_dimension_probability_for_test <- function(d_current, p) {
  if (d_current >= p) {
    return(0)
  }
  if (d_current == 1) {
    return(0.4)
  }
  0.2
}

remove_dimension_probability_for_test <- function(d_current, p) {
  if (d_current <= 1) {
    return(0)
  }
  if (d_current == p) {
    return(0.4)
  }
  0.2
}

add_centre_probability_for_test <- function(c_current) {
  if (c_current == 1) {
    return(0.4)
  }
  0.2
}

remove_centre_probability_for_test <- function(c_current) {
  if (c_current > 1) {
    return(0.2)
  }
  0
}

expect_acceptance_ratio <- function(R_old, n_old, R_new, n_new,
                                    tess_rows, dim_count, sigma_sq, sigma_mu,
                                    modification, omega, lambda_rate,
                                    p, expected) {
  observed <- acceptanceProbability(
    R_ijOld = R_old,
    n_ijOld = n_old,
    R_ijNew = R_new,
    n_ijNew = n_new,
    tess_j_star = matrix(0, nrow = tess_rows, ncol = dim_count),
    dim_j_star = seq_len(dim_count),
    SigmaSquared = sigma_sq,
    Modification = modification,
    SigmaSquaredMu = sigma_mu,
    Omega = omega,
    LambdaRate = lambda_rate,
    NumCovariates = p
  )

  expect_equal(observed, expected, tolerance = 1e-12)
}

test_that("tile add and remove ratios use old counts and boundary factors", {
  sigma_sq <- 1.7
  sigma_mu <- 0.4
  lambda_rate <- 3.2
  omega <- 2
  p <- 5

  check_ac <- function(k_old, R_old, n_old, R_new, n_new) {
    log_lik <- log_likelihood_ratio_for_test(
      R_old, n_old, R_new, n_new, sigma_sq, sigma_mu
    )
    k_new <- k_old + 1
    expected <- log_lik +
      log(lambda_rate) -
      2 * log(k_new) +
      log(remove_centre_probability_for_test(k_new)) -
      log(add_centre_probability_for_test(k_old)) +
      0.5 * log(sigma_sq)

    expect_acceptance_ratio(
      R_old, n_old, R_new, n_new,
      tess_rows = k_new, dim_count = 2,
      sigma_sq = sigma_sq,
      sigma_mu = sigma_mu,
      modification = "AC",
      omega = omega,
      lambda_rate = lambda_rate,
      p = p,
      expected = expected
    )
  }

  check_rc <- function(k_old, R_old, n_old, R_new, n_new) {
    log_lik <- log_likelihood_ratio_for_test(
      R_old, n_old, R_new, n_new, sigma_sq, sigma_mu
    )
    k_new <- k_old - 1
    expected <- log_lik +
      2 * log(k_old) -
      log(lambda_rate) +
      log(add_centre_probability_for_test(k_new)) -
      log(remove_centre_probability_for_test(k_old)) -
      0.5 * log(sigma_sq)

    expect_acceptance_ratio(
      R_old, n_old, R_new, n_new,
      tess_rows = k_new, dim_count = 2,
      sigma_sq = sigma_sq,
      sigma_mu = sigma_mu,
      modification = "RC",
      omega = omega,
      lambda_rate = lambda_rate,
      p = p,
      expected = expected
    )
  }

  check_ac(
    k_old = 3,
    R_old = c(0.5, -0.2, 0.1),
    n_old = c(2, 3, 4),
    R_new = c(0.3, -0.1, 0.2, 0.4),
    n_new = c(2, 2, 3, 2)
  )
  check_ac(
    k_old = 1,
    R_old = 0.5,
    n_old = 9,
    R_new = c(0.2, 0.3),
    n_new = c(4, 5)
  )
  check_rc(
    k_old = 4,
    R_old = c(0.3, -0.1, 0.2, 0.4),
    n_old = c(2, 2, 3, 2),
    R_new = c(0.5, -0.2, 0.1),
    n_new = c(2, 3, 4)
  )
  check_rc(
    k_old = 2,
    R_old = c(0.2, 0.3),
    n_old = c(4, 5),
    R_new = 0.5,
    n_new = 9
  )
})

test_that("dimension add and remove ratios use old counts and boundary factors", {
  sigma_sq <- 1.2
  sigma_mu <- 0.3
  lambda_rate <- 2.5
  R_old <- c(0.4, -0.2)
  n_old <- c(4, 5)
  R_new <- c(0.1, 0.1)
  n_new <- c(3, 6)

  check_ad <- function(d_old, p, omega) {
    prob <- omega / p
    d_new <- d_old + 1
    log_lik <- log_likelihood_ratio_for_test(
      R_old, n_old, R_new, n_new, sigma_sq, sigma_mu
    )
    expected <- log_lik +
      2 * log(p - d_old) -
      log(d_old) -
      log(d_new) +
      log(prob) -
      log1p(-prob) +
      log(remove_dimension_probability_for_test(d_new, p)) -
      log(add_dimension_probability_for_test(d_old, p))

    expect_acceptance_ratio(
      R_old, n_old, R_new, n_new,
      tess_rows = length(n_new), dim_count = d_new,
      sigma_sq = sigma_sq,
      sigma_mu = sigma_mu,
      modification = "AD",
      omega = omega,
      lambda_rate = lambda_rate,
      p = p,
      expected = expected
    )
  }

  check_rd <- function(d_old, p, omega) {
    prob <- omega / p
    d_new <- d_old - 1
    log_lik <- log_likelihood_ratio_for_test(
      R_old, n_old, R_new, n_new, sigma_sq, sigma_mu
    )
    expected <- log_lik +
      log(d_old) +
      log(d_new) -
      2 * log(p - d_new) +
      log1p(-prob) -
      log(prob) +
      log(add_dimension_probability_for_test(d_new, p)) -
      log(remove_dimension_probability_for_test(d_old, p))

    expect_acceptance_ratio(
      R_old, n_old, R_new, n_new,
      tess_rows = length(n_new), dim_count = d_new,
      sigma_sq = sigma_sq,
      sigma_mu = sigma_mu,
      modification = "RD",
      omega = omega,
      lambda_rate = lambda_rate,
      p = p,
      expected = expected
    )
  }

  check_ad(d_old = 3, p = 6, omega = 2)
  check_ad(d_old = 1, p = 5, omega = 2)
  check_ad(d_old = 4, p = 5, omega = 2)
  check_rd(d_old = 4, p = 6, omega = 2)
  check_rd(d_old = 2, p = 5, omega = 2)
  check_rd(d_old = 5, p = 5, omega = 2)
})

test_that("dimension boundary factors do not double apply when p is two", {
  sigma_sq <- 1.2
  sigma_mu <- 0.3
  lambda_rate <- 2.5
  omega <- 1
  p <- 2
  prob <- omega / p
  R_old <- c(0.4, -0.2)
  n_old <- c(4, 5)
  R_new <- c(0.1, 0.1)
  n_new <- c(3, 6)
  log_lik <- log_likelihood_ratio_for_test(
    R_old, n_old, R_new, n_new, sigma_sq, sigma_mu
  )

  expected_ad <- log_lik +
    2 * log(1) -
    log(1) -
    log(2) +
    log(prob) -
    log1p(-prob)

  expect_acceptance_ratio(
    R_old, n_old, R_new, n_new,
    tess_rows = length(n_new), dim_count = 2,
    sigma_sq = sigma_sq,
    sigma_mu = sigma_mu,
    modification = "AD",
    omega = omega,
    lambda_rate = lambda_rate,
    p = p,
    expected = expected_ad
  )

  expected_rd <- log_lik +
    log(2) +
    log(1) -
    2 * log(1) +
    log1p(-prob) -
    log(prob)

  expect_acceptance_ratio(
    R_old, n_old, R_new, n_new,
    tess_rows = length(n_new), dim_count = 1,
    sigma_sq = sigma_sq,
    sigma_mu = sigma_mu,
    modification = "RD",
    omega = omega,
    lambda_rate = lambda_rate,
    p = p,
    expected = expected_rd
  )
})
