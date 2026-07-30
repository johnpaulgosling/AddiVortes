# Structural-move prior check.
#
# With the likelihood switched off, AD/RD/AC/RC should sample the cell-count
# and covariate-count priors from section 2.3.2. The log acceptance ratios here
# must stay equivalent to those in src/addi_vortes_code.cpp
# (log_acceptance_components).

selection_probs <- function(b, d, p) {
  c(
    AD = if (d != p) (if (d == 1) 0.4 else 0.2) else 0,
    RD = if (d > 1) (if (d == p) 0.4 else 0.2) else 0,
    AC = if (b == 1) 0.4 else 0.2,
    RC = if (b > 1) 0.2 else 0,
    Change = if (d == p) 0.2 else 0.1,
    Swap = if (d < p) 0.1 else 0
  )
}

# Corrected structure ratios + Appendix B boundary factors, matching C++.
log_acc_structure <- function(move, b, d, p, omega, lambda_c) {
  nC_new <- b + (move == "AC") - (move == "RC")
  d_new <- d + (move == "AD") - (move == "RD")

  a <- switch(
    move,
    AD = {
      acc <- log(p - d_new + 1) - log(d_new - 1) + log(omega) - log(p - omega)
      if (d_new == 2) acc <- acc - log(2)
      if (d_new == p) acc <- acc + log(2)
      acc
    },
    RD = {
      acc <- log(d_new) - log(p - d_new) + log(p - omega) - log(omega)
      if (d_new == (p - 1)) acc <- acc - log(2)
      if (d_new == 1) acc <- acc + log(2)
      acc
    },
    AC = {
      acc <- log(lambda_c) - log(nC_new - 1)
      if (nC_new == 2) acc <- acc - log(2)
      acc
    },
    RC = {
      acc <- log(nC_new) - log(lambda_c)
      if (nC_new == 1) acc <- acc + log(2)
      acc
    },
    0
  )
  a
}

run_structural_chain <- function(n_iter, p, omega, lambda_c, b_max = 60,
                                 seed = 1) {
  set.seed(seed)
  b <- 1L
  d <- 1L
  bs <- integer(n_iter)
  ds <- integer(n_iter)
  for (i in seq_len(n_iter)) {
    q <- selection_probs(b, d, p)
    move <- sample(names(q), 1L, prob = q)
    ok <- switch(
      move,
      AD = d < p,
      RD = d > 1,
      AC = b < b_max,
      RC = b > 1,
      TRUE
    )
    if (ok && log(runif(1)) < log_acc_structure(move, b, d, p, omega, lambda_c)) {
      b <- b + (move == "AC") - (move == "RC")
      d <- d + (move == "AD") - (move == "RD")
    }
    bs[i] <- b
    ds[i] <- d
  }
  list(b = bs, d = ds)
}

test_that("structural moves recover cell and covariate count priors", {
  skip_on_cran()

  p <- 10
  omega <- 3
  lambda_c <- 25
  n_iter <- 80000

  chain <- run_structural_chain(n_iter, p, omega, lambda_c, seed = 42)

  # b - 1 ~ Poisson(lambda_c); d - 1 ~ Binomial(p - 1, omega / p)
  expect_lt(abs(mean(chain$b) - (1 + lambda_c)), 1.0)
  expect_lt(abs(mean(chain$d) - (1 + (p - 1) * omega / p)), 0.15)
})
