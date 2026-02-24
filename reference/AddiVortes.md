# AddiVortes

The AddiVortes function is a Bayesian nonparametric regression model
that uses a tessellation to model the relationship between the
covariates and the output values. The model uses a backfitting algorithm
to sample from the posterior distribution of the output values for each
tessellation. The function returns the RMSE value for the test samples.

For spherical data, it is assumed that the final spherical dimension is
the polar angle: i.e. that with range 0, 2\*pi.

## Usage

``` r
AddiVortes(
  y,
  x,
  m = 200,
  totalMCMCIter = 1200,
  mcmcBurnIn = 200,
  nu = 6,
  q = 0.85,
  k = 3,
  sd = 0.8,
  Omega = min(3, ncol(x)),
  LambdaRate = 25,
  InitialSigma = "Linear",
  thinning = 1,
  metric = "E",
  showProgress = interactive()
)
```

## Arguments

- y:

  A vector of the output values.

- x:

  A matrix of the covariates.

- m:

  The number of tessellations.

- totalMCMCIter:

  The number of iterations.

- mcmcBurnIn:

  The number of burn in iterations.

- nu:

  The degrees of freedom.

- q:

  The quantile.

- k:

  The number of centres.

- sd:

  The standard deviation used in centre proposals.

- Omega:

  Omega/(number of covariates) is the prior probability of adding a
  dimension.

- LambdaRate:

  The rate of the Poisson distribution for the number of centres.

- InitialSigma:

  The method used to calculate the initial variance.

- thinning:

  The thinning rate.

- metric:

  Either "E" (Euclidean, default) or "S" (Spherical).

- showProgress:

  Logical; if TRUE, progress bars and messages are shown during fitting.

## Value

An AddiVortes object containing the posterior samples of the
tessellations, dimensions and predictions.

## Examples

``` r
# \donttest{
# Simple example with simulated data
set.seed(123)
x <- matrix(rnorm(50), 10, 5)
y <- rnorm(10)
# Fit model with reduced iterations for quick example
fit <- AddiVortes(y, x, m = 5, totalMCMCIter = 50, mcmcBurnIn = 10)
# }
```
