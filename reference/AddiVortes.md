# AddiVortes

The AddiVortes function is a Bayesian nonparametric regression model
that uses a tessellation to model the relationship between the
covariates and the output values. The model uses a backfitting algorithm
to sample from the posterior distribution of the output values for each
tessellation. The function returns the RMSE value for the test samples.

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
  Omega = 3,
  LambdaRate = 25,
  IntialSigma = "Linear",
  thinning = 1,
  showProgress = TRUE
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

  The standard deviation.

- Omega:

  Omega/(number of covariates) is the prior probability of adding a
  dimension.

- LambdaRate:

  The rate of the Poisson distribution for the number of centres.

- IntialSigma:

  The method used to calculate the initial variance.

- thinning:

  The thinning rate.

- showProgress:

  Logical; if TRUE (default), progress bars and messages are shown
  during fitting.

## Value

An AddiVortesFit object containing the posterior samples of the
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
#> Fitting AddiVortes model to input data...
#> Input dimensions: 10 observations, 5 covariates
#> Model configuration: 5 tessellations, 50 total iterations (10 burn-in)
#> 
#> Phase 1: Burn-in sampling (10 iterations)
#>   |                                                          |                                                  |   0%  |                                                          |=====                                             |  10%  |                                                          |==========                                        |  20%  |                                                          |===============                                   |  30%  |                                                          |====================                              |  40%  |                                                          |=========================                         |  50%  |                                                          |==============================                    |  60%  |                                                          |===================================               |  70%  |                                                          |========================================          |  80%  |                                                          |=============================================     |  90%  |                                                          |==================================================| 100%
#> 
#> Phase 2: Posterior sampling (40 iterations)
#>   |                                                          |                                                  |   0%  |                                                          |==                                                |   5%  |                                                          |====                                              |   8%  |                                                          |=====                                             |  10%  |                                                          |======                                            |  12%  |                                                          |========                                          |  15%  |                                                          |=========                                         |  18%  |                                                          |==========                                        |  20%  |                                                          |===========                                       |  22%  |                                                          |============                                      |  25%  |                                                          |==============                                    |  28%  |                                                          |===============                                   |  30%  |                                                          |================                                  |  32%  |                                                          |==================                                |  35%  |                                                          |===================                               |  38%  |                                                          |====================                              |  40%  |                                                          |=====================                             |  42%  |                                                          |======================                            |  45%  |                                                          |========================                          |  48%  |                                                          |=========================                         |  50%  |                                                          |==========================                        |  52%  |                                                          |============================                      |  55%  |                                                          |=============================                     |  58%  |                                                          |==============================                    |  60%  |                                                          |===============================                   |  62%  |                                                          |================================                  |  65%  |                                                          |==================================                |  68%  |                                                          |===================================               |  70%  |                                                          |====================================              |  72%  |                                                          |======================================            |  75%  |                                                          |=======================================           |  78%  |                                                          |========================================          |  80%  |                                                          |=========================================         |  82%  |                                                          |==========================================        |  85%  |                                                          |============================================      |  88%  |                                                          |=============================================     |  90%  |                                                          |==============================================    |  92%  |                                                          |================================================  |  95%  |                                                          |================================================= |  98%  |                                                          |==================================================| 100%
#> 
#> MCMC sampling completed.
#> 
# }
```
