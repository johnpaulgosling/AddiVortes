# Predict Method for AddiVortesFit

Predicts outcomes for new data using a fitted `AddiVortesFit` model
object. It can return mean predictions, quantiles and optionally
calculate the Root Mean Squared Error (RMSE) if true outcomes are
provided.

## Usage

``` r
# S3 method for class 'AddiVortesFit'
predict(
  object,
  newdata,
  type = c("response", "quantile"),
  quantiles = c(0.025, 0.975),
  interval = c("credible", "prediction"),
  showProgress = TRUE,
  parallel = TRUE,
  cores = NULL,
  ...
)
```

## Arguments

- object:

  An object of class `AddiVortesFit`, typically the result of a call to
  [`AddiVortes()`](https://johnpaulgosling.github.io/AddiVortes/reference/AddiVortes.md).

- newdata:

  A matrix of covariates for the new test set. The number of columns
  must match the original training data.

- type:

  The type of prediction required. The default `"response"` gives the
  mean prediction. The alternative `"quantile"` returns the quantiles
  specified by the `quantiles` argument.

- quantiles:

  A numeric vector of probabilities with values in \[0, 1\] to compute
  for the predictions when `type = "quantile"`.

- interval:

  The type of interval calculation. The default `"credible"` accounts
  only for uncertainty in the mean (similar to lm's confidence
  interval). The alternative `"prediction"` also includes the model's
  error variance, producing wider intervals (similar to lm's prediction
  interval).

- showProgress:

  Logical; if TRUE (default), a progress bar is shown during prediction.

- parallel:

  Logical; if TRUE (default), predictions are computed in parallel.

- cores:

  The number of CPU cores to use for parallel processing. If NULL
  (default), it defaults to one less than the total number of available
  cores.

- ...:

  Further arguments passed to or from other methods (currently unused).

## Value

If `type = "response"`, a numeric vector of mean predictions. If
`type = "quantile"`, a matrix where each row corresponds to an
observation in `newdata` and each column to a quantile.

## Details

This function relies on the internal helper function
`applyScaling_internal` being available in the environment, which is
used by the main `AddiVortes` function.

When `interval = "prediction"` and `type = "quantile"`, the function
samples additional Gaussian noise with variance equal to the sampled
sigma squared from the posterior. This accounts for the inherent
variability in individual predictions, not just uncertainty in the mean
function. The noise is added in the scaled space before unscaling
predictions.

## Examples

``` r
# \donttest{
# Fit a model
set.seed(123)
X <- matrix(rnorm(100), 20, 5)
Y <- rnorm(20)
fit <- AddiVortes(Y, X, m = 5, totalMCMCIter = 50, mcmcBurnIn = 10)
#> Fitting AddiVortes model to input data...
#> Input dimensions: 20 observations, 5 covariates
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

# New data for prediction
X_new <- matrix(rnorm(25), 5, 5)

# Mean predictions
pred_mean <- predict(fit, X_new, type = "response")
#> Generating predictions for 5 observations using 40 posterior samples...
#> 
#> Prediction generation completed.
#> 

# Credible intervals (uncertainty in mean only)
pred_conf <- predict(fit, X_new, type = "quantile", 
                    interval = "credible",
                    quantiles = c(0.025, 0.975))
#> Generating predictions for 5 observations using 40 posterior samples...
#> 
#> Prediction generation completed.
#> 

# Prediction intervals (includes error variance)
pred_pred <- predict(fit, X_new, type = "quantile",
                    interval = "prediction",
                    quantiles = c(0.025, 0.975))
#> Generating predictions for 5 observations using 40 posterior samples...
#> 
#> Prediction generation completed.
#> 

# Prediction intervals are wider than credible intervals
mean(pred_pred[, 2] - pred_pred[, 1]) > mean(pred_conf[, 2] - pred_conf[, 1])
#> [1] TRUE
# }
```
