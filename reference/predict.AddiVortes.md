# Predict Method for AddiVortes

Predicts outcomes for new data using a fitted `AddiVortes` model object.
Regression fits return means or quantiles of the response.
Classification fits return class probabilities, class labels,
latent-scale values, or quantiles of the class probabilities.

## Usage

``` r
# S3 method for class 'AddiVortes'
predict(
  object,
  newdata,
  type = c("response", "quantile", "class", "link"),
  quantiles = c(0.025, 0.975),
  interval = c("credible", "prediction"),
  showProgress = interactive(),
  ...
)
```

## Arguments

- object:

  An object of class `AddiVortes`, typically the result of a call to
  [`AddiVortes()`](https://johnpaulgosling.github.io/AddiVortes/reference/AddiVortes.md).

- newdata:

  A matrix of covariates for the new test set. The number of columns
  must match the original training data.

- type:

  The type of prediction required. The default `"response"` gives the
  mean prediction (class probabilities for classification). `"quantile"`
  returns the quantiles specified by `quantiles`. `"class"` returns
  predicted class labels (classification only). `"link"` returns the
  latent sum of tessellations \\G(x)\\ on the model scale before any
  response unscaling.

- quantiles:

  A numeric vector of probabilities to compute for the predictions when
  `type = "quantile"`.

- interval:

  The type of interval calculation. The default `"credible"` accounts
  only for uncertainty in the mean (similar to `lm`'s confidence
  interval). The alternative `"prediction"` also includes the model's
  error variance, producing wider intervals (similar to `lm`'s
  prediction interval). Not used for classification models.

- showProgress:

  Logical; if TRUE, a progress bar is shown during prediction.

- ...:

  Further arguments passed to or from other methods (currently unused).

## Value

If `type = "response"`, a numeric vector of mean predictions for
regression or binary classification, or an \\n \times K\\ probability
matrix for multinomial classification. If `type = "quantile"`, a matrix
of quantiles (binary/regression) or a named list of such matrices
(multinomial). If `type = "class"`, a factor of predicted labels. If
`type = "link"`, the latent function \\G(x)\\ on the model scale before
response unscaling.

## Details

This function relies on the internal helper function
`applyScaling_internal` being available in the environment, which is
used by the main `AddiVortes` function.

Predictions traverse all retained draws and tessellations in a single
C++ call, avoiding repeated R/C++ boundary crossings per tessellation.

When `interval = "prediction"` and `type = "quantile"`, the function
samples additional Gaussian noise with variance equal to the sampled
sigma squared from the posterior. This accounts for the inherent
variability in individual predictions, not just uncertainty in the mean
function. The noise is added in the scaled space before unscaling
predictions. Classification uses a probit link with residual variance
fixed at 1, so prediction intervals are not defined; use
`type = "quantile"` for credible intervals on probabilities.

For regression, `"response"` unscales predictions back to the original
response units, while `"link"` returns the posterior mean of the latent
scaled function \\G(x)\\.

For binary classification, `"response"` is the posterior mean of
\\\Phi(G^{(s)}(x))\\. For multinomial classification, class
probabilities are estimated from independent \\N(G, I)\\ latents, with
the first class as the reference.

## Examples

``` r
# \donttest{
# Fit a model
set.seed(123)
X <- matrix(rnorm(100), 20, 5)
Y <- rnorm(20)
fit <- AddiVortes(Y, X, m = 5, totalMCMCIter = 50, mcmcBurnIn = 10)

# New data for prediction
X_new <- matrix(rnorm(25), 5, 5)

# Mean predictions
pred_mean <- predict(fit, X_new, type = "response")

# Credible intervals (uncertainty in mean only)
pred_conf <- predict(fit, X_new,
  type = "quantile",
  interval = "credible",
  quantiles = c(0.025, 0.975)
)

# Prediction intervals (includes error variance)
pred_pred <- predict(fit, X_new,
  type = "quantile",
  interval = "prediction",
  quantiles = c(0.025, 0.975)
)

# Prediction intervals are wider than credible intervals
mean(pred_pred[, 2] - pred_pred[, 1]) > mean(pred_conf[, 2] - pred_conf[, 1])
#> [1] TRUE
# }
```
