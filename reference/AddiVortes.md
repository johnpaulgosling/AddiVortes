# Fit an AddiVortes regression or classification model

The AddiVortes model is a Bayesian nonparametric model that uses
additive Voronoi tessellations to relate covariates to a response. For a
numeric response the model is Gaussian regression. For a classification
response it uses a probit link with Albert-Chib latent variables
(binary) or independent multinomial probit latents (three or more
classes). The task is chosen automatically from `y`.

The function can handle multiple types of covariates, including
continuous, spherical and categorical. Categorical covariates are
automatically detected. By default (`cat.onehot = TRUE`) they are
one-hot encoded, with the first level of each categorical variable used
as the reference category; the `catScaling` parameter then controls the
weight of categorical differences in distance calculations. Setting
`cat.onehot = FALSE` instead keeps each categorical covariate as a
single integer-coded column and uses Eskin distance (Eskin et al.,
2002). For spherical covariates, the function assumes that the final
spherical dimension corresponds to the polar angle, which has a range of
0 to 2\*pi. The `metric` parameter can be used to specify the type of
each covariate (Euclidean, Spherical, or Categorical), and the `members`
parameter can indicate membership of covariates into different subspaces
when using multiple spheres in covariate space.

## Usage

``` r
AddiVortes(
  y,
  x,
  m = 200,
  totalMCMCIter = 2000,
  mcmcBurnIn = 500,
  nu = 6,
  q = 0.85,
  k = 3,
  sd = 0.8,
  Omega = min(3, ncol(x)),
  LambdaRate = 5,
  InitialSigma = "Linear",
  thinning = 1,
  metric = "E",
  members = NULL,
  catScaling = 1,
  cat.onehot = TRUE,
  showProgress = interactive()
)
```

## Arguments

- y:

  A vector of response values. Numeric `y` is treated as regression,
  except when it has exactly two unique values in `{0, 1}`, which is
  binary classification. Factor, character and logical vectors are
  treated as classification: two levels give a binary probit model and
  three or more levels give a multinomial probit model. The first factor
  level (or 0 for numeric 0/1 responses) is the reference class. Missing
  values are not allowed.

- x:

  A matrix or data frame of the covariates. Character and factor columns
  are treated as categorical variables and automatically converted to
  d-1 binary indicator variables via one-hot encoding (with the first
  level as reference).

- m:

  The number of tessellations. For multinomial classification this is
  the number of tessellations **per latent dimension** (there are
  \\K-1\\ latents for \\K\\ classes).

- totalMCMCIter:

  The number of MCMC iterations. Default `2000`.

- mcmcBurnIn:

  The number of burn-in iterations. Default `500`.

- nu:

  The degrees of freedom for the inverse-gamma prior on the residual
  variance. Ignored for classification, where the latent residual
  variance is fixed at 1.

- q:

  The quantile used to set the inverse-gamma prior on the residual
  variance. Ignored for classification.

- k:

  Prior scale for tessellation output values. For regression,
  \\\sigma\_\mu = 0.5/(k\sqrt{m})\\ on the scaled response. For
  classification, \\\sigma\_\mu = 3/(k\sqrt{m})\\ on the latent probit
  scale.

- sd:

  The standard deviation used in centre proposals.

- Omega:

  Omega/(number of covariates) is the prior probability of adding a
  dimension.

- LambdaRate:

  The rate of the Poisson distribution for the number of centres.

- InitialSigma:

  The method used to calculate the initial residual variance for
  regression (`"Linear"` or `"Naive"`). Ignored for classification.

- thinning:

  The thinning rate.

- metric:

  Either "E" (Euclidean, default), "S" (Spherical), or "C"
  (Categorical).

- members:

  If needed, indicates membership of covariates into different subspaces
  (needed if using multiple spheres in covariate space). Default `NULL`.

- catScaling:

  Numeric scalar controlling the scale of binary indicator variables
  created from categorical covariates. Each binary indicator takes
  values 0 (reference level) or `catScaling` (non-reference level). The
  default value of 1 matches the range of continuous covariates, which
  are normalised to `[-0.5, 0.5]` (range = 1) during fitting, so
  categorical differences receive comparable weight to continuous
  differences in the distance calculations. Increase above 1 to give
  categorical differences more weight; decrease below 1 to give them
  less weight. Binary indicator columns are named `<colname>_<level>`
  (e.g. a column `grp` with levels `"A"`, `"B"`, `"C"` produces columns
  `grp_B` and `grp_C`, with `"A"` as the reference level).

- cat.onehot:

  Should categorical covariates be one-hot encoded? Default `TRUE`. When
  `TRUE`, each categorical covariate with *d* levels is expanded to *d*
  − 1 binary indicators and distances are Euclidean (weighted by
  `catScaling`). When `FALSE`, categories are kept as a single
  integer-coded column and mismatches use Eskin distance (Eskin et al.,
  2002), with squared cost \\2 / d^2\\ when levels differ and 0 when
  they match; `catScaling` is then ignored. See the categorical
  covariates vignette for a comparison and guidance on which to use.

- showProgress:

  Logical; if TRUE, a progress bar is shown during fitting.

## Value

An AddiVortes object containing the posterior samples of the
tessellations, dimensions and predictions, plus per-iteration trace
statistics used by
[`traceplots.AddiVortes()`](https://johnpaulgosling.github.io/AddiVortes/reference/traceplots.AddiVortes.md).
Classification fits also store `task`, `classLevels`, `nLatents` and
in-sample accuracy.

## References

Stone, A. and Gosling, J.P. (2025). AddiVortes: (Bayesian) additive
Voronoi tessellations. *Journal of Computational and Graphical
Statistics*.

Stone, A.J., Ogundimu, E. and Gosling, J.P. (2026). Binary AddiVortes:
(Bayesian) Additive Voronoi Tessellations for Binary Classification with
an application to Predicting Home Mortgage Application Outcomes.

Albert, J.H. and Chib, S. (1993). Bayesian analysis of binary and
polychotomous response data. *Journal of the American Statistical
Association*, 88(422), 669–679.

Kindo, B.P., Wang, H. and Peña, E.A. (2016). Multinomial probit Bayesian
additive regression trees. *Stat*, 5(1), 171–181.

## Examples

``` r
# \donttest{
# Simple example with simulated data
set.seed(123)
x <- matrix(rnorm(50), 10, 5)
y <- rnorm(10)
# Fit model with reduced iterations for quick example
fit <- AddiVortes(y, x, m = 5, totalMCMCIter = 50, mcmcBurnIn = 10)

# Larger example with categorical covariates (d=2 and d=3) and a test set
set.seed(456)
n_train <- 200
n_test <- 50
x_train <- data.frame(
  x1   = rnorm(n_train),
  x2   = runif(n_train),
  grp2 = sample(c("A", "B"), n_train, replace = TRUE),
  grp3 = sample(c("low", "mid", "high"), n_train, replace = TRUE)
)
y_train <- x_train$x1 + ifelse(x_train$grp2 == "B", 1, 0) + rnorm(n_train, sd = 0.5)

fit2 <- AddiVortes(y_train, x_train,
  m = 10, totalMCMCIter = 200, mcmcBurnIn = 50,
  catScaling = 1, showProgress = FALSE
)

x_test <- data.frame(
  x1   = rnorm(n_test),
  x2   = runif(n_test),
  grp2 = sample(c("A", "B"), n_test, replace = TRUE),
  grp3 = sample(c("low", "mid", "high"), n_test, replace = TRUE)
)
y_test <- x_test$x1 + ifelse(x_test$grp2 == "B", 1, 0) + rnorm(n_test, sd = 0.5)

preds <- predict(fit2, x_test, showProgress = FALSE)
test_rmse <- sqrt(mean((y_test - preds)^2))

# Binary classification is selected automatically from a 0/1 or factor y
set.seed(789)
x_clf <- matrix(runif(80), 40, 2)
y_clf <- factor(ifelse(x_clf[, 1] + x_clf[, 2] > 1, "yes", "no"))
fit_clf <- AddiVortes(y_clf, x_clf, m = 8, totalMCMCIter = 80,
                      mcmcBurnIn = 20, showProgress = FALSE)
p_clf <- predict(fit_clf, x_clf, type = "response", showProgress = FALSE)
cls_clf <- predict(fit_clf, x_clf, type = "class", showProgress = FALSE)
# }
```
