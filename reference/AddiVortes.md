# Fit AddiVortes regression model

The AddiVortes model is a Bayesian nonparametric regression model that
uses a tessellation to model the relationship between the covariates and
the output values. The model uses a backfitting algorithm to sample from
the posterior distribution of the output values for each tessellation.
Alongside fitting details and the posterior sample, the function returns
the RMSE value for the test samples.

The function can handle multiple types of covariates, including
continuous, spherical and categorical. Categorical covariates are
automatically detected and one-hot encoded, with the first level of each
categorical variable used as the reference category. The `catScaling`
parameter allows control over the weight of categorical differences in
distance calculations. For spherical covariates, the function assumes
that the final spherical dimension corresponds to the polar angle, which
has a range of 0 to 2\*pi. The `metric` parameter can be used to specify
the type of each covariate (Euclidean, Spherical, or Categorical), and
the `members` parameter can indicate membership of covariates into
different subspaces when using multiple spheres in covariate space.

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
  members = NULL,
  catScaling = 1,
  showProgress = interactive()
)
```

## Arguments

- y:

  A vector of the output values.

- x:

  A matrix or data frame of the covariates. Character and factor columns
  are treated as categorical variables and automatically converted to
  d-1 binary indicator variables via one-hot encoding (with the first
  level as reference).

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
# }
```
