# API reference

## `AddiVortesRegressor`

```python
from addivortes import AddiVortesRegressor
```

`AddiVortesRegressor` is the primary estimator class.

### Constructor

```python
AddiVortesRegressor(
    *,
    n_tessellations=200,
    total_mcmc_iter=1200,
    burn_in=200,
    thinning=1,
    nu=6.0,
    q=0.85,
    k=3,
    omega=None,
    lambda_rate=25.0,
    initial_sigma="linear",
    metric="euclidean",
    members=None,
    cat_scaling=1.0,
    random_state=None,
    verbose=False,
)
```

### Main methods

#### `fit(X, y)`

Fit the model to covariates `X` and response `y`.

- `X` may be a numpy array or pandas DataFrame.
- pandas object, string, and categorical columns are treated as categorical.
- `y` must be one-dimensional and numeric.

Returns the fitted estimator.

#### `predict(X, *, kind="response", quantiles=(0.025, 0.975), interval="credible", random_state=None)`

Predict from a fitted model.

- `kind="response"` returns posterior mean predictions with shape `(n_obs,)`.
- `kind="quantile"` returns posterior quantiles with shape
  `(n_obs, len(quantiles))`.
- `interval="credible"` summarizes uncertainty in the posterior mean.
- `interval="prediction"` includes Gaussian observation noise using posterior
  sigma samples.

#### `fit_predict(X, y, **predict_kwargs)`

Fit the model and return predictions for the training covariates.

#### `score(X, y)`

Return the coefficient of determination on test data.

#### `summary()`

Return a dictionary with fitted-model summary information, including:

- number of features,
- number of tessellations,
- number of posterior samples,
- in-sample RMSE,
- estimated lambda,
- average cells per tessellation.

### Fitted attributes

After `fit`, common fitted attributes include:

- `posterior_`
- `x_centres_`
- `x_ranges_`
- `y_centre_`
- `y_range_`
- `in_sample_rmse_`
- `cat_encoding_`
- `feature_names_in_`
- `n_features_in_`

## Alias

`AddiVortes` is an alias for `AddiVortesRegressor`:

```python
from addivortes import AddiVortes
```
