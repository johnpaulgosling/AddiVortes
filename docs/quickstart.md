# Quick start

## Numpy arrays

```python
import numpy as np
from addivortes import AddiVortesRegressor

rng = np.random.default_rng(123)
X = rng.normal(size=(100, 4))
y = X[:, 0] - 0.5 * X[:, 1] + rng.normal(scale=0.2, size=100)

model = AddiVortesRegressor(
    n_tessellations=25,
    total_mcmc_iter=300,
    burn_in=100,
    random_state=123,
)
model.fit(X, y)

mean_predictions = model.predict(X[:5])
credible_intervals = model.predict(
    X[:5],
    kind="quantile",
    interval="credible",
    quantiles=(0.025, 0.975),
)
prediction_intervals = model.predict(
    X[:5],
    kind="quantile",
    interval="prediction",
    quantiles=(0.025, 0.975),
)
```

## Pandas data frames and categorical covariates

```python
import numpy as np
import pandas as pd
from addivortes import AddiVortesRegressor

rng = np.random.default_rng(456)
X = pd.DataFrame(
    {
        "x1": rng.normal(size=80),
        "group": pd.Categorical(rng.choice(["a", "b", "c"], size=80)),
    }
)
y = X["x1"].to_numpy() + (X["group"].astype(str) == "b").to_numpy(dtype=float)

model = AddiVortesRegressor(
    n_tessellations=20,
    total_mcmc_iter=250,
    burn_in=75,
    cat_scaling=1.0,
    random_state=456,
)
model.fit(X, y)

new_X = pd.DataFrame({"x1": [0.0, 1.0], "group": ["b", "unknown"]})
predictions = model.predict(new_X)
```

Unseen categorical levels at prediction time are treated as the reference level
for that categorical covariate.

## Model summary

```python
summary = model.summary()
print(summary["posterior_samples"])
print(summary["in_sample_rmse"])
```
