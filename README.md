# AddiVortes

AddiVortes is a Python package for **Bayesian Additive Voronoi Tessellation**
regression. It provides a machine-learning style estimator for non-parametric
regression, spatial modeling, uncertainty quantification, and complex function
approximation.

The package exposes a Python API backed by a C++20 extension. It supports numpy
arrays and pandas data frames, including categorical covariates via one-hot
encoding.

## Features

- Pythonic estimator API with `fit`, `predict`, `fit_predict`, `score`, and
  sklearn-style parameter access.
- Bayesian posterior samples for mean predictions and credible intervals.
- Prediction intervals that include posterior error variance.
- Numeric and categorical covariate preprocessing.
- C++20 backend for the MCMC sampler and nearest-cell assignment.

## Installation

The package name is `addivortes` and it targets Python 3.10 or newer.

After a release has been published to PyPI:

```bash
python -m pip install addivortes
```

To install the current GitHub version with pip:

```bash
python -m pip install "git+https://github.com/johnpaulgosling/AddiVortes.git"
```

For local development from a source checkout:

```bash
python -m pip install -e ".[dev]"
```

## Quick start

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

predictions = model.predict(X[:5])
intervals = model.predict(X[:5], kind="quantile", quantiles=(0.025, 0.975))
```

## Data frames and categorical covariates

```python
import pandas as pd
from addivortes import AddiVortesRegressor

X = pd.DataFrame(
    {
        "x1": [0.1, 0.4, 0.2, 0.8],
        "group": pd.Categorical(["a", "b", "a", "c"]),
    }
)
y = [0.0, 1.1, 0.2, 0.9]

model = AddiVortesRegressor(
    n_tessellations=5,
    total_mcmc_iter=50,
    burn_in=10,
    random_state=123,
)
model.fit(X, y)
predictions = model.predict(X)
```

## Development

Run the Python test suite with:

```bash
python -m pytest
```

Build a wheel locally with:

```bash
python -m build
```

## Citation

If you use AddiVortes in research, please cite:

Stone, A. and Gosling, J.P. (2025). AddiVortes: (Bayesian) additive Voronoi
tessellations. Journal of Computational and Graphical Statistics.

## License

AddiVortes is distributed under the GPL-3.0-or-later license.
