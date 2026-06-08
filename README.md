# AddiVortes: Bayesian Additive Voronoi Tessellations

## Overview

AddiVortes implements the **Bayesian Additive Voronoi Tessellation** model for machine learning regression and non-parametric statistical modeling. This R package provides a flexible alternative to **BART (Bayesian Additive Regression Trees)**, using Voronoi tessellations instead of trees for spatial partitioning.

## Key Features

- **Machine Learning Regression**: Advanced Bayesian regression modeling for complex datasets
- **Alternative to BART**: Uses Voronoi tessellations instead of trees for more flexible spatial modeling
- **Spatial Data Analysis**: Excellent for geographic and spatial datasets
- **Non-parametric Modeling**: No assumptions about functional form
- **Bayesian Framework**: Full posterior inference with uncertainty quantification
- **Complex Function Approximation**: Captures non-linear relationships and interactions

## Applications

AddiVortes is particularly well-suited for:

- **Spatial regression** and geographic data analysis
- **Machine learning** tasks requiring interpretable models
- **Non-parametric regression** where the functional form is unknown
- **Bayesian modeling** with uncertainty quantification
- **Complex surface modeling** and function approximation
- **Alternative to BART** for researchers seeking different ensemble approaches

## Installation

You can install the latest version of AddiVortes from GitHub with:

```R
# install.packages("devtools")
devtools::install_github("johnpaulgosling/AddiVortes", 
                         build_vignettes = TRUE)
```

The Python package is named `addivortes` and targets Python 3.10 or newer.
For local development or installation from a source checkout:

```bash
python -m pip install .
```

## Quick Start

### R

```R
library(AddiVortes)

# Load your data
# X <- your_predictors
# y <- your_response

# Fit the AddiVortes model
# model <- AddiVortes(X, y)

# Make predictions
# predictions <- predict(model, newdata = X_test)
```

### Python

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

## Documentation

- [Getting Started Guide](https://johnpaulgosling.github.io/AddiVortes/articles/introduction.html)
- [Prediction Examples](https://johnpaulgosling.github.io/AddiVortes/articles/prediction.html)
- [Function Reference](https://johnpaulgosling.github.io/AddiVortes/reference/)

## Comparison with BART

While **BART (Bayesian Additive Regression Trees)** uses tree-based partitioning, **AddiVortes** uses Voronoi tessellations, which can provide:

- More natural spatial partitioning
- Flexible geometric boundaries
- Alternative ensemble approach for machine learning
- Enhanced performance on spatial data

## Cite Us

If you use this package in your research, please cite:

```R
citation("AddiVortes")
```

## References

Stone, A. and Gosling, J.P. (2025). AddiVortes: (Bayesian) additive Voronoi tessellations. Journal of Computational and Graphical Statistics.

## Keywords

Bayesian machine learning, BART alternative, Voronoi tessellation, spatial regression, non-parametric regression, ensemble methods, statistical modeling, R package
