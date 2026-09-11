# AddiVortes: (Bayesian) Additive Voronoi Tessellations

## Overview

AddiVortes implements the **Bayesian Additive Voronoi Tessellation** model for machine learning regression, classification and non-parametric statistical modelling. This R package provides a flexible alternative to **BART (Bayesian Additive Regression Trees)**, using Voronoi tessellations instead of trees for spatial partitioning.

## Key Features

- **Machine Learning Regression**: Advanced Bayesian regression modelling for complex datasets
- **Classification**: Binary and multi-category probit models, chosen automatically from the response
- **Alternative to BART**: Uses Voronoi tessellations instead of trees for more flexible spatial modelling
- **Spatial Data Analysis**: Excellent for geographic and spatial datasets
- **Non-parametric Modelling**: No assumptions about functional form
- **Bayesian Framework**: Full posterior inference with uncertainty quantification
- **Complex Function Approximation**: Captures non-linear relationships and interactions

## Applications

AddiVortes is particularly well-suited for:

- **Spatial regression** and geographic data analysis
- **Binary and multi-category classification** with posterior class probabilities
- **Machine learning** tasks requiring interpretable models
- **Non-parametric regression** where the functional form is unknown
- **Bayesian modelling** with uncertainty quantification
- **Complex surface modelling** and function approximation
- **Alternative to BART** for researchers seeking different ensemble approaches

## Installation

You can install the latest version of AddiVortes from GitHub with:

```R
pak::pak("johnpaulgosling/AddiVortes")
```

## Quick Start

```R
library(AddiVortes)

# Load your data
# X <- your_predictors
# y <- your_response

# Fit the AddiVortes model (regression or classification, from y)
# model <- AddiVortes(y, X)

# Make predictions
# predictions <- predict(model, newdata = X_test)
# For classification, type = "response" gives probabilities
# and type = "class" gives labels
```

## Documentation

Vignettes:

- [Machine Learning with AddiVortes](https://johnpaulgosling.github.io/AddiVortes/articles/introduction.html)
- [Bayesian Regression and Prediction](https://johnpaulgosling.github.io/AddiVortes/articles/prediction.html)
- [Modelling Spherical Data with AddiVortes](https://johnpaulgosling.github.io/AddiVortes/articles/spherical.html)
- [Using Categorical Covariates with AddiVortes](https://johnpaulgosling.github.io/AddiVortes/articles/categorical.html)
- [Classification with AddiVortes](https://johnpaulgosling.github.io/AddiVortes/articles/classification.html)

See also the [function reference](https://johnpaulgosling.github.io/AddiVortes/reference/).

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

Stone, A. and Gosling, J.P. (2025). AddiVortes: (Bayesian) additive Voronoi tessellations. *Journal of Computational and Graphical Statistics*, **34**, 859–71. [doi](https://doi.org/10.1080/10618600.2024.2414104)

## Keywords

Bayesian machine learning, BART alternative, Voronoi tessellation, spatial regression, non-parametric regression, ensemble methods, statistical modelling, R package
