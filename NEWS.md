# AddiVortes News

## AddiVortes 0.4.1

* Added parallel processing to predict function.
* Removed unnecessary square-root function from KD-tree search.

## AddiVortes 0.3.3

* Added warning for situation where covariates exceed observations.

## AddiVortes 0.3.2

* Improved progress bar output to give more information during MCMC sampling and prediction processes.
* Removed outputted progress bars from vignette examples to reduce clutter.

## AddiVortes 0.3.1

* Enhanced SEO and discoverability with comprehensive keyword optimization across package documentation
* Added formal Keywords field in DESCRIPTION for better search indexing
* Updated README.md with clear positioning as BART alternative
* Optimized vignette titles for machine learning and Bayesian regression search terms
* Enhanced package documentation with machine learning focus

## AddiVortes 0.3.0

* Implemented nearest neighbour search in C++.
* Removed dependency on FNN package.

## AddiVortes 0.2.5

### Package Features

* Added progress bars to MCMC sampling functions and subsequent prediction functions for better user feedback during long computations.

## AddiVortes 0.2.4

### CRAN Submission Preparation

* Removed compiled object files from package source
* Added examples to main `AddiVortes()` function
* Updated DESCRIPTION file for CRAN compliance:
  - Added Author and Maintainer fields  
  - Added R version dependency (>= 3.5.0)
  - Updated Date field
  - Fixed grammatical error in Description
* Updated .Rbuildignore with standard exclusions
* Enhanced documentation with examples for key functions

### Package Features

* Implements Bayesian Additive Voronoi Tessellation models
* Non-parametric regression with tessellation-based approach  
* Posterior sampling via backfitting algorithm
* Prediction and visualization methods for fitted models
* Comprehensive test suite and vignettes