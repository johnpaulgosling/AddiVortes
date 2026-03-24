# AddiVortes News

## AddiVortes 0.6.0

* Consolidated model fitting into a single C++ call (`addi_vortes_mcmc_cpp`).
  The complete MCMC loop — sigma sampling, tessellation proposals, cell-index
  assignment, residual aggregation, acceptance-probability evaluation, and mu
  sampling — now runs entirely in C++, eliminating the per-iteration R↔C++
  overhead that was present in previous versions.
* The `AddiVortes` class and all associated S3 methods are unchanged.

## AddiVortes 0.5.4

* Added data to package to help with vignettes and testing.

## AddiVortes 0.5.3

* Added vignette "Modelling Spherical Data with AddiVortes" (by Andy Iskauskas
  and John Paul Gosling) demonstrating the use of great-circle distance via
  `metric = "S"`, the coordinate convention (latitude/longitude in radians),
  and a comparison with the Euclidean metric on synthetic globe data.
* Added vignette "Using Categorical Covariates with AddiVortes" explaining
  automatic one-hot encoding (d categories → d-1 binary indicators), the
  `catScaling` parameter, the `catEncoding` metadata field, and handling of
  unseen category levels at prediction time.

## AddiVortes 0.5.2

* Added support for categorical covariates via automatic one-hot encoding.
  Character and factor columns are converted to d-1 binary indicator variables
  (first level as reference), controlled by the new `catScaling` parameter.
* `predict.AddiVortes()` now accepts `data.frame` inputs when the model was
  trained with categorical covariates, applying the stored encoding automatically.
* Unseen categories at prediction time are mapped to the reference level (all zeros).
* Encoding metadata stored in the `catEncoding` field of the fitted model object.
* Added comprehensive tests for categorical covariate support.
* Added second example demonstrating categorical covariates with a larger
  training (n=200) and test set (n=50).

## AddiVortes 0.5.1

* Reinstated C++ speed-ups for new distance metrics
* Restructured testthat tests

## AddiVortes 0.5.0

* Added compatibility for purely spherical input data and mixed Euclidean and Spherical inputs.

## AddiVortes 0.4.12

* Added comprehensive test suite for AddiVortes class methods syntax validation.
* Tests cover input validation and error handling for all S3 methods (print, summary, predict, plot).

## AddiVortes 0.4.11

* Renamed class `AddiVortesFit` to `AddiVortes` to maintain consistency with the package name.
* Updated all S3 methods (print, summary, predict, plot) to use the new class name.
* Updated constructor function from `new_AddiVortesFit` to `new_AddiVortes`.

## AddiVortes 0.4.10

* Fixed variance/standard deviation mismatch in tessellation proposal step.
* Fixed typo in variable names and added more camelCase.

## AddiVortes 0.4.9

* Improved legends in plot.AddiVortesFit().
* Fix typos in package description.

## AddiVortes 0.4.8

* Initial CRAN release of AddiVortes package.
* Fixed installation failure on r-devel-linux-x86_64-fedora-clang:
  - Added missing `<cstring>` header for `memcpy()` function in C++ code
  - Clang compiler requires explicit inclusion of standard library headers

## AddiVortes 0.4.7

* Fixed DESCRIPTION file for CRAN compliance:
  - Removed non-standard 'Keywords' field
  - Removed redundant 'Author' and 'Maintainer' fields (now auto-derived from Authors@R)
* Resolved R CMD check NOTEs for DESCRIPTION meta-information

## AddiVortes 0.4.6

* Linting and formatting improvements.
* Final preparations for CRAN submission.

## AddiVortes 0.4.5

* Fixed bug in tessellation proposal when number of covariates equals number of selected dimensions.
* Add Dimension (AD) modification now properly checks if all covariates are already selected before attempting to add a new one.
* Added comprehensive test suite for small covariate counts (1, 2, and 3 covariates).

## AddiVortes 0.4.4

* Fixed test automation bug stemming from too many cores being assumed available.


## AddiVortes 0.4.3

* Cleaned up tests folder for CRAN submission preparation.
* Removed `TestSuite.R`, `CodeProfiler.R`, and `TestHelper.R` from tests directory.
* Incorporated relevant tests from `TestSuite.R` into testthat framework.
* Added test for thinning parameter functionality.
* Only `testthat.R` remains in tests folder alongside the testthat directory per CRAN requirements.

## AddiVortes 0.4.2

* Added prediction interval support to `predict.AddiVortesFit()` function with new `interval` parameter.
* Fixed bug where `posteriorSigma` was not stored in model objects.
* Prediction intervals now available alongside confidence intervals, similar to `lm` predict function.
* Updated tests to use named access instead of numeric indexing for improved robustness.

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
