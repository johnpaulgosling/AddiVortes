# AddiVortes 0.2.4

## CRAN Submission Preparation

* Removed compiled object files from package source
* Added examples to main `AddiVortes()` function
* Updated DESCRIPTION file for CRAN compliance:
  - Added Author and Maintainer fields  
  - Added R version dependency (>= 3.5.0)
  - Updated Date field
  - Fixed grammatical error in Description
* Updated .Rbuildignore with standard exclusions
* Enhanced documentation with examples for key functions

## Package Features

* Implements Bayesian Additive Voronoi Tessellation models
* Non-parametric regression with tessellation-based approach  
* Posterior sampling via backfitting algorithm
* Prediction and visualization methods for fitted models
* Comprehensive test suite and vignettes