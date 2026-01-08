# Print Method for AddiVortesFit

Prints a summary of a fitted `AddiVortesFit` object, providing
information about the model structure, dimensions, and fit quality
similar to the output of a linear model summary.

## Usage

``` r
# S3 method for class 'AddiVortesFit'
print(x, ...)
```

## Arguments

- x:

  An object of class `AddiVortesFit`, typically the result of a call to
  [`AddiVortes()`](https://johnpaulgosling.github.io/AddiVortes/reference/AddiVortes.md).

- ...:

  Further arguments passed to or from other methods (currently unused).

## Value

The function is called for its side effect of printing model information
and returns the input object `x` invisibly.

## Details

The print method displays:

- The model formula representation

- Number of covariates and posterior samples

- Number of tessellations used

- In-sample RMSE

- Covariate scaling information
