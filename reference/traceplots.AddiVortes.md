# Trace Plot Diagnostics for AddiVortes

Displays four MCMC trace plots for a fitted `AddiVortes` object: the
average number of centres per tessellation, the standard deviation of
the number of centres per tessellation, the average number of dimensions
used per tessellation, and the error standard deviation.

## Usage

``` r
# S3 method for class 'AddiVortes'
traceplots(x, sigma_trace = NULL, ask = FALSE, ...)
```

## Arguments

- x:

  An object of class `AddiVortes`, typically the result of a call to
  [`AddiVortes()`](https://johnpaulgosling.github.io/AddiVortes/reference/AddiVortes.md).

- sigma_trace:

  An optional numeric vector of error standard deviation values from
  MCMC samples. If not provided, the method uses
  `sqrt(x$posteriorSigma)`.

- ask:

  Logical; if TRUE, the user is asked to press Enter before each plot.

- ...:

  Additional arguments passed to plotting functions.

## Value

This function is called for its side effect of creating plots and
returns `NULL` invisibly.

## Details

The four trace plots are:

1.  **Average Centres**: Average number of centres per tessellation.

2.  **Centre Count Standard Deviation**: Standard deviation of the
    number of centres per tessellation.

3.  **Average Dimensions**: Average number of active dimensions used per
    tessellation.

4.  **Error Standard Deviation**: MCMC trace for the error standard
    deviation.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming 'fit' is a trained AddiVortes object
traceplots(fit)

# With a custom error standard deviation trace
traceplots(fit, sigma_trace = my_sigma_samples)
} # }
```
