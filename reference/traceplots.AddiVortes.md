# Trace Plot Diagnostics for AddiVortes

Displays four MCMC trace plots for a fitted `AddiVortes` object: the
average number of centres per tessellation, the standard deviation of
the number of centres per tessellation, the average number of dimensions
used per tessellation, and the retained-state log-likelihood component
whose differences form the likelihood part of the acceptance ratio.

## Usage

``` r
# S3 method for class 'AddiVortes'
traceplots(x, ask = FALSE, ...)
```

## Arguments

- x:

  An object of class `AddiVortes`, typically the result of a call to
  [`AddiVortes()`](https://johnpaulgosling.github.io/AddiVortes/reference/AddiVortes.md).

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

4.  **Log Likelihood**: Average retained-state log-likelihood component
    at the end of each MCMC iteration.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming 'fit' is a trained AddiVortes object
traceplots(fit)
} # }
```
