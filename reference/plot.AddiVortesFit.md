# Plot Method for AddiVortesFit

Generates comprehensive diagnostic plots for a fitted `AddiVortesFit`
object. This function creates multiple diagnostic plots including
residuals, MCMC traces for sigma, and tessellation complexity over
iterations.

## Usage

``` r
# S3 method for class 'AddiVortesFit'
plot(
  x,
  x_train,
  y_train,
  sigma_trace = NULL,
  which = c(1, 2, 3),
  ask = FALSE,
  ...
)
```

## Arguments

- x:

  An object of class `AddiVortesFit`, typically the result of a call to
  [`AddiVortes()`](https://johnpaulgosling.github.io/AddiVortes/reference/AddiVortes.md).

- x_train:

  A matrix of the original training covariates.

- y_train:

  A numeric vector of the original training true outcomes.

- sigma_trace:

  An optional numeric vector of sigma values from MCMC samples. If not
  provided, the method will attempt to extract it from the model object.

- which:

  A numeric vector specifying which plots to generate: 1 = Residuals
  plot, 2 = Sigma trace, 3 = Tessellation complexity trace, 4 =
  Predicted vs Observed. Default is c(1, 2, 3).

- ask:

  Logical; if TRUE, the user is asked to press Enter before each plot.

- ...:

  Additional arguments passed to plotting functions.

## Value

This function is called for its side effect of creating plots and
returns `NULL` invisibly.

## Details

The function generates up to four diagnostic plots:

1.  **Residuals Plot**: Residuals vs fitted values with smoothed trend
    line

2.  **Sigma Trace**: MCMC trace plot for the error variance parameter

3.  **Tessellation Complexity**: Trace of average tessellation size over
    iterations

4.  **Predicted vs Observed**: Scatter plot with confidence intervals

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming 'fit' is a trained AddiVortesFit object
plot(fit, x_train = x_train_data, y_train = y_train_data)

# Show only specific plots
plot(fit, x_train = x_train_data, y_train = y_train_data, which = c(1, 3))

# With custom sigma trace
plot(fit, x_train = x_train_data, y_train = y_train_data, 
     sigma_trace = my_sigma_samples)
} # }
```
