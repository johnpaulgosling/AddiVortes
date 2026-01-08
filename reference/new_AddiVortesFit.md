# Create an AddiVortesFit Object

A constructor for the AddiVortesFit class.

## Usage

``` r
new_AddiVortesFit(
  posteriorTess,
  posteriorDim,
  posteriorSigma,
  posteriorPred,
  xCentres,
  xRanges,
  yCentre,
  yRange,
  inSampleRmse
)
```

## Arguments

- posteriorTess:

  A list of the posterior samples of the tessellations.

- posteriorDim:

  A list of the posterior samples of the dimensions.

- posteriorSigma:

  A list of the posterior samples of the error variance.

- posteriorPred:

  A list of the posterior samples of the predictions.

- xCentres:

  The centres of the covariates.

- xRanges:

  The ranges of the covariates.

- yCentre:

  The centre of the output values.

- yRange:

  The range of the output values.

- inSampleRmse:

  The in-sample RMSE.

## Value

An object of class AddiVortesFit.
