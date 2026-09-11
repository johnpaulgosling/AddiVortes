# Create an AddiVortes Object

A constructor for the AddiVortes class.

## Usage

``` r
new_AddiVortes(
  posteriorTess,
  posteriorDim,
  posteriorSigma,
  posteriorPred,
  xCentres,
  xRanges,
  yCentre,
  yRange,
  inSampleRmse,
  metric = "E",
  members = rep(1, length(xCentres)),
  metric_aug = "E",
  member_aug = rep(1, length(xCentres)),
  catEncoding = NULL,
  traceStats = NULL,
  task = "regression",
  classLevels = NULL,
  nLatents = 1L,
  mPerLatent = NA_integer_,
  inSampleAccuracy = NA_real_,
  inSampleBrier = NA_real_
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

- metric:

  The metric used for scaling covariates (default "E" for Euclidean).

- members:

  The membership vector for the covariates

- metric_aug:

  The augmented metric after categorical variables are converted to
  one-hot

- member_aug:

  The membership vector corresponding to metric_aug

- catEncoding:

  Optional list of categorical encoding metadata returned by
  `encodeCategories_internal`, or `NULL` if no categorical covariates
  were present.

- traceStats:

  Optional data frame of per-iteration MCMC trace statistics.

- task:

  The modelling task: `"regression"`, `"binary"` or `"multinomial"`.

- classLevels:

  Character vector of class labels for classification fits, or `NULL`
  for regression.

- nLatents:

  Number of latent probit dimensions. 1 for regression and binary
  classification; \\K-1\\ for \\K\\-class multinomial models.

- mPerLatent:

  Number of tessellations per latent ensemble.

- inSampleAccuracy:

  In-sample classification accuracy, or `NA` for regression.

- inSampleBrier:

  In-sample Brier score for binary classification, or `NA` otherwise.

## Value

An object of class AddiVortes.
