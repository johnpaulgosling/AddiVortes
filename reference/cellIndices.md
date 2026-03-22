# Assign Observations to Tessellation Cells

For a given tessellation, this function identifies which cell (centre)
each observation belongs to based on nearest neighbour classification.

## Usage

``` r
cellIndices(x, tess, dim, metric = "Euclidean")
```

## Arguments

- x:

  A numeric matrix of covariates where each row is an observation.

- tess:

  A numeric matrix representing the tessellation centres, where each row
  is a unique centre.

- dim:

  An integer vector specifying the column indices of `x` to be used for
  calculating distance.

- metric:

  Either "Euclidean" or "Spherical".

## Value

A numeric vector of integers where each element corresponds to a row in
`x` and its value is the row index of the nearest centre in `tess`.

## Details

It finds the closest tessellation centre for each observation (row) in
the covariate matrix, considering only the specified dimensions. This is
achieved using the k-nearest neighbour algorithm where k=1.
