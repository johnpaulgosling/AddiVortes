#' @title Assign Observations to Tessellation Cells
#' @description For a given tessellation, this function identifies which cell
#'   (centre) each observation belongs to based on nearest neighbour classification.
#'
#' @details It finds the closest tessellation centre for each observation (row) in
#'   the covariate matrix, considering only the specified dimensions. This is
#'   achieved using the k-nearest neighbour algorithm where k=1.
#'
#' @param x A numeric matrix of covariates where each row is an observation.
#' @param tess A numeric matrix representing the tessellation centres, where each
#'   row is a unique centre.
#' @param dim An integer vector specifying the column indices of `x` to be used
#'   for calculating distance.
#' @param metric Either "Euclidean" or "Spherical".
#' @param members A vector indicating covariate membership.
#'
#' @return A numeric vector of integers where each element corresponds to a row
#'   in `x` and its value is the row index of the nearest centre in `tess`.
#'
#' @export
cellIndices <- function(x, tess, dim, metric = "E", members) {
  n_tess <- nrow(tess)
  n_x <- nrow(x)

  if (n_tess == 1L) { # only 1 centre
    CellsForGivenTess <- rep.int(1L, n_x)
  } else { # multiple
    # Always place the tessellation centre coordinates at their global column
    # positions before computing distances. `tess` stores its columns in the
    # order given by `dim` (the active dimensions), whereas the distance code
    # (knnx_index_cpp) treats column i of `tess` as global covariate i. Guarding
    # this remap on `ncol(tess) != ncol(x)` misses the case where every
    # covariate is active (ncol(tess) == ncol(x)) but `dim` is a permutation of
    # the columns, which silently mismatches coordinates (and, for spherical
    # groups, the order-sensitive distance). Remapping unconditionally is a
    # no-op when `dim` is already the identity order.
    new_tess <- matrix(0, nrow = n_tess, ncol = ncol(x))
    new_tess[, dim] <- tess
    tess <- new_tess
    CellsForGivenTess <- knnx_index(
      tess,
      x,
      dim, metric, members
    )
  }
  return(CellsForGivenTess)
} # Implicit return of CellsForGivenTess
