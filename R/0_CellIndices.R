#' @title Assign Observations to Tessellation Cells (Brute Force)
#' @description For a given tessellation, this function identifies which cell
#'   (centre) each observation belongs to based on nearest neighbour classification
#'   using a brute force approach.
#'
#' @details It finds the closest tessellation centre for each observation (row) in
#'   the covariate matrix, considering only the specified dimensions. This is
#'   achieved using a brute force distance calculation which is faster than
#'   KD-tree for typical tessellation sizes.
#'
#' @param obs A numeric matrix where each row is an observation.
#' @param centers A numeric matrix where each row is a center.
#'
#' @return A numeric vector of integers where each element corresponds to a row
#'   in `obs` and its value is the row index of the nearest centre in `centers`.
#'
#' @keywords internal
#' @noRd
assign_bruteforce <- function(obs, centers) {
  obs_sq <- rowSums(obs^2)
  cn_sq  <- rowSums(centers^2)
  D <- outer(obs_sq, cn_sq, "+") - 2*(obs %*% t(centers))
  max.col(-D)
}

#' @title Assign Observations to Tessellation Cells
#' @description For a given tessellation, this function identifies which cell
#'   (centre) each observation belongs to based on nearest neighbour classification.
#'
#' @details It finds the closest tessellation centre for each observation (row) in
#'   the covariate matrix, considering only the specified dimensions. This is
#'   achieved using a brute force distance calculation.
#'
#' @param x A numeric matrix of covariates where each row is an observation.
#' @param tess A numeric matrix representing the tessellation centres, where each
#'   row is a unique centre.
#' @param dim An integer vector specifying the column indices of `x` to be used
#'   for calculating distance.
#'
#' @return A numeric vector of integers where each element corresponds to a row
#'   in `x` and its value is the row index of the nearest centre in `tess`.
#'
#' @keywords internal
#' @noRd
cellIndices <- function(x, tess, dim) {
  if (length(tess[, 1]) == 1) { # only 1 centre
    CellsForGivenTess <- rep(1, length(x[, 1]))
  } else { # multiple
    # Extract the relevant dimensions from observations
    obs_subset <- matrix(x[, dim], ncol = length(dim))
    # Use brute force assignment
    CellsForGivenTess <- assign_bruteforce(obs_subset, tess)
  }
} # Implicit return of CellsForGivenTess
