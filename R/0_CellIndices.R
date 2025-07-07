#' @title Cell Indices
#'
#' @description A function that gives the row (the centre) of the tessellation
#' that each observation falls within.
#'
#' @param x The covariate matrix.
#' @param tess The tessellation.
#' @param dim The dimensions of the tessellation.
#'
#' @return The row (the centre) of the tessellation that each observation falls 
#' within.
#' @keywords internal
#' @noRd
#'
cellIndices <- function(x, tess, dim) {
  if (length(tess[, 1]) == 1) { # only 1 centre
    CellsForGivenTess <- rep(1, length(x[, 1]))
  } else { # multiple
    CellsForGivenTess <- knnx.index(
      tess,
      matrix(x[, dim],
             ncol = length(dim)
      ), 1
    )
  }
} # Implicit return of CellsForGivenTess
