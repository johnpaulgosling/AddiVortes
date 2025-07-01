#' @title Cell Indices
#'
#' @description A function that gives the row (the centre) of the tessellation
#' that each observation falls within.
#'
#' @param x The covariate matrix.
#' @param Tess The tessellation.
#' @param Dim The dimensions of the tessellation.
#'
#' @return The row (the centre) of the tessellation that each observation falls within.
#' @keywords internal
#' @noRd
#'
cellIndices <- function(x, Tess, Dim) {
  if (length(Tess[, 1]) == 1) { # only 1 centre
    CellsForGivenTess <- rep(1, length(x[, 1]))
  } else { # multiple
    CellsForGivenTess <- knnx.index(Tess,
                                    matrix(x[, Dim],
                                           ncol = length(Dim)), 1)
  }
  return(CellsForGivenTess)
}
