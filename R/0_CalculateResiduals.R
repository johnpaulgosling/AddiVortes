#' @title Calculate Residuals
#'
#' @description A function that calculates the n-vector of partial residuals derived from a
#' fitting process that excludes the jth tessellation.
#'
#' @param y The outcome variable.
#' @param x The covariate matrix.
#' @param j The index of the tessellation.
#' @param SumOfAllTess The sum of all tessellations.
#' @param Tess The tessellation.
#' @param Dim The dimensions of the tessellation.
#' @param Pred The prediction set.
#' @param TessStar The proposed tessellation.
#' @param DimStar The proposed dimensions.
#' @param LastTessPred The last tessellation prediction.
#'
#' @return A list containing the residuals of the old tessellation,
#' the number of observations in each cell of the old tessellation,
#' the residuals of the new tessellation, the number of observations in each cell of the
#' new tessellation, the sum of all tessellations,
#' the indexes of the proposed tessellation, and the indexes of the original tessellation.
#'
#' @keywords internal
#' @noRd
#'
calculateResiduals <- function(y, x, j, SumOfAllTess, Tess, Dim,
                                Pred, TessStar, DimStar, LastTessPred) {
  if (j == 1) {
    indexes <- cellIndices(x, Tess[[j]], Dim[[j]])
    CurrentTessPred <- Pred[[j]][indexes]
    SumOfAllTess <- SumOfAllTess - CurrentTessPred
  } else {
    indexes <- cellIndices(x, Tess[[j]], Dim[[j]])
    CurrentTessPred <- Pred[[j]][indexes]
    SumOfAllTess <- SumOfAllTess + LastTessPred - CurrentTessPred
  }

  IndexesStar <- cellIndices(x, TessStar[[j]], DimStar[[j]])
  R_j <- y - SumOfAllTess

  # Initializing Sizes
  # Determine the number of levels/groups expected
  num_levels_old <- length(Pred[[j]])

  if (num_levels_old > 0) {
    # Use tabulate for efficient counting of integer occurrences (1 to num_levels_old)
    n_ijOld <- tabulate(indexes,
                        nbins = num_levels_old)
    # Use tapply to sum R_j based on the groups defined by 'indexes'
    R_ijOld <- tapply(R_j, factor(indexes,
                                  levels = 1:num_levels_old),
                      sum, na.rm = TRUE)
    # Ensure the output is a simple vector (tapply can return an array)
    R_ijOld <- as.vector(R_ijOld)
  } else {
    # Handle empty case
    n_ijOld <- numeric(0)
    R_ijOld <- numeric(0)
  }

  # Determine the number of levels/groups expected (using the length from the original loop)
  num_levels_new <- length(TessStar[[j]][, 1])

  if (num_levels_new > 0) {
    # Use tabulate for counts
    n_ijNew <- tabulate(IndexesStar,
                        nbins = num_levels_new)

    # Use tapply for sums
    R_ijNew <- tapply(R_j, factor(IndexesStar,
                                  levels = 1:num_levels_new),
                      sum, na.rm = TRUE)

    # Ensure vector output
    R_ijNew <- as.vector(R_ijNew)
  } else {
    # Handle empty case
    n_ijNew <- numeric(0)
    R_ijNew <- numeric(0)
  }

  # Return the results as a list
  return(list(R_ijOld, n_ijOld,
              R_ijNew, n_ijNew,
              SumOfAllTess, IndexesStar, indexes))
}
