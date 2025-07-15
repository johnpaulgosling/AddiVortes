#' @title Calculate Residuals
#'
#' @description A function that calculates the n-vector of partial residuals
#' derived from a fitting process that excludes the jth tessellation.
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
#' the residuals of the new tessellation, the number of observations in each
#' cell of the new tessellation, the sum of all tessellations,
#' the indexes of the proposed tessellation, and the indexes of the original
#' tessellation.
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
  
  # Determine the number of levels/groups expected for the old tessellation
  num_levels_old <- length(Pred[[j]])
  n_ijOld <- tabulate(indexes, nbins = num_levels_old)
  
  # --- OPTIMISATION START ---
  # Replaced tapply with a faster rowsum-based approach
  if (num_levels_old > 0) {
    R_ijOld_vec <- numeric(num_levels_old)
    res_rowsum_old <- rowsum(R_j, indexes)
    R_ijOld_vec[as.numeric(rownames(res_rowsum_old))] <- res_rowsum_old[,1]
    R_ijOld <- R_ijOld_vec
  } else {
    R_ijOld <- numeric(0)
  }
  # --- OPTIMISATION END ---
  
  
  # Determine the number of levels/groups expected for the new tessellation
  num_levels_new <- length(TessStar[[j]][, 1])
  n_ijNew <- tabulate(IndexesStar, nbins = num_levels_new)
  
  # --- OPTIMISATION START ---
  # Replaced tapply with a faster rowsum-based approach
  if (num_levels_new > 0) {
    R_ijNew_vec <- numeric(num_levels_new)
    res_rowsum_new <- rowsum(R_j, IndexesStar)
    R_ijNew_vec[as.numeric(rownames(res_rowsum_new))] <- res_rowsum_new[,1]
    R_ijNew <- R_ijNew_vec
  } else {
    R_ijNew <- numeric(0)
  }
  # --- OPTIMISATION END ---
  
  # Return the results as a list
  return(list(
    R_ijOld, n_ijOld,
    R_ijNew, n_ijNew,
    SumOfAllTess, IndexesStar, indexes
  ))
}