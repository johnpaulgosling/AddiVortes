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
calculateResiduals <- function(y, j, SumOfAllTess, Pred, TessStar, LastTessPred,
                               indexes, indexesStar) {
  # This function now receives 'indexes' and 'indexesStar' and does no
  # k-NN searching itself.
  
  if (j == 1) {
    CurrentTessPred <- Pred[[j]][indexes]
    SumOfAllTess <- SumOfAllTess - CurrentTessPred
  } else {
    CurrentTessPred <- Pred[[j]][indexes]
    SumOfAllTess <- SumOfAllTess + LastTessPred - CurrentTessPred
  }
  
  R_j <- y - SumOfAllTess
  
  # Determine the number of levels/groups expected for the old tessellation
  num_levels_old <- length(Pred[[j]])
  n_ijOld <- tabulate(indexes, nbins = num_levels_old)
  
  if (num_levels_old > 0) {
    R_ijOld_vec <- numeric(num_levels_old)
    res_rowsum_old <- rowsum(R_j, indexes)
    R_ijOld_vec[as.numeric(rownames(res_rowsum_old))] <- res_rowsum_old[, 1]
    R_ijOld <- R_ijOld_vec
  } else {
    R_ijOld <- numeric(0)
  }
  
  # Determine the number of levels/groups expected for the new tessellation
  num_levels_new <- length(TessStar[[j]][, 1])
  n_ijNew <- tabulate(indexesStar, nbins = num_levels_new)
  
  if (num_levels_new > 0) {
    R_ijNew_vec <- numeric(num_levels_new)
    res_rowsum_new <- rowsum(R_j, indexesStar)
    R_ijNew_vec[as.numeric(rownames(res_rowsum_new))] <- res_rowsum_new[, 1]
    R_ijNew <- R_ijNew_vec
  } else {
    R_ijNew <- numeric(0)
  }
  
  # Return the results as a list
  return(list(
    R_ijOld, n_ijOld,
    R_ijNew, n_ijNew,
    SumOfAllTess, indexesStar, indexes
  ))
}