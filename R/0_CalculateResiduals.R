#' @title Calculate Partial Residuals for a Tessellation
#' @description Calculates the vector of partial residuals by excluding the effect
#'   of the jth tessellation. It then aggregates these residuals and observation
#'   counts for both the current and the proposed new tessellation structures.
#'
#' @param y The numeric vector of the outcome or response variable.
#' @param j The index of the tessellation currently being updated.
#' @param SumOfAllTess A numeric vector representing the sum of predictions from
#'   all tessellations.
#' @param Pred A list where each element contains the predicted values for a
#'   tessellation.
#' @param lastTessPred The prediction from the jth tessellation in the previous
#'   iteration, used to correctly update the sum of all tessellation predictions.
#' @param indexes A numeric vector of integers mapping each observation to a
#'   cell in the current tessellation.
#' @param indexesStar A numeric vector of integers mapping each observation to a
#'   cell in the proposed new tessellation.
#' @param num_centres_new The number of centres (cells) in the proposed new
#'   tessellation.
#'
#' @return A list containing seven elements in the following order:
#'   \enumerate{
#'     \item A numeric vector of summed residuals for each cell in the current tessellation (`R_ijOld`).
#'     \item An integer vector of observation counts for each cell in the current tessellation (`n_ijOld`).
#'     \item A numeric vector of summed residuals for each cell in the proposed tessellation (`R_ijNew`).
#'     \item An integer vector of observation counts for each cell in the proposed tessellation (`n_ijNew`).
#'     \item The updated sum of all tessellation predictions (`SumOfAllTess`).
#'     \item The integer vector of cell indexes for the proposed tessellation (`indexesStar`).
#'     \item The integer vector of cell indexes for the current tessellation (`indexes`).
#'   }
#'
#' @keywords internal
#' @noRd
calculateResiduals <- function(y, j, SumOfAllTess, Pred, lastTessPred,
                               indexes, indexesStar, num_centres_new) {
  # This version is no longer passed the temporary TessStar list.
  # Instead, it receives the number of new centres directly.
  
  if (j == 1) {
    CurrentTessPred <- Pred[[j]][indexes]
    SumOfAllTess <- SumOfAllTess - CurrentTessPred
  } else {
    CurrentTessPred <- Pred[[j]][indexes]
    SumOfAllTess <- SumOfAllTess + lastTessPred - CurrentTessPred
  }
  R_j <- y - SumOfAllTess
  
  # --- Old Tessellation Logic (Unchanged) ---
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
  
  # --- New Tessellation Logic (Refactored) ---
  n_ijNew <- tabulate(indexesStar, nbins = num_centres_new)
  if (num_centres_new > 0) {
    R_ijNew_vec <- numeric(num_centres_new)
    res_rowsum_new <- rowsum(R_j, indexesStar)
    R_ijNew_vec[as.numeric(rownames(res_rowsum_new))] <- res_rowsum_new[, 1]
    R_ijNew <- R_ijNew_vec
  } else {
    R_ijNew <- numeric(0)
  }
  
  return(list(
    R_ijOld, n_ijOld,
    R_ijNew, n_ijNew,
    SumOfAllTess, indexesStar, indexes
  ))
}